#include "parallhull.h"

#include <math.h>
#include <immintrin.h>
#ifdef NON_MPI_MODE
    #include <time.h>
    #include <stdio.h>
    #include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#else
    #include <mpi.h>
#endif



#define USE_MANUAL_AVX_PIPELINE_OPTIMIZATION // without this takes removeCoveredPoints takes more than triple the time

#define HULL_ALLOC_ELEMS 1000


static void getExtremeCoordsPts(Data *pts, size_t ptIndices[4]);
static void extremeCoordsInit(Data *hull, Data *uncoveredPts, size_t ptIndices[4]);
static void removeCoveredPoints(Data *hull, Data *uncoveredPts, char *uncoveredCache, ProcThreadIDCombo *id);
static void findFarthestPts(Data *hull, Data *uncoveredPts, size_t *maxDistPtIndices);
static void addPtsToHull(Data *hull, Data *uncoveredPts, size_t **maxDistPtIndicesPtr, size_t **offsetCounterPtr, size_t *allocatedElemsCount, ProcThreadIDCombo *id);

Data quickhull (Data *d, ProcThreadIDCombo *id)
{
    int iterCount = 0;

    #ifdef NON_MPI_MODE
        struct timespec timeStruct;
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        double startTime = cvtTimespec2Double(timeStruct);
        double iterTime = startTime;
        double previousIterTime = 0;
    #else
        double startTime = MPI_Wtime();
        double iterTime = startTime;
        double previousIterTime = 0;
    #endif

    Data hull, uncoveredPts;
    uncoveredPts = *d;

    hull.n = 0;
    size_t allocatedElemsCount = HULL_ALLOC_ELEMS < uncoveredPts.n ? HULL_ALLOC_ELEMS+1 : uncoveredPts.n+1;
    hull.X = malloc(allocatedElemsCount * sizeof(float) + MALLOC_PADDING);
    if (hull.X == NULL)
        throwError("p[%2d] t[%3d] extremeCoordsInit: Failed to allocate initial memory for the hull coordinates", id->p, id->t);
    hull.Y = malloc(allocatedElemsCount * sizeof(float) + MALLOC_PADDING);
    if (hull.Y == NULL)
        throwError("p[%2d] t[%3d] extremeCoordsInit: Failed to allocate initial memory for the hull coordinates", id->p, id->t);
    size_t *offsetCounter = malloc(allocatedElemsCount * 2 * sizeof(size_t) + MALLOC_PADDING*2);
    if (offsetCounter == NULL)
        throwError("p[%2d] t[%3d] quickhull: Failed to allocate memory", id->p, id->t);
    size_t *maxDistPtIndices = &offsetCounter[allocatedElemsCount];

    // init
    size_t ptIndices[4];
    getExtremeCoordsPts(d, ptIndices);
    extremeCoordsInit(&hull, &uncoveredPts, ptIndices);
    
    while (uncoveredPts.n > 0)
    {
        if (allocatedElemsCount * sizeof(size_t) * 2 >= uncoveredPts.n)
            removeCoveredPoints(&hull, &uncoveredPts, (char*)offsetCounter, id);
        else
            removeCoveredPoints(&hull, &uncoveredPts, NULL, id);
        
        if (uncoveredPts.n == 0) break;

        iterCount++;   
        #ifdef NON_MPI_MODE 
            previousIterTime = iterTime;
            clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
            iterTime = cvtTimespec2Double(timeStruct);
        #else
            iterTime = MPI_Wtime();
        #endif

        LOG(LOG_LVL_TRACE, "p[%2d] t[%3d] quickhull: Iteration %5d lasted %.3es, %.3es from the begining. nUncovered=%.3e, hullSize=%ld", id->p, id->t, iterCount, iterTime-previousIterTime, iterTime-startTime, (float)uncoveredPts.n, hull.n);
        
        #ifdef QUICKHULL_STEP_DEBUG
            // show partial hull with gnuplot at each iteration and wait for user input to resume
            char plotTitle[200];
            sprintf(plotTitle, "p[%2d] t[%3d] Partial Hull: size=%lu, uncovered=%lu", id->p, id->t, hull.n, uncoveredPts.n);
            plotData(d, &hull, uncoveredPts.n, plotTitle);
            getchar();
        #endif

        findFarthestPts(&hull, &uncoveredPts, maxDistPtIndices);
        #ifdef DEBUG
            size_t oldNUncovered = uncoveredPts.n;
        #endif

        addPtsToHull(&hull, &uncoveredPts, &maxDistPtIndices, &offsetCounter, &allocatedElemsCount, id);

        #ifdef DEBUG
            if (hullConvexityCheck(&hull, id))
                throwError("p[%2d] t[%3d] quickhull: Hull is not convex. n=%ld, hullSize=%ld, nUncovered=%ld, oldNUncovered=%ld", id->p, id->t, d->n, hull.n, uncoveredPts.n, oldNUncovered);
        #endif
    }

    #ifdef QUICKHULL_STEP_DEBUG
        // show partial hull with gnuplot at each iteration and wait for user input to resume
        char plotTitle[200];
        sprintf(plotTitle, "p[%2d] t[%3d] quickhull: size=%lu, uncovered=%lu", id->p, id->t, hull.n, uncoveredPts.n);
        plotData(d, &hull, uncoveredPts.n, plotTitle);
        getchar();
    #endif

    #ifdef DEBUG
        LOG(LOG_LVL_DEBUG, "p[%2d] t[%3d] quickhull: DEBUG macro is defined! Now checking whether all points are actually inside the hull", id->p, id->t);
        if (finalCoverageCheck(&hull, d, id) != 0)
            throwError("p[%2d] t[%3d] quickhull: There are still %ld points that are not inside the hull", id->p, id->t, uncoveredPts.n);
    #endif

    free(offsetCounter);
    hull.X = realloc(hull.X, (hull.n+1) * sizeof(float) + MALLOC_PADDING);
    hull.Y = realloc(hull.Y, (hull.n+1) * sizeof(float) + MALLOC_PADDING);

    return hull;
}

static void getExtremeCoordsPts(Data *pts, size_t ptIndices[4])
{
    for (int i = 0; i < 4; i++)
        ptIndices[i] = 0;
    
    for (size_t i = 1; i < pts->n; i++)
    {
        if  ((pts->X[i] > pts->X[ptIndices[1]]) || 
            ((pts->X[i] == pts->X[ptIndices[1]]) && (pts->Y[i] > pts->Y[ptIndices[1]]))) //xMax rightmust upward if multiple choices
            ptIndices[1] = i;
        else if ((pts->X[i] < pts->X[ptIndices[3]]) || 
                ((pts->X[i] == pts->X[ptIndices[3]]) && (pts->Y[i] < pts->Y[ptIndices[3]]))) //xMin leftmost loward if multiple choices
            ptIndices[3] = i;

        if ((pts->Y[i] > pts->Y[ptIndices[2]]) || 
            ((pts->Y[i] == pts->Y[ptIndices[2]]) && (pts->X[i] < pts->X[ptIndices[2]]))) //yMax upmost leftward if multiple choices
            ptIndices[2] = i;
        else if ((pts->Y[i] < pts->Y[ptIndices[0]]) || 
                ((pts->Y[i] == pts->Y[ptIndices[0]]) && (pts->X[i] > pts->X[ptIndices[0]]))) //yMin lowmost rightward if multiple choices
            ptIndices[0] = i;
    }
}

static void extremeCoordsInit(Data *hull, Data *uncoveredPts, size_t ptIndices[4])
{    
    // check extreme coords for duplicates
    for (int i = 0; i < 3; i++)
    {
        if (ptIndices[i] == -1) break;

        for (int j = i+1; j < 4; j++)
        {
            if (ptIndices[i] == ptIndices[j])
            {
                ptIndices[j] = -1;
                for (int k = j; k < 3; k++)
                    swapElems(ptIndices[k], ptIndices[k+1])
            }
        }
    }
    
    for (int i = 0; (i < 4) && (ptIndices[i] != -1); i++)
    {
        hull->X[i] = uncoveredPts->X[ptIndices[i]];
        hull->Y[i] = uncoveredPts->Y[ptIndices[i]];
        hull->n++;
    }
    hull->X[hull->n] = hull->X[0];
    hull->Y[hull->n] = hull->Y[0];

    // remove points from uncovered set
    // sort first
    for (int i = 0; (i < 4) && (ptIndices[i] != -1); i++)
        for (int j =i+1; (j < 4) && (ptIndices[j] != -1); j++)
            if (ptIndices[j] > ptIndices[i])
                swapElems(ptIndices[i], ptIndices[j])
    
    // now remove the elements
    for (int i = 0; (i < 4) && (ptIndices[i] != -1); i++)
    {
        uncoveredPts->n--;
        swapElems(uncoveredPts->X[ptIndices[i]], uncoveredPts->X[uncoveredPts->n]);
        swapElems(uncoveredPts->Y[ptIndices[i]], uncoveredPts->Y[uncoveredPts->n]);
    }
}

static void removeCoveredPoints(Data *hull, Data *uncoveredPts, char *uncoveredCache, ProcThreadIDCombo *id)
{
    bool allocatedMem = false;
    if (uncoveredCache == NULL)
    {
        allocatedMem = true;
        uncoveredCache = malloc(uncoveredPts->n * sizeof(char) + MALLOC_PADDING);
        if (uncoveredCache == NULL)
            throwError("p[%2d] t[%3d] removeCoveredPoints: Failed to allocate memory for the uncoveredCache", id->p, id->t);
    }

    // zero cache
    for (size_t i = 0; i < uncoveredPts->n; i++)
        uncoveredCache[i] = 0;
    
    for (size_t h = 0; h < hull->n; h++)
    {
        __m256i shuffleMask = _mm256_setr_epi64x(0xFFFFFFFFFFFF0800, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFF0800FFFF , 0xFFFFFFFFFFFFFFFF);
        __m256d a, b, c;
        {
            __m128 hX0s = _mm_broadcast_ss(&hull->X[h]);
            __m128 hX1s = _mm_broadcast_ss(&hull->X[h+1]);
            __m128 hY0s = _mm_broadcast_ss(&hull->Y[h]);
            __m128 hY1s = _mm_broadcast_ss(&hull->Y[h+1]);
            __m256d hX0 = _mm256_cvtps_pd(hX0s);
            __m256d hX1 = _mm256_cvtps_pd(hX1s);
            __m256d hY0 = _mm256_cvtps_pd(hY0s);
            __m256d hY1 = _mm256_cvtps_pd(hY1s);
            a = _mm256_sub_pd(hX1, hX0);
            b = _mm256_sub_pd(hY0, hY1);
            c = _mm256_fmsub_pd(hX0, hY1, _mm256_mul_pd(hX1, hY0));
        }
        size_t l = 12;
        #ifdef USE_MANUAL_AVX_PIPELINE_OPTIMIZATION
        for (size_t i=0, j=4, k=8; l < uncoveredPts->n; i+=16, j+=16, k+=16, l+=16) // loop using inline
        {
            __m256d iDist, jDist, kDist, lDist;
            __m256d ptXi = _mm256_cvtps_pd(_mm_loadu_ps(&uncoveredPts->X[i]));
            __m256d ptYi = _mm256_cvtps_pd(_mm_loadu_ps(&uncoveredPts->Y[i]));
            __m256d ptXj = _mm256_cvtps_pd(_mm_loadu_ps(&uncoveredPts->X[j]));
            __m256d ptYj = _mm256_cvtps_pd(_mm_loadu_ps(&uncoveredPts->Y[j]));
            __m256d ptXk = _mm256_cvtps_pd(_mm_loadu_ps(&uncoveredPts->X[k]));
            __m256d ptYk = _mm256_cvtps_pd(_mm_loadu_ps(&uncoveredPts->Y[k]));
            __m256d ptXl = _mm256_cvtps_pd(_mm_loadu_ps(&uncoveredPts->X[l]));
            __m256d ptYl = _mm256_cvtps_pd(_mm_loadu_ps(&uncoveredPts->Y[l]));
        
            __m128i u128 = _mm_loadu_si128((__m128i_u*)&uncoveredCache[i]);

            iDist = _mm256_fmadd_pd(b, ptXi, c);
            jDist = _mm256_fmadd_pd(b, ptXj, c);
            kDist = _mm256_fmadd_pd(b, ptXk, c);
            lDist = _mm256_fmadd_pd(b, ptXl, c);

            iDist = _mm256_fmadd_pd(a, ptYi, iDist);
            jDist = _mm256_fmadd_pd(a, ptYj, jDist);
            kDist = _mm256_fmadd_pd(a, ptYk, kDist);
            lDist = _mm256_fmadd_pd(a, ptYl, lDist);

            __m256i iMask = _mm256_castpd_si256(_mm256_cmp_pd(iDist, _mm256_setzero_pd(), _CMP_LT_OQ));
            __m256i jMask = _mm256_castpd_si256(_mm256_cmp_pd(jDist, _mm256_setzero_pd(), _CMP_LT_OQ));
            __m256i kMask = _mm256_castpd_si256(_mm256_cmp_pd(kDist, _mm256_setzero_pd(), _CMP_LT_OQ));
            __m256i lMask =_mm256_castpd_si256( _mm256_cmp_pd(lDist, _mm256_setzero_pd(), _CMP_LT_OQ));

            iMask = _mm256_shuffle_epi8(iMask, shuffleMask);
            jMask = _mm256_shuffle_epi8(jMask, shuffleMask);
            kMask = _mm256_shuffle_epi8(kMask, shuffleMask);
            lMask = _mm256_shuffle_epi8(lMask, shuffleMask);

            __m128i mask;
            {
            __m128i iMaskLow = _mm256_castsi256_si128(iMask);
            __m128i jMaskLow = _mm256_castsi256_si128(jMask);
            __m128i kMaskLow = _mm256_castsi256_si128(kMask);
            __m128i lMaskLow = _mm256_castsi256_si128(lMask);

            __m128i iMaskHigh = _mm256_extractf128_si256(iMask, 1);
            __m128i jMaskHigh = _mm256_extractf128_si256(jMask, 1);
            __m128i kMaskHigh = _mm256_extractf128_si256(kMask, 1);
            __m128i lMaskHigh = _mm256_extractf128_si256(lMask, 1);

            iMaskLow = _mm_or_si128(iMaskLow, iMaskHigh);
            jMaskLow = _mm_or_si128(jMaskLow, jMaskHigh);
            kMaskLow = _mm_or_si128(kMaskLow, kMaskHigh);
            lMaskLow = _mm_or_si128(lMaskLow, lMaskHigh);

            jMaskLow = _mm_shuffle_epi32(jMaskLow, 0xF3);
            kMaskLow = _mm_shuffle_epi32(kMaskLow, 0xCF);
            lMaskLow = _mm_shuffle_epi32(lMaskLow, 0x3F);

            iMaskLow = _mm_or_si128(iMaskLow, jMaskLow);
            kMaskLow = _mm_or_si128(kMaskLow, lMaskLow);
            mask = _mm_or_si128(iMaskLow, kMaskLow);
            }

            mask = _mm_or_si128(mask, u128);

            _mm_storeu_si128((__m128i_u*)&uncoveredCache[i], mask);
        }// */
        #endif
        for (size_t i=l-12; i < uncoveredPts->n; i+=4)
        {
            __m256d iDist;
            __m256d ptX = _mm256_cvtps_pd(_mm_loadu_ps(&uncoveredPts->X[i]));
            __m256d ptY = _mm256_cvtps_pd(_mm_loadu_ps(&uncoveredPts->Y[i]));

            __m128i u128 = _mm_loadu_si128((__m128i_u*)&uncoveredCache[i]);

            iDist = _mm256_fmadd_pd(b, ptX, c);
            iDist = _mm256_fmadd_pd(a, ptY, iDist);

            __m256i iMask = _mm256_castpd_si256(_mm256_cmp_pd(iDist, _mm256_setzero_pd(), _CMP_LT_OQ));

            iMask = _mm256_shuffle_epi8(iMask, shuffleMask);

            __m128i mask;
            {
            __m128i iMaskLow = _mm256_castsi256_si128(iMask);
            __m128i iMaskHigh = _mm256_extractf128_si256(iMask, 1);
            iMaskLow = _mm_or_si128(iMaskLow, iMaskHigh);
            mask = iMaskLow;
            }

            mask = _mm_or_si128(mask, u128);

            _mm_storeu_si128((__m128i_u*)&uncoveredCache[i], mask);
        }// */
    }


    size_t i = 0;
    size_t j = uncoveredPts->n - 1;
    while (i <= j)
    {
        if (uncoveredCache[i] == 0)
        {
            while ((j > i) && (uncoveredCache[j] == 0))
                j--;

            if (i == j) break;

            swapElems(uncoveredCache[i], uncoveredCache[j])
            swapElems(uncoveredPts->X[i], uncoveredPts->X[j])
            swapElems(uncoveredPts->Y[i], uncoveredPts->Y[j])
            
            j--;
        }
        i++;
    }

    if (allocatedMem)
        free(uncoveredCache);

    uncoveredPts->n = i;
}

static void findFarthestPts(Data *hull, Data *uncoveredPts, size_t *maxDistPtIndices)
{
    for (size_t k = 0; k < hull->n; k+=4)
    {
        __m256d a, b, c;
        {
            __m256d hX0 = _mm256_cvtps_pd(_mm_loadu_ps(&hull->X[k]));
            __m256d hX1 = _mm256_cvtps_pd(_mm_loadu_ps(&hull->X[k+1]));
            __m256d hY0 = _mm256_cvtps_pd(_mm_loadu_ps(&hull->Y[k]));
            __m256d hY1 = _mm256_cvtps_pd(_mm_loadu_ps(&hull->Y[k+1]));
            a = _mm256_sub_pd(hX1, hX0);
            b = _mm256_sub_pd(hY0, hY1);
            c = _mm256_fmsub_pd(hX0, hY1, _mm256_mul_pd(hX1, hY0));
        }

        __m256d maxDist = _mm256_setzero_pd();
        __m256i maxDistPtsID = _mm256_set1_epi64x(-1);
        for (size_t i = 0; i < uncoveredPts->n; i++)
        {
            __m256d x = _mm256_cvtps_pd(_mm_broadcast_ss(&uncoveredPts->X[i]));
            __m256d y = _mm256_cvtps_pd(_mm_broadcast_ss(&uncoveredPts->Y[i]));

            __m256i ptsID = _mm256_set1_epi64x(i);

            __m256d dist = _mm256_fmadd_pd(b, x, c);
            dist = _mm256_fmadd_pd(a, y, dist);

            __m256d cmpMask = _mm256_cmp_pd(dist, maxDist, _CMP_LT_OQ);
            maxDist = _mm256_min_pd(maxDist, dist);
            maxDistPtsID = _mm256_blendv_epi8(maxDistPtsID, ptsID, _mm256_castpd_si256(cmpMask));
        }
        
        _mm256_storeu_si256((__m256i_u*)&maxDistPtIndices[k], maxDistPtsID);
    }
}

static void addPtsToHull(Data *hull, Data *uncoveredPts, size_t **maxDistPtIndicesPtr, size_t **offsetCounterPtr, size_t *allocatedElemsCount, ProcThreadIDCombo *id)
{
    size_t *offsetCounter = *offsetCounterPtr;
    size_t *maxDistPtIndices = *maxDistPtIndicesPtr;

    // compute offset to apply in hull in order to minimize the number of swaps
    offsetCounter[0] = 0; // first element always stays there
    for (size_t i = 1; i <= hull->n; i++)
    {
        if (maxDistPtIndices[i-1] != -1)
            offsetCounter[i] = offsetCounter[i-1] + 1;
        else
            offsetCounter[i] = offsetCounter[i-1];
    }
    
    // realloc memory if necessary
    bool reallocMemory = false;
    if (hull->n + offsetCounter[hull->n] > *allocatedElemsCount)
    {
        reallocMemory = true;
        *allocatedElemsCount *= 4;
        hull->X = realloc(hull->X, *allocatedElemsCount * sizeof(float) + MALLOC_PADDING);
        if (hull->X == NULL)
            throwError("p[%2d] t[%3d] addPtToHull: Failed to reallocate hull.X. uncoveredPts=%ld", id->p, id->t, uncoveredPts->n);
        hull->Y = realloc(hull->Y, *allocatedElemsCount * sizeof(float) + MALLOC_PADDING);
        if (hull->Y == NULL)
            throwError("p[%2d] t[%3d] addPtToHull: Failed to reallocate hull.Y. uncoveredPts=%ld", id->p, id->t, uncoveredPts->n);
    }

    // make space in hull.X and hull.Y to fit new points
    for (size_t i = hull->n; offsetCounter[i] != 0; i--)
    {
        size_t iOffset = i + offsetCounter[i];
        hull->X[iOffset] = hull->X[i];
        hull->Y[iOffset] = hull->Y[i];
        #ifdef DEBUG
            hull->X[i] = 0;
            hull->Y[i] = 0;
        #endif
    }
    // now hull.X and hull.Y have exact spacing to fit the whole of new points to add without needing to be moved again

    // put new points in the empty spaces left by the above procedure
    for (size_t i = 0; i < hull->n; i++)
    {
        if (maxDistPtIndices[i] != -1)
        {
            size_t iOffset = i + offsetCounter[i] + 1;
            hull->X[iOffset] = uncoveredPts->X[maxDistPtIndices[i]];
            hull->Y[iOffset] = uncoveredPts->Y[maxDistPtIndices[i]];
        }
    }

    // sort maxDistPtIndices to safely remove points added to the hull from the uncovered set use selection sort because i'm lazy and don't want to write quicksort
    // first move all -1 (if any) to the end of the array to avoid worthless moves(especially useful in reducing complexity towards the last iterations)
    size_t i = 0;
    size_t j = hull->n - 1;
    while (i <= j)
    {
        if (maxDistPtIndices[i] == -1)
        {
            while ((j > i) && (maxDistPtIndices[j] == -1))
                j--;
            
            if (j == i) break;

            swapElems(maxDistPtIndices[i], maxDistPtIndices[j])
            j--;
        }
        i++;
    }
    // selection sort here
    for (size_t k = 0; k < i; k++)
    {
        size_t maxk = k;
        for (size_t l = k+1; l < i; l++)
            if (maxDistPtIndices[maxk] < maxDistPtIndices[l])
                maxk = l;
        if (k != maxk)
            swapElems(maxDistPtIndices[maxk], maxDistPtIndices[k]);
    }
    
    // remove element from hull
    for (size_t k = 0; k < i; k++)
    {
        #ifdef DEBUG
            if (maxDistPtIndices[k] == -1)
                throwError("p[%2d] t[%3d] addPtToHull: Unable to remove all -1 from maxDistPtIndices before the sorting. uncoveredPts=%ld, i=%ld, j=%ld", id->p, id->t, uncoveredPts->n, i, j);
        #endif
        uncoveredPts->n--;
        swapElems(uncoveredPts->X[maxDistPtIndices[k]], uncoveredPts->X[uncoveredPts->n]);
        swapElems(uncoveredPts->Y[maxDistPtIndices[k]], uncoveredPts->Y[uncoveredPts->n]);
    }

    size_t addedElemsCount = offsetCounter[hull->n];
    if (reallocMemory)
    {
        free(offsetCounter);
        offsetCounter = malloc(*allocatedElemsCount * 2 * sizeof(size_t) + MALLOC_PADDING*2);
        if (offsetCounter == NULL)
            throwError("p[%2d] t[%3d] addPtToHull: Failed to reallocate memory. uncoveredPts=%ld", id->p, id->t, uncoveredPts->n);
        *maxDistPtIndicesPtr = &offsetCounter[*allocatedElemsCount];
        *offsetCounterPtr = offsetCounter;
    }

    hull->n += addedElemsCount;
}


#ifdef DEBUG
int hullConvexityCheck(Data *hull, ProcThreadIDCombo *id)
{
    int retval = 0;
    for (size_t i = 0; i < hull->n; i++)
    {
        size_t ip1 = i+1;
        if (ip1 == hull->n) ip1 = 0;
        
        double a = (double)hull->X[ip1] - hull->X[i];
        double b = (double)hull->Y[i] - hull->Y[ip1];
        double c = (double)hull->X[i] * hull->Y[ip1] - (double)hull->X[ip1] * hull->Y[i];

        for (size_t j = 0; j < hull->n; j++)
        {
            if ((j == i) || (j == ip1)) continue;

            double dist = a * hull->Y[j] + b * hull->X[j] + c;
            if (dist <= 0)
            {
                LOG(LOG_LVL_ERROR, "p[%2d] t[%3d] quikchull-hullConvexityCheck: Hull is not convex. pt[%ld] is not to the left of line between pt[%ld] and pt[%ld]", id->p, id->t, j, i, ip1);
                retval = 1;
            }
        }
    }
    return retval;
}

int finalCoverageCheck(Data *hull, Data *pts, ProcThreadIDCombo *id)
{
    Data p = *pts;

    removeCoveredPoints(hull, &p, NULL, id);

    for (size_t i = 0; i < p.n; i++)
    {
        for (size_t j = 0; j < hull->n; j++)
        {
            if ((p.X[i] == hull->X[j]) && (p.Y[i] == hull->Y[j])) // point is part of the hull
            {
                p.n--;
                swapElems(p.X[i], p.X[p.n])
                swapElems(p.Y[i], p.Y[p.n])
                i--;
                break;
            }
        }
    }
    
    for (size_t i = 0; i < p.n; i++)
        LOG(LOG_LVL_ERROR, "p[%2d] t[%3d] quikchull-finalCoverageCheck: Hull does not cover point %ld", id->p, id->t, i);
    
    if (p.n == 0) return 0;
    else return 1;
}

#endif