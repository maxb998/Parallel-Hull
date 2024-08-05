#include "parallhull.h"

#include <math.h>
#include <pthread.h>

#ifdef NON_MPI_MODE
    #include <time.h>
    #include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#else
    //#include <mpi.h>
#endif
#if defined(PARALLHULL_STEP_DEBUG) || defined(PARALLHULL_MERGE_OUTPUT_PLOT)
    #include <stdio.h>
#endif

#define MAX_THREADS 256

typedef struct {
    Data fullData;
    Data *hulls;
    size_t reducedProblemUB; // upper bound on the size of the problem on which quickhull will run.
    int nThreads;
    ProcThreadIDCombo id;
    int *finishRecord;
    size_t *dataSize;
} ThreadData;

static void *parallhullThread(void *arg);
static Data mergeHulls(Data *h1, Data *h2, ProcThreadIDCombo *id);
static inline bool findNextMergePoint(Data *mergedH, Data *mainH, Data *altH, size_t *mainHIndex, size_t *altHIndex, ProcThreadIDCombo *id);
#ifdef DEBUG
    static inline bool hullConvexityCheck(Data *h, ProcThreadIDCombo *id);
    static inline bool mergeHullCoverageCheck(Data *h0, Data *h1, Data *h2, ProcThreadIDCombo *id);
#endif

Data parallhullThreaded(Data *d, size_t reducedProblemUB, int procID, int nThreads)
{
    double startTime;
    #ifdef NON_MPI_MODE
        struct timespec timeStruct;
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        startTime = cvtTimespec2Double(timeStruct);
    #else
    #endif

    int finishRecord[MAX_THREADS];
    for (int i = 0; i < MAX_THREADS; i++)
        finishRecord[i] = 0;

    // decide how many elements each thread gets(at the beginning)
    size_t dataSize[MAX_THREADS];
    double n = d->n;
    n = ceil(n/nThreads);
    for (int i = 0; i < nThreads-1; i++)
        dataSize[i] = (size_t)n;
    dataSize[nThreads-1] = d->n - ((size_t)n * (nThreads-1)); // last thread might have less elements due to non-divisibility of d->n with nThreads
    
    pthread_t threads[MAX_THREADS];
    ThreadData ds[MAX_THREADS];
    Data hulls[MAX_THREADS];
    for (int i = 0; i < nThreads; i++)
    {
        ds[i].fullData = *d;
        ds[i].hulls = hulls;
        ds[i].reducedProblemUB = reducedProblemUB;
        ds[i].nThreads = nThreads;
        ds[i].id.p = procID;
        ds[i].id.t = i;
        ds[i].finishRecord = finishRecord;
        ds[i].dataSize = dataSize;
        pthread_create(&threads[i], NULL, parallhullThread, (void*)&ds[i]);
    }

    for (int i = nThreads - 1; i > 0; i--)
    {
        pthread_join(threads[i], NULL);
        if (hulls[i].X - d->X > d->n)
            free(hulls[i].X);
    }
    pthread_join(threads[0], NULL); // join outside so that final hull memory alloc is not free

    #ifdef DEBUG
    #endif
    
    double finishTime;
    #ifdef NON_MPI_MODE
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        finishTime = cvtTimespec2Double(timeStruct);
        LOG(LOG_LVL_NOTICE, "Exec time of ParallHull using threads is %lfs", finishTime - startTime);
    #else
    #endif

    return hulls[0];
}


static void *parallhullThread(void *arg)
{
    ThreadData *thData = (ThreadData*)arg;
    int thID = thData->id.t;
    Data rd = { .n=thData->dataSize[thID] };
    {
        Data fd = thData->fullData;
        size_t startElem = 0;
        for (int i = 0; i < thID; i++)
            startElem += thData->dataSize[i];
        rd.X = &fd.X[startElem];
        rd.Y = &fd.Y[startElem];
    }

    // P1: each thread works on its own data in the first part here
    if (rd.n > thData->reducedProblemUB)
    {
        size_t nParts = (size_t)ceil((double)rd.n / thData->reducedProblemUB);
        Data *hulls = malloc(nParts * sizeof(Data));
        if (hulls == NULL)
            throwError("p[%2d] t[%3d] parallhullThread: Could not allocate %ld elements of 8 bytes in memory on thread %d", thData->id.p, thID, nParts * sizeof(Data), thID);

        size_t avgPartSize = (size_t)(ceil((double)rd.n / nParts));

        // P1.1: sequentially compute quickhull on each partition generated using the specified upper bound on the size of the rrd(Reduced Reduced problem Data)
        for (size_t i = 0; i < nParts; i++)
        {
            LOG(LOG_LVL_TRACE, "p[%2d] t[%3d] parallhullThread: Solving reduced problem %ld", thData->id.p, thID, i);

            size_t startPos = avgPartSize * i;

            hulls[i].n = avgPartSize;
            hulls[i].X = &rd.X[startPos];
            hulls[i].Y = &rd.Y[startPos];

            if (i == nParts-1)
                hulls[i].n = rd.n - avgPartSize * (nParts-1);

            hulls[i].n = quickhull(&hulls[i], &thData->id);
        }

        LOG(LOG_LVL_INFO, "p[%2d] t[%3d] parallhullThread: Quickhull on subproblem/s done, now merging", thData->id.p, thID);

        // P1.2: sequentially merge the convex hulls generated in every rrd
        while (nParts > 1)
        {
            size_t halfNParts = nParts / 2;
            for (size_t i = 0; i < halfNParts; i++)
            {
                LOG(LOG_LVL_TRACE, "p[%2d] t[%3d] parallhullThread: Merging thread internal hulls %ld and %ld", thData->id.p, thID, i, i + halfNParts);
                Data h = mergeHulls(&hulls[i], &hulls[i+halfNParts], &thData->id);

                #ifdef DEBUG
                    hullConvexityCheck(&h, &thData->id);
                    mergeHullCoverageCheck(&h, &hulls[i], &hulls[i + halfNParts], &thData->id);
                #endif

                if (hulls[i].X - rd.X > rd.n) // free memory allocated by mergehulls
                    free(hulls[i].X);
                if (hulls[i+halfNParts].X - rd.X > rd.n)
                    free(hulls[i+halfNParts].X);
                
                hulls[i] = h;
            }
            if (nParts & 1UL)
            {
                hulls[halfNParts] = hulls[nParts-1];
                nParts++;
            }
            nParts /= 2;
        }
        
        thData->hulls[thID] = hulls[0];
        
        free(hulls);
    }
    else
    {
        thData->hulls[thID].X = rd.X;
        thData->hulls[thID].Y = rd.Y;
        thData->hulls[thID].n = quickhull(&rd, &thData->id);
    }

    LOG(LOG_LVL_INFO, "p[%2d] t[%3d] parallhullThread: Thread subproblem solved", thData->id.p, thID);

    thData->finishRecord[thID] = 1;

    // P2: thread merge their results with each other in a ordered manner
    int i = 1, thID2merge = thID + 1;
    while ((thID % i == 0) && (thID2merge < thData->nThreads))
    {
        while (thData->finishRecord[thID2merge] < i) // spinlock (the assumption here is that threads should take more or less the same amount of time to merge, and this kinds of keeps the cpu "warm")
            __builtin_ia32_pause();
        
        Data h = mergeHulls(&thData->hulls[thID], &thData->hulls[thID2merge], &thData->id);

        LOG(LOG_LVL_INFO, "p[%2d] t[%3d] parallhullThread: Merging hull with hull in thread %d", thData->id.p, thID, thID2merge);

        #ifdef DEBUG
            hullConvexityCheck(&h, &thData->id);
            mergeHullCoverageCheck(&h, &thData->hulls[thID], &thData->hulls[thID2merge], &thData->id);
        #endif

        if (thData->hulls[thID].X - thData->fullData.X > thData->fullData.n)
            free(thData->hulls[thID].X);
        
        thData->hulls[thID] = h;

        i = i << 1;
        thID2merge = thID + i;

        // set thread finished computation
        thData->finishRecord[thID] = i;
    }

    thData->finishRecord[thID] = 0x7FFFFFFF; // cannot stall spinlock anymore
    
    return NULL;
}

static Data mergeHulls(Data *h1, Data *h2, ProcThreadIDCombo *id)
{
    size_t h1Index = 0, h2Index = 0;

    Data h0;
    h0.n = 1;
    h0.X = malloc((h1->n + h2->n) * 2 * sizeof(float));
    if (h0.X == NULL)
        throwError("p[%2d] t[%3d] mergeHulls: Failed to allocate memory for merged hull", id->p, id->t);
    h0.Y = &h0.X[h1->n + h2->n];

    bool mainHullIsH1, h1IsFirstMain;

    {   // set first point of h0 as the lowest one
    if (h1->Y[0] < h2->Y[0])
    {
        h0.X[0] = h1->X[0];
        h0.Y[0] = h1->Y[0];
        h1Index = 1;
        mainHullIsH1 = true;
    }
    else if (h1->Y[0] > h2->Y[0])
    {
        h0.X[0] = h2->X[0];
        h0.Y[0] = h2->Y[0];
        h2Index = 1;
        mainHullIsH1 = false;
    }
    else if (h1->X[0] >= h2->X[0]) // last two cases are here in case of equality (look at the X coordinate and choose the rightmost) (likely a rare occurence)
    {
        h0.X[0] = h1->X[0];
        h0.Y[0] = h1->Y[0];
        h1Index = 1;
        mainHullIsH1 = true;
    }
    else
    {
        h0.X[0] = h2->X[0];
        h0.Y[0] = h2->Y[0];
        h2Index = 1;
        mainHullIsH1 = false;
    }
    h1IsFirstMain = mainHullIsH1;
    }

    //################## Merge Procedure ##########################

    while ((h1Index <= h1->n) || (h2Index <= h2->n))
    {
        #ifdef PARALLHULL_STEP_DEBUG
            // show partial hull with gnuplot at each iteration and wait for user input to resume
            char plotTitle[200];
            sprintf(plotTitle, "Merged Hull: size=%lu, h1.n=%ld, h2.n=%ld, h1Index=%ld, h2Index=%ld", h0.n, h1->n, h2->n, h1Index, h2Index);
            LOG(LOG_LVL_DEBUG, "%s", plotTitle);
            plotHullMergeStep(h1, h2, &h0, h1Index, h2Index, plotTitle, false);
            getchar();
        #endif

        if (mainHullIsH1)
        {
            if (findNextMergePoint(&h0, h1, h2, &h1Index, &h2Index, id))
                mainHullIsH1 = ! mainHullIsH1;
        }
        else
        {
            if (findNextMergePoint(&h0, h2, h1, &h2Index, &h1Index, id))
                mainHullIsH1 = !mainHullIsH1;
        }

        #ifdef DEBUG
            if (mainHullIsH1 && (h1Index >= h1->n+1))
            {
                plotHullMergeStep(h1, h2, &h0, h1Index, h2Index, "Plot of the error", false);
                throwError("p[%2d] t[%3d] mergeHull: Cannot add points from an hull that has already been read completely. h1Index=%ld, h1Size=%ld, h2Index=%ld, h2Size=%ld, h0Size=%ld, h1IsMain=%d", id->p, id->t, h1Index, h1->n, h2Index, h2->n, h0.n, mainHullIsH1);
            }
            if (!mainHullIsH1 && (h2Index >= h2->n+1))
            {
                plotHullMergeStep(h1, h2, &h0, h1Index, h2Index, "Plot of the error", false);
                throwError("p[%2d] t[%3d] mergeHull: Cannot add points from an hull that has already been read completely. h1Index=%ld, h1Size=%ld, h2Index=%ld, h2Size=%ld, h0Size=%ld, h1IsMain=%d", id->p, id->t, h1Index, h1->n, h2Index, h2->n, h0.n, mainHullIsH1);
            }
        #endif

        if (mainHullIsH1)
        {
            size_t swapIndex = h1Index;
            if (h1Index > h1->n) break;
            if (h1Index == h1->n)
            { 
                if (!h1IsFirstMain) swapIndex = 0; 
                else break;
            }
            h0.X[h0.n] = h1->X[swapIndex];
            h0.Y[h0.n] = h1->Y[swapIndex];
            h1Index++;
        }
        else
        {
            size_t swapIndex = h2Index;
            if (h2Index > h2->n) break;
            if (h2Index == h2->n)
            { 
                if (h1IsFirstMain) swapIndex = 0; 
                else break;
            }
            h0.X[h0.n] = h2->X[swapIndex];
            h0.Y[h0.n] = h2->Y[swapIndex];
            h2Index++;
        }
        h0.n++;
    }
    
    #ifdef PARALLHULL_MERGE_OUTPUT_PLOT
        // show partial hull with gnuplot at each iteration and wait for user input to resume
        char plotTitle[200];
        sprintf(plotTitle, "Merged Hull: size=%lu, h1.n=%ld, h2.n=%ld", h0.n, h1->n, h2->n);
        LOG(LOG_LVL_DEBUG, "%s", plotTitle);
        plotHullMergeStep(h1, h2, &h0, h1Index, h2Index, plotTitle, true);
        getchar();
    #endif

    return h0;
}

static inline bool findNextMergePoint(Data *mergedH, Data *mainH, Data *altH, size_t *mainHIndex, size_t *altHIndex, ProcThreadIDCombo *id)
{
    #ifdef DEBUG
        if (*mainHIndex > mainH->n+1)
            throwError("p[%2d] t[%3d] findNextMergePoint: mainHIndex > mainH.n", id->p, id->t);
        if (*altHIndex > altH->n+1)
            throwError("p[%2d] t[%3d] findNextMergePoint: altHIndex > altH.n", id->p, id->t);
        // if ((*mainHIndex == mainH->n) && (*altHIndex == altH->n))
        //     throwError("findNextMergePoint: mainHIndex == mainH->n && altHIndex == altH->n");
    #endif

    bool swapMainToAlt=false;

    if (*altHIndex > altH->n+1)
        return false;

    size_t lastIndex = mergedH->n - 1;

    {
        size_t mainHIndexFix = *mainHIndex < mainH->n ? *mainHIndex : *mainHIndex % mainH->n;
        double a,b,c;
        a = (double)mainH->X[mainHIndexFix] - mergedH->X[lastIndex];
        b = (double)mergedH->Y[lastIndex] - mainH->Y[mainHIndexFix];
        c = (double)mergedH->X[lastIndex] * mainH->Y[mainHIndexFix] - (double)mainH->X[mainHIndexFix] * mergedH->Y[lastIndex];
        
        // now find min dist point starting from 
        double previousDist = INFINITY;
        size_t altHIndexFixed = *altHIndex < altH->n ? *altHIndex : *altHIndex % altH->n;
        double currentDist = a * altH->Y[altHIndexFixed] + b * altH->X[altHIndexFixed] + c;

        while ((previousDist >= currentDist) && (*altHIndex <= altH->n) && (currentDist >= 0)) // while distance is decreasing and we are not looking through already looked data
        {
            (*altHIndex)++;
            altHIndexFixed = *altHIndex < altH->n ? *altHIndex : *altHIndex % altH->n;
            previousDist = currentDist;
            currentDist = a * altH->Y[altHIndexFixed] + b * altH->X[altHIndexFixed] + c;
        }

        if (currentDist < 0)
            swapMainToAlt = true;
        else 
            (*altHIndex)--;
    }

    // need to check if point in altH[altHIndex] is the correct one for the chosen hull or not. It might not belong to the convex hull and be covered by the line casted a successive point in the hull
    if ((swapMainToAlt) && (*altHIndex < altH->n))
    {
        while (*altHIndex < altH->n)
        {
            double a,b,c;
            a = (double)altH->X[*altHIndex] - mergedH->X[lastIndex];
            b = (double)mergedH->Y[lastIndex] - altH->Y[*altHIndex];
            c = (double)mergedH->X[lastIndex] * altH->Y[*altHIndex] - (double)altH->X[*altHIndex] * mergedH->Y[lastIndex];

            size_t altHIndexP1 = *altHIndex + 1;
            if (altHIndexP1 >= altH->n)
                altHIndexP1 = 0;
            double dist = a * altH->Y[altHIndexP1] + b * altH->X[altHIndexP1] + c;

            if (dist > 0) break; // ignore dist == 0 since it is on the same line just get the point farther
            
            (*altHIndex)++;
        }
    }

    return swapMainToAlt;
}

#ifdef DEBUG
static inline bool hullConvexityCheck(Data *h, ProcThreadIDCombo *id)
{
    bool retval = false;
    for (size_t i = 0; i < h->n; i++)
    {
        size_t ip1 = i+1;
        if (ip1 == h->n)
            ip1=0;
        
        double a = (double)h->X[ip1] - h->X[i];
        double b = (double)h->Y[i] - h->Y[ip1];
        double c = (double)h->X[i] * h->Y[ip1] - (double)h->X[ip1] * h->Y[i];

        for (size_t j = 0; j < h->n; j++)
        {
            if ((j == i) || (j == ip1))
                continue;
            double dist = a * h->Y[j] + b * h->X[j] + c;
            if (dist < 0)
            {
                LOG(LOG_LVL_ERROR, "p[%2d] t[%3d] hullConvexityCheck: Merged Hull is not convex. pt[%ld] is not to the left of line between pt[%ld] and pt[%ld]", id->p, id->t, j, i, ip1);
                retval = true;
            }
        }
    }
    return retval;
}
static inline bool mergeHullCoverageCheck(Data *h0, Data *h1, Data *h2, ProcThreadIDCombo *id)
{
    bool retval = false;
    for (size_t i = 0; i < h0->n; i++)
    {
        size_t ip1 = i+1;
        if (ip1 == h0->n)
            ip1=0;
        
        double a = (double)h0->X[ip1] - h0->X[i];
        double b = (double)h0->Y[i] - h0->Y[ip1];
        double c = (double)h0->X[i] * h0->Y[ip1] - (double)h0->X[ip1] * h0->Y[i];

        for (size_t j = 0; j < h1->n; j++)
        {
            if (((h1->X[j] != h0->X[i]) && (h1->Y[j] != h0->Y[i])) ||
                ((h1->X[j] != h0->X[ip1]) && (h1->X[j] != h0->Y[ip1])))
                continue;
            double dist = a * h1->Y[j] + b * h1->X[j] + c;
            if (dist < 0)
            {
                LOG(LOG_LVL_ERROR, "p[%2d] t[%3d] mergeHullCoverageCheck: Merged Hull does not contain point[%ld] of h1. It is not to the right of line between pt[%ld] and pt[%ld]", id->p, id->t, j, i, ip1);
            }
        }
        for (size_t j = 0; j < h2->n; j++)
        {
            if (((h2->X[j] != h0->X[i]) && (h2->Y[j] != h0->Y[i])) ||
                ((h2->X[j] != h0->X[ip1]) && (h2->X[j] != h0->Y[ip1])))
                continue;
            double dist = a * h2->Y[j] + b * h2->X[j] + c;
            if (dist < 0)
            {
                LOG(LOG_LVL_ERROR, "p[%2d] t[%3d] mergeHullCoverageCheck: Merged Hull does not contain point[%ld] of h2. It is not to the right of line between pt[%ld] and pt[%ld]", id->p, id->t, j, i, ip1);
                retval = true;
            }
        }
    }
    return retval;
}
#endif
