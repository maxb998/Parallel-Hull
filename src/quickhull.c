#include "parallhull.h"

#define SMARTMODE

#include <math.h>
#ifdef NON_MPI_MODE
    #include <time.h>
    #include <stdio.h>
    #include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#endif

typedef struct 
{
    size_t node;
    size_t anchor;
} SuccessorData;

static void getExtremeCoordsPts(Data *d, size_t ptIndices[4]);
static size_t setExtremeCoordsAsInit(Data *d, size_t ptIndices[4]);
static size_t removeCoveredPointsSimple(Data *d, size_t hullSize, size_t nUncovered);
#ifndef SMARTMODE
    static SuccessorData findFarthestPtSimple(Data *d, size_t hullSize, size_t nUncovered);
    static void addPtToHullSimple(Data *d, SuccessorData pt, size_t hullSize, size_t nUncovered);
#else
    static size_t removeCoveredPointsSmart(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    static void setupMaxDistArrays(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    static SuccessorData findFarthestPtSmart(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    static void addPtToHullSmart(Data *d, SuccessorData pt, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    static void adaptMinArrays(Data *d, SuccessorData pt, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    #ifdef DEBUG
        static bool integrityCheck(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors, ProcThreadIDCombo *id);
    #endif
#endif
#ifdef DEBUG
    static bool hullConvexityCheck(Data *h, size_t hullSize, ProcThreadIDCombo *id);
#endif

size_t quickhull(Data *d, ProcThreadIDCombo *id)
{
    double startTime, iterTime, previousIterTime;
    int iterCount = 0, shiftAmount = -1;

    #ifdef NON_MPI_MODE
        struct timespec timeStruct;
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        startTime = cvtTimespec2Double(timeStruct);
        iterTime = startTime;
        previousIterTime = 0;
    #else
    #endif

    // init
    size_t ptIndices[4];
    getExtremeCoordsPts(d, ptIndices);
    size_t hullSize = setExtremeCoordsAsInit(d, ptIndices);
    
    size_t nUncovered = d->n - hullSize;
    nUncovered = removeCoveredPointsSimple(d, hullSize, nUncovered);

    #ifdef SMARTMODE
        float *minDistsAllocPtr = malloc(nUncovered * (sizeof(float) + sizeof(size_t)));
        if (minDistsAllocPtr == NULL)
            throwError("Could not allocate %lu bytes of memory for the minDists structures", nUncovered * sizeof(float) + nUncovered + sizeof(size_t));
            
        float *minDists = minDistsAllocPtr;
        size_t *minDistAnchors = (size_t*)(&minDists[nUncovered]);
        setupMaxDistArrays(d, hullSize, nUncovered, minDists, minDistAnchors);
    #endif

    while (nUncovered > 0)
    {
        iterCount++;   
        #ifdef NON_MPI_MODE 
            previousIterTime = iterTime;
            clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
            iterTime = cvtTimespec2Double(timeStruct);
        #endif

        inlineLOG(LOG_LVL_TRACE, "p[%2d] t[%3d] quickhull: Iteration %5d took %.3es. nUncovered=%.3e, hullSize=%ld", id->p, id->t, iterCount, iterTime-previousIterTime, (float)nUncovered, hullSize);
        if ((iterCount >> shiftAmount) || (iterCount < 3))
        {
            LOG(LOG_LVL_TRACE, "p[%2d] t[%3d] quickhull: Iteration %5d completed after %.3es from begining. nUncovered=%.3e, hullSize=%ld", id->p, id->t, iterCount, iterTime-startTime, (float)nUncovered, hullSize);
            shiftAmount++;
        }
        
        #if defined(QUICKHULL_STEP_DEBUG) && defined(NON_MPI_MODE)
            // show partial hull with gnuplot at each iteration and wait for user input to resume
            char plotTitle[200];
            sprintf(plotTitle, "Partial Hull: size=%lu, uncovered=%lu", hullSize, nUncovered);
            Data hull = { .n=hullSize, .X=d->X, .Y=d->Y };
            plotData(d, &hull, nUncovered, plotTitle);
            getchar();
        #endif

        #ifndef SMARTMODE
            SuccessorData s = findFarthestPtSimple(d, hullSize, nUncovered);
            addPtToHullSimple(d, s, hullSize, nUncovered);
            hullSize++;
            nUncovered--;
            nUncovered = removeCoveredPointsSimple(d, hullSize, nUncovered);
        #else
            SuccessorData s = findFarthestPtSmart(d, hullSize, nUncovered, minDists, minDistAnchors);
            addPtToHullSmart(d, s, hullSize, nUncovered, minDists, minDistAnchors);
            minDists++; minDistAnchors++;
            hullSize++;
            #ifdef DEBUG
                size_t oldNUncovered = nUncovered;
            #endif

            nUncovered--;
            nUncovered = removeCoveredPointsSmart(d, hullSize, nUncovered, minDists, minDistAnchors);
            if (hullSize <= 3)
                setupMaxDistArrays(d, hullSize, nUncovered, minDists, minDistAnchors);
            else
                adaptMinArrays(d, s, hullSize, nUncovered, minDists, minDistAnchors);
            #ifdef DEBUG
                if (integrityCheck(d, hullSize, nUncovered, minDists, minDistAnchors, id))
                    throwError("p[%2d] t[%3d] quickhull: Error detected in values inside minDist* arrays. n=%ld, hullSize=%ld, nUncovered=%ld, oldNUncovered=%ld, s.node=%ld, s.anchor=%ld", id->p, id->t, d->n, hullSize, nUncovered, oldNUncovered, s.node, s.anchor);
            #endif
        #endif

        #ifdef DEBUG
            if (hullConvexityCheck(d, hullSize, id))
                throwError("p[%2d] t[%3d] quickhull: Hull is not convex. n=%ld, hullSize=%ld, nUncovered=%ld, oldNUncovered=%ld, s.node=%ld, s.anchor=%ld", id->p, id->t, d->n, hullSize, nUncovered, oldNUncovered, s.node, s.anchor);

        #endif
    }

    LOG(LOG_LVL_TRACE, "p[%2d] t[%3d] quickhull: Iteration %5d took %.3es. nUncovered=%.3e, hullSize=%ld", id->p, id->t, iterCount, iterTime-previousIterTime, (float)nUncovered, hullSize);

    #ifdef DEBUG
        LOG(LOG_LVL_DEBUG, "p[%2d] t[%3d] quickhull: DEBUG macro is defined! Now checking whether all points are actually inside the hull", id->p, id->t);
        nUncovered = d->n - hullSize;
        nUncovered = removeCoveredPointsSimple(d, hullSize, nUncovered);
        if (nUncovered != 0)
            throwError("p[%2d] t[%3d] quickhull: There are still %ld points that are not inside the hull", id->p, id->t, nUncovered);
    #endif

    #ifdef SMARTMODE
        free(minDistsAllocPtr);
    #endif
    return hullSize;
}

static void getExtremeCoordsPts(Data *d, size_t ptIndices[4])
{
    for (int i = 0; i < 4; i++)
        ptIndices[i] = 0;
    
    for (size_t i = 1; i < d->n; i++)
    {
        if  ((d->X[i] > d->X[ptIndices[1]]) || 
            ((d->X[i] == d->X[ptIndices[1]]) && (d->Y[i] > d->Y[ptIndices[1]]))) //xMax rightmust upward if multiple choices
            ptIndices[1] = i;
        else if ((d->X[i] < d->X[ptIndices[3]]) || 
                ((d->X[i] == d->X[ptIndices[3]]) && (d->Y[i] < d->Y[ptIndices[3]]))) //xMin leftmost loward if multiple choices
            ptIndices[3] = i;

        if ((d->Y[i] > d->Y[ptIndices[2]]) || 
            ((d->Y[i] == d->Y[ptIndices[2]]) && (d->X[i] < d->X[ptIndices[2]]))) //yMax upmost leftward if multiple choices
            ptIndices[2] = i;
        else if ((d->Y[i] < d->Y[ptIndices[0]]) || 
                ((d->Y[i] == d->Y[ptIndices[0]]) && (d->X[i] > d->X[ptIndices[0]]))) //yMin lowmost rightward if multiple choices
            ptIndices[0] = i;
    }
}

static size_t setExtremeCoordsAsInit(Data *d, size_t ptIndices[4])
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
    
    size_t hullSize = 0;
    for (int i = 0; (i < 4) && (ptIndices[i] != -1); i++)
    {
        // adaptation in case of swapping an element that is already in the top 4 positions
        for (int j = i+1; j < 4; j++)
            if (ptIndices[j] == i)
                ptIndices[j] = ptIndices[i];
        
        swapElems(d->X[i], d->X[ptIndices[i]]);
        swapElems(d->Y[i], d->Y[ptIndices[i]]);

        hullSize++;
    }

    return hullSize;
}

static size_t removeCoveredPointsSimple(Data *d, size_t hullSize, size_t nUncovered)
{
    size_t i = hullSize;
    size_t j = nUncovered + hullSize - 1; // position of last covered element in d->X and d->Y
    while (i <= j)
    {
        bool iIsCovered = true;
        for (size_t k1 = 0, k2 = hullSize-1; k1 < hullSize; k2 = k1++)
        {
            if ((d->Y[k1] < d->Y[i]) != (d->Y[k2] < d->Y[i]))
            {
                double a = (double)d->X[k1] - d->X[k2];
                double b = (double)d->Y[k2] - d->Y[k1];
                double c = (double)d->X[k2] * d->Y[k1] - (double)d->X[k1] * d->Y[k2];
                double dist = a * d->Y[i] + b * d->X[i] + c;
                if (dist >= 0)
                    iIsCovered = !iIsCovered;
            }
        }

        if (iIsCovered)
        {
            while(i < j)
            {
                bool jIsCovered = true;
                for (size_t k1 = 0, k2 = hullSize-1; k1 < hullSize; k2 = k1++)
                {
                    if ((d->Y[k1] < d->Y[j]) != (d->Y[k2] < d->Y[j]))
                    {
                        double a = (double)d->X[k1] - d->X[k2];
                        double b = (double)d->Y[k2] - d->Y[k1];
                        double c = (double)d->X[k2] * d->Y[k1] - (double)d->X[k1] * d->Y[k2];
                        double dist = a * d->Y[j] + b * d->X[j] + c;
                        if (dist >= 0)
                            jIsCovered = !jIsCovered;
                    }
                }
                if (jIsCovered)
                    j--;
                else break;
            }

            if (j == i) break; // all points up to i included are covered by the current hull

            swapElems(d->X[i], d->X[j]);
            swapElems(d->Y[i], d->Y[j]);
            j--;
            i++; // if j was covered there is no need to check i in the next iteration
        }
        else
            i++;
    }

    return i - hullSize;
}

#ifndef SMARTMODE

static SuccessorData findFarthestPtSimple(Data *d, size_t hullSize, size_t nUncovered)
{
    SuccessorData farPt = { .node=0, .anchor=-1 };
    float farPtDist = -1;

    for (size_t i = hullSize; i < hullSize + nUncovered; i++)
    {
        size_t closestAnchor = 0;
        float closestDist = INFINITY;
        for (size_t j = 0, k = hullSize-1; j < hullSize; k=j++)
        {
            float dist0 = d->X[j] - d->X[k];
            float dist1 = d->X[k] - d->X[i];
            float dist2 = d->X[i] - d->X[j];
            float yDiff0 = d->Y[j] - d->Y[k];
            float yDiff1 = d->Y[k] - d->Y[i];
            float yDiff2 = d->Y[i] - d->Y[j];
            dist0 = sqrtf(dist0 * dist0 + yDiff0 * yDiff0);
            dist1 = sqrtf(dist1 * dist1 + yDiff1 * yDiff1);
            dist2 = sqrtf(dist2 * dist2 + yDiff2 * yDiff2);
            float dist = dist1 + dist2 - dist0;
            if (closestDist > dist)
            {
                closestDist = dist;
                closestAnchor = k;
            }
        }
        if (farPtDist < closestDist)
        {
            farPtDist = closestDist;
            farPt.node = i;
            farPt.anchor = closestAnchor;
        }
    }

    #ifdef DEBUG
        if (((long)farPt.node < (long)hullSize) || ((long)farPt.node > (long)hullSize + nUncovered) || ((long)farPt.anchor > (long)hullSize) || ((long)farPt.anchor < 0))
            throwError("Invalid anchor(%ld) and node(%ld) @ hullsize=%ld, nUncovered=%ld", farPt.anchor, farPt.node, hullSize, nUncovered);
    #endif

    return farPt;
}

static void addPtToHullSimple(Data *d, SuccessorData pt, size_t hullSize, size_t nUncovered)
{
    swapElems(d->X[pt.node], d->X[hullSize]);
    swapElems(d->Y[pt.node], d->Y[hullSize]);
    
    size_t anchor = pt.anchor;
    if (hullSize == 2) //  need to check on which side the point must be added
    {
        double a = (double)d->X[1] - d->X[0];
        double b = (double)d->Y[0] - d->Y[1];
        double c = (double)d->X[0] * d->Y[1] - (double)d->X[1] * d->Y[0];
        double dist = a * d->Y[2] + b * d->X[2] + c;
        if (dist < 0)
            anchor = 0;
    }

    float xBak = d->X[hullSize], yBak = d->Y[hullSize];
    for (size_t i = hullSize-1; i > anchor; i--)
    {
        d->X[i+1] = d->X[i];
        d->Y[i+1] = d->Y[i];
    }

    d->X[anchor+1] = xBak;
    d->Y[anchor+1] = yBak;
}
#endif

#ifdef SMARTMODE
static SuccessorData findFarthestPtSmart(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors)
{
    SuccessorData farPt = { .node=hullSize, .anchor=minDistAnchors[0] };
    float farPtDist = minDists[0];

    for (size_t i = 1; i < nUncovered; i++)
    {
        if (farPtDist < minDists[i])
        {
            farPtDist = minDists[i];
            farPt.node = i + hullSize;
            farPt.anchor = minDistAnchors[i];
        }
    }

    #ifdef DEBUG
        if (((long)farPt.node < (long)hullSize) || ((long)farPt.node > (long)hullSize + nUncovered) || ((long)farPt.anchor > (long)hullSize) || ((long)farPt.anchor < 0))
            throwError("Invalid anchor(%ld) and node(%ld) @ hullsize=%ld, nUncovered=%ld", farPt.anchor, farPt.node, hullSize, nUncovered);
    #endif

    return farPt;
}

static size_t removeCoveredPointsSmart(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors)
{
    size_t i = hullSize;
    size_t j = nUncovered + hullSize - 1;
    while (i <= j)
    {
        bool iIsCovered = true;
        for (size_t k1 = 0, k2 = hullSize-1; k1 < hullSize; k2 = k1++)
        {
            if ((d->Y[k1] < d->Y[i]) != (d->Y[k2] < d->Y[i]))
            {
                double a = (double)d->X[k1] - d->X[k2];
                double b = (double)d->Y[k2] - d->Y[k1];
                double c = (double)d->X[k2] * d->Y[k1] - (double)d->X[k1] * d->Y[k2];
                double dist = a * d->Y[i] + b * d->X[i] + c;
                if (dist >= 0)
                    iIsCovered = !iIsCovered;
            }
        }

        if (iIsCovered)
        {
            while(i < j)
            {
                bool jIsCovered = true;
                for (size_t k1 = 0, k2 = hullSize-1; k1 < hullSize; k2 = k1++)
                {
                    if ((d->Y[k1] < d->Y[j]) != (d->Y[k2] < d->Y[j]))
                    {
                        double a = (double)d->X[k1] - d->X[k2];
                        double b = (double)d->Y[k2] - d->Y[k1];
                        double c = (double)d->X[k2] * d->Y[k1] - (double)d->X[k1] * d->Y[k2];
                        double dist = a * d->Y[j] + b * d->X[j] + c;
                        if (dist >= 0)
                            jIsCovered = !jIsCovered;
                    }
                }
                if (jIsCovered)
                    j--;
                else break;
            }
            
            if (j == i) break; // all points up to i included are covered by the current hull

            swapElems(d->X[i], d->X[j]);
            swapElems(d->Y[i], d->Y[j]);
            minDists[i - hullSize] = minDists[j - hullSize];
            minDistAnchors[i - hullSize] = minDistAnchors[j - hullSize];
            j--;
            i++;
        }
        else
            i++;
    }

    return i - hullSize;
}

static void setupMaxDistArrays(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors)
{
    for (size_t i = hullSize; i < hullSize + nUncovered; i++)
    {
        size_t i2 = i - hullSize;
        minDists[i2] = INFINITY;
        minDistAnchors[i2] = 0;
        for (size_t j = 0, k = hullSize-1; j < hullSize; k=j++)
        {
            float dist0 = d->X[j] - d->X[k];
            float dist1 = d->X[k] - d->X[i];
            float dist2 = d->X[i] - d->X[j];
            float yDiff0 = d->Y[j] - d->Y[k];
            float yDiff1 = d->Y[k] - d->Y[i];
            float yDiff2 = d->Y[i] - d->Y[j];
            dist0 = sqrtf(dist0 * dist0 + yDiff0 * yDiff0);
            dist1 = sqrtf(dist1 * dist1 + yDiff1 * yDiff1);
            dist2 = sqrtf(dist2 * dist2 + yDiff2 * yDiff2);
            float dist = dist1 + dist2 - dist0;
            if (minDists[i2] > dist)
            {
                minDists[i2] = dist;
                minDistAnchors[i2] = k;
            }
        }
        #ifdef DEBUG
            if (minDists[i2] == INFINITY)
                throwError("Could not find the closest edge for %d uncovered node", i2);
        #endif
    }
}

static void addPtToHullSmart(Data *d, SuccessorData pt, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors)
{
    // move the node closer to the hull
    swapElems(d->X[hullSize], d->X[pt.node]);
    swapElems(d->Y[hullSize], d->Y[pt.node]);
    minDists[pt.node-hullSize] = minDists[0];
    minDistAnchors[pt.node-hullSize] = minDistAnchors[0];

    size_t anchor = pt.anchor;
    if (hullSize == 2) //  need to check on which side the point must be added
    {
        double a = (double)d->X[1] - d->X[0];
        double b = (double)d->Y[0] - d->Y[1];
        double c = (double)d->X[0] * d->Y[1] - (double)d->X[1] * d->Y[0];
        double dist = a * d->Y[2] + b * d->X[2] + c;
        if (dist < 0)
            anchor = 0;
    }

    float xBak = d->X[hullSize], yBak = d->Y[hullSize];
    for (size_t i = hullSize-1; i > anchor; i--)
    {
        d->X[i+1] = d->X[i];
        d->Y[i+1] = d->Y[i];
    }

    d->X[anchor+1] = xBak;
    d->Y[anchor+1] = yBak;
}

static void adaptMinArrays(Data *d, SuccessorData pt, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors)
{
    // recompute min distance wherever it is necessary
    size_t j0 = pt.anchor, j1 = (pt.anchor+1) % hullSize, j2 = (pt.anchor+2) % hullSize;
    for (size_t i = 0; i < nUncovered; i++)
    {
        if (minDistAnchors[i] == pt.anchor)
        {
            size_t k = i + hullSize;
            minDists[i] = INFINITY;
            {
                float dist0 = d->X[j0] - d->X[j1];
                float dist1 = d->X[j1] - d->X[k];
                float dist2 = d->X[k] - d->X[j0];
                float yDiff0 = d->Y[j0] - d->Y[j1];
                float yDiff1 = d->Y[j1] - d->Y[k];
                float yDiff2 = d->Y[k] - d->Y[j0];
                dist0 = sqrtf(dist0 * dist0 + yDiff0 * yDiff0);
                dist1 = sqrtf(dist1 * dist1 + yDiff1 * yDiff1);
                dist2 = sqrtf(dist2 * dist2 + yDiff2 * yDiff2);
                minDists[i] = dist1 + dist2 - dist0;
            }
            {
                float dist0 = d->X[j1] - d->X[j2];
                float dist1 = d->X[j2] - d->X[k];
                float dist2 = d->X[k] - d->X[j1];
                float yDiff0 = d->Y[j1] - d->Y[j2];
                float yDiff1 = d->Y[j2] - d->Y[k];
                float yDiff2 = d->Y[k] - d->Y[j1];
                dist0 = sqrtf(dist0 * dist0 + yDiff0 * yDiff0);
                dist1 = sqrtf(dist1 * dist1 + yDiff1 * yDiff1);
                dist2 = sqrtf(dist2 * dist2 + yDiff2 * yDiff2);
                float dist = dist1 + dist2 - dist0;
                if (dist < minDists[i])
                {
                    minDists[i] = dist;
                    minDistAnchors[i] = j1;
                }
            }
        }
        else if (minDistAnchors[i] > pt.anchor)
            minDistAnchors[i]++;
    }
}

#ifdef DEBUG
static bool integrityCheck(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors, ProcThreadIDCombo *id)
{
    bool quit = false;
    for (size_t i = hullSize; i < hullSize + nUncovered; i++)
    {
        size_t closestAnchor = 0;
        float closestDist = INFINITY;
        for (size_t j = 0, k = hullSize-1; j < hullSize; k=j++)
        {
            float dist0 = d->X[j] - d->X[k];
            float dist1 = d->X[k] - d->X[i];
            float dist2 = d->X[i] - d->X[j];
            float yDiff0 = d->Y[j] - d->Y[k];
            float yDiff1 = d->Y[k] - d->Y[i];
            float yDiff2 = d->Y[i] - d->Y[j];
            dist0 = sqrtf(dist0 * dist0 + yDiff0 * yDiff0);
            dist1 = sqrtf(dist1 * dist1 + yDiff1 * yDiff1);
            dist2 = sqrtf(dist2 * dist2 + yDiff2 * yDiff2);
            float dist = dist1 + dist2 - dist0;
            if (closestDist > dist)
            {
                closestDist = dist;
                closestAnchor = k;
            }
        }
        
        size_t i2 = i - hullSize;
        if ((minDists[i2] != closestDist) && (minDistAnchors[i2] != closestAnchor))
        {
            quit = true;
            LOG(LOG_LVL_ERROR, "p[%2d] t[%3d] integrityCheck: @[%ld] minDists=%.3e instead of %.3e,\t minDistAnchors=%ld instead of %ld", id->p, id->t, i2, minDists[i2], closestDist, minDistAnchors[i2], closestAnchor);
        }
        else if (minDists[i2] != closestDist)
        {
            quit = true;
            LOG(LOG_LVL_ERROR, "p[%2d] t[%3d] integrityCheck: @[%ld] minDists=%.3e instead of %.3e", id->p, id->t, i2, minDists[i2], closestDist);
        }
        else if (minDistAnchors[i2] != closestAnchor)
        {
            quit = true;
            LOG(LOG_LVL_ERROR, "p[%2d] t[%3d] integrityCheck: @[%ld] minDistAnchors=%ld instead of %ld", id->p, id->t, i2, minDistAnchors[i2], closestAnchor);
        }
    }
    
    return quit;
}
#endif
#endif

#ifdef DEBUG
static bool hullConvexityCheck(Data *h, size_t hullSize, ProcThreadIDCombo *id)
{
    bool retval = false;
    for (size_t i = 0; i < hullSize; i++)
    {
        size_t ip1 = i+1;
        if (ip1 == hullSize)
            ip1=0;
        
        double a = (double)h->X[ip1] - h->X[i];
        double b = (double)h->Y[i] - h->Y[ip1];
        double c = (double)h->X[i] * h->Y[ip1] - (double)h->X[ip1] * h->Y[i];

        for (size_t j = 0; j < hullSize; j++)
        {
            if ((j == i) || (j == ip1))
                continue;
            double dist = a * h->Y[j] + b * h->X[j] + c;
            if (dist < 0)
            {
                LOG(LOG_LVL_ERROR, "p[%2d] t[%3d] hullConvexityCheck: Hull is not convex. pt[%ld] is not to the left of line between pt[%ld] and pt[%ld]", id->p, id->t, j, i, ip1);
                retval = true;
            }
        }
    }
    return retval;
}
#endif