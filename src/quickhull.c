#include "paralhull.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

typedef struct 
{
    size_t node;
    size_t anchor;
} SuccessorData;

static inline void getExtremeCoordsPts(Data *d, size_t ptIndices[4]);
static inline size_t setExtremeCoordsAsInit(Data *d, size_t ptIndices[4]);
static void removeConveredPoints(Data *d, size_t hullSize, size_t *nUncovered);
static inline SuccessorData findFarthestPt(Data *d, size_t hullSize, size_t nUncovered);
static inline void addPtToHull(Data *d, SuccessorData pt, size_t hullSize);

int quickhull(Data *d)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // init
    size_t ptIndices[4];
    getExtremeCoordsPts(d, ptIndices);
    size_t hullSize = setExtremeCoordsAsInit(d, ptIndices);
    
    size_t nUncovered = d->n - hullSize;
    removeConveredPoints(d, hullSize, &nUncovered);

    int iterCount = 0;
    double iterTime=startTime, previousIterTime;
    while (nUncovered > 0)
    {
        #ifdef LOCAL_DEBUG
            // show partial hull with gnuplot at each iteration and wait for user input to resume
            char plotTitle[200];
            sprintf(plotTitle, "Partial Hull: size=%d, uncovered=%d", hullSize, nUncovered);
            plotData(d, hullSize, nUncovered, "1920,1080", plotTitle);
            getchar();
        #endif

        SuccessorData s = findFarthestPt(d, hullSize, nUncovered);
        addPtToHull(d, s, hullSize);
        hullSize++;

        removeConveredPoints(d, hullSize, &nUncovered);

        previousIterTime = iterTime;
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        iterTime = cvtTimespec2Double(timeStruct);
        LOG(LOG_LVL_INFO, "iteration %d took %lfs", iterCount, iterTime-previousIterTime);

        iterCount++;
    }

    return hullSize;
}

static inline void getExtremeCoordsPts(Data *d, size_t ptIndices[4])
{
    for (int i = 0; i < 4; i++)
        ptIndices[i] = 0;
    
    for (size_t i = 1; i < d->n; i++)
    {
        if (d->X[i] > d->X[ptIndices[1]]) //xMax
            ptIndices[1] = i;
        else if (d->X[i] < d->X[ptIndices[3]]) //xMin
            ptIndices[3] = i;

        if (d->Y[i] > d->Y[ptIndices[2]]) //yMax
            ptIndices[2] = i;
        else if (d->Y[i] < d->Y[ptIndices[0]]) //yMin
            ptIndices[0] = i;
    }
}

static inline size_t setExtremeCoordsAsInit(Data *d, size_t ptIndices[4])
{    
    // check extreme coords for duplicates
    if (ptIndices[2] == ptIndices[0])
        ptIndices[2] = -1;
    if (ptIndices[3] == ptIndices[0])
        ptIndices[3] = -1;
    if (ptIndices[2] == ptIndices[1])
        ptIndices[2] = -1;
    if (ptIndices[3] == ptIndices[1])
        ptIndices[3] = -1;
    
    size_t hullSize = 0;
    for (int i = 0; i < 4; i++)
    {
        if (ptIndices[i] == -1)
            continue;

        int minSubIndex = 0;
        for (int j = 1; j < 4; j++)
            if (ptIndices[minSubIndex] > ptIndices[j] && ptIndices[j] != -1)
                minSubIndex = j;

        swapElems(d->X[minSubIndex], d->X[ptIndices[minSubIndex]]);
        swapElems(d->Y[minSubIndex], d->Y[ptIndices[minSubIndex]]);

        ptIndices[minSubIndex] = __INT_MAX__;

        hullSize++;
    }

    return hullSize;
}

static void removeConveredPoints(Data *d, size_t hullSize, size_t *nUncovered)
{
    size_t i = hullSize;
    size_t j = *nUncovered + hullSize - 1; // position of last covered element in d->X and d->Y
    while (i <= j)
    {
        bool iIsCovered = true;
        for (size_t k1 = 0, k2 = hullSize-1; k1 < hullSize; k2 = k1++)
        {
            if ((d->Y[k1] < d->Y[i]) != (d->Y[k2] < d->Y[i]))
            {
                float a = d->X[k1] - d->X[k2];
                float b = d->Y[k2] - d->Y[k1];
                float c = d->X[k2] * d->Y[k1] - d->X[k1] * d->Y[k2];
                if (a * d->Y[i] + b * d->X[i] + c > 0)
                    iIsCovered = !iIsCovered;
            }
        }

        if (iIsCovered)
        {
            swapElems(d->X[i], d->X[j]);
            swapElems(d->Y[i], d->Y[j]);
            j--;
        }
        else
            i++;
    }

    *nUncovered = i - hullSize;
}

static inline SuccessorData findFarthestPt(Data *d, size_t hullSize, size_t nUncovered)
{
    SuccessorData farthestPt = { .node=0, .anchor=-1 };
    float farthestPtDist = -1;

    for (size_t i = hullSize; i < hullSize + nUncovered; i++)
    {
        size_t closestAnchor = 0;
        float closestDist = INFINITY;
        for (size_t j = 0, k = hullSize-1; j < hullSize; k=j++)
        {
            float xDiff0 = d->X[j] - d->X[k];
            float xDiff1 = d->X[k] - d->X[i];
            float xDiff2 = d->X[i] - d->X[j];
            float yDiff0 = d->Y[j] - d->Y[k];
            float yDiff1 = d->Y[k] - d->Y[i];
            float yDiff2 = d->Y[i] - d->Y[j];
            float dist0 = sqrtf(xDiff0 * xDiff0 + yDiff0 * yDiff0);
            float dist1 = sqrtf(xDiff1 * xDiff1 + yDiff1 * yDiff1);
            float dist2 = sqrtf(xDiff2 * xDiff2 + yDiff2 * yDiff2);
            float dist = dist1 + dist2 - dist0;
            if (closestDist > dist)
            {
                closestDist = dist;
                closestAnchor = k;
            }
        }
        if (farthestPtDist < closestDist)
        {
            farthestPtDist = closestDist;
            farthestPt.node = i;
            farthestPt.anchor = closestAnchor;
        }
    }

    if (farthestPt.node == 0 || farthestPt.anchor == -1)
        throwError("Invalid anchor and node @ hullsize=%d, nUncovered=%d", hullSize, nUncovered);

    return farthestPt;
}

static inline void addPtToHull(Data *d, SuccessorData pt, size_t hullSize)
{
    swapElems(d->X[pt.node], d->X[hullSize]);
    swapElems(d->Y[pt.node], d->Y[hullSize]);

    float xBak = d->X[hullSize], yBak = d->Y[hullSize];
    for (size_t i = hullSize-1; i > pt.anchor; i--)
    {
        d->X[i+1] = d->X[i];
        d->Y[i+1] = d->Y[i];
    }

    d->X[pt.anchor+1] = xBak;
    d->Y[pt.anchor+1] = yBak;
}
