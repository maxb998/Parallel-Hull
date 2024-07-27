#include "paralhull.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

//#define SMARTMODE

typedef struct 
{
    size_t node;
    size_t anchor;
} SuccessorData;

static void getExtremeCoordsPts(Data *d, size_t ptIndices[4]);
static size_t setExtremeCoordsAsInit(Data *d, size_t ptIndices[4]);
static size_t removeConveredPointsNaive(Data *d, size_t hullSize, size_t nUncovered);
#ifndef SMARTMODE
    static SuccessorData findFarthestPtNaive(Data *d, size_t hullSize, size_t nUncovered);
    static void addPtToHullNaive(Data *d, SuccessorData pt, size_t hullSize, size_t nUncovered);
#else
    static size_t removeConveredPointsSmart(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    static void setupMaxDistArrays(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    static SuccessorData findFarthestPtSmart(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    static void addPtToHullSmart(Data *d, SuccessorData pt, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    static void adaptMinArrays(Data *d, SuccessorData pt, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    #ifdef DEBUG
        static bool integrityCheck(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors);
    #endif
#endif

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
    nUncovered = removeConveredPointsNaive(d, hullSize, nUncovered);

    #ifdef SMARTMODE
        float *minDists = malloc(nUncovered * (sizeof(float) + sizeof(size_t)));
        if (minDists == NULL)
            throwError("Could not allocate %lu bytes of memory for the minDists structures", nUncovered * sizeof(float) + nUncovered + sizeof(size_t));
        size_t *minDistAnchors = (size_t*)(&minDists[nUncovered]);
        setupMaxDistArrays(d, hullSize, nUncovered, minDists, minDistAnchors);
    #endif

    int iterCount = 0;
    double iterTime=startTime, previousIterTime;
    while (nUncovered > 0)
    {
        #ifdef LOCAL_DEBUG
            // show partial hull with gnuplot at each iteration and wait for user input to resume
            char plotTitle[200];
            sprintf(plotTitle, "Partial Hull: size=%lu, uncovered=%lu", hullSize, nUncovered);
            plotData(d, hullSize, nUncovered, "1920,1080", plotTitle);
            getchar();
        #endif

        #ifndef SMARTMODE
            SuccessorData s = findFarthestPtNaive(d, hullSize, nUncovered);
            addPtToHullNaive(d, s, hullSize, nUncovered);
            hullSize++;
            nUncovered = removeConveredPointsNaive(d, hullSize, nUncovered);
        #else
            SuccessorData s = findFarthestPtSmart(d, hullSize, nUncovered, minDists, minDistAnchors);
            addPtToHullSmart(d, s, hullSize, nUncovered, minDists, minDistAnchors);
            hullSize++;
            #ifdef DEBUG
                size_t oldNUncovered = nUncovered;
            #endif
            nUncovered--;
            nUncovered = removeConveredPointsSmart(d, hullSize, nUncovered, minDists, minDistAnchors);
            adaptMinArrays(d, s, hullSize, nUncovered, minDists, minDistAnchors);
            #ifdef DEBUG
                if (integrityCheck(d, hullSize, nUncovered, minDists, minDistAnchors))
                    throwError("Error detected in values inside minDist* arrays. n=%ld, hullSize=%ld, nUncovered=%ld, oldNUncovered=%ld, s.node=%ld, s.anchor=%ld", d->n, hullSize, nUncovered, oldNUncovered, s.node, s.anchor);
            #endif
        #endif

        previousIterTime = iterTime;
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        iterTime = cvtTimespec2Double(timeStruct);
        LOG(LOG_LVL_INFO, "iteration %d took %lfs", iterCount, iterTime-previousIterTime);
        LOG(LOG_LVL_TRACE, "nUncovered=%ld, hullSize=%ld", nUncovered, hullSize);

        iterCount++;
    }

    #ifdef SMARTMODE
        free(minDists);
    #endif
    return hullSize;
}

static void getExtremeCoordsPts(Data *d, size_t ptIndices[4])
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

static size_t setExtremeCoordsAsInit(Data *d, size_t ptIndices[4])
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

static size_t removeConveredPointsNaive(Data *d, size_t hullSize, size_t nUncovered)
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

    return i - hullSize;
}

#ifndef SMARTMODE
static SuccessorData findFarthestPtNaive(Data *d, size_t hullSize, size_t nUncovered)
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

static void addPtToHullNaive(Data *d, SuccessorData pt, size_t hullSize, size_t nUncovered)
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

static size_t removeConveredPointsSmart(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors)
{
    size_t i = hullSize;// i2 = 0;
    size_t j = nUncovered + hullSize - 1;// j2 = nUncovered; // position of last covered element in d->X and d->Y
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
            minDists[i - hullSize] = minDists[j - hullSize];
            minDistAnchors[i - hullSize] = minDistAnchors[j - hullSize];
            j--; //j2--;
        }
        else
        {
            i++; //i2++;
        }
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
    for (size_t i = 0; i < nUncovered-1; i++)
    {
        minDists[i] = minDists[i+1];
        minDistAnchors[i] = minDistAnchors[i+1];
    }

    float xBak = d->X[hullSize], yBak = d->Y[hullSize];
    for (size_t i = hullSize-1; i > pt.anchor; i--)
    {
        d->X[i+1] = d->X[i];
        d->Y[i+1] = d->Y[i];
    }

    d->X[pt.anchor+1] = xBak;
    d->Y[pt.anchor+1] = yBak;
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
                //minDistAnchors[i] = j; // not needed
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
static bool integrityCheck(Data *d, size_t hullSize, size_t nUncovered, float *minDists, size_t *minDistAnchors)
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
            LOG(LOG_LVL_ERROR, "@[%ld] minDists=%.3e instead of %.3e,\t minDistAnchors=%ld instead of %ld", i2, minDists[i2], closestDist, minDistAnchors[i2], closestAnchor);
        }
        else if (minDists[i2] != closestDist)
        {
            quit = true;
            LOG(LOG_LVL_ERROR, "@[%ld] minDists=%.3e instead of %.3e", i2, minDists[i2], closestDist);
        }
        else if (minDistAnchors[i2] != closestAnchor)
        {
            quit = true;
            LOG(LOG_LVL_ERROR, "@[%ld] minDistAnchors=%ld instead of %ld", i2, minDistAnchors[i2], closestAnchor);
        }
        // else
        //     LOG(LOG_LVL_INFO, "@[%ld] all ok", i2);
    }
    
    return quit;
}
#endif
#endif