#include "paralhull.h"

#include <math.h>
#include <stdio.h>

typedef struct 
{
    int node;
    int anchor;
} SuccessorData;

static inline void getExtremeCoordsPts(Data *d, int ptIndices[4]);
static inline int setExtremeCoordsAsInit(Data *d, int ptIndices[4]);
static void removeConveredPoints(Data *d, int hullSize, int *nUncovered);
static inline SuccessorData findFarthestPt(Data *d, int hullSize, int nUncovered);
static inline void addPtToHull(Data *d, SuccessorData pt, int hullSize);

int quickhull(Data *d)
{
    // init
    int ptIndices[4];
    getExtremeCoordsPts(d, ptIndices);
    int hullSize = setExtremeCoordsAsInit(d, ptIndices);
    
    int nUncovered = d->n - hullSize;
    removeConveredPoints(d, hullSize, &nUncovered);

    while (nUncovered > 0)
    {
        plotPartialData(d, hullSize, nUncovered, "1600,900");
        getchar();

        SuccessorData s = findFarthestPt(d, hullSize, nUncovered);
        addPtToHull(d, s, hullSize);
        hullSize++;

        removeConveredPoints(d, hullSize, &nUncovered);
    }

    return hullSize;
}

static inline void getExtremeCoordsPts(Data *d, int ptIndices[4])
{
    for (int i = 0; i < 4; i++)
        ptIndices[i] = 0;
    
    for (int i = 1; i < d->n; i++)
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

static inline int setExtremeCoordsAsInit(Data *d, int ptIndices[4])
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
    
    int hullSize = 0;
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

static void removeConveredPoints(Data *d, int hullSize, int *nUncovered)
{
    int i = hullSize;
    int j = *nUncovered + hullSize - 1; // position of last covered element in d->X and d->Y
    while (i < j)
    {
        bool iIsCovered = true;
        for (int k1 = 0, k2 = hullSize-1; k1 < hullSize; k2 = k1++)
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

static inline SuccessorData findFarthestPt(Data *d, int hullSize, int nUncovered)
{
    SuccessorData farthestPt = { .node=0, .anchor=-1 };
    float farthestPtDist = 0;

    for (int i = hullSize; i < hullSize + nUncovered; i++)
    {
        int closestAnchor = 0;
        float closestDist = INFINITY;
        for (int j = 0, k = hullSize-1; j < hullSize; k=j++)
        {
            register float xDiff1 = d->X[k] - d->X[i];
            register float xDiff2 = d->X[i] - d->X[j];
            register float yDiff1 = d->Y[k] - d->Y[i];
            register float yDiff2 = d->Y[i] - d->Y[j];
            register float dist1 = sqrtf(xDiff1 * xDiff1 + yDiff1 * yDiff1);
            register float dist2 = sqrtf(xDiff2 * xDiff2 + yDiff2 * yDiff2);
            register float dist = dist1 + dist2;
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
        throwError("Invalid anchor and node");

    return farthestPt;
}

static inline void addPtToHull(Data *d, SuccessorData pt, int hullSize)
{
    swapElems(d->X[pt.node], d->X[hullSize]);
    swapElems(d->Y[pt.node], d->Y[hullSize]);

    float xBak = d->X[hullSize], yBak = d->Y[hullSize];
    for (int i = hullSize-1; i > pt.anchor; i--)
    {
        d->X[i+1] = d->X[i];
        d->Y[i+1] = d->Y[i];
    }

    d->X[pt.anchor+1] = xBak;
    d->Y[pt.anchor+1] = yBak;
}
