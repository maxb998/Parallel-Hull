#include "paralhull.h"

static inline void getExtremeCoordsPts(Data *d, int ptIndices[4]);
static inline int setExtremeCoordsAsInit(Data *d, int ptIndices[4]);
static void removeConveredPoints(Data *d, int hullSize, int *nUncovered);

void quickhull(Data *d)
{
    // init
    int ptIndices[4];
    getExtremeCoordsPts(d, ptIndices);
    int hullSize = setExtremeCoordsAsInit(d, ptIndices);
    
    

}

static inline void getExtremeCoordsPts(Data *d, int ptIndices[4])
{
    for (int i = 0; i < 4; i++)
        ptIndices[i] == 0;
    
    for (int i = 1; i < d->n; i++)
    {
        if (d->X[i] > d->X[ptIndices[0]]) //xMax
            ptIndices[0] = i;
        else if (d->X[i] < d->X[ptIndices[1]]) //xMin
            ptIndices[1] = i;

        if (d->Y[i] > d->Y[ptIndices[2]]) //yMax
            ptIndices[2] = i;
        else if (d->Y[i] < d->Y[ptIndices[3]]) //yMin
            ptIndices[3] = i;
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
            if (ptIndices[minSubIndex] > ptIndices[j])
                minSubIndex = j;

        swapElems(d->Y[minSubIndex], d->Y[ptIndices[minSubIndex]]);
        swapElems(d->Y[minSubIndex], d->Y[ptIndices[minSubIndex]]);
        hullSize++;
    }

    return hullSize;
}

static void removeConveredPoints(Data *d, int hullSize, int *nUncovered)
{
    int i = hullSize;
    int j = d->n - *nUncovered + hullSize - 1;
    while (i < j)
    {
        bool iIsCovered = false;
        for (int k1 = 0, k2 = hullSize; k1 < hullSize; k2 = k1++)
            if (d->Y[k1] < d->Y[i] != (d->Y[k2] > d->Y[i]) && (d->X[i] < (d->X[k2]-d->X[k1]) * (d->Y[i]-d->Y[k1]) / (d->Y[k2]-d->Y[k1]) + d->X[k1]))
                iIsCovered = !iIsCovered;
        
        if (iIsCovered)
        {
            swapElems(d->X[i], d->X[j]);
            swapElems(d->Y[i], d->Y[j]);
            *nUncovered += 1;
            j--;
        }
        else
            i++;
    }
}