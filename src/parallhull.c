#include "parallhull.h"

#define MAX_THREADS 256

#ifdef LOCAL_MODE
    #include <time.h>
    #include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#endif
#include <math.h>
#include <pthread.h>
#include <immintrin.h>

typedef struct {
    Data fullData;
    size_t reducedProbUB; // upper bound on the size of the problem on which quickhull will run.
    int nThreads;
    int threadID;
    bool *finishRecord;
    size_t *dataSize;
} ThreadData;

static void *parallhullThread(void *arg);
static size_t mergeHulls(Data d, size_t h1Size, size_t h2Size, size_t rd1ElemsCount);


size_t parallhullThreaded(Data *d, size_t reducedProbUB, int nThreads)
{
    #ifdef LOCAL_MODE
        struct timespec timeStruct;
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        double startTime = cvtTimespec2Double(timeStruct);
    #endif

    bool finishRecord[MAX_THREADS];
    for (int i = 0; i < MAX_THREADS; i++)
        finishRecord[i] = false;

    // decide how many elements each thread gets(at the beginning)
    size_t dataSize[MAX_THREADS];
    double n = d->n;
    n = ceil(n/nThreads);
    for (int i = 0; i < nThreads-1; i++)
        dataSize[i] = (size_t)n;
    dataSize[nThreads-1] = d->n - ((size_t)n * (nThreads-1)); // last thread might have less elements due to non-divisibility of d->n with nThreads
    
    pthread_t threads[MAX_THREADS];
    ThreadData ds[MAX_THREADS];
    for (int i = 0; i < nThreads; i++)
    {
        ds[i].fullData = *d;
        ds[i].reducedProbUB = reducedProbUB;
        ds[i].nThreads = nThreads;
        ds[i].threadID = i;
        ds[i].finishRecord = finishRecord;
        ds[i].dataSize = dataSize;
        pthread_create(&threads[i], NULL, parallhullThread, (void*)&ds[i]);
    }
    
}


static void *parallhullThread(void *arg)
{
    ThreadData *thData = (ThreadData*)arg;
    Data rd = { .n=thData->dataSize[thData->threadID] };
    {
        Data fd = thData->fullData;
        size_t startElem = 0;
        for (int i = 0; i < thData->threadID; i++)
            startElem += thData->dataSize[i];
        rd.X = &fd.X[startElem];
        rd.Y = &fd.Y[startElem];
    }

    // P1: each thread works on its own data in the first part here
    size_t hullSize;
    if (rd.n > thData->reducedProbUB)
    {
        size_t nParts = (size_t)ceil((double)rd.n / thData->reducedProbUB);
        size_t *hullSizes = malloc(nParts * sizeof(size_t));
        if (hullSizes == NULL)
            throwError("Could not allocate %ld elements of 8 bytes in memory on thread %d", 2*nParts, thData->threadID);

        size_t avgPartSize = (size_t)(ceil((double)rd.n / nParts));

        // P1.1: sequentially compute quickhull on each partition generated using the specified upper bound on the size of the rrd(Reduced Reduced problem Data)
        for (size_t i = 0; i < nParts; i++)
        {
            size_t startPos = avgPartSize * i;
            Data rrd = { .n=avgPartSize, .X=&rd.X[startPos], .Y=&rd.Y[startPos]  };
            if (i == nParts-1)
                rrd.n = rd.n - avgPartSize * (nParts-1);
            hullSizes[i] = quickhull(&rrd);
        }

        // P1.2: sequentially merge the convex hulls generated in every rrd

        
        free(hullSizes);
    }
    else
        hullSize = quickhull(&rd);

    // P2: thread merge their results with each other in a ordered manner
}

static size_t mergeHulls(Data d, size_t h1Size, size_t h2Size, size_t rd1ElemsCount)
{
    // move all elements of hull2 close to hull1 without deleting any points, so swaps it is
    for (size_t i = 0; i < h2Size; i++)
    {
        swapElems(d.X[h1Size + i], d.X[rd1ElemsCount + i])
        swapElems(d.Y[h1Size + i], d.Y[rd1ElemsCount + i])
    } // rd1ElemsCount has now become useless

    // move hull that has the lowest point in position 1 if it is in position 2
    if (d.Y[0] > d.Y[h1Size])
    {
        float *h1BakX = malloc(h1Size * 2 * sizeof(float));
        if (h1BakX == NULL)
            throwError("Failed to allocate %ld bytes of memory for the X and Y backups swap", h1Size * 2 * sizeof(float));
        float *h1BakY = &h1BakX[h1Size];

        for (size_t i = 0; i < h1Size; i++)
        {
            h1BakX[i] = d.X[i];
            h1BakY[i] = d.Y[i];
        }
        for (size_t i = 0, j = h1Size; i < h2Size; i++, j++)
        {
            d.X[i] = d.X[j];
            d.Y[i] = d.Y[j];
        }
        for (size_t i = 0, j = h2Size; i < h1Size; i++, j++)
        {
            d.X[j] = h1BakX[i];
            d.Y[j] = h1BakY[i];
        }
        swapElems(h1Size, h2Size);
        free(h1BakX);   
    }

    //################## Actual merge procedure ##########################

    size_t h1Index = 1, h2Index = h1Size; // h1Index also tracks the merged hull build progress(size)
    for (size_t h2Index = h1Size; h2Index < h1Size + h2Size; h1Index++)
    {
        float a,b,c;
        size_t previousIndex = h1Index - 1;
        a = d.X[h1Index] - d.X[previousIndex];
        b = d.Y[previousIndex] - d.Y[h1Index];
        c = d.X[previousIndex] * d.Y[h1Index] - d.X[h1Index] * d.Y[previousIndex];
        
        // now find min dist point starting from 
        float previousDist = INFINITY;
        float currentDist = INFINITY;

        while (previousDist >= currentDist) // while distant is decreasing
        {
            currentDist = a * d.Y[h2Index] + b * d.X[h2Index] + c;
            if (currentDist < 0) // candidate point is not in the final hull
            {
                swapElems(d.X[h1Index], d.X[h2Index]);
                swapElems(d.Y[h1Index], d.Y[h2Index]);
                h2Index++;
                break;
            }
            previousDist = currentDist;
            h2Index++;
        }

        h1Index++;
    }
}

