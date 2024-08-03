#include "parallhull.h"

#define MAX_THREADS 256

#ifdef NON_MPI_MODE
    #include <time.h>
    #include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#else
    #include <mpi.h>
#endif
#include <math.h>
#include <pthread.h>

typedef struct {
    Data fullData;
    Data *hulls;
    size_t reducedProblemUB; // upper bound on the size of the problem on which quickhull will run.
    int nThreads;
    int threadID;
    int *finishRecord;
    size_t *dataSize;
} ThreadData;

static void *parallhullThread(void *arg);
static Data mergeHulls(Data *h1, Data *h2);
static inline bool findNextMergePoint(Data *mergedH, Data *mainH, Data *altH, size_t *mainHIndex, size_t *altHIndex);

Data parallhullThreaded(Data *d, size_t reducedProblemUB, int nThreads)
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
        ds[i].threadID = i;
        ds[i].finishRecord = finishRecord;
        ds[i].dataSize = dataSize;
        pthread_create(&threads[i], NULL, parallhullThread, (void*)&ds[i]);
    }

    pthread_join(threads[0], NULL); // join outside so that final hull memory alloc is not free
    for (int i = 1; i < nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        if (hulls[i].X - d->X > d->n)
            free(hulls[i].X);
    }
    
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
    if (rd.n > thData->reducedProblemUB)
    {
        size_t nParts = (size_t)ceil((double)rd.n / thData->reducedProblemUB);
        Data *hulls = malloc(nParts * sizeof(Data));
        if (hulls == NULL)
            throwError("Could not allocate %ld elements of 8 bytes in memory on thread %d", nParts * sizeof(Data), thData->threadID);

        size_t avgPartSize = (size_t)(ceil((double)rd.n / nParts));

        // P1.1: sequentially compute quickhull on each partition generated using the specified upper bound on the size of the rrd(Reduced Reduced problem Data)
        for (size_t i = 0; i < nParts; i++)
        {
            size_t startPos = avgPartSize * i;

            hulls[i].n = avgPartSize;
            hulls[i].X = &rd.X[startPos];
            hulls[i].Y = &rd.Y[startPos];

            if (i == nParts-1)
                hulls[i].n = rd.n - avgPartSize * (nParts-1);

            hulls[i].n = quickhull(&hulls[i], thData->threadID);
        }

        // P1.2: sequentially merge the convex hulls generated in every rrd
        while (nParts > 1)
        {
            for (size_t i = 0; i+1 < nParts; i += 2)
            {
                Data h = mergeHulls(&hulls[i], &hulls[i+1]);

                if (hulls[i].X - rd.X > rd.n) // free memory allocated by mergehulls
                    free(hulls[i].X);
                if (hulls[i+1].X - rd.X > rd.n)
                    free(hulls[i+1].X);
                
                hulls[i/2] = h;
            }
            if (nParts % 2 == 1)
            {
                hulls[nParts/2 + 1] = hulls[nParts-1];
                nParts++;
            }
            nParts /= 2;
        }
        
        thData->hulls[thData->threadID] = hulls[0];
        
        free(hulls);
    }
    else
    {
        thData->hulls[thData->threadID].X = rd.X;
        thData->hulls[thData->threadID].Y = rd.Y;
        thData->hulls[thData->threadID].n = quickhull(&rd, thData->threadID);
    }

    thData->finishRecord[thData->threadID] = 1;

    // P2: thread merge their results with each other in a ordered manner
    int threadIDmod = thData->threadID;
    int secThreadIndexOffset = 1;
    while (threadIDmod % 2 == 0)
    {
        while (thData->finishRecord[thData->threadID + secThreadIndexOffset] < secThreadIndexOffset) // spinlock
            __builtin_ia32_pause();
        
        Data h = mergeHulls(&thData->hulls[thData->threadID], &thData->hulls[thData->threadID + secThreadIndexOffset]);

        if (thData->hulls[thData->threadID].X - thData->fullData.X > thData->fullData.n)
            free(thData->hulls[thData->threadID].X);
        
        thData->hulls[thData->threadID] = h;

        secThreadIndexOffset *= 2;
        threadIDmod /= 2;
    }

    thData->finishRecord[thData->threadID] = secThreadIndexOffset;
    
    return NULL;
}

static Data mergeHulls(Data *h1, Data *h2)
{
    size_t h1Index = 0, h2Index = 0;

    Data h0;
    h0.n = 0;
    h0.X = malloc((h1->n + h2->n) * 2 * sizeof(float));
    if (h0.X == NULL)
        throwError("Failed to allocate memory for merged hull");

    bool mainHullIsH1;

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
    }

    //################## Merge Procedure ##########################

    while ((h1Index <= h1->n) || (h2Index <= h2->n))
    {
        if (mainHullIsH1)
            mainHullIsH1 = findNextMergePoint(&h0, h1, h2, &h1Index, &h2Index);
        else
            mainHullIsH1 = findNextMergePoint(&h0, h2, h1, &h2Index, &h1Index);

        #ifdef DEBUG
            if (mainHullIsH1 && (h1Index >= h1->n+1))
                throwError("mergeHull: Cannot add points from an hull that has already been read completely. h1Index=%ld, h1Size=%ld, h2Index=%ld, h2Size=%ld, h0Size=%ld", h1Index, h1->n, h2Index, h2->n, h0.n);
            if (!mainHullIsH1 && (h2Index >= h2->n+1))
                throwError("mergeHull: Cannot add points from an hull that has already been read completely. h1Index=%ld, h1Size=%ld, h2Index=%ld, h2Size=%ld, h0Size=%ld", h1Index, h1->n, h2Index, h2->n, h0.n);
        #endif

        if (mainHullIsH1)
        {
            h0.X[h0.n] = h1->X[h1Index];
            h0.Y[h0.n] = h1->Y[h1Index];
            h1Index++;
        }
        else
        {
            h0.X[h0.n] = h2->X[h2Index];
            h0.Y[h0.n] = h2->Y[h2Index];
            h2Index++;
        }
        h0.n++;
    }
    
    return h0;
}

static inline bool findNextMergePoint(Data *mergedH, Data *mainH, Data *altH, size_t *mainHIndex, size_t *altHIndex)
{
    #ifdef DEBUG
        if (*mainHIndex > mainH->n+1)
            throwError("findNextMergePoint: mainHIndex > mainH.n");
        if (*altHIndex > altH->n+1)
            throwError("findNextMergePoint: altHIndex > altH.n");
        if ((*mainHIndex == mainH->n) && (*altHIndex == altH->n))
            throwError("findNextMergePoint: mainHIndex == mainH->n && altHIndex == altH->n");
    #endif

    bool swapMainToAlt=false;

    size_t lastIndex = mergedH->n - 1;

    {
        double a,b,c;
        a = (double)mainH->X[*mainHIndex] - mergedH->X[lastIndex];
        b = (double)mergedH->Y[lastIndex] - mainH->Y[*mainHIndex];
        c = (double)mergedH->X[lastIndex] * mainH->Y[*mainHIndex] - (double)mainH->X[*mainHIndex] * mergedH->Y[lastIndex];
        
        // now find min dist point starting from 
        double previousDist = INFINITY;
        double currentDist = INFINITY;

        while ((previousDist >= currentDist) && (*altHIndex < altH->n)) // while distance is decreasing and we are not looking through already looked data
        {
            size_t altHIndexP1 = *altHIndex + 1;
            if (altHIndexP1 == altH->n)
                altHIndexP1 = 0;

            previousDist = currentDist;
            currentDist = a * altH->Y[altHIndexP1] + b * altH->X[altHIndexP1] + c;
            if (currentDist < 0) // candidate point is not in the final hull (because it is inside of it)
            {
                swapMainToAlt = true;
                break;
            }
            (*altHIndex)++;
        }
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
            if (altHIndexP1 == altH->n)
                altHIndexP1 = 0;
            double dist = a * altH->Y[altHIndexP1] + b * altH->X[altHIndexP1] + c;

            if (dist > 0) break; // ignore dist == 0 since it is on the same line just get the point farther
            
            (*altHIndex)++;
        }
    }

    return swapMainToAlt;
}

