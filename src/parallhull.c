#include "parallhull.h"

#include <math.h>
#include <pthread.h>

#ifdef NON_MPI_MODE
    #include <time.h>
    #include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#else
    #include <mpi.h>
#endif
#if defined(PARALLHULL_STEP_DEBUG) || defined(PARALLHULL_MERGE_OUTPUT_PLOT)
    #include <stdio.h>
#endif

#define MAX_THREADS 256

typedef enum {
    NPO_CONTINUE,
    NPO_SWAP
} NextPtOp;


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
static inline NextPtOp findNextMergePoint(Data *mergedH, Data *mainH, Data *altH, size_t *mainHIndex, size_t *altHIndex, ProcThreadIDCombo *id);
#ifdef DEBUG
    static inline bool mergeHullCoverageCheck(Data *h0, Data *h1, Data *h2, ProcThreadIDCombo *id);
#endif

Data parallhullThreaded(Data *d, size_t reducedProblemUB, int procID, int nThreads)
{
    #ifdef NON_MPI_MODE
        struct timespec timeStruct;
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        double startTime = cvtTimespec2Double(timeStruct);
    #else
        double startTime = MPI_Wtime();
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

    pthread_join(threads[0], NULL); // join outside so that final hull memory alloc is not free
    for (int i = 1; i < nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        free(hulls[i].X);
        free(hulls[i].Y);
    }

    #ifdef DEBUG
    #endif
    
    #ifdef NON_MPI_MODE
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        double finishTime = cvtTimespec2Double(timeStruct);
        LOG(LOG_LVL_NOTICE, "Exec time of ParallHull using threads is %lfs", finishTime - startTime);
    #else
        double finishTime = MPI_Wtime();
        LOG(LOG_LVL_NOTICE, "p[%2d] parallhull: finished hull computation in %lfs", procID, finishTime - startTime);
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

            Data pts = { .X=&rd.X[startPos], .Y=&rd.Y[startPos], .n=avgPartSize };

            if (i == nParts-1)
                pts.n = rd.n - avgPartSize * (nParts-1);

            hulls[i] = quickhull(&pts, &thData->id);
        }

        LOG(LOG_LVL_INFO, "p[%2d] t[%3d] parallhullThread: Quickhull on subproblem/s done, now merging", thData->id.p, thID);

        // P1.2: sequentially merge the convex hulls generated in every rrd
        while (nParts > 1)
        {
            size_t halfNParts = nParts / 2;
            for (size_t i = 0; i < halfNParts; i++)
            {
                LOG(LOG_LVL_TRACE, "p[%2d] t[%3d] parallhullThread: Merging thread internal hulls %ld(size=%ld) and %ld(size=%ld)", thData->id.p, thID, i, hulls[i].n, i + halfNParts, hulls[i+halfNParts].n);
                Data h = mergeHulls(&hulls[i], &hulls[i+halfNParts], &thData->id);

                #ifdef DEBUG
                    if (hullConvexityCheck(&h, &thData->id))
                    {
                        plotHullMergeStep(&hulls[i], &hulls[i+halfNParts], &h, 0, 0, "Plot of the error", false);
                        throwError("p[%2d] t[%3d] parallhullThread: Merged hull is not convex", thData->id.p, thID);
                    }
                    if (mergeHullCoverageCheck(&h, &hulls[i], &hulls[i+halfNParts], &thData->id))
                    {
                        plotHullMergeStep(&hulls[i], &hulls[i+halfNParts], &h, 0, 0, "Plot of the error", false);
                        throwError("p[%2d] t[%3d] parallhullThread: Merged Hull does not cover all the points in the hull", thData->id.p, thID);
                    }
                #endif

                free(hulls[i].X);
                free(hulls[i].Y);
                free(hulls[i+halfNParts].X);
                free(hulls[i+halfNParts].Y);
                
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
        thData->hulls[thID] = quickhull(&rd, &thData->id);

    LOG(LOG_LVL_INFO, "p[%2d] t[%3d] parallhullThread: Thread subproblem solved", thData->id.p, thID);

    thData->finishRecord[thID] = 1;

    // P2: thread merge their results with each other in a ordered manner
    int s = 0;
    int thID2merge = thID + 1;
    while ((((thID>>s) & 1) == 0) && (thID2merge < thData->nThreads))
    {
        while (thData->finishRecord[thID2merge] < s+1) // spinlock (the assumption here is that threads should take more or less the same amount of time to merge, and this kinds of keeps the cpu "warm")
            __builtin_ia32_pause();
        
        Data h = mergeHulls(&thData->hulls[thID], &thData->hulls[thID2merge], &thData->id);

        LOG(LOG_LVL_INFO, "p[%2d] t[%3d] parallhullThread: Merging hull with hull in thread %d. s=%d", thData->id.p, thID, thID2merge, s);

        #ifdef DEBUG
            if (hullConvexityCheck(&h, &thData->id))
            {
                plotHullMergeStep(&thData->hulls[thID], &thData->hulls[thID2merge], &h, 0, 0, "Plot of the error", false);
                throwError("p[%2d] t[%3d] parallhullThread: Merged hull is not convex", thData->id.p, thID);
            }
            if (mergeHullCoverageCheck(&h, &thData->hulls[thID], &thData->hulls[thID2merge], &thData->id))
            {
                plotHullMergeStep(&thData->hulls[thID], &thData->hulls[thID2merge], &h, 0, 0, "Plot of the error", false);
                throwError("p[%2d] t[%3d] parallhullThread: Merged Hull does not cover all the points in the hull", thData->id.p, thID);
            }
        #endif

        // each thread manages to free its own memory
        free(thData->hulls[thID].X);
        free(thData->hulls[thID].Y);
        
        thData->hulls[thID] = h;

        s++;
        thID2merge = thID + (1<<s);

        thData->finishRecord[thID] = s+1;
    }

    thData->finishRecord[thID] = 0x7FFFFFFF; // cannot stall spinlock anymore
    
    return NULL;
}

#ifndef NON_MPI_MODE
void mpiHullMerge(Data *h1, int rank, int nProcs)
{
    int s = 0;
    int rank2receive = rank + 1;
    int MPIErrCode;

    int modRank = rank;
    while ((((modRank>>s) & 1) == 0) && (rank2receive < nProcs))
    {
        Data h2;
        { // receive
            LOG(LOG_LVL_DEBUG, "p[%2d] mpiHullMerge: receiving data from rank %d", rank, rank2receive);

            MPIErrCode = MPI_Recv(&h2.n, 1, MPI_UNSIGNED_LONG, rank2receive, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (MPIErrCode)
                throwError("p[%2d] mpiHullMerge: Got error %d on receiving the number of elements of partial hull from p[%d]", rank, MPIErrCode, rank2receive);

            h2.X = malloc((h2.n + 1) * sizeof(float) + MALLOC_PADDING);
            h2.Y = malloc((h2.n + 1) * sizeof(float) + MALLOC_PADDING);
            if ((h2.X == NULL) || (h2.Y == NULL))
                throwError("p[%2d] mpiHullMerge: Failed to allocate memory for the hull to be received");

            MPIErrCode = MPI_Recv(h2.X, h2.n, MPI_FLOAT, rank2receive, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (MPIErrCode)
                throwError("p[%2d] mpiHullMerge: Got error %d on receiving the partial hull Xs from p[%d]", rank, MPIErrCode, rank2receive);
            h2.X[h2.n] = h2.X[0];
            h2.X[h2.n+1] = h2.X[1];

            MPIErrCode = MPI_Recv(h2.Y, h2.n, MPI_FLOAT, rank2receive, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (MPIErrCode)
                throwError("p[%2d] mpiHullMerge: Got error %d on receiving the partial hull Xs from p[%d]", rank, MPIErrCode, rank2receive);
            h2.Y[h2.n] = h2.Y[0];
            h2.Y[h2.n+1] = h2.Y[1];
        }
        
        ProcThreadIDCombo id = {.p=rank, .t=0};
        Data h0 = mergeHulls(h1, &h2, &id);

        LOG(LOG_LVL_INFO, "p[%2d] mpiHullMerge: Merging hull with hull in proc %d", rank, rank2receive);

        #ifdef DEBUG
            if (hullConvexityCheck(&h0, &id))
            {
                plotHullMergeStep(h1, &h2, &h0, 0, 0, "Plot of the error", false);
                throwError("p[%2d] mpiHullMerge: Merged hull is not convex", rank);
            }
            if (mergeHullCoverageCheck(&h0, h1, &h2, &id))
            {
                plotHullMergeStep(h1, &h2, &h0, 0, 0, "Plot of the error", false);
                throwError("p[%2d] mpiHullMerge: Merged Hull does not cover all the points in the hull", rank);
            }
        #endif

        // each thread manages to free its own memory
        free(h1->X);
        free(h1->Y);
        free(h2.X);
        free(h2.Y);

        *h1 = h0;

        s++;
        rank2receive = rank + (1<<s);
    }
    if (rank != 0) // send
    {
        int rank2send = rank & (0xFFFFFFFE<<s);
        while (rank2send == rank)
        {
            s++;
            rank2send = rank & (0xFFFFFFFE<<s);
        }
        LOG(LOG_LVL_DEBUG, "p[%2d] mpiHullMerge: sending data to rank %d with s=%d", rank, rank2send, s);

        MPIErrCode = MPI_Send(&h1->n, 1, MPI_UNSIGNED_LONG, rank2send, 0, MPI_COMM_WORLD);
        if (MPIErrCode)
            throwError("p[%2d] mpiHullMerge: Got error %d on sending the number of elements of partial hull from p[%d]", rank, MPIErrCode, rank2send);

        MPIErrCode = MPI_Send(h1->X, h1->n, MPI_FLOAT, rank2send, 1, MPI_COMM_WORLD);
        if (MPIErrCode)
            throwError("p[%2d] mpiHullMerge: Got error %d on sending the partial hull Xs from p[%d]", rank, MPIErrCode, rank2send);

        MPIErrCode = MPI_Send(h1->Y, h1->n, MPI_FLOAT, rank2send, 2, MPI_COMM_WORLD);
        if (MPIErrCode)
            throwError("p[%2d] mpiHullMerge: Got error %d on sending the partial hull Xs from p[%d]", rank, MPIErrCode, rank2send);
    }
}
#endif

static Data mergeHulls(Data *h1, Data *h2, ProcThreadIDCombo *id)
{
    Data h0;
    h0.n = 0;
    h0.X = malloc((h1->n + h2->n) * sizeof(float) + MALLOC_PADDING);
    h0.Y = malloc((h1->n + h2->n) * sizeof(float) + MALLOC_PADDING);
    if ((h0.X == NULL) || (h0.Y == NULL))
        throwError("p[%2d] t[%3d] mergeHulls: Failed to allocate memory for merged hull", id->p, id->t);

    Data *mainH, *altH;
    size_t mainHIndex = 0, altHindex = 0;

    // set an additional point at the end of h1 and h2 to avoid errors(memory space available in padding)
    h1->X[h1->n+1] = h1->X[1];
    h1->Y[h1->n+1] = h1->Y[1];
    h2->X[h2->n+1] = h2->X[1];
    h2->Y[h2->n+1] = h2->Y[1];
    h0.X[0] = NAN; h0.Y[0] = NAN; // needed for first loop iter

    {
    if (h1->Y[0] < h2->Y[0])
    {
        mainH = h1;
        altH = h2;
    }
    else if (h1->Y[0] > h2->Y[0])
    {
        mainH = h2;
        altH = h1;
    }
    else if (h1->X[0] >= h2->X[0]) // last two cases are here in case of equality (look at the X coordinate and choose the rightmost) (likely a rare occurence)
    {
        mainH = h1;
        altH = h2;
    }
    else
    {
        mainH = h2;
        altH = h1;
    }
    }

    while ((mainHIndex < mainH->n) || ((mainHIndex <= mainH->n) && ((h0.X[h0.n-1] != h0.X[0]) || (h0.Y[h0.n-1] != h0.Y[0]))))
    {
        #ifdef DEBUG
            if (mainHIndex > mainH->n)
            {
                plotHullMergeStep(mainH, altH, &h0, mainHIndex, altHindex, "Plot of the error", false);
                throwError("p[%2d] t[%3d] mergeHull: Cannot add points from an hull that has already been read completely. mainHIndex=%ld, mainH.n=%ld, altHIndex=%ld, altH.n=%ld, h0.n=%ld", id->p, id->t, mainHIndex, mainH->n, altHindex, altH->n, h0.n);
            }
        #endif

        #ifdef PARALLHULL_STEP_DEBUG
            // show partial hull with gnuplot at each iteration and wait for user input to resume
            char plotTitle[200];
            sprintf(plotTitle, "Merged Hull: h0.n=%lu, mainH.n=%ld, altH.n=%ld, mainHIndex=%ld, altHIndex=%ld", h0.n, mainH->n, altH->n, mainHIndex, altHindex);
            LOG(LOG_LVL_DEBUG, "%s", plotTitle);
            plotHullMergeStep(mainH, altH, &h0, mainHIndex, altHindex, plotTitle, false);
            getchar();
        #endif

        // add pt
        h0.X[h0.n] = mainH->X[mainHIndex];
        h0.Y[h0.n] = mainH->Y[mainHIndex];
        h0.n++;
        mainHIndex++;

        NextPtOp npo = findNextMergePoint(&h0, mainH, altH, &mainHIndex, &altHindex, id);

        if (npo == NPO_SWAP)
        {
            swapElems(mainH, altH)
            swapElems(mainHIndex, altHindex)
        }
    }
    h0.n--;

    h0.X = realloc(h0.X, (h0.n + 1) * sizeof(float) + MALLOC_PADDING);
    h0.Y = realloc(h0.Y, (h0.n + 1) * sizeof(float) + MALLOC_PADDING);

    return h0;
}

static inline NextPtOp findNextMergePoint(Data *mergedH, Data *mainH, Data *altH, size_t *mainHIndex, size_t *altHIndex, ProcThreadIDCombo *id)
{
    #ifdef DEBUG
        if (*mainHIndex > mainH->n+1)
            throwError("p[%2d] t[%3d] findNextMergePoint: mainHIndex > mainH.n", id->p, id->t);
        if (*altHIndex > altH->n+1)
            throwError("p[%2d] t[%3d] findNextMergePoint: altHIndex > altH.n", id->p, id->t);
    #endif

    if (*altHIndex > altH->n)
        return NPO_CONTINUE;

    NextPtOp retval = NPO_CONTINUE;
    size_t lastIndex = mergedH->n - 1;
    {
        double a,b,c;
        a = (double)mainH->X[*mainHIndex] - mergedH->X[lastIndex];
        b = (double)mergedH->Y[lastIndex] - mainH->Y[*mainHIndex];
        c = (double)mergedH->X[lastIndex] * mainH->Y[*mainHIndex] - (double)mainH->X[*mainHIndex] * mergedH->Y[lastIndex];
        
        // now find min dist point starting from 
        double previousDist = INFINITY;
        double currentDist = a * altH->Y[*altHIndex] + b * altH->X[*altHIndex] + c;

        while ((previousDist >= currentDist) && (*altHIndex < altH->n) && (currentDist >= 0)) // while distance is decreasing and we are not looking through already looked data
        {
            (*altHIndex)++;
            previousDist = currentDist;
            currentDist = a * altH->Y[*altHIndex] + b * altH->X[*altHIndex] + c;
        }

        if (currentDist < 0)
            retval = NPO_SWAP;
        else // (previousDist >= currentDist) || (*altHIndex <= altH->n)
            (*altHIndex)--;
    }

    // need to check if point in altH[altHIndex] is the correct one for the chosen hull or not. It might not belong to the convex hull and be covered by the line casted a successive point in the hull
    if (retval == NPO_SWAP)
    {
        while (*altHIndex < altH->n)
        {
            double a,b,c;
            a = (double)altH->X[*altHIndex] - mergedH->X[lastIndex];
            b = (double)mergedH->Y[lastIndex] - altH->Y[*altHIndex];
            c = (double)mergedH->X[lastIndex] * altH->Y[*altHIndex] - (double)altH->X[*altHIndex] * mergedH->Y[lastIndex];

            size_t altHIndexP1 = *altHIndex + 1;
            double dist = a * altH->Y[altHIndexP1] + b * altH->X[altHIndexP1] + c;

            if (dist > 0) break; // ignore dist == 0 since it is on the same line just get the point farther
            
            (*altHIndex) = altHIndexP1;
        }
    }

    return retval;
}

#ifdef DEBUG
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
            if (((h1->X[j] == h0->X[i]) && (h1->Y[j] == h0->Y[i])) ||
                ((h1->X[j] == h0->X[ip1]) && (h1->Y[j] == h0->Y[ip1])))
                continue;
            double dist = a * h1->Y[j] + b * h1->X[j] + c;
            if (dist < 0)
            {
                LOG(LOG_LVL_ERROR, "p[%2d] t[%3d] mergeHullCoverageCheck: Merged Hull does not contain point[%ld] of h1. It is not to the right of line between pt[%ld] and pt[%ld]", id->p, id->t, j, i, ip1);
                retval = true;
            }
        }
        for (size_t j = 0; j < h2->n; j++)
        {
            if (((h2->X[j] == h0->X[i]) && (h2->Y[j] == h0->Y[i])) ||
                ((h2->X[j] == h0->X[ip1]) && (h2->Y[j] == h0->Y[ip1])))
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
