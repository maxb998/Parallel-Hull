#include "parallhull.h"

#include <stdio.h>

#ifdef NON_MPI_MODE
    #include <time.h>
    #include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#else
    #include <mpi.h>
    #include <math.h>
#endif


#ifdef NON_MPI_MODE
int main (int argc, char *argv[])
{
    double startTime, fileReadTime, quickhullTime;
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    startTime = cvtTimespec2Double(timeStruct);

    Data d = {
        .n=0,
        .X=NULL,
        .Y=NULL
    };
    Params p = argParse(argc, argv);

    readFile(&d, &p);
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    fileReadTime = cvtTimespec2Double(timeStruct);
    LOG(LOG_LVL_DEBUG, "Check endianity of raw file content: X[0]=%f  X[1]=%f", d.X[0], d.X[1]);
    LOG(LOG_LVL_NOTICE, "File read in %lfs", fileReadTime - startTime);

    Data hull = parallhullThreaded(&d, -1, 0, p.nThreads);
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    quickhullTime = cvtTimespec2Double(timeStruct);

    LOG(LOG_LVL_NOTICE, "Parallhull finished in %lfs", quickhullTime - fileReadTime);
    LOG(LOG_LVL_INFO, "Final Hull size = %ld", hull.n);

    #ifdef DEBUG
        // check hull for duplicates before common errors
        for (size_t i = 0; i < hull.n; i++)
            for (size_t j = i+1; j < hull.n; j++)
                if ((hull.X[i] == hull.X[j]) && (hull.Y[i] == hull.Y[j]))
                    LOG(LOG_LVL_ERROR, "Final Hull contains a duplicate at position %ld and %ld", i , j);
        

        ProcThreadIDCombo fakeID = { .p=0, .t=0 };
        if (hullConvexityCheck(&hull, &fakeID))
            throwError("Final Hull is not convex");
        free(d.X);
        d.n = 0; d.X = NULL; d.Y = NULL;
        readFile(&d, &p);
        if (finalCoverageCheck(&hull, &d, &fakeID))
            throwError("Final Hull does not cover all points");
    #endif

    #ifdef GUI_OUTPUT
        if (d.n < 200000)
            plotData(&d, &hull, 0, "Complete Hull");
    #endif

    // saveHullPointsTxt(&hull, "hullPts.txt");

    free(d.X);
    d.X = NULL; d.Y = NULL;

    return EXIT_SUCCESS;
}

#else

int main (int argc, char *argv[])
{
    double startTime, initTime, fileReadTime, localHullTime, mergeTime;
    int MPIErrCode;
    startTime = MPI_Wtime();

    Data d = {
        .n=0,
        .X=NULL,
        .Y=NULL
    };
    Params p = argParse(argc, argv);

    int rank;
    MPIErrCode = MPI_Init(&argc, &argv);
    if (MPIErrCode != MPI_SUCCESS)
        throwError("MPI_init failed with code %d", MPIErrCode); 
    MPIErrCode = MPI_Comm_size(MPI_COMM_WORLD, &p.nProcs);
    if (MPIErrCode !=MPI_SUCCESS)
        throwError("MPI_Comm_size failed with code %d", MPIErrCode);
    MPIErrCode = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (MPIErrCode !=MPI_SUCCESS)
        throwError("MPI_Comm_rank failed with code %d", MPIErrCode);

    initTime = MPI_Wtime();
    LOG(LOG_LVL_NOTICE, "p[%d] MPI run with: nProcs = %2d \tnThreads = %3d\n", rank, p.nProcs, p.nThreads);
    LOG(LOG_LVL_NOTICE, "p[%d] MPI init took %lfs", rank, initTime - startTime);

    readFilePart(&d, &p, rank);

    if (rank == 0)
        LOG(LOG_LVL_DEBUG, "Check endianity of raw file content: X[0]=%f  X[1]=%f", d.X[0], d.X[1]);

    fileReadTime = MPI_Wtime();
    LOG(LOG_LVL_NOTICE, "p[%d] File read in %lfs", rank, fileReadTime - startTime);

    Data hull = parallhullThreaded(&d, -1, rank, p.nThreads);

    localHullTime = MPI_Wtime();
    LOG(LOG_LVL_NOTICE, "p[%d] Local quickhull finished in %lfs", rank, localHullTime - fileReadTime);
    
    mpiHullMerge(&hull, rank, p.nProcs);

    mergeTime = MPI_Wtime();
    LOG(LOG_LVL_NOTICE, "p[%d] Hull merge finished in %lfs", rank, mergeTime - localHullTime);

    #ifdef DEBUG
        if (rank == 0)
        {
            Data fullData;
            readFile(&fullData, &p);

            ProcThreadIDCombo id = { .p=0, .t=0 };

            if (finalCoverageCheck(&hull, &fullData, &id))
                throwError("final coverage check failed");
            LOG(LOG_LVL_NOTICE, "Final checks ok!");

            free(fullData.X);
        }
    #endif

    #if defined(GUI_OUTPUT)
        if ((rank == 0) && (d.n < 200000))
            plotData(&d, &hull, 0, "Complete Hull");
    #endif

    if (rank == 0)
    {
        LOG(LOG_LVL_NOTICE, "Total time taken: %lfs", mergeTime - startTime);
        LOG(LOG_LVL_NOTICE, "Total time taken(without init): %lfs", mergeTime - initTime);
        LOG(LOG_LVL_NOTICE, "Computation time taken: %lfs", mergeTime - fileReadTime);
    }

    free(d.X);
    free(hull.X);
    free(hull.Y);

    MPI_Finalize();
    return EXIT_SUCCESS;
}

#endif