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

    Data hull = parallhullThreaded(&d, -1L, p.nThreads);
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    quickhullTime = cvtTimespec2Double(timeStruct);
    LOG(LOG_LVL_NOTICE, "Quickhull finished in %lfs", quickhullTime - fileReadTime);

    #ifdef GUI_OUTPUT
        if (d.n < 200000)
            plotData(&d, &hull, 0, GNUPLOT_RES, "Complete Hull");
    #endif

    free(d.X);
    d.X = NULL; d.Y = NULL;

    return EXIT_SUCCESS;
}

#else

int main (int argc, char *argv[])
{
    double startTime, initTime, fileReadTime, quickhullTime;
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

    if (rank == 0)
    {
        initTime = MPI_Wtime();
        LOG(LOG_LVL_NOTICE, "MPI run with:\n\tnProcs = %d\n\tnThreads = %d\n", p.nProcs, p.nThreads);
        LOG(LOG_LVL_NOTICE, "MPI init took %lfs", initTime - startTime);
    }

    if (rank == 0)
    {
        readFile(&d, &p);
        fileReadTime = MPI_Wtime();

        LOG(LOG_LVL_DEBUG, "Check endianity of raw file content: X[0]=%f  X[1]=%f", d.X[0], d.X[1]);
        LOG(LOG_LVL_NOTICE, "File read in %lfs", fileReadTime - startTime);
    }

    MPIErrCode = MPI_Bcast(&d.n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (MPIErrCode != MPI_SUCCESS)
        throwError("MPI_Bcast failed with code %d", MPIErrCode);
    
    size_t n = d.n;
    size_t reducedSize = (size_t)ceil((double)n / p.nProcs);

    if (rank == 0)
        LOG(LOG_LVL_NOTICE, "Quickhull finished in %lfs", quickhullTime - fileReadTime);

    #if defined(GUI_OUTPUT)
        if ((rank == 0) && (d.n < 200000))
            plotData(&d, &hull, 0, GNUPLOT_RES, "Complete Hull");
    #endif

    free(d.X);
    d.X = NULL; d.Y = NULL;

    return EXIT_SUCCESS;
}

#endif