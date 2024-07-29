#include "paralhull.h"

#include <stdio.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


int main (int argc, char *argv[])
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    Data d = {
        .n=0,
        .X=NULL,
        .Y=NULL
    };
    Params p = argParse(argc, argv);

    readFile(&d, &p);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double fileReadTime = cvtTimespec2Double(timeStruct);
    LOG(LOG_LVL_DEBUG, "Check endianity of raw file content: X[0]=%f  X[1]=%f", d.X[0], d.X[1]);
    LOG(LOG_LVL_NOTICE, "File read in %lfs", fileReadTime - startTime);

    int hullSize = quickhull(&d);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double quickhullTime = cvtTimespec2Double(timeStruct); 
    LOG(LOG_LVL_NOTICE, "Quickhull finished in %lfs", quickhullTime - fileReadTime);

    #ifdef GUI_OUTPUT
        if (d.n < 10000)
            plotData(&d, hullSize, 0, GNUPLOT_RES, "Complete Hull");
    #endif

    free(d.X);
    d.X = NULL; d.Y = NULL;

    return EXIT_SUCCESS;
}
