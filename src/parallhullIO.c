#include "parallhull.h"

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#ifdef GUI_OUTPUT
const char * logLevelString [] = {
	"[\033[1;31mFATAL\033[0m]", // 0
	"[\033[0;35mERROR\033[0m]", // 1
	"[\033[0;33mWARN\033[0m]", // 2
	"[\033[0;36mNOTICE\033[0m]", // 3
	"[\033[0;34mINFO\033[0m]", // 4
	"[\033[0;32mDEBUG\033[0m]", // 5
	"[\033[0;90mTRACE\033[0m]"  // 6
};
#else
const char * logLevelString [] = {
	"[FATAL]", // 0
	"[ERROR]", // 1
	"[WARN]", // 2
	"[NOTICE]", // 3
	"[INFO]", // 4
	"[DEBUG]", // 5
	"[TRACE]"  // 6
};
#endif

static enum LogLevel LOG_LEVEL = LOG_LVL_INFO;
void setLogLevel(enum LogLevel lvl)
{
    LOG_LEVEL = lvl;
}

void LOG (enum LogLevel lvl, char * line, ...)
{
    // check log level
    if (lvl > LOG_LEVEL) return;

    // print log level
    printf("\r%s ", logLevelString[lvl]);
    fflush(stdout);

    // print passed message and values
    va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);
    
    // add new line at the end
    if (line[strlen(line)-1] != '\n')
        printf("\n");
}

void inlineLOG (enum LogLevel lvl, char * line, ...)
{
    // check log level
    if (lvl > LOG_LEVEL) return;

    // print log level
    printf("\r%s ", logLevelString[lvl]);
    fflush(stdout);

    // print passed message and values
    va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);
    
    fflush(stdout);
}

void throwError (char * line, ...)
{
    printf("%s ", logLevelString[0]);

    va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);

    printf("\n");

    exit(EXIT_FAILURE);
}

#ifdef NON_MPI_MODE
void readFile(Data *d, Params *p)
{
    FILE *fileptr = fopen(p->inputFile, "rb");
    if (fileptr == NULL)
        throwError("Could not read file %s", p->inputFile);

    fseek(fileptr, 0, SEEK_END);
    d->n = ftell(fileptr) / (2 * sizeof(float));
    rewind(fileptr);

    d->X = malloc(d->n * 2 * sizeof(float));
    if (d->X == NULL)
        throwError("Failed to allocate memory for points");
    d->Y = &d->X[d->n];

    fread(d->X, d->n * 2 * sizeof(float), 1, fileptr);

    fclose(fileptr);
}
#else
size_t readFile(Data *d, Params *p, int rank)
{
    FILE *fileptr = fopen(p->inputFile, "rb");
    if (fileptr == NULL)
        throwError("Could not read file %s", p->inputFile);

    if (fseek(fileptr, 0, SEEK_END))
        throwError("readFile: fseek END failed");
    size_t n = ftell(fileptr) / (2 * sizeof(float));

    size_t stdReducedSize = (size_t)ceil((double)n / p->nProcs);
    d->n = stdReducedSize;
    if (rank == p->nProcs-1)
        d->n = n - stdReducedSize * (p->nProcs-1);

    d->X = malloc(d->n * 2 * sizeof(float));
    if (d->X == NULL)
        throwError("readFile: Failed to allocate memory for points");
    d->Y = &d->X[d->n];

    // read X
    size_t offset = stdReducedSize * rank * sizeof(float);
    if (fseek(fileptr, offset, SEEK_SET));
        throwError("readFile: fseek SET X Failed");
    fread(d->X, d->n * sizeof(float), 1, fileptr);

    // read Y
    offset += n * sizeof(float);
    if (fseek(fileptr, offset, SEEK_SET));
        throwError("readFile: fseek SET Y Failed");
    fread(d->Y, d->n * sizeof(float), 1, fileptr);

    fclose(fileptr);

    return n;
}
#endif

void plotData(Data *points, Data *hull, int nUncovered, const char * title)
{
    // creating the pipeline for gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent 2>/dev/null ", "w");

    // gnuplot settings
    fprintf(gnuplotPipe, "set title \"%s\"\n", title);
    fprintf(gnuplotPipe, "set terminal qt size %s\n", GNUPLOT_RES);

    // set plot linestyles
    fprintf(gnuplotPipe, "set style line 1 linecolor rgb '%s' pt 7 pointsize %d\n", "red", 1);
    fprintf(gnuplotPipe, "set style line 2 linecolor rgb '%s' pt 7 pointsize %d\n", "green", 1);
    fprintf(gnuplotPipe, "set style line 3 linecolor rgb '%s' pt 7 pointsize %d\n", "blue", 1);
    fprintf(gnuplotPipe, "set style line 4 linecolor rgb '%s' pointsize 0\n", "black");

    // populating the plot
    
    fprintf(gnuplotPipe, "plot '-' with point linestyle 1, '-' with point linestyle 2, '-' with point linestyle 3, '-' with linespoint linestyle 4\n");

    // first plot only the points
    for (int i = hull->n + nUncovered; i < points->n; i++)
        fprintf(gnuplotPipe, "%f %f\n", points->X[i], points->Y[i]);
    fprintf(gnuplotPipe, "e\n");

    for (int i = hull->n; i < hull->n + nUncovered; i++)
        fprintf(gnuplotPipe, "%f %f\n", points->X[i], points->Y[i]);
    fprintf(gnuplotPipe, "e\n");

    for (int i = 0; i < hull->n; i++)
        fprintf(gnuplotPipe, "%f %f\n", hull->X[i], hull->Y[i]);
    fprintf(gnuplotPipe, "e\n");

    // second print the tour
    for (int i = 0; i < hull->n; i++)
        fprintf(gnuplotPipe, "%f %f\n", hull->X[i], hull->Y[i]);
    fprintf(gnuplotPipe, "%f %f\n", hull->X[0], hull->Y[0]);
    fprintf(gnuplotPipe, "e\n");

    // force write on stream
    fflush(gnuplotPipe);

    // close stream
    pclose(gnuplotPipe);
}

void plotHullMergeStep(Data *h1, Data *h2, Data *h0, size_t h1Index, size_t h2Index, const char * title, const bool closeH0)
{
    // creating the pipeline for gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent 2>/dev/null ", "w");

    // gnuplot settings
    fprintf(gnuplotPipe, "set title \"%s\"\n", title);
    fprintf(gnuplotPipe, "set terminal qt size %s\n", GNUPLOT_RES);

    // set plot linestyles
    fprintf(gnuplotPipe, "set style line 1 linecolor rgb '%s' pointsize %d\n", "forest-green", 2);
    fprintf(gnuplotPipe, "set style line 2 linecolor rgb '%s' pointsize %d\n", "green",      4);
    fprintf(gnuplotPipe, "set style line 3 linecolor rgb '%s' pointsize %d\n", "medium-blue",  2);
    fprintf(gnuplotPipe, "set style line 4 linecolor rgb '%s' pointsize %d\n", "cyan",       4);
    fprintf(gnuplotPipe, "set style line 5 linecolor rgb '%s' pointsize %d\n", "red",        1);

    // populating the plot
    fprintf(gnuplotPipe, "plot '-' with linespoint linestyle 1, '-' with point linestyle 2, '-' with linespoint linestyle 3, '-' with point linestyle 4, '-' with linespoint linestyle 5\n");

    // print h1
    for (int i = 0; i < h1->n; i++)
        fprintf(gnuplotPipe, "%f %f\n", h1->X[i], h1->Y[i]);
    fprintf(gnuplotPipe, "%f %f\ne\n", h1->X[0], h1->Y[0]);
    fprintf(gnuplotPipe, "%f %f\ne\n", h1->X[h1Index%h1->n], h1->Y[h1Index%h1->n]); // plot where are we pointing in h1

    // print h2
    for (int i = 0; i < h2->n; i++)
        fprintf(gnuplotPipe, "%f %f\n", h2->X[i], h2->Y[i]);
    fprintf(gnuplotPipe, "%f %f\ne\n", h2->X[0], h2->Y[0]);
    fprintf(gnuplotPipe, "%f %f\ne\n", h2->X[h2Index%h2->n], h2->Y[h2Index%h2->n]); // plot where are we pointing in h1

    // print h0
    for (int i = 0; i < h0->n; i++)
        fprintf(gnuplotPipe, "%f %f\n", h0->X[i], h0->Y[i]);
    if (closeH0)
        fprintf(gnuplotPipe, "%f %f\n", h0->X[0], h0->Y[0]);

    // force write on stream
    fflush(gnuplotPipe);

    // close stream
    pclose(gnuplotPipe);
}
