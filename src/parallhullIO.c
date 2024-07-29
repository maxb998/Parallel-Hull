#include "paralhull.h"

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

void plotData(Data *points, int hullSize, int nUncovered, const char * plotPixelSize, const char * title)
{
    // creating the pipeline for gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent 2>/dev/null ", "w");

    // gnuplot settings
    fprintf(gnuplotPipe, "set title \"%s\"\n", title);
    fprintf(gnuplotPipe, "set terminal qt size %s\n", plotPixelSize);

    // set plot linestyles
    fprintf(gnuplotPipe, "set style line 1 linecolor rgb '%s' pt 7 pointsize %d\n", "red", 1);
    fprintf(gnuplotPipe, "set style line 2 linecolor rgb '%s' pt 7 pointsize %d\n", "green", 1);
    fprintf(gnuplotPipe, "set style line 3 linecolor rgb '%s' pt 7 pointsize %d\n", "blue", 1);
    fprintf(gnuplotPipe, "set style line 4 linecolor rgb '%s' pointsize 0\n", "black");

    // populating the plot
    
    fprintf(gnuplotPipe, "plot '-' with point linestyle 1, '-' with point linestyle 2, '-' with point linestyle 3, '-' with linespoint linestyle 4\n");

    // first plot only the points
    for (int i = hullSize + nUncovered; i < points->n; i++)
        fprintf(gnuplotPipe, "%f %f\n", points->X[i], points->Y[i]);
    fprintf(gnuplotPipe, "e\n");

    for (int i = hullSize; i < hullSize + nUncovered; i++)
        fprintf(gnuplotPipe, "%f %f\n", points->X[i], points->Y[i]);
    fprintf(gnuplotPipe, "e\n");

    for (int i = 0; i < hullSize; i++)
        fprintf(gnuplotPipe, "%f %f\n", points->X[i], points->Y[i]);
    fprintf(gnuplotPipe, "e\n");

    // second print the tour
    for (int i = 0; i < hullSize; i++)
        fprintf(gnuplotPipe, "%f %f\n", points->X[i], points->Y[i]);
    fprintf(gnuplotPipe, "%f %f\n", points->X[0], points->Y[0]);
    fprintf(gnuplotPipe, "e\n");

    // force write on stream
    fflush(gnuplotPipe);

    // close stream
    pclose(gnuplotPipe);
}