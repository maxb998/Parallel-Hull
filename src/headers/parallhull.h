#include <stdlib.h>
#include <stdbool.h>

#define AVX_VEC_SIZE 8

// #define QUICKHULL_STEP_DEBUG // plots data useful for debug at each iteration of the quickhull algorithm
#define DEBUG
#define GUI_OUTPUT
// #define NON_MPI_MODE

#define GNUPLOT_RES "1600,900"

#define swapElems(elem1,elem2) { register typeof(elem1) swapVarTemp = elem1; elem1 = elem2; elem2 = swapVarTemp; }


// Convert timespec struct to a double time in seconds
#define cvtTimespec2Double(t) (double)t.tv_sec + (double)t.tv_nsec / 1000000000.0

enum LogLevel
{
	LOG_LVL_FATAL,
	LOG_LVL_ERROR,
	LOG_LVL_WARN,
	LOG_LVL_NOTICE,
	LOG_LVL_INFO,
	LOG_LVL_DEBUG,
	LOG_LVL_TRACE
};

typedef struct
{
    int nProcs;
    int procID;
    int nThreads;

    char inputFile[1000];
    enum LogLevel logLevel;
    
} Params;

typedef struct
{
    size_t n;
    float *X;
    float *Y;
} Data;

void setLogLevel(enum LogLevel lvl);
void LOG (enum LogLevel lvl, char * line, ...);
void inlineLOG (enum LogLevel lvl, char * line, ...);
void throwError (char * line, ...);

Params argParse(int argc, char *argv[]);

void readFile(Data *d, Params *p);

void plotData(Data *points, Data *hull, int nUncovered, const char * plotPixelSize, const char * title);

size_t quickhull(Data *d, int threadID);

Data parallhullThreaded(Data *d, size_t reducedProblemUB, int nThreads);