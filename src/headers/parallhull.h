#include <stdlib.h>
#include <stdbool.h>

#define AVX_VEC_SIZE 8

// #define QUICKHULL_STEP_DEBUG // plots data useful for debug at each iteration of the quickhull algorithm
// #define PARALLHULL_STEP_DEBUG
// #define PARALLHULL_MERGE_OUTPUT_PLOT
#define DEBUG
#define GUI_OUTPUT
#define NON_MPI_MODE

#define GNUPLOT_RES "1920,1080"

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

typedef struct{
    int p;
    int t;
} ProcThreadIDCombo;

void setLogLevel(enum LogLevel lvl);
void LOG (enum LogLevel lvl, char * line, ...);
void inlineLOG (enum LogLevel lvl, char * line, ...);
void throwError (char * line, ...);

Params argParse(int argc, char *argv[]);

#ifdef NON_MPI_MODE
void readFile(Data *d, Params *p);
#else
// returns full number of elements in file
size_t readFile(Data *d, Params *p, int rank);
#endif

void plotData(Data *points, Data *hull, int nUncovered, const char * title);

void plotHullMergeStep(Data *h1, Data *h2, Data *h0, size_t h1Index, size_t h2Index, const char * title, const bool closeH0);

size_t quickhull(Data *d, ProcThreadIDCombo *id);

Data parallhullThreaded(Data *d, size_t reducedProblemUB, int procID, int nThreads);