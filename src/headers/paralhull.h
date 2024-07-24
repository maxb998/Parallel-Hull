#include <stdlib.h>
#include <stdbool.h>

#define AVX_VEC_SIZE 8

#define swapElems(elem1,elem2) { register typeof(elem1) swapVarTemp = elem1; elem1 = elem2; elem2 = swapVarTemp; }

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
    bool gnuplot;
    
} Params;

typedef struct
{
    size_t n;
    float *X;
    float *Y;
} Data;

void setLogLevel(enum LogLevel lvl);
void LOG (enum LogLevel lvl, char * line, ...);
void throwError (char * line, ...);

Params argParse(int argc, char *argv[]);

Data emptyData();

void readFile(Data *d);