#include "paralhull.h"

#include <argp.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#define SUBOPT_BLANKSPACE " "
#define SUBOPT_LOG_ERROR "error"
#define SUBOPT_LOG_CRITICAL "critical"
#define SUBOPT_LOG_WARNING "warning"
#define SUBOPT_LOG_NOTICE "notice"
#define SUBOPT_LOG_INFO "info"
#define SUBOPT_LOG_DEBUG "debug"
#define SUBOPT_LOG_TRACE "trace"
#define LOG_LEVEL_DOC "\
Specify the log level (DEFAULT=log)\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_ERROR "\t\t: Show only error messages\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_CRITICAL "\t: Show critical messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_WARNING "\t: Show warning and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_NOTICE "\t: Show notice messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_INFO "\t\t: Show info messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_DEBUG "\t\t: Show debug messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_TRACE "\t\t: Show all messages\n"
static const char *logLevelStrings[] = { SUBOPT_LOG_ERROR, SUBOPT_LOG_CRITICAL, SUBOPT_LOG_WARNING, SUBOPT_LOG_NOTICE, SUBOPT_LOG_INFO, SUBOPT_LOG_DEBUG, SUBOPT_LOG_TRACE };
static const int loglvlsCount = sizeof(logLevelStrings)/sizeof(*logLevelStrings);

enum argpKeys{
    ARGP_FILE='f',
    ARGP_NTHREADS='j',
    ARGP_LOG_LEVEL='l'
};


const char * logLevelString [] = {
	"[\033[1;31mFATAL\033[0m]", // 0
	"[\033[0;35mERROR\033[0m]", // 1
	"[\033[0;33mWARN\033[0m]", // 2
	"[\033[0;36mNOTICE\033[0m]", // 3
	"[\033[0;34mINFO\033[0m]", // 4
	"[\033[0;32mDEBUG\033[0m]", // 5
	"[\033[0;90mTRACE\033[0m]"  // 6
};

error_t argpParser(int key, char *arg, struct argp_state *state);
static void parseEnumOption(char *arg, int *savePtr, const char **optionsSet, const int from, const int to, const char *optionName);
static int parseUint(char *arg, char expectedEndChr, const char *paramName);
static Params emptyParams();

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
    printf("%s ", logLevelString[lvl]);
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

Params argParse(int argc, char *argv[])
{
    static struct argp_option argpOptions[] = {
        { .name="file", .key=ARGP_FILE, .arg="FILENAME", .flags=0, .doc="Location of the file containing the points used calculate the hull\n", .group=1 },
        { .name="threads", .key=ARGP_NTHREADS, .arg="UINT", .flags=0, .doc="Number of threads to use\n", .group=1 },
        { .name="loglvl", .key=ARGP_LOG_LEVEL, .arg="STRING", .flags=0, .doc=LOG_LEVEL_DOC, .group=1 },
        { 0 }
    };

    static struct argp argpData = {
        .options = argpOptions,
        .parser=argpParser, 
    };

    Params p = emptyParams();
    argp_parse(&argpData, argc, argv, 0, 0, &p);

    return p;
}

error_t argpParser(int key, char *arg, struct argp_state *state)
{
    Params *p = state->input;

    switch (key)
    {
    case ARGP_FILE: // get the input filename
        if (access(arg, R_OK)) // check if file exists and is accessible
        {
            LOG(LOG_LVL_ERROR, "File \"%s\" cannot be accessed or does not exist", arg);
            return ARGP_ERR_UNKNOWN;
        }
        strncpy(p->inputFile, arg, strlen(arg));
        break;

    case ARGP_NTHREADS:
        p->nThreads = parseUint(arg, 0, "nThreads");
        break;

    case ARGP_LOG_LEVEL:
        parseEnumOption(arg, (int*)&p->logLevel, logLevelStrings, 0, loglvlsCount, "loglvl");
        setLogLevel(p->logLevel);
        break;
    
    case ARGP_KEY_END:
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static void parseEnumOption(char *arg, int *savePtr, const char **optionsSet, const int from, const int to, const char *optionName)
{
    for (int i = from; i < to; i++)
    {
        if (strcmp(arg, optionsSet[i]) == 0)
        {
            *savePtr = i;
            return;
        }
    }
    
    throwError("%s: argument not valid", optionName);
}

static int parseUint(char *arg, char expectedEndChr, const char *paramName)
{
    char *endPtr;
    long cvt = strtol(arg, &endPtr, 10);
    if (cvt < 0)
        throwError("The value specified as %s cannot be negative", paramName);
    if (*endPtr != expectedEndChr)
        throwError("There are extra character after the %s value or formatting is not correct. Check formats with --help", paramName);

    return (int)cvt;
}

static Params emptyParams()
{
    Params p = {
        .inputFile={0},
        .logLevel=LOG_LVL_INFO,
        .nProcs=-1,
        .nThreads=1,
        .procID=-1,
        .gnuplot=false
    };

    return p;
}

Data emptyData()
{
    Data d = {
        .n=0,
        .X=NULL,
        .Y=NULL
    };

    return d;
}

Data cloneData(Data *d)
{
    Data r = *d;
    r.X = malloc(r.n * 2 * sizeof(float));
    if (r.X == NULL)
        throwError("cloneData: Failed to allocate memory for points");
    r.Y = &r.X[r.n+1];

    for (int i = 0; i < d->n; i++)
    {
        r.X[i] = d->X[i];
        r.Y[i] = d->Y[i];
    }
    
    return r;
}

void readFile(Data *d, Params *p)
{
    char *r = p->inputFile;

    while (*r && !isdigit(*r))
        r++;

    // Extract number from the file name
    for (; *r && isdigit(*r); r++)
        d->n = d->n * 10 + (*r - '0');

    // Allocate memory for the points
    d->X = malloc(d->n * 2 * sizeof(float));
    if (d->X == NULL)
        throwError("Failed to allocate memory for points");
    d->Y = &d->X[d->n];

    FILE* file = fopen(p->inputFile, "r");
    if (file == NULL)
        throwError("Could not read file %s", p->inputFile);

    // Read each line to get the points
    char line[1000];
    for (size_t i = 0; i < d->n; i++)
    {
        if (fgets(line, sizeof(line), file) == NULL)
        {
            fclose(file);
            throwError("Failed to read coordinates for point %d", i + 1);
        }

        // Parse the coordinates
        if (sscanf(line, "%f,%f", &d->X[i], &d->Y[i]) != 2)
        {
            fclose(file);
            throwError("Invalid format for coordinates on line %d", i + 2);
        }
    }

    fclose(file);
}

void plotData(Data *points, int hullSize, const char * plotPixelSize, const char * pointColor, const char * tourPointColor, const int pointSize)
{

    // creating the pipeline for gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent 2>/dev/null", "w");

    // gnuplot settings
    fprintf(gnuplotPipe, "set title \"Convex Hull\"\n");
    fprintf(gnuplotPipe, "set terminal qt size %s\n", plotPixelSize);

    // set plot linestyles
    fprintf(gnuplotPipe, "set style line 1 linecolor rgb '%s' pt 7 pointsize %d\n", pointColor, pointSize);
    fprintf(gnuplotPipe, "set style line 2 linecolor rgb '%s' pointsize 0\n", tourPointColor);//, pointSize);

    // populating the plot
    
    fprintf(gnuplotPipe, "plot '-' with point linestyle 1, '-' with linespoint linestyle 2\n");

    // first plot only the points
    for (int i = 0; i < points->n; i++)
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

void plotPartialData(Data *points, int hullSize, int nUncovered, const char * plotPixelSize)
{

    // creating the pipeline for gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent 2>/dev/null ", "w");

    // gnuplot settings
    fprintf(gnuplotPipe, "set title \"Convex Hull\"\n");
    fprintf(gnuplotPipe, "set terminal qt size %s\n", plotPixelSize);

    // set plot linestyles
    fprintf(gnuplotPipe, "set style line 1 linecolor rgb '%s' pt 7 pointsize %d\n", "red", 2);
    fprintf(gnuplotPipe, "set style line 2 linecolor rgb '%s' pt 7 pointsize %d\n", "green", 2);
    fprintf(gnuplotPipe, "set style line 3 linecolor rgb '%s' pointsize 0\n", "black");

    // populating the plot
    
    fprintf(gnuplotPipe, "plot '-' with point linestyle 1, '-' with point linestyle 2, '-' with linespoint linestyle 3\n");

    // first plot only the points
    for (int i = 0; i < points->n; i++)
        fprintf(gnuplotPipe, "%f %f\n", points->X[i], points->Y[i]);
    fprintf(gnuplotPipe, "e\n");

    for (int i = 0; i < hullSize + nUncovered; i++)
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