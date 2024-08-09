#include "parallhull.h"

#include <argp.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>

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

error_t argpParser(int key, char *arg, struct argp_state *state);
static void parseEnumOption(char *arg, int *savePtr, const char **optionsSet, const int from, const int to, const char *optionName);
static int parseUint(char *arg, char expectedEndChr, const char *paramName);

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

    Params p = {
        .inputFile={0},
        .logLevel=LOG_LVL_INFO,
        .nProcs=-1,
        .nThreads=1,
        .procID=-1
    };
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
        strncpy(p->inputFile, arg, 1000);
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