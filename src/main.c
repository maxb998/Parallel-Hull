#include "paralhull.h"

#include <stdio.h>
#include <time.h>

#define SEPARATOR_STR "##############################################################################################################################\n"



int main (int argc, char *argv[])
{
    Data d = emptyData();
    Params p = argParse(argc, argv);

    readFile(&d, &p);

    int hullSize = quickhull(&d);

    plotData(&d, hullSize, "1600,900", "green", "black", 1);

    return EXIT_SUCCESS;
}
