#include "paralhull.h"

#include <stdio.h>
#include <time.h>

#define SEPARATOR_STR "##############################################################################################################################\n"



int main (int argc, char *argv[])
{
    Data d = emptyData();
    Params p = argParse(argc, argv);

    readFile(&d);


    

    return EXIT_SUCCESS;
}
