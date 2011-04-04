/******************************************************************************
 *
 * File:           main.c
 *
 * Created:        02/11/2006
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Standalone utility for generation of 2D orthogonal grids. 
 *
 * Description:    For more information, see README and gridgen.c.
 *                 
 * Revisions:      02/10/2006 PS Minor structural changes associated with
 *                            the introduction of libgridgen.a. Changed the
 *                            structure name "gridspec" to "gridgen".
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include "gridgen.h"

static void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "error: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    exit(1);
}

static void usage(void)
{
    printf("Usage: gridgen [options] <parameter file>\n");
    printf("       gridgen -a\n");
    printf("       gridgen -h\n");
    printf("Options:\n");
    printf("  -a -- describe the algorithm used\n");
    printf("  -h -- describe parameter file format\n");
    printf("  -v -- verbose / version\n");
    printf("  -V -- very verbose\n");
    printf("  -VV -- extremely verbose\n");

    exit(0);
}

static void parse_commandline(int argc, char* argv[], char** prmfname)
{
    int i;

    for (i = 1; i < argc; ++i) {
        if (argv[i][0] != '-')
            *prmfname = argv[i];
        else
            switch (argv[i][1]) {
            case 'a':
                gridgen_printhelpalg();
                exit(0);
                break;
            case 'h':
                gridgen_printhelpprm();
                exit(0);
                break;
            case 'p':
                i++;
                if (i >= argc)
                    quit("no file name found after (obsolete) -p\n");
                if (argv[i][0] == '-')
                    quit("do not begin parameter file name with '-'\n");
                *prmfname = argv[i];
                break;
            case 'v':
                gridgen_setverbose(1);
                break;
            case 'V':
                gridgen_setverbose((argv[i][2] == 0) ? 2 : 3);
                break;
            default:
                usage();
                break;
            }
    }

    if (argc == 2 && tolower(argv[1][1]) == 'v') {
        gridgen_printversion();
        exit(0);
    }

    if (*prmfname == NULL)
        usage();
}

int main(int argc, char* argv[])
{
    char* prmfname = NULL;

    parse_commandline(argc, argv, &prmfname);
    gridgen_generategrid(prmfname);

    return 0;
}
