/******************************************************************************
 *
 *  File:           subgrid.c
 *  
 *  Created         23/08/2002
 *  
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *  
 *  Purpose:        Writes a subgrid of an input grid
 *
 *  Revisions:      none.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include "gucommon.h"
#include "gridnodes.h"

static void version()
{
    printf("subgrid/libgu version %s\n", gu_version);
    exit(0);
}

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

static void usage()
{
    printf("Usage: subgrid <grid file> [-i <imin>:<imax>] [-j <jmin>:<jmax>] [-x|-y] [-v]\n");
    printf("Try \"subgrid -h\" for more information\n");

    exit(0);
}

static void info()
{
    printf("Usage: subgrid <grid file> [-i <imin>:<imax>] [-j <jmin>:<jmax>]\n");
    printf("                [-x|-y] [-v]\n");
    printf("Where:\n");
    printf("  <grid file> -- text file with node coordinates (see remarks below)\n");
    printf("    (use \"stdin\" or \"-\" for standard input)\n");
    printf("Options:\n");
    printf("  -v -- verbose / version\n");
    printf("  -x -- print X coordinates only\n");
    printf("  -y -- print Y coordinates only\n");
    printf("Description:\n");
    printf("  `subgrid' reads grid nodes in \"X Y\\n\" format and writes nodes for the\n");
    printf("   specified subgrid\n");
    printf("Remarks:\n");
    printf("  1. The input must contain header describing the node array dimension:\n");
    printf("     ## <nx> x <ny>\n");
    printf("  2. After the header, the input must contain (nx * ny) lines with X and Y\n");
    printf("     node coordinates.\n");

    exit(0);
}

static void parse_commandline(int argc, char* argv[], char** fname, COORDTYPE* ct, int* imin, int* imax, int* jmin, int* jmax)
{
    int i;

    if (argc < 2)
        usage();

    i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') {
            *fname = argv[i];
            i++;
        } else {
            switch (argv[i][1]) {
            case 0:
                *fname = argv[i];
                i++;
                break;
            case 'i':
                i++;
                if (i == argc)
                    quit("no range found after \"-i\"\n");
                sscanf(argv[i], "%d:%d", imin, imax);
                i++;
                break;
            case 'j':
                i++;
                if (i == argc)
                    quit("no range found after \"-j\"\n");
                sscanf(argv[i], "%d:%d", jmin, jmax);
                i++;
                break;
            case 'h':
                info();
                break;
            case 'x':
                i++;
                *ct = CT_X;
                break;
            case 'y':
                i++;
                *ct = CT_Y;
                break;
            case 'v':
                i++;
                gu_verbose = 1;
                break;
            default:
                usage();
                break;
            }
        }
    }

    if (gu_verbose && argc == 2)
        version();

    if (*fname == NULL)
        usage();
}

int main(int argc, char* argv[])
{
    char* fname = NULL;
    gridnodes* grid = NULL;
    gridnodes* subgrid = NULL;
    COORDTYPE ct = CT_XY;
    int imin = INT_MIN;
    int imax = INT_MAX;
    int jmin = INT_MIN;
    int jmax = INT_MAX;

    parse_commandline(argc, argv, &fname, &ct, &imin, &imax, &jmin, &jmax);
    grid = gridnodes_read(fname, NT_NONE);
    subgrid = gridnodes_subgrid(grid, imin, imax, jmin, jmax);
    gridnodes_write(subgrid, "stdout", ct);
    gridnodes_destroy(subgrid);

    gridnodes_destroy(grid);

    return 0;
}
