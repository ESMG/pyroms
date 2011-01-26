/******************************************************************************
 *
 *  File:           setbathy.c
 *  
 *  Created         23/08/2002
 *  
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *  
 *  Purpose:        In a given array of values associated with a grid,
 *                  set values corresponding to certain indices to
 *                  predefined values.
 *
 *  Revisions:      none.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include "gucommon.h"

#define BUFSIZE 10240

static void version()
{
    printf("setbathy/libgu version %s\n", gu_version);
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
    printf("Usage: setbathy <file> [-i <imin>:<imax>] [-j <jmin>:<jmax>]\n");
    printf("                [-d <nx>x<ny>] [-o <line No. offset>] [-v] -z <value>\n");
    printf("Description: in a text file containing <nx> x <ny> lines with some data\n");
    printf("  e.g. bathymetry values, `setbathy' sets values in specified index range to\n");
    printf("  a specified value. The first (or pointed) line of the file should contain\n");
    printf("  the header \"## <nx> x <ny>\" with grid dimensions; otherwhile one must\n");
    printf("  enter them by using \"-d\" option.\n");

    exit(0);
}

static void parse_commandline(int argc, char* argv[], char** fname, int* imin, int* imax, int* jmin, int* jmax, int* nx, int* ny, int* offset, char** z)
{
    int i;

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
            case 'd':
                i++;
                if (i == argc)
                    quit("no dimension found after \"-d\"\n");
                if (sscanf(argv[i], "%dx%d", nx, ny) != 2)
                    usage();
                i++;
                break;
            case 'o':
                i++;
                if (i == argc)
                    quit("no offset found after \"-d\"\n");
                *offset = atoi(argv[i]);
                i++;
                break;
            case 'v':
                i++;
                gu_verbose = 1;
                break;
            case 'z':
                i++;
                if (i == argc)
                    quit("no value found after \"-z\"\n");
                *z = argv[i];
                i++;
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
    FILE* f = NULL;
    int imin = INT_MIN;
    int imax = INT_MAX;
    int jmin = INT_MIN;
    int jmax = INT_MAX;
    int nx = -1;
    int ny = -1;
    char* z = NULL;
    int offset = 0;
    char buf[BUFSIZE];
    int i, j;

    parse_commandline(argc, argv, &fname, &imin, &imax, &jmin, &jmax, &nx, &ny, &offset, &z);

    if (strcmp(fname, "stdin") == 0 || strcmp(fname, "-") == 0) {
        fname = strdup("stdin");
        f = stdin;
    } else
        f = gu_fopen(fname, "r");

    if (gu_verbose && offset > 0)
        fprintf(stderr, "## skipping %d lines...\n", offset);
    for (i = 0; i < offset; ++i)
        if (fgets(buf, BUFSIZE, f) == NULL)
            quit("%s: could not read %d-th line\n", fname, i);
        else
            fprintf(stdout, "%s", buf);

    if (nx < 0) {
        if (gu_verbose)
            fprintf(stderr, "## reading header...\n");
        /*
         * get grid size 
         */
        if (fgets(buf, BUFSIZE, f) == NULL)
            quit("%s: empty input\n", fname);

        if (sscanf(buf, "## %d x %d", &nx, &ny) != 2)
            quit("%s: could not read grid size: expected header in \"## %%d x %%d\" format\n or \"-d\" option\n", fname);
    }

    if (gu_verbose)
        fprintf(stderr, "## %d x %d grid\n", nx, ny);

    if (nx < 1)
        quit("gridnodes_read(): nx = %d: invalid grid size\n", nx);
    if (ny < 1)
        quit("gridnodes_read(): ny = %d: invalid grid size\n", ny);
    if ((double) nx * (double) ny > (double) INT_MAX)
        quit("gridnodes_read(): grid size (%d x %d) is too big\n", nx, ny);

    for (j = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            if (fgets(buf, BUFSIZE, f) == NULL)
                quit("%s: could not read %d-th point (%d x %d points expected)\n", fname, j * nx + i + 1, nx, ny);
            if (i >= imin && i <= imax && j >= jmin && j <= jmax)
                fprintf(stdout, "%s\n", z);
            else
                fprintf(stdout, "%s", buf);
        }
    }

    while (fgets(buf, BUFSIZE, f) != NULL)
        fprintf(stdout, "%s", buf);

    if (f == stdin)
        free(fname);

    return 0;
}
