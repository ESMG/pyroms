/******************************************************************************
 *
 *  File:           insertgrid.c
 *  
 *  Created         22/04/2005
 *  
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *  
 *  Purpose:        Insert a smaller grid inside a larger one. The inserted
 *                  nodes overwrite existing ones. Pretty straightforward.
 *
 *  Revisions:      PS 06/07/2005 Modified to allow the following:
 *                    (i) arbitrary positions of the grids relative to each
 *                        other
 *                    (ii) merge of two grids: if second grid has an invalid
 *                        node while the first grid has a valid node, use the
 *                        node from the first grid; otherwhile use node from
 *                        the second grid
 *                    (iii) specifying the grid node type.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "nan.h"
#include "gucommon.h"
#include "gridnodes.h"

int merge = 0;

static void version()
{
    printf("insertgrid/libgu version %s\n", gu_version);
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
    printf("Usage: insertgrid <master grid file> <grid file> <i> <j> [-i {DD|CO|CE}] [-m] [-v]\n");
    printf("       insertgrid -h\n");
    printf("       insertgrid -v\n");
    printf("Try \"insertgrid -h\" for more information\n");

    exit(0);
}

static void info()
{
    printf("Usage: insertgrid <master grid file> <grid file> <i> <j> [-v]\n");
    printf("       insertgrid -h\n");
    printf("       insertgrid -v\n");
    printf("Where:\n");
    printf("  <master grid file> -- text file with node coordinates of the grid to be\n");
    printf("                   modified (see remarks below) (use \"stdin\" or \"-\" for\n");
    printf("                   standard input)\n");
    printf("  <grid file>   -- text file with nodes of the grid to be inserted\n");
    printf("  <i>, <j>      -- indices of location of node (i,j) = (0,0) of the slave\n");
    printf("                   grid in the master grid\n");
    printf("  -i {DD|CO|CE} -- grid node type (double density, corner or center)\n");
    printf("Options:\n");
    printf("  -h -- print this information\n");
    printf("  -m -- merge grids (default = insert)\n");
    printf("  -v -- verbose / version\n");
    printf("Description:\n");
    printf("  `insertgrid' reads grid nodes from two files and inserts nodes from one\n");
    printf("  grid into another grid at specified location in index space.\n");
    printf("Remarks:\n");
    printf("  1. The input grid files must contain header describing the node array\n");
    printf("     dimension: ## <nx> x <ny>\n");
    printf("  2. After the header, the input files must contain (nx * ny) lines with\n");
    printf("     X and Y node coordinates.\n");
    printf("  3. The output grid will be verified if only the node type is specified.\n");

    exit(0);
}

static void parse_commandline(int argc, char* argv[], char** fname1, char** fname2, int* i0, int* j0, NODETYPE* type)
{
    int i;

    *fname1 = NULL;
    *fname2 = NULL;
    *i0 = -1;
    *j0 = -1;
    *type = NT_NONE;

    for (i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            if (*fname1 == NULL)
                *fname1 = argv[i];
            else if (*fname2 == NULL)
                *fname2 = argv[i];
            else if (*i0 == -1)
                *i0 = atoi(argv[i]);
            else if (*j0 == -1)
                *j0 = atoi(argv[i]);
            else
                usage();
        } else {
            switch (argv[i][1]) {
            case 0:
                if (*fname1 == NULL)
                    *fname1 = argv[i];
                else if (*fname2 == NULL)
                    *fname2 = argv[i];
                else
                    usage();
                break;
            case 'h':
                info();
                break;
            case 'i':
                i++;
                if (i == argc)
                    quit("no node type found after \"-i\"\n");
                if (strcasecmp("dd", argv[i]) == 0)
                    *type = NT_DD;
                else if (strcasecmp("ce", argv[i]) == 0)
                    *type = NT_CEN;
                else if (strcasecmp("co", argv[i]) == 0)
                    *type = NT_COR;
                else
                    quit("input node type \"%s\" not recognised\n", argv[i]);
                break;
            case 'm':
                merge = 1;
                break;
            case 'v':
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

    if (*fname1 == NULL || *fname2 == NULL || *i0 == -1 || *j0 == -1)
        usage();
}

int main(int argc, char* argv[])
{
    char* fname1 = NULL;
    char* fname2 = NULL;
    gridnodes* grid1 = NULL;
    gridnodes* grid2 = NULL;
    gridnodes* grid3 = NULL;
    NODETYPE type = NT_NONE;
    int i0 = INT_MAX;
    int j0 = INT_MAX;
    int nx1 = 0;
    int ny1 = 0;
    int nx2 = 0;
    int ny2 = 0;
    int nx3 = 0;
    int ny3 = 0;
    double** gx1 = NULL;
    double** gy1 = NULL;
    double** gx2 = NULL;
    double** gy2 = NULL;
    double** gx3 = NULL;
    double** gy3 = NULL;
    int i1start, i2start, j1start, j2start;
    int i1, i2, j1, j2, i3, j3;

    parse_commandline(argc, argv, &fname1, &fname2, &i0, &j0, &type);
    grid1 = gridnodes_read(fname1, NT_NONE);
    grid2 = gridnodes_read(fname2, NT_NONE);

    nx1 = gridnodes_getnx(grid1);
    ny1 = gridnodes_getny(grid1);
    nx2 = gridnodes_getnx(grid2);
    ny2 = gridnodes_getny(grid2);

    gx1 = gridnodes_getx(grid1);
    gy1 = gridnodes_gety(grid1);
    gx2 = gridnodes_getx(grid2);
    gy2 = gridnodes_gety(grid2);

    nx3 = (i0 + nx2 > nx1) ? i0 + nx2 : nx1;
    ny3 = (j0 + ny2 > ny1) ? j0 + ny2 : ny1;
    grid3 = gridnodes_create(nx3, ny3, type);
    gx3 = gridnodes_getx(grid3);
    gy3 = gridnodes_gety(grid3);

    if (i0 < 0) {
        i1start = i0;
        i2start = 0;
    } else {
        i1start = 0;
        i2start = -i0;
    }
    if (j0 < 0) {
        j1start = j0;
        j2start = 0;
    } else {
        j1start = 0;
        j2start = -j0;
    }
    if (!merge) {
        for (j3 = 0, j1 = j1start, j2 = j2start; j3 < ny3; ++j3, ++j1, ++j2) {
            for (i3 = 0, i1 = i1start, i2 = i2start; i3 < nx3; ++i3, ++i1, ++i2) {
                int ok2 = (j2 >= 0 && j2 < ny2 && i2 >= 0 && i2 <= nx2);

                if (ok2) {
                    gx3[j3][i3] = gx2[j2][i2];
                    gy3[j3][i3] = gy2[j2][i2];
                } else {
                    int ok1 = (j1 >= 0 && j1 < ny1 && i1 >= 0 && i1 <= nx1);

                    if (ok1) {
                        gx3[j3][i3] = gx1[j1][i1];
                        gy3[j3][i3] = gy1[j1][i1];
                    } else {
                        gx3[j3][i3] = NaN;
                        gy3[j3][i3] = NaN;
                    }
                }
            }
        }
    } else {
        for (j3 = 0, j1 = j1start, j2 = j2start; j3 < ny3; ++j3, ++j1, ++j2) {
            for (i3 = 0, i1 = i1start, i2 = i2start; i3 < nx3; ++i3, ++i1, ++i2) {
                int ok2 = (j2 >= 0 && j2 < ny2 && i2 >= 0 && i2 <= nx2);
                int ok1 = (j1 >= 0 && j1 < ny1 && i1 >= 0 && i1 <= nx1);

                if (!ok2 && !ok1) {
                    gx3[j3][i3] = NaN;
                    gy3[j3][i3] = NaN;
                } else if (ok2 && !ok1) {
                    gx3[j3][i3] = gx2[j2][i2];
                    gy3[j3][i3] = gy2[j2][i2];
                } else if (!ok2 && ok1) {
                    gx3[j3][i3] = gx1[j1][i1];
                    gy3[j3][i3] = gy1[j1][i1];
                } else {
                    if (!isnan(gx2[j2][i2])) {
                        gx3[j3][i3] = gx2[j2][i2];
                        gy3[j3][i3] = gy2[j2][i2];
                    } else {
                        gx3[j3][i3] = gx1[j1][i1];
                        gy3[j3][i3] = gy1[j1][i1];
                    }
                }
            }
        }
    }

    gridnodes_validate(grid3);
    gridnodes_write(grid3, "stdout", CT_XY);

    gridnodes_destroy(grid1);
    gridnodes_destroy(grid2);
    gridnodes_destroy(grid3);

    return 0;
}
