/******************************************************************************
 *
 *  File:           getnodes.c
 *  
 *  Created         22/01/2002
 *  
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *  
 *  Purpose:        Post-processes gridgen output grid file:
 *                    (i) reads input double density grid nodes
 *                    (ii) filters out nodes for partially defined cells
 *                    (iii) writes either all or center nodes
 *                  Also, can convert grid nodes of one type into another
 *
 *  Revisions:      none.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "nan.h"
#include "gucommon.h"
#include "gridnodes.h"

static void version()
{
    printf("getnodes/libgu version %s\n", gu_version);
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
    printf("Usage: getnodes <grid file> [-i <node type>] [-o <node type>]\n");
    printf("                [-m <mask file>] [-x|-y] [-v]\n");
    printf("Try \"getnodes -h\" for more information\n");

    exit(0);
}

static void info()
{
    printf("Usage: getnodes <grid file> [-i <node type>] [-o <node type>]\n");
    printf("                [-m <mask file>] [-v] [-x|-y]\n");
    printf("Where:\n");
    printf("  <grid file> -- text file with node coordinates (see remarks below)\n");
    printf("    (use \"stdin\" or \"-\" for standard input)\n");
    printf("Options:\n");
    printf("  -i <node type> -- input node type\n");
    printf("  -m <mask file> -- text file with nce1 x nce2 lines containing \"0\" or \"1\"\n");
    printf("    (use \"stdin\" or \"-\" for standard input)\n");
    printf("  -o <node type> -- output node type\n");
    printf("  -p -- tweak a bi- or tri- polar grid to a representation that can be\n");
    printf("        mapped (xy <-> ij) by using gridmap structure (node coordinates\n");
    printf("        are assumed to be in degrees)\n");
    printf("  -v -- verbose / version\n");
    printf("  -x -- print X coordinates only\n");
    printf("  -y -- print Y coordinates only\n");
    printf("Node types:\n");
    printf("  DD -- double density nodes (default) \n");
    printf("  CE -- cell center nodes\n");
    printf("  CO -- cell corner nodes\n");
    printf("Description:\n");
    printf("  `getnodes' (i) reads input grid nodes in \"X Y\" format;\n");
    printf("   (ii) validates them (sets node values for empty cells to NaNs);\n");
    printf("   (iii) if necessary, converts the nodes to a specified node type;\n");
    printf("   (iv) prints the requested nodes.\n");
    printf("Remarks:\n");
    printf("  1. The grid file must contain header describing the node array dimension:\n");
    printf("     ## <nx> x <ny>\n");
    printf("     where for double density nodes nx = nce1 * 2 + 1, ny = nce2 * 2 + 1;\n");
    printf("     for corner nodes  nx = nce1 + 1, ny = nce2 + 1; and for center nodes\n");
    printf("     nx = nce1, ny = nce2.\n");
    printf("  2. After the header, the grid file must contain (nx * ny) lines with X and Y\n");
    printf("     node coordinates.\n");
    printf("  3. An empty or commented line in the input grid file as well as NaNs for\n");
    printf("     node coordinates indicate an invalid node.\n");
    printf("  4. An optional mask file is a file with nce1 x nce2 lines containing \"1\" for\n");
    printf("     valid cells and \"0\" for invalid cells.\n");
    printf("  5. A grid cell is valid if all four corner nodes are valid (not NaNs).\n");
    printf("     If a cell mask was specified, then a valid corner node must also have\n");
    printf("     at least one valid neigbour cell.\n");
    printf("  6. The grid (union of all valid grid cells) must be simpy connected.\n");

    exit(0);
}

static void parse_commandline(int argc, char* argv[], char** gridfname, char** maskfname, NODETYPE* ntin, NODETYPE* ntout, COORDTYPE* ct, int* tweaknpolar)
{
    int i;

    if (argc < 2)
        usage();

    i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') {
            *gridfname = argv[i];
            i++;
        } else {
            switch (argv[i][1]) {
            case 0:
                *gridfname = argv[i];
                i++;
                break;
            case 'i':
                i++;
                if (i == argc)
                    quit("no node type found after \"-i\"\n");
                if (strcasecmp("dd", argv[i]) == 0)
                    *ntin = NT_DD;
                else if (strcasecmp("ce", argv[i]) == 0)
                    *ntin = NT_CEN;
                else if (strcasecmp("co", argv[i]) == 0)
                    *ntin = NT_COR;
                else
                    quit("input node type \"%s\" not recognised\n", argv[i]);
                i++;
                break;
            case 'h':
                info();
                break;
            case 'm':
                i++;
                if (i == argc)
                    quit("no file name found after \"-m\"\n");
                *maskfname = argv[i];
                i++;
                break;
            case 'o':
                i++;
                if (i == argc)
                    quit("no node type found after \"-o\"\n");
                if (strcasecmp("dd", argv[i]) == 0)
                    *ntout = NT_DD;
                else if (strcasecmp("ce", argv[i]) == 0)
                    *ntout = NT_CEN;
                else if (strcasecmp("co", argv[i]) == 0)
                    *ntout = NT_COR;
                else
                    quit("input node type \"%s\" not recognised\n", argv[i]);
                i++;
                break;
            case 'p':
                i++;
                *tweaknpolar = 1;
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

    if (*gridfname == NULL)
        usage();
}

static void gridnodes_tweaknpolar(gridnodes* gn)
{
    int nx = gridnodes_getnx(gn);
    int ny = gridnodes_getny(gn);
    double** gx = gridnodes_getx(gn);
    double** gy = gridnodes_gety(gn);
    int i, j;

    for (i = 0; i < nx; ++i)
        for (j = 1; j < ny; ++j)
            if (fabs(gx[j][i] - gx[j - 1][i]) > 180.0)
                gx[j][i] += 360.0;

    for (j = 0; j < ny; ++j)
        for (i = 1; i < nx; ++i)
            if (fabs(gx[j][i] - gx[j][i - 1]) > 180.0) {
                gx[j][i - 1] = NaN;
                gy[j][i - 1] = NaN;
            }

    for (j = 0; j < ny; ++j)
        for (i = 1; i < nx; ++i)
            if (fabs(gy[j][i] - gy[j][i - 1]) > 90.0) {
                gx[j][i - 1] = NaN;
                gy[j][i - 1] = NaN;
            }
}

int main(int argc, char* argv[])
{
    char* gridfname = NULL;
    char* maskfname = NULL;
    NODETYPE ntin = NT_DD;
    NODETYPE ntout = NT_DD;
    COORDTYPE ct = CT_XY;
    int tweaknpolar = 0;
    gridnodes* gn = NULL;

    parse_commandline(argc, argv, &gridfname, &maskfname, &ntin, &ntout, &ct, &tweaknpolar);
    gn = gridnodes_read(gridfname, ntin);
    gridnodes_validate(gn);
    if (ntin != ntout) {
        gridnodes* gnnew = gridnodes_transform(gn, ntout);

        gridnodes_destroy(gn);
        gn = gnnew;
    }
    if (maskfname != NULL) {
        int nx = gridnodes_getnce1(gn);
        int ny = gridnodes_getnce2(gn);
        int** mask = gu_readmask(maskfname, nx, ny);

        gridnodes_applymask(gn, mask);
        gu_free2d(mask);
    }
    if (tweaknpolar)
        gridnodes_tweaknpolar(gn);
    if (gu_verbose)
        gridnodes_calcstats(gn);
    gridnodes_write(gn, "stdout", ct);
    gridnodes_destroy(gn);

    return 0;
}
