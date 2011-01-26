/******************************************************************************
 *
 *  File:           getbound.c
 *  
 *  Created         13/08/2002
 *  
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *  
 *  Purpose:        Extract boundary polygon from grid points
 *
 *  Revisions:      none.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "gucommon.h"
#include "gridnodes.h"
#include "poly.h"

static void version()
{
    printf("getbound/libgu version %s\n", gu_version);
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
    printf("Usage: getbound <gridfile> [-c] [-i {DD|CO}] [-r] [-v]\n");
    printf("Try \"getbound -h\" for more information\n");

    exit(0);
}

static void info()
{
    printf("Usage: getbound [options] <grid file>\n");
    printf("Where:\n");
    printf("  <grid file> -- text file with node coordinates (see remarks below)\n");
    printf("    (use \"stdin\" or \"-\" for standard input)\n");
    printf("Options:\n");
    printf("  -c -- compact the polygon (exclude vertices that lie in between other\n");
    printf("        vertices)\n");
    printf("  -i <node type> -- input node type\n");
    printf("  -r -- Write boundary polygon vertex indices rather than physical coordinates\n");
    printf("  -v -- Verbose / version\n");
    printf("Node types:\n");
    printf("  DD -- double density nodes (default) \n");
    printf("  CO -- cell corner nodes\n");
    printf("Description:\n");
    printf("  `getbound' reads grid nodes from standard input. After that, it builds the\n");
    printf("   boundary polygon and writes it to standard output\n");
    printf("Remarks:\n");
    printf("  1. The grid file must contain header describing the node array dimension:\n");
    printf("     ## <nx> x <ny>\n");
    printf("     where for double density nodes nx = nce1 * 2 + 1, ny = nce2 * 2 + 1;\n");
    printf("     for corner nodes  nx = nce1 + 1, ny = nce2 + 1; and for center nodes\n");
    printf("     nx = nce1, ny = nce2.\n");
    printf("  2. After the header, the grid file must contain (nx * ny) lines with X and\n");
    printf("     Y node coordinates.\n");
    printf("  3. An empty or commented line in the input grid file as well as NaNs for\n");
    printf("     node coordinates indicate an invalid node.\n");
    printf("  4. A grid cell is valid if all corner nodes are valid (not NaNs).\n");
    printf("  5. The grid (union of all valid grid cells) must be simpy connected.\n");

    exit(0);
}

static void parse_commandline(int argc, char* argv[], char** fname, int* compact, int* ij, NODETYPE* nt)
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
            case 'c':
                *compact = 1;
                i++;
                break;
            case 'i':
                i++;
                if (i == argc)
                    quit("no node type found after \"-i\"\n");
                if (strcasecmp("dd", argv[i]) == 0)
                    *nt = NT_DD;
                else if (strcasecmp("ce", argv[i]) == 0)
                    quit("cell centre node type is not supported by getbound\n");
                else if (strcasecmp("co", argv[i]) == 0)
                    *nt = NT_COR;
                else
                    quit("input node type \"%s\" not recognised\n", argv[i]);
                i++;
                break;
            case 'h':
                info();
                break;
            case 'r':
                i++;
                *ij = 1;
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
    int ij = 0;
    int compact = 0;
    NODETYPE nt = NT_DD;
    gridnodes* gn = NULL;
    poly* pl = NULL;

    parse_commandline(argc, argv, &fname, &compact, &ij, &nt);

    if (nt == NT_DD) {
        /*
         * read DD grid nodes 
         */
        gridnodes* gndd = gridnodes_read(fname, NT_DD);

        gridnodes_validate(gndd);
        /*
         * get corner grid nodes from DD grid nodes 
         */
        gn = gridnodes_transform(gndd, NT_COR);
        gridnodes_destroy(gndd);
    } else {
        gn = gridnodes_read(fname, NT_COR);
        gridnodes_validate(gn);
    }

    /*
     * build boundary polygon 
     */
    if (!ij)
        pl = poly_formbound(gridnodes_getnce1(gn), gridnodes_getnce2(gn), gridnodes_getx(gn), gridnodes_gety(gn));
    else
        pl = poly_formboundij(gridnodes_getnce1(gn), gridnodes_getnce2(gn), gridnodes_getx(gn));

    if (compact)
        poly_compact(pl, 1.0e-10);

    poly_write(pl, stdout);

    poly_destroy(pl);
    gridnodes_destroy(gn);

    return 0;
}
