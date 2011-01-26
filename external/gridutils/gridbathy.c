/******************************************************************************
 *
 *  File:           gridbathy.c
 *
 *  Created:        04/08/2000
 *
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *
 *  Purpose:        For a given scattered data calculates (interpolates/
 *                  approximates/averages) values in centres of a specified
 *                  grid by using either:
 *                    -- approximation with bivariate C1 cubic spline
 *                    -- natural neighbours interpolation
 *                    -- linear interpolation
 *                    -- averaging
 *
 *  Description:    See info().
 *
 *  Revisions:      PS 28 Aug 2003: added option "-a 4" to average within
 *                    the grid cells.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <nn.h>
#include <csa.h>
#include <errno.h>
#include "nan.h"
#include "gridnodes.h"
#include "gridmap.h"
#include "gucommon.h"
#include "gridaverager.h"

#define PPE_DEF 3
#define PPE_MAX 10
#define ZMIN_DEF (-DBL_MAX)
#define ZMAX_DEF DBL_MAX

enum { CSA = 0, NN_SIBSON = 1, NN_NONSIBSONIAN = 2, LINEAR = 3, AVERAGE = 4 } rule = CSA;

int linear = 0;
int indexspace = 0;

static void version()
{
    printf("gridbathy/libgu version %s\nlibnn version %s\nlibcsa version %s\n", gu_version, nn_version, csa_version);
    exit(0);
}

static void usage()
{
    printf("Usage: gridbathy -b <bathymetry file> -g <grid file>\n");
    printf("                 [-a {0|1|2|3}]\n");
    printf("                 [-c <i> <j>]\n");
    printf("                 [-i <node type>]\n");
    printf("                 [-m <mask file>]\n");
    printf("                 [-n <points per edge>]\n");
    printf("                 [-r <min depth> <max depth>]\n");
    printf("                 [-v]\n");
    printf("                 [-x]\n");
    printf("Try \"gridbathy -h\" for more information\n");

    exit(0);
}

static void info()
{
    printf("Usage: gridbathy -b <bathymetry file> -g <grid file>\n");
    printf("                 [-a {0|1|2|3}]\n");
    printf("                 [-c <i> <j>]\n");
    printf("                 [-i <node type>]\n");
    printf("                 [-m <mask file>]\n");
    printf("                 [-n <points per edge>]\n");
    printf("                 [-r <min depth> <max depth>]\n");
    printf("                 [-v]\n");
    printf("                 [-x]\n");
    printf("Where:\n");
    printf("  -a 0                 -- use bivariate cubic spline approximation (default)\n");
    printf("  -a 1                 -- use non-Sibsonian Natural Neighbours interpolation\n");
    printf("  -a 2                 -- use Sibson Natural Neighbours interpolation\n");
    printf("  -a 3                 -- use linear interpolation\n");
    printf("  -a 4                 -- use averaging\n");
    printf("  -b <data file>       -- three-column (X Y Z) data file (use \"stdin\" or \"-\"\n");
    printf("                          for standard input)\n");
    printf("  -c <i> <j>           -- estimate depth for this cell only\n");
    printf("  -g <grid file>       -- two-column (X Y) grid file (see remarks below) (use\n");
    printf("                          \"stdin\" or \"-\" for standard input)\n");
    printf("  -i <node type>       -- grid node type; possible types:\n");
    printf("                          DD -- double density (default) \n");
    printf("                          CO -- cell corner\n");
    printf("  -m <mask file>       -- text file with nce1 x nce2 lines containing \"0\" or \"1\"\n");
    printf("                          (use \"stdin\" or \"-\" for standard input)\n");
    printf("  -n <points per edge> -- number of points per cell edge (default = 3)\n");
    printf("  -r <min> <max>       -- depth range (default = -infty +infty)\n");
    printf("  -v                   -- verbose / version\n");
    printf("  -x                   -- conduct interpolation in index space\n");
    printf("Description:\n");
    printf(" `gridbathy' interpolates scatterred 2D data into a grid. It interpolates input\n");
    printf("  data in a specified number of points within each \"valid\" cell (see below\n");
    printf("  for details) of the input grid using one of available schemes. After that,\n");
    printf("  it discards invalid points, replaces the ones that shot outside the allowed\n");
    printf("  range by the allowed min/max values and averages over the cell.\n");
    printf("Remarks:\n");
    printf("  1. The grid file must contain header describing the node array dimension:\n");
    printf("     ## <nx> x <ny>\n");
    printf("     where for double density nodes nx = nce1 * 2 + 1, ny = nce2 * 2 + 1;\n");
    printf("     for corner nodes  nx = nce1 + 1, ny = nce2 + 1; and for center nodes\n");
    printf("     nx = nce1, ny = nce2.\n");
    printf("  2. An empty or commented line in the input grid file as well as NaNs for\n");
    printf("     node coordinates indicate an invalid node.\n");
    printf("  3. A grid cell is valid if all corner nodes are valid (not NaNs).\n");
    printf("  4. The grid (union of all valid grid cells) must be simply connected.\n");
    printf("  5. The number of interpolations per cell is defined by -n option and is\n");
    printf("     equal to 3x3 = 9 by default. E.g., \"-n 4\" specifies 4x4 = 16\n");
    printf("     interpolations per cell; \"-n 1\" is equivalent to interpolation in cell\n");
    printf("     centers only.\n");
    printf("  6. An optional mask file is a file with nce1 x nce2 lines containing \"1\" for\n");
    printf("     valid cells and \"0\" for invalid cells.\n");

    exit(0);
}

static void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);
    fprintf(stderr, "error: gridbathy: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    exit(1);
}

static int str2double(char* token, double* value)
{
    char* end = NULL;

    if (token == NULL) {
        *value = NaN;
        return 0;
    }

    *value = strtod(token, &end);

    if (end == token) {
        *value = NaN;
        return 0;
    }

    return 1;
}

static int str2int(char* token, int* value)
{
    char* end = NULL;

    if (token == NULL) {
        *value = INT_MIN;
        return 0;
    }

    *value = strtol(token, &end, 10);

    if (end == token) {
        *value = INT_MIN;
        return 0;
    }

    return 1;
}

static void parse_commandline(int argc, char* argv[], char** bathyfname, char** gridfname, char** maskfname, NODETYPE* nt, int* ppe, double* zmin, double* zmax, int* ii, int* jj)
{
    int i;

    if (argc <= 1)
        usage();

    i = 1;
    while (i < argc) {
        switch (argv[i][1]) {
        case 'a':
            i++;
            if (i >= argc)
                quit("no interpolation rule number found after -a\n");
            rule = atoi(argv[i]);
            if (rule < CSA || rule > AVERAGE)
                quit("invalid interpolation rule number\n");
            i++;
            break;
        case 'b':
            i++;
            if (i >= argc)
                quit("no file name found after -b\n");
            *bathyfname = argv[i];
            i++;
            break;
        case 'c':
            i++;
            if (i >= argc)
                quit("no cell <i> index found after -c\n");
            *ii = atof(argv[i]);
            i++;
            if (i >= argc)
                quit("no cell <j> index found after -c\n");
            *jj = atof(argv[i]);
            i++;
            break;
        case 'g':
            i++;
            if (i >= argc)
                quit("no file name found after -g\n");
            *gridfname = argv[i];
            i++;
            break;
        case 'h':
            info();
            break;
        case 'i':
            i++;
            if (i == argc)
                quit("no node type found after \"-i\"\n");
            if (strcasecmp("dd", argv[i]) == 0)
                *nt = NT_DD;
            else if (strcasecmp("ce", argv[i]) == 0)
                *nt = NT_CEN;
            else if (strcasecmp("co", argv[i]) == 0)
                *nt = NT_COR;
            else
                quit("input node type \"%s\" not recognised\n", argv[i]);
            i++;
            break;
        case 'm':
            i++;
            if (i == argc)
                quit("no file name found after \"-m\"\n");
            *maskfname = argv[i];
            i++;
            break;
        case 'n':
            i++;
            if (i >= argc)
                quit("no points per edge value found after -n\n");
            if (!str2int(argv[i], ppe))
                quit("error: could not read points per edge after -n\n");
            i++;
            break;
        case 'r':
            i++;
            if (i >= argc)
                quit("no min depth found after -r\n");
            if (!str2double(argv[i], zmin))
                quit("error: could not read min depth after -r\n");
            i++;
            if (i >= argc)
                quit("no max depth found after -r\n");
            if (!str2double(argv[i], zmax))
                quit("could not read max depth after -r\n");
            i++;
            break;
        case 'v':
            i++;
            gu_verbose = 1;
            break;
        case 'x':
            i++;
            indexspace = 1;
            break;
        default:
            usage();
            break;
        }
    }

    if (gu_verbose && argc == 2)
        version();
}

int main(int argc, char* argv[])
{
    char* bathyfname = NULL;
    int nbathy = -1;
    point* pbathy = NULL;

    char* gridfname = NULL;
    NODETYPE nt = NT_DD;
    gridnodes* gn = NULL;
    gridmap* gm = NULL;

    char* maskfname = NULL;
    int** mask = NULL;

    int ppe = PPE_DEF;
    double zmin = ZMIN_DEF;
    double zmax = ZMAX_DEF;

    delaunay* d = NULL;

    void (*interpolate_point) (void*, point *) = NULL;
    void* interpolator = NULL;

    int i = -1, j = -1;

    parse_commandline(argc, argv, &bathyfname, &gridfname, &maskfname, &nt, &ppe, &zmin, &zmax, &i, &j);

    /*
     * sanity check 
     */
    if (bathyfname == NULL)
        quit("no input bathymetry data specified\n");
    if (gridfname == NULL)
        quit("no input grid data specified\n");
    if (ppe <= 0 || ppe > PPE_MAX)
        quit("number of points per edge specified = %d greater than %d\n", ppe, PPE_MAX);
    if (zmin >= zmax)
        quit("min depth = %.3g > max depth = %.3g\n", zmin, zmax);
    if (nt != NT_DD && nt != NT_COR)
        quit("unsupported node type\n");

    /*
     * read bathymetry 
     */
    points_read(bathyfname, 3, &nbathy, &pbathy);

    if (gu_verbose)
        fprintf(stderr, "## %d input bathymetry values\n", nbathy);
    if (nbathy < 3)
        quit("less than 3 input bathymetry values\n");

    /*
     * read and validate grid 
     */
    gn = gridnodes_read(gridfname, nt);
    gridnodes_validate(gn);
    /*
     * read mask
     */
    if (maskfname != NULL) {
        int nx = gridnodes_getnce1(gn);
        int ny = gridnodes_getnce2(gn);

        mask = gu_readmask(maskfname, nx, ny);
    }
    /*
     * transform grid nodes to corner type 
     */
    if (nt != NT_COR) {
        gridnodes* newgn = gridnodes_transform(gn, NT_COR);

        gridnodes_destroy(gn);
        gn = newgn;
    }
    /*
     * build the grid map for physical <-> index space conversions
     */
    gm = gridmap_build(gridnodes_getnce1(gn), gridnodes_getnce2(gn), gridnodes_getx(gn), gridnodes_gety(gn));

    /*
     * convert bathymetry to index space if necessary 
     */
    if (indexspace) {
        point* newpbathy = malloc(nbathy * sizeof(point));
        int newnbathy = 0;
        int ii;

        for (ii = 0; ii < nbathy; ++ii) {
            point* p = &pbathy[ii];
            point* newp = &newpbathy[newnbathy];
            double ic, jc;

            if (gridmap_xy2fij(gm, p->x, p->y, &ic, &jc)) {
                newp->x = ic;
                newp->y = jc;
                newp->z = p->z;
                newnbathy++;
            }
        }

        free(pbathy);
        pbathy = newpbathy;
        nbathy = newnbathy;
    }

    /*
     * create interpolator 
     */
    if (rule == CSA) {          /* using libcsa */
        interpolator = csa_create();
        csa_addpoints(interpolator, nbathy, pbathy);
        csa_calculatespline(interpolator);
        interpolate_point = (void (*)(void*, point *)) csa_approximatepoint;
    } else if (rule == AVERAGE) {
        interpolator = ga_create(gm);
        ga_addpoints(interpolator, nbathy, pbathy);
        interpolate_point = (void (*)(void*, point *)) ga_getvalue;
        ppe = 1;
    } else {                    /* using libnn */
        /*
         * triangulate 
         */
        if (gu_verbose) {
            fprintf(stderr, "## triangulating...");
            fflush(stdout);
        }
        d = delaunay_build(nbathy, pbathy, 0, NULL, 0, NULL);
        if (gu_verbose) {
            fprintf(stderr, "done\n");
            fflush(stderr);
        }

        if (rule == NN_SIBSON || rule == NN_NONSIBSONIAN) {
            interpolator = nnpi_create(d);
            if (rule == NN_SIBSON)
                nn_rule = SIBSON;
            else
                nn_rule = NON_SIBSONIAN;
            interpolate_point = (void (*)(void*, point *)) nnpi_interpolate_point;
        } else if (rule == LINEAR) {
            interpolator = lpi_build(d);
            interpolate_point = (void (*)(void*, point *)) lpi_interpolate_point;
        }
    }

    /*
     * main cycle -- over grid cells 
     */
    {
        double** gx = gridnodes_getx(gn);
        int jmin, jmax, imin, imax;

        if (i < 0) {
            imin = 0;
            imax = gridnodes_getnce1(gn) - 1;
            jmin = 0;
            jmax = gridnodes_getnce2(gn) - 1;
        } else {
            if (gu_verbose)
                fprintf(stderr, "## calculating depth for cell (%d,%d)\n", i, j);
            imin = i;
            imax = i;
            jmin = j;
            jmax = j;
        }

        for (j = jmin; j <= jmax; ++j) {
            for (i = imin; i <= imax; ++i) {
                double sum = 0.0;
                int count = 0;
                int ii, jj;

                if ((mask != NULL && mask[j][i] == 0) || isnan(gx[j][i]) || isnan(gx[j + 1][i + 1]) || isnan(gx[j][i + 1]) || isnan(gx[j + 1][i + 1])) {
                    printf("NaN\n");
                    continue;
                }

                for (ii = 0; ii < ppe; ++ii) {
                    for (jj = 0; jj < ppe; ++jj) {
                        double fi = (double) i + 0.5 / (double) ppe * (1.0 + 2.0 * (double) ii);
                        double fj = (double) j + 0.5 / (double) ppe * (1.0 + 2.0 * (double) jj);
                        point p;

                        if (!indexspace)
                            gridmap_fij2xy(gm, fi, fj, &p.x, &p.y);
                        else {
                            p.x = fi;
                            p.y = fj;
                        }

                        interpolate_point(interpolator, &p);

                        if (isnan(p.z))
                            continue;
                        else if (p.z < zmin)
                            p.z = zmin;
                        else if (p.z > zmax)
                            p.z = zmax;

                        sum += p.z;
                        count++;
                    }
                }

                if (count == 0)
                    printf("NaN\n");
                else
                    printf("%.2f\n", sum / (double) count);
                fflush(stdout);
            }
        }
    }

    /*
     * clean up, just because 
     */
    if (rule == CSA)
        csa_destroy(interpolator);
    else if (rule == AVERAGE)
        ga_destroy(interpolator);
    else {
        if (rule == NN_SIBSON || rule == NN_NONSIBSONIAN)
            nnpi_destroy(interpolator);
        else if (rule == LINEAR)
            lpi_destroy(interpolator);
        delaunay_destroy(d);
    }
    if (mask != NULL)
        gu_free2d(mask);
    gridmap_destroy(gm);
    gridnodes_destroy(gn);
    free(pbathy);

    return 0;
}
