/******************************************************************************
 *
 * File:           nnbathy.c
 *
 * Created:        04/08/2000
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Interpolate scalar 2D data in specified points using
 *                 Natural Neighbours interpolation.
 *
 * Description:    See usage().
 *
 *                 The default version of nnbathy allocates the whole output
 *                 grid in memory before interpolating; the compilier flag
 *                 NN_SERIAL compiles a bit more sophisticated version of that
 *                 interpolates in output points on one-by-one basis.
 *
 * Revisions:      29/05/2006: introduced NN_SERIAL, see the description above
 *                 01/06/2006: moved flags and other command-line input into
 *                             structure "specs"
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <errno.h>
#include "config.h"
#include "nan.h"
#include "minell.h"
#include "nn.h"
#if defined(NN_SERIAL)
#include "preader.h"
#endif

#if !defined(NN_SERIAL)
#define NMAX 4096
#endif
#define STRBUFSIZE 64

typedef struct {
    int generate_points;
    int thin;
    int nointerp;
    int linear;
    int invariant;
    int square;
    char* fin;
    char* fout;
    int nx;
    int ny;
    int nxd;
    int nyd;
    double rmax;
    double wmin;
    double zoom;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
#if defined(NN_SERIAL)
    int npoints;
#endif
} specs;

static void version()
{
    printf("  nnbathy/nn version %s\n", nn_version);
    exit(0);
}

static void usage()
{
    printf("Usage: nnbathy -i <XYZ file>\n");
    printf("               -o <XY file> | -n <nx>x<ny> [-c|-s] [-z <zoom>]\n");
    printf("               [-x <xmin xmax>] [-xmin <xmin>] [-xmax <xmax>]\n");
    printf("               [-y <ymin ymax>] [-ymin <ymin>] [-ymax <ymax>]\n");
    printf("               [-v|-T <vertex id>|-V]\n");
    printf("               [-D [<nx>x<ny>]]\n");
    printf("               [-L <dist>]\n");
    printf("               [-N]\n");
    printf("               [-P alg={l|nn|ns}]\n");
    printf("               [-W <min weight>]\n");
#if defined(NN_SERIAL)
    printf("               [-%% [npoints]]\n");
#endif
    printf("Options:\n");
    printf("  -c               -- scale internally so that the enclosing minimal ellipse\n");
    printf("                      turns into a circle (this produces results invariant to\n");
    printf("                      affine transformations)\n");
    printf("  -i <XYZ file>    -- three-column file with points to interpolate from\n");
    printf("                      (use \"-i stdin\" or \"-i -\" for standard input)\n");
    printf("  -n <nx>x<ny>     -- generate <nx>x<ny> output rectangular grid\n");
    printf("  -o <XY file>     -- two-column file with points to interpolate in\n");
    printf("                      (use \"-o stdin\" or \"-o -\" for standard input)\n");
    printf("  -s               -- scale internally so that Xmax - Xmin = Ymax - Ymin\n");
    printf("  -x <xmin> <xmax> -- set Xmin and Xmax for the output grid\n");
    printf("  -xmin <xmin>     -- set Xmin for the output grid\n");
    printf("  -xmax <xmax>     -- set Xmin for the output grid\n");
    printf("  -y <ymin> <ymax> -- set Ymin and Ymax for the output grid\n");
    printf("  -ymin <ymin>     -- set Ymin for the output grid\n");
    printf("  -ymax <ymax>     -- set Ymin for the output grid\n");
    printf("  -v               -- verbose / version\n");
    printf("  -z <zoom>        -- zoom in (if <zoom> < 1) or out (<zoom> > 1) (activated\n");
    printf("                      only when used in conjunction with -n)\n");
    printf("  -D [<nx>x<ny>]   -- thin input data by averaging X, Y and Z values within\n");
    printf("                      every cell of the rectangular <nx>x<ny> grid (size\n");
    printf("                      optional with -n)\n");
    printf("  -L <dist>        -- thin input data by averaging X, Y and Z values within\n");
    printf("                      clusters of consequitive input points such that the\n");
    printf("                      sum of distances between points within each cluster\n");
    printf("                      does not exceed the specified maximum value\n");
    printf("  -N               -- do not interpolate, only pre-process\n");
    printf("  -P alg=<l|nn|ns> -- use the following algorithm:\n");
    printf("                        l -- linear interpolation\n");
    printf("                        nn -- Sibson interpolation (default)\n");
    printf("                        ns -- Non-Sibsonian interpolation\n");
    printf("  -T <vertex id>   -- verbose; in weights output print weights associated\n");
    printf("                      with this vertex only\n");
    printf("  -V               -- very verbose / version\n");
    printf("  -W <min weight>  -- restricts extrapolation by assigning minimal allowed\n");
    printf("                      weight for a vertex (normally \"-1\" or so; lower\n");
    printf("                      values correspond to lower reliability; \"0\" means\n");
    printf("                      no extrapolation)\n");
#if defined(NN_SERIAL)
    printf("  -%% [npoints]     -- print percent of the work done to standard error;\n");
    printf("                      npoints -- total number of points to be done (optional\n");
    printf("                      with -n)\n");
#endif
    printf("Description:\n");
    printf("  `nnbathy' interpolates scalar 2D data in specified points using Natural\n");
    printf("  Neighbours interpolation. The interpolated values are written to standard\n");
    printf("  output.\n");

    exit(0);
}

static void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);             /* just in case, to have exit message last */

    fprintf(stderr, "  error: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    exit(1);
}

static double str2double(char* token, char* option)
{
    char* end = NULL;
    double value = NaN;

    if (token != NULL)
        value = strtod(token, &end);

    if (token == NULL || end == token) {
        fprintf(stderr, "  error: command-line option \"%s\": could not convert \"%s\" to double\n", option, (token != NULL) ? token : "NULL");

        exit(1);
    }

    return value;
}

static specs* specs_create(void)
{
    specs* s = malloc(sizeof(specs));

    s->generate_points = 0;
    s->thin = 0;
    s->nointerp = 0;
    s->linear = 0;
    s->invariant = 0;
    s->square = 0;
    s->fin = NULL;
    s->fout = NULL;
    s->nx = -1;
    s->ny = -1;
    s->nxd = -1;
    s->nyd = -1;
    s->rmax = NaN;
    s->wmin = -DBL_MAX;
    s->zoom = 1.0;
    s->xmin = NaN;
    s->xmax = NaN;
    s->ymin = NaN;
    s->ymax = NaN;
#if defined(NN_SERIAL)
    s->npoints = 0;
#endif
    return s;
}

void specs_destroy(specs * s)
{
    free(s);
}

static void parse_commandline(int argc, char* argv[], specs * s)
{
    int i;

    if (argc < 2)
        usage();

    i = 1;
    while (i < argc) {
        if (argv[i][0] != '-')
            usage();

        switch (argv[i][1]) {
        case 'c':
            i++;
            s->square = 0;
            s->invariant = 1;
            break;
        case 'i':
            i++;
            if (i >= argc)
                quit("no file name found after -i\n");
            s->fin = argv[i];
            i++;
            break;
        case 'l':
            i++;
            s->linear = 1;
            break;
        case 'n':
            i++;
            s->fout = NULL;
            s->generate_points = 1;
            if (i >= argc)
                quit("no grid dimensions found after -n\n");
            if (sscanf(argv[i], "%dx%d", &s->nx, &s->ny) != 2)
                quit("could not read grid dimensions after \"-n\"\n");
#if defined(NN_SERIAL)
            if (s->nx <= 0 || s->ny <= 0)
#else
            if (s->nx <= 0 || s->nx > NMAX || s->ny <= 0 || s->ny > NMAX)
#endif
                quit("invalid size for output grid\n");
            i++;
            break;
        case 'o':
            i++;
            if (i >= argc)
                quit("no file name found after -o\n");
            s->fout = argv[i];
            i++;
            break;
        case 's':
            i++;
            s->square = 1;
            s->invariant = 0;
            break;
        case 'x':
            if (argv[i][2] == 0) {
                i++;
                if (i >= argc)
                    quit("no xmin value found after -x\n");
                s->xmin = str2double(argv[i], "-x");
                i++;
                if (i >= argc)
                    quit("no xmax value found after -x\n");
                s->xmax = str2double(argv[i], "-x");
                i++;
            } else if (strcmp(argv[i], "-xmin") == 0) {
                i++;
                if (i >= argc)
                    quit("no value found after -xmin\n");
                s->xmin = str2double(argv[i], "-xmin");
                i++;
            } else if (strcmp(argv[i], "-xmax") == 0) {
                i++;
                if (i >= argc)
                    quit("no value found after -xmax\n");
                s->xmax = str2double(argv[i], "-xmax");
                i++;
            } else
                usage();
            break;
        case 'y':
            if (argv[i][2] == 0) {
                i++;
                if (i >= argc)
                    quit("no ymin value found after -y\n");
                s->ymin = str2double(argv[i], "-y");
                i++;
                if (i >= argc)
                    quit("no ymax value found after -y\n");
                s->ymax = str2double(argv[i], "-y");
                i++;
            } else if (strcmp(argv[i], "-ymin") == 0) {
                i++;
                if (i >= argc)
                    quit("no value found after -ymin\n");
                s->ymin = str2double(argv[i], "-ymin");
                i++;
            } else if (strcmp(argv[i], "-ymax") == 0) {
                i++;
                if (i >= argc)
                    quit("no value found after -ymax\n");
                s->ymax = str2double(argv[i], "-ymax");
                i++;
            } else
                usage();
            break;
        case 'v':
            i++;
            nn_verbose = 1;
            break;
        case 'z':
            i++;
            if (i >= argc)
                quit("no zoom value found after -z\n");
            s->zoom = str2double(argv[i], "-z");
            i++;
            break;
        case 'D':
            i++;
            s->thin = 1;
            if (argc > i && argv[i][0] != '-') {
                if (sscanf(argv[i], "%dx%d", &s->nxd, &s->nyd) != 2)
                    quit("could not read grid dimensions after \"-D\"\n");
                if (s->nxd <= 0 || s->nyd <= 0)
                    quit("invalid value for nx = %d or ny = %d after -D option\n", s->nxd, s->nyd);
#if !defined(NN_SERIAL)
                if (s->nxd > NMAX || s->nyd > NMAX)
                    quit("too big value after -D option (expected < %d)\n", NMAX);
#endif
                i++;
            }
            break;
        case 'L':
            i++;
            s->thin = 2;
            if (i >= argc)
                quit("no value found after -L\n");
            s->rmax = str2double(argv[i], "-L");
            i++;
            break;
        case 'N':
            i++;
            s->nointerp = 1;
            break;
        case 'P':{
                char delim[] = "=";
                char prmstr[STRBUFSIZE] = "";
                char* token;

                i++;
                if (i >= argc)
                    quit("no input found after -P\n");

                if (strlen(argv[i]) >= STRBUFSIZE)
                    quit("could not interpret \"%s\" after -P option\n", argv[i]);

                strcpy(prmstr, argv[i]);
                token = strtok(prmstr, delim);
                if (token == NULL)
                    quit("could not interpret \"%s\" after -P option\n", argv[i]);

                if (strcmp(token, "alg") == 0) {
                    token = strtok(NULL, delim);
                    if (token == NULL)
                        quit("could not interpret \"%s\" after -P option\n", argv[i]);

                    if (strcmp(token, "nn") == 0) {
                        nn_rule = SIBSON;
                        s->linear = 0;
                    } else if (strcmp(token, "ns") == 0) {
                        nn_rule = NON_SIBSONIAN;
                        s->linear = 0;
                    } else if (strcmp(token, "l") == 0) {
                        s->linear = 1;
                    } else
                        usage();
                }

                i++;
                break;
            }
        case 'W':
            i++;
            if (i >= argc)
                quit("no minimal allowed weight found after -W\n");
            s->wmin = str2double(argv[i], "-W");
            i++;
            break;
        case 'T':
            i++;
            if (i >= argc)
                quit("no vertex id found after -T\n");
            nn_test_vertice = atoi(argv[i]);
            nn_verbose = 1;
            i++;
            break;
        case 'V':
            i++;
            nn_verbose = 2;
            break;
#if defined(NN_SERIAL)
        case '%':
            i++;
            if (i < argc && argv[i][0] != '-') {
                s->npoints = atoi(argv[i]);
                i++;
            } else
                s->npoints = 1;
            break;
#endif
        default:
            usage();
            break;
        }
    }

    if (nn_verbose && argc == 2)
        version();

    if (s->thin) {
        if (s->nxd == -1)
            s->nxd = s->nx;
        if (s->nyd == -1)
            s->nyd = s->ny;
        if (s->nxd <= 0 || s->nyd <= 0)
            quit("invalid grid size for thinning\n");
    }
#if defined(NN_SERIAL)
    if (s->npoints == 1) {
        if (s->nx <= 0)
            s->npoints = 0;
        else
            s->npoints = s->nx * s->ny;
    }
#endif
}

static void points_write(int n, point* points)
{
    int i;

    for (i = 0; i < n; ++i) {
        point* p = &points[i];

        if (isnan(p->z))
            printf("%.15g %.15g NaN\n", p->x, p->y);
        else
            printf("%.15g %.15g %.15g\n", p->x, p->y, p->z);
    }
}

#if !defined(NN_SERIAL)
/* A simpler version of nnbathy that allocates the whole output grid in memory
 */
int main(int argc, char* argv[])
{
    specs* s = specs_create();
    int nin = 0;
    point* pin = NULL;
    minell* me = NULL;
    int nout = 0;
    point* pout = NULL;
    double k = NaN;

    parse_commandline(argc, argv, s);

    if (s->fin == NULL)
        quit("no input data\n");

    if (!s->generate_points && s->fout == NULL && !s->nointerp)
        quit("no output grid specified\n");

    points_read(s->fin, 3, &nin, &pin);

    if (nin < 3)
        return 0;

    if (s->thin == 1)
        points_thingrid(&nin, &pin, s->nxd, s->nyd);
    else if (s->thin == 2)
        points_thinlin(&nin, &pin, s->rmax);

    if (s->nointerp) {
        points_write(nin, pin);
        specs_destroy(s);
        free(pin);
        return 0;
    }

    if (s->generate_points) {
        /*
         * points_getrange() only writes the proper values to those arguments
         * which do not point to NaNs 
         */
        points_getrange(nin, pin, s->zoom, &s->xmin, &s->xmax, &s->ymin, &s->ymax);
        points_generate(s->xmin, s->xmax, s->ymin, s->ymax, s->nx, s->ny, &nout, &pout);
    } else
        points_read(s->fout, 2, &nout, &pout);

    if (s->invariant) {
        me = minell_build(nin, pin);
        minell_scalepoints(me, nin, pin);
        minell_scalepoints(me, nout, pout);
    } else if (s->square) {
        k = points_scaletosquare(nin, pin);
        points_scale(nout, pout, k);
    }

    if (s->linear)
        lpi_interpolate_points(nin, pin, nout, pout);
    else
        nnpi_interpolate_points(nin, pin, s->wmin, nout, pout);

    if (s->invariant)
        minell_rescalepoints(me, nout, pout);
    else if (s->square)
        points_scale(nout, pout, 1.0 / k);

    points_write(nout, pout);

    if (me != NULL)
        minell_destroy(me);
    specs_destroy(s);
    free(pin);
    free(pout);

    return 0;
}
#else                           /* NN_SERIAL */
/* A version of nnbathy that interpolates output points serially. Can save a
 * bit of memory for large output grids.
 */
int main(int argc, char* argv[])
{
    specs* s = specs_create();
    int nin = 0;
    point* pin = NULL;
    minell* me = NULL;
    point* pout = NULL;
    double k = NaN;
    preader* pr = NULL;
    delaunay* d = NULL;
    void* interpolator = NULL;
    int ndone = 0;

    parse_commandline(argc, argv, s);

    if (s->fin == NULL)
        quit("no input data\n");

    if (!s->generate_points && s->fout == NULL && !s->nointerp)
        quit("no output grid specified\n");

    points_read(s->fin, 3, &nin, &pin);

    if (nin < 3)
        return 0;

    if (s->thin == 1)
        points_thingrid(&nin, &pin, s->nxd, s->nyd);
    else if (s->thin == 2)
        points_thinlin(&nin, &pin, s->rmax);

    if (s->nointerp) {
        points_write(nin, pin);
        specs_destroy(s);
        free(pin);
        return 0;
    }

    if (s->generate_points) {
        points_getrange(nin, pin, s->zoom, &s->xmin, &s->xmax, &s->ymin, &s->ymax);
        pr = preader_create1(s->xmin, s->xmax, s->ymin, s->ymax, s->nx, s->ny);
    } else
        pr = preader_create2(s->fout);

    if (s->invariant) {
        me = minell_build(nin, pin);
        minell_scalepoints(me, nin, pin);
    } else if (s->square)
        k = points_scaletosquare(nin, pin);

    d = delaunay_build(nin, pin, 0, NULL, 0, NULL);

    if (s->linear)
        interpolator = lpi_build(d);
    else {
        interpolator = nnpi_create(d);
        nnpi_setwmin(interpolator, s->wmin);
    }

    while ((pout = preader_getpoint(pr)) != NULL) {
        if (s->invariant)
            minell_scalepoints(me, 1, pout);
        else if (s->square)
            points_scale(1, pout, k);

        if (s->linear)
            lpi_interpolate_point(interpolator, pout);
        else
            nnpi_interpolate_point(interpolator, pout);

        if (s->invariant)
            minell_rescalepoints(me, 1, pout);
        else if (s->square)
            points_scale(1, pout, 1.0 / k);

        points_write(1, pout);
        ndone++;

        if (s->npoints > 0)
            fprintf(stderr, "  %5.2f%% done\r", 100.0 * ndone / s->npoints);
    }
    if (s->npoints > 0)
        fprintf(stderr, "                \r");

    if (me != NULL)
        minell_destroy(me);
    if (s->linear)
        lpi_destroy(interpolator);
    else
        nnpi_destroy(interpolator);
    delaunay_destroy(d);
    preader_destroy(pr);
    specs_destroy(s);
    free(pin);

    return 0;
}
#endif
