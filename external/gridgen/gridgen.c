/******************************************************************************
 *
 * File:           gridgen.c
 *
 * Created:        19/03/2001
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Code for generation of 2D orthogonal grids.
 *
 * Description:    This program is based on the CRDT algorithm as described
 *                 in:
 *                 Tobin D. Driscoll and Stephen A. Vavasis, "Numerical
 *                 conformal mapping using cross-ratios and Delaunay
 *                 triangulation", SIAM J. Sci. Comp., 1998, 19(6), 1783-1803.
 *                 See http://www.math.udel.edu/~driscoll/pubs/crdt.pdf.
 *                 
 * Revisions:      02/11/2006 PS Minor structural changes associated with
 *                            the introduction of libgridgen.a. Changed the
 *                            structure name "gridspec" to "gridgen".
 *                 15/02/2007 PS Introduced a new function gridgen_create2().
 *                            Some associated structural changes.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include "gridgen.h"
#include "ode.h"
#include "istack.h"
#include "delaunay.h"
#include "vertlist.h"
#include "swcr.h"
#include "broyden.h"
#include "geom.h"
#include "issimplepoly.h"
#include "hash.h"
#include "version.h"
#include "config.h"
#include "nan.h"

#define BUFSIZE 1024
#define NMIN 3                  /* minimal number of vertices */
#define NMAX 10000              /* maximal number of vertices */
#define M_2PI (M_PI * 2.0)
#define EPS 1.0e-6
#define EPS_DEF 1.0e-12         /* default precision */
#define EPS_MAX 1.0e-3
#define EPS_MIN 1.0e-14
#define NNODES_MIN 5
#define NNODES_MAX 20
#define NNODES_DEF 10           /* default number of nodes */
#define N_DEF 25
#define N_MIN 2                 /* min. number of nodes in grid edge */
#define N_MAX 10001             /* max. number of nodes in grid edge */
#define BIGDOUBLE 1.0e+6
#define THIN_DEF 1              /* thin input vertices by default */
#define CHECKSIMPLEPOLY_DEF 1   /* check input for self-intersections */
#define NEWTON_DEF 1            /* newton method by default */
#define M 4                     /* number of stored iterations in generalised
                                 * Broyden update */
#define NPPE_DEF 3              /* number of points per internal edge in
                                 * polygon images */
#define ZZERO 0.0 + 0.0 * I

int gg_verbose = 0;
int tr_verbose = 0;

/* The diagonal goes from vids[0] to vids[2];
 * tids[0] corresponds to vids[0], vids[1] and vids[2];
 * tids[1] corresponds to vids[4], vids[2], vids[3].
 */
typedef struct {
    int id;
    int vids[4];
    int tids[2];

    /*
     * quadrilaterals which share with this one one of the triangles 
     */
    int nneighbours;
    int* neighbours;
} quadrilateral;

typedef struct {
    char* prmfname;
    char* datafname;
    char* sigmafname;
    char* rectfname;
    FILE* out;
#if !defined (GRIDGEN_STANDALONE) && defined(HAVE_GRIDNODES_H)
    gridnodes* gn;
#endif
    FILE* fsigma;
    FILE* frect;
    int* nrect;
    double** xrect;
    double** yrect;

    double xmin;
    double ymin;
    double xmax;
    double ymax;

    int thin;
    int checksimplepoly;
    int reversed;

    int nx;
    int ny;

    int ngridpoints;
    zdouble* gridpoints;

    /*
     * specifications 
     */
    double eps;
    int nnodes;
    int newton;

    /*
     * z 
     */
    vertlist* vertices;
    delaunay* d;
    int nquadrilaterals;
    quadrilateral* quadrilaterals;
    int lastqi;

    int nz;
    zdouble* zs;

    /*
     * z<->w 
     */
    double* betas;
    swcr* sc;
    double* cis;
    int* nsigmas;
    double** sigmas;
    zdouble* ws;
    zdouble* As;
    zdouble* Bs;

    int ncorners;               /* number of corners in the mapped polygonal; 
                                 * 4 by default */

    /*
     * z1 
     */
    zdouble* rzs;
    int* rvids;

    double newxmin;
    double newymin;
    double newxmax;
    double newymax;

    /*
     * w<->z1 
     */
    double* newbetas;
    swcr* newsc;
    zdouble* newAs;
    zdouble* newBs;
    zdouble* newzs;
    zdouble* newrzs;

    /*
     * What we have is a number of mappings, each working best in a certain
     * quadrilateral of the complex plane of the original polygon. It looses
     * a bit of precision for the neighbour quadrilaterals (but still works 
     * OK for them), and can deteriorate substantially or completely for points
     * in distant quadrilaterals.
     *
     * Now, let one map a point from the complex plane of the image to the 
     * complex plane of the original. For best results, one needs to find the
     * best mapping to be used. To find it, one needs to know which
     * quadrilateral in the original plane the inverse transformation of a 
     * point in the image belongs to.
     *
     * The problem is that one does not know the images of the quadrilaterals, 
     * only the vertice images.
     *
     * To get coarse images of the quadrilaterals, we put gg->nppe points on 
     * the diagonal of each quadrilateral, and map them to the image region. 
     * The original quadrilaterals would then effectively transform into a
     * number of non-intersecting polygons in the image plane.
     *
     * (Each interior edge of the triangulation is a diagonal. By mapping
     * all diagonals, we map all interior edges of all quadrilaterals. In
     * this way, we ensure that each internal edge is mapped only once.
     * We then need to be able to get effectively the diagonal id for a given
     * edge of a quadrilateral. This is done via hashtable that allows to get 
     * the diagonal id from the edge vertex ids.)
     *
     * After that, almost any point in the image space can be easily mapped to
     * a certain quadrilateral. In rare cases when we will get a neighbour of
     * the right quadrilateral instead of the quadrilateral itself (because of
     * approximating the boundary with a polyline) this still will be good
     * enough to map the point.
     *
     * The mapping is controled by parameter `nppe'. To switch it off, set
     * nppe = 0 in the parameter file. In this case, coarser images of 
     * quadrilaterals will be used that are formed by straight lines
     * connecting the quadrilateral vertex images.
     */

    /*
     * quasi-images of quadrilaterals
     */
    int nppe;                   /* number of points per internal edge */
    int nppq;                   /* number of points allocated per
                                 * quadrilateral (= nppe * 4 + 5) */
    int* nqivertices;           /* number of vertices in this image
                                 * [nquadrilaterals] */
    zdouble* qivertices;        /* vertices [nquadrilaterals * nppq] */

    /*
     * id of quadrilateral containing the last valid grid point
     */
    int lastqid;

    /*
     * temporal stuff -- some output of interest to me. Ignore it. 
     */
    int diagonal;
} gridgen;

void gridgen_setverbose(int verbose)
{
    if (verbose <= 0) {
        gg_verbose = 0;
        tr_verbose = 0;
    } else if (verbose == 1) {
        gg_verbose = 1;
        tr_verbose = 0;
    } else if (verbose == 2) {
        gg_verbose = 2;
        tr_verbose = 1;
    } else if (verbose >= 3) {
        gg_verbose = 2;
        tr_verbose = 2;
    }
}

void gridgen_printversion(void)
{
    printf("gridgen version %s\n", gridgen_version);
}

void gridgen_printhelpalg(void)
{
    printf("  `gridgen' generates quasi-orthogonal grids for polygonal regions by CRDT\n");
    printf("algorithm. It is based on the Schwarz-Christoffel formula that allows to\n");
    printf("transform conformally a unit circumcircle into a simply connected polygon with\n");
    printf("specified vertex angles.\n");
    printf("  In the transformation defined by Schwarz-Christoffel formula, each vertex of\n");
    printf("the polygon corresponds to a point on the circumcircle (\"prevertex\"). By\n");
    printf("specifying different angles, it is possible to transform the same set of\n");
    printf("prevertices into different polygons. A conjunction of one forward and one\n");
    printf("backward transformation with the same prevertices but different vertex angles in\n");
    printf("the Schwarz-Christoffel formula yields a single conformal mapping between the\n");
    printf("original polygon and a new (\"image\") polygon.\n");
    printf("  The underlying idea for using conformal mappings is that these mappings do\n");
    printf("preserve local angles. Therefore, if the image polygon is set to be a rectangle\n");
    printf("or a set of intersecting rectangles, one can easily put an orthogonal grid on it\n");
    printf("and ...provided he can do the inverse transform... this orthogonal grid would\n");
    printf("yield a (quasi)-orthogonal grid in the original polygon.\n");
    printf("  The main obstacle to this idyllic picture is a potential loss of precision on\n");
    printf("the way, when a number of vertices degenerate into a single prevertex\n");
    printf("(\"crowding\"). The effect of crowding is particularly strong for polygons with\n");
    printf("long (thin) arms, and this is where the CRDT algorithm comes into play.\n");
    printf("  The essence of the CRDT algorithm is in splitting the original polygon into a\n");
    printf("number of intersecting quadrilaterals and generating not one but a number of\n");
    printf("conformally equivalent transformations (\"embeddings\") of the original polygon\n");
    printf("to the unit circumcircle, each of them attributed with a certain quadrilateral\n");
    printf("and working numerically best for this and adjacent quadrilaterals.\n");
    printf("  All these details are transparent to a user of `gridgen', who has to define\n");
    printf("only the original and the image polygons. The image polygon is defined in a\n");
    printf("semi-implicit manner, by specifying new angles (but not edge lengths) for each\n");
    printf("vertex. In the simplest case, a user must specify which four vertices are to\n");
    printf("become the vertices of the image rectangle by marking them with \"1\" in the\n");
    printf("input file specifying the original polygon. After that, `gridgen' will calculate\n");
    printf("images of the nodes of a rectangular grid with specified parameters (or images\n");
    printf("of a custom set of points) within the image polygon in the original polygon, and\n");
    printf("this will complete the task.\n");
    printf("  There are 4 main stages in generating a quasi-orthogonal grid:\n");
    printf("  -- preprocessing of the input polygon by removing duplicate vertices and\n");
    printf("     inserting new ones when necessary (to avoid thin triangles);\n");
    printf("  -- calculation of the conformal mappings that transfer the input polygon\n");
    printf("     into a unit circumcircle;\n");
    printf("  -- calculation of conformal mappings of a unit circumcircle into a polygon\n");
    printf("     that transfer the prevertices into vertices of a polygon with specified\n");
    printf("     vertex angles;\n");
    printf("  -- inverse mapping of specified points (grid nodes) from the image polygon\n");
    printf("     into the original one.\n");
}

void gridgen_printhelpprm(void)
{
    printf("  `gridgen': Generates grid for a polygonal region by CRDT algorithm.\n");
    printf("  Parameter file format:\n");
    printf("    input <input polygon data file>\n");
    printf("    [output <output grid data file>]\n");
    printf("    [grid <input grid data file>]\n");
    printf("    [nx <number of nodes in X direction>] (default = 25)\n");
    printf("    [ny <number of nodes in Y direction>] (25)\n");
    printf("    [nnodes <number of nodes in Gauss-Jacobi quadrature>] (12)\n");
    printf("    [precision <precision>] (1e-10)\n");
    printf("    [thin {0|1}] (1)\n");
    printf("    [checksimplepoly {0|1}] (1)\n");
    printf("    [newton {0|1}] (1)\n");
    printf("    [sigmas <intermediate backup file>]\n");
    printf("    [rectangle <output image polygon file>]\n");
    printf("    [nppe <number of points per internal edge>] (3)\n");
    printf("  Input polygon data file format:\n");
    printf("    <x_1> <y_1> [<beta_1>[*]]\n");
    printf("    <x_2> <y_2> [<beta_2>[*]]\n");
    printf("    ...\n");
    printf("    <x_n> <y_n> [<beta_n>[*]]\n");
    printf("  Input grid data file format:\n");
    printf("    <x_1> <y_1>\n");
    printf("    ...\n");
    printf("    <x_n> <y_n>\n");
    printf("    where 0 <= x_i, y_i <= 1\n");
    printf("  Remarks:\n");
    printf("    1. The input polygon must be simply connected.\n");
    printf("    2. Values `beta_i' in the input polygon data file (see above) define angles\n");
    printf("       of the image region corners:\n");
    printf("       beta_i = (1 - angle_i / pi) * 2; 0 by default,\n");
    printf("       where angle_i -- value of the interior angle at the i_th vertex. To get\n");
    printf("       a valid polygon, the sum of beta_i must be equal to 4. The image of the\n");
    printf("       very first vertex marked by \"*\" is placed in point (0,1), while the\n");
    printf("       edge connecting it with the next vertex is set directed towards\n");
    printf("       -i*infty. If not specified, beta_i is set to 0.\n");
    printf("    3. A uniform grid is specified by `nx' and `ny' parameters. This may be\n");
    printf("       overridden by specifying a custom grid from an input grid file. Such a\n");
    printf("       file should contain node coordinates scaled to 1x1 square.\n");
    printf("    4. The `nnodes' parameter affects precision and computation time; it is\n");
    printf("       advised to be set to the same value as -log10(<precision>) or slightly\n");
    printf("       higher.\n");
    printf("    5. If `newton 1' specified, the nonlinear system for sigmas is solved using\n");
    printf("       Gauss-Newton solver with Broyden update; otherwhile simple iterations\n");
    printf("       are used.\n");
    printf("    6. If no output file specified, the results are written to the standard\n");
    printf("       output.\n");
    printf("    7. If `sigmas' specified, the sigmas (an intermediate solution of nonlinear\n");
    printf("       system) will be read from/saved to a binary file. This may save a fair\n");
    printf("       bit of time in the case of multiple runs of the program.\n");
    printf("    8. If `rectangle' specified, the image polygon vertex coordinates will be\n");
    printf("       written to the corresponding text file.\n");
    printf("    9. Parameter `nppe' controls the coarseness of mapping of quadrilaterals.\n");
    printf("       Large values can take some time; small values can lead to failing\n");
    printf("       to map some points in difficult cases when quadrilateral images are\n");
    printf("       strongly distorted.\n");
    printf("  Acknowledgments. This program uses the following public code/algorithms:\n");
    printf("    1. CRDT algorithm by Tobin D. Driscoll and Stephen A. Vavasis -- for\n");
    printf("       conformal mapping.\n");
    printf("    2. `Triangle' by Jonathan Richard Shewchuk -- for Delaunay triangulation.\n");
    printf("    3. `SCPACK' by Lloyd N. Trefethen -- for Schwarz-Christoffel transform.\n");
    printf("    4. `DOPRI8' by E. Hairer, S.P. Norsett, G. Wanner -- for solving ODEs.\n");
    printf("    5. Shamos-Hoey algorithm implementation by softSurfer -- for testing a\n");
    printf("       polyline on self-intersections.\n");
    printf("\n");
    printf("  Good luck!\n");
    printf("\n");
    printf("  Pavel Sakov\n");
    printf("  April-May 2001 - January 2006\n");
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

static FILE* gg_fopen(char* fname, char* mode)
{
    FILE* f = fopen(fname, mode);

    if (f == NULL)
        quit("could not open \"%s\" in \"%s\" mode: %s\n", fname, mode, strerror(errno));

    return f;
}

/* This is modified skipToKeyEnd() from `sjwlib'.
 */
static int key_find(char* fname, FILE* fp, char* key)
{
    char buf[BUFSIZE];
    int len = strlen(key);
    int fpos;
    char* s;
    int rewound = 0;

    do {
        fpos = ftell(fp);
        s = fgets(buf, BUFSIZE, fp);
        if (s == NULL) {
            if (!rewound) {
                fpos = 0;
                if (fseek(fp, fpos, 0))
                    quit("%s: could not rewind: %s\n", fname, strerror(errno));
                s = fgets(buf, BUFSIZE, fp);
                rewound = 1;
            } else
                return 0;
        }

        if (s == NULL)
            break;

        while (isspace(s[0]))
            s++;

        if (isspace(s[len]))
            s[len] = 0;
    } while (strcmp(key, s) != 0);

    if (s == NULL)
        return 0;

    if (fseek(fp, fpos + (s - buf) + len, 0))
        quit("%s: no characters found after \"%s\"\n", fname, key);

    return 1;
}

/* This is modified readkeyprm_s() from `sjwlib'.
 */
static int prm_read(char* fname, FILE* fp, char* key, char* p)
{
    char line[BUFSIZE];
    char* s = line;
    char* r = p;

    if (!key_find(fname, fp, key))
        return 0;

    if (fgets(line, BUFSIZE, fp) != line)
        quit("%s: key \"%s\": read failed: %s\n", fname, key, strerror(errno));

    for (s = line; *s && isspace(*s); s++);

    for (r = p; *s && (*s != '\n'); *r++ = *s++);

    *r = 0;

    if (gg_verbose)
        fprintf(stderr, "-> %s = \"%s\"\n", key, p);

    return 1;
}

static gridgen* gridgen_init(void)
{
    gridgen* gg = malloc(sizeof(gridgen));

    gg->prmfname = NULL;
    gg->datafname = NULL;
    gg->sigmafname = NULL;
    gg->rectfname = NULL;
    gg->out = NULL;
#if !defined (GRIDGEN_STANDALONE) && defined(HAVE_GRIDNODES_H)
    gg->gn = NULL;
#endif
    gg->fsigma = NULL;
    gg->frect = NULL;
    gg->nrect = NULL;
    gg->xrect = NULL;
    gg->yrect = NULL;
    gg->xmin = DBL_MAX;
    gg->xmax = -DBL_MAX;
    gg->ymin = DBL_MAX;
    gg->ymax = -DBL_MAX;
    gg->thin = THIN_DEF;
    gg->checksimplepoly = CHECKSIMPLEPOLY_DEF;
    gg->reversed = 0;
    gg->nx = -1;
    gg->ny = -1;
    gg->ngridpoints = 0;
    gg->gridpoints = NULL;
    gg->eps = EPS_DEF;
    gg->nnodes = NNODES_DEF;
    gg->newton = NEWTON_DEF;
    gg->vertices = NULL;
    gg->d = NULL;
    gg->nquadrilaterals = 0;
    gg->lastqi = 0;
    gg->nz = 0;
    gg->zs = NULL;
    gg->betas = NULL;
    gg->sc = NULL;
    gg->nsigmas = NULL;
    gg->sigmas = NULL;
    gg->ws = NULL;
    gg->As = NULL;
    gg->Bs = NULL;
    gg->ncorners = -1;
    gg->rzs = NULL;
    gg->rvids = NULL;
    gg->newxmin = DBL_MAX;
    gg->newxmax = -DBL_MAX;
    gg->newymin = DBL_MAX;
    gg->newymax = -DBL_MAX;
    gg->newbetas = NULL;
    gg->newsc = NULL;
    gg->newAs = NULL;
    gg->newBs = NULL;
    gg->newzs = NULL;
    gg->newrzs = NULL;
    gg->nppe = NPPE_DEF;
    gg->nppq = gg->nppe * 4 + 5;
    gg->nqivertices = NULL;
    gg->qivertices = NULL;
    gg->lastqid = -1;
    gg->diagonal = -1;

    return gg;
}

static gridgen* gridgen_create(char* prmfname)
{
    gridgen* gg = gridgen_init();
    char buf[BUFSIZE];
    FILE* prm = NULL;
    FILE* data = NULL;
    int nold;

    prm = gg_fopen(prmfname, "r");
    gg->prmfname = prmfname;

    gg->vertices = vertlist_create();

    if (!prm_read(prmfname, prm, "input", buf))
        quit("%s: key \"input\" not found\n", prmfname);
    gg->datafname = strdup(buf);
    data = gg_fopen(buf, "r");

    if (gg_verbose)
        fprintf(stderr, "reading:\n");
    vertlist_read(gg->vertices, data);

    if (gg->vertices->n == 0)
        quit("%s: no XY points found\n", buf);
    if (gg->vertices->n < NMIN)
        quit("%s: %d points found -- too few (should be >= %d)\n", buf, gg->vertices->n, NMIN);
    if (gg->vertices->n > NMAX)
        quit("%s: %d points found -- too many (should be <= %d)\n", buf, gg->vertices->n, NMAX);
    if (gg_verbose)
        fprintf(stderr, "  %d vertices after read\n", gg->vertices->n);

    if (gg_verbose > 1)
        vertlist_print(gg->vertices, stderr);

    gg->reversed = 0;
    if (vertlist_area(gg->vertices) < 0.0) {
        if (gg_verbose)
            fprintf(stderr, "  clockwise, reversing\n");
        vertlist_change_order(&gg->vertices);
        gg->reversed = 1;
    } else if (gg_verbose)
        fprintf(stderr, "  counterclockwise ok\n");

    vertlist_find_minmax(gg->vertices, &gg->xmin, &gg->xmax, &gg->ymin, &gg->ymax);

    nold = gg->vertices->n;
    if (prm_read(prmfname, prm, "thin", buf))
        gg->thin = atoi(buf);
    if (gg->thin) {
        if (gg_verbose)
            fprintf(stderr, "thinning:\n");
        vertlist_thin(gg->vertices, (gg->xmax - gg->xmin) / BIGDOUBLE, (gg->ymax - gg->ymin) / BIGDOUBLE);

        if (gg->vertices->n < NMIN)
            quit("%s: %d vertices after thinning -- too little (should be >= %d)\n", buf, gg->vertices->n, NMIN);

        if (gg_verbose)
            fprintf(stderr, "  %d vertices after thinning\n", gg->vertices->n);
    }

    if (prm_read(prmfname, prm, "checksimplepoly", buf))
        gg->checksimplepoly = atoi(buf);
    if (gg->checksimplepoly) {
        double* x = NULL;
        double* y = NULL;

        if (gg_verbose)
            fprintf(stderr, "checking for self-intersections:\n");
        vertlist_toxy(gg->vertices, &x, &y);
        if (!issimplepolygon(gg->vertices->n, x, y)) {
            fprintf(stderr, "  Beware that code for testing a polyline on self-intersections is still in beta\n");
            fprintf(stderr, "  stage. You may try to skip the test by adding parameter \"checksimplepoly 0\"\n");
            fprintf(stderr, "  to your paremeter file. (If true, you should get a segfault during\n");
            fprintf(stderr, "  triangulation.)\n");
            quit("%s: not a simple polyline\n", gg->datafname);
        }
        free(x);
        free(y);
    }

    if (gg_verbose > 1 && gg->vertices->n != nold)
        vertlist_print(gg->vertices, stderr);

    if (prm_read(prmfname, prm, "nnodes", buf))
        gg->nnodes = atoi(buf);
    if (gg->nnodes < NNODES_MIN)
        gg->nnodes = NNODES_MIN;
    if (gg->nnodes > NNODES_MAX)
        gg->nnodes = NNODES_MAX;
    if (gg_verbose)
        fprintf(stderr, "nnodes = %d\n", gg->nnodes);

    if (prm_read(prmfname, prm, "newton", buf))
        gg->newton = atoi(buf);

    if (prm_read(prmfname, prm, "precision", buf))
        gg->eps = atof(buf);
    if (gg->eps < EPS_MIN)
        gg->eps = EPS_MIN;
    if (gg->eps > EPS_MAX)
        gg->eps = EPS_MAX;
    if (gg_verbose)
        fprintf(stderr, "precision = %3g\n", gg->eps);

    if (prm_read(prmfname, prm, "grid", buf)) {
        FILE* gridfile = gg_fopen(buf, "r");
        vertlist* l = vertlist_create();

        if (gridfile == NULL)
            quit("could not open \"%s\": %s\n", buf, strerror(errno));

        vertlist_read(l, gridfile);
        gg->gridpoints = vertlist_tozdouble(l);
        gg->ngridpoints = l->n;

        fclose(gridfile);
        vertlist_destroy(l);
    }

    if (gg->gridpoints == NULL) {
        if (prm_read(prmfname, prm, "nx", buf))
            gg->nx = atoi(buf);
        else
            gg->nx = N_DEF;
        if (prm_read(prmfname, prm, "ny", buf))
            gg->ny = atoi(buf);
        else
            gg->ny = N_DEF;
        if (gg->nx < N_MIN)
            gg->nx = N_MIN;
        if (gg->nx > N_MAX)
            gg->nx = N_MAX;
        if (gg->ny < N_MIN)
            gg->ny = N_MIN;
        if (gg->ny > N_MAX)
            gg->ny = N_MAX;

        if (gg_verbose)
            fprintf(stderr, "going to generate %dx%d grid\n", gg->nx, gg->ny);
    }

    if (prm_read(prmfname, prm, "output", buf))
        gg->out = gg_fopen(buf, "w");
    else
        gg->out = stdout;

    if (prm_read(prmfname, prm, "sigmas", buf)) {
        gg->sigmafname = strdup(buf);
        gg->fsigma = fopen(buf, "r");
    }

    if (prm_read(prmfname, prm, "rectangle", buf)) {
        gg->rectfname = strdup(buf);
        gg->frect = gg_fopen(buf, "w");
    }

    if (prm_read(prmfname, prm, "ppe", buf))
        gg->nppe = atoi(buf);
    gg->nppq = gg->nppe * 4 + 5;

    fclose(data);
    fclose(prm);

    return gg;
}

#if !defined (GRIDGEN_STANDALONE) && defined(HAVE_GRIDNODES_H)

static gridgen* gridgen_create2(int nbdry, double xbdry[], double ybdry[], double beta[], int ul, int nx, int ny, int ngrid, double xgrid[], double ygrid[], int nnodes, int newton, double precision, int checksimplepoly, int thin, int nppe, int verbose, int* nsigmas, double** sigmas, int* nrect, double** xrect, double** yrect)
{
    gridgen* gg = gridgen_init();
    int nold;

    gridgen_setverbose(verbose);

    if (gg_verbose)
        fprintf(stderr, "reading:\n");
    gg->vertices = vertlist_create2(nbdry, xbdry, ybdry, beta);
    vertlist_setfirstnode(gg->vertices, ul);

    if (gg->vertices->n < NMIN)
        quit("%d points found -- too few (should be >= %d)\n", gg->vertices->n, NMIN);
    if (gg->vertices->n > NMAX)
        quit("%d points found -- too many (should be <= %d)\n", gg->vertices->n, NMAX);
    if (gg_verbose)
        fprintf(stderr, "  %d vertices after read\n", gg->vertices->n);

    if (gg_verbose > 1)
        vertlist_print(gg->vertices, stderr);

    gg->reversed = 0;
    if (vertlist_area(gg->vertices) < 0.0) {
        if (gg_verbose)
            fprintf(stderr, "  clockwise, reversing\n");
        vertlist_change_order(&gg->vertices);
        gg->reversed = 1;
    } else if (gg_verbose)
        fprintf(stderr, "  counterclockwise ok\n");

    vertlist_find_minmax(gg->vertices, &gg->xmin, &gg->xmax, &gg->ymin, &gg->ymax);

    nold = gg->vertices->n;
    gg->thin = thin;
    if (gg->thin) {
        if (gg_verbose)
            fprintf(stderr, "thinning:\n");
        vertlist_thin(gg->vertices, (gg->xmax - gg->xmin) / BIGDOUBLE, (gg->ymax - gg->ymin) / BIGDOUBLE);

        if (gg->vertices->n < NMIN)
            quit("%d vertices after thinning -- too little (should be >= %d)\n", gg->vertices->n, NMIN);

        if (gg_verbose)
            fprintf(stderr, "  %d vertices after thinning\n", gg->vertices->n);
    }

    gg->checksimplepoly = checksimplepoly;
    if (gg->checksimplepoly) {
        double* x = NULL;
        double* y = NULL;

        if (gg_verbose)
            fprintf(stderr, "checking for self-intersections:\n");
        vertlist_toxy(gg->vertices, &x, &y);
        if (!issimplepolygon(gg->vertices->n, x, y)) {
            fprintf(stderr, "  Beware that code for testing a polyline on self-intersections is still in beta\n");
            fprintf(stderr, "  stage. You may try to skip the test by adding parameter \"checksimplepoly 0\"\n");
            fprintf(stderr, "  to your paremeter file. (If true, you should get a segfault during\n");
            fprintf(stderr, "  triangulation.)\n");
            quit("%s: not a simple polyline\n", gg->datafname);
        }
        free(x);
        free(y);
    }

    if (gg_verbose > 1 && gg->vertices->n != nold)
        vertlist_print(gg->vertices, stderr);

    gg->nnodes = nnodes;
    if (gg->nnodes < NNODES_MIN)
        gg->nnodes = NNODES_MIN;
    if (gg->nnodes > NNODES_MAX)
        gg->nnodes = NNODES_MAX;
    if (gg_verbose)
        fprintf(stderr, "nnodes = %d\n", gg->nnodes);

    gg->newton = newton;

    gg->eps = precision;
    if (gg_verbose)
        fprintf(stderr, "precision = %3g\n", gg->eps);

    if (ngrid > 0) {
        int i;

        gg->ngridpoints = ngrid;
        gg->gridpoints = malloc(ngrid * sizeof(zdouble));

        for (i = 0; i < ngrid; ++i)
            gg->gridpoints[i] = xgrid[i] + I * ygrid[i];
    }

    if (gg->gridpoints == NULL) {
        gg->nx = nx;
        gg->ny = ny;
        if (gg->nx < N_MIN)
            gg->nx = N_MIN;
        if (gg->nx > N_MAX)
            gg->nx = N_MAX;
        if (gg->ny < N_MIN)
            gg->ny = N_MIN;
        if (gg->ny > N_MAX)
            gg->ny = N_MAX;

        if (gg_verbose)
            fprintf(stderr, "going to generate %dx%d grid\n", gg->nx, gg->ny);
    }

    if (gg->gridpoints == NULL)
        gg->gn = gridnodes_create(gg->nx, gg->ny, NT_NONE);
    else
        gg->gn = gridnodes_create(gg->ngridpoints, 1, NT_NONE);

    /*
     * Storage coordinates of the working stuff
     */
    gg->nsigmas = nsigmas;
    gg->sigmas = sigmas;
    gg->nrect = nrect;
    gg->xrect = xrect;
    gg->yrect = yrect;

    return gg;
}

#endif

static void get_rpoints(gridgen* gg)
{
    vertlist* l = gg->vertices;
    vertnode* now = l->first;
    int n = l->n;
    double sum, mult;
    int i, count;

    if (gg_verbose) {
        fprintf(stderr, "getting marked vertices:\n");
        fflush(stderr);
    }

    /*
     * count corner markers 
     */
    for (i = 0, count = 0, sum = 0.0; i < n; ++i) {
        point* p = &now->p;

        if (p->z != 0) {
            if (fabs(p->z) > 2.0)
                quit("vertex %d (%.15g,%.15g): |beta| = %.15g > 2\n", i, p->x, p->y, fabs(p->z));
            count++;
            sum += p->z;
        }
        now = now->next;
    }

    if (gg_verbose)
        fprintf(stderr, "  image region: %d corners\n", count);

    if (count < 3)
        quit("less than 3 corners in the image region\n", count);
    if (fabs(fabs(sum) - 4.0) > EPS)
        quit("sum of beta_i = %.15f != 4\n", sum);

    /*
     * normalize sum of image region corner angles if a small discrepancy
     * exists 
     */
    mult = sum / 4.0;
    if (mult != 1.0)
        for (i = 0; i < n; ++i) {
            now->p.z *= mult;
            now = now->next;
        }

    gg->ncorners = count;
    gg->rzs = malloc(sizeof(zdouble) * count);
    gg->rvids = malloc(sizeof(int) * count);
    gg->newrzs = malloc(sizeof(zdouble) * count);

    /*
     * some extra vertices may be inserted, so we will get `rvids' later 
     */
    for (i = 0; i < gg->ncorners; ++i)
        gg->rvids[i] = -1;

    for (i = 0, count = 0; i < n; ++i) {
        point* p = &now->p;

        if (p->z != 0.0) {
            gg->rzs[count] = p->x + I * p->y;
            count++;
        }
        now = now->next;
    }

    if (gg_verbose > 1) {
        fprintf(stderr, "image region corners:\n");
        for (i = 0; i < gg->ncorners; ++i)
            fprintf(stderr, "  %d: (%.10g,%.10g)\n", i, creal(gg->rzs[i]), cimag(gg->rzs[i]));
    }
}

static double* calculate_betas(vertlist* l)
{
    int n = l->n;
    double* betas = malloc(n * sizeof(double));
    vertnode* now = l->first;
    int i;

    point* pprev = &now->prev->p;
    point* pnow = &now->p;
    point* pnext = &now->next->p;

    double dx1 = pnow->x - pprev->x;
    double dy1 = pnow->y - pprev->y;
    double dx2 = pnext->x - pnow->x;
    double dy2 = pnext->y - pnow->y;

    if (gg_verbose) {
        fprintf(stderr, "calculating betas:\n");
        fflush(stderr);
    }

    for (i = 0; i < n; ++i) {

        double arg = carg(-(dx1 + I * dy1) / (dx2 + I * dy2));

        if (arg < 0)
            arg += M_2PI;
        betas[i] = arg / M_PI - 1.0;

        if (gg_verbose > 1)
            fprintf(stderr, "  %d: %f\n", i, betas[i]);

        now = now->next;

        pprev = pnow;
        pnow = pnext;
        pnext = &now->next->p;
        dx1 = dx2;
        dy1 = dy2;
        dx2 = pnext->x - pnow->x;
        dy2 = pnext->y - pnow->y;
    }

    return betas;
}

static quadrilateral* set_quadrilaterals(delaunay* d)
{
    quadrilateral* out = malloc((d->npoints - 3) * sizeof(quadrilateral));
    istack* diags = istack_create();
    istack* neighbours = istack_create();
    int* edges = d->edges;
    int maxdiff = d->npoints - 1;
    int n = 0;
    int i, ii;

    if (gg_verbose) {
        fprintf(stderr, "getting quadrilaterals:\n");
        fflush(stderr);
    }

    for (i = 0, ii = 0; i < d->nedges; ++i, ii += 2) {
        int diff = abs(edges[ii] - edges[ii + 1]);

        if (diff != 1 && diff != maxdiff) {
            istack_push(diags, edges[ii]);
            istack_push(diags, edges[ii + 1]);
            n++;
        }
    }

    if (gg_verbose)
        fprintf(stderr, "  %d diagonals\n", n);

    for (i = 0; i < n; ++i) {
        quadrilateral* q = &out[i];
        int vid1 = istack_pop(diags);
        int vid2 = istack_pop(diags);
        int i1, i2;
        triangle* triangles[2];
        int nt = 0;

        q->id = i;
        for (i1 = 0; i1 < d->n_point_triangles[vid1]; ++i1) {
            int tid1 = d->point_triangles[vid1][i1];

            for (i2 = 0; i2 < d->n_point_triangles[vid2]; ++i2) {
                int tid2 = d->point_triangles[vid2][i2];

                if (tid1 == tid2) {
                    triangles[nt] = &d->triangles[tid1];
                    q->tids[nt] = tid1;
                    nt++;
                }
            }
        }
        q->vids[0] = vid1;
        q->vids[2] = vid2;
        for (ii = 0; ii < 3; ++ii) {
            int vid = triangles[0]->vids[ii];

            if (vid != vid1 && vid != vid2)
                q->vids[1] = vid;
        }
        for (ii = 0; ii < 3; ++ii) {
            int vid = triangles[1]->vids[ii];

            if (vid != vid1 && vid != vid2)
                q->vids[3] = vid;
        }

        for (ii = 0; triangles[0]->vids[ii] != vid1; ++ii);
        if (triangles[0]->vids[(ii + 1) % 3] == vid2) {
            int tmp = q->vids[1];

            q->vids[1] = q->vids[3];
            q->vids[3] = tmp;
            tmp = q->tids[0];
            q->tids[0] = q->tids[1];
            q->tids[1] = tmp;
        }
    }

    /*
     * set neighbours 
     */
    n = d->npoints - 3;
    for (i = 0; i < n; ++i) {
        quadrilateral* q = &out[i];
        int t0 = q->tids[0];
        int t1 = q->tids[1];
        int ii;

        istack_reset(neighbours);

        for (ii = 0; ii < n; ++ii) {
            quadrilateral* qq = &out[ii];
            int tt0 = qq->tids[0];
            int tt1 = qq->tids[1];

            if (ii == i)
                continue;

            if (t0 == tt0 || t0 == tt1 || t1 == tt0 || t1 == tt1)
                istack_push(neighbours, ii);
        }

        q->nneighbours = neighbours->n;
        q->neighbours = malloc(q->nneighbours * sizeof(int));
        memcpy(q->neighbours, neighbours->v, q->nneighbours * sizeof(int));
    }

    if (gg_verbose > 1)
        for (i = 0; i < n; ++i) {
            quadrilateral* q = &out[i];
            int ii;

            fprintf(stderr, "  %d:  (%d,%d,%d,%d) (%d,%d)", i, q->vids[0], q->vids[1], q->vids[2], q->vids[3], q->tids[0], q->tids[1]);
            for (ii = 0; ii < q->nneighbours; ++ii)
                fprintf(stderr, (ii == 0) ? " (%d" : ",%d", q->neighbours[ii]);
            fprintf(stderr, ")\n");
        }

    istack_destroy(neighbours);
    istack_destroy(diags);

    return out;
}

static double* calculate_cis(int nquads, quadrilateral* quads, delaunay* d)
{
    double* cis = malloc(nquads * sizeof(double));
    point* points = d->points;
    int i;

    if (gg_verbose) {
        fprintf(stderr, "calculating log|ro|:\n");
        fflush(stderr);
    }

    for (i = 0; i < nquads; ++i) {
        quadrilateral* q = &quads[i];
        point* p1 = &points[q->vids[0]];
        point* p2 = &points[q->vids[1]];
        point* p3 = &points[q->vids[2]];
        point* p4 = &points[q->vids[3]];
        zdouble z1 = p1->x + I * p1->y;
        zdouble z2 = p2->x + I * p2->y;
        zdouble z3 = p3->x + I * p3->y;
        zdouble z4 = p4->x + I * p4->y;

        cis[i] = log(cabs((z4 - z1) * (z2 - z3) / (z3 - z4) / (z1 - z2)));
        if (gg_verbose > 1)
            fprintf(stderr, "  %d: %.10g\n", i, cis[i]);
    }

    return cis;
}

static void process_quadri(quadrilateral qs[], int qindex, int do_second_triangle, double x[], zdouble w[], int done[], int* count)
{
    quadrilateral* q = &qs[qindex];
    double ro = -exp(x[qindex]);
    int i;

    if (do_second_triangle) {
        zdouble* pa = &w[q->vids[0]];
        zdouble* pb = &w[q->vids[1]];
        zdouble* pc = &w[q->vids[2]];

        if (*pc != *pb) {
            zdouble h = ro * (*pb - *pa) / (*pc - *pb);

            w[q->vids[3]] = (h * *pc + *pa) / (h + 1.0);
        } else
            /*
             * OK, crowding occurs here. The principal point with CRDT is
             * that it occurs for quadrilaterals far from the starting one
             * for this (local) mapping and can be safely ignored. 
             */
            w[q->vids[3]] = *pc;
    } else {
        zdouble* pa = &w[q->vids[2]];
        zdouble* pb = &w[q->vids[3]];
        zdouble* pc = &w[q->vids[0]];

        if (*pc != *pb) {
            zdouble h = ro * (*pb - *pa) / (*pc - *pb);

            w[q->vids[1]] = (h * *pc + *pa) / (h + 1.0);
        } else
            /*
             * see the comment above 
             */
            w[q->vids[1]] = *pc;
    }
    done[qindex] = 1;
    (*count)++;

    for (i = 0; i < q->nneighbours; ++i) {
        int nid = q->neighbours[i];
        quadrilateral* q0 = &qs[nid];

        if (done[nid])
            continue;
        if (q0->tids[0] == q->tids[0] || q0->tids[0] == q->tids[1])
            process_quadri(qs, nid, 1, x, w, done, count);
        else
            process_quadri(qs, nid, 0, x, w, done, count);
    }
}

static void calculate_fi(int nq, quadrilateral qs[], int qindex, double x[], zdouble w[])
{
    int count = 0;
    int* done = calloc(nq, sizeof(int));

    /*
     * startup 
     */
    double ro = exp(x[qindex]);
    double alpha = atan(sqrt(ro));
    double sin_alpha = sin(alpha);
    double cos_alpha = cos(alpha);
    quadrilateral* q = &qs[qindex];

    w[q->vids[0]] = cos_alpha + I * sin_alpha;
    w[q->vids[1]] = -cos_alpha + I * sin_alpha;
    w[q->vids[2]] = -cos_alpha - I * sin_alpha;

    /*
     * go 
     */
    process_quadri(qs, qindex, 1, x, w, done, &count);  /* recursive */

    free(done);
}

static zdouble integrate(swcr* sc, int vid, zdouble w[])
{
    return sc_w2z(sc, w[vid], vid, ZZERO, ZZERO, -1, 1.0 + 0.0 * I, w);
}

static void calculate_dzetas(int nq, quadrilateral qs[], int nz, zdouble zs[], int qindex, swcr* sc, double x[], zdouble w[], zdouble* A, zdouble* B, zdouble dzetas[])
{
    quadrilateral* q = &qs[qindex];
    int* vids = q->vids;
    zdouble z0 = zs[vids[0]];
    zdouble z2 = zs[vids[2]];
    int i;

    calculate_fi(nq, qs, qindex, x, w);

    for (i = 0; i < 4; ++i) {
        int vid = vids[i];

        dzetas[i] = integrate(sc, vid, w);
    }

    *B = (z2 - z0) / (dzetas[2] - dzetas[0]);
    *A = z0 - *B * dzetas[0];
}

static void F(double x[], double f[], void* p)
{
    gridgen* gg = (gridgen*) p;
    int nq = gg->nquadrilaterals;
    int nz = gg->nz;
    zdouble dzetas[4];
    int i;

    for (i = 0; i < nq; ++i) {
        zdouble dzeta0, dzeta1, dzeta2, dzeta3;

        calculate_dzetas(nq, gg->quadrilaterals, nz, gg->zs, i, gg->sc, x, &gg->ws[i * nz], &gg->As[i], &gg->Bs[i], dzetas);

        dzeta0 = dzetas[0];
        dzeta1 = dzetas[1];
        dzeta2 = dzetas[2];
        dzeta3 = dzetas[3];
        f[i] = log(cabs((dzeta3 - dzeta0) * (dzeta1 - dzeta2) / (dzeta2 - dzeta3) / (dzeta0 - dzeta1))) - gg->cis[i];
        if (isnan(f[i])) {
            if (gg_verbose && sc_issingular(gg->sc)) {
                fprintf(stderr, "  You have encountered a (very rare) case when a pre-vertice in an arising\n");
                fprintf(stderr, "  Schwarz-Christoffel integral coincides with the point of crowding. You may\n");
                fprintf(stderr, "  try to overcome this by either:\n");
                fprintf(stderr, "    -- deleting the file with the initial sigma values (the one pointed by\n");
                fprintf(stderr, "       \"sigmas\" entry in the parameter file);\n");
                fprintf(stderr, "    -- changing \"newton 1\" to \"newton 0\" or vice versa in the parameter file;\n");
                fprintf(stderr, "    -- slightly moving the boundary polygon vertices.\n");
            }
            quit("F(): NaN detected\n");
        }
    }
}

static void find_sigmas(gridgen* gg, func F)
{
    int n = gg->nquadrilaterals;
    double* x = malloc(n * sizeof(double));
    double* f = malloc(n * sizeof(double));
    double* w = NULL;
    double error = DBL_MAX;
    double error_prev;
    int count = 0;
    int nlastbad = 0;

    if (gg_verbose) {
        fprintf(stderr, "solving for sigmas:\n");
        fflush(stderr);
    }

    if (gg->sigmas != NULL && *gg->nsigmas == n && gg->sigmas != NULL && *gg->sigmas != NULL)
        memcpy(x, *gg->sigmas, n * sizeof(double));
    else if (gg->fsigma == NULL || (int) fread(x, sizeof(double), n, gg->fsigma) != n)
        memcpy(x, gg->cis, n * sizeof(double));

    gg->ws = calloc(n * gg->vertices->n, sizeof(zdouble));
    gg->As = calloc(gg->vertices->n, sizeof(zdouble));
    gg->Bs = calloc(gg->vertices->n, sizeof(zdouble));

    if (gg->newton == 1) {
        double* ww;
        int i;

        w = calloc(n * n, sizeof(double));

        for (i = 0, ww = w; i < n; ++i, ww += n + 1)
            ww[0] = -1.0;
    } else if (gg->newton > 1)
        w = calloc(n * M * 2, sizeof(double));

    if (gg_verbose == 1)
        fprintf(stderr, "  ");

    do {
        int i;

        error_prev = error;
        error = 0.0;

        if (gg->newton) {
            /*
             * newton method with broyden update
             */
            if (count == 0)
                F(x, f, gg);
            else
                broyden_update(F, n, x, f, w, gg);
        } else {
            /*
             * simple iterations 
             */
            F(x, f, gg);
            for (i = 0; i < n; ++i)
                x[i] -= f[i];
        }

        for (i = 0; i < n; ++i)
            error += f[i] * f[i];
        error = sqrt(error);

        if (error < error_prev)
            nlastbad = 0;
        else
            nlastbad++;

        if (gg_verbose < 2 && gg_verbose + nlastbad > 2) {
            fprintf(stderr, "\n  error[%d] = %.3g, error[%d] = %.3g", count, error_prev, count + 1, error);
            fflush(stderr);
        } else if (gg_verbose > 1) {
            fprintf(stderr, "  %d: error = %.3g\n", count, error);
            fflush(stderr);
        } else if (gg_verbose) {
            fprintf(stderr, ".");
            fflush(stderr);
        }

        count++;
    } while (error > gg->eps && (!gg->newton || error_prev > gg->eps));

    if (gg_verbose)
        fprintf(stderr, "\n");
    if (gg_verbose > 1) {
        int i;

        fprintf(stderr, "  sigmas:\n");
        for (i = 0; i < n; ++i)
            fprintf(stderr, "  %d: %.15g\n", i, x[i]);
        fflush(stderr);
    }

    if (count > 1 && gg->sigmafname != NULL) {
        if (gg->fsigma != NULL)
            fclose(gg->fsigma);
        gg->fsigma = gg_fopen(gg->sigmafname, "w");
        if (gg_verbose)
            fprintf(stderr, "  saving sigmas to \"%s\":\n", gg->sigmafname);
        if ((int) fwrite(x, sizeof(double), n, gg->fsigma) != n)
            quit("could not save sigmas to \"%s\": %s\n", gg->sigmafname, strerror(errno));
        fclose(gg->fsigma);
        gg->fsigma = NULL;
    }

    if (w != NULL)
        free(w);
    free(f);

    if (gg->sigmas != NULL) {
        if (*gg->sigmas != NULL)
            free(*gg->sigmas);
        *gg->nsigmas = n;
        *gg->sigmas = x;
    } else
        free(x);
}

static void quadrilaterals_destroy(int n, quadrilateral* quads)
{
    int i;

    if (quads == NULL)
        return;
    for (i = 0; i < n; ++i)
        free(quads[i].neighbours);
    free(quads);
}

static void gridgen_destroy(gridgen* gg)
{
    if (gg == NULL)
        return;
    if (gg->zs != NULL)
        free(gg->zs);
    if (gg->ws != NULL)
        free(gg->ws);
    if (gg->As != NULL)
        free(gg->As);
    if (gg->Bs != NULL)
        free(gg->Bs);
    sc_destroy(gg->sc);
    if (gg->newbetas != NULL)
        free(gg->newbetas);
    if (gg->betas != NULL)
        free(gg->betas);
    quadrilaterals_destroy(gg->nquadrilaterals, gg->quadrilaterals);
    if (gg->cis != NULL)
        free(gg->cis);
    delaunay_destroy(gg->d);
    vertlist_destroy(gg->vertices);
    if (gg->datafname != NULL)
        free(gg->datafname);
    sc_destroy(gg->newsc);
    if (gg->newzs != NULL)
        free(gg->newzs);
    if (gg->newAs != NULL)
        free(gg->newAs);
    if (gg->Bs != NULL)
        free(gg->newBs);
    if (gg->gridpoints != NULL)
        free(gg->gridpoints);
    if (gg->out != stdout && gg->out != NULL)
        fclose(gg->out);
    if (gg->sigmafname != NULL)
        free(gg->sigmafname);
    if (gg->fsigma != NULL)
        fclose(gg->fsigma);
    if (gg->rectfname != NULL)
        free(gg->rectfname);
    if (gg->frect != NULL)
        fclose(gg->frect);
    if (gg->rzs != NULL)
        free(gg->rzs);
    if (gg->rvids != NULL)
        free(gg->rvids);
    if (gg->newrzs != NULL)
        free(gg->newrzs);
    if (gg->nqivertices != NULL) {
        free(gg->nqivertices);
        free(gg->qivertices);
    }
    free(gg);
}

static void gridgen_check(gridgen* gg)
{
    int n = gg->vertices->n;
    double* betas = gg->betas;
    double sum = 0.0;
    int i;

    if (gg_verbose) {
        fprintf(stderr, "checking input:\n");
        fflush(stderr);
    }

    for (i = 0; i < n; ++i)
        sum += betas[i];

    if (fabs(sum + 2.0) > EPS)
        quit("gridgen_check(): sum of betas = %.15f (should be 2)\n", sum);

    if (gg_verbose > 1)
        fprintf(stderr, "  sum of betas = 2.0 + %.3g\n", sum + 2.0);
}

static void get_rvids(gridgen* gg)
{
    vertlist* l = gg->vertices;
    vertnode* now = l->first;
    int n = l->n;
    int count, i;

    if (gg_verbose) {
        fprintf(stderr, "getting image region corner vertex indices:\n");
        fflush(stderr);
    }

    for (i = 0, count = 0; i < n; ++i) {
        point* p = &now->p;

        if (p->z != 0.0 && !isnan(p->z)) {
            gg->rvids[count] = i;
            count++;
        }
        now = now->next;
    }

    if (gg_verbose > 1) {
        zdouble* rzs = gg->rzs;

        for (i = 0; i < gg->ncorners; ++i)
            fprintf(stderr, "  %d: (%.10g,%.10g) %d\n", i, creal(rzs[i]), cimag(rzs[i]), gg->rvids[i]);
        fflush(stderr);
    }
}

static void set_newbetas(gridgen* gg)
{
    point* points = gg->d->points;
    int n = gg->ncorners;
    int i;

    if (gg_verbose) {
        fprintf(stderr, "setting new betas:\n");
        fflush(stderr);
    }

    gg->newbetas = calloc(gg->d->npoints, sizeof(double));

    for (i = 0; i < n; ++i) {
        int index = gg->rvids[i];

        if (points[index].z != 0)
            gg->newbetas[index] = -points[index].z * 0.5;
    }
}

static int get_first_quadrilateral(gridgen* gg)
{
    int qid = -1;
    istack* qids = istack_create();
    int n = gg->nquadrilaterals;
    int vid = gg->rvids[0];
    int vid1 = (vid + 1) % gg->d->npoints;
    int i;

    for (i = 0; i < n; ++i) {
        int* vids = gg->quadrilaterals[i].vids;

        if (vid == vids[0] || vid == vids[1] || vid == vids[2] || vid == vids[3])
            istack_push(qids, i);
    }

    assert(qids->n > 0);

    if (qids->n == 1)
        qid = qids->v[0];
    else
        for (i = 0; i < qids->n; ++i) {
            int* vids = gg->quadrilaterals[qids->v[i]].vids;

            if (vid1 == vids[0] || vid1 == vids[1] || vid1 == vids[2] || vid1 == vids[3]) {
                qid = qids->v[i];
                break;
            }
        }

    assert(qid >= 0);

    istack_destroy(qids);

    return qid;
}

/* Finds index (0..3) of a vertex within quadrilateral "newqid" not 
 * contained in the quadrilateral "qid".
 */
static int find_newvertice(int nq, quadrilateral qs[], int qid, int pqid)
{
    int* pvids = qs[pqid].vids;
    int* vids = qs[qid].vids;
    int i, j;

    for (i = 0; i < 4; ++i) {
        int vid = vids[i];

        for (j = 0; j < 4; ++j)
            if (vid == pvids[j])
                break;
        if (j == 4)
            return vid;
    }
    return -1;
}

static void calculate_zs(int nq, quadrilateral qs[], int nz, zdouble zs[], int qindex, int vid, swcr* sc, zdouble w[], int* count, int done[], zdouble A[], zdouble B[])
{
    quadrilateral* q = &qs[qindex];
    int* vids = q->vids;
    zdouble z0 = zs[vids[0]];
    zdouble z2 = zs[vids[2]];
    zdouble dzeta0, dzeta2, dzeta;
    int i;

    if (*count) {
        dzeta0 = integrate(sc, vids[0], &w[qindex * nz]);
        dzeta2 = integrate(sc, vids[2], &w[qindex * nz]);
        dzeta = integrate(sc, vid, &w[qindex * nz]);

        B[qindex] = (z2 - z0) / (dzeta2 - dzeta0);
        A[qindex] = z0 - B[qindex] * dzeta0;

        zs[vid] = A[qindex] + B[qindex] * dzeta;
    }

    (*count)++;
    done[qindex] = 1;

    for (i = 0; i < q->nneighbours; ++i) {
        int nid = q->neighbours[i];

        if (done[nid])
            continue;
        {
            int vid = find_newvertice(nq, qs, nid, qindex);

            calculate_zs(nq, qs, nz, zs, nid, vid, sc, w, count, done, A, B);
        }
    }
}

static void find_zminmax(int n, zdouble zs[], double* xmin, double* xmax, double* ymin, double* ymax)
{
    int i;

    *xmin = DBL_MAX;
    *xmax = -DBL_MAX;
    *ymin = DBL_MAX;
    *ymax = -DBL_MAX;

    for (i = 0; i < n; ++i) {
        zdouble z = zs[i];
        double x = creal(z);
        double y = cimag(z);

        if (x < *xmin)
            *xmin = x;
        if (x > *xmax)
            *xmax = x;
        if (y < *ymin)
            *ymin = y;
        if (y > *ymax)
            *ymax = y;
    }
}

static void align_newzs(int n, zdouble rzs[], int rvids[], int nz, zdouble zs[], double betas[], double eps)
{
    double eps5 = eps * 5.0;
    double beta;
    int i;

    for (i = 0; i < n; ++i) {
        if (fabs(creal(rzs[i])) < eps5)
            rzs[i] = I * cimag(rzs[i]);
        else if (fabs(1.0 - creal(rzs[i])) < eps5)
            rzs[i] = 1.0 + I * cimag(rzs[i]);
        if (fabs(cimag(rzs[i])) < eps5)
            rzs[i] = creal(rzs[i]);
        else if (fabs(1.0 - cimag(rzs[i])) < eps5)
            rzs[i] = creal(rzs[i]) + I * 1.0;
    }

    for (i = 0; i < n; ++i)
        zs[rvids[i]] = rzs[i];

    for (i = 0, beta = 0.0; i < n; ++i) {
        int vid0 = rvids[i];
        int vid1 = rvids[(i + 1) % n];
        zdouble z = zs[vid0];
        int j, jj;

        if (fabs(fmod(beta, 0.5)) < eps) {
            int k = rint(fabs(beta / 0.5));

            if (k % 2 == 0) {
                double x = creal(z);

                for (j = 0, jj = (vid0 + 1) % nz; j < (vid1 + nz - vid0) % nz; ++j, jj = (jj + 1) % nz)
                    zs[jj] = x + I * cimag(zs[jj]);
            } else {
                double y = cimag(z);

                for (j = 0, jj = (vid0 + 1) % nz; j < (vid1 + nz - vid0) % nz; ++j, jj = (jj + 1) % nz)
                    zs[jj] = creal(zs[jj]) + I * y;
            }
        }
        beta += betas[vid1];
    }

    for (i = 0; i < n; ++i)
        rzs[i] = zs[rvids[i]];
}

static void write_zarray(int n, zdouble* z, FILE* f)
{
    int i;

    for (i = 0; i < n; ++i)
        fprintf(f, "%f %f\n", creal(z[i]), cimag(z[i]));
    fflush(f);
}

static void calculate_newzs(gridgen* gg)
{
    int nz = gg->nz;
    int nq = gg->nquadrilaterals;
    int* done = calloc(nq, sizeof(int));
    int count = 0;
    int vid0 = gg->rvids[0];
    int qid0;
    quadrilateral* q0 = NULL;
    int i;

    if (gg_verbose) {
        fprintf(stderr, "calculating image region:\n");
        fflush(stderr);
    }

    gg->newzs = calloc(nz, sizeof(zdouble));
    gg->newAs = malloc(nq * sizeof(zdouble));
    gg->newBs = malloc(nq * sizeof(zdouble));

    qid0 = get_first_quadrilateral(gg);
    q0 = &gg->quadrilaterals[qid0];
    if (gg_verbose > 1) {
        int* vids = q0->vids;

        fprintf(stderr, "  start from quadrilateral %d (%d,%d,%d,%d)\n", qid0, vids[0], vids[1], vids[2], vids[3]);
    }

    /*
     * startup 
     */
    {
        zdouble* A = gg->newAs;
        zdouble* B = gg->newBs;
        int* vids = q0->vids;
        int vid1 = (vid0 + 1) % nz;
        int qvid0 = -1;
        int qvid1 = -1;
        zdouble dzetas[4];

        for (i = 0; i < 4; ++i) {
            int vid = vids[i];

            dzetas[i] = integrate(gg->newsc, vid, &gg->ws[qid0 * nz]);
        }
        for (i = 0; i < 4; ++i) {
            if (vids[i] == vid0)
                qvid0 = i;
            if (vids[i] == vid1)
                qvid1 = i;
        }
        B[qid0] = -I * 1.0 / (dzetas[qvid1] - dzetas[qvid0]);
        A[qid0] = 1.0 * I - B[qid0] * dzetas[qvid0];

        for (i = 0; i < 4; ++i)
            gg->newzs[vids[i]] = A[qid0] + B[qid0] * dzetas[i];
    }

    calculate_zs(nq, gg->quadrilaterals, nz, gg->newzs, qid0, -1, gg->newsc, gg->ws, &count, done, gg->newAs, gg->newBs);       /* recursive 
                                                                                                                                 */

    for (i = 0; i < gg->ncorners; ++i) {
        int vid = gg->rvids[i];

        gg->newrzs[i] = gg->newzs[vid];
    }

    if (gg_verbose > 1) {
        fprintf(stderr, "image vertices:\n");
        for (i = 0; i < nz; ++i)
            fprintf(stderr, "  %d: (%.15g, %.15g)\n", i, creal(gg->newzs[i]), cimag(gg->newzs[i]));
        fprintf(stderr, "image corners:\n");
        for (i = 0; i < gg->ncorners; ++i)
            fprintf(stderr, "  %d: (%.15g, %.15g)\n", i, creal(gg->newrzs[i]), cimag(gg->newrzs[i]));
        fprintf(stderr, "aligning:\n");
    }

    align_newzs(gg->ncorners, gg->newrzs, gg->rvids, gg->nz, gg->newzs, gg->newbetas, gg->eps);

    if (gg_verbose > 1) {
        fprintf(stderr, "image vertices:\n");
        for (i = 0; i < nz; ++i)
            fprintf(stderr, "  %d: (%.15g, %.15g)\n", i, creal(gg->newzs[i]), cimag(gg->newzs[i]));
        fprintf(stderr, "image corners:\n");
        for (i = 0; i < gg->ncorners; ++i)
            fprintf(stderr, "  %d: (%.15g, %.15g)\n", i, creal(gg->newrzs[i]), cimag(gg->newrzs[i]));
    }

    find_zminmax(gg->ncorners, gg->newrzs, &gg->newxmin, &gg->newxmax, &gg->newymin, &gg->newymax);
    if (gg_verbose)
        fprintf(stderr, "  conformal modulus = %.5g\n", (gg->newxmax - gg->newxmin) / (gg->newymax - gg->newymin));

    if (gg->frect != NULL) {
        if (gg_verbose)
            fprintf(stderr, "saving image region:\n");
        write_zarray(gg->ncorners, gg->newrzs, gg->frect);
    }
    if (gg->nrect != NULL) {
        *gg->nrect = gg->ncorners;
        if (*gg->xrect != NULL)
            free(*gg->xrect);
        *gg->xrect = malloc(sizeof(double) * gg->ncorners);
        if (*gg->yrect != NULL)
            free(*gg->yrect);
        *gg->yrect = malloc(sizeof(double) * gg->ncorners);

        for (i = 0; i < gg->ncorners; ++i) {
            (*gg->xrect)[i] = creal(gg->newrzs[i]);
            (*gg->yrect)[i] = cimag(gg->newrzs[i]);
        }
    }

    free(done);
}

static int z2q_simple(gridgen* gg, zdouble z, zdouble zs[], double eps)
{
    int nq = gg->nquadrilaterals;
    int i, qi;

    for (i = 0, qi = gg->lastqi; i < nq; ++i, qi = (qi + 1) % nq) {
        quadrilateral* q = &gg->quadrilaterals[qi];

        if (in_quadrilateral(&z, zs, q->vids, eps)) {
            gg->lastqi = qi;
            return qi;
        }
    }

    return -1;
}

static int z2q(gridgen* gg, zdouble z, zdouble zs[], double eps)
{
    if (gg->nppe > 0) {
        int nq = gg->nquadrilaterals;
        int i, qi;

        for (i = 0, qi = gg->lastqi; i < nq; ++i, qi = (qi + 1) % nq) {
            if (in_poly(z, gg->nqivertices[qi], &gg->qivertices[qi * gg->nppq], eps)) {
                gg->lastqi = qi;
                return qi;
            }
        }
    } else
        return z2q_simple(gg, z, zs, eps);

    return -1;
}

static void output_point(gridgen* gg, double x, double y)
{
#if !defined(GRIDGEN_STANDALONE) && defined(HAVE_GRIDNODES_H)
    if (gg->gn == NULL) {
#endif
        if (!isnan(x))
            fprintf(gg->out, "%.10g %.10g\n", x, y);
        else
            fprintf(gg->out, "NaN NaN\n");
        fflush(gg->out);
#if !defined(GRIDGEN_STANDALONE) && defined(HAVE_GRIDNODES_H)
    } else
        gridnodes_readnextpoint(gg->gn, x, y);
#endif
}

static void map_point(gridgen* gg, zdouble z)
{
    zdouble zz = NaN;
    double eps = gg->eps;
    quadrilateral* qs = gg->quadrilaterals;
    int nz = gg->nz;
    zdouble* zs = gg->newzs;
    zdouble* ws = gg->ws;

    int status;
    zdouble w;

    int qid = z2q(gg, z, zs, eps);
    int vid = -1;
    int vidmin = -1;
    double distmin = DBL_MAX;
    int* vids;
    int i;

    if (qid < 0)
        quit("map_point(): could not map z = (%.10g,%.10g) to a quadrilateral\n", creal(z), cimag(z));

    vids = qs[qid].vids;

    /*
     * find the closest vertex
     */
    for (i = 0; i < 4; ++i) {
        double dist = cabs(z - zs[vids[i]]);

        if (dist < eps) {
            w = gg->ws[qid * nz + vids[i]];
            vid = vids[i];
            status = 0;
            goto know_w;
        }
        if (dist < distmin) {
            vidmin = vids[i];
            distmin = dist;
        }
    }

    /*
     * (Note that because we eliminated the possibility of z being the nearest
     * vertex, w can not be a transform of a polygon vertex. The only danger
     * is that because of approximate nature of z2q() we could get a wrong
     * quadrilateral.)
     */
    w = sc_z2w(gg->newsc, z, zs[vidmin], ws[nz * qid + vidmin], gg->newAs[qid], ZZERO, gg->newAs[qid], ZZERO, gg->newBs[qid], &ws[qid * nz], eps * 10.0, &status);

  know_w:
    if (status < 2) {
        zz = sc_w2z(gg->sc, w, vid, ZZERO, gg->As[qid], -1, gg->Bs[qid], &ws[qid * nz]);
        if (status != 0)
            (void) z2q_simple(gg, zz, gg->zs, eps);
        output_point(gg, creal(zz), cimag(zz));
    } else {
        output_point(gg, NaN, NaN);
    }
    if (gg_verbose && gg->out != stdout) {
        if (gg_verbose > 1)
            fprintf(stderr, "(%.10g,%.10g)(%d)->(%.10f,%.10f)->(%.10g,%.10g)\n", creal(z), cimag(z), qid, creal(w), cimag(w), creal(zz), cimag(zz));
        else {
            if (status == 0)
                fprintf(stderr, ".");
            else if (status == 1)
                fprintf(stderr, "o");
            else
                fprintf(stderr, "-");
        }
        fflush(stderr);
    }
}

static void align_onboundaries(gridgen* gg, zdouble* z)
{
    double eps = gg->eps;
    double x = creal(z[0]);
    double y = cimag(z[0]);

    if (fabs(x - gg->newxmin) < eps)
        x = gg->newxmin;
    else if (fabs(x - gg->newxmax) < eps)
        x = gg->newxmax;
    if (fabs(y - gg->newymin) < eps)
        y = gg->newymin;
    else if (fabs(y - gg->newymax) < eps)
        y = gg->newymax;

    z[0] = x + I * y;
}

static void map_quadrilaterals(gridgen* gg)
{
    int nppe = gg->nppe;
    int nq = gg->nquadrilaterals;

    /*
     * (there is 1-to-1 correspondence between inner edges of the
     * triangulation and the diagonals of the quadrilaterals) 
     */
    zdouble* zm_all = malloc(nppe * nq * sizeof(zdouble));

    int maxdiff = gg->nz - 1;
    double eps = gg->eps;
    int nz = gg->nz;
    zdouble* zs = gg->zs;
    zdouble* ws = gg->ws;
    hashtable* e2q = NULL;      /* internal edge to quadrilateral */

    int i, j;

    if (gg_verbose)
        fprintf(stderr, "mapping quadrilaterals (nppe = %d):\n", gg->nppe);
    if (gg_verbose == 1)
        fprintf(stderr, "  ");

    for (i = 0; i < nq; ++i) {
        quadrilateral* q = &gg->quadrilaterals[i];
        zdouble* zm = &zm_all[nppe * i];
        zdouble z0, z1;

        if (gg_verbose > 1)
            fprintf(stderr, "  %d:\n", i);
        /*
         * end points of a diagonal 
         */
        z0 = gg->zs[q->vids[0]];
        z1 = gg->zs[q->vids[2]];

        for (j = 0; j < nppe; ++j) {
            zdouble z = (z0 * (nppe - j) + z1 * (j + 1)) / (nppe + 1);
            int vid = (j < nppe / 2) ? q->vids[0] : q->vids[2];
            int status = 0;

            zdouble w = sc_z2w(gg->sc, z, zs[vid], ws[nz * i + vid], gg->As[i], ZZERO, gg->As[i], ZZERO, gg->Bs[i], &ws[i * nz], eps * 10.0, &status);

            assert(status < 2); /* I can see no reason why the mapping could
                                 * fail here. Bail out. */

            zm[j] = sc_w2z(gg->newsc, w, -1, ZZERO, gg->newAs[i], -1, gg->newBs[i], &ws[i * nz]);

            if (gg_verbose == 1) {
                if (status == 0)
                    fprintf(stderr, ".");
                else if (status == 1)
                    fprintf(stderr, "o");
                else
                    fprintf(stderr, "-");
            }
        }

        if (gg_verbose > 1) {
            z0 = gg->newzs[q->vids[0]];
            z1 = gg->newzs[q->vids[1]];

            fprintf(stderr, "    (%.15g %.15g)\n", creal(z0), cimag(z0));
            for (j = 0; j < nppe; ++j)
                fprintf(stderr, "    (%.15g %.15g)\n", creal(zm[j]), cimag(zm[j]));
            fprintf(stderr, "    (%.15g %.15g)\n", creal(z1), cimag(z1));
        }
    }
    if (gg_verbose == 1)
        fprintf(stderr, "\n");

    /*
     * a bit of muscle building
     */
    e2q = ht_create_i2(gg->nquadrilaterals);
    for (i = 0; i < nq; ++i) {
        quadrilateral* q = &gg->quadrilaterals[i];
        int key[2];

        if (q->vids[0] < q->vids[2]) {
            key[0] = q->vids[0];
            key[1] = q->vids[2];
        } else {
            key[0] = q->vids[2];
            key[1] = q->vids[0];
        }

        /*
         * make sure that the edge generates a unique key
         */
        assert(ht_insert(e2q, &key, q) == NULL);
    }

    /*
     * storing the found images
     */
    gg->nqivertices = calloc(nq, sizeof(int));
    gg->qivertices = malloc(nq * gg->nppq * sizeof(zdouble));
    for (i = 0; i < nq; ++i) {
        quadrilateral* q = &gg->quadrilaterals[i];
        zdouble* qivertices = &gg->qivertices[i * gg->nppq];
        zdouble* zm;
        int edge;

        for (edge = 0; edge <4; ++edge) {
            int vid0 = q->vids[edge];
            int vid1 = q->vids[(edge +1) %4];
            int diff = abs(vid0 - vid1);

            /*
             * store the vertex
             */
            qivertices[gg->nqivertices[i]] = gg->newzs[vid0];
            gg->nqivertices[i]++;

            /*
             * if necessary, store the rest 
             */
            if (diff != 1 && diff != maxdiff) {
                int key[2];
                quadrilateral* q1;
                int qid;
                int j;

                if (vid0 < vid1) {
                    key[0] = vid0;
                    key[1] = vid1;
                } else {
                    key[0] = vid1;
                    key[1] = vid0;
                }
                q1 = ht_find(e2q, &key);
                assert(q1 != NULL);     /* found in the table */
                qid = q1->id;

                zm = &zm_all[nppe * qid];
                if (vid0 == gg->quadrilaterals[qid].vids[0]) {
                    for (j = 0; j < nppe; ++j) {
                        qivertices[gg->nqivertices[i]] = zm[j];
                        gg->nqivertices[i]++;
                    }
                } else {
                    for (j = nppe - 1; j >= 0; --j) {
                        qivertices[gg->nqivertices[i]] = zm[j];
                        gg->nqivertices[i]++;
                    }
                }
            }
        }
        /*
         * close up the polygon
         */
        qivertices[gg->nqivertices[i]] = gg->newzs[q->vids[0]];
        gg->nqivertices[i]++;
    }

    ht_destroy(e2q);
    free(zm_all);
}

static void generate_grid(gridgen* gg)
{
    int count = 0;
    int i, j;

    if (gg_verbose) {
        fprintf(stderr, "generating grid:\n");
        fflush(stderr);
    }

    ode_silent = 1;
    ode_stoponnan = 0;

    if (gg_verbose == 1 && gg->out != stdout)
        fprintf(stderr, "  ");

    if (gg->gridpoints == NULL) {
        int nx = gg->nx;
        int ny = gg->ny;
        double dzx = (gg->newxmax - gg->newxmin) / (nx - 1);
        double dzy = (gg->newymax - gg->newymin) / (ny - 1);
        zdouble z0 = gg->newxmin + I * gg->newymin;

        if (gg->out != NULL)
            fprintf(gg->out, "## %d x %d\n", nx, ny);

        for (j = 0; j < ny; ++j) {
            zdouble zy = z0 + I * dzy * j;

            for (i = 0; i < nx; ++i) {
                zdouble z = zy + dzx * i;

                align_onboundaries(gg, &z);
                if (in_poly(z, gg->ncorners, gg->newrzs, gg->eps)) {
                    map_point(gg, z);
                    count++;
                } else
                    output_point(gg, NaN, NaN);
            }
        }
    } else {
        double dzx = gg->newxmax - gg->newxmin;
        double dzy = gg->newymax - gg->newymin;
        zdouble z0 = gg->newxmin + I * gg->newymin;
        int n = gg->ngridpoints;
        zdouble* zs = (zdouble*) gg->gridpoints;

        for (i = 0; i < n; ++i) {
            zdouble z = z0 + creal(zs[i]) * dzx + I * cimag(zs[i]) * dzy;

            align_onboundaries(gg, &z);
            if (in_poly(z, gg->ncorners, gg->newrzs, gg->eps)) {
                map_point(gg, z);
                count++;
            } else
                output_point(gg, NaN, NaN);
        }
    }

    if (gg_verbose)
        fprintf(stderr, " (%d nodes)\n", count);
}

void gridgen_generategrid(char* prmfname)
{
    gridgen* gg = NULL;

    gg = gridgen_create(prmfname);

    /*
     * get turn points 
     */
    get_rpoints(gg);

    /*
     * check that it triangulates OK 
     */
    gg->d = vertlist_triangulate(gg->vertices, stderr);
    delaunay_destroy(gg->d);

    /*
     * insert extra vertices -- phase 1 
     */
    vertlist_process_phase1(gg->vertices);

    /*
     * see how many triangles are there now 
     */
    if (gg_verbose) {
        gg->d = vertlist_triangulate(gg->vertices, stderr);
        delaunay_destroy(gg->d);
    }

    /*
     * insert extra vertices -- phase 2
     */
    vertlist_process_phase2(gg->vertices);      /* recursive */

    /*
     * final triangulation 
     */
    gg->d = vertlist_triangulate(gg->vertices, stderr);

    /*
     * calculate vertex angles 
     */
    gg->betas = calculate_betas(gg->vertices);

    /*
     * check sum of angles 
     */
    gridgen_check(gg);

    /*
     * get quadrilaterals from the triangulation 
     */
    gg->quadrilaterals = set_quadrilaterals(gg->d);
    gg->nquadrilaterals = gg->d->npoints - 3;

    gg->nz = gg->vertices->n;
    gg->zs = vertlist_tozdouble(gg->vertices);

    /*
     * get corner vertex indices for image polygon 
     */
    get_rvids(gg);

    /*
     * calculating log|ro| for each quadrilateral 
     */
    gg->cis = calculate_cis(gg->nquadrilaterals, gg->quadrilaterals, gg->d);

    /*
     * create Schwarz-Christoffel transform 
     */
    gg->sc = sc_create(gg->vertices->n, gg->nnodes, gg->betas);

    /*
     * solve nonlinear system 
     */
    find_sigmas(gg, F);

    /*
     * set vertex angles for image polygon 
     */
    set_newbetas(gg);

    /*
     * create Schwarz-Christoffel transform for the image 
     */
    gg->newsc = sc_create(gg->vertices->n, gg->nnodes, gg->newbetas);

    /*
     * calculate image polygon 
     */
    calculate_newzs(gg);

    if (gg->nppe > 0)
        map_quadrilaterals(gg);

    generate_grid(gg);

    gridgen_destroy(gg);
}

#if !defined (GRIDGEN_STANDALONE) && defined(HAVE_GRIDNODES_H)

gridnodes* gridgen_generategrid2(int nbdry, double xbdry[], double ybdry[], double beta[], int ul, int nx, int ny, int ngrid, double xgrid[], double ygrid[], int nnodes, int newton, double precision, int checksimplepoly, int thin, int nppe, int verbose, int* nsigmas, double** sigmas, int* nrect, double** xrect, double** yrect)
{
    gridgen* gg = NULL;
    gridnodes* gn = NULL;

    gg = gridgen_create2(nbdry, xbdry, ybdry, beta, ul, nx, ny, ngrid, xgrid, ygrid, nnodes, newton, precision, checksimplepoly, thin, nppe, verbose, nsigmas, sigmas, nrect, xrect, yrect);

    /*
     * get turn points 
     */
    get_rpoints(gg);

    /*
     * check that it triangulates OK 
     */
    gg->d = vertlist_triangulate(gg->vertices, stderr);
    delaunay_destroy(gg->d);

    /*
     * insert extra vertices -- phase 1 
     */
    vertlist_process_phase1(gg->vertices);

    /*
     * see how many triangles are there now 
     */
    if (gg_verbose) {
        gg->d = vertlist_triangulate(gg->vertices, stderr);
        delaunay_destroy(gg->d);
    }

    /*
     * insert extra vertices -- phase 2
     */
    vertlist_process_phase2(gg->vertices);      /* recursive */

    /*
     * final triangulation 
     */
    gg->d = vertlist_triangulate(gg->vertices, stderr);

    /*
     * calculate vertex angles 
     */
    gg->betas = calculate_betas(gg->vertices);

    /*
     * check sum of angles 
     */
    gridgen_check(gg);

    /*
     * get quadrilaterals from the triangulation 
     */
    gg->quadrilaterals = set_quadrilaterals(gg->d);
    gg->nquadrilaterals = gg->d->npoints - 3;

    gg->nz = gg->vertices->n;
    gg->zs = vertlist_tozdouble(gg->vertices);

    /*
     * get corner vertex indices for image polygon 
     */
    get_rvids(gg);

    /*
     * calculating log|ro| for each quadrilateral 
     */
    gg->cis = calculate_cis(gg->nquadrilaterals, gg->quadrilaterals, gg->d);

    /*
     * create Schwarz-Christoffel transform 
     */
    gg->sc = sc_create(gg->vertices->n, gg->nnodes, gg->betas);

    /*
     * solve nonlinear system 
     */
    find_sigmas(gg, F);

    /*
     * set vertex angles for image polygon 
     */
    set_newbetas(gg);

    /*
     * create Schwarz-Christoffel transform for the image 
     */
    gg->newsc = sc_create(gg->vertices->n, gg->nnodes, gg->newbetas);

    /*
     * calculate image polygon 
     */
    calculate_newzs(gg);

    if (gg->nppe > 0)
        map_quadrilaterals(gg);

    generate_grid(gg);

    gn = gg->gn;

    gridgen_destroy(gg);

    return gn;
}

#endif
