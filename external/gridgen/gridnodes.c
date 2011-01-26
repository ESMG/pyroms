/******************************************************************************
 *
 *  File:           gridnodes.c
 *  
 *  Created         22/01/2002
 *  
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *  
 *  Purpose:        Handling of grid node arrays
 *  Revisions:      14 Feb 2007 PS: added gridnodes_readnextpoint()
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>
#include <errno.h>
#include <assert.h>
#include "nan.h"
#include "gridnodes.h"
#include "gridmap.h"
#include "gucommon.h"

#define BUFSIZE 10240

/* Deviation from Orthogonality = 90 - theta
 * Aspect Ratio = max(dx,dy) / min(dx,dy)
 */
typedef struct {
    double mdo;                 /* maximum deviation from orthogonality */
    int imdo;
    int jmdo;
    double ado;                 /* average deviation from orthogonality */
    double mar;                 /* maximum aspect ratio */
    double aar;                 /* average aspect ratio */
} gridstats;

struct gridnodes {
    int nx;
    int ny;
    double** gx;
    double** gy;
    NODETYPE type;
    int validated;
    gridstats* stats;
    int nextpoint;
};

char* nodetype2str[] = {
    "not specified",
    "double density",
    "center",
    "corner"
};

/* Constructor. Reads double density grid nodes into arrays of X and Y
 * coordinates.
 * @param fname File name with grid nodes; can be "stdin"
 * @param type Node type
 * @return Gridnodes structure
 */
gridnodes* gridnodes_read(char* fname, NODETYPE type)
{
    gridnodes* gn = malloc(sizeof(gridnodes));
    FILE* f = NULL;
    int count;
    char buf[BUFSIZE];
    double* xx;
    double* yy;
    int i, j;

    if (gu_verbose)
        fprintf(stderr, "## grid input: reading from \"%s\"\n", fname);

    if (strcasecmp(fname, "stdin") == 0 || strcmp(fname, "-") == 0)
        f = stdin;
    else
        f = gu_fopen(fname, "r");

    gn->type = type;

    /*
     * get grid size 
     */
    if (fgets(buf, BUFSIZE, f) == NULL)
        gu_quit("%s: empty file\n", fname);

    if (sscanf(buf, "## %d x %d", &gn->nx, &gn->ny) != 2)
        gu_quit("%s: could not read grid size: expected header in \"## %%d x %%d\" format\n", fname);

    if (gu_verbose)
        fprintf(stderr, "##   %d x %d grid\n", gn->nx, gn->ny);

    if (gn->nx < 1)
        gu_quit("gridnodes_read(): nx = %d: invalid grid size\n", gn->nx);
    if (gn->ny < 1)
        gu_quit("gridnodes_read(): ny = %d: invalid grid size\n", gn->ny);
    if ((double) gn->nx * (double) gn->ny > (double) INT_MAX)
        gu_quit("gridnodes_read(): grid size (%d x %d) is too big\n", gn->nx, gn->ny);

    if (type == NT_DD) {
        if (gn->nx % 2 == 0)
            gu_quit("gridnodes_read(): nx = %d must be odd for double density grid nodes\n", gn->nx);
        if (gn->ny % 2 == 0)
            gu_quit("gridnodes_read(): ny = %d must be odd for double density grid nodes\n", gn->ny);
    }

    /*
     * allocate storage 
     */
    gn->gx = gu_alloc2d(gn->nx, gn->ny, sizeof(double));
    gn->gy = gu_alloc2d(gn->nx, gn->ny, sizeof(double));

    /*
     * read grid nodes 
     */
    for (j = 0, xx = gn->gx[0], yy = gn->gy[0], count = 0; j < gn->ny; ++j) {
        for (i = 0; i < gn->nx; ++i, ++xx, ++yy) {
            if (fgets(buf, BUFSIZE, f) == NULL)
                gu_quit("%s: could not read %d-th point (%d x %d points expected)\n", fname, j * gn->nx + i + 1, gn->nx, gn->ny);
            if (sscanf(buf, "%lf %lf", xx, yy) != 2) {
                *xx = NaN;
                *yy = NaN;
                continue;
            }
            if (!isnan(*xx))
                count++;
        }
    }

    if (gu_verbose) {
        fprintf(stderr, "##   %d non-empty grid nodes (%.1f%%)\n", count, 100.0 * count / gn->nx / gn->ny);
        fflush(stderr);
    }

    if (f != stdin)
        fclose(f);

    gn->validated = 0;
    gn->stats = NULL;

    return gn;
}

/* Constructor. Creates an empty grid.
 *
 * @param nx Number of columns
 * @param nx Number of rows
 * @param type Node type
 * @return Gridnodes structure
 */
gridnodes* gridnodes_create(int nx, int ny, NODETYPE type)
{
    gridnodes* gn = malloc(sizeof(gridnodes));

    gn->nx = nx;
    gn->ny = ny;
    gn->gx = gu_alloc2d(nx, ny, sizeof(double));
    gn->gy = gu_alloc2d(nx, ny, sizeof(double));
    gn->type = type;
    gn->validated = 0;
    gn->stats = NULL;
    gn->nextpoint = 0;

    return gn;
}

/* Allows to fill the created gridnodes object point by point.
 *
 * @param gn Grid nodes
 * @param x X coordinate
 * @param y Y coordinate
 */
void gridnodes_readnextpoint(gridnodes* gn, double x, double y)
{
    int j = gn->nextpoint / gn->nx;
    int i = gn->nextpoint % gn->nx;

    gn->gx[j][i] = x;
    gn->gy[j][i] = y;
    gn->nextpoint = (gn->nextpoint + 1) % (gn->nx * gn->ny);
}

/* Destructor.
 * @param gn Grid nodes
 */
void gridnodes_destroy(gridnodes* gn)
{
    if (gn->stats != NULL)
        free(gn->stats);
    gu_free2d(gn->gx);
    gu_free2d(gn->gy);
    free(gn);
}

/*
 * A grid generator calculates particular nodes in double density grid just as
 * a set of points. Therefore, there may be cells with only some of the
 * bounding double-density nodes valid (not NaNs), and such cells must be
 * marked as "non-valid".
 *
 * Following is a set of rules used for marking the cells in DD grid:
 *
 * 1. Gridwise, a cell is defined by its four corner nodes: a cell is valid if
 *    and only if all four corner cells are valid (not NaNs).
 * 2. A modeller may use the centre node to mark a cell as non-valid within
 *    an application (e.g. to mark it as a land cell or an outside cell).
 * 3. A modeler may use edge nodes (those in between corner cells) to mark
 *    an edge as non-valid within an application (e.g. to mark it as a 
 *    non-penetrable edge).
 *
 * The following validation procedures are supposed to exclude non-valid cells
 * from the grid.
 *
 * Note: theoretically, it is possible that one or more of the four edge nodes
 * or the centre node will come from the grid generator as NaNs, while all
 * four corner nodes come as valid nodes. For now, to restore these values from
 * the corner values, use gridnodes_transform() procedure to transform the grid
 * from NT_COR to NT_DD type (transform to NT_COR first if you start from NT_DD
 * grid).
 */

/* Validates double density grid nodes.
 * @param gn Grid nodes
 */
static void gridnodes_validate_dd(gridnodes* gn)
{
    int nx = gn->nx;
    int ny = gn->ny;
    double** x = gn->gx;
    double** y = gn->gy;
    int count;
    int i, j;

    /*
     * corner nodes:
     *
     * mark as non-valid if one or more of the other three corner nodes in each
     * of the four adjacent cells are non-valid; otherwise leave as is
     */
    for (j = 0; j < ny; j += 2) {
        for (i = 0; i < nx; i += 2) {
            if (j > 0 && i > 0 && !isnan(x[j - 2][i - 2]) && !isnan(x[j - 2][i]) && !isnan(x[j][i - 2]))
                continue;
            if (j + 2 < ny && i > 0 && !isnan(x[j + 2][i - 2]) && !isnan(x[j + 2][i]) && !isnan(x[j][i - 2]))
                continue;
            if (j > 0 && i + 2 < nx && !isnan(x[j - 2][i + 2]) && !isnan(x[j - 2][i]) && !isnan(x[j][i + 2]))
                continue;
            if (j + 2 < ny && i + 2 < nx && !isnan(x[j + 2][i + 2]) && !isnan(x[j + 2][i]) && !isnan(x[j][i + 2]))
                continue;
            x[j][i] = NaN;
            y[j][i] = NaN;
        }
    }

    /*
     * centre nodes:
     *
     * mark as non-valid if one or more of the four adjacent corner cells is
     * non-valid; otherwise leave as is
     *
     * (note that it is possible for a valid cell to have a non-valid edge
     * node; this can be a way to mark an impenetrable face)
     */
    for (j = 1; j < ny; j += 2) {
        for (i = 1; i < nx; i += 2) {
            if (isnan(x[j][i]))
                continue;
            if (isnan(x[j - 1][i - 1]) || isnan(x[j - 1][i + 1]) || isnan(x[j + 1][i - 1]) || isnan(x[j + 1][i + 1])) {
                x[j][i] = NaN;
                y[j][i] = NaN;
            }
        }
    }

    /*
     * edge nodes:
     *
     * mark as non-valid if one or more of the two adjacent corner cells are
     * non-valid; otherwise leave as is
     */
    for (j = 0; j < ny; j += 2) {
        for (i = 1; i < nx; i += 2) {
            if (isnan(x[j][i]))
                continue;
            if ((i - 1 >= 0 && isnan(x[j][i - 1])) || (i + 1 < nx && isnan(x[j][i + 1]))) {
                x[j][i] = NaN;
                y[j][i] = NaN;
            }
        }
    }
    for (j = 1; j < ny; j += 2) {
        for (i = 0; i < nx; i += 2) {
            if (isnan(x[j][i]))
                continue;
            if ((j - 1 >= 0 && isnan(x[j - 1][i])) || (j + 1 < ny && isnan(x[j + 1][i]))) {
                x[j][i] = NaN;
                y[j][i] = NaN;
            }
        }
    }

    /*
     * count valid cells 
     */
    if (gu_verbose) {
        for (j = 1, count = 0; j < ny; j += 2)
            for (i = 1; i < nx; i += 2)
                if (!isnan(x[j][i]))
                    count++;
        fprintf(stderr, "##   %d valid cells (%.1f%%)\n", count, count * 100.0 / (nx / 2) / (ny / 2));
    }
}

/* Validates corner grid nodes.
 * @param gn Grid nodes
 */
static void gridnodes_validate_cor(gridnodes* gn)
{
    int nx = gn->nx;
    int ny = gn->ny;
    double** x = gn->gx;
    double** y = gn->gy;
    int count;
    int i, j;

    for (j = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            if (i > 0 && j > 0 && !isnan(x[j - 1][i - 1]) && !isnan(x[j - 1][i]) && !isnan(x[j][i - 1]))
                continue;
            if (i > 0 && j < ny - 1 && !isnan(x[j][i - 1]) && !isnan(x[j + 1][i - 1]) && !isnan(x[j + 1][i]))
                continue;
            if (i < nx - 1 && j > 0 && !isnan(x[j - 1][i]) && !isnan(x[j - 1][i + 1]) && !isnan(x[j][i + 1]))
                continue;
            if (i < nx - 1 && j < ny - 1 && !isnan(x[j + 1][i]) && !isnan(x[j][i + 1]) && !isnan(x[j + 1][i + 1]))
                continue;
            /*
             * the node does not belong to any valid cell around it 
             */
            x[j][i] = NaN;
            y[j][i] = NaN;
        }
    }

    if (gu_verbose) {
        ny--;
        nx--;
        for (j = 0, count = 0; j < ny; ++j)
            for (i = 0; i < nx; ++i)
                if (!isnan(x[j + 1][i + 1]) && !isnan(x[j + 1][i]) && !isnan(x[j][i + 1]) && !isnan(x[j][i]))
                    count++;
        fprintf(stderr, "##   %d valid cells (%.1f%%)\n", count, count * 100.0 / nx / ny);
    }
}

/* Validates grid nodes.
 * Sets all nodes not belonging to any valid cells to NaNs.
 * @param gn Grid nodes
 */
void gridnodes_validate(gridnodes* gn)
{
    if (gu_verbose)
        fprintf(stderr, "## grid validation:\n");
    if (gn->type == NT_DD)
        gridnodes_validate_dd(gn);
    else if (gn->type == NT_COR)
        gridnodes_validate_cor(gn);
    else if (gu_verbose)
        fprintf(stderr, "## gridnodes_validate(): nothing to do for nodes of \"%s\" type\n", nodetype2str[gn->type]);

    gn->validated = 1;

    if (gu_verbose)
        fflush(stdout);
}

/* Makes a deep copy of grid nodes.
 * @param old Source grid nodes
 * @return Destination grid nodes
 */
gridnodes* gridnodes_copy(gridnodes* old)
{
    gridnodes* new = malloc(sizeof(gridnodes));

    new->nx = old->nx;
    new->ny = old->ny;
    new->type = old->type;
    new->gx = gu_alloc2d(old->nx, old->ny, sizeof(double));
    new->gy = gu_alloc2d(old->nx, old->ny, sizeof(double));
    memcpy(&new->gx[0][0], &old->gx[0][0], old->nx * old->ny * sizeof(double));
    memcpy(&new->gy[0][0], &old->gy[0][0], old->nx * old->ny * sizeof(double));

    return new;
}

/* Makes a deep copy of subgrid nodes with indices [imin:imax][jmin:jmax].
 * @param old Source grid nodes
 * @param imin Minimal i index
 * @param imax Maximal i index 
 * @param jmin Minimal j index
 * @param jmax Maximal j index 
 * @return Destination grid nodes
 */
gridnodes* gridnodes_subgrid(gridnodes* gn, int imin, int imax, int jmin, int jmax)
{
    gridnodes* new = NULL;
    int i, j, ii, jj;

    if (imin < 0)
        imin = 0;
    if (imax >= gn->nx)
        imax = gn->nx - 1;
    if (jmin < 0)
        jmin = 0;
    if (jmax >= gn->ny)
        jmax = gn->ny - 1;

    if (imin == 0 && imax == gn->nx - 1 && jmin == 0 && jmax == gn->ny - 1)
        return gridnodes_copy(gn);

    new = malloc(sizeof(gridnodes));

    new->nx = imax - imin + 1;
    new->ny = jmax - jmin + 1;

    new->gx = gu_alloc2d(new->nx, new->ny, sizeof(double));
    new->gy = gu_alloc2d(new->nx, new->ny, sizeof(double));

    for (j = jmin, jj = 0; j <= jmax; ++j, ++jj) {
        for (i = imin, ii = 0; i <= imax; ++i, ++ii) {
            new->gx[jj][ii] = gn->gx[j][i];
            new->gy[jj][ii] = gn->gy[j][i];
        }
    }

    return new;
}

/* Transforms grid nodes of one type into grid nodes of another type.
 * @param gn Grid nodes
 * @param type Type of new grid nodes
 * @return New grid nodes
 */
gridnodes* gridnodes_transform(gridnodes* gn, NODETYPE type)
{
    gridnodes* gn1 = NULL;
    int i, j, i1, j1;

    if (!gn->validated)
        gridnodes_validate(gn);

    if (gn->type == type || type == NT_NONE)
        return gridnodes_copy(gn);

    gn1 = malloc(sizeof(gridnodes));
    gn1->type = type;
    gn1->stats = NULL;

    if (gn->type == NT_DD) {
        if (type == NT_COR) {
            gn1->nx = gn->nx / 2 + 1;
            gn1->ny = gn->ny / 2 + 1;
            gn1->gx = gu_alloc2d(gn1->nx, gn1->ny, sizeof(double));
            gn1->gy = gu_alloc2d(gn1->nx, gn1->ny, sizeof(double));
            for (j = 0, j1 = 0; j < gn->ny; j += 2, ++j1) {
                for (i = 0, i1 = 0; i < gn->nx; i += 2, ++i1) {
                    gn1->gx[j1][i1] = gn->gx[j][i];
                    gn1->gy[j1][i1] = gn->gy[j][i];
                }
            }
        } else if (type == NT_CEN) {
            gn1->nx = gn->nx / 2;
            gn1->ny = gn->ny / 2;
            gn1->gx = gu_alloc2d(gn1->nx, gn1->ny, sizeof(double));
            gn1->gy = gu_alloc2d(gn1->nx, gn1->ny, sizeof(double));
            for (j = 1, j1 = 0; j < gn->ny; j += 2, ++j1) {
                for (i = 1, i1 = 0; i < gn->nx; i += 2, ++i1) {
                    gn1->gx[j1][i1] = gn->gx[j][i];
                    gn1->gy[j1][i1] = gn->gy[j][i];
                }
            }
        }
    } else if (gn->type == NT_COR) {
        if (type == NT_CEN) {
            gridmap* map = gridmap_build(gn->nx - 1, gn->ny - 1, gn->gx, gn->gy);

            gn1->nx = gn->nx - 1;
            gn1->ny = gn->ny - 1;
            gn1->gx = gu_alloc2d(gn1->nx, gn1->ny, sizeof(double));
            gn1->gy = gu_alloc2d(gn1->nx, gn1->ny, sizeof(double));

            /*
             * this may take a while 
             */
            for (i = 0; i < gn1->nx; ++i)
                for (j = 0; j < gn1->ny; ++j)
                    gridmap_fij2xy(map, i + 0.5, j + 0.5, &gn1->gx[j][i], &gn1->gy[j][i]);

            gridmap_destroy(map);
        } else if (type == NT_DD) {
            gridmap* map = gridmap_build(gn->nx - 1, gn->ny - 1, gn->gx, gn->gy);

            gn1->nx = gn->nx * 2 - 1;
            gn1->ny = gn->ny * 2 - 1;
            gn1->gx = gu_alloc2d(gn1->nx, gn1->ny, sizeof(double));
            gn1->gy = gu_alloc2d(gn1->nx, gn1->ny, sizeof(double));

            /*
             * this may take a while 
             */
            for (i = 0; i < gn1->nx; ++i)
                for (j = 0; j < gn1->ny; ++j)
                    gridmap_fij2xy(map, i / 2.0, j / 2.0, &gn1->gx[j][i], &gn1->gy[j][i]);

            gridmap_destroy(map);
        }
    } else if (gn->type == NT_CEN) {
        gu_quit("transforming nodes of type \"%s\" into \"%s\" is not supported\n", nodetype2str[gn->type], nodetype2str[type]);
    }

    gn1->validated = 1;         /* an internally generated grid is supposed
                                 * to be OK */

    return gn1;
}

void gridnodes_applymask(gridnodes* gn, int** mask)
{
    if (gn->type == NT_DD) {
        double** x = gn->gx;
        double** y = gn->gy;
        int nx = gn->nx;
        int ny = gn->ny;
        int i, j, ii, jj;

        for (j = 1, jj = 0; j < ny; j += 2, ++jj) {
            for (i = 1, ii = 0; i < nx; i += 2, ++ii) {
                if (mask[jj][ii] == 0) {
                    x[j][i] = NaN;
                    y[j][i] = NaN;
                }
            }
        }

        /*
         * corner nodes:
         *
         * mark as non-valid if all four center nodes of the adjacent cells are
         * non-valid
         */
        for (j = 0; j < ny; j += 2) {
            for (i = 0; i < nx; i += 2) {
                if (j > 0 && i > 0 && !isnan(x[j - 1][i - 1]))
                    continue;
                if (j + 1 < ny && i > 0 && !isnan(x[j + 1][i - 1]))
                    continue;
                if (j > 0 && i + 1 < nx && !isnan(x[j - 1][i + 1]))
                    continue;
                if (j + 1 < ny && i + 1 < nx && !isnan(x[j + 1][i + 1]))
                    continue;
                x[j][i] = NaN;
                y[j][i] = NaN;
            }
        }

        gridnodes_validate(gn);

    } else if (gn->type == NT_CEN) {
        double** x = gn->gx;
        double** y = gn->gy;
        int nx = gn->nx;
        int ny = gn->ny;
        int i, j;

        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                if (mask[j][i] == 0) {
                    x[j][i] = NaN;
                    y[j][i] = NaN;
                }
            }
        }
    } else if (gn->type == NT_COR)
        gu_quit("gridnodes_applymask(): applying mask to nodes of type \"%s\" is not supported\n", nodetype2str[NT_COR]);
    else if (gn->type == NT_NONE)
        gu_quit("gridnodes_applymask(): nodes type not specified\n");
}

/* Writes grid nodes into a file.
 * @param gn Grid nodes
 * @param fname File name; can be "stdout"
 * @param ctype Output coordinate type: CT_XY for XY, CT_X for X and CT_Y for Y
 */
void gridnodes_write(gridnodes* gn, char* fname, COORDTYPE ctype)
{
    FILE* f = NULL;
    double** x = gn->gx;
    double** y = gn->gy;
    int nx = gn->nx;
    int ny = gn->ny;
    int count = 0;
    int i, j;

    if (gu_verbose) {
        fprintf(stderr, "## grid output: writing to \"%s\"\n", fname);
        fprintf(stderr, "##   %d x %d grid\n", gn->nx, gn->ny);
    }

    if (!strcasecmp(fname, "stdout")) {
        f = stdout;
    } else
        f = gu_fopen(fname, "w");

    fprintf(f, "## %d x %d\n", nx, ny);

    if (ctype == CT_XY) {
        for (j = 0; j < ny; ++j) {
            for (i = 0; i < nx; ++i) {
                if (isnan(x[j][i]))
                    fprintf(f, "NaN NaN\n");
                else {
                    fprintf(f, "%.15g %.15g\n", x[j][i], y[j][i]);
                    count++;
                }
            }
        }
    } else if (ctype == CT_X) {
        for (j = 0; j < ny; ++j) {
            for (i = 0; i < nx; ++i) {
                if (isnan(x[j][i]))
                    fprintf(f, "NaN\n");
                else {
                    fprintf(f, "%.15g\n", x[j][i]);
                    count++;
                }
            }
        }
    } else if (ctype == CT_Y) {
        for (j = 0; j < ny; ++j) {
            for (i = 0; i < nx; ++i) {
                if (isnan(x[j][i]))
                    fprintf(f, "NaN\n");
                else {
                    fprintf(f, "%.15g\n", y[j][i]);
                    count++;
                }
            }
        }
    }

    if (gu_verbose)
        fprintf(stderr, "##   %d non-empty grid nodes (%.1f%%)\n", count, 100.0 * count / gn->nx / gn->ny);

    if (f != stdout)
        fclose(f);
    else
        fflush(stdout);
}

static void gridstats_init(gridstats* gs)
{
    gs->mdo = NaN;
    gs->imdo = -1;
    gs->jmdo = -1;
    gs->ado = NaN;
    gs->mar = NaN;
    gs->aar = NaN;
}

static double dtheta(double x1, double y1, double x2, double y2)
{
    double cos_th = (x1 * x2 + y1 * y2) / hypot(x1, y1) / hypot(x2, y2) / 2.0;
    double th = acos(cos_th) * 180.0 / M_PI;
    double dth = fabs(90.0 - th);

    return dth;
}

static void gridnodes_calcstats_cor(gridnodes* gn)
{
    int nx = gn->nx - 1;
    int ny = gn->ny - 1;
    double** x = gn->gx;
    double** y = gn->gy;
    double dor_sum = 0.0;
    double dor_max = 0.0;
    double ar_max = 1.0;
    double ar_sum = 0.0;
    int ncell = 0;
    int imdo = -1;
    int jmdo = -1;
    int i, j;

    if (gn->stats == NULL)
        gn->stats = malloc(sizeof(gridstats));

    gridstats_init(gn->stats);

    for (j = 0; j < ny; ++j)
        for (i = 0; i < nx; ++i)
            if (!isnan(x[j + 1][i + 1]) && !isnan(x[j + 1][i]) && !isnan(x[j][i + 1]) && !isnan(x[j][i])) {
                double dor[4];
                double ar;
                int k;

                dor[0] = dtheta(x[j][i + 1] - x[j][i], y[j][i + 1] - y[j][i], x[j + 1][i] - x[j][i], y[j + 1][i] - y[j][i]);
                dor[1] = dtheta(x[j + 1][i + 1] - x[j][i + 1], y[j + 1][i + 1] - y[j][i + 1], x[j][i + 1] - x[j][i], y[j][i + 1] - y[j][i]);
                dor[2] = dtheta(x[j + 1][i] - x[j + 1][i + 1], y[j + 1][i] - y[j + 1][i + 1], x[j + 1][i + 1] - x[j][i + 1], y[j + 1][i + 1] - y[j][i + 1]);
                dor[3] = dtheta(x[j + 1][i] - x[j + 1][i + 1], y[j + 1][i] - y[j + 1][i + 1], x[j + 1][i] - x[j][i], y[j + 1][i] - y[j][i]);

                for (k = 0; k < 4; ++k) {
                    if (dor[k] > dor_max) {
                        dor_max = dor[k];
                        imdo = i;
                        jmdo = j;
                    }
                    dor_sum += dor[k];
                }

                ar = hypot(x[j][i] + x[j][i + 1] - x[j + 1][i] - x[j + 1][i + 1], y[j][i] + y[j][i + 1] - y[j + 1][i] - y[j + 1][i + 1]) / hypot(x[j][i] + x[j + 1][i] - x[j][i + 1] - x[j + 1][i + 1], y[j][i] + y[j + 1][i] - y[j][i + 1] - y[j + 1][i + 1]);
                if (ar < 1.0)
                    ar = 1.0 / ar;
                if (ar > ar_max)
                    ar_max = ar;
                ar_sum += ar;

                ncell++;
            }

    gn->stats->mdo = dor_max;
    gn->stats->imdo = imdo;
    gn->stats->jmdo = jmdo;
    gn->stats->ado = dor_sum / ncell / 4.0;
    gn->stats->mar = ar_max;
    gn->stats->aar = ar_sum / ncell;
}

void gridnodes_calcstats(gridnodes* gn)
{
    gridnodes* gn1;

    if (gn->stats == NULL)
        gn->stats = malloc(sizeof(gridstats));

    gridstats_init(gn->stats);

    if (gn->type == NT_CEN || gn->type == NT_NONE) {
        if (gu_verbose)
            fprintf(stderr, "## gridnodes_calcstats(): do not know what to do with nodes of \"%s\" type\n", nodetype2str[gn->type]);
        return;
    }

    if (!gn->validated)
        gridnodes_validate(gn);

    gn1 = (gn->type == NT_COR) ? gn : gridnodes_transform(gn, NT_COR);
    gridnodes_calcstats_cor(gn1);

    if (gn->type != NT_COR) {
        memcpy(gn->stats, gn1->stats, sizeof(gridstats));
        gridnodes_destroy(gn1);
    }

    if (gu_verbose) {
        fprintf(stderr, "## maximum deviation from orthogonality = %.3f deg, in cell (%d,%d)\n", gn->stats->mdo, gn->stats->imdo, gn->stats->jmdo);
        fprintf(stderr, "## mean deviation from orthogonality = %.3f deg\n", gn->stats->ado);
        fprintf(stderr, "## maximum aspect ratio = %.3f\n", gn->stats->mar);
        fprintf(stderr, "## mean aspect ratio = %.3f\n", gn->stats->aar);
    }
}

int gridnodes_getnx(gridnodes* gn)
{
    return gn->nx;
}

int gridnodes_getny(gridnodes* gn)
{
    return gn->ny;
}

double** gridnodes_getx(gridnodes* gn)
{
    return gn->gx;
}

double** gridnodes_gety(gridnodes* gn)
{
    return gn->gy;
}

int gridnodes_getnce1(gridnodes* gn)
{
    if (gn->type == NT_DD)
        return (gn->nx - 1) / 2;
    else if (gn->type == NT_COR)
        return gn->nx - 1;
    else if (gn->type == NT_CEN)
        return gn->nx;
    else
        gu_quit("gridnodes_getnce1(): node type not specified\n");

    return 0;
}

int gridnodes_getnce2(gridnodes* gn)
{
    if (gn->type == NT_DD)
        return (gn->ny - 1) / 2;
    else if (gn->type == NT_COR)
        return gn->ny - 1;
    else if (gn->type == NT_CEN)
        return gn->ny;
    else
        gu_quit("gridnodes_getnce2(): node type not specified\n");

    return 0;
}
