/******************************************************************************
 *
 * File:           gridmap.c
 *  
 * Created:        Fri Feb 19 11:36:24 EST 1992 (as xytoij.c)
 *  
 * Author:         Daniel Delbourgo/Stephen Walker
 *                 CSIRO Division of Oceanography
 *
 * Purpose:        Calculates transformations between physical and index
 *                 space within a numerical grid
 *
 * Revisions:      110898 JRW
 *                 Added IJtoXY conversion and
 *                 fractional XYtoIJ and IJtoXY conversion
 *
 *                 2000 Pavel Sakov
 *                 Added branch calculation to handle both right- and
 *                 left-handed grids.
 *
 *                 April 2002 Pavel Sakov
 *                 Major mods to handle topologically non-rectangular grids
 *                 Renamed to gridmap.c
 *
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "poly.h"
#include "gridmap.h"
#include "gucommon.h"

#define EPS 1.0e-8
#define EPS_ZERO 1.0e-5
#define EPS_COMPACT 1.0e-10

typedef struct subgrid {
    gridmap* gmap;              /* gridf map this subgrid belongs to */
    poly* bound;                /* boundary polygon */
    int mini;                   /* minimal i index within the subgrid */
    int maxi;                   /* maximal i index within the subgrid */
    int minj;                   /* minimal j index within the subgrid */
    int maxj;                   /* maximal j index within the subgrid */
    struct subgrid* half1;      /* child 1 */
    struct subgrid* half2;      /* child 2 */
} subgrid;

struct gridmap {
    poly* bound;                /* boundary polygon */
    subgrid* trunk;             /* binary tree trunk */
    int nleaves;                /* for debugging purposes */
    int nce1;                   /* number of cells in e1 direction */
    int nce2;                   /* number of cells in e2 direction */
    double** gx;                /* reference to array of X coords
                                 * [nce2+1][nce2+1] */
    double** gy;                /* reference to array of Y coords
                                 * [nce2+1][nce2+1] */
    int sign;                   /* sqrt branch */
};

/* Creates a subgrid.
 * @param gm Grid map
 * @param pl Boundary polygon for the subgrid
 * @param mini Minimal i index within the subgrid
 * @param maxi Maximal i index within the subgrid
 * @param minj Minimal j index within the subgrid
 * @param maxj Maximal j index within the subgrid
 * @return Subgrid
 */
static subgrid* subgrid_create(gridmap* gm, poly* pl, int i1, int i2, int j1, int j2)
{
    subgrid* l = malloc(sizeof(subgrid));

    double** gx = gm->gx;
    double** gy = gm->gy;
    int n = pl->n;
    double x, y;
    int i = i1;                 /* to eliminate warning */
    int j, ii;

    l->bound = pl;
    l->gmap = gm;
    l->mini = INT_MAX;
    l->maxi = INT_MIN;
    l->minj = INT_MAX;
    l->maxj = INT_MIN;
    l->half1 = NULL;
    l->half2 = NULL;

    if (n == 0)
        return l;

    x = pl->x[0];
    y = pl->y[0];

    for (j = j1; j <= j2; ++j)
        for (i = i1; i <= i2; ++i)
            if (x == gx[j][i] && y == gy[j][i])
                goto aftersearch;

  aftersearch:

    if (j > j2)
        gu_quit("subgrid_create(): boundary vertex not in the grid\n");

    l->mini = i;
    l->maxi = i;
    l->minj = j;
    l->maxj = j;

    for (ii = 1; ii < n; ++ii) {
        x = pl->x[ii];
        y = pl->y[ii];

        if (i > i1 && x == gx[j][i - 1] && y == gy[j][i - 1])
            i--;
        else if (i < i2 && x == gx[j][i + 1] && y == gy[j][i + 1])
            i++;
        else if (j > j1 && x == gx[j - 1][i] && y == gy[j - 1][i])
            j--;
        else if (j < j2 && x == gx[j + 1][i] && y == gy[j + 1][i])
            j++;
        else if (x == gx[j][i] && y == gy[j][i])
            continue;
        else
            gu_quit("subgrid_create(): boundary vertex not in the grid\n");

        if (l->mini > i)
            l->mini = i;
        if (l->maxi < i)
            l->maxi = i;
        if (l->minj > j)
            l->minj = j;
        if (l->maxj < j)
            l->maxj = j;
    }

    return l;
}

/* Destroys a subgrid.
 * @param l Subgrid
 */
static void subgrid_destroy(subgrid* l)
{
    if (l->half1 != NULL)
        subgrid_destroy(l->half1);
    if (l->half2 != NULL)
        subgrid_destroy(l->half2);
    poly_destroy(l->bound);
    free(l);
}

/* Cuts boundary polygon in two. 
 * The cut goes either horizontally ([fixed][changes]) or vertically 
 * ([changes][fixed]) in index space; the physical nodes are given by
 * input double arrays; first two intersections of the cutting polyline
 * with the polygon are used to form the new polygons.
 * @param pl Original polygon
 * @param gx Array of x cell corner coordinates
 * @param gy Array of y cell corner coordinates
 * @param horiz flag: 1 for horizontal cut; 0 otherwise
 * @param index Value of "fixed" index
 * @param start Start value of "variable" index
 * @param end End value of "variable" index
 * @param pl1 Output polygon 1
 * @param pl1 Output polygon 2
 */
static void cut_boundary(poly* pl, double** gx, double** gy, int horiz, int index, int start, int end, poly** pl1, poly** pl2)
{
    int n = pl->n;
    int i = -1;
    int i1 = -1;                /* array index of the first intersection */
    int i2 = -1;                /* array index of the second intersection */
    int ii1 = -1;               /* polygon index of the first intersection */
    int ii2 = -1;               /* polygon index of the second intersection */
    int tmp = -1;

    /*
     * if the polygon has been explicitely closed, ignore the last point 
     */
    if (poly_isclosed(pl, 1.0e-15))
        n--;

    if (horiz) {                /* horizontal cut */
        /*
         * find first intersection 
         */
        for (i = start; i < end; ++i) {
            ii1 = poly_findindex(pl, gx[index][i], gy[index][i]);
            if (ii1 < 0)
                /*
                 * this node does not belong to the boundary 
                 */
                continue;
            /*
             * this node belongs to the boundary polygon To accept it as a
             * valid intersection, we must ensure that the very next node is 
             * either inside the polygon or belongs to the boundary but is
             * not the next boundary node. 
             */
            tmp = poly_findindex(pl, gx[index][i + 1], gy[index][i + 1]);
            if ((tmp < 0 && poly_containspoint(pl, gx[index][i + 1], gy[index][i + 1])) || (tmp >= 0 && abs(tmp - ii1) > 1 && abs(tmp - ii1) < n - 1))
                break;
        }

        /*
         * how the for-cycle ended 
         */
        if (i < end)            /* ok */
            i1 = i;
        else                    /* no intersection found */
            return;

        /*
         * find second intersection start from the node next to the first
         * intersection 
         */
        for (i = i1 + 1; i <= end; ++i) {
            ii2 = poly_findindex(pl, gx[index][i], gy[index][i]);
            if (ii2 >= 0)
                /*
                 * this node must be inside the boundary polygon -- skip 
                 */
                break;
        }

        if (ii2 < 0)
            /*
             * no intersection found 
             */
            return;
        else
            /*
             * ok 
             */
            i2 = i;

        /*
         * we found all necessary details, now form the new polygons 
         */
        *pl1 = poly_create();
        *pl2 = poly_create();

        /*
         * add the portion of perimeter 
         */
        for (i = ii1; i != ii2; i = (i + 1) % n)
            poly_addpoint(*pl1, pl->x[i], pl->y[i]);
        /*
         * add the cutting section 
         */
        for (i = i2; i > i1; --i)
            poly_addpoint(*pl1, gx[index][i], gy[index][i]);

        /*
         * add the portion of perimeter 
         */
        for (i = ii2; i != ii1; i = (i + 1) % n)
            poly_addpoint(*pl2, pl->x[i], pl->y[i]);
        /*
         * add the cutting section 
         */
        for (i = i1; i < i2; ++i)
            poly_addpoint(*pl2, gx[index][i], gy[index][i]);

    } else {                    /* vertical cut */
        for (i = start; i < end; ++i) {
            ii1 = poly_findindex(pl, gx[i][index], gy[i][index]);
            if (ii1 < 0)
                continue;
            tmp = poly_findindex(pl, gx[i + 1][index], gy[i + 1][index]);
            if ((tmp < 0 && poly_containspoint(pl, gx[i + 1][index], gy[i + 1][index])) || (tmp >= 0 && abs(tmp - ii1) > 1 && abs(tmp - ii1) < n - 1))
                break;
        }

        if (i < end)
            i1 = i;
        else
            return;

        for (i = i1 + 1; i <= end; ++i) {
            ii2 = poly_findindex(pl, gx[i][index], gy[i][index]);
            if (ii2 >= 0)
                break;
        }

        if (ii2 < 0)
            return;
        else
            i2 = i;

        *pl1 = poly_create();
        *pl2 = poly_create();

        for (i = ii1; i != ii2; i = (i + 1) % n)
            poly_addpoint(*pl1, pl->x[i], pl->y[i]);
        for (i = i2; i > i1; --i)
            poly_addpoint(*pl1, gx[i][index], gy[i][index]);

        for (i = ii2; i != ii1; i = (i + 1) % n)
            poly_addpoint(*pl2, pl->x[i], pl->y[i]);
        for (i = i1; i < i2; ++i)
            poly_addpoint(*pl2, gx[i][index], gy[i][index]);
        /*
         * There used to be closure of the polylines here:
         *        poly_close(*pl1);
         *        poly_close(*pl2);
         * While this does not harm, it causes the boundary of a quadrilateral
         * to contain 5 points rather than 4, which caused allocation for 8
         * points... Anyway, considering the way search in poly_containspoint()
         * works, this is not necessary.
         */
    }
}

/* Divides a subgrid in two.
 * @param sg The subgrid to divide
 * @param subgrid1 Output subgrid 1
 * @param subgrid2 Output subgrid 2
 * @param gx Array of x cell corner coordinates
 * @param gy Array of y cell corner coordinates
 */
static void subgrid_divide(subgrid* sg, subgrid** sg1, subgrid** sg2)
{
    poly* pl1 = NULL;
    poly* pl2 = NULL;
    gridmap* gm = sg->gmap;
    int index;

    if ((sg->maxi <= sg->mini + 1) && (sg->maxj <= sg->minj + 1)) {
        *sg1 = *sg2 = NULL;
        return;
    }

    if (sg->maxi - sg->mini > sg->maxj - sg->minj) {
        /*
         * divide "vertically" 
         */
        index = (sg->mini + sg->maxi) / 2;
        cut_boundary(sg->bound, gm->gx, gm->gy, 0, index, sg->minj, sg->maxj, &pl1, &pl2);
    } else {
        /*
         * divide "horizontally" 
         */
        index = (sg->minj + sg->maxj) / 2;
        cut_boundary(sg->bound, gm->gx, gm->gy, 1, index, sg->mini, sg->maxi, &pl1, &pl2);
    }

    if (pl1 == NULL || pl2 == NULL)
        gu_quit("dividesubgrid(): could not cut the boundary\n");

    *sg1 = subgrid_create(gm, pl1, sg->mini, sg->maxi, sg->minj, sg->maxj);
    *sg2 = subgrid_create(gm, pl2, sg->mini, sg->maxi, sg->minj, sg->maxj);
}

static void gridmap_subdivide(gridmap* gm, subgrid* sg)
{
    subgrid* sg1 = NULL;
    subgrid* sg2 = NULL;

    subgrid_divide(sg, &sg1, &sg2);

    if (sg1 != NULL) {
        sg->half1 = sg1;
        ++(gm->nleaves);
        gridmap_subdivide(gm, sg1);
    }
    if (sg2 != NULL) {
        gridmap_subdivide(gm, sg2);
        sg->half2 = sg2;
        ++(gm->nleaves);
    }
    poly_compact(sg->bound, EPS_COMPACT);
}

/* Builds a grid map structure to facilitate conversion from coordinate
 * to index space.
 *
 * @param gx array of X coordinates (of size (nce1+1)*(nce2+1))
 * @param gy array of Y coordinates (of size (nce1+1)*(nce2+1))
 * @param nce1 number of cells in e1 direction
 * @param nce2 number of cells in e2 direction
 * @return a map tree to be use by xy2ij
 */
gridmap* gridmap_build(int nce1, int nce2, double** gx, double** gy)
{
    gridmap* gm = malloc(sizeof(gridmap));
    poly* bound;
    subgrid* trunk;

    gm->nce1 = nce1;
    gm->nce2 = nce2;
    gm->gx = gx;
    gm->gy = gy;
    gm->sign = 0;

    bound = poly_formbound(nce1, nce2, gx, gy);
    trunk = subgrid_create(gm, bound, 0, nce1, 0, nce2);

    gm->bound = bound;
    gm->trunk = trunk;
    gm->nleaves = 1;

    gridmap_subdivide(gm, trunk);       /* recursive */

    return gm;
}

/* Destroys a grid map.
 * @param gm Grid map
 */
void gridmap_destroy(gridmap* gm)
{
    subgrid_destroy(gm->trunk);
    free(gm);
}

/* Calculates indices (i,j) of a grid cell containing point (x,y).
 *
 * @param gm Grid map
 * @param x X coordinate
 * @param y Y coordinate
 * @param i pointer to returned I indice value
 * @param j pointer to returned J indice value
 * @return 1 if successful, 0 otherwhile
 */
int gridmap_xy2ij(gridmap* gm, double x, double y, int* i, int* j)
{
    subgrid* sg = gm->trunk;

    /*
     * check if point is in grid outline 
     */
    sg = gm->trunk;
    if (!poly_containspoint(sg->bound, x, y))
        return 0;

    /*
     * do the full search 
     */
    while (sg->half1 != NULL) {
        /*
         * Test on the point being within the boundary polyline is the most
         * expensive part of the mapping; therefore, perform it in a branch
         * that contains a smaller (number of points-wise) polyline.
         */
        if (sg->half1->bound->n <= sg->half2->bound->n)
            sg = (poly_containspoint(sg->half1->bound, x, y)) ? sg->half1 : sg->half2;
        else
            sg = (poly_containspoint(sg->half2->bound, x, y)) ? sg->half2 : sg->half1;
    }

    *i = sg->mini;
    *j = sg->minj;

    return 1;
}

/* Calculates (x,y) coordinates for a point within the grid with
 * specified fractional indices (i,j).
 *
 * The transformation used to compute the coords is a forward
 * tetragonal bilinear texture mapping.
 *
 * @param gm a tree structure returned from xytoij_init
 * @param fi I indice value
 * @param fj J indice value
 * @param x Pointer to returned X coordinate
 * @param y Pointer to returned Y coordinate
 * @return non-zero if successful
 */
int gridmap_fij2xy(gridmap* gm, double fi, double fj, double* x, double* y)
{
    int status = 1;
    int i, j;
    double u, v;
    double** gx = gm->gx;
    double** gy = gm->gy;
    double a, b, c, d, e, f, g, h;

    /*
     * Trim I to range 0 to nce1 
     */
    if (fi < 0) {
        fi = 0.0;
        status = 0;
    }

    if (fi > gm->nce1) {
        fi = gm->nce1 - EPS;
        status = 0;
    }

    /*
     * Trim J to range 0 to nce2 
     */
    if (fj < 0) {
        fj = 0.0;
        status = 0;
    }

    if (fj > gm->nce2) {
        fj = gm->nce2 - EPS;
        status = 0;
    }

    i = (int) fi;
    j = (int) fj;

    u = fi - i;
    v = fj - j;

    if (u == 0.0 && v == 0.0) {
        *x = gx[j][i];
        *y = gy[j][i];
    } else if (u == 0.0) {
        *x = gx[j + 1][i] * v + gx[j][i] * (1.0 - v);
        *y = gy[j + 1][i] * v + gy[j][i] * (1.0 - v);
    } else if (v == 0.0) {
        *x = gx[j][i + 1] * u + gx[j][i] * (1.0 - u);
        *y = gy[j][i + 1] * u + gy[j][i] * (1.0 - u);
    } else {
        a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
        b = gx[j][i + 1] - gx[j][i];
        c = gx[j + 1][i] - gx[j][i];
        d = gx[j][i];
        e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
        f = gy[j][i + 1] - gy[j][i];
        g = gy[j + 1][i] - gy[j][i];
        h = gy[j][i];

        *x = a * u * v + b * u + c * v + d;
        *y = e * u * v + f * u + g * v + h;
    }

    return status;
}

/* Calculates the branch of sqrt() to be taken in gridmap_xy2fij(). Has to be
 * called only once for a grid.
 * 
 * @param gm Grid map
 * @param x X coordinate
 * @param y Y coordinate
 * @return 1 or -1 if successful; 0 otherwhile
 */
static int calc_branch(gridmap* gm, double x, double y)
{
    int i, j;
    double** gx = gm->gx;
    double** gy = gm->gy;
    int sign = 1;
    double error[2];

    /*
     * normally one tries xytoij() before calling calc_branch() 
     */
    if (gridmap_xy2ij(gm, x, y, &i, &j) == 0)
        return 0;               /* failed */

    {
        double a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
        double b = gx[j][i + 1] - gx[j][i];
        double c = gx[j + 1][i] - gx[j][i];
        double d = gx[j][i];
        double e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
        double f = gy[j][i + 1] - gy[j][i];
        double g = gy[j + 1][i] - gy[j][i];
        double h = gy[j][i];

        double A = a * f - b * e;

        double B, C;
        int k;

        /*
         * normally one checks A before calling calc_branch() 
         */
        if (fabs(A) < EPS_ZERO)
            return 0;           /* failed */

        B = e * x - a * y + a * h - d * e + c * f - b * g;
        C = g * x - c * y + c * h - d * g;

        for (k = 0; k < 2; ++k) {
            double u = (-B + sign * sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
            double v_denom = a * u + c;
            double v = (fabs(v_denom) < EPS_ZERO) ? (y - f * u - h) / (e * u + g) : (x - b * u - d) / v_denom;

            error[k] = 0.0;

            if (u < 0.0)
                error[k] -= u;
            else if (u > 1.0)
                error[k] += (u - 1.0);
            if (v < 0.0)
                error[k] -= v;
            else if (v > 1.0)
                error[k] += (v - 1.0);

            sign = -1;
        }
    }

    if (error[0] < error[1])
        return 1;
    return -1;
}

/* Calculates (x,y) coordinates for a point within a numerical grid specified
 * by fractional indices (i,j).
 *
 * The transformation used to compute the indices is an inverse
 * tetragonal bilinear texture mapping.
 *
 * At the moment we assume that either there is only one grid in the model or
 * all the grids use a uniform branch in xy2fij() convertion.
 *
 * @param gm Grid map
 * @param x X coordinate
 * @param y Y coordinate
 * @param x Pointer to returned fractional I index
 * @param y Pointer to returned fractional J index
 * @return 1 if successful, 0 otherwise
 */
int gridmap_xy2fij(gridmap* gm, double x, double y, double* fi, double* fj)
{
    int i, j;
    double** gx = gm->gx;
    double** gy = gm->gy;

    if (gridmap_xy2ij(gm, x, y, &i, &j) == 0)
        return 0;               /* failed */

    {
        double a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
        double b = gx[j][i + 1] - gx[j][i];
        double c = gx[j + 1][i] - gx[j][i];
        double d = gx[j][i];
        double e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
        double f = gy[j][i + 1] - gy[j][i];
        double g = gy[j + 1][i] - gy[j][i];
        double h = gy[j][i];

        double A = a * f - b * e;
        double B = e * x - a * y + a * h - d * e + c * f - b * g;
        double C = g * x - c * y + c * h - d * g;

        double u, v, d1, d2;

        if (fabs(A) < EPS_ZERO)
            u = -C / B * (1.0 + A * C / B / B);
        else {
            if (gm->sign == 0) {
                gm->sign = calc_branch(gm, x, y);
                if (gm->sign == 0)
                    return 0;   /* failed */
            }
            u = (-B + gm->sign * sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
        }
        d1 = a * u + c;
        d2 = e * u + g;
        v = (fabs(d2) > fabs(d1)) ? (y - f * u - h) / d2 : (x - b * u - d) / d1;

        if (u < 0.0)
            u = 0.0;
        else if (u >= 1.0)
            u = 1.0 - EPS;
        if (v < 0.0)
            v = 0.0;
        else if (v >= 1.0)
            v = 1.0 - EPS;

        *fi = i + u;
        *fj = j + v;
    }

    return 1;
}

int gridmap_getnce1(gridmap* gm)
{
    return gm->nce1;
}

int gridmap_getnce2(gridmap* gm)
{
    return gm->nce2;
}

void gridmap_getextent(gridmap* gm, double* xmin, double* xmax, double* ymin, double* ymax)
{
    extent* e = &gm->trunk->bound->e;

    *xmin = e->xmin;
    *xmax = e->xmax;
    *ymin = e->ymin;
    *ymax = e->ymax;
}
