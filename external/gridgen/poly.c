/******************************************************************************
 *
 * File:           poly.c
 *
 * Created:        early 2002
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *  
 * Purpose:        Library routines for polylines
 *
 * Revisions:      none
 *
 *****************************************************************************/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include "poly.h"

#define POLY_NSTART 4
#define POLY_MAXLINELEN 2048

/* Clears extent.
 * @param e Extent
 */
static void extent_clear(extent* e)
{
    e->xmin = DBL_MAX;
    e->ymin = DBL_MAX;
    e->xmax = -DBL_MAX;
    e->ymax = -DBL_MAX;
}

/* Checks whether a point belongs to extent.
 * @param e Extent
 * @param x X coordinate
 * @param y Y coordinate
 * @return 1 for yes, 0 for no
 */
static int extent_containspoint(extent* e, double x, double y)
{
    if (x >= e->xmin && x <= e->xmax && y >= e->ymin && y <= e->ymax)
        return 1;

    return 0;
}

/* Updates extent to include a specified point.
 * @param e Extent
 * @param x X coordinate
 * @param y Y coordinate
 */
static void extent_update(extent* e, double x, double y)
{
    if (x < e->xmin)
        e->xmin = x;
    if (y < e->ymin)
        e->ymin = y;
    if (x > e->xmax)
        e->xmax = x;
    if (y > e->ymax)
        e->ymax = y;
}

/* Tests whether two points coincide. A threshold distance is used to
 * define the tolerance within which two points may be considered the same.
 *
 * @param p1 first point
 * @param p2 second point
 * @param eps threshold distance
 * @return non-zero if successful
 */
static int point_equal(double x1, double y1, double x2, double y2, double eps)
{
    return ((fabs(x1 - x2) <= eps) && (fabs(y1 - y2) <= eps));
}

/* Re-calculates extent of a polyline.
 * @param pl Polyline
 * @param x X coordinate
 * @param y Y coordinate
 */
static void poly_recalcextent(poly* pl)
{
    double* xs = pl->x;
    double* ys = pl->y;
    extent* e = &pl->e;
    int n = pl->n;
    int i;

    extent_clear(e);
    for (i = 0; i < n; ++i)
        extent_update(e, xs[i], ys[i]);
}

/* Appends a point to the tail of a polyline.
 * @param pl Polyline
 * @param x X coordinate
 * @param y Y coordinate
 */
void poly_addpoint(poly* pl, double x, double y)
{
    if (isnan(x) || isnan(y)) {
        fflush(stdout);
        fprintf(stderr, "error: poly_addpoint(): NaN detected\n");
        exit(1);
    }

    if (pl->n == pl->nallocated) {
        pl->x = realloc(pl->x, pl->nallocated * sizeof(double) * 2);
        pl->y = realloc(pl->y, pl->nallocated * sizeof(double) * 2);
        pl->nallocated *= 2;
    }

    pl->x[pl->n] = x;
    pl->y[pl->n] = y;
    pl->n++;
    extent_update(&pl->e, x, y);
}

/* Adds a point to a polyline at a given position.
 * No action for invalid index.
 * @param pl Polyline
 * @param index Index
 * @param x X coordinate
 * @param y Y coordinate
 */
void poly_addpointat(poly* pl, int index, double x, double y)
{
    if (index > pl->n - 1) {
        poly_addpoint(pl, x, y);
        return;
    }

    if (pl->n == pl->nallocated) {
        pl->x = realloc(pl->x, pl->nallocated * sizeof(double) * 2);
        pl->y = realloc(pl->y, pl->nallocated * sizeof(double) * 2);
        pl->nallocated *= 2;
    }

    memmove(&pl->x[index + 1], &pl->x[index], (pl->n - index) * sizeof(double));
    memmove(&pl->y[index + 1], &pl->x[index], (pl->n - index) * sizeof(double));
    pl->x[index] = x;
    pl->y[index] = y;
    pl->n++;
    extent_update(&pl->e, x, y);
}

/* Appends one polyline to another, like strcat().
 * @param pl1 Destination polyline
 * @param pl2 Source polyline
 */
void poly_append(poly* pl1, poly* pl2)
{
    int n = pl1->n + pl2->n;
    int sizechanged = 0;

    while (n < pl1->nallocated) {
        pl1->nallocated *= 2;
        sizechanged = 1;
    }

    if (sizechanged) {
        pl1->x = realloc(pl1->x, pl1->nallocated * sizeof(double));
        pl1->y = realloc(pl1->y, pl1->nallocated * sizeof(double));
    }

    memcpy(&pl1->x[pl1->n], pl2->x, pl2->n * sizeof(double));
    memcpy(&pl1->y[pl1->n], pl2->y, pl2->n * sizeof(double));
}

/* Computes the area of a polygon. 
 * For an open polyline, implies it being closed. The area is positive for
 * counterclockwise polygon, negative otherwise.
 * @param pl Polyline
 */
double poly_area(poly* pl)
{
    double area = 0.0;
    double* x = pl->x;
    double* y = pl->y;
    int n = pl->n;
    int i;

    for (i = 0; i < n; ++i) {
        int i1 = (i + 1) % n;

        area += (x[i1] - x[i]) * (y[i1] + y[i]);
    }

    return area / 2.0;
}

/* Clears a polyline. 
 * Does not deallocate memory; use poly_destroy() for that.
 * @param pl Polyline
 */
void poly_clear(poly* pl)
{
    pl->n = 0;
    extent_clear(&pl->e);
}

/* Makes a deep copy of a polyline.
 * @param pl Polyline
 * @return Polyline
 */
poly* poly_copy(poly* pl)
{
    poly* pl1 = malloc(sizeof(poly));

    pl1->n = pl->n;
    pl1->nallocated = pl->nallocated;
    pl1->e.xmin = pl->e.xmin;
    pl1->e.xmax = pl->e.xmax;
    pl1->e.ymin = pl->e.ymin;
    pl1->e.ymax = pl->e.ymax;

    pl1->x = malloc(pl1->nallocated * sizeof(double));
    pl1->y = malloc(pl1->nallocated * sizeof(double));
    memcpy(pl1->x, pl->x, pl->n * sizeof(double));
    memcpy(pl1->y, pl->y, pl->n * sizeof(double));

    return pl1;
}

/* Closes a polyline by adding the first point to the tail if necessary.
 * @param pl Polyline
 * @return Polyline
 */
void poly_close(poly* pl)
{
    if (!poly_isclosed(pl, 0.0))
        poly_addpoint(pl, pl->x[0], pl->y[0]);
}

/* Tests whether a point is inside a polygon.
 * The polyline is assumed to be closed: an extra line segment from the end 
 * point back to the start point is assumed if necessary.
 * @param pl Polyline
 * @param x X coordinate
 * @param y Y coordinate
 * @return 1 for yes, 0 for no
 */
int poly_containspoint(poly* pl, double x, double y)
{
    double* xs = pl->x;
    double* ys = pl->y;
    int n = pl->n;
    int hits;
    int i;

    if (n <= 1)
        return 0;
    if (!extent_containspoint(&pl->e, x, y))
        return 0;

    for (i = 0, hits = 0; i < n; ++i) {
        int i1 = (i + 1) % n;
        double x1 = xs[i] - x;
        double y1 = ys[i] - y;
        double x2 = xs[i1] - x;
        double y2 = ys[i1] - y;

        if (y1 == 0.0 && y2 == 0.0) {
            if (x1 * x2 <= 0)
                return 1;
        } else if (y1 == 0.0) {
            if (x1 == 0.0)
                return 1;
            if (x1 > 0.0)
                hits += (y2 > 0.0) ? 1 : -1;
        } else if (y2 == 0.0) {
            if (x2 == 0.0)
                return 1;
            if (x2 > 0.0)
                hits += (y1 < 0.0) ? 1 : -1;
        } else if (y1 * y2 < 0.0) {
            if (x1 > 0.0 && x2 > 0.0)
                hits += 2;
            else if (x1 * x2 <= 0.0) {
                double xx = x1 - (x2 - x1) * y1 / (y2 - y1);

                if (xx == 0)
                    return 1;
                if (xx > 0.0)
                    hits += 2;
            }
        }
    }

    if ((hits / 2) % 2)
        return 1;

    return 0;
}

/* Constructor.
 * @return Polyline
 */
poly* poly_create()
{
    poly* pl = malloc(sizeof(poly));
    pl->x = malloc(POLY_NSTART * sizeof(double));
    pl->y = malloc(POLY_NSTART * sizeof(double));
    pl->n = 0;
    pl->nallocated = POLY_NSTART;
    extent_clear(&pl->e);

    return pl;
}

/* Removes a point at a given position from a polyline.
 * @param pl Polyline
 * @param index Position
 */
void poly_deletepoint(poly* pl, int index)
{
    double xx = pl->x[index];
    double yy = pl->y[index];
    extent* e = &pl->e;

    memmove(&pl->x[index], &pl->x[index + 1], (pl->n - index - 1) * sizeof(double));
    memmove(&pl->y[index], &pl->y[index + 1], (pl->n - index - 1) * sizeof(double));
    pl->n--;

    if (e->xmin == xx || e->xmax == xx || e->ymin == yy || e->ymax == yy) {
        double* x = pl->x;
        double* y = pl->y;
        int i;

        extent_clear(e);
        for (i = 0; i < pl->n; ++i)
            extent_update(e, x[i], y[i]);
    }
}

/* Destructor.
 * @param pl Polyline
 */
void poly_destroy(poly* pl)
{
    free(pl->x);
    free(pl->y);
    free(pl);
}

/* Finds index of a point within polyline.
 * @param pl Polyline
 * @param y Y coordinate
 * @return Index if found; -1 otherwise
 */
int poly_findindex(poly* pl, double x, double y)
{
    double* xs = pl->x;
    double* ys = pl->y;
    int n = pl->n;
    int i;

    for (i = 0; i < n; ++i) {
        if (x == xs[i] && y == ys[i])
            return i;
    }

    return -1;
}

/* Forms polyline boundary around a grid.
 * Note: supposed to handle a grid of corner nodes only.
 * @param nce1 Number of cells in X direction
 * @param nce2 Number of cells in Y direction
 * @param x X coordinates of grid nodes
 * @param y Y coordinates of grid nodes
 * @return Boundary polyline
 */
poly* poly_formbound(int nce1, int nce2, double** x, double** y)
{
    poly* pl = poly_create();
    int direction = 0;          /* 0 - down, 1 - right, 2 - up, 3 -left */
    int iinc[] = { 0, 1, 0, -1 };
    int jinc[] = { 1, 0, -1, 0 };
    int i, j, istart, jstart;

    for (j = 0; j <= nce2; ++j)
        for (i = 0; i <= nce1; ++i)
            if (!isnan(x[j][i]))
                goto ok;

  ok:
    if (j > nce2)
        return pl;

    istart = i;
    jstart = j;
    poly_addpoint(pl, x[j][i], y[j][i]);

    do {
        int direction_stop = (direction + 2) % 4;
        int inext, jnext;

        direction = (direction + 3) % 4;
        inext = i + iinc[direction];
        jnext = j + jinc[direction];

        while ((inext < 0 || jnext < 0 || inext > nce1 || jnext > nce2 || isnan(x[jnext][inext])) && direction != direction_stop) {
            direction = (direction + 1) % 4;
            inext = i + iinc[direction];
            jnext = j + jinc[direction];
        }

        if (direction == direction_stop)
            break;

        i = inext;
        j = jnext;

        poly_addpoint(pl, x[j][i], y[j][i]);

    } while (j != jstart || i != istart);

    poly_close(pl);

    return (pl);
}

/* Forms polyline boundary around a grid in index space.
 * @param nce1 Number of cells in X direction
 * @param nce2 Number of cells in Y direction
 * @param x X coordinates of grid nodes
 * @return Boundary polyline
 */
poly* poly_formboundij(int nce1, int nce2, double** x)
{
    poly* pl = poly_create();
    int direction = 0;          /* 0 - down, 1 - right, 2 - up, 3 -left */
    int iinc[] = { 0, 1, 0, -1 };
    int jinc[] = { 1, 0, -1, 0 };
    int i, j, istart, jstart;

    for (j = 0; j <= nce2; ++j)
        for (i = 0; i <= nce1; ++i)
            if (!isnan(x[j][i]))
                goto ok;

  ok:
    if (j > nce2)
        return pl;

    istart = i;
    jstart = j;
    poly_addpoint(pl, i, j);

    do {
        int direction_stop = (direction + 2) % 4;
        int inext, jnext;

        direction = (direction + 3) % 4;
        inext = i + iinc[direction];
        jnext = j + jinc[direction];

        while ((inext < 0 || jnext < 0 || inext > nce1 || jnext > nce2 || isnan(x[jnext][inext])) && direction != direction_stop) {
            direction = (direction + 1) % 4;
            inext = i + iinc[direction];
            jnext = j + jinc[direction];
        }

        if (direction == direction_stop)
            break;

        i = inext;
        j = jnext;

        poly_addpoint(pl, i, j);

    } while (j != jstart || i != istart);

    poly_close(pl);

    return (pl);
}

/* Checks whether the polyline is closed.
 * @param pl Polyline
 * @param eps Distance tolerance
 * @return 1 for yes, 0 for no
 */
int poly_isclosed(poly* pl, double eps)
{
    return pl->n > 1 && point_equal(pl->x[0], pl->y[0], pl->x[pl->n - 1], pl->y[pl->n - 1], eps);
}

static int iscommentline(char* line)
{
    return (line[0] == '#');
}

static int isblankline(char* line)
{
    while (*line && isspace(*line))
        line++;
    if (*line)
        return 0;
    return 1;
}

static int nextline(char* line, int n, FILE* fp)
{
    char* s;

    do {
        s = fgets(line, n, fp);
    } while (s && (iscommentline(line) || isblankline(line)));

    if (s == NULL)
        return 0;

    return 1;
}

/* Reads points from a file stream and appends them to a polyline.
 * @param pl Polyline
 * @param fp File handle
 * @return The number of points in the polyline
 */
int poly_read(poly* pl, FILE* fp)
{
    double x, y;
    char buf[POLY_MAXLINELEN];

    /*
     * skip comments and blank lines 
     */
    if (nextline(buf, POLY_MAXLINELEN, fp) == 0)
        return 0;

    while (sscanf(buf, "%lf %lf", &x, &y) == 2) {
        poly_addpoint(pl, x, y);
        if (fgets(buf, POLY_MAXLINELEN, fp) == NULL)
            break;
    }

    return (pl->n);
}

/* Resamples a polyline including only points more than the threshold
 ** distance apart.
 * @param pl Polyline
 * @param eps Threshold distance
 */
void poly_resample(poly* pl, double eps)
{
    double* xs = pl->x;
    double* ys = pl->y;
    int n = pl->n;
    int i, nn;

    if (pl->n <= 1)
        return;

    for (i = 1, nn = 0; i < n; ++i) {
        if (!point_equal(xs[nn], ys[nn], xs[i], ys[i], eps)) {
            nn++;
            xs[nn] = xs[i];
            ys[nn] = ys[i];
        }
    }

    if (pl->n != nn + 1) {
        pl->n = nn + 1;
        poly_recalcextent(pl);
    }
}

/* Reverses points in a polyline.
 * @param pl Polyline
 */
void poly_reverse(poly* pl)
{
    double* x = pl->x;
    double* y = pl->y;
    int n = pl->n;
    int n2 = n / 2;
    int i, j;

    if (pl->n <= 1)
        return;

    for (i = 0, j = pl->n - 1; i < n2; ++i, --j) {
        double tmp = x[i];

        x[i] = x[j];
        x[j] = tmp;
        tmp = y[i];
        y[i] = y[j];
        y[j] = tmp;
    }
}

/* Removes spikes from a polyline.
 * @param pl Polyline
 * @param maxdist Maximal allowed distance between two adjacent points
 */
void poly_despike(poly* pl, double maxdist)
{
    double* xs = pl->x;
    double* ys = pl->y;
    int n = pl->n;
    int i, nn;

    if (pl->n <= 1)
        return;

    for (i = 1, nn = 0; i < n; ++i) {
        if (hypot(xs[nn] - xs[i], ys[nn] - ys[i]) < maxdist) {
            nn++;
            xs[nn] = xs[i];
            ys[nn] = ys[i];
        }
    }

    if (pl->n != nn + 1) {
        pl->n = nn + 1;
        poly_recalcextent(pl);
    }
}

/* Writes polyline to a file stream.
 * @param pl Polyline
 * @param fp File handle
 */
void poly_write(poly* pl, FILE* fp)
{
    int i;

    fprintf(fp, "## %d\n", pl->n);
    for (i = 0; i < pl->n; ++i)
        fprintf(fp, "%.15g %.15g\n", pl->x[i], pl->y[i]);
}

static int ononeline(double* x, double* y, int i1, int i2, int i3, double eps)
{
    return fabs((x[i1] - x[i2]) * (y[i3] - y[i2]) - (x[i3] - x[i2]) * (y[i1] - y[i2])) <= eps;
}

/* Deletes redundant nodes.
 * @param pl Polyline
 * @param eps A small number used in tests on two points being the same or
 *            three points belonging to one line.
 */
void poly_compact(poly* pl, double eps)
{
    int n = pl->n;
    double* x = pl->x;
    double* y = pl->y;
    int* ids = NULL;
    int nnew;
    int ileft, imiddle, iright;
    int i;

    if (n <= 4)
        return;                 /* do not bother */

    ids = malloc(pl->n * sizeof(int));
    nnew = 0;

    for (i = 0, imiddle = 0, ileft = n - 1; i < n - 1; ++i)
        if (!point_equal(x[ileft], y[ileft], x[imiddle], y[imiddle], eps))
            break;
        else
            ileft--;

    for (i = 0, iright = 1; i < n; ++i) {
        if (!point_equal(x[iright], y[iright], x[imiddle], y[imiddle], eps)) {
            if (!ononeline(x, y, ileft, imiddle, iright, eps)) {
                ids[nnew++] = imiddle;
                ileft = imiddle;
            }
            imiddle = iright;
        }
        iright = (iright + 1) % n;
    }

    if (nnew != n) {
        for (i = 0; i < nnew; ++i) {
            x[i] = x[ids[i]];
            y[i] = y[ids[i]];
        }
        pl->n = nnew;
    }

    free(ids);
    pl->x = realloc(pl->x, sizeof(double) * pl->n);
    pl->y = realloc(pl->y, sizeof(double) * pl->n);
    pl->nallocated = pl->n;
}
