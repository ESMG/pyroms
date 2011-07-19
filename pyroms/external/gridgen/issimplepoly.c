/******************************************************************************
 *
 * File:        issimplepoly.c
 *
 * Created:     24/02/2003
 *
 * Author:      Pavel Sakov
 *              CSIRO Marine Research
 *
 * Description: Testing a polygon on self-intersections by using Shamos-Hoey 
 *              ("sweep-line") algorithm.
 *
 *              Beware that a polygon with duplicated vertices is considered
 *              self-intersecting.
 *
 *              The original code was borrowed from 
 *              http://geometryalgorithms.com/Archive/algorithm_0108/
 *              /algorithm_0108.htm.
 *              It contained 
 *              "Copyright 2001, softSurfer (www.softsurfer.com)".
 *
 *              The current implementation is pretty far away from this code.
 *
 * Revisions:   PS 180106 Corrected slseg_compare() -- thanks to
 *              Dr. Juergen Zangers for a number of torture tests.
 *
 *              PS 180106 Got rid of using the avltree structure. The
 *              code became cleaner and supposedly more efficient.
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <limits.h>
#include "issimplepoly.h"

static int verbose = 0;
static int issimplepoly_closed = 1;

typedef enum { FIRST = 0, SECOND = 1 } ORDER;

typedef struct {
    int eid;                    /* polygon edge id */
    ORDER order;
    double* x;                  /* X coordinate of the endpoint */
    double* y;                  /* Y coordinate of the endpoint */
} endpoint;

static int point_compare(double* x1, double* y1, double* x2, double* y2)
{
    if (*x1 > *x2)
        return 1;
    if (*x1 < *x2)
        return -1;
    if (*y1 > *y2)
        return 1;
    if (*y1 < *y2)
        return -1;
    return 0;
}

static int endpoint_compare(const void* p1, const void* p2)
{
    endpoint* e1 = *(endpoint**) p1;
    endpoint* e2 = *(endpoint**) p2;

    return point_compare(e1->x, e1->y, e2->x, e2->y);
}

typedef struct {
    int n;                      /* total number of endpoints */
    int next;                   /* index of next endpoint on queue */
    endpoint* data;             /* array of all endpoints */
    endpoint** sorted;          /* sorted list of endpoint pointers */
} endpointqueue;

static endpointqueue* eq_create(int n, double x[], double y[])
{
    endpointqueue* eq = malloc(sizeof(endpointqueue));
    int nsegs = n;
    int i;

    if (!issimplepoly_closed)
        nsegs = n - 1;

    eq->n = 2 * nsegs;          /* two vertex endpoints for each edge */
    eq->next = 0;

    eq->data = malloc(eq->n * sizeof(endpoint));
    for (i = 0; i < nsegs; ++i) {
        endpoint* e1 = &eq->data[i * 2];
        endpoint* e2 = &eq->data[i * 2 + 1];

        e1->eid = i;
        e2->eid = i;
        e1->x = &x[i];
        e1->y = &y[i];
        e2->x = &x[(i + 1) % n];
        e2->y = &y[(i + 1) % n];

        if (point_compare(e1->x, e1->y, e2->x, e2->y) > 0) {
            e1->order = SECOND;
            e2->order = FIRST;
        } else {
            e1->order = FIRST;
            e2->order = SECOND;
        }
    }

    eq->sorted = malloc(eq->n * sizeof(endpoint*));
    for (i = 0; i < eq->n; ++i)
        eq->sorted[i] = &eq->data[i];

    qsort(eq->sorted, eq->n, sizeof(endpoint*), endpoint_compare);

    return eq;
}

static void eq_destroy(endpointqueue * eq)
{
    free(eq->data);
    free(eq->sorted);
    free(eq);
}

static endpoint* eq_getnext(endpointqueue * eq)
{
    if (eq->next >= eq->n)
        return NULL;
    else
        return eq->sorted[eq->next++];
}

struct slseg;
typedef struct slseg slseg;

struct sweepline;
typedef struct sweepline sweepline;

struct slseg {
    int eid;                    /* polygon edge id */
    double lxy[2];              /* left vertex */
    double rxy[2];              /* right vertex */
    sweepline* sl;
};

struct segnode;
typedef struct segnode segnode;

struct segnode {
    slseg* s;
    segnode* next;
};

struct sweepline {
    int n;                      /* number of vertices in polygon */
    double* x;                  /* x coordinates [n] */
    double* y;                  /* y coordinates [n] */
    segnode* firstnode;
    double pos;                 /* position */
};

static slseg* slseg_create(sweepline* sl, endpoint* e)
{
    int n = sl->n;
    double* x = sl->x;
    double* y = sl->y;
    slseg* s = malloc(sizeof(slseg));
    int i = e->eid;
    int i1 = (e->eid + 1) % n;

    s->eid = i;
    if (point_compare(&x[i], &y[i], &x[i1], &y[i1]) < 0) {
        s->lxy[0] = x[i];
        s->lxy[1] = y[i];
        s->rxy[0] = x[i1];
        s->rxy[1] = y[i1];
    } else {
        s->rxy[0] = x[i];
        s->rxy[1] = y[i];
        s->lxy[0] = x[i1];
        s->lxy[1] = y[i1];
    }

    s->sl = sl;

    return s;
}

static void slseg_destroy(slseg* s)
{
    free(s);
}

static sweepline* sl_create(int n, double x[], double y[])
{
    sweepline* sl = malloc(sizeof(sweepline));

    sl->n = n;
    sl->x = x;
    sl->y = y;
    sl->firstnode = NULL;
    sl->pos = -DBL_MAX;

    return sl;
}

static void sl_destroy(sweepline* sl)
{
    segnode* sn = sl->firstnode;

    while (sn != NULL) {
        segnode* next = sn->next;

        slseg_destroy(sn->s);
        free(sn);
        sn = next;
    }

    free(sl);
}

static int slseg_compare(void* p1, void* p2)
{
    slseg* s1 = (slseg*) p1;
    slseg* s2 = (slseg*) p2;
    double x = s1->sl->pos;
    double* xy11 = s1->lxy;
    double* xy12 = s1->rxy;
    double* xy21 = s2->lxy;
    double* xy22 = s2->rxy;
    double y1, y2;

    if (xy11[0] == x) {
        if (xy12[0] == x)
            y1 = (xy11[1] + xy12[1]) / 2.0;
        else
            y1 = xy11[1];
    } else if (xy12[0] == x)
        y1 = xy12[1];
    else
        y1 = (xy11[1] * (xy12[0] - x) + xy12[1] * (x - xy11[0])) / (xy12[0] - xy11[0]);

    if (xy21[0] == x) {
        if (xy22[0] == x)
            y2 = (xy21[1] + xy22[1]) / 2.0;
        else
            y2 = xy21[1];
    } else if (xy22[0] == x)
        y2 = xy22[1];
    else
        y2 = (xy21[1] * (xy22[0] - x) + xy22[1] * (x - xy21[0])) / (xy22[0] - xy21[0]);

    if (y1 > y2)
        return 1;
    if (y1 < y2)
        return -1;

    if (s1->eid > s2->eid)
        return 1;
    if (s1->eid < s2->eid)
        return -1;

    return 0;
}

static slseg* sl_add(sweepline* sl, endpoint* e, slseg** above, slseg** below)
{
    segnode* sn = sl->firstnode;
    segnode* last = NULL;
    segnode* new = NULL;
    slseg* s = slseg_create(sl, e);

    sl->pos = *e->x;
    *above = NULL;
    *below = NULL;
    while (sn != NULL) {
        slseg* ss = sn->s;
        int res = slseg_compare(ss, s);

        assert(res != 0);
        if (res > 0) {
            if (*above == NULL || slseg_compare(ss, *above) < 0)
                *above = ss;
        } else {
            if (*below == NULL || slseg_compare(ss, *below) > 0)
                *below = ss;
        }
        last = sn;
        sn = sn->next;
    }

    new = malloc(sizeof(segnode));
    new->s = s;
    new->next = NULL;
    if (sl->firstnode == NULL)
        sl->firstnode = new;
    if (last != NULL)
        last->next = new;

    return s;
}

static void sl_remove(sweepline* sl, endpoint* e, slseg** above, slseg** below)
{
    slseg* s = slseg_create(sl, e);
    segnode* sn = sl->firstnode;
    segnode* prev = NULL;
    int removed = 0;

    sl->pos = *e->x;
    *above = NULL;
    *below = NULL;
    while (sn != NULL) {
        slseg* ss = sn->s;
        segnode* next = sn->next;
        int res = slseg_compare(ss, s);

        if (res > 0) {
            if (*above == NULL || slseg_compare(ss, *above) < 0)
                *above = ss;
        } else if (res < 0) {
            if (*below == NULL || slseg_compare(ss, *below) > 0)
                *below = ss;
        } else {
            if (prev == NULL)
                sl->firstnode = sn->next;
            else
                prev->next = sn->next;
            slseg_destroy(sn->s);
            free(sn);

            removed = 1;
        }
        prev = sn;
        sn = next;
    }
    assert(removed);
    slseg_destroy(s);
}

static double vectmult(double* xy1, double* xy2, double* xy3)
{
    return (xy2[0] - xy1[0]) * (xy3[1] - xy1[1]) - (xy3[0] - xy1[0]) * (xy2[1] - xy1[1]);
}

static int _sl_intersect(sweepline* sl, slseg* s1, slseg* s2)
{
    int adjacent = 0;
    int i1 = -1, i2 = -1, i3 = -1;
    double lsign, rsign;

    if (s1 == NULL || s2 == NULL)
        return 0;

    if ((s1->eid + 1) % sl->n == s2->eid) {
        adjacent = 1;
        i1 = s1->eid;
        i2 = s2->eid;
        i3 = (i2 + 1) % sl->n;
    } else if ((s2->eid + 1) % sl->n == s1->eid) {
        adjacent = 1;
        i1 = s2->eid;
        i2 = s1->eid;
        i3 = (i2 + 1) % sl->n;
    }
    if (adjacent) {
        double xy1[2], xy2[2], xy3[2];

        xy1[0] = 0.0;
        xy1[1] = 0.0;
        xy2[0] = sl->x[i2] - sl->x[i1];
        xy2[1] = sl->y[i2] - sl->y[i1];
        xy3[0] = sl->x[i3] - sl->x[i1];
        xy3[1] = sl->y[i3] - sl->y[i1];

        if (fabs(vectmult(xy1, xy2, xy3)) != 0.0)
            return 0;
        if (xy1[0] != xy2[0])
            return (xy3[0] - xy2[0]) * (xy2[0] - xy1[0]) < 0.0;
        if (xy1[1] != xy2[1])
            return (xy3[1] - xy2[1]) * (xy2[1] - xy1[1]) < 0.0;
        return 1;
    }

    lsign = vectmult(s1->lxy, s1->rxy, s2->lxy);
    rsign = vectmult(s1->lxy, s1->rxy, s2->rxy);

    if (lsign * rsign > 0.0)
        return 0;

    lsign = vectmult(s2->lxy, s2->rxy, s1->lxy);
    rsign = vectmult(s2->lxy, s2->rxy, s1->rxy);

    if (lsign * rsign > 0.0)
        return 0;

    return 1;
}

static int sl_intersect(sweepline* sl, slseg* s1, slseg* s2)
{
    if (_sl_intersect(sl, s1, s2)) {
        if (verbose)
            printf("  edge %d [(%.7g,%.7g),(%.7g,%.7g)] intersects with edge %d [(%.7g,%.7g),(%.7g,%.7g)]\n", s1->eid, s1->lxy[0], s1->lxy[1], s1->rxy[0], s1->rxy[1], s2->eid, s2->lxy[0], s2->lxy[1], s2->rxy[0], s2->rxy[1]);
        return 1;
    }
    return 0;
}

static void print_endpoint(endpoint* e, int n, double x[], double y[], slseg* above, slseg* below)
{
    int static count = 0;

    printf("    %d: (%.7g,%.7g): %s of [(%.7g,%.7g),(%.7g,%.7g)]", count, *e->x, *e->y, (e->order == FIRST) ? "left endpoint" : "right endpoint", x[e->eid], y[e->eid], x[(e->eid + 1) % n], y[(e->eid + 1) % n]);
    if (above == NULL)
        printf(", above = NULL");
    else
        printf(", above = [(%.7g,%.7g),(%.7g,%.7g)]", above->lxy[0], above->lxy[1], above->rxy[0], above->rxy[1]);
    if (below == NULL)
        printf(", below = NULL");
    else
        printf(", below = [(%.7g,%.7g),(%.7g,%.7g)]", below->lxy[0], below->lxy[1], below->rxy[0], below->rxy[1]);
    printf("\n");
    count++;
}

int issimplepolygon(int n, double x[], double y[])
{
    endpointqueue* eq = eq_create(n, x, y);
    sweepline* sl = sl_create(n, x, y);
    endpoint* e;                /* the current endpoint */
    slseg* s;                   /* the current segment */
    slseg* above;
    slseg* below;
    int result = 1;

    while ((e = eq_getnext(eq)) != NULL) {
        if (e->order == FIRST) {
            s = sl_add(sl, e, &above, &below);
            if (verbose > 1)
                print_endpoint(e, n, x, y, above, below);
            if (sl_intersect(sl, s, above) || sl_intersect(sl, s, below)) {
                result = 0;
                break;
            }
        } else {
            sl_remove(sl, e, &above, &below);
            if (verbose > 1)
                print_endpoint(e, n, x, y, above, below);
            if (sl_intersect(sl, above, below)) {
                result = 0;
                break;
            }
        }
    }

    eq_destroy(eq);
    sl_destroy(sl);

    return result;
}

#if defined(TEST_SI)

#define BUFSIZE 1024
#define NALLOC_START 100

#include <stdarg.h>
#include <errno.h>

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

static void usage(char* programname)
{
    printf("Usage: %s [-o] [-v|-V] <polyline file>\n", programname);
    printf("Options:\n");
    printf("  -o -- consider the polyline open\n");
    printf("  -v -- verbose\n");
    printf("  -V -- more verbose\n");
    exit(0);
}

static void xy_read(FILE* f, int* n, double** x, double** y)
{
    char buf[BUFSIZE];
    int nallocated = 0;

    *n = 0;
    *x = NULL;
    *y = NULL;

    while (fgets(buf, BUFSIZE, f) != NULL) {
        char seps[] = " ,;\t";
        double xx, yy;
        char* token;

        if (buf[0] == '#')
            continue;

        if ((token = strtok(buf, seps)) == NULL)
            continue;
        xx = atof(token);

        if ((token = strtok(NULL, seps)) == NULL)
            continue;
        yy = atof(token);

        if (*n == nallocated) {
            nallocated = (nallocated == 0) ? NALLOC_START : nallocated * 2;
            *x = realloc(*x, nallocated * sizeof(double));
            *y = realloc(*y, nallocated * sizeof(double));
        }

        (*x)[*n] = xx;
        (*y)[*n] = yy;
        (*n)++;
    }

    if (*n > 0)
        if ((*x)[*n - 1] == (*x)[0] && (*y)[*n - 1] == (*y)[0])
            (*n)--;
}

static void parse_commandline(int argc, char* argv[], char** fname)
{
    int i;

    *fname = NULL;

    i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') {
            if (*fname == NULL) {
                *fname = argv[i];
                i++;
            } else
                usage(argv[0]);
        } else {
            switch (argv[i][1]) {
            case 'o':
                issimplepoly_closed = 0;
                i++;
                break;
            case 'v':
                verbose = 1;
                i++;
                break;
            case 'V':
                verbose = 2;
                i++;
                break;
            default:
                usage(argv[0]);
                break;
            }
        }
    }
}

int main(int argc, char* argv[])
{
    char* fname = argv[1];
    FILE* f = NULL;
    int n = 0;
    double* x = NULL;
    double* y = NULL;

    parse_commandline(argc, argv, &fname);

    if (fname == NULL) {
        printf("error: no input data\n");
        usage(argv[0]);
    }

    f = fopen(fname, "r");
    if (f == NULL)
        quit("could not open \"%s\" for read: %s\n", fname, strerror(errno));
    xy_read(f, &n, &x, &y);
    if (verbose > 1)
        printf("  %d points\n", n);
    fclose(f);

    if (!issimplepolygon(n, x, y))
        printf("  is NOT a simple polygon\n");
    else
        printf("  IS a simple polygon\n");

    if (x != NULL)
        free(x);
    if (y != NULL)
        free(y);

    return 0;
}

#endif                          /* TEST_SI */
