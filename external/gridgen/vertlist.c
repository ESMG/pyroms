/******************************************************************************
 *
 * File:           vertlist.c
 *
 * Created:        19/03/2001
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Handling list of polygon vertices
 *
 * Description:    None
 *
 * Revisions:      None
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "config.h"
#include "delaunay.h"
#include "istack.h"
#include "geom.h"
#include "vertlist.h"
#include "nan.h"

#define SQRT2x3 (M_SQRT2 * 3.0)

void vertlist_init(vertlist* l)
{
    l->first = NULL;
    l->n = 0;
}

vertlist* vertlist_create()
{
    vertlist* l = malloc(sizeof(vertlist));

    vertlist_init(l);

    return l;
}

void vertlist_add(vertlist* l, point* p)
{
    vertnode* new = malloc(sizeof(vertnode));

    new->p.x = p->x;
    new->p.y = p->y;
    new->p.z = p->z;
    new->protected = 0;

    if (l->n == 0) {
        l->first = new;
        new->prev = new;
    } else {
        vertnode* last = l->first->prev;

        last->next = new;
        l->first->prev = new;
        new->prev = last;
    }
    new->next = l->first;
    l->n++;
}

vertlist* vertlist_create2(int n, double x[], double y[], double z[])
{
    vertlist* l = vertlist_create();
    int i;

    for (i = 0; i < n; ++i) {
        point p = { x[i], y[i], z[i] };

        vertlist_add(l, &p);
    }

    return l;
}

void vertlist_delete(vertlist* l, vertnode* v)
{
    vertnode* prev = v->prev;
    vertnode* next = v->next;

    if (l->first == v)
        l->first = v->next;

    free(v);

    if (prev != next) {
        prev->next = next;
        next->prev = prev;
        l->n--;
    } else
        vertlist_init(l);
}

void vertlist_insert_after(vertlist* l, vertnode* v, point* p)
{
    vertnode* new = malloc(sizeof(vertnode));

    new->p.x = p->x;
    new->p.y = p->y;
    new->p.z = p->z;
    new->protected = 0;
    new->prev = v;
    new->next = v->next;
    v->next = new;
    new->next->prev = new;
    l->n++;
}

void vertlist_clear(vertlist* l)
{
    vertnode* first = l->first;
    vertnode* now = l->first;

    do {
        vertnode* next = now->next;

        free(now);
        now = next;
    } while (now != first);

    vertlist_init(l);
}

void vertlist_destroy(vertlist* l)
{
    if (l == NULL)
        return;
    vertlist_clear(l);
    free(l);
}

void vertlist_find_minmax(vertlist* l, double* xmin, double* xmax, double* ymin, double* ymax)
{
    vertnode* first = l->first;
    vertnode* now = l->first;

    *xmin = DBL_MAX;
    *xmax = -DBL_MAX;
    *ymin = DBL_MAX;
    *ymax = -DBL_MAX;

    do {
        point* p = &now->p;

        if (p->x < *xmin)
            *xmin = p->x;
        if (p->x > *xmax)
            *xmax = p->x;
        if (p->y < *ymin)
            *ymin = p->y;
        if (p->y > *ymax)
            *ymax = p->y;
        now = now->next;
    } while (now != first);
}

void vertlist_thin(vertlist* l, double dxmin, double dymin)
{
    vertnode* now = l->first;
    int n = l->n;
    int i = 0;

    if (now == NULL)
        return;

    do {
        vertnode* next = now->next;
        point* p = &now->p;
        point* p1 = &next->p;

        if (fabs(p->x - p1->x) < dxmin && fabs(p->y - p1->y) < dymin) {
            if (p1->z != 0)
                p->z = p1->z;
            vertlist_delete(l, next);
            if (gg_verbose > 1)
                fprintf(stderr, "  vertex %d deleted\n", i);
        } else
            now = next;

        i++;
    } while (i < n);
}

void vertlist_print(vertlist* l, FILE* f)
{
    int i = 0;
    vertnode* first = l->first;
    vertnode* now = l->first;

    fprintf(stderr, "\n");

    if (now == NULL) {
        fprintf(stderr, "  <no vertices>\n\n");
        return;
    }

    do {
        point* p = &now->p;

        fprintf(f, "  %d:  (%.8g,%.8g)\n", i, p->x, p->y);
        ++i;
        now = now->next;
    } while (now != first);
    fprintf(f, "\n");

    fflush(f);
}

void vertlist_read(vertlist* l, FILE* f)
{
    char buf[BUFSIZE];
    int firstnode_index = 0;

    while (fgets(buf, BUFSIZE, f) != NULL) {
        char seps[] = " ,;\t";
        char* token;
        point p;

        if (buf[0] == '#')
            continue;

        if ((token = strtok(buf, seps)) == NULL)
            continue;
        p.x = atof(token);

        if ((token = strtok(NULL, seps)) == NULL)
            continue;
        p.y = atof(token);

        if ((token = strtok(NULL, seps)) == NULL)
            p.z = 0;
        else {
            p.z = atof(token);
            if (strchr(token, '*') != NULL)
                firstnode_index = l->n;
        }

        vertlist_add(l, &p);
    }

    vertlist_setfirstnode(l, firstnode_index);
}

point* vertlist_topoint(vertlist* l)
{
    point* out;
    vertnode* now = l->first;
    int n = l->n;
    int i;

    if (n == 0)
        return NULL;

    out = malloc(n * sizeof(point));

    for (i = 0; i < n; ++i) {
        point* pin = &now->p;
        point* pout = &out[i];

        pout->x = pin->x;
        pout->y = pin->y;
        pout->z = pin->z;
        now = now->next;
    }

    return out;
}

double vertlist_area(vertlist* l)
{
    vertnode* now = l->first;
    int n = l->n;
    double area = 0.0;
    int i;

    if (n < 3)
        return 0.0;

    for (i = 0; i < n; ++i) {
        point* p = &now->p;
        point* p1 = &now->next->p;

        area += (p->x - p1->x) * (p1->y + p->y);

        now = now->next;
    }

    return area / 2.0;
}

void vertlist_change_order(vertlist** l)
{
    vertlist* lnew = vertlist_create();
    vertnode* now = (*l)->first;
    vertnode* first = (*l)->first;

    do {
        point* p = &now->p;

        vertlist_add(lnew, p);
        now = now->prev;
    } while (now != first);

    vertlist_destroy(*l);
    (*l) = lnew;
}

static zdouble p2z(point* p)
{
    return p->x + I * p->y;
}

zdouble* vertlist_tozdouble(vertlist* l)
{
    int n = l->n;
    zdouble* zs = malloc(n * sizeof(zdouble));
    vertnode* now = l->first;
    int i;

    for (i = 0; i < n; ++i) {
        zs[i] = p2z(&now->p);
        now = now->next;
    }

    return zs;
}

void vertlist_toxy(vertlist* l, double** x, double** y)
{
    int n = l->n;
    vertnode* now = l->first;
    int i;

    *x = malloc(n * 2 * sizeof(double));
    *y = malloc(n * 2 * sizeof(double));
    for (i = 0; i < n; ++i) {
        point* p = &now->p;

        (*x)[i] = p->x;
        (*y)[i] = p->y;

        now = now->next;
    }
}

static double vertnode_angle_cos(vertnode* v)
{
    point* p = &v->p;
    point* p1 = &v->next->p;
    point* p2 = &v->prev->p;
    double dx1 = p1->x - p->x;
    double dx2 = p2->x - p->x;
    double dy1 = p1->y - p->y;
    double dy2 = p2->y - p->y;

    return (dx1 * dx2 + dy1 * dy2) / sqrt((dx1 * dx1 + dy1 * dy1) * (dx2 * dx2 + dy2 * dy2));
}

void vertlist_process_phase1(vertlist* l)
{
    vertnode* first = l->first;
    vertnode* now = l->first;
    int n = 0;

    if (gg_verbose) {
        fprintf(stderr, "phase 1 split:\n");
        fflush(stderr);
    }

    do {
        if (vertnode_angle_cos(now) > M_SQRT1_2) {
            point pnew;
            point* p = &now->p;
            point* p1 = &now->next->p;
            point* p2 = &now->prev->p;
            double dx1 = p1->x - p->x;
            double dx2 = p2->x - p->x;
            double dy1 = p1->y - p->y;
            double dy2 = p2->y - p->y;
            double l0, l1, l2;
            double k1, k2;

            if (gg_verbose) {
                if (n == 0)
                    fprintf(stderr, "  ");
                fprintf(stderr, ".");
                fflush(stderr);
            }

            l1 = point_point_distance(p, p1);
            l2 = point_point_distance(p, p2);
            l0 = (l1 < l2) ? l1 : l2;
            k1 = l0 / l1 / 2.0;
            k2 = l0 / l2 / 2.0;

            pnew.x = p->x + dx2 * k2;
            pnew.y = p->y + dy2 * k2;
            pnew.z = NaN;

            vertlist_insert_after(l, now->prev, &pnew);

            pnew.x = p->x + dx1 * k1;
            pnew.y = p->y + dy1 * k1;
            pnew.z = NaN;

            vertlist_insert_after(l, now, &pnew);

            now->protected = 1;
            n++;
        }
        now = now->next;
    } while (now != first);

    if (gg_verbose)
        fprintf(stderr, "  %d vertices after phase 1 split\n", l->n);
    if (n > 0 && gg_verbose > 1)
        vertlist_print(l, stderr);
}

delaunay* vertlist_triangulate(vertlist* l, FILE* log)
{
    delaunay* d;
    point* vertices = vertlist_topoint(l);
    int n = l->n;
    int* edges = malloc(n * 2 * sizeof(int));
    int i, j;

    if (gg_verbose && log != NULL) {
        fprintf(log, "triangulating:\n");
        fflush(log);
    }

    for (i = 0, j = 0; i < n; ++i) {
        edges[j] = i;
        j++;
        edges[j] = (i + 1) % n;
        j++;
    }

    d = delaunay_build(n, vertices, n, edges, 0, NULL);

    free(edges);
    free(vertices);

    if (gg_verbose && log != NULL) {
        fprintf(log, "  %d triangles\n", d->ntriangles);
        fprintf(log, "  %d edges\n", d->nedges);
    }

    return d;
}

void vertlist_process_phase2(vertlist* l)
{
    static int times_called = 0;
    vertnode* first = l->first;
    vertnode* now = l->first;
    delaunay* d = vertlist_triangulate(l, NULL);
    int n = l->n;
    point* vertices = d->points;
    int i = 0;
    istack* tocheck = istack_create();

    if (gg_verbose) {
        if (times_called == 0)
            fprintf(stderr, "phase 2 split:\n  ");
        if (gg_verbose == 1)
            fprintf(stderr, ".");
        fflush(stderr);
    }

    times_called++;

    do {
        if (!now->protected && !now->next->protected) {
            int i1 = (i + 1) % n;
            point* p = &now->p;
            point* p1 = &now->next->p;
            double edgelength = point_point_distance(p, p1);
            double mindist = DBL_MAX;
            int ii, count;

            istack_reset(tocheck);

            for (ii = i, count = 0; count < 2; ii = (ii + 1) % n, ++count) {
                int t;

                for (t = 0; t < d->n_point_triangles[i]; ++t) {
                    int tid = d->point_triangles[i][t];
                    int v;

                    for (v = 0; v < 3; ++v) {
                        int vertex = d->triangles[tid].vids[v];

                        if (vertex != i && vertex != i1 && !istack_contains(tocheck, vertex))
                            istack_push(tocheck, vertex);
                    }
                }
            }

            while (tocheck->n > 0) {
                int vertex = istack_pop(tocheck);
                double dist = point_edge_distance(&vertices[vertex], p, p1);

                if (dist < mindist)
                    mindist = dist;
            }

            if (mindist < edgelength / SQRT2x3) {
                double dx = p1->x - p->x;
                double dy = p1->y - p->y;
                point pnew;

                pnew.x = p->x + dx / 3.0;
                pnew.y = p->y + dy / 3.0;
                pnew.z = NaN;

                vertlist_insert_after(l, now, &pnew);

                pnew.x = p->x + dx * 2.0 / 3.0;
                pnew.y = p->y + dy * 2.0 / 3.0;

                vertlist_insert_after(l, now->next, &pnew);

                vertlist_process_phase2(l);

                break;
            }
        }
        ++i;
        now = now->next;
    } while (now != first);

    times_called--;

    if (gg_verbose && times_called == 0)
        fprintf(stderr, "\n  %d vertices after phase 2 split\n", l->n);

    delaunay_destroy(d);
    istack_destroy(tocheck);
}

void vertlist_setfirstnode(vertlist* l, int index)
{
    vertnode* first = l->first;
    int i;

    for (i = 0; i < index; ++i)
        first = first->next;

    l->first = first;
}
