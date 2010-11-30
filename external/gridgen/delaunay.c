/******************************************************************************
 *
 * File:           delaunay.c
 *
 * Created:        04/08/2000
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Delaunay triangulation - a wrapper to triangulate().
 *
 * Description:    This is a stripped version from the Natural Neigbours 
 *                 interpolation library `nn' available at
 *                 http://www.marine.csiro.au/~sakov/nn.tar.gz
 *
 * Revisions:      None
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "triangle.h"
#include "delaunay.h"

extern int tr_verbose;

static void tio_init(struct triangulateio* tio)
{
    tio->pointlist = NULL;
    tio->pointattributelist = NULL;
    tio->pointmarkerlist = NULL;
    tio->numberofpoints = 0;
    tio->numberofpointattributes = 0;
    tio->trianglelist = NULL;
    tio->triangleattributelist = NULL;
    tio->trianglearealist = NULL;
    tio->neighborlist = NULL;
    tio->numberoftriangles = 0;
    tio->numberofcorners = 0;
    tio->numberoftriangleattributes = 0;
    tio->segmentlist = 0;
    tio->segmentmarkerlist = NULL;
    tio->numberofsegments = 0;
    tio->holelist = NULL;
    tio->numberofholes = 0;
    tio->regionlist = NULL;
    tio->numberofregions = 0;
    tio->edgelist = NULL;
    tio->edgemarkerlist = NULL;
    tio->normlist = NULL;
    tio->numberofedges = 0;
}

static void tio_destroy(struct triangulateio* tio)
{
    if (tio->pointlist != NULL)
        free(tio->pointlist);
    if (tio->pointattributelist != NULL)
        free(tio->pointattributelist);
    if (tio->pointmarkerlist != NULL)
        free(tio->pointmarkerlist);
    if (tio->trianglelist != NULL)
        free(tio->trianglelist);
    if (tio->triangleattributelist != NULL)
        free(tio->triangleattributelist);
    if (tio->trianglearealist != NULL)
        free(tio->trianglearealist);
    if (tio->neighborlist != NULL)
        free(tio->neighborlist);
    if (tio->segmentlist != NULL)
        free(tio->segmentlist);
    if (tio->segmentmarkerlist != NULL)
        free(tio->segmentmarkerlist);
    if (tio->holelist != NULL)
        free(tio->holelist);
    if (tio->regionlist != NULL)
        free(tio->regionlist);
    if (tio->edgelist != NULL)
        free(tio->edgelist);
    if (tio->edgemarkerlist != NULL)
        free(tio->edgemarkerlist);
    if (tio->normlist != NULL)
        free(tio->normlist);
}

/* Builds Delaunay triangulation of the given array of points.
 *
 * @param np Number of points
 * @param points Array of points [np] (input)
 * @param ns Number of forced segments
 * @param segments Array of (forced) segment endpoint indices [2*ns]
 * @param nh Number of holes
 * @param holes Array of hole (x,y) coordinates [2*nh]
 * @return Delaunay triangulation structure with triangulation results
 */
delaunay* delaunay_build(int np, point points[], int ns, int segments[], int nh, double holes[])
{
    delaunay* d = malloc(sizeof(delaunay));
    struct triangulateio tio_in;
    struct triangulateio tio_out;
    char cmd[64] = "eznC";
    int i, j;

    assert(sizeof(REAL) == sizeof(double));

    tio_init(&tio_in);

    if (np == 0) {
        free(d);
        return NULL;
    }

    tio_in.pointlist = malloc(np * 2 * sizeof(double));
    tio_in.numberofpoints = np;
    for (i = 0, j = 0; i < np; ++i) {
        tio_in.pointlist[j++] = points[i].x;
        tio_in.pointlist[j++] = points[i].y;
    }

    if (ns > 0) {
        tio_in.segmentlist = malloc(ns * 2 * sizeof(int));
        tio_in.numberofsegments = ns;
        memcpy(tio_in.segmentlist, segments, ns * 2 * sizeof(int));
    }

    if (nh > 0) {
        tio_in.holelist = malloc(nh * 2 * sizeof(double));
        tio_in.numberofholes = nh;
        memcpy(tio_in.holelist, holes, nh * 2 * sizeof(double));
    }

    tio_init(&tio_out);

    if (!tr_verbose)
        strcat(cmd, "Q");
    else if (tr_verbose > 1)
        strcat(cmd, "VV");
    if (ns != 0)
        strcat(cmd, "p");

    if (tr_verbose)
        fflush(stderr);

    /*
     * climax 
     */
    triangulate(cmd, &tio_in, &tio_out, NULL);

    if (tr_verbose)
        fflush(stderr);

    /*
     * I assume that all input points appear in tio_out in the same order as 
     * they were written to tio_in. I have seen no exceptions so far, even
     * if duplicate points were presented. Just to be reasonably sure, let
     * us make a couple of checks. 
     */
    assert(tio_out.numberofpoints == np);
    assert(tio_out.pointlist[2 * np - 2] == points[np - 1].x && tio_out.pointlist[2 * np - 1] == points[np - 1].y);

    d->npoints = np;
    d->points = malloc(np * sizeof(point));
    for (i = 0, j = 0; i < np; ++i) {
        point* p = &d->points[i];

        p->x = tio_out.pointlist[j++];
        p->y = tio_out.pointlist[j++];
        p->z = points[i].z;
    }
    if (tr_verbose) {
        fprintf(stderr, "input:\n");
        for (i = 0, j = 0; i < np; ++i) {
            point* p = &d->points[i];

            fprintf(stderr, "  %d: %15.7g %15.7g %15.7g\n", i, p->x, p->y, p->z);
        }
    }

    d->ntriangles = tio_out.numberoftriangles;
    d->triangles = malloc(d->ntriangles * sizeof(triangle));

    if (tr_verbose)
        fprintf(stderr, "triangles:\n");
    for (i = 0; i < d->ntriangles; ++i) {
        int offset = i * 3;
        triangle* t = &d->triangles[i];

        t->vids[0] = tio_out.trianglelist[offset];
        t->vids[1] = tio_out.trianglelist[offset + 1];
        t->vids[2] = tio_out.trianglelist[offset + 2];

        if (tr_verbose)
            fprintf(stderr, "  %d: (%d,%d,%d)\n", i, t->vids[0], t->vids[1], t->vids[2]);
    }

    d->n_point_triangles = calloc(d->npoints, sizeof(int));
    for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];

        for (j = 0; j < 3; ++j)
            d->n_point_triangles[t->vids[j]]++;
    }
    d->point_triangles = malloc(d->npoints * sizeof(int*));
    for (i = 0; i < d->npoints; ++i) {
        if (d->n_point_triangles[i] > 0)
            d->point_triangles[i] = malloc(d->n_point_triangles[i] * sizeof(int));
        else
            d->point_triangles[i] = NULL;
        d->n_point_triangles[i] = 0;
    }
    for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];

        for (j = 0; j < 3; ++j) {
            int vid = t->vids[j];

            d->point_triangles[vid][d->n_point_triangles[vid]] = i;
            d->n_point_triangles[vid]++;
        }
    }

    if (tio_out.edgelist != NULL) {
        d->nedges = tio_out.numberofedges;
        d->edges = malloc(d->nedges * 2 * sizeof(int));
        memcpy(d->edges, tio_out.edgelist, d->nedges * 2 * sizeof(int));
    } else {
        d->nedges = 0;
        d->edges = NULL;
    }

    tio_destroy(&tio_in);
    tio_destroy(&tio_out);

    return d;
}

/* Releases memory engaged in the Delaunay triangulation structure.
 *
 * @param d Structure to be destroyed
 */
void delaunay_destroy(delaunay* d)
{
    int i;

    for (i = 0; i < d->npoints; ++i)
        if (d->point_triangles[i] != NULL)
            free(d->point_triangles[i]);
    if (d->nedges > 0)
        free(d->edges);
    free(d->point_triangles);
    free(d->n_point_triangles);
    free(d->triangles);
    free(d->points);
    free(d);
}
