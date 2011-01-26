/******************************************************************************
 *
 * File:           gridaverager.c
 *  
 * Created:        Thu Aug 28 10:00:07 EST 2003
 *  
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *  
 * Purpose:        Calculate average value of a 2D field a cell of a grid
 *
 * Revisions:      None
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "nan.h"
#include "gridmap.h"
#include "gucommon.h"
#include "gridaverager.h"

struct gridaverager {
    gridmap* gm;
    double** v;
    int** n;
};

gridaverager* ga_create(gridmap* gm)
{
    gridaverager* ga = malloc(sizeof(gridaverager));

    ga->gm = gm;
    ga->v = gu_alloc2d(gridmap_getnce1(gm), gridmap_getnce2(gm), sizeof(double));
    ga->n = gu_alloc2d(gridmap_getnce1(gm), gridmap_getnce2(gm), sizeof(int));

    return ga;
}

void ga_destroy(gridaverager* ga)
{
    gu_free2d(ga->v);
    gu_free2d(ga->n);
    free(ga);
}

void ga_addpoints(gridaverager* ga, int n, point points[])
{
    gridmap* gm = ga->gm;
    int ii;

    for (ii = 0; ii < n; ++ii) {
        point* p = &points[ii];
        int i, j;

        if (gridmap_xy2ij(gm, p->x, p->y, &i, &j) != 0) {
            ga->v[j][i] += p->z;
            ga->n[j][i]++;
        }
    }
}

void ga_getvalue(gridaverager* ga, point * p)
{
    int i, j, n;

    if (gridmap_xy2ij(ga->gm, p->x, p->y, &i, &j) != 0 && (n = ga->n[j][i]) != 0)
        p->z = ga->v[j][i] / (double) n;
    else
        p->z = NaN;
}
