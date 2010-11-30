/******************************************************************************
 *
 * File:           geom.c
 *
 * Created:        12/09/2003
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        A few geometrical utilities
 *
 * Revisions:      None
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "geom.h"

double point_point_distance(point* p1, point* p2)
{
    return hypot(p1->x - p2->x, p1->y - p2->y);
}

double point_edge_distance(point* p, point* p1, point* p2)
{
    double dx10 = p->x - p1->x;
    double dy10 = p->y - p1->y;
    double dx20 = p->x - p2->x;
    double dy20 = p->y - p2->y;
    double dx12 = p2->x - p1->x;
    double dy12 = p2->y - p1->y;

    if (dx10 * dx12 + dy10 * dy12 < 0.0)
        return hypot(dx10, dy10);

    if (dx20 * dx12 + dy20 * dy12 > 0.0)
        return hypot(dx20, dy20);

    return fabs(dx10 * dy20 - dx20 * dy10) / hypot(dx12, dy12);
}

int on_left_side(zdouble* z, zdouble* z0, zdouble* z1, double eps)
{
    double delta = cimag((*z1 - *z) * conj(*z0 - *z));

    if (delta > eps)
        return 1;
    else if (delta < -eps)
        return -1;
    return 0;
}

int in_triangle(zdouble* z, zdouble* z0, zdouble* z1, zdouble* z2, double eps)
{
    int l0 = on_left_side(z, z0, z1, eps);
    int l1 = on_left_side(z, z1, z2, eps);
    int l2 = on_left_side(z, z2, z0, eps);
    int sum = l0 + l1 + l2;

    if (sum == 3)
        return 1;
    if (l0 < 0 || l1 < 0 || l2 < 0)
        return 0;
    if (sum > 0)
        return 1;

    {
        double x = creal(*z);
        double y = cimag(*z);
        double xmin = DBL_MAX;
        double xmax = -DBL_MAX;
        double ymin = DBL_MAX;
        double ymax = -DBL_MAX;
        double xx[3], yy[3];
        int i;

        xx[0] = creal(*z0);
        xx[1] = creal(*z1);
        xx[2] = creal(*z2);
        yy[0] = cimag(*z0);
        yy[1] = cimag(*z1);
        yy[2] = cimag(*z2);

        for (i = 0; i < 3; ++i) {
            if (xx[i] < xmin)
                xmin = xx[i];
            if (xx[i] > xmax)
                xmax = xx[i];
            if (yy[i] < ymin)
                ymin = yy[i];
            if (yy[i] > ymax)
                ymax = yy[i];
        }
        if (x > xmax || x < xmin || y > ymax || y < ymin)
            return 0;
    }
    return 1;
}

int in_quadrilateral(zdouble* z, zdouble zs[], int vids[], double eps)
{
    return (in_triangle(z, &zs[vids[0]], &zs[vids[1]], &zs[vids[2]], eps) || in_triangle(z, &zs[vids[2]], &zs[vids[3]], &zs[vids[0]], eps));
}

int in_poly(zdouble z, int n, zdouble zs[], double eps)
{
    double x = creal(z);
    double y = cimag(z);
    int cnum;
    int i;

    for (i = 0, cnum = 0; i < n; ++i) {
        int i1 = (i + 1) % n;
        double x1 = creal(zs[i]) - x;
        double y1 = cimag(zs[i]) - y;
        double x2 = creal(zs[i1]) - x;
        double y2 = cimag(zs[i1]) - y;

        if (y1 == 0.0 && y2 == 0.0) {
            if (x1 * x2 <= 0)
                return 1;
        } else if (y1 == 0.0) {
            if (x1 == 0.0)
                return 1;
            if (x1 > 0.0)
                cnum += (y2 > 0.0) ? 1 : -1;
        } else if (y2 == 0.0) {
            if (x2 == 0.0)
                return 1;
            if (x2 > 0.0)
                cnum += (y1 < 0.0) ? 1 : -1;
        } else if (y1 * y2 < 0.0) {
            if (x1 > 0.0 && x2 > 0.0)
                cnum += 2;
            else if (x1 * x2 <= 0.0) {
                double xx = x1 - (x2 - x1) * y1 / (y2 - y1);

                if (xx == 0)
                    return 1;
                if (xx > 0.0)
                    cnum += 2;
            }
        }
    }

    if ((cnum / 2) % 2)
        return 1;

    for (i = 0; i < n; ++i) {
        point p;
        point p1;
        point p2;

        p.x = creal(z);
        p.y = cimag(z);
        p1.x = creal(zs[i]);
        p1.y = cimag(zs[i]);
        p2.x = creal(zs[(i + 1) % n]);
        p2.y = cimag(zs[(i + 1) % n]);

        if (point_edge_distance(&p, &p1, &p2) < eps)
            return 1;
    }

    return 0;
}
