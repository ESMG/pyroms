/******************************************************************************
 *
 * File:           geom.h
 *
 * Created:        12/09/2003
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        A header for few geometrical utilities
 *
 * Revisions:      None
 *
 *****************************************************************************/

#if !defined(_GEOM_H)
#define _GEOM_H

#include "config.h"

#if !defined(_POINT_STRUCT)
#define _POINT_STRUCT
typedef struct {
    double x;
    double y;
    double z;
} point;
#endif

double point_point_distance(point* p1, point* p2);
double point_edge_distance(point* p, point* p1, point* p2);
int on_left_side(zdouble* z, zdouble* z0, zdouble* z1, double eps);
int in_triangle(zdouble* z, zdouble* z0, zdouble* z1, zdouble* z2, double eps);
int in_quadrilateral(zdouble* z, zdouble zs[], int vids[], double eps);
int in_poly(zdouble z, int n, zdouble zs[], double eps);

#endif
