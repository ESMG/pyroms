/******************************************************************************
 *
 * File:           delaunay.h
 *
 * Created:        04/08/2000
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Header for delaunay triangulation wrapper
 *
 * Description:    None
 *
 * Revisions:      None
 *
 *****************************************************************************/

#if !defined(_DELAUNAY_H)
#define _DELAUNAY_H

#if !defined(_POINT_STRUCT)
#define _POINT_STRUCT
typedef struct {
    double x;
    double y;
    double z;
} point;
#endif

typedef struct {
    int vids[3];
} triangle;

typedef struct {
    int npoints;
    point* points;

    int ntriangles;
    triangle* triangles;

    int* n_point_triangles;     /* n_point_triangles[i] is number of
                                 * triangles i-th point belongs to */
    int** point_triangles;      /* point_triangles[i][j] is index of j-th
                                 * triangle i-th pont belongs to */

    int nedges;
    int* edges;                 /* n-th edge is formed by points[edges[n*2]]
                                 * and points[edges[n*2+1]] */
} delaunay;

delaunay* delaunay_build(int np, point points[], int ns, int segments[], int nh, double holes[]);
void delaunay_destroy(delaunay* d);

#endif
