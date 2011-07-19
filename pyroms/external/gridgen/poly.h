/******************************************************************************
 *
 * File:           polyline.h
 *  
 * Created:        Early 2002
 *  
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *  
 * Purpose:        Library routines for polylines
 *
 *****************************************************************************/

typedef struct {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
} extent;

typedef struct {
    int n;                      /* number of points */
    int nallocated;             /* number of allocated points */
    extent e;                   /* bounding rectangle */
    double* x;                  /* array of x coordinates [n] */
    double* y;                  /* array of y coordinates [n] */
} poly;

poly* poly_create();
void poly_destroy(poly* pl);

void poly_addpoint(poly* pl, double x, double y);
void poly_addpointat(poly* pl, int index, double x, double y);
void poly_append(poly* pl1, poly* pl2);
double poly_area(poly* pl);
void poly_clear(poly* pl);
void poly_close(poly* pl);
int poly_containspoint(poly* pl, double x, double y);
poly* poly_copy(poly* pl);
void poly_deletepoint(poly* pl, int index);
void poly_despike(poly* pl, double maxdist);
int poly_findindex(poly* pl, double x, double y);
poly* poly_formbound(int nce1, int nce2, double** x, double** y);
poly* poly_formboundij(int nce1, int nce2, double** x);
int poly_isclosed(poly* pl, double eps);
int poly_read(poly* pl, FILE* fp);
void poly_resample(poly* pl, double eps);
void poly_reverse(poly* pl);
void poly_write(poly* pl, FILE* fp);
void poly_compact(poly* pl, double eps);
