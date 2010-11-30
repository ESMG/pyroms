/******************************************************************************
 *
 * File:           vertlist.h
 *
 * Created:        19/03/2001
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Header file for vertlist.c
 *
 * Revisions:      None
 *
 *****************************************************************************/

#include "config.h"

extern int gg_verbose;

#define BUFSIZE 1024

#if !defined(_VERTLIST_H)
#define _VERTLIST_H

#if !defined(_POINT_STRUCT)
#define _POINT_STRUCT
typedef struct {
    double x;
    double y;
    double z;
} point;
#endif

struct vertnode;
typedef struct vertnode vertnode;

struct vertnode {
    vertnode* prev;
    vertnode* next;
    point p;
    int protected;
};

typedef struct {
    vertnode* first;
    int n;
} vertlist;

void vertlist_add(vertlist* l, point* p);
double vertlist_area(vertlist* l);
void vertlist_change_order(vertlist** l);
void vertlist_clear(vertlist* l);
vertlist* vertlist_create(void);
vertlist* vertlist_create2(int n, double x[], double y[], double z[]);
void vertlist_thin(vertlist* l, double dxmin, double dymin);
void vertlist_delete(vertlist* l, vertnode* v);
void vertlist_destroy(vertlist* l);
void vertlist_find_minmax(vertlist* l, double* xmin, double* xmax, double* ymin, double* ymax);
void vertlist_init(vertlist* l);
void vertlist_insert_after(vertlist* l, vertnode* v, point* p);
void vertlist_print(vertlist* l, FILE* f);
void vertlist_read(vertlist* l, FILE* f);
point* vertlist_topoint(vertlist* l);
zdouble* vertlist_tozdouble(vertlist* l);
void vertlist_toxy(vertlist* l, double** x, double** y);
void vertlist_process_phase1(vertlist* l);
void vertlist_process_phase2(vertlist* l);
delaunay* vertlist_triangulate(vertlist* l, FILE* log);
void vertlist_setfirstnode(vertlist* l, int index);

#endif
