/******************************************************************************
 *
 * File:           gridmap.h
 *  
 * Created:        Thu Jan 23 14:00:00 EST 1997
 *  
 * Author:         Daniel Delbourgo/Stephen Walker/Jason Waring
 *                 CSIRO Marine Research
 *  
 * Purpose:        Calculates transformations between physical and index
 *                 space within a numerical grid
 *
 * Revisions:      110898 JRW
 *                 Added IJtoXY conversion and
 *                 fractional XYtoIJ and IJtoXY conversion
 *
 *                 2000 Pavel Sakov
 *                 Added branch calculation to handle both right- and
 *                 left-handed grids (calc_branch()).
 *
 *                 April 2002 Pavel Sakov
 *                 Major mods to handle topologically non-rectangular grids
 *
 *****************************************************************************/

#if !defined(_GRIDMAP_H)
#define _GRIDMAP_H

struct gridmap;
typedef struct gridmap gridmap;

gridmap* gridmap_build(int nce1, int nce2, double** gx, double** gy);
void gridmap_destroy(gridmap* gm);
int gridmap_fij2xy(gridmap* gm, double fi, double fj, double* x, double* y);
int gridmap_xy2ij(gridmap* gm, double x, double y, int* i, int* j);
int gridmap_xy2fij(gridmap* map, double x, double y, double* fi, double* fj);
int gridmap_getnce1(gridmap* gm);
int gridmap_getnce2(gridmap* gm);
void gridmap_getextent(gridmap* gm, double* xmin, double* xmax, double* ymin, double* ymax);

#endif
