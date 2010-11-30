/******************************************************************************
 *
 * File:           gridaverager.h
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

#if !defined(_GRIDAVERAGER_H)
#define _GRIDAVERAGER_H

struct gridaverager;
typedef struct gridaverager gridaverager;

#if !defined(_POINT_STRUCT)
#define _POINT_STRUCT
typedef struct {
    double x;
    double y;
    double z;
} point;
#endif

gridaverager* ga_create(gridmap* gm);
void ga_destroy(gridaverager* ga);
void ga_addpoints(gridaverager* ga, int n, point points[]);
void ga_getvalue(gridaverager* ga, point * p);

#endif
