/******************************************************************************
 *
 * File:        issimplepoly.h
 *
 * Created:     24/02/2003
 *
 * Author:      Pavel Sakov
 *              CSIRO Marine Research
 *
 * Description: Testing a polygon on self-intersections by using Shamos-Hoey 
 *              ("sweep-line") algorithm.
 *
 *****************************************************************************/

#if !defined(_ISSIMPLEPOLY_H)
#define _ISSIMPLEPOLY_H

int issimplepolygon(int n, double x[], double y[]);

#endif
