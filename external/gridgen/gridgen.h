/******************************************************************************
 *
 * File:           gridgen.h
 *
 * Created:        2/11/2006
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Top-level header file for the gridgen library.
 *                 For more information, see README and gridgen.c.
 *
 * Description:    Contains prototypes of functions for running gridgen.
 *                 To use gridgen from your C code, #include this header
 *                 and link with libgridgen.a.
 *                 
 * Revisions:       15/02/2007 PS Introduced a new function gridgen_create2().
 *
 *****************************************************************************/

#if !defined(_GRIDGEN_H)
#define _GRIDGEN_H

#include "config.h"
#if defined(HAVE_GRIDNODES_H)
#include "gridnodes.h"
#endif

void gridgen_setverbose(int verbose);
void gridgen_printversion(void);
void gridgen_printhelpalg(void);
void gridgen_printhelpprm(void);
void gridgen_generategrid(char* prm);

#if defined(HAVE_GRIDNODES_H)

gridnodes* gridgen_generategrid2(int nbdry, double xbdry[], double ybdry[], double beta[], int ul, int nx, int ny, int ngrid, double xgrid[], double ygrid[], int nnodes, int newton, double precision, int checksimplepoly, int thin, int nppe, int verbose, int* nsigmas, double** sigmas, int* nrect, double** xrect, double** yrect);

#endif

#endif                          /* _GRIDGEN_H */
