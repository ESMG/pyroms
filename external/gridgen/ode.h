/******************************************************************************
 *
 * File:           ode.h
 *
 * Created:        24/02/2000
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Created         27.06.2001
 *
 * Purpose:        Header file for ODE integration.
 *
 * Revisions:      None.
 *
 *****************************************************************************/

#if !defined(_ODE_H)
#define _ODE_H

typedef void (*fluxfn) (double t, double* y, double* y1, void* p);

typedef void (*intout) (int nstep, int n, double x, double* y, void* p);

typedef void (*intfinal) (int n, double x, double* y, double* y1, int naccpt, int nrejct, int nfcn, void* p);

typedef int (*integrator) (fluxfn calc, int n, double x, double* y0, double xend, double eps, double hmax, double* h, intout out, intfinal final, void* p);

int dopri8(fluxfn calc, int n, double x, double* y, double xend, double eps, double hmax, double* h0, intout out, intfinal final, void* custom_data);

extern int ode_silent;          /* default = 0 */
extern int ode_stoponnan;       /* default = 1 */

#endif
