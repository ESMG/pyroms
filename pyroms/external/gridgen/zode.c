/******************************************************************************
 *
 * File:           zode.c
 *
 * Created         4.4.2001
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Laboratories
 *
 * Purpose:        Solving complex ODEs
 *
 * Description:    Wrapper to "ode.c" for solving complex ODEs
 *
 * Revisions:      none
 *	    
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ode.h"
#include "zode.h"

typedef struct {
    int nz;
    zfunc zcalc;
    zdouble z0;
    zdouble z1;
    void* p;
    void* work;
} calcstruct;

static void calc(double x, double* y, double* y1, void* p)
{
    calcstruct* zs = (calcstruct*) p;
    zdouble* f = (zdouble*) zs->work;
    zdouble* f1 = &((zdouble*) zs->work)[zs->nz];
    zdouble z = (1.0 - x) * zs->z0 + x * zs->z1;
    int i, ii;

    for (i = 0, ii = 0; i < zs->nz; ++i, ii += 2) {
        f[i] = y[ii] + I * y[ii + 1];
        f1[i] = y1[ii] + I * y1[ii + 1];
    }

    zs->zcalc(z, f, f1, zs->p);

    for (i = 0, ii = 0; i < zs->nz; ++i, ii += 2) {
        y[ii] = creal(f[i]);
        y[ii + 1] = cimag(f[i]);
        y1[ii] = creal(f1[i]);
        y1[ii + 1] = cimag(f1[i]);
    }
}

/* Complex ODE integrator.
 * @param ode Real ODE integrator used
 * @param zcalc Derivative function
 * @param nz Dimension
 * @param z Initial argument value
 * @param f Function value (input/output)
 * @param zend Final argument value
 * @param Initial relative step guess (0 < h < 1)
 * @param eps Local precision
 * @param p Custom data; will be passed to `zcalc'
 * @return Passed from `ode'; 1 on success; 0 otherwhile
 */
int zode(integrator ode, zfunc zcalc, int nz, zdouble z, zdouble f[], zdouble zend, double* h, double eps, void* p)
{
    double* y = (double*) f;
    int status = -1;
    calcstruct cs;

    cs.nz = nz;
    cs.zcalc = zcalc;
    cs.z0 = z;
    cs.z1 = zend;
    cs.p = p;
    cs.work = malloc(sizeof(zdouble) * nz * 2);

    status = ode(calc, nz * 2, 0.0, y, 1.0, eps, 1.0, h, NULL, NULL, &cs);

    free(cs.work);

    return status;
}
