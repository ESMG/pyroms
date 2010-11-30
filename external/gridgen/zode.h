/******************************************************************************
 *
 * File:           zode.h
 *
 * Created:        24/02/2000
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Complex ODE integration
 *
 * Revisions:      None
 *
 * Description:    Header file for complex ODE integration
 *  
 *****************************************************************************/

#if !defined(_ZODE_H)
#define _ZODE_H

#include "ode.h"
#include "config.h"

typedef void (*zfunc) (zdouble z, zdouble f[], zdouble f1[], void* p);

/** Complex ODE integrator.
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
int zode(integrator ode, zfunc zcalc, int nz, zdouble z, zdouble f[], zdouble zend, double* h, double eps, void* p);

#endif
