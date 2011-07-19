/******************************************************************************
 *
 * File:           nonl.h
 *
 * Created:        24/02/2000
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Header file for nonlinear solver.
 *
 * Revisions:      None.
 *
 * Description:    Header file for nonlinear solver.
 *  
 *****************************************************************************/

#if !defined(_NONL_H)
#define _NONL_H

typedef void (*func) (double* x, double* f, void* p);

/** Makes one iteration of the Gauss-Newton nonlinear solver with Broyden
 * update.
 * @param F Function 
 * @param n System dimension
 * @param x Argument vector [n] (input/output)
 * @param f Function values for `x' [n] (input/output)
 * @param W Negative inverse Jacobian approximation [n][n] (input/output)
 * @param p Custom data; will be passed to `F'
 */
void newton_iteration(func F, int n, double* x, double* f, double* W, void* p);

#endif
