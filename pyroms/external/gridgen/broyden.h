/******************************************************************************
 *
 * File:           broyden.h
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

#if !defined(_BROYDEN_H)
#define _BROYDEN_H

typedef void (*func) (double* x, double* f, void* p);

/* Makes one iteration of the Gauss-Newton nonlinear solver with Broyden
 * update.
 * @param F Function 
 * @param n System dimension
 * @param x Argument [n] (input/output)
 * @param f F(x) [n] (input/output)
 * @param W Negative inverse Jacobian approximation [n^2] (input/output)
 * @param p Custom data; will be passed to `F'
 *
 * Broyden method:
 *   xnew = x + W f
 *   fnew = F(xnew)
 *   Wnew = W - (W fnew) (s^T W) / s^T W (fnew - f),
 * where
 *   s = xnew - x (= Wf)
 */
void broyden_update(func F, int n, double* x, double* f, double* W, void* custom);

#endif
