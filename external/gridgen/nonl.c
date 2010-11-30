/******************************************************************************
 *
 * File:           nonl.c
 *
 * Created:        24/02/2000
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Nonlinear solver.
 *
 * Revisions:      None.
 *
 * Description:    None.
 *  
 *****************************************************************************/

#include <stdlib.h>
#include "nonl.h"

/* Makes one iteration of the Gauss-Newton nonlinear solver with Broyden
 * update.
 * @param F Function 
 * @param n System dimension
 * @param x Argument vector [n] (input/output)
 * @param f Function values for `x' [n] (input/output)
 * @param W Negative inverse Jacobian approximation [n][n] (input/output)
 * @param p Custom data; will be passed to `F'
 *
 * Broyden method:
 *
 * xnew = x + Wf(x)
 * Wnew = W - (p + W y) p^T W / (p^T W y), where
 * p = xnew - x
 * y = fnew - f
 */
void newton_iteration(func F, int n, double* x, double* f, double* W, void* custom)
{
    double* fnew = malloc(n * sizeof(double));
    double* p = calloc(n, sizeof(double));
    double* y = NULL;
    double* pwy = calloc(n, sizeof(double));
    double* DW = calloc(n * n, sizeof(double));
    double denom;
    double* wij;
    double* dwij;
    int i, j, k;

    for (j = 0, wij = W; j < n; ++j) {
        for (i = 0; i < n; ++i, ++wij)
            p[j] += wij[0] * f[i];
        x[j] += p[j];
    }

    F(x, fnew, custom);

    y = f;
    for (i = 0; i < n; ++i)
        y[i] = fnew[i] - f[i];

    for (j = 0, wij = W, denom = 0.0; j < n; ++j) {
        double wy = 0.0;

        for (i = 0; i < n; ++i, ++wij)
            wy += wij[0] * y[i];
        denom += p[j] * wy;
        pwy[j] = p[j] + wy;
    }

    for (j = 0, dwij = DW; j < n; ++j) {
        double pwyj = pwy[j];

        for (i = 0; i < n; ++i, ++dwij)
            for (k = 0, wij = &W[i]; k < n; ++k, wij += n)
                dwij[0] += pwyj * p[k] * wij[0];
    }

    for (j = 0, wij = W, dwij = DW; j < n; ++j)
        for (i = 0; i < n; ++i, ++wij, ++dwij)
            wij[0] -= dwij[0] / denom;

    for (i = 0; i < n; ++i)
        f[i] = fnew[i];

    free(DW);
    free(p);
    free(fnew);
    free(pwy);
}
