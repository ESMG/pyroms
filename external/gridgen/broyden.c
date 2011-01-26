/******************************************************************************
 *
 * File:           broyden.c
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
#include <assert.h>
#include "broyden.h"

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
void broyden_update(func F, int n, double* x, double* f, double* W, void* custom)
{
    double* fnew = calloc(n * 4, sizeof(double));
    double* s = &fnew[n];
    double* stw = &fnew[n * 2];
    double* wf1 = &fnew[n * 3];
    double denom;
    double* wij;
    int i, j;

    for (j = 0, wij = W; j < n; ++j) {
        for (i = 0; i < n; ++i, ++wij)
            s[j] += wij[0] * f[i];
        x[j] += s[j];
    }

    F(x, fnew, custom);

    for (j = 0, wij = W; j < n; ++j) {
        double sj = s[j];

        for (i = 0; i < n; ++i, ++wij)
            stw[i] += sj * wij[0];
    }

    for (i = 0, denom = 0.0; i < n; ++i)
        denom += stw[i] * (fnew[i] - f[i]);

    for (j = 0, wij = W; j < n; ++j)
        for (i = 0; i < n; ++i, ++wij)
            wf1[j] += wij[0] * fnew[i];

    for (j = 0, wij = W; j < n; ++j) {
        wf1[j] /= denom;
        for (i = 0; i < n; ++i, ++wij)
            wij[0] -= wf1[j] * stw[i];
    }

    for (i = 0; i < n; ++i)
        f[i] = fnew[i];

    free(fnew);
}

#if defined(TEST_BROYDEN)

#include <math.h>
#include <stdio.h>

#define N 5
#define COUNT_MAX 100
double d[] = { 3.0, 2.0, 1.5, 1.0, 0.5 };
double c = 0.01;

static void F(double* x, double* f, void* p)
{
    int i;

    for (i = 0; i < 5; ++i)
        f[i] = -x[i] * (d[i] + c * x[i] * x[i]);
}

int simple = 1;

int main(int argc, char* argv[])
{
    int count = 0;
    double x[N] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    double y[N];
    double w[N * N];
    double error;
    int i;

    for (i = 0; i < N * N; ++i)
        w[i] = 0.0;
    for (i = 0; i < N; ++i)
        w[i * N + i] = 1.0;

    printf("  iteration: log10(error)\n");
    do {
        if (count == 0)
            F(x, y, NULL);
        else
            broyden_update(F, N, x, y, w, NULL);
        for (i = 0, error = 0.0; i < N; ++i)
            error += y[i] * y[i];
        error = log10(error) / 2.0;
        printf("  %d: %.3g\n", count, error);
        count++;
    } while (error > -7.0 && count < COUNT_MAX);

    return 0;
}
#endif                          /* TEST_BROYDEN */
