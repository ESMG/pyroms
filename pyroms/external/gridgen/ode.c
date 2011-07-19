/******************************************************************************
 *
 * File:           ode.c
 *
 * Created         27.06.2001
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Laboratories
 *
 * Purpose:        ODE integration.
 *
 * Description:    Contains dopri8 ODE integrator ported from FORTRAN77
 *                 -- see E.Hairer, S.P. Norsett, G. Wanner, "Solving Ordinary
 *                 Differential Equations I. Nonstiff Problems", 
 *                 Springer-Verlag, 1987.
 *
 * Revisions:      none
 *	    
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "ode.h"

#if !defined(max)
#define max(x,y) ((x)>(y) ? (x) : (y) )
#endif

#if !defined(min)
#define min(x,y) ((x)<(y) ? (x) : (y) )
#endif

int ode_silent = 0;
int ode_stoponnan = 1;
int ode_stoponlong = 1;

static void stats(int nfcn, int nstep, int naccpt, int nrejct);

/* dopri8():
 *
 * Numerical solution of a system of first order ordinary differential
 * equations Y'=F(X,Y). This is an embedded Runge-Kutta method of order7(8)
 * due to Dormand & Prince (with stepsize control).
 * 
 * Ported from FORTRAN77 code. E.Hairer, S.P. Norsett, G. Wanner, "Solving
 * Ordinary Differential Equations I. Nonstiff Problems", Springer-Verlag,
 * 1987.
 * 
 * Excellent for intermediate precisions (1e-8 to 1e-13) and smooth functions.
 *
 * The call to `dopri8()' normally returns 1; 0 otherwise.
 * The integration results are stored in the function vector passed. The
 * procedure also stores the last valid stepsize in the initial stepsize
 * guess value passed so it could be used at the next call if necessary.
 *
 * NOTICE. The derivative arrays supplied as arguments to `calc' may contain
 * some noise, it is a duty of the derivative calculation routine to set ALL
 * of the return values.
 */

/* Smallest number satisfying 1.0 + UROUND > 1.0. */
#define UROUND 1.11e-16

/* Maximal number of steps. */
#define NMAX 2000

/* Coefficients. */
#define C2    (1.0 / 18.0)
#define C3    (1.0 / 12.0)
#define C4    (1.0 / 8.0)
#define C5    (5.0 / 16.0)
#define C6    (3.0 / 8.0)
#define C7    (59.0 / 400.0)
#define C8    (93.0 / 200.0)
#define C9    (5490023248.0 / 9719169821.0)
#define C10   (13.0 / 20.0)
#define C11   (1201146811.0 / 1299019798.0)
#define C12   (1.0)
#define C13   (1.0)
#define A21   (C2)
#define A31   (1.0 / 48.0)
#define A32   (1.0 / 16.0)
#define A41   (1.0 / 32.0)
#define A43   (3.0 / 32.0)
#define A51   (5.0 / 16.0)
#define A53   (-75.0 / 64.0)
#define A54   (-A53)
#define A61   (3.0 / 80.0)
#define A64   (3.0 / 16.0)
#define A65   (3.0 / 20.0)
#define A71   (29443841.0 / 614563906.0)
#define A74   (77736538.0 / 692538347.0)
#define A75   (-28693883.0 / 1125.0e+6)
#define A76   (23124283.0 / 18.0e+8)
#define A81   (16016141.0 / 946692911.0)
#define A84   (61564180.0 / 158732637.0)
#define A85   (22789713.0 / 633445777.0)
#define A86   (545815736.0 / 2771057229.0)
#define A87   (-180193667.0 / 1043307555.0)
#define A91   (39632708.0 / 573591083.0)
#define A94   (-433636366.0 / 683701615.0)
#define A95   (-421739975.0 / 2616292301.0)
#define A96   (100302831.0 / 723423059.0)
#define A97   (790204164.0 / 839813087.0)
#define A98   (800635310.0 / 3783071287.0)
#define A101  (246121993.0 / 1340847787.0)
#define A104  (-37695042795.0 / 15268766246.0)
#define A105  (-309121744.0 / 1061227803.0)
#define A106  (-12992083.0 / 490766935.0)
#define A107  (6005943493.0 / 2108947869.0)
#define A108  (393006217.0 / 1396673457.0)
#define A109  (123872331.0 / 1001029789.0)
#define A111  (-1028468189.0 / 846180014.0)
#define A114  (8478235783.0 / 508512852.0)
#define A115  (1311729495.0 / 1432422823.0)
#define A116  (-10304129995.0 / 1701304382.0)
#define A117  (-48777925059.0 / 3047939560.0)
#define A118  (15336726248.0 / 1032824649.0)
#define A119  (-45442868181.0 / 3398467696.0)
#define A1110 (3065993473.0 / 597172653.0)
#define A121  (185892177.0 / 718116043.0)
#define A124  (-3185094517.0 / 667107341.0)
#define A125  (-477755414.0 / 1098053517.0)
#define A126  (-703635378.0 / 230739211.0)
#define A127  (5731566787.0 / 1027545527.0)
#define A128  (5232866602.0 / 850066563.0)
#define A129  (-4093664535.0 / 808688257.0)
#define A1210 (3962137247.0 / 1805957418.0)
#define A1211 (65686358.0 / 487910083.0)
#define A131  (403863854.0 / 491063109.0)
#define A134  (-5068492393.0 / 434740067.0)
#define A135  (-411421997.0 / 543043805.0)
#define A136  (652783627.0 / 914296604.0)
#define A137  (11173962825.0 / 925320556.0)
#define A138  (-13158990841.0 / 6184727034.0)
#define A139  (3936647629.0 / 1978049680.0)
#define A1310 (-160528059.0 / 685178525.0)
#define A1311 (248638103.0 / 1413531060.0)
#define B1    (14005451.0 / 335480064.0)
#define B6    (-59238493.0 / 1068277825.0)
#define B7    (181606767.0 / 758867731.0)
#define B8    (561292985.0 / 797845732.0)
#define B9    (-1041891430.0 / 1371343529.0)
#define B10   (760417239.0 / 1151165299.0)
#define B11   (118820643.0 / 751138087.0)
#define B12   (-528747749.0 / 2220607170.0)
#define B13   (1.0 / 4.0)
#define BH1   (13451932.0 / 455176623.0)
#define BH6   (-808719846.0 / 976000145.0)
#define BH7   (1757004468.0 / 5645159321.0)
#define BH8   (656045339.0 / 265891186.0)
#define BH9   (-3867574721.0 / 1518517206.0)
#define BH10  (465885868.0 / 322736535.0)
#define BH11  (53011238.0 / 667516719.0)
#define BH12  (2.0 / 45.0)

int dopri8(fluxfn calc,         /* function computing the first derivatives */
           int n,               /* dimension of the system */
           double x,            /* initial X-value */
           double* y,           /* Y-values */
           double xend,         /* final X-value */
           double eps,          /* local tolerance */
           double hmax,         /* maximal stepsize */
           double* h0,          /* initial stepsize guess */
           intout out,          /* output procedure */
           intfinal final,      /* final procedure */
           void* custom_data)
{
    /*
     * number of function evaluations 
     */
    int nfcn = 0;

    /*
     * number of computed steps 
     */
    int nstep = 0;

    /*
     * number of accepted steps 
     */
    int naccpt = 0;

    /*
     * Number of rejected steps. Please note that the number of rejected
     * steps does not increase when Dopri8 has to reduce the initial step
     * size, while the number of computed steps increases. Hence, it is
     * possible to have nstep > naccpt + nrejct. 
     */
    int nrejct = 0;

    /*
     * direction from x to xmax 
     */
    double posneg = (xend > x) ? 1.0 : -1.0;

    /*
     * if the step has been rejected 
     */
    int reject = 0;

    /*
     * work arrays 
     */
    double* k = malloc(sizeof(double) * n * 8);

    double* k1 = k;
    double* k2 = &k[n];
    double* k3 = &k[n * 2];
    double* k4 = &k[n * 3];
    double* k5 = &k[n * 4];
    double* k6 = &k[n * 5];
    double* k7 = &k[n * 6];
    double* y1 = &k[n * 7];

    double h = *h0;

    /*
     * flag -- whether the step has been truncated to match the end point 
     */
    int trunc = 0;              /* init to eliminate gcc warning */
    double xph, err, fac, hnew;
    int i;

    eps = max(eps, 13.0 * UROUND);
    h = min(max(1.0e-10, fabs(h)), hmax) * posneg;

    if (out != NULL)
        out(1, n, x, y, custom_data);

    /*
     * main cycle 
     */
    while (((x - xend) * posneg + UROUND) <= 0.0) {
        if (ode_stoponlong && nstep > NMAX) {
            if (!ode_silent) {
                fprintf(stderr, "dopri8(): Could not proceed at x = %.5g: nstep > NMAX.\n", x);
                stats(nfcn, nstep, naccpt, nrejct);
            }
            free(k);
            return 0;
        }
        if ((x + 0.03 * h) == x) {
            if (!ode_silent) {
                fprintf(stderr, "dopri8(): Could not proceed at x = %.5g : (x + 0.03 * step) == x.\n", x);
                stats(nfcn, nstep, naccpt, nrejct);
            }
            free(k);
            return 0;
        }
        if (!reject) {
            if (((x + h - xend) * posneg) > 0.0) {
                h = xend - x;
                trunc = 1;
            } else
                trunc = 0;

            /*
             * First nine stages. 
             */
            calc(x, y, k1, custom_data);
        }
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * A21 * k1[i];

        calc(x + C2 * h, y1, k2, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);

        calc(x + C3 * h, y1, k3, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (A41 * k1[i] + A43 * k3[i]);

        calc(x + C4 * h, y1, k4, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (A51 * k1[i] + A53 * k3[i] + A54 * k4[i]);

        calc(x + C5 * h, y1, k5, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (A61 * k1[i] + A64 * k4[i] + A65 * k5[i]);

        calc(x + C6 * h, y1, k6, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (A71 * k1[i] + A74 * k4[i] + A75 * k5[i] + A76 * k6[i]);

        calc(x + C7 * h, y1, k7, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (A81 * k1[i] + A84 * k4[i] + A85 * k5[i] + A86 * k6[i] + A87 * k7[i]);

        calc(x + C8 * h, y1, k2, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (A91 * k1[i] + A94 * k4[i] + A95 * k5[i] + A96 * k6[i] + A97 * k7[i] + A98 * k2[i]);

        calc(x + C9 * h, y1, k3, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (A101 * k1[i] + A104 * k4[i] + A105 * k5[i] + A106 * k6[i] + A107 * k7[i] + A108 * k2[i] + A109 * k3[i]);

        /*
         * compute intermediate sums 
         */
        for (i = 0; i < n; ++i) {
            double y11s = A111 * k1[i] + A114 * k4[i] + A115 * k5[i] + A116 * k6[i] + A117 * k7[i] + A118 * k2[i] + A119 * k3[i];
            double y12s = A121 * k1[i] + A124 * k4[i] + A125 * k5[i] + A126 * k6[i] + A127 * k7[i] + A128 * k2[i] + A129 * k3[i];

            k4[i] = A131 * k1[i] + A134 * k4[i] + A135 * k5[i] + A136 * k6[i] + A137 * k7[i] + A138 * k2[i] + A139 * k3[i];
            k5[i] = B1 * k1[i] + B6 * k6[i] + B7 * k7[i] + B8 * k2[i] + B9 * k3[i];
            k6[i] = BH1 * k1[i] + BH6 * k6[i] + BH7 * k7[i] + BH8 * k2[i] + BH9 * k3[i];
            k2[i] = y11s;
            k3[i] = y12s;
        }

        /*
         * last 4 stages 
         */
        calc(x + C10 * h, y1, k7, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (k2[i] + A1110 * k7[i]);

        calc(x + C11 * h, y1, k2, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (k3[i] + A1210 * k7[i] + A1211 * k2[i]);

        xph = x + h;
        calc(xph, y1, k3, custom_data);
        for (i = 0; i < n; ++i)
            y1[i] = y[i] + h * (k4[i] + A1310 * k7[i] + A1311 * k2[i]);

        calc(xph, y1, k4, custom_data);
        for (i = 0; i < n; ++i) {
            k5[i] = h * (k5[i] + B10 * k7[i] + B11 * k2[i] + B12 * k3[i] + B13 * k4[i]);
            k6[i] = k5[i] - h * (k6[i] + BH10 * k7[i] + BH11 * k2[i] + BH12 * k3[i]);
            k5[i] += y[i];
        }

        nfcn += 13;

        /*
         * error estimation 
         */
        err = 0.0;
        for (i = 0; i < n; ++i) {
            double denom = max(max(1.0e-6, fabs(k5[i])), max(fabs(y[i]), 2.0 * UROUND / eps));

            denom = k6[i] / denom;
            err += denom * denom;
        }
        err = sqrt(err / n);

        if (isnan(err)) {
            if (ode_stoponnan) {
                if (!ode_silent)
                    fprintf(stderr, "dopri8: Stopped: NaN detected\n");
                free(k);
                return 0;
            }
        }

        /*
         * New step size. We require 0.333 <= hnew / h <= 6.0. 
         */
        if (isnan(err))
            fac = 3.0;
        else
            fac = max(1.0 / 6.0, min(3.0, pow(err / eps, 0.125) / 0.9));
        hnew = h / fac;

        if (err < eps) {
            naccpt++;
            for (i = 0; i < n; ++i)
                y[i] = k5[i];

            x = xph;
            if (out != NULL)
                out(naccpt + 1, n, x, y, custom_data);
            if (fabs(hnew) > hmax)
                hnew = posneg * hmax;
            if (reject)
                hnew = posneg * min(fabs(hnew), fabs(h));
            reject = 0;
            h = hnew;
            /*
             * store back the step size if this step has not been truncated
             * to match the end point or it was the first step 
             */
            if (naccpt == 1 || !trunc)
                *h0 = fabs(hnew);
        } else {
            reject = 1;
            h = hnew;
            if (naccpt >= 1)
                nrejct++;
            nfcn--;
        }

        nstep++;
    }

    if (final != NULL)
        final(n, x, y, k1, naccpt, nrejct, nfcn, custom_data);

    free(k);
    return 1;
}

static void stats(int nfcn, int nstep, int naccpt, int nrejct)
{
    fprintf(stderr, "  Number of function evaluations = %d.\n", nfcn);
    fprintf(stderr, "  Number of computed steps = %d.\n", nstep);
    fprintf(stderr, "  Number of accepted steps = %d.\n", naccpt);
    fprintf(stderr, "  Number of rejected steps = %d.\n", nrejct);
}
