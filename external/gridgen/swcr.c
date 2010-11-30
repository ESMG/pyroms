/******************************************************************************
 *
 * File:           swcr.c
 *
 * Created:        19/03/2001
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Calculation of Schwarz-Christoffel transform
 *
 * Description:    Calculates Schwarz-Christoffel transform as described in:
 *                 Lloyd N. Trefethen, "Numerical Computation of the
 *                 Schwartz-Christoffel transformation", SIAM J. Sci. Stat.
 *                 Comput., 1980, 1(1), 82-102.
 *                 A fair amount of mission-critical code has been borrowed 
 *                 from SCPACK (ported from FORTRAN), 
 *                 see http://www.netlib.org/conformal.
 *
 * Revisions:      None
 *
 *****************************************************************************/

#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "zode.h"
#include "swcr.h"
#include "nan.h"

#define NEWTON_NITER_MAX 20
#define INTEGRATION_NINT_MAX 1000
#define WABS_OK_MAX 1.1
#define WABS_NEWTON_MAX 2.0
#define ZZERO 0.0 + 0.0 * I

#define min(x,y) ((x) < (y) ? (x) : (y))

#if HAVE_TGAMMA
double tgamma(double);
#elif HAVE_LGAMMA
double lgamma(double);
#endif

struct swcr {
    int n;
    int nq;
    double* betas;
    double* nodes;
    double* weights;
    double* work;
    int singular;               /* flag */
};

static void quit(char* format, ...)
{
    va_list args;

    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    exit(1);
}

/* Borrowed from "sclibdbl.f" (SCPACK).
 * This procedure supplies the coefficients a(j), b(j) of the recurrence
 * relation
 *
 *   b[j] p[j](x) = (x - a[j]) p[j-1](x) - b[j-1] p[j-2](x)
 *
 * for the various classical (normalized) orthogonal polynomials, and the
 * zero-th moment
 *
 *   muzero = integral w(x) dx
 *
 * of the given polynomial's weight function w(x). Since the polynomials are
 * ortho-normalized, the tridiagonal matrix is guaranteed to be symmetric.
 */
static double class(int n, double alpha, double beta, double* b, double* a)
{
    int nm1 = n - 1;
    double ab = alpha + beta;
    double abi = ab + 2.0;
    double abi2;
    double bbaa = beta * beta - alpha * alpha;
    double a1 = alpha + 1;
    double b1 = beta + 1;

#if HAVE_TGAMMA
    double muzero = pow(2.0, ab + 1.0) * tgamma(a1) * tgamma(b1) / tgamma(abi);
#elif HAVE_LGAMMA
    double muzero = pow(2.0, ab + 1.0) * exp(lgamma(a1) + lgamma(b1) - lgamma(abi));
#endif
    int i;

    a[0] = (beta - alpha) / abi;
    b[0] = sqrt(4.0 * a1 * b1 / (ab + 3.0) / abi / abi);

    for (i = 1; i < nm1; ++i) {
        int i1 = i + 1;

        abi = 2.0 * i1 + ab;
        abi2 = abi * abi;
        a[i] = bbaa / (abi - 2.0) / abi;
        b[i] = sqrt(4.0 * i1 * (alpha + i1) * (beta + i1) * (ab + i1) / (abi2 - 1.0) / abi2);
    }

    abi += 2.0;
    a[nm1] = bbaa / (abi - 2.0) / abi;
    b[nm1] = 0.0;

    return muzero;
}

/* Borrowed from "sclibdbl.f" (SCPACK).
 * This is a modified version of the eispack routine imtql2(). It finds the
 * eigenvalues and first components of the eigenvectors of a symmetric
 * tridiagonal matrix by the implicit ql method.
 */
static void imtql2(int n, double* d, double* e, double* z)
{
    int nm1 = n - 1;
    double b, c, f, g, p, r, s;
    int i, l;

    if (n <= 1)
        return;

    e[nm1] = 0.0;

    for (l = 0; l < n; ++l) {
        int j = 0;
        int m, mml;

        while (1) {
            for (m = l; m < nm1; ++m)
                if (fabs(e[m]) < UROUND * (fabs(d[m]) + fabs(d[m + 1])))
                    break;
            p = d[l];

            if (m == l)
                goto continue_l;

            if (j == 30)
                quit("error: imtql2(): number of iterations > 30\n");

            j++;

            g = (d[l + 1] - p) / e[l] / 2.0;
            r = sqrt(g * g + 1.0);
            g = d[m] - p + e[l] / (g + copysign(r, g));
            s = 1.0;
            c = 1.0;
            p = 0.0;
            mml = m - l;

            for (i = m - 1; i >= l; --i) {
                f = s * e[i];
                b = c * e[i];
                if (fabs(f) >= fabs(g)) {
                    c = g / f;
                    r = sqrt(c * c + 1.0);
                    e[i + 1] = f * r;
                    s = 1.0 / r;
                    c *= s;
                } else {
                    s = f / g;
                    r = sqrt(s * s + 1.0);
                    e[i + 1] = g * r;
                    c = 1.0 / r;
                    s *= c;
                }
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
                f = z[i + 1];
                z[i + 1] = s * z[i] + c * f;
                z[i] = c * z[i] - s * f;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }                       /* while (1) */
      continue_l:
        ;
    }

    for (i = 0; i < nm1; ++i) {
        int k = i;
        int j;

        p = d[i];

        for (j = i + 1; j < n; ++j) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }

        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            p = z[i];
            z[i] = z[k];
            z[k] = p;
        }
    }
}

/* Borrowed from "sclibdbl.f" (SCPACK).
 * This routine computes the nodes t[] and weights w[] for Gauss-Jacobi
 * quadrature formulas. These are used when one wishes to approximate
 *
 *   integral (from a to b) f(x) w(x) dx
 *
 * by
 *
 *    n
 *   sum w[j] f(t[j])
 *   j=1
 *
 * (w(x) and w[j] have no connection with each other), where w(x) is the weight
 * function,
 *
 *   w(x) = (1-x)**alpha * (1+x)**beta
 *
 * on (-1, 1); alpha, beta > -1.
 *
 * @param n     The number of points used for the quadrature rule
 * @param alpha Parameter of the weight function
 * @param beta  Parameter of the weight function
 * @param b     Work array [n]
 * @param t     Nodes [n] (output)
 * @param w     Weights [n] (output)
 *
 * Subroutines required: class, imtql2
 *
 * Reference:
 *
 * The routine has been adapted from the more general routine gaussq by Golub
 * and Welsch. See Golub, G. H., and Welsch, J. H., "Calculation of Gaussian
 * quadrature rules", Mathematics of computation, 23 (April 1969), pp. 221-230.
 */
static void gaussj(int n, double alpha, double beta, double* b, double* t, double* w)
{
    double mu0 = class(n, alpha, beta, b, t);
    int i;

    w[0] = 1.0;
    imtql2(n, t, b, w);

    for (i = 0; i < n; ++i)
        w[i] = mu0 * w[i] * w[i];
}

/* Creates Schwarz-Christoffel transform structure.
 * @param n Number of prevertices
 * @param nq Number of nodes in Gauss-Jacobi quadrature formulas
 * @param betas Array of angle parameters [n]: betas[i] = alpha[i]/pi - 1, 
 *              where alpha[i] -- interior angle at i-th vertex
 * @return Pointer to Schwarz-Christoffel transform structure
 */
swcr* sc_create(int n, int nq, double betas[])
{
    swcr* sc;
    int i, ii;

    if (nq <= 0)
        return NULL;

    sc = malloc(sizeof(swcr));

    sc->n = n;
    sc->nq = nq;
    sc->betas = betas;
    sc->nodes = malloc((n + 1) * nq * sizeof(double));
    sc->weights = calloc((n + 1) * nq, sizeof(double));
    sc->work = malloc(nq * sizeof(double));
    sc->singular = 0;

    for (i = 0; i < n; ++i) {
        ii = nq * i;
        if (betas[i] > -1.0)
            gaussj(nq, 0.0, betas[i], sc->work, &sc->nodes[ii], &sc->weights[ii]);
    }
    ii = nq * n;
    gaussj(nq, 0.0, 0.0, sc->work, &sc->nodes[ii], &sc->weights[ii]);

    return sc;
}

/* Destroys Schwarz-Christoffel transform structure.
 * @param sc Structure to be destroyed
 */
void sc_destroy(swcr* sc)
{
    if (sc == NULL)
        return;
    free(sc->nodes);
    free(sc->weights);
    free(sc->work);
    free(sc);
}

/* Borrowed from "scpdbl.f" (SCPACK). Originally named zprod().
 * Computes the integrand
 *
 *     n
 *   prod (1-w/ws[i])**beta(i),
 *    i=1
 *
 * taking argument only (not modulus) for term i = k.
 * This is the innermost subroutine in `gridgen' calculations. The complex log
 * calculation below may account for as much as half of the total execution
 * time.
 */
static zdouble sc_prod(swcr* sc, zdouble w, int k, zdouble ws[])
{
    int n = sc->n;
    double* betas = sc->betas;
    zdouble wsum = 0.0;
    int i;

    for (i = 0; i < n; ++i) {
        if (betas[i] != 0.0) {
            zdouble ztmp = 1.0 - w / ws[i];

            if (i == k)
                ztmp /= cabs(ztmp);
            wsum += betas[i] * clog(ztmp);
        }
    }

    return cexp(wsum);
}

/* Borrowed from "scpdbl.f" (SCPACK).
 * Determines the distance from z0 to the nearest singularity w[] other than
 * w[i0].
 *
 * @param  z0 Distance from this point
 * @param  i0 Index to be excluded
 * @param  n  Dimension of `w'
 * @param  w  Complex array [n]
 * @return Distance
 */
static double dist(zdouble z0, int i0, int n, zdouble w[])
{
    int i;
    double d = DBL_MAX;

    for (i = 0; i < n; ++i) {
        double dd;

        if (i == i0)
            continue;
        dd = min(d, cabs(z0 - w[i]));
        if (dd < d)
            d = dd;
    }

    return d;
}

/* Borrowed from "scpdbl.f" (SCPACK). Originally named zqsum().
 * Computes the integral of sc_prod() from `z0' to `z1' by applying a one-sided
 * Gauss-Jacobi formula with possible singularity at w[i0].
 */
static zdouble sc_sum(swcr* sc, zdouble z0, zdouble z1, int i0, zdouble w[])
{
    zdouble zsum = ZZERO;
    zdouble zh = (z1 - z0) / 2.0;
    zdouble zc = (z1 + z0) / 2.0;
    int nq = sc->nq;
    int k = (i0 < 0) ? sc->n : i0;
    double* nodes = &sc->nodes[nq * k];
    double* weights = &sc->weights[nq * k];
    int i;

    for (i = 0; i < nq; ++i)
        zsum += weights[i] * sc_prod(sc, zc + zh * nodes[i], i0, w);
    zsum *= zh;

    if (cabs(zh) != 0.0 && i0 >= 0)
        zsum *= pow(cabs(zh), sc->betas[k]);

    return zsum;
}

/* Borrowed from "scpdbl.f" (SCPACK). Originally named zquad1().
 * Computes the complex line integral of sc_prod() from z0 to z1 along a
 * straight line segment within the unit disk. Compound one-sided Gauss-Jacobi
 * quadrature is used, using dist() to determine the distance to the nearest
 * singularity w[].
 */
static zdouble sc_integrate_(swcr* sc, zdouble z0, zdouble z1, int i0, zdouble w[])
{
    zdouble zsum, z00, z11;
    double r;
    int count = 0;

    r = min(1.0, dist(z0, i0, sc->n, w) / cabs(z0 - z1));

    if (r == 0.0)
        sc->singular = 1;

    z00 = z0 + r * (z1 - z0);
    zsum = sc_sum(sc, z0, z00, i0, w);

    while (r != 1.0) {
        /*
         * This check is a hack: sometimes (very rarely) the cycle hangs.
         * Need to sort this out. 
         */
        if (count == INTEGRATION_NINT_MAX)
            break;

        r = min(1.0, dist(z00, -1, sc->n, w) / cabs(z00 - z1));
        z11 = z00 + r * (z1 - z00);
        zsum += sc_sum(sc, z00, z11, -1, w);
        z00 = z11;
        count++;
    }

    return zsum;
}

/* Borrowed from "scpdbl.f" (SCPACK). Originally named zquad().
 * Computes the complex line integral of sc_prod() from `w0' to `w1' along a
 * straight line segment within the unit disk. sc_integrate_() is called twice,
 * once for each half of this integral.
 */
static zdouble sc_integrate(swcr* sc, zdouble w0, int i0, zdouble w1, int i1, zdouble w[])
{
    zdouble wmid = (w0 + w1) / 2.0;

    if (cabs(w0) > WABS_OK_MAX || cabs(w1) > WABS_OK_MAX)
        quit("error: sc_integrate(): w outside unit disk\n");

    if (cabs(w0 - w1) == 0.0)
        return 0.0;

    return sc_integrate_(sc, w1, wmid, i1, w) - sc_integrate_(sc, w0, wmid, i0, w);
}

/* Schwarz-Christoffel transform.
 *
 *                    w       n          beta[j]
 * z(w) = z0 + B * integral prod(1-w/w[j])      dw
 *                    w0     j=1     
 *
 * @param sc Schwarz-Christoffel transform structure
 * @param w Point to be mapped
 * @param k Index of prevertex coinciding with `w' if such exists; -1
 *          otherwhile
 * @param w0 Point: w0 -> z0
 * @param z0 Point: w0 -> z0
 * @param k0 Index of prevertex coinciding with `w0' if such exists; -1
 *           otherwhile
 * @param B Integral multiple
 * @param wi Prevertices [sc->n]
 * @param eps Allowable overshoot beyond unit disk
 * @return Image of `w'
 */
zdouble sc_w2z(swcr* sc, zdouble w, int k, zdouble w0, zdouble z0, int k0, zdouble B, zdouble wi[])
{
    return z0 + B * sc_integrate(sc, w, k, w0, k0, wi);
}

typedef struct {
    swcr* sc;
    zdouble B;
    zdouble* w;
} odestruct;

static void f(zdouble z, zdouble f[], zdouble f1[], void* p)
{
    odestruct* os = (odestruct*) p;
    int n = os->sc->n;
    double* betas = os->sc->betas;
    zdouble* w = os->w;
    zdouble f0 = f[0];
    zdouble wsum = 0.0;
    int i;

    for (i = 0; i < n; ++i)
        if (betas[i] != 0.0)
            wsum += betas[i] * clog(1.0 - f0 / w[i]);

    f1[0] = cexp(-wsum) / os->B;
}

/* Inverse Schwarz-Christoffel transform.
 *
 * Assumes that neither z0 not z coincide with any of the polygon image 
 * vertices (singular points). (In `gridgen', it is easier to check this
 * outside `sc_z2w'. BUT, perhaps, a more robust solution exists, as the
 * check may give wrong results because of the numerical precision-related 
 * errors.)
 *
 * @param sc Schwarz-Christoffel transform structure
 * @param z Point to be mapped
 * @param z1 Point: z1 -> w1 (used in startup)
 * @param w1 Point: z1 -> w1 (used in startup)
 * @param z2 Point: z2 -> w2 (used in startup)
 * @param w2 Point: z2 -> w2 (used in startup)
 * @param z0 Point: z0 -> w0 (used in integrations)
 * @param w0 Point: z0 -> w0 (used in integrations)
 * @param B Integral multiple
 * @param wi Prevertices [sc->n]
 * @param eps Precision
 * @param status Status: 0 if found by secant method; 1 if by ODE solving; 2 if fails
 * @return Image of `z'
 */
zdouble sc_z2w(swcr* sc, zdouble z, zdouble z1, zdouble w1, zdouble z2, zdouble w2, zdouble z0, zdouble w0, zdouble B, zdouble wi[], double eps, int* status)
{
    int count = 0;
    int ok = 1;
    zdouble w = w2;

    *status = 0;

    /*
     * Start with secant method first. It will occasionally fail because the 
     * function structure may be too fine, or because of the branch cuts
     * of complex logs starting from prevertices w_i. 
     */
    if (cabs(w1 - w2) > eps) {
        while (cabs(w1 - w2) > eps) {
            double wabs;

            w = (w2 * (z - z1) - w1 * (z - z2)) / (z2 - z1);
            wabs = cabs(w);

            if (wabs > WABS_NEWTON_MAX || count >= NEWTON_NITER_MAX) {
                ok = 0;
                break;
            }
            /*
             * My understanding is that w(z) is analytical inside unit disk, 
             * has branch points in vertice images on the unit
             * circumcircle and branch cuts from these to infinity. Hence,
             * if w shoots outside unit disk, it may be reasonably used in
             * some circumstances but not others. Therefore, we allow
             * overshoots by making WABS_OK_MAX > 1. BUT, because of that,
             * if unlucky, the secant method may converge to some value
             * outside the unit disk. Such values will be discarded and ODE
             * method tried. 
             */
            else if (wabs > WABS_OK_MAX)
                w /= wabs;
            w1 = w2;
            z1 = z2;
            w2 = w;
            z2 = sc_w2z(sc, w, -1, w0, z0, -1, B, wi);
            count++;
        }

        /*
         * (see the comment above) 
         */
        if (cabs(w) > 1.0 + eps)
            ok = 0;

        if (ok) {
            zdouble w = (w2 * (z - z1) - w1 * (z - z2)) / (z2 - z1);

            if (cabs(w - w2) < eps)
                return w;
            return w2;
        }

        *status = 1;
    }

    /*
     * So, the secant method has failed, lets try solving the ODE. This is
     * pretty robust although much slower. In very rare cases can fail in
     * vicinity of prevertices because of numerical errors. 
     */
    {
        static double h = 1.0;  /* try to solve ODE in one step */
        odestruct os;

        os.sc = sc;
        os.B = B / (z - z0);
        os.w = wi;

        if (!zode(dopri8, f, 1, z0, &w0, z, &h, eps / 10.0, &os)) {
            *status = 2;
            return NaN;
        }
    }

    if (cabs(w0) > 1.0)
        w0 /= cabs(w0);

    return w0;
}

int sc_issingular(swcr* sc)
{
    return sc->singular;
}
