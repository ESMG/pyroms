/******************************************************************************
 *
 * File:           swcr.h
 *
 * Created:        19/03/2001
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Header file for Schwarz-Christoffel transform
 *
 * Revisions:      None
 *
 *****************************************************************************/

#if !defined _SWCR_H
#define _SWCR_H

#include "config.h"

#define UROUND 1.11e-16

struct swcr;
typedef struct swcr swcr;

/** Creates Schwarz-Christoffel transform structure.
 *
 * @param n Number of prevertices
 * @param nq Number of nodes in Gauss-Jacobi quadrature formulas
 * @param betas Array of angle parameters [n]: betas[i] = alpha[i]/pi - 1, 
 *              where alpha[i] -- interior angle at i-th vertex
 * @return Pointer to Schwarz-Christoffel transform structure
 */
swcr* sc_create(int n, int nq, double betas[]);

/** Destroys Schwarz-Christoffel transform structure.
 * @param sc Structure to be destroyed
 */
void sc_destroy(swcr* sc);

/** Schwarz-Christoffel transform.
 *
 *                    w       n         beta[j]
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
 * @return Image of `w'
 */
zdouble sc_w2z(swcr* sc, zdouble w, int k, zdouble w0, zdouble z0, int k0, zdouble B, zdouble wi[]);

/** Inverse Schwarz-Christoffel transform.
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
 * @param status Status: 0 if found by secant method; 1 if by ODE solving; 2
 *               if fails
 * @return Image of `z'
 */
zdouble sc_z2w(swcr* sc, zdouble z, zdouble z1, zdouble w1, zdouble z2, zdouble w2, zdouble z0, zdouble w0, zdouble B, zdouble wi[], double eps, int* status);

int sc_issingular(swcr* sc);

#endif
