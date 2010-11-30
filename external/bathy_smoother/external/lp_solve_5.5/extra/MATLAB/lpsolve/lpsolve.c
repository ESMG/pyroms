/*
    Author:        Peter Notebaert
    Contact:       lpsolve@peno.be
    License terms: LGPL.

    Release notes:
      interface to the lp_solve 5.5 toolkit.

      This file and the MATLAB interface code is originally developed by
      Jeffrey C. Kantor (jeff@control.cheg.nd.edu, cchen@darwin.helion.nd.edu)
      for lp_solve version 2.3
      The original name was lpmex.c


      It is completely revised and redesigned by Peter Notebaert for lp_solve version 5.5

      lpsolve.c needs to be linked with hash.c and the lp_solve library.

      The hash function is needed to have a quick access from the provided functionname to the
      implementation of the function.

    Change history:
     v5.1.0.0
      - First implementation
     v5.5.0.1  12 nov 2005
      - set_outputfile(lp0, ""); added to disable messages on console
      - For routines that return an integer,
        CreateDoubleMatrix and SetDoubleMatrix replaced by CreateLongMatrix and SetLongMatrix
     v5.5.0.2  19 nov 2005
      - For routines that need an lp handle, this handle can now also be the model name.
      - New call get_handle to get handle from model name.
     v5.5.0.4  25 aug 2006
      - New routine guess_basis added.
*/

#define driverVERSION "5.5.0.4"

#include <stdio.h>
#include <signal.h>
#include <string.h>

#include "lpsolvecaller.h"

/* Declare a global lp record */

#define LPSTEP 100

static lprec   **lp;
static lprec   *lp0;
static hashtable *cmdhash, *handlehash;
static int     h;
static int     lp_last;
static int     result;

/* Some working storage */

static short   initialized = FALSE;
static short   interrupted = FALSE;
static char    *cmd; /* command */
static char    *errmsg; /* error message */

static Double  *dpr, *dpr0;
static Long    *ipr, *ipr0;

#define buf    errmsg

#define cmdsz     50
#define errmsgsz 200
#define bufsz    errmsgsz

typedef void (impl_routine)();

static void impl_set_obj_fn();

static jmp_buf exit_mark;

void exitnow()
{
  longjmp(exit_mark, -1);
}

static void Check_nrhs(char *functionname, int nrhs0, int nrhs)
{
	if (nrhs - 1 != nrhs0) {
                sprintf(errmsg, "%s requires %d argument%s.", functionname, nrhs0, (nrhs0 == 1) ? "" : "s");
		ErrMsgTxt(errmsg);
	}
}


/* callback function for lp_solve. Messages are reported via this routine and printed in the application */

static void __WINAPI mylog(lprec *lp, void *userhandle, char *buf)
{
  	Printf("%s", buf);
}


static int __WINAPI myabort(lprec *lp, void *userhandle)
{
        return(interrupted);
}


/* put lp on list */

static int create_handle(lprec *lp0, char *err)
{
        int i;

        if (lp0 == NULL)
        	ErrMsgTxt(err);
	for (i = 0; (i <= lp_last) && (lp[i] != NULL); i++);
	if (i > lp_last) {
	  	i = ++lp_last;
                if ((i % LPSTEP) == 0) {
                        if (i == 0)
                                lp = (lprec **) malloc(LPSTEP * sizeof(*lp));
                        else
                                lp = (lprec **) realloc(lp, (i + LPSTEP) * sizeof(*lp));
                        memset(lp + i, 0, LPSTEP * sizeof(*lp));
                }
        }
        lp[i] = lp0;

        putlogfunc(lp0, mylog, NULL);
        set_outputfile(lp0, "");
        putabortfunc(lp0, myabort, NULL);

        return(i);
}

/* lp validation test */

static int handle_valid(int handle)
{
        if (handle < 0 || handle > lp_last || lp[handle] == NULL)
	   return(0);
        else
           return(1);
}

/* free lp from list */

static void delete_handle(int handle)
{
        if (handle_valid(handle)) {
                char *name;
                lprec *lp0 = lp[handle];

                name = get_lp_name(lp0);
                if ((name != NULL) && (*name) && (strcmp(name, "Unnamed")))
                        drophash(name, NULL, handlehash);
                delete_lp(lp0);
                lp[handle] = NULL;
        }
}


static void set_handlename(char *name, int h)
{
        if (*name) {
        	if (handlehash == NULL)
        		handlehash = create_hash_table(100, 0);
                else {
                        char *oldname;

                        oldname = get_lp_name(lp0);
                        if ((oldname != NULL) && (*oldname) && (strcmp(oldname, "Unnamed")))
                                drophash(oldname, NULL, handlehash);

                }
                if (findhash(name, handlehash) == NULL)
                	puthash(name, h, NULL, handlehash);
        }
}


/* An exit function to clean up the lp structure.
   called on exit */

void ExitFcn()
{
	int	i;

        if (initialized) {
        	for (i = 0; i <= lp_last; i++)
                        delete_handle(i);
                free_hash_table(cmdhash);
                if (handlehash != NULL)
                	free_hash_table(handlehash);
                free(cmd);
        	free(errmsg);
                exit_lpsolve_lib();
#               if defined DEBUG
		        Printf("Terminated\n");
#               endif
        }
}


#if defined DEBUG || defined DEMO
/* xxlpsolve('demo') */
/* This demo is not a necessary part of the program */

static void impl_demo()
{
        int i, h;

	/* Printf("%.15f", GetRealScalar(prhs, 1)); */

        Check_nrhs(cmd, 0, nrhs);

        h = create_handle(make_lp(0, 4), "make_lp failed");
        lp0 = lp[h];

        Printf("min: 2 C1 +3 C2 -2 C3 +3 C4;\n");
	str_set_obj_fn(lp0, "2 3 -2 3");

        Printf("3 C1 +2 C2 +2 C3 +1 C4 <= 4.0;\n");
        str_add_constraint(lp0, "3 2 2 1", LE, 4.0);

        Printf("0 C1 +4 C2 +3 C3 +1 C4 >= 3.0;\n");
	str_add_constraint(lp0, "0 4 3 1", GE, 3.0);

        Printf("solve: %d\n", solve(lp0));
        Printf("obj: %f\n", get_objective(lp0));
        for (i = 1; i <= get_Ncolumns(lp0); i++)
        	Printf("C%d: %f\n", i, get_var_primalresult(lp0, get_Nrows(lp0) + i));
        write_lp(lp0, "a.lp");

        delete_handle(h);
}
#endif


/* since always a matrix vector is provided to both add_column and add_columnex, always call the more
   performant sparse version of the two routines */
/* return = xxlpsolve('add_column', lp, [column]) */
/* return = xxlpsolve('add_columnex', lp, [column]) */

static void impl_add_column()
{
        int m, count;
        int	*index;
	REAL	*vec;

        Check_nrhs(cmd, 2, nrhs);
        m = get_Nrows(lp0);
        vec = (REAL *) matCalloc(1 + m, sizeof(*vec));
        index = (int *) matCalloc(1 + m, sizeof(*index));
	count = GetRealSparseVector(prhs, 2, vec, index, 0, 1 + m, 0);
	result = add_columnex(lp0, count, vec, index);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        matFree(index);
        matFree(vec);
}


/* since always a matrix vector is provided to both add_constraint and add_constraintex, always call the more
   performant sparse version of the two routines */
/* return = xxlpsolve('add_constraint', lp, [row], constr_type, rh) */
/* return = xxlpsolve('add_constraintex', lp, [row], constr_type, rh) */

static void impl_add_constraint()
{
        int type, n, count;
        int	*index;
	REAL	*vec, value;

        Check_nrhs(cmd, 4, nrhs);
	type = (int) GetRealScalar(prhs, 3);
	if ((type != LE) && (type != EQ) && (type != GE))
		ErrMsgTxt("invalid constraint type.");
        value = GetRealScalar(prhs, 4);
        n = get_Ncolumns(lp0);
        vec = (REAL *) matCalloc(n, sizeof(*vec));
        index = (int *) matCalloc(n, sizeof(*index));
	count = GetRealSparseVector(prhs, 2, vec, index, 1, n, 0);
	result = add_constraintex(lp0, count, vec, index, type, value);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        matFree(index);
        matFree(vec);
}


/* return = xxlpsolve('add_SOS', lp, name, sostype, priority, [sosvars], [weights]) */

static void impl_add_SOS()
{
        int n, *sosvars;
        REAL *weights;
        int count1, count2;

        Check_nrhs(cmd, 6, nrhs);
        GetString(prhs, 2, buf, bufsz, TRUE);
        n = get_Ncolumns(lp0);
        sosvars = (int *) matCalloc(n, sizeof(*sosvars));
        weights = (REAL *) matCalloc(n, sizeof(*weights));
	count1 = GetIntVector(prhs, 5, sosvars, 0, n, FALSE);
        count2 = GetRealVector(prhs, 6, weights, 0, n, FALSE);
        if (count1 != count2) {
          matFree(weights);
          matFree(sosvars);
          ErrMsgTxt("add_SOS: sosvars and weights vector must have same size.");
        }
	result = add_SOS(lp0, buf, (int) GetRealScalar(prhs, 3), (int) GetRealScalar(prhs, 4), count1, sosvars,weights);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        matFree(weights);
        matFree(sosvars);
}


/* return = xxlpsolve('column_in_lp', lp, [column]) */

static void impl_column_in_lp()
{
        int n;
	REAL	*vec;

        Check_nrhs(cmd, 2, nrhs);
        n = get_Nrows(lp0);
        vec = (REAL *) matCalloc(1 + n, sizeof(*vec));
        GetRealVector(prhs, 2, vec, 0, 1 + n, TRUE);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = column_in_lp(lp0, vec);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        matFree(vec);
}


/* xxlpsolve('copy_lp', lp) */

static void impl_copy_lp()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = create_handle(copy_lp(lp0), "copy_lp failed");
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('default_basis', lp) */

static void impl_default_basis()
{
        Check_nrhs(cmd, 1, nrhs);
        default_basis(lp0);
}


/* return = xxlpsolve('del_column', lp, column) */

static void impl_del_column()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = del_column(lp0, (int) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('del_constraint', lp, del_row) */

static void impl_del_constraint()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = del_constraint(lp0, (int) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('delete_lp', lp) */
/* xxlpsolve('free_lp', lp) */

static void impl_delete_lp()
{
        Check_nrhs(cmd, 1, nrhs);
        delete_handle(h);
}


/* xxlpsolve('dualize_lp', lp) */

static void impl_dualize_lp()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = dualize_lp(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_anti_degen', lp) */
static void impl_get_anti_degen()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_anti_degen(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* [bascolumn] = xxlpsolve('get_basis', lp {, nonbasic}) */

static void impl_get_basis()
{
        int n, i, *bascolumn, *bascolumn0;
        MYBOOL nonbasic;

        if (nrhs == 1+1)
                n = 1;
        else
                n = 2;
        Check_nrhs(cmd, n, nrhs);
        if (n == 1)
                nonbasic = 0;
        else
        	nonbasic = (MYBOOL) GetRealScalar(prhs, 2);
        n = get_Nrows(lp0) + ((nonbasic) ? get_Ncolumns(lp0) : 0);
        bascolumn0 = bascolumn = (int *) matCalloc(1 + n, sizeof(*bascolumn));
        if(get_basis(lp0, bascolumn, nonbasic)) {
                ipr0 = ipr = CreateLongMatrix(n, 1, plhs, 0);
                for (i = 0; i < n; i++)
                  *(ipr++) = *(++bascolumn);
        }
        else {
                n = 0;
                ipr0 = ipr = CreateLongMatrix(n, 1, plhs, 0);
        }
        SetLongMatrix(ipr0, n, 1, plhs, 0, TRUE);
        matFree(bascolumn0);
}


/* return = xxlpsolve('get_basiscrash', lp) */

static void impl_get_basiscrash()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_basiscrash(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_bb_depthlimit', lp) */

static void impl_get_bb_depthlimit()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_bb_depthlimit(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_bb_floorfirst', lp) */

static void impl_get_bb_floorfirst()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_bb_floorfirst(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_bb_rule', lp) */

static void impl_get_bb_rule()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_bb_rule(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_bounds_tighter', lp) */

static void impl_get_bounds_tighter()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_bounds_tighter(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_break_at_value', lp) */

static void impl_get_break_at_value()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_break_at_value(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* name = xxlpsolve('get_col_name', lp, column) */
/* [names] = xxlpsolve('get_col_name', lp) */

static void impl_get_col_name()
{
        char *name;

        if (nrhs == 1+1) {
                int n, i;
                char **names;

                n = get_Ncolumns(lp0);

                names = (char **) matCalloc(n, sizeof(*names));
                for (i = 0; i < n; i++)
                        names[i] = strdup(get_col_name(lp0, i + 1));
                CreateString(names, n, plhs, 0);
                for (i = 0; i < n; i++)
                        free(names[i]);
                matFree(names);
        }
        else {
        	Check_nrhs(cmd, 2, nrhs);
                name = get_col_name(lp0, (int) GetRealScalar(prhs, 2));
                CreateString(&name, 1, plhs, 0);
        }
}


/* [column, return] = xxlpsolve('get_column', lp, col_nr) */
/* [column, return] = xxlpsolve('get_columnex', lp, col_nr) */

static void impl_get_column()
{
        int col;

        Check_nrhs(cmd, 2, nrhs);
	col = (int) GetRealScalar(prhs, 2);
        dpr = CreateDoubleMatrix(1 + get_Nrows(lp0), 1, plhs, 0);
	result = get_column(lp0, col, dpr);
        SetDoubleMatrix(dpr, 1 + get_Nrows(lp0), 1, plhs, 0, TRUE);
        if (nlhs > 1) {
                ipr = CreateLongMatrix(1, 1, plhs, 1);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 1, TRUE);
        }
}


/* return = xxlpsolve('get_constr_type', lp, row) */
/* [constr_type] = xxlpsolve('get_constr_type', lp) */

static void impl_get_constr_type()
{
	if (nrhs == 1+1) {
                int i, m;

                Check_nrhs(cmd, 1, nrhs);
                m = get_Nrows(lp0);
                ipr0 = ipr = CreateLongMatrix(m, 1, plhs, 0);
                for (i = 1; i <= m; i++)
                        *(ipr++) = get_constr_type(lp0, i);
                SetLongMatrix(ipr0, m, 1, plhs, 0, TRUE);
        }
        else {
        	Check_nrhs(cmd, 2, nrhs);
                ipr = CreateLongMatrix(1, 1, plhs, 0);
        	*ipr = get_constr_type(lp0, (int) GetRealScalar(prhs, 2));
                SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('get_constr_value', lp, row, [primsolution]) */

static void impl_get_constr_value()
{
        int n, count;
        int	*nzindex;
	REAL	*primsolution;

        if (nrhs == 1+2) {
                Check_nrhs(cmd, 2, nrhs);
                primsolution = NULL;
                nzindex = NULL;
                count = 0;
        }
        else {
                Check_nrhs(cmd, 3, nrhs);
                n = get_Ncolumns(lp0);
                if (n == 0)
                        n = 1;
                primsolution = (REAL *) matCalloc(n, sizeof(*primsolution));
                nzindex = (int *) matCalloc(n, sizeof(*nzindex));
        	count = GetRealSparseVector(prhs, 3, primsolution, nzindex, 1, n, 0);
        }
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_constr_value(lp0, (int) GetRealScalar(prhs, 2), count, primsolution, nzindex);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
        if (nzindex != NULL)
                matFree(nzindex);
        if (primsolution != NULL)
                matFree(primsolution);
}


/* [constr, return] = xxlpsolve('get_constraints', lp) */

static void impl_get_constraints()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(get_Nrows(lp0), 1, plhs, 0);
	result = get_constraints(lp0, dpr);
        SetDoubleMatrix(dpr, get_Nrows(lp0), 1, plhs, 0, TRUE);
        if (nlhs > 1) {
                ipr = CreateLongMatrix(1, 1, plhs, 1);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 1, TRUE);
        }
}


/* [duals, return] = xxlpsolve('get_dual_solution', lp) */

static void impl_get_dual_solution()
{
        int i;
	REAL	*vec = NULL;

        Check_nrhs(cmd, 1, nrhs);
	result = get_ptr_dual_solution(lp0, &vec);
        if ((!result) || (vec == NULL))
                ErrMsgTxt("get_dual_solution: sensitivity unknown.");
        i = get_Nrows(lp0) + get_Ncolumns(lp0);
        dpr = CreateDoubleMatrix(i, 1, plhs, 0);
        memcpy(dpr, vec + 1, i * sizeof(*vec));
        SetDoubleMatrix(dpr, i, 1, plhs, 0, TRUE);
        if (nlhs > 1) {
                ipr = CreateLongMatrix(1, 1, plhs, 1);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 1, TRUE);
        }
}


/* return = xxlpsolve('get_epsb', lp) */

static void impl_get_epsb()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_epsb(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_epsd', lp) */

static void impl_get_epsd()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_epsd(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_epsel', lp) */

static void impl_get_epsel()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_epsel(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_epsint', lp) */

static void impl_get_epsint()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_epsint(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_epsperturb', lp) */

static void impl_get_epsperturb()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_epsperturb(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_epspivot', lp) */

static void impl_get_epspivot()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_epspivot(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_improve', lp) */

static void impl_get_improve()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_improve(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_infinite', lp) */

static void impl_get_infinite()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_infinite(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_lowbo', lp, column) */
/* [return] = xxlpsolve('get_lowbo', lp) */

static void impl_get_lowbo()
{
	if (nrhs == 1+1) {
                int n, i;

                Check_nrhs(cmd, 1, nrhs);
                n = get_Ncolumns(lp0);
                dpr0 = dpr = CreateDoubleMatrix(n, 1, plhs, 0);
                for (i = 1; i <= n; i++)
                  *(dpr++) = get_lowbo(lp0, i);
                SetDoubleMatrix(dpr0, n, 1, plhs, 0, TRUE);
        }
        else {
        	Check_nrhs(cmd, 2, nrhs);
                dpr = CreateDoubleMatrix(1, 1, plhs, 0);
		*dpr = get_lowbo(lp0, (int) GetRealScalar(prhs, 2));
                SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('get_lp_index', lp, orig_index) */

static void impl_get_lp_index()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_lp_index(lp0, (int) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* name = xxlpsolve('get_lp_name', lp) */

static void impl_get_lp_name()
{
        char *name;

        Check_nrhs(cmd, 1, nrhs);
        name = get_lp_name(lp0);
        CreateString(&name, 1, plhs, 0);
}


/* value = xxlpsolve('get_mat', lp, row, col) */
/* [matrix, return] = xxlpsolve('get_mat', lp[, sparse]) */

static void impl_get_mat()
{
	if ((nrhs == 1+1) || (nrhs == 1+2)) {
                int m, n, i, sparse = FALSE;
        	REAL	*vec;

                Check_nrhs(cmd, nrhs - 1, nrhs);
                m = get_Nrows(lp0);
                n = get_Ncolumns(lp0);
                vec = (REAL *) matCalloc(1 + m, sizeof(*vec));
                result = TRUE;
                if (nrhs == 1+2)
                        sparse = (int) GetRealScalar(prhs, 2);
                if (!sparse) { /* return the matrix in dense format */
                        dpr0 = dpr = CreateDoubleMatrix(m, n, plhs, 0);
                        for (i = 1; (i <= n) && (result); i++) {
                		result = get_column(lp0, i, vec);
                        	memcpy(dpr, vec + 1, m * sizeof(*vec));
                                dpr += m;
                        }
                }
                else {         /* return the matrix in sparse format */
                        int nz = 0;

                        dpr0 = dpr = CreateDoubleSparseMatrix(m, n, plhs, 0);
                        for (i = 1; (i <= n) && (result); i++) {
                		result = get_column(lp0, i, vec);
                                SetColumnDoubleSparseMatrix(plhs, 0, m, n, dpr, i, vec + 1, NULL, m, &nz);
                        }
                }
                SetDoubleMatrix(dpr0, m, n, plhs, 0, TRUE);
                matFree(vec);
                if (nlhs > 1) {
                        ipr = CreateLongMatrix(1, 1, plhs, 1);
                        *ipr = result;
                        SetLongMatrix(ipr, 1, 1, plhs, 1, TRUE);
                }
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
                dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        	*dpr = get_mat(lp0, (int) GetRealScalar(prhs, 2), (int) GetRealScalar(prhs, 3));
                SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('get_max_level', lp) */

static void impl_get_max_level()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_max_level(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_maxpivot', lp) */

static void impl_get_maxpivot()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_maxpivot(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_mip_gap', lp, absolute) */

static void impl_get_mip_gap()
{
        Check_nrhs(cmd, 2, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_mip_gap(lp0, (MYBOOL) GetRealScalar(prhs, 2));
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_nameindex', lp, name, isrow) */

static void impl_get_nameindex()
{
        Check_nrhs(cmd, 3, nrhs);
        GetString(prhs, 2, buf, bufsz, TRUE);
        result = get_nameindex(lp0, buf, (MYBOOL) GetRealScalar(prhs, 3));
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_Ncolumns', lp) */

static void impl_get_Ncolumns()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_Ncolumns(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_negrange', lp) */

static void impl_get_negrange()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_negrange(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_nonzeros', lp) */

static void impl_get_nonzeros()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_nonzeros(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_Norig_columns', lp) */

static void impl_get_Norig_columns()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_Norig_columns(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_Norig_rows', lp) */

static void impl_get_Norig_rows()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_Norig_rows(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_Nrows', lp) */

static void impl_get_Nrows()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_Nrows(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_obj_bound', lp) */

static void impl_get_obj_bound()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_obj_bound(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* [row_vec, return] = xxlpsolve('get_obj_fn', lp) */
/* [row_vec, return] = xxlpsolve('get_obj_fun', lp) */

static void impl_get_obj_fn()
{
        int n;
	REAL	*vec;

        Check_nrhs(cmd, 1, nrhs);
        n = get_Ncolumns(lp0);
        dpr = CreateDoubleMatrix(1, n, plhs, 0);
        vec = (REAL *) matCalloc(1 + n, sizeof(*vec));
	result = get_row(lp0, 0, vec);
        memcpy(dpr, vec + 1, n * sizeof(*vec));
        SetDoubleMatrix(dpr, 1, n, plhs, 0, TRUE);
        matFree(vec);
        if (nlhs > 1) {
                ipr = CreateLongMatrix(1, 1, plhs, 1);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 1, TRUE);
        }
}


/* return = xxlpsolve('get_objective', lp) */

static void impl_get_objective()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_objective(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* name = xxlpsolve('get_objective_name', lp) */

static void impl_get_objective_name()
{
        char *name;

        Check_nrhs(cmd, 1, nrhs);
        name = get_row_name(lp0, 0);
        CreateString(&name, 1, plhs, 0);
}


/* return = xxlpsolve('get_orig_index', lp, lp_index) */

static void impl_get_orig_index()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_orig_index(lp0, (int) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* name = xxlpsolve('get_origcol_name', lp, column) */
/* [names] = xxlpsolve('get_origcol_name', lp) */

static void impl_get_origcol_name()
{
        char *name;

	if (nrhs == 1+1) {
                int n, i;
        	char **names;

                n = get_Ncolumns(lp0);

                names = (char **) matCalloc(n, sizeof(*names));
                for (i = 0; i < n; i++)
                        names[i] = strdup(get_origcol_name(lp0, i + 1));
                CreateString(names, n, plhs, 0);
                for (i = 0; i < n; i++)
                        free(names[i]);
                matFree(names);
        }
        else {
	        Check_nrhs(cmd, 2, nrhs);
                name = get_origcol_name(lp0, (int) GetRealScalar(prhs, 2));
                CreateString(&name, 1, plhs, 0);
        }
}


/* name = xxlpsolve('get_origrow_name', lp, row) */
/* [names] = xxlpsolve('get_origrow_name', lp) */

static void impl_get_origrow_name()
{
        char *name;

	if (nrhs == 1+1) {
                int m, i;
        	char **names;

                m = get_Nrows(lp0);

                names = (char **) matCalloc(m, sizeof(*names));
                for (i = 0; i < m; i++)
                        names[i] = strdup(get_origrow_name(lp0, i + 1));
                CreateString(names, m, plhs, 0);
                for (i = 0; i < m; i++)
                        free(names[i]);
                matFree(names);
        }
        else {
        	Check_nrhs(cmd, 2, nrhs);
                name = get_origrow_name(lp0, (int) GetRealScalar(prhs, 2));
                CreateString(&name, 1, plhs, 0);
        }
}


/* return = xxlpsolve('get_pivoting', lp) */

static void impl_get_pivoting()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_pivoting(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_presolve', lp) */

static void impl_get_presolve()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_presolve(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_presolveloops', lp) */

static void impl_get_presolveloops()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_presolveloops(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* [pv, return] = xxlpsolve('get_primal_solution', lp) */

static void impl_get_primal_solution()
{
        int i;

        Check_nrhs(cmd, 1, nrhs);
        i = 1 + get_Nrows(lp0) + get_Ncolumns(lp0);
        dpr = CreateDoubleMatrix(i, 1, plhs, 0);
	result = get_primal_solution(lp0, dpr);
        SetDoubleMatrix(dpr, i, 1, plhs, 0, TRUE);
        if (nlhs > 1) {
                ipr = CreateLongMatrix(1, 1, plhs, 1);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 1, TRUE);
        }
}


/* return = xxlpsolve('get_print_sol', lp) */

static void impl_get_print_sol()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_print_sol(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_rh', lp, row) */
/* [rh] = xxlpsolve('get_rh', lp) */

static void impl_get_rh()
{
	if (nrhs == 1+1) {
                int m, i;

                Check_nrhs(cmd, 1, nrhs);
                m = get_Nrows(lp0);
                dpr0 = dpr = CreateDoubleMatrix(1 + m, 1, plhs, 0);
                for (i = 0; i <= m; i++)
                  *(dpr++) = get_rh(lp0, i);
                SetDoubleMatrix(dpr0, 1 + m, 1, plhs, 0, TRUE);
        }
        else {
                Check_nrhs(cmd, 2, nrhs);
                dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        	*dpr = get_rh(lp0, (int) GetRealScalar(prhs, 2));
                SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('get_rh_range', lp, row) */
/* [rh_ranges] = xxlpsolve('get_rh_range', lp) */

static void impl_get_rh_range()
{
	if (nrhs == 1+1) {
                int m, i;

                Check_nrhs(cmd, 1, nrhs);
                m = get_Nrows(lp0);
                dpr0 = dpr = CreateDoubleMatrix(m, 1, plhs, 0);
                for (i = 1; i <= m; i++)
                  *(dpr++) = get_rh_range(lp0, i);
                SetDoubleMatrix(dpr0, m, 1, plhs, 0, TRUE);
        }
        else {
	        Check_nrhs(cmd, 2, nrhs);
                dpr = CreateDoubleMatrix(1, 1, plhs, 0);
		*dpr = get_rh_range(lp0, (int) GetRealScalar(prhs, 2));
                SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
        }
}


/* [row, return] = xxlpsolve('get_row', lp, row_nr) */
/* [row, return] = xxlpsolve('get_rowex', lp, row_nr) */

static void impl_get_row()
{
        int n;
	REAL	*vec;

        Check_nrhs(cmd, 2, nrhs);
        n = get_Ncolumns(lp0);
        dpr = CreateDoubleMatrix(1, n, plhs, 0);
        vec = (REAL *) matCalloc(1 + n, sizeof(*vec));
	result = get_row(lp0, (int) GetRealScalar(prhs, 2), vec);
        memcpy(dpr, vec + 1, n * sizeof(*vec));
        SetDoubleMatrix(dpr, 1, n, plhs, 0, TRUE);
        matFree(vec);
        if (nlhs > 1) {
                ipr = CreateLongMatrix(1, 1, plhs, 1);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 1, TRUE);
        }
}


/* name = xxlpsolve('get_row_name', lp, row) */
/* [names] = xxlpsolve('get_row_name', lp) */

static void impl_get_row_name()
{
        char *name;

        if (nrhs == 1+1) {
                int m, i;
        	char **names;

                m = get_Nrows(lp0);

                names = (char **) matCalloc(m, sizeof(*names));
                for (i = 0; i < m; i++)
                        names[i] = strdup(get_row_name(lp0, i + 1));
                CreateString(names, m, plhs, 0);
                for (i = 0; i < m; i++)
                        free(names[i]);
                matFree(names);
        }
        else {
                Check_nrhs(cmd, 2, nrhs);
                name = get_row_name(lp0, (int) GetRealScalar(prhs, 2));
                CreateString(&name, 1, plhs, 0);
        }
}


/* return = xxlpsolve('get_scalelimit', lp) */

static void impl_get_scalelimit()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_scalelimit(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_scaling', lp) */

static void impl_get_scaling()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_scaling(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* [objfrom, objtill, objfromvalue, objtillvalue, return] = xxlpsolve('get_sensitivity_obj', lp) */
/* [objfrom, objtill, objfromvalue, objtillvalue, return] = xxlpsolve('get_sensitivity_objex', lp) */

static void impl_get_sensitivity_objex()
{
        int n;
	REAL	*objfrom = NULL, *objtill = NULL, *objfromvalue, *objtillvalue;

        Check_nrhs(cmd, 1, nrhs);
	result = get_ptr_sensitivity_obj(lp0, &objfrom, &objtill);
        if ((!result) || (objfrom == NULL) || (objtill == NULL))
                ErrMsgTxt("get_sensitivity_obj: sensitivity unknown.");
        n = get_Ncolumns(lp0);
        objfrom = CreateDoubleMatrix(1, n, plhs, 0);
        if (nlhs > 1)
        	objtill = CreateDoubleMatrix(1, n, plhs, 1);
        else
                objtill = NULL;
        if (nlhs > 2)
        	objfromvalue = CreateDoubleMatrix(1, n, plhs, 2);
        else
                objfromvalue = NULL;
        if (nlhs > 3) {
        	objtillvalue = CreateDoubleMatrix(1, n, plhs, 3);
                memset(objtillvalue, 0, n * sizeof(*objtillvalue));
        }
        else
                objtillvalue = NULL;
	result = get_sensitivity_objex(lp0, objfrom, objtill, objfromvalue, NULL /* objtillvalue */);
        SetDoubleMatrix(objfrom, 1, n, plhs, 0, TRUE);
        SetDoubleMatrix(objtill, 1, n, plhs, 1, TRUE);
        SetDoubleMatrix(objfromvalue, 1, n, plhs, 2, TRUE);
        SetDoubleMatrix(objtillvalue, 1, n, plhs, 3, TRUE);
        if (nlhs > 4) {
                ipr = CreateLongMatrix(1, 1, plhs, 4);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 4, TRUE);
        }
}


/* [duals, dualsfrom, dualstill, return] = xxlpsolve('get_sensitivity_rhs', lp) */
/* [duals, dualsfrom, dualstill, return] = xxlpsolve('get_sensitivity_rhsex', lp) */

static void impl_get_sensitivity_rhsex()
{
        int i;
        REAL *duals = NULL, *dualsfrom = NULL, *dualstill = NULL;

        Check_nrhs(cmd, 1, nrhs);
        result = get_ptr_sensitivity_rhs(lp0, &duals, &dualsfrom, &dualstill);
        if ((!result) || (duals == NULL) || (dualsfrom == NULL) || (dualstill == NULL))
                ErrMsgTxt("get_sensitivity_rhs: sensitivity unknown.");
        i = get_Nrows(lp0) + get_Ncolumns(lp0);
        duals = CreateDoubleMatrix(i, 1, plhs, 0);
        if (nlhs > 1)
        	dualsfrom = CreateDoubleMatrix(i, 1, plhs, 1);
        else
		dualsfrom = NULL;
        if (nlhs > 2)
        	dualstill = CreateDoubleMatrix(i, 1, plhs, 2);
        else
                dualstill = NULL;
	result = get_sensitivity_rhs(lp0, duals, dualsfrom, dualstill);
        SetDoubleMatrix(duals, i, 1, plhs, 0, TRUE);
        SetDoubleMatrix(dualsfrom, i, 1, plhs, 1, TRUE);
        SetDoubleMatrix(dualstill, i, 1, plhs, 2, TRUE);
        if (nlhs > 3) {
                ipr = CreateLongMatrix(1, 1, plhs, 3);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 3, TRUE);
        }
}


/* return = xxlpsolve('get_simplextype', lp) */

static void impl_get_simplextype()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_simplextype(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* [obj, x, duals, return] = xxlpsolve('get_solution', lp) */

static void impl_get_solution()
{
        int m, n;
        REAL *duals;

        Check_nrhs(cmd, 1, nrhs);

        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_objective(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);

        if (nlhs > 1) {
	        n = get_Ncolumns(lp0);
        	dpr = CreateDoubleMatrix(n, 1, plhs, 1);
	        result = get_variables(lp0, dpr);
                SetDoubleMatrix(dpr, n, 1, plhs, 1, TRUE);
        }

        if (nlhs > 2) {
                m = get_Nrows(lp0);
        	dpr = CreateDoubleMatrix(m, 1, plhs, 2);
		result &= get_ptr_dual_solution(lp0, &duals);
                memcpy(dpr, duals + 1, m * sizeof(*dpr));
                SetDoubleMatrix(dpr, m, 1, plhs, 2, TRUE);
        }

        if (nlhs > 3) {
                ipr = CreateLongMatrix(1, 1, plhs, 3);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 3, TRUE);
        }
}


/* return = xxlpsolve('get_solutioncount', lp) */

static void impl_get_solutioncount()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_solutioncount(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_solutionlimit', lp) */

static void impl_get_solutionlimit()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_solutionlimit(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_status', lp) */

static void impl_get_status()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_status(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_statustext', lp, statuscode) */

static void impl_get_statustext()
{
        char *name;

        Check_nrhs(cmd, 2, nrhs);
        name = get_statustext(lp0, (int) GetRealScalar(prhs, 2));
	CreateString(&name, 1, plhs, 0);
}


/* return = xxlpsolve('get_timeout', lp) */

static void impl_get_timeout()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = get_timeout(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_total_iter', lp) */

static void impl_get_total_iter()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = (double) get_total_iter(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_total_nodes', lp) */

static void impl_get_total_nodes()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = (double) get_total_nodes(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_upbo', lp, column) */
/* [upbo] = xxlpsolve('get_upbo', lp) */

static void impl_get_upbo()
{
	if (nrhs == 1+1) {
                int n, i;

                Check_nrhs(cmd, 1, nrhs);
                n = get_Ncolumns(lp0);
                dpr0 = dpr = CreateDoubleMatrix(n, 1, plhs, 0);
                for (i = 1; i <= n; i++)
                  *(dpr++) = get_upbo(lp0, i);
                SetDoubleMatrix(dpr0, n, 1, plhs, 0, TRUE);
        }
        else {
	        Check_nrhs(cmd, 2, nrhs);
                dpr = CreateDoubleMatrix(1, 1, plhs, 0);
		*dpr = get_upbo(lp0, (int) GetRealScalar(prhs, 2));
                SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('get_var_branch', lp, column) */
/* [var_branch] = xxlpsolve('get_var_branch', lp) */

static void impl_get_var_branch()
{
	if (nrhs == 1+1) {
                int n, i;

                Check_nrhs(cmd, 1, nrhs);
                n = get_Ncolumns(lp0);
                ipr0 = ipr = CreateLongMatrix(n, 1, plhs, 0);
                for (i = 1; i <= n; i++)
                  *(ipr++) = get_var_branch(lp0, i);
                SetLongMatrix(ipr0, n, 1, plhs, 0, TRUE);
        }
        else {
                Check_nrhs(cmd, 2, nrhs);
                ipr = CreateLongMatrix(1, 1, plhs, 0);
        	*ipr = get_var_branch(lp0, (int) GetRealScalar(prhs, 2));
                SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('get_var_dualresult', lp, index) */

static void impl_get_var_dualresult()
{
        Check_nrhs(cmd, 2, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
	*dpr = get_var_dualresult(lp0, (int) GetRealScalar(prhs, 2));
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_var_primalresult', lp, index) */

static void impl_get_var_primalresult()
{
        Check_nrhs(cmd, 2, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
	*dpr = get_var_primalresult(lp0, (int) GetRealScalar(prhs, 2));
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_var_priority', lp, column) */
/* [var_priority] = xxlpsolve('get_var_priority', lp) */

static void impl_get_var_priority()
{
	if (nrhs == 1+1) {
                int n, i;

                Check_nrhs(cmd, 1, nrhs);
                n = get_Ncolumns(lp0);
                ipr0 = ipr = CreateLongMatrix(n, 1, plhs, 0);
                for (i = 1; i <= n; i++)
                  *(ipr++) = get_var_priority(lp0, i);
                SetLongMatrix(ipr0, n, 1, plhs, 0, TRUE);
        }
        else {
        	Check_nrhs(cmd, 2, nrhs);
                ipr = CreateLongMatrix(1, 1, plhs, 0);
		*ipr = get_var_priority(lp0, (int) GetRealScalar(prhs, 2));
                SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        }
}


/* [var, return] = xxlpsolve('get_variables', lp) */

static void impl_get_variables()
{
        int n;

        Check_nrhs(cmd, 1, nrhs);
        n = get_Ncolumns(lp0);
        dpr = CreateDoubleMatrix(n, 1, plhs, 0);
        result = get_variables(lp0, dpr);
        SetDoubleMatrix(dpr, n, 1, plhs, 0, TRUE);
        if (nlhs > 1) {
                ipr = CreateLongMatrix(1, 1, plhs, 1);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 1, TRUE);
        }
}


/* return = xxlpsolve('get_verbose', lp) */

static void impl_get_verbose()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = get_verbose(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('get_working_objective', lp) */

static void impl_get_working_objective()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = get_working_objective(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* [basisvector, return] = xxlpsolve('guess_basis', lp, [guessvector]) */

static void impl_guess_basis()
{
        int i, n, m;
        REAL *guessvector;
        int *basisvector, *basisvector0;

        Check_nrhs(cmd, 2, nrhs);

        n = get_Ncolumns(lp0);
        m = get_Nrows(lp0);
        guessvector = (REAL *) matCalloc(1 + n, sizeof(REAL));
        basisvector0 = basisvector = (int *) matCalloc(1 + n + m, sizeof(int));
	GetRealVector(prhs, 2, guessvector, 1, n, TRUE);
        result = guess_basis(lp0, guessvector, basisvector);
        ipr0 = ipr = CreateLongMatrix(n + m, 1, plhs, 0);
        for (i = 0; i < n + m; i++)
          *(ipr++) = *(++basisvector);
        SetLongMatrix(ipr0, n + m, 1, plhs, 0, TRUE);
        matFree(basisvector0);
        matFree(guessvector);
        if (nlhs > 1) {
                ipr = CreateLongMatrix(1, 1, plhs, 1);
                *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 1, TRUE);
        }
}


/* return = xxlpsolve('has_BFP', lp) */

static void impl_has_BFP()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = has_BFP(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('has_XLI', lp) */

static void impl_has_XLI()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = has_XLI(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_add_rowmode', lp) */

static void impl_is_add_rowmode()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = is_add_rowmode(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_anti_degen', lp, testmask) */

static void impl_is_anti_degen()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = is_anti_degen(lp0, (int) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_binary', lp, column) */
/* [binary] = xxlpsolve('is_binary', lp) */

static void impl_is_binary()
{
        if (nrhs == 1+1) {
                int n, i;

                Check_nrhs(cmd, 1, nrhs);
                n = get_Ncolumns(lp0);
                ipr0 = ipr = CreateLongMatrix(n, 1, plhs, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_binary(lp0, i);
                SetLongMatrix(ipr0, n, 1, plhs, 0, TRUE);
        }
        else {
        	Check_nrhs(cmd, 2, nrhs);
                ipr = CreateLongMatrix(1, 1, plhs, 0);
		*ipr = is_binary(lp0, (int) GetRealScalar(prhs, 2));
                SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('is_break_at_first', lp) */

static void impl_is_break_at_first()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = is_break_at_first(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_constr_type', lp, row, mask) */

static void impl_is_constr_type()
{
        Check_nrhs(cmd, 3, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = is_constr_type(lp0, (int) GetRealScalar(prhs, 2), (int) GetRealScalar(prhs, 3));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_debug', lp) */

static void impl_is_debug()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = is_debug(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_feasible', lp, [values] {, threshold}) */

static void impl_is_feasible()
{
        int i, n;
	REAL	*vec, threshold;

        if (nrhs == 2+1)
                n = 2;
        else
                n = 3;
        Check_nrhs(cmd, n, nrhs);
        i = get_Nrows(lp0) + get_Ncolumns(lp0);
        vec = (REAL *) matCalloc(1 + i, sizeof(REAL));
	GetRealVector(prhs, 2, vec, 1, i, TRUE);
        if (n == 2)
                threshold = get_epsint(lp0);
        else
                threshold = GetRealScalar(prhs, 3);
        result = is_feasible(lp0, vec, threshold);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
	matFree(vec);
}


/* return = xxlpsolve('is_free', lp, column) */
/* return = xxlpsolve('is_unbounded', lp, column) */
/* [free] = xxlpsolve('is_free', lp) */
/* [free] = xxlpsolve('is_unbounded', lp) */

static void impl_is_free()
{
        if (nrhs == 1+1) {
                int n, i;

                Check_nrhs(cmd, 1, nrhs);
                n = get_Ncolumns(lp0);
                ipr0 = ipr = CreateLongMatrix(n, 1, plhs, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_unbounded(lp0, i);
                SetLongMatrix(ipr0, n, 1, plhs, 0, TRUE);
        }
        else {
        	Check_nrhs(cmd, 2, nrhs);
                ipr = CreateLongMatrix(1, 1, plhs, 0);
		*ipr = is_unbounded(lp0, (int) GetRealScalar(prhs, 2));
                SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('is_infinite', lp, value) */

static void impl_is_infinite()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = is_infinite(lp0, GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_int', lp, column) */
/* [int] = xxlpsolve('is_int', lp) */

static void impl_is_int()
{
        if (nrhs == 1+1) {
                int n, i;

                Check_nrhs(cmd, 1, nrhs);
                n = get_Ncolumns(lp0);
                ipr0 = ipr = CreateLongMatrix(n, 1, plhs, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_int(lp0, i);
                SetLongMatrix(ipr0, n, 1, plhs, 0, TRUE);
        }
        else {
                Check_nrhs(cmd, 2, nrhs);
                ipr = CreateLongMatrix(1, 1, plhs, 0);
        	*ipr = is_int(lp0, (int) GetRealScalar(prhs, 2));
                SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('is_integerscaling', lp) */

static void impl_is_integerscaling()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = is_integerscaling(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_maxim', lp) */

static void impl_is_maxim()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = is_maxim(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_nativeBFP', lp) */

static void impl_is_nativeBFP()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = is_nativeBFP(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_nativeXLI', lp) */

static void impl_is_nativeXLI()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = is_nativeXLI(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_negative', lp, column) */
/* [negative] = xxlpsolve('is_negative', lp) */

static void impl_is_negative()
{
        if (nrhs == 1+1) {
                int n, i;

                Check_nrhs(cmd, 1, nrhs);
                n = get_Ncolumns(lp0);
                ipr0 = ipr = CreateLongMatrix(n, 1, plhs, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_negative(lp0, i);
                SetLongMatrix(ipr0, n, 1, plhs, 0, TRUE);
        }
        else {
	        Check_nrhs(cmd, 2, nrhs);
                ipr = CreateLongMatrix(1, 1, plhs, 0);
		*ipr = is_negative(lp0, (int) GetRealScalar(prhs, 2));
                SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('is_piv_mode', lp, testmask) */

static void impl_is_piv_mode()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = is_piv_mode(lp0, (int) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_piv_rule', lp, rule) */

static void impl_is_piv_rule()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = is_piv_rule(lp0, (int) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_presolve', lp, testmask) */

static void impl_is_presolve()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = is_presolve(lp0, (int) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_scalemode', lp, testmask) */

static void impl_is_scalemode()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = is_scalemode(lp0, (int) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_scaletype', lp, scaletype) */

static void impl_is_scaletype()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = is_scaletype(lp0, (int) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_semicont', lp, column) */
/* [semicont] = xxlpsolve('is_semicont', lp) */

static void impl_is_semicont()
{
        if (nrhs == 1+1) {
                int n, i;

                Check_nrhs(cmd, 1, nrhs);
                n = get_Ncolumns(lp0);
                ipr0 = ipr = CreateLongMatrix(n, 1, plhs, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_semicont(lp0, i);
                SetLongMatrix(ipr0, n, 1, plhs, 0, TRUE);
        }
        else {
        	Check_nrhs(cmd, 2, nrhs);
                ipr = CreateLongMatrix(1, 1, plhs, 0);
		*ipr = is_semicont(lp0, (int) GetRealScalar(prhs, 2));
                SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('is_SOS_var', lp, column) */
/* [SOS_var] = xxlpsolve('is_SOS_var', lp) */

static void impl_is_SOS_var()
{
        if (nrhs == 1+1) {
                int n, i;

                Check_nrhs(cmd, 1, nrhs);
                n = get_Ncolumns(lp0);
                ipr0 = ipr = CreateLongMatrix(n, 1, plhs, 0);
                for (i = 1; i <= n; i++)
        		*(ipr++) = is_SOS_var(lp0, i);
                SetLongMatrix(ipr0, n, 1, plhs, 0, TRUE);
        }
        else {
        	Check_nrhs(cmd, 2, nrhs);
                ipr = CreateLongMatrix(1, 1, plhs, 0);
		*ipr = is_SOS_var(lp0, (int) GetRealScalar(prhs, 2));
                SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        }
}


/* return = xxlpsolve('is_trace', lp) */

static void impl_is_trace()
{
        Check_nrhs(cmd, 1, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = is_nativeXLI(lp0);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('is_use_names', lp, isrow) */

static void impl_is_use_names()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = (Long) is_use_names(lp0, (MYBOOL) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* versionstring = xxlpsolve('lp_solve_version') */

static void impl_lp_solve_version()
{
        int majorversion, minorversion, release, build;

        Check_nrhs(cmd, 0, nrhs);
        lp_solve_version(&majorversion, &minorversion, &release, &build);
        sprintf(buf, "%d.%d.%d.%d", majorversion, minorversion, release, build);
	CreateString(&buf, 1, plhs, 0);
}


/* lp = xxlpsolve('make_lp', rows, columns) */

static void impl_make_lp()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = create_handle(make_lp((int) GetRealScalar(prhs, 1), (int) GetRealScalar(prhs, 2)), "make_lp failed");
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* lp = xxlpsolve('resize_lp', lp, rows, columns) */

static void impl_resize_lp()
{
        Check_nrhs(cmd, 3, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = resize_lp(lp0, (int) GetRealScalar(prhs, 2), (int) GetRealScalar(prhs, 3));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('print_constraints', lp {, columns}) */

static void impl_print_constraints()
{
        int n, columns;

        if (nrhs == 1+1)
        	n = 1;
        else
        	n = 2;
        Check_nrhs(cmd, n, nrhs);
        if (n == 1)
                columns = 1;
        else
                columns = (int) GetRealScalar(prhs, 2);
	print_constraints(lp0, columns);
}


/* return = xxlpsolve('print_debugdump', lp, filename) */

static void impl_print_debugdump()
{
        char filename[260];

        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, filename, sizeof(filename), TRUE);
        result = print_debugdump(lp0, filename);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('print_duals', lp) */

static void impl_print_duals()
{
        Check_nrhs(cmd, 1, nrhs);
	print_duals(lp0);
}


/* xxlpsolve('print_lp', lp) */

static void impl_print_lp()
{
        Check_nrhs(cmd, 1, nrhs);
	print_lp(lp0);
}


/* xxlpsolve('print_objective', lp) */

static void impl_print_objective()
{
        Check_nrhs(cmd, 1, nrhs);
	print_objective(lp0);
}


/* xxlpsolve('print_scales', lp) */

static void impl_print_scales()
{
        Check_nrhs(cmd, 1, nrhs);
	print_scales(lp0);
}


/* xxlpsolve('print_solution', lp {, columns}) */

static void impl_print_solution()
{
        int n, columns;

        if (nrhs == 1+1)
        	n = 1;
        else
        	n = 2;
        Check_nrhs(cmd, n, nrhs);
        if (n == 1)
                columns = 1;
        else
                columns = (int) GetRealScalar(prhs, 2);
	print_solution(lp0, columns);
}


/* xxlpsolve('print_str', lp, str) */

static void impl_print_str()
{
        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, buf, bufsz, TRUE);
        print_str(lp0, buf);
}


/* xxlpsolve('print_tableau', lp) */

static void impl_print_tableau()
{
        Check_nrhs(cmd, 1, nrhs);
	print_tableau(lp0);
}


/* [handle_vec] = xxlpsolve('print_handle') */
/* print all used handles */

static void impl_print_handle()
{
        int i, j, k;

	j = 0;
        for (i = 0; i <= lp_last; i++)
	  if (lp[i] != NULL)
	     j++;
        if (j)
          k = 1;
        else
          k = 0;
        ipr0 = ipr = CreateLongMatrix(j, k, plhs, 0);
	for (i = 0; i <= lp_last; i++)
	  if (lp[i] != NULL)
            *(ipr++) = i;
        SetLongMatrix(ipr0, j, k, plhs, 0, TRUE);
}


/* [handle_vec] = xxlpsolve('get_handle', 'name') */
/* get handle from model name */

static void impl_get_handle()
{
        hashelem *hp;

        Check_nrhs(cmd, 1, nrhs);
        GetString(prhs, 1, buf, bufsz, TRUE);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        if (handlehash != NULL)
        	hp = findhash(buf, handlehash);
        else
                hp = NULL;
        if (hp == NULL)
                *ipr = -1;
        else
        	*ipr = hp->index;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* [return, info] = xxlpsolve('read_basis', lp, filename) */

static void impl_read_basis()
{
        char filename[260];
#       define info filename

        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, filename, sizeof(filename), TRUE);
        result = read_basis(lp0, filename, (nlhs > 1) ? info : NULL);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        if (nlhs > 1) {
                char *ptr = info;
                CreateString(&ptr, 1, plhs, 1);
        }
#       undef info
}


/* lp = xxlpsolve('read_freeMPS', filename {, verbose}) */

static void impl_read_freeMPS()
{
        int n, verbose;
        char filename[260];

        if (nrhs == 1+1)
        	n = 1;
        else
        	n = 2;
        Check_nrhs(cmd, n, nrhs);
        if (n >= 2)
        	verbose = (int) GetRealScalar(prhs, 2);
        else
                verbose = NORMAL;
	GetString(prhs, 1, filename, sizeof(filename), TRUE);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = create_handle(read_freeMPS(filename, verbose), "read_freeMPS can't read file.");
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* lp = xxlpsolve('read_lp_file', filename {, verbose {, lp_name}}) */
/* lp = xxlpsolve('read_lp', filename {, verbose {, lp_name}}) */
/* lp = xxlpsolve('read_LP', filename {, verbose {, lp_name}}) */

static void impl_read_LP()
{
        int n, verbose;
        char filename[260], lp_name[50];

        if (nrhs == 1+1)
        	n = 1;
        else if (nrhs == 1+2)
        	n = 2;
        else
        	n = 3;
        Check_nrhs(cmd, n, nrhs);
	GetString(prhs, 1, filename, sizeof(filename), TRUE);
        if (n >= 2)
        	verbose = (int) GetRealScalar(prhs, 2);
        else
                verbose = NORMAL;
        if (n >= 3)
        	GetString(prhs, 3, lp_name, sizeof(lp_name), TRUE);
        else
                *lp_name = 0;
        lp0 = read_LP(filename, verbose, lp_name);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = create_handle(lp0, "read_LP can't read file.");
        set_handlename(lp_name, (int) *ipr);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* lp = xxlpsolve('read_mps', filename {, verbose}) */
/* lp = xxlpsolve('read_MPS', filename {, verbose}) */

static void impl_read_MPS()
{
        int n, verbose;
        char filename[260];

        if (nrhs == 1+1)
        	n = 1;
        else
        	n = 2;
        Check_nrhs(cmd, n, nrhs);
	GetString(prhs, 1, filename, sizeof(filename), TRUE);
        if (n >= 2)
        	verbose = (int) GetRealScalar(prhs, 2);
        else
                verbose = NORMAL;
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = create_handle(read_MPS(filename, verbose), "read_MPS can't read file.");
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* lp = xxlpsolve('read_params', lp, filename {, options}) */

static void impl_read_params()
{
        int n;
        char filename[260], options[50];

        if (nrhs == 1+2)
        	n = 2;
        else
        	n = 3;
        Check_nrhs(cmd, n, nrhs);
	GetString(prhs, 2, filename, sizeof(filename), TRUE);
        if (n >= 3)
        	GetString(prhs, 3, options, sizeof(options), TRUE);
        else
                *options = 0;
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = read_params(lp0, filename, options);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* lp = xxlpsolve('read_XLI', xliname, modelname {, dataname {, options {, verbose}}} */

static void impl_read_XLI()
{
        int n, verbose;
        char xliname[260], modelname[260], dataname[260], options[260];

        if (nrhs == 1+2)
        	n = 2;
        else if (nrhs == 1+3)
        	n = 3;
        else if (nrhs == 1+4)
        	n = 4;
        else
        	n = 5;
        Check_nrhs(cmd, n, nrhs);
	GetString(prhs, 1, xliname, sizeof(xliname), TRUE);
        GetString(prhs, 2, modelname, sizeof(modelname), TRUE);
        if (n >= 3)
        	GetString(prhs, 3, dataname, sizeof(dataname), TRUE);
        else
                *dataname = 0;
        if (n >= 4)
        	GetString(prhs, 4, options, sizeof(options), TRUE);
        else
                *options = 0;
        if (n >= 5)
        	verbose = (int) GetRealScalar(prhs, 5);
        else
                verbose = NORMAL;
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = create_handle(read_XLI(xliname, modelname, (*dataname) ? dataname : NULL, options, verbose), "read_XLI can't read file.");
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('reset_params', lp) */

static void impl_reset_params()
{
        Check_nrhs(cmd, 1, nrhs);
        reset_params(lp0);
}


/* return = xxlpsolve('set_add_rowmode', lp, turnon) */

static void impl_set_add_rowmode()
{
        Check_nrhs(cmd, 2, nrhs);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = set_add_rowmode(lp0, (MYBOOL) GetRealScalar(prhs, 2));
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('set_anti_degen', lp, anti_degen) */

static void impl_set_anti_degen()
{
        Check_nrhs(cmd, 2, nrhs);
	set_anti_degen(lp0, (int) GetRealScalar(prhs, 2));
}


/* return = xxlpsolve('set_basis', lp, [bascolumn], nonbasic) */

static void impl_set_basis()
{
        int i, *bascolumn;
        MYBOOL nonbasic;

        Check_nrhs(cmd, 3, nrhs);
        nonbasic = (MYBOOL) GetRealScalar(prhs, 3);
        i = get_Nrows(lp0) + ((nonbasic) ? get_Ncolumns(lp0) : 0);
        bascolumn = (int *) matCalloc(1 + i, sizeof(*bascolumn));
	GetIntVector(prhs, 2, bascolumn, 1, i, TRUE);
	result = set_basis(lp0, bascolumn, nonbasic);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
	matFree(bascolumn);
}


/* xxlpsolve('set_basiscrash', lp, mode) */

static void impl_set_basiscrash()
{
        Check_nrhs(cmd, 2, nrhs);
	set_basiscrash(lp0, (int) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_basisvar', lp, basisPos, enteringCol) */

static void impl_set_basisvar()
{
        Check_nrhs(cmd, 3, nrhs);
	set_basisvar(lp0, (int) GetRealScalar(prhs, 2), (int) GetRealScalar(prhs, 3));
}


/* xxlpsolve('set_bb_depthlimit', lp, bb_maxlevel) */

static void impl_set_bb_depthlimit()
{
        Check_nrhs(cmd, 2, nrhs);
	set_bb_depthlimit(lp0, (int) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_bb_floorfirst', lp, bb_floorfirst) */

static void impl_set_bb_floorfirst()
{
        Check_nrhs(cmd, 2, nrhs);
	set_bb_floorfirst(lp0, (int) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_bb_rule', lp, bb_rule) */

static void impl_set_bb_rule()
{
        Check_nrhs(cmd, 2, nrhs);
	set_bb_rule(lp0, (int) GetRealScalar(prhs, 2));
}


/* return = xxlpsolve('set_BFP', lp, filename) */

static void impl_set_BFP()
{
        char filename[260];

        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, filename, sizeof(filename), TRUE);
        result = set_BFP(lp0, filename);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('set_binary', lp, column, must_be_bin) */
/* return = xxlpsolve('set_binary', lp, [must_be_bin]) */

static void impl_set_binary()
{
        if (nrhs == 1+2) {
                int i, n, *vec;

                Check_nrhs(cmd, 2, nrhs);
                n = get_Ncolumns(lp0);
                vec = (int *) matCalloc(n, sizeof(*vec));
        	GetIntVector(prhs, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_binary(lp0, i + 1, (MYBOOL) vec[i]);
                matFree(vec);
        }
        else {
                Check_nrhs(cmd, 3, nrhs);
        	result = set_binary(lp0, (int) GetRealScalar(prhs, 2), (MYBOOL) GetRealScalar(prhs, 3));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('set_bounds', lp, column, lower, upper) */
/* return = xxlpsolve('set_bounds', lp, [lower], [upper]) */

static void impl_set_bounds()
{
        if (nrhs == 1+3) {
                int i, n;
                REAL	*lower, *upper;

                Check_nrhs(cmd, 3, nrhs);
                n = get_Ncolumns(lp0);
                lower = (REAL *) matCalloc(n, sizeof(REAL));
                upper = (REAL *) matCalloc(n, sizeof(REAL));
        	GetRealVector(prhs, 2, lower, 0, n, TRUE);
                GetRealVector(prhs, 3, upper, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_bounds(lp0, i + 1, lower[i], upper[i]);
                matFree(upper);
                matFree(lower);
        }
        else {
	        Check_nrhs(cmd, 4, nrhs);
		result = set_bounds(lp0, (int) GetRealScalar(prhs, 2), GetRealScalar(prhs, 3), GetRealScalar(prhs, 4));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('set_bounds_tighter', lp, tighten) */

static void impl_set_bounds_tighter()
{
        Check_nrhs(cmd, 2, nrhs);
	set_bounds_tighter(lp0, (MYBOOL) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_break_at_first', lp, break_at_first) */

static void impl_set_break_at_first()
{
        Check_nrhs(cmd, 2, nrhs);
	set_break_at_first(lp0, (MYBOOL) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_break_at_value', lp, break_at_value) */

static void impl_set_break_at_value()
{
        Check_nrhs(cmd, 2, nrhs);
	set_break_at_value(lp0, GetRealScalar(prhs, 2));
}


/* return = xxlpsolve('set_col_name', lp, column, name) */
/* return = xxlpsolve('set_col_name', lp, [names]) */

static void impl_set_col_name()
{
        if (nrhs == 1+2) {
                int n, i;
                strArray pa;

                Check_nrhs(cmd, 2, nrhs);
                n = get_Ncolumns(lp0);
                pa = GetCellCharItems(prhs, 2, n);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++) {
                	GetCellString(pa, i, buf, bufsz);
                        result = set_col_name(lp0, i + 1, buf);
                }
                FreeCellCharItems(pa, n);
        }
        else {
                Check_nrhs(cmd, 3, nrhs);
                GetString(prhs, 3, buf, bufsz, TRUE);
                result = set_col_name(lp0, (int) GetRealScalar(prhs, 2), buf);
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* since always a matrix vector is provided to both set_column and set_columnex, always call the more
   performant sparse version of the two routines */
/* return = xxlpsolve('set_column', lp, col_no, [column]) */
/* return = xxlpsolve('set_columnex', lp, col_no, [column]) */

static void impl_set_column()
{
        int m, count;
        int	*index;
	REAL	*vec;

        Check_nrhs(cmd, 3, nrhs);
        m = get_Nrows(lp0);
        vec = (REAL *) matCalloc(1 + m, sizeof(*vec));
        index = (int *) matCalloc(1 + m, sizeof(*index));
	count = GetRealSparseVector(prhs, 3, vec, index, 0, 1 + m, 0);
	result = set_columnex(lp0, (int) GetRealScalar(prhs, 2), count, vec, index);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        matFree(index);
        matFree(vec);
}


/* return = xxlpsolve('set_constr_type', lp, row, con_type) */
/* return = xxlpsolve('set_constr_type', lp, [con_type]) */

static void impl_set_constr_type()
{
        if (nrhs == 1+2) {
                int i, m, *vec;

                Check_nrhs(cmd, 2, nrhs);
                m = get_Nrows(lp0);
                vec = (int *) matCalloc(m, sizeof(*vec));
        	GetIntVector(prhs, 2, vec, 0, m, TRUE);
                result = TRUE;
                for (i = 0; (i < m) && (result); i++)
                        result = set_constr_type(lp0, i + 1, vec[i]);
                matFree(vec);
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
		result = set_constr_type(lp0, (int) GetRealScalar(prhs, 2), (int) GetRealScalar(prhs, 3));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('set_debug', lp, debug) */

static void impl_set_debug()
{
        Check_nrhs(cmd, 2, nrhs);
	set_debug(lp0, (MYBOOL) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_epsb', lp, epsb) */

static void impl_set_epsb()
{
        Check_nrhs(cmd, 2, nrhs);
	set_epsb(lp0, GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_epsd', lp, epsd) */

static void impl_set_epsd()
{
        Check_nrhs(cmd, 2, nrhs);
	set_epsd(lp0, GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_epsel', lp, epsel) */

static void impl_set_epsel()
{
        Check_nrhs(cmd, 2, nrhs);
	set_epsel(lp0, GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_epsint', lp, epsint) */

static void impl_set_epsint()
{
        Check_nrhs(cmd, 2, nrhs);
	set_epsint(lp0, GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_epslevel', lp, epslevel) */

static void impl_set_epslevel()
{
        Check_nrhs(cmd, 2, nrhs);
	set_epslevel(lp0, (int) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_epsperturb', lp, epsperturb) */

static void impl_set_epsperturb()
{
        Check_nrhs(cmd, 2, nrhs);
	set_epsperturb(lp0, GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_epspivot', lp, epspivot) */

static void impl_set_epspivot()
{
        Check_nrhs(cmd, 2, nrhs);
	set_epspivot(lp0, GetRealScalar(prhs, 2));
}


/* return = xxlpsolve('set_free', lp, column) */
/* return = xxlpsolve('set_unbounded', lp, column) */

static void impl_set_free()
{
        Check_nrhs(cmd, 2, nrhs);
	result = set_unbounded(lp0, (int) GetRealScalar(prhs, 2));
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('set_improve', lp, improve) */

static void impl_set_improve()
{
        Check_nrhs(cmd, 2, nrhs);
	set_improve(lp0, (int) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_infinite', lp, infinite) */

static void impl_set_infinite()
{
        Check_nrhs(cmd, 2, nrhs);
	set_infinite(lp0, GetRealScalar(prhs, 2));
}


/* return = xxlpsolve('set_int', lp, column, must_be_int) */
/* return = xxlpsolve('set_int', lp, [must_be_int]) */

static void impl_set_int()
{
        if (nrhs == 1+2) {
                int i, n, *vec;

                Check_nrhs(cmd, 2, nrhs);
                n = get_Ncolumns(lp0);
                vec = (int *) matCalloc(n, sizeof(*vec));
        	GetIntVector(prhs, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_int(lp0, i + 1, (MYBOOL) vec[i]);
                matFree(vec);
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
		result = set_int(lp0, (int) GetRealScalar(prhs, 2), (MYBOOL) GetRealScalar(prhs, 3));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('set_lowbo', lp, column, value) */
/* return = xxlpsolve('set_lowbo', lp, [values]) */

static void impl_set_lowbo()
{
        if (nrhs == 1+2) {
                int i, n;
                REAL *vec;

                Check_nrhs(cmd, 2, nrhs);
                n = get_Ncolumns(lp0);
                vec = (REAL *) matCalloc(n, sizeof(*vec));
        	GetRealVector(prhs, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_lowbo(lp0, i + 1, vec[i]);
                matFree(vec);
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
		result = set_lowbo(lp0, (int) GetRealScalar(prhs, 2), GetRealScalar(prhs, 3));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('set_lp_name', lp, name) */

static void impl_set_lp_name()
{
        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, buf, bufsz, TRUE);
        set_handlename(buf, h);
        result = set_lp_name(lp0, buf);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('set_mat', lp, [matrix]) */
/* return = xxlpsolve('set_mat', lp, row, column, value) */

static void impl_set_mat()
{
        if (nrhs == 1+2) {
                int m, n, j, *index, *index1, count;
                REAL *obj, *obj1, *vec, *vec1, a;
                rMatrix mat;

                mat = GetpMatrix(prhs, 2);
		/* Called with a matrix argument */
		m = GetM(mat);
		n = GetN(mat);
		if ((get_Nrows(lp0) != m) || (get_Ncolumns(lp0) != n))
			ErrMsgTxt("Invalid matrix dimension.");

                obj = obj1 = (REAL *) matCalloc(1 + n, sizeof(*obj));
                result = get_row(lp0, 0, obj);
                vec = (REAL *) matCalloc(1 + m, sizeof(*vec));
                index = (int *) matCalloc(1 + m, sizeof(*index));
                for (j = 1; (j <= n) && (result); j++) {
                        vec1 = vec;
                        index1 = index;
                        count = 0;
                        if ((a = (*(++obj1))) != 0.0) {
                                *(vec1++) = a;
                                *(index1++) = 0;
                                count ++;
                        }
			count += GetRealSparseVector2(prhs, mat, 2, vec1, index1, 1, m, j);
                        result = set_columnex(lp0, j, count, vec, index);
                }
                matFree(index);
                matFree(vec);
                matFree(obj);
                Check_nrhs(cmd, 2, nrhs);
        }
        else { /* called with a single matrix element */
        	Check_nrhs(cmd, 4, nrhs);
		result = set_mat(lp0, (int) GetRealScalar(prhs, 2), (int) GetRealScalar(prhs, 3), GetRealScalar(prhs, 4));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('set_maxim', lp) */

static void impl_set_maxim()
{
        Check_nrhs(cmd, 1, nrhs);
        set_maxim(lp0);
}


/* xxlpsolve('set_maxpivot', max_num_inv) */

static void impl_set_maxpivot()
{
        Check_nrhs(cmd, 2, nrhs);
	set_maxpivot(lp0, (int) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_minim', lp) */

static void impl_set_minim()
{
        Check_nrhs(cmd, 1, nrhs);
        set_minim(lp0);
}


/* xxlpsolve('set_mip_gap', lp, absolute, mip_gap) */

static void impl_set_mip_gap()
{
        Check_nrhs(cmd, 3, nrhs);
	set_mip_gap(lp0, (MYBOOL) GetRealScalar(prhs, 2), GetRealScalar(prhs, 3));
}


/* xxlpsolve('set_negrange', negrange) */

static void impl_set_negrange()
{
        Check_nrhs(cmd, 2, nrhs);
	set_negrange(lp0, GetRealScalar(prhs, 2));
}


/* return = xxlpsolve('set_obj', lp, column, value) */
/* return = xxlpsolve('set_obj', lp, [values]) */

static void impl_set_obj()
{
        if (nrhs == 1+2) {
                impl_set_obj_fn();
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
		result = set_obj(lp0, (int) GetRealScalar(prhs, 2), GetRealScalar(prhs, 3));
                ipr = CreateLongMatrix(1, 1, plhs, 0);
	        *ipr = result;
                SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        }
}


/* xxlpsolve('set_obj_bound', lp, obj_bound) */

static void impl_set_obj_bound()
{
        Check_nrhs(cmd, 2, nrhs);
	set_obj_bound(lp0, GetRealScalar(prhs, 2));
}


/* since always a matrix vector is provided to both set_obj_fn and set_obj_fnex, always call the more
   performant sparse version of the two routines */
/* return = xxlpsolve('set_obj_fn', lp, [row]) */
/* return = xxlpsolve('set_obj_fnex', lp, [row]) */

static void impl_set_obj_fn()
{
        int n, count;
        int	*index;
	REAL	*vec;

        Check_nrhs(cmd, 2, nrhs);
        n = get_Ncolumns(lp0);
        vec = (REAL *) matCalloc(1 + n, sizeof(*vec));
        index = (int *) matCalloc(1 + n, sizeof(*index));
	count = GetRealSparseVector(prhs, 2, vec, index, 1, n, 0);
	result = set_obj_fnex(lp0, count, vec, index);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        matFree(index);
        matFree(vec);
}


/* return = xxlpsolve('set_outputfile', lp, filename) */

static void impl_set_outputfile()
{
        char filename[260];

        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, filename, sizeof(filename), TRUE);
        result = set_outputfile(lp0, (*filename) ? filename : NULL);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('set_pivoting', lp, pivoting) */

static void impl_set_pivoting()
{
        Check_nrhs(cmd, 2, nrhs);
	set_pivoting(lp0, (int) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_preferdual', lp, dodual) */

static void impl_set_preferdual()
{
        Check_nrhs(cmd, 2, nrhs);
	set_preferdual(lp0, (MYBOOL) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_presolve', lp, do_presolve {, maxloops}) */

static void impl_set_presolve()
{
        int maxloops;

        if (nrhs == 1+2) {
                Check_nrhs(cmd, 2, nrhs);
                maxloops = get_presolveloops(lp0);
        }
        else {
                Check_nrhs(cmd, 3, nrhs);
                maxloops = (int) GetRealScalar(prhs, 3);
        }
	set_presolve(lp0, (int) GetRealScalar(prhs, 2), maxloops);
}


/* xxlpsolve('set_print_sol', lp, print_sol) */

static void impl_set_print_sol()
{
        Check_nrhs(cmd, 2, nrhs);
	set_print_sol(lp0, (int) GetRealScalar(prhs, 2));
}


/* return = xxlpsolve('set_rh', lp, row, value) */
/* return = xxlpsolve('set_rh', lp, [values]) */

static void impl_set_rh()
{
        if (nrhs == 1+2) {
                int i, m;
                REAL *vec;

                Check_nrhs(cmd, 2, nrhs);
                m = get_Nrows(lp0);
                vec = (REAL *) matCalloc(1 + m, sizeof(*vec));
        	GetRealVector(prhs, 2, vec, 0, 1 + m, TRUE);
                result = TRUE;
                for (i = 0; (i <= m) && (result); i++)
                        result = set_rh(lp0, i, vec[i]);
                matFree(vec);
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
		result = set_rh(lp0, (int) GetRealScalar(prhs, 2), GetRealScalar(prhs, 3));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('set_rh_range', lp, row, deltavalue) */
/* return = xxlpsolve('set_rh_range', lp, [deltavalues]) */

static void impl_set_rh_range()
{
        if (nrhs == 1+2) {
                int i, m;
                REAL *vec;

                Check_nrhs(cmd, 2, nrhs);
                m = get_Nrows(lp0);
                vec = (REAL *) matCalloc(1 + m, sizeof(*vec));
        	GetRealVector(prhs, 2, vec, 0, 1 + m, TRUE);
                result = TRUE;
                for (i = 0; (i < m) && (result); i++)
                        result = set_rh_range(lp0, i + 1, vec[i]);
                matFree(vec);
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
		result = set_rh_range(lp0, (int) GetRealScalar(prhs, 2), GetRealScalar(prhs, 3));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('set_rh_vec', lp, [rh]) */

static void impl_set_rh_vec()
{
        int m;
	REAL	*vec;

        Check_nrhs(cmd, 2, nrhs);
        m = get_Nrows(lp0);
	vec = (REAL *) matCalloc(1 + m, sizeof(REAL));
	GetRealVector(prhs, 2, vec, 1, m, TRUE);
	set_rh_vec(lp0, vec);
        matFree(vec);
}


/* since always a matrix vector is provided to both set_row and set_rowex, always call the more
   performant sparse version of the two routines */
/* return = xxlpsolve('set_row', lp, row_no, [row]) */
/* return = xxlpsolve('set_rowex', lp, row_no, [row]) */

static void impl_set_row()
{
        int n, count;
        int	*index;
	REAL	*vec;

        Check_nrhs(cmd, 3, nrhs);
        n = get_Ncolumns(lp0);
        vec = (REAL *) matCalloc(1 + n, sizeof(*vec));
        index = (int *) matCalloc(1 + n, sizeof(*index));
	count = GetRealSparseVector(prhs, 3, vec, index, 1, n, 0);
	result = set_rowex(lp0, (int) GetRealScalar(prhs, 2), count, vec, index);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        matFree(index);
        matFree(vec);
}


/* return = xxlpsolve('set_row_name', lp, row, name) */
/* return = xxlpsolve('set_row_name', lp, [names]) */

static void impl_set_row_name()
{
        if (nrhs == 1+2) {
                int m, i;
                strArray pa;

                Check_nrhs(cmd, 2, nrhs);
                m = get_Nrows(lp0);
                pa = GetCellCharItems(prhs, 2, m);
                result = TRUE;
                for (i = 0; (i < m) && (result); i++) {
                	GetCellString(pa, i, buf, bufsz);
                        result = set_row_name(lp0, i + 1, buf);
                }
                FreeCellCharItems(pa, m);
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
	        GetString(prhs, 3, buf, bufsz, TRUE);
	        result = set_row_name(lp0, (int) GetRealScalar(prhs, 2), buf);
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('set_scalelimit', lp, scalelimit) */

static void impl_set_scalelimit()
{
        Check_nrhs(cmd, 2, nrhs);
	set_scalelimit(lp0, GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_scaling', lp, scalemode) */

static void impl_set_scaling()
{
        Check_nrhs(cmd, 2, nrhs);
	set_scaling(lp0, (int) GetRealScalar(prhs, 2));
}


/* return = xxlpsolve('set_semicont', lp, column, must_be_sc) */
/* return = xxlpsolve('set_semicont', lp, [must_be_sc]) */

static void impl_set_semicont()
{
        if (nrhs == 1+2) {
                int i, n, *vec;

                Check_nrhs(cmd, 2, nrhs);
                n = get_Ncolumns(lp0);
                vec = (int *) matCalloc(n, sizeof(*vec));
        	GetIntVector(prhs, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_semicont(lp0, i + 1, (MYBOOL) vec[i]);
                matFree(vec);
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
		result = set_semicont(lp0, (int) GetRealScalar(prhs, 2), (MYBOOL) GetRealScalar(prhs, 3));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('set_sense', lp, maximize) */

static void impl_set_sense()
{
        Check_nrhs(cmd, 2, nrhs);
	set_sense(lp0, (MYBOOL) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_simplextype', lp, simplextype) */

static void impl_set_simplextype()
{
        Check_nrhs(cmd, 2, nrhs);
	set_simplextype(lp0, (int) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_solutionlimit', lp, simplextype) */

static void impl_set_solutionlimit()
{
        Check_nrhs(cmd, 2, nrhs);
	set_solutionlimit(lp0, (int) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_timeout', lp, sectimeout) */

static void impl_set_timeout()
{
        Check_nrhs(cmd, 2, nrhs);
	set_timeout(lp0, (long) GetRealScalar(prhs, 2));
}


/* xxlpsolve('set_trace', lp, trace) */

static void impl_set_trace()
{
        Check_nrhs(cmd, 2, nrhs);
	set_trace(lp0, (MYBOOL) GetRealScalar(prhs, 2));
}


/* return = xxlpsolve('set_upbo', lp, column, value) */
/* return = xxlpsolve('set_upbo', lp, [values]) */

static void impl_set_upbo()
{
        if (nrhs == 1+2) {
                int i, n;
                REAL *vec;

                Check_nrhs(cmd, 2, nrhs);
                n = get_Ncolumns(lp0);
                vec = (REAL *) matCalloc(n, sizeof(*vec));
        	GetRealVector(prhs, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_upbo(lp0, i + 1, vec[i]);
                matFree(vec);
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
		result = set_upbo(lp0, (int) GetRealScalar(prhs, 2), GetRealScalar(prhs, 3));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('set_use_names', lp, isrow, use_names) */

static void impl_set_use_names()
{
        Check_nrhs(cmd, 3, nrhs);
        set_use_names(lp0, (MYBOOL) GetRealScalar(prhs, 2), (MYBOOL) GetRealScalar(prhs, 3));
}


/* return = xxlpsolve('set_var_branch', lp, column, branch_mode) */
/* return = xxlpsolve('set_var_branch', lp, [branch_mode]) */

static void impl_set_var_branch()
{
        if (nrhs == 1+2) {
                int i, n, *vec;

                Check_nrhs(cmd, 2, nrhs);
                n = get_Ncolumns(lp0);
                vec = (int *) matCalloc(n, sizeof(*vec));
        	GetIntVector(prhs, 2, vec, 0, n, TRUE);
                result = TRUE;
                for (i = 0; (i < n) && (result); i++)
                        result = set_var_branch(lp0, i + 1, vec[i]);
                matFree(vec);
        }
        else {
	        Check_nrhs(cmd, 3, nrhs);
		result = set_var_branch(lp0, (int) GetRealScalar(prhs, 2), (int) GetRealScalar(prhs, 3));
        }
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('set_var_weights', lp, [weights]) */

static void impl_set_var_weights()
{
        int n;
	REAL	*vec;

        Check_nrhs(cmd, 2, nrhs);
        n = get_Ncolumns(lp0);
	vec = (REAL *) matCalloc(n, sizeof(REAL));
	GetRealVector(prhs, 2, vec, 0, n, TRUE);
	result = set_var_weights(lp0, vec);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        matFree(vec);
}


/* xxlpsolve('set_verbose', lp, verbose) */

static void impl_set_verbose()
{
        Check_nrhs(cmd, 2, nrhs);
	set_verbose(lp0, (int) GetRealScalar(prhs, 2));
}


/* return = xxlpsolve('set_XLI', lp, filename) */

static void impl_set_XLI()
{
        char filename[260];

        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, filename, sizeof(filename), TRUE);
        result = set_XLI(lp0, filename);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* result = xxlpsolve('solve', lp) */

static void impl_solve()
{
        /* int m, n, i */;

        Check_nrhs(cmd, 1, nrhs);
	result = solve(lp0);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
        switch (result) {
        case OPTIMAL:
        case SUBOPTIMAL:
        case PROCBREAK:
        case FEASFOUND:
/*
          if (get_verbose(lp0) >= DETAILED) {
	    Printf("Branch & Bound depth: %d\n", get_max_level(lp0));
	    Printf("Nodes processed: %.0f\n", (double) get_total_nodes(lp0));
	    Printf("Simplex pivots: %.0f\n", (double) get_total_iter(lp0));
            Printf("Number of equal solutions: %d\n\n", get_solutioncount(lp0));
          }
          if (get_verbose(lp0) >= NORMAL) {
            Printf("Value of objective function: %f\n\n", get_objective(lp0));
            Printf("Actual values of the variables:\n");
            m = get_Nrows(lp0);
            n = get_Ncolumns(lp0);
            for (i = 1; i <= n; i++)
              Printf("%s: %f\n", get_col_name(lp0, i), get_var_primalresult(lp0, m + i));
          }
*/
          break;
        case NOMEMORY:
          if (get_verbose(lp0) >= NORMAL)
          	Printf("Out of memory\n");
          break;
        case INFEASIBLE:
          if (get_verbose(lp0) >= NORMAL)
	  	Printf("This problem is infeasible\n");
          break;
        case UNBOUNDED:
          if (get_verbose(lp0) >= NORMAL)
          	Printf("This problem is unbounded\n");
          break;
        case PROCFAIL:
          if (get_verbose(lp0) >= NORMAL)
          	Printf("The B&B routine failed\n");
          break;
        case TIMEOUT:
          if (get_verbose(lp0) >= NORMAL)
          	Printf("Timeout\n");
          break;
        case USERABORT:
          if (get_verbose(lp0) >= NORMAL)
          	Printf("User aborted\n");
          break;
        case DEGENERATE:
          if (get_verbose(lp0) >= NORMAL)
          	Printf("This problem is degenerative\n");
          break;
        case NUMFAILURE:
          if (get_verbose(lp0) >= NORMAL)
          	Printf("Numerical failure encountered\n");
          break;
        case NOFEASFOUND:
          if (get_verbose(lp0) >= NORMAL)
          	Printf("No feasible branch and bound solution found\n");
          break;
        default:
          if (get_verbose(lp0) >= NORMAL)
          	Printf("lp_solve failed\n");
          break;
        }
}


/* return = xxlpsolve('time_elapsed', lp) */

static void impl_time_elapsed()
{
        Check_nrhs(cmd, 1, nrhs);
        dpr = CreateDoubleMatrix(1, 1, plhs, 0);
        *dpr = time_elapsed(lp0);
        SetDoubleMatrix(dpr, 1, 1, plhs, 0, TRUE);
}


/* xxlpsolve('unscale', lp) */

static void impl_unscale()
{
        Check_nrhs(cmd, 1, nrhs);
        unscale(lp0);
}


/* return = xxlpsolve('write_basis', lp, filename) */

static void impl_write_basis()
{
        char filename[260];

        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, filename, sizeof(filename), TRUE);
        result = write_basis(lp0, filename);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('write_freemps', lp, filename) */
/* return = xxlpsolve('write_freeMPS', lp, filename) */

static void impl_write_freemps()
{
        char filename[260];

        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, filename, sizeof(filename), TRUE);
        result = write_freemps(lp0, filename);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('write_lp', lp, filename) */
/* return = xxlpsolve('write_LP', lp, filename) */

static void impl_write_lp()
{
        char filename[260];

        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, filename, sizeof(filename), TRUE);
        result = write_lp(lp0, filename);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('write_mps', lp, filename) */
/* return = xxlpsolve('write_MPS', lp, filename) */

static void impl_write_mps()
{
        char filename[260];

        Check_nrhs(cmd, 2, nrhs);
        GetString(prhs, 2, filename, sizeof(filename), TRUE);
        result = write_mps(lp0, filename);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* lp = xxlpsolve('write_params', lp, filename {, options}) */

static void impl_write_params()
{
        int n;
        char filename[260], options[50];

        if (nrhs == 1+2)
        	n = 2;
        else
        	n = 3;
        Check_nrhs(cmd, n, nrhs);
	GetString(prhs, 2, filename, sizeof(filename), TRUE);
        if (n >= 3)
        	GetString(prhs, 3, options, sizeof(options), TRUE);
        else
                *options = 0;
        ipr = CreateLongMatrix(1, 1, plhs, 0);
        *ipr = write_params(lp0, filename, options);
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


/* return = xxlpsolve('write_XLI', lp, filename {, options {, results}}) */

static void impl_write_XLI()
{
        char filename[260], options[50];
        int n;
        MYBOOL results;

        if (nrhs == 1+2)
                n = 2;
        else if (nrhs == 1+3)
                n = 3;
        else
                n = 4;

        Check_nrhs(cmd, n, nrhs);
        GetString(prhs, 2, filename, sizeof(filename), TRUE);
        if (n >= 3)
        	GetString(prhs, 3, options, sizeof(options), TRUE);
        else
                *options = 0;
        if (n >= 4)
                results = (MYBOOL) GetRealScalar(prhs, 4);
        else
                results = FALSE;
        result = write_XLI(lp0, filename, options, results);
        ipr = CreateLongMatrix(1, 1, plhs, 0);
	*ipr = result;
        SetLongMatrix(ipr, 1, 1, plhs, 0, TRUE);
}


static struct {
        char *cmd;
        impl_routine *routine;
        int needshandle;
} routines[] = {
  { "add_column", impl_add_column, TRUE },
  { "add_columnex", impl_add_column, TRUE },
  { "add_constraint", impl_add_constraint, TRUE },
  { "add_constraintex", impl_add_constraint, TRUE },
  { "add_SOS", impl_add_SOS, TRUE },
  { "column_in_lp", impl_column_in_lp, TRUE },
  { "copy_lp", impl_copy_lp, TRUE },
  { "default_basis", impl_default_basis, TRUE },
  { "del_column", impl_del_column, TRUE },
  { "del_constraint", impl_del_constraint, TRUE },
  { "delete_lp", impl_delete_lp, TRUE },
  { "dualize_lp", impl_dualize_lp, TRUE },
  { "free_lp", impl_delete_lp, TRUE },
  { "get_anti_degen", impl_get_anti_degen, TRUE },
  { "get_basis", impl_get_basis, TRUE },
  { "get_basiscrash", impl_get_basiscrash, TRUE },
  { "get_bb_depthlimit", impl_get_bb_depthlimit, TRUE },
  { "get_bb_floorfirst", impl_get_bb_floorfirst, TRUE },
  { "get_bb_rule", impl_get_bb_rule, TRUE },
  { "get_bounds_tighter", impl_get_bounds_tighter, TRUE },
  { "get_break_at_value", impl_get_break_at_value, TRUE },
  { "get_col_name", impl_get_col_name, TRUE },
  { "get_column", impl_get_column, TRUE },
  { "get_columnex", impl_get_column, TRUE },
  { "get_constr_type", impl_get_constr_type, TRUE },
  { "get_constr_value", impl_get_constr_value, TRUE },
  { "get_constraints", impl_get_constraints, TRUE },
  { "get_dual_solution", impl_get_dual_solution, TRUE },
  { "get_epsb", impl_get_epsb, TRUE },
  { "get_epsd", impl_get_epsd, TRUE },
  { "get_epsel", impl_get_epsel, TRUE },
  { "get_epsint", impl_get_epsint, TRUE },
  { "get_epsperturb", impl_get_epsperturb, TRUE },
  { "get_epspivot", impl_get_epspivot, TRUE },
  { "get_improve", impl_get_improve, TRUE },
  { "get_infinite", impl_get_infinite, TRUE },
  { "get_lowbo", impl_get_lowbo, TRUE },
  { "get_lp_index", impl_get_lp_index, TRUE },
  { "get_lp_name", impl_get_lp_name, TRUE },
  { "get_mat", impl_get_mat, TRUE },
  { "get_max_level", impl_get_max_level, TRUE },
  { "get_maxpivot", impl_get_maxpivot, TRUE },
  { "get_mip_gap", impl_get_mip_gap, TRUE },
  { "get_nameindex", impl_get_nameindex, TRUE },
  { "get_Ncolumns", impl_get_Ncolumns, TRUE },
  { "get_negrange", impl_get_negrange, TRUE },
  { "get_nonzeros", impl_get_nonzeros, TRUE },
  { "get_Norig_columns", impl_get_Norig_columns, TRUE },
  { "get_Norig_rows", impl_get_Norig_rows, TRUE },
  { "get_Nrows", impl_get_Nrows, TRUE },
  { "get_obj_bound", impl_get_obj_bound, TRUE },
  { "get_objective", impl_get_objective, TRUE },
  { "get_orig_index", impl_get_orig_index, TRUE },
  { "get_origcol_name", impl_get_origcol_name, TRUE },
  { "get_origrow_name", impl_get_origrow_name, TRUE },
  { "get_pivoting", impl_get_pivoting, TRUE },
  { "get_presolve", impl_get_presolve, TRUE },
  { "get_presolveloops", impl_get_presolveloops, TRUE },
  { "get_primal_solution", impl_get_primal_solution, TRUE },
  { "get_print_sol", impl_get_print_sol, TRUE },
  { "get_rh", impl_get_rh, TRUE },
  { "get_rh_range", impl_get_rh_range, TRUE },
  { "get_row", impl_get_row, TRUE },
  { "get_rowex", impl_get_row, TRUE },
  { "get_row_name", impl_get_row_name, TRUE },
  { "get_scalelimit", impl_get_scalelimit, TRUE },
  { "get_scaling", impl_get_scaling, TRUE },
  { "get_sensitivity_obj", impl_get_sensitivity_objex, TRUE },
  { "get_sensitivity_objex", impl_get_sensitivity_objex, TRUE },
  { "get_sensitivity_rhs", impl_get_sensitivity_rhsex, TRUE },
  { "get_sensitivity_rhsex", impl_get_sensitivity_rhsex, TRUE },
  { "get_simplextype", impl_get_simplextype, TRUE },
  { "get_solution", impl_get_solution, TRUE },
  { "get_solutioncount", impl_get_solutioncount, TRUE },
  { "get_solutionlimit", impl_get_solutionlimit, TRUE },
  { "get_status", impl_get_status, TRUE },
  { "get_statustext", impl_get_statustext, TRUE },
  { "get_timeout", impl_get_timeout, TRUE },
  { "get_total_iter", impl_get_total_iter, TRUE },
  { "get_total_nodes", impl_get_total_nodes, TRUE },
  { "get_upbo", impl_get_upbo, TRUE },
  { "get_var_branch", impl_get_var_branch, TRUE },
  { "get_var_dualresult", impl_get_var_dualresult, TRUE },
  { "get_var_primalresult", impl_get_var_primalresult, TRUE },
  { "get_var_priority", impl_get_var_priority, TRUE },
  { "get_variables", impl_get_variables, TRUE },
  { "get_verbose", impl_get_verbose, TRUE },
  { "get_working_objective", impl_get_working_objective, TRUE },
  { "guess_basis", impl_guess_basis, TRUE },
  { "has_BFP", impl_has_BFP, TRUE },
  { "has_XLI", impl_has_XLI, TRUE },
  { "is_add_rowmode", impl_is_add_rowmode, TRUE },
  { "is_anti_degen", impl_is_anti_degen, TRUE },
  { "is_binary", impl_is_binary, TRUE },
  { "is_break_at_first", impl_is_break_at_first, TRUE },
  { "is_constr_type", impl_is_constr_type, TRUE },
  { "is_debug", impl_is_debug, TRUE },
  { "is_feasible", impl_is_feasible, TRUE },
  { "is_free", impl_is_free, TRUE },
  { "is_infinite", impl_is_infinite, TRUE },
  { "is_int", impl_is_int, TRUE },
  { "is_integerscaling", impl_is_integerscaling, TRUE },
  { "is_maxim", impl_is_maxim, TRUE },
  { "is_nativeBFP", impl_is_nativeBFP, TRUE },
  { "is_nativeXLI", impl_is_nativeXLI, TRUE },
  { "is_negative", impl_is_negative, TRUE },
  { "is_piv_mode", impl_is_piv_mode, TRUE },
  { "is_piv_rule", impl_is_piv_rule, TRUE },
  { "is_presolve", impl_is_presolve, TRUE },
  { "is_scalemode", impl_is_scalemode, TRUE },
  { "is_scaletype", impl_is_scaletype, TRUE },
  { "is_semicont", impl_is_semicont, TRUE },
  { "is_SOS_var", impl_is_SOS_var, TRUE },
  { "is_trace", impl_is_trace, TRUE },
  { "is_unbounded", impl_is_free, TRUE },
  { "is_use_names", impl_is_use_names, TRUE },
  { "lp_solve_version", impl_lp_solve_version, FALSE },
  { "make_lp", impl_make_lp, FALSE },
  { "print_constraints", impl_print_constraints, TRUE },
  { "print_debugdump", impl_print_debugdump, TRUE },
  { "print_duals", impl_print_duals, TRUE },
  { "print_lp", impl_print_lp, TRUE },
  { "print_objective", impl_print_objective, TRUE },
  { "print_scales", impl_print_scales, TRUE },
  { "print_solution", impl_print_solution, TRUE },
  { "print_str", impl_print_str, TRUE },
  { "print_tableau", impl_print_tableau, TRUE },
  { "read_basis", impl_read_basis, TRUE },
  { "read_freemps", impl_read_freeMPS, FALSE },
  { "read_freeMPS", impl_read_freeMPS, FALSE },
  { "read_lp", impl_read_LP, FALSE },
  { "read_LP", impl_read_LP, FALSE },
  { "read_mps", impl_read_MPS, FALSE },
  { "read_MPS", impl_read_MPS, FALSE },
  { "read_params", impl_read_params, TRUE },
  { "read_XLI", impl_read_XLI, FALSE },
  { "reset_params", impl_reset_params, TRUE },
  { "resize_lp", impl_resize_lp, TRUE },
  { "set_add_rowmode", impl_set_add_rowmode, TRUE },
  { "set_anti_degen", impl_set_anti_degen, TRUE },
  { "set_basis", impl_set_basis, TRUE },
  { "set_basiscrash", impl_set_basiscrash, TRUE },
  { "set_basisvar", impl_set_basisvar, TRUE },
  { "set_bb_depthlimit", impl_set_bb_depthlimit, TRUE },
  { "set_bb_floorfirst", impl_set_bb_floorfirst, TRUE },
  { "set_bb_rule", impl_set_bb_rule, TRUE },
  { "set_BFP", impl_set_BFP, TRUE },
  { "set_binary", impl_set_binary, TRUE },
  { "set_bounds", impl_set_bounds, TRUE },
  { "set_bounds_tighter", impl_set_bounds_tighter, TRUE },
  { "set_break_at_first", impl_set_break_at_first, TRUE },
  { "set_break_at_value", impl_set_break_at_value, TRUE },
  { "set_col_name", impl_set_col_name, TRUE },
  { "set_column", impl_set_column, TRUE },
  { "set_columnex", impl_set_column, TRUE },
  { "set_constr_type", impl_set_constr_type, TRUE },
  { "set_debug", impl_set_debug, TRUE },
  { "set_epsb", impl_set_epsb, TRUE },
  { "set_epsd", impl_set_epsd, TRUE },
  { "set_epsel", impl_set_epsel, TRUE },
  { "set_epsint", impl_set_epsint, TRUE },
  { "set_epslevel", impl_set_epslevel, TRUE },
  { "set_epsperturb", impl_set_epsperturb, TRUE },
  { "set_epspivot", impl_set_epspivot, TRUE },
  { "set_free", impl_set_free, TRUE },
  { "set_improve", impl_set_improve, TRUE },
  { "set_infinite", impl_set_infinite, TRUE },
  { "set_int", impl_set_int, TRUE },
  { "set_lowbo", impl_set_lowbo, TRUE },
  { "set_lp_name", impl_set_lp_name, TRUE },
  { "set_mat", impl_set_mat, TRUE },
  { "set_maxim", impl_set_maxim, TRUE },
  { "set_maxpivot", impl_set_maxpivot, TRUE },
  { "set_minim", impl_set_minim, TRUE },
  { "set_mip_gap", impl_set_mip_gap, TRUE },
  { "set_negrange", impl_set_negrange, TRUE },
  { "set_obj", impl_set_obj, TRUE },
  { "set_obj_bound", impl_set_obj_bound, TRUE },
  { "set_obj_fn", impl_set_obj_fn, TRUE },
  { "set_obj_fnex", impl_set_obj_fn, TRUE },
  { "set_outputfile", impl_set_outputfile, TRUE },
  { "set_pivoting", impl_set_pivoting, TRUE },
  { "set_preferdual", impl_set_preferdual, TRUE },
  { "set_presolve", impl_set_presolve, TRUE },
  { "set_print_sol", impl_set_print_sol, TRUE },
  { "set_rh", impl_set_rh, TRUE },
  { "set_rh_range", impl_set_rh_range, TRUE },
  { "set_rh_vec", impl_set_rh_vec, TRUE },
  { "set_row", impl_set_row, TRUE },
  { "set_rowex", impl_set_row, TRUE },
  { "set_row_name", impl_set_row_name, TRUE },
  { "set_scalelimit", impl_set_scalelimit, TRUE },
  { "set_scaling", impl_set_scaling, TRUE },
  { "set_semicont", impl_set_semicont, TRUE },
  { "set_sense", impl_set_sense, TRUE },
  { "set_simplextype", impl_set_simplextype, TRUE },
  { "set_solutionlimit", impl_set_solutionlimit, TRUE },
  { "set_timeout", impl_set_timeout, TRUE },
  { "set_trace", impl_set_trace, TRUE },
  { "set_unbounded", impl_set_free, TRUE },
  { "set_upbo", impl_set_upbo, TRUE },
  { "set_use_names", impl_set_use_names, TRUE },
  { "set_var_branch", impl_set_var_branch, TRUE },
  { "set_var_weights", impl_set_var_weights, TRUE },
  { "set_verbose", impl_set_verbose, TRUE },
  { "set_XLI", impl_set_XLI, TRUE },
  { "solve", impl_solve, TRUE },
  { "time_elapsed", impl_time_elapsed, TRUE },
  { "unscale", impl_unscale, TRUE },
  { "write_basis", impl_write_basis, TRUE },
  { "write_freemps", impl_write_freemps, TRUE },
  { "write_freeMPS", impl_write_freemps, TRUE },
  { "write_lp", impl_write_lp, TRUE },
  { "write_LP", impl_write_lp, TRUE },
  { "write_mps", impl_write_mps, TRUE },
  { "write_MPS", impl_write_mps, TRUE },
  { "write_params", impl_write_params, TRUE },
  { "write_XLI", impl_write_XLI, TRUE },

/* extra routines */

#if defined DEBUG || defined DEMO
  { "demo", impl_demo, FALSE },
#endif
  { "get_col_names", impl_get_col_name, TRUE },
  { "get_constr_types", impl_get_constr_type, TRUE },
  { "get_int", impl_is_int, TRUE },
  { "get_no_cols", impl_get_Ncolumns, TRUE },
  { "get_no_rows", impl_get_Nrows, TRUE },
  { "get_objective_name", impl_get_objective_name, TRUE },
  { "get_obj_fn", impl_get_obj_fn, TRUE },
  { "get_obj_fun", impl_get_obj_fn, TRUE },
  { "get_problem_name", impl_get_lp_name, TRUE },
  { "get_reduced_costs", impl_get_dual_solution, TRUE },
  { "get_row_names", impl_get_row_name, TRUE },
  { "mat_elm", impl_get_mat, TRUE },
  { "print_handle", impl_print_handle, FALSE },
  { "read_lp_file", impl_read_LP, FALSE },
  { "get_handle", impl_get_handle, FALSE },
};

#if defined WIN32
#  define signalconvention __cdecl
#else
#  define signalconvention
#endif

static void signalconvention SIGINT_func(int sig)
{
        interrupted = TRUE;
}

static void mainloop()
{
	int	i;
        hashelem *hp;

        interrupted = FALSE;
        signal(SIGINT, SIGINT_func);

        if (setjmp(exit_mark) == 0) {

        	if (!initialized) {
        		/* Register the Exit Function */
                        registerExitFcn();

        		/* Allocate a string array to store command */

        		cmd = (char *) calloc(cmdsz, sizeof(*cmd));

        		/* Allocate a string array to store error message */

        		errmsg = (char *) calloc(errmsgsz, sizeof(*errmsg));

                        /* create hashtable of all callbable commands to find them back quickly */
                        cmdhash = create_hash_table(sizeof(routines) / sizeof(*routines), 0);
                        for (i = 0; i < (int) (sizeof(routines)/sizeof(*routines)); i++)
                          	puthash(routines[i].cmd, i, NULL, cmdhash);

        		/* Initialise the lp array, and pointer to the last lp */

        		lp_last = -1;

                        if (!init_lpsolve_lib())
                                ErrMsgTxt("Failed to initialise lpsolve library.\n");

                        handlehash = NULL;

                	initialized = TRUE;
#                       if defined DEBUG
                        	Printf("Initialised\n");
#                       endif
        	}

        	/* Get the first argument as a string matrix */

        	if (nrhs < 1) {
        		int majorversion, minorversion, release, build;

        		/* ErrMsgTxt("At least one command is required."); */
                        lp_solve_version(&majorversion, &minorversion, &release, &build);
                        Printf(strdrivername "  " caller " Interface version " driverVERSION "\n" \
                               "using lpsolve version %d.%d.%d.%d\n\n" \
                               "Usage: [ret1, ret2, ...] = " strdrivername "(%sfunctionname%s, arg1, arg2, ...)\n",
                               majorversion, minorversion, release, build, quotechar, quotechar);
                        return;
                }

        	GetString(prhs, 0, cmd, cmdsz, TRUE);

#               if defined DEBUG
               		Printf("%s\n", cmd);
#               endif

        	/* Now call the required part of the lp toolkit */

                hp = findhash(cmd, cmdhash);
                if (hp == NULL) {
                        strcpy(errmsg, cmd);
                        strncat(errmsg,": Unimplemented.", errmsgsz);
                        ErrMsgTxt(errmsg);
                }
                i = hp->index;

                if (routines[i].needshandle) {
                        if (nrhs < 2)
        		        ErrMsgTxt("An lp handle is required.");

                        if (GetString(prhs, 1, buf, bufsz, FALSE)) {
                                if (handlehash != NULL)
                                	hp = findhash(buf, handlehash);
                                else
                                        hp = NULL;
                                if (hp == NULL) {
                                        char name[bufsz];

                                        strcpy(name, buf);
                                        sprintf(buf, "Invalid model name: %s", name);
                                        ErrMsgTxt(buf);
                                }
                                h = hp->index;
                        }
                        else {
                                h = (int) GetRealScalar(prhs, 1);
                        }

                	if ((!handle_valid(h)) || ((lp0 = lp[h]) == NULL)) {
                        	strcpy(errmsg, cmd);
                        	strncat(errmsg, ": Invalid lp handle.", errmsgsz);
                		ErrMsgTxt(errmsg);
                	}
                }
                BEGIN_INTERRUPT_IMMEDIATELY_IN_FOREIGN_CODE;
                routines[i].routine();
                END_INTERRUPT_IMMEDIATELY_IN_FOREIGN_CODE;
        }
}

callerPrototype(drivername)
{
        publicargs();

        mainloop();

        ExitcallerPrototype;
}
