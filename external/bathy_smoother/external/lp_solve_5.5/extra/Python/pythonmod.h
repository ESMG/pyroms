#include "Python.h"

#include "lp_lib.h"

#define quotechar "'"
#define drivername lpsolve
#define strdrivername "lpsolve"
#define caller "Python"
#define matCalloc calloc
#define matFree free

#define MatrixEl matrix
#define pMatrix MatrixObject *
#define rMatrix PyObject *
#define strArray char **

#define putlogfunc put_logfunc
#define putabortfunc put_abortfunc

#define init_lpsolve_lib() TRUE
#define exit_lpsolve_lib()

#define callerPrototype(callername) PyObject *lpsolve(LprecObject *self, PyObject *args)

#define publicargs() setargs(self, args)

#define registerExitFcn()

#define ExitcallerPrototype if(lhs.type == -1) return(NULL); if(lhs.PyObject == NULL) { Py_INCREF(Py_None); return Py_None; } else { return lhs.PyObject; }

#define BEGIN_INTERRUPT_IMMEDIATELY_IN_FOREIGN_CODE
#define END_INTERRUPT_IMMEDIATELY_IN_FOREIGN_CODE

typedef struct {
    int      type;
    PyObject *PyObject;
} MatrixObject;

typedef struct {
    PyObject_HEAD
    PyObject    *x_attr;        /* Attributes dictionary */
    lprec       *mylprec;       /* the real lprec */
} LprecObject;

extern LprecObject *self;
extern PyObject *args;
extern int nrhs;
extern MatrixObject lhs;

#define prhs args
#define plhs &lhs
#define nlhs 99

#define Double double
#define Long long

extern void setargs(LprecObject *self, PyObject *args);
extern int ErrMsgTxt(char *str);
extern void Printf(char *format, ...);
extern PyObject *GetpMatrix(PyObject *pm, int element);
extern int GetM(PyObject *pm);
extern int GetN(PyObject *pm);
extern Double GetRealScalar(PyObject *pm, int element);
extern int GetIntVector(PyObject *pm, int element, int *vec, int start, int len, int exactcount);
extern int GetRealVector(PyObject *pm, int element, Double *vec, int start, int len, int exactcount);
extern int GetRealSparseVector(PyObject *pm, int element, Double *vec, int *index, int start, int len, int col);
extern int GetString(PyObject *pm, int element, char *buf, int size, int ShowError);
extern strArray GetCellCharItems(PyObject *pm, int element, int len);
extern void GetCellString(char **pa, int element, char *buf, int len);
extern void FreeCellCharItems(strArray pa, int len);
extern double *CreateDoubleMatrix(int m, int n, pMatrix lhs, int element);
extern double *CreateDoubleSparseMatrix(int m, int n, pMatrix lhs, int element);
extern void SetColumnDoubleSparseMatrix(pMatrix lhs, int element, int m, int n, double *mat, int column, Double *arry, int *index, int size, int *nz);
extern void SetDoubleMatrix(double *mat, int m, int n, pMatrix lhs, int element, int freemat);
extern long *CreateLongMatrix(long m, long n, pMatrix lhs, long element);
extern void SetLongMatrix(long *mat, int m, int n, pMatrix lhs, int element, int freemat);
extern void CreateString(char **str, int m, pMatrix lhs, int element);
#define GetRealSparseVector2(prhs, mat, element, vec, index, start, len, col) GetRealSparseVector(prhs, element, vec, index, start, len, col)
