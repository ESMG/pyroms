#include "mex.h"

#include "lp_lib.h"

#if defined WIN32 || defined _WIN32 || defined MSDOS || defined DOS
# define LPSOLVEAPIEXPLICIT

# include "lp_explicit.h"
#else

#if 0
/* static linking */
# define init_lpsolve_lib() TRUE
# define exit_lpsolve_lib()
# define putlogfunc   put_logfunc
# define putabortfunc put_abortfunc

#else
/* dynamic linking */
#if defined driverVERSION

# define LPSOLVEAPIFROMLIB

# include "lp_explicit.h"

static hlpsolve hlpsolve_ = NULL;

#endif

# define init_lpsolve_lib() ((hlpsolve_ != NULL) || (((hlpsolve_ = open_lpsolve_lib(NULL)) != NULL) && init_lpsolve(hlpsolve_)))
# define exit_lpsolve_lib() { if (hlpsolve_ != NULL) close_lpsolve_lib(hlpsolve_); hlpsolve_ = NULL; }
# define putlogfunc put_logfunc
# define putabortfunc put_abortfunc
#endif

#endif

#define quotechar "'"
#define ErrMsgTxt mexErrMsgTxt
#define drivername mxlpsolve
#define strdrivername "mxlpsolve"
#define caller "MATLAB"
#define IsNumeric mxIsNumeric
#define IsComplex mxIsComplex
#define IsSparse mxIsSparse
#define GetM mxGetM
#define GetN mxGetN
#define GetpMatrix(pm, element) ((pm)[element])
#define GetPr mxGetPr
#define matCalloc mxCalloc
#define matRealloc mxRealloc
#define matFree mxFree
#define GetCellString(ppm, element, buf, size) GetString(ppm, element, buf, size, TRUE)
#define Printf mexPrintf

#define MatrixEl mxArray *
#define pMatrix MatrixEl *
#define rMatrix MatrixEl
#define strArray pMatrix

#define callerPrototype(callername) void mexFunction(int nlhs, mxArray *plhs0[], int nrhs, const mxArray *prhs0[])

#define publicargs() \
        mxArray **plhs, **prhs; \
\
        plhs = (mxArray **) plhs0; \
        prhs = (mxArray **) prhs0; \
        setargs(nlhs, nrhs, plhs, prhs)

#define registerExitFcn() if (mexAtExit(ExitFcn)) ErrMsgTxt("Failed to register exit function.\n")

#define ExitcallerPrototype return

#define BEGIN_INTERRUPT_IMMEDIATELY_IN_FOREIGN_CODE
#define END_INTERRUPT_IMMEDIATELY_IN_FOREIGN_CODE

extern int nlhs;
extern int nrhs;
extern pMatrix plhs;
extern pMatrix prhs;

#define Double double
#define Long Double

extern void setargs(int	nlhs0, int nrhs0, pMatrix plhs0, pMatrix prhs0);
extern Double GetRealScalar(pMatrix ppm, int element);
extern int GetIntVector(pMatrix ppm, int element, int *vec, int start, int len, int exactcount);
extern int GetRealVector(pMatrix ppm, int element, Double *vec, int start, int len, int exactcount);
extern int GetRealSparseVector(pMatrix ppm, int element, Double *vec, int *index, int start, int len, int col);
extern int GetString(pMatrix ppm, int element, char *buf, int size, int ShowError);
extern strArray GetCellCharItems(pMatrix ppm, int element, int len);
extern void FreeCellCharItems(strArray pa, int len);
extern Double *CreateDoubleMatrix(int m, int n, pMatrix plhs, int element);
extern Double *CreateDoubleSparseMatrix(int m, int n, pMatrix plhs, int element);
extern void SetColumnDoubleSparseMatrix(pMatrix plhs, int element, int m, int n, double *mat, int column, Double *arry, int *index, int size, int *nz);
#define CreateLongMatrix CreateDoubleMatrix
extern void CreateString(char **str, int n, pMatrix plhs, int element);
#define SetDoubleMatrix(mat, m, n, plhs, element, freemat)
#define SetLongMatrix SetDoubleMatrix
#define GetRealSparseVector2(prhs, mat, element, vec, index, start, len, col) GetRealSparseVector(prhs, element, vec, index, start, len, col)
