#include "lpsolvecaller.h"

#define maxArgs 10

#if 1
  /* represent arrays via lists */
# define PyArray_Size PyList_Size
# define PyArray_New  PyList_New
# define PyArray_SET_ITEM PyList_SET_ITEM
# define PyArray_Resize(Array, Size) { /* Printf("PyList_Size(*Array) was %d, resizing to %d\n", PyList_Size(*Array), Size); */ while ((PyList_Size(*Array) < Size) && (PyList_Append(*Array, Py_None) == 0)); /* Printf("PyList_Size(*Array) becomes %d\n", PyList_Size(*Array)); */ }
#else
  /* represent arrays via tuples */
# define PyArray_Size PyTuple_Size
# define PyArray_New  PyTuple_New
# define PyArray_SET_ITEM PyTuple_SET_ITEM
# define PyArray_Resize _PyTuple_Resize
#endif

LprecObject *self;
PyObject *args;

PyObject *Lprec_ErrorObject;
int Lprec_errorflag;
int nrhs;

MatrixObject lhs;

extern PyObject *lpsolve(LprecObject *self, PyObject *args);

/* List of functions defined in the module */

/* This information is shown in Python via the command help(lpsolve) */

static char lpsolve_doc[] =
"lpsolve('functionname', arg1, arg2, ...) ->  execute lpsolve functionname with args";

static PyMethodDef lpsolve_methods[] = {
    {strdrivername,   (PyCFunction) lpsolve,       METH_VARARGS, lpsolve_doc},
    {NULL,              NULL}           /* sentinel */
};


/* Initialization function for the module (*must* be called initxx) */

DL_EXPORT(void)
    initlpsolve55(void)
{
    PyObject *m, *d;

    /* Create the module and add the functions */
    m = Py_InitModule("lpsolve55", lpsolve_methods);

    /* Add some symbolic constants to the module */
    d = PyModule_GetDict(m);
#if PYTHON_API_VERSION >= 1007
    Lprec_ErrorObject = PyErr_NewException("lpsolve.error", NULL, NULL);
#else
    Lprec_ErrorObject = Py_BuildValue("s", "lpsolve.error");
#endif
    PyDict_SetItemString(d, "error", Lprec_ErrorObject);
    PyDict_SetItemString(d, "LE", PyInt_FromLong(LE));
    PyDict_SetItemString(d, "EQ", PyInt_FromLong(EQ));
    PyDict_SetItemString(d, "GE", PyInt_FromLong(GE));
    PyDict_SetItemString(d, "FR", PyInt_FromLong(FR));
    PyDict_SetItemString(d, "SCALE_NONE", PyInt_FromLong(SCALE_NONE));
    PyDict_SetItemString(d, "SCALE_EXTREME", PyInt_FromLong(SCALE_EXTREME));
    PyDict_SetItemString(d, "SCALE_RANGE", PyInt_FromLong(SCALE_RANGE));
    PyDict_SetItemString(d, "SCALE_MEAN", PyInt_FromLong(SCALE_MEAN));
    PyDict_SetItemString(d, "SCALE_GEOMETRIC", PyInt_FromLong(SCALE_GEOMETRIC));
    PyDict_SetItemString(d, "SCALE_CURTISREID", PyInt_FromLong(SCALE_CURTISREID));
    PyDict_SetItemString(d, "SCALE_QUADRATIC", PyInt_FromLong(SCALE_QUADRATIC));
    PyDict_SetItemString(d, "SCALE_LOGARITHMIC", PyInt_FromLong(SCALE_LOGARITHMIC));
    PyDict_SetItemString(d, "SCALE_USERWEIGHT", PyInt_FromLong(SCALE_USERWEIGHT));
    PyDict_SetItemString(d, "SCALE_POWER2", PyInt_FromLong(SCALE_POWER2));
    PyDict_SetItemString(d, "SCALE_EQUILIBRATE", PyInt_FromLong(SCALE_EQUILIBRATE));
    PyDict_SetItemString(d, "SCALE_INTEGERS", PyInt_FromLong(SCALE_INTEGERS));
    PyDict_SetItemString(d, "SCALE_DYNUPDATE", PyInt_FromLong(SCALE_DYNUPDATE));
    PyDict_SetItemString(d, "SCALE_ROWSONLY", PyInt_FromLong(SCALE_ROWSONLY));
    PyDict_SetItemString(d, "SCALE_COLSONLY", PyInt_FromLong(SCALE_COLSONLY));
    PyDict_SetItemString(d, "IMPROVE_NONE", PyInt_FromLong(IMPROVE_NONE));
    PyDict_SetItemString(d, "IMPROVE_SOLUTION", PyInt_FromLong(IMPROVE_SOLUTION));
    PyDict_SetItemString(d, "IMPROVE_DUALFEAS", PyInt_FromLong(IMPROVE_DUALFEAS));
    PyDict_SetItemString(d, "IMPROVE_THETAGAP", PyInt_FromLong(IMPROVE_THETAGAP));
    PyDict_SetItemString(d, "IMPROVE_BBSIMPLEX", PyInt_FromLong(IMPROVE_BBSIMPLEX));
    PyDict_SetItemString(d, "PRICER_FIRSTINDEX", PyInt_FromLong(PRICER_FIRSTINDEX));
    PyDict_SetItemString(d, "PRICER_DANTZIG", PyInt_FromLong(PRICER_DANTZIG));
    PyDict_SetItemString(d, "PRICER_DEVEX", PyInt_FromLong(PRICER_DEVEX));
    PyDict_SetItemString(d, "PRICER_STEEPESTEDGE", PyInt_FromLong(PRICER_STEEPESTEDGE));
    PyDict_SetItemString(d, "PRICE_PRIMALFALLBACK", PyInt_FromLong(PRICE_PRIMALFALLBACK));
    PyDict_SetItemString(d, "PRICE_MULTIPLE", PyInt_FromLong(PRICE_MULTIPLE));
    PyDict_SetItemString(d, "PRICE_PARTIAL", PyInt_FromLong(PRICE_PARTIAL));
    PyDict_SetItemString(d, "PRICE_ADAPTIVE", PyInt_FromLong(PRICE_ADAPTIVE));
    PyDict_SetItemString(d, "PRICE_RANDOMIZE", PyInt_FromLong(PRICE_RANDOMIZE));
    PyDict_SetItemString(d, "PRICE_AUTOPARTIAL", PyInt_FromLong(PRICE_AUTOPARTIAL));
    PyDict_SetItemString(d, "PRICE_LOOPLEFT", PyInt_FromLong(PRICE_LOOPLEFT));
    PyDict_SetItemString(d, "PRICE_LOOPALTERNATE", PyInt_FromLong(PRICE_LOOPALTERNATE));
    PyDict_SetItemString(d, "PRICE_HARRISTWOPASS", PyInt_FromLong(PRICE_HARRISTWOPASS));
    PyDict_SetItemString(d, "PRICE_TRUENORMINIT", PyInt_FromLong(PRICE_TRUENORMINIT));
    PyDict_SetItemString(d, "PRESOLVE_NONE", PyInt_FromLong(PRESOLVE_NONE));
    PyDict_SetItemString(d, "PRESOLVE_ROWS", PyInt_FromLong(PRESOLVE_ROWS));
    PyDict_SetItemString(d, "PRESOLVE_COLS", PyInt_FromLong(PRESOLVE_COLS));
    PyDict_SetItemString(d, "PRESOLVE_LINDEP", PyInt_FromLong(PRESOLVE_LINDEP));
    PyDict_SetItemString(d, "PRESOLVE_SOS", PyInt_FromLong(PRESOLVE_SOS));
    PyDict_SetItemString(d, "PRESOLVE_REDUCEMIP", PyInt_FromLong(PRESOLVE_REDUCEMIP));
    PyDict_SetItemString(d, "PRESOLVE_KNAPSACK", PyInt_FromLong(PRESOLVE_KNAPSACK));
    PyDict_SetItemString(d, "PRESOLVE_ELIMEQ2", PyInt_FromLong(PRESOLVE_ELIMEQ2));
    PyDict_SetItemString(d, "PRESOLVE_IMPLIEDFREE", PyInt_FromLong(PRESOLVE_IMPLIEDFREE));
    PyDict_SetItemString(d, "PRESOLVE_REDUCEGCD", PyInt_FromLong(PRESOLVE_REDUCEGCD));
    PyDict_SetItemString(d, "PRESOLVE_PROBEFIX", PyInt_FromLong(PRESOLVE_PROBEFIX));
    PyDict_SetItemString(d, "PRESOLVE_PROBEREDUCE", PyInt_FromLong(PRESOLVE_PROBEREDUCE));
    PyDict_SetItemString(d, "PRESOLVE_ROWDOMINATE", PyInt_FromLong(PRESOLVE_ROWDOMINATE));
    PyDict_SetItemString(d, "PRESOLVE_COLDOMINATE", PyInt_FromLong(PRESOLVE_COLDOMINATE));
    PyDict_SetItemString(d, "PRESOLVE_MERGEROWS", PyInt_FromLong(PRESOLVE_MERGEROWS));
    PyDict_SetItemString(d, "PRESOLVE_IMPLIEDSLK", PyInt_FromLong(PRESOLVE_IMPLIEDSLK));
    PyDict_SetItemString(d, "PRESOLVE_COLFIXDUAL", PyInt_FromLong(PRESOLVE_COLFIXDUAL));
    PyDict_SetItemString(d, "PRESOLVE_BOUNDS", PyInt_FromLong(PRESOLVE_BOUNDS));
    PyDict_SetItemString(d, "PRESOLVE_DUALS", PyInt_FromLong(PRESOLVE_DUALS));
    PyDict_SetItemString(d, "PRESOLVE_SENSDUALS", PyInt_FromLong(PRESOLVE_SENSDUALS));
    PyDict_SetItemString(d, "ANTIDEGEN_NONE", PyInt_FromLong(ANTIDEGEN_NONE));
    PyDict_SetItemString(d, "ANTIDEGEN_FIXEDVARS", PyInt_FromLong(ANTIDEGEN_FIXEDVARS));
    PyDict_SetItemString(d, "ANTIDEGEN_COLUMNCHECK", PyInt_FromLong(ANTIDEGEN_COLUMNCHECK));
    PyDict_SetItemString(d, "ANTIDEGEN_STALLING", PyInt_FromLong(ANTIDEGEN_STALLING));
    PyDict_SetItemString(d, "ANTIDEGEN_NUMFAILURE", PyInt_FromLong(ANTIDEGEN_NUMFAILURE));
    PyDict_SetItemString(d, "ANTIDEGEN_LOSTFEAS", PyInt_FromLong(ANTIDEGEN_LOSTFEAS));
    PyDict_SetItemString(d, "ANTIDEGEN_INFEASIBLE", PyInt_FromLong(ANTIDEGEN_INFEASIBLE));
    PyDict_SetItemString(d, "ANTIDEGEN_DYNAMIC", PyInt_FromLong(ANTIDEGEN_DYNAMIC));
    PyDict_SetItemString(d, "ANTIDEGEN_DURINGBB", PyInt_FromLong(ANTIDEGEN_DURINGBB));
    PyDict_SetItemString(d, "ANTIDEGEN_RHSPERTURB", PyInt_FromLong(ANTIDEGEN_RHSPERTURB));
    PyDict_SetItemString(d, "ANTIDEGEN_BOUNDFLIP", PyInt_FromLong(ANTIDEGEN_BOUNDFLIP));
    PyDict_SetItemString(d, "CRASH_NONE", PyInt_FromLong(CRASH_NONE));
    PyDict_SetItemString(d, "CRASH_MOSTFEASIBLE", PyInt_FromLong(CRASH_MOSTFEASIBLE));
    PyDict_SetItemString(d, "CRASH_LEASTDEGENERATE", PyInt_FromLong(CRASH_LEASTDEGENERATE));
    PyDict_SetItemString(d, "SIMPLEX_PRIMAL_PRIMAL", PyInt_FromLong(SIMPLEX_PRIMAL_PRIMAL));
    PyDict_SetItemString(d, "SIMPLEX_DUAL_PRIMAL", PyInt_FromLong(SIMPLEX_DUAL_PRIMAL));
    PyDict_SetItemString(d, "SIMPLEX_PRIMAL_DUAL", PyInt_FromLong(SIMPLEX_PRIMAL_DUAL));
    PyDict_SetItemString(d, "SIMPLEX_DUAL_DUAL", PyInt_FromLong(SIMPLEX_DUAL_DUAL));
    PyDict_SetItemString(d, "NODE_FIRSTSELECT", PyInt_FromLong(NODE_FIRSTSELECT));
    PyDict_SetItemString(d, "NODE_GAPSELECT", PyInt_FromLong(NODE_GAPSELECT));
    PyDict_SetItemString(d, "NODE_RANGESELECT", PyInt_FromLong(NODE_RANGESELECT));
    PyDict_SetItemString(d, "NODE_FRACTIONSELECT", PyInt_FromLong(NODE_FRACTIONSELECT));
    PyDict_SetItemString(d, "NODE_PSEUDOCOSTSELECT", PyInt_FromLong(NODE_PSEUDOCOSTSELECT));
    PyDict_SetItemString(d, "NODE_PSEUDONONINTSELECT", PyInt_FromLong(NODE_PSEUDONONINTSELECT));
    PyDict_SetItemString(d, "NODE_PSEUDORATIOSELECT", PyInt_FromLong(NODE_PSEUDORATIOSELECT));
    PyDict_SetItemString(d, "NODE_USERSELECT", PyInt_FromLong(NODE_USERSELECT));
    PyDict_SetItemString(d, "NODE_WEIGHTREVERSEMODE", PyInt_FromLong(NODE_WEIGHTREVERSEMODE));
    PyDict_SetItemString(d, "NODE_BRANCHREVERSEMODE", PyInt_FromLong(NODE_BRANCHREVERSEMODE));
    PyDict_SetItemString(d, "NODE_GREEDYMODE", PyInt_FromLong(NODE_GREEDYMODE));
    PyDict_SetItemString(d, "NODE_PSEUDOCOSTMODE", PyInt_FromLong(NODE_PSEUDOCOSTMODE));
    PyDict_SetItemString(d, "NODE_DEPTHFIRSTMODE", PyInt_FromLong(NODE_DEPTHFIRSTMODE));
    PyDict_SetItemString(d, "NODE_RANDOMIZEMODE", PyInt_FromLong(NODE_RANDOMIZEMODE));
    PyDict_SetItemString(d, "NODE_GUBMODE", PyInt_FromLong(NODE_GUBMODE));
    PyDict_SetItemString(d, "NODE_DYNAMICMODE", PyInt_FromLong(NODE_DYNAMICMODE));
    PyDict_SetItemString(d, "NODE_RESTARTMODE", PyInt_FromLong(NODE_RESTARTMODE));
    PyDict_SetItemString(d, "NODE_BREADTHFIRSTMODE", PyInt_FromLong(NODE_BREADTHFIRSTMODE));
    PyDict_SetItemString(d, "NODE_AUTOORDER", PyInt_FromLong(NODE_AUTOORDER));
    PyDict_SetItemString(d, "NODE_RCOSTFIXING", PyInt_FromLong(NODE_RCOSTFIXING));
    PyDict_SetItemString(d, "NODE_STRONGINIT", PyInt_FromLong(NODE_STRONGINIT));
    PyDict_SetItemString(d, "NOMEMORY", PyInt_FromLong(NOMEMORY));
    PyDict_SetItemString(d, "OPTIMAL", PyInt_FromLong(OPTIMAL));
    PyDict_SetItemString(d, "SUBOPTIMAL", PyInt_FromLong(SUBOPTIMAL));
    PyDict_SetItemString(d, "INFEASIBLE", PyInt_FromLong(INFEASIBLE));
    PyDict_SetItemString(d, "UNBOUNDED", PyInt_FromLong(UNBOUNDED));
    PyDict_SetItemString(d, "DEGENERATE", PyInt_FromLong(DEGENERATE));
    PyDict_SetItemString(d, "NUMFAILURE", PyInt_FromLong(NUMFAILURE));
    PyDict_SetItemString(d, "USERABORT", PyInt_FromLong(USERABORT));
    PyDict_SetItemString(d, "TIMEOUT", PyInt_FromLong(TIMEOUT));
    PyDict_SetItemString(d, "PROCFAIL", PyInt_FromLong(PROCFAIL));
    PyDict_SetItemString(d, "PROCBREAK", PyInt_FromLong(PROCBREAK));
    PyDict_SetItemString(d, "FEASFOUND", PyInt_FromLong(FEASFOUND));
    PyDict_SetItemString(d, "NOFEASFOUND", PyInt_FromLong(NOFEASFOUND));
    PyDict_SetItemString(d, "BRANCH_CEILING", PyInt_FromLong(BRANCH_CEILING));
    PyDict_SetItemString(d, "BRANCH_FLOOR", PyInt_FromLong(BRANCH_FLOOR));
    PyDict_SetItemString(d, "BRANCH_AUTOMATIC", PyInt_FromLong(BRANCH_AUTOMATIC));
    PyDict_SetItemString(d, "MSG_PRESOLVE", PyInt_FromLong(MSG_PRESOLVE));
    PyDict_SetItemString(d, "MSG_LPFEASIBLE", PyInt_FromLong(MSG_LPFEASIBLE));
    PyDict_SetItemString(d, "MSG_LPOPTIMAL", PyInt_FromLong(MSG_LPOPTIMAL));
    PyDict_SetItemString(d, "MSG_MILPEQUAL", PyInt_FromLong(MSG_MILPEQUAL));
    PyDict_SetItemString(d, "MSG_MILPFEASIBLE", PyInt_FromLong(MSG_MILPFEASIBLE));
    PyDict_SetItemString(d, "MSG_MILPBETTER", PyInt_FromLong(MSG_MILPBETTER));
    PyDict_SetItemString(d, "NEUTRAL", PyInt_FromLong(NEUTRAL));
    PyDict_SetItemString(d, "CRITICAL", PyInt_FromLong(CRITICAL));
    PyDict_SetItemString(d, "SEVERE", PyInt_FromLong(SEVERE));
    PyDict_SetItemString(d, "IMPORTANT", PyInt_FromLong(IMPORTANT));
    PyDict_SetItemString(d, "NORMAL", PyInt_FromLong(NORMAL));
    PyDict_SetItemString(d, "DETAILED", PyInt_FromLong(DETAILED));
    PyDict_SetItemString(d, "FULL", PyInt_FromLong(FULL));
    PyDict_SetItemString(d, "Infinite", PyFloat_FromDouble(DEF_INFINITE));
}

void Printf(char *format, ...)
{
        va_list ap;
        static char buf[4096];

	va_start(ap, format);
	vsnprintf(buf, sizeof(buf), format, ap);
        /* printf("%s", buf); */
        PySys_WriteStdout("%s", buf);
	va_end(ap);
}

int ErrMsgTxt(char *str)
{
        PyErr_SetString(Lprec_ErrorObject, str);
	/* Printf("%s\n", str); */
        lhs.type = -1;
        exitnow();

        return(0);
}

void setargs(LprecObject *self0, PyObject *args0)
{
        PyObject *ar[maxArgs];

        self = self0;
        args = args0;

        for(nrhs = (sizeof(ar) / sizeof(*ar)) - 1; nrhs >= 0; nrhs--)
                ar[nrhs] = NULL;
        PyArg_UnpackTuple(args, strdrivername, 0, sizeof(ar)/sizeof(*ar), &ar[0], &ar[1], &ar[2], &ar[3], &ar[4], &ar[5], &ar[6], &ar[7], &ar[8], &ar[9]);
        for(nrhs = (sizeof(ar) / sizeof(*ar)) - 1; (nrhs >= 0) && (ar[nrhs] == NULL); nrhs--);
        nrhs++;
        lhs.type = 0;
        lhs.PyObject = NULL;
}

PyObject *GetpMatrix(PyObject *pm, int element)
{
        PyObject *ar[maxArgs];
        int i;
        char OK = TRUE;

        for(i = (sizeof(ar) / sizeof(*ar)) - 1; i >= 0; i--)
                ar[i] = NULL;
        PyArg_UnpackTuple(pm, strdrivername, 0, sizeof(ar)/sizeof(*ar), &ar[0], &ar[1], &ar[2], &ar[3], &ar[4], &ar[5], &ar[6], &ar[7], &ar[8], &ar[9]);
        if((element >= sizeof(ar) / sizeof(*ar)) ||
           (ar[element] == NULL)) {
                PyErr_Clear();
                return(NULL);
        }
        return(ar[element]);
}

int GetM(PyObject *pm)
{
        int m;

        if (PyNumber_Check(pm)) {
                m = 1;
        }
        else {
            	m = PyObject_Length(pm);
        }
        return(m);
}

int GetN(PyObject *pm)
{
        int n;
        PyObject *item;

        if (PyNumber_Check(pm)) {
                n = 1;
        }
        else {
                item = PySequence_GetItem(pm, 0);
                if (item != NULL) {
                        if (PyNumber_Check(item)) {
                                n = 1;
                        }
                        else {
                                n = PyObject_Length(item);
                        }
                        Py_DECREF(item);
                }
                else
                        n = 0;
        }
        return(n);
}

/* Function to get a real scalar with error checking */

Double GetRealScalar(PyObject *pm, int element)
{
        PyObject *item;
        Double a;
        char OK = TRUE;

        item = GetpMatrix(pm, element);
        if (item == NULL)
                OK = FALSE;
        else {
                if (!PyNumber_Check(item))
                        OK = FALSE;
                else
        		a = PyFloat_AsDouble(item);
        }

        if (!OK) {
		ErrMsgTxt("Expecting a scalar argument.");
        }

        return(a);
}

#define GetVector(pm, element, vec, cast, start, len, exactcount, ret) \
{ \
        PyObject *vector, *item; \
        int k, m, n, count = 0; \
        char OK = TRUE, isvector; \
 \
        vector = GetpMatrix(pm, element); \
        if (vector == NULL) \
                OK = FALSE; \
        else { \
                if (PyNumber_Check(vector)) { \
                        m = 1; \
                        isvector = FALSE; \
                } \
                else { \
                    	m = PyObject_Length(vector); \
                        isvector = TRUE; \
                } \
                n = 1; \
 \
        	if ( !((m == 1) || (n == 1)) || \
                     ((m == 1) && (((exactcount) && (len != n)) || ((!exactcount) && (n > len)))) || \
                     ((n == 1) && (((exactcount) && (len != m)) || ((!exactcount) && (m > len))))) { \
                        PyErr_Clear(); \
                        ErrMsgTxt("invalid vector."); \
        	} \
 \
                if (n == 1) \
                        len = m; \
                else \
                        len = n; \
                vec += start; \
 \
                k = 0; \
        	for (; k < m; k++) { \
        	    Lprec_errorflag = 0; \
                    if (isvector) \
        	    	item = PySequence_GetItem(vector, k); \
                    else \
                        item = vector; \
        	    if ((item == NULL) || (!PyNumber_Check(item))) { \
                        if ((isvector) && (item != NULL)) \
            		    Py_DECREF(item); \
                        ErrMsgTxt("invalid vector (non-numerical item)."); \
        	    } \
                    *(vec++) = (cast) PyFloat_AsDouble(item); \
                    if (isvector) \
        	    	Py_DECREF(item); \
        	    if (Lprec_errorflag) \
                        ErrMsgTxt("invalid vector."); \
        	} \
        } \
 \
        count = len; \
 \
        ret = count; \
}

/* Functions to get len elements from a Python vector. */

int GetIntVector(PyObject *pm, int element, int *vec, int start, int len, int exactcount)
{
        int ret;

        GetVector(pm, element, vec, int, start, len, exactcount, ret)

        return(ret);
}

int GetRealVector(PyObject *pm, int element, Double *vec, int start, int len, int exactcount)
{
        int ret;

        GetVector(pm, element, vec, Double, start, len, exactcount, ret)

        return(ret);
}


/* Function to get max len elements from a Python vector.
   Python doesn't know sparse matrixes, so a matrix is converted to sparse */

int GetRealSparseVector(PyObject *pm, int element, Double *vec, int *index, int start, int len, int col)
{
        PyObject *vector, *item, *item1, *item2;
        int k, m, n, count = 0;
        Double a;
        char OK = TRUE, isvector;

        vector = GetpMatrix(pm, element);
        if (vector == NULL)
                OK = FALSE;
        else {
                if (PyNumber_Check(vector)) {
                        m = 1;
                        isvector = FALSE;
                }
                else {
                    	m = PyObject_Length(vector);
                        isvector = TRUE;
                }
                if ((col == 0) || (!isvector))
                	n = 1;
                else
                        n = col;

        	if (  ((col == 0) && (((m != 1) && (n != 1)) || ((m == 1) && (n > len)) || ((n == 1) && (m > len)))) ||
                      ((col != 0) && ((m > len) || (col > n))) ) {
                        PyErr_Clear();
        		/* Printf("1: m=%d, n=%d, col=%d, len=%d\n", m,n,col,len); */
        		ErrMsgTxt("invalid vector.");
        	}

                if ((((n == 1) || (col != 0)) && (m != len)) || ((col == 0) && (m == 1) && (n != len))) {
        		/* Printf("2: m=%d, n=%d, col=%d, len=%d\n", m,n,col,len); */
                        ErrMsgTxt("invalid vector.");
                }

        	for (k = 0; k < m; k++) {
        	    Lprec_errorflag = 0;
                    item1 = item2 = NULL;
                    if (isvector)
        	    	item = item1 = PySequence_GetItem(vector, k);
                    else
                        item = vector;
                    if ((item != NULL) && (col) && (isvector)) {
                        if (!PyNumber_Check(item)) {
                            n = PyObject_Length(item);
                            if (col <= n)
                            	item = item2 = PySequence_GetItem(item, col - 1);
                        }
                        else
                            n = 1;
                    }
                    else
                        n = 1;
                    if (col > n) {
                            if (item2 != NULL)
                	    	Py_DECREF(item2);
                            if (item1 != NULL)
                	    	Py_DECREF(item1);
                            PyErr_Clear();
                            ErrMsgTxt("invalid vector.");
                    }

        	    if ((item == NULL) || (!PyNumber_Check(item))) {
                        if (item2 != NULL)
            		    Py_DECREF(item2);
                        if (item1 != NULL)
            		    Py_DECREF(item1);
                        PyErr_Clear();
                        ErrMsgTxt("invalid vector (non-numerical item).");
        	    }
                    a = PyFloat_AsDouble(item);
                    if (a) {
                        *(vec++) = a;
                        *(index++) = start + k;
                        count++;
                    }
                    if (item2 != NULL)
        	    	Py_DECREF(item2);
                    if (item1 != NULL)
        	    	Py_DECREF(item1);
        	    if (Lprec_errorflag)
                        ErrMsgTxt("invalid vector.");
        	}
        }

        return(count);
}


int GetString(PyObject *pm, int element, char *buf, int size, int ShowError)
{
        PyObject *item;
        int size1;
        char *ptr = NULL;

        item = GetpMatrix(pm, element);
        if((item == NULL) ||
           (PyString_AsStringAndSize(item, &ptr, &size1) != 0) ||
           (ptr == NULL)) {
                PyErr_Clear();
                if(ShowError)
                        ErrMsgTxt("Expecting a character element.");
                return(FALSE);
        }
        if(size1 < size)
                size = size1;
        memcpy(buf, ptr, size);
        buf[size] = 0;
        return(TRUE);
}

strArray GetCellCharItems(PyObject *pm, int element, int len)
{
        int m, i, size1;
        PyObject *vector, *item;
        char **pa0, **pa, *str, isvector, *ptr;

        vector = GetpMatrix(pm, element);
        if(vector == NULL) {
                PyErr_Clear();
                ErrMsgTxt("Expecting a character array.");
        }
        if (PyString_Check(vector)) {
                m = 1;
                isvector = FALSE;
        }
        else {
            	m = PyObject_Length(vector);
                if (m == -1) {
                	PyErr_Clear();
                	ErrMsgTxt("Expecting a character array.");
                }
                isvector = TRUE;
        }

        if (!(m == len))
                ErrMsgTxt("invalid vector.");

        pa = pa0 = (char **) matCalloc(len, sizeof(*pa));
        for (i = 0; i < len; i++) {
    		Lprec_errorflag = 0;
                if (isvector)
    		    item = PySequence_GetItem(vector, i);
                else
                    item = vector;
    		if ((item == NULL) || (!PyString_Check(item))) {
                    PyErr_Clear();
                    if ((isvector) && (item != NULL))
        		    Py_DECREF(item);
                    FreeCellCharItems(pa, i);
                    ErrMsgTxt("invalid vector (non-string item).");
    		}
           	if ((PyString_AsStringAndSize(item, &ptr, &size1) != 0) ||
                    (ptr == NULL)) {
                	PyErr_Clear();
                        if (isvector)
            			Py_DECREF(item);
                        FreeCellCharItems(pa, i);
                	ErrMsgTxt("Expecting a character element.");
                }

                str = pa[i] = (char *) matCalloc(size1 + 1, sizeof(*str));
                memcpy(str, ptr, size1);
                str[size1] = 0;
                if (isvector)
    			Py_DECREF(item);
    		if (Lprec_errorflag) {
                    FreeCellCharItems(pa, i + 1);
                    ErrMsgTxt("invalid vector.");
                }
        }
        return(pa0);
}

void GetCellString(char **pa, int element, char *buf, int len)
{
        strncpy(buf, pa[element], len);
}

void FreeCellCharItems(strArray pa, int len)
{
        int i;

        for (i = 0; i < len; i++)
                free(pa[i]);
	matFree(pa);
}

static void setlhs(pMatrix lhs, long element, PyObject *PyObject1)
{
        int size;
        PyObject *PyObject2;

        if(element == 0) {
                lhs->type = 1;
                lhs->PyObject = PyObject1;
        }
        else {
                if (lhs->type == 2) {
                	size = PyArray_Size(lhs->PyObject);
                        if(size == -1)
                        	PyErr_Clear();
                }
                else
                        size = -1;
                if(size == -1) {
                        PyObject2 = lhs->PyObject;
                        lhs->type = 2;
                        lhs->PyObject = PyArray_New(1 + element);
                        if(PyObject2 != NULL)
                                PyArray_SET_ITEM(lhs->PyObject, 0, PyObject2);
                }
                else if (element >= size) {
                        PyArray_Resize(&(lhs->PyObject), 1 + element);
                }
                PyArray_SET_ITEM(lhs->PyObject, element, PyObject1);
        }
}

double *CreateDoubleMatrix(int m, int n, pMatrix lhs, int element)
{
        return((double *) malloc(m * n * sizeof(double)));
}

void SetDoubleMatrix(double *mat, int m, int n, pMatrix lhs, int element, int freemat)
{
        if (mat != NULL) {
                if (m * n == 1)
                        setlhs(lhs, element, PyFloat_FromDouble(*mat));
                else {
                        PyObject *PyObject1, *PyObject2;
                        int i, j, len, size;
                        double *mat1;

                        if (m == 1) {
                                len = m;
                                size = n;
                        }
                        else {
                                len = n;
                                size = m;
                        }

                        PyObject1 = PyArray_New(size);
                        mat1 = mat;
                        for (i = 0; i < size; i++) {
                                if (len == 1) {
                    			PyArray_SET_ITEM(PyObject1, i, PyFloat_FromDouble(*(mat1++)));
                                }
                                else {
                                        PyObject2 = PyArray_New(len);
                                        mat1 = mat + i;
                                        for (j = 0; j < len; j++, mat1 += size)
                                                PyArray_SET_ITEM(PyObject2, j, PyFloat_FromDouble(*mat1));
                                        PyArray_SET_ITEM(PyObject1, i, PyObject2);
                                }
                        }
                        setlhs(lhs, element, PyObject1);
                }

        	if (freemat)
        		free(mat);
        }
}

long *CreateLongMatrix(long m, long n, pMatrix lhs, long element)
{
        return((long *) malloc(m * n * sizeof(long)));
}

void SetLongMatrix(long *mat, int m, int n, pMatrix lhs, int element, int freemat)
{
        if (mat != NULL) {
                if (m * n == 1)
                        setlhs(lhs, element, PyLong_FromLong(*mat));
                else {
                        PyObject *PyObject1, *PyObject2;
                        int i, j, len, size;
                        long *mat1;

                        if (m == 1) {
                                len = m;
                                size = n;
                        }
                        else {
                                len = n;
                                size = m;
                        }
                        PyObject1 = PyArray_New(size);
                        mat1 = mat;
                        for (i = 0; i < size; i++) {
                                if (len == 1) {
                    			PyArray_SET_ITEM(PyObject1, i, PyLong_FromLong(*(mat1++)));
                                }
                                else {
                                        PyObject2 = PyArray_New(len);
                                        mat1 = mat + i;
                                        for (j = 0; j < len; j++, mat1 += size)
                                                PyArray_SET_ITEM(PyObject2, j, PyLong_FromLong(*mat1));
                                        PyArray_SET_ITEM(PyObject1, i, PyObject2);
                                }
                        }
                        setlhs(lhs, element, PyObject1);
                }

        	if (freemat)
        		free(mat);
        }
}

double *CreateDoubleSparseMatrix(int m, int n, pMatrix lhs, int element)
{
        return((double *) malloc(m * n * sizeof(double)));
}

void SetColumnDoubleSparseMatrix(pMatrix lhs, int element, int m, int n, double *mat, int column, double *arry, int *index, int size, int *nz)
{
        double *sr = mat + (column - 1) * m, a;
        int ii, i, j = -1;

        for (i = 0; (i < size); i++) {
                a = arry[i];
                if (a) {
                        if (index == NULL)
                                ii = i;
                        else
                                ii = index[i] - 1;

                        while (++j < ii)
                                sr[j] = 0.0;

                        sr[ii] = a;
                }
        }

        while (++j < m)
                sr[j] = 0.0;

        *nz += m;
}

void CreateString(char **str, int m, pMatrix lhs, int element)
{
        if(m == 1)
                setlhs(lhs, element, PyString_FromString(*str));
        else {
                PyObject *PyObject1;
                int i, len = m;

                PyObject1 = PyArray_New(len);
                for (i = 0; i < len; i++)
            		PyArray_SET_ITEM(PyObject1, i, PyString_FromString(*(str++)));
                setlhs(lhs, element, PyObject1);
        }
}
