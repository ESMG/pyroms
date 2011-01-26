#include "matlab.h"

int nlhs;
int nrhs;
pMatrix plhs;
pMatrix prhs;

void setargs(int nlhs0, int nrhs0, pMatrix plhs0, pMatrix prhs0)
{
        nlhs = nlhs0;
        nrhs = nrhs0;
        plhs = plhs0;
        prhs = prhs0;
}

/* Function to get a real scalar with error checking */

Double GetRealScalar(pMatrix ppm, int element)
{
        rMatrix pm = ppm[element];

	if ((GetM(pm) == 1) && (GetN(pm) == 1)
	  && (IsNumeric(pm)) && (!IsComplex(pm)) ) {
		return(mxGetScalar(pm));
	} else {
		ErrMsgTxt("Expecting a scalar argument.");
	}
        return(0.0);
}

#define GetVector(pm, vec, cast, start, len, exactcount, ret) \
{ \
	int	j, k, k1, m, n, count = 0; \
	int	*ir, *jc; \
	double	*pr, *pr0; \
        cast    *vec0; \
 \
	m = GetM(pm); \
	n = GetN(pm); \
 \
	if ( !((m == 1) || (n == 1)) || \
             ((m == 1) && (((exactcount) && (len != n)) || ((!exactcount) && (n > len)))) || \
             ((n == 1) && (((exactcount) && (len != m)) || ((!exactcount) && (m > len)))) || \
	     !IsNumeric(pm) || IsComplex(pm)) { \
		ErrMsgTxt("invalid vector."); \
	} \
 \
	pr = GetPr(pm); \
 \
	if (!IsSparse(pm)) { \
                if (n == 1) \
                        len = m; \
                else \
                        len = n; \
                vec += start; \
		for (k = 0; k < len; k++, pr++, vec++) { \
			*vec = (cast) *pr; \
		} \
                count = len; \
	} else if (IsSparse(pm)) { \
		jc = mxGetJc(pm); \
		ir = mxGetIr(pm); \
                pr0 = pr; \
                vec0 = vec; \
                for (j = 0; j < n; j++) { \
                        k = jc[j]; \
                        k1 = jc[j + 1]; \
                        pr = pr0 + k; \
                        vec = vec0; \
                        vec += start + j * m; \
                        for (; k < k1; k++, pr++) { \
				vec[ir[k]] = (cast) *pr; \
                                count++; \
			} \
		} \
	} else { \
		ErrMsgTxt("Can't figure out this matrix."); \
	} \
 \
        ret = count; \
}

/* Functions to get len elements from a MATLAB vector. Matrix
   can be either full or sparse. Elements are stored in indices
   start..start+n-1  Errors out if the MATLAB vector is not length len */

int GetIntVector(pMatrix ppm, int element, int *vec, int start, int len, int exactcount)
{
        int ret;
        rMatrix pm = ppm[element];

        GetVector(pm, vec, int, start, len, exactcount, ret)

        return(ret);
}

int GetRealVector(pMatrix ppm, int element, Double *vec, int start, int len, int exactcount)
{
        int ret;
        rMatrix pm = ppm[element];

        GetVector(pm, vec, Double, start, len, exactcount, ret)

        return(ret);
}


/* Function to get max len elements from a MATLAB sparse vector. Matrix
   can be either full or sparse. Elements are stored in indices
   start..start+n-1  Errors out if the MATLAB vector is longer than length len */

int GetRealSparseVector(pMatrix ppm, int element, Double *vec, int *index, int start, int len, int col)
{
	int	j, k, k1, m, n, start1, count = 0;
	int	*ir, *jc;
	double	*pr, *pr0;
        rMatrix pm = ppm[element];

	m = GetM(pm);
	n = GetN(pm);

	if (  ((col == 0) && (((m != 1) && (n != 1)) || ((m == 1) && (n > len)) || ((n == 1) && (m > len)))) ||
              ((col != 0) && ((m > len) || (col > n))) ||
	      !IsNumeric(pm) ||
              IsComplex(pm)  ) {
		ErrMsgTxt("invalid vector.");
	}

	pr = GetPr(pm);

	if (!IsSparse(pm)) {
                if ((((n == 1) || (col != 0)) && (m != len)) || ((col == 0) && (m == 1) && (n != len)))
                        ErrMsgTxt("invalid vector.");

                if (col)
                	pr += (col - 1) * m;
                for (k = 0; k < len; k++, pr++) {
                        if (*pr) {
				*(vec++) = *pr;
                        	*(index++) = start + k;
                        	count++;
                        }
		}
	} else if (IsSparse(pm)) {
                int j1, j2;

		jc = mxGetJc(pm);
		ir = mxGetIr(pm);
                pr0 = pr;
                if (col == 0) {
                        j1 = 0;
                        j2 = n;
                }
                else {
                        j1 = col - 1;
                        j2 = col;
                }
		for (j = j1; j < j2; j++) {
                        k = jc[j];
                        k1 = jc[j + 1];
                        pr = pr0 + k;
                        start1 = start;
                        if (col == 0)
                        	start1 += j * m;
                        for (; k < k1; k++, pr++, vec++, index++) {
                                *vec = *pr;
                                *index = start1 + ir[k];
                                count++;
			}
		}
	} else {
		ErrMsgTxt("Can't figure out this matrix.");
	}

        return(count);
}


int GetString(pMatrix ppm, int element, char *buf, int size, int ShowError)
{
        rMatrix pm = ppm[element];

	if (!mxIsChar(pm)) {
                if (ShowError)
                	ErrMsgTxt("Expecting a character element.");
                return(FALSE);
        }
	mxGetString(pm, buf, size);
        return(TRUE);
}

strArray GetCellCharItems(pMatrix ppm, int element, int len)
{
        int m, n, i;
        rMatrix pm = ppm[element];
        pMatrix pa0, **pa;

        if (!mxIsCell(pm))
                ErrMsgTxt("Expecting a cell argument.");

        m = GetM(pm);
        n = GetN(pm);
        if (!(((m == 1) && (n == len)) || ((n == 1) && (m == len))))
                ErrMsgTxt("invalid vector.");
        pa = pa0 = (pMatrix) matCalloc(len, sizeof(*pa));
        for (i = 0; i < len; i++) {
        	*pa = mxGetCell(pm, i);
        	if (!mxIsChar(*pa))
                        break;
                pa++;
        }
        if (i < len) {
                matFree(pa0);
        	ErrMsgTxt("Expecting a character cell element.");
        }
        return(pa0);
}


void FreeCellCharItems(strArray pa, int len)
{
	matFree(pa);
}

double *CreateDoubleMatrix(int m, int n, rMatrix plhs[], int element)
{
	plhs[element] = mxCreateDoubleMatrix(m, n, mxREAL);
	return(GetPr(plhs[element]));
}

double *CreateDoubleSparseMatrix(int m, int n, rMatrix plhs[], int element)
{
        int nzmax;

        if (m > n)
                nzmax = m;
        else
                nzmax = n;
        if (nzmax < 10)
                nzmax = 10;

        plhs[element] = mxCreateSparse(m, n, nzmax, mxREAL);

        return(NULL);
}

void SetColumnDoubleSparseMatrix(rMatrix plhs[], int element, int m, int n, double *mat, int column, double *arry, int *index, int size, int *nz)
{
        int *jcs, *irs, i, ii, nzmax;
        double *sr, a;

        jcs = mxGetJc(plhs[element]);
        irs = mxGetIr(plhs[element]);
        sr  = mxGetPr(plhs[element]);
        if (m > n)
                nzmax = m;
        else
                nzmax = n;
        if (nzmax < 10)
                nzmax = 10;

        jcs[column - 1] = *nz;
        for (i = 0; (i < size); i++) {
                a = arry[i];
                if (a) {
                        if (index == NULL)
                                ii = i;
                        else
                                ii = index[i] - 1;
                        if ((*nz != 0) && (((*nz) % nzmax) == 0)) {
                                mxSetNzmax(plhs[element], *nz + nzmax);
                                mxSetPr(plhs[element], matRealloc(sr, (*nz + nzmax) * sizeof(double)));
                                mxSetIr(plhs[element], matRealloc(irs, (*nz + nzmax) * sizeof(int)));
                                jcs = mxGetJc(plhs[element]);
                                irs = mxGetIr(plhs[element]);
                                sr  = mxGetPr(plhs[element]);
                        }
                        sr[*nz] = a;
                        irs[*nz] = ii;
                        (*nz)++;
                }
        }
        jcs[column] = *nz;
}

void CreateString(char **str, int n, rMatrix plhs[], int element)
{
        rMatrix pa;

        if (n == 1) {
         	pa = mxCreateString(str[0]);
                if (plhs != NULL)
                	plhs[element] = pa;
        }
        else {
                int i;
        	rMatrix pm;

                pm = plhs[element] = mxCreateCellMatrix(1, n);
                for (i = 0; i < n; i++) {
                        pa = mxCreateString(str[i]);
                	mxSetCell(pm, i, pa);
        	}
        }
}
