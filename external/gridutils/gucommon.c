/******************************************************************************
 *
 *  File:           gucommon.c
 *  
 *  Created         22/01/2002
 *  
 *  Author:         Pavel Sakov
 *                  CSIRO Marine Research
 *  
 *  Purpose:        Some common stuff for grid utilities
 *  Revisions:      none
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <assert.h>
#include <errno.h>
#include "version.h"

#define BUFSIZE 10240

int gu_verbose = 0;

void gu_quit(char* format, ...)
{
    va_list args;

    fflush(stdout);
    fprintf(stderr, "\nerror: gridutils: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    exit(1);
}

FILE* gu_fopen(const char* path, const char* mode)
{
    FILE* f = NULL;

    f = fopen(path, mode);
    if (f == NULL)
        gu_quit("%s: could not open for \"%s\" : %s\n", path, mode, strerror(errno));

    return f;
}

/* Allocates n1xn2 matrix of something. Note that it will be accessed as 
 * [n2][n1].
 * @param n1 Number of columns
 * @param n2 Number of rows
 * @param unitsize Size of one matrix element in bytes
 * @return Matrix
 */
void* gu_alloc2d(int n1, int n2, size_t unitsize)
{
    size_t size;
    char* p;
    char** pp;
    int i;

    if (n1 <= 0 || n2 <= 0)
        gu_quit("alloc2d(): invalid size (n1 = %d, n2 = %d)\n", n1, n2);

    size = n1 * n2;
    if ((p = calloc(size, unitsize)) == NULL)
        gu_quit("gu_alloc2d(): %s\n", strerror(errno));

    size = n2 * sizeof(void*);
    if ((pp = malloc(size)) == NULL)
        gu_quit("gu_alloc2d(): %s\n", strerror(errno));
    for (i = 0; i < n2; i++)
        pp[i] = &p[i * n1 * unitsize];

    return pp;
}

/* Destroys a matrix.
 * @param pp Matrix
 */
void gu_free2d(void* pp)
{
    void* p;

    p = ((void**) pp)[0];
    free(pp);
    free(p);
}

int** gu_readmask(char* fname, int nx, int ny)
{
    int** v = NULL;
    FILE* f = NULL;
    char buf[BUFSIZE];
    int count;
    int i, j;

    if (strcasecmp(fname, "stdin") == 0 || strcmp(fname, "-") == 0)
        f = stdin;
    else
        f = gu_fopen(fname, "r");

    v = gu_alloc2d(nx, ny, sizeof(int));

    for (j = 0, count = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            if (fgets(buf, BUFSIZE, f) == NULL)
                gu_quit("%s: could not read %d-th mask value (%d x %d values expected)\n", fname, j * nx + i + 1, nx, ny);
            buf[strlen(buf) - 1] = 0;
            if (strcmp(buf, "0") == 0)
                v[j][i] = 0;
            else if (strcmp(buf, "1") == 0) {
                v[j][i] = 1;
                count++;
            } else
                gu_quit("%s: could not interpret %d-th mask value = \"%s\" (expected \"0\" or \"1\"\n", fname, j * nx + i + 1, buf);
        }
    }

    if (f != stdin)
        fclose(f);

    if (gu_verbose) {
        int n = nx * ny;

        fprintf(stderr, "## mask: %d valid cells (%.1f%%), %d masked cells (%.1f%%)\n", count, 100.0 * count / n, n - count, 100.0 * (n - count) / n);
        fflush(stderr);
    }

    return v;
}
