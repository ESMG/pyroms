/* config.h.  Generated from config.h.in by configure.  */
#define HAVE_TGAMMA 1
/* #undef HAVE_LGAMMA */
#define HAVE_COMPLEX 1
/* #undef HAVE__COMPLEX__ */
/* #undef NOCOMPLEX */
#define HAVE_GRIDNODES_H 1

#if defined(HAVE_COMPLEX)
 #include <complex.h>
 #define zdouble _Complex double
#elif defined(HAVE__COMPLEX__)
 #include <complex.h>
 #define zdouble __complex__ double
#elif defined(NOCOMPLEX)
 #include "c99-min/complex.h"
 #define zdouble _Complex double
#endif
