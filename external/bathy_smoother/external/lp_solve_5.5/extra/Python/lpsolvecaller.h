#include <setjmp.h>

extern void exitnow();

#if defined MATLAB
# include "matlab.h"
#elif defined OMATRIX
# include "omatrix.h"
#elif defined SCILAB
# include "scilab.h"
#elif defined OCTAVE
# include "octave.h"
#elif defined PYTHON
# include "pythonmod.h"
#else
# error "Unknown platform. Either MATLAB, OMATRIX, SCILAB, OCTAVE, PYTHON must be defined."
#endif
