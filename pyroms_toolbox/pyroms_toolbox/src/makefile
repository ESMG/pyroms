#!/bin/csh

F2PY = f2py
FCOMPILER = pg 
F2PY_LIB = 
F2PY_INCLUDE = 
F2PY_F90_FLAG= '-r8'
SRCDIR  = .
EXEDIR  = ..

all: $(SRCDIR)/_average.so

_average.so: $(SRCDIR)/average.f90
	$(F2PY) --fcompiler=$(FCOMPILER) --f90flags=$(F2PY_F90_FLAG) \
	$(F2PY_LIB) $(F2PY_INCLUDE) -m _average \
	-c average.f90

install:
	/bin/mv *.so ..

clean:
	/bin/rm *.so
