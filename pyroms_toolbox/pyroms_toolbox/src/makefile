#!/bin/csh

F2PY = f2py
FCOMPILER = gfortran
F2PY_LIB = 
F2PY_INCLUDE = 
F2PY_F90_FLAG= '-fdefault-real-8'
SRCDIR  = .
EXEDIR  = ..

all: $(SRCDIR)/_average.so $(SRCDIR)/_move_runoff.so $(SRCDIR)/creep.so \
       $(SRCDIR)/_move_river_t.so

_average.so: $(SRCDIR)/average.f90
	$(F2PY) --fcompiler=$(FCOMPILER) --f90flags=$(F2PY_F90_FLAG) \
	$(F2PY_LIB) $(F2PY_INCLUDE) -m _average \
	-c average.f90

_move_runoff.so: $(SRCDIR)/move_runoff.f90
	$(F2PY) --fcompiler=$(FCOMPILER) --f90flags=$(F2PY_F90_FLAG) \
	$(F2PY_LIB) $(F2PY_INCLUDE) -m _move_runoff \
	-c move_runoff.f90

_move_river_t.so: $(SRCDIR)/move_river_t.f90
	$(F2PY) --fcompiler=$(FCOMPILER) --f90flags=$(F2PY_F90_FLAG) \
	$(F2PY_LIB) $(F2PY_INCLUDE) -m _move_river_t \
	-c move_river_t.f90

creep.so: $(SRCDIR)/creeping_sea.f90
	$(F2PY) -c -m creep creeping_sea.f90

install:
	/bin/mv *.so ..

clean:
	/bin/rm *.so
