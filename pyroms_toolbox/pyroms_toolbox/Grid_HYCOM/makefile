# Makefile for my utils
# R. Dussin, 2013

F2PY=f2py

all: creep.so

creep.so: creeping_sea.f90
	$(F2PY) -c -m creep creeping_sea.f90

clean:
	@rm *.so
