#!/bin/csh
#
# Makefile for interpolation code
#
# CVS:$Id: makefile,v 1.7 2000/04/19 21:46:44 pwjones Exp $
#
#===============================================================================
# SVN $Id$
# SVN $URL$
#===============================================================================

NC_CONFIG = nf-config
LIBDIR = $(shell $(NC_CONFIG) --prefix)/lib
INCDIR = $(shell $(NC_CONFIG) --prefix)/include
SRCDIR  = .
PREFIX  = /usr/local

COMPILE = gfortran
FLAGS = -g -fdefault-real-8 -ffixed-form -O2 -fPIC
LIB  =  $(shell $(NC_CONFIG) --flibs)
INCLUDE = -I$(INCDIR)
F2PY = f2py
F2PY_LIB = -L${LIBDIR} -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl
F2PY_INCLUDE = -I${INCDIR}
F2PY_F90_FLAG= '-g -fdefault-real-8 -ffixed-form -O2 -fPIC'

OBJSET  = \
	kinds_mod.o \
	constants.o \
	iounits.o \
	netcdf.o \
	grids.o \
	remap_vars.o \
	remap_distwgt.o \
	remap_conserv.o \
	remap_bilinear.o \
	remap_bicubic.o \
	timers.o \
	remap_write.o \
	scrip.o

OBJTEST  = \
	kinds_mod.o \
	constants.o \
	iounits.o \
	netcdf.o \
	grids.o \
	timers.o \
	remap_vars.o \
	remap_read.o \
	remap.o

OBJF2PY  = \
        kinds_mod.o \
        constants.o \
        iounits.o \
        netcdf.o \
        grids.o \
        remap_vars.o \
        remap_distwgt.o \
        remap_conserv.o \
        remap_bilinear.o \
        remap_bicubic.o \
        remap_read.o \
        remap.o \
        timers.o \
        remap_write.o \


all: $(SRCDIR)/scrip $(SRCDIR)/scrip_test

$(SRCDIR)/scrip: $(OBJSET)
	$(COMPILE) $(FLAGS) $(OBJSET) $(LIB) $(INCLUDE) -o $(SRCDIR)/scrip

$(SRCDIR)/scrip_test: $(OBJTEST) scrip_test.o
	$(COMPILE) $(FLAGS) $(OBJTEST) scrip_test.o $(LIB) \
	-o $(SRCDIR)/scrip_test

$(SRCDIR)/scrip_test_repeat: $(OBJTEST) scrip_test_repeat.o
	$(COMPILE) $(FLAGS) $(OBJTEST) scrip_test_repeat.o $(LIB) \
	-o $(SRCDIR)/scrip_test_repeat

kinds_mod.o: $(SRCDIR)/kinds_mod.f
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/kinds_mod.f

constants.o: $(SRCDIR)/constants.f kinds_mod.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/constants.f

iounits.o: $(SRCDIR)/iounits.f kinds_mod.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/iounits.f

netcdf.o: $(SRCDIR)/netcdf.f kinds_mod.o constants.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/netcdf.f

grids.o: $(SRCDIR)/grids.f kinds_mod.o constants.o iounits.o netcdf.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/grids.f

remap_vars.o: $(SRCDIR)/remap_vars.f kinds_mod.o constants.o grids.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/remap_vars.f

remap_conserv.o: $(SRCDIR)/remap_conserv.f kinds_mod.o constants.o \
		timers.o remap_vars.o grids.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/remap_conserv.f

remap_distwgt.o: $(SRCDIR)/remap_distwgt.f kinds_mod.o constants.o \
		remap_vars.o grids.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/remap_distwgt.f

remap_bilinear.o: $(SRCDIR)/remap_bilinear.f kinds_mod.o constants.o \
		remap_vars.o grids.o timers.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/remap_bilinear.f

remap_bicubic.o: $(SRCDIR)/remap_bicubic.f kinds_mod.o constants.o \
		remap_vars.o grids.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/remap_bicubic.f

timers.o: $(SRCDIR)/timers.f kinds_mod.o constants.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/timers.f

remap_write.o: $(SRCDIR)/remap_write.f kinds_mod.o constants.o \
		netcdf.o remap_vars.o grids.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/remap_write.f

remap_read.o: $(SRCDIR)/remap_read.f kinds_mod.o constants.o netcdf.o \
		remap_vars.o grids.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/remap_read.f

remap.o: $(SRCDIR)/remap.f kinds_mod.o constants.o
	$(COMPILE) $(FLAGS) -c $(SRCDIR)/remap.f

scrip.o: $(SRCDIR)/scrip.f kinds_mod.o constants.o iounits.o timers.o \
		remap_vars.o grids.o remap_conserv.o remap_distwgt.o \
		remap_bilinear.o remap_bicubic.o remap_write.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/scrip.f

scrip_test.o: $(SRCDIR)/scrip_test.f kinds_mod.o constants.o iounits.o \
		netcdf.o remap_vars.o grids.o remap.o remap_read.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/scrip_test.f

scrip_test_repeat.o: $(SRCDIR)/scrip_test_repeat.f kinds_mod.o \
		constants.o netcdf.o iounits.o remap_vars.o grids.o
	$(COMPILE) $(FLAGS) $(INCLUDE) -c $(SRCDIR)/scrip_test_repeat.f


f2py: $(SRCDIR)/scrip.so

$(SRCDIR)/scrip.so: $(SRCDIR)/pyscrip.f90 $(OBJF2PY)
	$(F2PY) --fcompiler=gnu95 --f90flags=$(F2PY_F90_FLAG) \
	$(F2PY_LIB) $(F2PY_INCLUDE) -m scrip \
	-c $(OBJF2PY) pyscrip.f90

install:
	/bin/cp $(SRCDIR)/scrip $(PREFIX)/bin
	/bin/cp $(SRCDIR)/scrip_test $(PREFIX)/bin
	/bin/cp $(SRCDIR)/scrip*.so $(PREFIX)/lib

clean:
	/bin/rm *.o *.mod *.so scrip scrip_test
