#!/bin/csh
#
# Makefile for SCRIP interpolation code adapted for PyROMS (https://github.com/ESMG/pyroms).
#
# This file is based on one marked v 1.7 2000/04/19 21:46:44 pwjones.
#
# It has subsequently been incorporated into PyROMS and modified, see the history in
#   https://github.com/ESMG/pyroms/blob/python3/pyroms/external/scrip/source/makefile
#
# Recent history:
#   MGH (mark.hadf@gmail.com) 2019-10-16):
#     Renamed the netCDF-related macros and made more use of Gnu make "?=" assignments to
#     allow macros to be overridden by environment variables.
#   MGH (mark.hadf@gmail.com) 2019-10-31):
#     - Back to ":=" and "=" assignments.
#     - Renamed several more variables for greater clarity.
#     - Omitted the SCRIP and SCRIP_TEST binary executable targets.
#     - Added the DEVELOP variable to allow support for editable and non-editable
#       installations of pyroms.

SRCDIR = .
PREFIX = /usr/local

# Specify the directory into which the SCRIP extension will be copied.
# The default value depends on whether pyroms is installed in editable mode (DEVELOP=1)
ifdef DEVELOP
    SCRIP_EXT_DIR := ../../../pyroms
else
    SCRIP_EXT_DIR := $(shell python -c "import site; print(site.getsitepackages()[0])")/pyroms
endif

# The nf-config command should make linking with netCDF trivial. Unfortunately it is
# not available in Conda installations.
ifdef USE_NF_CONFIG
    NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
    NETCDF_LIBDIR ?= $(shell $(NF_CONFIG) --prefix)/lib
    NETCDF_LIBS  ?=  $(shell $(NF_CONFIG) --flibs)
else
    NETCDF_INCDIR ?= $(PREFIX)/include
    NETCDF_LIBDIR ?= $(PREFIX)/lib
    NETCDF_LIBS ?= -lnetcdf -lnetcdff -ljpeg -lmfhdf -ldf -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran
endif

FC = gfortran
F2PY = f2py
FFLAGS = -g -fdefault-real-8 -ffixed-form -O2 -fPIC
F2PY_F90_FLAGS = '-g -fdefault-real-8 -ffixed-form -O2 -fPIC'
INCLUDE = -I$(SRCDIR) -I$(NETCDF_INCDIR)
LIB = -L$(SRCDIR) -L$(NETCDF_LIBDIR) $(NETCDF_LIBS)

OBJ_SCRIP_EXT  = \
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


all: scrip_ext

kinds_mod.o: $(SRCDIR)/kinds_mod.f
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/kinds_mod.f

constants.o: $(SRCDIR)/constants.f kinds_mod.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/constants.f

iounits.o: $(SRCDIR)/iounits.f kinds_mod.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/iounits.f

netcdf.o: $(SRCDIR)/netcdf.f kinds_mod.o constants.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/netcdf.f

grids.o: $(SRCDIR)/grids.f kinds_mod.o constants.o iounits.o netcdf.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/grids.f

remap_vars.o: $(SRCDIR)/remap_vars.f kinds_mod.o constants.o grids.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/remap_vars.f

remap_conserv.o: $(SRCDIR)/remap_conserv.f kinds_mod.o constants.o \
		timers.o remap_vars.o grids.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/remap_conserv.f

remap_distwgt.o: $(SRCDIR)/remap_distwgt.f kinds_mod.o constants.o \
		remap_vars.o grids.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/remap_distwgt.f

remap_bilinear.o: $(SRCDIR)/remap_bilinear.f kinds_mod.o constants.o \
		remap_vars.o remap_distwgt.o grids.o timers.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/remap_bilinear.f

remap_bicubic.o: $(SRCDIR)/remap_bicubic.f kinds_mod.o constants.o \
		remap_vars.o grids.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/remap_bicubic.f

timers.o: $(SRCDIR)/timers.f kinds_mod.o constants.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/timers.f

remap_write.o: $(SRCDIR)/remap_write.f kinds_mod.o constants.o \
		netcdf.o remap_vars.o grids.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/remap_write.f

remap_read.o: $(SRCDIR)/remap_read.f kinds_mod.o constants.o netcdf.o \
		remap_vars.o grids.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/remap_read.f

remap.o: $(SRCDIR)/remap.f kinds_mod.o constants.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $(SRCDIR)/remap.f

scrip_ext: $(SRCDIR)/pyscrip.f90 $(OBJ_SCRIP_EXT)
	$(F2PY) --fcompiler=gnu95 --f90exec=$(FC) --f90flags=$(F2PY_F90_FLAGS) \
	$(LIB) $(INCLUDE) -m scrip -c $(OBJ_SCRIP_EXT) $(SRCDIR)/pyscrip.f90

install: all
	cp -vf scrip*.so $(SCRIP_EXT_DIR)

uninstall:
	rm -vf $(SCRIP_EXT_DIR)/scrip*.so

clean:
	rm -vf *.o *.mod *.so
