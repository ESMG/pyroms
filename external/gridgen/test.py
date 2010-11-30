from numpy import *
from ctypes import *

pydll.LoadLibrary("./libgridgen.so")
libgridgen = CDLL("./libgridgen.so")
libgridgen.gridgen_generategrid2.restype = c_void_p
libgridgen.gridnodes_getx.restype = POINTER(POINTER(c_double))
libgridgen.gridnodes_gety.restype = POINTER(POINTER(c_double))
libgridgen.gridnodes_getnce1.restype = c_int
libgridgen.gridnodes_getnce2.restype = c_int
libgridgen.gridmap_build.restype = c_void_p

NULL = POINTER(c_int)()

nsigmas = c_int(0)
sigmas = c_void_p(0)
nrect = c_int(0)
xrect =  c_void_p(0)
yrect = c_void_p(0)

NBDRY = 13
nbdry = c_int(NBDRY)
bdry_t = c_double * NBDRY
xbdry = bdry_t(0, 7, 7, 6, 4, 2, 1, 1, 4, 7, 7, 0, 0)
ybdry = bdry_t(0, 0, 2, 2, 8, 2, 2, 7, 11, 7, 12, 12, 0)
beta = bdry_t(0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0)
ul = c_int(2)
nx = c_int(10)
ny = c_int(5)
ngrid = c_int(0)
xgrid = POINTER(c_double)()
ygrid = POINTER(c_double)()
nnodes = c_int(12)
newton = c_int(1)
precision = c_double(1.0e-11)
checksimplepoly = c_int(1)
thin = c_int(1)
nppe = c_int(3)
verbose = c_int(1)

gn = libgridgen.gridgen_generategrid2(nbdry, xbdry, ybdry, beta, ul, nx, ny, ngrid, xgrid, ygrid, nnodes, newton, precision, checksimplepoly, thin, nppe, verbose, byref(nsigmas), byref(sigmas), byref(nrect), byref(xrect), byref(yrect))
#
# One may use generated grid now, e.g access node coordinates:
#
x = libgridgen.gridnodes_getx(gn)
y = libgridgen.gridnodes_gety(gn)
nx = libgridgen.gridnodes_getnx(gn);
ny = libgridgen.gridnodes_getny(gn);
for j in range(ny):
    for i in range(nx):
        print '  [', j, '][', i, ']: x =', x[j][i], ', y=', y[j][i]
#
# or print to a file:
#
libgridgen.gridnodes_write(gn, "test.txt", c_int(2))
#
# or map between index and physical space:
#
gm = libgridgen.gridmap_build(nx - 1, ny - 1, x, y);
xx = c_double();
yy = c_double();
fi = c_double(0.9);
fj = c_double(0.9);
libgridgen.gridmap_fij2xy(gm, fi, fj, byref(xx), byref(yy));
print '  i =', fi.value, ', j =', fj.value, ' -> x =', xx.value, ', y=', yy.value
libgridgen.gridmap_xy2fij(gm, xx, yy, byref(fi), byref(fj));
print '  x =', xx.value, ', y=', yy.value, ' ->  i =', fi.value, ', j =', fj.value
#
# Clean-up
#
libgridgen.gridmap_destroy(gm)
libgridgen.gridnodes_destroy(gn)
#
# (One may analyse/draw image stored in xrect, yrect now)
#

#
# This time gridgen should use the calculated values of sigmas
# (there should be only one dot in a log after "solving for sigmas:")
#
nx = 100
ny = 20
gn = libgridgen.gridgen_generategrid2(nbdry, xbdry, ybdry, beta, ul, nx, ny, ngrid, xgrid, ygrid, nnodes, newton, precision, checksimplepoly, thin, nppe, verbose, byref(nsigmas), byref(sigmas), byref(nrect), byref(xrect), byref(yrect))
libgridgen.gridnodes_write(gn, "test-1.txt", c_int(2))
libgridgen.gridnodes_destroy(gn)


nx = 150
ny = 25
gn = libgridgen.gridgen_generategrid2(nbdry, xbdry, ybdry, beta, ul, nx, ny, ngrid, xgrid, ygrid, nnodes, newton, precision, checksimplepoly, thin, nppe, verbose, NULL, NULL, NULL, NULL)
libgridgen.gridnodes_write(gn, "test-2.txt", c_int(2))
libgridgen.gridnodes_destroy(gn)
