# encoding: utf-8
'''
PYROMS is a toolkit for working with ROMS ocean models

pyroms is based on the python/numpy/matplotlib scientific python suite.
NetCDF I/O is based on the NetCDF4-python package. The toolkit contains
general modeling tools for dealing with arrays, diagnosing standard
properties, curvilinear grid generation, and interpolation.
'''

from . import cf
from . import vgrid
from . import extern
from . import hgrid
from . import grid
from . import io
from . import sta_hgrid
from . import sta_grid
from . import tools
from . import remapping
from . import utility

__authors__ = ['Frederic Castruccio (frederic@marine.rutgers.edu)']

__version__ = '0.1.0'
