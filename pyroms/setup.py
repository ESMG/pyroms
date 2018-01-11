#!/usr/bin/env python

"""pyroms is a suite of tools for working with ROMS.

Requires:
    NumPy (http://numpy.scipy.org)
    matplotlib with the Basemap toolkit (http://matplotlib.sourceforge.net)
    netCDF4 (http://code.google.com/p/netcdf4-python/)

Contains:
    hgrid  -  Tools for dealing with curvilinear grids
      BoundaryInteractor
      _Focus_x
      _Focus_y
       Focus
      CGrid
      CGrid_geo
      Gridgen
      edit_mask_mesh
      get_position_from_map

    vgrid  -  Various vertical coordinates
      s_coordinate
      z_r
      z_w
      z_coordinate

    grid  -  ROMS Grid object
      ROMS_Grid
      ROMS_gridinfo
      print_ROMS_gridinfo
      list_ROMS_gridid
      get_ROMS_hgrid
      get_ROMS_vgrid
      get_ROMS_grid
      write_ROMS_grid

    io  -  wrapper for netCDF4
      Dataset
      MFDataset

    cf  -  CF compliant files tools
      time

    tools  -  Tools specific to the Regional Ocean Modeling System
      roms2z
      z2roms
      zslice
      sslice
      islice
      jslice
      isoslice
      section
      lonslice
      latslice
      section_transport

    utility  -  Some basic tools
      get_lonlat
      get_ij
      roms_varlist
      get_roms_var
      get_bottom
      get_surface
"""

classifiers = """\
Development Status :: beta
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: MIT
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries :: Python Modules
"""

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pyroms', parent_package, top_path)
    config.add_extension('_interp', sources = ['pyroms/src/interp.f'])
    config.add_extension('_obs_interp', sources = ['pyroms/src/obs_interp.f'])
    config.add_extension('_remapping', sources = ['pyroms/src/remapping.f'])
    config.add_extension('_remapping_fast', sources = ['pyroms/src/remapping_fast.f'])
    config.add_extension('_remapping_fast_weighted', sources = ['pyroms/src/remapping_fast_weighted.f'])
    config.add_extension('_iso', sources = ['pyroms/src/iso.f'])
    return config

doclines = __doc__.split("\n")

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "pyroms",
          version = '0.1.0',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          url = "https://github.com/ESMG/",
          packages = ['pyroms',
                      'pyroms.remapping',
                      'pyroms.extern'],
          license = 'BSD',
          platforms = ["any"],
          )
    setup(**configuration(top_path='').todict())
