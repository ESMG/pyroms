"""
bathy_smoother is a suite of tools for working with ROMS bathymetry.
(ripped from matlab script LP_bathymetry)    

Requires:
    NumPy (http://numpy.scipy.org)
    lpsolve (http://lpsolve.sourceforge.net/)

Contains:
    bathy_smoothing - Tools for smoothing the bathymetry

    bathy_tools - Various bathymetry tools

    LP_bathy_smoothing - Tools for smoothing the bathymetry using LP

    LP_bathy_tools - LP tools for smoothing the bathymetry

    LP_tools - Various LP tools

"""

from numpy.distutils.core import Extension

doclines = __doc__.split("\n")

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "bathy_smoother",
          version = '0.1',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          author = "Frederic Castruccio",
          author_email = "frederic@marine.rutgers.edu",
          packages = ['bathy_smoother'],
          license = 'BSD',
          platforms = ["any"],
          )
