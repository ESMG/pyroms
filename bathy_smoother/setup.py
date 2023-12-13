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

doclines = __doc__.split("\n")

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('bathy_smoother',parent_package,top_path,
         package_path='bathy_smoother')
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True)
#                       quiet=True)
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = 'bathy_smoother',
          version = '0.2.0',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          url = 'https://github.com/ESMG/pyroms',
          license = 'BSD',
          platforms = ["any"],
          configuration=configuration,
          )
