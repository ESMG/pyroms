#!/usr/bin/env python

"""
pyroms_toolbox is a suite of tools for working with ROMS.

Requires:
    pyroms (https://github.com/ESMG/pyroms)

Contains:
    many things...

"""

from numpy.distutils.core import Extension

average = Extension(name = '_average',
                    sources = ['pyroms_toolbox/src/average.f90'])

creep = Extension(name = 'creep',
                  sources = ['pyroms_toolbox/src/creeping_sea.f90'])

move_river = Extension(name = '_move_river_t',
                       sources = ['pyroms_toolbox/src/move_river_t.f90'])

move_runoff = Extension(name = '_move_runoff',
                        sources = ['pyroms_toolbox/src/move_runoff.f90'])


doclines = __doc__.split("\n")

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None,parent_package,top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True)
#                       quiet=True)
    config.add_subpackage('pyroms_toolbox')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "pyroms_toolbox",
          version = '0.2.0',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          author = "ESMG",
          url = 'https://github.com/ESMG/pyroms',
          license = 'BSD',
          platforms = ["any"],
          ext_modules=[average, creep, move_river, move_runoff],
          configuration=configuration,
          )
