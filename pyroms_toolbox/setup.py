#!/usr/bin/env python

"""
pyroms_toolbox is a suite of tools for working with ROMS.

Requires:
    pyroms (https://github.com/kshedstrom/pyroms)

Contains:
    many things...

"""

doclines = __doc__.split("\n")

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pyroms_toolbox',parent_package,top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True)
#                       quiet=True)
    config.add_subpackage('pyroms_toolbox.BGrid_GFDL')
    config.add_subpackage('pyroms_toolbox.BGrid_POP')
    config.add_subpackage('pyroms_toolbox.BGrid_SODA')
    config.add_subpackage('pyroms_toolbox.seawater')
    config.add_subpackage('pyroms_toolbox.SODA')
    config.add_library('_average', sources=['src/average.f90']),
    config.add_library('_move_runoff', sources=['src/move_runoff.f90']),
    config.add_extension('_average',
          sources = ['src/average.f90'],
          libraries = ['_average']
          )
    config.add_extension('_move_runoff',
          sources = ['src/move_runoff.f90'],
          libraries = ['_move_runoff']
          )
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "pyroms_toolbox",
          version = '0.1',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          author = "Frederic Castruccio",
          author_email = "frederic@marine.rutgers.edu",
          url = 'https://github.com/kshedstrom/pyroms',
          license = 'BSD',
          platforms = ["any"],
          configuration=configuration,
          )
