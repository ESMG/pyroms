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

    define_macros = []
    define_macros.append(('YY_NEVER_INTERACTIVE', None))
    define_macros.append(('PARSER_LP', None))
    define_macros.append(('INVERSE_ACTIVE', 'INVERSE_LUSOL'))
    define_macros.append(('RoleIsExternalInvEngine', None))

    sources = ['external/lp_solve_5.5/lp_MDO.c',
               'external/lp_solve_5.5/shared/*.c',
               'external/lp_solve_5.5/colamd/*.c',
               'external/lp_solve_5.5/bfp/bfp_LUSOL/lp_LUSOL.c',
               'external/lp_solve_5.5/bfp/bfp_LUSOL/LUSOL/lusol.c',
               'external/lp_solve_5.5/ini.c',
               'external/lp_solve_5.5/fortify.c',
               'external/lp_solve_5.5/lp_rlp.c',
               'external/lp_solve_5.5/lp_crash.c',
               'external/lp_solve_5.5/lp_Hash.c',
               'external/lp_solve_5.5/lp_lib.c',
               'external/lp_solve_5.5/lp_wlp.c',
               'external/lp_solve_5.5/lp_matrix.c',
               'external/lp_solve_5.5/lp_mipbb.c',
               'external/lp_solve_5.5/lp_MPS.c',
               'external/lp_solve_5.5/lp_params.c',
               'external/lp_solve_5.5/lp_presolve.c',
               'external/lp_solve_5.5/lp_price.c',
               'external/lp_solve_5.5/lp_pricePSE.c',
               'external/lp_solve_5.5/lp_report.c',
               'external/lp_solve_5.5/lp_scale.c',
               'external/lp_solve_5.5/lp_simplex.c',
               'external/lp_solve_5.5/lp_SOS.c',
               'external/lp_solve_5.5/lp_utils.c',
               'external/lp_solve_5.5/yacc_read.c'],
    inc_dirs = ['external/lp_solve_5.5',
                'external/lp_solve_5.5/bfp',
                'external/lp_solve_5.5/bfp/bfp_LUSOL',
                'external/lp_solve_5.5/bfp/bfp_LUSOL/LUSOL',
                'external/lp_solve_5.5/colamd',
                'external/lp_solve_5.5/shared'],
    config.add_library('lpsolve55',
          sources = sources,
          include_dirs = inc_dirs,
          macros=define_macros)
    config.add_extension('lpsolve55',
          sources = sources,
          include_dirs = inc_dirs,
          libraries = ['lpsolve55']
          )
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = '',
          version = '0.1',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          url = 'https://github.com/ESMG/pyroms',
          license = 'BSD',
          platforms = ["any"],
          configuration=configuration,
          )
