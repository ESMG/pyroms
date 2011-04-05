from distutils.core import setup, Extension
from os import getenv
windir = getenv('windir')
if windir == None:
  WIN32 = 'NOWIN32'
else:
  WIN32 = 'WIN32'

lpsolve55 = Extension(name = 'lpsolve55',
                sources = ['lpsolve.c', 'hash.c', 'pythonmod.c'],
                define_macros=[('PYTHON', '1'), (WIN32, '1'), ('NODEBUG', '1')],
                include_dirs=['../..'],
                library_dirs=['/u1/uaf/kate/lib/'],
                libraries=['lpsolve55'])

setup (name = "lpsolve55",
       version = "5.5.0.4",
       description = "Linear Program Solver, Interface to lpsolve",
       author = "Peter Notebaert",
       author_email = "lpsolve@peno.be",
       url = "http://www.peno.be/",
       py_modules=['lp_solve', 'lp_maker'],
       ext_modules = [lpsolve55],
       )
