from numpy.distutils.fcompiler import get_default_fcompiler
import sys
import os

fcompile = get_default_fcompiler(os.name,sys.platform,True)

install_dir = sys.argv[1]

os.system("python setup.py config_fc --fcompiler="+str(fcompile)+" build_clib --fcompiler="+str(fcompile)+" build_ext --fcompiler="+str(fcompile)+" install --prefix='"+str(install_dir)+"'")
