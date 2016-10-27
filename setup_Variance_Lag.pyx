from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = "VarLag",	
  ext_modules = cythonize("Variance_Lag.pyx"),
  include_dirs=[numpy.get_include()]
)
