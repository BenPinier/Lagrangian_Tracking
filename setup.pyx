from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = "Lag",	
  ext_modules = cythonize("Lag.pyx"),
  include_dirs=[numpy.get_include()]
)
