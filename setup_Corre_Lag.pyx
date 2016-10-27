from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = "Correlation_Vel_lag",	
  ext_modules = cythonize("Correlation_Vel_lag.pyx"),
  include_dirs=[numpy.get_include()]
)
