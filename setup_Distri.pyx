from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = "Distribution_vel",	
  ext_modules = cythonize("Distribution_vel.pyx"),
  include_dirs=[numpy.get_include()],
)
