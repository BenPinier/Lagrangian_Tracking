from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = "Pos2Flu",	
  ext_modules = cythonize("Pos2Flu.pyx"),
  include_dirs=[numpy.get_include()]
)
