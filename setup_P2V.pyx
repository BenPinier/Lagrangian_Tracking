from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = "Pos2Vel",	
  ext_modules = cythonize("Pos2Vel.pyx"),
  include_dirs=[numpy.get_include()],
)
