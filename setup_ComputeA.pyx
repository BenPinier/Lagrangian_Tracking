from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = "Compute_Azz",	
  ext_modules = cythonize("Compute_A.pyx"),
  include_dirs=[numpy.get_include()]
)
