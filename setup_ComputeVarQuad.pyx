from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = "Compute_VarQuad",	
  ext_modules = cythonize("Compute_VarQuad.pyx"),
  include_dirs=[numpy.get_include()]
)
