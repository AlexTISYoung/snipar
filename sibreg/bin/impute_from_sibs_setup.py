"""Setup for just compiling the cython code

Just run "python sibreg/bin/impute_from_sibs.py build_ext --inplace"
"""

from distutils.core import setup
from Cython.Build import cythonize
import numpy
setup(name='sib imputation',
      ext_modules=cythonize("sibreg/bin/impute_from_sibs.pyx", annotate = True, language='c++'),
      include_dirs=[numpy.get_include()]
      )