from distutils.core import setup
from Cython.Build import cythonize
import numpy
setup(name='sib imputation',
      ext_modules=cythonize("sibreg/bin/impute_from_sibs.pyx", annotate = True, language='c++', include_path=[numpy.get_include()]),
      include_dirs=[numpy.get_include()]
      )