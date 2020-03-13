from distutils.core import setup
from Cython.Build import cythonize

setup(name='sib imputation',
      ext_modules=cythonize("sibreg/bin/cython_impute_from_sibs.pyx", annotate = True, language='c++'))