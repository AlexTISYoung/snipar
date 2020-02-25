from distutils.core import setup
from Cython.Build import cythonize

setup(name='sib imputation',
      ext_modules=cythonize("tests/test_impute_from_sibs.pyx", language='c++'))
      