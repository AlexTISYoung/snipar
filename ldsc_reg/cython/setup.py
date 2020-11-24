from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    name='logll',
    ext_modules=cythonize("sib_ldsc_cython.pyx"),
    include_dirs=[np.get_include()]
)