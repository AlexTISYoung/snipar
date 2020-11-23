from setuptools import setup
from Cython.Build import cythonize

setup(
    name='logll',
    ext_modules=cythonize("sib_ldsc_cython.pyx")
)