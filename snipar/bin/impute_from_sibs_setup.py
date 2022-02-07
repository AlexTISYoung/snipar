"""Setup for just compiling the cython code

Just run "python snipar/bin/impute_from_sibs.py build_ext --inplace"
"""

from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "snipar.bin.impute_from_sibs",
        ["snipar/bin/impute_from_sibs.pyx"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
        language="c++"
    )
]

setup(
    name='sib imputation',
    ext_modules=cythonize(ext_modules),
    include_dirs=[numpy.get_include()]

)
