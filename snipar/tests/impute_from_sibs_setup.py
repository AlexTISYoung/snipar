from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "snipar.tests.test_impute_from_sibs",
        ["snipar/tests/test_impute_from_sibs.pyx"],
        extra_compile_args=['-Xpreprocessor', '-fopenmp'],
        extra_link_args=['-Xpreprocessor', '-fopenmp'],
        language = "c++"
    )
]

setup(
    name='test sib imputation',
    ext_modules=cythonize(ext_modules),
    include_dirs=[numpy.get_include()]

)
