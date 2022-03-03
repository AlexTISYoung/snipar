from setuptools import setup, find_packages, Extension
from setuptools import dist

class MyExt(Extension):
    def __init__(self, *args, **kwargs):
        self.__include_dirs = []
        super().__init__(*args, **kwargs)

    @property
    def include_dirs(self):
        import numpy
        return self.__include_dirs + [numpy.get_include()]

    @include_dirs.setter
    def include_dirs(self, dirs):
        self.__include_dirs = dirs

setup(name='snipar',
      # version='1.2.0a1',
      # description='Functions for performing robust GWAS using sibpairs in a random effects model',
      # url='http://github.com/alexTISYoung/sibreg',
      # download_url='https://github.com/AlexTISYoung/hlmm/archive/1.2.0a1.tar.gz',
      # author='Alexander I. Young',
      # author_email='alextisyoung@gmail.com',
      # license='MIT',
      include_package_data=True,
      scripts=['gwas.py', 'pgs.py', 'impute.py', 'ibd.py','correlate.py'],
      classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 4 - Beta',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: MIT License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 3.7',
      ],
      keywords='statistics genetics',
      packages=['snipar', 'snipar.imputation','snipar.read','snipar.tests','snipar.example'],
      setup_requires=['numpy==1.19.3', 'Cython==0.29.21'],
      install_requires=[
            'bgen_reader==4.0.7',
            'pandas==1.1.1',
            'numpy==1.19.3',
            'Cython==0.29.25',
            'scipy==1.7.1',
            'pysnptools==0.5.3',
            'networkx==2.2',
            'h5py==2.10.0',
            'pooch==1.5.1',
            'numba==0.50.0',
            'gitpython==3.1.24',
            'scikit-learn==1.0.2',
            'statsmodels==0.13.1',
            ],
      test_suite="snipar/tests",
      zip_safe=False,
      ext_modules=[MyExt("snipar.imputation.impute_from_sibs",
			     ["snipar/imputation/impute_from_sibs.pyx"],
			     language='c++',
			     extra_compile_args=['-fopenmp'],
			     extra_link_args=['-fopenmp'],
			     ),
                  MyExt("snipar.tests.test_impute_from_sibs",
			     ["snipar/tests/test_impute_from_sibs.pyx"],
			     language='c++',
			     extra_compile_args=['-fopenmp'],
			     extra_link_args=['-fopenmp'],
			     ),
                  ],
)
