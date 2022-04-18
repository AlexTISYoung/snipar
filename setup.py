from setuptools import setup, find_packages, Extension
from setuptools import dist
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()
class MyExt(Extension):
    def __init__(self, *args, **kwargs):
        self.__include_dirs = ["."]
        super().__init__(*args, **kwargs)

    @property
    def include_dirs(self):
        import numpy
        return self.__include_dirs + [numpy.get_include()]

    @include_dirs.setter
    def include_dirs(self, dirs):
        self.__include_dirs = dirs

setup(name='snipar',
      version='0.0.8',
      description='Library and command line scripts for inferring identity-by-descent (IBD) segments shared between siblings, imputing missing parental genotypes, and for performing family based genome-wide association and polygenic score analyses.',
      long_description = long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/alexTISYoung/snipar',
      download_url='https://pypi.org/project/snipar/files',
      author='Alexander I. Young, Moeen Nehzati',
      author_email='alextisyoung@gmail.com',
      license='MIT',
      include_package_data=True,
      package_data={'': ['*.pxd', '*.pyx']},
      scripts=['snipar/scripts/gwas.py', 'snipar/scripts/pgs.py', 'snipar/scripts/impute.py', 'snipar/scripts/ibd.py','snipar/scripts/correlate.py','snipar/example/snipar_example_data.py'],
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

      ],
      python_requires='>=3.7',
      keywords='statistics genetics',
      packages=['snipar', 'snipar.imputation', 'snipar.read', 'snipar.tests', 'snipar.example', 'snipar.scripts'],
      install_requires=[
            'numpy==1.21.1',
            'scipy==1.7.1',
            'bgen_reader==4.0.7',
            'pandas==1.1.4',
            'pysnptools==0.5.3',
            'networkx==2.5',
            'h5py==3.6.0',
            'pooch==1.5.1',
            'numba==0.55.0',
            'gitpython==3.1.24',
            'scikit-learn==1.0.2',
            'statsmodels==0.13.2',
            ],
      test_suite="snipar/tests",
      zip_safe=False,
      ext_modules=[
                MyExt("snipar.imputation.impute_from_sibs",
                    ["snipar/imputation/impute_from_sibs.pyx"],
                    language='c++',
                    extra_compile_args=['-fopenmp'],
                    extra_link_args=['-fopenmp'],
			     ),
                MyExt("snipar.tests.test_impute_from_sibs",
                    ["snipar/tests/test_impute_from_sibs.pyx"],
                    language='c++',
			    ),
                  ],
)
