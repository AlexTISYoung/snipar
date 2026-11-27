import os
import sys
from pathlib import Path

import numpy
from setuptools import Extension, setup


this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

# The native imputation loop uses OpenMP locks as well as parallel loops.
if sys.platform == "win32":
    openmp_options = dict(extra_compile_args=["/openmp"])
elif sys.platform == "darwin":
    prefix = os.environ.get("SNIPAR_OPENMP_PREFIX")
    if not prefix:
        for candidate in (Path(sys.prefix), Path("/opt/homebrew/opt/libomp"),
                          Path("/usr/local/opt/libomp")):
            if (candidate / "include/omp.h").is_file():
                prefix = str(candidate)
                break
    openmp_options = dict(extra_compile_args=["-Xpreprocessor", "-fopenmp"],
                          libraries=["omp"])
    if prefix:
        openmp_options.update(library_dirs=[str(Path(prefix) / "lib")],
                              runtime_library_dirs=[str(Path(prefix) / "lib")])
else:
    openmp_options = dict(extra_compile_args=["-fopenmp"], extra_link_args=["-fopenmp"])

include_dirs = [numpy.get_include()]
if sys.platform == "darwin" and prefix:
    include_dirs.append(str(Path(prefix) / "include"))

install_requires = [
    "numpy>=2.1.2,<3",
    "scipy>=1.14.1,<2",
    # cbgen wheels cover these Python/platform combinations. Elsewhere the
    # native bgen reader supplies the same score, GWAS and imputation inputs.
    'bgen-reader>=4.0.9,<5; python_version < "3.13" and '
    '(python_version < "3.12" or sys_platform != "linux" or platform_machine != "aarch64")',
    'bgen>=1.11,<2; python_version >= "3.13" or '
    '(python_version >= "3.12" and sys_platform == "linux" and platform_machine == "aarch64")',
    "pandas>=2.2.3,<3",
    "pysnptools>=1.0.1,<2",
    "networkx>=3.3,<4",
    "h5py>=3.12,<4",
    'numba>=0.61,<0.63; sys_platform == "darwin" and platform_machine == "x86_64"',
    'numba>=0.61,<1; sys_platform != "darwin" or platform_machine != "x86_64"',
    "statsmodels>=0.14.4,<1",
    "numdifftools>=0.9.41,<1",
    "tqdm>=4.66.5",
]

setup(name='snipar',
      version='0.0.24',
      description='Library and command line scripts for inferring identity-by-descent (IBD) segments shared between siblings, imputing missing parental genotypes, performing family based genome-wide association and polygenic score analyses.',
      long_description = long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/AlexTISYoung/snipar',
      download_url='https://pypi.org/project/snipar/files',
      author='Alexander Strudwick Young, Moeen Nehzati, Junming Guan',
      author_email='alextisyoung@gmail.com',
      license='MIT',
      include_package_data=True,
      package_data={'': ['*.pxd', '*.pyx']},
      scripts=['snipar/scripts/gwas.py', 'snipar/scripts/pgs.py', 'snipar/scripts/pgs.R', 'snipar/scripts/impute.py', 'snipar/scripts/ibd.py','snipar/scripts/correlate.py','snipar/example/snipar_example_data.py','snipar/scripts/simulate.py'],
      classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 4 - Beta',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      python_requires='>=3.11.5,<3.15',
      keywords='statistics genetics',
      packages=['snipar', 'snipar.imputation', 'snipar.read', 'snipar.tests', 'snipar.example', 'snipar.scripts'],
      install_requires=install_requires,
      extras_require={'legacy-bgen': ['bgen-reader>=4.0.9,<5']},
      zip_safe=False,
      ext_modules=[
                Extension("snipar.imputation.impute_from_sibs",
                    ["snipar/imputation/impute_from_sibs.pyx"],
                    language='c++',
                    include_dirs=include_dirs,
                    **openmp_options,
			     ),
                Extension("snipar.tests.test_impute_from_sibs",
                    ["snipar/tests/test_impute_from_sibs.pyx"],
                    include_dirs=[numpy.get_include()],
                    language='c++',
			    ),
                  ],
)
