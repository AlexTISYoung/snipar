from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(name='sibreg',
      version='1.2.0a1',
      description='Functions for performing robust GWAS using sibpairs in a random effects model',
      url='http://github.com/alexTISYoung/sibreg',
      download_url='https://github.com/AlexTISYoung/hlmm/archive/1.2.0a1.tar.gz',
      author='Alexander I. Young',
      author_email='alextisyoung@gmail.com',
      license='MIT',
      scripts=['bin/sGWAS.py'],
      classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 3 - Alpha',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: MIT License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 2.7',
      ],
      keywords='statistics genetics',
      packages=['sibreg'],
      install_requires=[
            'numpy',
            'scipy',
            'pysnptools',
            'panda',
            'cython'
        ],
      extras_require={
            'test': ['numdifftools'],
      },
      zip_safe=False,
      ext_modules=[cythonize("sibreg/bin/impute_from_sibs.pyx", language='c++'),
                   cythonize("tests/test_impute_from_sibs.pyx", language='c++')
                  ]
)