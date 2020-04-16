import os
os.system("python tests/impute_from_sibs_setup.py build_ext --inplace")
os.system("python sibreg/bin/impute_from_sibs_setup.py build_ext --inplace")
os.system("cp impute_from_sibs.cpython-37m-x86_64-linux-gnu.so sibreg/bin/")
from tests.test_generated import *
from tests.test_sibreg import *
# from tests.test_impute_from_sibs import *
from tests.test_commandline import *
from tests.test_pedigree_creation import *

if __name__ == '__main__':
    unittest.main()
