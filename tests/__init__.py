import os
os.system("python tests/impute_from_sibs_setup.py build_ext --inplace")

from test_generated import *
from test_sibreg import *
from test_impute_from_sibs import *
if __name__ == '__main__':
    unittest.main()
