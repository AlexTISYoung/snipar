import time
import logging
import unittest
import os
import warnings
from tempfile import TemporaryDirectory

tests_root = os.path.dirname(__file__)
_test_directory = TemporaryDirectory(prefix="snipar-tests-")
output_root = _test_directory.name

class SniparTest(unittest.TestCase):
    p_value_threshold = 0.01
    subsample_snp = 50
    log = False
    @classmethod
    def setUpClass(cls):        
        logging.basicConfig(level=logging.ERROR, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
        warnings.filterwarnings("ignore")

    @classmethod
    def tearDownClass(cls):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print('%s: %.3f' % (self.id(), t))