import time
import logging
import unittest
import os
import warnings

tests_root = os.path.dirname(__file__)
output_root = os.path.join(tests_root, "tmp")
if not os.path.exists(output_root):
    os.mkdir(output_root)

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