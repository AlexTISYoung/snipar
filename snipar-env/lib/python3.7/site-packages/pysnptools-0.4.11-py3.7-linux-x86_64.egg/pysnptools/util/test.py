import os
import numpy as np
import logging
import doctest
import unittest

class TestUtilTools(unittest.TestCase):
    def test_util(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        import pysnptools.util.pheno
        import pysnptools.util.mapreduce1.examples
        from pysnptools.util.mapreduce1 import map_reduce #Needed to work around thread local variable issue
        for mod in [
                    pysnptools.util.pheno,
                    pysnptools.util.mapreduce1.examples,
                    pysnptools.util.mapreduce1.mapreduce,
                    pysnptools.util.mapreduce1.runner.localmultiproc,
                    pysnptools.util.mapreduce1.runner.localmultithread,
                    ]:
            result = doctest.testmod(mod,optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
            assert result.failed == 0, "failed doc test: " + __file__
        os.chdir(old_dir)

    def test_util_mapreduce1_testmod(self):
        import pysnptools.util.mapreduce1
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/mapreduce1")
        for mod in [
                    pysnptools.util.mapreduce1.runner.local,
                    pysnptools.util.mapreduce1.runner.localinparts,
                    pysnptools.util.mapreduce1.runner.localmultiproc,
                    pysnptools.util.mapreduce1.runner.localmultithread,
                    ]:
            result = doctest.testmod(mod,optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
            assert result.failed == 0, "failed doc test: " + __file__
        os.chdir(old_dir)


def getTestSuite():
    """
    set up composite test suite
    """

    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestUtilTools))
    return test_suite

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    r.run(suites)
