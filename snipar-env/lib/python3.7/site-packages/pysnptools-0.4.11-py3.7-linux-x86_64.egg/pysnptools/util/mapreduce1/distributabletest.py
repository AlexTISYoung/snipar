from __future__ import absolute_import
from __future__ import print_function
from pysnptools.util.mapreduce1.runner import *
import logging
import unittest
import io

class DistributableTest(object) : #implements IDistributable
    '''
    This is a class for distributing any testing.
    '''
    def __init__(self,test_or_suite,tempdirectory=None):
        self._test_or_suite = test_or_suite
        self.__tempdirectory = tempdirectory
        self.name = str(test_or_suite)

    @staticmethod
    def deep_iter(test_or_suite):
        try:
            for sub in test_or_suite:
                for subsub in DistributableTest.deep_iter(sub):
                    yield subsub
        except TypeError as detail:
            yield test_or_suite

    def do_work(self, test):

        ###############test_result = unittest.TestResult()
        r = unittest.TextTestRunner(failfast=False)
        logging.info(test)
        test.setUpClass()
        test_result = r.run(test)
        ##############test.run(result=test_result)
        test.tearDownClass()
        return test_result

 #start of IDistributable interface--------------------------------------

    _work_count = None
    @property
    def work_count(self):
        if self._work_count is None:
            self._work_count = 0
            for test in self.deep_iter(self._test_or_suite):
                self._work_count += 1
        return self._work_count


    def work_sequence(self):
        for i, test in enumerate(self.deep_iter(self._test_or_suite)):
            yield lambda test=test : self.do_work(test)  # the 'test=test' is need to get around a strangeness in Python

    def reduce(self, result_sequence):
        fp = io.StringIO() if sys.version_info >= (3,0) else io.BytesIO()
        error_count = 0
        failure_count = 0
        test_result_list = []
        for test_result in result_sequence:
            test_result_list.append(test_result)
            for error in test_result.errors:
                error_count += 1
                fp.write("Error # {0} - {1}\n{2}\n============\n".format(error_count, error[0],error[1]))
            for failure in test_result.failures:
                failure_count += 1
                fp.write("Error # {0} - {1}\n{2}\n============\n".format(failure_count, failure[0],failure[1]))

        fp.write("\nFinal Error # {0}. Final Failure # {1}\n".format(error_count,failure_count))
        s = fp.getvalue()
        fp.close()

        print(s)
        return s

    @property
    def tempdirectory(self):
        return self.__tempdirectory

    #optional override -- the str name of the instance is used by the cluster as the job name
    def __str__(self):
        return "{0}(...)".format(self.__class__.__name__)
 #end of IDistributable interface---------------------------------------

