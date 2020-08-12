'''
Runs a distributable job locally in one process. Returns the value of the job.

See SamplePi.py for examples.
'''

from __future__ import absolute_import
from __future__ import print_function
from pysnptools.util.mapreduce1.runner import Runner, _run_all_in_memory
import os, sys
import logging

class Local(Runner):
    '''
    A :class:`.Runner` that runs a :func:`.map_reduce` locally. To save memory, it will feed the results of the mapper to the reducer as those results are computed.

    **Constructor:**
        :Parameters: * **mkl_num_threads** (*number*) -- (default None) Limit on the number threads used by the NumPy MKL library.
        :Parameters: * **logging_handler** (*stream*) --  (default stdout) Where to output logging messages.
        
        :Example:

        >>> from six.moves import range #Python 2 & 3 compatibility
        >>> from pysnptools.util.mapreduce1 import map_reduce
        >>> from pysnptools.util.mapreduce1.runner import Local
        >>> def holder1(n,runner):
        ...     def mapper1(x):
        ...         return x*x
        ...     def reducer1(sequence):
        ...        return sum(sequence)
        ...     return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        >>> holder1(100,Local())
        328350

    '''
    def __init__(self, mkl_num_threads = None, logging_handler=logging.StreamHandler(sys.stdout)):
        logger = logging.getLogger()
        if not logger.handlers:
            logger.setLevel(logging.INFO)
        for h in list(logger.handlers):
            logger.removeHandler(h)
        logger.addHandler(logging_handler)
        if logger.level == logging.NOTSET:
            logger.setLevel(logging.INFO)
        
        if mkl_num_threads != None:
            os.environ['MKL_NUM_THREADS'] = str(mkl_num_threads)

    def run(self, distributable):
        _JustCheckExists().input(distributable)
        result = _run_all_in_memory(distributable)
        _JustCheckExists().output(distributable)
        return result

class _JustCheckExists(object): #Implements ICopier

    def __init__(self,doPrintOutputNames=False):
        self.doPrintOutputNames = doPrintOutputNames
    
    def input(self,item):
        if isinstance(item, str):
            if not os.path.exists(item): raise Exception("Missing input file '{0}'".format(item))
        elif hasattr(item,"copyinputs"):
            item.copyinputs(self)
        # else -- do nothing

    def output(self,item):
        if isinstance(item, str):
            if not os.path.exists(item): raise Exception("Missing output file '{0}'".format(item))
            if self.doPrintOutputNames:
                print(item)
        elif hasattr(item,"copyoutputs"):
            item.copyoutputs(self)
        # else -- do nothing

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
