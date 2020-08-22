from __future__ import absolute_import
from pysnptools.util.mapreduce1.runner import Runner,_JustCheckExists, _run_all_in_memory, _shape_to_desired_workcount, _work_sequence_for_one_index
import os
import logging
from six.moves import range
try:
    import dill as pickle
except:
    logging.warning("Can't import dill, so won't be able to clusterize lambda expressions. If you try, you'll get this error 'Can't pickle <type 'function'>: attribute lookup __builtin__.function failed'")
    import six.moves.cPickle as pickle
import subprocess, sys, os.path
import threading
import pysnptools.util as util
from six.moves.queue import PriorityQueue

class LocalMultiThread(Runner):
    '''
    A :class:`.Runner` that runs a :func:`.map_reduce` as multiple threads on a single machine.

    Note that Python has problems running some programs efficiently on multiple threads. (Search 'python global interpreter lock' for details.)
    If this :class:`.Runner` doesn't work, consider :class:`LocalMultiProc` instead.

    **Constructor:**
        :Parameters: * **taskcount** (*number*) -- The number of threads to run on.
        :Parameters: * **mkl_num_threads** (*number*) -- (default None) Limit on the number threads used by the NumPy MKL library.
        :Parameters: * **just_one_process** (*bool*) -- (default False) Divide the work for multiple threads, but sequentially on one thread. Can be useful for debugging.
        
        :Example:

        >>> from pysnptools.util.mapreduce1 import map_reduce
        >>> from pysnptools.util.mapreduce1.runner import LocalMultiThread
        >>> from six.moves import range #Python 2 & 3 compatibility
        >>> def holder1(n,runner):
        ...     def mapper1(x):
        ...         return x*x
        ...     def reducer1(sequence):
        ...        return sum(sequence)
        ...     return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        >>> holder1(100,LocalMultiThread(4))
        328350

    '''

    def __init__(self, taskcount, mkl_num_threads = None, just_one_process = False,):
        if not 0 < taskcount: raise Exception("Expect taskcount to be positive")

        self.taskcount = taskcount
        self.just_one_process = just_one_process
        if mkl_num_threads != None:
            os.environ['MKL_NUM_THREADS'] = str(mkl_num_threads)

    def _result_sequence(self,thread_list,priority_queue,shaped_distributable):
        for thread in thread_list:
            if not self.just_one_process:
                thread.join()
            result_sequence = priority_queue.get()[1]
            for result in result_sequence:
                yield result

    def run(self, distributable):
        _JustCheckExists().input(distributable)

        priority_queue = PriorityQueue()
        thread_list = []
        shaped_distributable = _shape_to_desired_workcount(distributable, self.taskcount)
        for taskindex in range(self.taskcount):
            def _target(taskindex=taskindex):
                result_list = []
                for work in _work_sequence_for_one_index(shaped_distributable, self.taskcount, taskindex):
                    result_list.append(_run_all_in_memory(work))
                priority_queue.put((taskindex,result_list))
            if not self.just_one_process:
                thread = threading.Thread(target=_target,name=str(taskindex))
                thread_list.append(thread)
                thread.start()
            else:
                thread_list.append(None)
                _target()
        
        result_sequence = self._result_sequence(thread_list, priority_queue,shaped_distributable)
        result = shaped_distributable.reduce(result_sequence)

        _JustCheckExists().output(distributable)
        return result

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
