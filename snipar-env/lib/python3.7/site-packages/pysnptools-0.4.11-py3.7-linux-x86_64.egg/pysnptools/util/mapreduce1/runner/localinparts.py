from pysnptools.util.mapreduce1.runner import Runner,_JustCheckExists, _run_one_task
import os, sys
import logging
import pysnptools.util as pstutil
try:
    import dill as pickle
except:
    logging.warning("Can't import dill, so won't be able to clusterize lambda expressions. If you try, you'll get this error 'Can't pickle <type 'function'>: attribute lookup __builtin__.function failed'")
    import cPickle as pickle

class LocalInParts(Runner):
    '''
    A :class:`.Runner` that runs one piece of a :func:`.map_reduce` job locally. Partial results are saved to disk.
    Clustering runners and :class:`LocalMultiProc` use this runner internally.

    **Constructor:**
        :Parameters: * **taskindex** (*number*) -- Which piece of work to run. When 0 to **taskcount**-1, does map work. When **taskcount**, does the reduce work. 
                     * **taskcount** (*number*) -- The number of pieces into which to divide the work.
                     * **mkl_num_threads** (*number*) -- (default None) Limit on the number threads used by the NumPy MKL library.
                     * **result_file** (*string*) -- (default None) Where to pickle the final results. If no file is given, the final results are returned, but not saved to a file.
                     * **result_dir** (*string*) -- (default None) The directory for any result_file. Defaults to the current working directory.
                     * **temp_dir** (*string*) -- (default None) The directory for partial results. Defaults to the **result_dir**/.working_directory.{map_reduce's Name}.
                     * **logging_handler** (*stream*) --  (default stdout) Where to output logging messages.
        
        :Example:

        >>> from pysnptools.util.mapreduce1 import map_reduce
        >>> from pysnptools.util.mapreduce1.runner import LocalInParts
        >>> from six.moves import range #Python 2 & 3 compatibility
        >>> def holder1(n,runner):
        ...     def mapper1(x):
        ...         return x*x
        ...     def reducer1(sequence):
        ...        return sum(sequence)
        ...     return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        >>> holder1(100,LocalInParts(0,4)) #Run part 0 of 4 and save partial results to disk as '0.4.p'.
        >>> holder1(100,LocalInParts(1,4)) #Run part 1 of 4 and save partial results to disk as '1.4.p'.
        >>> holder1(100,LocalInParts(2,4)) #Run part 2 of 4 and save partial results to disk as '2.4.p'.
        >>> holder1(100,LocalInParts(3,4)) #Run part 3 of 4 and save partial results to disk as '3.4.p'.
        >>> holder1(100,LocalInParts(4,4)) #Read the all partial results and then apply the reduce function and return the result.
        328350

    '''

    def __init__(self, taskindex, taskcount, mkl_num_threads=None, result_file=None, run_dir=".", temp_dir=None, logging_handler=logging.StreamHandler(sys.stdout)):
        logger = logging.getLogger()
        if not logger.handlers:
            logger.setLevel(logging.INFO)
        for h in list(logger.handlers):
            logger.removeHandler(h)
        if logger.level == logging.NOTSET:
            logger.setLevel(logging.INFO)
        logger.addHandler(logging_handler)

        self.temp_dir = temp_dir
        self.run_dir = run_dir
        self.result_file = os.path.join(run_dir,result_file) if result_file else None
        self.taskindex = taskindex
        self.taskcount = taskcount
        if mkl_num_threads != None:
            os.environ['MKL_NUM_THREADS'] = str(mkl_num_threads)


    def run(self, distributable):
        if self.temp_dir is not None:
            tempdir = self.temp_dir
        else:
            tempdir = os.path.join(self.run_dir,distributable.tempdirectory)
        tempdir = os.path.realpath(tempdir)
        if self.taskindex != self.taskcount:
            _JustCheckExists().input(distributable)
            return _run_one_task(distributable, self.taskindex, self.taskcount, tempdir)
        else:
            result = _run_one_task(distributable, self.taskindex, self.taskcount, tempdir)
            if self.result_file is not None:
                pstutil.create_directory_if_necessary(self.result_file)
                with open(self.result_file, mode='wb') as f:
                    pickle.dump(result, f, pickle.HIGHEST_PROTOCOL)

            return result

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
