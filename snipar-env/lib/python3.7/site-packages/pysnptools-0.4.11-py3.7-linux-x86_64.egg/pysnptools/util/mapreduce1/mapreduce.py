from __future__ import absolute_import
import logging
from pysnptools.util.mapreduce1.runner import _run_all_in_memory, Local
from contextlib import contextmanager
import threading
from six.moves import range

def _identity(x):
    return x

class _MapReduce(object): #implements IDistributable
    """
    class to run distributed map using the idistributable back-end
    """

    def __init__(self, input_seq, mapper, nested, reducer, input_files=None, output_files=None, name=None):

        self.input_seq = input_seq
        self.mapper = mapper
        self.nested = nested
        if (self.mapper is not _identity) and (self.nested is not None):
            raise Exception("'mapper' and 'nested' should not both be set")
        self.reducer = reducer
        self.name = name

        if input_files is None:
            self.input_files = []
        else:
            self.input_files = input_files

        if output_files is None:
            self.output_files = []
        else:
            self.output_files = output_files


    def __str__(self):
        return "map_reduce(name='{0}',...)".format(self.name)
    def __repr__(self):
        return str(self)


#start of IDistributable interface--------------------------------------
    @property
    def work_count(self):
        return len(self.input_seq)

    def work_sequence_range(self, start, stop):
        for i in range(start,stop):
            input_arg = self.input_seq[i]
            if self.nested is None:
                #logging.debug("\nrandom access executing %i" % i)
                with _MapReduce._dyn_vars(is_in_nested=False):
                    yield lambda i=i, input_arg=input_arg: self.dowork(i, input_arg)   # the 'i=i',etc is need to get around a strangeness in Python
            else:
                assert self.nested is not None, "real assert"
                with _MapReduce._dyn_vars(is_in_nested=True):
                    dist = self.nested(input_arg)
                    yield dist

    def work_sequence(self):
        for i, input_arg in enumerate(self.input_seq):
            if self.nested is None:
                #logging.debug("\nexecuting %i" % i)
                with _MapReduce._dyn_vars(is_in_nested=False):
                    yield lambda i=i, input_arg=input_arg: self.dowork(i, input_arg)  # the 'i=i',etc is need to get around a strangeness in Python
            else:
                assert self.nested is not None, "real assert"
                with _MapReduce._dyn_vars(is_in_nested=True):
                    dist = self.nested(input_arg)
                    yield dist


    def reduce(self, output_seq):
        '''
        '''

        return self.reducer(output_seq)


    #optional override -- the str name of the instance is used by the cluster as the job name
    def __str__(self):
        if self.name is None:
            return "map_reduce()"
        else:
            return self.name
 #end of IDistributable interface---------------------------------------

    def dowork(self, i, input_arg):
        #logging.info("{0}, {1}".format(len(train_snp_idx), len(test_snp_idx)))
        #logging.debug("\nexecuting {0}".format(input_arg))
        work = lambda : self.mapper(input_arg)
        result = _run_all_in_memory(work)
        return result

   
    # required by IDistributable
    @property
    def tempdirectory(self):
        return ".work_directory.{0}".format(self.name)
        

    def copyinputs(self, copier):
        for fn in self.input_files:
            copier.input(fn)

    def copyoutputs(self,copier):
        for fn in self.output_files:
            copier.output(fn)


    dyn = threading.local()

    # from short example in http://stackoverflow.com/questions/2001138/how-to-create-dynamical-scoped-variables-in-python999
    @staticmethod
    @contextmanager
    def _dyn_vars(**new):
        old = {}
        for name, value in new.items():
            old[name] = getattr(_MapReduce.dyn, name, None)
            setattr(_MapReduce.dyn, name, value)
        yield
        for name, value in old.items():
            setattr(_MapReduce.dyn, name, value)


def _is_in_nested():
    return hasattr(_MapReduce.dyn,"is_in_nested") and _MapReduce.dyn.is_in_nested
    
def map_reduce(input_seq, mapper=_identity, reducer=list, input_files=None, output_files=None, name=None, runner=None, nested=None):
    """
    Runs a function on sequence of inputs and runs a second function on the results. Can be nested and clusterized.

    :param input_seq: a sequence of inputs. The sequence must support the len function and be indexable. e.g. a list, xrange(100)
    :type input_seq: a sequence

    :param mapper: A function to apply to each set of inputs (optional). Defaults to the identity function.
    :type mapper: a function

    :param reducer: A function to turn the results from the mapper to a single value (optional). Defaults to creating a list of the results.
    :type reducer: a function that takes a sequence

    :param input_files: An optional list that tells what input files are needed. The list can contain the names of files (strings), None (ignored), or
        objects such as :class:`.SnpReader`'s that can self-report their input files.
    :type input_files: a list

    :param output_files: An optional list that tells what output files will be produced. The list can contain the names of files (strings), None (ignored), or
        objects such as :class:`.SnpReader`'s that can self-report their output files.
    :type output_files: a list

    :param name: A name to be displayed if this work is done on a cluster.
    :type name: a string

    :param runner: a :class:`.Runner`, optional: Tells how to run locally, multi-processor, or on a cluster.
        If not given, the function is run locally.
    :type runner: :class:`.Runner`

    :param nested: a mapper function that is itself a map_reduce. Some runners can efficiently clusterize such nested mappers. 
    :type nested: a function

    :rtype: The results from the reducer.


    :Example:

    Square the numbers 0 to 99 and report their sum, locally:

        >>> from pysnptools.util.mapreduce1 import map_reduce
        >>> from six.moves import range #Python 2 & 3 compatibility
        >>> map_reduce(range(100), 
        ...        mapper=lambda x: x*x,
        ...        reducer=sum)
        328350

    Compute it again, this time run on four processors:

        >>> from pysnptools.util.mapreduce1.runner import LocalMultiProc
        >>> from six.moves import range #Python 2 & 3 compatibility
        >>> map_reduce(range(100),
        ...        mapper=lambda x: x*x,
        ...        reducer=sum,
        ...        runner=LocalMultiProc(4))
        328350

    Compute it using named functions, again using four processors:

        >>> def holder1(n,runner):
        ...     def mapper1(x):
        ...         return x*x
        ...     def reducer1(sequence):
        ...        return sum(sequence)
        ...     return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        >>> holder1(100,LocalMultiProc(4))
        328350

    """

    dist = _MapReduce(input_seq, mapper=mapper, nested=nested, reducer=reducer, input_files=input_files, output_files=output_files,name=name)
    if runner is None and _is_in_nested():
        return dist

    if runner is None:
        runner = Local()

    result = runner.run(dist)
    return result
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    from pysnptools.util.mapreduce1 import map_reduce #Needed to work around thread local variable issue
    import doctest
    doctest.testmod()
