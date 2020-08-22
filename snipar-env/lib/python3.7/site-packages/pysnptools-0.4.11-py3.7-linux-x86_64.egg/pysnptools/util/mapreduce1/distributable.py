'''
A distributable job is a map/reduce (aka scatter/rather aka split/tabulate) job defined by a class that implements the
IDistributable interface (which Python doesn't check). FastLmmSet is an example of an distributable job.
Such distributable jobs can by run by IRunners.

See SamplePi.py for examples.

'''


import os
import logging
try:
    import dill as pickle
except:
    logging.warning("Can't import dill, so won't be able to clusterize lambda expressions. If you try, you'll get this error 'Can't pickle <type 'function'>: attribute lookup __builtin__.function failed'")
    import cPickle as pickle
import subprocess, sys, os.path
from pysnptools.util.mapreduce1.runner import *

class IDistributable(object):
    @property
    def work_count(self):
        raise NotImplementedError( "Should have implemented this" )
        return # Tells how many work items there are

    def work_sequence_range(self, start, stop):
        # Enumerates a sequence of work items starting at 'start' and ending right before 'end'
        # for example by 'yield'ing lambda (aka function pointers) that can be evaluated with no parameters
        #         the result of evaluating a work item can be anything, but it must be picklable
        # Can also yield IDistributable objects (thus the map-reduce can be nested).
        # If you like, you can define 'work_sequence_range' in terms of 'work_sequence'
        #         import itertools
        #         return islice(self.work_sequence(),start,stop)
        raise NotImplementedError( "Should have implemented this" )


    def work_sequence(self):
        # Enumerates a sequence of work items
        # for example by 'yield'ing lambda (aka function pointers) that can be evaluated with no parameters
        #         the result of evaluating a work item can be anything, but it must be picklable
        # Can also yield IDistributable objects (thus the map-reduce can be nested).
        # If you like, you can define 'work_sequence' in terms of 'work_sequence_range'
        #        return self.work_sequence_range(0,self.work_count)
        raise NotImplementedError( "Should have implemented this" )

    def reduce(self, result_sequence):
        # given a sequence of results and tabulates them
        # The tabulation may be written a file. Some value (perhaps 'None') is returned.
        raise NotImplementedError( "Should have implemented this" )

    @property
    def tempdirectory(self):
        raise NotImplementedError( "Should have implemented this" )
        return # the name of a directory where intermediate results could be stored. Does not need to exist; it will be created if necessary.

    #optional
    def __str__(self):
        raise NotImplementedError( "Should have implemented this" )
        return # a name for this job (is used by the cluster)

if __name__ == "__main__":
    '''
    From the command-line run a pickled distributable job.
    '''
    if len(sys.argv) != 3:
        print("Usage:python distributable.py <distributable.p> <runner>\n")
        sys.exit(1)

    distributablep_filename = sys.argv[1]
    if not os.path.exists(distributablep_filename): raise Exception(distributablep_filename + " does not exist")

    with open(distributablep_filename, mode='rb') as f:
        try:
            distributable = pickle.load(f)
        except AttributeError as e:
            raise AttributeError("An AttributeError when loading the pickle file is often caused by having the __main__ in the same file as a needed class. One possible fix is to add an import statement in the __main__ for the class. [Original message: '{0}'".format(e))
            

    runner_string = sys.argv[2]
    exec("runner = " + runner_string)
    runner.run(distributable)
