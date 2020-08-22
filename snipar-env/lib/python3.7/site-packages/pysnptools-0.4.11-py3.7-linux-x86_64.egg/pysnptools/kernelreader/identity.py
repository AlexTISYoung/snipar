import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.kernelreader import KernelReader
from pysnptools.pstreader import PstData
from pysnptools.pstreader import PstReader

class Identity(KernelReader):
    '''
    A :class:`.KernelReader` for that represents an identity matrix. No memory for the values is allocated until :meth:`Identity.read` is called.

    See :class:`.KernelReader` for general examples of using KernelReaders.

    **Constructor:**
        :Parameters: * **iid** (an array of strings) -- The :attr:`KernelReader.iid` information

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.kernelreader import Identity
        >>> identity = Identity(iid=[['fam0','iid0'],['fam0','iid1']])
        >>> print(identity.iid_count)
        2
        >>> print(identity.read().val) # '...' is possible space character
        [[...1.  0.]
         [...0.  1.]]

        >>> identity = Identity(iid=[['fam0','iid0'],['fam0','iid1'],['fam0','iid2']],test=[['fam0','iid1'],['fam0','iid3']])
        >>> print((identity.iid0_count, identity.iid1_count))
        (3, 2)
        >>> print(identity.read().val) # '...' is possible space character
        [[...0.  0.]
         [...1.  0.]
         [...0.  0.]]

    '''
    def __init__(self, iid, test=None): #!!! add docs and test for test
        super(Identity, self).__init__()

        if test is None:
            test = iid

        if len(iid)>0:
            self._row0 = iid
        else:
            self._row0 = self._empty

        if len(test)>0:
            self._row1 = test
        else:
            self._row1 = self._empty

    _empty = np.empty([0,2],dtype='str')

    def __repr__(self): 
        return "{0}({1}x{2})".format(self.__class__.__name__, self.row_count, self.col_count)


    @property
    def row(self):
        return self._row0

    @property
    def col(self):
        return self._row1

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        if row_index_or_none is None and col_index_or_none is None and self._row0 is self._row1: #read all of a square ID
            return np.identity(self.row_count,dtype=dtype)
        else: #Non-square
            #!!! This is also less efficient than it could be because it create a big identity matrix and then slices it.

            #In about O(col_count + row_count) fill in zeros
            big = np.zeros([self.row_count,self.col_count],dtype=dtype)
            common = set([PstReader._makekey(x) for x in self.row]) & set([PstReader._makekey(x) for x in self.col])
            big[self.row_to_index(common),self.col_to_index(common)] = 1.0
            val, shares_memory = self._apply_sparray_or_slice_to_val(big, row_index_or_none, col_index_or_none, order, dtype, force_python_only)
            return val

    def __getitem__(self, iid_indexer_and_snp_indexer):
        if isinstance(iid_indexer_and_snp_indexer,tuple):
            row_index_or_none, col_index_or_none = iid_indexer_and_snp_indexer
        else:
            row_index_or_none = iid_indexer_and_snp_indexer
            col_index_or_none = row_index_or_none

        #If iid0 will be iid1 then make 'is' equal
        iid=self.iid0[row_index_or_none]
        test=self.iid1[col_index_or_none]
        if len(iid) == len(test) and np.array_equal(iid,test):
            result = Identity(iid)
        else: 
            result = Identity(iid=iid,test=test)
        return result


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
