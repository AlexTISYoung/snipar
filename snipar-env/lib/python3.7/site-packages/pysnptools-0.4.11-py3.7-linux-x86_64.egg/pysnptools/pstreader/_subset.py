import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.pstreader import PstReader
from pysnptools.pstreader import PstData

class _PstSubset(PstReader):

    def __init__(self, internal, row_indexer, col_indexer):
        '''
        an indexer can be:
             an integer i (same as [i])
             a slice
             a list of integers (including negatives)
             a list of Booleans
        '''
        super(_PstSubset, self).__init__()
        self._ran_once = False

        self._internal = internal
        self._row_indexer = PstReader._make_sparray_or_slice(row_indexer)
        self._col_indexer = PstReader._make_sparray_or_slice(col_indexer)


    def __repr__(self):
        s = "{0}[{1},{2}]".format(self._internal,_PstSubset.static_nice_string(self,self._row_indexer),_PstSubset.static_nice_string(self,self._col_indexer))
        return s

    def copyinputs(self, copier):
        self._internal.copyinputs(copier)

    @property
    def row(self):
        self._run_once()
        return self._row

    @property
    def col(self):
        self._run_once()
        return self._col

    @property
    def row_property(self):
        self._run_once()
        return self._row_property

    @property
    def col_property(self):
        self._run_once()
        return self._col_property

    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = True
    def _read(self, row_indexer, col_indexer, order, dtype, force_python_only, view_ok):
        self._run_once()

        if hasattr(self._internal,'_read_accepts_slices'):
            assert self._internal._read_accepts_slices, "If an object has the _read_accepts_slices attribute, it must have value 'True'"
            composed_row_index_or_none = _PstSubset.compose_indexer_with_indexer(self._internal.row_count, self._row_indexer, self.row_count, row_indexer)
            composed_col_index_or_none = _PstSubset.compose_indexer_with_indexer(self._internal.col_count, self._col_indexer, self.col_count, col_indexer)
            val = self._internal._read(composed_row_index_or_none, composed_col_index_or_none, order, dtype, force_python_only, view_ok)
            return val
        else:
            row_index_or_none = PstReader._make_sparray_from_sparray_or_slice(self.row_count, row_indexer)
            composed_row_index_or_none = _PstSubset.compose_indexer_with_index_or_none(self._internal.row_count, self._row_indexer, self.row_count, row_index_or_none)
            col_index_or_none = PstReader._make_sparray_from_sparray_or_slice(self.col_count, col_indexer)
            composed_col_index_or_none = _PstSubset.compose_indexer_with_index_or_none(self._internal.col_count, self._col_indexer, self.col_count, col_index_or_none)
            val = self._internal._read(composed_row_index_or_none, composed_col_index_or_none, order, dtype, force_python_only, view_ok)
            return val

    def _run_once(self):
        if self._ran_once:
            return

        self._ran_once = True
        self._row = self._internal.row[self._row_indexer] #!!! would be nice if these four calls were replaced by calls to, e.g. "def _col_index(self,col_indexer)", etc that could be overridden (by for example Pairs)
        self._col = self._internal.col[self._col_indexer]
        if self._row.dtype == self._col.dtype and np.array_equal(self._row,self._col): #When an object is square, keep the row and col the same object.
            self._col = self._row
        self._row_property = self._internal.row_property[self._row_indexer]
        self._col_property = self._internal.col_property[self._col_indexer]

    _slice_format = {(False,False,False):":",
                     (False,False,True):"::{2}",
                     (False,True,False):":{1}",
                     (False,True,True):":{1}:{2}",
                     (True,False,False):"{0}:",
                     (True,False,True):"{0}::{2}",
                     (True,True,False):"{0}:{1}",
                     (True,True,True):"{0}:{1}:{2}"}

    @staticmethod
    def static_nice_string(self, some_slice):
        if isinstance(some_slice,slice):
            return _PstSubset._slice_format[(some_slice.start is not None, some_slice.stop is not None, some_slice.step is not None)].format(some_slice.start, some_slice.stop, some_slice.step)
        elif len(some_slice) == 1:
            return str(some_slice[0])
        elif len(some_slice) < 10:
            return "[{0}]".format(",".join([str(i) for i in some_slice]))
        else:
            return "[{0},...]".format(",".join([str(i) for i in some_slice[:10]]))


    #!!commented out because doesn't guarantee that the shortcut will return with the dtype and order requested.
    #                  Also, didn't handle stacked do-nothing subsets
    #def read(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False):
    #    if view_ok and hasattr(self._internal,"val") and _PstSubset._is_all_slice(self._row_indexer) and _PstSubset._is_all_slice(self._col_indexer):
    #        return self._internal
    #    else:
    #        return PstReader.read(self, order, dtype, force_python_only, view_ok)


    @staticmethod
    def compose_indexer_with_index_or_none(countA, indexerA, countB, index_or_noneB):
        if _PstSubset._is_all_slice(indexerA):
            return index_or_noneB

        indexA = PstReader._make_sparray_from_sparray_or_slice(countA, indexerA)

        if _PstSubset._is_all_slice(index_or_noneB):
            return indexA

        indexAB = indexA[index_or_noneB]

        return indexAB


    @staticmethod
    def compose_indexer_with_indexer(countA, indexerA, countB, indexerB):
        if _PstSubset._is_all_slice(indexerA):
            return indexerB

        if _PstSubset._is_all_slice(indexerB):
            return indexerA

        indexA = PstReader._make_sparray_from_sparray_or_slice(countA, indexerA)
        indexB = PstReader._make_sparray_from_sparray_or_slice(countB, indexerB)

        indexAB = indexA[indexB]

        return indexAB
