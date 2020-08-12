import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.pstreader import PstReader
from pysnptools.pstreader.pstdata import PstData
import pysnptools.util as pstutil
import warnings

class _OneShot(PstReader):
    '''
    A utility class that makes it easier to define readers that read all the data and into memory before doing any subsetting.
    '''

    def _read_pstdata():
        raise NotImplementedError("{0} needs to define its own _read_pstdata".format(self.__class__.__name__))

    def __init__(self):
        super(_OneShot, self).__init__()

        self._ran_once = False

    def __repr__(self): 
        if hasattr(self,"filename"):
            return "{0}('{1}')".format(self.__class__.__name__,self.filename)
        else:
            return "{0}".format(self.__class__.__name__)

    def copyinputs(self, copier):
        if hasattr(self,"filename"):
            # doesn't need to self._run_once()
            copier.input(self.filename)

    def _run_once(self):
        if (self._ran_once):
            return
        self._ran_once = True
        self._data = self._read_pstdata()

    @property
    def row(self):
        self._run_once()
        return self._data._row

    @property
    def col(self):
        self._run_once()
        return self._data._col

    @property
    def row_property(self):
        self._run_once()
        return self._data._row_property

    @property
    def col_property(self):
        self._run_once()
        return self._data._col_property

    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = True
    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        self._run_once()
        val = self._data._read(row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok)
        return val

