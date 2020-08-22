import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.pstreader import PstReader
from pysnptools.pstreader.pstdata import PstData
import pysnptools.util as pstutil
import warnings

class PstNpz(PstReader):
    '''
    A :class:`.PstReader` for reading \*.pst.npz files from disk.

    See :class:`.PstReader` for general examples of using PstReaders.

    The general NPZ format is described in http://docs.scipy.org/doc/numpy/reference/generated/numpy.savez.html. The PstNpz format stores
    val, row, col, row_property, and col_property information in NPZ format.
   
    **Constructor:**
        :Parameters: * **filename** (*string*) -- The PstNpz file to read.

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.pstreader import PstNpz
        >>> data_on_disk = PstNpz('../../tests/datasets/little.pst.npz')
        >>> print(data_on_disk.row_count)
        300

    **Methods beyond** :class:`.NpzReader`

    '''


    def __init__(self, filename):
        '''
        filename    : string of the name of the npz file.
        '''
        super(PstNpz, self).__init__() #We don't send "file name" up because we know about super doesn't want it.
        self._ran_once = False

        self._filename = filename

    def __repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self._filename)

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

    def _run_once(self):
        if (self._ran_once):
            return
        self._ran_once = True

        with np.load(self._filename,allow_pickle=True) as data: #!! similar code in epistasis
            if len(data.keys()) == 2 and 'arr_0' in data.keys(): #for backwards compatibility
                self._row = data['arr_0']
                self._col = self._row
                self._row_property = np.empty((len(self._row),0))
                self._col_property = np.empty((len(self._col),0))
            else:
                self._row = data['row']
                self._col = data['col']
                if self._row.dtype == self._col.dtype and np.array_equal(self._row, self._col): #If it's square, mark it so by making the col and row the same object
                    self._col = self._row
                self._row_property = data['row_property']
                self._col_property = data['col_property']

        return self

    def copyinputs(self, copier):
        # doesn't need to self._run_once()
        copier.input(self._filename)

    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = True
    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        # 'view_ok' doesn't mean anything here because we are always ready fresh from disk.
        #!! could use mmap so only rows of interest are loaded.
        self._run_once()

        #np.load does the right thing and doesn't load 'val' into memory until accessed here.
        with np.load(self._filename,allow_pickle=True) as data: #!! similar code in epistasis
            if len(data.keys()) == 2 and  'arr_1' in data.keys(): #for backwards compatibility
               val = data['arr_1']
            else:
               val = data['val']

        val, _ = self._apply_sparray_or_slice_to_val(val, row_index_or_none, col_index_or_none, order, dtype, force_python_only)
        return val

    @staticmethod
    def write(filename, pstdata):
        """Writes a :class:`PstData` to PstNpz format and returns the :class:`.PstNpz`.

        :param filename: the name of the file to create
        :type filename: string
        :param pstdata: The in-memory data that should be written to disk.
        :type pstdata: :class:`PstData`
        :rtype: :class:`.PstNpz`

        >>> from pysnptools.pstreader import PstData, PstNpz
        >>> import pysnptools.util as pstutil
        >>> data1 = PstData(row=['a','b','c'],col=['y','z'],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        >>> pstutil.create_directory_if_necessary("tempdir/tiny.pst.npz")
        >>> PstNpz.write("tempdir/tiny.pst.npz",data1)          # Write data in PstNz format
        PstNpz('tempdir/tiny.pst.npz')
        """
        np.savez(filename, row=pstdata.row, col=pstdata.col, row_property=pstdata.row_property, col_property=pstdata.col_property,val=pstdata.val)
        logging.debug("Done writing " + filename)
        return PstNpz(filename)


if __name__ == "__main__":
    import doctest
    logging.basicConfig(level=logging.INFO)
    
    doctest.testmod()

    #from pysnptools.snpreader.dat import Dat
    #snpreader = Dat(r'../tests/datasets/all_chr.maf0.001.N300.dat')
    #snp_matrix = snpreader.read()
    #print(len(snp_matrix['sid']))
    #snp_matrix = snpreader[:,:].read()
    #print(len(snp_matrix['sid']))
    #sid_index_list = snpreader.sid_to_index(['23_9','23_2'])
    #snp_matrix = snpreader[:,sid_index_list].read()
    #print(",".join(snp_matrix['sid']))
    #snp_matrix = snpreader[:,0:10].read()
    #print(",".join(snp_matrix['sid']))

    #print(snpreader.iid_count)
    #print(snpreader.sid_count)
    #print(len(snpreader.pos))

    #snpreader2 = snpreader[::-1,4]
    #print(snpreader.iid_count)
    #print(snpreader2.sid_count)
    #print(len(snpreader2.pos))

    #snp_matrix = snpreader2.read()
    #print(len(snp_matrix['iid']))
    #print(len(snp_matrix['sid']))

    #snp_matrix = snpreader2[5,:].read()
    #print(len(snp_matrix['iid']))
    #print(len(snp_matrix['sid']))

    #iid_index_list = snpreader2.iid_to_index(snpreader2.iid[::2])
    #snp_matrix = snpreader2[iid_index_list,::3].read()
    #print(len(snp_matrix['iid']))
    #print(len(snp_matrix['sid']))

    #snp_matrix = snpreader[[4,5],:].read()
    #print(len(snp_matrix['iid']))
    #print(len(snp_matrix['sid']))

    #print(snpreader2)
    #print(snpreader[::-1,4])
    #print(snpreader2[iid_index_list,::3])
    #print(snpreader[:,sid_index_list])
    #print(snpreader2[5,:])
    #print(snpreader[[4,5],:])
