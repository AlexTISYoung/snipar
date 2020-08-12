try:
    import h5py
except:
    pass

import logging
import numpy as np
from pysnptools.pstreader import PstReader
from pysnptools.pstreader.pstdata import PstData
import warnings
from six.moves import range

class PstHdf5(PstReader):
    '''
    A :class:`.PstReader` for reading \*.pst.hdf5 files from disk.

    See :class:`.PstReader` for general examples of using PstReaders.

    The general HDF5 format is described in http://www.hdfgroup.org/HDF5/. The PstHdf5 format stores
    val, row, col, row_property, and col_property information in Hdf5 format.
   
    **Constructor:**
        :Parameters: * **filename** (*string*) -- The PstHdf5 file to read.

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.pstreader import PstHdf5
        >>> on_disk = PstHdf5('../examples/toydata.iidmajor.snp.hdf5') # PstHdf5 can load .pst.hdf5, .snp.hdf5, and kernel.hdf5
        >>> print(on_disk.row_count)
        500

    **Methods beyond** :class:`.PstReader`
    '''

    def __init__(self, filename):
        super(PstHdf5, self).__init__() #We know PstReader doesn't want the file name

        self._block_size = 5000

        self._ran_once = False
        self._h5 = None

        self.filename=filename


    def __repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self.filename) #!!LATER print non-default values, too

    def copyinputs(self, copier):
        copier.input(self.filename)

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

    def _find_vocab(self):
        vocab_list = [['row','col','val','row_property','col_property'],['iid','sid','val',None,'pos'],['iid','rs','snps',None,'pos']]

        for vocab in vocab_list:
            if all((key is None or key in self._h5) for key in vocab):
                return vocab
        raise Exception("Don't know how to read HDF5 with these keys: " + ",".join(self._h5.iterkeys()))


    def _run_once(self):
        if self._ran_once:
            return
        try:
            self._h5 = h5py.File(self.filename, "r")
        except IOError as e:
            raise IOError("Missing or unopenable file '{0}' -- Native error message: {1}".format(self.filename,e))

        row_key,col_key,val_key,row_property_key,col_property_key = self._find_vocab()

        self._row = PstData._fixup_input(self._h5[row_key])
        self._col = PstData._fixup_input(self._h5[col_key])
        if self._row.dtype == self._col.dtype and np.array_equal(self._row,self._col):  #If it's square, mark it so by making the col and row the same object
            self._row = self._col
        self._row_property = PstData._fixup_input(self._h5[row_property_key] if row_property_key else None,count=len(self._row))  #Extra "if ... else" for backwards compatibility.
        self._col_property = PstData._fixup_input(self._h5[col_property_key],count=len(self._col))
        self.val_in_file = self._h5[val_key]

        self.is_col_major = None
        if "col-major" in self.val_in_file.attrs:
            self.is_col_major = self.val_in_file.attrs["col-major"]
        elif "SNP-major" in self.val_in_file.attrs:
            self.is_col_major = self.val_in_file.attrs["SNP-major"]
        assert self.is_col_major is not None, "In Hdf5 the 'val' matrix must have a Boolean 'col-major' (or 'SNP-major') attribute"

        S_original = len(self._col)
        N_original = len(self._row)
        if self.is_col_major:
            if not self.val_in_file.shape == (S_original, N_original) : raise Exception("In Hdf5, the val matrix dimensions don't match those of 'row' and 'col'")
        else:
            if not self.val_in_file.shape == (N_original, S_original) : raise Exception("In Hdf5, the val matrix dimensions don't match those of 'row' and 'col'")

        self._ran_once = True


    @staticmethod
    def _is_sorted_without_repeats(list):
        if len(list) < 2:
            return True
        for i in range(1,len(list)):
            if not list[i-1] < list[i]:
                return False
        return True
    

    def __del__(self):
        if self._h5 != None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
            self._h5.close()

    def _read_direct(self, val, val_order, selection=np.s_[:,:]):
        if self.is_col_major:
            selection = tuple(reversed(selection))

        if val_order == "F":
            self.val_in_file.read_direct(val.T,selection)
        else:
            assert val_order == "C", "real assert"
            self.val_in_file.read_direct(val,selection)

    def _create_block(self, block_size, order, dtype):
        matches_order = self.is_col_major == (order =="F")
        opposite_order = "C" if order == "F" else "F"
        if matches_order:
            return np.empty([len(self._row),block_size], dtype=dtype, order=order), order
        else:
            return np.empty([len(self._row),block_size], dtype=dtype, order=opposite_order), opposite_order

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        self._run_once()

        assert order in ['F','C','A'], "Expect order to be 'F', 'C' or 'A'"

        if order == "A":
            order = "F" if self.is_col_major else "C"

        opposite_order = "C" if order == "F" else "F"


        if row_index_or_none is not None:
            row_index_count = len(row_index_or_none)
            row_index_list = row_index_or_none
            row_is_sorted = PstHdf5._is_sorted_without_repeats(row_index_list)
            if isinstance(row_index_list,np.ndarray):
                row_index_list = row_index_list.tolist()
        else:
            row_index_count = self.row_count
            row_index_list = range(self.row_count)
            row_is_sorted = True

        if col_index_or_none is not None:
            col_index_count = len(col_index_or_none)
            col_index_list = col_index_or_none
            if isinstance(col_index_list,np.ndarray):
                col_index_list = col_index_list.tolist()
        else:
            col_index_count = self.col_count
            col_index_list = list(range(self.col_count))
        #Check if snps and iids indexes are in order and in range
        col_are_sorted = PstHdf5._is_sorted_without_repeats(col_index_list)

        val = np.empty([row_index_count, col_index_count], dtype=dtype, order=order)

        matches_order = self.is_col_major == (order=="F")
        is_simple = not force_python_only and row_is_sorted and col_are_sorted and matches_order #If 'is_simple' may be able to use a faster reader

        # case 0 -- zero elements in val
        if row_index_count == 0 or col_index_count == 0:
            pass

        # case 1 - all cols & all rows requested
        elif is_simple and col_index_count == self.col_count and row_index_count == self.row_count:
            self._read_direct(val, order)

        # case 2 - some cols and all rows
        elif is_simple and row_index_count == self.row_count:
            self._read_direct(val, order, np.s_[:,col_index_list])

        # case 3 all cols and some row
        elif is_simple and col_index_count == self.col_count:
            self._read_direct(val, order, np.s_[row_index_list,:])

        # case 4 some cols and some rows -- use blocks
        else:
            block_size = min(self._block_size, col_index_count)
            block, block_order = self._create_block(block_size, order, dtype)

            if not col_are_sorted:
                col_index_index_list = np.argsort(col_index_list).tolist()
                col_index_list_sorted = list(np.array(col_index_list)[col_index_index_list])
            else:
                col_index_index_list = np.arange(col_index_count)
                col_index_list_sorted = col_index_list

            for start in range(0, col_index_count, block_size):
                #print start
                stop = min(start+block_size,col_index_count)
                if stop-start < block_size:  #On the last loop, the buffer might be too big, so make it smaller
                    block, block_order = self._create_block(stop-start, order, dtype)
                col_index_list_forblock = col_index_list_sorted[start:stop]
                col_index_index_list_forblock = col_index_index_list[start:stop]
                self._read_direct(block, block_order, np.s_[:,col_index_list_forblock])
                val[:,col_index_index_list_forblock] = block[row_index_list,:]

        #!!LATER does this test work when the size is 1 x 1 and order if F? iid_index_or_none=[0], sid_index_or_none=[1000] (based on test_blocking_hdf5)
        has_right_order = (order=="C" and val.flags["C_CONTIGUOUS"]) or (order=="F" and val.flags["F_CONTIGUOUS"])
        assert val.shape == (row_index_count, col_index_count) and val.dtype == dtype and has_right_order
        return val




    @staticmethod
    def write(filename, pstdata, hdf5_dtype=None, col_major=True):
        """Writes a :class:`PstData` to PstHdf5 format and returns the :class:`.PstHdf5`.

        :param filename: the name of the file to create
        :type filename: string
        :param pstdata: The in-memory data that should be written to disk.
        :type pstdata: :class:`PstData`
        :param hdf5_dtype: None (use the .val's dtype) or a Hdf5 dtype, e.g. 'f8','f4',etc.
        :type hdf5_dtype: string
        :param col_major: Tells if vals should be stored on disk in col_major (default) or row_major format.
        :type col_major: bool
        :rtype: :class:`.PstHdf5`

        >>> import numpy as np
        >>> from pysnptools.pstreader import PstData, PstHdf5
        >>> import pysnptools.util as pstutil
        >>> data1 = PstData(row=['a','b','c'],col=['y','z'],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        >>> pstutil.create_directory_if_necessary("tempdir/tiny.pst.hdf5")
        >>> PstHdf5.write("tempdir/tiny.pst.hdf5",data1)          # Write data in PstHdf5 format
        PstHdf5('tempdir/tiny.pst.hdf5')
        """

        assert hdf5_dtype is None or (isinstance(hdf5_dtype, str) and len(hdf5_dtype) == 2 and  hdf5_dtype[0] == 'f'), "Expect hdf5_dtype to be None or to start with 'f', e.g. 'f4' for single, 'f8' for double"

        val = (pstdata.val.T) if col_major else pstdata.val

        def any_u_to_a(possible_unicode):
            #If it's any kind of string, encode it as ascii
            if possible_unicode.dtype.char == 'S' or possible_unicode.dtype.char == 'U': #not using np.issubdtype because of future warning
                return np.array(possible_unicode,dtype='S')
            else: #Otherwise, just leave it.
                return possible_unicode

        with h5py.File(filename, "w") as h5:
            h5.create_dataset('row', data=any_u_to_a(pstdata.row))
            h5.create_dataset('col', data=any_u_to_a(pstdata.col))
            h5.create_dataset('row_property', data=any_u_to_a(pstdata.row_property))
            h5.create_dataset('col_property', data=any_u_to_a(pstdata.col_property))
            h5.create_dataset('val', data=val,dtype=hdf5_dtype,shuffle=True)#compression="gzip", doesn't seem to work with Anaconda
            h5['val'].attrs["col-major"] = col_major

        return PstHdf5(filename)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
