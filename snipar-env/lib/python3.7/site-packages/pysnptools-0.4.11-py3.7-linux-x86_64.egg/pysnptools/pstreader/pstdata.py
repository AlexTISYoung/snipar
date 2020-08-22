from __future__ import print_function

import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.pstreader import PstReader


def _default_empty_creator(count):
    return np.empty([count or 0, 0],dtype='str')

def _default_empty_creator_val(row_count,col_count):
    return np.empty([row_count,col_count],dtype='str')

class PstData(PstReader):
    '''A :class:`.PstReader` for holding values in-memory, along with related row and col information.
    It is usually created by calling the :meth:`.PstReader.read` method on another :class:`.PstReader`, for example, :class:`.PstNpz`. It can also be constructed directly.

    See :class:`.PstReader` for details and examples.

    **Constructor:**
        :Parameters: * **row** (an array of anything) -- The :attr:`.PstReader.row` information
                     * **col** (an array of anything) -- The :attr:`.PstReader.col` information
                     * **val** (a 2-D array of floats) -- The values
                     * **row_property** (optional, an array of anything) -- Additional information associated with each row.
                     * **col_property** (optional, an array of anything) -- Additional information associated with each col.
                     * **name** (optional, string) -- Information to be display about the origin of this data
                     * **copyinputs_function** (optional, function) -- *Used internally by optional clustering code*

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.pstreader import PstData
        >>> pstdata = PstData(row=[['fam0','iid0'],['fam0','iid1']], col=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]])
        >>> print(pstdata.val[0,1], pstdata.row_count, pstdata.col_count)
        2.0 2 3

    **Equality:**

        Two PstData objects are equal if their five arrays (:attr:`.PstData.val`, :attr:`.PstReader.row`, :attr:`.PstReader.col`, :attr:`.PstReader.row_property`, and :attr:`.PstReader.col_property`) arrays are equal.
        (Their 'name' does not need to be the same).  
        If either :attr:`.PstData.val` contains NaN, the objects will not be equal. However, :meth:`.PstData.allclose` can be used to treat NaN as regular values.
        Any NaN's in the other four arrays are treated as regular values.

        :Example:

        >>> import numpy as np
        >>> from pysnptools.pstreader import PstData
        >>> pstdata1 = PstData(row=[['fam0','iid0'],['fam0','iid1']], col=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]])
        >>> pstdata2 = PstData(row=[['fam0','iid0'],['fam0','iid1']], col=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]])
        >>> print(pstdata1 == pstdata2) #True, because all the arrays have the same values.
        True
        >>> print(pstdata1.val is pstdata2.val) #False, because the two arrays have different memory.
        False
        >>> pstdata3 = PstData(row=['a','b'], col=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]])
        >>> pstdata4 = PstData(row=[['fam0','iid0'],['fam0','iid1']], col=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]])
        >>> print(pstdata3 == pstdata4) #False, because the rows have different ids.
        False
        >>> pstdata5 = PstData(row=[['fam0','iid0'],['fam0','iid1']], col=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,np.nan]])
        >>> pstdata6 = PstData(row=[['fam0','iid0'],['fam0','iid1']], col=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,np.nan]])
        >>> print(pstdata5 == pstdata6) #False, because the val's contain NaN
        False
        >>> print(pstdata5.allclose(pstdata6)) #True, if we consider the NaN as regular values, all the arrays have the same values.
        True

    **Methods beyond** :class:`.PstReader`


    '''
    def __init__(self, row, col, val, row_property=None, col_property=None, name=None, parent_string=None, copyinputs_function=None):
        super(PstData, self).__init__()

        self._val = None

        self._row = PstData._fixup_input(row)
        self._col = PstData._fixup_input(col)
        if self._row.dtype == self._col.dtype and np.array_equal(self._row, self._col): #If it's square, mark it so by making the col and row the same object
            self._col = self._row
        self._row_property = PstData._fixup_input(row_property,count=len(self._row))
        self._col_property = PstData._fixup_input(col_property,count=len(self._col))


        self._val = PstData._fixup_input_val(val,row_count=len(self._row),col_count=len(self._col))

        name = name or parent_string or ""
        if parent_string is not None:
            warnings.warn("'parent_string' is deprecated. Use 'name'", DeprecationWarning)
        self._name = name

    def __eq__(a,b):
        return a.allclose(b,equal_nan=False)

    def allclose(self,value,equal_nan=True):
        '''
        :param value: Other object with which to compare.
        :type value: :class:`PstData`
        :param equal_nan: (Default: True) Tells if NaN in :attr:`.PstData.val` should be treated as regular values when testing equality.
        :type equal_nan: bool

        >>> import numpy as np
        >>> pstdata5 = PstData(row=[['fam0','iid0'],['fam0','iid1']], col=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,np.nan]])
        >>> pstdata6 = PstData(row=[['fam0','iid0'],['fam0','iid1']], col=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,np.nan]])
        >>> print(pstdata5.allclose(pstdata6)) #True, if we consider the NaN as regular values, all the arrays have the same values.
        True
        >>> print(pstdata5.allclose(pstdata6,equal_nan=False)) #False, if we consider the NaN as special values, all the arrays are not equal.
        False

        '''
        try:
            return (PstData._allclose(self.row,value.row,equal_nan=True) and
                    PstData._allclose(self.col,value.col,equal_nan=True) and
                    PstData._allclose(self.row_property,value.row_property,equal_nan=True) and
                    PstData._allclose(self.col_property,value.col_property,equal_nan=True) and
                    np.allclose(self.val,value.val,equal_nan=equal_nan))
        except:
            return False

    @staticmethod
    def _allclose(a,b,equal_nan=True):
        if not equal_nan:
            return np.array_equal(a,b)
        if a.dtype == 'float' and b.dtype == 'float':
            return np.allclose(a,b,equal_nan=equal_nan)
        return np.array_equal(a,b)

    @staticmethod
    def _fixup_input(input,count=None, empty_creator=_default_empty_creator,dtype=None):
        if input is None or len(input)==0:
            input = empty_creator(count)
        elif not isinstance(input,np.ndarray):
            input = np.array(input,dtype=dtype)

        assert count is None or len(input) == count, "Expect length of {0} for input {1}".format(count,input)

        return input

    @staticmethod
    def _fixup_input_val(input,row_count,col_count,empty_creator=_default_empty_creator_val):
        if input is None:
            assert row_count == 0 or col_count == 0, "If val is None, either row_count or col_count must be 0"
            input = _default_empty_creator_val(row_count, col_count)
        elif not isinstance(input,np.ndarray or (input.dtype not in [np.float32,np.float64])):
            input = np.array(input,dtype=np.float64)

        assert len(input.shape)==2, "Expect val to be two dimensional."
        assert input.shape[0] == row_count, "Expect number of rows ({0}) in val to match the number of row names given ({1})".format(input.shape[0], row_count)
        assert input.shape[1] == col_count, "Expect number of columns ({0}) in val to match the number of column names given ({1})".format(input.shape[1], col_count)

        return input



    def __repr__(self):
        if self._name == "":
            return "{0}()".format(self.__class__.__name__)
        else:
            return "{0}({1})".format(self.__class__.__name__,self._name)

    def copyinputs(self, copier):
        pass

    @property
    def row(self):
        return self._row

    @property
    def col(self):
        return self._col

    @property
    def row_property(self):
        return self._row_property
    
    @property
    def col_property(self):
        return self._col_property

    def _get_val(self):
        return self._val

    def _set_val(self, new_value):
        self._val = new_value

    val = property(_get_val,_set_val)
    """The 2D NumPy array of floats that represents the values.

    >>> from pysnptools.pstreader import PstNpz
    >>> pstdata = PstNpz('../examples/toydata.pst.npz')[:5,:].read() #read data for first 5 rows
    >>> print(pstdata.val[4,100]) #print one of the values
    2.0
    """

    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = True
    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        val, shares_memory = self._apply_sparray_or_slice_to_val(self.val, row_index_or_none, col_index_or_none, order, dtype, force_python_only)
        if shares_memory and not view_ok:
            val = val.copy(order='K')
        return val


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
