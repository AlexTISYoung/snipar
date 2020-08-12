import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.snpreader import SnpReader
from pysnptools.standardizer import Unit
from pysnptools.standardizer import Identity
from pysnptools.pstreader import PstData
import warnings
import time

class SnpData(PstData,SnpReader):
    """A :class:`.SnpReader` for holding SNP values (or similar values) in-memory, along with related iid and sid information.
    It is usually created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`. It can also be constructed directly.

    See :class:`.SnpReader` for details and examples.

    **Constructor:**
        :Parameters: * **iid** (an array of string pair) -- The :attr:`.SnpReader.iid` information.
                     * **sid** (an array of strings) -- The :attr:`.SnpReader.sid` information.
                     * **val** (a 2-D array of floats) -- The SNP values
                     * **pos** (optional, an array of strings) -- The :attr:`.SnpReader.pos` information
                     * **name** (optional, string) -- Information to be display about the origin of this data
                     * **copyinputs_function** (optional, function) -- *Used internally by optional clustering code*

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import SnpData
        >>> snpdata = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]])
        >>> print((snpdata.val[0,1], snpdata.iid_count, snpdata.sid_count))
        (2.0, 2, 3)

    **Equality:**

        Two SnpData objects are equal if their four arrays (:attr:`.SnpData.val`, :attr:`SnpReader.iid`, :attr:`.SnpReader.sid`, and :attr:`.SnpReader.pos`)
        are 'array_equal'. (Their 'name' does not need to be the same).
        If either :attr:`.SnpData.val` contains NaN, the objects will not be equal. However, :meth:`.SnpData.allclose` can be used to treat NaN as
        regular values.

        :Example:

        >>> import numpy as np
        >>> from pysnptools.snpreader import SnpData
        >>> snpdata1 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> snpdata2 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> print(snpdata1 == snpdata2) #True, because all the arrays have the same values.
        True
        >>> print(snpdata1.val is snpdata2.val) #False, because the two arrays have different memory.
        False
        >>> snpdata3 = SnpData(iid=[['a','0'],['b','0']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> snpdata4 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> print(snpdata3 == snpdata4) #False, because the iids are different.
        False
        >>> snpdata5 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,np.nan]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> snpdata6 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,np.nan]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> print(snpdata5 == snpdata6) #False, because the val's contain NaN
        False
        >>> print(snpdata5.allclose(snpdata6)) #True, if we consider the NaN as regular values, all the arrays have the same values.
        True

    **Methods beyond** :class:`.SnpReader`
    """

    def __init__(self, iid, sid, val, pos=None, name=None, parent_string=None, copyinputs_function=None):

        #We don't have a 'super(SnpData, self).__init__()' here because SnpData takes full responsibility for initializing both its superclasses

        self._val = None

        if parent_string is not None:
            warnings.warn("'parent_string' is deprecated. Use 'name'", DeprecationWarning)
        self._row = PstData._fixup_input(iid,empty_creator=lambda ignore:np.empty([0,2],dtype='str'),dtype='str')
        self._col = PstData._fixup_input(sid,empty_creator=lambda ignore:np.empty([0],dtype='str'),dtype='str')
        self._row_property = PstData._fixup_input(None,count=len(self._row),empty_creator=lambda count:np.empty([count,0],dtype='str'),dtype='str')
        self._col_property = PstData._fixup_input(pos,count=len(self._col),empty_creator=lambda count:np.full([count, 3], np.nan))
        self._val = PstData._fixup_input_val(val,row_count=len(self._row),col_count=len(self._col),empty_creator=lambda row_count,col_count:np.empty([row_count,col_count],dtype=np.float64))
        self._assert_iid_sid_pos()
        self._name = name or parent_string or ""
        self._std_string_list = []

    val = property(PstData._get_val,PstData._set_val)
    """The 2D NumPy array of floats that represents the values of the SNPs.

    >>> from pysnptools.snpreader import Bed
    >>> snpdata = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)[:5,:].read() #read data for first 5 iids
    >>> print(snpdata.val[4,100]) #print one of the SNP values
    2.0
    """

    def allclose(self,value,equal_nan=True):
        '''
        :param value: Other object with which to compare.
        :type value: :class:`SnpData`
        :param equal_nan: (Default: True) Tells if NaN in :attr:`.SnpData.val` should be treated as regular values when testing equality.
        :type equal_nan: bool

        >>> import numpy as np
        >>> snpdata5 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,np.nan]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> snpdata6 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,np.nan]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> print(snpdata5.allclose(snpdata6)) #True, if we consider the NaN as regular values, all the arrays have the same values.
        True
        >>> print(snpdata5.allclose(snpdata6,equal_nan=False)) #False, if we consider the NaN as special values, all the arrays are not equal.
        False

        '''
        return PstData.allclose(self,value,equal_nan=equal_nan)

    def train_standardizer(self, apply_in_place, standardizer=Unit(), force_python_only=False):
        """
        .. deprecated:: 0.2.23
           Use :meth:`standardize` with return_trained=True instead.
        """
        warnings.warn("train_standardizer is deprecated. standardize(...,return_trained=True,...) instead", DeprecationWarning)
        assert apply_in_place, "code assumes apply_in_place"
        self._std_string_list.append(str(standardizer))
        _, trained_standardizer = standardizer.standardize(self, return_trained=True, force_python_only=force_python_only)
        return trained_standardizer

    #LATER should there be a single warning if Unit() finds and imputes NaNs?
    def standardize(self, standardizer=Unit(), block_size=None, return_trained=False, force_python_only=False):
        """Does in-place standardization of the in-memory
        SNP data. By default, it applies 'Unit' standardization, that is: the values for each SNP will have mean zero and standard deviation 1.0.
        NaN values are then filled with zero, the mean (consequently, if there are NaN values, the final standard deviation will not be zero.
        Note that, for efficiency, this method works in-place, actually changing values in the ndarray. Although it works in place, for convenience
        it also returns the SnpData.

        :param standardizer: optional -- Specify standardization to be applied. 
             Any :class:`.Standardizer` may be used. Some choices include :class:`.Unit` (default, makes values for each SNP have mean zero and
             standard deviation 1.0) and :class:`.Beta`.
        :type standardizer: :class:`.Standardizer`

        :param block_size: optional -- Deprecated.
        :type block_size: None

        :param return_trained: If true, returns a second value containing a constant :class:`.Standardizer` trained on this data.
        :type return_trained: bool

        :param force_python_only: optional -- If true, will use pure Python instead of faster C++ libraries.
        :type force_python_only: bool

        :rtype: :class:`.SnpData` (standardizes in place, but for convenience, returns 'self')

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify some data on disk in Bed format
        >>> snpdata1 = snp_on_disk.read() # read all SNP values into memory
        >>> print(snpdata1) # Prints the specification for this SnpData
        SnpData(Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False))
        >>> print(snpdata1.val[0,0])
        2.0
        >>> snpdata1.standardize() # standardize changes the values in snpdata1.val and changes the specification.
        SnpData(Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False),Unit())
        >>> print('{0:.6f}'.format(snpdata1.val[0,0]))
        0.229416
        >>> snpdata2 = snp_on_disk.read().standardize() # Read and standardize in one expression with only one ndarray allocated.
        >>> print('{0:.6f}'.format(snpdata2.val[0,0]))
        0.229416
        """
        self._std_string_list.append(str(standardizer))
        _, trained = standardizer.standardize(self, return_trained=True, force_python_only=force_python_only)
        if return_trained:
            return self, trained
        else:
            return self

    def _read_kernel(train, standardizer, block_size=None, order='A', dtype=np.float64, force_python_only=False, view_ok=False, return_trained=False):
        '''
        The method creates a kernel for the in-memory SNP data. It handles these cases
                * No standardization is needed & everything is in memory  OR uses the FROM-DISK method
        '''
        from pysnptools.pstreader import PstReader


        #Just do a 'python' dot, if no standardization is needed and everything is the right type
        if isinstance(standardizer,Identity) and train.val.dtype == dtype:
            ts = time.time()
            #is_worth_logging = train.val.shape[0] * train.val.shape[1] * test.val.shape[0] > 1e9
            #if is_worth_logging: logging.info("  _read_kernel about to multiply train{0} x test{1}".format(train.val.shape,test.val.shape))
            if order == 'F': #numpy's 'dot' always returns 'C' order
                K = (train.val.dot(train.val.T)).T
            else:
                K = train.val.dot(train.val.T)
            assert PstReader._array_properties_are_ok(K,order,dtype), "internal error: K is not of the expected order or dtype"
            #if is_worth_logging: logging.info("  _read_kernel took %.2f seconds" % (time.time()-ts))
            if return_trained:
                return K, standardizer
            else:
                return K
        else: #Do things the more general SnpReader way.
            return SnpReader._read_kernel(train, standardizer, block_size=block_size, order=order, dtype=dtype, force_python_only=force_python_only,view_ok=view_ok, return_trained=return_trained)

    def __repr__(self):
        if self._name == "":
            if len(self._std_string_list) > 0:
                s = "{0}({1})".format(self.__class__.__name__,",".join(self._std_string_list))
            else:
                s = "{0}()".format(self.__class__.__name__)
        else:
            if len(self._std_string_list) > 0:
                s = "{0}({1},{2})".format(self.__class__.__name__,self._name,",".join(self._std_string_list))
            else:
                s = "{0}({1})".format(self.__class__.__name__,self._name)
        return s


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
