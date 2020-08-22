import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.kernelreader import KernelReader
from pysnptools.pstreader import PstData
from pysnptools.kernelstandardizer import Identity as KS_Identity
from pysnptools.kernelstandardizer import DiagKtoN

class KernelData(KernelReader,PstData):
    """A :class:`.KernelReader` for holding kernel values in-memory, along with related iid information.
    It is usually created by calling the :meth:`.SnpReader.read_kernel` method on a :class:`.SnpReader`, for example, :class:`.Bed`.
    It can also be created by calling the :meth:`.KernelReader.read` method on :class:`.SnpKernel`, for example, :class:`.KernelNpz`.
    It can also be constructed directly.

    See :class:`.KernelReader` for details and examples.

    **Constructor:**
        :Parameters: * **iid** (an array of string pairs) -- The :attr:`KernelReader.iid` information.
                     * **iid0** (an array of string pairs) -- The :attr:`KernelReader.iid0` information.
                     * **iid1** (an array of strings pairs) -- The :attr:`KernelReader.iid1` information.
                     * **val** (a 2-D array of floats) -- The SNP values
                     * **name** (optional, string) -- Information to be display about the origin of this data

        If *iid* is provided, don't provide *iid0* and *iid1*. Likewise, if *iid0* and *iid1* are provided, don't provide *iid*.

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.kernelreader import KernelData
        >>> kerneldata = KernelData(iid=[['fam0','iid0'],['fam0','iid1']], val=[[1.,.5],[.5,1.]])
        >>> print((kerneldata.val[0,1], kerneldata.iid_count))
        (0.5, 2)

    **Equality:**

        Two KernelData objects are equal if their three arrays (:attr:`.iid0`, :attr:`.iid1`, and :attr:`.KernelData.val`) are 'array_equal'.
        (Their 'string' does not need to be the same).
        If either :attr:`.KernelData.val` contains NaN, the objects will not be equal. However, :meth:`.KernelData.allclose` can be used to treat NaN as
        regular values.


        :Example:

        >>> import numpy as np
        >>> from pysnptools.kernelreader import KernelData
        >>> kerneldata1 = KernelData(iid=[['fam0','iid0'],['fam0','iid1']], val=[[1.,.5],[.5,1.]])
        >>> kerneldata2 = KernelData(iid=[['fam0','iid0'],['fam0','iid1']], val=[[1.,.5],[.5,1.]])
        >>> print(kerneldata1 == kerneldata2) #True, because all the arrays have the same values.
        True
        >>> print(kerneldata1.val is kerneldata2.val) #False, because the two arrays have different memory.
        False
        >>> kerneldata3 = KernelData(iid=[['a','0'],['b','0']], val=[[1.,.5],[.5,1.]])
        >>> kerneldata4 = KernelData(iid=[['fam0','iid0'],['fam0','iid1']], val=[[1.,.5],[.5,1.]])
        >>> print(kerneldata3 == kerneldata4) #False, because the iids are different.
        False
        >>> kerneldata5 = KernelData(iid=[['fam0','iid0'],['fam0','iid1']], val=[[1.,.5],[.5,np.nan]])
        >>> kerneldata6 = KernelData(iid=[['fam0','iid0'],['fam0','iid1']], val=[[1.,.5],[.5,np.nan]])
        >>> print(kerneldata5 == kerneldata6) #False, because the val's contain NaN
        False
        >>> print(kerneldata5.allclose(kerneldata6)) #True, if we consider the NaN as regular values, all the arrays have the same values.
        True


    **Methods beyond** :class:`.KernelReader`
    """

    def __init__(self, iid=None, iid0=None, iid1=None, val=None, name=None, parent_string=None): #!!!autodoc doesn't generate good doc for this constructor
        #We don't have a 'super(KernelData, self).__init__()' here because KernelData takes full responsiblity for initializing both its superclasses

        self._val = None

        #!!why does SnpData __init__ have a copy_inputs, but KernelData doesn't?
        assert (iid is None) != (iid0 is None and iid1 is None), "Either 'iid' or both 'iid0' 'iid1' must be provided."
        assert name is None or parent_string is None, "Can't set both 'name' and the deprecated 'parent_string'"
        if parent_string is not None:
            warnings.warn("'parent_string' is deprecated. Use 'name'", DeprecationWarning)

        if iid is not None:
            self._row = PstData._fixup_input(iid,empty_creator=lambda ignore:np.empty([0,2],dtype='str'),dtype='str')
            self._col = self._row
        else:
            self._row = PstData._fixup_input(iid0,empty_creator=lambda ignore:np.empty([0,2],dtype='str'),dtype='str')
            self._col = PstData._fixup_input(iid1,empty_creator=lambda ignore:np.empty([0,2],dtype='str'),dtype='str')
        self._row_property = PstData._fixup_input(None,count=len(self._row),empty_creator=lambda count:np.empty([count,0],dtype='str'),dtype='str')
        self._col_property = PstData._fixup_input(None,count=len(self._col),empty_creator=lambda count:np.empty([count,0],dtype='str'),dtype='str')
        self._val = PstData._fixup_input_val(val,row_count=len(self._row),col_count=len(self._col),empty_creator=lambda row_count,col_count:np.empty([row_count,col_count],dtype=np.float64))

        self._assert_iid0_iid1() 
        self._name = name or parent_string or ""
        self._std_string_list = []


    val = property(PstData._get_val,PstData._set_val)
    """The 2D NumPy array of floats that represents the values of the kernel.

    >>> from pysnptools.kernelreader import KernelData
    >>> kerneldata = KernelData(iid=[['fam0','iid0'],['fam0','iid1']], val=[[1.,.5],[.5,1.]])
    >>> print((kerneldata.val[0,1], kerneldata.iid_count))
    (0.5, 2)
    """

    def allclose(self, value,equal_nan=True):
        '''
        :param value: Other object with which to compare.
        :type value: :class:`KernelData`
        :param equal_nan: (Default: True) Tells if NaN in :attr:`.KernelData.val` should be treated as regular values when testing equality.
        :type equal_nan: bool

        >>> import numpy as np
        >>> kerneldata5 = KernelData(iid=[['fam0','iid0'],['fam0','iid1']], val=[[1.,.5],[.5,np.nan]])
        >>> kerneldata6 = KernelData(iid=[['fam0','iid0'],['fam0','iid1']], val=[[1.,.5],[.5,np.nan]])
        >>> print(kerneldata5.allclose(kerneldata6)) #True, if we consider the NaN as regular values, all the arrays have the same values.
        True
        >>> print(kerneldata5.allclose(kerneldata6,equal_nan=False)) #False, if we consider the NaN as special values, all the arrays are not equal.
        False

        '''
        return PstData.allclose(self,value,equal_nan=equal_nan)


    #!! SnpData.standardize() changes the str to help show that the data has been standardized. Should this to that too?
    def standardize(self, standardizer=DiagKtoN(), return_trained=False, force_python_only=False):
        """Does in-place standardization of the in-memory
        kernel data. The method multiples the values with a scalar factor such that the diagonal sums to iid_count. Although it works in place, for convenience
        it also returns the KernelData.

        :rtype: :class:`.KernelData` (standardizes in place, but for convenience, returns 'self')

        >>> from pysnptools.kernelreader import KernelNpz
        >>> import numpy as np
        >>> kernel_on_disk = KernelNpz('../examples/toydata.kernel.npz')
        >>> kerneldata1 = kernel_on_disk.read() # read all kernel values into memory
        >>> print(np.diag(kerneldata1.val).sum())
        5000000.0
        >>> kerneldata1.standardize() # standardize changes the values in kerneldata1.val
        KernelData(KernelNpz('../examples/toydata.kernel.npz'))
        >>> print(np.diag(kerneldata1.val).sum())
        500.0
        >>> kerneldata2 = kernel_on_disk.read().standardize() # Read and standardize in one expression with only one ndarray allocated.
        >>> print(np.diag(kerneldata2.val).sum())
        500.0
        """
        return standardizer.standardize(self, return_trained=return_trained, force_python_only=force_python_only)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
