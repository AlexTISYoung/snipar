from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import subprocess, sys
import os.path
from itertools import *
import pandas as pd
import logging
import time
import pysnptools.util as pstutil
from pysnptools.pstreader import PstReader
from pysnptools.kernelstandardizer import DiagKtoN

class KernelReader(PstReader):
    """A KernelReader is one of three things:

    * A class such as :class:`KernelNpz` for you to specify data in a file. For example,

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.kernelreader import KernelNpz
        >>> 
        >>> kernel_on_disk = KernelNpz('../examples/toydata.kernel.npz')
        >>> print(kernel_on_disk) # prints specification for reading from file
        KernelNpz('../examples/toydata.kernel.npz')
        >>> kernel_on_disk.iid_count # prints the number of iids (but doesn't read any kernel values)
        500

    * A :class:`.KernelData` class that holds kernel data in memory, typically after computing from a SnpReader or reading it from a KernelReader:

        >>> # Compute kernel from a SnpReader
        >>> from pysnptools.snpreader import Bed
        >>> from pysnptools.standardizer import Unit
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)
        >>> kerneldata1 = snp_on_disk.read_kernel(Unit()) #reads the SNP values and computes the kernel
        >>> type(kerneldata1.val).__name__ # The val property is an ndarray of kernel values
        'ndarray'
        >>> print(kerneldata1) # prints the specification of the in-memory kernel information
        KernelData(SnpKernel(Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False),standardizer=Unit()))
        >>> kerneldata1.iid_count #prints the number of iids (number of individuals) in this in-memory data
        300
        >>> # Read kernel from a KernelReader
        >>> kernel_on_disk = KernelNpz('../examples/toydata.kernel.npz')
        >>> kerneldata2 = kernel_on_disk.read() #reads the kernel values
        >>> print(kerneldata2) # prints the specification of the in-memory kernel information
        KernelData(KernelNpz('../examples/toydata.kernel.npz'))
        >>> kerneldata2.iid_count #prints the number of iids (number of individuals) in this in-memory data
        500


    * A subset of any KernelReader, specified with "[ *iid_index* ]" (or  specified with "[ *iid0_index* , *iid1_index* ]"), to read only some kernel values. It can
      also be used to re-order the values.

        >>> kernel_on_disk = KernelNpz('../examples/toydata.kernel.npz')
        >>> subset_on_disk1 = kernel_on_disk[[3,4]] # specification for a subset of the data on disk. No kernel values are read yet.
        >>> print(subset_on_disk1.iid_count) # prints the number of iids in this subset (but still doesn't read any kernel values)
        2
        >>> print(subset_on_disk1) #prints a specification of 'subset_on_disk1'
        KernelNpz('../examples/toydata.kernel.npz')[[3,4],[3,4]]
        >>> kerneldata_subset = subset_on_disk1.read() # efficiently (if possible) reads the specified subset of values from the disk
        >>> print(kerneldata_subset) # prints the specification of the in-memory kernel information
        KernelData(KernelNpz('../examples/toydata.kernel.npz')[[3,4],[3,4]])
        >>> print((int(kerneldata_subset.val.shape[0]), int(kerneldata_subset.val.shape[1]))) # The dimensions of the ndarray of kernel values
        (2, 2)
        >>> subset_on_disk2 = kernel_on_disk[[3,4],::2] # specification for a subset of the data on disk. No kernel values are read yet.
        >>> print((subset_on_disk2.iid0_count, subset_on_disk2.iid1_count))
        (2, 250)


    The KernelReaders Classes

        ================================== ================== ====================== ================== ====================
        *Class*                            *Format*           *Random Access*        *Suffixes*         *Write* method?
        :class:`.KernelData`               in-memory          Yes                    *n/a*              *n/a*
        :class:`.KernelNpz`                binary             No                     .kernel.npz        Yes
        :class:`.KernelHdf5`               binary             Yes                    .kernel.hdf5       Yes
        :class:`.Identity`                 *n/a*              Yes                    *n/a*              No
        :class:`.SnpKernel`                depends            depends                *n/a*              No
        ================================== ================== ====================== ================== ====================
    
  
    Methods & Properties:

        Every KernelReader, such as :class:`.KernelNpz` and :class:`.KernelData`, when square has these properties: :attr:`iid`, :attr:`iid_count`, 
        and these methods: :meth:`read`, and :meth:`iid_to_index`. A square kernel is one that has the same iid list for both its rows and columns.

        More generally, KernelReaders can have one iid list for its rows and a different iid list for its columns, so these properties and methods are also defined: :attr:`iid0`, :attr:`iid1`, :attr:`iid0_count`, 
        :attr:`iid1_count`, :meth:`iid0_to_index`, and :meth:`iid1_to_index`.
       
        See below for details.

        :class:`.KernelData` is a KernelReader so it supports the above properties and methods. In addition, it supports property :attr:`.KernelData.val`, method :meth:`.KernelData.standardize`, and equality testing.
        See below for details.

        Some of the classes, such as :class:`.KernelNpz`, also provide a static :meth:`KernelNpz.write` method for writing :class:`.KernelData`.

        >>> # create a kernel from a Bed file and write to KernelNpz format
        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import Bed
        >>> from pysnptools.standardizer import Unit
        >>> import pysnptools.util as pstutil
        
        >>> kerneldata = Bed('../examples/toydata.bed',count_A1=False).read_kernel(Unit())     # Create a kernel from the data in the Bed file
        >>> pstutil.create_directory_if_necessary("tempdir/toydata.kernel.npz")
        >>> KernelNpz.write("tempdir/toydata.kernel.npz",kerneldata)      # Write data in KernelNpz format
        KernelNpz('tempdir/toydata.kernel.npz')

    iids:

        Individual are identified with an iid, which is a ndarray of two strings: a family ID and a case ID. For example:

        >>> kernel_on_disk = KernelNpz('../examples/toydata.kernel.npz')
        >>> print(kernel_on_disk.iid[:3]) # print the first three iids
        [['per0' 'per0']
         ['per1' 'per1']
         ['per2' 'per2']]
        >>> print(kernel_on_disk.iid_to_index([['per2','per2'],['per1','per1']])) #Find the indexes for two iids.
        [2 1]

    :class:`.KernelReader` is a kind of :class:`.PstReader`. See the documentation for :class:`.PstReader` to learn about:
        
        * When Data is Read
        * When Data is Re-Read and Copied
        * Avoiding Unwanted ndarray Allocations
        * Creating Subsetting PstReaders with Indexing

    The :meth:`read` Method
  
        By default the :meth:`read` returns a ndarray of scipy.float64 laid out in memory in F-contiguous order (iid0-index varies the fastest). You may, instead,
        ask for scipy.float32 or for C-contiguous order or any order. See :meth:`read` for details.

    The :meth:`.KernelData.standardize` Method
        The :meth:`.KernelData.standardize` method, available only on :class:`.KernelData`, does in-place standardization of the in-memory
        kernel data. The method multiples the values with a scalar factor such that the diagonal sums to iid_count. Although it works in place, for convenience
        it also returns itself. See :meth:`.KernelData.standardize` for details.

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
   
    Details of Methods & Properties:
    """
    def __init__(self, *args, **kwargs):
        super(KernelReader, self).__init__(*args, **kwargs)

    @property
    def iid(self):
        """A ndarray of the iids. Each iid is a ndarray of two strings (a family ID and a case ID) that identifies an individual.
        Assumes the kernel is square, so will throw an exception if the row iids are different from the column iids.

        :rtype: ndarray (length :attr:`.iid_count`) of ndarray (length 2) of strings

        This property (to the degree practical) reads only iid and sid data from the disk, not kernel value data. Moreover, the iid data is read from file only once.

        :Example:

        >>> from pysnptools.kernelreader import KernelNpz
        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> kernel_on_disk = KernelNpz('../examples/toydata.kernel.npz')
        >>> print(kernel_on_disk.iid[:3]) # print the first three iids
        [['per0' 'per0']
         ['per1' 'per1']
         ['per2' 'per2']]
        """
        assert self.iid0 is self.iid1, "When 'iid' is used, iid0 must be the same as iid1"
        return self.iid0

    @property
    def iid0(self):
        """
        A ndarray of the row iids. See :attr:`.iid`
        """
        return self.row

    @property
    def iid1(self):
        """
        A ndarray of the column iids. See :attr:`.iid`
        """
        return self.col

    @property
    def iid_count(self):
        """number of iids
        Assumes the kernel is square, so will throw an exception if the row iids are different from the column iids.

        :rtype: integer

        This property (to the degree practical) reads only iid data from the disk, not kernel value data. Moreover, the iid data is read from file only once.
        """
        assert self.iid0 is self.iid1, "When 'iid_count' is used, iid0 must be the same as iid1"
        return self.iid0_count

    @property
    def iid0_count(self):
        """number of row iids. See :attr:`iid_count`

        :rtype: integer
        """
        return self.row_count

    @property
    def iid1_count(self):
        """number of column iids. See :attr:`iid_count`

        :rtype: integer
        """
        return self.col_count


    @property
    def row_property(self):
        """
        Defined for compatibility with :class:`.PstReader`. Will always be empty.
        """
        if not hasattr(self,"_row_property"):
            self._row_property = np.empty((self.row_count,0))
        return self._row_property

    @property
    def col_property(self):
        """
        Defined for compatibility with :class:`.PstReader`. Will always be empty.
        """
        if not hasattr(self,"_col_property"):
            self._col_property = np.empty((self.col_count,0))
        return self._col_property



    #!!check that views always return contiguous memory by default
    def read(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False):
        """Reads the kernel values and returns a :class:`.KernelData` (with :attr:`.KernelData.val` property containing a new ndarray of the kernel values).

        :param order: {'F' (default), 'C', 'A'}, optional -- Specify the order of the ndarray. If order is 'F' (default),
            then the array will be in F-contiguous order (iid0-index varies the fastest).
            If order is 'C', then the returned array will be in C-contiguous order (iid1-index varies the fastest).
            If order is 'A', then the :attr:`.KernelData.val`
            ndarray may be in any order (either C-, Fortran-contiguous, or even discontiguous).
        :type order: string or None

        :param dtype: {scipy.float64 (default), scipy.float32}, optional -- The data-type for the :attr:`.KernelData.val` ndarray.
        :type dtype: data-type

        :param force_python_only: optional -- If False (default), may use outside library code. If True, requests that the read
            be done without outside library code.
        :type force_python_only: bool


        :param view_ok: optional -- If False (default), allocates new memory for the :attr:`.KernelData.val`'s ndarray. If True,
            if practical and reading from a :class:`KernelData`, will return a new 
            :class:`KernelData` with a ndarray shares memory with the original :class:`KernelData`.
            Typically, you'll also wish to use "order='A'" to increase the chance that sharing will be possible.
            Use these parameters with care because any change to either ndarray (for example, via :meth:`.KernelData.standardize`) will effect
            the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
            share memory and so it may ignore your suggestion and allocate a new ndarray anyway.
        :type view_ok: bool

        :rtype: :class:`.KernelData`

        Calling the method again causes the kernel values to be re-read and creates a new in-memory :class:`.KernelData` with a new ndarray of kernel values.

        If you request the values for only a subset of the sids or iids, (to the degree practical) only that subset will be read from disk.

        :Example:

        >>> from pysnptools.kernelreader import KernelNpz
        >>> kernel_on_disk = KernelNpz('../examples/toydata.kernel.npz')
        >>> kerneldata1 = kernel_on_disk.read() # Read all the kernel data returning a KernelData instance
        >>> print(type(kerneldata1.val).__name__) # The KernelData instance contains a ndarray of the data.
        ndarray
        >>> subset_kerneldata = kernel_on_disk[::2].read() # From the disk, read kernel values for every other iid
        >>> print('{0:.6f}'.format(subset_kerneldata.val[0,0])) # Print the first kernel value in the subset
        9923.069928
        >>> subsub_kerneldata = subset_kerneldata[:10].read(order='A',view_ok=True) # Create an in-memory subset of the subset with kernel values for the first ten iids. Share memory if practical.
        >>> import numpy as np
        >>> #print(np.may_share_memory(subset_kerneldata.val, subsub_kerneldata.val)) # Do the two ndarray's share memory? They could. Currently they won't.       
        """
        val = self._read(None, None, order, dtype, force_python_only, view_ok)
        from pysnptools.kernelreader import KernelData
        ret = KernelData(iid0=self.iid0, iid1=self.iid1, val=val, name=str(self))
        return ret

    def iid_to_index(self, list):
        """Takes a list of iids and returns a list of index numbers.
        Assumes the kernel is square, so will throw an exception if the row iids are different from the column iids.

        :param list: list of iids
        :type order: list of list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only iid from the disk, not kernel value data. Moreover, the iid data is read from file only once.

        :Example:

        >>> from pysnptools.kernelreader import KernelNpz
        >>> kernel_on_disk = KernelNpz('../examples/toydata.kernel.npz')
        >>> print(kernel_on_disk.iid_to_index([['per2','per2'],['per1','per1']])) #Find the indexes for two iids.
        [2 1]
        """
        assert self.iid0 is self.iid1, "When 'iid_to_index' is used, iid0 must be the same as iid1"
        return self.iid0_to_index(list)

    def iid0_to_index(self, list):
        """Takes a list of row iids and returns a list of index numbers. See :attr:`iid_to_index`
        """
        return self.row_to_index(list)

    @staticmethod
    def _makekey(item):
        return tuple(str(i) for i in item)


    def iid1_to_index(self, list):
        """Takes a list of column iids and returns a list of index numbers. See :attr:`iid_to_index`
        """
        return self.col_to_index(list)

    def __getitem__(self, iid_indexer_and_snp_indexer):
        from pysnptools.kernelreader._subset import _KernelSubset
        if isinstance(iid_indexer_and_snp_indexer,tuple):
            iid0_indexer, iid1_indexer = iid_indexer_and_snp_indexer
        else:
            iid0_indexer = iid_indexer_and_snp_indexer
            iid1_indexer = iid0_indexer

        return _KernelSubset(self, iid0_indexer, iid1_indexer)

    def _assert_iid0_iid1(self):
        assert self._row.dtype.type is np.str_ and len(self._row.shape)==2 and self._row.shape[1]==2, "iid0 should be dtype str, have two dimensions, and the second dimension should be size 2"
        assert self._col.dtype.type is np.str_ and len(self._col.shape)==2 and self._col.shape[1]==2, "iid1 should be dtype str have two dimensions, and the second dimension should be size 2"

    def _read_with_standardizing(self, to_kerneldata, snp_standardizer=None, kernel_standardizer=DiagKtoN(), return_trained=False):
        assert to_kerneldata, "When working with non-SnpKernels, to_kerneldata must be 'True'"
        kernel, kernel_trained = self.read().standardize(kernel_standardizer,return_trained=True)

        if return_trained:
            return kernel, None, kernel_trained
        else:
            return kernel


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
    print("done")