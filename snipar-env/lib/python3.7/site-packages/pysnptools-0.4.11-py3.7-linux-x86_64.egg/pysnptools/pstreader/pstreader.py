import numpy as np
import subprocess, sys
import os.path
from itertools import *
import pandas as pd
import logging
import time
import pysnptools.util as pstutil
import numbers
try:
    from builtins import int
except:
    pass

class PstReader(object):
    """A PstReader is one of three things:

    * A :class:`.PstData` class that holds matrix data in memory:

        >>> import numpy as np
        >>> from pysnptools.pstreader import PstData
        >>> data1 = PstData(row=['a','b','c'],col=['y','z'],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        >>> type(data1.val).__name__ # The val property is an ndarray of values
        'ndarray'
        >>> data1.row_count #prints the number of rows in this in-memory data
        3

    * A class such as :class:`.PstNpz` for you to specify data in file. For example,

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> print(on_disk) # prints specification for reading from file
        PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> on_disk.col_count # prints the number of columns (but doesn't read any matrix values)
        1015

    * A subset of any PstReader, specified with "[ *row_index* , *col_index* ]", to read just some values.

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> subset_on_disk = on_disk[[3,4],::2] # specification for a subset of the data on disk. No values are read yet.
        >>> print(subset_on_disk.col_count) # prints the number of columns in this subset (but still doesn't read any values)
        508
        >>> print(subset_on_disk) #prints a specification of 'subset_on_disk'
        PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz')[[3,4],::2]
        >>> data_subset = subset_on_disk.read() # efficiently (if practical) reads the specified subset of values from the disk
        >>> print(data_subset) # prints the specification of the in-memory information
        PstData(PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz')[[3,4],::2])
        >>> int(data_subset.val.shape[0]),int(data_subset.val.shape[1]) # The dimensions of the ndarray of values
        (2, 508)

    The PstReaders Classes

        ================================== ================================== ====================== =====================
        *Class*                            *Format*                           *Random Access*        *Suffixes*
        :class:`.PstData`                  in-memory floats                   Yes                    *n/a*  
        :class:`.PstNpz`                   binary, floats                     No                     .pst.npz,.snp.npz,
                                                                                                     .kernel.npz
        :class:`.PstHdf5`                  binary, floats                     Yes                    .pst.hdf5,.snp.hdf5,
                                                                                                     .kernel.hdf5
        :class:`.PstMemMap`                mem-mapped floats                  Yes                    .pst.memmap,.snp.memmap
        various :class:`.SnpReader`        varies                             varies                 varies
        various :class:`.KernelReader`     varies                             varies                 varies
        ================================== ================================== ====================== =====================

        
            A :class:`.SnpReader` and :class:`.KernelReader` are each a kind of :class:`.PstReader`. They have some restrictions summarized here:
            ================================== =============== ============ ============ ==================== ====================
            *Class*                            *val type*      *row type*   *col type*   *row_property type*  *col_property type*
            :class:`.PstReader`                float           any          any          any                  any     
            :class:`.SnpReader`                float           str,str      str          none                 float,float,float
            :class:`.KernelReader`             float           str,str      str,str      none                 none
            ================================== =============== ============ ============ ==================== ====================

            For convenience, they allow additional ways to access rows and columns.

            ================================== =============== ============ ============ ==================== ====================
            *Class*                            *val name*      *row name*   *col name*   *row_property name*  *col_property name*
            :class:`.PstReader`                val             row          col          row_property         col_property
            :class:`.SnpReader`                val             iid,row      sid,col      none                 col_property,pos
            :class:`.KernelReader`             val             iid0,row,iid iid1,col,iid none                 none
            ================================== =============== ============ ============ ==================== ====================

            :Note: A :attr:`.KernelReader.iid` may be used when :attr:`.KernelReader.iid0` is equal to :attr:`.KernelReader.iid1`
  
    Methods & Properties:

        Every PstReader, such as :class:`.PstNpz` and :class:`.PstData`, has these properties: :attr:`row`, :attr:`row_count`, :attr:`col`, :attr:`col_count`,
        :attr:`row_property`, :attr:`col_property` and these methods: :meth:`read`, :meth:`row_to_index`, :meth:`col_to_index`. See below for details.

        :class:`.PstData` is a PstReader so it supports the above properties and methods. In addition, it supports property :attr:`.PstData.val` and equality testing.
        See below for details.

    Rows and Cols:

        Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.pstreader import PstHdf5
        >>> from pysnptools.util import print2 #print bytes strings and Unicode strings the same
        >>> on_disk = PstHdf5('../examples/toydata.iidmajor.snp.hdf5') # PstHdf5 can load .pst.hdf5, .snp.hdf5, and kernel.hdf5
        >>> print2(on_disk.row[:3]) # print the first three rows
        [['per0' 'per0']
         ['per1' 'per1']
         ['per2' 'per2']]
        >>> print2(on_disk.col[:7]) # print the first seven columns
        ['null_0' 'null_1' 'null_2' 'null_3' 'null_4' 'null_5' 'null_6']
        >>> print2(on_disk.row_to_index([[b'per2', b'per2'],[b'per1', b'per1']])) #Find the indexes for two rows.
        [2 1]
        
    When Matrix Data is Read:

        Matrix data can be enormous so we generally avoid reading it to the degree practical. Specifically,
        
        * Constructing and printing a PstReader causes no file reading. For example, these commands read no data:

            >>> on_disk = PstHdf5('../examples/toydata.iidmajor.snp.hdf5') # Construct a PstHdf5 PstReader. No data is read.
            >>> print(on_disk) # Print the PstHdf5 PstReader specification. No data is read.
            PstHdf5('../examples/toydata.iidmajor.snp.hdf5')
            >>> subset_on_disk = on_disk[[3,4],::2] # Construct a subsetting PstReader. No data is read.
            >>> print(subset_on_disk) # print the subset PstReader. No data is read.
            PstHdf5('../examples/toydata.iidmajor.snp.hdf5')[[3,4],::2]

        * Properties and methods related to the rows and columns (to the degree practical) read only row and col data from the disk,
          not value data. Moreover, the row and col data is read from file only once. Consider these commands:

            >>> on_disk = PstHdf5('../examples/toydata.iidmajor.snp.hdf5') # Construct a PstHdf5 PstReader. No data is read.
            >>> print2(on_disk.col[:7]) # without reading any values data from disk, read the row and col data from disk, cache it, and then print the first seven cols.
            ['null_0' 'null_1' 'null_2' 'null_3' 'null_4' 'null_5' 'null_6']
            >>> print(on_disk.col_to_index([b'null_7',b'null_2'])) #use the cached col information to find the indexes of 'null_7' and 'null_2'. (No data is read from disk.)
            [7 2]

        * The only method that reads values from file (to the degree practical) is :meth:`read`. For example:

            >>> on_disk = PstHdf5('../examples/toydata.iidmajor.snp.hdf5') # Construct a PstHdf5 PstReader. No data is read.
            >>> data1 = on_disk.read() #read all the values from disk, creating a new PstData instance that keeps these values in memory
            >>> print(data1.val[0,2]) # print the value for the row with index 0 and the col with index 2. (No data is read from disk.)
            2.0

        * If you request the values for only a subset of the rows or columns, (to the degree practical) only that subset will be read from disk.
          for example:

            >>> on_disk = PstHdf5('../../tests/datasets/all_chr.maf0.001.N300.snp.hdf5')[[3,4],::2] # Construct a subsetting PstReader. No data is read.
            >>> data_subset = subset_on_disk.read() # from disk, read the values for the rows with index 3 and 4 AND cols with even numbered indexes.
            >>> print(data_subset.val[0,2]) # print the value with subset row index 0 and col index 2 (corresponding to row index 3 and col index 4 in the full data). No data is read from disk.
            1.0

    When Matrix Data is Re-Read and Copied:

        Every time you call a PstReader's :meth:`read` method, the PstReader re-reads the value data and returns a new in-memory :class:`.PstData`
        (with :attr:`.PstData.val` property containing a new ndarray of the values).

        Here is an example of what not to do, because it causes all the matrix value data to be read twice.

            >>> on_disk = PstHdf5('../examples/toydata.iidmajor.snp.hdf5') # Construct a PstHdf5 PstReader. No data is read.
            >>> # Not recommended because it inefficiently reads all the values twice.
            >>> print(on_disk.read().val[0,2]) # read all values into a new PstData, print a value
            2.0
            >>> print(on_disk.read().val[0,3]) # read all values (again) into a second new PstData, print a value
            2.0

        Here are two efficient alternatives. First, if all values can all fit in memory, read them once into a :class:`PstData` and then
        access that :class:`PstData` multiple times.

            >>> on_disk = PstHdf5('../examples/toydata.iidmajor.snp.hdf5') # Construct a PstHdf5 PstReader. No data is read.
            >>> data1 = on_disk.read() # read all values into a new PstData
            >>> print(data1.val[0,2]) # print a value from data1's in-memory ndarray
            2.0
            >>> print(data1.val[0,3]) # print another value from data1's in-memory ndarray.
            2.0

        Second, if the value data is too large to fit in memory, use subsetting to read only the values of interest from disk.
       
            >>> on_disk = PstHdf5('../examples/toydata.iidmajor.snp.hdf5') # Construct a PstHdf5 PstReader. No data is read.
            >>> print(on_disk[0,2].read().val[0,0]) #Define the subset of data and read only that subset from disk.
            2.0
            >>> print(on_disk[0,3].read().val[0,0]) #Define a second subset of data and read only that subset from disk.
            2.0

        Because the in-memory :class:`.PstData` class is a kind of PstReader, you may read from it, too.
        Doing so create a new :class:`.PstData` instance containing a copy of the matrix values in a new ndarray.

            >>> on_disk = PstHdf5('../examples/toydata.iidmajor.snp.hdf5') # Construct a PstHdf5 PstReader. No data is read.
            >>> data1 = on_disk.read() # read all matrix values from disk into a new PstData
            >>> print(data1.val is data1.val) # Do the in-memory SNP values use the same memory as themselves? Yes
            True
            >>> data2 = data1.read() # copy all the matrix values into a new ndarray in a new PstData
            >>> print(data2.val is data1.val) # Do the two ndarrays of in-memory matrix values use the same memory?
            False


    Avoiding Unwanted ndarray Allocations

        You may want a subset of matrix values from an in-memory :class:`PstData` and you may know that this subset and the original :class:`PstData`
        can safely share the memory of the ndarray of matrix values. For this case, the :meth:`read` has optional parameters called view_ok and order. If you override 
        the defaults of "view_ok=False,order='F'" with "view_ok=True,order='A', the :meth:`read` will, if practical, return a new 
        :class:`PstData` with a ndarray shares memory with the original ndarray.
        Use these parameters with care because any change to either ndarray will effect
        the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
        share memory and so it may ignore your suggestion and allocate a new ndarray anyway.

            >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Construct a PstNpz PstReader. No data is read.
            >>> data1 = on_disk.read() # read all data from disk into a PstData with a new ndarray
            >>> column01 = data1[:,0:1].read(view_ok=True,order='A') #create PstData with the data from just the first two columns. Sharing memory is OK. The memory may be laid out in any order (that is col-major and row-major are both OK).
            >>> import numpy as np
            >>> #print(np.may_share_memory(data1.val, column01.val)) # Do the two ndarray's share memory? They could (but currently they won't)
            >>> column201 = data1[:,[2,0,1]].read(view_ok=True,order='A') #create PstData with the data from three columns, permuted. Sharing memory is OK.
            >>> print(np.may_share_memory(data1.val, column201.val)) # Do the two ndarray's share memory? No, ndarray decided that this indexing was too complex for sharing.
            False

    Creating Subsetting PstReaders with Indexing

        You often don't want to read the matrix values for all rows and cols. You can use indexing to create a subsetting PstReader that
        will read only the matrix values of interest.

        PstReaders support the indexing formats supported by ndarray plus two generalizations. Here are examples of indexing with an array
        of indexes, with slicing, and with an array of Booleans.

            >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify some data on disk in PstNpz format
            >>> subset_reader_1 = on_disk[[3,-1],:] #index with an array of indexes (negatives count from end)
            >>> print((subset_reader_1.row_count, subset_reader_1.col_count))
            (2, 1015)
            >>> data1 = subset_reader_1.read() # read just the two rows of interest from the disk
            >>> subset_reader_2 = on_disk[:,:0:-2] #index with a slice
            >>> print((subset_reader_2.row_count, subset_reader_2.col_count))
            (300, 507)
            >>> boolindexes = [s.startswith(b'23_') for s in on_disk.col] # create a Boolean index of cols that start '23_'
            >>> subset_reader_3 = on_disk[:,boolindexes] #index with array of Booleans
            >>> print((subset_reader_3.row_count, subset_reader_3.col_count))
            (300, 24)

        The first generalization over what ndarray offers is full indexing on both the row dimension and the col dimension, in other words,
        full multidimensional indexing. For example,

            >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify some data on disk in PstNpz format
            >>> subset_reader_4 = on_disk[[3,4],:0:-2] # index on two dimensions at once
            >>> print((subset_reader_4.row_count, subset_reader_4.col_count))
            (2, 507)

        The second generalization is indexing on a single integer index.

            >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify some data on disk in PstNpz format
            >>> subset_reader_5 = on_disk[5,:] #index with single integer
            >>> print((subset_reader_5.row_count, subset_reader_5.col_count))
            (1, 1015)

        Indexing is also useful when you have matrix values in memory via a :class:`PstData` index and want to copy a subset of those values.
        While you could instead index directly on the `.PstData.val` ndarray, by indexing on the :class:`PstData` instance you
        also get row and col information.

            >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify some data on disk in PstNpz format
            >>> data1 = on_disk.read() # read all matrix values into memory
            >>> print2(data1.col[:9]) # print the first 9 cols
            ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4']

            >>> data_subset = data1[:,::2].read(view_ok=True,order='A') # create a copy or view with every other col
            >>> print2(data_subset.col[:9]) # print the first 9 cols in the subset
            ['1_12' '1_10' '1_28' '1_36' '1_4' '1_11' '1_32' '1_9' '1_17']


        You can apply indexing on top of indexing to specify subsets of subsets of data to read. In this example, 
        only the column values for every 16th col is actually read from the disk.

            >>> # These are just PstReaders, nothing is read from disk yet
            >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify some data on disk in PstNpz format
            >>> half_reader = on_disk[:,::2] # a reader for half the cols
            >>> quarter_reader = half_reader[:,::2] # a reader for half of half the cols
            >>> sixteenth_reader = quarter_reader[:,::2][:,::2] # a reader for half of half of half of half the cols
            >>> print(sixteenth_reader) #Print the specification of this reader
            PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz')[:,::2][:,::2][:,::2][:,::2]
            >>> # Now we read from disk. Only values for one col in every 16 will be read.
            >>> data_sixteenth = sixteenth_reader.read()
            >>> print(data_sixteenth.val[0,3])
            2.0

    The :meth:`read` Method
  
        By default the :meth:`read` returns a ndarray of scipy.float64 laid out in memory in F-contiguous order (row-index varies the fastest). You may, instead,
        ask for scipy.float32 or for C-contiguous order or any order. See :meth:`read` for details.

    Details of Methods & Properties:
    """

    def __init__(self, *args, **kwargs): #Ignore any inputs because our parent is 'object'
        super(PstReader, self).__init__()

    @property
    def row(self):
        """A ndarray of the row ids. Each id can be anything, for example, a string, an array of two strings, a number, etc.

        :rtype: ndarray (length :attr:`.row_count`)

        This property (to the degree practical) reads only row and col data from the disk, not matrix value data.
        Moreover, the row and col data is read from file only once.

        :Example:

        >>> import numpy as np
        >>> from pysnptools.pstreader import PstData
        >>> data1 = PstData(row=['a','b','c'],col=['y','z'],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        >>> print(data1.row[:2]) # print the first two row ids
        ['a' 'b']
        >>> data2 = PstData(row=[[1,0],[2,0],[2,1]],col=['y','z'],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        >>> print(data2.row[:2]) # print the first two row ids
        [[1 0]
         [2 0]]
        """
        raise NotImplementedError

    @property
    def row_count(self):
        """number of rows

        :rtype: integer

        This property (to the degree practical) reads only row and col data from the disk, not matrix value data. Moreover, the row and col data is read from file only once.
        """
        return len(self.row)

    @property
    def col(self):
        """A ndarray of the cols id. Each id can be anything, for example, a string, an array of two strings, a number, etc.

        :rtype: ndarray (length :attr:`.col_count`)

        This property (to the degree practical) reads only row and col data from the disk, not matrix value data.
        Moreover, the row and col data is read from file only once.

        :Example:

        >>> import numpy as np
        >>> from pysnptools.pstreader import PstData
        >>> data1 = PstData(row=['a','b','c'],col=['y','z'],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        >>> print(data1.col[-1]) # print the last column id.
        z
        >>> data2 = PstData(row=['a','b','c'],col=[[1,0],[2,0]],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        >>> print(data2.col[-1]) # print the last column id.
        [2 0]
        """
        raise NotImplementedError

    @property
    def col_count(self):
        """number of cols

        :rtype: integer

        This property (to the degree practical) reads only row and col data from the disk, not matrix value data. Moreover, the row and col data is read from file only once.

        """
        return len(self.col)

    @property
    def shape(self):
        """number of rows and number of cols

        :rtype: tuple of two integers

        This property (to the degree practical) reads only row and col data from the disk, not matrix value data. Moreover, the row and col data is read from file only once.

        """
        return (len(self.row),len(self.col))

    @property
    def row_property(self):
        """A ndarray of the additional information for each row. Each element is a ndarray.

        :rtype: 1- or 2-dimensional ndarray (length :attr:`.row_count`)

        This property (to the degree practical) reads only row and col data from the disk, not matrix value data. Moreover, the row and col data is read from file only once.

        :Example:

        >>> import numpy as np
        >>> from pysnptools.pstreader import PstData
        >>> data1 = PstData(row=['a','b','c'],col=['y','z'],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        >>> print(data1.row_property[:2]) # print the first two row property values
        ['A' 'B']
        """
        raise NotImplementedError

    @property
    def col_property(self):
        """A ndarray of the additional information for each col. Each element is a ndarray.

        :rtype: 1- or 2-dimensional ndarray (length :attr:`.col_count`)

        This property (to the degree practical) reads only row and col data from the disk, not matrix value data. Moreover, the row and col data is read from file only once.

        :Example:

        :Example:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> print(on_disk.col_property[:3]) # print column information for the first three cols: #The '...' is an optional space
        [[...1.          0.00800801  0.        ]
         [...1.          0.023023    1.        ]
         [...1.          0.0700701   4.        ]]
        """
        raise NotImplementedError


    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        raise NotImplementedError



    #!!check that views always return contiguous memory by default
    def read(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False):
        """Reads the matrix values and returns a :class:`.PstData` (with :attr:`.PstData.val` property containing a new ndarray of the matrix values).

        :param order: {'F' (default), 'C', 'A'}, optional -- Specify the order of the ndarray. If order is 'F' (default),
            then the array will be in F-contiguous order (row-index varies the fastest).
            If order is 'C', then the returned array will be in C-contiguous order (col-index varies the fastest).
            If order is 'A', then the :attr:`.PstData.val`
            ndarray may be in any order (either C-, Fortran-contiguous).
        :type order: string or None

        :param dtype: {scipy.float64 (default), scipy.float32}, optional -- The data-type for the :attr:`.PstData.val` ndarray.
        :type dtype: data-type

        :param force_python_only: optional -- If False (default), may use outside library code. If True, requests that the read
            be done without outside library code.
        :type force_python_only: bool


        :param view_ok: optional -- If False (default), allocates new memory for the :attr:`.PstData.val`'s ndarray. If True,
            if practical and reading from a :class:`PstData`, will return a new 
            :class:`PstData` with a ndarray shares memory with the original :class:`PstData`.
            Typically, you'll also wish to use "order='A'" to increase the chance that sharing will be possible.
            Use these parameters with care because any change to either ndarray will effect
            the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
            share memory and so it may ignore your suggestion and allocate a new ndarray anyway.
        :type view_ok: bool

        :rtype: :class:`.PstData`

        Calling the method again causes the matrix values to be re-read and creates a new in-memory :class:`.PstData` with a new ndarray of matrix values.

        If you request the values for only a subset of the sids or iids, (to the degree practical) only that subset will be read from disk.

        :Example:

        >>> from pysnptools.pstreader import PstHdf5
        >>> on_disk = PstHdf5('../examples/toydata.iidmajor.snp.hdf5') # Specify matrix data on disk
        >>> pstdata1 = on_disk.read() # Read all the matrix data returning a PstData instance
        >>> print(type(pstdata1.val).__name__) # The PstData instance contains a ndarray of the data.
        ndarray
        >>> subset_pstdata = on_disk[:,::2].read() # From the disk, read matrix values for every other sid
        >>> print(subset_pstdata.val[0,0]) # Print the first matrix value in the subset
        1.0
        >>> subsub_pstdata = subset_pstdata[:10,:].read(order='A',view_ok=True) # Create an in-memory subset of the subset with matrix values for the first ten iids. Share memory if practical.
        >>> import numpy as np
        >>> # print(np.may_share_memory(subset_snpdata.val, subsub_snpdata.val)) # Do the two ndarray's share memory? They could. Currently they won't.       
        """
        val = self._read(None, None, order, dtype, force_python_only, view_ok)
        from pysnptools.pstreader import PstData
        ret = PstData(self.row, self.col, val, row_property=self.row_property, col_property=self.col_property, name=str(self))
        return ret

    def row_to_index(self, list):
        """Takes a list of row ids and returns a list of index numbers

        :param list: list of rows
        :type order: list

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only row and col data from the disk, not matrix value data. Moreover, the row and col data is read from file only once.

        :Example:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify matrix data on disk
        >>> print(on_disk.row_to_index([[b'POP1',b'44'],[b'POP1',b'12']])) #Find the indexes for two rows.
        [2 1]
        """
        if not hasattr(self, "_row_to_index"):
            self._row_to_index = {}
            for index, item in enumerate(self.row):
                key = self._makekey(item)
                if key in self._row_to_index:
                   raise Exception("Expect row to appear in data only once. ({0})".format(key))
                self._row_to_index[key] = index
        index = np.fromiter((self._row_to_index[PstReader._makekey(item1)] for item1 in list),np.int)
        return index

    def col_to_index(self, list):
        """Takes a list of column ds and returns a list of index numbers

        :param list: list of cols
        :type list: list

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only row and col data from the disk, not matrix value data. Moreover, the row and col data is read from file only once.

        :Example:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify matrix data on disk
        >>> print(on_disk.col_to_index([b'1_10',b'1_13'])) #Find the indexes for two cols.
        [2 9]
        """
        if not hasattr(self, "_col_to_index"):
            logging.debug("Creating _col_to_index")
            col_set = None
            try:
                col_set = set(self.col)
            except:
                pass

            if col_set is not None:
                assert len(col_set) == self.col_count, "Expect col to appear in data only once."
                self._col_to_index = {item : index for index, item in enumerate(self.col)}
            else:
                col_list = [PstReader._makekey(item) for item in self.col]
                assert len(set(col_list)) == self.col_count, "Expect col to appear in data only once."
                self._col_to_index = {item : index for index, item in enumerate(col_list)}
            logging.debug("Finished creating _col_to_index")
        index = np.fromiter((self._col_to_index[PstReader._makekey(item1)] for item1 in list),np.int)
        return index

    @staticmethod
    def _makekey(item):
        if isinstance(item,str):
            return item
        if isinstance(item,(numbers.Integral,float)): #return quickly from known items
            return item
        try:
            hash(item)
            return item
        except:
            pass

        return tuple((PstReader._makekey(subitem) for subitem in item))

    def __getitem__(self, row_indexer_and_col_indexer):
        from pysnptools.pstreader._subset import _PstSubset
        row_indexer, col_indexer = row_indexer_and_col_indexer
        return _PstSubset(self, row_indexer, col_indexer)


    def copyinputs(self, copier):
        raise NotImplementedError

    @staticmethod
    def _is_all_slice(index_or_none):
        if index_or_none is None:
            return True
        return  isinstance(index_or_none,slice) and index_or_none == slice(None)

    @staticmethod
    def _make_sparray_or_slice(indexer):
        if indexer is None:
            return slice(None)

        if isinstance(indexer,np.ndarray):
            return PstReader._process_ndarray(indexer)

        if isinstance(indexer, slice):
            return indexer

        if np.isscalar(indexer):
            assert isinstance(indexer, numbers.Integral), "Expect scalar indexes to be integers"
            return np.array([indexer])

        return PstReader._process_ndarray(np.array(indexer))

    @staticmethod
    def _process_ndarray(indexer):
        if len(indexer)==0: # If it's zero length, the type is unreliable and unneeded.
            return np.zeros((0),dtype=np.integer)
        if indexer.dtype == bool:
            return np.arange(len(indexer),dtype=np.integer)[indexer]
        assert np.issubdtype(indexer.dtype, np.integer), "Indexer of unknown type"
        return indexer


    @staticmethod
    def _make_sparray_from_sparray_or_slice(count, indexer):
        if isinstance(indexer,slice):
            return np.arange(*indexer.indices(count))
        if isinstance(indexer,np.ndarray) and np.issubdtype(indexer.dtype, np.integer) and np.any(indexer<0):
            return np.arange(count)[indexer]
        return indexer

    @staticmethod
    def _array_properties_are_ok(val, order, dtype):
        if val.dtype != dtype:
            return False
        if order is 'F':
            return val.flags['F_CONTIGUOUS']
        elif order is 'C':
            return val.flags['C_CONTIGUOUS']

        return True

    def _apply_sparray_or_slice_to_val(self, val, row_indexer_or_none, col_indexer_or_none, order, dtype, force_python_only):
        if (PstReader._is_all_slice(row_indexer_or_none) and PstReader._is_all_slice(col_indexer_or_none)  and not force_python_only and 
                (order == 'A' or (order == 'F' and val.flags['F_CONTIGUOUS']) or (order == 'C' and val.flags['C_CONTIGUOUS'])) and
                (dtype is None or  val.dtype == dtype)):
            return val, True

        row_indexer = PstReader._make_sparray_or_slice(row_indexer_or_none)
        col_indexer = PstReader._make_sparray_or_slice(col_indexer_or_none)
        if not force_python_only:
            row_index = PstReader._make_sparray_from_sparray_or_slice(self.row_count, row_indexer)
            col_index = PstReader._make_sparray_from_sparray_or_slice(self.col_count, col_indexer)
            sub_val = pstutil.sub_matrix(val, row_index, col_index, order=order, dtype=dtype)
            return sub_val, False


        if PstReader._is_all_slice(row_indexer) or PstReader._is_all_slice(col_indexer):
            sub_val = val[row_indexer, col_indexer] #!!is this faster than the C++?
        else: 
            row_index = PstReader._make_sparray_from_sparray_or_slice(self.row_count, row_indexer)
            col_index = PstReader._make_sparray_from_sparray_or_slice(self.col_count, col_indexer)
            #See http://stackoverflow.com/questions/21349133/numpy-array-integer-indexing-in-more-than-one-dimension
            sub_val = val[row_index.reshape(-1,1), col_index]

        assert len(sub_val.shape)==2, "Expect result of subsetting to be 2 dimensional"

        if not PstReader._array_properties_are_ok(sub_val, order, dtype):
            if order is None:
                order = "K"
            if dtype is None:
                dtype = sub_val.dtype
            sub_val = sub_val.astype(dtype, order, copy=True)

        shares_memory =  np.may_share_memory(val, sub_val)
        assert(PstReader._array_properties_are_ok(sub_val, order, dtype))
        return sub_val, shares_memory

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
    print("done")