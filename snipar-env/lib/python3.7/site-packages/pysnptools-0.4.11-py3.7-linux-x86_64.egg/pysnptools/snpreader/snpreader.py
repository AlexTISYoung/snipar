import numpy as np
import subprocess, sys
import os.path
from itertools import *
import pandas as pd
import logging
import time
import pysnptools.util as pstutil
from pysnptools.pstreader import PstReader
import warnings
import pysnptools.standardizer as stdizer
from six.moves import range

#!!why do the examples use ../tests/datasets instead of "examples"?
class SnpReader(PstReader):
    """A SnpReader is one of three things:

    * A class such as :class:`.Bed` for you to specify data in a file. For example,

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)
        >>> print(snp_on_disk) # prints specification for reading from file
        Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)
        >>> snp_on_disk.sid_count # prints the number of SNPS (but doesn't read any SNP values)
        1015

    * A :class:`.SnpData` class that holds SNP data in memory, typically after reading it from disk:

        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)
        >>> snpdata1 = snp_on_disk.read() #reads the SNP values
        >>> type(snpdata1.val).__name__ # The val property is an ndarray of SNP values
        'ndarray'
        >>> print(snpdata1) # prints the specification of the in-memory SNP information
        SnpData(Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False))
        >>> snpdata1.iid_count #prints the number of iids (number of individuals) in this in-memory data
        300

    * A subset of any SnpReader, specified with "[ *iid_index* , *sid_index* ]", to read only some SNP values. It can
      also be used to re-order the values.

        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)
        >>> subset_on_disk = snp_on_disk[[3,4],::2] # specification for a subset of the data on disk. No SNP values are read yet.
        >>> print(subset_on_disk.sid_count) # prints the number of sids in this subset (but still doesn't read any SNP values)
        508
        >>> print(subset_on_disk) #prints a specification of 'subset_on_disk'
        Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)[[3,4],::2]
        >>> snpdata_subset = subset_on_disk.read() # efficiently reads the specified subset of values from the disk
        >>> print(snpdata_subset) # prints the specification of the in-memory SNP information
        SnpData(Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)[[3,4],::2])
        >>> print((int(snpdata_subset.val.shape[0]), int(snpdata_subset.val.shape[1]))) # The dimensions of the ndarray of SNP values
        (2, 508)


    The SnpReaders Classes

        ========================= =================== ====================== ================== ======================
        *Class*                   *Format*            *Random Access*        *Suffixes*         *Write* method?
        :class:`.SnpData`         in-memory floats    Yes                    *n/a*              *n/a*              
        :class:`.Bed`             binary, 0,1,2       Yes (by sid)           .bed/.bim/.fam     Yes
        :class:`.Pheno`           text, floats        No                     .txt,.phe          Yes
        :class:`.Dat`             text, floats        No                     .dat/.map/.fam     Yes
        :class:`.Ped`             text, 0,1,2         No                     .ped/.map          Yes
        :class:`.Dense`           text, 0,1,2         No                     .dense.txt         Yes
        :class:`.SnpNpz`          binary, floats      No                     .snp.npz           Yes
        :class:`.SnpHdf5`         binary, floats      Yes (by sid or iid)    .snp.hdf5          Yes
        :class:`.SnpMemMap`       mem-mapped floats   Yes                    .snp.memmap        Yes              
        :class:`.SnpGen`          generated values    Yes (by sid)           *n/a*              *n/a*              
        ========================= =================== ====================== ================== ======================
    
  
    Methods & Properties:

        Every SnpReader, such as :class:`.Bed` and :class:`.SnpData`, has these properties: :attr:`iid`, :attr:`iid_count`, :attr:`sid`, :attr:`sid_count`,
        :attr:`pos` and these methods: :meth:`read`, :meth:`iid_to_index`, :meth:`sid_to_index`, :meth:`read_kernel`. See below for details.

        :class:`.SnpData` is a SnpReader so it supports the above properties and methods. In addition, it supports property :attr:`.SnpData.val`, method :meth:`.SnpData.standardize`, and equality testing.
        See below for details.

        Many of the classes, such as :class:`.Bed`, also provide a static :meth:`Bed.write` method for writing :class:`.SnpData` to disk.

        >>> # read from Pheno, write to Bed
        >>> from pysnptools.snpreader import Pheno, Bed
        >>> import pysnptools.util as pstutil
        
        >>> snpdata = Pheno('../examples/toydata.phe').read() # Read data from Pheno format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata.bed")
        >>> Bed.write("tempdir/toydata.bed",snpdata,count_A1=False)   # Write data in Bed format
        Bed('tempdir/toydata.bed',count_A1=False)


    iids and sids:

        Individuals are identified with an iid, which is a ndarray of two strings: a family ID and a case ID.
        SNP locations are identified with an sid string in a ndarray. For example:

        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)
        >>> print(snp_on_disk.iid[:3]) # print the first three iids
        [['POP1' '0']
         ['POP1' '12']
         ['POP1' '44']]
        >>> print(snp_on_disk.sid[:9]) # print the first nine sids
        ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4']
        >>> print(snp_on_disk.iid_to_index([['POP1','44'],['POP1','12']])) #Find the indexes for two iids.
        [2 1]
        
    Selecting and Reordering Individuals and SNPs

        You often don't want to read the SNP values for all iids and sids. You can use indexing to create a subsetting SnpReader that
        will read only the SNP values of interest.

        SnpReaders support the indexing formats supported by ndarray plus two generalizations. Here are examples of indexing with an array
        of indexes, with slicing, and with an array of Booleans.

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify some data on disk in Bed format
            >>> subset_snpreader_1 = snp_on_disk[[3,-1],:] #index with an array of indexes. Negatives count from end.
            >>> print((subset_snpreader_1.iid_count, subset_snpreader_1.sid_count))
            (2, 1015)
            >>> snpdata1 = subset_snpreader_1.read() # read just the two rows of interest from the disk
            >>> subset_snpreader_2 = snp_on_disk[:,:0:-2] #index with a slice
            >>> print((subset_snpreader_2.iid_count, subset_snpreader_2.sid_count))
            (300, 507)
            >>> boolindexes = [s.startswith('23_') for s in snp_on_disk.sid] # create a Boolean index of sids that start '23_'
            >>> subset_snpreader_3 = snp_on_disk[:,boolindexes] #index with array of Booleans
            >>> print((subset_snpreader_3.iid_count, subset_snpreader_3.sid_count))
            (300, 24)

        The first generalization over what ndarray offers is full indexing on both the iid dimension and the sid dimension, in other words,
        full multidimensional indexing. For example,

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify some data on disk in Bed format
            >>> subset_snpreader_4 = snp_on_disk[[3,4],:0:-2] # index on two dimensions at once
            >>> print((subset_snpreader_4.iid_count, subset_snpreader_4.sid_count))
            (2, 507)

        The second generalization is indexing on a single integer index.

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify some data on disk in Bed format
            >>> subset_snpreader_5 = snp_on_disk[5,:] #index with single integer
            >>> print((subset_snpreader_5.iid_count, subset_snpreader_5.sid_count))
            (1, 1015)

        Indexing is also useful when you have SNP values in memory via a :class:`SnpData` index and want to copy a subset of those values.
        While you could instead index directly on the `.SnpData.val` ndarray, by indexing on the :class:`SnpData` instance you
        also get iid and cid information.

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify some data on disk in Bed format
            >>> snpdata1 = snp_on_disk.read() # read all SNP values into memory
            >>> print(snpdata1.sid[:9]) # print the first 9 sids
            ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4']
            >>> snpdata_subset = snpdata1[:,::2].read(view_ok=True,order='A') # create a copy or view with every other sid
            >>> print(snpdata_subset.sid[:9])# print the first 9 sids in the subset
            ['1_12' '1_10' '1_28' '1_36' '1_4' '1_11' '1_32' '1_9' '1_17']


        You can apply indexing on top of indexing to specify subsets of subsets of data to read. In this example, 
        only the SNP values for every 16th sid is actually read from the disk.

            >>> # These are just SnpReaders, nothing is read from disk yet
            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify some data on disk in Bed format
            >>> half_snpreader = snp_on_disk[:,::2] # a reader for half the sids
            >>> quarter_snpreader = half_snpreader[:,::2] # a reader for half of half the sids
            >>> sixteenth_snpreader = quarter_snpreader[:,::2][:,::2] # a reader for half of half of half of half the sids
            >>> print(sixteenth_snpreader) #Print the specification of this reader
            Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)[:,::2][:,::2][:,::2][:,::2]
            >>> # Now we read from disk. Only values for one sid in every 16 will be read.
            >>> snpdata_sixteenth = sixteenth_snpreader.read()
            >>> print(snpdata_sixteenth.val[0,3])
            2.0
        
    When Data is Read:

        SNP data can be enormous so we generally avoid reading it to the degree practical. Specifically,
        
        * Constructing and printing a SnpReader causes no file reading. For example, these commands read no data:

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Construct a Bed SnpReader. No data is read.
            >>> print(snp_on_disk) # Print the Bed SnpReader specification. No data is read.
            Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)
            >>> subset_on_disk = snp_on_disk[[3,4],::2] # Construct a subsetting SnpReader. No data is read.
            >>> print(subset_on_disk) # print the subset SnpReader. No data is read.
            Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)[[3,4],::2]

        * Properties and methods related to the iids and sids (to the degree practical) read just some iid and sid data from the disk,
          not SNP value data. Moreover, the iid and sid data is read from file only once. Consider these commands:

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Construct a Bed SnpReader. No data is read.
            >>> print(snp_on_disk.sid[:9]) # without reading any SNP values data from disk, read the sid and iid data from disk, cache it, and then print the first nine sids.
            ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4']
            >>> print(snp_on_disk.sid_to_index(['1_10','1_13'])) #use the cached sid information to find the indexes of '1_10' and '1_13'. (No data is read from disk.)
            [2 9]

        * The only methods that read SNP values from file are :meth:`read` and :meth:`read_kernel` (to the degree practical). For example:

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Construct a Bed SnpReader. No data is read.
            >>> snpdata1 = snp_on_disk.read() #read all the SNP values from disk, creating a new SnpData instance that keeps these values in memory
            >>> print(snpdata1.val[0,2])# print the SNP value for the iid with index 0 and the sid with index 2. (No data is read from disk.)
            1.0

        * If you request the values for only a subset of the iids or sids, (to the degree practical) only that subset will be read from disk.
          for example:

            >>> subset_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)[[3,4],::2] # Construct a subsetting SnpReader. No data is read.
            >>> snpdata_subset = subset_on_disk.read() # from disk, read the SNP values for the iids with index 3 and 4 AND sids with even numbered indexes.
            >>> print(snpdata_subset.val[0,2]) # print the SNP value with subset iid index 0 and sid index 2 (corresponding to iid index 3 and sid index 4 in the full data). No data is read from disk.
            2.0

    When Data is Re-Read and Copied:

        Every time you call a SnpReader's :meth:`read` method, the SnpReader re-reads the SNP value data and returns a new in-memory :class:`.SnpData`
        (with :attr:`.SnpData.val` property containing a new ndarray of the SNP values). Likewise, when you call the :meth:`kernel` method, the SnpReader re-reads
        the data and returns a new kernel ndarray.

        Here is an example of what not to do, because it causes all the SNP value data to be read twice.

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Construct a Bed SnpReader. No data is read.
            >>> # The following is not recommended because it inefficiently reads all the SNP values twice.
            >>> print(snp_on_disk.read().val[0,2]) # read all values into a new SnpData, print a SNP value
            1.0
            >>> print(snp_on_disk.read().val[0,3]) # read all values (again) into a second new SnpData, print a SNP value
            2.0

        Here are two efficient alternatives. First, if all SNP values can all fit in memory, read them once into a :class:`SnpData` and then
        access that :class:`SnpData` multiple times.

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Construct a Bed SnpReader. No data is read.
            >>> snpdata1 = snp_on_disk.read() # read all values into a new SnpData
            >>> print(snpdata1.val[0,2]) # print a SNP value from snpdata1's in-memory ndarray
            1.0
            >>> print(snpdata1.val[0,3]) # print another SNP value from snpdata1's in-memory ndarray.
            2.0

        Second, if the SNP value data is too large to fit in memory, use subsetting to read only the SNP values of interest from disk.
       
            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Construct a Bed SnpReader. No data is read.
            >>> print(snp_on_disk[0,2].read().val[0,0]) #Define the subset of data and read only that subset from disk.
            1.0
            >>> print(snp_on_disk[0,3].read().val[0,0]) #Define a second subset of data and read only that subset from disk.
            2.0

        Because the in-memory :class:`.SnpData` class is a kind of SnpReader, you may read from it, too.
        Doing so create a new :class:`.SnpData` instance containing a copy of the SNP values in a new ndarray.

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Construct a Bed SnpReader. No data is read.
            >>> snpdata1 = snp_on_disk.read() # read all SNP values from disk into a new SnpData
            >>> print(snpdata1.val is snpdata1.val) # Do the in-memory SNP values use the same memory as themselves? Yes
            True
            >>> snpdata2 = snpdata1.read() # copy all the SNP values into a new ndarray in a new SnpData
            >>> print(snpdata2.val is snpdata1.val) # Do the two ndarrays of in-memory SNP values use the same memory?
            False


    Avoiding Unwanted ndarray Allocations

        You may want a subset of SNPs values from an in-memory :class:`SnpData` and you may know that this subset and the original :class:`SnpData`
        can safely share the memory of the ndarray of SNP values. For this case, the :meth:`read` has optional parameters called view_ok and order. If you override 
        the defaults of "view_ok=False,order='F'" with "view_ok=True,order='A', the :meth:`read` will, if practical, return a new 
        :class:`SnpData` with a ndarray shares memory with the original ndarray.
        Use these parameters with care because any change to either ndarray (for example, via :meth:`.SnpData.standardize`) will effect
        the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
        share memory and so it may ignore your suggestion and allocate a new ndarray anyway.

            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Construct a Bed SnpReader. No data is read.
            >>> snpdata1 = snp_on_disk.read() # read all data from disk into a SnpData with a new ndarray
            >>> column01 = snpdata1[:,0:1].read(view_ok=True,order='A') #create SnpData with the data from just the first two SNPs. Sharing memory is OK. The memory may be laid out in any order (that is sid-major and iid-major are both OK).
            >>> import numpy as np
            >>> #print np.may_share_memory(snpdata1.val, column01.val) # Do the two ndarray's share memory? They could (but currently they won't)
            >>> column201 = snpdata1[:,[2,0,1]].read(view_ok=True,order='A') #create SnpData with the data from three SNPs, permuted. Sharing memory is OK.
            >>> print(np.may_share_memory(snpdata1.val, column201.val)) # Do the two ndarray's share memory? No, ndarray decided that this indexing was too complex for sharing.
            False

    The :meth:`read` Method
  
        By default the :meth:`read` returns a ndarray of scipy.float64 laid out in memory in F-contiguous order (iid-index varies the fastest). You may, instead,
        ask for scipy.float32 or for C-contiguous order or any order. See :meth:`read` for details.

    The :meth:`.SnpData.standardize` Method
        The :meth:`.SnpData.standardize` method, available only on :class:`.SnpData`, does in-place standardization of the in-memory
        SNP data. By default, it applies 'Unit' standardization, that is: the values for each SNP will have mean zero and standard deviation 1.0.
        NaN values are then filled with zeros (consequently, if there are NaN values, the final standard deviation will not be zero).
        Note that, for efficiently, this method works in-place, actually changing values in the ndarray. Although it works in place, for convenience
        it also returns itself. See :meth:`.SnpData.standardize` for options and details.

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
   
    The :meth:`read_kernel` Method

        The :meth:`read_kernel` method, available on any SnpReader, returns a :class:`KernelData`. The :attr:`val` property of the :class:`KernelData` is
        an ndarray of the (possibility standardized) SNP values multiplied with their transposed selves. When applied to an read-from-disk SnpReader, such as :class:`.Bed`,
        the method can save memory by reading (and standardizing) the data in blocks. See :meth:`read_kernel` for details.


            >>> from pysnptools.standardizer import Unit
            >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify some data on disk in Bed format
            >>> kerneldata1 = snp_on_disk.read_kernel(Unit()) #Create an in-memory kernel from the snp data on disk.
            >>> print('{0:.6f}'.format(kerneldata1.val[0,0]))
            901.421836
            >>> kerneldata2 = snp_on_disk.read_kernel(Unit(),block_size=10) #Create an in-memory kernel from the snp data on disk, but only read 10 SNPS at a time from the disk.
            >>> print('{0:.6f}'.format(kerneldata2.val[0,0]))
            901.421836


    Details of Methods & Properties:
    """

    def __init__(self, *args, **kwargs):
        super(SnpReader, self).__init__(*args, **kwargs)

    @property
    def iid(self):
        """A ndarray of the iids. Each iid is a ndarray of two strings (a family ID and a case ID) that identifies an individual.

        :rtype: ndarray of strings with shape [:attr:`.iid_count`,2]

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)
        >>> print(snp_on_disk.iid[:3]) # print the first three iids
        [['POP1' '0']
         ['POP1' '12']
         ['POP1' '44']]
        """
        return self.row

    @property
    def iid_count(self):
        """number of iids

        :rtype: integer

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.
        """
        return self.row_count

    @property
    def sid(self):
        """A ndarray of the sids. Each sid is a string that identifies a SNP.

        :rtype: ndarray (length :attr:`.sid_count`) of strings

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)
        >>> print(snp_on_disk.sid[:9]) # print the first nine sids
        ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4']

        """
        return self.col

    @property
    def sid_count(self):
        """number of sids

        :rtype: integer

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        """
        return self.col_count

    #!!document that chr must not be X,Y,M only numbers (as per the PLINK BED format)
    #!!Also what about telling the ref and alt allele? Also, what about tri and quad alleles, etc?
    @property
    def pos(self):
        """A ndarray of the position information for each sid. Each element is a ndarray of three scipy.numbers's (chromosome, genetic distance, basepair distance).

        :rtype: ndarray of float64 with shape [:attr:`.sid_count`, 3]

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)
        >>> print(snp_on_disk.pos[:3,]) # print position information for the first three sids: #The '...' is for possible space char
        [[...1.          0.00800801  0.        ]
         [...1.          0.023023    1.        ]
         [...1.          0.0700701   4.        ]]
        """
        return self.col_property

    @property
    def row_property(self):
        """Defined as a zero-width array for compatibility with :class:`PstReader`, but not used.
        """
        if not hasattr(self,'_row_property'):
            self._row_property = np.empty((self.row_count,0))
        return self._row_property


    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok):
        raise NotImplementedError
    
    #!!check that views always return contiguous memory by default
    def read(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False):
        """Reads the SNP values and returns a :class:`.SnpData` (with :attr:`.SnpData.val` property containing a new ndarray of the SNP values).

        :param order: {'F' (default), 'C', 'A'}, optional -- Specify the order of the ndarray. If order is 'F' (default),
            then the array will be in F-contiguous order (iid-index varies the fastest).
            If order is 'C', then the returned array will be in C-contiguous order (sid-index varies the fastest).
            If order is 'A', then the :attr:`.SnpData.val`
            ndarray may be in any order (either C-, Fortran-contiguous).
        :type order: string or None

        :param dtype: {scipy.float64 (default), scipy.float32}, optional -- The data-type for the :attr:`.SnpData.val` ndarray.
        :type dtype: data-type

        :param force_python_only: optional -- If False (default), may use outside library code. If True, requests that the read
            be done without outside library code.
        :type force_python_only: bool

        :param view_ok: optional -- If False (default), allocates new memory for the :attr:`.SnpData.val`'s ndarray. If True,
            if practical and reading from a :class:`SnpData`, will return a new 
            :class:`SnpData` with a ndarray shares memory with the original :class:`SnpData`.
            Typically, you'll also wish to use "order='A'" to increase the chance that sharing will be possible.
            Use these parameters with care because any change to either ndarray (for example, via :meth:`.SnpData.standardize`) will effect
            the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
            share memory and so it may ignore your suggestion and allocate a new ndarray anyway.
        :type view_ok: bool

        :rtype: :class:`.SnpData`

        Calling the method again causes the SNP values to be re-read and creates a new in-memory :class:`.SnpData` with a new ndarray of SNP values.

        If you request the values for only a subset of the sids or iids, (to the degree practical) only that subset will be read from disk.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300.bed',count_A1=False) # Specify SNP data on disk
        >>> snpdata1 = snp_on_disk.read() # Read all the SNP data returning a SnpData instance
        >>> print(type(snpdata1.val).__name__) # The SnpData instance contains a ndarray of the data.
        ndarray
        >>> subset_snpdata = snp_on_disk[:,::2].read() # From the disk, read SNP values for every other sid
        >>> print(subset_snpdata.val[0,0]) # Print the first SNP value in the subset
        2.0
        >>> subsub_snpdata = subset_snpdata[:10,:].read(order='A',view_ok=True) # Create an in-memory subset of the subset with SNP values for the first ten iids. Share memory if practical.
        >>> import numpy as np
        >>> # print np.may_share_memory(subset_snpdata.val, subsub_snpdata.val) # Do the two ndarray's share memory? They could. Currently they won't.       
        """
        val = self._read(None, None, order, dtype, force_python_only, view_ok)
        from pysnptools.snpreader import SnpData
        ret = SnpData(self.iid,self.sid,val,pos=self.pos,name=str(self))
        return ret

    def iid_to_index(self, list):
        """Takes a list of iids and returns a list of index numbers

        :param list: list of iids
        :type order: list of list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify SNP data on disk
        >>> print(snp_on_disk.iid_to_index([['POP1','44'],['POP1','12']])) #Find the indexes for two iids.
        [2 1]
        """
        return self.row_to_index(list)

    def sid_to_index(self, list):
        """Takes a list of sids and returns a list of index numbers

        :param list: list of sids
        :type list: list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify SNP data on disk
        >>> print(snp_on_disk.sid_to_index(['1_10','1_13'])) #Find the indexes for two sids.
        [2 9]
        """
        return self.col_to_index(list)

    def __getitem__(self, iid_indexer_and_snp_indexer):
        import os
        from pysnptools.snpreader._subset import _SnpSubset
        iid_indexer, snp_indexer = iid_indexer_and_snp_indexer
        return _SnpSubset(self, iid_indexer, snp_indexer)

    def read_kernel(self, standardizer=None, block_size=None, order='A', dtype=np.float64, force_python_only=False, view_ok=False):
        """Returns a :class:`KernelData` such that the :meth:`KernelData.val` property will be a ndarray of the standardized SNP values multiplied with their transposed selves.

        :param standardizer: -- (required) Specify standardization to be applied before the matrix multiply. Any :class:`.Standardizer` may be used. Some choices include :class:`Standardizer.Identity` 
            (do nothing), :class:`.Unit` (make values for each SNP have mean zero and standard deviation 1.0) and :class:`Beta`.
        :type standardizer: :class:`.Standardizer`

        :param block_size: optional -- Default of None (meaning to load all). Suggested number of sids to read into memory at a time.
        :type block_size: int or None

        :rtype: class:`KernelData`

        Calling the method again causes the SNP values to be re-read and allocates a new class:`KernelData`.

        When applied to an read-from-disk SnpReader, such as :class:`.Bed`, the method can save memory by reading (and standardizing) the data in blocks.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> from pysnptools.standardizer import Unit
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify SNP data on disk
        >>> kerneldata1 = snp_on_disk.read_kernel(Unit())
        >>> print((int(kerneldata1.iid_count), '{0:.6f}'.format(kerneldata1.val[0,0])))
        (300, '901.421836')
        """
        assert standardizer is not None, "'standardizer' must be provided"

        from pysnptools.kernelreader import SnpKernel
        snpkernel = SnpKernel(self,standardizer=standardizer,block_size=block_size)
        kerneldata = snpkernel.read(order, dtype, force_python_only, view_ok)
        return kerneldata

    def kernel(self, standardizer, allowlowrank=False, block_size=10000, blocksize=None):
        """ .. Warning:: Deprecated. Use :meth:`read_kernel` instead.

        Returns a ndarray of size iid_count x iid_count. The returned array has the value of the standardized SNP values multiplied with their transposed selves.

        :param standardizer: -- Specify standardization to be applied before the matrix multiply. Any :class:`.Standardizer` may be used. Some choices include :class:`Standardizer.Identity` 
            (do nothing), :class:`.Unit` (make values for each SNP have mean zero and standard deviation 1.0), :class:`Beta`.
        :type standardizer: :class:`.Standardizer`

        :param block_size: optional -- Default of 10000. None means to load all. Suggested number of sids to read into memory at a time.
        :type block_size: int or None

        :rtype: ndarray of size :attr:`.iid_count` x :attr:`.iid_count`

        Calling the method again causes the SNP values to be re-read and allocates a new ndarray.

        When applied to an read-from-disk SnpReader, such as :class:`.Bed`, the method can save memory by reading (and standardizing) the data in blocks.

        If you request the values for only a subset of the sids or iids, (to the degree practical) only that subset will be read from disk.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> from pysnptools.standardizer import Unit
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False) # Specify SNP data on disk
        >>> kernel = snp_on_disk.kernel(Unit())
        >>> print(((int(kernel.shape[0]),int(kernel.shape[1])), '{0:.6f}'.format(kernel[0,0])))
        ((300, 300), '901.421836')
        """        #print "entering kernel with {0},{1},{2}".format(self, standardizer, blocksize)
        warnings.warn(".kernel(...) is deprecated. Use '.read_kernel(...).val", DeprecationWarning)
        if blocksize is not None:
            warnings.warn(".kernel(...blocksize...) is deprecated. Use '.kernel(...block_size=...)", DeprecationWarning)
            block_size = blocksize
        return self._read_kernel(standardizer, block_size=block_size)

    @staticmethod
    def _as_snpdata(snpreader, standardizer, force_python_only, dtype):
        '''
        Like 'read' except (1) won't read if already a snpdata and (2) returns the standardizer
        '''
        from pysnptools.snpreader import SnpData
        if isinstance(snpreader,SnpData) and snpreader.val.dtype==dtype and isinstance(standardizer,stdizer.Identity):
            return snpreader, stdizer.Identity()
        else:
            return snpreader.read(order='A',dtype=dtype).standardize(standardizer,return_trained=True,force_python_only=force_python_only)
    
    def _read_kernel(self, standardizer, block_size=None, order='A', dtype=np.float64, force_python_only=False, view_ok=False, return_trained=False):
        #Do all-at-once (not in blocks) if 1. No block size is given or 2. The #ofSNPs < Min(block_size,iid_count)
        if block_size is None or (self.sid_count <= block_size or self.sid_count <= self.iid_count):
            train_data,trained_standardizer  = SnpReader._as_snpdata(self,standardizer=standardizer,dtype=dtype,force_python_only=force_python_only)
            kernel = train_data._read_kernel(stdizer.Identity(), order=order,dtype=dtype,force_python_only=force_python_only,view_ok=False)
            if return_trained:
                return kernel, trained_standardizer
            else:
                return kernel

        else: #Do in blocks
            #Set the default order to 'C' because with kernels any order is fine and the Python .dot method likes 'C' best.
            if order=='A':
                order = 'C'
            t0 = time.time()
            K = np.zeros([self.iid_count,self.iid_count],dtype=dtype,order=order)
            trained_standardizer_list = []

            logging.info("reading {0} SNPs in blocks of {1} and adding up kernels (for {2} individuals)".format(self.sid_count, block_size, self.iid_count))

            ct = 0
            ts = time.time()

            for start in range(0, self.sid_count, block_size):
                ct += block_size
                train_data,trained_standardizer = SnpReader._as_snpdata(self[:,start:start+block_size],standardizer=standardizer,dtype=dtype,force_python_only=force_python_only)
                trained_standardizer_list.append(trained_standardizer)
                K += train_data._read_kernel(stdizer.Identity(),block_size=None,order=order,dtype=dtype,force_python_only=force_python_only,view_ok=False)
                if ct % block_size==0:
                    diff = time.time()-ts
                    if diff > 1: logging.info("read %s SNPs in %.2f seconds" % (ct, diff))

            t1 = time.time()
            logging.info("%.2f seconds elapsed" % (t1-t0))

            if return_trained:
                return K, standardizer._merge_trained(trained_standardizer_list) #turns this into a single standardizer
            else:
                return K

    def copyinputs(self, copier):
        raise NotImplementedError

    def _assert_iid_sid_pos(self):
        assert self._row.dtype.type is np.str_ and len(self._row.shape)==2 and self._row.shape[1]==2, "iid should be dtype str, have two dimensions, and the second dimension should be size 2"
        assert self._col.dtype.type is np.str_ and len(self._col.shape)==1, "sid should be of dtype of str and one dimensional"

    @staticmethod
    def _name_of_other_file(filename,remove_suffix,add_suffix):
        if filename.lower().endswith(remove_suffix.lower()):
            filename = filename[0:-1-len(remove_suffix)]
        return filename+"."+add_suffix

    @staticmethod
    def _write_fam(snpdata, basefilename, remove_suffix):
        famfile = SnpReader._name_of_other_file(basefilename, remove_suffix, "fam")

        with open(famfile,"w") as fam_filepointer:
            for iid_row in snpdata.iid:
                fam_filepointer.write("{0} {1} 0 0 0 0\n".format(iid_row[0],iid_row[1]))


    @staticmethod
    def _write_map_or_bim(snpdata, basefilename, remove_suffix, add_suffix):
        mapfile = SnpReader._name_of_other_file(basefilename, remove_suffix, add_suffix)

        with open(mapfile,"w") as map_filepointer:
            for sid_index, sid in enumerate(snpdata.sid):
                posrow = snpdata.pos[sid_index]
                map_filepointer.write("%r\t%s\t%r\t%r\tA\tC\n" % (posrow[0], sid, posrow[1], posrow[2]))


    @staticmethod
    def _read_fam(basefilename, remove_suffix):
        famfile = SnpReader._name_of_other_file(basefilename, remove_suffix, "fam")

        logging.info("Loading fam file {0}".format(famfile))
        if os.path.getsize(famfile)>0:
            iid = np.loadtxt(famfile, dtype = 'str',usecols=(0,1),comments=None)
        else:
            iid = np.empty((0,2), dtype = 'str')
        if len(iid.shape) == 1: #When empty or just one item, make sure the result is (x,2)
            iid = iid.reshape((len(iid)//2,2))
        return iid


    @staticmethod
    def _read_map_or_bim( basefilename, remove_suffix, add_suffix):
        mapfile = SnpReader._name_of_other_file(basefilename, remove_suffix, add_suffix)

        logging.info("Loading {0} file {1}".format(add_suffix, mapfile))
        if os.path.getsize(mapfile) == 0: #If the map/bim file is empty, return empty arrays
            sid = np.array([],dtype='str')
            pos = np.array([[]],dtype=int).reshape(0,3)
            return sid,pos
        else:
            fields = pd.read_csv(mapfile,delimiter = '\t',usecols = (0,1,2,3),header=None,index_col=False,comment=None)
            sid = np.array(fields[1].tolist(),dtype='str')
            pos = fields[[0,2,3]].values
            return sid,pos


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc t
    print("done")