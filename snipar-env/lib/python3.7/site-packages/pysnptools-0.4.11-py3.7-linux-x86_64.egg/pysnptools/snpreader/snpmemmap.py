from __future__ import print_function

import logging
import os
import numpy as np
import unittest
import doctest
import pysnptools.util as pstutil
from pysnptools.pstreader import PstData
from pysnptools.pstreader import PstMemMap
from pysnptools.snpreader import SnpReader, SnpData


class SnpMemMap(PstMemMap,SnpData):
    '''
    A :class:`.SnpData` that keeps its data in a memory-mapped file. This allows data large than fits in main memory.

    See :class:`.SnpData` for general examples of using SnpData.

    **Constructor:**
        :Parameters: **filename** (*string*) -- The *\*.snp.memmap* file to read.
        
        Also see :meth:`.SnpMemMap.empty` and :meth:`.SnpMemMap.write`.

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import SnpMemMap
        >>> snp_mem_map = SnpMemMap('../examples/tiny.snp.memmap')
        >>> print(snp_mem_map.val[0,1], snp_mem_map.iid_count, snp_mem_map.sid_count)
        2.0 2 3

    **Methods inherited from** :class:`.SnpData`

        :meth:`.SnpData.allclose`, :meth:`.SnpData.standardize`

    **Methods beyond** :class:`.SnpReader`

    '''

    def __init__(self, *args, **kwargs):
        super(SnpMemMap, self).__init__(*args, **kwargs)

    val = property(PstMemMap._get_val,PstMemMap._set_val)
    """The 2D NumPy memmap array of floats that represents the values.

    >>> from pysnptools.snpreader import SnpMemMap
    >>> snp_mem_map = SnpMemMap('../examples/tiny.snp.memmap')
    >>> print(snp_mem_map.val[0,1])
    2.0
    """

    @property
    def offset(self):
        '''The byte position in the file where the memory-mapped values start.
       
        (The disk space before this is used to store :attr:`.SnpReader.iid`, etc. information.
        This property is useful when interfacing with, for example, external Fortran and C matrix libraries.)
        
        '''
        self._run_once()
        return self._offset

    @property
    def filename(self):
        '''The name of the memory-mapped file
        '''
        #Don't need '_run_once'
        return self._filename

    @staticmethod
    def empty(iid, sid, filename, pos=None,order="F",dtype=np.float64):
        '''Create an empty :class:`.SnpMemMap` on disk.

        :param iid: The :attr:`.SnpReader.iid` information
        :type iid: an array of string pairs

        :param sid: The :attr:`.SnpReader.sid` information
        :type sid: an array of strings

        :param filename: name of memory-mapped file to create
        :type filename: string

        :param pos: optional -- The additional :attr:`.SnpReader.pos` information associated with each sid. Default: None
        :type pos: an array of numeric triples

        :param order: {'F' (default), 'C'}, optional -- Specify the order of the ndarray.
        :type order: string or None

        :param dtype: {scipy.float64 (default), scipy.float32}, optional -- The data-type for the :attr:`.SnpMemMap.val` ndarray.
        :type dtype: data-type

        :rtype: :class:`.SnpMemMap`

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.snpreader import SnpMemMap
        >>> filename = "tempdir/tiny.snp.memmap"
        >>> pstutil.create_directory_if_necessary(filename)
        >>> snp_mem_map = SnpMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename,order="F",dtype=np.float64)
        >>> snp_mem_map.val[:,:] = [[0.,2.,0.],[0.,1.,2.]]
        >>> snp_mem_map.flush()

        '''

        self = SnpMemMap(filename)
        self._empty_inner(row=iid, col=sid, filename=filename, row_property=None, col_property=pos,order=order,dtype=dtype)
        return self

    def flush(self):
        '''Flush :attr:`.SnpMemMap.val` to disk and close the file. (If values or properties are accessed again, the file will be reopened.)

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.snpreader import SnpMemMap
        >>> filename = "tempdir/tiny.snp.memmap"
        >>> pstutil.create_directory_if_necessary(filename)
        >>> snp_mem_map = SnpMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename,order="F",dtype=np.float64)
        >>> snp_mem_map.val[:,:] = [[0.,2.,0.],[0.,1.,2.]]
        >>> snp_mem_map.flush()

        '''
        if self._ran_once:
            self.val.flush()
            del self._val
            self._ran_once = False


    @staticmethod
    def write(filename, snpdata):
        """Writes a :class:`SnpData` to :class:`SnpMemMap` format.

        :param filename: the name of the file to create
        :type filename: string
        :param snpdata: The in-memory data that should be written to disk.
        :type snpdata: :class:`SnpData`
        :rtype: :class:`.SnpMemMap`

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.snpreader import SnpData, SnpMemMap
        >>> data1 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],val= [[0.,2.,0.],[0.,1.,2.]])
        >>> pstutil.create_directory_if_necessary("tempdir/tiny.pst.memmap")
        >>> SnpMemMap.write("tempdir/tiny.snp.memmap",data1)      # Write data1 in SnpMemMap format
        SnpMemMap('tempdir/tiny.snp.memmap')
        """

        #We write iid and sid in ascii for compatibility between Python 2 and Python 3 formats.
        row_ascii = np.array(snpdata.row,dtype='S') #!!!avoid this copy when not needed
        col_ascii = np.array(snpdata.col,dtype='S') #!!!avoid this copy when not needed
        self = PstMemMap.empty(row_ascii, col_ascii, filename, row_property=snpdata.row_property, col_property=snpdata.col_property,order=PstMemMap._order(snpdata),dtype=snpdata.val.dtype)
        self.val[:,:] = snpdata.val
        self.flush()
        logging.debug("Done writing " + filename)
        return SnpMemMap(filename)



    def _run_once(self):
            if (self._ran_once):
                return
            row_ascii,col_ascii,val,row_property,col_property = self._run_once_inner()
            row = np.array(row_ascii,dtype='str') #!!!avoid this copy when not needed
            col = np.array(col_ascii,dtype='str') #!!!avoid this copy when not needed

            SnpData.__init__(self,iid=row,sid=col,val=val,pos=col_property,name="np.memmap('{0}')".format(self._filename))

class TestSnpMemMap(unittest.TestCase):     

    def test1(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        filename2 = "tempdir/tiny.snp.memmap"
        pstutil.create_directory_if_necessary(filename2)
        snpreader2 = SnpMemMap.empty(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],filename=filename2,order="F",dtype=np.float64)
        assert isinstance(snpreader2.val,np.memmap)
        snpreader2.val[:,:] = [[0.,2.,0.],[0.,1.,2.]]
        assert np.array_equal(snpreader2[[1],[1]].read(view_ok=True).val,np.array([[1.]]))
        snpreader2.flush()
        assert isinstance(snpreader2.val,np.memmap)
        assert np.array_equal(snpreader2[[1],[1]].read(view_ok=True).val,np.array([[1.]]))
        snpreader2.flush()

        snpreader3 = SnpMemMap(filename2)
        assert np.array_equal(snpreader3[[1],[1]].read(view_ok=True).val,np.array([[1.]]))
        assert isinstance(snpreader3.val,np.memmap)

        logging.info("in TestSnpMemMap test1")
        snpreader = SnpMemMap('../examples/tiny.snp.memmap')
        assert snpreader.iid_count == 2
        assert snpreader.sid_count == 3
        assert isinstance(snpreader.val,np.memmap)

        snpdata = snpreader.read(view_ok=True)
        assert isinstance(snpdata.val,np.memmap)
        os.chdir(old_dir)


def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpMemMap))
    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    r.run(suites)

    result = doctest.testmod()
    assert result.failed == 0, "failed doc test: " + __file__
