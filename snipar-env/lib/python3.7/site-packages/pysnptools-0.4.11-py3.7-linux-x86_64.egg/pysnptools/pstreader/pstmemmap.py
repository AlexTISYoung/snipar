from __future__ import absolute_import
from __future__ import print_function
import logging
import numpy as np
import unittest
import os
import doctest
import pysnptools.util as pstutil
from pysnptools.pstreader import PstReader, PstData

class PstMemMap(PstData):
    '''
    A :class:`.PstData` that keeps its data in a memory-mapped file. This allows data large than fits in main memory.

    See :class:`.PstData` for general examples of using PstData.

    **Constructor:**
        :Parameters: **filename** (*string*) -- The *\*.pst.memmap* file to read.
        
        Also see :meth:`.PstMemMap.empty` and :meth:`.PstMemMap.write`.

        :Example:

        >>> from pysnptools.pstreader import PstMemMap
        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> pst_mem_map = PstMemMap('../examples/tiny.pst.memmap')
        >>> print(pst_mem_map.val[0,1], pst_mem_map.row_count, pst_mem_map.col_count)
        2.0 3 2

    **Methods beyond** :class:`PstReader`

    '''


    def __init__(self, filename):
        PstReader.__init__(self)
        self._ran_once = False
        self._filename = filename

    def __repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self._filename)

    def _get_val(self):
        self._run_once()
        return self._val

    def _set_val(self, new_value):
        self._run_once()
        if self._val is new_value:
            return
        raise Exception("PstMemMap val's cannot be set to a different array")

    val = property(_get_val,_set_val)
    """The 2D NumPy memmap array of floats that represents the values.

    >>> from pysnptools.pstreader import PstMemMap
    >>> pst_mem_map = PstMemMap('../examples/tiny.pst.memmap')
    >>> print(pst_mem_map.val[0,1])
    2.0
    """


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

    @property
    def offset(self):
        '''The byte position in the file where the memory-mapped values start.
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
    def empty(row, col, filename, row_property=None, col_property=None,order="F",dtype=np.float64):
        '''Create an empty :class:`.PstMemMap` on disk.

        :param row: The :attr:`.PstReader.row` information
        :type row: an array of anything

        :param col: The :attr:`.PstReader.col` information
        :type col: an array of anything

        :param filename: name of memory-mapped file to create
        :type filename: string

        :param row_property: optional -- The additional :attr:`.PstReader.row_property` information associated with each row. Default: None
        :type row_property: an array of anything

        :param col_property: optional -- The additional :attr:`.PstReader.col_property` information associated with each col. Default: None
        :type col_property: an array of anything

        :param order: {'F' (default), 'C'}, optional -- Specify the order of the ndarray.
        :type order: string or None

        :param dtype: {scipy.float64 (default), scipy.float32}, optional -- The data-type for the :attr:`.PstMemMap.val` ndarray.
        :type dtype: data-type

        :rtype: :class:`.PstMemMap`

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.pstreader import PstMemMap
        >>> filename = "tempdir/tiny.pst.memmap"
        >>> pstutil.create_directory_if_necessary(filename)
        >>> pst_mem_map = PstMemMap.empty(row=['a','b','c'],col=['y','z'],filename=filename,row_property=['A','B','C'],order="F",dtype=np.float64)
        >>> pst_mem_map.val[:,:] = [[1,2],[3,4],[np.nan,6]]
        >>> pst_mem_map.flush()

        '''

        self = PstMemMap(filename)
        self._empty_inner(row, col, filename, row_property, col_property,order,dtype)
        return self


    def _empty_inner(self, row, col, filename, row_property, col_property,order,dtype):
        self._ran_once = True
        self._dtype = dtype
        self._order = order

        row = PstData._fixup_input(row)
        col = PstData._fixup_input(col)
        row_property = PstData._fixup_input(row_property,count=len(row))
        col_property = PstData._fixup_input(col_property,count=len(col))

        with open(filename,'wb') as fp:
            np.save(fp, row)
            np.save(fp, col)
            np.save(fp, row_property)
            np.save(fp, col_property)
            np.save(fp, np.array([self._dtype]))
            np.save(fp, np.array([self._order]))
            self._offset = fp.tell()

        logging.info("About to start allocating memmap '{0}'".format(filename))
        val = np.memmap(filename, offset=self._offset, dtype=dtype, mode="r+", order=order, shape=(len(row),len(col)))
        logging.info("Finished allocating memmap '{0}'. Size is {1}".format(filename,os.path.getsize(filename)))
        PstData.__init__(self,row,col,val,row_property,col_property,name="np.memmap('{0}')".format(filename))

    def _run_once(self):
        if (self._ran_once):
            return
        row,col,val,row_property,col_property = self._run_once_inner()
        PstData.__init__(self,row,col,val,row_property,col_property,name="np.memmap('{0}')".format(self._filename))


    def _run_once_inner(self):
        self._ran_once = True

        logging.debug("np.load('{0}')".format(self._filename))
        with open(self._filename,'rb') as fp:
            row = np.load(fp,allow_pickle=True)
            col = np.load(fp,allow_pickle=True)
            row_property = np.load(fp,allow_pickle=True)
            col_property = np.load(fp,allow_pickle=True)
            self._dtype = np.load(fp,allow_pickle=True)[0]
            self._order = np.load(fp,allow_pickle=True)[0]
            self._offset = fp.tell()
        val = np.memmap(self._filename, offset=self._offset, dtype=self._dtype, mode='r', order=self._order, shape=(len(row),len(col)))
        return row,col,val,row_property,col_property

        

    def copyinputs(self, copier):
        # doesn't need to self._run_once()
        copier.input(self._filename)

    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = True
    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        val, shares_memory = self._apply_sparray_or_slice_to_val(self.val, row_index_or_none, col_index_or_none, order, dtype, force_python_only)
        if not shares_memory and view_ok:
            logging.warn("Read from {0} required copy".format(self)) #!!!keep this warning?
        if shares_memory and not view_ok:
            val = val.copy(order='K')
        return val

    @staticmethod
    def _order(pstdata):
        if pstdata.val.flags['F_CONTIGUOUS']:
            return "F"
        if pstdata.val.flags['C_CONTIGUOUS']:
            return "C"
        raise Exception("Don't know order of PstData's value")


    def flush(self):
        '''Flush :attr:`.PstMemMap.val` to disk and close the file. (If values or properties are accessed again, the file will be reopened.)

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.pstreader import PstMemMap
        >>> filename = "tempdir/tiny.pst.memmap"
        >>> pstutil.create_directory_if_necessary(filename)
        >>> pst_mem_map = PstMemMap.empty(row=['a','b','c'],col=['y','z'],filename=filename,row_property=['A','B','C'],order="F",dtype=np.float64)
        >>> pst_mem_map.val[:,:] = [[1,2],[3,4],[np.nan,6]]
        >>> pst_mem_map.flush()

        '''
        if self._ran_once:
            self.val.flush()
            del self._val
            self._ran_once = False

    @staticmethod
    def write(filename, pstdata):
        """Writes a :class:`PstData` to :class:`PstMemMap` format and returns the :class:`.PstMemMap`.

        :param filename: the name of the file to create
        :type filename: string
        :param pstdata: The in-memory data that should be written to disk.
        :type pstdata: :class:`PstData`
        :rtype: :class:`.PstMemMap`

        >>> import pysnptools.util as pstutil
        >>> from pysnptools.pstreader import PstData, PstMemMap
        >>> data1 = PstData(row=['a','b','c'],col=['y','z'],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        >>> pstutil.create_directory_if_necessary("tempdir/tiny.pst.memmap")
        >>> PstMemMap.write("tempdir/tiny.pst.memmap",data1)      # Write data1 in PstMemMap format
        PstMemMap('tempdir/tiny.pst.memmap')
        """

        self = PstMemMap.empty(pstdata.row, pstdata.col, filename, row_property=pstdata.row_property, col_property=pstdata.col_property,order=PstMemMap._order(pstdata),dtype=pstdata.val.dtype)
        self.val[:,:] = pstdata.val
        self.flush()
        logging.debug("Done writing " + filename)

        return self


class TestPstMemMap(unittest.TestCase):     

    def test1(self):
        logging.info("in TestPstMemMap test1")
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        filename2 = "tempdir/tiny.pst.memmap"
        pstutil.create_directory_if_necessary(filename2)
        pstreader2 = PstMemMap.empty(row=['a','b','c'],col=['y','z'],filename=filename2,row_property=['A','B','C'],order="F",dtype=np.float64)
        assert isinstance(pstreader2.val,np.memmap)
        pstreader2.val[:,:] = [[1,2],[3,4],[np.nan,6]]
        assert np.array_equal(pstreader2[[0],[0]].read(view_ok=True).val,np.array([[1.]]))
        pstreader2.flush()
        assert isinstance(pstreader2.val,np.memmap)
        assert np.array_equal(pstreader2[[0],[0]].read(view_ok=True).val,np.array([[1.]]))
        pstreader2.flush()

        pstreader3 = PstMemMap(filename2)
        assert np.array_equal(pstreader3[[0],[0]].read(view_ok=True).val,np.array([[1.]]))
        assert isinstance(pstreader3.val,np.memmap)

        pstreader = PstMemMap('../examples/tiny.pst.memmap')
        assert pstreader.row_count == 3
        assert pstreader.col_count == 2
        assert isinstance(pstreader.val,np.memmap)

        pstdata = pstreader.read(view_ok=True)
        assert isinstance(pstdata.val,np.memmap)
        os.chdir(old_dir)

    def test2(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        filename = "tempdir/x.pst.memmap"
        pstutil.create_directory_if_necessary(filename)

        a = PstMemMap.empty(row=['a','b','c'],col=['y','z'],filename=filename,row_property=['A','B','C'],order="F",dtype=np.float64)
        b = PstData(row=['a','b','c'],col=['y','z'],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        pointer1, read_only_flag = a.val.__array_interface__['data']
        a.val+=1
        a.val+=b.val
        pointer2, read_only_flag = a.val.__array_interface__['data']
        assert pointer1==pointer2
        os.chdir(old_dir)

    def test3(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        filename = "tempdir/x.pst.memmap"
        pstutil.create_directory_if_necessary(filename)

        a = PstMemMap.empty(row=['a','b','c'],col=['y','z'],filename=filename,row_property=['A','B','C'],order="F",dtype=np.float64)
        pstdata = a.read(order='C',view_ok=True)
        os.chdir(old_dir)





def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPstMemMap))
    return test_suite



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    if True:
        from pysnptools.pstreader import PstMemMap
        fn = '../examples/tiny.pst.memmap'
        os.getcwd()
        print((os.path.exists(fn)))
        pst_mem_map = PstMemMap(fn)
        print((pst_mem_map.val[0,1]))


    if False:
        a=np.ndarray([2,3])
        pointer, read_only_flag = a.__array_interface__['data']
        print(pointer)
        a*=2
        pointer, read_only_flag = a.__array_interface__['data']
        print(pointer)
        a = PstMemMap.empty(row=['a','b','c'],col=['y','z'],filename=r'c:\deldir\a.memmap',row_property=['A','B','C'],order="F",dtype=np.float64)
        b = PstData(row=['a','b','c'],col=['y','z'],val=[[1,2],[3,4],[np.nan,6]],row_property=['A','B','C'])
        pointer, read_only_flag = a.val.__array_interface__['data']
        print(pointer)
        a.val+=1
        a.val+=b.val
        pointer, read_only_flag = a.val.__array_interface__['data']
        print(pointer)


    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    r.run(suites)

    result = doctest.testmod()
    assert result.failed == 0, "failed doc test: " + __file__
