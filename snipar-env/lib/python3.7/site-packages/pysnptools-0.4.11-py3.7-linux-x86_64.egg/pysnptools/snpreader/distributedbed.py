from __future__ import absolute_import
import os
import shutil
import numpy as np
import logging
import unittest
from pysnptools.snpreader import _MergeSIDs
from pysnptools.snpreader import SnpReader, Bed
from pysnptools.pstreader import PstReader
from pysnptools.util import _multiopen
from pysnptools.snpreader import _snps_fixup
from pysnptools.util import log_in_place
from pysnptools.util.mapreduce1 import map_reduce
from pysnptools.util.filecache import FileCache
import six
from six.moves import range

class DistributedBed(SnpReader):
    '''
    A class that implements the :class:`SnpReader` interface. It stores :class:`.Bed`-like data in pieces on storage. When you request data, it retrieves only the needed pieces.

    **Constructor:**
        :Parameters: **storage** (string or :class:`.FileCache`) -- Tells where the DistirubtedBed data is stored.
                      A string can be given and will be interpreted as the path to a directory.
                      A :class:`.FileCache` instance can be given, which provides a method to specify cluster-distributed storage.
        :type storage: string or :class:`.FileCache`

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import DistributedBed
        >>> data_on_disk = DistributedBed('../examples/toydataSkip10.distributedbed')
        >>> print((data_on_disk.iid_count, data_on_disk.sid_count))
        (500, 1000)

    '''
    def __init__(self, storage):
        super(DistributedBed, self).__init__()

        self._ran_once = False
        self._storage = FileCache._fixup(storage)
        self._merge = None


    def __repr__(self): 
        return "{0}({1})".format(self.__class__.__name__,self._storage)

    @property
    def row(self):
        """*same as* :attr:`iid`
        """
        self._run_once()
        return self._merge.row

    @property
    def col(self):
        """*same as* :attr:`sid`
        """
        self._run_once()
        return self._merge.col

    @property
    def col_property(self):
        """*same as* :attr:`pos`
        """
        self._run_once()
        return self._merge.col_property

    def _run_once(self):
        if self._ran_once:
            return
        self._ran_once = True

        _metadatanpz = "metadata.npz"
        with self._storage.open_read(_metadatanpz) as handle_metadatanpz_file_name:
            #self._file_dict["metadatanpz"] = handle_metadatanpz
            _reader_name_listnpz = "reader_name_list.npz"
            with self._storage.open_read(_reader_name_listnpz) as handle_reader_name_listnpz_file_name:
                reader_name_list = np.array(np.load(handle_reader_name_listnpz_file_name)['reader_name_list'],dtype='str')
                #self._file_dict["reader_name_listnpz"] = handle_reader_name_listnpz

                reader_list = [_Distributed1Bed(reader_name,self._storage) for reader_name in reader_name_list]

                self._merge = _MergeSIDs(reader_list,cache_file=handle_metadatanpz_file_name,skip_check=True)

                for reader in reader_list:
                    reader._row = self._merge.row
            

    #def __del__(self):
    #    for handle in self._file_dict.itervalues():
    #        handle.close()

    def copyinputs(self, copier):
        pass

    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok):
        self._run_once()
        return self._merge._read(iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok)

    #!!! in future could make a default for piece_per_chrom_count that made each piece some GB size
    @staticmethod
    def write(storage, snpreader, piece_per_chrom_count=1, updater=None, runner=None): #!!! might want to set pieces_per_chrom such that it is a certain size
        '''
        Uploads from any :class:`.Bed`-like data to cluster storage for efficient retrieval later.
        If some of the contents already exists in storage, it skips uploading that part of the contents. (To avoid this behavior,
        clear the storage.)

        :param storage: Tells where to store SNP data.
                      A string can be given and will be interpreted as the path of a local directory to use for storage. (The local
                      directory will **not** be automatically erased and so must be user managed.) 
                      A :class:`.FileCache` instance can be given, which provides a
                      method to specify cluster-distributed storage. (:class:`.FileCache`'s will **not** be automatically erased and must be user managed.)
                      If `None`, the storage will be in an automatically-erasing temporary directory. (If the TEMP environment variable is set, Python places the temp directory under it.)
                      
        :type storage: string or :class:`.FileCache` or None.

        :param snpreader: A :class:`.Bed` or other :class:`.SnpReader` with values of 0,1,2, or missing.
            (Note that this differs from most other `write` methods that take a :class:`.SnpData`)
        :type snpreader: :class:`.SnpReader`

        :param piece_per_chrom_count: The number of pieces in which to store the data from each chromosome. Data is split across
            SNPs. For exmple, if `piece_per_chrom_count` is set to 100 and 22 chromosomes are uploaded, then data will be stored in 2200 pieces. Later, when data is requested
            only the pieces necessary for the request will be downloaded to local storage.
        :type piece_per_chrom_count: A number

        :param updater: A single argument function to write logging message to, for example, the function created by :func:`.log_in_place`.
        :type updater: A function or lambda

        :param runner: a :class:`.Runner`, optional: Tells how to run.
            (Note that :class:`.Local` and :class:`.LocalMultProc` are good options.)
            If not given, the function is run locally.
        :type runner: :class:`.Runner`

        :rtype: DistributedBed

        >>> from pysnptools.snpreader import DistributedBed, Bed
        >>> import shutil
        >>> directory = 'tempdir/toydataSkip10.distributedbed'
        >>> if os.path.exists(directory):
        ...     shutil.rmtree(directory)
        >>> snpreader = Bed('../examples/toydata.bed',count_A1=False)[:,::10]  # Read every 10 snps from Bed format
        >>> DistributedBed.write(directory,snpreader,piece_per_chrom_count=5)  # Write data in DistributedBed format
        DistributedBed(LocalCache('tempdir/toydataSkip10.distributedbed'))


        '''
        from pysnptools.util import _file_transfer_reporter
        from pysnptools.util.filecache import FileCache

        count_A1 = True #Make all these's the same for reading and writing so that nothing will change.
        snpreader = _snps_fixup(snpreader, count_A1=count_A1)

        storage = FileCache._fixup(storage)

        chrom_set = sorted(set(snpreader.pos[:,0]))
        for chrom in chrom_set:
            assert chrom==chrom and chrom==int(chrom), "DistributedBed.write expects all chromosomes to be integers (not '{0}')".format(chrom)
        with _file_transfer_reporter("DistributedBed.write", size=0, updater=updater) as updater2:
            def mapper_closure(chrom):
                chrom_reader = snpreader[:,snpreader.pos[:,0]==chrom]
                def nested_closure(piece_per_chrom_index):
                    start = chrom_reader.sid_count * piece_per_chrom_index // piece_per_chrom_count
                    stop = chrom_reader.sid_count * (piece_per_chrom_index+1) // piece_per_chrom_count
                    piece_reader = chrom_reader[:,start:stop]
                    _piece_name_list = ["chrom{0}.piece{1}of{2}.{3}".format(int(chrom),piece_per_chrom_index,piece_per_chrom_count,suffix) for suffix in ['bim','fam','bed']]
                    exist_list = [storage.file_exists(_piece_name) for _piece_name in _piece_name_list]
                    if sum(exist_list) < 3: #If all three of the BIM/FAM/BED files are already there, then skip the upload, otherwise do the upload
                        for i in range(3): #If one or two of BIM/FAM/BED are there, remove them
                            if exist_list[i]:
                                storage.remove(_piece_name_list[i])
                        _Distributed1Bed.write(_piece_name_list[-1],storage,piece_reader.read(),count_A1=count_A1,updater=updater2)
                    return _piece_name_list[-1]
                return map_reduce(range(piece_per_chrom_count),
                    mapper=nested_closure,
                    )
            list_list_pair = map_reduce(chrom_set,
                nested = mapper_closure,
                runner=runner,
                )                

        reader_name_list = []
        reader_list = []
        for chrom_result in list_list_pair:
            for _piece_name in chrom_result:
                reader_name_list.append(_piece_name)
                reader_list.append(_Distributed1Bed(_piece_name,storage))
                

        _metadatanpz = "metadata.npz"
        with storage.open_write(_metadatanpz) as local_metadatanpz:
            _reader_name_listnpz = "reader_name_list.npz"
            with storage.open_write(_reader_name_listnpz) as local_reader_name_listnpz:
                reader_name_list_ascii = np.array(reader_name_list,dtype='S')
                np.savez(local_reader_name_listnpz,reader_name_list=reader_name_list_ascii)
                if os.path.exists(local_metadatanpz):
                    os.remove(local_metadatanpz)
                _MergeSIDs(reader_list,cache_file=local_metadatanpz,skip_check=True)

        return DistributedBed(storage)

class _Distributed1Bed(SnpReader):
    '''
    An atomic set of bed/bim/fam files stored somewhere. Can answer metadata questions without downloading the *.bed file.
    But does download the whole *.bed file when any SNP value is requested.
    '''
    def __init__(self,path,storage):
        super(_Distributed1Bed, self).__init__()

        self._ran_once = False
        self._file_dict = {}

        self._storage = storage
        self.path = path
        self.local = None

    def __repr__(self): 
        return "{0}('{1}','{2}')".format(self.__class__.__name__,self.path,self._storage)


    @property
    def row(self):
        """*same as* :attr:`iid`
        """
        if not hasattr(self,"_row"):
            _fam = SnpReader._name_of_other_file(self.path,remove_suffix="bed", add_suffix="fam")
            local_fam = self._storage.open_read(_fam)
            self._row = SnpReader._read_fam(local_fam.__enter__(),remove_suffix="fam")
            self._file_dict["fam"] = local_fam
        return self._row

    @property
    def col(self):
        """*same as* :attr:`sid`
        """
        if not hasattr(self,"_col"):
            _bim = SnpReader._name_of_other_file(self.path,remove_suffix="bed", add_suffix="bim")
            local_bim = self._storage.open_read(_bim)
            self._col, self._col_property = SnpReader._read_map_or_bim(local_bim.__enter__(),remove_suffix="bim", add_suffix="bim")
            self._file_dict["bim"] = local_bim
        return self._col

    @property
    def col_property(self):
        """*same as* :attr:`pos`
        """
        if not hasattr(self,"_col"):
            self.col #get col info
        return self._col_property

    def _run_once(self):
        if self._ran_once:
            return
        self._ran_once = True
        self.row # get row info
        self.col # get col info

        _bed = SnpReader._name_of_other_file(self.path,remove_suffix="bed", add_suffix="bed")
        local_bed = self._storage.open_read(_bed)
        self.local = Bed(local_bed.__enter__(),count_A1=True,iid=self.row,sid=self.col,pos=self.col_property,skip_format_check=True)
        self._file_dict["bed"] = local_bed

    def __del__(self):
        for handle in six.itervalues(self._file_dict):
            handle.__exit__(None,None,None)
        self._file_dict = {}

    def copyinputs(self, copier):
        pass

    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok):
        self._run_once()
        return self.local._read(iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok)
    
    @staticmethod
    def write(path, storage, snpdata, count_A1=True, updater=None):
        file_list = [SnpReader._name_of_other_file(path,remove_suffix="bed", add_suffix=new_suffix) for new_suffix in ["bim","fam","bed"]] #'bed' should be last
        with _multiopen(lambda file_name:storage.open_write(file_name,updater=updater),file_list) as local_file_name_list:
            Bed.write(local_file_name_list[-1],snpdata,count_A1=count_A1)

        return _Distributed1Bed(path,storage)


class TestDistributedBed(unittest.TestCase):     

    def test1(self):
        logging.info("in TestDistributedBed test1")
        from pysnptools.snpreader import SnpGen, DistributedBed
        snpgen = SnpGen(seed=0,iid_count=100,sid_count=100)

        temp_dir = 'tempdir/distributed_bed_test1'
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        distributed_bed = DistributedBed.write(temp_dir,snpgen,piece_per_chrom_count=2)
        snpdata = distributed_bed.read()

        ref1 = DistributedBed(os.path.dirname(os.path.realpath(__file__))+'/../../tests/datasets/distributed_bed_test1').read()
        assert(snpdata.allclose(ref1,equal_nan=True))

        ref2 = Bed(os.path.dirname(os.path.realpath(__file__))+'/../../tests/datasets/distributed_bed_test1_X',count_A1=False).read()
        assert(snpdata.allclose(ref2,equal_nan=True))


def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDistributedBed))
    return test_suite


if __name__ == "__main__":
    import doctest
    logging.basicConfig(level=logging.INFO)

    if False:
        from pysnptools.snpreader import DistributedBed, Bed
        import shutil
        directory = 'tempdir/toydataSkip10.distributedbed'
        if os.path.exists(directory):
            shutil.rmtree(directory)
        snpreader = Bed('../examples/toydata.bed',count_A1=False)[:,::10]  # Read every 10 snps from Bed format
        DistributedBed.write(directory,snpreader,piece_per_chrom_count=5)  # Write data in DistributedBed format

    result = doctest.testmod()
    assert result.failed == 0, "failed doc test: " + __file__


    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    r.run(suites)


