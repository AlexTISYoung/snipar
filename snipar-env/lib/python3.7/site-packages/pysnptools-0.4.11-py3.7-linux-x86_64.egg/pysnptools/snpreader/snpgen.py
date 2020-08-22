from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import logging
import os
import unittest
import doctest
from pysnptools.snpreader import SnpReader
from six.moves import range

class SnpGen(SnpReader):
    '''
    A :class:`.SnpReader` that generates deterministic SNP data on the fly.

    See :class:`.SnpReader` for general examples of using SnpReader.

    **Constructor:**
        :Parameters: * **seed** (*number*) -- The random seed that (along with *sid_batch_size*) determines the SNP values.
        :Parameters: * **iid_count** (*number*) --  the number of iids (number of individuals)
        :Parameters: * **sid_count** (*number*) --  the number of sids (number of SNPs)
        :Parameters: * **chrom_count** (*number*) --  the number of chromosomes to generate (must be 22 or fewer)
        :Parameters: * **cache_file** (*string*) -- (default None) If provided, tells where to cache the common iid, sid, and pos information. Using it can save time.
        :Parameters: * **sid_batch_size** (*number*) -- (default 1000) Tells how many SNP to generate at once. The default value is usually good.
        
        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import SnpGen
        >>> #Prepare to generate data for 1000 individuals and 1,000,000 SNPs
        >>> snp_gen = SnpGen(seed=332,iid_count=1000,sid_count=1000*1000)
        >>> print(snp_gen.iid_count,snp_gen.sid_count)
        1000 1000000
        >>> snp_data = snp_gen[:,200*1000:201*1000].read() #Generate for all users and for SNPs 200K to 201K
        >>> print(snp_data.val[1,1], snp_data.iid_count, snp_data.sid_count)
        0.0 1000 1000


    :Also See: :func:`.snp_gen`

    '''

    def __init__(self, seed, iid_count, sid_count, chrom_count=22, cache_file=None, sid_batch_size=1000):
        self._ran_once = False
        self._cache_file = cache_file

        self._seed = seed
        self._iid_count = iid_count
        self._sid_count = sid_count
        self._chrom_count = chrom_count
        self._sid_batch_size = sid_batch_size


        if cache_file is not None:
            if not os.path.exists(cache_file):
                import pysnptools.util as pstutil
                pstutil.create_directory_if_necessary(cache_file)
                self._run_once()
                np.savez(cache_file, _row=self._row, _col=self._col, _col_property=self._col_property)
            else:
                with np.load(cache_file,allow_pickle=True) as data:
                    self._row = data['_row']
                    assert len(self._row) == iid_count, "The iid in the cache file has a different length than iid_count"
                    self._col = data['_col']
                    assert len(self._col) == sid_count, "The sid in the cache file has a different length than sid_count"
                    self._col_property = data['_col_property']
                    assert len(self._col_property) == sid_count, "The pos in the cache file has a different length than sid_count"
                    self._ran_once = True


    def __repr__(self): 
        return "{0}(seed={1},iid_count={2},sid_count={3},chrom_count={4},sid_batch_size={5},cache_file={6})".format(self.__class__.__name__,self._seed,self._iid_count,self._sid_count,self._chrom_count,self._sid_batch_size,self._cache_file)

    @property
    def row(self):
        self._run_once()
        return self._row

    @property
    def col(self):
        self._run_once()
        return self._col

    @property
    def col_property(self):
        self._run_once()
        return self._col_property

    _chrom_size = np.array([263,255,214,203,194,183,171,155,145,144,144,143,114,109,106,98,92,85,67,72,50,56],dtype=np.int64)*int(1e6) #The approximate size of human chromosomes in base pairs

    def _run_once(self):
        if (self._ran_once):
            return
        self._ran_once = True
        self._row = np.array([('0','iid_{0}'.format(i)) for i in range(self._iid_count)])
        self._col = np.array(['sid_{0}'.format(i) for i in range(self._sid_count)])
        self._col_property = np.zeros(((self._sid_count),3))

        chrom_total = SnpGen._chrom_size[:self._chrom_count].sum()
        step = chrom_total // self._sid_count
        start = 0
        chrom_size_so_far = 0
        for chrom_index in range(self._chrom_count):
            chrom_size_so_before = chrom_size_so_far
            chrom_size_so_far += SnpGen._chrom_size[chrom_index]
            stop = chrom_size_so_far * self._sid_count // chrom_total
            self._col_property[start:stop,0] = chrom_index+1
            self._col_property[start:stop,2] = np.arange(0,stop-start)*step+chrom_size_so_before+1
            #print(chrom_index+1,start,stop,self._sid_count)
            start = stop

    def copyinputs(self, copier):
        self._run_once()
        copier.input(self._cache_file)

    # Most _read's support only indexlists or None, but this one supports Slices, too.
    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        self._run_once()
        import pysnptools.util as pstutil

        row_index_count = len(row_index_or_none) if row_index_or_none is not None else self._iid_count # turn to a count of the index positions e.g. all of them
        col_index = col_index_or_none if col_index_or_none is not None else np.arange(self._sid_count) # turn to an array of index positions, e.g. 0,1,200,2200,10
        batch_index = col_index // self._sid_batch_size  #find the batch index of each index position, e.g. 0,0,0,2,0
        val = np.empty((row_index_count,len(col_index))) #allocate memory for result
        list_batch_index = list(set(batch_index))
        for i in list_batch_index:  #for each distinct batch index, generate snps
            logging.info("working on snpgen batch {0} of {1}".format(i,len(list_batch_index))) #!!!why does this produce messages like 'working on snpgen batch 8 of 2'?
            start = i*self._sid_batch_size  #e.g. 0 (then 2000)
            stop = start + self._sid_batch_size #e.g. 1000, then 3000
            batch_val = self._get_val2(start,stop) # generate whole batch
            a = (batch_index==i) #e.g. [True,True,True,False,True], then [False,False,False,True,False]
            b = col_index[a]-start #e.g.  0,1,200,10, then 200
            val[:,a] = batch_val[:,b] if row_index_or_none is None else pstutil.sub_matrix(batch_val, row_index_or_none, b)

        return val

    def _get_dist(self):
        w = np.array([-0.6482249 , -8.49790398])
    
        def dist_fit(x):
            log_y = w[0] * np.log(x) + w[1]
            return np.exp(log_y)
    
        x_sample = np.logspace(np.log10(.1/self._iid_count),np.log10(.5),100,base=10) #discretize from 1e-7 (assuming iid_count=1e6) to .5, logarithmically
        y_sample = np.array([dist_fit(x) for x in x_sample])                    #Find the relative weight of each point
        dist = y_sample/y_sample.sum()
        return x_sample, dist



    @staticmethod
    def _get_sid(sid_start, sid_stop):
        sid = ["sid_{0}".format(i) for i in range(sid_start,sid_stop)]
        return sid

    def _get_val2(self, sid_start, sid_stop):
        x_sample, dist = self._get_dist() # The discrete distribution of minor allele frequencies (based on curve fitting to real data)
        val = SnpGen._get_val(x_sample,dist,self._iid_count,sid_start,sid_stop, self._seed)
        return val

    @staticmethod
    def _get_val(x_sample,dist,iid_count,sid_start,sid_stop,seed):
        missing_rate = .218
        sid_batch_size_a = sid_stop-sid_start
        sid_batch_size_b = sid_stop-sid_start

        np.random.seed(seed+sid_start)
        val_a = np.empty((iid_count,sid_batch_size_a),order="F")
        sid_so_far = 0
    
        while sid_so_far < sid_batch_size_a:
            logging.debug(sid_so_far)
            sid_index_to_freq = np.random.choice(x_sample,size=sid_batch_size_b,replace=True,p=dist) #For each sid, pick its minor allele freq
            val_b = (np.random.rand(iid_count,2,sid_batch_size_b) < sid_index_to_freq).sum(axis=1).astype(float) #Sample each allele and then sum the minor alleles #!!!slowest part
            missing = np.random.rand(iid_count,sid_batch_size_b)<missing_rate
            val_b[missing] = np.nan

            sid_index_to_minor_allele_count = np.nansum(val_b,axis=0) #Find the # of minor alleles for each allele. By design this will be mostly 0's
    
            #Keep only the sids with a minor allele count between 1 and iid_count.
            legal = (sid_index_to_minor_allele_count != 0) * (sid_index_to_minor_allele_count <= iid_count)
            legal_count = legal.sum()
            val_a[:,sid_so_far:sid_so_far+legal_count]=val_b[:,legal][:,:sid_batch_size_a-sid_so_far] #!!! takes about 1/3 of the time. Faster if both were order="F"?
            sid_so_far += legal_count
    
        return val_a


class TestSnpGen(unittest.TestCase):     

    def test1(self):
        logging.info("in TestSnpGen test1")
        seed = 0
        snpgen = SnpGen(seed=seed,iid_count=1000,sid_count=1000*1000)
        snpdata = snpgen[:,[0,1,200,2200,10]].read()
        snpdata2 = snpgen[:,[0,1,200,2200,10]].read()
        assert(snpdata.allclose(snpdata2))

        from pysnptools.snpreader import Bed
        ref = Bed( os.path.dirname(os.path.realpath(__file__))+'/../../tests/datasets/snpgen.bed',count_A1=False).read()
        assert(snpdata.allclose(ref,equal_nan=True))

        cache_file = 'tempdir/cache_file_test1.npz'
        os.remove(cache_file) if os.path.exists(cache_file) else None
        snpgen3 = SnpGen(seed=seed,iid_count=1000,sid_count=1000*1000,cache_file=cache_file)
        snpdata3 = snpgen3[::10,[0,1,200,2200,10]].read()
        assert(snpdata3.allclose(snpdata2[::10,:].read()))
        snpgen4 = SnpGen(seed=seed,iid_count=1000,sid_count=1000*1000,cache_file=cache_file)
        snpdata4 = snpgen4[::10,[0,1,200,2200,10]].read()
        assert(snpdata4.allclose(snpdata2[::10,:].read()))

def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpGen))
    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    r.run(suites)

    result = doctest.testmod()
    assert result.failed == 0, "failed doc test: " + __file__

