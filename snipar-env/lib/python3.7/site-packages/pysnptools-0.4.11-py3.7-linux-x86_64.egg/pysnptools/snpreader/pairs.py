import os
import numpy as np
import logging
from pysnptools.snpreader import SnpReader

class _Pairs(SnpReader):
    '''
    Experimental.
    '''

    #!!!see fastlmm\association\epistasis.py for code that allows ranges of snps to be specified when making pairs
    def __init__(self, snpreader0,snpreader1=None, do_standardize=True,sid_materialize_limit=1000*1000,_include_single_times_single=False): 
        #!!! could add option to change snp separator and another to encode chrom, etc in the snp name
        super(Pairs, self).__init__()
        self._ran_once = False
        self.snpreader0 = snpreader0
        self.snpreader1 = snpreader1
        self.do_standardize = do_standardize
        self.sid_materialize_limit = sid_materialize_limit
        self._include_single_times_single=_include_single_times_single

    def __repr__(self):
        part2 = '' if self.snpreader1 is None else ',{0}'.format(self.snpreader1)
        return "{0}({1}{2})".format(self.__class__.__name__,self.snpreader0,part2)

    @property
    def row(self):
        """*same as* :attr:`iid`
        """
        if not hasattr(self,"_row"):
            self._row = self.snpreader0.row
            assert self.snpreader1 is None or np.array_equal(self._row,self.snpreader1.row), "Expect snpreaders to have the same iids in the same order"
        return self._row

    @property
    def col(self):
        """*same as* :attr:`sid`
        """
        if not hasattr(self,"_col"):
            assert self.col_count < self.sid_materialize_limit, '{:,} is too many sids to materialize'.format(self.col_count)
            if self.snpreader1 is None:
                col_0 = self.snpreader0.col
                col_list = []
                self.index0_list = [] #!!!should be _index0_list, etc
                self.index1_list = []
                for index0 in xrange(self.snpreader0.col_count):
                    #logging.info("index0={0} of {1}".format(index0,self.snpreader0.col_count))#!!!
                    start1 = index0 if self._include_single_times_single else index0+1
                    for index1 in xrange(start1,self.snpreader0.col_count):
                        col_list.append('{0},{1}'.format(col_0[index0],col_0[index1]))
                        self.index0_list.append(index0)
                        self.index1_list.append(index1)
                self._col = np.array(col_list)
                self.index0_list = np.array(self.index0_list)
                self.index1_list = np.array(self.index1_list)
            else:
                col_0 = self.snpreader0.col
                col_1 = self.snpreader1.col
                assert len(set(col_0)&set(col_1))==0,"Pairs currently requires two snpreaders to have no sids in common"
                col_list = []
                self.index0_list = []
                self.index1_list = []
                for index0 in xrange(self.snpreader0.col_count):
                    #logging.info("index0={0} of {1}".format(index0,self.snpreader0.col_count))#!!!
                    for index1 in xrange(self.snpreader1.col_count):
                        col_list.append('{0},{1}'.format(col_0[index0],col_1[index1]))
                        self.index0_list.append(index0)
                        self.index1_list.append(index1)
                self._col = np.array(col_list)
                self.index0_list = np.array(self.index0_list)
                self.index1_list = np.array(self.index1_list)

        assert self.col_count == len(self._col), "real assert"
        return self._col

    @property
    def col_count(self):
        n0 = self.snpreader0.col_count
        if self.snpreader1 is None:
            if self._include_single_times_single:
                return (n0*n0+n0)//2
            else:
                return (n0*n0-n0)//2
        else:
            return n0*self.snpreader1.col_count

    @property
    def col_property(self):
        """*same as* :attr:`pos`
        """
        if not hasattr(self,"_col_property"):
            self._col_property = np.zeros([self.sid_count,3],dtype=np.int64)
        return self._col_property

    def copyinputs(self, copier):
        # doesn't need to self._run_once() because only uses original inputs !!!is this true?
        self.snpreader0.copyinputs(copier)
        if self.snpreader1 is not None:
            self.snpreader1.copyinputs(copier)

    def _run_once(self):
        if self._ran_once:
            return

        self._ran_once = True
        self.col

    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok):
        self._run_once()

        iid_count_in = self.iid_count
        sid_count_in = self.sid_count

        if iid_index_or_none is not None:
            iid_count_out = len(iid_index_or_none)
            iid_index_out = iid_index_or_none
        else:
            iid_count_out = iid_count_in
            iid_index_out = range(iid_count_in)

        if sid_index_or_none is not None:
            sid_count_out = len(sid_index_or_none)
            sid_index_out = sid_index_or_none
        else:
            sid_count_out = sid_count_in
            sid_index_out = range(sid_count_in)

        
        sid_index_inner_0 = self.index0_list[sid_index_out] #Find the sid_index of the left snps of interest
        sid_index_inner_1 = self.index1_list[sid_index_out] #Find the sid_index of the right snps of interest
        if self.snpreader1 is None:
            sid_index_inner_01 = np.unique(np.r_[sid_index_inner_0,sid_index_inner_1]) #Index of every snp of interest
            inner_01 = self.snpreader0[iid_index_or_none,sid_index_inner_01].read(order=order, dtype=dtype, force_python_only=force_python_only, view_ok=True) #read every val of interest
            val_inner_01 = inner_01.standardize().val if self.do_standardize else inner_01.val

            sid_index_inner_01_reverse = {v:i for i,v in enumerate(sid_index_inner_01)} #Dictionary of snp_index to position in sid_index_inner_01
            sid_index_inner_0_in_val = np.array([sid_index_inner_01_reverse[i] for i in sid_index_inner_0])  #Replace snp_index0 with column # in val_inner_01
            sid_index_inner_1_in_val = np.array([sid_index_inner_01_reverse[i] for i in sid_index_inner_1])  #Replace snp_index1 with column # in val_inner_01
            val_inner_0 = val_inner_01[:,sid_index_inner_0_in_val] #Extract the vals for the left snps of interest
            val_inner_1 = val_inner_01[:,sid_index_inner_1_in_val]#Extract the vals for the right snps of interest
        else:
            inner_0 = self.snpreader0[iid_index_or_none,sid_index_inner_0].read(order=order, dtype=dtype, force_python_only=force_python_only, view_ok=True) #read every val of interest
            inner_1 = self.snpreader1[iid_index_or_none,sid_index_inner_1].read(order=order, dtype=dtype, force_python_only=force_python_only, view_ok=True) #read every val of interest
            val_inner_0 = inner_0.standardize().val if self.do_standardize else inner_0.val
            val_inner_1 = inner_1.standardize().val if self.do_standardize else inner_1.val
        val = val_inner_0*val_inner_1 #Element multiplication creates the vals for the pairs
        return val



def split_on_sids(snpreader,part_count):
    sid_count = snpreader.sid_count
    start = 0
    for part_index in xrange(1,part_count+1):
        end = part_index*sid_count//part_count
        yield snpreader[:,start:end]
        start=end

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    from pysnptools.snpreader import Pheno, Bed
    import pysnptools.util as pstutil

    data_file = 'd:\OneDrive\programs\epiCornell\syndata.bed'
    if False:
        from pysnptools.snpreader import SnpData
        import numpy as np
        bed1 = Bed("../../tests/datasets/synth/all")
        print(bed1.iid_count, bed1.sid_count, bed1.iid_count * bed1.sid_count)
        #goal 1500 individuals x 27000 SNP
        snpdata1 = bed1.read()
        iid = bed1.iid
        sid = ['sid{0}'.format(i) for i in xrange(27000)]
        val = np.tile(snpdata1.val,(3,6))[:,:27000].copy()
        #snpdata = Pheno('pysnptools/examples/toydata.phe').read()         # Read data from Pheno format
        snpdata2 = SnpData(iid, sid, val)
        print(snpdata2.iid_count, snpdata2.sid_count, snpdata2.iid_count * snpdata2.sid_count)
        Bed.write(snpdata2,data_file,count_A1=False)

    synbed = Bed(data_file)
    print(synbed.iid_count, synbed.sid_count, synbed.iid_count * synbed.sid_count)

    part_count = 1000
    part_list = list(split_on_sids(synbed,part_count))

    pairs00 = _Pairs(part_list[0])
    from fastlmm.association import single_snp
    pheno_fn = r"d:\OneDrive\programs\epiCornell\pheno.txt"
    cov_fn = r"d:\OneDrive\programs\epiCornell\cov.txt"
    results_df = single_snp(pairs00, K0=synbed, pheno=pheno_fn, covar=cov_fn, leave_out_one_chrom=False, count_A1=True)

    if False:
        for i,synbed_part_i in enumerate(synbed_part_list):
            for j,synbed_part_j in enumerate(synbed_part_list):
                if j<i:
                    continue #not break
                print("Looking at pair {0},{1}".format(i,j))
                pairs = _Pairs(synbed_part_i) if i==j else _Pairs(synbed_part_i,synbed_part_j)
                #print(pairs.iid)
                print('{:,}'.format(pairs.sid_count))
                #print(pairs.sid)
                #print(pairs.pos)
                #print(pairs.row_property)
                snpdata = pairs.read()#
                #print(snpdata.val)

    import datetime
    from pysnptools.kernelreader import SnpKernel
    from pysnptools.standardizer import Unit
    from pysnptools.util.mapreduce1.runner import LocalMultiProc
    from pysnptools.util.mapreduce1 import map_reduce
    #runner=None
    runner = LocalMultiProc(1,just_one_process=False)

    part_pair_count = (part_count*part_count+part_count)//2
    part_pair_index = -1
    print("part_pair_count={0:,}".format(part_pair_count))

    K0 = SnpKernel(synbed,standardizer=Unit()).read() #Precompute the similarity

    start_time = datetime.datetime.now()
    for i,part_i in enumerate(part_list):
        def mapper1(j):
            #from fastlmm.association import single_snp
            #from pysnptools.snpreader import Pairs
            #print('Z')
            #part_j = part_list[j]
            #print('A')
            print("Looking at pair {0},{1} which is {2} of {3}".format(i,j,part_pair_index+j+1,part_pair_count))
            #pairs = Pairs(part_i) if i==j else Pairs(part_i,part_j)
            #result_df_ij = single_snp(pairs, K0=K0, pheno=pheno_fn, covar=cov_fn, leave_out_one_chrom=False, count_A1=True)
            #print(result_df_ij[:1])
            #return result_df_ij

        result_df_i = map_reduce(range(i,part_count),
                                 mapper=mapper1,
                                 reducer=lambda result_j_list:pd.concat(result_j_list),
                                 runner=runner,
                                 name='js')
        part_pair_index+=(part_count-i)
        time_so_far = datetime.datetime.now()-start_time
        total_time_estimate = time_so_far*part_pair_count/(part_pair_index+1)
        print(total_time_estimate)

