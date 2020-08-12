from __future__ import absolute_import
from __future__ import print_function
import os
import logging
import numpy as np
import scipy as sp
from six.moves import range
import sys
from contextlib import contextmanager
import time
import datetime
from pysnptools.util.intrangeset import IntRangeSet

def _testtest(data, iididx):
    return (data[0][iididx],data[1][iididx])

def intersect_apply(data_list, sort_by_dataset=True, intersect_before_standardize=True, is_test=False): #!!!add doc and tests for is_test
    """Intersects and sorts the iids from a list of datasets, returning new version of the datasets with all the same iids in the same order.

    :param data_list: list of datasets
    :type data_list: list
    :param sort_by_dataset: optional, If True (default), the iids are ordered according to the first non-None dataset.
        If False, the order is arbitrary, but consistent.
    :type sort_by_dataset: bool

    :param intersect_before_standardize: optional. Special code for :class:`.SnpKernel`, the class that postpones computing a kernel from SNP data. 
        If True (default), :class:`.SnpKernel` will remove any iids before SNP standardization before computing the kernel.
        If False, SNPs will be standardized with all iids, then the kernel will be computed, then any iids will be removed.
    :type intersect_before_standardize: bool

    :rtype: list of datasets

    Here are the dataset formats understood and what is returned for each.

    ============================================== ================================================================
    Dataset Format                                 What is Returned
    ============================================== ================================================================
    None                                           None
    A :class:`.SnpReader`                          A new subsetting :class:`.SnpReader` with adjusted iid
    A :class:`.KernelReader`                       A new subsetting :class:`.KernelReader` with adjusted iid
    A dictionary with ['iid'] and ['vals'] keys    The same dictionary but with the iid and vals values adjusted
    Tuple of the form (val ndarray, iid list)      A new tuple with the val ndarray and iid list adjusted
    ============================================== ================================================================
    
    If the iids in all the datasets are already the same and in the same order, then the datasets are returned without change.

    Notice that only dictionaries are processed in-place. Inputting a :class:`.SnpReader` and :class:`.KernelReader` returns a new class of the same type (unless its iids
    are already ok). Inputting a tuple returns a new tuple (unless its iids are already ok).

    :Example:

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> from pysnptools.snpreader import Bed, Pheno
    >>> from pysnptools.kernelreader import SnpKernel
    >>> from pysnptools.standardizer import Unit
    >>> #Create five datasets in different formats
    >>> ignore_in = None
    >>> kernel_in = SnpKernel(Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False),Unit()) # Create a kernel from a Bed file
    >>> pheno_in = Pheno('../../tests/datasets/phenSynthFrom22.23.N300.randcidorder.txt',missing="")
    >>> cov = Pheno('../../tests/datasets/all_chr.maf0.001.covariates.N300.txt',missing="").read()
    >>> cov_as_tuple_in = (cov.val,cov.iid) #We could do cov directly, but as an example we make it a tuple.
    >>>
    >>> # Create five new datasets with consistent iids
    >>> ignore_out, kernel_out, pheno_out, cov_as_tuple_out = intersect_apply([ignore_in, kernel_in, pheno_in, cov_as_tuple_in])
    >>> # Print the first five iids from each dataset
    >>> print(ignore_out, kernel_out.iid[:5], pheno_out.iid[:5], cov_as_tuple_out[1][:5])
    None [['POP1' '0']
     ['POP1' '12']
     ['POP1' '44']
     ['POP1' '58']
     ['POP1' '65']] [['POP1' '0']
     ['POP1' '12']
     ['POP1' '44']
     ['POP1' '58']
     ['POP1' '65']] [['POP1' '0']
     ['POP1' '12']
     ['POP1' '44']
     ['POP1' '58']
     ['POP1' '65']]
    """

    iid_list = []
    reindex_list = []

    #LATER this doesn't cover non-square kernel readers. What should it do when there are two different iid lists?
    from pysnptools.kernelreader import SnpKernel
    from pysnptools.kernelreader import Identity as IdentityKernel
    for data in data_list:
        if data is None:
            iid = None
            reindex = lambda data, iididx : None
        elif intersect_before_standardize and isinstance(data,SnpKernel):
            iid = data.iid1 if is_test else data.iid0
            reindex = lambda data, iididx,is_test=is_test : _reindex_snpkernel(data, iididx, is_test)
        elif intersect_before_standardize and isinstance(data,IdentityKernel):
            iid = data.iid1 if is_test else data.iid0
            reindex = lambda data, iididx,is_test=is_test : _reindex_identitykernel(data, iididx, is_test)
        else:
            try: #pheno dictionary
                iid = data['iid'] 
                reindex = lambda data, iididx : _reindex_phen_dict(data, iididx)
            except:
                if hasattr(data,'iid1') and is_test: #test kernel
                    iid = data.iid1
                    reindex = lambda data, iididx : data[:,iididx]
                else:
                    try:
                        iid = data.iid
                        try:
                            if iid is data.col: #If the 'iid' shares memory with the 'col', then it's a square-kernel-like thing and should be processed with just once index
                                reindex = lambda data, iididx : data[iididx]
                            else:
                                reindex = lambda data, iididx : data[iididx,:]
                        except:
                            reindex = lambda data, iididx : data[iididx,:]
                    except AttributeError: #tuple of (val,iid)
                        iid = data[1]
                        reindex = lambda data, iididx : _testtest(data,iididx)

        iid_list.append(iid)
        reindex_list.append(reindex)

    if len(iid_list) == 0: raise Exception("Expect a least one input item")

    if _all_same(iid_list):
        logging.debug("iids match up across {0} data sets".format(len(iid_list)))
        return data_list
    else:
        logging.debug("iids do not match up, so intersecting the data over individuals")            
        indarr = intersect_ids(iid_list)
        assert indarr.shape[0] > 0, "no individuals remain after intersection, check that ids match in files"

        if sort_by_dataset:
            #Look for first non-None iid
            for i, iid in enumerate(iid_list):
                if iid is not None:
                    #sort the indexes so that SNPs ids in their original order (and 
                    #therefore we have to move things around in memory the least amount)
                    sortind=np.argsort(indarr[:,i])
                    indarr=indarr[sortind]
                    break

        #!!! for the case in which some of the data items don't need to change, can we avoid calling reindex? Alternatively, should the _read code notice that all the iids are the same and in the same order                    
        data_out_list = []
        for i in range(indarr.shape[1]):
            data = data_list[i]
            iididx = indarr[:,i]
            reindex = reindex_list[i]
            new_data = reindex(data, iididx)
            data_out_list.append(new_data)

        return data_out_list

def _reindex_phen_dict(phen_dict, iididx):
    if len(phen_dict['vals'].shape)==1:
        phen_dict['vals'] = phen_dict['vals'][iididx]
    else:
        phen_dict['vals'] = phen_dict['vals'][iididx,:]
    phen_dict['iid'] = phen_dict['iid'][iididx]
    return phen_dict

def _reindex_snpkernel(snpkernel, iididx, is_test=False):
    from pysnptools.kernelreader import SnpKernel
    if not is_test:
        new_reader = snpkernel.snpreader[iididx,:]
        reference = new_reader
        result = SnpKernel(new_reader,snpkernel.standardizer,block_size=snpkernel.block_size)
    else:
        new_reader = snpkernel.test[iididx,:]
        result = SnpKernel(snpkernel.snpreader,snpkernel.standardizer,block_size=snpkernel.block_size)
    return result

def _reindex_identitykernel(identitykernel, iididx, is_test=False):
    from pysnptools.kernelreader import Identity as IdentityKernel
    if not is_test:
        iid = identitykernel.iid[iididx]
        result = IdentityKernel(iid)
    else:
        result = IdentityKernel(identitykernel.iid0,test=identitykernel.iid1[iididx])
    return result

def _all_same(iids_list):
    for i in range(len(iids_list)-1):
        iidA = iids_list[i]
        iidB = iids_list[i+1]
        if not np.array_equal(iidA,iidB):
            return False
    return True

def intersect_ids(idslist):
    '''
    .. deprecated::
       Use :func:`.intersect_apply` instead.    
    
    Takes a list of 2d string arrays of family and case ids.
    These are intersected.

    :rtype: indarr, an array of size N x L, where N is the number of
        individuals in the intersection, and L is the number of lists in idslist, and which
        contains the index to use (in order) such that all people will be identical and in order
        across all data sets.

    If one of the lists=None, it is ignored (but still has values reported in indarr, all equal to -1),
    '''
    id2ind={}    
    L=len(idslist)
    observed=sp.zeros(L,dtype='bool')
    first = True
    for l, id_list in enumerate(idslist):
        if id_list is not None:
            observed[l]=1
            if first:
                first = False
                for i in range(id_list.shape[0]):
                    id=(id_list[i,0], id_list[i,1])
                    entry=sp.zeros(L)*sp.nan #id_list to contain the index for this id, for all lists provided
                    entry[l]=i                 #index for the first one
                    id2ind[id]=entry
            else:
                for i in range(id_list.shape[0]):
                    id=(id_list[i,0], id_list[i,1])
                    if id in id2ind:
                        id2ind[id][l]=i

    indarr=sp.array(list(id2ind.values()),dtype='float')  #need float because may contain NaNs
    indarr[:,~observed]=-1                          #replace all Nan's from empty lists to -1
    inan = sp.isnan(indarr).any(1)                  #find any rows that contain at least one Nan
    indarr=indarr[~inan]                            #keep only rows that are not NaN
    indarr=sp.array(indarr,dtype='int')             #convert to int so can slice 
    indarr = indarr[indarr[:,0].argsort()]          #make sure Py2 and Py3 return results in the same arbitrary order
    return indarr   


def sub_matrix(val, row_index_list, col_index_list, order='A', dtype=sp.float64):
    """
    Efficiently creates a sub-matrix from a 2-D ndarray.

    :param val: The ndarray from which to copy.
    :type val: ndarray
    :param row_index_list: Tells which rows (and in which order) to copy into the sub-matrix
    :type row_index_list: list of integers
    :param col_index_list: Tells which columns (and in which order) to copy into the sub-matrix
    :type col_index_list: list of integers
    :param order: {'A' (default), 'C', 'F'}, optional -- Specify the order of the sub-matrix to create.
        If order is 'C', then the returned sub-matrix will be in C-contiguous order (first index varies the fastest).
        If order is 'F', then the array will be in F-contiguous order (second index varies the fastest).
        If order is 'A', then sub-matrix may be in any order F or C.
    :type order: string or None
    :param dtype: {scipy.float64 (default), scipy.float32}, optional -- The data-type for sub-matrix created.
    :type dtype: data-type

    :rtype: ndarray

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> import numpy as np
    >>> import pysnptools.util as pstutil
    >>> np.random.seed(0) # set seed so that results are deterministic
    >>> matrix = np.random.rand(12,7) # create a 12 x 7 ndarray
    >>> submatrix = pstutil.sub_matrix(matrix,[0,2,11],[6,5,4,3,2,1,0])
    >>> print(int(submatrix.shape[0]),int(submatrix.shape[1]))
    3 7
    >>> print(matrix[2,0] == submatrix[1,6]) #The row # 2 is now #1, the column #0 is now #6.
    True

    Note: Behind the scenes, for performance, this function selects and then calls one of 16 C++ helper functions.
    """
    from pysnptools.snpreader import wrap_matrix_subset

    if order == 'A':
        if val.flags['F_CONTIGUOUS']:
            effective_order = 'F'
        else:
            effective_order = 'C'
    else:
        effective_order = order

    sub_val = sp.empty((len(row_index_list), len(col_index_list)),dtype=dtype,order=effective_order)

    logging.debug("About to call cython matrixSubset")

    iid_count, sid_count = val.shape

    if val.flags['F_CONTIGUOUS']:
        if val.dtype ==  sp.float64:
            if dtype == sp.float64:
                if effective_order=="F":
                    wrap_matrix_subset.matrixSubsetDoubleFToDoubleFAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                elif effective_order=="C":
                    wrap_matrix_subset.matrixSubsetDoubleFToDoubleCAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                else:
                    raise Exception("order '{0}' not known, only 'F' and 'C'".format(effective_order));
            elif dtype == sp.float32:
                if effective_order=="F":
                    wrap_matrix_subset.matrixSubsetDoubleFToSingleFAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                elif effective_order=="C":
                    wrap_matrix_subset.matrixSubsetDoubleFToSingleCAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                else:
                    raise Exception("dtype '{0}' not known, only float64 and float32".format(dtype))
        elif val.dtype ==  sp.float32:
            if dtype == sp.float64:
                if effective_order=="F":
                    wrap_matrix_subset.matrixSubsetSingleFToDoubleFAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                elif effective_order=="C":
                    wrap_matrix_subset.matrixSubsetSingleFToDoubleCAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                else:
                    raise Exception("order '{0}' not known, only 'F' and 'C'".format(effective_order));
            elif dtype == sp.float32:
                if effective_order=="F":
                    wrap_matrix_subset.matrixSubsetSingleFToSingleFAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                elif effective_order=="C":
                    wrap_matrix_subset.matrixSubsetSingleFToSingleCAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                else:
                    raise Exception("dtype '{0}' not known, only float64 and float32".format(dtype))
        else:
            raise Exception("input dtype '{0}' not known, only float64 and float32".format(val.dtype))
    elif val.flags['C_CONTIGUOUS']:
        if val.dtype ==  sp.float64:
            if dtype == sp.float64:
                if effective_order=="F":
                    wrap_matrix_subset.matrixSubsetDoubleCToDoubleFAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                elif effective_order=="C":
                    wrap_matrix_subset.matrixSubsetDoubleCToDoubleCAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                else:
                    raise Exception("order '{0}' not known, only 'F' and 'C'".format(effective_order));
            elif dtype == sp.float32:
                if effective_order=="F":
                    wrap_matrix_subset.matrixSubsetDoubleCToSingleFAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                elif effective_order=="C":
                    wrap_matrix_subset.matrixSubsetDoubleCToSingleCAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                else:
                    raise Exception("dtype '{0}' not known, only float64 and float32".format(dtype))
        elif val.dtype ==  sp.float32:
            if dtype == sp.float64:
                if effective_order=="F":
                    wrap_matrix_subset.matrixSubsetSingleCToDoubleFAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                elif effective_order=="C":
                    wrap_matrix_subset.matrixSubsetSingleCToDoubleCAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                else:
                    raise Exception("order '{0}' not known, only 'F' and 'C'".format(effective_order));
            elif dtype == sp.float32:
                if effective_order=="F":
                    wrap_matrix_subset.matrixSubsetSingleCToSingleFAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                elif effective_order=="C":
                    wrap_matrix_subset.matrixSubsetSingleCToSingleCAAA(val, iid_count, sid_count, row_index_list, col_index_list, sub_val)
                else:
                    raise Exception("dtype '{0}' not known, only float64 and float32".format(dtype))
        else:
            raise Exception("input dtype '{0}' not known, only float64 and float32".format(val.dtype))
    else:
        raise Exception("input order must be 'F' or 'C'");


    logging.debug("Back from cython matrixSubset")
    return sub_val


def create_directory_if_necessary(name, isfile=True, robust=False):
    '''
    Create a directory for a file if the directory doesn't already exist.

    :param name: file or directory name
    :type name: string
    :param isfile: If True (default), the name is a file, otherwise it is a directory.
    :type isfile: bool
    :param robust: If False (default), will try once to create the directory. If True, will try 25 times.
    :type robust: bool

    '''
    if isfile:
        directory_name = os.path.dirname(name)
    else:
        directory_name = name

    if directory_name != "":
        if not robust:
            try:
                os.makedirs(directory_name)
            except OSError as e:
                if not os.path.isdir(directory_name):
                    raise Exception("not valid path: '{0}'. (Working directory is '{1}'".format(directory_name,os.getcwd()))
        else:
            is_ok = False
            time_to_sleep = 10.0
            for i in range(25):
                try:
                    os.makedirs(directory_name)
                    is_ok = True
                    break
                except OSError as e:
                    if not os.path.isdir(directory_name):
                        time_to_sleep *= 1.1
                        warnings.warn("creating directory robust=True, try#{0},time={3} error: not valid path: '{1}'. (Working directory is '{2}'".format(i, directory_name,os.getcwd(),int(time_to_sleep)))
                        time.sleep(int(time_to_sleep)) ; #make random?
                    else:
                        is_ok = True
                        break
            if not is_ok:
                raise Exception("not valid path: '{0}'. (Working directory is '{1}'".format(directory_name,os.getcwd()))


def weighted_mean(ys, weights):
    '''
    :param ys: The ndarray of values
    :type ys: ndarray
    :param weights: the weight of each value (unnormalized is fine)
    :type weights: ndarray
    :rtype: the weight mean


    >>> ys = np.array([103.664086,89.80645161,83.86888046,90.54141176])
    >>> weights = np.array([2.340862423,4.982888433,0.17522245,0.098562628])
    >>> round(weighted_mean(ys, weights),5)
    93.9487
    '''
    mean = ys.dot(weights)/weights.sum()
    return mean


def weighted_simple_linear_regression(xs, ys, weights):
    '''
    :param xs: The ndarray of independent values
    :type xs: ndarray
    :param ys: The ndarray of dependent values
    :type ys: ndarray
    :param weights: the weight of each case (unnormalized is fine)
    :type weights: ndarray
    :rtype: slope, intercept, xmean, ymean

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> xs = np.array([53.8329911,57.49486653,60.07392197,60.21081451])
    >>> ys = np.array([103.664086,89.80645161,83.86888046,90.54141176])
    >>> weights = np.array([2.340862423,4.982888433,0.17522245,0.098562628])
    >>> slope, intercept, xmean, ymean = weighted_simple_linear_regression(xs, ys, weights)
    >>> print(round(slope,5), round(intercept,5), round(xmean,5), round(ymean,5))
    -3.52643 293.05586 56.46133 93.9487

    '''
    xmean = weighted_mean(xs,weights)
    xs_less_mean = xs - xmean
    ymean = weighted_mean(ys,weights)
    ys_less_mean = ys - ymean
    weighted_xs_less_mean = xs_less_mean * weights
    slope = ys_less_mean.dot(weighted_xs_less_mean)/xs_less_mean.dot(weighted_xs_less_mean)
    intercept = ymean - xmean * slope
    return slope, intercept, xmean, ymean

def print2(arg):
    '''
    Make printing under Python3 look the same as under Python2.
    '''
    s = str(arg).replace("b'","'")
    print(s)

def format_delta(delta_seconds):
    '''
    Format a time delta nicely.

    :param delta_seconds: The number of seconds
    :type delta_seconds: number
    :rtype: string

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> from pysnptools.util import format_delta
    >>> print(format_delta(86403.5))
    1 day, 0:00:03.500000
    '''
    return datetime.timedelta(seconds=delta_seconds)

def _mbps(size, delta):
    return size * 8 / delta / 1e6

def _mbps_str(t0, size, total=0):
    delta = time.time()-t0
    mbps0 = _mbps(size, delta) if delta > 0 else 0
    percent = float(size)/total if total > 0 else 1
    left = delta/percent*(1-percent) if total > 0 else 0
    return "{0:0.2f}, {1:0.1f}%, left={2}".format(mbps0,percent*100,format_delta(left))


@contextmanager
def log_in_place(name, level, time_lambda=time.time, show_log_diffs=False):
    '''
        Create an one-argument function to write messages to. If the logging level
        is met, messages will appear. All messages will be on the same line and
        will include time.

        :Example:

        .. code-block:: python

            from pysnptools.util import log_in_place
            from six.moves import range
            import logging
            import time
            logging.basicConfig(level=logging.INFO)
         
            with log_in_place("counting", logging.INFO) as updater:
            for i in range(100):
                    updater(i)
                    time.sleep(.1) #typically, some work -- not a delay -- goes here

        Outputs 100 messages on the same line, ending with something like "counting -- time=0:00:09.99, 99"
    '''
    #!!! what if logging messages aren't suppose to go to stdout?
    t_wait = time_lambda()
    last_len = [0] #We have to make this an array so that the value is by reference.
    last_message_hash = [None]
    line_end = '\r'
    every_printed = [False] #Don't print the final newline if nothing is ever printed

    def writer(message):
        if logging.getLogger().level > level:
            return
        time_str = str(datetime.timedelta(seconds=time_lambda()-t_wait))
        if '.' in time_str:
            time_str = time_str[:time_str.index('.')+3] #Time to the 1/100th of a sec
        s = "{0} -- time={1}, {2}".format(name,time_str,message)
        if show_log_diffs:
            message_hash = hash(message)
            if message_hash !=  last_message_hash[0] and last_message_hash[0] is not None:
                sys.stdout.write('\n')
            last_message_hash[0] = message_hash
        sys.stdout.write("{0}{1}\r".format(s," "*max(0,last_len[0]-len(s)))) #Pad with spaces to cover up previous message
        every_printed[0] = True
        last_len[0] = len(s)

    yield writer

    if logging.getLogger().level > level:
        return
    if every_printed[0]:
        sys.stdout.write("\n")                

@contextmanager
def _file_transfer_reporter(name,size=None,updater=None):
    '''
    If an update is given, we use that. Otherwise, we create our own.
    '''
    if updater is None:
        bytes_so_far = [0] #We have to make this an array so that the value is by reference.
        t0 = time.time()

        with log_in_place(name, logging.INFO) as writer:
            def updater2(byte_size):
                bytes_so_far[0] += byte_size
                writer("Mbps={0}".format(_mbps_str(t0,bytes_so_far[0],size)))
            yield updater2
    else:
        yield updater


@contextmanager
def _multiopen(open_lambda, input_list):
    handle_list = [open_lambda(input) for input in input_list]        # Open the related inputs
    list_to_yield = [handle.__enter__() for handle in handle_list]    # Get the list to yield
    yield list_to_yield                                               # yield it
    for handle in handle_list:                                        # Close them
        handle.__exit__(None,None,None)

def _datestamp(appendrandom=False):
    import datetime
    now = datetime.datetime.now()
    s = str(now)[:19].replace(" ","_").replace(":","_")
    if appendrandom:
        import random
        s += "_" + str(random.random())[2:]
    return s

from pysnptools.util.generate import snp_gen

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
    import doctest
    doctest.testmod()
