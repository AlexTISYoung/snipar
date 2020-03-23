# distutils: language = c++
#!/well/kong/users/wiw765/anaconda2/bin/python
import numpy as np
import pandas as pd
import logging
import time
from libcpp.map cimport map as cmap
from libcpp.string cimport string as cstring
from libcpp.pair cimport pair as cpair
from cpython cimport array
import array
cimport numpy as cnp
from libcpp.vector cimport vector
import cython
from libc.math cimport isnan
import h5py
from datetime import datetime

cdef float nan_float = np.nan
#TODO move readind IBD and pedigree to outside prepare_data

cdef cmap[cpair[cstring, cstring], vector[int]] dict_to_cmap(dict the_dict):
    #Converts a dictionary of (str, str)->int[] to cmap[cpair[cstring, cstring], vector[int]]
    cdef cpair[cstring,cstring] map_key
    cdef vector[int] map_val
    cdef cpair[cpair[cstring,cstring], vector[int]] map_element
    cdef cmap[cpair[cstring, cstring], vector[int]] c_dict
    for key,val in the_dict.items():
        map_key = key
        map_val = val
        map_element = (map_key, map_val)
        c_dict.insert(map_element)
    return c_dict
#TODO use scipy interface

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef float impute_snp(int snp, 
                      cnp.ndarray[cnp.int_t, ndim=2] snp_ibd0,
                      cnp.ndarray[cnp.int_t, ndim=2] snp_ibd1,
                      cnp.ndarray[cnp.int_t, ndim=2] snp_ibd2,
                      float f,
                      cnp.ndarray[cnp.float_t, ndim=2] bed,
                      int len_snp_ibd0,
                      int len_snp_ibd1,
                      int len_snp_ibd2):
    cdef float result = nan_float
    cdef float additive
    cdef float sibsum = 0
    cdef int sib1, sib2, pair_index
    if len_snp_ibd0 > 0:
        #if there is any ibd state0 we have observed all of the parents' genotypes,
        #therefore we can discard other ibd statuses
        result = 0        
        for pair_index in range(len_snp_ibd0):
            sib1 = snp_ibd0[pair_index, 0]
            sib2 = snp_ibd0[pair_index, 1]
            result += (bed[sib1, snp]+bed[sib2, snp])
        result = result/len_snp_ibd0

    elif len_snp_ibd1 > 0:
        #Because ibd2 is similar to having just one individual, we can discard ibd2s
        result = 0
        for pair_index in range(len_snp_ibd1):
            sib1 = snp_ibd1[pair_index, 0]
            sib2 = snp_ibd1[pair_index, 1]
            sibsum = (bed[sib1, snp]+bed[sib2, snp])
            additive = 0
            if sibsum==0:
                additive = f
            elif sibsum==1:
                additive = 1+f
            elif sibsum==2:
                additive = 1+2*f
            elif sibsum==3:
                additive = 2+f
            elif sibsum==4:
                additive = 3+f
            result += additive
        result = result/len_snp_ibd1
    
    elif len_snp_ibd2 > 0:
        #As ibd2 simillar to having one individual, we divide the sum of the pair by two
        result = 0
        for pair_index in range(len_snp_ibd2):
            sib1 = snp_ibd2[pair_index, 0]
            sib2 = snp_ibd2[pair_index, 1]
            result += (bed[sib1, snp]+bed[sib2, snp])/2. + 2*f
        result = result/len_snp_ibd2

    return result
        
cdef int get_IBD_type(cstring id1,
                      cstring id2,
                      int loc,
                      cmap[cpair[cstring, cstring], vector[int]]& ibd_dict):
    # TODO use get
    #the value for ibd_dict is like this: [start1, end1, ibd_type1, start2, end2, ibd_type2,...]
    cdef int result = 0
    cdef int index
    cdef cpair[cstring, cstring] key1
    cdef cpair[cstring, cstring] key2    
    cdef vector[int] segments
    key1.first = id1
    key1.second = id2
    key2.first = id2
    key2.second = id1

    if ibd_dict.count(key1) > 0:
        segments = ibd_dict[key1]

    elif ibd_dict.count(key2) > 0:
        segments = ibd_dict[key2]

    for index in range(segments.size()//3):
        if segments[3*index] <= loc <= segments[3*index+1]:
            result = segments[3*index+2]
            break

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def impute(sibships, iid_to_bed_index,  gts, ibd, pos, hdf5_output_dict, output_address = None):
    logging.info("imputing data ...")
    #converting python obejcts to c
    #sibships
    cdef int max_sibs = np.max(sibships["sib_count"])
    cdef int max_ibd_pairs = max_sibs*(max_sibs-1)/2
    cdef int number_of_fams = sibships.shape[0]
    cdef cnp.ndarray[cnp.double_t, ndim=1]freqs = np.nanmean(gts,axis=0)/2.0
    cdef vector[vector[cstring]] fams
    for fam in range(number_of_fams):
        fams.push_back(sibships["IID"].iloc[fam])
    cdef cnp.ndarray[cnp.int_t, ndim=1] sib_count = sibships["sib_count"].values    
    #iid_to_bed_index
    cdef cmap[cstring, int] c_iid_to_bed_index = iid_to_bed_index
    #gts
    cdef cnp.ndarray[cnp.float_t, ndim=2] c_gts = gts
    cdef int number_of_snps = c_gts.shape[1]
    #ibd
    cdef cmap[cpair[cstring, cstring], vector[int]] c_ibd = dict_to_cmap(ibd)
    #pos
    cdef cnp.ndarray[cnp.int_t, ndim=1] c_pos = pos
    
    cdef int len_snp_ibd0 = 0
    cdef int len_snp_ibd1 = 0
    cdef int len_snp_ibd2 = 0
    cdef cnp.ndarray[cnp.int_t, ndim=2] snp_ibd0 = np.zeros([max_ibd_pairs, 2], dtype=np.int)
    cdef cnp.ndarray[cnp.int_t, ndim=2] snp_ibd1 = np.zeros([max_ibd_pairs, 2], dtype=np.int)
    cdef cnp.ndarray[cnp.int_t, ndim=2] snp_ibd2 = np.zeros([max_ibd_pairs, 2], dtype=np.int)
    
    cdef int i, j, loc, ibd_type, sib1_index, sib2_index, progress
    cdef cstring sib1_id, sib2_id
    cdef cnp.ndarray[cnp.int_t, ndim=1] sibs_index = np.zeros(max_sibs).astype(int)
    cdef cnp.ndarray[cnp.double_t, ndim=2] imputed_par_gts = np.zeros((number_of_fams, number_of_snps))
    progress = -1
    for index in range(number_of_fams):
        if (index*100)//number_of_fams > progress:
            progress = (index*100)//number_of_fams
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print("["+ current_time +"] imputation progress:"+str(progress)+"% ["+ ("_"*(progress//2)) + ((50-progress//2-1)*"=") +"]")
        for i in range(sib_count[index]):
            sibs_index[i] = c_iid_to_bed_index[fams[index][i]]        
        for snp in range(number_of_snps):
            len_snp_ibd0 = 0
            len_snp_ibd1 = 0
            len_snp_ibd2 = 0
            loc = c_pos[snp]
            for i in range(1, sib_count[index]):
                for j in range(i):
                    sib1_index = sibs_index[i]
                    sib2_index = sibs_index[j]
                    sib1_id = fams[index][i]
                    sib2_id = fams[index][j]
                    sib1_gene_isnan = isnan(c_gts[sib1_index, snp])
                    sib2_gene_isnan = isnan(c_gts[sib2_index, snp])
                    ibd_type = get_IBD_type(sib1_id, sib2_id, loc, c_ibd)
                    if sib1_gene_isnan  and sib2_gene_isnan:
                        continue
                    #if one sib is nan, create a ibd2 pair consisting of the other sib
                    elif not sib1_gene_isnan  and sib2_gene_isnan:
                        snp_ibd2[len_snp_ibd2,0] = sib1_index
                        snp_ibd2[len_snp_ibd2,1] = sib1_index
                        len_snp_ibd2 += 1

                    elif sib1_gene_isnan  and not sib2_gene_isnan:
                        snp_ibd2[len_snp_ibd2,0] = sib2_index
                        snp_ibd2[len_snp_ibd2,1] = sib2_index
                        len_snp_ibd2 += 1

                    elif not sib1_gene_isnan and not sib2_gene_isnan:
                        if ibd_type == 2:
                            snp_ibd2[len_snp_ibd2,0] = sib1_index
                            snp_ibd2[len_snp_ibd2,1] = sib2_index
                            len_snp_ibd2 += 1
                        if ibd_type == 1:
                            snp_ibd1[len_snp_ibd1,0] = sib1_index
                            snp_ibd1[len_snp_ibd1,1] = sib2_index
                            len_snp_ibd1 += 1
                        if ibd_type == 0:
                            snp_ibd0[len_snp_ibd0,0] = sib1_index
                            snp_ibd0[len_snp_ibd0,1] = sib2_index
                            len_snp_ibd0 += 1
                    
            imputed_par_gts[index][snp] = impute_snp(snp, snp_ibd0, snp_ibd1, snp_ibd2, freqs[snp], c_gts, len_snp_ibd0, len_snp_ibd1, len_snp_ibd2)

    if output_address is not None:
        logging.info("Writing the results as a hdf5 file to "+output_address)
        with h5py.File(output_address,'w') as f:
            imputed_par_gts = imputed_par_gts
            f.create_dataset('imputed_par_gts',(number_of_fams, number_of_snps),dtype = 'f',chunks = True, compression = 'gzip', compression_opts=9, data = imputed_par_gts)
            f['families'] = sibships["FID"].values.tolist()
            for key, value in hdf5_output_dict.items():
                f[key] = value

    return sibships["FID"].values.tolist(), imputed_par_gts