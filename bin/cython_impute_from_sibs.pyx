#!/well/kong/users/wiw765/anaconda2/bin/python
import numpy as np
import pandas as pd
import logging
from pysnptools.snpreader import Bed
import json
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

cdef float nan_float = np.nan
def prepare_data(ped_address, genotypes_address, ibd_address, chr, start=None, end=None, bim_address = None):
    logging.info("initializing data")
    logging.info("loading and filtering pedigree file ...")
    #keeping individuals with no parents
    ped = pd.read_csv(ped_address, sep = " ")
    #TODO find parentless people smartly
    # ped["has_father"] = ~ ped["FATHER_ID"].str.endswith("_P")
    # ped["has_mother"] = ~ ped["MOTHER_ID"].str.endswith("_M")
    ped["has_father"] = ~ped["FATHER_ID"].isin(ped["IID"])
    ped["has_mother"] = ~ped["MOTHER_ID"].isin(ped["IID"])
    ped = ped[~(ped["has_mother"] | ped["has_father"])]
    #TODO handle sibs with one parent
    ped_ids =  set(ped["IID"].tolist())
    #finding siblings in each family
    sibships = ped.groupby(["FID", "FATHER_ID", "MOTHER_ID"]).agg({'IID':lambda x: list(x)}).reset_index()
    sibships["sib_count"] = sibships["IID"].apply(len)
    sibships = sibships[sibships["sib_count"]>1]
    logging.info("loading bim file ...")
    if bim_address is None:
        bim_file = genotypes_address+'.bim'
    else:
        bim_file = bim_address
    bim = pd.read_csv(bim_file, sep = "\t", header=None, names=["Chr", "id", "morgans", "coordinate", "allele1", "allele2"])
    #tODO what if people are related but do not have ibd on chrom
    logging.info("loading and transforming ibd file ...")
    ibd = pd.read_csv(ibd_address, sep = "\t").astype(str)
    #Adding location of start and end of each 
    ibd = ibd[ibd["Chr"] == str(chr)][["ID1", "ID2", "IBDType", "StartSNP", "StopSNP"]]
    ibd["IBDType"] = ibd["IBDType"].apply(lambda x: 2 if x=="IBD2" else 1)
    temp = bim[["id", "coordinate"]].rename(columns = {"id":"StartSNP","coordinate":"StartSNPLoc"})
    ibd= ibd.merge(temp, on="StartSNP")
    temp = bim[["id", "coordinate"]].rename(columns = {"id":"StopSNP","coordinate":"StopSNPLoc"})
    ibd = ibd.merge(temp, on="StopSNP")
    # ibd['segment'] = ibd[['StartSNPLoc', 'StopSNPLoc', "IBDType"]].apply(list, axis=1)
    ibd['segment'] = ibd[['StartSNPLoc', 'StopSNPLoc', "IBDType"]].values.tolist()
    def create_seg_list(x):
        elements = list(x)
        result = []
        for el in elements:
            result = result+el
        return result
    ibd = ibd.groupby(["ID1", "ID2"]).agg({'segment':lambda x:create_seg_list(x)}).to_dict()["segment"]
    #TODO orders
    logging.info("loading bed file ...")
    gts_f = Bed(genotypes_address+".bed")
    ids_in_ped = [(id in ped_ids) for id in gts_f.iid[:,1]]
    gts_ids = gts_f.iid[ids_in_ped]
    if end is not None:        
        gts = gts_f[ids_in_ped , start:end].read().val.astype(float)
        pos = gts_f.pos[start:end, 2]
        sid = gts_f.sid[start:end]
    else:
        gts = gts_f[ [(id in ped_ids) for id in gts_f.iid[:,1]], :].read().val.astype(float)
        pos = gts_f.pos[:, 2]
        sid = gts_f.sid
    iid_to_bed_index = {i:index for index, i in enumerate(gts_ids[:,1])}
    logging.info("initializing data done ...")
    return sibships, iid_to_bed_index, gts, ibd, pos


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

# #TODO use scipy interface

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
    #TODO check division by two
    #check for NANs
    cdef float result = nan_float

    #TODO throw exceptions instead
    #TODO what happens with nan
    cdef float additive
    cdef float sibsum = 0
    cdef int sib1, sib2, pair_index
    if len_snp_ibd0 > 0:
        result = 0        
        for pair_index in range(len_snp_ibd0):
            sib1 = snp_ibd0[pair_index, 0]
            sib2 = snp_ibd0[pair_index, 1]
            result += (bed[sib1, snp]+bed[sib2, snp])
        result = result/len_snp_ibd0

    elif len_snp_ibd1 > 0:
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
    # todo use get
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

    return result

@cython.wraparound(False)
@cython.boundscheck(False)
def impute(sibships, iid_to_bed_index,  gts, ibd, pos):
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
    
    cdef int i, j, loc, ibd_type, sib1_index, sib2_index
    cdef cstring sib1_id, sib2_id
    cdef cnp.ndarray[cnp.int_t, ndim=1] sibs_index = np.zeros(max_sibs).astype(int)
    cdef cnp.ndarray[cnp.double_t, ndim=2] imputed_par_gts = np.zeros((number_of_fams, number_of_snps))

    for index in range(number_of_fams):        
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

    return sibships["FID"].values.tolist(), imputed_par_gts



# python bin/impute_from_sibs_setup.py build_ext --inplace; mv cython_impute_from_sibs.so bin/; python bin/impute_runner.py --king --bim test_data/filtered_ukb_chr22.bim --start 100 --end 200 22 test_data/IBD.segments.gz test_data/filtered_ukb_chr22 test_data/pedigree test_data/parent_imputed_chr22
# python bin/impute_from_sibs_setup.py build_ext --inplace; python bin/impute_runner.py --bim test_data/filtered_ukb_chr22.bim --start 100 --end 200 22 test_data/IBD.segments.gz test_data/filtered_ukb_chr22 test_data/pedigree test_data/parent_imputed_chr22