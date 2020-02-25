from libcpp.map cimport map as cmap
from libcpp.string cimport string as cstring
from libcpp.pair cimport pair as cpair
cimport numpy as cnp
from libcpp.vector cimport vector


cdef int get_IBD_type(cstring id1,
                      cstring id2,
                      int loc,
                      cmap[cpair[cstring, cstring], vector[int]]& ibd_dict)

cdef float impute_snp(int snp, 
                      cnp.ndarray[cnp.int_t, ndim=2] snp_ibd0,
                      cnp.ndarray[cnp.int_t, ndim=2] snp_ibd1,
                      cnp.ndarray[cnp.int_t, ndim=2] snp_ibd2,
                      float f,
                      cnp.ndarray[cnp.int_t, ndim=2] bed,
                      int len_snp_ibd0,
                      int len_snp_ibd1,
                      int len_snp_ibd2)

cdef cmap[cpair[cstring, cstring], vector[int]] dict_to_cmap(dict the_dict)