from libcpp.map cimport map as cmap
from libcpp.string cimport string as cstring
from libcpp.pair cimport pair as cpair
cimport numpy as cnp
from libcpp.vector cimport vector


cdef int get_IBD_type(cstring id1,
                      cstring id2,
                      int loc,
                      cmap[cpair[cstring, cstring], vector[int]]& ibd_dict) nogil

cdef float impute_snp_from_offsprings(int snp, 
                                    int[:, :] snp_ibd0,
                                    int[:, :] snp_ibd1,
                                    int[:, :] snp_ibd2,
                                    float f,
                                    double[:, :] bed,
                                    int len_snp_ibd0,
                                    int len_snp_ibd1,
                                    int len_snp_ibd2) nogil


cdef float impute_snp_from_parent_offsprings(int snp,
                                            int parent,
                                            int[:, :] snp_ibd0,
                                            int[:, :] snp_ibd1,
                                            int[:, :] snp_ibd2,
                                            float f,
                                            double[:, :] bed,
                                            int len_snp_ibd0,
                                            int len_snp_ibd1,
                                            int len_snp_ibd2,
                                            ) nogil

cdef cmap[cpair[cstring, cstring], vector[int]] dict_to_cmap(dict the_dict)

cdef char is_possible_child(int child, int parent) nogil