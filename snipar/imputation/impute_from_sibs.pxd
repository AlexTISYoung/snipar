from libcpp.map cimport map as cmap
from libcpp.string cimport string as cstring
from libcpp.pair cimport pair as cpair
cimport numpy as cnp
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free

cdef int get_IBD_type(cstring id1,
                      cstring id2,
                      int loc,
                      cmap[cpair[cstring, cstring], vector[int]]& ibd_dict) nogil

cdef cpair[double, bint] impute_snp_from_offsprings(int snp,
                      int[:] sib_indexes,
                      int sib_count,
                      int[:, :] snp_ibd0,
                      int[:, :] snp_ibd1,
                      int[:, :] snp_ibd2,
                      float f,
                      double[:] parent_genotype_prob,
                      signed char[:, :, :] phased_gts,
                      signed char[:, :] unphased_gts,
                      int[:, :, :] sib_hap_IBDs,                      
                      int len_snp_ibd0,
                      int len_snp_ibd1,
                      int len_snp_ibd2,
                      bint use_backup
                      ) nogil


cdef cpair[double, cpair[int, bint]] impute_snp_from_parent_offsprings(int snp,
                      int parent,
                      int[:] sib_indexes,
                      int sib_count,
                      int[:, :] snp_ibd0,
                      int[:, :] snp_ibd1,
                      int[:, :] snp_ibd2,
                      float f,
                      double[:] parent_genotype_prob,
                      signed char[:, :, :] phased_gts,
                      signed char[:, :] unphased_gts,
                      int[:, :, :] sib_hap_IBDs,
                      int[:, :, :] parent_offspring_hap_IBDs,
                      int len_snp_ibd0,
                      int len_snp_ibd1,
                      int len_snp_ibd2,
                      bint use_backup
                      ) nogil

cdef cmap[cpair[cstring, cstring], vector[int]] dict_to_cmap(dict the_dict)

cdef bint is_possible_child(int child, int parent) nogil

cdef void get_IBD(signed char[:] hap1,
                  signed char[:] hap2,
                  int length,
                  int half_window,
                  double threshold,
                  int[:] agreement_count,
                  double[:] agreement_percentage,
                  int[:] agreement) nogil
