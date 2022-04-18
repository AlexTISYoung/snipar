"""Contains functions in cython for doing the parent sum imputation from offsprings and parents(if they are observed).

Functions
----------
    get_probability_of_both_parents_conditioned_on_offsprings
    get_probability_of_one_parent_conditioned_on_offsprings_and_parent
    get_IBD
    get_hap_index
    is_possible_child
    dict_to_cmap
    impute_snp_from_offsprings
    impute_snp_from_parent_offsprings
    get_IBD_type
    impute
"""
# distutils: language = c++
import numpy as np
import logging
from libcpp.map cimport map as cmap
from libcpp.string cimport string as cstring
from libcpp.pair cimport pair as cpair
cimport numpy as cnp
from libcpp.vector cimport vector
import cython
import h5py
from cython.parallel import prange
cimport openmp
from snipar.config import nan_integer as python_integer_nan
from libc.stdio cimport printf
cdef float nan_float = np.nan
cdef int nan_integer = python_integer_nan
#prob_offspring_on_parent[i, j, k] shows the probability the offspring having the genotype k when the parents are i and j
cdef double[:,:,:] prob_offspring_on_parent = np.array(
[[[1.0, 0.0, 0.0],
 [0.5, 0.5, 0.0],
 [0.0, 1.0, 0.0],],

[[0.5, 0.5, 0.0],
 [0.25, 0.5, 0.25],
 [0.0, 0.5, 0.5],],
 
[[0.0, 1.0, 0.0],
 [0.0, 0.5, 0.5],
 [0.0, 0.0, 1.0],],]
 )


cdef double get_probability_of_both_parents_conditioned_on_offsprings(int snp,
                                                             int gp1,
                                                             int gp2,
                                                             int[:] sib_indexes,
                                                             int sib_count,
                                                             signed char [:,:] unphased_gts,
                                                             double[:] parent_genotype_prob) nogil:
    """Returns the probability of both parents being gp1 and gp2 given the offsprings specified in sib_indexes at SNP snp

        The function takes advantage of conditional independence of offsprings on parents. P(gp1, gp2| gs1, g2, ...) = P(gs1, gs2, ...|gp1, gp2)*P(gp1, gp2)/P(gs1, gs2, ...)=
        = P(gp1)P(gp2)/P(gs1, gs2, ...)*Prod[P(gsi|gp1, gp2)]
    Args:
        snp : int
            Index of the snp where the probability is computed.

        gp1 : int
            Genotype of parent 1

        gp2 : int
            Genotype of parent 2

        sib_indexes : int[:]
            Determines the gts index for each sibling from index of each sibling between offsprings

        sib_count : int
            Number of the offspring the parents have. The gts indexes of offsprings should be the first sib_count elements of sib_indexes.

        unphased_gts : signed char [:,:]
            A two-dimensional array containing genotypes for all individuals and SNPs respectively.

        parent_genotype_prob : double[:]
            An array with three elements, probability of parental genotype being 0, 1, 2 respectively.

    Returns:
        double
            Probability of parents being gp1 and gp2 given offsprings at SNP snp
    """
    cdef double numerator = parent_genotype_prob[gp1]*parent_genotype_prob[gp2]
    cdef double denumerator = 0.
    cdef int flag, gs, _gp1, _gp2, index
    cdef double tmp
    flag = 0
    for index in range(sib_count):
        gs = unphased_gts[sib_indexes[index], snp]
        if not (gs == nan_integer):
            numerator = numerator*prob_offspring_on_parent[gp1, gp2, gs]
            flag = 1
    if flag==0:
        return nan_float
    if numerator == 0:
        return 0
    for _gp1 in range(3):
        for _gp2 in range(3):
            tmp = 1
            for index in range(sib_count):
                gs = unphased_gts[sib_indexes[index], snp]
                tmp = tmp*prob_offspring_on_parent[_gp1, _gp2, gs]
            denumerator = denumerator+tmp*parent_genotype_prob[_gp1]*parent_genotype_prob[_gp2]
    if denumerator==0:
        return nan_float
    return numerator/denumerator


cdef double get_probability_of_one_parent_conditioned_on_offsprings_and_parent(int snp,
                                                                        int known_gp,
                                                                        int unknown_gp,
                                                                        int[:] sib_indexes,
                                                                        int sib_count,
                                                                        signed char [:,:] unphased_gts,
                                                                        double[:] parent_genotype_prob) nogil:
    """Returns the probability of a parents being unknown_gp given the other parent being known_gp the offsprings specified in sib_indexes at SNP snp

        The function takes advantage of conditional independence of offsprings on parents. P(unknown_gp| known_gp, gs1, gs2, ...) = P(unknown_gp, known_gp, gs1, gs2, ...)/P(known_gp, gs1, gs2, ...) = 
        = P(gs1, gs2, ...|known_gp, unknown_gp)*P(known_gp, unknown_gp)/P(known_gp, gs1, g2, ...) = 
        = Prod[gsi|known_gp, unknown_gp]*P(known_gp)P(unknown_gp)/Sum_unknown[P(unknown, known_gp, gs1, g2, ...)] = 
        = Prod[gsi|known_gp, unknown_gp]*P(known_gp)P(unknown_gp)/Sum_unknown[Prod[gsi|known_gp, unknown_gp]*P(known_gp)P(unknown_gp)] 

    Args:
        snp : int
            Index of the snp where the probability is computed.

        known_gp : int
            Genotype of the known parent(used in the condition)

        unknown_gp : int
            Genotype of the unknown parent(the one you are after its probability)

        sib_indexes : int[:]
            Determines the gts index for each sibling from index of each sibling between offsprings

        sib_count : int
            Number of the offspring the parents have. The gts indexes of offsprings should be the first sib_count elements of sib_indexes.

        unphased_gts : signed char [:,:]
            A two-dimensional array containing genotypes for all individuals and SNPs respectively.

        parent_genotype_prob : double[:]
            An array with three elements, probability of parental genotype being 0, 1, 2 respectively.

    Returns:
        double
            Probability of a parent being unknown_gp given offsprings and known_gp at SNP snp
    """
    cdef double numerator = parent_genotype_prob[known_gp]*parent_genotype_prob[unknown_gp]
    cdef double denumerator = 0
    cdef int flag, gs, _unknown_gp, index
    cdef double tmp
    flag = 0
    for index in range(sib_count):
        gs = unphased_gts[sib_indexes[index], snp]
        if not (gs == nan_integer):
            numerator = numerator*prob_offspring_on_parent[known_gp, unknown_gp, gs]
            flag = 1
    if flag==0:
        return nan_float
    if numerator==0:
        return 0.
    for _unknown_gp in range(3):
        tmp = 1
        for index in range(sib_count):
            gs = unphased_gts[sib_indexes[index], snp]
            if not (gs == nan_integer):
                tmp = tmp*prob_offspring_on_parent[known_gp, _unknown_gp, gs]
        denumerator = denumerator+tmp*parent_genotype_prob[known_gp]*parent_genotype_prob[_unknown_gp]
    if denumerator==0:
        return nan_float
    return numerator/denumerator

cdef extern from * nogil:
    r"""
    #include <omp.h>
    #include <stdio.h>  
    #include <string.h>
    static omp_lock_t cnt_lock;
    //the counter for progress
    static int cnt = 0;
    //Want to log the progress across different threads so made a lock for the variable cnt
    void reset(){
        omp_init_lock(&cnt_lock);
        cnt = 0;
    }
    void destroy(){
        omp_destroy_lock(&cnt_lock);
    }
    void report(int mod, char* chromosomes, int total){
        //writes the progress once every mod
        time_t now;
        char* text;
        omp_set_lock(&cnt_lock);
        cnt++;
        if(cnt%mod == 0){
            now = time(NULL);
            text = ctime(&now);
            text[strlen(text)-1] = 0;
            printf("%s INFO impute with chromosome %s: progress is %d\% \n", text, chromosomes, (100*cnt)/total);
        }
        omp_unset_lock(&cnt_lock);
    }
    """
    void reset()
    void destroy()
    void report(int mod, char* pre_message_info, int total)

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef void get_IBD(signed char[:] hap1,
                  signed char[:] hap2,
                  int length,
                  int half_window,
                  double threshold,
                  int[:] agreement_count,
                  double[:] agreement_percentage,
                  int[:] agreement) nogil:
    """Inferes IBD status between two haplotypes. For the location i, it checks [i-half_window, i+half_window] size, if they are the same on more than threshold portiona of the locations, it's IBD.
    Agreement, agreement_count, agreement percentage will contain the inffered IBDs, number of similarities in the window and its percentage for all the locations respectively.

    Args:
        hap1 : signed char[:]
            First haplotype

        hap2 : signed char[:]
            Second haplotype

        length : int
            Length of the haplotypes

        half_window : int
            For each location i, the IBD inference is restricted to [i-half_window, i+half_window] segment

        threshold : double
            We have an IBD segment if agreement_percentage is more than threshold

        agreement_count : int[:]
            For each location i, it's number of the time haplotypes agree with each other in the [i-half_window, i+half_window] window

        agreement_percentage : double[:]
            For each location i, it's the ratio of agreement between haplotypes in [i-half_window, i+half_window] window

        agreement : int[:]
            For each location i, it's the IBD status between haplotypes"""
    cdef int i
    cdef int first, last
    agreement_count[0] = 0
    last = min(half_window, length-1)
    for i in range(last+1):
        agreement_count[0] += (hap1[i] == hap2[i])
    agreement_percentage[0] = agreement_count[0]/<double>(last+1)

    for i in range(1, length):
        agreement_count[i] = agreement_count[i-1]
        last = i+half_window
        first = i-half_window-1
        if 0 <= first:
            agreement_count[i] = agreement_count[i-1]-(hap1[first] == hap2[first])
        if last < length:
            agreement_count[i] = agreement_count[i]+(hap1[last] == hap2[last])
        agreement_percentage[i] = agreement_count[i]/<double>(min(last+1, length) - max(0, first+1))
    
    for i in range(length):
        agreement[i] = (agreement_percentage[i]>threshold)
        if hap1[i] != hap2[i]:
            agreement[i] = 0


cdef int get_hap_index(int i, int j) nogil:
    """Maps an unordered pair of integers to a single integer. Mapping is unique and continous.
    Args:
        i : int
            
        j : int

    Returns:
        int"""
    if i > j:
        return i*(i-1)/2+j
    return j*(j-1)/2+i

cdef bint is_possible_child(int child, int parent) nogil:
    """Checks whether the child genotype is a possible offspring of the parent genotype. Returns False if any of them are nan
        Args:
            child : int
            parent : int

        Returns: bint
    """
    if parent == nan_integer or child == nan_integer:
        return False
    
    if (parent == 2) and (2 >= child > 0):
        return True

    if parent == 1 and (2 >= child >= 0):
        return True
    
    if (parent == 0) and (0 <= child < 2):
        return True

    return False

cdef cmap[cpair[cstring, cstring], vector[int]] dict_to_cmap(dict the_dict):
    """ Converts a (str,str)->list[int] dictionary to cmap[cpair[cstring, cstring], vector[int]]

    Args:
        the_dict : (str,str)->list[int]

    Returns:
        cmap[cpair[cstring, cstring], vector[int]]
    """
    cdef cpair[cstring,cstring] map_key
    cdef vector[int] map_val
    cdef cpair[cpair[cstring,cstring], vector[int]] map_element
    cdef cmap[cpair[cstring, cstring], vector[int]] c_dict
    for key,val in the_dict.items():
        map_key.first = key[0].encode('ASCII')
        map_key.second = key[1].encode('ASCII')
        map_val = val
        map_element = (map_key, map_val)
        c_dict.insert(map_element)
    return c_dict

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
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
                      bint use_backup,
                      ) nogil:
    """Imputes the parent sum divided by two for a single SNP from offsprings and returns the imputed value
        If phased_gts is not NULL, it tries to do the imputation with that first and if there is a problem, it falls back to unphased data.

    Args:
        snp : int
            index of each sibling between offsprings

        sib_count : int
            Number of the offspring the parents have. The gts indexes of offsprings should be the first sib_count elements of sib_indexes.
        
        sib_indexes: int[:]
            Determines the gts index for each sibling from index of each sibling between offsprings

        snp_ibd0 : int[:,:]
            List of sib pairs that are ibd0 in this SNP. It is assumed that there are len_snp_ibd0 sib pairs is this list.

        snp_ibd1 : int[:,:]
            List of sib pairs that are ibd1 in this SNP. It is assumed that there are len_snp_ibd1 sib pairs is this list

        snp_ibd2 : int[:,:]
            List of sib pairs that are ibd2 in this SNP. It is assumed that there are len_snp_ibd2 sib pairs is this list

        f : float
            Minimum allele frequency for the SNP.
        
        parent_genotype_prob : double[:]
            An array with three elements, probability of parental genotype being 0, 1, 2 respectively.
        
        phased_gts : signed char[:,:,:]
            A three-dimensional array containing genotypes for all individuals, SNPs and, haplotypes respectively.

        unphased_gts : signed char[:,:]
            A two-dimensional array containing genotypes for all individuals and SNPs respectively.
        
        sib_hap_IBDs: int[:, :, :]
            The IBD statuses of haplotypes. For each pair of siblings, first index is obtained by get_hap_index.
            Second index determines which haplotype pair(0 is 0-0, 1 is 0-1, 2 is 1-0, 3 is 1-1),
            third index the location of interest on the haplotypes

        len_snp_ibd0 : int
            The number of sibling pairs in snp_ibd0.

        len_snp_ibd1 : int
            The number of sibling pairs in snp_ibd1.

        len_snp_ibd2 : int
            The number of sibling pairs in snp_ibd2.
        
        use_backup : bint
            Whether it should use backup bayesian imputation where there is no ibd infomation available.

    Returns:
        cpair[double, bint]
            First value is imputed parent sum divided by two. NAN if all the children are NAN in this SNP. Second value is whether the imputation has been done using backup imputation.

    """

    cdef double result = nan_float
    cdef float additive
    cdef int sibsum = 0
    cdef int counter, sib1, sib2, pair_index, sib_index1, sib_index2, hap_index, h00, h01, h10, h11, gs10, gs11, gs20, gs21, gp1, gp2, gs1, gs2
    cdef cpair[double, bint] return_val
    #unless otherwise stated, it'll be false
    return_val.first = nan_float
    return_val.second = False


    if phased_gts != None:
        # The only time that having phased data matters is when we have IBD1
        if len_snp_ibd0==0 and len_snp_ibd1>0:
            result = 0
            counter = 0
            for pair_index in range(len_snp_ibd1):
                sib1 = snp_ibd1[pair_index, 0]
                sib2 = snp_ibd1[pair_index, 1]
                sib_index1 = sib_indexes[sib1]
                sib_index2 = sib_indexes[sib2]
                hap_index = get_hap_index(sib1, sib2)
                h00 = sib_hap_IBDs[hap_index, 0, snp]
                h01 = sib_hap_IBDs[hap_index, 1, snp]
                h10 = sib_hap_IBDs[hap_index, 2, snp]
                h11 = sib_hap_IBDs[hap_index, 3, snp]
                
                gs10 = phased_gts[sib_index1, snp, 0]
                gs11 = phased_gts[sib_index1, snp, 1]
                gs20 = phased_gts[sib_index2, snp, 0]
                gs21 = phased_gts[sib_index2, snp, 1]
                #checks whether inferred haplotype IBDs are consistend with the given IBD status
                if h00+h01+h10+h11 == 1:
                    # From the four observed alleles two are shared. So the imputation result is sum of f and the three distinct alleles divided by two.
                    if h00==1:
                        result += (f + gs10 + gs11 + gs21)/2
                    if h01==1:
                        result += (f + gs10 + gs11 + gs20)/2
                    if h10==1:
                        result += (f + gs11 + gs10 + gs21)/2
                    if h11==1:
                        result += (f + gs11 + gs10 + gs20)/2
                    counter += 1
            #TODO make sure things won't fall apart at the end of ibd segment
            if counter>0:
                result = result/counter
                return_val.first = result
                return return_val

    if use_backup and (len_snp_ibd0 == len_snp_ibd1 == len_snp_ibd2 == 0):
        result = 0
        for gp1 in range(3):
            for gp2 in range(3):
                result += (gp1+gp2)/2*get_probability_of_both_parents_conditioned_on_offsprings(snp, gp1, gp2, sib_indexes, sib_count, unphased_gts, parent_genotype_prob)
        if 0. <= result <= 2.:
            return_val.first = result
        else:            
            return_val.first = nan_float
        return_val.second = True
    
    elif len_snp_ibd0 > 0:
        #if there is any ibd state0 we have observed all of the parents' genotypes,
        #therefore we can discard other ibd statuses
        result = 0
        for pair_index in range(len_snp_ibd0):
            sib1 = sib_indexes[snp_ibd0[pair_index, 0]]
            sib2 = sib_indexes[snp_ibd0[pair_index, 1]]
            result += (unphased_gts[sib1, snp]+unphased_gts[sib2, snp])/2.
        return_val.first = result/len_snp_ibd0
    
    elif len_snp_ibd1 > 0:
        #Because ibd2 is similar to having just one individual, we can discard ibd2s
        result = 0
        for pair_index in range(len_snp_ibd1):
            sib1 = sib_indexes[snp_ibd1[pair_index, 0]]
            sib2 = sib_indexes[snp_ibd1[pair_index, 1]]
            gs1 = unphased_gts[sib1, snp]
            gs2 = unphased_gts[sib2, snp]
            if (gs1 == 0 and gs2 == 2) or (gs1 == 2 and gs2 == 0):
                break
            sibsum = gs1 + gs2
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
            result += additive/2.
        else:
            return_val.first = result/len_snp_ibd1

    elif len_snp_ibd2 > 0:
        #As ibd2 simillar to having one individual, we divide snp sum of the pair by two
        result = 0
        for pair_index in range(len_snp_ibd2):
            sib1 = sib_indexes[snp_ibd2[pair_index, 0]]
            sib2 = sib_indexes[snp_ibd2[pair_index, 1]]
            gs1 = unphased_gts[sib1, snp]
            gs2 = unphased_gts[sib2, snp]
            if gs1 != gs2:
                break
            result += (gs1+gs2)/4. + f
        else:
            return_val.first = result/len_snp_ibd2
    
    return return_val

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
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
                      bint use_backup,
                      ) nogil:
    """Imputes the missing parent for a single SNP from the other parent and offsprings and returns the imputed value, number of mendelian errors and whether backup linear imputation was used.
    
    If returns Nan if there are no sibling pairs that can be children of the existing parent.

    Args:
        snp : int
            The SNP index

        parent : int
            The index of parent's row in the bed matrix

        sib_indexes: int[:]
            Determines the gts index for each sibling from index of each sibling between offsprings

        sib_count : int
            Number of the offspring the parents have. The gts indexes of offsprings should be the first sib_count elements of sib_indexes.
        
        snp_ibd0 : cnp.ndarray[cnp.int_t, ndim=2]
            List of sib pairs that are ibd0 in this SNP. It is assumed that there are len_snp_ibd0 sib pairs is this list.

        snp_ibd1 : cnp.ndarray[cnp.int_t, ndim=2]
            List of sib pairs that are ibd1 in this SNP. It is assumed that there are len_snp_ibd1 sib pairs is this list

        snp_ibd2 : cnp.ndarray[cnp.int_t, ndim=2]
            List of sib pairs that are ibd2 in this SNP. It is assumed that there are len_snp_ibd2 sib pairs is this list

        f : float
            Minimum allele frequency for the SNP.

        phased_gts : signed char[:,:,:]
            A three-dimensional array containing genotypes for all individuals, SNPs and, haplotypes respectively.

        unphased_gts : signed char[:,:]
            A two-dimensional array containing genotypes for all individuals and SNPs respectively.

        sib_hap_IBDs: int[:, :, :]
            The IBD statuses of haplotypes. For each pair of siblings, first index is obtained by get_hap_index.
            Second index determines which haplotype pair(0 is 0-0, 1 is 0-1, 2 is 1-0, 3 is 1-1),
            third index the location of interest on the haplotypes

        parent_offspring_hap_IBDs: int[:, :, :]
            The IBD statuses of haplotypes. For each pair of parent offspring, first index is obtained by sib index between siblings.
            Second index determines which haplotype pair(0 is 0-0, 1 is 0-1, 2 is 1-0, 3 is 1-1),
            third index the location of interest on the haplotypes

        len_snp_ibd0 : int
            The number of sibling pairs in snp_ibd0.

        len_snp_ibd1 : int
            The number of sibling pairs in snp_ibd1.

        len_snp_ibd2 : int
            The number of sibling pairs in snp_ibd2.
        
        use_backup : bint
            Whether it should use backup bayesian imputation where there is no ibd infomation available.

    Returns:
        cpair[double, cpair[int, bint]]
            First value is imputed missing parent. NAN if all the children are NAN in this SNP or there is a mendelian error. Second value is number of mendelian errors and whether the imputation has been done using backup imputation accordingly.

    """    

    cdef double result
    cdef double additive
    cdef int gs1, gs2, gp1, i, dummy_gp
    cdef double sibsum = 0
    cdef int sib1, sib2, pair_index, counter, sib_index1, sib_index2, hap_index
    cdef int sibs_h00, sibs_h01, sibs_h10, sibs_h11, sibship_shared_allele_sib1, sibship_shared_allele_sib2
    cdef int parent_sib1_h00, parent_sib1_h01, parent_sib1_h10, parent_sib1_h11, parent_offspring1_shared_allele_parent, parent_offspring1_shared_allele_offspring
    cdef int parent_sib2_h00, parent_sib2_h01, parent_sib2_h10, parent_sib2_h11, parent_offspring2_shared_allele_parent, parent_offspring2_shared_allele_offspring
    cdef int gp = unphased_gts[parent, snp]
    cdef bint is_backup = False
    cdef int mendelian_error_count = 0
    cdef cpair[double, cpair[int, bint]] return_val
    #this is default
    return_val.first = nan_float
    return_val.second.first = 0
    return_val.second.second = False
    for i in range(sib_count):
        gs1 = unphased_gts[sib_indexes[i], snp]
        if not is_possible_child(<int> gs1, <int> gp):
            mendelian_error_count += 1
    return_val.second.first = mendelian_error_count

    if mendelian_error_count>0:
        return return_val    

    if phased_gts != None:
        #having phased data does not matter with IBD state 0
        if len_snp_ibd0==0 and len_snp_ibd1>0:
            result = 0
            counter = 0
            for pair_index in range(len_snp_ibd1):
                sib1 = snp_ibd1[pair_index, 0]
                sib2 = snp_ibd1[pair_index, 1]
                sib_index1 = sib_indexes[sib1]
                sib_index2 = sib_indexes[sib2]
                hap_index = get_hap_index(sib1, sib2)
                sibs_h00 = sib_hap_IBDs[hap_index, 0, snp]
                sibs_h01 = sib_hap_IBDs[hap_index, 1, snp]
                sibs_h10 = sib_hap_IBDs[hap_index, 2, snp]
                sibs_h11 = sib_hap_IBDs[hap_index, 3, snp]
                sibship_shared_allele_sib1 = sibs_h10 + sibs_h11
                sibship_shared_allele_sib2 = sibs_h01 + sibs_h11
                #checks inferred haplotype IBDs are consistent with the given IBD status
                if sibs_h00 + sibs_h10 + sibs_h01 + sibs_h11 != 1:
                    continue

                parent_sib1_h00 = parent_offspring_hap_IBDs[sib1, 0, snp]
                parent_sib1_h01 = parent_offspring_hap_IBDs[sib1, 1, snp]
                parent_sib1_h10 = parent_offspring_hap_IBDs[sib1, 2, snp]
                parent_sib1_h11 = parent_offspring_hap_IBDs[sib1, 3, snp]
                parent_offspring1_shared_allele_parent = parent_sib1_h10 + parent_sib1_h11
                parent_offspring1_shared_allele_offspring = parent_sib1_h01 + parent_sib1_h11
                #checks inferred haplotype IBDs are consistent with the natural IBD status
                if parent_sib1_h00 + parent_sib1_h10 + parent_sib1_h01 + parent_sib1_h11 != 1:
                    continue

                parent_sib2_h00 = parent_offspring_hap_IBDs[sib2, 0, snp]
                parent_sib2_h01 = parent_offspring_hap_IBDs[sib2, 1, snp]
                parent_sib2_h10 = parent_offspring_hap_IBDs[sib2, 2, snp]
                parent_sib2_h11 = parent_offspring_hap_IBDs[sib2, 3, snp]
                parent_offspring2_shared_allele_parent = parent_sib2_h10 + parent_sib2_h11
                parent_offspring2_shared_allele_offspring = parent_sib2_h01 + parent_sib2_h11
                #checks inferred haplotype IBDs are consistent with the natural IBD status
                if parent_sib2_h00 + parent_sib2_h10 + parent_sib2_h01 + parent_sib2_h11 != 1:
                    continue

                if parent_offspring1_shared_allele_offspring == sibship_shared_allele_sib1 and parent_offspring2_shared_allele_offspring == sibship_shared_allele_sib2:
                    #if the allele shared between offspring is also shared between those and the existing parent
                    result += phased_gts[sib_index1, snp, 1-parent_offspring1_shared_allele_offspring]+phased_gts[sib_index2, snp, 1-parent_offspring2_shared_allele_offspring]
                    counter+=1

                elif parent_offspring1_shared_allele_offspring != sibship_shared_allele_sib1 and parent_offspring2_shared_allele_offspring != sibship_shared_allele_sib2:
                    #if the allele shared between offspring is not shared between those and the existing parent
                    result += phased_gts[sib_index2, snp, sibship_shared_allele_sib1]+f
                    counter+=1
                # else:TODO
                    # printf("here is the bug")

            if counter > 0:
                return_val.first = result/counter
                return return_val

        elif len_snp_ibd0==0 and len_snp_ibd1==0 and len_snp_ibd2>0:
            result = 0
            counter = 0
            for pair_index in range(len_snp_ibd2):
                sib1 = snp_ibd2[pair_index, 0]
                sib2 = snp_ibd2[pair_index, 1]
                sib_index1 = sib_indexes[sib1]
                sib_index2 = sib_indexes[sib2]
                parent_sib1_h00 = parent_offspring_hap_IBDs[sib1, 0, snp]
                parent_sib1_h01 = parent_offspring_hap_IBDs[sib1, 1, snp]
                parent_sib1_h10 = parent_offspring_hap_IBDs[sib1, 2, snp]
                parent_sib1_h11 = parent_offspring_hap_IBDs[sib1, 3, snp]
                parent_offspring1_shared_allele_parent = parent_sib1_h10 + parent_sib1_h11
                parent_offspring1_shared_allele_offspring = parent_sib1_h01 + parent_sib1_h11
                #checks inferred haplotype IBDs are consistent with the natural IBD status
                if parent_sib1_h00 + parent_sib1_h10 + parent_sib1_h01 + parent_sib1_h11 != 1:
                    continue
                result += phased_gts[sib_index1, snp, 1-parent_offspring1_shared_allele_offspring]+f
                counter += 1
            if counter > 0:
                result = result/counter
                if 0. <= result <= 2.:
                    return_val.first = result
                else:            
                    return_val.first = nan_float
                return return_val

    result = nan_float
    counter = 0
    if use_backup and (len_snp_ibd0 == len_snp_ibd1 == len_snp_ibd2 == 0):
        result = 0
        for dummy_gp in range(3):
            #TODO this line has changed check if it fixes
            result += (dummy_gp)*get_probability_of_one_parent_conditioned_on_offsprings_and_parent(snp, gp, dummy_gp, sib_indexes, sib_count, unphased_gts, parent_genotype_prob)
        return_val.first = result
        return_val.second.second = True

    elif len_snp_ibd0 > 0:
        #if there is any ibd state0 we have observed all of the parents' genotypes,
        #therefore we can discard other ibd statuses
        result = 0
        for pair_index in range(len_snp_ibd0):
            sib1 = snp_ibd0[pair_index, 0]
            gs1 = unphased_gts[sib_indexes[sib1], snp]
            sib2 = snp_ibd0[pair_index, 1]
            gs2 = unphased_gts[sib_indexes[sib2], snp]
            additive = (gs1 + gs2 - gp)
            if additive>2 or additive<0:
                break
            result += additive
        else:
            return_val.first = result/len_snp_ibd0

    elif len_snp_ibd1 > 0:
        #Because ibd2 is similar to having just one individual, we can discard ibd2s
        result = 0
        for pair_index in range(len_snp_ibd1):
            sib1 = snp_ibd1[pair_index, 0]
            sib2 = snp_ibd1[pair_index, 1]
            gs1 = unphased_gts[sib_indexes[sib1], snp]
            gs2 = unphased_gts[sib_indexes[sib2], snp]            
            additive = 0
            if gp == 0 and (gs1 == 0 and gs2 == 0):
                additive = 0.5*f*(1-f)/((1-f)**2 + 0.5*f*(1-f))
            
            elif gp == 0 and ((gs1 == 0 and gs2 == 1) or (gs1 == 1 and gs2 == 0)):
                additive = 1

            elif gp == 0 and (gs1 == 1 and gs2 == 1):
                additive = (0.5*f*(1-f) + 2*f**2)/(0.5*f*(1-f)+f**2)

            elif gp == 1 and (gs1 == 0 and gs2 == 0):
                additive = 0
            
            elif gp == 1 and ((gs1 == 0 and gs2 == 1) or (gs1 == 1 and gs2 == 0)):
                additive = f*(1-f)/(0.5*(1-f)**2 + f*(1-f))

            elif gp == 1 and (gs1 == 1 and gs2 == 1):
                additive = 0.5*f**2/(0.25*f**2 + 0.25*(1-f)**2)

            elif gp == 1 and ((gs1 == 1 and gs2 == 2) or (gs1 == 2 and gs2 == 1)):
                additive = f*(1-f)/(f*(1-f) + 0.5*f**2) + f**2/(f*(1-f) + 0.5*f**2)
            
            elif gp == 1 and (gs1 == 2 and gs2 == 2):
                additive = 2
            
            elif gp == 2 and (gs1 == 1 and gs2 == 1):
                additive = 0.5*f*(1-f)/(0.5*f*(1-f)+(1-f)**2)

            elif gp == 2 and ((gs1 == 1 and gs2 == 2) or (gs1 == 2 and gs2 == 1)):
                additive = 1

            elif gp == 2 and (gs1 == 2 and gs2 == 2):
                additive = 0.5*f*(1-f)/(0.5*f*(1-f) + f**2) + 2*f**2/(0.5*f*(1-f) + f**2)
            else:
                break
            result += additive
        else:
            return_val.first = result/len_snp_ibd1        

    elif len_snp_ibd2 > 0:
        #As ibd2 simillar to having one individual, we dividsnpe the sum of the pair by two
        result = 0
        for pair_index in range(len_snp_ibd2):
            sib1 = snp_ibd2[pair_index, 0]
            sib2 = snp_ibd2[pair_index, 1]
            gs1 = unphased_gts[sib_indexes[sib1], snp]
            gs2 = unphased_gts[sib_indexes[sib2], snp]
            additive = 0
            if gs1 != gs2:
                break
            if gp == 0 and gs1 == 0:
                additive = f*(1-f)/((1-f)**2 + f*(1-f))

            elif gp == 0 and gs1 == 1:
                additive = (f*(1-f) + 2*(f**2))/(f*(1-f) + f**2)

            elif gp == 1 and gs1 == 0:
                additive = 0.5*f*(1-f)/(0.5*f*(1-f) + 0.5*(1-f)**2)

            elif gp == 1 and gs1 == 1:
                additive = (f*(1-f) + f**2)/(0.5*(1-f)**2 + f*(1-f) + 0.5*f**2)

            elif gp == 1 and gs1 == 2:
                additive = (0.5*f*(1-f) + f**2)/(0.5*f*(1-f) + 0.5*f**2)

            elif gp == 2 and gs1 == 1:
                additive = f*(1-f)/((1-f)**2 + f*(1-f))

            elif gp == 2 and gs1 == 2:
                additive = (f*(1-f) + 2*f**2)/(f*(1-f) + f**2)
            
            else:
                break
            result += additive
        else:
            return_val.first = result/len_snp_ibd2

    return return_val

cdef int get_IBD_type(cstring id1,
                      cstring id2,
                      int loc,
                      cmap[cpair[cstring, cstring], vector[int]]& ibd_dict) nogil:
    """Returns the IBD status of individuals with id1 and id2 in the SNP located at loc. Returns nan_integer if ambiguous.

    Args:
        id1 : cstring
            IID of individual 1

        id2 : cstring
            IID of individual 2

        loc : int
            Location of the SNP

        ibd_dict : cmap[cpair[cstring, cstring], vector[int]]
            A dictionary containing flattened IBD segments for each pair of related individuals.
            Each segment consists of three integers, start, end, and IBD_status (start and end are inclusive).
            Values are lists of integers in this fashion: [start0, end0, ibd_status0, start2, end2, ibd_status2, ...]
            Sibreg.imputation.preprocess_data.prepare_data can be used to create this.

    Returns:
        int
            the IBD status of individuals with id1 and id2 in the SNP located at loc. nan_integer if ambiguous.

    """

    #the value for ibd_dict is like this: [start1, end1, ibd_type1, start2, end2, ibd_type2,...]
    cdef int result = nan_integer
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
def impute(sibships, iid_to_bed_index,  phased_gts, unphased_gts, ibd, pos, hdf5_output_dict, chromosome, freqs, output_address = None, threads = None, output_compression = None, output_compression_opts = None, half_window=50, ibd_threshold = 0.999, silent_progress=False, use_backup=False):
    """Does the parent sum imputation for families in sibships and all the SNPs in unphased_gts and returns the results.

        Inputs and outputs of this function are ascii bytes instead of strings. It writes result of the imputation to the output_address.

    Args:
        sibships : pandas.Dataframe
            A pandas DataFrame with columns ['FID', 'FATHER_ID', 'MOTHER_ID', 'IID'] where IID columns is a list of the IIDs of individuals in that family.
            It only contains families with more than one child. The parental sum is computed for all these families.

        iid_to_bed_index : str->int
            A dictionary mapping IIDs of people to their location in the bed file.

        phased_gts : numpy.array[signed char]
            Numpy array containing the phased genotype data. Axes are individulas and SNPS respectively.
            It's elements should be 0 or 1 except NaN values which should be equal to nan_integer specified in the config.

        unphased_gts : numpy.array[signed char]
            Numpy array containing the unphased genotype data from a bed file. Axes are individulas, SNPS and haplotype number respectively.
            It's elements should be 0 or 1 except NaN values which should be equal to nan_integer specified in the config.


        ibd : pandas.Dataframe
            A pandas DataFrame with columns "ID1", "ID2", 'segment'. The segments column is a list of IBD segments between ID1 and ID2.
            Each segment consists of a start, an end, and an IBD status. The segment list is flattened meaning it's like [start0, end0, ibd_status0, start1, end1, ibd_status1, ...]

        pos : numpy.array
            A numpy array with the position of each SNP in the order of appearance in phased and unphased gts.
        
        hdf5_output_dict : dict
            Other key values to be added to the HDF5 output. Usually contains:                
                'bim_columns' : Columns of the resulting bim file
                'bim_values' : Contents of the resulting bim file
                'pedigree' : pedigree table Its columns are has_father, has_mother, single_parent respectively.
                'non_duplicates' : Indexes of the unique snps. Imputation is restricted to them.
                'standard_f' : Whether the allele frequencies are just population average instead of MAFs estimated using PCs

        chromosome: str
            Name of the chromosome(s) that's going to be imputed. Only used for logging purposes.

        freqs: list[float]
            A two-dimensional array containing estimated fs for all individuals and SNPs respectively.
        
        output_address : str, optional
            If presented, the results would be written to this address in HDF5 format.
            Aside from all the key, value pairs inside hdf5_output_dict, the following are also written to the file.
                'imputed_par_gts' : imputed genotypes
                'pos' : the position of SNPs(in the order of appearance in genotypes)                
                'families' : family ids of the imputed parents(in the order of appearance in genotypes)
                'parental_status' : a numpy array where each row shows the family status of the family of the corresponding row in families.
                'sib_ratio_backup' : An array with the size of number of snps. Show the ratio of backup imputation among offspring imputations in each snp.
                'parent_ratio_backup' : An array with the size of number of snps. Show the ratio of backup imputation among parent-offspring imputations in each snp.
                'mendelian_error_ratio' : Ratio of mendelian errors among parent-offspring pairs for each snp
                'estimated_genotyping_error' : estimated for each snp using mendelian_error_ratio and maf
                'ratio_ibd0' : ratio of families with offsprings in ibd0 to all the fams.
        
        threads : int, optional
            Specifies the Number of threads to be used. If None there will be only one thread.

        output_compression : str
            Optional compression algorithm used in writing the output as an hdf5 file. It can be either gzip or lzf. None means no compression.

        output_compression_opts : int
            Additional settings for the optional compression algorithm. Take a look at the create_dataset function of h5py library for more information. None means no compression setting.
        
        half_window : int, optional
            For each location i, the IBD inference for the haplotypes is restricted to [i-half_window, i+half_window].

        ibd_threshold : float, optional
            Minimum ratio of agreement between haplotypes for declaring IBD.
        
        silent_progress : boolean, optional
            Whether it should log the percentage of imputation's progress
        
        use_backup : boolean, optional
            Whether it should use backup imputation where there is no ibd infomation available. It's false by default.

    Returns:
        tuple(list, numpy.array)
            The second element is imputed parental genotypes and the first element is family ids of the imputed parents(in the order of appearance in the first element).                        
    """
    logging.info("with chromosome " + str(chromosome)+": " + "imputing ...")
    if sibships.empty:
        logging.warning("with chromosome " + str(chromosome)+": " + "Error: No families to be imputed")
        return [], np.array()

    cdef int number_of_threads = 1
    if threads is not None:
        number_of_threads = threads
    logging.info("with chromosome " + str(chromosome)+": " + "imputing data ...")
    #converting python obejcts to c
    #sibships
    cdef int max_sibs = np.max(sibships["sib_count"])
    cdef int max_ibd_pairs = max_sibs*(max_sibs-1)//2
    cdef int number_of_fams = sibships.shape[0]
    cdef cnp.ndarray[cnp.double_t, ndim=2] c_freqs = freqs
    cdef double f, fvars
    cdef int p, q
    sibships["parent"] = sibships["FATHER_ID"]
    sibships.loc[sibships["has_father"], "parent"] = sibships["FATHER_ID"][sibships["has_father"]]
    sibships.loc[sibships["has_mother"], "parent"] = sibships["MOTHER_ID"][sibships["has_mother"]]
    cdef vector[cstring] parents
    cdef vector[vector[cstring]] fams
    for fam in range(number_of_fams):
        fams.push_back(sibships["IID"].iloc[fam])
        parents.push_back(sibships["parent"].iloc[fam])
    cdef int[:] sib_count = sibships["sib_count"].values.astype("i")
    cdef cnp.ndarray[cnp.uint8_t, ndim=1] single_parent = sibships["single_parent"].astype('uint8').values    
    #iid_to_bed_index
    cdef cmap[cstring, int] c_iid_to_bed_index = iid_to_bed_index
    #unphased_gts
    cdef signed char[:, :] c_unphased_gts = unphased_gts
    cdef signed char[:, :, :] c_phased_gts = phased_gts
    cdef int number_of_snps = c_unphased_gts.shape[1]
    #ibd
    cdef cmap[cpair[cstring, cstring], vector[int]] c_ibd = dict_to_cmap(ibd)
    #pos
    cdef cnp.ndarray[cnp.int_t, ndim=1] c_pos = pos
    cdef int len_snp_ibd0 = 0
    cdef int len_snp_ibd1 = 0
    cdef int len_snp_ibd2 = 0
    cdef int[:,:,:] snp_ibd0 = np.ones([number_of_threads, max_ibd_pairs, 2], dtype=np.dtype("i"))
    cdef int[:,:,:] snp_ibd1 = np.ones([number_of_threads, max_ibd_pairs, 2], dtype=np.dtype("i"))
    cdef int[:,:,:] snp_ibd2 = np.ones([number_of_threads, max_ibd_pairs, 2], dtype=np.dtype("i"))
    cdef int i, j, loc, ibd_type, sib1_index, sib2_index, progress, where
    cdef cstring sib1_id, sib2_id
    cdef int[:, :] sibs_index = np.zeros((number_of_threads, max_sibs)).astype("i")
    cdef double[:, :] parent_genotype_prob = np.zeros((number_of_threads, 3))
    cdef double[:,:] imputed_par_gts = np.zeros((number_of_fams, number_of_snps))
    imputed_par_gts[:] = nan_float
    cdef int snp, this_thread, sib1_gene_isnan, sib2_gene_isnan, index
    byte_chromosome = chromosome.encode("ASCII")
    cdef char* chromosome_c = byte_chromosome
    cdef int mod = (number_of_fams+1)//100
    #For hap_ibds, axes denote thread, individual pair, haplotypes and SNPs.
    cdef int [:,:,:,:] sib_hap_IBDs = np.ones([number_of_threads, max_ibd_pairs, 4, number_of_snps], dtype=np.dtype("i"))
    cdef int [:,:,:,:] parent_offspring_hap_IBDs = np.ones([number_of_threads, max_sibs, 4, number_of_snps], dtype=np.dtype("i"))
    cdef double[:, :] agreement_percentages = np.zeros((number_of_threads, number_of_snps))
    cdef int[:, :] agreement_counts = np.ones([number_of_threads, number_of_snps], dtype=np.dtype("i"))
    cdef int half_window_c = half_window
    cdef float ibd_threshold_c = ibd_threshold
    cdef long[:] counter_ibd0 = np.zeros(number_of_snps).astype(long)
    cdef long[:] counter_nonnan_input = np.zeros(number_of_snps).astype(long)
    cdef bint is_backup
    cdef cpair[double, bint] o_result
    cdef cpair[double, cpair[int, bint]] po_result
    cdef long[:] sib_backup_count = np.zeros(number_of_snps).astype(long)
    cdef long[:] single_parent_backup_count = np.zeros(number_of_snps).astype(long)
    cdef long[:] single_parent_mendelian_error_count = np.zeros(number_of_snps).astype(long)
    cdef double[:] single_parent_fvars = np.zeros(number_of_snps)
    cdef int mendelian_error_count
    cdef bint c_silent_progress = silent_progress
    cdef bint c_use_backup = use_backup
    cdef int[:, :] sib_is_nan = np.zeros((number_of_threads, max_sibs)).astype("i")
    cdef bint has_non_nan_offspring = False
    cdef int[:] number_of_non_nan_offspring_per_snp = np.zeros(number_of_snps).astype("i")
    reset()
    logging.info("with chromosome " + str(chromosome)+": " + "using "+str(number_of_threads)+" threads")
    for index in prange(number_of_fams, nogil = True, num_threads = number_of_threads):
        if c_silent_progress == 0:
            report(mod, chromosome_c, number_of_fams)
        this_thread = openmp.omp_get_thread_num()
        for i in range(sib_count[index]):
            sibs_index[this_thread, i] = c_iid_to_bed_index[fams[index][i]]
        if c_phased_gts != None:
            # First fills hap_ibds
            for i in range(1, sib_count[index]):
                for j in range(i):
                    where = get_hap_index(i, j)
                    get_IBD(c_phased_gts[sibs_index[this_thread, i],:,0], c_phased_gts[sibs_index[this_thread, j],:,0], number_of_snps, half_window_c, ibd_threshold_c, agreement_counts[this_thread, :], agreement_percentages[this_thread, :], sib_hap_IBDs[this_thread, where, 0, :])
                    get_IBD(c_phased_gts[sibs_index[this_thread, i],:,0], c_phased_gts[sibs_index[this_thread, j],:,1], number_of_snps, half_window_c, ibd_threshold_c, agreement_counts[this_thread, :], agreement_percentages[this_thread, :], sib_hap_IBDs[this_thread, where, 1, :])
                    get_IBD(c_phased_gts[sibs_index[this_thread, i],:,1], c_phased_gts[sibs_index[this_thread, j],:,0], number_of_snps, half_window_c, ibd_threshold_c, agreement_counts[this_thread, :], agreement_percentages[this_thread, :], sib_hap_IBDs[this_thread, where, 2, :])
                    get_IBD(c_phased_gts[sibs_index[this_thread, i],:,1], c_phased_gts[sibs_index[this_thread, j],:,1], number_of_snps, half_window_c, ibd_threshold_c, agreement_counts[this_thread, :], agreement_percentages[this_thread, :], sib_hap_IBDs[this_thread, where, 3, :])

            if single_parent[index]:
                for i in range(0, sib_count[index]):
                    get_IBD(c_phased_gts[c_iid_to_bed_index[parents[index]],:,0], c_phased_gts[sibs_index[this_thread, i],:,0], number_of_snps, half_window_c, ibd_threshold_c, agreement_counts[this_thread, :], agreement_percentages[this_thread, :], parent_offspring_hap_IBDs[this_thread, i, 0, :])
                    get_IBD(c_phased_gts[c_iid_to_bed_index[parents[index]],:,0], c_phased_gts[sibs_index[this_thread, i],:,1], number_of_snps, half_window_c, ibd_threshold_c, agreement_counts[this_thread, :], agreement_percentages[this_thread, :], parent_offspring_hap_IBDs[this_thread, i, 1, :])
                    get_IBD(c_phased_gts[c_iid_to_bed_index[parents[index]],:,1], c_phased_gts[sibs_index[this_thread, i],:,0], number_of_snps, half_window_c, ibd_threshold_c, agreement_counts[this_thread, :], agreement_percentages[this_thread, :], parent_offspring_hap_IBDs[this_thread, i, 2, :])
                    get_IBD(c_phased_gts[c_iid_to_bed_index[parents[index]],:,1], c_phased_gts[sibs_index[this_thread, i],:,1], number_of_snps, half_window_c, ibd_threshold_c, agreement_counts[this_thread, :], agreement_percentages[this_thread, :], parent_offspring_hap_IBDs[this_thread, i, 3, :])
        snp = 0
        while snp < number_of_snps:
            len_snp_ibd0 = 0
            len_snp_ibd1 = 0
            len_snp_ibd2 = 0
            loc = c_pos[snp]
            has_non_nan_offspring = True
            for i in range(sib_count[index]):
                sib1_index = sibs_index[this_thread, i]
                sib_is_nan[this_thread, i] = (c_unphased_gts[sib1_index, snp] == nan_integer)
                if sib_is_nan[this_thread, i] > 0:
                    has_non_nan_offspring = True
            if has_non_nan_offspring:
                number_of_non_nan_offspring_per_snp[snp] += 1
                if sib_count[index] > 1:
                    for i in range(1, sib_count[index]):
                        for j in range(i):
                            sib1_id = fams[index][i]
                            sib2_id = fams[index][j]
                            sib1_gene_isnan = sib_is_nan[this_thread, i]
                            sib2_gene_isnan = sib_is_nan[this_thread, j]
                            ibd_type = get_IBD_type(sib1_id, sib2_id, loc, c_ibd)
                            if sib1_gene_isnan  and sib2_gene_isnan:
                                continue
                            
                            elif not sib1_gene_isnan  and sib2_gene_isnan:
                                snp_ibd2[this_thread, len_snp_ibd2,0] = i
                                snp_ibd2[this_thread, len_snp_ibd2,1] = i
                                len_snp_ibd2 = len_snp_ibd2+1

                            elif sib1_gene_isnan  and not sib2_gene_isnan:
                                snp_ibd2[this_thread, len_snp_ibd2,0] = j
                                snp_ibd2[this_thread, len_snp_ibd2,1] = j
                                len_snp_ibd2 = len_snp_ibd2 + 1

                            elif not sib1_gene_isnan and not sib2_gene_isnan:
                                if ibd_type == 2:
                                    snp_ibd2[this_thread, len_snp_ibd2,0] = i
                                    snp_ibd2[this_thread, len_snp_ibd2,1] = j
                                    len_snp_ibd2 = len_snp_ibd2 + 1
                                if ibd_type == 1:
                                    snp_ibd1[this_thread, len_snp_ibd1,0] = i
                                    snp_ibd1[this_thread, len_snp_ibd1,1] = j
                                    len_snp_ibd1 = len_snp_ibd1 + 1
                                if ibd_type == 0:
                                    snp_ibd0[this_thread, len_snp_ibd0,0] = i
                                    snp_ibd0[this_thread, len_snp_ibd0,1] = j
                                    len_snp_ibd0 = len_snp_ibd0 + 1
                    if len_snp_ibd0>0:
                        counter_ibd0[snp] = counter_ibd0[snp]+1
                else :
                    sib1_index = sibs_index[this_thread, 0]
                    if not (c_unphased_gts[sib1_index, snp] == nan_integer):
                        snp_ibd2[this_thread, len_snp_ibd2,0] = 0
                        snp_ibd2[this_thread, len_snp_ibd2,1] = 0
                        len_snp_ibd2 = len_snp_ibd2 + 1
                f = 0.
                fvars = 0.
                for p in range(sib_count[index]):
                    q = sibs_index[this_thread, p]
                    fvars = fvars + c_freqs[q, snp]*(1-c_freqs[q, snp])
                    f = f + c_freqs[q, snp]
                f = f/sib_count[index]
                parent_genotype_prob[this_thread, 0] = (1-f)**2
                parent_genotype_prob[this_thread, 1] = 2*f*(1-f)
                parent_genotype_prob[this_thread, 2] = f**2
                if single_parent[index] and c_unphased_gts[c_iid_to_bed_index[parents[index]], snp] != nan_integer:
                    #2 time MAF of offspring is MAF of sum of the parents. That minus the existing parent results in MAF of the missing parent.
                    f = 2*f-c_freqs[c_iid_to_bed_index[parents[index]], snp]
                    if f>1.:
                        f=1.
                    elif f<0.:
                        f=0.
                    po_result = impute_snp_from_parent_offsprings(snp,
                                                                                    c_iid_to_bed_index[parents[index]],
                                                                                    sibs_index[this_thread, :],
                                                                                    sib_count[index],
                                                                                    snp_ibd0[this_thread,:,:],
                                                                                    snp_ibd1[this_thread,:,:],
                                                                                    snp_ibd2[this_thread,:,:],
                                                                                    f,
                                                                                    parent_genotype_prob[this_thread, :],
                                                                                    c_phased_gts,
                                                                                    c_unphased_gts,
                                                                                    sib_hap_IBDs[this_thread,:,:,:],
                                                                                    parent_offspring_hap_IBDs[this_thread,:,:,:],
                                                                                    len_snp_ibd0,
                                                                                    len_snp_ibd1,
                                                                                    len_snp_ibd2,
                                                                                    c_use_backup,
                                                                                    )
                    imputed_par_gts[index, snp] = po_result.first
                    mendelian_error_count = po_result.second.first
                    single_parent_mendelian_error_count[snp] += mendelian_error_count
                    single_parent_fvars[snp] += fvars
                    is_backup = po_result.second.second
                    single_parent_backup_count[snp] += is_backup
                else:
                    o_result = impute_snp_from_offsprings(snp,
                                                                            sibs_index[this_thread, :],
                                                                            sib_count[index],
                                                                            snp_ibd0[this_thread,:,:],
                                                                            snp_ibd1[this_thread,:,:],
                                                                            snp_ibd2[this_thread,:,:],
                                                                            f,
                                                                            parent_genotype_prob[this_thread, :],
                                                                            c_phased_gts,
                                                                            c_unphased_gts,
                                                                            sib_hap_IBDs[this_thread,:,:,:],
                                                                            len_snp_ibd0,
                                                                            len_snp_ibd1,
                                                                            len_snp_ibd2,
                                                                            c_use_backup,
                                                                            )
                    imputed_par_gts[index, snp] = o_result.first
                    is_backup = o_result.second
                    sib_backup_count[snp] += is_backup
            snp = snp+1
    destroy()
    number_of_po_pairs = sum(sibships[sibships["single_parent"]]["sib_count"])
    mendelian_error_ratio = np.array([c/number_of_po_pairs if c!=0 else 0 for c in single_parent_mendelian_error_count])
    estimated_genotyping_error = np.array(single_parent_mendelian_error_count) / np.array(single_parent_fvars)
    multi_sib_fams = sum(sibships["sib_count"]>1)
    single_parent_fams = np.sum(single_parent)
    no_parent_fams = number_of_fams - single_parent_fams
    multi_sib_fams_ratio = multi_sib_fams/number_of_fams
    sib_ratio_backup = np.array([b/no_parent_fams  if no_parent_fams>0 else 0. for b in sib_backup_count])
    parent_ratio_backup = np.array([b/single_parent_fams  if single_parent_fams>0 else 0. for b in single_parent_backup_count])
    ratio_ibd0 = np.array([<float>counter_ibd0[i]/number_of_non_nan_offspring_per_snp[i] if number_of_non_nan_offspring_per_snp[i]>0 else 0 for i in range(number_of_snps)])
    total_ratio_ibd0 = (<float>np.sum(counter_ibd0))/np.sum(number_of_non_nan_offspring_per_snp) if np.sum(number_of_non_nan_offspring_per_snp)>0 else 0
    expected_total_ibd0 = 0
    for c in sib_count:
        expected_total_ibd0 += (1-0.5**(c-1))**2
    expected_total_ratio_ibd0 = expected_total_ibd0/number_of_fams
    logging.info(f"with chromosome {chromosome} :total number of fams is {number_of_fams}")
    logging.info(f"with chromosome {chromosome} :more than one offspring is genotyped for {multi_sib_fams} which is {multi_sib_fams_ratio:.2f} of total fams")
    logging.info(f"with chromosome {chromosome} :IBD0 state observed for {total_ratio_ibd0*100:.2f}%")
    logging.info(f"with chromosome {chromosome} :Expected to observe IBD0 state for {expected_total_ratio_ibd0*100:.2f}%")
    if use_backup:
        logging.info(f"with chromosome {chromosome} :ratio of backup imputation among sibs is {np.mean(sib_ratio_backup):.2f}")
        logging.info(f"with chromosome {chromosome} :ratio of backup imputation among parent-offsprings is {np.mean(parent_ratio_backup):.2f}")
    if total_ratio_ibd0>0.5:
        logging.warning("with chromosome " + str(chromosome)+": ibd0 ratio is too high")
    output_nan_count = np.sum(np.isnan(imputed_par_gts))
    nan_ratio = output_nan_count/imputed_par_gts.size
    logging.info("with chromosome " + str(chromosome)+f": total number of nans: {output_nan_count}, ratio of nan snps to all snps is {nan_ratio:.2f}")
    if nan_ratio > 0.01:
        logging.warning("Too much inconsistencies in the data")
    
    hdf5_output_dict["mendelian_error_ratio"] = mendelian_error_ratio
    hdf5_output_dict["estimated_genotyping_error"] = estimated_genotyping_error
    hdf5_output_dict["ratio_ibd0"] = ratio_ibd0
    hdf5_output_dict["sib_ratio_backup"] = sib_ratio_backup
    hdf5_output_dict["parent_ratio_backup"] = parent_ratio_backup
    hdf5_output_dict['families'] = np.array(sibships["FID"].values, dtype='S')
    hdf5_output_dict['parental_status'] = sibships[["has_father", "has_mother", "single_parent"]]
    hdf5_output_dict['pos'] = pos
    hdf5_output_dict['imputed_par_gts'] = imputed_par_gts
    if output_address is not None:
        logging.info("with chromosome " + str(chromosome)+": " + "Writing the results as a hdf5 file to "+output_address + ".hdf5")
        with h5py.File(output_address+".hdf5",'w') as file:                        
            for key, val in hdf5_output_dict.items():
                if key=='imputed_par_gts':
                    file.create_dataset(key, val.shape, dtype = 'float16', chunks = True, compression = output_compression, compression_opts=output_compression_opts, data = val)
                else:
                    file[key] = val
    return sibships["FID"].values.tolist(), np.array(imputed_par_gts)