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
from cython.operator cimport dereference as deref
cimport openmp
from snipar.config import nan_integer as python_integer_nan
from libc.stdio cimport printf
from libc.math cimport isnan
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
                if gs != nan_integer:
                    tmp *= prob_offspring_on_parent[_gp1, _gp2, gs]
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
    #include <atomic>
    #include <stdio.h>

    class SniparProgress {
        omp_lock_t lock;
        std::atomic<int> count{0};
        int last_reported = 0;
    public:
        SniparProgress() { omp_init_lock(&lock); }
        ~SniparProgress() { omp_destroy_lock(&lock); }
        void report(int interval, const char* chromosome, int total) {
            const int complete = 1 + count.fetch_add(1, std::memory_order_relaxed);
            if (complete % interval == 0 || complete == total) {
                omp_set_lock(&lock);
                if (complete > last_reported) {
                    printf("INFO impute with chromosome %s: progress is %d%%\n",
                           chromosome, (int)(100LL * complete / total));
                    last_reported = complete;
                }
                omp_unset_lock(&lock);
            }
        }
    };
    """
    cdef cppclass SniparProgress:
        void report(int interval, const char* chromosome, int total) noexcept

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
                  int[:] agreement) noexcept nogil:
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
    if length == 0:
        return
    agreement_count[0] = 0
    last = min(half_window, length-1)
    for i in range(last+1):
        agreement_count[0] += (hap1[i] == hap2[i] and hap1[i] != nan_integer)
    agreement_percentage[0] = agreement_count[0]/<double>(last+1)

    for i in range(1, length):
        agreement_count[i] = agreement_count[i-1]
        last = i+half_window
        first = i-half_window-1
        if 0 <= first:
            agreement_count[i] = agreement_count[i-1]-(hap1[first] == hap2[first] and hap1[first] != nan_integer)
        if last < length:
            agreement_count[i] = agreement_count[i]+(hap1[last] == hap2[last] and hap1[last] != nan_integer)
        agreement_percentage[i] = agreement_count[i]/<double>(min(last+1, length) - max(0, first+1))
    
    for i in range(length):
        agreement[i] = (agreement_percentage[i]>threshold)
        if hap1[i] != hap2[i] or hap1[i] == nan_integer:
            agreement[i] = 0


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef void fill_haplotype_ibd(signed char[:] hap1, signed char[:] hap2,
                             int length, int half_window, double threshold,
                             int[:] agreement) noexcept nogil:
    """Find agreement using a rolling count and one output write per marker."""
    cdef int i, first, last = min(half_window, length - 1), count = 0
    for i in range(last + 1):
        count += hap1[i] == hap2[i] and hap1[i] != nan_integer
    for i in range(length):
        last = i + half_window
        first = i - half_window
        agreement[i] = (hap1[i] == hap2[i] and hap1[i] != nan_integer and
                        count / <double>(min(last + 1, length) - max(0, first)) > threshold)
        if first >= 0:
            count -= hap1[first] == hap2[first] and hap1[first] != nan_integer
        if last + 1 < length:
            count += hap1[last + 1] == hap2[last + 1] and hap1[last + 1] != nan_integer



cdef int get_hap_index(int i, int j) nogil:
    """Maps an unordered pair of integers to a single integer. Mapping is unique and continuous.
    Args:
        i : int
            
        j : int

    Returns:
        int"""
    if i > j:
        return i*(i-1)//2+j
    return j*(j-1)//2+i

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
cdef cpair[double, cpair[bint, float]] impute_snp_from_offsprings(int snp,
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
                      ) noexcept nogil:
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
            Whether it should use backup bayesian imputation where there is no ibd information available.

    Returns:
        cpair[double, cpair[bint, float]]
            First value is imputed parent sum divided by two. NAN if all the children are NAN in this SNP. Second value is a pair: first value showing whether the imputation has been done using backup imputation and the second value showing number of observed parental alleles if it's known.

    """

    cdef double result = nan_float
    cdef float additive
    cdef int sibsum = 0
    cdef int counter, sib1, sib2, pair_index, sib_index1, sib_index2, hap_index, h00, h01, h10, h11, gs10, gs11, gs20, gs21, gp1, gp2, gs1, gs2
    cdef cpair[double, cpair[bint, float]] return_val
    #[imputed, [backup, num_observed_parental alleles]]
    #unless otherwise stated, it'll be false
    return_val.first = nan_float
    return_val.second.first = False
    return_val.second.second = nan_float


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
                #checks whether inferred haplotype IBDs are consistent with the given IBD status
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
                return_val.second.second = 3
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
        return_val.second.first = True
        return_val.second.second = nan_float
    
    elif len_snp_ibd0 > 0:
        #if there is any ibd state0 we have observed all of the parents' genotypes,
        #therefore we can discard other ibd statuses        
        result = 0
        for pair_index in range(len_snp_ibd0):
            sib1 = sib_indexes[snp_ibd0[pair_index, 0]]
            sib2 = sib_indexes[snp_ibd0[pair_index, 1]]
            result += (unphased_gts[sib1, snp]+unphased_gts[sib2, snp])/2.
        return_val.first = result/len_snp_ibd0
        return_val.second.second = 4
    
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
            return_val.second.second = 3

    elif len_snp_ibd2 > 0:
        num_observed_parental_alleles = 2
        #As ibd2 similar to having one individual, we divide snp sum of the pair by two
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
            return_val.second.second = 2
    
    return return_val

@cython.cdivision(True)
cdef inline double impute_parent_child(float f, int known, int child) noexcept nogil:
    """Expected genotype of the missing parent for a feasible parent-child pair."""
    if known == 1:
        if child == 0:
            return 0.5*f*(1-f)/(0.5*f*(1-f) + 0.5*(1-f)**2)
        if child == 1:
            return (f*(1-f) + f**2)/(0.5*(1-f)**2 + f*(1-f) + 0.5*f**2)
        return (0.5*f*(1-f) + f**2)/(0.5*f*(1-f) + 0.5*f**2)
    if child == (known // 2):
        return f*(1-f)/((1-f)**2 + f*(1-f))
    return (f*(1-f) + 2*f**2)/(f*(1-f) + f**2)


@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int transmitted_haplotype(int row, int marker, int parent_dosage,
                                      int child, int snp,
                                      const signed char[:, :] parent_origin,
                                      signed char[:, :, :] phased_gts,
                                      signed char[:, :] unphased_gts) noexcept nogil:
    if row < 0:
        return -1
    cdef int origin = parent_origin[row, marker]
    if not 0 <= origin <= 1:
        return -1
    cdef int observed = phased_gts[child, snp, origin]
    cdef int missing = phased_gts[child, snp, 1 - origin]
    if (0 <= observed <= 1 and 0 <= missing <= 1 and
            0 <= parent_dosage - observed <= 1 and
            observed + missing == unphased_gts[child, snp]):
        return origin
    return -1


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef cpair[double, float] impute_known_origin(
        int snp, int marker, int parent, int parent_dosage, float f,
        size_t family_start, int[:] sib_indexes, const int[:] origin_rows,
        const signed char[:, :] parent_origin, const signed char[:, :] parent_haplotype,
        signed char[:, :, :] phased_gts, signed char[:, :] unphased_gts,
        int[:, :, :] sib_hap_IBDs, int[:, :] pairs, int n_pairs,
        bint shares_one) noexcept nogil:
    """Apply Supplementary Note §5.2 to offspring with identified transmissions."""
    cdef int pair, first, second, child1, child2, row1, row2, origin1, origin2
    cdef int missing1, missing2, observed1, observed2, p1, p2, hap_pair
    cdef int shared_observed, shared_missing, same_parent_copy
    cdef int observed_sum = 0, observed_count = 0, complete_sum = 0, complete_count = 0
    cdef cpair[double, float] result
    result.second = 0
    for pair in range(n_pairs):
        first = pairs[pair, 0]
        second = pairs[pair, 1]
        child1 = sib_indexes[first]
        child2 = sib_indexes[second]
        row1 = origin_rows[family_start + first]
        row2 = origin_rows[family_start + second]
        origin1 = transmitted_haplotype(row1, marker, parent_dosage, child1, snp,
                                        parent_origin, phased_gts, unphased_gts)
        origin2 = transmitted_haplotype(row2, marker, parent_dosage, child2, snp,
                                        parent_origin, phased_gts, unphased_gts)
        if not shares_one:
            if (origin1 < 0 and origin2 < 0 or
                    unphased_gts[child1, snp] != unphased_gts[child2, snp]):
                continue
            if origin1 >= 0:
                missing1 = phased_gts[child1, snp, 1 - origin1]
                if origin2 >= 0 and missing1 != phased_gts[child2, snp, 1 - origin2]:
                    continue
            else:
                missing1 = phased_gts[child2, snp, 1 - origin2]
            observed_sum += missing1
            observed_count += 1
            continue

        if origin1 < 0 or origin2 < 0:
            continue
        observed1 = phased_gts[child1, snp, origin1]
        observed2 = phased_gts[child2, snp, origin2]
        missing1 = phased_gts[child1, snp, 1 - origin1]
        missing2 = phased_gts[child2, snp, 1 - origin2]
        same_parent_copy = -1
        if parent_dosage == 1:
            # A heterozygous parent's allele identifies its transmitted homolog.
            same_parent_copy = observed1 == observed2
        elif parent_haplotype != None:
            p1 = parent_haplotype[row1, marker]
            p2 = parent_haplotype[row2, marker]
            if (0 <= p1 <= 1 and 0 <= p2 <= 1 and
                    phased_gts[parent, snp, p1] == observed1 and
                    phased_gts[parent, snp, p2] == observed2):
                same_parent_copy = p1 == p2
        if same_parent_copy < 0:
            if missing1 != missing2:
                same_parent_copy = 1
            else:
                # With identical alleles, use sibling sharing along the two
                # possible same-parent alignments of the input haplotypes.
                hap_pair = get_hap_index(first, second)
                shared_observed = sib_hap_IBDs[hap_pair, 2 * origin1 + origin2, snp]
                shared_missing = sib_hap_IBDs[hap_pair, 3 - (2 * origin1 + origin2), snp]
                if shared_observed == shared_missing:
                    continue
                same_parent_copy = shared_observed
        if same_parent_copy:
            complete_sum += missing1 + missing2
            complete_count += 1
        elif missing1 == missing2:
            observed_sum += missing1
            observed_count += 1

    if complete_count:
        result.first = <double> complete_sum / complete_count
        result.second = 4
    elif observed_count:
        result.first = <double> observed_sum / observed_count + f
        result.second = 3
    return result


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef cpair[cpair[double, float], cpair[int, bint]] impute_snp_from_parent_offsprings(int snp,
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
                      int origin_snp,
                      size_t family_start,
                      const int[:] origin_rows,
                      const signed char[:, :] parent_origin,
                      const signed char[:, :] parent_haplotype,
                      ) noexcept nogil:
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
        
        snp_ibd0 : cnp.ndarray[cnp.int64_t, ndim=2]
            List of sib pairs that are ibd0 in this SNP. It is assumed that there are len_snp_ibd0 sib pairs is this list.

        snp_ibd1 : cnp.ndarray[cnp.int64_t, ndim=2]
            List of sib pairs that are ibd1 in this SNP. It is assumed that there are len_snp_ibd1 sib pairs is this list

        snp_ibd2 : cnp.ndarray[cnp.int64_t, ndim=2]
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
            Whether it should use backup bayesian imputation where there is no ibd information available.

    Returns:
        cpair[cpair[double, float], cpair[int, bint]]
            First value is a pair of imputed missing parent and number of observed parental alleles. NAN if all the children are NAN in this SNP or there is a mendelian error. Second value is number of mendelian errors and whether the imputation has been done using backup imputation accordingly.

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
    cdef float num_observed_parental_alleles = nan_float
    cdef cpair[cpair[double, float], cpair[int, bint]] return_val
    cdef cpair[double, float] origin_result
    #this is default
    return_val.first.first = nan_float
    return_val.first.second = nan_float
    return_val.second.first = 0
    return_val.second.second = False
    for i in range(sib_count):
        gs1 = unphased_gts[sib_indexes[i], snp]
        if not is_possible_child(<int> gs1, <int> gp):
            mendelian_error_count += 1
    return_val.second.first = mendelian_error_count

    if mendelian_error_count>0:
        return return_val    

    if parent_origin != None and len_snp_ibd0 == 0:
        if len_snp_ibd1 > 0:
            origin_result = impute_known_origin(
                snp, origin_snp, parent, gp, f, family_start, sib_indexes, origin_rows,
                parent_origin, parent_haplotype, phased_gts, unphased_gts, sib_hap_IBDs,
                snp_ibd1, len_snp_ibd1, True)
        elif len_snp_ibd2 > 0:
            origin_result = impute_known_origin(
                snp, origin_snp, parent, gp, f, family_start, sib_indexes, origin_rows,
                parent_origin, parent_haplotype, phased_gts, unphased_gts, sib_hap_IBDs,
                snp_ibd2, len_snp_ibd2, False)
        else:
            origin_result.second = 0
        if origin_result.second:
            return_val.first = origin_result
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
                    #if the allele shared between offsprings is also shared between those and the existing parent
                    #TODO shouldn't num_observed_parental_alleles be averaged, both in the case and next elif
                    num_observed_parental_alleles = 4
                    result += phased_gts[sib_index1, snp, 1-parent_offspring1_shared_allele_offspring]+phased_gts[sib_index2, snp, 1-parent_offspring2_shared_allele_offspring]
                    counter+=1

                elif parent_offspring1_shared_allele_offspring != sibship_shared_allele_sib1 and parent_offspring2_shared_allele_offspring != sibship_shared_allele_sib2:
                    #if the allele shared between offsprings is not shared between those and the existing parent
                    num_observed_parental_alleles = 3
                    result += phased_gts[sib_index2, snp, sibship_shared_allele_sib1]+f
                    counter+=1
                # else:TODO
                    # printf("here is the bug")

            if counter > 0:
                return_val.first.first = result/counter
                return_val.first.second = num_observed_parental_alleles
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
                    return_val.first.first = result
                    return_val.first.second = 3
                else:            
                    return_val.first.first = nan_float
                    return_val.first.second = nan_float
                return return_val

    result = nan_float
    counter = 0
    if use_backup and (len_snp_ibd0 == len_snp_ibd1 == len_snp_ibd2 == 0):
        result = 0
        for dummy_gp in range(3):
            #TODO this line has changed check if it fixes
            result += (dummy_gp)*get_probability_of_one_parent_conditioned_on_offsprings_and_parent(snp, gp, dummy_gp, sib_indexes, sib_count, unphased_gts, parent_genotype_prob)
        return_val.first.first = result
        return_val.first.second = nan_float
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
            return_val.first.first = result/len_snp_ibd0
            return_val.first.second = 4

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
            return_val.first.first = result/len_snp_ibd1
            return_val.first.second = nan_float

    elif len_snp_ibd2 > 0:
        #As ibd2 similar to having one individual, we dividsnpe the sum of the pair by two
        result = 0
        for pair_index in range(len_snp_ibd2):
            sib1 = snp_ibd2[pair_index, 0]
            sib2 = snp_ibd2[pair_index, 1]
            gs1 = unphased_gts[sib_indexes[sib1], snp]
            gs2 = unphased_gts[sib_indexes[sib2], snp]
            if gs1 != gs2:
                break
            result += impute_parent_child(f, gp, gs1)
        else:
            return_val.first.first = result/len_snp_ibd2
            return_val.first.second = 3

    return return_val

@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int posterior_mode(const float[:, :, ::1] probabilities, int row, int marker) noexcept nogil:
    cdef int state, selected = -1
    cdef float best = 0
    for state in range(3):
        if probabilities[row, marker, state] > best:
            selected = state
            best = probabilities[row, marker, state]
    return selected


@cython.cdivision(True)
cdef inline float expected_parent_alleles_ibd1(float f, int parent, int first, int second) noexcept nogil:
    # IBD1 shares either the known parent's homolog (four observed alleles) or
    # the missing parent's homolog (three). Condition both possibilities on
    # the observed offspring dosages and the missing parent's allele frequency.
    cdef int copy, known, left, right
    cdef double shared_known = 0, shared_missing = 0
    for copy in range(2):
        known = parent > copy
        left, right = first - known, second - known
        if 0 <= left <= 1 and 0 <= right <= 1:
            shared_known += (f if left else 1 - f) * (f if right else 1 - f)
        right = second - (parent - known)
        if 0 <= left <= 1 and left == right:
            shared_missing += f if left else 1 - f
    if shared_known + shared_missing > 0:
        return 3 + shared_known / (shared_known + shared_missing)
    return nan_float


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef cpair[double, cpair[float, bint]] impute_pair_posterior(
        int snp, int out_snp, int row, const float[:, :, ::1] probabilities,
        int parent, int[:] siblings, int[:, :] ibd0, int[:, :] ibd1, int[:, :] ibd2,
        float f, double[:] genotype_probabilities, signed char[:, :, :] phased,
        signed char[:, :] dosages, int[:, :, :] sibling_haplotypes,
        int[:, :, :] duo_haplotypes, size_t family_start, int[:] origin_rows,
        const signed char[:, :] parent_origin, const signed char[:, :] parent_haplotype,
        int hard_state, double hard_dosage, float hard_observed) noexcept nogil:
    cdef cpair[double, cpair[float, bint]] result
    cdef cpair[double, cpair[bint, float]] without_parent
    cdef cpair[cpair[double, float], cpair[int, bint]] with_parent
    cdef double dosage, total_dosage = 0, total_observed = 0, total_weight = 0
    cdef float weight, observed
    cdef int state
    ibd0[0, 0], ibd0[0, 1] = 1, 0
    ibd1[0, 0], ibd1[0, 1] = 1, 0
    ibd2[0, 0], ibd2[0, 1] = 1, 0
    for state in range(3):
        weight = probabilities[row, out_snp, state]
        if weight <= 0:
            continue
        if state == hard_state:
            dosage = hard_dosage
            observed = hard_observed
        elif parent < 0:
            without_parent = impute_snp_from_offsprings(
                snp, siblings, 2, ibd0, ibd1, ibd2, f, genotype_probabilities, phased,
                dosages, sibling_haplotypes, state == 0, state == 1, state == 2, False)
            dosage = without_parent.first
            observed = without_parent.second.second
        else:
            with_parent = impute_snp_from_parent_offsprings(
                snp, parent, siblings, 2, ibd0, ibd1, ibd2, f, genotype_probabilities, phased,
                dosages, sibling_haplotypes, duo_haplotypes, state == 0, state == 1, state == 2,
                False, out_snp, family_start, origin_rows, parent_origin, parent_haplotype)
            dosage = with_parent.first.first
            observed = with_parent.first.second
        if not isnan(dosage):
            if state == 1 and parent >= 0 and isnan(observed):
                observed = expected_parent_alleles_ibd1(
                    f, dosages[parent, snp], dosages[siblings[0], snp], dosages[siblings[1], snp])
            total_weight += weight
            total_dosage += weight * dosage
            total_observed += weight * observed
    result.first = hard_dosage
    result.second.first = hard_observed
    result.second.second = total_weight > 0
    if total_weight > 0:
        result.first = total_dosage / total_weight
        result.second.first = total_observed / total_weight
    return result


cdef const vector[int]* find_ibd_segments(
        const cstring& id1, const cstring& id2,
        const cmap[cpair[cstring, cstring], vector[int]]& ibd_dict) noexcept nogil:
    cdef cmap[cpair[cstring, cstring], vector[int]].const_iterator it
    it = ibd_dict.find(cpair[cstring, cstring](id1, id2))
    if it == ibd_dict.end():
        it = ibd_dict.find(cpair[cstring, cstring](id2, id1))
    if it == ibd_dict.end():
        return NULL
    return &deref(it).second

cdef int ibd_at_position(const vector[int]* segments, int loc) noexcept nogil:
    cdef size_t index
    if segments != NULL:
        for index in range(0, deref(segments).size(), 3):
            if deref(segments)[index] <= loc <= deref(segments)[index + 1]:
                return deref(segments)[index + 2]
    return nan_integer

cdef int get_IBD_type(cstring id1, cstring id2, int loc,
                      cmap[cpair[cstring, cstring], vector[int]]& ibd_dict) nogil:
    """Find the IBD state at an inclusive physical coordinate for either pair orientation."""
    return ibd_at_position(find_ibd_segments(id1, id2, ibd_dict), loc)

cdef class ImputationFamilies:
    """Compact, reusable sample and IBD segment indexes for one chromosome."""

    cdef vector[int] parents
    cdef vector[int] children
    cdef vector[size_t] offspring_offsets
    cdef vector[size_t] family_pair_offsets
    cdef vector[size_t] segment_offsets
    cdef vector[int] segments
    cdef vector[unsigned char] sorted_segments
    cdef vector[int] posterior_rows
    cdef int n_posterior_rows
    cdef int max_sibs
    cdef int n_samples
    cdef object family_ids
    cdef object parental_status
    cdef object single_parent
    cdef long long n_parent_pairs
    cdef int n_single_parents
    cdef int n_multi_sib
    cdef double expected_ibd0

    def __init__(self, sibships, iid_to_bed_index, ibd, int n_samples, *, posterior_index=None):
        cdef int i, j, offset, sample_index, segment_start, segment_end, state
        cdef object first, second, segment
        cdef bint sorted_pair
        cdef int previous_end
        self.n_samples = n_samples
        self.n_posterior_rows = 0 if posterior_index is None else len(posterior_index)
        self.family_ids = np.asarray(sibships["FID"].values, dtype='S')
        self.parental_status = sibships[["has_father", "has_mother", "single_parent"]].to_numpy(dtype=bool)
        self.single_parent = self.parental_status[:, 2].astype(np.uint8)
        self.max_sibs = 0
        self.offspring_offsets.push_back(0)
        self.family_pair_offsets.push_back(0)
        self.segment_offsets.push_back(0)
        for father, mother, iids, has_father, has_mother in sibships[
                ["FATHER_ID", "MOTHER_ID", "IID", "has_father", "has_mother"]
            ].itertuples(index=False, name=None):
            if not len(iids):
                raise ValueError("Each family must contain at least one offspring")
            self.max_sibs = max(self.max_sibs, len(iids))
            self.n_multi_sib += len(iids) > 1
            self.expected_ibd0 += (1 - 0.5 ** (len(iids) - 1)) ** 2
            if has_father != has_mother:
                self.n_parent_pairs += len(iids)
                self.n_single_parents += 1
            for iid in iids:
                sample_index = iid_to_bed_index.get(iid, -1)
                if not 0 <= sample_index < n_samples:
                    raise ValueError(f"Missing genotype for sibling {iid!r}")
                self.children.push_back(sample_index)
            self.offspring_offsets.push_back(self.children.size())
            sample_index = -1
            if has_father or has_mother:
                parent = father if has_father else mother
                sample_index = iid_to_bed_index.get(parent, -1)
                if not 0 <= sample_index < n_samples:
                    raise ValueError(f"Missing genotype for parent {parent!r}")
            self.parents.push_back(sample_index)
            # Match the preferred pair orientation from the original lookup.
            labels = [iid.decode() for iid in iids]
            for i in range(1, len(labels)):
                first = labels[i]
                for j in range(i):
                    second = labels[j]
                    if posterior_index is not None:
                        key = (iids[i], iids[j]) if iids[i] <= iids[j] else (iids[j], iids[i])
                        sample_index = posterior_index.get(key, -1)
                        if not -1 <= sample_index < self.n_posterior_rows:
                            raise ValueError("Sibling IBD posterior row is out of bounds")
                        self.posterior_rows.push_back(sample_index)
                    segment = ibd.get((first, second))
                    if segment is None:
                        segment = ibd.get((second, first), ())
                    if len(segment) % 3:
                        raise ValueError("Each IBD segment must contain a start, end and state")
                    sorted_pair = True
                    previous_end = -1
                    for offset in range(0, len(segment), 3):
                        segment_start = segment[offset]
                        segment_end = segment[offset + 1]
                        state = segment[offset + 2]
                        if segment_start < previous_end or segment_end < segment_start:
                            sorted_pair = False
                        previous_end = segment_end
                        self.segments.push_back(segment_start)
                        self.segments.push_back(segment_end)
                        self.segments.push_back(state)
                    self.sorted_segments.push_back(sorted_pair)
                    self.segment_offsets.push_back(self.segments.size())
            self.family_pair_offsets.push_back(self.sorted_segments.size())


@cython.wraparound(False)
@cython.boundscheck(False)
cdef int segment_state(ImputationFamilies families, size_t pair, int loc,
                       bint positions_sorted, size_t* cursor) noexcept nogil:
    cdef size_t end = families.segment_offsets[pair + 1]
    cdef size_t j
    if positions_sorted and families.sorted_segments[pair]:
        j = cursor[0]
        while j < end and families.segments[j + 1] < loc:
            j += 3
        cursor[0] = j
        if j < end and families.segments[j] <= loc:
            return families.segments[j + 2]
    else:
        for j in range(families.segment_offsets[pair], end, 3):
            if families.segments[j] <= loc <= families.segments[j + 1]:
                return families.segments[j + 2]
    return nan_integer


@cython.wraparound(False)
@cython.boundscheck(False)
cdef bint valid_dosages(signed char[:, :] dosages) noexcept nogil:
    cdef Py_ssize_t row, col
    cdef signed char call
    for row in range(dosages.shape[0]):
        for col in range(dosages.shape[1]):
            call = dosages[row, col]
            if call != nan_integer and not 0 <= call <= 2:
                return False
    return True


@cython.wraparound(False)
@cython.boundscheck(False)
def decode_phased_genotypes(const float[:, :, :] probabilities, int n_rows,
                            const Py_ssize_t[:] output_rows, const unsigned char[:] flipped):
    """Decode BGEN probabilities directly into compact hard calls in BED allele order."""
    cdef Py_ssize_t marker, begin, end, source, target, haplotype
    cdef Py_ssize_t n_samples = probabilities.shape[0]
    cdef Py_ssize_t n_markers = probabilities.shape[1]
    cdef signed char call
    cdef float threshold = 0.99
    if probabilities.shape[2] != 4 or n_samples != len(output_rows) or n_markers != len(flipped):
        raise ValueError("BGEN imputation requires phased biallelic SNPs")
    if n_samples and (np.min(output_rows) < 0 or np.max(output_rows) >= n_rows):
        raise ValueError("BGEN sample rows do not match the prepared genotype data")
    calls = np.full((n_rows, n_markers, 2), nan_integer, dtype=np.int8)
    cdef signed char[:, :, ::1] result = calls
    with nogil:
        for begin in range(0, n_samples, 128):
            end = min(begin + 128, n_samples)
            for marker in range(n_markers):
                for source in range(begin, end):
                    target = output_rows[source]
                    for haplotype in range(2):
                        call = nan_integer
                        if probabilities[source, marker, 2 * haplotype + 1] > threshold:
                            call = 0
                        elif probabilities[source, marker, 2 * haplotype] > threshold:
                            call = 1
                        if call >= 0 and flipped[marker]:
                            call = 1 - call
                        result[target, marker, haplotype] = call
    return calls


@cython.wraparound(False)
@cython.boundscheck(False)
def impute(sibships, iid_to_bed_index,  phased_gts, unphased_gts, ibd, pos, hdf5_output_dict, chromosome, freqs, output_address = None, threads = None, output_compression = None, output_compression_opts = None, half_window=50, ibd_threshold = 0.999, silent_progress=False, use_backup=False, *,
           prepared=None, return_genotypes=True, output_start=0, output_stop=None,
           parent_origin=None, parent_origin_index=None, parent_haplotype=None,
           ibd_posteriors=None, posterior_hard_fallback=False):
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
                'MAF_*' : info about the MAF estimator if MAF estimator is used.

        chromosome: str
            Name of the chromosome(s) that's going to be imputed. Only used for logging purposes.

        freqs: list[float]
            A two-dimensional array containing estimated fs for all individuals and SNPs respectively.
        
        output_address : str, optional
            If presented, the results would be written to this address in HDF5 format.
            Aside from all the key, value pairs inside hdf5_output_dict, the following are also written to the file.
                'imputed_par_gts' : imputed parental genotypes. It's the imputed missing parent if only one parent is missing and the imputed average of the both parents if both are missing.
                'pos' : the position of SNPs(in the order of appearance in genotypes)                
                'families' : family ids of the imputed parents(in the order of appearance in genotypes)
                'parental_status' : a numpy array where each row shows the family status of the family of the corresponding row in families. Columns are has_father, has_mother and, single_parent.
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
            Whether it should use backup imputation where there is no ibd information available. It's false by default.

        parent_origin : numpy.ndarray, optional
            Signed byte array with one row per inferred duo and one column per output SNP.
            Values 0 and 1 select the offspring haplotype inherited from the observed
            parent, and -1 selects the original imputation fallback.

        parent_origin_index : dict, optional
            Maps (observed parent genotype row, offspring genotype row) to an origin row.

        parent_haplotype : numpy.ndarray, optional
            Matching signed byte array for the observed parent's transmitted haplotype,
            with -1 for uncertain calls.

    Returns:
        tuple(list, numpy.array)
            The second element is imputed parental genotypes and the first element is family ids of the imputed parents(in the order of appearance in the first element).                        
    """
    logging.info("with chromosome " + str(chromosome)+": " + "imputing ...")
    if sibships.empty:
        logging.warning("with chromosome " + str(chromosome)+": " + "Error: No families to be imputed")
        return [], np.array()

    cdef int number_of_threads = 1 if threads is None else threads
    if number_of_threads < 1:
        raise ValueError("threads must be positive")
    if half_window < 0:
        raise ValueError("half_window must be non-negative")
    if unphased_gts.ndim != 2:
        raise ValueError("Genotypes must be two-dimensional")
    cdef int first = output_start
    cdef int stop = unphased_gts.shape[1] if output_stop is None else output_stop
    cdef int output_width = stop - first
    if not 0 <= first < stop <= unphased_gts.shape[1]:
        raise ValueError("Imputation requires a nonempty output slice")
    if freqs.shape != (unphased_gts.shape[0], output_width):
        raise ValueError("Genotypes and frequencies must have the same two-dimensional shape in the output slice")
    if (phased_gts is not None and
            phased_gts.shape != (*unphased_gts.shape, 2)):
        raise ValueError("Phased genotype dimensions must agree with unphased genotypes")
    if len(pos) != unphased_gts.shape[1] or not len(pos):
        raise ValueError("Provide exactly one position per SNP and at least one SNP")
    logging.info("with chromosome " + str(chromosome)+": " + "imputing data ...")
    cdef signed char[:, :] c_unphased_gts = unphased_gts
    cdef signed char[:, :, :] c_phased_gts = phased_gts
    cdef bint all_calls_valid
    with nogil:
        all_calls_valid = valid_dosages(c_unphased_gts)
    if not all_calls_valid:
        raise ValueError("Unphased genotypes must contain hard calls 0, 1 or 2, or the missing sentinel")
    cdef ImputationFamilies families
    if prepared is None:
        families = ImputationFamilies(sibships, iid_to_bed_index, ibd, unphased_gts.shape[0])
    else:
        families = prepared
        if families.n_samples != unphased_gts.shape[0]:
            raise ValueError("Prepared family indexes do not match the genotype sample count")
    cdef int max_sibs = families.max_sibs
    cdef int max_ibd_pairs = max(1, max_sibs*(max_sibs-1)//2)
    cdef int number_of_fams = families.parents.size()
    cdef const float[:, :, ::1] c_posteriors = ibd_posteriors
    cdef bint c_posterior_fallback = posterior_hard_fallback
    if ibd_posteriors is not None and (
            ibd_posteriors.shape != (families.n_posterior_rows, output_width, 3) or
            families.posterior_rows.size() != families.sorted_segments.size()):
        raise ValueError("IBD posteriors require a matching prepared family index and one column per output SNP")
    cdef const signed char[:, :] c_parent_origin = parent_origin
    cdef const signed char[:, :] c_parent_haplotype = parent_haplotype
    cdef int[:] origin_rows = None
    cdef Py_ssize_t fam, offspring, origin_row
    if parent_origin is not None:
        if phased_gts is None or parent_origin_index is None or parent_origin.shape[1] != output_width:
            raise ValueError("Parent origin requires phased genotypes, a duo index, and one column per output SNP")
        if parent_haplotype is not None and parent_haplotype.shape != parent_origin.shape:
            raise ValueError("Parent haplotypes must match the parent origin array shape")
        origin_rows = np.full(families.children.size(), -1, dtype=np.intc)
        for fam in range(number_of_fams):
            if families.single_parent[fam]:
                for offspring in range(families.offspring_offsets[fam], families.offspring_offsets[fam + 1]):
                    origin_row = parent_origin_index.get((families.parents[fam], families.children[offspring]), -1)
                    if not -1 <= origin_row < parent_origin.shape[0]:
                        raise ValueError("Parent origin duo index is out of bounds")
                    origin_rows[offspring] = origin_row
    elif parent_origin_index is not None or parent_haplotype is not None:
        raise ValueError("Parent origin array required with a duo index or parent haplotypes")
    cdef const double[:, :] c_freqs = freqs
    cdef int number_of_snps = c_unphased_gts.shape[1]
    cdef const unsigned char[:] single_parent = families.single_parent
    cdef double f, fvars
    cdef int p, q, nsibs
    cdef int i, j, loc, ibd_type, sib1_index, sib2_index, where, parent_index
    cdef size_t family_start, pair
    cdef size_t[:, :] segment_cursors = np.empty((number_of_threads, max_ibd_pairs), dtype=np.uintp)
    cdef bint positions_sorted = np.all(np.asarray(pos)[1:] >= np.asarray(pos)[:len(pos)-1])
    cdef cnp.ndarray[cnp.int64_t, ndim=1] c_pos = np.asarray(pos, dtype=np.int64)
    cdef int len_snp_ibd0 = 0
    cdef int len_snp_ibd1 = 0
    cdef int len_snp_ibd2 = 0
    cdef int[:,:,:] snp_ibd0 = np.ones([number_of_threads, max_ibd_pairs, 2], dtype=np.dtype("i"))
    cdef int[:,:,:] snp_ibd1 = np.ones([number_of_threads, max_ibd_pairs, 2], dtype=np.dtype("i"))
    cdef int[:,:,:] snp_ibd2 = np.ones([number_of_threads, max_ibd_pairs, 2], dtype=np.dtype("i"))
    cdef int[:, :] sibs_index = np.zeros((number_of_threads, max_sibs)).astype("i")
    cdef double[:, :] parent_genotype_prob = np.zeros((number_of_threads, 3))
    imputed = np.empty((number_of_fams, output_width))
    observed = np.empty((number_of_fams, output_width), dtype=np.float32)
    cdef double[:, ::1] imputed_par_gts = imputed
    cdef float[:, ::1] num_observed_parental_alleles = observed
    cdef double[:, ::1] robust = None
    cdef float[:, ::1] expected_observed = None
    cdef int[:, ::1] posterior_counts = None
    if ibd_posteriors is not None:
        robust = np.empty_like(imputed)
        expected_observed = np.empty_like(observed)
        posterior_counts = np.zeros((number_of_threads, output_width), dtype=np.intc)
    cdef int pair_posterior_row, hard_state, called_parent
    cdef cpair[double, cpair[float, bint]] posterior_result
    cdef long long[:] missing_output = np.zeros(number_of_threads, dtype=np.int64)
    cdef int snp, out_snp, this_thread, sib1_gene_isnan, sib2_gene_isnan, index
    byte_chromosome = chromosome.encode("ASCII")
    cdef char* chromosome_c = byte_chromosome
    cdef int mod = max(1, (number_of_fams + 1) // 100)
    cdef SniparProgress progress_counter
    #For hap_ibds, axes denote thread, individual pair, haplotypes and SNPs.
    cdef int [:,:,:,:] sib_hap_IBDs = None if phased_gts is None else np.empty([number_of_threads, max_ibd_pairs, 4, number_of_snps], dtype=np.dtype("i"))
    cdef int [:,:,:,:] parent_offspring_hap_IBDs = None if phased_gts is None else np.empty([number_of_threads, max_sibs, 4, number_of_snps], dtype=np.dtype("i"))
    cdef int half_window_c = half_window
    cdef float ibd_threshold_c = ibd_threshold
    cdef long[:, :] counter_ibd0 = np.zeros((number_of_threads, output_width), dtype="l")
    cdef bint is_backup
    cdef cpair[double,cpair[bint,float]] o_result
    cdef cpair[cpair[double,float],cpair[int,bint]] po_result
    cdef long[:, :] sib_backup_count = np.zeros((number_of_threads, output_width), dtype="l")
    cdef long[:, :] single_parent_backup_count = np.zeros((number_of_threads, output_width), dtype="l")
    cdef long[:, :] single_parent_mendelian_error_count = np.zeros((number_of_threads, output_width), dtype="l")
    cdef double[:, :] single_parent_fvars = np.zeros((number_of_threads, output_width))
    cdef int mendelian_error_count
    cdef bint c_silent_progress = silent_progress
    cdef bint c_use_backup = use_backup
    cdef int[:, :] sib_is_nan = np.zeros((number_of_threads, max_sibs)).astype("i")
    cdef bint has_non_nan_offspring = False
    cdef long[:, :] number_of_non_nan_offspring_per_snp = np.zeros((number_of_threads, output_width), dtype="l")
    cdef bint population_freq = c_freqs.strides[0] == 0
    cdef bint fast_parent_child
    cdef int child, known
    cdef double child_f, inferred_f, estimate
    cdef double[:, :, ::1] singleton = None
    if population_freq and c_phased_gts is None and families.n_single_parents:
        singleton = np.empty((output_width, 3, 3))
        with nogil:
            for out_snp in range(output_width):
                inferred_f = c_freqs[0, out_snp]
                if inferred_f > 1:
                    inferred_f = 1
                elif inferred_f < 0:
                    inferred_f = 0
                for known in range(3):
                    for child in range(3):
                        singleton[out_snp, known, child] = (
                            impute_parent_child(inferred_f, known, child)
                            if is_possible_child(child, known) else nan_float)
    logging.info("with chromosome " + str(chromosome)+": " + "using "+str(number_of_threads)+" threads")
    for index in prange(number_of_fams, nogil=True, num_threads=number_of_threads, schedule="static"):
        if c_silent_progress == 0:
            progress_counter.report(mod, chromosome_c, number_of_fams)
        this_thread = openmp.omp_get_thread_num()
        parent_index = families.parents[index]
        family_start = families.offspring_offsets[index]
        nsibs = families.offspring_offsets[index + 1] - family_start
        fast_parent_child = nsibs == 1 and c_phased_gts is None and single_parent[index]
        for i in range(nsibs):
            sibs_index[this_thread, i] = families.children[family_start + i]
        for i in range(families.family_pair_offsets[index + 1] - families.family_pair_offsets[index]):
            segment_cursors[this_thread, i] = families.segment_offsets[families.family_pair_offsets[index] + i]
        if c_phased_gts != None:
            # First fills hap_ibds
            for i in range(1, nsibs):
                for j in range(i):
                    where = get_hap_index(i, j)
                    fill_haplotype_ibd(c_phased_gts[sibs_index[this_thread, i],:,0], c_phased_gts[sibs_index[this_thread, j],:,0], number_of_snps, half_window_c, ibd_threshold_c, sib_hap_IBDs[this_thread, where, 0, :])
                    fill_haplotype_ibd(c_phased_gts[sibs_index[this_thread, i],:,0], c_phased_gts[sibs_index[this_thread, j],:,1], number_of_snps, half_window_c, ibd_threshold_c, sib_hap_IBDs[this_thread, where, 1, :])
                    fill_haplotype_ibd(c_phased_gts[sibs_index[this_thread, i],:,1], c_phased_gts[sibs_index[this_thread, j],:,0], number_of_snps, half_window_c, ibd_threshold_c, sib_hap_IBDs[this_thread, where, 2, :])
                    fill_haplotype_ibd(c_phased_gts[sibs_index[this_thread, i],:,1], c_phased_gts[sibs_index[this_thread, j],:,1], number_of_snps, half_window_c, ibd_threshold_c, sib_hap_IBDs[this_thread, where, 3, :])

            if single_parent[index]:
                for i in range(0, nsibs):
                    fill_haplotype_ibd(c_phased_gts[parent_index,:,0], c_phased_gts[sibs_index[this_thread, i],:,0], number_of_snps, half_window_c, ibd_threshold_c, parent_offspring_hap_IBDs[this_thread, i, 0, :])
                    fill_haplotype_ibd(c_phased_gts[parent_index,:,0], c_phased_gts[sibs_index[this_thread, i],:,1], number_of_snps, half_window_c, ibd_threshold_c, parent_offspring_hap_IBDs[this_thread, i, 1, :])
                    fill_haplotype_ibd(c_phased_gts[parent_index,:,1], c_phased_gts[sibs_index[this_thread, i],:,0], number_of_snps, half_window_c, ibd_threshold_c, parent_offspring_hap_IBDs[this_thread, i, 2, :])
                    fill_haplotype_ibd(c_phased_gts[parent_index,:,1], c_phased_gts[sibs_index[this_thread, i],:,1], number_of_snps, half_window_c, ibd_threshold_c, parent_offspring_hap_IBDs[this_thread, i, 3, :])
        for snp in range(first, stop):
            out_snp = snp - first
            if fast_parent_child:
                child = c_unphased_gts[sibs_index[this_thread, 0], snp]
                known = c_unphased_gts[parent_index, snp]
                if child != nan_integer and known != nan_integer:
                    number_of_non_nan_offspring_per_snp[this_thread, out_snp] += 1
                    child_f = c_freqs[sibs_index[this_thread, 0], out_snp]
                    single_parent_fvars[this_thread, out_snp] += child_f * (1 - child_f)
                    if is_possible_child(child, known):
                        if population_freq:
                            estimate = singleton[out_snp, known, child]
                        else:
                            inferred_f = 2 * child_f - c_freqs[parent_index, out_snp]
                            if inferred_f > 1:
                                inferred_f = 1
                            elif inferred_f < 0:
                                inferred_f = 0
                            estimate = impute_parent_child(inferred_f, known, child)
                        num_observed_parental_alleles[index, out_snp] = 3
                    else:
                        estimate = nan_float
                        num_observed_parental_alleles[index, out_snp] = nan_float
                        single_parent_mendelian_error_count[this_thread, out_snp] += 1
                    imputed_par_gts[index, out_snp] = estimate
                    if c_posteriors is not None:
                        robust[index, out_snp] = estimate
                        expected_observed[index, out_snp] = num_observed_parental_alleles[index, out_snp]
                    missing_output[this_thread] += isnan(estimate)
                    continue
            len_snp_ibd0 = 0
            len_snp_ibd1 = 0
            len_snp_ibd2 = 0
            loc = c_pos[snp]
            pair_posterior_row = -1
            hard_state = -1
            called_parent = -1
            has_non_nan_offspring = False
            for i in range(nsibs):
                sib1_index = sibs_index[this_thread, i]
                sib_is_nan[this_thread, i] = (c_unphased_gts[sib1_index, snp] == nan_integer)
                if not sib_is_nan[this_thread, i]:
                    has_non_nan_offspring = True
            if has_non_nan_offspring:
                number_of_non_nan_offspring_per_snp[this_thread, out_snp] += 1
                if nsibs > 1:
                    for i in range(1, nsibs):
                        for j in range(i):
                            sib1_gene_isnan = sib_is_nan[this_thread, i]
                            sib2_gene_isnan = sib_is_nan[this_thread, j]
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

                            else:
                                where = get_hap_index(i, j)
                                pair = families.family_pair_offsets[index] + where
                                ibd_type = segment_state(families, pair, loc, positions_sorted,
                                                         &segment_cursors[this_thread, where])
                                if c_posteriors is not None:
                                    if c_posterior_fallback and ibd_type < 0 and families.posterior_rows[pair] >= 0:
                                        ibd_type = posterior_mode(c_posteriors, families.posterior_rows[pair], out_snp)
                                    if nsibs == 2:
                                        pair_posterior_row = families.posterior_rows[pair]
                                        hard_state = ibd_type
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
                        counter_ibd0[this_thread, out_snp] = counter_ibd0[this_thread, out_snp]+1
                else :
                    sib1_index = sibs_index[this_thread, 0]
                    if not (c_unphased_gts[sib1_index, snp] == nan_integer):
                        snp_ibd2[this_thread, len_snp_ibd2,0] = 0
                        snp_ibd2[this_thread, len_snp_ibd2,1] = 0
                        len_snp_ibd2 = len_snp_ibd2 + 1
                f = 0.
                fvars = 0.
                for p in range(nsibs):
                    q = sibs_index[this_thread, p]
                    fvars = fvars + c_freqs[q, out_snp]*(1-c_freqs[q, out_snp])
                    f = f + c_freqs[q, out_snp]
                f = f/nsibs
                parent_genotype_prob[this_thread, 0] = (1-f)**2
                parent_genotype_prob[this_thread, 1] = 2*f*(1-f)
                parent_genotype_prob[this_thread, 2] = f**2
                if single_parent[index] and c_unphased_gts[parent_index, snp] != nan_integer:
                    called_parent = parent_index
                    #2 time MAF of offspring is MAF of sum of the parents. That minus the existing parent results in MAF of the missing parent.
                    f = 2*f-c_freqs[parent_index, out_snp]
                    if f>1.:
                        f=1.
                    elif f<0.:
                        f=0.
                    po_result = impute_snp_from_parent_offsprings(snp,
                                                                                    parent_index,
                                                                                    sibs_index[this_thread, :],
                                                                                    nsibs,
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
                                                                                    out_snp,
                                                                                    family_start,
                                                                                    origin_rows,
                                                                                    c_parent_origin,
                                                                                    c_parent_haplotype,
                                                                                    )
                    imputed_par_gts[index, out_snp] = po_result.first.first
                    num_observed_parental_alleles[index, out_snp] = po_result.first.second
                    mendelian_error_count = po_result.second.first
                    single_parent_mendelian_error_count[this_thread, out_snp] += mendelian_error_count
                    single_parent_fvars[this_thread, out_snp] += fvars
                    is_backup = po_result.second.second
                    single_parent_backup_count[this_thread, out_snp] += is_backup
                else:
                    o_result = impute_snp_from_offsprings(snp,
                                                                            sibs_index[this_thread, :],
                                                                            nsibs,
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
                    imputed_par_gts[index, out_snp] = o_result.first
                    is_backup = o_result.second.first
                    num_observed_parental_alleles[index, out_snp] = o_result.second.second
                    sib_backup_count[this_thread, out_snp] += is_backup
            else:
                imputed_par_gts[index, out_snp] = nan_float
                num_observed_parental_alleles[index, out_snp] = nan_float
            if c_posteriors is not None:
                estimate = imputed_par_gts[index, out_snp]
                robust[index, out_snp] = estimate
                expected_observed[index, out_snp] = num_observed_parental_alleles[index, out_snp]
                if pair_posterior_row >= 0:
                    posterior_result = impute_pair_posterior(
                        snp, out_snp, pair_posterior_row, c_posteriors, called_parent,
                        sibs_index[this_thread, :], snp_ibd0[this_thread, :, :], snp_ibd1[this_thread, :, :],
                        snp_ibd2[this_thread, :, :], f, parent_genotype_prob[this_thread, :], c_phased_gts,
                        c_unphased_gts, sib_hap_IBDs[this_thread, :, :, :], parent_offspring_hap_IBDs[this_thread, :, :, :],
                        family_start, origin_rows, c_parent_origin, c_parent_haplotype, hard_state,
                        estimate, expected_observed[index, out_snp])
                    imputed_par_gts[index, out_snp] = posterior_result.first
                    expected_observed[index, out_snp] = posterior_result.second.first
                    posterior_counts[this_thread, out_snp] += posterior_result.second.second
            missing_output[this_thread] += isnan(imputed_par_gts[index, out_snp])
    number_of_po_pairs = families.n_parent_pairs
    mendelian_errors = np.sum(single_parent_mendelian_error_count, axis=0)
    fvar_totals = np.sum(single_parent_fvars, axis=0)
    mendelian_error_ratio = mendelian_errors / max(1, number_of_po_pairs)
    with np.errstate(invalid="ignore", divide="ignore"):
        estimated_genotyping_error = mendelian_errors / fvar_totals
    multi_sib_fams = families.n_multi_sib
    single_parent_fams = families.n_single_parents
    no_parent_fams = number_of_fams - single_parent_fams
    multi_sib_fams_ratio = multi_sib_fams / number_of_fams
    sib_ratio_backup = np.sum(sib_backup_count, axis=0) / max(1, no_parent_fams)
    parent_ratio_backup = np.sum(single_parent_backup_count, axis=0) / max(1, single_parent_fams)
    ibd0_counts = np.sum(counter_ibd0, axis=0)
    observed_counts = np.sum(number_of_non_nan_offspring_per_snp, axis=0)
    ratio_ibd0 = ibd0_counts / np.maximum(1, observed_counts)
    total_ratio_ibd0 = ibd0_counts.sum() / max(1, observed_counts.sum())
    expected_total_ratio_ibd0 = families.expected_ibd0 / number_of_fams
    logging.info(f"with chromosome {chromosome} :total number of fams is {number_of_fams}")
    logging.info(f"with chromosome {chromosome} :more than one offspring is genotyped for {multi_sib_fams} which is {multi_sib_fams_ratio:.2f} of total fams")
    logging.info(f"with chromosome {chromosome} :IBD0 state observed for {total_ratio_ibd0*100:.2f}%")
    logging.info(f"with chromosome {chromosome} :Expected to observe IBD0 state for {expected_total_ratio_ibd0*100:.2f}%")
    if use_backup:
        logging.info(f"with chromosome {chromosome} :ratio of backup imputation among sibs is {np.mean(sib_ratio_backup):.2f}")
        logging.info(f"with chromosome {chromosome} :ratio of backup imputation among parent-offsprings is {np.mean(parent_ratio_backup):.2f}")
    if total_ratio_ibd0>0.5:
        logging.warning("with chromosome " + str(chromosome)+": ibd0 ratio is too high")
    output_nan_count = np.sum(missing_output)
    nan_ratio = output_nan_count/imputed.size
    logging.info("with chromosome " + str(chromosome)+f": total number of nans: {output_nan_count}, ratio of nan snps to all snps is {nan_ratio:.2f}")
    if nan_ratio > 0.01:
        logging.warning("Too much inconsistencies in the data")
    
    hdf5_output_dict["mendelian_error_ratio"] = mendelian_error_ratio
    hdf5_output_dict["estimated_genotyping_error"] = estimated_genotyping_error
    hdf5_output_dict["ratio_ibd0"] = ratio_ibd0
    hdf5_output_dict["sib_ratio_backup"] = sib_ratio_backup
    hdf5_output_dict["parent_ratio_backup"] = parent_ratio_backup
    hdf5_output_dict['families'] = families.family_ids
    hdf5_output_dict['parental_status'] = families.parental_status
    hdf5_output_dict['pos'] = pos[first:stop]
    hdf5_output_dict['imputed_par_gts'] = imputed
    hdf5_output_dict['num_observed_parental_alleles'] = observed
    if ibd_posteriors is not None:
        hdf5_output_dict['robust_imputed_par_gts'] = np.asarray(robust)
        hdf5_output_dict['num_expected_observed_parental_alleles'] = np.asarray(expected_observed)
        hdf5_output_dict['posterior_num_sibships'] = np.sum(posterior_counts, axis=0)
    if output_address is not None:
        logging.info("with chromosome " + str(chromosome)+": " + "Writing the results as a hdf5 file to "+output_address + ".hdf5")
        with h5py.File(output_address+".hdf5",'w') as file:                        
            for key, val in hdf5_output_dict.items():
                if key in {'imputed_par_gts', 'num_observed_parental_alleles',
                           'robust_imputed_par_gts', 'num_expected_observed_parental_alleles'}:
                    file.create_dataset(key, val.shape, dtype = 'float16', chunks = True, compression = output_compression, compression_opts=output_compression_opts, data = val)
                else:
                    file[key] = val
    if return_genotypes:
        return sibships["FID"].values.tolist(), imputed
    return None
