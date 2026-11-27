import numpy as np
import numpy.ma as ma
from numba import njit
from snipar.utilities import make_id_dict
from snipar.types import SparseGRMRepr, Ids
from scipy.sparse import csc_matrix, tril
from scipy.sparse.linalg import splu
# import logging

# logger = logging.getLogger(__name__)
class gtarray(object):
    """Define a genotype or PGS array that stores individual IDs, family IDs, and SNP information.

    Args:
        garray : :class:`~numpy:numpy.array`
            2 or 3 dimensional numpy array of genotypes/PGS values. First dimension is individuals. For a 2 dimensional array, the second dimension is SNPs or PGS values.
            For a 3 dimensional array, the second dimension indexes the individual and his/her relatives' genotypes (for example: proband, paternal, and maternal); and
            the third dimension is the SNPs.
        ids : :class:`~numpy:numpy.array`
            vector of individual IDs of same length as first dimension of garray
        sid : :class:`~numpy:numpy.array`
            vector of SNP ids, equal in length, L, to last dimension of array
        alleles : :class:`~numpy:numpy.array`
            [L x 2] matrix of ref and alt alleles for the SNPs. L must match size of sid
        pos : :class:`~numpy:numpy.array`
            vector of SNP positions; must match size of sid
        chrom : :class:`~numpy:numpy.array`
            vector of SNP chromosomes; must match size of sid
        map : :class:`~numpy:numpy.array`
            vector of SNP chromosomes; must match size of sid
        fams : :class:`~numpy:numpy.array`
            vector of family IDs; must match size of ids
        par_status : :class:`~numpy:numpy.array'
             [N x 2] numpy matrix that records whether parents have observed or imputed genotypes/PGS, where N matches size of ids.
             The first column is for the father of that individual; the second column is for the mother of that individual.
             If the parent is neither observed nor imputed, the value is -1; if observed, 0; and if imputed, 1.
        num_obs_par_al : :class:`~numpy:numpy.array'
             2 dimensioanal numpy array, with each entry holding the number of observed parental alleles for each individual on each locus.

    Returns:
        G : :class:`snipar.gtarray`

    """
    def __init__(self, garray, ids, sid=None, alleles=None, pos=None, chrom=None, map=None, error_probs=None, fams=None, par_status=None, num_obs_par_al=None, ped=None, standard_f=None):
        if type(garray) == np.ndarray or type(garray) == np.ma.core.MaskedArray:
            if type(garray) == np.ndarray:
                mask = np.isnan(garray)
                self.has_NAs = bool(mask.any())
                self.gts = ma.array(garray, mask=mask if self.has_NAs else ma.nomask)
            else:
                self.gts = garray
                self.has_NAs = bool(np.any(garray.mask))
            self.shape = garray.shape
            self.ndim = garray.ndim
            self.dtype = garray.dtype
            self.freqs = None
        else:
            raise ValueError('Genotypes must be a numpy ndarray')
        if self.ndim not in (2, 3):
            raise ValueError('Genotypes must be a 2- or 3-dimensional array')
        if garray.shape[0] == ids.shape[0]:
            self.ids = ids
            self.id_dict = make_id_dict(ids)
        else:
            raise ValueError('Shape of genotypes and ids does not match')
        self.snp_index = self.ndim - 1
        self.sid = sid
        self.sid_dict = None
        if sid is not None:
            if sid.shape != (garray.shape[-1],):
                raise ValueError('Shape of SNP ids (sid) does not match shape of genotype array')
            self.sid_dict = make_id_dict(sid)
        if alleles is not None:
            if self.sid is not None:
                if alleles.shape[0] == self.sid.shape[0]:
                    self.alleles = alleles
                else:
                    raise ValueError('Size of alleles does not match size of genotypes')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.alleles = None
        if pos is not None:
            if self.sid is not None:
                if pos.shape[0] == self.sid.shape[0]:
                    self.pos = pos
                else:
                    raise ValueError('Size of position vector does not match size of genotypes')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.pos = None
        if chrom is not None:
            if self.sid is not None:
                if chrom.shape[0] == self.sid.shape[0]:
                    self.chrom = chrom
                else:
                    raise ValueError('Size of map does not match number of SNPs')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.chrom = None
        if map is not None:
            if self.sid is not None:
                if map.shape[0] == self.sid.shape[0]:
                    self.map = map
                else:
                    raise ValueError('Size of map does not match number of SNPs')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.map = None
        if error_probs is not None:
            if self.sid is not None:
                if error_probs.shape[0] == self.sid.shape[0]:
                    self.error_probs = error_probs
                else:
                    raise ValueError('Size of map does not match number of SNPs')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.error_probs = None
        if fams is not None:
            if fams.shape[0] == ids.shape[0] and fams.ndim==1:
                self.fams = fams
            else:
                raise ValueError('Fams not of same length as IDs')
        else:
            self.fams = None

        if par_status is not None:
            if par_status.shape == (ids.shape[0], 2):
                self.par_status = par_status
            else:
                raise ValueError('Incompatible par status array')
        else:
            self.par_status = None
        
        if ped is not None:
            self.ped = ped
        else:
            self.ped = None
        
        if num_obs_par_al is not None:
            if num_obs_par_al.shape == (self.shape[0], self.shape[-1]):
                self.num_obs_par_al = num_obs_par_al
            else:
                raise ValueError('Incompatible par status array')
        else:
            self.num_obs_par_al = None
        
        if standard_f is not None:
            if self.sid is not None and standard_f.shape == self.sid.shape:
                self.standard_f = standard_f
            else:
                raise ValueError('Size of standard_f does not match number of SNPs')
        else:
            self.standard_f = None

        self.mean_normalised = False

        self.info = None
        self.variances = None

    def compute_freqs(self):
        """
        Computes the frequencies of the SNPs. Stored in self.freqs.
        """
        proband = self.gts if self.ndim == 2 else self.gts[:, 0, :]
        self.freqs = ma.filled(ma.mean(proband, axis=0, dtype=np.float64) / 2.0, np.nan)

    def filter(self, filter_pass):
        if self.freqs is not None:
            self.freqs = self.freqs[filter_pass]
        self.gts = self.gts[..., filter_pass]
        self.shape = self.gts.shape
        if self.has_NAs:
            self.has_NAs = bool(np.any(self.gts.mask))
            if not self.has_NAs:
                self.gts.shrink_mask()
        if self.info is not None:
            self.info = self.info[filter_pass]
            self.variances = self.variances[filter_pass]
        if self.sid is not None:
            self.sid = self.sid[filter_pass]
            self.sid_dict = make_id_dict(self.sid)
        if self.pos is not None:
            self.pos = self.pos[filter_pass]
        if self.alleles is not None:
            self.alleles = self.alleles[filter_pass]
        if self.chrom is not None:
            self.chrom = self.chrom[filter_pass]
        if self.map is not None:
            self.map = self.map[filter_pass]
        if self.error_probs is not None:
            self.error_probs = self.error_probs[filter_pass]
        if self.num_obs_par_al is not None:
            self.num_obs_par_al = self.num_obs_par_al[:, filter_pass]
        if self.standard_f is not None:
            self.standard_f = self.standard_f[filter_pass]

    def filter_maf(self, min_maf = 0.01, verbose=False):
        """
        Filter SNPs based on having minor allele frequency (MAF) greater than min_maf, and have % missing observations less than max_missing.
        """
        if self.freqs is None:
            self.compute_freqs()
        freqs_pass = np.logical_and(self.freqs > min_maf, self.freqs < (1 - min_maf))
        if verbose:
            print(str(self.freqs.shape[0] - np.sum(freqs_pass)) + ' SNPs with MAF<' + str(min_maf))
        if not np.all(freqs_pass):
            self.filter(freqs_pass)
        return freqs_pass

    def filter_missingness(self, max_missing = 5, verbose=False):
        mask = self.gts.mask
        if np.ndim(mask) == 0:
            missingness_pass = np.full(self.shape[-1], 100 * mask <= max_missing)
        else:
            missingness = mask.mean(axis=tuple(range(self.snp_index)))
            missingness_pass = 100 * missingness <= max_missing
        if verbose:
            print(str(self.shape[-1] - np.sum(missingness_pass)) + ' SNPs with missingness >' + str(max_missing) + '%')
        if not np.all(missingness_pass):
            self.filter(missingness_pass)
        return missingness_pass

    def compute_info(self):
        if self.freqs is None:
            self.compute_freqs()
        if self.ndim == 2:
            self.variances = np.var(self.gts, axis=0)
        elif self.ndim==3:
            self.variances = np.var(self.gts[:,0,:], axis=0)
        denominator = 2.0 * self.freqs * (1 - self.freqs)
        self.info = np.divide(ma.filled(self.variances, np.nan), denominator,
                              out=np.full(denominator.shape, np.nan), where=denominator > 0)

    def filter_info(self, min_info = 0.99, verbose=False):
        if self.info is None:
            self.compute_info()
        info_pass = self.info > min_info
        if verbose:
            print(str(self.info.shape[0] - np.sum(info_pass)) + ' SNPs with INFO <' + str(min_info))
        self.filter(info_pass)

    def filter_ids(self,keep_ids, verbose=False):
        """
        Keep only individuals with ids given by keep_ids
        """
        keep_ids = np.asarray(keep_ids)
        in_ids = np.array([x in self.id_dict for x in keep_ids], dtype=bool)
        n_filtered = np.sum(in_ids)
        if n_filtered==0:
            raise(ValueError('No individuals would be left after filtering'))
        else:
            if verbose:
                print('After filtering, '+str(n_filtered)+' individuals remain')
            indices = np.array([self.id_dict[x] for x in keep_ids[in_ids]])
            self.gts = self.gts[indices]
            self.ids = self.ids[indices]
            self.id_dict = make_id_dict(self.ids)
            self.shape = self.gts.shape
            if self.fams is not None:
                self.fams = self.fams[indices]
            if self.par_status is not None:
                self.par_status = self.par_status[indices, :]
            if self.num_obs_par_al is not None:
                self.num_obs_par_al = self.num_obs_par_al[indices, :]
            for field in ('complete_trios_inds', 'sibs_inds'):
                if hasattr(self, field):
                    setattr(self, field, np.flatnonzero(np.isin(indices, getattr(self, field))))
            self.freqs = self.info = self.variances = None
            if self.has_NAs:
                self.has_NAs = bool(np.any(self.gts.mask))
                if not self.has_NAs:
                    self.gts.shrink_mask()
            self.mean_normalised = False

    def mean_normalise(self):
        """
        This normalises the SNPs/PGS columns to have mean-zero.
        """
        if not self.mean_normalised:
            self.gts = self.gts - ma.mean(self.gts, axis=0)
            self.dtype = self.gts.dtype
            self.mean_normalised = True

    def scale(self):
        """
        This normalises the SNPs/PGS columns to have variance 1.
        """
        self.gts = self.gts / ma.std(self.gts, axis=0)
        self.dtype = self.gts.dtype
        self.has_NAs = bool(np.any(self.gts.mask))

    def fill_NAs(self):
        """
        This normalises the SNP columns to have mean-zero, then fills in NA values with zero.
        """
        if not self.mean_normalised:
            self.mean_normalise()
        mask = self.gts.mask
        if np.ndim(mask) == 0:
            self.has_NAs = False
            return np.zeros(self.shape[1:], dtype=int)
        NAs = mask.sum(axis=0)
        if NAs.any():
            self.gts[mask] = 0
        self.gts.mask = ma.nomask
        self.has_NAs = False
        return NAs

    def add(self,garray):
        """
        Adds another gtarray of the same dimension to this array and returns the sum. It matches IDs before summing.
        """
        if type(garray)==gtarray:
            pass
        else:
            raise ValueError('Must add to another gtarray')

        if not self.gts.ndim == garray.gts.ndim:
            raise ValueError('Arrays must have same number of dimensions')

        if self.gts.ndim == 2:
            if not self.gts.shape[1] == garray.gts.shape[1]:
                raise ValueError('Arrays must have same dimensions (apart from first)')

        if self.gts.ndim == 3:
            if not self.gts.shape[1:3] == garray.gts.shape[1:3]:
                raise ValueError('Arrays must have same dimensions (apart from first)')

        # Match IDs
        common_ids = [iid for iid in self.id_dict if iid in garray.id_dict]
        if len(common_ids) == 0:
            raise ValueError('No IDs in common')
        self_index = np.array([self.id_dict[x] for x in common_ids])
        other_index = np.array([garray.id_dict[x] for x in common_ids])

        # Out
        if self.ids.ndim == 1:
            ids_out = self.ids[self_index]
        else:
            ids_out = self.ids[self_index, :]

        if self.gts.ndim ==2:
            add_gts = self.gts[self_index, :]+garray.gts[other_index, :]
        else:
            add_gts = self.gts[self_index, :, :] + garray.gts[other_index, :, :]
        if self.ped is None:
            ped = garray.ped
        else:
            ped = self.ped
        fams = None if self.fams is None else self.fams[self_index]
        par_status = None if self.par_status is None else self.par_status[self_index]
        return gtarray(add_gts, ids_out, self.sid, alleles=self.alleles,
                       fams=fams, par_status=par_status, ped=ped)

    def diagonalise(self,inv_root):
        """
        This will transform the genotype array based on the inverse square root of the phenotypic covariance matrix
        from the family based linear mixed model.
        """
        if self.fams is None:
            raise(ValueError('Family labels needed for diagonalization'))
        if not self.mean_normalised:
            self.mean_normalise()
        groups = {}
        for i, fam in enumerate(self.fams):
            groups.setdefault(fam, []).append(i)
        fam_indices = {fam: np.asarray(idx) for fam, idx in groups.items()}
        for fam in fam_indices:
            famsize = fam_indices[fam].shape[0]
            if self.ndim == 2:
                if famsize == 1:
                    self.gts[fam_indices[fam], :] = inv_root[1]*self.gts[fam_indices[fam],:]
                else:
                    self.gts[fam_indices[fam],:] = inv_root[famsize].dot(self.gts[fam_indices[fam],:])
            elif self.ndim == 3:
                if famsize == 1:
                    self.gts[fam_indices[fam], : , :] = inv_root[1]*self.gts[fam_indices[fam], : , :]
                else:
                    for j in range(self.shape[1]):
                        self.gts[fam_indices[fam],j, :] = inv_root[famsize].dot(self.gts[fam_indices[fam],j, :])
        self.fam_indices = fam_indices


@njit
def _is_whole(n):
    return n % 1 == 0


@njit
def _lin_imp(freq, g, gp):
    res = (4 * freq + 2 * g - gp) / 3.
    return min(2, max(res, 0))


@njit
def _po_imp(freq, g, gp):
    if g == 0:
        return freq
    elif g == 1:
        if gp == 0:
            return 1 + freq
        elif gp == 1:
            return 0.5 + freq
        else:  # gp == 2
            return freq
    else:  # g == 2
        return 1 + freq


@njit
def _impute_unrel_pat_gts(gts: np.ndarray, freqs: np.ndarray, parsum: bool) -> np.ndarray:
    if gts.shape[2] == 0:
        return gts
    for i in range(gts.shape[0]):
        if gts[i, 1, 0] != -1:
            continue
        for j in range(gts.shape[2]):
            estimate = 0.5 * gts[i, 0, j] + freqs[j]
            gts[i, 1, j] = 2 * estimate if parsum else estimate
            if not parsum:
                gts[i, 2, j] = estimate
    return gts


def _impute_unrel_pat_gts_cond_gau(G: gtarray, ped: np.ndarray,
                                   grm: SparseGRMRepr, unrelated_inds: Ids, parsum: bool) -> np.ndarray:
    freqs = G.freqs
    id_dict = G.id_dict
    gts = G.gts.astype(np.result_type(G.dtype, np.float32), copy=False).filled(np.nan)
    if not len(unrelated_inds):
        return gts
    unrelated_inds_dict = make_id_dict(unrelated_inds)
    R12_data, R12_row, R12_col = [], [], []
    for i in range(len(grm[0])):
        row = grm[1][i]
        col = grm[2][i]
        row_in = row in unrelated_inds_dict
        col_in = col in unrelated_inds_dict
        if not row_in and not col_in:
            continue
        elif row_in and col_in:
            if row == col:
                continue
            R12_row.append(unrelated_inds_dict[row])
            R12_col.append(col)
            R12_data.append(grm[0][i])
            R12_row.append(unrelated_inds_dict[col])
            R12_col.append(row)
            R12_data.append(grm[0][i])
        else:
            if row_in and not col_in:
                ind1, ind2 = row, col
            else:
                ind2, ind1 = row, col
            R12_row.append(unrelated_inds_dict[ind1])
            R12_col.append(ind2)
            R12_data.append(grm[0][i])

    R12_col += [idx for idx in unrelated_inds]
    R12_row += [unrelated_inds_dict[idx] for idx in unrelated_inds]
    R12_data += [0.5 for k in range(len(unrelated_inds))]
    if ped is not None:
        for iid, father, mother in ped[:, 1:4]:
            child_index = id_dict.get(iid)
            if child_index is None:
                continue
            for parent in (father, mother):
                parent_index = id_dict.get(parent)
                if parent_index in unrelated_inds_dict:
                    R12_row.append(unrelated_inds_dict[parent_index])
                    R12_col.append(child_index)
                    R12_data.append(0.25)

    R12_row = np.array(R12_row, dtype='uint32')
    R12_col = np.array(R12_col, dtype='uint32')
    R12_data = np.array(R12_data)
    R12 = csc_matrix((R12_data, (R12_row, R12_col)), shape=(len(unrelated_inds), G.shape[0]))
    R22 = csc_matrix((grm[0], (grm[1], grm[2])), shape=(G.shape[0], G.shape[0]))
    R22_triu: csc_matrix = tril(R22, k=-1, format='csc')
    R22 = R22 + R22_triu.T
  
    means = 2 * freqs
    lu = splu(R22)
    for start in range(0, G.shape[-1], 256):
        block = slice(start, start + 256)
        proband = gts[:, 0, block] - means[block]
        proband[np.isnan(proband)] = 0
        estimate = np.clip(means[block] + R12 @ lu.solve(proband.astype(float, copy=False)), 0, 2)
        gts[unrelated_inds, 1, block] = 2 * estimate if parsum else estimate
        if not parsum:
            gts[unrelated_inds, 2, block] = estimate
    return gts


@njit
def _impute_missing(gts: np.ndarray, freqs: np.ndarray) -> np.ndarray:
    parents_separate = gts.shape[1] == 3
    for j in range(gts.shape[0]):
        for i in range(gts.shape[2]):
            child = gts[j, 0, i]
            first = gts[j, 1, i]
            no_child, no_first = np.isnan(child), np.isnan(first)
            if not parents_separate:
                if no_child:
                    gts[j, 0, i] = 2 * freqs[i] if no_first else first / 2
                if no_first:
                    gts[j, 1, i] = 4 * freqs[i] if no_child else child + 2 * freqs[i]
                continue
            second = gts[j, 2, i]
            no_second = np.isnan(second)
            if not (no_child or no_first or no_second):
                continue
            freq = freqs[i]
            if no_child:
                if no_first and no_second:
                    gts[j, :, i] = 2 * freq
                elif no_first:
                    gts[j, 0, i] = second / 2 + freq
                    gts[j, 1, i] = 2 * freq
                elif no_second:
                    gts[j, 0, i] = first / 2 + freq
                    gts[j, 2, i] = 2 * freq
                else:
                    gts[j, 0, i] = (first + second) / 2
            elif no_first and no_second:
                gts[j, 1:3, i] = child / 2 + freq
            elif no_first:
                if _is_whole(child) and _is_whole(second):
                    gts[j, 1, i] = _po_imp(freq, child, second)
                else:
                    gts[j, 1, i] = _lin_imp(freq, child, second)
            else:
                if _is_whole(child) and _is_whole(first):
                    gts[j, 2, i] = _po_imp(freq, child, first)
                else:
                    gts[j, 2, i] = _lin_imp(freq, child, first)
    return gts


def impute_unrel_par_gts(G: gtarray, sib: bool = False, parsum: bool = False, ped: np.ndarray = None, grm: SparseGRMRepr = None, unrelated_inds: Ids = None, true_par=None, par_ids=None) -> gtarray:
    """Impute parental genotypes of unrelated individuals.

    Args:
        G (gtarray): gtarray object holding genetic data.
        sib (bool, optional): Whether sib effect is modelled. Defaults to False.
        parsum (bool, optional): Whether sum of parental genotypes is modelled. Defaults to False.

    Returns:
        gtarray: gtarray object with imputed parental genotypes.
    """ 
    if sib:
        # logger.warning(
        #     'No need to impute parental genotypes of unrelated individuals if sib effect is modelled.'
        # )
        # print('WARNING: No need to impute parental genotypes of unrelated individuals if sib effect is modelled.')
        return G
    if G.ndim != 3 or G.shape[1] != (2 if parsum else 3):
        raise ValueError('Genotype matrix does not match the parental genotype design.')
    if G.freqs is None:
        G.compute_freqs()
    if grm is not None and unrelated_inds is not None:
        # logger.info(
        #     'Imputing parental geotypes using conditional gaussian.'
        # )
        # print('Imputing parental geotypes using conditional gaussian')
        gts = _impute_unrel_pat_gts_cond_gau(G, parsum=parsum, ped=ped, grm=grm, unrelated_inds=unrelated_inds)
    elif grm is None and unrelated_inds is None:
        gts = G.gts.astype(np.result_type(G.dtype, np.float32), copy=False).filled(np.nan)
        gts = _impute_unrel_pat_gts(gts, G.freqs, parsum=parsum)
    else:
        raise TypeError('One of grm and unrelated_inds is None.')
    mask = np.isnan(gts)
    G.has_NAs = bool(mask.any())
    G.gts = ma.array(gts, mask=mask if G.has_NAs else ma.nomask)
    G.dtype = G.gts.dtype
    # logger.info(
    #     'Done imputing parental geotypes of unrelated individuals.'
    # )
    # print('Done imputing parental geotypes of unrelated individuals.')
    return G


def impute_missing(G: gtarray) -> gtarray:
    """Imputing missing entries of the genotype design matrix.

    Args:
        G (gtarray): genotype design matrix.

    Returns:
        gtarray: genotype design matrix with NAs imputed.
    """
    if G.freqs is None:
        G.compute_freqs()
    if G.ndim != 3 or G.shape[1] not in (2, 3):
        raise ValueError('Missing-genotype imputation requires a proband and two parents or a parental sum')
    gts = _impute_missing(G.gts.astype(np.result_type(G.dtype, np.float32), copy=False).filled(np.nan), G.freqs)
    if np.isnan(gts).any():
        raise ValueError('Genotypes remain missing after imputation')
    G.gts = ma.array(gts)
    G.dtype = G.gts.dtype
    G.has_NAs = False
    # logger.info(
    #     'Done imputing missing genotypes.'
    # )
    # print('Done imputing missing genotypes.')
    return G
