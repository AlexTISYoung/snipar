from snipar.gtarray import gtarray
import numpy as np
from snipar.read import bed, bgen, get_gts_matrix
from snipar.read.variants import NoMatchingVariantsError, bgen_variant_ids
import h5py
from snipar.utilities import *
from scipy.optimize import fmin_l_bfgs_b
from numba import njit, prange
import numpy.ma as ma
import snipar.slmm as slmm
from pysnptools.snpreader import Bed
import numdifftools as nd
import warnings
from pathlib import Path
from contextlib import nullcontext
from snipar.pgs_am import *
class pgs(object):
    """Define a polygenic score based on a set of SNPs with weights and ref/alt allele pairs.

    Args:
        snp_ids : :class:`~numpy:numpy.array`
            [L] vector of SNP ids
        weights : :class:`~numpy:numpy.array`
            [L] vector of weights of equal length to snp_ids
        alleles : :class:`~numpy:numpy.array`
            [L x 2] matrix of ref and alt alleles for the SNPs. L must match size of snp_ids

    Returns:
        pgs : :class:`snipar.pgs`

    """
    def __init__(self,snp_ids,weights,alleles):
        snp_ids = np.asarray(snp_ids)
        weights = np.asarray(weights)
        alleles = np.asarray(alleles)
        if snp_ids.ndim != 1 or weights.shape != snp_ids.shape or alleles.shape != (snp_ids.size, 2):
            raise ValueError('PGS SNP IDs, weights and allele pairs have incompatible shapes')
        if not np.isfinite(weights).all():
            raise ValueError('PGS weights must be finite')
        self.snp_ids = snp_ids
        self.snp_dict = make_id_dict(snp_ids)
        if len(self.snp_dict) != len(snp_ids):
            raise ValueError('PGS weights must have unique SNP IDs')
        self.weights = weights
        self.alleles = alleles
    
    def remove_zeros(self):
        zero_weight = self.weights==0
        n_zero = np.sum(zero_weight)
        if n_zero>0:
            print('Filtering '+str(n_zero)+' SNPs with zero weight')
            self.snp_ids = self.snp_ids[~zero_weight]
            self.snp_dict = make_id_dict(self.snp_ids)
            self.weights = self.weights[~zero_weight]
            self.alleles = self.alleles[~zero_weight]
        else:
            print('No zero weight SNPs found')

    def compute(self, garray, cols=None):
        """Compute polygenic score values from a given genotype array. Finds the SNPs in the genotype array
        that have weights in the pgs and matching alleles, and computes the PGS based on these SNPs and the
        weights after allele-matching.


        Args:
            garray : :class:`sbreg.gtarray`
                genotype array to compute PGS values for
            cols : :class:`numpy:numpy.array`
                names to give the columns in the output gtarray

        Returns:
            pg : :class:`snipar.gtarray`
                2d gtarray with PGS values. If a 3d gtarray is input, then each column corresponds to
                the second dimension on the input gtarray (for example, individual, paternal, maternal PGS).
                If a 2d gtarray is input, then there will be only one column in the output gtarray. The
                names given in 'cols' are stored in 'sid' attribute of the output.

        """
        if not isinstance(garray, gtarray):
            raise ValueError('Must be of gtarray class')
        if garray.alleles is None:
            raise ValueError('Alleles of genotype matrix must be provided')
        weight_indices = np.fromiter((self.snp_dict.get(snp, -1) for snp in garray.sid),
                                     dtype=np.intp, count=len(garray.sid))
        variant_indices = np.flatnonzero(weight_indices >= 0)
        if not variant_indices.size:
            print('No overlap between PGS SNPs and genotype SNPs')
            return None

        weight_indices = weight_indices[variant_indices]
        score_alleles = self.alleles[weight_indices]
        genotype_alleles = garray.alleles[variant_indices]
        forward = np.all(score_alleles == genotype_alleles, axis=1)
        reverse = np.all(score_alleles == genotype_alleles[:, ::-1], axis=1) & ~forward
        matched = forward | reverse
        if not np.all(matched):
            print(f'Removing {np.sum(~matched)} SNPs due to allele mismatch between genotypes and PGS alleles')
        if not np.any(matched):
            return None

        indices = variant_indices[matched]
        weights = self.weights[weight_indices[matched]] * np.where(reverse[matched], -1, 1)
        genotypes = garray.gts
        n_cols = 1 if garray.ndim == 2 else genotypes.shape[1]
        if cols is None:
            cols = ['proband'] + [f'relative_{i}' for i in range(1, n_cols)]
        cols = np.asarray(cols)
        if cols.shape != (n_cols,):
            raise ValueError('Provide one PGS column name per genotype column')

        # Limit each working array to one generation. Missing genotypes
        # contribute zero after centering by the proband allele dosage.
        genotypes = genotypes[:, None, :] if garray.ndim == 2 else genotypes
        if indices.size != genotypes.shape[-1]:
            proband = genotypes[:, 0, indices]
        else:
            proband = genotypes[:, 0, :]
        complete = not np.any(ma.getmask(genotypes))
        means = proband.data.mean(axis=0) if complete else proband.mean(axis=0)
        if complete:
            score_mean = means @ weights
        dtype = genotypes.dtype if garray.ndim == 3 and genotypes.dtype.kind == 'f' else np.float64
        result = np.empty((genotypes.shape[0], n_cols), dtype=dtype)
        for column, name in enumerate(cols):
            if column == 0:
                values = proband
            elif indices.size == genotypes.shape[-1]:
                values = genotypes[:, column, :]
            else:
                values = genotypes[:, column, indices]
            scale = 2 if name == 'parental' else 1
            if complete:
                # Center the projected score, avoiding a second sample x SNP
                # matrix and the full float64 cast used by mixed-dtype dot.
                result[:, column] = np.einsum('ij,j->i', values.data, weights, dtype=float) - scale * score_mean
            else:
                result[:, column] = ma.dot(values - scale * means, weights).filled(np.nan)
        return pgarray(result, garray.ids, sid=cols, fams=garray.fams,
                       par_status=garray.par_status, ped=garray.ped)


def read_weights(weights, SNP='SNP', beta_col='b', A1='A1', A2='A2', sep=None):
    table = np.loadtxt(weights, dtype=str, delimiter=sep, ndmin=2)
    if table.shape[0] < 2:
        raise ValueError('Weights file must contain a header and at least one variant')
    header = table[0].tolist()
    columns = [SNP, beta_col, A1, A2]
    if any(header.count(name) != 1 for name in columns):
        raise ValueError('Weight file must contain each required column exactly once: ' + ', '.join(columns))
    snp, beta, a1, a2 = [header.index(name) for name in columns]
    print(f'Read weights for {table.shape[0] - 1} variants')
    return pgs(table[1:, snp], table[1:, beta].astype(float), table[1:, [a1, a2]])


def read_pgs(pgs_file):
    table = np.loadtxt(pgs_file, dtype=str, comments=None, ndmin=2)
    if not table.size or not np.array_equal(table[0, :2], ['FID', 'IID']):
        raise ValueError('First two columns of PGS file must be FID and IID')
    has_parents = np.array_equal(table[0, 2:4], ['FATHER_ID', 'MOTHER_ID'])
    n_metadata = 4 if has_parents else 2
    header = table[0, n_metadata:]
    if not header.size or table.shape[0] < 2:
        raise ValueError('PGS file must contain at least one score column and one observation')
    if np.unique(header).size != header.size:
        raise ValueError('PGS column names must be unique')
    ids = table[1:, 1]
    if np.unique(ids).size != ids.size:
        raise ValueError('PGS individual IDs (IID) must be unique')
    ped = table[1:, :4].copy() if has_parents else None
    par_status = (ped[:, 2:4] == 'NA').astype(int) if has_parents else None
    values = table[1:, n_metadata:]
    values[np.isin(values, ['NA', 'N/A'])] = 'nan'
    return pgarray(values.astype(float), ids, sid=header, fams=table[1:, 0],
                   par_status=par_status, ped=ped)


def compute(pgs, bedfile=None, bgenfile=None, par_gts_f=None, ped=None, sib=False, compute_controls=False, verbose=True, batch_size=None):
    """Compute a polygenic score (PGS) for the individuals with observed genotypes and observed/imputed parental genotypes.

    Args:
        par_gts_f : :class:`str`
            path to HDF5 file with imputed parental genotypes
        gts_f : :class:`str`
            path to bed file with observed genotypes
        pgs : :class:`snipar.pgs`
            the PGS, defined by the weights for a set of SNPs and the alleles of those SNPs
        sib : :class:`bool`
            Compute the PGS for genotyped individuals with at least one genotyped sibling and observed/imputed parental genotypes. Default False.
        compute_controls : :class:`bool`
            Compute polygenic scores for control families (families with observed parental genotypes set to missing). Default False.

    Returns:
        pg : :class:`snipar.gtarray`
            Return the polygenic score as a genotype array with columns: individual's PGS, mean of their siblings' PGS, observed/imputed paternal PGS,
            observed/imputed maternal PGS

    """
    if (bedfile is None) == (bgenfile is None):
        raise ValueError('Provide one BED or BGEN file')
    if compute_controls and par_gts_f is None:
        raise ValueError('Control PGSs require an imputed parental genotype file')
    if batch_size is not None and (int(batch_size) != batch_size or batch_size < 1):
        raise ValueError('PGS batch size must be a positive integer')

    source = nullcontext(Bed(bedfile, count_A1=True)) if bedfile is not None else open_bgen(bgenfile, verbose=False)
    matched = None
    with source as observed:
        snp_ids = observed.sid if bedfile is not None else bgen_variant_ids(observed)[0]
        pgs_snp_set = set(pgs.snp_ids)
        snps = np.array([snp for snp in snp_ids if snp in pgs_snp_set], dtype=str)
        _, indices, counts = np.unique(snps, return_index=True, return_counts=True)
        snps = snps[np.sort(indices[counts == 1])]
        if not snps.size:
            print('No unique variants in common between weights file and observed genotypes')
            return None
        if par_gts_f is not None:
            try:
                with h5py.File(par_gts_f, 'r') as imputed:
                    if bedfile is not None:
                        matched = bed._match_variants(observed, imputed, Path(bedfile).with_suffix('.bim'), snps)
                    else:
                        matched = bgen._match_variants(observed, imputed, snps)
            except NoMatchingVariantsError as exc:
                print(exc)
                return None
            snps = matched.sid

    n_snps = snps.size
    batch_size = n_snps if batch_size is None else min(int(batch_size), n_snps)
    scores = None
    for start in range(0, n_snps, batch_size):
        batch = slice(start, start + batch_size)
        genotypes = get_gts_matrix(bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f,
                                   ped=ped, snp_ids=snps[batch], sib=sib, compute_controls=compute_controls,
                                   verbose=verbose and start == 0,
                                   matched_variants=None if matched is None else matched.subset(batch))
        for index, array in enumerate(genotypes if compute_controls else (genotypes,)):
            cols = ['proband', 'sibling'] if sib else ['proband']
            cols += ['paternal', 'maternal'] if array.shape[1] == len(cols) + 2 else ['parental']
            score = pgs.compute(array, cols)
            if score is None:
                continue
            if scores is None:
                scores = [None] * (4 if compute_controls else 1)
            scores[index] = score if scores[index] is None else scores[index].add(score)
    if scores is None:
        return None
    if compute_controls and any(score is None for score in scores):
        raise ValueError('No matching weighted variants for one or more control groups')
    return scores if compute_controls else scores[0]


class pgarray(gtarray):

    def add(self, garray):
        """Sum scores across independent variant sets, matching by individual ID."""
        if not isinstance(garray, pgarray):
            raise ValueError('Must add to another pgarray')
        if self.gts.shape[1:] != garray.gts.shape[1:] or not np.array_equal(self.sid, garray.sid):
            raise ValueError('PGS arrays must have the same score columns')
        if np.array_equal(self.ids, garray.ids):
            left = right = slice(None)
        else:
            common = [iid for iid in self.id_dict if iid in garray.id_dict]
            if not common:
                raise ValueError('No IDs in common')
            left = [self.id_dict[iid] for iid in common]
            right = [garray.id_dict[iid] for iid in common]
        a, b = self.gts[left], garray.gts[right]
        combined = ma.array(a.filled(0) + b.filled(0), mask=ma.getmask(a) & ma.getmask(b))
        fams = None if self.fams is None else self.fams[left]
        par_status = None if self.par_status is None else self.par_status[left]
        ped = self.ped if self.ped is not None else garray.ped
        return pgarray(combined, self.ids[left], self.sid, alleles=self.alleles,
                       fams=fams, par_status=par_status, ped=ped)

    def filter_bpg(self):
        if self.par_status is None:
            raise(ValueError('Parental genotype status unknown so cannot restrict to both parents genotyped sample'))
        else:
            bpg_ids = self.ids[np.sum(self.par_status==0, axis=1)==2]
            self.filter_ids(bpg_ids)

    def estimate_r(self, return_se=True, parents_only=False):
        if self.fams is None or self.par_status is None:
            raise ValueError('Family IDs and parental genotype status are needed to estimate parental correlation')
        family_ids, inverse, sizes = np.unique(self.fams, return_inverse=True, return_counts=True)
        family_sizes = dict(zip(family_ids, sizes))
        if not len(family_ids):
            return np.nan, np.nan, family_sizes
        # Group once, keeping the first member in the original order.
        order = np.argsort(inverse, kind='stable')
        first = order[np.concatenate(([0], np.cumsum(sizes[:-1])))]
        scores = ma.asarray(self.gts, dtype=float).filled(np.nan)
        parents = np.full((len(family_ids), 2), np.nan)
        for i, name in enumerate(('paternal', 'maternal')):
            observed = self.par_status[first, i] == 0
            if name in self.sid_dict:
                parents[observed, i] = scores[first[observed], self.sid_dict[name]]
        has_both = np.isfinite(parents).all(axis=1)
        n_both = has_both.sum()
        r_par = np.nan
        if n_both >= 2 and np.all(parents[has_both].std(axis=0) > 0):
            r_par = np.clip(np.corrcoef(parents[has_both].T)[0, 1], -1, 1)
        if parents_only:
            se = (1 - r_par**2) / np.sqrt(n_both - 1) if n_both > 2 else np.nan
            return r_par, se, family_sizes

        proband = scores[:, self.sid_dict['proband']]
        usable = np.isfinite(proband)
        proband = _standardize(proband[usable])
        usable_indices = inverse[usable][np.isfinite(proband)]
        proband = proband[np.isfinite(proband)]
        n_sib = np.bincount(usable_indices, minlength=len(family_ids))
        sib_sum = np.bincount(usable_indices, weights=proband, minlength=len(family_ids))
        sib_sq = np.bincount(usable_indices, weights=proband**2, minlength=len(family_ids))
        for i in range(2):
            observed = np.isfinite(parents[:, i])
            parents[observed, i] = _standardize(parents[observed, i])
        n_par = np.isfinite(parents).sum(axis=1)
        par_sum = np.nansum(parents, axis=1)
        par_sq = np.nansum(parents**2, axis=1)
        informative = (n_sib + n_par) >= 2
        if not np.any(informative):
            return np.nan, np.nan, family_sizes

        sib_mean = sib_sum / np.sqrt(np.maximum(n_sib, 1))
        par_mean = par_sum / np.sqrt(np.maximum(n_par, 1))
        counts, group = np.unique(np.column_stack((n_sib[informative], n_par[informative])),
                                  axis=0, return_inverse=True)
        moments = np.column_stack([np.bincount(group, weights=values[informative]) for values in (
            sib_mean**2, par_mean**2, sib_mean * par_mean,
            np.maximum(sib_sq - sib_mean**2, 0), np.maximum(par_sq - par_mean**2, 0))])
        group_sizes = np.bincount(group)
        n_sib_pairs = ((n_sib * (n_sib - 1)) // 2).sum()
        r_sib = (np.sum(sib_sum**2 - sib_sq) / n_sib_pairs - 1) if n_sib_pairs else np.nan
        init_values = np.array([r_par, r_sib])
        init_weights = np.array([n_both, n_sib_pairs])
        finite = np.isfinite(init_values)
        r_init = np.average(init_values[finite], weights=init_weights[finite]) if np.any(finite) else 0.
        r_init = np.clip(r_init, -.999, .999)

        print('Finding MLE for correlation between parents scores')
        def objective(x):
            return _pgs_profile_loglik(x[0], counts, group_sizes, moments)

        fit, _, info = fmin_l_bfgs_b(objective, x0=[r_init], approx_grad=True,
                                    bounds=[(-.999, .999)])
        if info['warnflag'] == 0:
            r = fit[0]
            curvature = nd.Hessian(objective, step=min(1e-3, (1 - abs(r)) / 4))(fit).item()
            se = np.sqrt(2 / curvature) if curvature > 0 else np.nan
        else:
            print('Could not find MLE for correlation. Returning the observed family correlations.')
            r = r_init
            se = (1 - r**2) / np.sqrt(n_both + n_sib_pairs - 1) if n_both + n_sib_pairs > 1 else np.nan
        return r, se, family_sizes

    def am_adj(self):
        print('Estimating correlation between maternal and paternal PGSs assuming equilibrium')
        r, r_se, fsizes = self.estimate_r()
        print(f'Estimated correlation between maternal and paternal PGSs: {r:.4f} S.E.={r_se:.4f}')
        if np.isfinite(r_se) and r > 1.5 * r_se:
            family_sizes = np.array([fsizes[fam] for fam in self.fams])
            n_parents = (self.par_status == 0).sum(axis=1)
            no_parents = n_parents == 0
            print('Adjusting imputed PGSs for assortative mating')
            if 'parental' in self.sid_dict and 'paternal' not in self.sid_dict:
                self.gts[no_parents, self.sid_dict['parental']] *= npg_am_adj(r, family_sizes[no_parents])
            else:
                father, mother = (self.sid_dict[name] for name in ('paternal', 'maternal'))
                factor = npg_am_adj(r, family_sizes[no_parents])
                self.gts[no_parents, father] *= factor
                self.gts[no_parents, mother] *= factor
                for index, (imputed, observed) in enumerate(((father, mother), (mother, father))):
                    use = (n_parents == 1) & (self.par_status[:, index] == 1)
                    self.gts[use, imputed] = opg_am_adj(
                        self.gts[use, imputed], self.gts[use, observed], r, family_sizes[use])
        else:
            print('Signal to noise ratio for correlation estimate too small, so not performing assortative mating adjustment.')
        return r, r_se

    def scale(self, sf=None):
        if sf is None:
            if 'proband' in self.sid:
                proband_index = np.where(self.sid=='proband')[0][0]
            else:
                raise(ValueError('Cannot scale as no proband column found and no scale factor provided'))
            sf = np.std(self.gts[:, proband_index])
        if ma.is_masked(sf) or not np.isfinite(sf) or sf <= 0:
            raise ValueError('PGS scale factor must be finite and positive')
        # Rescale by observed proband PGS
        self.gts = self.gts / sf
        return sf

    def write(self, filename, scale=False):
        if scale:
            self.scale()
        ####### Write PGS to file ########  
        if self.fams is None:
            raise ValueError('Family IDs are required to write a PGS file')
        metadata = [self.fams, self.ids]
        header = ['FID', 'IID']
        if self.ped is not None and self.par_status is not None:
            ped_dict = make_id_dict(self.ped, 1)
            parent_ids = self.ped[[ped_dict[iid] for iid in self.ids], 2:4].astype(object)
            parent_ids[self.par_status != 0] = 'NA'
            metadata.extend(parent_ids.T)
            header.extend(['FATHER_ID', 'MOTHER_ID'])
        values = ma.asarray(self.gts, dtype=float).filled(np.nan)
        pg_out = np.vstack((header + self.sid.tolist(), np.column_stack((*metadata, values))))
        print('Writing PGS to ' + str(filename))
        np.savetxt(filename, pg_out, fmt='%s')
        return pg_out
    
    def compute_grandpar(self, r):
        # Check we have required information
        if self.ped is None:
            raise(ValueError('Pedigree needed for grandparental model'))
        if self.par_status is None:
            raise(ValueError('Parental genotype status needed to compute grandparents'))
        # Check pgs columns
        if 'paternal' in self.sid:
            paternal_index = np.where(self.sid=='paternal')[0][0]
        else:
            raise(ValueError('No paternal PGS column found'))
        if 'maternal' in self.sid:
            maternal_index = np.where(self.sid=='maternal')[0][0]
        else:
            raise(ValueError('No maternal PGS column found'))
        # Find individuals with both parents genotyped
        bpg = np.sum(self.par_status==0,axis=1)==2
        n_bpg = np.sum(bpg)
        print('Found '+str(n_bpg)+' individuals with both parents genotyped for which grandparental scores will be computed')
        ## Create grandparental PGS matrix
        gpar = np.zeros((self.gts.shape[0],4))
        gpar[:] = np.nan
        bpg_ped = self.ped[:0]
        if n_bpg>0:
            # Find their parents' IDs
            ped_dict = make_id_dict(self.ped,1)
            bpg_ped = self.ped[[ped_dict[x] for x in self.ids[bpg]],:]
            # Find the mean of the mothers and fathers
            parent_means = np.mean(self.gts[bpg,:],axis=0)[[paternal_index, maternal_index]]
            ## Fill
            for i in range(bpg_ped.shape[0]):
                i_index = self.id_dict[bpg_ped[i,1]]
                ## Paternal grandparental scores
                # Find imputed grandparents
                if bpg_ped[i,2] in self.id_dict:
                    gpar[i_index, 0:2] = self.gts[self.id_dict[bpg_ped[i,2]],[paternal_index,maternal_index]]
                else:
                    # linear imputation from father if no imputed
                    gpar[i_index, 0:2] = parent_means[0]+(1+r)*(self.gts[i_index, paternal_index]-parent_means[0])/2.0
                ## Maternal grandparental scores
                # Find imputed grandparents
                if bpg_ped[i,3] in self.id_dict:
                    gpar[i_index, 2:4] = self.gts[self.id_dict[bpg_ped[i,3]],[paternal_index,maternal_index]]
                else:
                    # linear imputation from mother if no imputed
                    gpar[i_index, 2:4] = parent_means[1]+(1+r)*(self.gts[i_index, maternal_index]-parent_means[1])/2.0
        ## Append to PGS matrix
        self.gts = ma.column_stack((self.gts, ma.masked_invalid(gpar)))
        self.sid = np.hstack((self.sid,np.array(['gpp','gmp','gpm','gmm'])))
        self.sid_dict = make_id_dict(self.sid)
        self.shape = self.gts.shape
        return bpg_ped
        
def _standardize(values):
    if values.size < 2:
        return np.full_like(values, np.nan)
    scale = values.std()
    return (values - values.mean()) / scale if scale > 0 else np.full_like(values, np.nan)


def opg_am_adj(pgi_imp, pgi_obs, r, n):
    decay = np.exp2(-np.asarray(n, dtype=float))
    coefficient = r * decay / (1 + r - r * decay)
    return coefficient * pgi_obs + (1 + coefficient) * pgi_imp


def npg_am_adj(r, n):
    return (1 + r) / (1 + (1 - np.exp2(1 - np.asarray(n, dtype=float))) * r)


@njit(cache=True)
def _pgs_profile_loglik(r, counts, group_sizes, moments):
    # Sibling contrasts have variance (1-r)/2 and the parent contrast
    # (when present) has variance 1-r. Only the normalized sibling and
    # parent means have a joint covariance, of dimension at most two.
    k = (1 + r) / 2
    residual = (1 - r) / 2
    total = 0.
    for group in range(len(counts)):
        ns, np_ = counts[group]
        n = group_sizes[group]
        sib_sq, parent_sq, cross, sib_contrast, parent_contrast = moments[group]
        if ns > 1:
            total += n * (ns - 1) * np.log(residual) + sib_contrast / residual
        if np_ == 2:
            parent_variance = 1 + r
            total += n * np.log(1 - r) + parent_contrast / (1 - r)
        else:
            parent_variance = 1.
        if np_ > 0:
            total += n * np.log(parent_variance) + parent_sq / parent_variance
        if ns > 0:
            if np_ > 0:
                slope = np.sqrt(ns * np_) * k / parent_variance
                conditional_variance = residual if np_ == 2 else residual * (1 + ns * k)
                sib_sq = sib_sq - 2 * slope * cross + slope**2 * parent_sq
            else:
                conditional_variance = residual + ns * k
            total += n * np.log(conditional_variance) + sib_sq / conditional_variance
    return total


@njit
def pgs_corr_matrix(r,is_sib_fam,is_parent_fam):
    n_sib = np.sum(is_sib_fam)
    n_par = np.sum(is_parent_fam)
    # Full matrix
    R = np.zeros((n_sib+n_par,n_sib+n_par), dtype=np.float64)
    np.fill_diagonal(R, 1.0)
    # sib submatrix
    r_sib = (1+r)/2
    R_sib = r_sib*np.ones((n_sib,n_sib),dtype=np.float64)
    np.fill_diagonal(R_sib,1.0)
    R[0:n_sib,0:n_sib] = R_sib
    # sib to parent
    if n_par>0:
        R[0:n_sib,n_sib:(n_sib+n_par)] = r_sib
        R[n_sib:(n_sib+n_par),0:n_sib] = r_sib
    # parent-to-parent
    if n_par>1:
        R[n_sib+n_par-1,n_sib+n_par-2] = r
        R[n_sib+n_par-2,n_sib+n_par-1] = r
    # return 
    return R

@njit
def pgs_corr_likelihood_fam(r,pg_fam,is_sib_fam,is_parent_fam):
    sib_or_parent = np.logical_or(is_sib_fam,is_parent_fam)
    R = pgs_corr_matrix(r,is_sib_fam,is_parent_fam)
    slogdet_R = np.linalg.slogdet(R)
    pg_vec = pg_fam[sib_or_parent].reshape((np.sum(sib_or_parent),1))
    L_fam = slogdet_R[0]*slogdet_R[1]+pg_vec.T @ np.linalg.inv(R) @ pg_vec
    return L_fam[0,0]

@njit(parallel=True)
def pgs_corr_likelihood(r,pgs_fam,is_sib,is_parent):
    L = 0.0
    for i in prange(pgs_fam.shape[0]):
        L += pgs_corr_likelihood_fam(r, pgs_fam[i,:], is_sib[i,:], is_parent[i,:])
    return L

def pgs_cor_lik(r, *args):
    pgs_fam, is_sib, is_parent = args
    return pgs_corr_likelihood(r[0], pgs_fam, is_sib, is_parent)

def simulate_r_inf(r, nfam, nsib, npar):
    is_sib = np.zeros((nfam,nsib+npar),dtype=np.bool_)
    is_sib[:,0:nsib] = True
    is_parent = np.zeros((nfam,nsib+npar),dtype=np.bool_)
    is_parent[:,nsib:(nsib+npar)] = True
    R = pgs_corr_matrix(r,is_sib[0,:],is_parent[0,:])
    pgs_fam = np.random.multivariate_normal(np.zeros((nsib+npar)), R, size=nfam)
    # correlation between parents
    if npar==2:
        bpg = np.ones((nfam),dtype=bool)
    else:
        bpg = np.zeros((nfam),dtype=bool)
    n_bpg = np.sum(bpg)
    if n_bpg>0:
        r_bpg = np.corrcoef(pgs_fam[:,nsib],pgs_fam[:,nsib+1])[0,1]
    else:
        r_bpg = 0
    # Correlation from sibs
    if nsib>1:
        sib_2 = np.ones((nfam),dtype=bool)
    else:
        sib_2 = np.zeros((nfam),dtype=bool)
    n_sib_2 = np.sum(sib_2)
    if n_sib_2>0:
        r_sib = 2*np.corrcoef(pgs_fam[sib_2,0],pgs_fam[sib_2,1])[0,1]-1
    else:
        r_sib = 0
    # Initialize at weighted average
    r_init = n_bpg*r_bpg+n_sib_2*r_sib
    r_init = r_init/(n_bpg+n_sib_2)
    # Find MLE
    optimized = fmin_l_bfgs_b(func=pgs_cor_lik, x0=r_init,
                            args=(pgs_fam, is_sib, is_parent),
                            approx_grad=True,
                            bounds=[(-0.999,0.999)])
    return np.array([r_init,optimized[0][0]])

def fit_pgs_model(y, pg, ngen, ibdrel_path=None, covariates=None, fit_sib=False, parsum=False, gparsum=False, outprefix=None, sparse_thresh=0.025, return_fit=False, sample_ids=None, write_effects=True):
    if ngen not in (1, 2, 3):
        raise ValueError('ngen must be 1, 2, or 3')
    if not np.array_equal(y.ids, pg.ids):
        pg.filter_ids(y.ids)
        y.filter_ids(pg.ids)
    if covariates is not None and not np.array_equal(pg.ids, covariates.ids):
        covariates.filter_ids(pg.ids)
        if not np.array_equal(pg.ids, covariates.ids):
            raise ValueError('Covariates are required for every PGS individual in the analysis')

    pg_cols = ['proband']
    score_sums = {}
    if ngen == 1:
        print('Fitting 1 generation model (proband only)')
    else:
        if fit_sib:
            if 'sibling' not in pg.sid:
                raise ValueError('Sibling PGS not found (use --fit_sib when calculating PGS)')
            pg_cols.append('sibling')
        if parsum:
            pg_cols.append('parental')
            if {'paternal', 'maternal'}.issubset(pg.sid):
                score_sums['parental'] = ('paternal', 'maternal')
            elif 'parental' not in pg.sid:
                raise ValueError('Maternal and paternal PGS not found so cannot sum (--parsum option given)')
        else:
            pg_cols.extend(('paternal', 'maternal'))
        if ngen == 2:
            print('Fitting 2 generation model (proband and observed/imputed parents)')
        else:
            print('Fitting 3 generation model: observed proband and observed parents, and observed/imputed grandparents')
            if gparsum:
                pg_cols.extend(('gp', 'gm'))
                if {'gpp', 'gmp', 'gpm', 'gmm'}.issubset(pg.sid):
                    score_sums.update(gp=('gpp', 'gmp'), gm=('gpm', 'gmm'))
                elif not {'gp', 'gm'}.issubset(pg.sid):
                    raise ValueError('Grandparental PGSs not found so cannot sum (--gparsum option given)')
            else:
                pg_cols.extend(('gpp', 'gpm', 'gmp', 'gmm'))
    fit = make_and_fit_model(y, pg, pg_cols, ibdrel_path=ibdrel_path,
                             covariates=covariates, sparse_thresh=sparse_thresh,
                             return_fit=True, sample_ids=sample_ids, score_sums=score_sums)
    alpha, alpha_cols = fit['alpha'], fit['cols']
    # Save to file
    if outprefix is not None:
        write_estimates(str(outprefix)+'.'+str(ngen), alpha, alpha_cols, write_effects=False)
        if write_effects:
            write_pgs_effects(outprefix, {ngen: (alpha, alpha_cols)})
    return fit if return_fit else (alpha, alpha_cols)

def make_and_fit_model(y, pg, pg_cols, ibdrel_path=None, covariates=None, sparse_thresh=0.025,
                       return_fit=False, sample_ids=None, score_sums=None):
    score_sums = {} if score_sums is None else score_sums
    score_columns = []
    for name in pg_cols:
        sources = score_sums.get(name, (name,))
        if any(source not in pg.sid_dict for source in sources):
            raise ValueError('Required PGS column not found: ' + name)
        col = ma.asarray(pg.gts[:, pg.sid_dict[sources[0]]], dtype=float).filled(np.nan)
        for source in sources[1:]:
            col = col + ma.asarray(pg.gts[:, pg.sid_dict[source]], dtype=float).filled(np.nan)
        score_columns.append(col)
    X = np.column_stack(score_columns)
    if covariates is not None:
        X = np.hstack((ma.asarray(covariates.gts, dtype=float).filled(np.nan), X))
        X_cols = np.hstack((['intercept'], covariates.sid, pg_cols))
    else:
        X_cols = np.hstack((['intercept'], pg_cols))
    response = ma.asarray(y.gts[:, 0], dtype=float).filled(np.nan)
    no_NA = np.all(np.isfinite(X), axis=1) & np.isfinite(response)
    if sample_ids is not None:
        no_NA &= np.isin(pg.ids, sample_ids)
    n = int(no_NA.sum())
    print(f'Sample size: {n}')
    if n <= X.shape[1] + 1:
        raise ValueError('Too few complete PGS and phenotype observations for this model')
    ## Read GRM
    id_dict = make_id_dict(pg.ids[no_NA])
    varcomp_lst = make_grms(pg.fams[no_NA], ibdrel_path=ibdrel_path, id_dict=id_dict, keep=pg.ids[no_NA], sparse_thresh=sparse_thresh)
    ## Define LMM
    slmm_model = slmm.LinearMixedModel(response[no_NA],
                    varcomp_arr_lst=varcomp_lst, covar_X=np.array(X[no_NA,:],dtype=float), add_intercept=True)
    # Optimize Lmm
    print('Optimizing linear mixed model')
    slmm_model.scipy_optimize()
    # Print variance components
    varcomps = list(i for i in slmm_model.varcomps)
    vcomp_index=0
    print('Variance components (after projecting out PGS and covariates):')
    if ibdrel_path is not None:
        print(f'GRM variance component: {varcomps[vcomp_index]:.3f}')
        print(f'Variance explained by GRM: {varcomps[vcomp_index] / np.sum(varcomps) * 100:.1f}%')
        vcomp_index+=1
    print(f'Sibling variance component: {varcomps[vcomp_index]:.3f}')
    print(f'Variance explained by sibling component: {varcomps[vcomp_index] / np.sum(varcomps) * 100:.1f}%')
    print(f'Total variance: {np.sum(varcomps):.3f}')
    # Print implied sibling correlation
    if ibdrel_path is None:
        print(f'Implied sibling correlation: {varcomps[0] / np.sum(varcomps):.3f}')
    else:
        print(f'Implied sibling correlation: {(0.5*varcomps[0] + varcomps[1]) / np.sum(varcomps):.3f}')
    # Compute estimates
    ZT_Vinv_Z_imp = slmm_model.Z.T @ slmm_model.Vinv_Z
    alpha = [np.linalg.solve(ZT_Vinv_Z_imp, slmm_model.Z.T @ slmm_model.Vinv_y), np.linalg.inv(ZT_Vinv_Z_imp)]
    if return_fit:
        return dict(alpha=alpha, cols=X_cols, ids=pg.ids[no_NA].copy(),
                    model=slmm_model)
    return alpha, X_cols


def pgs_model_cross_covariance(population_fit, direct_fit):
    """Sampling covariance between two GLS fits on the same individuals.

    Let A_j = (X_j' V_j^-1 X_j)^-1 X_j' V_j^-1. The residual from
    the full, two-generation model is shared by both fits; omitted parental
    contributions in the population residual are assumed orthogonal to it.
    Therefore Cov(b_1, b_2) = A_1 V_2 A_2'. Since V_2 A_2' = X_2 C_2,
    this needs only sparse solves already performed during fitting. It also
    works when separately estimated variance components are not proportional.

    The returned matrix has population-model rows and direct-model columns.
    Variance-component estimation uncertainty is ignored, as in the existing
    GLS standard errors. The model-specific marginal variances remain C_1
    and C_2; geometric-mean scaling of residual variances is not appropriate.
    """
    if not np.array_equal(population_fit['ids'], direct_fit['ids']):
        raise ValueError('PGS ratio models must use the same individuals in the same order')
    model_1, model_2 = population_fit['model'], direct_fit['model']
    if not np.array_equal(model_1.y, model_2.y):
        raise ValueError('PGS ratio models must use the same phenotype and scale')
    covariance_1, covariance_2 = population_fit['alpha'][1], direct_fit['alpha'][1]
    return covariance_1 @ (model_1.Vinv_Z.T @ model_2.Z) @ covariance_2


def direct_population_ratio(population_fit, direct_fit, outprefix=None, write_effects=True):
    """Estimate the 2-gen/1-gen proband coefficient ratio and delta-method SE.

    Uses the marginal variance from each fitted model and their shared-error
    covariance. A weak denominator is flagged because the delta approximation
    can be inaccurate near zero. Substantially non-positive-semidefinite
    fitted covariance estimates are flagged and produce an undefined SE.
    With outprefix, writes the cross-model covariance matrix and one effects table
    containing both models and the ratio, unless write_effects is False.
    """
    cross = pgs_model_cross_covariance(population_fit, direct_fit)
    pop_index = list(population_fit['cols']).index('proband')
    direct_index = list(direct_fit['cols']).index('proband')
    population = float(population_fit['alpha'][0][pop_index])
    direct = float(direct_fit['alpha'][0][direct_index])
    pop_var = float(population_fit['alpha'][1][pop_index, pop_index])
    direct_var = float(direct_fit['alpha'][1][direct_index, direct_index])
    covariance = float(cross[pop_index, direct_index])
    joint = np.array([[pop_var, covariance], [covariance, direct_var]])
    result = dict(ratio=np.nan, se=np.nan, direct=direct, population=population,
                  direct_se=np.sqrt(direct_var), population_se=np.sqrt(pop_var),
                  covariance=covariance,
                  sampling_correlation=covariance / np.sqrt(pop_var * direct_var),
                  n=len(population_fit['ids']), status='ok',
                  method='GLS_shared_full_model_residual_delta')
    if not np.all(np.isfinite(joint)) or min(pop_var, direct_var) <= 0:
        result['status'] = 'invalid_covariance'
    elif np.linalg.eigvalsh(joint)[0] < -1e-10 * max(pop_var, direct_var):
        result['status'] = 'inconsistent_marginal_covariance'
    if not np.isfinite(direct):
        result['status'] = 'invalid_direct_coefficient'
    elif population == 0 or not np.isfinite(population):
        result['status'] = 'undefined_zero_denominator'
    else:
        result['ratio'] = direct / population
        if result['status'] == 'ok':
            gradient = np.array([-direct / population ** 2, 1.0 / population])
            variance = float(gradient @ joint @ gradient)
            result['se'] = np.sqrt(max(0.0, variance))
            if abs(population) < 1.96 * np.sqrt(pop_var):
                result['status'] = 'weak_denominator'
    if result['status'] != 'ok':
        warnings.warn('Direct/population PGS ratio: ' + result['status'], RuntimeWarning)
    print('Direct/population PGS ratio: {} (SE {}); {}'.format(
        result['ratio'], result['se'], result['status']))
    if outprefix is not None:
        outprefix = str(outprefix)
        Path(outprefix).parent.mkdir(parents=True, exist_ok=True)
        cross_out = np.vstack((
            np.hstack((['term'], direct_fit['cols'])),
            np.column_stack((population_fit['cols'], cross))))
        np.savetxt(outprefix + '.1_2.cross_vcov.txt', cross_out, fmt='%s', delimiter='\t')
        if write_effects:
            write_pgs_effects(outprefix, {
                1: (population_fit['alpha'], population_fit['cols']),
                2: (direct_fit['alpha'], direct_fit['cols']),
            }, ratio=result)
    return result


def fit_pgs_model_pair(y, pg, write_effects=True, **kwargs):
    """Fit both generations on the full model's complete cases and their ratio.

    The full model is fitted first so absent parental scores cannot change the
    sample used for the population coefficient. Existing return values are
    retained; output combines both models and the ratio in one effects table.
    """
    direct_fit = fit_pgs_model(y, pg, 2, return_fit=True, write_effects=False, **kwargs)
    population_kwargs = dict(kwargs, sample_ids=direct_fit['ids'])
    population_fit = fit_pgs_model(y, pg, 1, return_fit=True, write_effects=False,
                                    **population_kwargs)
    ratio = direct_population_ratio(population_fit, direct_fit,
                                    outprefix=kwargs.get('outprefix'),
                                    write_effects=write_effects)
    return population_fit, direct_fit, ratio

def write_pgs_effects(outprefix, models, ratio=None, adjusted=None):
    """Write all requested coefficient estimates and SEs to one table.

    models maps generation numbers to the (alpha, column_names) return values
    of fit_pgs_model. Optional adjusted is the (estimates, SEs) pair returned
    by am_adj_2gen. Covariance matrices are written separately by the fitters.
    """
    rows = []
    for generation in sorted(models):
        alpha, cols = models[generation]
        for name, estimate, se in zip(cols, alpha[0], np.sqrt(np.diag(alpha[1]))):
            statistic = '{}.{}'.format(generation, name) if generation != '' else str(name)
            rows.append((statistic, estimate, se))
    if ratio is not None:
        rows.append(('direct_population_ratio', ratio['ratio'], ratio['se']))
    if adjusted is not None:
        estimates, ses = adjusted
        for name in ('delta', 'alpha', 'rk', 'k', 'r', 'h2_eq', 'rho',
                     'alpha_delta', 'v_eta_delta'):
            rows.append(('am_adj.' + name, estimates[name], ses[name]))
    outprefix = str(outprefix)
    Path(outprefix).parent.mkdir(parents=True, exist_ok=True)
    print('Saving effect estimates to ' + outprefix + '.effects.txt')
    with open(outprefix + '.effects.txt', 'w') as output:
        output.write('statistic\test\tSE\n')
        for name, estimate, se in rows:
            values = [str(float(value)) if np.isfinite(value) else 'NA'
                      for value in (estimate, se)]
            output.write('\t'.join([name] + values) + '\n')


def write_estimates(outprefix, alpha, cols, write_effects=True):
    """Write a model covariance, optionally with its coefficient table."""
    outprefix = str(outprefix)
    Path(outprefix).parent.mkdir(parents=True, exist_ok=True)
    if write_effects:
        write_pgs_effects(outprefix, {'': (alpha, cols)})
    cols = cols.reshape((cols.shape[0],1))
    vcols = np.zeros((1,cols.shape[0]+1),dtype=cols.dtype)
    vcols[0,1:vcols.shape[1]] = cols[:,0]
    alpha_cov_out = np.vstack((vcols, np.hstack((cols,alpha[1]))))
    print('Saving sampling variance-covariance matrix to '+outprefix+ '.vcov.txt')
    np.savetxt(outprefix+ '.vcov.txt',alpha_cov_out, fmt='%s')

def make_grms(fams, ibdrel_path=None, id_dict=None, keep=None, sparse_thresh=0.05):
    if ibdrel_path is not None:
        grm_data, grm_row_ind, grm_col_ind = slmm.build_ibdrel_arr(
            ibdrel_path, id_dict=id_dict, keep=keep, thres=sparse_thresh)
    ## Build sparse sib GRM
    sib_data, sib_row_ind, sib_col_ind = slmm.build_sib_arr(fams)
    ## GRM list
    if 'grm_data' in locals():
        varcomp_lst = (
            (grm_data, grm_row_ind, grm_col_ind),
            (sib_data, sib_row_ind, sib_col_ind),
        )
    else:
        varcomp_lst = (
            (sib_data, sib_row_ind, sib_col_ind),
        )
    return varcomp_lst

######## Assortative mating adjustment in two-generation model ###########

def h2f_parse(h2f_str):
    h2f_split = h2f_str.split(',')
    if not len(h2f_split)==2:
        raise(ValueError('Invalid h2f input. Please use estimate,se format, such as 0.5,0.01.'))
    h2f = float(h2f_split[0])
    h2f_se = float(h2f_split[1].split('\n')[0])
    return h2f, h2f_se

def am_adj_2gen_calc(delta, delta_se, ab, ab_se, r_delta_ab, h2f, h2f_se, rk, rk_se, is_beta=False, verbose=True):
    # Get estimates
    if is_beta:
        estimates = v_eta_delta(delta, delta_se, h2f, h2f_se, rk, beta=ab)
        estimates['beta'] = ab
    else:
        estimates = v_eta_delta(delta, delta_se, h2f, h2f_se, rk, alpha=ab)
        estimates['alpha'] = ab
    estimates['rk'] = rk
    estimates['delta'] = delta
    
    # Get ses
    ses = {'rk':rk_se,
           'k':k_se(delta, delta_se, h2f, h2f_se, rk, rk_se),
           'r':r_se(delta, delta_se, h2f, h2f_se, rk, rk_se),
           'h2_eq':h2eq_se(delta, delta_se, h2f, h2f_se, rk, rk_se),
           'rho':rho_se(delta, delta_se, h2f, h2f_se, rk, rk_se),
           'delta':delta_se,
           'v_eta_delta':se_v_eta_delta(delta, delta_se, ab, ab_se, r_delta_ab, h2f, h2f_se, rk, rk_se, is_beta=is_beta)}
    if is_beta:
        ses['alpha_delta'] = se_alpha_from_beta(delta, delta_se, ab, ab_se, r_delta_ab, h2f, h2f_se, rk, rk_se)
        ses['beta'] = ab_se
        ses['r_delta_beta'] = r_delta_ab
    else:
        ses['alpha_delta'] = se_alpha_from_alpha(delta, delta_se, ab, ab_se, r_delta_ab, h2f, h2f_se, rk, rk_se)
        ses['alpha'] = ab_se
        ses['r_delta_alpha'] = r_delta_ab
    # Print
    if verbose:
        for par in ['rk','k','r','h2_eq','rho','alpha_delta','v_eta_delta']:
            print(par+': '+str(round(estimates[par],4))+' ('+str(round(ses[par],4))+')')
    # Return
    return estimates, ses

def write_2gen_adj_ests(estimates,ses, outprefix=''):
    pars = np.array(['delta','alpha','rk','k','r','h2_eq','rho','alpha_delta','v_eta_delta'])
    outarray = np.zeros((pars.shape[0],2))
    parcount = 0
    for par in pars:
        outarray[parcount,:] = [estimates[par],ses[par]]
        parcount += 1
    outarray = np.hstack((pars.reshape((pars.shape[0],1)),outarray))
    outarray = np.vstack((np.array(['parameter','estimate','SE']).reshape((1,3)),outarray))
    np.savetxt(outprefix+'.am_adj_pars.txt',outarray,fmt='%s')
    return outarray

def am_adj_2gen(estimates, estimate_cols, h2f, h2f_se, rk=None, rk_se=None, pg=None, y_std=1, pg_std=1):
    """
    Adjust 2-generation model results for assortative mating (assuming equilibrium)
    """
    ##### Check we have r, or a pg object to estimate r with #####
    if pg is None:
        if rk is None:
            raise(ValueError('Need estimate of correlation between parents PGIs if not providing pg object to use to estimate'))
        if rk_se is None:
            raise(ValueError('Need standard error of correlation between parents PGIs if not providing pg object to use to estimate')) 
    ##### Find estimates and their sampling variance-covariance #####
    if 'proband' in estimate_cols:
        proband_index = np.where(estimate_cols=='proband')[0][0]
    else:
        raise(ValueError('Proband coefficient (direct effect) estimate needed for adjustment'))
    # If separate maternal and paternal NTC estimates, take average
    if 'paternal' in estimate_cols and 'maternal' in estimate_cols:
        paternal_index = np.where(estimate_cols=='paternal')[0][0]
        maternal_index = np.where(estimate_cols=='maternal')[0][0]
        est_cols = [proband_index, paternal_index, maternal_index]
        estimate = estimates[0][est_cols].reshape((3,1))
        estimate_cov = estimates[1][np.ix_(est_cols,est_cols)]
        A =  np.array([[1,0,0],[0,0.5,0.5]])
        estimate = A @ estimate
        estimate_cov = A @ estimate_cov @ A.T
    elif 'parental' in estimate_cols:
        est_cols = [proband_index, np.where(estimate_cols=='parental')[0][0]]
        estimate = estimates[0][est_cols].reshape((2,1))
        estimate_cov = estimates[1][np.ix_(est_cols,est_cols)]
    else:
        raise(ValueError('Need parental NTC estimate(s) for adjustment')) 
    # Adjust for non-normalized pgs/phenotype
    estimate = estimate*pg_std/y_std
    estimate_cov = estimate_cov*((pg_std/y_std)**2)
    ##### Estimate correlation between parents' PGIs ######
    if rk is None:
        rk, rk_se, _ = pg.estimate_r(parents_only=True)
    if rk_se is None or not np.isfinite(rk_se) or rk_se < 0 or not np.isfinite(rk):
        raise ValueError('Provide a finite parental PGS correlation and SE with --rk when they cannot be estimated')
    ##### Estimate equilibrium quantities ######
    print('Computing equilibrium adjusted quantities')
    adj_estimates, adj_ses = am_adj_2gen_calc(estimate[0,0], np.sqrt(estimate_cov[0,0]), 
                                              estimate[1,0], np.sqrt(estimate_cov[1,1]),
                                              estimate_cov[0,1]/np.sqrt(estimate_cov[0,0]*estimate_cov[1,1]), 
                                              h2f, h2f_se, rk, rk_se, verbose=True)
    return adj_estimates, adj_ses
