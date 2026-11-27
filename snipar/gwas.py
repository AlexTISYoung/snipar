import gzip
import math
import time

import h5py
import numpy as np
from pysnptools.snpreader import Bed
from scipy.special import log_ndtr
import snipar.read as read
from snipar.read.variants import MatchedVariants, bgen_variant_ids, match_variants
import snipar.slmm as slmm
from snipar.gtarray import impute_missing, impute_unrel_par_gts
from snipar.utilities import encode_str_array, make_id_dict, open_bgen, outfile_name
from snipar.preprocess import find_par_gts
from typing import List
from multiprocessing import Pool, RawArray
from functools import partial


def find_common_ind_ids(obsfiles, impfiles, pheno_ids, from_chr=None, covar=None, keep=None, include_unrel=False):
    """Return samples eligible on every chromosome with observed and imputed genotypes."""
    chromosomes = [from_chr] if from_chr is not None else [chrom for chrom in obsfiles if chrom in impfiles]
    ids = np.asarray(pheno_ids)
    if keep is not None:
        ids = ids[np.isin(ids, np.loadtxt(keep, dtype=str, usecols=0, ndmin=1))]
    if covar is not None:
        ids = ids[np.isin(ids, covar.ids)]
    if not chromosomes or not ids.size:
        raise ValueError('No common ids among chromosome files.')
    for chromosome in sorted(chromosomes):
        ped, imp_fams = read.build_ped_from_par_gts(impfiles[chromosome])
        ids, fam_labels = read.get_ids_with_par(obsfiles[chromosome], ped, imp_fams, ids,
                                               include_unrel=include_unrel)
    return ids, fam_labels


def write_output(chrom, snp_ids, pos, alleles, outfile, parsum, sib, sib_diff, alpha, alpha_ses, alpha_cov, sigmas, freqs, robust=False, standard_gwas=False, trios_sibs=False, trios_only=False):
    """
    Write fitted SNP effects and other parameters to output HDF5 file.
    """
    print('Writing output to ' + outfile)
    outcols = ['population' if standard_gwas else 'direct']
    if not (robust or sib_diff or standard_gwas or (trios_sibs and not trios_only)):
        if sib:
            outcols.append('sib')
        outcols.extend(['avg_NTC'] if parsum else ['paternal', 'maternal'])
    with h5py.File(outfile, 'w') as output:
        output['bim'] = encode_str_array(np.column_stack((chrom, snp_ids, pos, alleles)))
        output['estimate_cols'] = encode_str_array(outcols)
        for name, data in (('estimate', alpha), ('estimate_ses', alpha_ses),
                           ('estimate_covariance', alpha_cov)):
            output.create_dataset(name, data=data, dtype='f', shuffle=True,
                                  compression='gzip', compression_opts=1)
        for i, sigma in enumerate(sigmas):
            output[f'sigma_{i}'] = sigma
        output['freqs'] = freqs


def outarray_effect(est, ses, freqs, vy):
    with np.errstate(invalid='ignore', divide='ignore'):
        N_effective = vy / (2 * freqs * (1 - freqs) * np.square(ses))
        Z = est / ses
        # P(chi-square(1) >= z²) equals twice the normal lower tail at -|z|.
        P = -(math.log(2) + log_ndtr(-np.abs(Z).astype(float))) / math.log(10)
    return np.column_stack((np.round(N_effective), np.round(est, 6), np.round(ses, 6),
                            np.round(Z, 6), np.round(P, 6)))


def write_txt_output(chrom, snp_ids, pos, alleles, outfile, parsum, sib, sib_diff, alpha, alpha_cov, sigmas, freqs, robust=False, standard_gwas=False, trios_sibs=False, trios_only=False):
    effects = ['population' if standard_gwas else 'direct']
    corr_pairs = []
    transform = np.eye(1)
    if not (robust or sib_diff or standard_gwas or (trios_sibs and not trios_only)):
        if sib:
            effects.append('sib')
        if not parsum:
            effects.extend(['paternal', 'maternal'])
        effects.extend(['avg_NTC', 'population'])
        transform = np.eye(len(effects), alpha.shape[1])
        avg_ntc = effects.index('avg_NTC')
        if not parsum:
            transform[avg_ntc, -2:] = 0.5
        transform[-1] = transform[avg_ntc]
        transform[-1, 0] = 1
        if sib:
            transform[-1, 1] = 0.5
        corr_pairs = [(0, avg_ntc), (0, avg_ntc + 1)]
        if sib:
            corr_pairs.append((0, 1))
        if not parsum:
            pat, mat = effects.index('paternal'), effects.index('maternal')
            corr_pairs.extend([(0, pat), (0, mat), (pat, mat)])
    header = ['chromosome', 'SNP', 'pos', 'A1', 'A2', 'freq']
    header.extend(f'{effect}_{field}' for effect in effects for field in ('N', 'Beta', 'SE', 'Z', 'log10_P'))
    header.extend(f'r_{effects[i]}_{effects[j]}' for i, j in corr_pairs)
    vy = sum(sigmas)
    print('Writing text output to ' + outfile)
    open_output = partial(gzip.open, compresslevel=1) if outfile.endswith('.gz') else open
    with open_output(outfile, 'wt') as output:
        output.write(' '.join(header) + '\n')
        for start in range(0, len(snp_ids), 8192):
            rows = slice(start, start + 8192)
            data = [chrom[rows], snp_ids[rows], pos[rows], alleles[rows], np.round(freqs[rows], 3)]
            if corr_pairs:
                betas = alpha[rows] @ transform.T
                covariance = transform @ alpha_cov[rows] @ transform.T
            else:
                betas, covariance = alpha[rows], alpha_cov[rows]
            with np.errstate(invalid='ignore', divide='ignore'):
                se = np.sqrt(np.diagonal(covariance, axis1=1, axis2=2))
                data.extend(outarray_effect(betas[:, i], se[:, i], freqs[rows], vy)
                            for i in range(len(effects)))
                data.extend(np.round(covariance[:, i, j] / (se[:, i] * se[:, j]), 6)
                            for i, j in corr_pairs)
            # Keep numeric columns numeric until formatting each row.
            blocks = [np.asarray(block).reshape(len(betas), -1) for block in data]
            output.writelines(' '.join(str(value) for part in row for value in part) + '\n'
                              for row in zip(*blocks))


def split_batches(nsnp: int, num_workers: int, batch_size = None) -> List[np.ndarray]:
    """
    Partition SNP indices [0, nsnp) into batches:
      - each batch has length <= batch_size
      - prefer the number of batches to be a multiple of num_workers
      - if that would force empty batches (nsnp too small), choose the largest feasible
        multiple; if none >= minimal needed, fall back to the minimal (non-multiple),
        and ALWAYS drop empties.
    """
    if batch_size is None:
        batch_size = nsnp
    if nsnp <= 0:
        return []
    if num_workers < 1:
        raise ValueError("num_workers must be >= 1")
    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")
    def round_up_multiple(x: int, m: int) -> int:
        return ((x + m - 1) // m) * m
    def round_down_multiple(x: int, m: int) -> int:
        return (x // m) * m
    # Minimal batches to satisfy the max-size constraint
    min_batches = math.ceil(nsnp / batch_size)
    # Candidate: nearest multiple of workers at or above min_batches
    cand = round_up_multiple(max(min_batches, num_workers), num_workers)
    if cand <= nsnp:
        total_batches = cand  # multiple of workers, no empties
    else:
        # Try the largest multiple <= nsnp that still respects min_batches
        down = round_down_multiple(nsnp, num_workers)
        if down >= max(min_batches, 1):
            total_batches = down
        else:
            # Fall back to minimal feasible (may NOT be a multiple)
            total_batches = min_batches
    # Balanced positive sizes summing to nsnp
    base = nsnp // total_batches
    rem  = nsnp % total_batches
    sizes = [base + 1] * rem + [base] * (total_batches - rem)
    # Build batches and DROP any accidental empties (defensive)
    idx = np.arange(nsnp, dtype=np.int64)
    batches, start = [], 0
    for s in sizes:
        if s <= 0:
            continue
        batches.append(idx[start:start + s])
        start += s
    # Sanity checks
    assert sum(len(b) for b in batches) == nsnp
    if batches:
        assert max(len(b) for b in batches) <= batch_size
    return batches

_var_dict = {}

def _init_worker(y_, n_vararr_, varcomps_, covar_, covar_shape, batch_processor=None, state=None, **kwargs):
    worker = _var_dict if state is None else state
    worker.clear()
    worker['y_'] = y_
    worker['varcomps_'] = varcomps_
    worker['batch_processor'] = batch_processor
    for i in range(n_vararr_):
        worker[f'varcomp_data_{i}_'] = kwargs[f'varcomp_data_{i}_']
        worker[f'varcomp_row_ind_{i}_'] = kwargs[f'varcomp_row_ind_{i}_']
        worker[f'varcomp_col_ind_{i}_'] = kwargs[f'varcomp_col_ind_{i}_']
    if covar_ is not None and covar_shape is not None:
        worker['covar_'] = covar_
        worker['covar_shape'] = covar_shape


def _run_worker_batch(batch):
    return _var_dict['batch_processor'](batch)


def process_batch(snp_ids, ped, n_vararr, imp_fams, pheno_ids=None, bedfile=None, bgenfile=None, par_gts_f=None,
                  parsum=False, fit_sib=False, sib_diff=False, max_missing=5, min_maf=0.01, verbose=False, 
                  print_sample_info=False, impute_unrel=False, standard_gwas=False, robust=False, trios_sibs=False,
                  add_jitter=False, worker_state=None):
    """Fit a batch given SNP IDs or their precomputed observed/imputed alignment."""
    worker = _var_dict if worker_state is None else worker_state
    matched_variants = None
    if isinstance(snp_ids, MatchedVariants):
        matched_variants, snp_ids = snp_ids, None
    ####### Construct family based genotype matrix #######
    G = read.get_gts_matrix(ped=ped, imp_fams=imp_fams, bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, snp_ids=snp_ids, 
                            ids=pheno_ids, parsum=parsum, sib=fit_sib, sib_diff=sib_diff,
                            include_unrel=impute_unrel, robust=robust, trios_sibs=trios_sibs,
                            verbose=False, print_sample_info=print_sample_info, matched_variants=matched_variants)
    if not np.array_equal(G.ids, pheno_ids):
        G.filter_ids(pheno_ids)
        if not np.array_equal(G.ids, pheno_ids):
            raise ValueError('Genotype batch does not contain every individual in the analysis cohort')
    if not impute_unrel and np.any(G.fams == ''):
        raise ValueError('No family label from pedigree for some individuals')
    G.compute_freqs()
    if impute_unrel and G.sid.size:
        G = impute_unrel_par_gts(G, sib=fit_sib, parsum=parsum, ped=ped)
    G.filter_maf(min_maf, verbose=verbose)
    G.filter_missingness(max_missing, verbose=verbose)
    if G.sid.size == 0:
        return G.freqs, G.sid, None, None, None
    if G.has_NAs:
        if sib_diff or fit_sib:
            G.fill_NAs()
        else:
            G = impute_missing(G)
    if 'model' not in worker:
        y = np.frombuffer(worker['y_'], dtype='float')
        varcomps = worker['varcomps_']
        # no need to provide residual var arr
        varcomp_arr_lst = tuple(
            (
                np.frombuffer(worker[f'varcomp_data_{i}_'], dtype='float'),
                np.frombuffer(worker[f'varcomp_row_ind_{i}_'], dtype='uint32'),
                np.frombuffer(worker[f'varcomp_col_ind_{i}_'], dtype='uint32'),
            ) for i in range(n_vararr)
        )
        covar = np.frombuffer(worker['covar_'], dtype='float').reshape(
            worker['covar_shape']) \
                if 'covar_' in worker else None
        if n_vararr == len(varcomps) == 2:
            # A leading GRM may be provided solely for imputation.
            varcomp_arr_lst = varcomp_arr_lst[1:]
        worker['model'] = slmm.LinearMixedModel(y, varcomp_arr_lst=varcomp_arr_lst,
                                                    varcomps=varcomps, covar_X=covar,
                                                    add_intercept=True, add_jitter=add_jitter)
    model = worker['model']
    if robust:
        G_sib = read.get_gts_matrix(ped=ped, imp_fams=imp_fams, bedfile=bedfile, bgenfile=bgenfile,
                                    snp_ids=G.sid, ids=pheno_ids, sib_diff=True, match_snp_ids=True)
        alpha, alpha_cov, alpha_ses = model.robust_est(G, G_sib)
    elif sib_diff:
        alpha, alpha_cov, alpha_ses = model.sib_diff_est(G.gts.data)
    elif trios_sibs:
        if G.sibs_inds.size == 0:
            genotypes = G.gts.data
            if parsum:
                genotypes = np.stack((genotypes[:, 0, :], genotypes[:, 1:, :].sum(axis=1)), axis=1)
            alpha, alpha_cov, alpha_ses = model.fit_snps_eff(genotypes)
        elif G.complete_trios_inds.size == 0:
            alpha, alpha_cov, alpha_ses = model.sib_diff_est(G.gts.data[:, :2, :])
        else:
            alpha, alpha_cov, alpha_ses = model.trios_sibs_est(G.gts.data, G.complete_trios_inds, G.sibs_inds)
    else:
        alpha, alpha_cov, alpha_ses = model.fit_snps_eff(G.gts.data, standard_gwas=standard_gwas)

    return G.freqs, G.sid, alpha, alpha_cov, alpha_ses


def process_chromosome(chrom_out, y, varcomp_lst,
                       ped, imp_fams, sigmas, outprefix, covariates=None, bedfile=None, bgenfile=None, par_gts_f=None,
                       fit_sib=False, sib_diff=False, parsum=False, impute_unrel=False, max_missing=5, min_maf=0.01, 
                       no_hdf5_out=False, no_txt_out=False, cpus=1, debug=False, robust=False, trios_sibs=False, trios_only=False,
                       add_jitter=False, standard_gwas=False, batch_size=None):
    ######## Check for bed/bgen #######
    if bedfile is None and bgenfile is None:
        raise(ValueError('Must supply either bed or bgen file with observed genotypes'))
    if bedfile is not None and bgenfile is not None:
        raise(ValueError('Both --bed and --bgen specified. Please specify one only'))
    if cpus < 1:
        raise ValueError('cpus must be positive')
    aliases = None
    if bedfile is not None:
        bed = Bed(bedfile,count_A1 = True)
        gts_id_dict = make_id_dict(bed.iid,1)
        snp_ids = bed.sid
        pos = np.array(bed.pos[:,2],dtype=int)
        alleles = np.loadtxt(bedfile.split('.bed')[0]+'.bim',dtype=str,usecols=(4,5),ndmin=2)
        chrom = np.array(bed.pos[:,0],dtype=int)
    elif bgenfile is not None:
        with open_bgen(bgenfile, verbose=False) as bgen:
            gts_id_dict = make_id_dict(bgen.samples)
            snp_ids, aliases = bgen_variant_ids(bgen)
            snp_ids = snp_ids.copy()
            if aliases is not None:
                aliases = aliases.copy()
            pos = np.array(bgen.positions, dtype=int)
            alleles = np.array([allele_ids.split(',') for allele_ids in bgen.allele_ids])
            chrom = np.array(bgen.chromosomes, dtype=str)
        if np.any(chrom == ''):
            chrom = np.where(chrom == '', str(chrom_out), chrom)
    # Check for observed parents if not using parsum
    if not parsum and not sib_diff:
        par_status, gt_indices, fam_labels = find_par_gts(y.ids, ped, gts_id_dict)
        parcount = np.sum(par_status==0,axis=1)
        if np.sum(parcount>0)==0 and not trios_sibs:
            # logger.warning('No individuals with genotyped parents found. Using sum of imputed maternal and paternal genotypes to prevent collinearity.')
            print('WARNING: no individuals with genotyped parents found. Using sum of imputed maternal and paternal genotypes to prevent collinearity.')
            parsum = True
        elif 100 > np.sum(parcount>0) > 0:
            # logger.warning('Warning: low number of individuals with observed parental genotypes. Consider using the --parsum argument to prevent issues due to collinearity.')
            print('WARNING: less than 100 individuals with observed parental genotypes. If the total sample size is small, consider using the --parsum argument to prevent issues due to collinearity.')
    ####### Compute batches #######
    print('Found '+str(snp_ids.shape[0])+' SNPs')
    # Align to original file columns before removing duplicate observed IDs from the output.
    matched = None
    if par_gts_f is not None and not sib_diff:
        with h5py.File(par_gts_f, 'r') as imputation:
            try:
                matched = match_variants(snp_ids, chrom, pos, alleles, imputation, aliases=aliases)
            except ValueError as error:
                label = chrom_out or ','.join(np.unique(chrom).astype(str))
                location = f'chromosome {label}, ' if label and label != '0' else ''
                raise ValueError(f'{error} ({location}{bedfile or bgenfile}, {par_gts_f})') from error
        print(f'Matched {matched.sid.size} SNPs between observed and imputed genotypes')
    # Remove duplicates
    unique_snps, counts = np.unique(snp_ids, return_counts=True)
    if np.any(counts > 1):
        non_duplicate = set(unique_snps[counts == 1])
        print('Removing '+str(np.sum(counts>1))+' duplicate SNP ids')
        not_duplicated = np.array([x in non_duplicate for x in snp_ids])
        snp_ids = snp_ids[not_duplicated]
        pos = pos[not_duplicated]
        chrom = chrom[not_duplicated]
        alleles = alleles[not_duplicated,:]
    snp_dict = make_id_dict(snp_ids)
    if robust or sib_diff or standard_gwas or (trios_sibs and not trios_only):
        alpha_dim = 1
    else:
        alpha_dim = 2
        if fit_sib:
            alpha_dim += 1
        if not parsum:
            alpha_dim += 1
    # Keep untested SNPs in output in original file order.
    alpha = np.full((snp_ids.size, alpha_dim), np.nan, dtype=np.float32)
    alpha_cov = np.full((snp_ids.size, alpha_dim, alpha_dim), np.nan, dtype=np.float32)
    alpha_ses = np.full_like(alpha, np.nan)
    freqs = np.full(snp_ids.size, np.nan, dtype=np.float32)
    variant_count = snp_ids.size if matched is None else matched.sid.size
    indices = (slice(idx[0], idx[-1] + 1) for idx in split_batches(variant_count, cpus, batch_size))
    batches = [snp_ids[idx] for idx in indices] if matched is None else [matched.subset(idx) for idx in indices]
    print(f'Using {len(batches)} batch' if len(batches) == 1 else f'Using {len(batches)} batches')

    def worker_array(values, dtype):
        values = np.ascontiguousarray(values, dtype=dtype)
        if cpus == 1:
            return values
        buffer = RawArray(np.ctypeslib.as_ctypes_type(values.dtype), values.size)
        np.copyto(np.frombuffer(buffer, dtype=values.dtype), values.ravel())
        return buffer

    y_ = worker_array(y.gts.data.reshape(-1), 'float64')
    covar_ = worker_array(covariates.gts.data, 'float64') if covariates is not None else None
    covar_shape = covariates.shape if covariates is not None else None
    varcomp_dict = {
        f'varcomp_{name}_{i}_': worker_array(values, dtype)
        for i, component in enumerate(varcomp_lst)
        for name, values, dtype in zip(('data', 'row_ind', 'col_ind'), component,
                                       ('float64', 'uint32', 'uint32'))
    }
    process_batch_ = partial(process_batch, pheno_ids=y.ids, n_vararr=len(varcomp_lst),
                             bedfile=bedfile, bgenfile=bgenfile,
                             par_gts_f=par_gts_f, ped=ped, imp_fams=imp_fams,
                             parsum=parsum, fit_sib=fit_sib, sib_diff=sib_diff, standard_gwas=standard_gwas,
                             max_missing=max_missing, min_maf=min_maf,
                             impute_unrel=impute_unrel,
                             robust=robust,
                             trios_sibs=trios_sibs,
                             add_jitter=add_jitter)
    _init_worker_ = partial(_init_worker, batch_processor=process_batch_, **varcomp_dict)
    worker_args = (y_, len(varcomp_lst), sigmas, covar_, covar_shape)
    if debug and batches:
        state = {}
        _init_worker_(*worker_args, state=state)
        process_batch_(batches[0], worker_state=state)
        raise SystemExit('Debug finished.')
    start = time.time()

    def store_results(results):
        for i, (batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses) in enumerate(results, 1):
            if len(batch_snps):
                batch_indices = np.array([snp_dict[x] for x in batch_snps])
                alpha[batch_indices] = batch_alpha
                alpha_cov[batch_indices] = batch_alpha_cov
                alpha_ses[batch_indices] = batch_alpha_ses
                freqs[batch_indices] = batch_freqs
            print(f'Done batch {i} out of {len(batches)} #snps:{len(batch_snps)}')

    if cpus == 1:
        state = {}
        _init_worker_(*worker_args, state=state)
        store_results(map(partial(process_batch_, worker_state=state), batches))
    elif batches:
        with Pool(processes=min(cpus, len(batches)), initializer=_init_worker_, initargs=worker_args) as pool:
            store_results(pool.imap_unordered(_run_worker_batch, batches, chunksize=1))
    print('Time used for inference: ', f'{time.time() - start:.2f}s')
    if not no_hdf5_out:
        if chrom_out==0:
            hdf5_outfile = outfile_name(outprefix, '.sumstats.hdf5')
        else:
            hdf5_outfile = outfile_name(outprefix, '.sumstats.hdf5', chrom=chrom_out)
        write_output(chrom, snp_ids, pos, alleles, hdf5_outfile, parsum, fit_sib, sib_diff, alpha, alpha_ses, alpha_cov,
                     sigmas, freqs, robust=robust, standard_gwas=standard_gwas, trios_sibs=trios_sibs, trios_only=trios_only)
    if not no_txt_out:
        if chrom_out==0:
            txt_outfile = outfile_name(outprefix, '.sumstats.gz')
        else:
            txt_outfile = outfile_name(outprefix, '.sumstats.gz', chrom=chrom_out)
        write_txt_output(chrom, snp_ids, pos, alleles, txt_outfile, parsum, fit_sib, sib_diff, alpha, alpha_cov,
                     sigmas, freqs, robust=robust, standard_gwas=standard_gwas, trios_sibs=trios_sibs, trios_only=trios_only)
