import gzip
import math
import time

import h5py
import numpy as np
import pandas as pd
from scipy.special import log_ndtr
import snipar.read as read
from snipar.read.variants import MatchedVariants
import snipar.slmm as slmm
from snipar.gtarray import impute_missing, impute_unrel_par_gts
from snipar.utilities import encode_str_array, make_id_dict, outfile_name
from typing import List
from functools import partial
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, wait
from contextlib import ExitStack
from itertools import islice
from multiprocessing import RawArray, get_context
from multiprocessing.util import Finalize
from threading import Lock, local
from typing import NamedTuple

from snipar.resources import MIB, available_cpus, limit_native_threads, plan_batches, worker_environment


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
    index = np.arange(max(nsnp, 0), dtype=np.int64)
    return [index[batch] for batch in _batch_slices(nsnp, num_workers, batch_size)]


def _number_of_batches(nsnp, num_workers, batch_size):
    if nsnp <= 0:
        return 0
    if batch_size is None:
        batch_size = nsnp
    if num_workers < 1:
        raise ValueError('num_workers must be >= 1')
    if batch_size < 1:
        raise ValueError('batch_size must be >= 1')
    needed = (nsnp + batch_size - 1) // batch_size
    up = ((needed + num_workers - 1) // num_workers) * num_workers
    if up <= nsnp:
        return up
    down = (nsnp // num_workers) * num_workers
    return down if down >= needed else needed


def _batch_slices(nsnp, num_workers, batch_size):
    total = _number_of_batches(nsnp, num_workers, batch_size)
    if total:
        size, extra = divmod(nsnp, total)
        start = 0
        for i in range(total):
            stop = start + size + (i < extra)
            yield slice(start, stop)
            start = stop


def _bounded_results(executor, task, batches, workers):
    batches = iter(batches)
    pending = {executor.submit(task, batch) for batch in islice(batches, workers)}
    while pending:
        completed, pending = wait(pending, return_when=FIRST_COMPLETED)
        for future in completed:
            yield future.result()
            batch = next(batches, None)
            if batch is not None:
                pending.add(executor.submit(task, batch))


class _SharedArray(NamedTuple):
    memory: object
    shape: tuple
    dtype: np.dtype

    def view(self):
        array = np.frombuffer(self.memory, dtype=self.dtype).reshape(self.shape)
        array.flags.writeable = False
        return array


def _share_array(array):
    if array is None:
        return None
    if array.dtype.hasobject:
        raise ValueError('GWAS worker arrays must have a fixed-width data type')
    memory = RawArray('b', array.nbytes)
    np.copyto(np.frombuffer(memory, dtype=array.dtype).reshape(array.shape), array)
    return _SharedArray(memory, array.shape, array.dtype)


def _share_plans(plans):
    if plans is not None:
        return tuple(type(plan)(*(_share_array(field) if isinstance(field, np.ndarray) else field
                                   for field in plan)) for plan in plans)


def _restore_plans(plans):
    if plans is not None:
        return tuple(type(plan)(*(field.view() if isinstance(field, _SharedArray) else field
                                   for field in plan)) for plan in plans)


def _fit_model(y, components, sigmas, covariates, family_labels, add_jitter, solver='lu'):
    if len(components) == len(sigmas) == 2:
        # A leading GRM may be provided solely for imputation.
        components = components[1:]
    return slmm.LinearMixedModel(y, components, sigmas, covariates, add_intercept=True,
                                 add_jitter=add_jitter, family_labels=family_labels,
                                 snp_solver=solver)


def _readers(stack, options, plans, sibling_plans):
    reader = stack.enter_context(read.prepare_gts_matrix(**options, sample_plans=plans))
    sibling = None
    if sibling_plans is not None:
        sibling_opts = dict(options, par_gts_f=None, sib_diff=True, robust=False)
        sibling = stack.enter_context(read.prepare_gts_matrix(**sibling_opts, sample_plans=sibling_plans))
    return reader, sibling


def _start_process_worker(options, plans, sibling_plans, y, components, sigmas, covar,
                          family_labels, add_jitter, threads, processor, solver, model_ids, variants):
    from numba import config, set_num_threads

    set_num_threads(min(threads, config.NUMBA_NUM_THREADS))
    stack = ExitStack()
    try:
        stack.enter_context(limit_native_threads(threads))
        reader, sibling = _readers(stack, options, _restore_plans(plans), _restore_plans(sibling_plans))
        arrays = tuple(tuple(a.view() for a in component) for component in components)
        model = _fit_model(y.view(), arrays, sigmas, None if covar is None else covar.view(),
                           None if family_labels is None else family_labels.view(), add_jitter, solver)
    except BaseException:
        stack.close()
        raise
    Finalize(None, stack.close, exitpriority=1)
    _var_dict.clear()
    _var_dict.update(model=model, reader=reader, sibling_reader=sibling, batch_processor=processor,
                     model_ids=model_ids.view(), variants=MatchedVariants(*(a.view() for a in variants)))


_var_dict = {}

def _run_worker_batch(batch):
    return _var_dict['batch_processor'](_var_dict['variants'].subset(batch))


def process_batch(snp_ids, ped, n_vararr, imp_fams, pheno_ids=None, bedfile=None, bgenfile=None, par_gts_f=None,
                  parsum=False, fit_sib=False, sib_diff=False, max_missing=5, min_maf=0.01, verbose=False, 
                  print_sample_info=False, impute_unrel=False, standard_gwas=False, robust=False, trios_sibs=False,
                  add_jitter=False, worker_state=None):
    """Fit a batch given SNP IDs or their precomputed observed/imputed alignment."""
    worker = _var_dict if worker_state is None else worker_state
    matched_variants = None
    if isinstance(snp_ids, MatchedVariants):
        matched_variants, snp_ids = snp_ids, None
    # Prepared readers keep aligned row plans and handles for the worker's lifetime.
    reader = worker.get('reader')
    if reader is None:
        G = read.get_gts_matrix(ped=ped, imp_fams=imp_fams, bedfile=bedfile, bgenfile=bgenfile,
                                par_gts_f=par_gts_f, snp_ids=snp_ids, ids=pheno_ids, parsum=parsum,
                                sib=fit_sib, sib_diff=sib_diff, include_unrel=impute_unrel,
                                robust=robust, trios_sibs=trios_sibs, verbose=False,
                                print_sample_info=print_sample_info, matched_variants=matched_variants)
    else:
        G = reader.read(matched_variants)
        pheno_ids = worker['model_ids']
    if not np.array_equal(G.ids, pheno_ids):
        G.filter_ids(pheno_ids)
        if not np.array_equal(G.ids, pheno_ids):
            raise ValueError('Genotype batch does not contain every individual in the analysis cohort')
    if reader is None and not impute_unrel and np.any(G.fams == ''):
        raise ValueError('No family label from pedigree for some individuals')
    G.compute_freqs()
    if impute_unrel and G.sid.size:
        G = impute_unrel_par_gts(G, sib=fit_sib, parsum=parsum, ped=ped)
    maf_pass = G.filter_maf(min_maf, verbose=verbose)
    missing_pass = G.filter_missingness(max_missing, verbose=verbose)
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
        sibling_reader = worker.get('sibling_reader')
        if sibling_reader is None:
            G_sib = read.get_gts_matrix(ped=ped, imp_fams=imp_fams, bedfile=bedfile, bgenfile=bgenfile,
                                        snp_ids=G.sid, ids=pheno_ids, sib_diff=True, match_snp_ids=True)
        else:
            tested = matched_variants.subset(maf_pass).subset(missing_pass)
            G_sib = sibling_reader.read(tested)
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
                       add_jitter=False, standard_gwas=False, batch_size=None, *,
                       threads=1, family_labels=None, solver='lu', inference_model=None):
    if (bedfile is None) == (bgenfile is None):
        raise ValueError('Must supply exactly one BED or BGEN genotype file')
    if cpus < 1 or threads < 1:
        raise ValueError('cpus and threads must be positive')
    threads = min(threads, available_cpus())
    cpus = min(cpus, max(1, available_cpus() // threads))
    options = dict(bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, parsum=parsum,
                   sib=fit_sib, sib_diff=sib_diff, include_unrel=impute_unrel, robust=robust,
                   trios_sibs=trios_sibs, num_threads=threads)
    with ExitStack() as stack:
        stack.enter_context(limit_native_threads(threads))
        ids = None if y is None else y.ids
        reader = stack.enter_context(read.prepare_gts_matrix(ped=ped, imp_fams=imp_fams, ids=ids,
                                                           **options))
        observed = reader.observed_variants
        print(f'Found {observed.sid.size} SNPs')
        try:
            matched = reader.match_variants(match_snp_ids=True)
        except ValueError as error:
            label = chrom_out or ','.join(np.unique(observed.chromosome).astype(str))
            location = f'chromosome {label}, ' if label and label != '0' else ''
            raise ValueError(f'{error} ({location}{bedfile or bgenfile}, {par_gts_f})') from error
        # BGEN string metadata can arrive as Python objects. Fixed-width
        # strings can be shared across spawned workers without object pointers.
        matched = MatchedVariants(*(a.astype(str) if a.dtype.hasobject else a for a in matched))
        if par_gts_f is not None and not sib_diff:
            print(f'Matched {matched.sid.size} SNPs between observed and imputed genotypes')
        _, index, counts = np.unique(observed.sid, return_inverse=True, return_counts=True)
        retain = counts[index] == 1
        if not np.all(retain):
            print(f'Removing {np.count_nonzero(counts > 1)} duplicate SNP ids')
            observed = observed.subset(retain)
        chrom, snp_ids, pos, alleles = observed[:4]
        if bgenfile is not None:
            chrom = np.where(chrom == '', str(chrom_out), chrom)
        snp_dict = make_id_dict(snp_ids)

        plans = reader.sample_plans
        trios_only = trios_sibs and (plans[0].trio_count == len(plans[0].ids))
        if not impute_unrel and np.any(plans[0].fams == ''):
            raise ValueError('No family label from pedigree for some individuals')
        if reader.parsum and not parsum:
            print('No observed parental genotypes found. Fitting the sum of parental genotypes.')
        parsum = reader.parsum
        options['parsum'] = parsum
        sibling_reader = sibling_plans = None
        if robust:
            sibling_opts = dict(options, par_gts_f=None, sib_diff=True, robust=False)
            sibling_reader = stack.enter_context(read.prepare_gts_matrix(ped=ped, ids=ids, **sibling_opts))
            sibling_plans = sibling_reader.sample_plans

        if robust or sib_diff or standard_gwas or (trios_sibs and not trios_only):
            alpha_dim = 1
        else:
            alpha_dim = 2 + fit_sib + (not parsum)
        # Untested SNPs remain present in the original file order.
        alpha = np.full((snp_ids.size, alpha_dim), np.nan, dtype=np.float32)
        alpha_cov = np.full((snp_ids.size, alpha_dim, alpha_dim), np.nan, dtype=np.float32)
        alpha_ses = np.full_like(alpha, np.nan)
        freqs = np.full(snp_ids.size, np.nan, dtype=np.float32)

        covar = None if covariates is None else covariates.gts.data
        response = y.gts.data.reshape(-1)
        n = len(response)
        structured = family_labels is not None and not len(varcomp_lst)
        use_threads = structured and bedfile is not None and not robust and (not trios_sibs or trios_only)
        # Allow also for filtering/imputation copies and the dense RHS blocks of a sparse solve.
        per_variant = reader.bytes_per_variant + (10 if structured else 26) * reader.n_components * n
        reader_fixed = reader.fixed_read_bytes
        if sibling_reader is not None:
            per_variant += sibling_reader.bytes_per_variant
            reader_fixed += sibling_reader.fixed_read_bytes
        model_bytes = n * (64 + 64 * robust + (0 if covar is None else 24 * covar.shape[1]))
        model_bytes += sum(len(component[0]) * 72 for component in varcomp_lst)
        worker_overhead = 0 if use_threads or cpus == 1 else 96 * MIB
        shared_bytes = 0 if use_threads or cpus == 1 else sum(
            array.nbytes for plan in (*plans, *(sibling_plans or ())) for array in plan
            if isinstance(array, np.ndarray))
        if shared_bytes:
            shared_bytes += response.nbytes + ids.nbytes + (0 if covar is None else covar.nbytes)
            shared_bytes += sum(sum(array.nbytes for array in component) for component in varcomp_lst)
            shared_bytes += sum(array.nbytes for array in matched)
        if batch_size is None:
            batch_size = max(1, (256 * MIB) // per_variant)
        try:
            workers, width = plan_batches(
                len(matched.sid), per_variant, workers=cpus, batch_size=batch_size,
                fixed_bytes_per_worker=reader_fixed + worker_overhead + (0 if use_threads else model_bytes),
                fixed_bytes=shared_bytes + (model_bytes if use_threads and inference_model is None else 0),
            )
        except MemoryError:
            if not shared_bytes:
                raise
            # A spawn pool's shared serialization buffers can make an otherwise
            # feasible chromosome exceed its budget. Use the in-process reader.
            workers, width = plan_batches(len(matched.sid), per_variant, batch_size=batch_size,
                                         fixed_bytes_per_worker=reader_fixed + model_bytes)
        batch_size = width
        total = _number_of_batches(len(matched.sid), workers, batch_size)
        batches = _batch_slices(len(matched.sid), workers, batch_size)
        print(f"Using {total} {'batch' if total == 1 else 'batches'} with up to {batch_size} variants "
              f'and {workers} {"worker" if workers == 1 else "workers"}')

        processor = partial(process_batch, ped=None, n_vararr=len(varcomp_lst), imp_fams=None,
                            parsum=parsum, fit_sib=fit_sib, sib_diff=sib_diff, standard_gwas=standard_gwas,
                            max_missing=max_missing, min_maf=min_maf, impute_unrel=impute_unrel,
                            robust=robust, trios_sibs=trios_sibs)

        def store_results(results):
            for i, (batch_freqs, batch_snps, batch_alpha, batch_alpha_cov, batch_alpha_ses) in enumerate(results, 1):
                if len(batch_snps):
                    indices = np.fromiter((snp_dict[s] for s in batch_snps), np.intp, count=len(batch_snps))
                    alpha[indices] = batch_alpha
                    alpha_cov[indices] = batch_alpha_cov
                    alpha_ses[indices] = batch_alpha_ses
                    freqs[indices] = batch_freqs
                if total <= 50 or i == 1 or i == total or i % (total // 10) == 0:
                    print(f'Done batch {i} out of {total} #snps:{len(batch_snps)}')

        start = time.perf_counter()
        if workers == 1 or debug or use_threads:
            model = inference_model
            if model is None:
                model = _fit_model(response, varcomp_lst, sigmas, covar, family_labels, add_jitter, solver)
            if workers == 1 or debug:
                state = dict(model=model, reader=reader, sibling_reader=sibling_reader, model_ids=ids)
                def task(batch):
                    return processor(matched.subset(batch), worker_state=state)

                if debug:
                    task(next(batches))
                    raise SystemExit('Debug finished.')
                store_results(map(task, batches))
            else:
                # Populate shared read-only model state before any thread enters a solve.
                model._family_snp_terms
                state = local()
                lock = Lock()

                def start_thread():
                    from numba import config, set_num_threads

                    set_num_threads(min(threads, config.NUMBA_NUM_THREADS))
                    # Each worker owns its BED/HDF5 handles, using the same numeric sample plan.
                    with lock:
                        private_reader, _ = _readers(stack, options, plans, None)
                    state.values = dict(model=model, reader=private_reader, model_ids=ids)

                def run_thread(batch):
                    return processor(matched.subset(batch), worker_state=state.values)

                with ThreadPoolExecutor(max_workers=workers, initializer=start_thread) as pool:
                    store_results(_bounded_results(pool, run_thread, batches, workers))
        else:
            # Spawn keeps native thread pools and open HDF5 handles out of child state.
            labels = family_labels
            if structured:
                labels, _ = pd.factorize(family_labels, sort=False)
            shared = tuple(tuple(_share_array(np.asarray(array)) for array in component)
                           for component in varcomp_lst)
            initializer_args = (options, _share_plans(plans), _share_plans(sibling_plans),
                                _share_array(response), shared, sigmas, _share_array(covar),
                                _share_array(labels), add_jitter, threads, processor, solver, _share_array(ids),
                                tuple(_share_array(array) for array in matched))
            with worker_environment(threads), ProcessPoolExecutor(
                max_workers=workers, mp_context=get_context('spawn'), initializer=_start_process_worker,
                initargs=initializer_args,
            ) as pool:
                store_results(_bounded_results(pool, _run_worker_batch, batches, workers))
    print('Time used for inference: ', f'{time.perf_counter() - start:.2f}s')
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
