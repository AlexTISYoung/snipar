import snipar.preprocess as preprocess
from snipar.utilities import open_bgen
import snipar.read.bed as bed
import snipar.read.bgen as bgen
import snipar.read.phenotype as phenotype
import h5py
import numpy as np
import pandas as pd
from snipar.utilities import make_id_dict
from typing import Tuple
from contextlib import contextmanager, ExitStack
from functools import cached_property

from bed_reader import open_bed
from numba import config, get_num_threads, set_num_threads

from snipar.gtarray import gtarray
from snipar.read.family import IdIndex, _sample_plan, _assemble, _assemble_parallel, _read_imputed
from snipar.read.variants import MatchedVariants, NoMatchingVariantsError, bgen_variant_ids
from snipar.read.variants import match_variants as align_imputed_variants
from snipar.read.variants import observed_variant_indices

def build_ped_from_par_gts(par_gts_f, *, include_controls=False):
    """Read the pedigree and family-row identifiers stored with imputed genotypes."""
    with h5py.File(par_gts_f, 'r') as source:
        ped = source['pedigree'].asstr(encoding='utf-8')[1:].astype(str)
        families = source['families'].asstr(encoding='utf-8')[:].astype(str)
    if not include_controls:
        ped = ped[~np.char.startswith(ped[:, 0], '_')]
    return ped, families


def get_gts_matrix(ped=None, imp_fams=None, bedfile=None, bgenfile=None, par_gts_f=None, snp_ids = None, ids = None, parsum=False, sib = False, sib_diff = False, compute_controls = False, include_unrel=False, robust=False, trios_sibs=False, match_snp_ids=False, verbose=False, print_sample_info=False, matched_variants=None):
    """Reads observed and imputed genotypes and constructs a family based genotype matrix for the individuals with
    observed/imputed parental genotypes, and if sib=True, at least one genotyped sibling.

    Args:
        par_gts_f : :class:`str`
            path to HDF5 file with imputed parental genotypes
        gts_f : :class:`str`
            path to bed file with observed genotypes
        snp_ids : :class:`numpy.ndarray`
            If provided, only obtains the subset of SNPs specified that are present in both imputed and observed genotypes
        ids : :class:`numpy.ndarray`
            If provided, only obtains the ids with observed genotypes and imputed/observed parental genotypes (and observed sibling genotypes if sib=True)
        sib : :class:`bool`
            Retrieve genotypes for individuals with at least one genotyped sibling along with the average of their siblings' genotypes and observed/imputed parental genotypes. Default False.
        compute_controls : :class:`bool`
            Compute polygenic scores for control families (families with observed parental genotypes set to missing). Default False.
        parsum : :class:`bool`
            Return the sum of maternal and paternal observed/imputed genotypes rather than separate maternal/paternal genotypes. Default False.

    Returns:
        G : :class:`snipar.gtarray`
            Genotype array for the subset of genotyped individuals with complete imputed/obsereved parental genotypes. The array is [N x k x L], where
            N is the number of individuals; k depends on whether sib=True and whether parsum=True; and  L is the number of SNPs. If sib=False and parsum=False,
            then k=3 and this axis indexes individual's genotypes, individual's father's imputed/observed genotypes, individual's mother's imputed/observed genotypes.
            If sib=True and parsum=False, then k=4, and this axis indexes the individual, the sibling, the paternal, and maternal genotypes in that order. If parsum=True and sib=False,
            then k=2, and this axis indexes the individual and sum of paternal and maternal genotypes; etc.
            If compute_controls=True, then a list is returned, where the first element is as above, and the following elements give equivalent genotyping arrays for control families where the mother has been set
            to missing, the father has been set to missing, and both parents have been set to missing.

    """
    with prepare_gts_matrix(ped=ped, imp_fams=imp_fams, bedfile=bedfile, bgenfile=bgenfile,
                            par_gts_f=par_gts_f, ids=ids, parsum=parsum, sib=sib, sib_diff=sib_diff,
                            compute_controls=compute_controls, include_unrel=include_unrel, robust=robust,
                            trios_sibs=trios_sibs, verbose=verbose or print_sample_info) as reader:
        variants = matched_variants
        if variants is None:
            variants = reader.match_variants(snp_ids, match_snp_ids=match_snp_ids)
        return reader.read(variants)


@contextmanager
def prepare_gts_matrix(ped=None, imp_fams=None, bedfile=None, bgenfile=None, par_gts_f=None,
                       ids=None, parsum=False, sib=False, sib_diff=False, compute_controls=False,
                       include_unrel=False, robust=False, trios_sibs=False, verbose=False,
                       num_threads=None, sample_plans=None, *, _split_pedigree=True):
    """Open a chromosome and reuse its sample and variant alignments across reads.

    Call ``match_variants`` once and pass slices of its result to ``read``.
    Independent workers for the same chromosome can share ``sample_plans`` and
    open their own file handles. Plans refer to row positions on these files.
    """
    if (bedfile is None) == (bgenfile is None):
        raise ValueError('Must provide exactly one BED or BGEN file')
    if num_threads is not None and num_threads < 1:
        raise ValueError('Genotype reader threads must be positive')
    if sum((sib, sib_diff, trios_sibs)) > 1:
        raise ValueError('Select one sibling or combined trio-and-sibling family design')
    with ExitStack() as stack:
        source = par_gts_f
        if source is not None and not isinstance(source, h5py.Group):
            source = stack.enter_context(h5py.File(source, 'r'))
        if bedfile is not None:
            filename = str(bedfile)
            if not filename.lower().endswith('.bed'):
                filename += '.bed'
            observed = stack.enter_context(open_bed(
                filename, count_A1=True, num_threads=num_threads,
                iid_count=None if sample_plans is None else sample_plans[0].file_n_observed,
            ))
        else:
            filename = None
            observed = stack.enter_context(open_bgen(bgenfile, verbose=verbose))
        plan_args = (ped, imp_fams, ids, compute_controls, _split_pedigree)
        plan_options = dict(parsum=parsum, sib=sib, sib_diff=sib_diff, include_unrel=include_unrel,
                            trios_sibs=trios_sibs, verbose=verbose)
        yield _PreparedFamilyReader(observed, source, sample_plans, plan_args, plan_options,
                                    filename, robust, num_threads)


class _PreparedFamilyReader:
    def __init__(self, observed, imputed, plans, plan_args, options, bedfile, robust, num_threads):
        self._observed, self._imputation = observed, imputed
        self._plans, self._plan_args, self._plan_options = plans, plan_args, options
        self._bedfile = bedfile
        self._is_bed = bedfile is not None
        self._robust, self._num_threads = robust, num_threads
        sib_diff, trios_sibs = options['sib_diff'], options['trios_sibs']
        if plans is not None:
            sib_diff, trios_sibs = plans[0].sib_diff, plans[0].trios_sibs
        self._uses_imputation = imputed is not None and not (sib_diff or trios_sibs)
        if robust and not self._uses_imputation:
            raise ValueError('Observed parental allele counts require accompanying imputation files')

    @cached_property
    def _planning(self):
        ped, imp_fams, ids, compute_controls, split_pedigree = self._plan_args
        if self._plans is not None:
            return self._plans, [ped] * len(self._plans)
        source = self._imputation
        if ped is None or compute_controls:
            if source is None:
                reason = 'Control genotypes require an imputation file with control families' if compute_controls else (
                    'Must provide a pedigree or an imputed parental genotypes file')
                raise ValueError(reason)
            ped = source['pedigree'].asstr(encoding='utf-8')[1:].astype(str)
        if not self._uses_imputation:
            imp_fams = None
        elif imp_fams is None:
            imp_fams = source['families'].asstr(encoding='utf-8')[:].astype(str)
        pedigrees = [ped]
        if split_pedigree:
            noncontrols = ~np.char.startswith(ped[:, 0], '_')
            if not noncontrols.all():
                pedigrees[0] = ped[noncontrols]
            if compute_controls:
                for prefix in ('_p_', '_m_', '_o_'):
                    group = np.char.startswith(ped[:, 0], prefix)
                    if not np.any(group):
                        raise ValueError('Imputation file must contain paternal, maternal, and both-parent control families')
                    pedigrees.append(ped[group])
        observed_ids = (bed._read_iids(self._bedfile, self._observed) if self._is_bed
                        else np.asarray(self._observed.samples, dtype=str))
        lookup = IdIndex(observed_ids)
        plans = tuple(_sample_plan(group, observed_ids, lookup, imp_fams, ids, **self._plan_options)
                      for group in pedigrees)
        return plans, pedigrees

    @property
    def sample_plans(self):
        return self._planning[0]

    @cached_property
    def group_shapes(self):
        return tuple((len(p.ids), p.n_components) for p in self.sample_plans)

    @property
    def n_samples(self):
        return self.group_shapes[0][0]

    @property
    def n_components(self):
        return self.group_shapes[0][1]

    @cached_property
    def n_observed(self):
        return max(len(p.observed) for p in self.sample_plans)

    @cached_property
    def n_imputed(self):
        return max(len(p.imputed) for p in self.sample_plans)

    @property
    def parsum(self):
        return self.sample_plans[0].parsum

    @cached_property
    def bytes_per_variant(self):
        # Retained group matrices, masks, decoding buffers, and robust-count
        # conversion. The group inputs are read one at a time.
        decode_bytes = 1 if self._is_bed else 24
        retained = sum(n * (5 * k + (2 if self._robust else 0)) for n, k in self.group_shapes)
        return retained + max(len(p.observed) * decode_bytes + len(p.imputed) * (8 if self._robust else 4)
                              + (4 * len(p.ids) if self._robust else 0) for p in self.sample_plans)

    @cached_property
    def fixed_read_bytes(self):
        size = (8 << 20) if self._uses_imputation else 0
        if not self._is_bed:
            reader = self._observed
            n_file = reader.nsamples if hasattr(reader, 'nsamples') else len(reader.samples)
            size += 32 * n_file * (self._num_threads or 1)
        return size

    @cached_property
    def _imputed_gts(self):
        return self._imputation['imputed_par_gts'] if self._uses_imputation else None

    @cached_property
    def _num_al(self):
        return self._imputation['num_observed_parental_alleles'] if self._robust else None

    @cached_property
    def _metadata(self):
        source = self._observed
        if self._is_bed:
            sid, aliases = source.sid, None
            chrom = source.chromosome
            try:
                chrom = np.where(chrom == 'X', '23', chrom)
                for label, code in (('Y', '24'), ('XY', '25'), ('MT', '26')):
                    chrom = np.where(chrom == label, code, chrom)
                chrom = chrom.astype(float).astype(int)
            except ValueError:
                chrom = source.chromosome
            positions = source.bp_position
            alleles = np.column_stack((source.allele_1, source.allele_2))
        else:
            sid, aliases = bgen_variant_ids(source)
            sid = sid.copy()
            aliases = None if aliases is None else aliases.copy()
            chrom, positions = np.array(source.chromosomes), np.array(source.positions)
            alleles = np.array([pair.split(',') for pair in source.allele_ids])
        indices = np.arange(len(sid), dtype=np.intp)
        false = np.broadcast_to(False, len(sid))
        absent = np.broadcast_to(-1, len(sid))
        return MatchedVariants(chrom, sid, positions, alleles, false, absent, indices), aliases

    @property
    def observed_variants(self):
        """Full observed metadata, including duplicate IDs, in file order."""
        return self._metadata[0]

    def match_variants(self, snp_ids=None, *, match_snp_ids=False, start=0, end=None):
        observed = self.observed_variants
        if self._uses_imputation:
            return align_imputed_variants(observed.sid, observed.chromosome, observed.pos, observed.alleles,
                                          self._imputation, snp_ids, start, end, self._metadata[1])
        rows = observed_variant_indices(observed.sid, snp_ids, match_snp_ids)
        if not rows.size:
            raise NoMatchingVariantsError('No SNPs found in observed genotypes')
        return observed.subset(rows)

    def _read_group(self, variants, plan, pedigree, order):
        n_variants = len(variants.sid)
        samples = plan.observed
        if samples[-1] - samples[0] + 1 == len(samples):
            samples = slice(int(samples[0]), int(samples[-1]) + 1)
        if self._is_bed:
            observed = self._observed.read((samples, variants.observed), dtype=np.int8, order=order)
        else:
            observed = bgen._first_allele_dosage(self._observed, samples, variants.observed,
                                                num_threads=self._num_threads, order=order)
        imputed = np.empty((0, n_variants), np.float32)
        if plan.imputed.size:
            if self._imputed_gts is None or np.any(variants.imputed < 0):
                raise ValueError('Selected variants have no aligned imputed genotypes')
            imputed = _read_imputed(self._imputed_gts, plan.imputed, variants.imputed)
        num_al = _read_imputed(self._num_al, plan.imputed, variants.imputed) if self._robust else imputed[:0]
        counts = np.empty((len(plan.ids) if self._robust else 0, n_variants), np.float32)
        output = np.empty((len(plan.ids), plan.n_components, n_variants), np.float32, order=order)
        params = (output, observed, imputed, plan.genotype_rows, plan.parent_status, plan.sibling_members,
                  plan.member_offsets, plan.sibling_targets, plan.target_offsets, plan.parsum and not plan.trios_sibs,
                  plan.sib, plan.sib_diff, plan.trios_sibs, variants.flipped, self._robust, counts, num_al)
        if self._num_threads is not None and self._num_threads > 1 and output.size >= 1_000_000:
            previous = get_num_threads()
            set_num_threads(min(self._num_threads, config.NUMBA_NUM_THREADS))
            try:
                _assemble_parallel(*params)
            finally:
                set_num_threads(previous)
        else:
            _assemble(*params)
        result = gtarray(output, plan.ids, variants.sid, chrom=variants.chromosome, pos=variants.pos,
                         alleles=variants.alleles, fams=plan.fams, ped=pedigree,
                         par_status=None if plan.sib_diff else plan.parent_status,
                         num_obs_par_al=counts.astype(np.float16) if self._robust else None)
        if plan.trios_sibs:
            result.complete_trios_inds = np.arange(plan.trio_count)
            result.sibs_inds = np.arange(plan.trio_count, len(plan.ids))
        return result

    def read(self, variants, *, order='C'):
        """Build family-genotype matrices for a block of matched variants."""
        plans, pedigrees = self._planning
        results = [self._read_group(variants, plan, ped, order) for plan, ped in zip(plans, pedigrees)]
        return results[0] if len(results) == 1 else results


def _exclude_duplicate_twins(ids, ibdrel_path):
    if ibdrel_path is None:
        return ids
    with open(ibdrel_path) as source:
        if 'InfType' not in source.readline().split():
            return ids
    remove = set()
    n_pairs = 0
    with pd.read_csv(ibdrel_path, sep=r'\s+', usecols=['ID1', 'InfType'], dtype=str,
                     keep_default_na=False, chunksize=200_000) as chunks:
        for chunk in chunks:
            pairs = chunk.loc[chunk['InfType'] == 'Dup/MZ', 'ID1']
            n_pairs += len(pairs)
            remove.update(pairs)
    if not n_pairs:
        return ids
    keep = ~np.isin(ids, list(remove))
    print(f'WARNING: {n_pairs} pairs of Duplicates or MZ twins in pedigree. '
          f'Removed {np.count_nonzero(~keep)} individuals.')
    return ids[keep]


def _observed_ids(filename):
    filename = str(filename)
    if filename.endswith('.bed'):
        with open_bed(filename, count_A1=True) as reader:
            return bed._read_iids(filename, reader)
    if filename.endswith('.bgen'):
        with open_bgen(filename) as reader:
            return np.array(reader.samples, dtype=str)
    raise ValueError(f'Unknown filetype for observed genotypes file: {filename}')


def _fill_unrelated_labels(fam_labels):
    missing = fam_labels == ''
    if not np.any(missing):
        return fam_labels
    # Use labels that cannot collide with genuine pedigree families.
    used = set(fam_labels[~missing])
    replacements = []
    index = 0
    for _ in range(np.count_nonzero(missing)):
        label = f'_not_{index}_'
        while label in used:
            index += 1
            label = f'_not_{index}_'
        replacements.append(label)
        index += 1
    replacements = np.asarray(replacements)
    result = fam_labels.astype(np.promote_types(fam_labels.dtype, replacements.dtype))
    result[missing] = replacements
    return result


def get_ids_with_par(gts_f: str,
                     ped: np.ndarray,
                     imp_fams: np.ndarray,
                     ids: np.ndarray = None,
                     sib: bool = False,
                     include_unrel: bool = False,
                     ibdrel_path: str = None,
                     return_info: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Find ids with observed/imputed parents and family labels.
    """
    gts_ids = _observed_ids(gts_f)
    ids = _exclude_duplicate_twins(gts_ids if ids is None else np.asarray(ids), ibdrel_path)
    if return_info:
        ids, observed_indices, imp_indices, _ = preprocess.get_indices_given_ped(
            ped, gts_ids, imp_fams=imp_fams, ids=ids, sib=sib, include_unrel=include_unrel,
        )
        gts_id_dict = make_id_dict(gts_ids[observed_indices])
        if imp_fams is not None:
            imp_fams = imp_fams[imp_indices]
        par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict, imp_fams)
    else:
        ids, fam_labels, par_status, _ = _sample_plan(
            ped, gts_ids, IdIndex(gts_ids), imp_fams, ids, sib=sib, include_unrel=include_unrel, select_only=True,
        )

    fam_labels = _fill_unrelated_labels(fam_labels)

    missing = np.count_nonzero(par_status < 0, axis=1)
    with_both = np.count_nonzero(missing == 0)
    imputed = np.count_nonzero(par_status == 1, axis=1)
    print(f'{with_both} individuals with phenotype observations and complete observed/imputed genotype observations.')
    print(f'{np.count_nonzero(imputed == 2)} individuals with imputed but no observed parental genotypes.')
    print(f'{np.count_nonzero((missing == 0) & (imputed == 1))} individuals with one imputed and one observed parental genotypes.')
    print(f'{np.count_nonzero((missing == 0) & (imputed == 0))} individuals with both parents observed.')

    if include_unrel:
        unrelated_inds = np.flatnonzero(missing == 2)
        print(f'{len(unrelated_inds)} samples without imputed or observed parental genotypes will be included through linear imputation.')
    one_missing = np.count_nonzero(missing == 1)
    if one_missing:
        raise RuntimeError(f'{one_missing} individuals have one missing parent and one non-missing parent.')
    
    if return_info:
        if include_unrel:
            return ids, fam_labels, unrelated_inds
        return ids, fam_labels, par_status, gt_indices
    return ids, fam_labels


def get_ids_with_sibs(gts_f: str,
                      ped: np.ndarray, 
                      ids: np.ndarray = None,
                      ibdrel_path: str = None,
                      return_info: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Find ids with sibs and family labels.
    """
    gts_ids = _observed_ids(gts_f)
    ids = _exclude_duplicate_twins(gts_ids if ids is None else np.asarray(ids), ibdrel_path)
    if not return_info:
        ids, fam_labels, _, _ = _sample_plan(ped, gts_ids, IdIndex(gts_ids), None, ids,
                                              sib_diff=True, select_only=True)
        print(f'Using {len(ids)} samples for the sib-difference estimator')
        return ids, _fill_unrelated_labels(fam_labels)
    ids, observed_indices = preprocess.get_indices_given_ped_sibs(ped, gts_ids, ids=ids, verbose=False)
    gts_ids = gts_ids[observed_indices]


    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    gt_indices, fam_labels = preprocess.find_gts(ids, ped, gts_id_dict)
    # Find which individuals can be used
    none_missing = gt_indices >= 0
    N = np.sum(none_missing)
    if N == 0:
        raise ValueError(
            'No individuals with phenotype observations and genotypes.')
    # Take those that can be used
    # gt_indices = gt_indices[none_missing]
    ids = ids[none_missing]
    fam_labels = _fill_unrelated_labels(fam_labels[none_missing])

    print(f'Using {len(ids)} samples for the sib-difference estimator')
    return ids, fam_labels, gt_indices

def get_ids_with_trios_sibs(gts_f: str,
                            ped: np.ndarray, 
                            ids: np.ndarray = None,
                            ibdrel_path: str = None,
                            return_info: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Find ids with sibs and both parents, and their family labels.
    """
    gts_ids = _observed_ids(gts_f)
    ids = _exclude_duplicate_twins(gts_ids if ids is None else np.asarray(ids), ibdrel_path)
    if not return_info:
        ids, fam_labels, _, n_trios = _sample_plan(
            ped, gts_ids, IdIndex(gts_ids), None, ids, trios_sibs=True, verbose=True, select_only=True,
        )
        return ids, _fill_unrelated_labels(fam_labels), n_trios == len(ids)
    ids, observed_indices, _, _ = preprocess.get_indices_given_ped_trios_sibs(
        ped, gts_ids, ids=ids, verbose=False,
    )
    gts_ids = gts_ids[observed_indices]


    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict)
    # Find which individuals can be used
    none_missing = gt_indices[:, 0] >= 0
    N = np.sum(none_missing)
    if N == 0:
        raise ValueError(
            'No individuals with phenotype observations and genotypes.')
    # Take those that can be used
    gt_indices = gt_indices[none_missing]
    ids = ids[none_missing]
    fam_labels = _fill_unrelated_labels(fam_labels[none_missing])

    return ids, fam_labels, gt_indices
    
def get_ids_with_trios(gts_f: str,
                       ped: np.ndarray, 
                       ids: np.ndarray = None,
                       ibdrel_path: str = None,
                       return_info: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Remove ids with sibs only. Used for the unified estimator when only pedigree file is supplied.
    """
    gts_ids = _observed_ids(gts_f)
    ids = _exclude_duplicate_twins(gts_ids if ids is None else np.asarray(ids), ibdrel_path)
    ids, observed_indices = preprocess.get_indices_given_ped_trios(ped, gts_ids, ids=ids, verbose=False)
    gts_ids = gts_ids[observed_indices]


    gts_id_dict = make_id_dict(gts_ids)

    # Find indices in reduced data
    par_status, gt_indices, fam_labels = preprocess.find_par_gts(ids, ped, gts_id_dict)
    # TODO: might suffice to only check if the first column == 0?
    # Find which individuals can be used
    none_missing = np.sum(gt_indices == -1, axis=1)
    none_missing = none_missing == 0
    N = sum(none_missing)
    if N == 0:
        raise ValueError(
            'No individuals with phenotype observations, genotypes and complete parental genotypes.')
    # Take those that can be used
    gt_indices = gt_indices[none_missing]
    ids = ids[none_missing]
    fam_labels = _fill_unrelated_labels(fam_labels[none_missing])

    if return_info:
        return ids, fam_labels, gt_indices
    return ids, fam_labels
