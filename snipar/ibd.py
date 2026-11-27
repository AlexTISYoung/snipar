from contextlib import ExitStack
import gzip
from pathlib import Path
from shutil import disk_usage
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from numba import get_num_threads, njit, prange
from snipar.map import decode_map_from_pos, get_map_positions
from snipar.ld import compute_ld_scores_packed, pack_ld_hardcalls
import numpy as np
from snipar.observed import allele_stats, open_observed
from snipar.read.family import IdIndex
from snipar.resources import plan_batches, get_memory_budget
from snipar.utilities import outfile_name

####### Transition Matrix ######
@njit
def transition_matrix(cM):
    """Compute probabilities of transitioning between IBD states as a function of genetic distance.
    Args:
        cM : :class:`float`
            genetic distance in centiMorgans (cM)
    Returns:
        log(P) : :class:`~numpy:numpy.array`
            the natural logarithm (element-wise) of the matrix of transition probabilities
    """
    P = np.identity(3)
    r_array = np.array([[-1.0,1.0,0.0],[0.5,-1.0,0.5],[0.0,1.0,-1.0]],dtype=np.float64)
    r = 1.0 - np.power((1.0 + np.exp(-cM / 50)) / 2.0, 4)
    P += r_array*r
    return np.log(P)

#
@njit
def p_ibd_0(f):
    """Compute Joint-PMF for sibling pair genotypes given IBD0.
    Args:
        f : :class:`float`
            allele frequency
    Returns:
        P : :class:`~numpy:numpy.array`
            matrix of probabilities
    """
    P_vec = np.array([(1.0-f)**2.0,2.0*f*(1-f),f**2.0]).reshape((3,1))
    return P_vec @ P_vec.T

@njit
def p_ibd_1(f):
    """Compute Joint-PMF for sibling pair genotypes given IBD1.
    Args:
        f : :class:`float`
            allele frequency
    Returns:
        P : :class:`~numpy:numpy.array`
            matrix of probabilities
    """
    q = 1 - f
    return np.array([[q**3, f*q*q, 0.0],
                     [f*q*q, f*q, f*f*q],
                     [0.0, f*f*q, f**3]])

@njit
def p_ibd_2(f):
    """Compute Joint-PMF for sibling pair genotypes given IBD2.
    Args:
        f : :class:`float`
            allele frequency
    Returns:
        P : :class:`~numpy:numpy.array`
            matrix of probabilities
    """
    P = np.zeros((3,3),dtype=np.float64)
    fdiag = np.array([(1.0-f)**2.0,2.0*f*(1.0-f),f**2.0])
    np.fill_diagonal(P,fdiag.reshape((3,1)))
    return P

@njit
def p_obs_given_IBD(g1_obs,g2_obs,f,p):
    """Compute Joint-PMF for sibling pair genotypes given IBD0.
    Args:
        g1_obs : :class:`integer`
            observed genotype for sibling 1
        g2_obs : :class:`integer`
            observed genotype for sibling 2
        f : :class:`float`
            allele frequency
        p : :class:'float'
            genotyping error probability
    Returns:
        log(P) : :class:`~numpy:numpy.array`
            vector giving log-probabilities of observing g1_obs,g2_obs give IBD 0,1,2
    """
    P_out = np.zeros((3))
    # Get probabilities of observed genotypes given true
    err_vectors = np.array([[1.0 - p, p / 2.0, 0.0],
                            [p, 1.0 - p, p],
                            [0.0,   p / 2.0, 1.0 - p]])
    # Get probabilities of true genotypes given IBD
    probs = np.zeros((3, 3, 3))
    probs[0, ...] = p_ibd_0(f)
    probs[1, ...] = p_ibd_1(f)
    probs[2, ...] = p_ibd_2(f)
    # Compute probabilities conditional on IBD
    for i in range(3):
        P_out[i] = err_vectors[g1_obs, ...].T @ probs[i,...] @ err_vectors[g2_obs, ...]
    # Return
    return np.log(P_out)

@njit
def make_dynamic(g1, g2, freqs, map, weights, error_probs):
    """Make state-matrix and pointer matrix for a sibling pair by dynamic programming
    Args:
        g1 : :class:`~numpy:numpy.array`
            integer vector of first sibling's genotypes
        g2 : :class:`~numpy:numpy.array`
            integer vector of first sibling's genotypes
        freqs : :class:`~numpy:numpy.array`
            floating point vector of allele frequencies
        map : :class:`~numpy:numpy.array`
            floating point vector of genetic positions in cM
        weights : :class:`~numpy:numpy.array`
            floating point vector of SNP weights (usually inverse LD-scores)
        p : :class:'float'
            genotyping error probability
    Returns:
        state_matrix : :class:`~numpy:numpy.array`
            matrix where each column gives the prob of max prob path to that state, where each row is IBD 0,1,2
        pointers : :class:`~numpy:numpy.array`
            integer vector giving the pointer to which state from previous position lead to max at this position
    """
    state_matrix = np.zeros((3, g1.shape[0]), dtype=np.float64)
    pointers = np.zeros((3, g1.shape[0]), dtype=np.int8)
    if not g1.size:
        return state_matrix, pointers
    # Check for nans
    not_nan = np.logical_not(np.logical_or(np.isnan(g1), np.isnan(g2)))
    # Initialise
    state_matrix[:,0] = np.log(np.array([0.25, 0.5, 0.25],dtype=np.float64))
    if not_nan[0]:
        state_matrix[:, 0] += weights[0]*p_obs_given_IBD(np.int8(g1[0]), np.int8(g2[0]), freqs[0], error_probs[0])
    # Compute
    for l in range(1, g1.shape[0]):
        if not_nan[l]:
            probs = weights[l]*p_obs_given_IBD(np.int8(g1[l]), np.int8(g2[l]), freqs[l], error_probs[l])
        else:
            probs = np.zeros((3))
        tmatrix = transition_matrix(map[l]-map[l-1])
        for i in range(3):
            tprobs = tmatrix[:, i]+state_matrix[:, l-1]
            state_matrix[i, l] = np.max(tprobs)+probs[i]
            pointers[i, l] = np.argmax(tprobs)
    return state_matrix, pointers

@njit
def viterbi(state_matrix,pointers):
    """Get viterbi path from state_matrix and pointers output of make_dynamic
    Args:
        state_matrix : :class:`~numpy:numpy.array`
            matrix where each column gives the prob of max prob path to that state, where each row is IBD 0,1,2
        pointers : :class:`~numpy:numpy.array`
            integer vector giving the pointer to which state from previous position lead to max at this position
    Returns:
            path : :class:`~numpy:numpy.array`
            integer vector giving the Viterbi path through the IBD states
    """
    path = np.zeros(state_matrix.shape[1],dtype=np.int8)
    if not path.size:
        return path
    path[path.shape[0]-1] = np.argmax(state_matrix[:, state_matrix.shape[1]-1])
    for i in range(1,path.shape[0]):
        path[path.shape[0]-(i+1)] = pointers[path[path.shape[0]-i], path.shape[0]-i]
    return path

@njit(cache=True)
def _chromosome_hmm(freqs, genetic_map, weights, error_probs):
    n = len(freqs)
    transitions = np.empty((n, 3, 3))
    emissions = np.empty((n, 3, 3, 3))
    for marker in range(n):
        if marker:
            transitions[marker] = transition_matrix(genetic_map[marker] - genetic_map[marker - 1])
        f = freqs[marker]
        q = 1 - f
        e = error_probs[marker]
        # Conditional on the shared haplotype's allele, the unshared alleles
        # of the two siblings are independent draws from the population.
        error = np.array([[1-e, e/2, 0], [e, 1-e, e], [0, e/2, 1-e]])
        for x in range(3):
            x0 = q * error[x, 0] + f * error[x, 1]
            x1 = q * error[x, 1] + f * error[x, 2]
            for y in range(3):
                y0 = q * error[y, 0] + f * error[y, 1]
                y1 = q * error[y, 1] + f * error[y, 2]
                if weights[marker] == 0:
                    emissions[marker, x, y] = 0
                    continue
                probabilities = (
                    (q*x0 + f*x1) * (q*y0 + f*y1),
                    q*x0*y0 + f*x1*y1,
                    q*q*error[x, 0]*error[y, 0] + 2*q*f*error[x, 1]*error[y, 1]
                    + f*f*error[x, 2]*error[y, 2],
                )
                for state in range(3):
                    emissions[marker, x, y, state] = weights[marker] * np.log(probabilities[state])
    return transitions, emissions


@njit(cache=True)
def infer_ibd(sibpairs, gts, freqs, map, weights, error_probs):
    """Return unsmoothed paths for an in-memory dosage matrix."""
    if not gts.shape[1]:
        return np.empty((len(sibpairs), 0), dtype=np.int8)
    transitions, emissions = _chromosome_hmm(freqs, map, weights, error_probs)
    return _infer_prepared(sibpairs, gts, np.arange(gts.shape[1]), transitions, emissions, 0)


@njit(inline='always')
def _call_at(genotypes, row, column, layout):
    if layout == 0:
        value = genotypes[row, column]
        return int(value) if np.isfinite(value) and value >= 0 else -1
    if layout == 1:
        call = (int(genotypes[row, column // 4]) >> (2 * (column % 4))) & 3
        return -1 if call == 3 else int(call)
    call = (int(genotypes[column, row // 4]) >> (2 * (row % 4))) & 3
    return (2, -1, 1, 0)[call]


@njit(cache=True, parallel=True)
def _infer_prepared(sibpairs, gts, markers, transitions, emissions, layout):
    n = markers.size
    ibd = np.empty((len(sibpairs), n), dtype=np.int8)
    for i in prange(len(sibpairs)):
        row1, row2 = sibpairs[i, 0], sibpairs[i, 1]
        previous = np.log(np.array([.25, .5, .25]))
        current = np.empty(3)
        pointers = np.empty((n, 3), dtype=np.int8)
        call1 = _call_at(gts, row1, markers[0], layout)
        call2 = _call_at(gts, row2, markers[0], layout)
        if min(call1, call2) >= 0:
            previous += emissions[0, call1, call2]
        for marker in range(1, n):
            call1 = _call_at(gts, row1, markers[marker], layout)
            call2 = _call_at(gts, row2, markers[marker], layout)
            observed = min(call1, call2) >= 0
            for dest in range(3):
                best = 0
                maximum = previous[0] + transitions[marker, 0, dest]
                for source in range(1, 3):
                    probability = previous[source] + transitions[marker, source, dest]
                    if probability > maximum:
                        best, maximum = source, probability
                pointers[marker, dest] = best
                current[dest] = maximum
                if observed:
                    current[dest] += emissions[marker, call1, call2, dest]
            previous, current = current, previous
        state = np.argmax(previous)
        ibd[i, n - 1] = state
        for marker in range(n - 1, 0, -1):
            state = pointers[marker, state]
            ibd[i, marker - 1] = state
    return ibd


@njit(cache=True, parallel=True)
def _pack_calls(values, out):
    for row in prange(values.shape[0]):
        for byte in range(out.shape[1]):
            word = 0
            for shift in range(min(4, values.shape[1] - 4 * byte)):
                call = values[row, 4 * byte + shift]
                code = int(call) if np.isfinite(call) and call >= 0 else 3
                word |= code << (2 * shift)
            out[row, byte] = word


class segment(object):
    def __init__(self,start_index,end_index,start_bp,end_bp,start_snp,end_snp,length,state):
        self.start = start_index
        self.end = end_index
        self.start_bp = start_bp
        self.end_bp = end_bp
        self.start_snp = start_snp
        self.end_snp = end_snp
        self.length = length
        self.state = state
    def to_text(self,id1,id2,chr,end=False):
        seg_txt = str(id1)+'\t'+str(id2)+'\t'+str(self.state)+'\t'+str(chr)+'\t'+str(self.start_bp)+'\t'+str(self.end_bp)+'\t'+str(self.start_snp)+'\t'+str(self.end_snp)+'\t'+str(self.length)
        if end:
            return seg_txt
        else:
            return seg_txt+'\n'

def find_segments(path, map, snps, pos):
    if not len(path):
        return []
    changes = np.flatnonzero(path[1:] != path[:-1]) + 1
    boundaries = np.concatenate(([0], changes, [len(path)]))
    segments = []
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        last = end - 1
        segments.append(segment(start, end, pos[start], pos[last], snps[start], snps[last],
                                map[last] - map[start], path[start]))
    return segments

def smooth_segments(path, map, snps, pos, min_length):
    segments = find_segments(path, map, snps, pos)
    if len(segments) <= 1:
        return path, segments
    changed = False
    for i, current in enumerate(segments):
        if current.length >= min_length:
            continue
        if i == 0:
            replacement = segments[1].state
        elif i == len(segments) - 1:
            replacement = segments[i - 1].state
        elif segments[i - 1].state == segments[i + 1].state:
            replacement = segments[i - 1].state
        else:
            continue
        # Segment indices are half-open, physical coordinates remain inclusive.
        path[current.start:current.end] = replacement
        changed = True
    return path, find_segments(path, map, snps, pos) if changed else segments

def smooth_ibd(ibd,map,snps,pos,min_length):
    allsegs = []
    for i in range(ibd.shape[0]):
        ibd_path, segments = smooth_segments(ibd[i,:],map,snps,pos,min_length)
        ibd[i,:] = ibd_path
        allsegs.append(segments)
    return ibd, allsegs

_SEGMENT_HEADER = "ID1\tID2\tIBDType\tChr\tstart_coordinate\tstop_coordinate\tstartSNP\tstopSNP\tlength\n"


def write_segs(sibpairs, allsegs, chr, outfile):
    with gzip.open(outfile, "wt", compresslevel=3) as handle:
        handle.write(_SEGMENT_HEADER)
        for i, ((id1, id2), segments) in enumerate(zip(sibpairs, allsegs)):
            lines = [s.to_text(id1, id2, chr, end=True) for s in segments]
            if lines:
                handle.write("\n".join(lines))
                if i != len(sibpairs) - 1:
                    handle.write("\n")


def write_segs_from_matrix(ibd, sibpairs, snps, pos, map, chrom, outfile, *, return_segments=True):
    if return_segments:
        allsegs = [find_segments(path, map, snps, pos) for path in ibd]
        write_segs(sibpairs, allsegs, chrom, outfile)
        return allsegs
    # Simulated chromosomes can contain millions of short segments. Stream
    # one sibling pair at a time for callers that only need the file.
    snps, positions = list(snps), [str(p) for p in pos]
    with gzip.open(outfile, "wt", compresslevel=3) as handle:
        handle.write(_SEGMENT_HEADER)
        for i, ((id1, id2), path) in enumerate(zip(sibpairs, ibd)):
            if not len(path):
                continue
            changes = np.flatnonzero(path[1:] != path[:-1]) + 1
            starts = np.concatenate(([0], changes))
            ends = np.concatenate((changes - 1, [len(path) - 1]))
            states = path[starts].tolist()
            lengths = (map[ends] - map[starts]).tolist()
            prefix = f"{id1}\t{id2}\t"
            lines = [
                f"{prefix}{state}\t{chrom}\t{positions[a]}\t{positions[b]}\t{snps[a]}\t{snps[b]}\t{length}"
                for a, b, state, length in zip(starts.tolist(), ends.tolist(), states, lengths)
            ]
            handle.write("\n".join(lines))
            if i != len(sibpairs) - 1:
                handle.write("\n")

def _workspace(stack, scratch, name, shape, dtype, memory, order='C'):
    if np.prod(shape, dtype=np.int64) * np.dtype(dtype).itemsize <= memory:
        return np.empty(shape, dtype=dtype, order=order)
    result = np.memmap(scratch / name, mode='w+', dtype=dtype, shape=shape, order=order)
    stack.callback(result._mmap.close)
    return result


def _ibd_map(source, chrom, positions, mapfile):
    if mapfile is not None:
        return get_map_positions(mapfile, SimpleNamespace(pos=positions))
    if source.bed:
        # Keep the map in float64. The BED reader exposes a float32 cM
        # field, which can move a pair across the exact LD window boundary.
        genetic_map = np.loadtxt(Path(source.reader.location).with_suffix('.bim'),
                                 dtype=np.float64, usecols=2, ndmin=1)
        if np.any(genetic_map != 0) and (len(genetic_map) == 1 or np.any(genetic_map != genetic_map[0])):
            if (not np.isfinite(genetic_map).all() or genetic_map[0] < 0 or genetic_map[-1] > 5000
                    or np.any(np.diff(genetic_map) < 0)):
                raise ValueError('BIM genetic map must be finite, nondecreasing, and within [0, 5000] cM')
            return genetic_map
    return decode_map_from_pos(chrom, positions)


def _write_ibd_batches(pairs, row_pairs, genotypes, layout, markers, transitions, emissions,
                       snps, positions, genetic_map, min_length, chrom, outprefix, ibdmatrix,
                       pair_batch):
    with ExitStack() as stack:
        output = stack.enter_context(gzip.open(outfile_name(outprefix, '.ibd.segments.gz', chrom),
                                               'wt', compresslevel=3))
        output.write(_SEGMENT_HEADER)
        matrix = None
        if ibdmatrix:
            matrix = stack.enter_context(gzip.open(outfile_name(outprefix, '.ibdmatrix.gz', chrom),
                                                   'wt', compresslevel=3))
            matrix.write('sib1 sib2 ' + ' '.join(snps) + '\n')
        for start in range(0, len(pairs), pair_batch):
            end = min(start + pair_batch, len(pairs))
            paths = _infer_prepared(row_pairs[start:end], genotypes, markers, transitions, emissions, layout)
            for (id1, id2), path in zip(pairs[start:end], paths):
                path, segments = smooth_segments(path, genetic_map, snps, positions, min_length)
                output.writelines(seg.to_text(id1, id2, chrom) for seg in segments)
                if matrix is not None:
                    matrix.write(f'{id1} {id2} ' + ' '.join(path.astype(str)) + '\n')
            del paths


def infer_ibd_chr(sibpairs, error_prob, error_probs, outprefix, bedfile=None, bgenfile=None, chrom=None,
                  min_length=0.01, mapfile=None, ibdmatrix=False, ld_out=False, min_maf=0.01,
                  max_missing=5, max_error=0.01, num_threads=None):
    """Stream observed markers and decode IBD in memory-bounded sibling-pair blocks."""
    threads = get_num_threads() if num_threads is None else num_threads
    # Marker decoding and Numba computation execute in separate phases.
    with ExitStack() as stack:
        source = stack.enter_context(open_observed(bedfile, bgenfile, num_threads=threads))
        sids = source.sid
        reader = source.reader
        positions = reader.bp_position if source.bed else np.asarray(reader.positions)
        if chrom is None:
            chromosomes = np.unique(reader.chromosome if source.bed else reader.chromosomes)
            if len(chromosomes) != 1:
                raise ValueError('Provide observed genotypes for exactly one chromosome')
            chrom = chromosomes[0] or 0
        row_pairs = IdIndex(source.ids).lookup(sibpairs)
        present = (row_pairs >= 0).all(axis=1)
        if not np.any(present):
            raise ValueError('No genotyped sibling pairs found')
        pairs, row_pairs = sibpairs[present], row_pairs[present]
        rows = np.unique(row_pairs)
        n, m = len(rows), len(sids)
        if not m:
            raise ValueError('No observed variants on chromosome')
        print(f'Inferring IBD for {len(pairs)} sibling pairs on chromosome {chrom}')
        ld_rows = np.random.randint(0, n, size=min(n, 1000))
        memory = get_memory_budget()
        current_dir = Path(stack.enter_context(TemporaryDirectory(prefix='snipar-ibd-')))
        packed_shape = (n, (m + 3) // 4)
        ld_shape = (m, 3, (len(ld_rows) + 63) // 64)
        packed_bytes = n * packed_shape[1]
        ld_bytes = np.prod(ld_shape, dtype=np.int64) * np.dtype(np.uint64).itemsize
        in_memory_packed = packed_bytes <= memory // 5
        in_memory_ld = ld_bytes <= memory // 5
        needed = (0 if in_memory_packed else packed_bytes) + (0 if in_memory_ld else ld_bytes)
        layout = 1
        available_disk = disk_usage(current_dir).free
        if needed > available_disk:
            if source.bed and (in_memory_ld or ld_bytes < available_disk):
                # BED is already packed. This fallback rereads its variant-major
                # pages during tracing but does not require extra genotype disk.
                with open(reader.location, 'rb') as stream:
                    if stream.read(3) != b'\x6c\x1b\x01':
                        raise ValueError('IBD requires a variant-major PLINK BED file')
                genotypes = np.memmap(reader.location, mode='r', dtype=np.uint8, offset=3,
                                      shape=(m, (len(source.ids) + 3) // 4))
                stack.callback(genotypes._mmap.close)
                layout = 2
            else:
                raise OSError(f'IBD scratch needs {needed:,} bytes, with {available_disk:,} bytes available')
        if layout == 1:
            genotypes = _workspace(stack, current_dir, 'genotypes', packed_shape, np.uint8, memory // 5)
            row_pairs = np.searchsorted(rows, row_pairs)
        ld_genotypes = _workspace(stack, current_dir, 'ld', ld_shape, np.uint64, memory // 5)
        fixed = source.fixed_bytes + (packed_bytes if layout == 1 and in_memory_packed else 0)
        fixed += ld_bytes if in_memory_ld else 0
        _, variant_batch = plan_batches(m, n * (2 if source.bed else 24), fixed_bytes=fixed,
                                        memory_limit=memory)
        if layout == 1 and variant_batch < m:
            if variant_batch < 4:
                raise MemoryError('IBD packing needs enough memory to decode four variants')
            variant_batch -= variant_batch % 4
        freqs = np.empty(m)
        missing = np.empty(m)
        for start in range(0, m, variant_batch):
            end = min(start + variant_batch, m)
            values = source.read(rows, np.arange(start, end), require_calls=True)
            total, called = allele_stats(values)
            freqs[start:end] = np.where(called > 0, total / np.maximum(2 * called, 1), np.nan)
            missing[start:end] = 100 * (n - called) / n
            pack_ld_hardcalls(values[ld_rows, :], ld_genotypes[start:end])
            if layout == 1:
                _pack_calls(values, genotypes[:, start // 4:(end + 3) // 4])
            del values
        p_error = np.full(m, error_prob, dtype=float)
        if error_probs is not None:
            index = IdIndex(error_probs.sid).lookup(sids)
            valid = index >= 0
            p_error[valid] = error_probs.error_ests[index[valid]]
        markers = np.flatnonzero((freqs > min_maf) & (freqs < 1 - min_maf)
                                 & (missing <= max_missing) & (p_error <= max_error))
        if not len(markers):
            raise ValueError('No SNPs remain after MAF, missingness and genotyping error filtering')
        if mapfile is None and source.bed:
            genetic_map = _ibd_map(source, chrom, positions, None)[markers]
        else:
            genetic_map = _ibd_map(source, chrom, positions[markers], mapfile)
        keep = np.isfinite(genetic_map)
        markers, genetic_map = markers[keep], genetic_map[keep]
        if not len(markers) or np.any(np.diff(genetic_map) < 0):
            raise ValueError('SNPs must have increasing genetic map positions')
        print('Computing LD weights')
        # Keep the original variant-major LD cache on disk and pass the retained
        # source column indices to the compiled kernel.
        ld = compute_ld_scores_packed(ld_genotypes, genetic_map, 1., markers)
        weights = 1 / np.where(np.isfinite(ld) & (ld > 0), ld, 1)
        snps, pos = sids[markers], positions[markers]
        transitions, emissions = _chromosome_hmm(freqs[markers], genetic_map, weights, p_error[markers])
        _, pair_batch = plan_batches(len(pairs), len(markers) + 1, batch_size=10000,
                                     fixed_bytes=threads * (3 * len(markers) + 128))
        print(f'Tracing IBD in blocks of up to {pair_batch} sibling pairs')
        _write_ibd_batches(pairs, row_pairs, genotypes, layout, markers, transitions, emissions,
                           snps, pos, genetic_map, min_length, chrom, outprefix, ibdmatrix, pair_batch)
        if ld_out:
            filename = outfile_name(outprefix, '.l2.ldscore.gz', chrom)
            with gzip.open(filename, 'wt', compresslevel=3) as output:
                output.write('CHR SNP BP L2\n')
                output.writelines(f'{chrom} {sid} {bp} {score}\n' for sid, bp, score in zip(snps, pos, ld))
