import gzip
from numba import njit, prange
from snipar.map import decode_map_from_pos, get_map_positions
from snipar.ld import compute_ld_scores
import numpy as np
from snipar.read.bed import read_sibs_from_bed
from snipar.read.bgen import read_sibs_from_bgen
from snipar.utilities import make_id_dict
from snipar.utilities import outfile_name
from snipar.utilities import open_bgen

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


@njit(cache=True, parallel=True)
def infer_ibd(sibpairs, gts, freqs, map, weights, error_probs):
    n = gts.shape[1]
    ibd = np.empty((len(sibpairs), n), dtype=np.int8)
    if not n:
        return ibd
    transitions, emissions = _chromosome_hmm(freqs, map, weights, error_probs)
    for i in prange(len(sibpairs)):
        g1, g2 = gts[sibpairs[i, 0]], gts[sibpairs[i, 1]]
        previous = np.log(np.array([.25, .5, .25]))
        current = np.empty(3)
        pointers = np.empty((n, 3), dtype=np.int8)
        if not (np.isnan(g1[0]) or np.isnan(g2[0])):
            previous += emissions[0, int(g1[0]), int(g2[0])]
        for marker in range(1, n):
            observed = not (np.isnan(g1[marker]) or np.isnan(g2[marker]))
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
                    current[dest] += emissions[marker, int(g1[marker]), int(g2[marker]), dest]
            previous, current = current, previous
        state = np.argmax(previous)
        ibd[i, n - 1] = state
        for marker in range(n - 1, 0, -1):
            state = pointers[marker, state]
            ibd[i, marker - 1] = state
    return ibd

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

def infer_ibd_chr(sibpairs, error_prob, error_probs, outprefix, bedfile=None, bgenfile=None, chrom=None, min_length=0.01, mapfile=None, ibdmatrix=False, ld_out=False, min_maf=0.01, max_missing=5, max_error=0.01):
    if bedfile is None and bgenfile is None:
        raise(ValueError('Must provide either bed file or bgenfile'))
    if bedfile is not None and bgenfile is not None:
        raise(ValueError('Provide either bed file or bgen file. Not both.'))
    if bedfile is not None:
        ## Read bed
        print('Reading genotypes from ' + bedfile)
        bimfile = bedfile.removesuffix(".bed") + '.bim'
        # Determine chromosome
        if chrom is None:
            chrom = np.loadtxt(bimfile, usecols=0, dtype=str)
            chrom = np.unique(chrom)
            if chrom.shape[0] > 1:
                raise (ValueError('More than 1 chromosome in input bedfile'))
            else:
                chrom = chrom[0]
        # Read sibling genotypes from bed file
        gts = read_sibs_from_bed(bedfile, sibpairs)
        print('Inferring IBD for chromosome ' + str(chrom))
    elif bgenfile is not None:
        ## Read bed
        print('Reading genotypes from ' + bgenfile)
        # Determine chromosome
        if chrom is None:
            bgen = open_bgen(bgenfile,verbose=False)
            chrom = bgen.chromosomes
            chrom = np.unique(chrom)
            if chrom.shape[0] > 1:
                raise (ValueError('More than 1 chromosome in input bgenfile'))
            else:
                chrom = chrom[0]
                if chrom=='':
                    chrom = 0
        # Read sibling genotypes from bed file
        gts = read_sibs_from_bgen(bgenfile, sibpairs)
        print('Inferring IBD for chromosome ' + str(chrom))
    # Calculate allele frequencies
    print('Calculating allele frequencies')
    gts.compute_freqs()
    # Check which sibling pairs have genotypes
    sibpair_indices = np.zeros((sibpairs.shape), dtype=bool)
    sibpair_indices[:, 0] = np.array([x in gts.id_dict for x in sibpairs[:, 0]])
    sibpair_indices[:, 1] = np.array([x in gts.id_dict for x in sibpairs[:, 1]])
    sibpairs = sibpairs[np.sum(sibpair_indices, axis=1) == 2, :]
    if sibpairs.shape[0] == 0:
        raise (ValueError('No genotyped sibling pairs found'))
    print(str(np.sum(sibpairs.shape[0])) + ' sibpairs have genotypes')
    # Find indices of sibpairs
    sibpair_indices = np.zeros((sibpairs.shape), dtype=int)
    sibpair_indices[:, 0] = np.array([gts.id_dict[x] for x in sibpairs[:, 0]])
    sibpair_indices[:, 1] = np.array([gts.id_dict[x] for x in sibpairs[:, 1]])
    # Filtering on MAF, LD score, and genotyping error
    # Find error probabilities
    p_error = np.zeros((gts.sid.shape[0]))
    p_error[:] = error_prob
    if error_probs is not None:
        in_error_probs = np.array([x in error_probs.sid_dict for x in gts.sid])
        error_index = np.array([error_probs.sid_dict[x] for x in gts.sid[in_error_probs]])
        p_error[in_error_probs] = error_probs.error_ests[error_index]
    gts.error_probs = p_error
    # Filter
    print('Before filtering on MAF, missingness, and genotyping error, there were ' + str(gts.shape[1]) + ' SNPs')
    gts.filter_maf(min_maf)
    gts.filter_missingness(max_missing)
    gts.filter(gts.error_probs <= max_error)
    if not gts.shape[1]:
        raise ValueError("No SNPs remain after MAF, missingness and genotyping error filtering")
    print('After filtering, there are ' + str(gts.shape[1]) + ' SNPs')
    # Read map file
    if mapfile is None and bedfile is not None:
        print('Separate genetic map not provided, so attempting to read map from ' + bimfile)
        map = np.loadtxt(bimfile, usecols=2, ndmin=1)
        map_snp_dict = make_id_dict(np.loadtxt(bimfile, usecols=1, dtype=str, ndmin=1))
        # Check for NAs
        if np.all(map == 0) or (map.size > 1 and np.var(map) == 0):
            print('Map information not found in bim file.')
            print('Using default map (decode sex averaged map on GRCh38 coordinates)')
            gts.map = decode_map_from_pos(chrom, gts.pos)
            pc_mapped = str(round(100*(1-np.mean(np.isnan(gts.map))),2))
            print('Found map positions for '+str(pc_mapped)+'% of SNPs')
            gts.filter(~np.isnan(gts.map))
        else:
            if not np.isfinite(map).all():
                raise (ValueError('Map cannot have NAs'))
            if np.min(map) < 0:
                raise (ValueError('Map file cannot have negative values'))
            # Check ordering
            ordered_map = np.sort(map)
            if np.array_equal(map, ordered_map):
                pass
            else:
                raise (ValueError('Map not monotonic. Please make sure input is ordered correctly'))
            # Check scale
            if np.max(map) > 5000:
                raise (ValueError('Maximum value of map too large'))
            gts.filter(np.array([x in map_snp_dict for x in gts.sid]))
            gts.map = map[[map_snp_dict[x] for x in gts.sid]]
    elif mapfile is None and bgenfile is not None:
        print('Map file not provided.')
        print('Using default map (decode sex averaged map on GRCh38 coordinates)')
        gts.map = decode_map_from_pos(chrom, gts.pos)
        pc_mapped = 100*(1-np.mean(np.isnan(gts.map)))
        if pc_mapped < 50:
            print('Warning: map positions not found for the majority of SNPs. Consider providing a genetic map using --map')
        print('Found map positions for '+str(round(pc_mapped,2))+'% of SNPs')
        gts.filter(~np.isnan(gts.map))
    else:
        print('Reading map from ' + str(mapfile))
        gts.map = get_map_positions(mapfile, gts)
    if not gts.shape[1]:
        raise ValueError("No SNPs have genetic map positions")
    if np.any(np.diff(gts.map) < 0):
        raise ValueError("SNPs must be in increasing genetic map order")
    print('Read map')
    # Weights
    print('Computing LD weights')
    ld_n = np.min([1000, gts.shape[0]])
    ld_sample = np.random.randint(0, gts.shape[0], size=ld_n)
    ld = compute_ld_scores(np.array(gts.gts[ld_sample,:], dtype=np.float64), gts.map, max_dist=1)
    gts.weights = 1 / np.where(np.isfinite(ld) & (ld > 0), ld, 1)
    # IBD
    print('Inferring IBD')
    ibd = infer_ibd(sibpair_indices, np.array(gts.gts,dtype=np.float64), gts.freqs, gts.map, gts.weights, gts.error_probs)
    ibd, allsegs = smooth_ibd(ibd, gts.map, gts.sid, gts.pos, min_length)
    ## Write output
    # Write segments
    segs_outfile = outfile_name(outprefix,'.ibd.segments.gz', chrom)
    print('Writing segments to ' + segs_outfile)
    write_segs(sibpairs, allsegs, chrom, segs_outfile)
    # Write matrix
    if ibdmatrix:
        outfile = outfile_name(outprefix,'.ibdmatrix.gz', chrom)
        print('Writing matrix output to ' + str(outfile))
        ibd = np.vstack(
            (np.column_stack((np.array(['sib1', 'sib2']).reshape((1, 2)), gts.sid.reshape(1, gts.shape[1]))),
             np.column_stack((sibpairs, ibd))))
        np.savetxt(outfile, ibd, fmt='%s')
    if ld_out:
        ld_outfile = outfile_name(outprefix,'.l2.ldscore.gz', chrom)
        print('Writing LD-scores to '+ld_outfile)
        ld_out = np.vstack((np.array(['CHR', 'SNP', 'BP', 'L2']).reshape((1,4)),np.vstack((np.array([chrom for x in gts.sid]), gts.sid, gts.pos, ld)).T))
        np.savetxt(ld_outfile, ld_out, fmt='%s')