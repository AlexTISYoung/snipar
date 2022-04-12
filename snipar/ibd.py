import gzip
from numba import njit, prange
from snipar.map import *
from snipar.ld import compute_ld_scores
import numpy as np
from snipar.read.bed import read_sibs_from_bed
from snipar.read.bgen import read_sibs_from_bgen
from snipar.utilities import make_id_dict
from snipar.utilities import outfile_name
from bgen_reader import open_bgen

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
    fmatrix = np.array([[0.0,1.0,0.0],[1.0,1.0,2.0],[0.0,2.0,3.0]])
    minus_fmatrix = np.array([[3.0,2.0,0.0],[2.0,1.0,1.0],[0.0,1.0,0.0]])
    P = np.exp(fmatrix*np.log(f)+minus_fmatrix*np.log(1-f))
    P[0,2] = 0
    P[2,0] = 0
    return P

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
    path[path.shape[0]-1] = np.argmax(state_matrix[:, state_matrix.shape[1]-1])
    for i in range(1,path.shape[0]):
        path[path.shape[0]-(i+1)] = pointers[path[path.shape[0]-i], path.shape[0]-i]
    return path

@njit(parallel=True)
def infer_ibd(sibpairs, gts, freqs, map, weights, error_probs):
    ibd = np.zeros((sibpairs.shape[0], gts.shape[1]), dtype=np.int8)
    for i in prange(sibpairs.shape[0]):
        sibpair = sibpairs[i, :]
        state_matrix, pointers = make_dynamic(gts[sibpair[0], :], gts[sibpair[1], :], freqs, map, weights, error_probs)
        ibd[i, ...] = viterbi(state_matrix, pointers)
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

def find_segments(path,map,snps,pos):
    segments = []
    ibd_start = path[0]
    ibd_start_index = 0
    for i in range(1,path.shape[0]):
        if not path[i] == ibd_start:
            segments.append(segment(ibd_start_index,i,
                                    pos[ibd_start_index],pos[i-1],
                                    snps[ibd_start_index],snps[i-1],
                                    map[i-1]-map[ibd_start_index],ibd_start))
            ibd_start_index = i
            ibd_start = path[i]
    segments.append(segment(ibd_start_index,i,
                            pos[ibd_start_index],pos[i],
                            snps[ibd_start_index],snps[i],
                            map[i]-map[ibd_start_index],ibd_start))
    return segments

def smooth_segments(path,map,snps,pos,min_length):
    # Smooth path
    segments = find_segments(path,map,snps,pos)
    if len(segments)>1:
        for i in range(len(segments)):
            if segments[i].length < min_length:
                if i==0:
                    if segments[i].start == segments[i].end:
                        path[segments[i].start] = segments[i+1].state
                    else:
                        path[segments[i].start:segments[i].end] = segments[i+1].state
                elif i<(len(segments)-1):
                    if segments[i-1].state == segments[i+1].state:
                        if segments[i].start == segments[i].end:
                            path[segments[i].start] = segments[i+1].state
                        else:
                            path[segments[i].start:segments[i].end] = segments[i+1].state
                else:
                    if segments[i].start == segments[i].end:
                        path[segments[i].start] = segments[i-1].state
                    else:
                        path[segments[i].start:segments[i].end] = segments[i - 1].state
        segments = find_segments(path,map,snps,pos)
    return path, segments

def smooth_ibd(ibd,map,snps,pos,min_length):
    allsegs = []
    for i in range(ibd.shape[0]):
        ibd_path, segments = smooth_segments(ibd[i,:],map,snps,pos,min_length)
        ibd[i,:] = ibd_path
        allsegs.append(segments)
    return ibd, allsegs

def write_segs(sibpairs,allsegs,chr,outfile):
    seg_out = gzip.open(outfile,'wb')
    # Header
    seg_out.write('ID1\tID2\tIBDType\tChr\tstart_coordinate\tstop_coordinate\tstartSNP\tstopSNP\tlength\n'.encode())
    # Write segment lines
    for i in range(0,sibpairs.shape[0]):
        nseg = len(allsegs[i])
        for j in range(nseg):
            if i == (sibpairs.shape[0]-1) and j == (nseg-1):
                seg_out.write(allsegs[i][j].to_text(sibpairs[i, 0], sibpairs[i, 1], chr, end=True).encode())
            else:
                seg_out.write(allsegs[i][j].to_text(sibpairs[i, 0], sibpairs[i, 1], chr, end=False).encode())
    seg_out.close()

def write_segs_from_matrix(ibd,sibpairs,snps,pos,map,chrom,outfile):
    # Get segments
    allsegs = []
    for i in range(sibpairs.shape[0]):
        allsegs.append(find_segments(ibd[i,:],map,snps,pos))
    # Write segments
    write_segs(sibpairs,allsegs,chrom,outfile)
    return allsegs

def infer_ibd_chr(sibpairs, error_prob, error_probs, outprefix, bedfile=None, bgenfile=None, chrom=None, min_length=0.01, mapfile=None, ibdmatrix=False, ld_out=False, min_maf=0.01, max_missing=5, max_error=0.01):
    if bedfile is None and bgenfile is None:
        raise(ValueError('Must provide either bed file or bgenfile'))
    if bedfile is not None and bgenfile is not None:
        raise(ValueError('Provide either bed file or bgen file. Not both.'))
    if bedfile is not None:
        ## Read bed
        print('Reading genotypes from ' + bedfile)
        bimfile = bedfile.split('.bed')[0] + '.bim'
        # Determine chromosome
        if chrom is None:
            chrom = np.loadtxt(bimfile, usecols=0, dtype=str)
            chrom = np.unique(chrom)
            if chrom.shape[0] > 1:
                raise (ValueError('More than 1 chromosome in input bedfile'))
            else:
                chrom = chrom[0]
        print('Inferring IBD for chromosome ' + str(chrom))
        # Read sibling genotypes from bed file
        gts = read_sibs_from_bed(bedfile, sibpairs)
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
        print('Inferring IBD for chromosome ' + str(chrom))
        # Read sibling genotypes from bed file
        gts = read_sibs_from_bgen(bgenfile, sibpairs)
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
    gts.filter(gts.error_probs < max_error)
    print('After filtering, there are ' + str(gts.shape[1]) + ' SNPs')
    # Read map file
    if mapfile is None and bedfile is not None:
        print('Separate genetic map not provided, so attempting to read map from ' + bimfile)
        map = np.loadtxt(bimfile, usecols=2)
        map_snp_dict = make_id_dict(np.loadtxt(bimfile, usecols=1, dtype=str))
        # Check for NAs
        if np.var(map) == 0:
            print('Map information not found in bim file.')
            print('Using default map (decode sex averaged map on Hg19 coordinates)')
            gts.map = decode_map_from_pos(chrom, gts.pos)
            pc_mapped = str(round(100*(1-np.mean(np.isnan(gts.map))),2))
            print('Found map positions for '+str(pc_mapped)+'% of SNPs')
            gts.filter(~np.isnan(gts.map))
        else:
            if np.sum(np.isnan(map)) > 0:
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
        print('Using default map (decode sex averaged map on Hg19 coordinates)')
        gts.map = decode_map_from_pos(chrom, gts.pos)
        pc_mapped = 100*(1-np.mean(np.isnan(gts.map)))
        if pc_mapped < 50:
            print('Warning: map positions not found for the majority of SNPs. Consider providing a genetic map using --map')
        print('Found map positions for '+str(round(pc_mapped,2))+'% of SNPs')
        gts.filter(~np.isnan(gts.map))
    else:
        print('Reading map from ' + str(mapfile))
        gts.map = get_map_positions(mapfile, gts)
    print('Read map')
    # Weights
    print('Computing LD weights')
    ld = compute_ld_scores(np.array(gts.gts, dtype=np.float_), gts.map, max_dist=1)
    gts.weights = np.power(ld, -1)
    # IBD
    print('Inferring IBD')
    ibd = infer_ibd(sibpair_indices, np.array(gts.gts,dtype=np.float_), gts.freqs, gts.map, gts.weights, gts.error_probs)
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
        ibd = np.row_stack(
            (np.column_stack((np.array(['sib1', 'sib2']).reshape((1, 2)), gts.sid.reshape(1, gts.shape[1]))),
             np.column_stack((sibpairs, ibd))))
        np.savetxt(outfile, ibd, fmt='%s')
    if ld_out:
        ld_outfile = outfile_name(outprefix,'.l2.ldscore.gz', chrom)
        print('Writing LD-scores to '+ld_outfile)
        ld_out = np.vstack((np.array(['CHR', 'SNP', 'BP', 'L2']).reshape((1,4)),np.vstack((np.array([chrom for x in gts.sid]), gts.sid, gts.pos, ld)).T))
        np.savetxt(ld_outfile, ld_out, fmt='%s')