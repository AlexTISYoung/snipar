from pysnptools.snpreader import Bed
from sibreg.sibreg import *
import argparse, gzip
from numba import njit, prange
from sibreg.bin.preprocess_data import create_pedigree

def get_sibpairs_from_ped(ped):
    parent_missing = np.array([ped[i,2]=='0' or ped[i,3]=='0' for i in range(ped.shape[0])])
    #print('Removing '+str(np.sum(parent_missing))+' rows from pedigree due to missing parent(s)')
    ped = ped[np.logical_not(parent_missing),:]
    # Find unique parent-pairs
    parent_pairs = np.array([ped[i,2]+ped[i,3] for i in range(ped.shape[0])])
    unique_pairs, sib_counts = np.unique(parent_pairs, return_counts=True)
    # Parent pairs with more than one offspring
    sib_parent_pairs = unique_pairs[sib_counts>1]
    fam_sizes = np.array([x*(x-1)/2 for x in sib_counts[sib_counts>1]])
    npairs = int(np.sum(fam_sizes))
    if npairs==0:
        raise(ValueError('No sibling pairs found'))
    # Find all sibling pairs
    sibpairs = np.zeros((npairs,2),dtype=ped.dtype)
    paircount = 0
    for ppair in sib_parent_pairs:
        sib_indices = np.where(parent_pairs==ppair)[0]
        for i in range(0,sib_indices.shape[0]-1):
            for j in range(i+1,sib_indices.shape[0]):
                sibpairs[paircount,:] = np.array([ped[sib_indices[i],1],ped[sib_indices[j],1]])
                paircount += 1
    return sibpairs

def get_sibpairs_from_king(kinfile):
    kin_header = np.array(open(kinfile,'r').readline().split('\t'))
    inf_type_index = np.where(np.array([x[0:7]=='InfType' for x in kin_header]))[0][0]
    id1_index = np.where(np.array(kin_header)=='ID1')[0][0]
    id2_index = np.where(np.array(kin_header)=='ID2')[0][0]
    sibpairs = np.loadtxt(kinfile,dtype=str,skiprows=1,usecols=(id1_index,id2_index,inf_type_index))
    sibpairs = sibpairs[sibpairs[:,2]=='FS',0:2]
    return sibpairs

def estimate_genotyping_error_rate(bedfile,ped):
    # Read bed
    print('Reading genotypes from ' + bedfile)
    bed = Bed(bedfile, count_A1=True)
    ids = bed.iid
    id_dict = make_id_dict(ids, 1)
    ## Find parent-offspring pairs
    # genotyped individuals
    genotyped = np.array([x in id_dict for x in ped[:, 1]])
    ped = ped[genotyped, :]
    # with genotyped father
    father_genotyped = np.array([x in id_dict for x in ped[:, 2]])
    # with genotyped mother
    mother_genotyped = np.array([x in id_dict for x in ped[:, 3]])
    # either
    opg = np.logical_or(father_genotyped, mother_genotyped)
    opg_ped = ped[opg, :]
    # number of pairs
    npair = np.sum(father_genotyped) + np.sum(mother_genotyped)
    if npair == 0:
        raise(ValueError('No parent-offspring pairs for genotype error probability estimation'))
    print(str(npair)+' parent-offspring pairs found')
    if npair*bed.sid.shape[0] < 10**5:
        print('Warning: limited information for estimation of genotyping error probability.')
    ## Read genotypes
    all_ids = np.unique(np.hstack((opg_ped[:, 1],
                                   ped[father_genotyped, 2],
                                   ped[mother_genotyped, 3])))
    all_ids_indices = np.sort(np.array([id_dict[x] for x in all_ids]))
    gts = bed[all_ids_indices, :].read().val
    id_dict = make_id_dict(ids[all_ids_indices, :], 1)
    gts = ma.array(gts, mask=np.isnan(gts))
    ## Get indices
    pair_indices = np.zeros((npair,2),dtype=int)
    pair_count = 0
    for i in range(opg_ped.shape[0]):
        o_index = id_dict[opg_ped[i,1]]
        if opg_ped[i, 2] in id_dict:
            pair_indices[pair_count,:] = np.array([o_index,id_dict[opg_ped[i,2]]])
            pair_count += 1
        if opg_ped[i, 3] in id_dict:
            pair_indices[pair_count,:] = np.array([o_index,id_dict[opg_ped[i,3]]])
            pair_count += 1
    ## Count Mendelian errors
    ME = count_ME(gts,pair_indices)
    # Compute allele frequencies
    freqs = ma.mean(gts, axis=0) / 2.0
    # Estimate error probability
    sum_het = np.sum(freqs * (1 - freqs))
    p_error = ME / (npair * sum_het)
    return p_error

@njit(parallel=True)
def count_ME(gts,pair_indices):
    ME = 0
    # Count Mendelian errors
    for i in prange(pair_indices.shape[0]):
        ME += np.sum(np.abs(gts[pair_indices[i,0],:]-gts[pair_indices[i,1],:])>1)
    return ME

def read_sibs_from_bed(bedfile,sibpairs):
    bed = Bed(bedfile, count_A1=True)
    ids = bed.iid
    id_dict = make_id_dict(ids, 1)
    # Find sibpairs in bed
    in_bed = np.vstack((np.array([x in id_dict for x in sibpairs[:,0]]),
                        np.array([x in id_dict for x in sibpairs[:, 1]]))).T
    both_in_bed = np.sum(in_bed,axis=1)==2
    # Remove pairs without both in bedfile
    if np.sum(both_in_bed)<sibpairs.shape[0]:
        print(str(sibpairs.shape[0]-np.sum(both_in_bed))+' sibpairs do not both have genotypes')
        sibpairs = sibpairs[both_in_bed,:]
    # Find indices of sibpairs
    sibindices = np.sort(np.array([id_dict[x] for x in sibpairs.flatten()]))
    gts = np.zeros((sibindices.shape[0],bed.sid.shape[0]),dtype=np.float32)
    gts[:] = bed[sibindices,:].read().val
    return gtarray(garray = gts, ids = ids[sibindices, 1], sid = bed.sid, pos = np.array(bed.pos[:,2],dtype=int))

# Read header of mapfile
def get_map_positions(mapfile,gts,min_map_prop = 0.5):
    map_file = open(mapfile,'r')
    map_header = map_file.readline()
    map_header = np.array(map_header.split(' '))
    map_header[len(map_header)-1] = map_header[len(map_header)-1].split('\n')[0]
    map_file.close()
    if 'pposition' in map_header and 'gposition' in map_header:
        bp_pos = np.loadtxt(mapfile,usecols = np.where(map_header=='pposition')[0][0], dtype=int, skiprows =1)
        pos_dict = make_id_dict(bp_pos)
        cm_pos = np.loadtxt(mapfile,usecols = np.where(map_header=='gposition')[0][0], dtype=float, skiprows =1)
        # Check for NAs
        if np.sum(np.isnan(cm_pos)) > 0:
            raise (ValueError('Map cannot have NAs'))
        if np.min(cm_pos) < 0:
            raise (ValueError('Map file cannot have negative values'))
        if np.var(cm_pos) == 0:
            raise (ValueError('Map file has no variation'))
        # Check ordering
        ordered_map = np.sort(cm_pos)
        if np.array_equal(cm_pos, ordered_map):
            pass
        else:
            raise (ValueError('Map not monotonic. Please make sure input is ordered correctly'))
        # Check scale
        if np.max(cm_pos) > 5000:
            raise (ValueError('Maximum value of map too large'))
        # Find positions of SNPs in map file
        map = np.zeros((gts.shape[1]),dtype=float)
        map[:] = np.nan
        in_map = np.array([x in pos_dict for x in gts.pos])
        # Check if we have at least 50% of SNPs in map
        prop_in_map = np.mean(in_map)
        if prop_in_map < min_map_prop:
            raise(ValueError('Only '+str(round(100*prop_in_map))+'% of SNPs have genetic positions in '+mapfile+'. Need at least '+str(round(100*min_map_prop))+'%'))
        print('Found genetic map positions for '+str(round(100*prop_in_map))+'% of SNPs in '+mapfile)
        # Fill in map values
        map[in_map] = cm_pos[[pos_dict[x] for x in gts.pos[in_map]]]
        # Linearly interpolate map
        if prop_in_map < 1:
            print('Linearly interpolating genetic map for SNPs not in input map')
            map = np.interp(gts.pos, gts.pos[in_map], map[in_map])
        return map
    else:
        raise(ValueError('Map file must contain columns pposition and gposition'))

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
    probs = np.zeros((3,3,3))
    probs[0,...] = p_ibd_0(f)
    probs[1,...] = p_ibd_1(f)
    probs[2,...] = p_ibd_2(f)
    # Compute probabilities conditional on IBD
    for i in range(3):
        P_out[i] = err_vectors[g1_obs, ...].T @ probs[i,...] @ err_vectors[g2_obs, ...]
    # Return
    return np.log(P_out)

@njit
def make_dynamic(g1,g2,freqs,map,weights,p):
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
    state_matrix = np.zeros((3,g1.shape[0]),dtype=np.float64)
    pointers = np.zeros((3,g1.shape[0]),dtype=np.int8)
    # Check for nans
    not_nan = np.logical_not(np.logical_or(np.isnan(g1),np.isnan(g2)))
    # Initialise
    state_matrix[:,0] = np.log(np.array([0.25,0.5,0.25],dtype=np.float64))
    if not_nan[0]:
        state_matrix[:,0] += weights[0]*p_obs_given_IBD(np.int8(g1[0]),np.int8(g2[0]),freqs[0],p)
    # Compute
    for l in range(1,g1.shape[0]):
        if not_nan[l]:
            probs = weights[l]*p_obs_given_IBD(np.int8(g1[l]),np.int8(g2[l]),freqs[l],p)
        else:
            probs = np.zeros((3))
        tmatrix = transition_matrix(map[l]-map[l-1])
        for i in range(3):
            tprobs = tmatrix[:,i]+state_matrix[:,l-1]
            state_matrix[i,l] = np.max(tprobs)+probs[i]
            pointers[i,l] = np.argmax(tprobs)
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
    path[path.shape[0]-1] = np.argmax(state_matrix[:,state_matrix.shape[1]-1])
    for i in range(1,path.shape[0]):
        path[path.shape[0]-(i+1)] = pointers[path[path.shape[0]-i],path.shape[0]-i]
    return path

@njit(parallel=True)
def infer_ibd(sibpairs,gts,freqs,map,weights,p):
    ibd = np.zeros((sibpairs.shape[0],gts.shape[1]),dtype=np.int8)
    for i in prange(sibpairs.shape[0]):
        sibpair = sibpairs[i,:]
        state_matrix, pointers = make_dynamic(gts[sibpair[0],:],gts[sibpair[1],:],freqs,map,weights,p)
        ibd[i,...] = viterbi(state_matrix,pointers)
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

def write_segs_from_matrix(ibd,bimfile,outfile):
    # Get IBD
    sibpairs = ibd[1:ibd.shape[0],0:2]
    snps = ibd[0,2:ibd.shape[1]]
    ibd = np.array(ibd[1:ibd.shape[0],2:ibd.shape[1]],dtype=np.int8)
    # Get map and chr from bim
    bim = np.loadtxt(bimfile, dtype=str)
    bim_dict = make_id_dict(bim,1)
    bim_indices = np.array([bim_dict[x] for x in snps])
    map = np.array(bim[bim_indices,2],dtype=float)
    chr = np.unique(bim[bim_indices,0])
    if chr.shape[0]>1:
        raise(ValueError('Must be one chromosome only'))
    else:
        chr = chr[0]
    # Get segments
    allsegs = []
    for i in range(sibpairs.shape[0]):
        allsegs.append(find_segments(ibd[i,:], map))
    # Write segments
    write_segs(sibpairs,allsegs,chr,snps,outfile)
    return allsegs

parser = argparse.ArgumentParser()
parser.add_argument('--bed',
                    type=str,help='Address of the unphased genotypes in .bed format.')
parser.add_argument('--king',
                    type=str,
                    default = None,
                    help='Address of the king file')
parser.add_argument('--agesex',
                    type=str,
                    default=None,
                    help='Address of file with age and sex information')
parser.add_argument('--pedigree',
                    type=str,
                    default=None,
                    help='Address of pedigree file')
parser.add_argument('--ldscores',type=str,help='Path to directory with LD scores',default=None)
parser.add_argument('--map',type=str,default=None)
parser.add_argument('--outprefix',
                    type=str,
                    default = 'ibd',
                    help="Writes the result of IBD inference to outprefix.ibd.segments.gz")
parser.add_argument('--chrom',type=str,default=None,help='Chromosome (only input one chromosome at a time)')
parser.add_argument('--p_error',type=float,help='Probability of genotyping error',default=None)
parser.add_argument('--min_length',type=float,help='Smooth segments with length less than min_length (cM)',
                    default=0.01)
parser.add_argument('--min_maf',type=float,help='Minimum minor allele frequency',default=0.01)
parser.add_argument('--ibdmatrix',action='store_true',default=False,help='Whether to output a matrix of SNP IBD states')
args = parser.parse_args()

min_length = args.min_length
bedfile = args.bed+'.bed'
kinfile = args.king
ldfile = args.ldscores
min_maf = args.min_maf
mapfile = args.map
outprefix = args.outprefix

#### Find sibling pairs ####
if args.pedigree is not None:
    print('Reading pedigree from '+str(args.pedigree))
    ped = np.loadtxt(args.pedigree,dtype=str)
    if ped.shape[1] < 4:
        raise(ValueError('Not enough columns in pedigree file'))
    elif ped.shape[1] > 4:
        print('Warning: pedigree file has more than 4 columns. The first four columns only will be used')
    # Remove rows with missing parents
    sibpairs = get_sibpairs_from_ped(ped)
elif kinfile is not None:
    print('Reading relationships from '+str(kinfile))
    sibpairs = get_sibpairs_from_king(kinfile)
else:
    raise(ValueError('Must provide either KING kinship file or pedigree'))

if sibpairs.shape[0]==0:
    raise(ValueError('No sibling pairs found'))
print('Found '+str(sibpairs.shape[0])+' full sibling pairs')

#### Get genotyping error probability ####
if args.p_error is None:
    print('No genotyping error probability provided. Will attempt to estimate from parent-offspring pairs.')
    if args.pedigree is None:
        if args.agesex is not None:
            print('Constructing pedigree from '+str(kinfile)+' and age and sex information from '+str(args.agesex))
            ped = np.array(create_pedigree(kinfile,args.agesex),dtype=str)
        else:
            raise(ValueError('Must provide age and sex information (--agesex) in addition to KING kinship file, if estimating genotyping error probability'))
    else:
        ped = np.loadtxt(args.pedigree, dtype=str)
    p = estimate_genotyping_error_rate(bedfile, ped)
    print('Estimated genotyping error probability: '+str(round(p,6)))
    if p > 0.05:
        print('Warning: high genotyping error rate detected. Check pedigree and/or genotype data.')

## Read bed
print('Reading genotypes from '+bedfile)
# Determine chromosome
bimfile = bedfile.split('.bed')[0]+'.bim'
chr = np.loadtxt(bimfile,usecols=0,dtype=str)
chr = np.unique(chr)
if chr.shape[0]>1:
    raise(ValueError('More than 1 chromosome in input bedfile'))
else:
    chr = chr[0]
print('Inferring IBD for chromosome '+str(chr))
# Read sibling genotypes from bed file
gts = read_sibs_from_bed(bedfile,sibpairs)

# Calculate allele frequencies
print('Calculating allele frequencies')
gts.compute_freqs()

# Check which sibling pairs have genotypes
sibpair_indices = np.zeros((sibpairs.shape),dtype=bool)
sibpair_indices[:,0] = np.array([x in gts.id_dict for x in sibpairs[:,0]])
sibpair_indices[:,1] = np.array([x in gts.id_dict for x in sibpairs[:,1]])
sibpairs = sibpairs[np.sum(sibpair_indices,axis=1)==2,:]
if sibpairs.shape[0]==0:
    raise(ValueError('No genotyped sibling pairs found'))
print(str(np.sum(sibpairs.shape[0]))+' sibpairs have genotypes in '+bedfile)
# Find indices of sibpairs
sibpair_indices = np.zeros((sibpairs.shape),dtype=int)
sibpair_indices[:,0] = np.array([gts.id_dict[x] for x in sibpairs[:,0]])
sibpair_indices[:,1] = np.array([gts.id_dict[x] for x in sibpairs[:,1]])

# LD scores
print('Reading LD scores from '+ldfile)
ld = np.loadtxt(ldfile,usecols=3,skiprows=1)
ldsnps = np.loadtxt(ldfile,usecols=1,skiprows=1,dtype=str)
ld_dict = make_id_dict(ldsnps)
in_ld_dict = np.array([x in ld_dict for x in gts.sid])

# Filtering on MAF and LD score
print('Before filtering on MAF and having an LD score, there were '+str(gts.shape[1])+' SNPs')
freq_pass = np.logical_and(np.logical_and(gts.freqs > min_maf, gts.freqs < 1-min_maf),in_ld_dict)
gts.filter(freq_pass)
print('After filtering on MAF and having an LD score, there are '+str(gts.shape[1])+' SNPs')

# Read map file
if mapfile is None:
    print('Separate genetic map not provided, so attempting to read map from '+bimfile)
    map = np.loadtxt(bimfile, usecols=2)
    # Check for NAs
    if np.sum(np.isnan(map))>0:
        raise(ValueError('Map cannot have NAs'))
    if np.min(map)<0:
        raise(ValueError('Map file cannot have negative values'))
    if np.var(map)==0:
        raise(ValueError('Map file has no variation'))
    # Check ordering
    ordered_map = np.sort(map)
    if np.array_equal(map,ordered_map):
        pass
    else:
        raise(ValueError('Map not monotonic. Please make sure input is ordered correctly'))
    # Check scale
    if np.max(map)>5000:
        raise(ValueError('Maximum value of map too large'))
    gts.map = map[freq_pass]
else:
    print('Reading map from '+str(mapfile))
    gts.map = get_map_positions(mapfile,gts)
print('Read map')

# Weights
gts.weights = np.power(ld[np.array([ld_dict[x] for x in gts.sid])],-1)

# IBD
print('Inferring IBD')
ibd = infer_ibd(sibpair_indices,np.array(gts.gts),gts.freqs,gts.map,gts.weights,p)
ibd, allsegs = smooth_ibd(ibd,gts.map,gts.sid,gts.pos,min_length)

## Write output
# Write segments
segs_outfile = outprefix+'.ibd.segments.gz'
print('Writing segments to '+segs_outfile)
write_segs(sibpairs,allsegs,chr,segs_outfile)
# Write matrix
if args.ibdmatrix:
    outfile = outprefix+'.ibd.gz'
    print('Writing matrix output to '+str(outfile))
    ibd = np.row_stack((np.column_stack((np.array(['sib1','sib2']).reshape((1,2)),gts.sid.reshape(1,gts.shape[1]))),
                        np.column_stack((sibpairs,ibd))))
    np.savetxt(outfile,ibd,fmt='%s')