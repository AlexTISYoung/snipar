from pysnptools.snpreader import Bed
from sibreg.sibreg import *
import argparse, gzip
from numba import njit, jit, prange

@njit
def transition_matrix(cM):
    P = np.identity(3)
    r_array = np.array([[-1.0,1.0,0.0],[0.5,-1.0,0.5],[0.0,1.0,-1.0]],dtype=np.float64)
    r = 1.0 - np.power((1.0 + np.exp(-cM / 50)) / 2.0, 4)
    P += r_array*r
    return np.log(P)

@njit
def p_ibd_0(f):
    P_vec = np.array([(1.0-f)**2.0,2.0*f*(1-f),f**2.0]).reshape((3,1))
    return P_vec @ P_vec.T

@njit
def p_ibd_1(f):
    fmatrix = np.array([[0.0,1.0,0.0],[1.0,1.0,2.0],[0.0,2.0,3.0]])
    minus_fmatrix = np.array([[3.0,2.0,0.0],[2.0,1.0,1.0],[0.0,1.0,0.0]])
    P = np.exp(fmatrix*np.log(f)+minus_fmatrix*np.log(1-f))
    P[0,2] = 0
    P[2,0] = 0
    return P

@njit
def p_ibd_2(f):
    P = np.zeros((3,3),dtype=np.float64)
    fdiag = np.array([(1.0-f)**2.0,2.0*f*(1.0-f),f**2.0])
    np.fill_diagonal(P,fdiag.reshape((3,1)))
    return P

@njit
def p_obs_given_IBD(g1_obs,g2_obs,f,p):
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
    state_matrix = np.zeros((3,g1.shape[0]),dtype=np.float64)
    pointers = np.zeros((3,g1.shape[0]),dtype=np.int8)
    # Initialise
    state_matrix[:,0] = np.log(np.array([0.25,0.5,0.25],dtype=np.float64))+weights[0]*p_obs_given_IBD(g1[0],g2[0],freqs[0],p)
    # Compute
    for l in range(1,g1.shape[0]):
        probs = weights[l]*p_obs_given_IBD(g1[l],g2[l],freqs[l],p)
        tmatrix = transition_matrix(map[l]-map[l-1])
        for i in range(3):
            tprobs = tmatrix[:,i]+state_matrix[:,l-1]
            state_matrix[i,l] = np.max(tprobs)+probs[i]
            pointers[i,l] = np.argmax(tprobs)
    return state_matrix, pointers

@njit
def viterbi(state_matrix,pointers):
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
    def __init__(self,start_index,end_index,length,state):
        self.start = start_index
        self.end = end_index
        self.length = length
        self.state = state
    def to_text(self,id1,id2,chr,snps,end=False):
        start_snp = snps[self.start]
        stop_snp = snps[self.end-1]
        seg_txt = str(id1)+'\t'+str(id2)+'\t'+str(self.state)+'\t'+str(chr)+'\t'+start_snp+'\t'+stop_snp+'\t'+str(self.length)
        if end:
            return seg_txt
        else:
            return seg_txt+'\n'

def find_segments(path,map):
    segments = []
    ibd_start = path[0]
    ibd_start_index = 0
    for i in range(1,path.shape[0]):
        if not path[i] == ibd_start:
            segments.append(segment(ibd_start_index,i,map[i-1]-map[ibd_start_index],ibd_start))
            ibd_start_index = i
            ibd_start = path[i]
    segments.append(segment(ibd_start_index,i,map[i]-map[ibd_start_index],ibd_start))
    return segments

def smooth_segments(path,map,p_length):
    # Smooth path
    segments = find_segments(path,map)
    if len(segments)>1:
        for i in range(len(segments)):
            length_prob = 1-np.exp(-segments[i].length/100.0)
            if length_prob<p_length:
                #print('Segment prob: '+str(length_prob))
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
        segments = find_segments(path,map)
    return path, segments

def smooth_ibd(ibd,map,p_length):
    allsegs = []
    for i in range(ibd.shape[0]):
        ibd_path, segments = smooth_segments(ibd[i,:],map,p_length)
        ibd[i,:] = ibd_path
        allsegs.append(segments)
    return ibd, allsegs

def write_segs(sibpairs,allsegs,chr,snps,outfile):
    seg_out = gzip.open(outfile,'wb')
    # Header
    seg_out.write('ID1\tID2\tIBDType\tChr\tStartSNP\tStopSNP\tLength\n'.encode())
    # Remove zero IBD segments
    for i in range(0,sibpairs.shape[0]):
        for seg in allsegs[i]:
            if seg.state == 0:
                allsegs[i].remove(seg)
    # Write segment lines
    for i in range(0,sibpairs.shape[0]):
        nseg = len(allsegs[i])
        for j in range(nseg):
            if i == (sibpairs.shape[0]-1) and j == (nseg-1):
                seg_out.write(allsegs[i][j].to_text(sibpairs[i, 0], sibpairs[i, 1], chr, snps, end=True).encode())
            else:
                seg_out.write(allsegs[i][j].to_text(sibpairs[i, 0], sibpairs[i, 1], chr, snps, end=False).encode())
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
                    type=str,help='Address of the unphased genotypes in .bed format. If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of [from_chr, to_chr) for each chromosome(from_chr and to_chr are two optional parameters for this script). Should be supplemented with from_chr and to_chr.')
parser.add_argument('--king',
                    type=str,
                    default = None,
                    help='Address of the king file')
parser.add_argument('--ldscores',type=str,help='Path to directory with LD scores',default=None)
parser.add_argument('--map',type=str,default=None)
parser.add_argument('--outprefix',
                    type=str,
                    default = "parent_imputed",
                    help="Writes the result of imputation for chromosome i to outprefix{i}")
parser.add_argument('--min_p_seg',type=float,help='Smooth short segments with probability of length smaller than that less than min_p_seg',
                    default=1e-3)
parser.add_argument('--p_error',type=float,help='Probability of genotyping error',default=1e-4)
parser.add_argument('--min_maf',type=float,help='Minimum minor allele frequency',default=1e-4)
args = parser.parse_args()

p = args.p_error
p_length = args.min_p_seg
bedfile = args.bed+'.bed'
kinfile = args.king
ldfile = args.ldscores
min_maf = args.min_maf
mapfile = args.map
outprefix = args.outprefix

#p = 1e-4
#p_length = 1e-3
#bedfile = 'bedfiles/ukb_hap_chr22.bed'
#kinfile = 'bedfiles/hap.kin0'
#ldfile = 'bedfiles/ldscores/22.l2.ldscore.gz'
#min_maf = 0.01
#mapfile = None
#outprefix = 'ibd/chr_22'

kin_header = np.array(open(kinfile,'r').readline().split('\t'))
inf_type_index = np.where(np.array([x[0:7]=='InfType' for x in kin_header]))[0][0]
id1_index = np.where(np.array(kin_header)=='ID1')[0][0]
id2_index = np.where(np.array(kin_header)=='ID2')[0][0]
sibpairs = np.loadtxt(kinfile,dtype=str,skiprows=1,usecols=(id1_index,id2_index,inf_type_index))
sibpairs = sibpairs[sibpairs[:,2]=='FS',0:2]
print('Found '+str(sibpairs.shape[0])+' full sibling pairs in '+kinfile)

# Read bed
print('Reading genotypes from '+bedfile)
bed = Bed(bedfile,count_A1=True)
ids = bed.iid
id_dict = make_id_dict(ids,1)
sibindices = np.sort(np.array([id_dict[x] for x in sibpairs.flatten()]))
gts = bed.read().val
freqs = np.mean(gts,axis=0)/2.0
gts = gts[sibindices,:]
gts = np.array(gts,dtype=np.int8)
id_dict = make_id_dict(ids[sibindices,:],1)
# Check which sibling pairs have genotypes
sibpair_indices = np.zeros((sibpairs.shape),dtype=bool)
sibpair_indices[:,0] = np.array([x in id_dict for x in sibpairs[:,0]])
sibpair_indices[:,1] = np.array([x in id_dict for x in sibpairs[:,1]])
sibpairs = sibpairs[np.sum(sibpair_indices,axis=1)==2,:]
print(str(np.sum(sibpairs.shape[0]))+' sibpairs have genotypes in '+bedfile)
# Find indices of sibpairs
sibpair_indices = np.zeros((sibpairs.shape),dtype=int)
sibpair_indices[:,0] = np.array([id_dict[x] for x in sibpairs[:,0]])
sibpair_indices[:,1] = np.array([id_dict[x] for x in sibpairs[:,1]])
# LD scores
print('Reading LD scores from '+ldfile)
ld = np.loadtxt(ldfile,usecols=3,skiprows=1)
ldsnps = np.loadtxt(ldfile,usecols=1,skiprows=1,dtype=str)
ld_dict = make_id_dict(ldsnps)
in_ld_dict = np.array([x in ld_dict for x in bed.sid])
# Filtering on MAF
print('Before filtering on MAF and having an LD score, there were '+str(gts.shape[1])+' SNPs')
freq_pass = np.logical_and(np.logical_and(freqs > min_maf, freqs < 1-min_maf),in_ld_dict)
gts = gts[:,freq_pass]
freqs = freqs[freq_pass]
print('After filtering on MAF and having an LD score, there are '+str(gts.shape[1])+' SNPs')
if mapfile is None:
    bimfile = bedfile.split('.bed')[0]+'.bim'
    print('Separate genetic map not provided, so attempting to read map from '+bimfile)
    map = np.loadtxt(bimfile, usecols=2)
    chr = np.loadtxt(bimfile,usecols=0,dtype=str)
    chr = np.unique(chr)
    if chr.shape[0]>1:
        raise(ValueError('More than 1 chromosome in input bedfile'))
    else:
        chr = chr[0]
    # Check for NAs
    if np.sum(np.isnan(map))>0:
        raise(ValueError('Map cannot have NAs'))
    if np.min(map)<0:
        raise(ValueError('Map file cannot have negative values'))
    # Check ordering
    ordered_map = np.sort(map)
    if np.array_equal(map,ordered_map):
        pass
    else:
        raise(ValueError('Map not monotonic. Please make sure input is ordered correctly'))
    # Check scale
    if np.max(map)>5000:
        raise(ValueError('Maximum value of map too large'))
    map = map[freq_pass]
    print('Read map')
else:
    map = np.loadtxt(mapfile)
# Weights
weights = np.power(ld[np.array([ld_dict[x] for x in bed.sid[freq_pass]])],-1)
# IBD
print('Inferring IBD')
ibd = infer_ibd(sibpair_indices,gts,freqs,map,weights,p)
ibd, allsegs = smooth_ibd(ibd,map,p_length)
## Write output
snps = bed.sid[freq_pass]
# Write segments
segs_outfile = outprefix+'.ibd.segments.gz'
print('Writing segments to '+segs_outfile)
write_segs(sibpairs,allsegs,chr,snps,segs_outfile)
outfile = outprefix+'.ibd.gz'
print('Writing matrix output to '+str(outfile))
ibd = np.row_stack((np.column_stack((np.array(['sib1','sib2']).reshape((1,2)),snps.reshape(1,gts.shape[1]))),
                    np.column_stack((sibpairs,ibd))))
np.savetxt(outfile,ibd,fmt='%s')
