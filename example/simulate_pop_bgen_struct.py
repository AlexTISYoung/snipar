from bgen_reader import open_bgen
import numpy as np
from numba import njit, prange, set_num_threads
from numba import config as numba_config
from sibreg.sibreg import *
import gzip, h5py, os

@njit
def simulate_recombinations(map):
    map_start = map[0]
    map_end = map[map.shape[0]-1]
    map_length = map_end-map_start
    n_recomb = np.random.poisson(map_length/100)
    recomb_points = map_start+np.sort(np.random.uniform(0,1,n_recomb))*map_length
    return n_recomb,recomb_points


@njit
def meiosis(map,n=1):
    # Recomb vector
    recomb_vector = np.zeros((n,map.shape[0]), dtype=np.bool_)
    # Do recombinations
    for r in range(n):
        n_recomb, recomb_points = simulate_recombinations(map)
        # Initialize
        last_hap = np.bool_(np.random.binomial(1,0.5,1)[0])
        recomb_vector[r,:] = np.bool_(last_hap)
        # Recombine
        if n_recomb>0:
            for i in range(n_recomb):
                first_snp = np.argmax(map>recomb_points[i])
                recomb_vector[r,first_snp:recomb_vector.shape[1]] = ~recomb_vector[r,first_snp:recomb_vector.shape[1]]
    # Return
    return recomb_vector


@njit(parallel=True)
def produce_next_gen(father_indices,mother_indices,gts,map):
    ngen = np.zeros((father_indices.shape[0],2,gts.shape[1],2),dtype=np.bool_)
    ibd = np.zeros((father_indices.shape[0],gts.shape[1]),dtype=np.int_)
    for i in prange(father_indices.shape[0]):
        # recombinations
        recomb_i = meiosis(map,n=4)
        # sib haplotypes and ibd
        for j in range(ibd.shape[1]):
            # paternal sib 1
            if recomb_i[0,j]:
                ngen[i, 0, j, 0] =  gts[father_indices[i], j, 0]
            else:
                ngen[i, 0, j, 0] = gts[father_indices[i], j, 1]
            # paternal sib 2
            if recomb_i[1,j]:
                ngen[i, 1, j, 0] =  gts[father_indices[i], j, 0]
            else:
                ngen[i, 1, j, 0] = gts[father_indices[i], j, 1]
            # maternal sib 1
            if recomb_i[2,j]:
                ngen[i, 0, j, 1] =  gts[mother_indices[i], j, 0]
            else:
                ngen[i, 0, j, 1] = gts[mother_indices[i], j, 1]
            # maternal sib 2
            if recomb_i[3,j]:
                ngen[i, 1, j, 1] =  gts[mother_indices[i], j, 0]
            else:
                ngen[i, 1, j, 1] = gts[mother_indices[i], j, 1]
            ibd[i,j] = int(recomb_i[0,j]==recomb_i[1,j])+int(recomb_i[2,j]==recomb_i[3,j])
    return ngen, ibd

@njit
def random_mating_indices(nfam):
    return np.random.choice(np.array([x for x in range(nfam)],dtype=np.int_),size=nfam,replace=False)

def am_indices(yp,ym,r):
    v = np.sqrt(np.var(yp)*np.var(ym))
    s2 = (1/r-1)*v
    zp = yp+np.sqrt(s2)*np.random.randn(yp.shape[0])
    zm = ym+np.sqrt(s2)*np.random.randn(ym.shape[0])
    rank_p = np.argsort(zp)
    rank_m = np.argsort(zm)
    return rank_p, rank_m

def structured_random_mating(ids,pop_membership):
    if ids.shape[0] == pop_membership.shape[0]:
        pass
    else:
        raise(ValueError('IDs and pop membership must be of same size'))
    pops = np.unique(pop_membership)
    id_dict = make_id_dict(ids)
    father_indices = []
    mother_indices = []
    ngen_pop_membership = []
    # Match within pop
    for pop in pops:
        ids_in_pop = np.random.permutation(ids[np.where(pop_membership==pop)[0]])
        npairs = np.int(np.floor(ids_in_pop.shape[0]/2.0))
        father_indices += [id_dict[x] for x in ids_in_pop[0:npairs]]
        mother_indices += [id_dict[x] for x in ids_in_pop[npairs:(2*npairs)]]
        ngen_pop_membership += [pop for x in range(npairs)]
    return np.array(father_indices,dtype=np.int_), np.array(mother_indices,dtype=np.int_), np.array(ngen_pop_membership)

def compute_genetic_component(haps,causal,a):
    causal = set(causal)
    G_m = np.zeros((haps[0].shape[0]))
    G_p = np.zeros((haps[0].shape[0]))
    snp_count = 0
    for chr in range(len(haps)):
        snp_index = np.array([snp_count+x for x in range(haps[chr].shape[2])])
        in_causal = np.array([snp_index[x] in causal for x in range(haps[chr].shape[2])])
        causal_gts = np.sum(haps[chr][:,:,in_causal,:],axis=3)
        G_p += causal_gts[:, 0, :].dot(a[snp_index[in_causal]])
        G_m += causal_gts[:, 1, :].dot(a[snp_index[in_causal]])
        snp_count += haps[chr].shape[2]
    return G_p, G_m

def compute_phenotype(haps,causal,a,sigma2):
    G_p, G_m = compute_genetic_component(haps,causal,a)
    Y_p = G_p+np.sqrt(sigma2)*np.random.randn(G_p.shape[0])
    Y_m = G_m+np.sqrt(sigma2)*np.random.randn(G_m.shape[0])
    return G_p, G_m, Y_p, Y_m

def compute_phenotype_struct(haps, pop_membership, ncausal, h2, pop_var):
    # Simulate effects
    nsnp_chr = np.array([x.shape[2] for x in haps])
    nsnp = np.sum(nsnp_chr)
    if ncausal > nsnp:
        raise (ValueError('Not enough SNPs to simulate phenotype with ' + str(ncausal) + ' causal SNPs'))
    a = np.zeros((nsnp))
    causal = np.sort(np.random.choice(np.arange(0, nsnp), ncausal, replace=False))
    a[causal] = np.random.randn(ncausal)
    G_p, G_m = compute_genetic_component(haps, causal, a)
    scale_fctr = np.sqrt(h2 / np.var(np.hstack((G_p, G_m))))
    a = a * scale_fctr
    # Simulate phenotype
    G_p, G_m = compute_genetic_component(haps,causal,a)
    # Add in pop specific environmental effect
    pops = np.unique(pop_membership,return_inverse=True)
    pop_effects = np.random.randn(pops[0].shape[0])
    pop_comp = pop_effects[pops[1]]
    pop_comp = np.sqrt(pop_var/np.var(pop_comp))*pop_comp
    Y_p = G_p+pop_comp+np.sqrt(1-h2-pop_var)*np.random.randn(G_p.shape[0])
    Y_m = G_m+pop_comp+np.sqrt(1-h2-pop_var)*np.random.randn(G_m.shape[0])
    return a, causal, G_p, G_m, Y_p, Y_m

def compute_phenotype_vert(haps,causal,a,sigma2,beta_vert,Y_p,Y_m):
    G_males, G_females = compute_genetic_component(haps, causal, a)
    Y_males = G_males + beta_vert*(Y_p+Y_m)+np.sqrt(sigma2) * np.random.randn(G_males.shape[0])
    Y_females = G_females + beta_vert*(Y_p+Y_m)+np.sqrt(sigma2) * np.random.randn(G_females.shape[0])
    return G_males, G_females, Y_males, Y_females

def compute_phenotype_indirect(haps,old_haps,father_indices,mother_indices,causal,a,b,sigma2):
    G_males, G_females = compute_genetic_component(haps,causal,a)
    G_p, G_m = compute_genetic_component(old_haps, causal, b)
    Y_males = G_males+G_p[father_indices]+G_m[mother_indices]+np.sqrt(sigma2)*np.random.randn(G_males.shape[0])
    Y_females = G_females+G_p[father_indices]+G_m[mother_indices]+np.sqrt(sigma2)*np.random.randn(G_males.shape[0])
    return G_males, G_females, Y_males, Y_females

def write_segs_from_matrix(ibd,sibpairs,snps,bimfile,outfile):
    # Get map and chr from bim
    bim = np.loadtxt(bimfile, dtype=str)
    pos = np.loadtxt(bimfile,usecols=3,dtype=int)
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
        allsegs.append(find_segments(ibd[i,:],map,snps,pos))
    # Write segments
    write_segs(sibpairs,allsegs,chr,outfile)
    return allsegs

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

class segment(object):
    def __init__(self, start_index, end_index, start_bp, end_bp, start_snp, end_snp, length, state):
        self.start = start_index
        self.end = end_index
        self.start_bp = start_bp
        self.end_bp = end_bp
        self.start_snp = start_snp
        self.end_snp = end_snp
        self.length = length
        self.state = state
    def to_text(self, id1, id2, chr, end=False):
        seg_txt = str(id1) + '\t' + str(id2) + '\t' + str(self.state) + '\t' + str(chr) + '\t' + str(
            self.start_bp) + '\t' + str(self.end_bp) + '\t' + str(self.start_snp) + '\t' + str(
            self.end_snp) + '\t' + str(self.length)
        if end:
            return seg_txt
        else:
            return seg_txt + '\n'

ngen_random = 1
ncausal = 1500
h2 = 0.5
chr_start = 20
chr_end = 22
structure_file = '/disk/genetics/ukb/alextisyoung/phenotypes/assessment_centre.txt'
pop_var = 0.3

simname = 'h2_' + str(h2)+'_pop_var_'+str(pop_var)

outdir = '/disk/genetics/ukb/alextisyoung/haplotypes/hm3/simulations/'+simname
os.mkdir(outdir)

# Read genotypes
haps = []
maps = []
snp_ids = []
alleles = []
for chr in range(chr_start,chr_end):
    print('Reading in chromosome '+str(chr))
    # Read map
    map = np.loadtxt('bedfiles/chr_'+str(chr)+'.bim', dtype=str)
    maps.append(np.array(map[:, 2], dtype=float))
    # Read bgen
    bgen = open_bgen('QCed/chr_'+str(chr)+'.bgen', verbose=True)
    # Snp
    snp_ids.append(bgen.ids)
    # Alleles
    alleles.append(np.array([x.split(',') for x in bgen.allele_ids]))
    # Read genotypes
    gts = bgen.read(([x for x in range(bgen.samples.shape[0])], [x for x in range(bgen.ids.shape[0])]), np.bool_)[:, :,
          np.array([0, 2])]
    ngen = np.zeros(gts.shape,dtype=np.bool_)
    ngen[:] = gts
    haps.append(ngen)

ids = bgen.samples
ids = np.array([x.split('_')[0] for x in ids])

if len(structure_file)>0:
    print('Loading population memberships')
    pop_membership = np.loadtxt(structure_file,skiprows=1,dtype=str)
    membership_dict = make_id_dict(pop_membership[:,0])
    pop_membership = pop_membership[[membership_dict[x] for x in ids],1]
else:
    pop_membership = np.ones((ids.shape[0]),dtype=str)

pops = np.unique(pop_membership)
n_pops = pops.shape[0]
print('Using '+str(n_pops)+' sub populations')

# Produce next generation
father_indices, mother_indices, pop_membership = structured_random_mating(ids,pop_membership)
# Generate haplotpyes of new generation
new_haps = []
ibd = []
for chr in range(0,len(haps)):
    print('Chromosome '+str(chr_start+chr))
    new_haps_chr, ibd_chr = produce_next_gen(father_indices,mother_indices,haps[chr],maps[chr])
    new_haps.append(new_haps_chr)
    ibd.append(ibd_chr)

print('Simulating structured phenotype')
a, causal, G_males, G_females, Y_males, Y_females = compute_phenotype_struct(new_haps, pop_membership, ncausal, h2, pop_var)
print('Sibling genotypic correlation: ' + str(round(np.corrcoef(G_males, G_females)[0, 1], 4)))
print('Sibling phenotypic correlation: ' + str(round(np.corrcoef(Y_males, Y_females)[0, 1], 4)))

# Final offspring generation
print('Saving output to file')
print('Saving pedigree')
nfam = G_males.shape[0]
# Produce pedigree
ped = np.zeros((nfam*2,6),dtype='U30')
for i in range(0,nfam):
    ped[(2*i):(2*(i+1)),0] = i
    ped[(2 * i):(2 * (i + 1)), 1] = np.array([str(i)+'_0',str(i)+'_1'],dtype=str)
    ped[(2 * i):(2 * (i + 1)),2] = 'M_'+str(father_indices[i])
    ped[(2 * i):(2 * (i + 1)), 3] = 'F_' + str(mother_indices[i])
    ped[(2 * i):(2 * (i + 1)), 4] = np.array(['0','1'])
    ped[(2 * i):(2 * (i + 1)), 5] = np.array([Y_males[i],Y_females[i]],dtype=ped.dtype)

sibpairs = ped[:,1].reshape((int(ped.shape[0]/2),2))
ped = np.vstack((np.array(['FID','IID','FATHER_ID','MOTHER_ID','SEX','PHENO']),ped))
np.savetxt(outdir+'/sibs.ped',ped[:,0:4],fmt='%s')
np.savetxt(outdir+'/sibs.fam',ped[1:ped.shape[0],:],fmt='%s')

## Save to HDF5 file
print('Saving offspring genotypes to '+outdir+'/genotypes.hdf5')
hf = h5py.File(outdir+'/genotypes.hdf5','w')
# save pedigree
hf['ped'] = encode_str_array(ped)
# save offspring genotypes
chr_count = 0
for chr in range(chr_start,chr_end):
    print('Writing genotypes for chromosome '+str(chr))
    bimfile = 'bedfiles/chr_' + str(chr) + '.bim'
    bim = np.loadtxt(bimfile, dtype=str)
    hf['bim_chr_'+str(chr)] = encode_str_array(bim)
    gts_chr = np.sum(new_haps[chr_count],axis=3,dtype=np.uint8)
    hf['chr_'+str(chr)] = gts_chr.reshape((gts_chr.shape[0]*2,gts_chr.shape[2]),order='C')
    chr_count += 1

hf.close()

# Write IBD segments
snp_count = 0
chr_count = 0

causal_out = np.zeros((a.shape[0],7),dtype='U30')
os.mkdir(outdir+'/ibd/')
for chr in range(chr_start,chr_end):
    print('Writing IBD segments for chromosome '+str(chr))
    bimfile = 'bedfiles/chr_'+str(chr)+'.bim'
    bim = np.loadtxt(bimfile, dtype=str)
    # Segments
    segs = write_segs_from_matrix(ibd[chr_count],sibpairs,bim[:,1],bimfile,outdir+'/ibd/chr_'+str(chr)+'.segments.gz')
    # Causal effects
    a_chr = a[snp_count:(snp_count + bim.shape[0])]
    causal_out[snp_count:(snp_count + bim.shape[0]),:] = np.hstack((bim,a_chr.reshape((bim.shape[0],1))))
    snp_count += bim.shape[0]
    # count
    chr_count += 1

np.savetxt(outdir+'/causal_effects.txt',causal_out,fmt='%s')