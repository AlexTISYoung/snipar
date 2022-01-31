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
def produce_next_gen(father_indices,mother_indices,males,females,map):
    ngen = np.zeros((father_indices.shape[0],2,males.shape[1],2),dtype=np.bool_)
    ibd = np.zeros((father_indices.shape[0],males.shape[1]),dtype=np.int_)
    for i in prange(father_indices.shape[0]):
        # recombinations
        recomb_i = meiosis(map,n=4)
        # sib haplotypes and ibd
        for j in range(ibd.shape[1]):
            # paternal sib 1
            if recomb_i[0,j]:
                ngen[i, 0, j, 0] =  males[father_indices[i], j, 0]
            else:
                ngen[i, 0, j, 0] = males[father_indices[i], j, 1]
            # paternal sib 2
            if recomb_i[1,j]:
                ngen[i, 1, j, 0] =  males[father_indices[i], j, 0]
            else:
                ngen[i, 1, j, 0] = males[father_indices[i], j, 1]
            # maternal sib 1
            if recomb_i[2,j]:
                ngen[i, 0, j, 1] =  females[mother_indices[i], j, 0]
            else:
                ngen[i, 0, j, 1] = females[mother_indices[i], j, 1]
            # maternal sib 2
            if recomb_i[3,j]:
                ngen[i, 1, j, 1] =  females[mother_indices[i], j, 0]
            else:
                ngen[i, 1, j, 1] = females[mother_indices[i], j, 1]
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
    unmatched = []
    # Match within pop
    for pop in pops:
        ids_in_pop = np.random.permutation(ids[np.where(pop_membership==pop)[0]])
        npairs = np.int(np.floor(ids_in_pop.shape[0]/2.0))
        father_indices += [id_dict[x] for x in ids_in_pop[0:npairs]]
        mother_indices += [id_dict[x] for x in ids_in_pop[npairs:(2*npairs)]]
        if ids_in_pop.shape[0]%2 == 1:
            unmatched.append(id_dict[ids_in_pop[2*npairs]])
    # Match unmatched
    unmatched = np.random.permutation(unmatched)
    n_unmatched_pairs = np.int(np.floor(len(unmatched)/2.0))
    father_indices += list(unmatched[0:n_unmatched_pairs])
    mother_indices += list(unmatched[n_unmatched_pairs:(2*n_unmatched_pairs)])
    if len(father_indices)==len(mother_indices) and 2*len(father_indices)==ids.shape[0]:
        pass
    else:
        raise(ValueError('Error in matching'))
    return np.array(father_indices), np.array(mother_indices)

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
ngen_am = 0
r_y = 0.5
ncausal = 1500
h2 = 0.5
beta_vert = 0
h2_total = 0.75
r_dir_alpha = 0.0
chr_start = 20
chr_end = 22
structure_file = '/disk/genetics/ukb/alextisyoung/phenotypes/assessment_centre.txt'
pop_var = 0.3

if h2_total>0:
    simname = 'h2_' + str(h2_total)+'_r_dir_alpha_'+str(r_dir_alpha)
else:
    simname = 'h2_'+str(h2)

if ngen_am>0:
    simname += '_amgen_'+str(ngen_am)+'_ry_'+str(r_y)

if np.abs(beta_vert)>0:
    simname += '_vert_'+str(beta_vert)

if len(structure_file)>0:
    simname += '_pop_var_'+str(pop_var)

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
    nfam = int(np.floor(gts.shape[0] / 2))
    ngen = np.zeros((nfam, 2, gts.shape[1], 2), dtype=np.bool_)
    ngen[:, 0, :, :] = gts[0:nfam, :, :]
    ngen[:, 1, :, :] = gts[nfam:(2 * nfam), :, :]
    del gts
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

# Simulate population
total_matings = ngen_random+ngen_am
V = np.zeros((total_matings+1,2))
a_count = 0
# Produce next generation
old_haps = haps
for gen in range(0,total_matings):
    # Simulate phenotype for AM
    if gen==0 and h2_total==0:
        print('Simulating phenotype')
        # Simulate phenotype
        nsnp_chr = np.array([x.shape[2] for x in old_haps])
        nsnp = np.sum(nsnp_chr)
        if ncausal>nsnp:
            raise(ValueError('Not enough SNPs to simulate phenotype with '+str(ncausal)+' causal SNPs'))
        a = np.zeros((nsnp))
        causal = np.sort(np.random.choice(np.arange(0,nsnp),ncausal,replace=False))
        a[causal] = np.random.randn(ncausal)
        G_p, G_m = compute_genetic_component(old_haps,causal,a)
        scale_fctr = np.sqrt(h2/np.var(np.hstack((G_p,G_m))))
        a = a*scale_fctr
    # Compute parental phenotypes
    if np.abs(beta_vert) > 0 and gen>0:
        print('Computing parental phenotypes')
        G_p, G_m, Y_p, Y_m = compute_phenotype_vert(old_haps, causal, a, 1-h2, beta_vert, Y_p[father_indices], Y_m[mother_indices])
    elif h2_total==0:
        print('Computing parental phenotypes')
        G_p, G_m, Y_p, Y_m = compute_phenotype(old_haps, causal, a, 1 - h2)
    # Record variance components
    if gen>0 or h2_total==0:
        V[a_count, :] = np.array([np.var(np.hstack((G_p, G_m))), np.var(np.hstack((Y_p, Y_m)))])
    print('Genetic variance: ' + str(round(V[a_count, 0], 4)))
    print('Phenotypic variance: ' + str(round(V[a_count, 1], 4)))
    print('Heritability: ' + str(round(V[a_count, 0] / V[a_count, 1], 4)))
    a_count += 1
    ## Match parents
    print('Mating ' + str(gen + 1))
    # Random mating
    if gen<ngen_random:
        if n_pops == 1:
            father_indices = random_mating_indices(nfam)
            mother_indices = random_mating_indices(nfam)
        else:
            father_indices, mother_indices = structured_random_mating(ids,pop_membership)
    # Assortative mating
    if gen>=ngen_random:
        # Compute parental phenotypes
        print('Computing parental phenotypes')
        # Match assortatively
        print('Matching assortatively')
        father_indices, mother_indices = am_indices(Y_p, Y_m, 0.5)
        # Print variances
        print('Parental phenotypic correlation: '+str(round(np.corrcoef(Y_p[father_indices],Y_m[mother_indices])[0,1],4)))
        print('Parental genotypic correlation: '+str(round(np.corrcoef(G_p[father_indices],G_m[mother_indices])[0,1],4)))
    # Generate haplotpyes of new generation
    new_haps = []
    ibd = []
    for chr in range(0,len(haps)):
        print('Chromosome '+str(chr_start+chr))
        new_haps_chr, ibd_chr = produce_next_gen(father_indices,mother_indices,old_haps[chr][:,0,:,:],old_haps[chr][:,1,:,:],maps[chr])
        new_haps.append(new_haps_chr)
        ibd.append(ibd_chr)
    # Compute indirect effect component
    if h2_total>0:
        if gen==0:
            print('Simulating indirect genetic effects')
            nsnp_chr = np.array([x.shape[2] for x in old_haps])
            nsnp = np.sum(nsnp_chr)
            if ncausal > nsnp:
                raise (ValueError('Not enough SNPs to simulate phenotype with ' + str(ncausal) + ' causal SNPs'))
            ab = np.zeros((nsnp,2))
            causal = np.sort(np.random.choice(np.arange(0, nsnp), ncausal, replace=False))
            ab[causal,:] = np.random.multivariate_normal(np.zeros((2)),
                                                       np.array([[1,r_dir_alpha],[r_dir_alpha,1]]),
                                                       size=ncausal)
            G_males, G_females, Y_males, Y_females = compute_phenotype_indirect(new_haps,old_haps,father_indices,mother_indices,causal,ab[:,0],ab[:,1],0)
            scale_fctr = np.sqrt(h2_total / np.var(np.hstack((Y_males, Y_females))))
            ab = ab*scale_fctr
        print('Computing parental phenotype')
        G_p, G_m, Y_p, Y_m = compute_phenotype_indirect(new_haps,old_haps,father_indices,mother_indices,causal,ab[:,0],ab[:,1],1-h2_total)
    if gen<(total_matings-1):
        old_haps = new_haps

print('Computing final generation phenotypes')
if np.abs(beta_vert)>0:
    G_males, G_females, Y_males, Y_females = compute_phenotype_vert(new_haps, causal, a, 1 - h2, beta_vert, Y_p[father_indices], Y_m[mother_indices])
elif h2_total>0:
    G_males, G_females, Y_males, Y_females = compute_phenotype_indirect(new_haps, old_haps, father_indices, mother_indices, causal,
                                                    ab[:, 0], ab[:, 1], 1 - h2_total)
else:
    G_males, G_females, Y_males, Y_females = compute_phenotype(new_haps, causal, a, 1 - h2)
print('Sibling genotypic correlation: ' + str(round(np.corrcoef(G_males, G_females)[0, 1], 4)))
print('Sibling phenotypic correlation: ' + str(round(np.corrcoef(Y_males, Y_females)[0, 1], 4)))
# Final offspring generation
V[a_count,:] = np.array([np.var(np.hstack((G_males,G_females))),np.var(np.hstack((Y_males, Y_females)))])

print('Saving output to file')
# Save variance
vcf = outdir+'/VCs.txt'
np.savetxt(outdir+'/VCs.txt',V)
print('Variance components saved to '+str(vcf))
print('Saving pedigree')
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

if h2_total>0:
    causal_out = np.zeros((ab.shape[0], 8), dtype='U30')
else:
    causal_out = np.zeros((a.shape[0],7),dtype='U30')
os.mkdir(outdir+'/ibd/')
for chr in range(chr_start,chr_end):
    print('Writing IBD segments for chromosome '+str(chr))
    bimfile = 'bedfiles/chr_'+str(chr)+'.bim'
    bim = np.loadtxt(bimfile, dtype=str)
    # Segments
    segs = write_segs_from_matrix(ibd[chr_count],sibpairs,bim[:,1],bimfile,outdir+'/ibd/chr_'+str(chr)+'.segments.gz')
    # Causal effects
    if h2_total==0:
        a_chr = a[snp_count:(snp_count + bim.shape[0])]
        causal_out[snp_count:(snp_count + bim.shape[0]),:] = np.hstack((bim,a_chr.reshape((bim.shape[0],1))))
    else:
        ab_chr = ab[snp_count:(snp_count + bim.shape[0]),:]
        causal_out[snp_count:(snp_count + bim.shape[0]),:] = np.hstack((bim,ab_chr.reshape((bim.shape[0],2))))
    snp_count += bim.shape[0]
    # count
    chr_count += 1

np.savetxt(outdir+'/causal_effects.txt',causal_out,fmt='%s')