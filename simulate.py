import argparse
from bgen_reader import open_bgen
import numpy as np
from numba import set_num_threads
from numba import config as numba_config
import gzip, h5py, os
import snipar.preprocess as preprocess

parser = argparse.ArgumentParser()
parser.add_argument('h2_direct',type=float,help='Heritability due to direct effects in first generation',default=None)
parser.add_argument('outprefix',type=str,help='Prefix for simulation output files')
parser.add_argument('--bgenfiles', type=str,
                    help='Address of genotype files in .bgen format (without .bgen suffix). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of 1-22.',
                    default=None)
parser.add_argument('--n_random',type=int,help='Number of generations of random mating',default=1)
parser.add_argument('--n_am',type=int,help='Number of generations of assortative mating',default=0)
parser.add_argument('--r_par',type=float,help='Phenotypic correlation of parents (for assortative mating)',default=None)
parser.add_argument('--n_causal',type=int,help='Number of causal loci',default=None)
parser.add_argument('--beta_vert',type=float,help='Vertical transmission coefficient',default=0)
parser.add_argument('--h2_total',type=float,help='Total variance explained by direct effects and indirect genetic effects from parents',default=None)
parser.add_argument('--r_dir_indir',type=float,help='Correlation between direct and indirect genetic effects',default=None)
args=parser.parse_args()

if args.beta_vert>0 and args.h2_total is not None:
    raise(ValueError('Cannot simulate both indirect effects and vertical transmission separately. Choose one'))

if args.ngen_random>=0:
    ngen_random = args.ngen_random
else:
    raise(ValueError('Number of generations cannot be negative'))

if args.n_am>=0:
    ngen_am = args.n_am
else:
    raise(ValueError('Number of generations cannot be negative'))

if (args.r_par**2)<=1:
    r_y = args.r_par
else:
    raise(ValueError('Parental correlation must be between -1 and 1'))

if 0 <= args.h2_direct <= 1:
    h2 = args.h2_direct
else:
    raise(ValueError('Heritability must be between 0 and 1'))

if args.h2_total is not None:
    if args.r_dir_indir is None:
        raise(ValueError('Must specify correlation between direct and indirect genetic effects'))
    else:
        if (args.r_dir_indir**2)<=1:
            r_dir_alpha = args.r_dir_indir
        else:
            raise(ValueError('Correlation between direct and indirect effects must be between -1 and 1'))
    if 0 <= args.h2_total <= 1:
        h2_total = args.h2_total
    else:
        raise(ValueError('Heritability must be between 0 and 1'))


beta_vert = args.beta_vert
ncausal = args.n_causal

bgenfiles = preprocess.parse_obsfiles(args.bgenfiles, obsformat='bgen')

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