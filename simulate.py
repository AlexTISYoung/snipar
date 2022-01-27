import argparse
from bgen_reader import open_bgen
import numpy as np
from numba import set_num_threads
from numba import config as numba_config
import gzip, h5py, os
import snipar.preprocess as preprocess
from snipar.gtarray import gtarray


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

if args.beta_vert > 0 and args.h2_total is not None:
    raise(ValueError('Cannot simulate both indirect effects and vertical transmission separately. Choose one'))

if args.ngen_random >= 0:
    ngen_random = args.ngen_random
else:
    raise(ValueError('Number of generations cannot be negative'))

if args.n_am >= 0:
    ngen_am = args.n_am
else:
    raise(ValueError('Number of generations cannot be negative'))

if (args.r_par**2) <= 1:
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
        if (args.r_dir_indir**2) <= 1:
            r_dir_alpha = args.r_dir_indir
        else:
            raise(ValueError('Correlation between direct and indirect effects must be between -1 and 1'))
    if 0 <= args.h2_total <= 1:
        h2_total = args.h2_total
    else:
        raise(ValueError('Heritability must be between 0 and 1'))


beta_vert = args.beta_vert
ncausal = args.n_causal

haps = read_haps(args.bgenfiles)



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
    ped[(2 * i):(2 * (i + 1)),2] = 'P_'+str(i)
    ped[(2 * i):(2 * (i + 1)), 3] = 'M_' + str(i)
    ped[(2 * i):(2 * (i + 1)), 4] = np.array(['0','1'])
    ped[(2 * i):(2 * (i + 1)), 5] = np.array([Y_males[i],Y_females[i]],dtype=ped.dtype)

sibpairs = ped[:,1].reshape((int(ped.shape[0]/2),2))
ped = np.vstack((np.array(['FID','IID','FATHER_ID','MOTHER_ID','SEX','PHENO']),ped))
np.savetxt(outdir+'/sibs.ped',ped[:,0:4],fmt='%s')
np.savetxt(outdir+'/sibs.fam',ped[1:ped.shape[0],:],fmt='%s')

## Genotyping errors
p_error = 0.01

@njit(parallel=True)
def introduce_errors(haps,p_error):
    for i in prange(haps.shape[0]):
        for j in range(haps.shape[1]):
            for k in range(haps.shape[2]):
                for l in range(haps.shape[3]):
                    if np.random.binomial(n=1,p=p_error,size=1)==1:
                        haps[i,j,k,l] = np.logical_not(haps[i,j,k,l])
    return haps

print('Simulating genotyping errors with error probability '+str(p_error))
for chr in range(len(new_haps)):
    new_haps[chr] = introduce_errors(new_haps[chr])
    old_haps[chr] = introduce_errors(old_haps[chr])


## Save to HDF5 file
print('Saving genotypes to '+outdir+'/genotypes.hdf5')
hf = h5py.File(outdir+'/genotypes.hdf5','w')
# save pedigree
hf['ped'] = encode_str_array(ped)
# save offspring genotypes
chr_count = 0
for chr in range(chr_start,chr_end):
    print('Writing genotypes for chromosome '+str(chr))
    bimfile = 'bedfiles/chr_' + str(chr) + '.bim'
    bim = np.loadtxt(bimfile, dtype=str)
    # Sibling genotypes
    hf['bim_chr_'+str(chr)] = encode_str_array(bim)
    gts_chr = np.sum(new_haps[chr_count],axis=3,dtype=np.uint8)
    hf['chr_'+str(chr)] = gts_chr.reshape((gts_chr.shape[0]*2,gts_chr.shape[2]),order='C')
    # Imputed parental genotypes
    print('Imputing parental genotypes and saving')
    freqs = np.mean(gts_chr, axis=(0, 1)) / 2.0
    # Phased
    phased_imp = impute_all_fams_phased(new_haps[chr_count],freqs,ibd[chr_count])
    hf['phased_imp_chr'+str(chr)] = phased_imp
    # Unphased
    ibd[chr_count] = np.sum(ibd[chr_count],axis=2)
    imp = impute_all_fams(gts_chr, freqs, ibd[chr_count])
    hf['imp_chr_'+str(chr)] = imp
    # Parental genotypes
    print('Saving true parental genotypes')
    par_gts_chr = np.zeros((old_haps[chr_count].shape[0],2,old_haps[chr_count].shape[2]),dtype=np.uint8)
    par_gts_chr[:,0,:] = np.sum(old_haps[chr_count][father_indices,0,:,:],axis=2)
    par_gts_chr[:,1,:] = np.sum(old_haps[chr_count][mother_indices,1,:,:],axis=2)
    hf['par_chr_'+str(chr)] = par_gts_chr
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