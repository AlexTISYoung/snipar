#!/usr/bin/env python
from bgen_reader import open_bgen
import numpy as np
import h5py, argparse
from snipar.ibd import write_segs_from_matrix
from snipar.map import decode_map_from_pos
from snipar.utilities import *
from snipar.simulate import *
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bgen',
                    type=str,help='Address of the phased genotypes in .bgen format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
    parser.add_argument('--chr_range',
                    type=parseNumRange,
                    nargs='*',
                    action=NumRangeAction,
                    help='number of the chromosomes to be imputed. Should be a series of ranges with x-y format or integers.', default=None)
    parser.add_argument('n_causal',type=int,help='Number of causal loci')
    parser.add_argument('outprefix',type=str,help='Prefix for simulation output files')
    parser.add_argument('--n_random',type=int,help='Number of generations of random mating',default=1)
    parser.add_argument('--n_am',type=int,help='Number of generations of assortative mating',default=0)
    parser.add_argument('--r_par',type=float,help='Phenotypic correlation of parents (for assortative mating)',default=None)
    parser.add_argument('--h2_direct',type=float,help='Heritability due to direct effects in first generation',default=None)
    parser.add_argument('--h2_total',type=float,help='Total variance explained by direct effects and indirect genetic effects from parents',default=None)
    parser.add_argument('--r_dir_indir',type=float,help='Correlation between direct and indirect genetic effects',default=None)
    parser.add_argument('--beta_vert',type=float,help='Vertical transmission coefficient',default=0)
    parser.add_argument('--save_par_gts',action='store_true',help='Save the genotypes of the parents of the final generation',default=False)
    parser.add_argument('--impute',action='store_true',help='Impute parental genotypes from phased sibling genotypes & IBD',default=False)
    parser.add_argument('--unphased_impute',action='store_true',help='Impute parental genotypes from unphased sibling genotypes & IBD',default=False)
    args=parser.parse_args()

    if args.beta_vert > 0 and args.h2_total is not None:
        raise(ValueError('Cannot simulate both indirect effects and vertical transmission separately. Choose one'))
    if args.beta_vert > 0 and args.h2_direct is None:
        raise(ValueError('Must provide h2_direct if simulating vertical transmission'))
    if args.h2_direct is None and args.h2_total is None:
        raise(ValueError('Must provide one of h2_direct or h2_total'))
    if args.h2_direct is not None and args.h2_total is not None:
        raise(ValueError('Must provide one of h2_direct or h2_total'))

    if args.n_random >= 0:
        ngen_random = args.n_random
    else:
        raise(ValueError('Number of generations cannot be negative'))

    if args.n_am > 0:
        if args.r_par is not None:
            if (args.r_par**2) <= 1:
                r_y = args.r_par
            else:
                raise(ValueError('Parental correlation must be between -1 and 1'))
        else:
            raise(ValueError('Must specify parental phenotypic correlation for assortative mating'))

    ngen_am = args.n_am

    if args.h2_direct is not None:
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
    else:
        h2_total = 0


    beta_vert = args.beta_vert
    ncausal = args.n_causal

    bgenfiles, chroms = parse_obsfiles(args.bgen, obsformat='bgen', chromosomes=args.chr_range)

    # Read genotypes
    haps = []
    maps = []
    snp_ids = []
    alleles = []
    positions = []
    for i in range(bgenfiles.shape[0]):
        print('Reading in chromosome '+str(chroms[i]))
        # Read bgen
        bgen = open_bgen(bgenfiles[i], verbose=True)
        # Read map
        map = decode_map_from_pos(chroms[i], bgen.positions)
        not_nan = np.logical_not(np.isnan(map))
        maps.append(map[not_nan])
        positions.append(bgen.positions[not_nan])
        # Snp
        snp_ids.append(bgen.rsids[not_nan])
        # Alleles
        alleles.append(np.array([x.split(',') for x in bgen.allele_ids[not_nan]]))
        # Read genotypes
        gts = bgen.read(([x for x in range(bgen.samples.shape[0])], [x for x in range(bgen.ids.shape[0])]), np.bool_)[:, :,
            np.array([0, 2])]
        gts = gts[:,not_nan]
        nfam = int(np.floor(gts.shape[0] / 2))
        ngen = np.zeros((nfam, 2, gts.shape[1], 2), dtype=np.bool_)
        ngen[:, 0, :, :] = gts[0:nfam, :, :]
        ngen[:, 1, :, :] = gts[nfam:(2 * nfam), :, :]
        del gts
        haps.append(ngen)

    # Simulate population
    total_matings = ngen_random+ngen_am
    if h2_total == 0:
        V = np.zeros((total_matings+1,2))
    else:
        V = np.zeros((total_matings,2))
    a_count = 0
    # Produce next generation
    old_haps = haps
    for gen in range(0, total_matings):
        # Simulate phenotype for AM
        if gen == 0 and h2_total == 0:
            print('Simulating phenotype')
            # Simulate phenotype
            nsnp_chr = np.array([x.shape[2] for x in old_haps])
            nsnp = np.sum(nsnp_chr)
            if ncausal > nsnp:
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
            father_indices = random_mating_indices(nfam)
            mother_indices = random_mating_indices(nfam)
        # Assortative mating
        if gen>=ngen_random:
            # Compute parental phenotypes
            print('Computing parental phenotypes')
            # Match assortatively
            print('Matching assortatively')
            father_indices, mother_indices = am_indices(Y_p, Y_m, r_y)
            # Print variances
            print('Parental phenotypic correlation: '+str(round(np.corrcoef(Y_p[father_indices],Y_m[mother_indices])[0,1],4)))
            print('Parental genotypic correlation: '+str(round(np.corrcoef(G_p[father_indices],G_m[mother_indices])[0,1],4)))
        # Generate haplotpyes of new generation
        new_haps = []
        ibd = []
        for chr in range(0,len(haps)):
            print('Chromosome '+str(chroms[0]+chr))
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
    vcf = args.outprefix+'VCs.txt'
    np.savetxt(vcf, V)
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
    np.savetxt(args.outprefix+'sibs.ped',ped[:,0:4],fmt='%s')
    np.savetxt(args.outprefix+'sibs.fam',ped[1:ped.shape[0],:],fmt='%s')

    ## Save to HDF5 file
    print('Saving genotypes to '+args.outprefix+'genotypes.hdf5')
    hf = h5py.File(args.outprefix+'genotypes.hdf5','w')
    # save pedigree
    hf['ped'] = encode_str_array(ped)
    # save offspring genotypes
    for i in range(len(haps)):
        print('Writing genotypes for chromosome '+str(chroms[i]))
        bim_i = np.vstack((snp_ids[i], maps[i], positions[i], alleles[i][:,0],alleles[i][:,1])).T
        hf['chr_'+str(chroms[i])+'_bim'] = encode_str_array(bim_i)
        gts_chr = np.sum(new_haps[i], axis=3, dtype=np.uint8)
        hf['chr_'+str(chroms[i])] = gts_chr.reshape((gts_chr.shape[0]*2,gts_chr.shape[2]),order='C')
        # # Imputed parental genotypes
        if args.impute or args.unphased_impute:
            print('Imputing parental genotypes and saving')
            freqs = np.mean(gts_chr, axis=(0, 1)) / 2.0
            # Phased
            if args.impute:
                phased_imp = impute_all_fams_phased(new_haps[i],freqs,ibd[i])
                hf['imp_chr'+str(chroms[i])] = phased_imp
            # Unphased
            if args.unphased_impute:
                ibd[i] = np.sum(ibd[i],axis=2)
                imp = impute_all_fams(gts_chr, freqs, ibd[i])
                hf['unphased_imp_chr_'+str(chroms[i])] = imp
            # Parental genotypes
        if args.save_par_gts:
            print('Saving true parental genotypes')
            par_gts_chr = np.zeros((old_haps[i].shape[0],2,old_haps[i].shape[2]),dtype=np.uint8)
            par_gts_chr[:,0,:] = np.sum(old_haps[i][father_indices,0,:,:],axis=2)
            par_gts_chr[:,1,:] = np.sum(old_haps[i][mother_indices,1,:,:],axis=2)
            hf['par_chr_'+str(chroms[i])] = par_gts_chr

    hf.close()

    # Write IBD segments
    snp_count = 0

    if h2_total>0:
        causal_out = np.zeros((ab.shape[0], 5), dtype='U30')
    else:
        causal_out = np.zeros((a.shape[0],4),dtype='U30')
    for i in range(len(haps)):
        print('Writing IBD segments for chromosome '+str(chroms[i]))
        # Segments
        if not args.unphased_impute:
            ibd[i] = np.sum(ibd[i],axis=2)
        segs = write_segs_from_matrix(ibd[i], sibpairs,
                                    snp_ids[i], positions[i],maps[i],chroms[i],
                                    args.outprefix+'chr_'+str(chroms[i])+'.segments.gz')
        # Causal effects
        if h2_total==0:
            a_chr = a[snp_count:(snp_count + snp_ids[i].shape[0])]
            causal_out[snp_count:(snp_count + snp_ids[i].shape[0]),:] = np.vstack((snp_ids[i],alleles[i][:,0],alleles[i][:,1],
                                                                                    a_chr)).T
        else:
            ab_chr = ab[snp_count:(snp_count + snp_ids[i].shape[0]),:]
            causal_out[snp_count:(snp_count + snp_ids[i].shape[0]),:] = np.vstack((snp_ids[i],alleles[i][:,0],alleles[i][:,1],
                                                                                    ab_chr[:,0],ab_chr[:,1])).T
        snp_count += snp_ids[i].shape[0]
        # count

    np.savetxt(args.outprefix+'causal_effects.txt',causal_out,fmt='%s')