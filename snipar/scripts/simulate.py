#!/usr/bin/env python
"""Simulates genotype-phenotype data using forward simulation. Phenotypes can be affected 
by direct genetic effects, indirect genetic effects (vertical transmission), and assortative mating.

Args:
@parser@

Results:
    genotype data in .bed format; full pedigree including phenotype and genetic components for all generations
        
"""
import numpy as np
import h5py, argparse
from snipar.ibd import write_segs_from_matrix
from snipar.map import decode_map_from_pos
from snipar.utilities import *
from snipar.simulate import *
from pysnptools.snpreader import SnpData,Bed
from snipar.utilities import get_parser_doc

parser = argparse.ArgumentParser()
parser.add_argument('n_causal',type=int,help='Number of causal loci')
parser.add_argument('h2',type=float,help='Heritability due to direct effects in first generation',default=None)
parser.add_argument('outprefix',type=str,help='Prefix for simulation output files')
parser.add_argument('--bgen',
                type=str,help='Address of the phased genotypes in .bgen format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--chr_range',
                type=parseNumRange,
                nargs='*',
                action=NumRangeAction,
                help='number of the chromosomes to be imputed. Should be a series of ranges with x-y format or integers.', default=None)
parser.add_argument('--nfam',type=int,help='Number of families to simulate. If inputting bgen and not given, will be one half of samples in bgen',default=None)
parser.add_argument('--min_maf',type=float,help='Minimum minor allele frequency for simulated genotyped, which will be simulted from density proportional to 1/x',default=0.05)
parser.add_argument('--maf',type=float,help='Minor allele frequency for simulated genotypes (not needed when providing bgen files)',default=None)
parser.add_argument('--n_random',type=int,help='Number of generations of random mating',default=0)
parser.add_argument('--n_am',type=int,help='Number of generations of assortative mating',default=0)
parser.add_argument('--r_par',type=float,help='Phenotypic correlation of parents (for assortative mating)',default=None)
parser.add_argument('--v_indir',type=float,help='Variance explained by parental indirect genetic effects as a fraction of the heritability, e.g 0.5',default=0)
parser.add_argument('--r_dir_indir',type=float,help='Correlation between direct and indirect genetic effects',default=None)
parser.add_argument('--beta_vert',type=float,help='Vertical transmission coefficient',default=0)
parser.add_argument('--save_par_gts',action='store_true',help='Save the genotypes of the parents of the final generation',default=False)
parser.add_argument('--impute',action='store_true',help='Impute parental genotypes from phased sibling genotypes & IBD',default=False)
parser.add_argument('--unphased_impute',action='store_true',help='Impute parental genotypes from unphased sibling genotypes & IBD',default=False)
__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))
def main(args):
    """"Calling this function with args is equivalent to running this script from commandline with the same arguments.
    Args:
        args: list
            list of all the desired options and arguments. The possible values are all the values you can pass this script from commandline.
    """
    if args.beta_vert > 0:
        if args.v_indir >0:
            raise(ValueError('Cannot simulate both indirect effects and vertical transmission separately. Choose one'))
        if args.beta_vert > np.power(2*(1+args.r_par),-0.5):
            raise(ValueError('Vertical transmission parameter would lead to infinite variance. Reduce to below 1/sqrt(2(1+r_par))'))
        if args.n_random == 0:
            print('Adding an additional initial generation of random-mating in order to initialize vertical transmission model')
            args.n_random = 1
    print('Simulating an initial generation by random-mating')
    print('Followed by '+str(args.n_random)+' generations of random-mating')
    print('Followed by '+str(args.n_am)+' generations of assortative mating')
    # Check parental correlation
    if args.n_am > 0:
        if args.r_par is not None:
            if (args.r_par**2) <= 1:
                r_y = args.r_par
            else:
                raise(ValueError('Parental correlation must be between -1 and 1'))
        else:
            raise(ValueError('Must specify parental phenotypic correlation for assortative mating'))
    # Check heritability
    if 0 <= args.h2 <= 1:
        h2 = args.h2
    else:
        raise(ValueError('Heritability must be between 0 and 1'))
    # Check v_indirect
    if not args.v_indir==0:
        if args.r_dir_indir is None:
            raise(ValueError('Must specify correlation between direct and indirect genetic effects'))
        else:
            if (args.r_dir_indir**2) > 1:
                raise(ValueError('Correlation between direct and indirect effects must be between -1 and 1'))
        if args.v_indir < 0:
            raise(ValueError('Heritability must be between 0 and 1'))

    #### Obtain haplotypes for base generation ###
    if args.bgen is None:
        if args.nfam is None:
            raise(ValueError('Must specify number of families to simulate'))
        unlinked = True
        haps, maps, snp_ids, alleles, positions, chroms = simulate_first_gen(args.nfam, args.n_causal, maf=None,min_maf=0.05)
    else:
        unlinked = False
        haps, maps, snp_ids, alleles, positions, chroms = haps_from_bgen(args.bgen, args.chr_range)
    ######################## Perform simulation ##################
    new_haps, haps, father_indices, mother_indices, ibd, ped, a, V = forward_sim(haps, maps, args.n_random, args.n_am, unlinked, args.n_causal, args.h2, 
                                                                    v_indirect=args.v_indir, r_direct_indirect=args.r_dir_indir, r_y=args.r_par, beta_vert=args.beta_vert)
    ##############################################################
    ### Save variance components ###
    print('Saving variance components to '+args.outprefix+'VCs.txt')
    np.savetxt(args.outprefix+'VCs.txt',V,fmt='%s')
    ### Save pedigree and fam files
    print('Writing pedigree and fam files')
    n_last = ped[ped.shape[0]-1,0].split('_')[0]
    last_gen = [x.split('_')[0]==n_last for x in ped[:,0]]
    np.savetxt(args.outprefix+'pedigree.txt',ped,fmt='%s')
    #np.savetxt(args.outprefix+'sibs.fam',ped[last_gen,:],fmt='%s')
    phen_out = ped[last_gen,:]
    np.savetxt(args.outprefix+'phenotype.txt',phen_out[:,[0,1,5]],fmt='%s')
    ## Save to HDF5 file
    print('Saving genotypes')
    # save offspring genotypes
    for i in range(len(new_haps)):
        print('Writing genotypes for chromosome '+str(chroms[i]))
        bim_i = np.vstack((np.repeat(chroms[i],snp_ids[i].shape[0]),snp_ids[i], maps[i], positions[i], alleles[i][:,0],alleles[i][:,1])).T
        gts_chr = SnpData(iid=ped[last_gen,0:2],sid=snp_ids[i],pos=bim_i[:,[0,2,3]],
                          val=np.sum(new_haps[i], axis=3, dtype=np.uint8).reshape((new_haps[i].shape[0]*2,new_haps[i].shape[2])))
        Bed.write(args.outprefix+'chr_'+str(chroms[i])+'.bed',gts_chr,count_A1=True,_require_float32_64=False)
        np.savetxt(args.outprefix+'chr_'+str(chroms[i])+'.bim',bim_i,fmt='%s',delimiter='\t')
        if args.save_par_gts:
            par_gen = [x.split('_')[0]==str(int(n_last)-1) for x in ped[:,0]]
            par_gts_chr = SnpData(iid=ped[par_gen,0:2],sid=snp_ids[i],pos=bim_i[:,[0,2,3]],
                          val=np.sum(haps[i], axis=3, dtype=np.uint8).reshape((haps[i].shape[0]*2,haps[i].shape[2])))
            Bed.write(args.outprefix+'chr_'+str(chroms[i])+'_par.bed',par_gts_chr,count_A1=True,_require_float32_64=False) 
            del par_gts_chr
            np.savetxt(args.outprefix+'chr_'+str(chroms[i])+'_par.bim',bim_i,fmt='%s',delimiter='\t')
        # # Imputed parental genotypes
        if args.impute or args.unphased_impute:
            print('Imputing parental genotypes and saving')
            freqs = np.mean(gts_chr.val, axis=0) / 2.0
            imp_ped = ped[last_gen,0:4]
            imp_ped = np.hstack((imp_ped,np.zeros((imp_ped.shape[0],2),dtype=bool)))
            # Phased
            if args.impute:
                hf = h5py.File(args.outprefix+'phased_impute_chr_'+str(chroms[i])+'.hdf5','w')
                phased_imp = impute_all_fams_phased(new_haps[i],freqs,ibd[i])
                hf['imputed_par_gts'] = phased_imp
                del phased_imp
                hf['bim_values'] = encode_str_array(bim_i)
                hf['bim_columns'] = encode_str_array(np.array(['rsid','map','position','allele1','allele2']))
                hf['pedigree'] = encode_str_array(imp_ped)
                hf['families'] = encode_str_array(imp_ped[0::2,0])
                hf.close()
            # Unphased
            if args.unphased_impute:
                hf = h5py.File(args.outprefix+'unphased_impute_chr_'+str(chroms[i])+'.hdf5','w')
                ibd[i] = np.sum(ibd[i],axis=2)
                imp = impute_all_fams(gts_chr, freqs, ibd[i])
                hf['imputed_par_gts'] = imp
                del imp
                hf['bim_values'] = encode_str_array(bim_i)
                hf['bim_columns'] = encode_str_array(np.array(['rsid','map','position','allele1','allele2']))
                hf['pedigree'] = encode_str_array(imp_ped)
                hf['families'] = encode_str_array(imp_ped[0::2,0])
                hf.close()
        del gts_chr
    #### Write IBD segments
    snp_count = 0
    sibpairs = ped[last_gen,1]
    sibpairs = sibpairs.reshape((int(sibpairs.shape[0]/2),2))
    if args.v_indir==0:
        causal_out = np.zeros((a.shape[0],5),dtype='U30')
    else:
        causal_out = np.zeros((a.shape[0],5),dtype='U30')
    for i in range(len(haps)):
        print('Writing IBD segments for chromosome '+str(chroms[i]))
        # Segments
        if not args.unphased_impute:
            ibd[i] = np.sum(ibd[i],axis=2)
        segs = write_segs_from_matrix(ibd[i], sibpairs,
                                    snp_ids[i], positions[i],maps[i],chroms[i],
                                    args.outprefix+'chr_'+str(chroms[i])+'.segments.gz')
        # Causal effects
        if args.v_indir==0:
            a_chr = a[snp_count:(snp_count + snp_ids[i].shape[0])]
            a_chr_v1 = a_chr+np.random.normal(0,np.std(a_chr),a_chr.shape)
            causal_out[snp_count:(snp_count + snp_ids[i].shape[0]),:] = np.vstack((snp_ids[i],alleles[i][:,0],alleles[i][:,1],
                                                                                    a_chr,a_chr_v1)).T
            if i==0:
                causal_out = np.vstack((np.array(['SNP','A1','A2','direct','direct_v1']).reshape((1,5)),causal_out))
        else:
            a_chr = a[snp_count:(snp_count + snp_ids[i].shape[0]),:]
            causal_out[snp_count:(snp_count + snp_ids[i].shape[0]),:] = np.vstack((snp_ids[i],alleles[i][:,0],alleles[i][:,1],
                                                                                    a_chr[:,0],a_chr[:,1])).T
            if i==0:
                causal_out = np.vstack((np.array(['SNP','A1','A2','direct','indirect']).reshape((1,5)),causal_out))
        snp_count += snp_ids[i].shape[0]
        # count
    np.savetxt(args.outprefix+'causal_effects.txt',causal_out,fmt='%s')
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)
