#!/usr/bin/env python
"""Simulates genotype-phenotype data using forward simulation. Phenotypes can be affected 
by direct genetic effects, indirect genetic effects (vertical transmission), assortative mating,
or multigeneration Gaussian stabilizing selection on breeders.

**Usage**

    ``simulate.py 100 0.5 selected_ --nfam 200 --n_stabilizing 3 --selection_variance 2 --save_par_gts``

Args:
@parser@

Results:
    Genotypes and phenotypes
        Final-generation genotypes in PLINK format and ``phenotype.txt`` with the offspring phenotypes.
        With ``--save_par_gts``, the simulation also writes the parents’ PLINK files.

    Pedigree and genetic effects
        ``pedigree.txt`` records phenotypes and genetic components across generations.
        ``causal_effects.txt`` contains the simulated locus effects.
        With ``--n_stabilizing`` greater than zero, ``selection_metadata.json`` records the selection settings and outcomes.
        All filenames begin with the chosen ``outprefix``.
        
"""
import numpy as np
import h5py, argparse, json
from snipar.ibd import write_segs_from_matrix
from snipar.utilities import *
from snipar.simulate import *
from pysnptools.snpreader import SnpData,Bed
from snipar.numrange import parseNumRange, NumRangeAction
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
parser.add_argument('--n_stabilizing', type=int, default=0,
                    help='Generations of Gaussian viability selection; final offspring are sampled before selection')
parser.add_argument('--selection_variance', type=float, default=None,
                    help='Gaussian fitness variance in baseline phenotype units; inf gives a matched neutral control')
parser.add_argument('--selection_optimum', type=float, default=0.0,
                    help='Fixed optimum relative to the initial mean genetic value')
parser.add_argument('--selection_nfam', type=int, default=None,
                    help='Breeding-pool pairs before selection; defaults to four times --nfam')
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
    selection = args.n_stabilizing != 0
    if selection:
        if args.n_stabilizing < 1 or args.selection_variance is None or args.nfam is None:
            raise ValueError('Selection requires positive --n_stabilizing, --selection_variance and --nfam')
        if args.n_random or args.n_am or args.v_indir or args.beta_vert:
            raise ValueError('Selection currently supports direct effects only, without additional random/assortative generations')
        if args.selection_nfam is not None and args.selection_nfam < args.nfam:
            raise ValueError("The selection breeding pool must contain at least --nfam pairs")
        if args.bgen and args.selection_nfam is not None:
            raise ValueError('The breeding pool is determined by the input haplotypes when using --bgen')
        stabilizing_selection_weights(np.zeros(1), args.selection_variance, args.selection_optimum)
    elif args.selection_variance is not None or args.selection_nfam is not None:
        raise ValueError('Selection parameters require --n_stabilizing')
    if args.n_causal < 1 or (args.nfam is not None and args.nfam < 1):
        raise ValueError("Number of causal loci and families must be positive")
    if args.n_random < 0 or args.n_am < 0:
        raise ValueError("Numbers of mating generations must be non-negative")
    if args.r_par is None and args.n_am:
        raise ValueError("Must specify parental phenotypic correlation for assortative mating")
    if args.r_par is not None and not -1 <= args.r_par <= 1:
        raise ValueError("Parental correlation must be between -1 and 1")
    if not 0 <= args.h2 <= 1:
        raise ValueError("Heritability must be between 0 and 1")
    if not np.isfinite(args.beta_vert):
        raise ValueError("Vertical transmission coefficient must be finite")
    if args.beta_vert:
        if args.v_indir:
            raise ValueError("Cannot simulate both indirect effects and vertical transmission separately. Choose one")
        parental_correlation = args.r_par if args.n_am else 0
        if 2 * (1 + parental_correlation) * args.beta_vert**2 >= 1:
            raise ValueError("Vertical transmission parameter would lead to infinite variance.")
        if args.n_random == 0:
            print("Adding an additional initial generation of random-mating in order to initialize vertical transmission model")
            args.n_random = 1
    print("Simulating an initial generation by random-mating")
    print(f"Followed by {args.n_random} generations of random-mating")
    print(f"Followed by {args.n_am} generations of assortative mating")
    if not np.isfinite(args.v_indir) or args.v_indir < 0:
        raise ValueError("Indirect genetic variance must be finite and non-negative")
    if args.v_indir:
        if args.r_dir_indir is None:
            raise ValueError("Must specify correlation between direct and indirect genetic effects")
        if not -1 <= args.r_dir_indir <= 1:
            raise ValueError("Correlation between direct and indirect effects must be between -1 and 1")

    #### Obtain haplotypes for base generation ###
    if args.bgen is None:
        if args.nfam is None:
            raise(ValueError('Must specify number of families to simulate'))
        unlinked = True
        n_breeding = args.nfam
        if selection:
            n_breeding = args.selection_nfam if args.selection_nfam is not None else 4 * args.nfam
        haps, maps, snp_ids, alleles, positions, chroms = simulate_first_gen(
            n_breeding, args.n_causal, maf=args.maf, min_maf=args.min_maf)
    else:
        unlinked = False
        haps, maps, snp_ids, alleles, positions, chroms = haps_from_bgen(args.bgen, args.chr_range, None if selection else args.nfam)
    ######################## Perform simulation ##################
    if selection:
        new_haps, haps, father_indices, mother_indices, ibd, ped, a, V, selection_metadata = forward_sim_stabilizing_selection(
            haps, maps, args.n_stabilizing, args.n_causal, args.h2, args.selection_variance,
            args.nfam, unlinked=unlinked, optimum=args.selection_optimum)
        if not np.isfinite(selection_metadata['selection_variance']):
            selection_metadata['selection_variance'] = None
            selection_metadata['neutral_control'] = True
        with open(args.outprefix + 'selection_metadata.json', 'w') as handle:
            json.dump(selection_metadata, handle, indent=2, allow_nan=False)
    else:
        new_haps, haps, father_indices, mother_indices, ibd, ped, a, V = forward_sim(haps, maps, args.n_random, args.n_am, unlinked, args.n_causal, args.h2,
                                                                        v_indirect=args.v_indir, r_direct_indirect=args.r_dir_indir, r_y=args.r_par, beta_vert=args.beta_vert)
    ##############################################################
    ### Save variance components ###
    print('Saving variance components to '+args.outprefix+'VCs.txt')
    np.savetxt(args.outprefix+'VCs.txt',V,fmt='%s')
    ### Save pedigree and fam files
    print('Writing pedigree and last generation relatedness and phenotype files')
    n_last = ped[ped.shape[0]-1,0].split('_')[0]
    last_gen = [x.split('_')[0]==n_last for x in ped[:,0]]
    np.savetxt(args.outprefix+'pedigree.txt',ped,fmt='%s')
    rel = (selection_family_relationships(ped, parent_siblings=args.n_stabilizing > 1)
           if selection else find_sibs_and_cousins(ped))
    np.savetxt(args.outprefix+'last_gen_rel.seg',rel,fmt='%s',delimiter='\t')
    #np.savetxt(args.outprefix+'sibs.fam',ped[last_gen,:],fmt='%s')
    phen_out = ped[last_gen,:]
    phen_out = np.vstack((np.array(['FID','IID','phenotype']).reshape((1,3)),phen_out[:,[0,1,5]]))
    np.savetxt(args.outprefix+'phenotype.txt',phen_out,fmt='%s')
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
                hf['bim_columns'] = encode_str_array(np.array(['chrom','rsid','map','position','allele1','allele2']))
                hf['pedigree'] = encode_str_array(imp_ped)
                hf['families'] = encode_str_array(imp_ped[0::2,0])
                hf.close()
            # Unphased
            if args.unphased_impute:
                ibd[i] = np.sum(ibd[i],axis=2)
                family_dosages = gts_chr.val.reshape(new_haps[i].shape[:3])
                imp = impute_all_fams(family_dosages, freqs, ibd[i])
                del family_dosages
                hf = h5py.File(args.outprefix+'unphased_impute_chr_'+str(chroms[i])+'.hdf5','w')
                hf['imputed_par_gts'] = imp
                del imp
                hf['bim_values'] = encode_str_array(bim_i)
                hf['bim_columns'] = encode_str_array(np.array(['chrom','rsid','map','position','allele1','allele2']))
                hf['pedigree'] = encode_str_array(imp_ped)
                hf['families'] = encode_str_array(imp_ped[0::2,0])
                hf.close()
        del gts_chr
    #### Write IBD segments
    sibpairs = ped[last_gen,1]
    sibpairs = sibpairs.reshape((int(sibpairs.shape[0]/2),2))
    for i in range(len(haps)):
        print('Writing IBD segments for chromosome '+str(chroms[i]))
        # Segments
        if not args.unphased_impute:
            ibd[i] = np.sum(ibd[i],axis=2)
        write_segs_from_matrix(ibd[i], sibpairs,
                                    snp_ids[i], positions[i],maps[i],chroms[i],
                                    args.outprefix+'chr_'+str(chroms[i])+'.segments.gz', return_segments=False)
    header = 'SNP A1 A2 direct ' + ('direct_v1' if args.v_indir == 0 else 'indirect')
    with open(args.outprefix+'causal_effects.txt', 'w') as output:
        output.write(header + '\n')
        start = 0
        for ids, chrom_alleles in zip(snp_ids, alleles):
            end = start + len(ids)
            effects = a[start:end]
            if args.v_indir == 0:
                noisy = effects + np.random.normal(0, np.std(effects), effects.shape)
                effects = np.column_stack((effects, noisy))
            np.savetxt(output, np.column_stack((ids, chrom_alleles.astype(str), effects)), fmt='%s')
            start = end
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)
