#!/usr/bin/env python
"""Infers direct effects, non-transmitted coefficients (NTCs), and population effects of a PGS on a phenotype.

Minimally: the script requires observed genotypes on  individuals and their parents, and/or 
parental genotypes imputed by snipar's impute.py script, along with a SNP weights file.

Args:
@parser@

Results:
    PGS file
        Output when inputting observed and imputed genotype files and a weights file. 
        A file with PGS values for each individual and their parents, with suffix .pgs.txt.
        Also includes sibling PGS if --fit_sib is specified, and grandparental PGS if --grandpar is specified.
    PGS effect estimates
        Output when inputting a phenotype file.
        A file with suffix effects.txt containing estimates of the PGS effects and their standard errors,
        and a file with suffix vcov.txt containing the sampling variance-covariance matrix of the effect estimates
        
"""
import argparse
import numpy as np
import snipar.pgs as pgs
import snipar.read as read
from snipar.utilities import *
from snipar.slmm import build_ibdrel_arr, build_sib_arr
from snipar.utilities import get_parser_doc
######### Command line arguments #########
parser=argparse.ArgumentParser()
parser.add_argument('out',type=str,help='Prefix for computed PGS file and/or regression results files')
parser.add_argument('--bgen',
                    type=str,help='Address of the phased genotypes in .bgen format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--bed',
                    type=str,help='Address of the unphased genotypes in .bed format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--imp', type=str, help='Address of hdf5 files with imputed parental genotypes (without .hdf5 suffix). If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range (chr_range is an optional parameters for this script).', default = None)
parser.add_argument('--chr_range',
                    type=parseNumRange,
                    nargs='*',
                    action=NumRangeAction,
                    help='number of the chromosomes to be imputed. Should be a series of ranges with x-y format or integers.', default=None)
parser.add_argument('--pedigree',type=str,help='Address of pedigree file. Must be provided if not providing imputed parental genotypes.',default=None)
parser.add_argument('--weights',type=str,help='Location of the PGS allele weights', default = None)
parser.add_argument('--SNP',type=str,help='Name of column in weights file with SNP IDs',default='SNP')
parser.add_argument('--beta_col',type=str,help='Name of column with betas/weights for each SNP',default='b')
parser.add_argument('--A1',type=str,help='Name of column with allele beta/weights are given with respect to',default='A1')
parser.add_argument('--A2',type=str,help='Name of column with alternative allele',default='A2')
parser.add_argument('--sep',type=str,help='Column separator in weights file. If not provided, an attempt to determine this will be made.',default=None)
parser.add_argument('--phenofile',type=str,help='Location of the phenotype file',default = None)
parser.add_argument('--pgs', type=str, help='Location of the pre-computed PGS file', default=None)
parser.add_argument('--covar',type=str,help='Path to file with covariates: plain text file with columns FID, IID, covar1, covar2, ..', default=None)
parser.add_argument('--fit_sib', action='store_true', default=False, help='Fit indirect effects from siblings')
parser.add_argument('--parsum',action='store_true',default = False, help='Use the sum of maternal and paternal PGS in the regression (useful when imputed from sibling data alone)')
parser.add_argument('--grandpar',action='store_true',default=False,help='Calculate imputed/observed grandparental PGS for individuals with both parents genotyped')
parser.add_argument('--gparsum',action='store_true',default = False, help='Use the sum of maternal grandparents and the sum of paternal grandparents in the regression (useful when no grandparents genotyped)')
parser.add_argument('--gen_models',
                    type=parseNumRange,
                    nargs='*',
                    action=NumRangeAction,
                    help='Which multi-generational models should be fit. Default fits 1 and 2 generation models. Specify a range by, for example, 1-3, where 3 fits a model with parental and grandparental scores', default='1-2')
parser.add_argument('--h2f',type=str,help='Provide heritability estimate in form h2f,h2f_SE (e.g. 0.5,0.01) from MZ-DZ comparison, RDR, or sibling realized relatedness. If provided when also fitting 2 generation model, will adjust results for assortative mating assuming equilibrium.',default=None)
parser.add_argument('--rk',type=str,help='Provide estimate of the correlation between parents PGIs in the form rk,rk_SE (e.g 0.1,0.01). If provided with h2f, will use for adjusting estimates for assortative mating.',default=None)
parser.add_argument('--bpg',action='store_true', default=False, help='Restrict sample to those with both parents genotyped')    
parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                    default=1)
parser.add_argument('--ibdrel_path', type=str,
                    help='Path to KING IBD segment inference output (without .seg prefix).', default=None)
parser.add_argument('--sparse_thresh', type=float,
                help='Threshold of GRM/IBD sparsity', default=0.05)
parser.add_argument('--scale_phen',action='store_true',help='Scale the phenotype to have variance 1',default=False)
parser.add_argument('--scale_pgs',action='store_true',help='Scale the PGS to have variance 1 among the phenotyped individuals',default=False)
parser.add_argument('--compute_controls', action='store_true', default=False,
                    help='Compute PGS for control families (default False)')
parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)',default='NA')
parser.add_argument('--no_am_adj',action='store_true',help='Do not adjust imputed parental PGSs for assortative mating',default=False)
parser.add_argument('--force_am_adj',action='store_true',help='Force assortative mating adjustment even when estimated correlation is noisy/not significant',default=False)
parser.add_argument('--batch_size',type=int,help='Batch size for reading in SNPs (default 10000)',default=10000)
__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))
def main(args):
    """"Calling this function with args is equivalent to running this script from commandline with the same arguments.
    Args:
        args: list
            list of all the desired options and arguments. The possible values are all the values you can pass this script from commandline.
    """
    if args.weights is not None:
        if args.bed is None and args.bgen is None:
            raise ValueError('Weights provided but no observed genotypes provided')
        if args.bed is not None and args.bgen is not None:
            raise ValueError('Provide only one of --bedfiles and --bgenfiles')
        print('Computing PGS from weights file')
        ####### Read PGS #######
        p = pgs.read_weights(args.weights, SNP=args.SNP, beta_col = args.beta_col, A1=args.A1, A2=args.A2, sep=args.sep)
        # Remove zeros
        p.remove_zeros()

        ###### Compute PGS ########
        # Find observed and imputed files
        if args.imp is None:
            print('Warning: no imputed parental genotypes provided. Will compute PGS only for individuals with both parents genotyped.')
            if args.bed is not None:
                bedfiles, chroms = parse_obsfiles(args.bed, 'bed', chromosomes=args.chr_range)
                bgenfiles = [None for x in range(chroms.shape[0])]
            elif args.bgen is not None:
                bgenfiles, chroms = parse_obsfiles(args.bgen, 'bgen', chromosomes=args.chr_range)
                bedfiles = [None for x in range(chroms.shape[0])]
            pargts_list = [None for x in range(chroms.shape[0])]
        else:
            if args.bed is not None:
                bedfiles, pargts_list, chroms = parse_filelist(args.bed, args.imp, 'bed', chromosomes=args.chr_range)
                bgenfiles = [None for x in range(chroms.shape[0])]
            elif args.bgen is not None:
                bgenfiles, pargts_list, chroms = parse_filelist(args.bgen, args.imp, 'bgen', chromosomes=args.chr_range)
                bedfiles = [None for x in range(chroms.shape[0])]
        if chroms.shape[0]==0:
            raise(ValueError('No input genotype files found'))
        # Get pedigree if no imputed parental genotypes provided
        if args.imp is None:
            if args.pedigree is None:
                raise(ValueError('Must provide pedigree if not providing imputed parental genotypes'))
            print('Reading pedigree from '+str(args.pedigree))
            ped = np.loadtxt(args.pedigree,dtype=str)
            if ped.shape[1] < 4:
                raise(ValueError('Not enough columns in pedigree file'))
            elif ped.shape[1] > 4:
                print('Warning: pedigree file has more than 4 columns. The first four columns only will be used')
        else:
            ped = None
        print('Computing PGS')
        if args.bed is not None:
            print('Observed genotypes file: '+bedfiles[0])
        if args.bgen is not None:
            print('Observed genotypes file: '+bgenfiles[0])
        if args.imp is not None:
            print('Imputed genotypes file: '+pargts_list[0])
        pg = pgs.compute(p, bedfile=bedfiles[0], bgenfile=bgenfiles[0], par_gts_f=pargts_list[0], ped=ped, sib=args.fit_sib, compute_controls=args.compute_controls, batch_size=args.batch_size)
        for i in range(1,chroms.shape[0]):
            if args.bed is not None:
                print('Observed genotypes file: '+bedfiles[i])
            if args.bgen is not None:
                print('Observed genotypes file: '+bgenfiles[i])
            if args.imp is not None:
                print('Imputed genotypes file: '+pargts_list[i])
            if args.compute_controls:
                pg_i = pgs.compute(p, bedfile=bedfiles[i], bgenfile=bgenfiles[i], par_gts_f=pargts_list[i], ped=ped, sib=args.fit_sib, compute_controls=args.compute_controls, batch_size=args.batch_size)
                if pg_i[0] is not None:
                    if pg is None:
                        pg = pg_i
                    else:
                        pg = [pg[x].add(pg_i[x]) for x in range(0, len(pg))]
            else:
                pg_i = pgs.compute(p, bedfile=bedfiles[i], bgenfile=bgenfiles[i], par_gts_f=pargts_list[i], ped=ped, sib=args.fit_sib, compute_controls=args.compute_controls, batch_size=args.batch_size)
                if pg_i is not None:
                    if pg is not None:
                        pg = pg.add(pg_i)
                    else:
                        pg = pg_i
        print('PGS computed')
        ####### Assortative mating adjustment #######
        if not args.no_am_adj:
            if args.compute_controls:
                r_am = pg[0].am_adj()
            else:
                r_am = pg.am_adj()
        else:
            r_am = 0
        ####### Compute grandparental PGSs #######
        if args.grandpar:
            if args.compute_controls:
                pg[0].compute_grandpar(r_am)
            else:
                pg.compute_grandpar(r_am)
        ####### Write PGS to file ########
        if args.compute_controls:
            pg[0].write(args.out + '.pgs.txt', scale=args.scale_pgs)
            pg[1].write(args.out + '.pgs.control_paternal.txt', scale=args.scale_pgs)
            pg[2].write(args.out + '.pgs.control_maternal.txt', scale=args.scale_pgs)
            pg[3].write(args.out + '.pgs.control_sibling.txt', scale=args.scale_pgs)
        else:
            pg.write(args.out + '.pgs.txt', scale=args.scale_pgs)
    elif args.pgs is not None:
        if args.phenofile is None:
            raise ValueError('Pre-computed PGS provided but no phenotype provided')
        print('Reading PGS from '+args.pgs)
        pg = pgs.read_pgs(args.pgs)
        if not args.parsum:
            if 'parental' in pg.sid and 'paternal' not in pg.sid:
                print('Using sum of paternal and maternal PGS as no paternal & maternal PGS values found')
                args.parsum = True
    else:
        raise ValueError('Weights or PGS must be provided')

    if args.phenofile is not None:
        print('Fitting PGS for '+str(args.phenofile))
        # Read phenotype
        y = read.phenotype.read_phenotype(args.phenofile, missing_char=args.missing_char, phen_index=args.phen_index)
        print('Number of non-missing phenotype observations: ' + str(y.shape[0]))
        if args.covar is not None:
            print('Reading covariates')
            covariates = read.phenotype.read_covariates(args.covar, pheno_ids=y.ids, missing_char=args.missing_char)
            # Match to pheno ids
            covariates.filter_ids(y.ids)
        else:
            covariates = None
        # Restrict sample to both parents genotyped
        if args.bpg:
            print('Restricting to individuals with both parents genotyped')
            pg.filter_bpg()
        # Remove individuals without phenotype observations from PGS
        # and match IDs
        pg.filter_ids(y.ids)
        y.filter_ids(pg.ids)
        if covariates is not None:
            covariates.filter_ids(pg.ids)
        print('Sample size of individuals with complete phenotype and PGS observations: '+str(y.shape[0]))
        # Scale
        if args.scale_phen:
            y.scale()
        if args.scale_pgs:
            pg.scale()
        ## Estimate models
        if '1' in args.gen_models:
            print('Estimating population effect (1 generation model)')
            alpha_1 = pgs.fit_pgs_model(y, pg, 1, ibdrel_path=args.ibdrel_path, covariates=covariates, fit_sib=args.fit_sib, parsum=args.parsum, gparsum=args.gparsum, outprefix=args.out, sparse_thresh=args.sparse_thresh)
        if '2' in args.gen_models:
            print('Estimating direct effect and parental NTCs (2 generation model)')
            alpha_2 = pgs.fit_pgs_model(y, pg, 2, ibdrel_path=args.ibdrel_path, covariates=covariates, fit_sib=args.fit_sib, parsum=args.parsum, gparsum=args.gparsum, outprefix=args.out, sparse_thresh=args.sparse_thresh)
            if args.h2f is not None:
                print('Adjusting two-generation model results and heritability estimate for assortative mating')
                h2f, h2f_se = pgs.h2f_parse(args.h2f)
                if args.rk is None:
                    adj_estimates, adj_ses = pgs.am_adj_2gen(alpha_2[0], alpha_2[1], h2f, h2f_se, pg=pg, y_std=np.std(y.gts[:,0]), pg_std=np.std(pg.gts[:,np.where(pg.sid=='proband')[0][0]]))
                else:
                    rk, rk_se = pgs.h2f_parse(args.rk)
                    adj_estimates, adj_ses = pgs.am_adj_2gen(alpha_2[0], alpha_2[1], h2f, h2f_se, rk=rk, rk_se=rk_se, y_std=np.std(y.gts[:,0]), pg_std=np.std(pg.gts[:,np.where(pg.sid=='proband')[0][0]]))
                pgs.write_2gen_adj_ests(adj_estimates, adj_ses, outprefix=args.out)
        if '3' in args.gen_models:
            print('Estimating direct effect and parental IGEs and grandparental coefficients (3 generation model)')
            alpha_3 = pgs.fit_pgs_model(y, pg, 3, ibdrel_path=args.ibdrel_path, covariates=covariates, fit_sib=args.fit_sib, parsum=args.parsum, gparsum=args.gparsum, outprefix=args.out, sparse_thresh=args.sparse_thresh)
    return None
if __name__ == "__main__":
    args=parser.parse_args()
    main(args)