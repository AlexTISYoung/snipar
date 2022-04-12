#!/usr/bin/env python
import argparse
import numpy as np
import snipar.pgs as pgs
from snipar.gtarray import gtarray
import snipar.read as read
import snipar.lmm as lmm
from snipar.utilities import *

######### Command line arguments #########
if __name__ == '__main__':
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
    parser.add_argument('--SNP',type=str,help='Name of column in weights file with SNP IDs',default='sid')
    parser.add_argument('--beta_col',type=str,help='Name of column with betas/weights for each SNP',default='ldpred_beta')
    parser.add_argument('--A1',type=str,help='Name of column with allele beta/weights are given with respect to',default='nt1')
    parser.add_argument('--A2',type=str,help='Name of column with alternative allele',default='nt2')
    parser.add_argument('--sep',type=str,help='Column separator in weights file. If not provided, an attempt to determine this will be made.',default=None)
    parser.add_argument('--phenofile',type=str,help='Location of the phenotype file',default = None)
    parser.add_argument('--pgs', type=str, help='Location of the pre-computed PGS file', default=None)
    parser.add_argument('--fit_sib', action='store_true', default=False, help='Fit indirect effects from siblings')
    parser.add_argument('--parsum',action='store_true',default = False, help='Use the sum of maternal and paternal PGS in the regression (useful when imputed from sibling data alone)')
    parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                        default=1)
    parser.add_argument('--scale_phen',action='store_true',help='Scale the phenotype to have variance 1',default=False)
    parser.add_argument('--scale_pgs',action='store_true',help='Scale the PGS to have variance 1 among the phenotyped individuals',default=False)
    parser.add_argument('--compute_controls', action='store_true', default=False,
                        help='Compute PGS for control families (default False)')
    parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)',default='NA')
    args=parser.parse_args()

    if args.weights is not None:
        if args.bed is None and args.bgen is None:
            raise ValueError('Weights provided but no observed genotypes provided')
        if args.bed is not None and args.bgen is not None:
            raise ValueError('Provide only one of --bedfiles and --bgenfiles')
        print('Computing PGS from weights file')
        ####### Read PGS #######
        if args.sep is None:
            weights = np.loadtxt(args.weights,dtype=str)
        else:
            weights = np.loadtxt(args.weights,dtype=args.sep)
        colnames = weights[0,:]
        weights = weights[1:weights.shape[0],:]
        print('Read weights for '+str(weights.shape[0])+' variants')
        beta = np.array(weights[:,np.where(colnames == args.beta_col)[0][0]],dtype=np.float64)
        allele_indices = np.array([np.where(colnames==args.A1)[0][0],np.where(colnames==args.A2)[0][0]])
        p = pgs.pgs(weights[:,np.where(colnames==args.SNP)[0][0]],
                beta,
                weights[:,allele_indices])

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
        pg = pgs.compute(p, bedfile=bedfiles[0], bgenfile=bgenfiles[0], par_gts_f=pargts_list[0], ped=ped, sib=args.fit_sib, compute_controls=args.compute_controls)
        for i in range(1,chroms.shape[0]):
            if args.compute_controls:
                pg_i = pgs.compute(p, bedfile=bedfiles[i], bgenfile=bgenfiles[i], par_gts_f=pargts_list[i], ped=ped, sib=args.fit_sib, compute_controls=args.compute_controls)
                pg = [pg[x].add(pg_i[x]) for x in range(0, len(pg))]
            else:
                pg = pg.add(pgs.compute(p, bedfile=bedfiles[i], bgenfile=bgenfiles[i], par_gts_f=pargts_list[i], ped=ped, sib=args.fit_sib, compute_controls=args.compute_controls))
        print('PGS computed')
        ####### Write PGS to file ########
        if args.compute_controls:
            pgs.write(pg[0], args.out + '.pgs.txt', scale_PGS=args.scale_pgs)
            pgs.write(pg[1],args.out + '.pgs.control_paternal.txt', scale_PGS=args.scale_pgs)
            pgs.write(pg[2], args.out + '.pgs.control_maternal.txt', scale_PGS=args.scale_pgs)
            pgs.write(pg[3],args.out + '.pgs.control_sibling.txt', scale_PGS=args.scale_pgs)
        else:
            pgs.write(pg, args.out + '.pgs.txt', scale_PGS=args.scale_pgs)
    elif args.pgs is not None:
        if args.phenofile is None:
            raise ValueError('Pre-computed PGS provided but no phenotype provided')
        print('Reading PGS from '+args.pgs)
        pgs_f = open(args.pgs,'r')
        pgs_header = pgs_f.readline().split(' ')
        pgs_header[len(pgs_header)-1] = pgs_header[len(pgs_header)-1].split('\n')[0]
        ncols = len(pgs_header)
        pgs_cols = tuple([x for x in range(2,ncols)])
        pg = gtarray(np.loadtxt(args.pgs,usecols = pgs_cols, skiprows=1),
                     np.loadtxt(args.pgs,usecols = 1, dtype=str, skiprows=1),
                     sid=np.array(pgs_header[2:ncols]),
                     fams=np.loadtxt(args.pgs,usecols = 0, dtype=str, skiprows=1))
    else:
        raise ValueError('Weights or PGS must be provided')

    if args.phenofile is not None:
        print('Fitting PGS for '+str(args.phenofile))
        # Read phenotype
        y = read.phenotype.read_phenotype(args.phenofile, missing_char=args.missing_char, phen_index=args.phen_index)
        print('Number of non-missing phenotype observations: ' + str(y.shape[0]))
        # Remove individuals without phenotype observations from PGS
        pg.filter_ids(y.ids)
        y.filter_ids(pg.ids)
        print('Final sample size of individuals with complete phenotype and PGS observations: '+str(y.shape[0]))
        # Parental sum
        if args.parsum:
            if 'maternal' in pg.sid and 'paternal' in pg.sid:
                parcols = np.sort(np.array([np.where(pg.sid=='maternal')[0][0],np.where(pg.sid=='paternal')[0][0]]))
                trans_matrix = np.identity(pg.gts.shape[1])
                trans_matrix[:,parcols[0]] += trans_matrix[:,parcols[1]]
                trans_matrix = np.delete(trans_matrix,parcols[1],1)
                pg.gts = pg.gts.dot(trans_matrix)
                pg.sid = np.delete(pg.sid,parcols[1])
                pg.sid[parcols[0]] = 'parental'
            else:
                raise(ValueError('Maternal and paternal PGS not found so cannot sum (--parsum option given)'))
        # Scale
        if args.scale_phen:
            y.scale()
        if args.scale_pgs:
            pg.scale()
        # Estimate effects
        print('Estimating direct effects and NTCs')
        alpha_imp = lmm.fit_model(y.gts[:,0], pg.gts, pg.fams, add_intercept=True, return_model=False, return_vcomps=False)
        # Estimate population effect
        print('Estimating population effect')
        alpha_proband = lmm.fit_model(y.gts[:,0], pg.gts[:, 0], pg.fams, add_intercept=True, return_model=False, return_vcomps=False)
        # Get print out for fixed mean effects
        alpha_out = np.zeros((pg.sid.shape[0]+1, 2))
        alpha_out[0:pg.sid.shape[0], 0] = alpha_imp[0][1:(1+pg.sid.shape[0])]
        alpha_out[0:pg.sid.shape[0], 1] = np.sqrt(np.diag(alpha_imp[1])[1:(1+pg.sid.shape[0])])
        alpha_out[pg.sid.shape[0],0] = alpha_proband[0][1]
        alpha_out[pg.sid.shape[0],1] = np.sqrt(np.diag(alpha_proband[1])[1])
        print('Saving estimates to '+args.out+ '.effects.txt')
        outcols = np.hstack((pg.sid,np.array(['population']))).reshape((pg.sid.shape[0]+1,1))
        np.savetxt(args.out + '.effects.txt',
                   np.hstack((outcols, np.array(alpha_out, dtype='S'))),
                   delimiter='\t', fmt='%s')
        print('Saving sampling covariance matrix of estimates to ' + args.out + '.vcov.txt')
        np.savetxt(args.out + '.vcov.txt', alpha_imp[1][1:(1+pg.sid.shape[0]),1:(1+pg.sid.shape[0])])
