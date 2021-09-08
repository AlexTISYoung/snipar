from sibreg.sibreg import *
import argparse
from os import path

def pgs_write(pg,filename,scale_PGS = False):
    if scale_PGS:
        # Rescale by observed proband PGS
        pg.gts = pg.gts / np.std(pg.gts[:, 0])
    ####### Write PGS to file ########
    pg_out = np.column_stack((pg.fams,pg.ids,pg.gts))
    pg_header = np.column_stack((np.array(['FID','IID']).reshape(1,2),pg.sid.reshape(1,pg.sid.shape[0])))
    pg_out = np.row_stack((pg_header,pg_out))
    print('Writing PGS to ' + filename)
    np.savetxt(filename, pg_out, fmt='%s')
    return None

def parse_filelist(obsfiles, impfiles, obsformat):
    obs_files = []
    imp_files = []
    if '~' in obsfiles and impfiles:
        bed_ixes = obsfiles.split('~')
        imp_ixes = impfiles.split('~')
        for i in range(1,23):
            obsfile = bed_ixes[0]+str(i)+bed_ixes[1]+'.'+obsformat
            impfile = imp_ixes[0]+str(i)+imp_ixes[1]+'.hdf5'
            if path.exists(impfile) and path.exists(obsfile):
                obs_files.append(obsfile)
                imp_files.append(impfile)
        print(str(len(imp_files))+' matched observed and imputed genotype files found')
    else:
            obs_files = [obsfiles+'.'+obsformat]
            imp_files = [impfiles+'.hdf5']
    return np.array(obs_files), np.array(imp_files)


######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('outprefix',type=str,help='Prefix for computed PGS file and/or regression results files')
    parser.add_argument('--bedfiles',type=str,help='Address of observed genotype files in .bed format (without .bed suffix). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of 1-22.', default = None)
    parser.add_argument('--bgenfiles',type=str,help='Address of observed genotype files in .bgen format (without .bgen suffix). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of 1-22.', default = None)
    parser.add_argument('--impfiles', type=str, help='Address of hdf5 files with imputed parental genotypes (without .hdf5 suffix). If there is a ~ in the address, ~ is replaced by the chromosome numbers in the range of 1-22.', default = None)
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
        if args.bedfiles is None and args.bgenfiles is None:
            raise ValueError('Weights provided but no observed genotypes provided')
        if args.bedfiles is not None and args.bgenfiles is not None:
            raise ValueError('Provide only one of --bedfiles and --bgenfiles')
        if args.impfiles is None:
            raise ValueError('Weights provided but no imputed parental genotypes provided')
        print('Computing PGS from weights file')
        ####### Read PGS #######
        if args.sep is None:
            weights = np.loadtxt(args.weights,dtype=str)
        else:
            weights = np.loadtxt(args.weight,dtype=args.sep)
        colnames = weights[0,:]
        weights = weights[1:weights.shape[0],:]
        print('Read weights for '+str(weights.shape[0])+' variants')
        beta = np.array(weights[:,np.where(colnames==args.beta_col)[0][0]],dtype=np.float64)
        allele_indices = np.array([np.where(colnames==args.A1)[0][0],np.where(colnames==args.A2)[0][0]])
        p = pgs(weights[:,np.where(colnames==args.SNP)[0][0]],
                beta,
                weights[:,allele_indices])

        ###### Compute PGS ########
        G_list = []
        if args.bedfiles is not None:
            gts_list, pargts_list = parse_filelist(args.bedfiles, args.impfiles, 'bed')
        elif args.bgenfiles is not None:
            gts_list, pargts_list = parse_filelist(args.bgenfiles, args.impfiles, 'bgen')
        if gts_list.shape[0]==0:
            raise(ValueError('No input genotype files found'))
        if not gts_list.shape[0] == pargts_list.shape[0]:
            raise ValueError('Lists of imputed and observed genotype files not of same length')
        print('Computing PGS')
        print('Using '+str(pargts_list[0])+' and '+str(gts_list[0]))
        pg = compute_pgs(pargts_list[0],gts_list[0], p, sib=args.fit_sib, compute_controls=args.compute_controls)
        for i in range(1,gts_list.shape[0]):
            print('Using ' + str(pargts_list[i]) + ' and ' + str(gts_list[i]))
            if args.compute_controls:
                pg_i = compute_pgs(pargts_list[i],gts_list[i], p, sib=args.fit_sib, compute_controls=args.compute_controls)
                pg = [pg[x].add(pg_i[x]) for x in range(0,len(pg))]
            else:
                pg = pg.add(compute_pgs(pargts_list[i],gts_list[i], p, sib=args.fit_sib, compute_controls=args.compute_controls))
        print('PGS computed')
        ####### Write PGS to file ########
        if args.compute_controls:
            pgs_write(pg[0], args.outprefix + '.pgs.txt', scale_PGS=args.scale_pgs)
            pgs_write(pg[1],args.outprefix + '.pgs.control_paternal.txt', scale_PGS=args.scale_pgs)
            pgs_write(pg[2], args.outprefix + '.pgs.control_maternal.txt', scale_PGS=args.scale_pgs)
            pgs_write(pg[3],args.outprefix + '.pgs.control_sibling.txt', scale_PGS=args.scale_pgs)
        else:
            pgs_write(pg, args.outprefix + '.pgs.txt', scale_PGS=args.scale_pgs)
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
        y, pheno_ids = read_phenotype(args.phenofile, missing_char=args.missing_char, phen_index=args.phen_index)
        print('Number of non-missing phenotype observations: ' + str(y.shape[0]))
        # Remove individuals without phenotype observations from PGS
        pg.filter_ids(pheno_ids)
        # Match phenotype to PGS
        y = match_phenotype(pg, y, pheno_ids)
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
            y = y/np.std(y)
        if args.scale_pgs:
            pg.gts = pg.gts / np.std(pg.gts[:, 0])
        # Estimate effects
        print('Estimating direct and indirect/parental effects')
        alpha_imp = fit_sibreg_model(y, pg.gts, pg.fams, add_intercept=True, return_model=False, return_vcomps=False)
        # Estimate population effect
        print('Estimating population effect')
        alpha_proband = fit_sibreg_model(y, pg.gts[:, 0], pg.fams, add_intercept=True, return_model=False, return_vcomps=False)
        # Get print out for fixed mean effects
        alpha_out = np.zeros((pg.sid.shape[0]+1, 2))
        alpha_out[0:pg.sid.shape[0], 0] = alpha_imp[0][1:(1+pg.sid.shape[0])]
        alpha_out[0:pg.sid.shape[0], 1] = np.sqrt(np.diag(alpha_imp[1])[1:(1+pg.sid.shape[0])])
        alpha_out[pg.sid.shape[0],0] = alpha_proband[0][1]
        alpha_out[pg.sid.shape[0],1] = np.sqrt(np.diag(alpha_proband[1])[1])
        print('Saving estimates to '+args.outprefix+ '.pgs_effects.txt')
        outcols = np.hstack((pg.sid,np.array(['population']))).reshape((pg.sid.shape[0]+1,1))
        np.savetxt(args.outprefix + '.pgs_effects.txt',
                   np.hstack((outcols, np.array(alpha_out, dtype='S'))),
                   delimiter='\t', fmt='%s')
        print('Saving sampling covariance matrix of estimates to ' + args.outprefix + '.pgs_vcov.txt')
        np.savetxt(args.outprefix + '.pgs_vcov.txt', alpha_imp[1][1:(1+pg.sid.shape[0]),1:(1+pg.sid.shape[0])])
