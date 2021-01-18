#!/well/kong/users/wiw765/anaconda3/bin/python
from pysnptools.snpreader import Bed, Pheno
from sibreg.sibreg import *
import h5py, argparse, code
import pandas as pd

def pgs_write(pg,pgs_out,name):
    pg.mean_normalise()
    # Rescale by observed proband PGS
    pg.gts = pg.gts / np.std(pg.gts[:, 0])
    print('PGS computed')
    ####### Write PGS to file ########
    pgs_out[name] = pg.gts
    pgs_out[name+'_ids'] = encode_str_array(pg.ids)
    pgs_out[name+'_cols'] = encode_str_array(pg.sid)
    pgs_out[name+'_fams'] = encode_str_array(pg.fams)
    return None


######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('outprefix',type=str,help='Location to output association statistic hdf5 file')
    parser.add_argument('--gts_list',type=str,help='File with list of bed files of observed genotypes', default = None)
    parser.add_argument('--pargts_list', type=str, help='File with list of imputed parental genotype HDF5 files', default = None)
    parser.add_argument('--weights',type=str,help='Location of the PGS allele weights', default = None)
    parser.add_argument('--phenofile',type=str,help='Location of the phenotype file',default = None)
    parser.add_argument('--pgs', type=str, help='Location of the pre-computed PGS file', default=None)
    parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                        default=1)
    parser.add_argument('--scale_phen',action='store_true',help='Scale the phenotype to have variance 1 in the analysis sample',default=False)
    parser.add_argument('--compute_controls', action='store_true', default=False,
                        help='Compute PGS for control families (default False)')
    parser.add_argument('--fit_sib',action='store_true',default=False,help='Fit indirect effects from siblings')
    parser.add_argument('--ped',type=str,help='Path to pedigree file. By default uses pedigree in imputed parental genotype HDF5 file',default=None)
    parser.add_argument('--tau_init',type=float,help='Initial value for ratio between shared family environmental variance and residual variance',
                        default=1)
    parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)',default=0.01)
    parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)',default='NA')
    parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)',default=5)
    parser.add_argument('--append',action='store_true',default=False,help='Append results to existing output file with given outprefix (default overwrites existing')
    parser.add_argument('--no_covariate_estimates',action='store_true',default=False,help='Suppress output of covariate effect estimates')
    args=parser.parse_args()

    if args.weights is not None:
        if args.gts_list is None:
            raise ValueError('Weights provided but no observed genotypes provided')
        if args.pargts_list is None:
            raise ValueError('Weights provided but no imputed parental genotypes provided')
        print('Computing PGS from weights file')
        ####### Read PGS #######
        weights = pd.read_csv(args.weights, delimiter='\t')
        w_sid = np.array(weights.loc[:, 'sid'], dtype='U')
        print('Read weights for '+str(w_sid.shape[0])+' SNPs')
        beta = np.array(weights.loc[:, 'ldpred_beta'])
        a1 = np.array(weights.loc[:, 'nt1'], dtype='U')
        a2 = np.array(weights.loc[:, 'nt2'], dtype='U')
        p = pgs(w_sid,beta,np.vstack((a1,a2)).T)

        ###### Compute PGS ########
        G_list = []
        gts_list = np.loadtxt(args.gts_list,dtype='U')
        pargts_list = np.loadtxt(args.pargts_list,dtype='U')
        if not gts_list.shape[0] == pargts_list.shape[0]:
            raise ValueError('Lists of imputed and observed genotype files not of same length')
        print('Computing PGS')
        print('Using '+str(pargts_list[0])+' and '+str(gts_list[0]))
        G = get_gts_matrix(pargts_list[0],gts_list[0],p.snp_ids,sib = args.fit_sib)
        pg = compute_pgs(pargts_list[0],gts_list[0],p, sib = args.fit_sib, compute_controls = args.compute_controls)
        for i in range(1,gts_list.shape[0]):
            print('Using ' + str(pargts_list[i]) + ' and ' + str(gts_list[i]))
            if args.compute_controls:
                pg_i = compute_pgs(pargts_list[i],gts_list[i],p, sib = args.fit_sib,compute_controls = args.compute_controls)
                pg = [pg[x].add(pg_i[x]) for x in range(0,len(pg))]
            else:
                pg = pg.add(compute_pgs(pargts_list[i],gts_list[i],p, sib = args.fit_sib,compute_controls = args.compute_controls))
        print('PGS computed')
        ####### Write PGS to file ########
        print('Writing PGS to ' + args.outprefix + '.pgs.hdf5')
        pgs_out = h5py.File(args.outprefix + '.pgs.hdf5', 'w')
        if args.compute_controls:
            pgs_write(pg[0], pgs_out, 'pgs')
            pgs_write(pg[1],pgs_out,'_p_')
            pgs_write(pg[2], pgs_out, '_m_')
            pgs_write(pg[3],pgs_out,'_o_')
        else:
            pgs_write(pg,pgs_out,'pgs')
        pgs_out.close()
    elif args.pgs is not None:
        if args.phenofile is None:
            raise ValueError('Pre-computed PGS provided but no phenotype provided')
        print('Reading PGS from '+args.pgs)
        pgs_f = h5py.File(args.pgs, 'r')
        pg = gtarray(np.array(pgs_f['pgs']),
                     convert_str_array(np.array(pgs_f['ids'])),
                     sid=convert_str_array(np.array(pgs_f['cols'])),
                     fams=convert_str_array(np.array(pgs_f['fams'])))
        print('Normalising PGS')
        pg.mean_normalise()
        pgs_f.close()
    else:
        raise ValueError('Weights or PGS must be provided')

    if args.phenofile is not None:
        print('Fitting PGS for '+str(args.phenofile))
        pheno = Pheno(args.phenofile, missing=args.missing_char).read()
        # pheno = Pheno('phenotypes/eduyears_resid.ped', missing='NA').read()
        y = np.array(pheno.val)
        pheno_ids = np.array(pheno.iid)[:,1]
        if y.ndim == 1:
            pass
        elif y.ndim == 2:
            y = y[:, args.phen_index - 1]
        else:
            raise ValueError('Incorrect dimensions of phenotype array')
        # Remove y NAs
        y_not_nan = np.logical_not(np.isnan(y))
        if np.sum(y_not_nan) < y.shape[0]:
            y = y[y_not_nan]
            pheno_ids = pheno_ids[y_not_nan]
        print('Number of non-missing phenotype observations: ' + str(y.shape[0]))
        in_pgs = np.array([x in pg.id_dict for x in pheno_ids])
        y = y[in_pgs]
        pheno_ids = pheno_ids[in_pgs]
        gt_indices = np.array([pg.id_dict[x] for x in pheno_ids])
        print('Final sample size: '+str(gt_indices.shape[0]))
        # Scale
        if args.scale_phen:
            y = y/np.std(y)
        pg.gts = pg.gts / np.std(pg.gts[gt_indices, 0])
        # Estimate effects
        alpha_imp = get_alpha_mle(y,pg.gts[gt_indices,:] , pg.fams[gt_indices], add_intercept = True)
        # Estimate proband only model
        alpha_proband = get_alpha_mle(y, pg.gts[gt_indices, 0], pg.fams[gt_indices], add_intercept=True)
        # Get print out for fixed mean effects
        alpha_out = np.zeros((pg.sid.shape[0]+1, 2))
        alpha_out[0:pg.sid.shape[0], 0] = alpha_imp[0][1:(1+pg.sid.shape[0])]
        alpha_out[0:pg.sid.shape[0], 1] = np.sqrt(np.diag(alpha_imp[1])[1:(1+pg.sid.shape[0])])
        alpha_out[pg.sid.shape[0],0] = alpha_proband[0][1]
        alpha_out[pg.sid.shape[0],1] = np.sqrt(np.diag(alpha_proband[1])[1])
        print('Saving estimates to '+args.outprefix+ '.pgs_effects.txt')
        outcols = np.hstack((pg.sid,np.array(['associative']))).reshape((pg.sid.shape[0]+1,1))
        np.savetxt(args.outprefix + '.pgs_effects.txt',
                   np.hstack((outcols, np.array(alpha_out, dtype='S20'))),
                   delimiter='\t', fmt='%s')
        print('Saving sampling covariance matrix of estimates to ' + args.outprefix + '.pgs_vcov.txt')
        np.savetxt(args.outprefix + '.pgs_vcov.txt', alpha_imp[1][1:(1+pg.sid.shape[0]),1:(1+pg.sid.shape[0])])
