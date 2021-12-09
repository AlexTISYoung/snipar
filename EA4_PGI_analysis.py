####################################################################################################################################################################################################################
# Script for estimating direct, paternal, maternal, and population effect PGIs using a linear mixed model
# Takes as input the PGI file with columns FID, IID, PGI
# The phenotype file with columns FID, IID, phenotype_1, phenotype_2, ..
# The pedigree file for the sample
# The hsq file output by GCTA REML which gives the variance components from a linear mixed model fit in GCTA with the PGI as a fixed effect
# The model has three GRMs: the SNP GRM with off-diagonal elements greater than a threshold (0.025) set to zero (snp_grm_lower)
# the SNP GRM with off-diagonal elements less than a threshold (0.025) set to zero (snp_grm_upper)
# a GRM with 1 on the off-diagonal IFF the pair are full-siblings, and zero otherwise (sib_grm)
# The hsq file should contain the variance component estimates for those three GRMs in that order: snp_grm_lower, snp_grm_upper, sib_grm
# The GRMs are used to construct the residual covariance model 
# The script will compute generalized least-squares estimates from a sibling-difference model if the --sibdiff option is used.
# The script will compute generalized least-squares estimates from a trio model if the --trios option is used
# If both the --sibdiff and --trios options are given, the script will meta-analyse the estimates from both, accounting for covariance between the estimates using the residual covariance model from the GRMs
# The script will estimate the correlation between maternal and paternal PGIs and between sibling PGIs
# The script will output the estimates of: direct, paternal, maternal, population, population-direct, and maternal-paternal
# The population and population-direct estimates are adjusted for the correlation between maternal and paternal and between sibling PGIs due to non-random mating
# The script will also output the variance-covariance matrix of these estimates (effects in _effects.txt, variance-covariance in _vcov.txt)
####################################################################################################################################################################################################################

from sibreg.sibreg import *
import argparse, code
from pysnptools.snpreader import Pheno
import numpy as np

######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('pgs', type=str, help='Location of the PGS file')
    parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
    parser.add_argument('ped',type=str,help='Location of pedigree file with FID giving sibships')
    parser.add_argument('hsq',type=str,help='Location of GCTA output with variance components')
    parser.add_argument('outprefix',type=str,help='Location to output association statistic hdf5 file')
    parser.add_argument('snp_grm_lower',type=str,help='Location of SNP grm for distant relatives')
    parser.add_argument('snp_grm_upper',type=str,help='Location of SNP grm for close relatives')
    parser.add_argument('sib_grm',type=str,help='Location of grm of sibling relations') 
    parser.add_argument('--vcomp',type=str,help='Variance components (SNP h^2, pedigree, sib, residual -- space separated)',default=None)
    parser.add_argument('--scale_phen',action='store_true',help='Scale phenotype to have variance 1',default=False)
    parser.add_argument('--sibdiff',action='store_true',default = False,help='Fit sibling difference in PGS model')
    parser.add_argument('--trios',action = 'store_true',default = False,help='Fit model with individuals with both parents genotyped')
    parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                        default=1)
    parser.add_argument('--tau_init',type=float,help='Initial value for ratio between shared family environmental variance and residual variance',
                        default=1)
    parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)',default='NA')
    args=parser.parse_args()

    if args.trios or args.sibdiff:
        pass
    else:
        raise ValueError('Must do at least one analysis from --trios or --sibdiff')

    # Read PGS
    print('Reading '+str(args.pgs))
    f = open(args.pgs, 'r')
    cols = f.readline()
    if len(cols.split('\t')) > len(cols.split(' ')):
        cols = np.array(cols.split('\t'))
        delim = '\t'
    else:
        cols = np.array(cols.split(' '))
        delim = ' '
    if cols[0] == 'FID' and cols[1]== 'IID':
        pass
    else:
        raise ValueError('First two columns of PGS must be FID, IID')
    f.close()
    ids = np.loadtxt(args.pgs, dtype='U', usecols=(0,1), delimiter=delim, skiprows=1)
    pgs_vals = np.loadtxt(args.pgs, usecols=tuple([x for x in range(2, cols.shape[0])]),delimiter=delim, skiprows=1)
    pg = gtarray(pgs_vals.reshape((pgs_vals.shape[0],1)), ids[:, 1], sid=cols[2:cols.shape[0]], fams=ids[:, 0])
    print('Normalising PGS to have mean zero and variance 1')
    pg.mean_normalise()
    pg.scale()

    # Read phenotype
    print('Reading '+str(args.phenofile))
    pheno = Pheno(args.phenofile, missing=args.missing_char).read()
    # pheno = Pheno('phenotypes/eduyears_resid.ped', missing='NA').read()
    y = np.array(pheno.val)
    pheno_ids = np.array(pheno.iid)[:, 1]
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
    y = y-np.mean(y)
    if args.scale_phen:
        y = y/np.std(y)
    print('Number of non-missing phenotype observations: ' + str(y.shape[0]))

    # Read pedigree
    print('Reading '+str(args.ped))
    ped = np.loadtxt(args.ped,dtype='U')
    controls = np.array([x[0]=='_' for x in ped[:,0]])
    ped = ped[np.logical_not(controls),:]
    ped[:,4] = np.array(['True' for x in range(ped.shape[0])])
    ped[:,5] = np.array(['True' for x in range(ped.shape[0])])
    ped_dict = make_id_dict(ped,1)

    # Construct phenotypic covariance matrix
    print('Constructing phenotypic covariance matrix')
    vcomps = np.zeros((4))
    if args.vcomp is None:
        hsq = open(args.hsq,'r')
        #hsq = open('varcomps/10.hsq', 'r')
        hsq_line = hsq.readline()
        while len(hsq_line)>0:
            hsq_line = hsq_line.split('\t')
            if hsq_line[0]=='V(G1)':
                vcomps[0] = float(hsq_line[1])
            if hsq_line[0]=='V(G2)':
                vcomps[1] = float(hsq_line[1])
            if hsq_line[0]=='V(G3)':
                vcomps[2] = float(hsq_line[1])
            if hsq_line[0]=='V(e)':
                vcomps[3] = float(hsq_line[1])
            hsq_line = hsq.readline()
    else:
        vcomps = np.array([float(x) for x in args.vcomp.split(' ')])

    R_sib_ids = np.loadtxt(args.sib_grm+'.grm.id', dtype=str)[:, 1]
    R_sib_dict = make_id_dict(R_sib_ids)
    R_ids = np.loadtxt(args.snp_grm_lower+'.grm.id', dtype=str)[:, 1]
    R_ids_dict = make_id_dict(R_ids)
    common_ids = np.intersect1d(R_sib_ids,R_ids)
    R = np.zeros((common_ids.shape[0],common_ids.shape[0]),dtype=np.float32)
    R_sib = np.zeros((R_sib_ids.shape[0],R_sib_ids.shape[0]),dtype=np.float32)
    R_sib[np.tril_indices(R_sib.shape[0])] =vcomps[2]*np.fromfile(args.sib_grm+'.grm.bin',dtype=np.float32)
    sib_indices = np.array([R_sib_dict[x] for x in common_ids])
    R += R_sib[np.ix_(sib_indices,sib_indices)]
    del R_sib
    R_snp = np.zeros((R_ids.shape[0],R_ids.shape[0]),dtype=np.float32)
    R_snp[np.tril_indices(R_snp.shape[0])] = vcomps[0]*np.fromfile(args.snp_grm_lower+'.grm.bin',dtype=np.float32)
    R_snp[np.tril_indices(R_snp.shape[0])] += vcomps[1]*np.fromfile(args.snp_grm_upper+'.grm.bin',dtype=np.float32)
    R_indices = np.array([R_ids_dict[x] for x in common_ids])
    R += R_snp[np.ix_(R_indices,R_indices)]
    del R_snp
    R += R.T
    np.fill_diagonal(R,np.diag(R)/2.0+vcomps[3])
    R_dict = make_id_dict(common_ids)


    if args.trios:
        print('Analysing individuals with both parents genotyped')
        par_status, gt_indices, fam_labels = find_par_gts(pg.ids, ped, pg.fams, pg.id_dict)
        # Identify individuals with all observed
        all_obs = np.min(gt_indices, axis=1) > 0
        n = np.sum(all_obs)
        if n == 0:
            raise ValueError('No individuals with both parents genotyped')
        print(str(n) + ' individuals with both parents genotyped')
        gt_indices = gt_indices[all_obs, :]
        bpg_ids = pg.ids[all_obs]
        fam_labels = fam_labels[all_obs]
        # Make array
        G = np.ones((n, 4), dtype=np.float32)
        G[:, 1] = pg.gts[gt_indices[:, 0], 0]
        G[:, 2] = pg.gts[gt_indices[:, 1], 0]
        G[:, 3] = pg.gts[gt_indices[:, 2], 0]
        # Compute correlation between parental PGS
        r_trio = np.cov(G[:,2:4].T)[0,1]
        print('Correlation between maternal and paternal PGS: '+str(r_trio))
        # Estimate
        bpg_id_dict = make_id_dict(bpg_ids)
        in_bpg = np.array([x in bpg_id_dict for x in pheno_ids])
        bpg_indices = np.array([bpg_id_dict[x] for x in pheno_ids[in_bpg]])
        R_trio_indices = np.array([R_dict[x] for x in pheno_ids[in_bpg]])
        R_trio_inv = np.linalg.inv(R[np.ix_(R_trio_indices,R_trio_indices)])
        XTX_trio = np.dot(G[bpg_indices,:].T,R_trio_inv.dot(G[bpg_indices,:]))
        XTY_trio = np.dot(G[bpg_indices,:].T,R_trio_inv.dot(y[in_bpg]))
        alpha_bpg = np.linalg.solve(XTX_trio,XTY_trio)
        alpha_bpg_cov = np.linalg.inv(XTX_trio)

    if args.sibdiff:
        print('Analysing with sibling difference method')
        fam_means, fam_sizes, fam_sums = get_fam_means(pg.ids, ped, pg.gts, pg.ids, remove_proband=False, return_famsizes = True)
        inv_fsize = np.mean(1/np.array(fam_sizes,dtype=np.float64))
        print('Mean of inverse of fam sizes: '+str(inv_fsize))
        # Remove overlap with trios
        if args.trios:
            in_bpg = np.array([x in bpg_id_dict for x in fam_means.ids])
            n_overlap = np.sum(in_bpg)
            if n_overlap == fam_means.ids.shape[0]:
                raise ValueError('No sibships without both parents genotyped')
            else:
                print('Removing '+str(n_overlap)+' individuals with both parents genotyped from sib difference analysis')
                fam_means = gtarray(fam_means.gts[np.logical_not(in_bpg),:],fam_means.ids[np.logical_not(in_bpg)])
        print('Found '+str(fam_means.ids.shape[0])+' individuals with genotyped siblings')
        G_sdiff = np.ones((fam_means.gts.shape[0],3),dtype = np.float32)
        pg_indices = np.array([pg.id_dict[x] for x in fam_means.ids])
        G_sdiff[:,1] = pg.gts[pg_indices,0]
        G_sdiff[:,2] = fam_means.gts[:,0]
        # Compute correlation between siblings PGS
        cov_sib = np.var(fam_sums[:,0]/fam_sizes)
        r_sib = (2*cov_sib-1-inv_fsize)/(1-inv_fsize)
        print('Correlation between siblings PGS: '+str((1+r_sib)/2.0))
        # Compute deviation from family mean
        G_sdiff[:,1] = G_sdiff[:,1] - G_sdiff[:,2]
        fam_labels = np.array([ped[ped_dict[x],0] for x in fam_means.ids])
        # Match with phenotype
        in_fam_means = np.array([x in fam_means.id_dict for x in pheno_ids])
        fam_means_indices = np.array([fam_means.id_dict[x] for x in pheno_ids[in_fam_means]])
        print('Estimating model using sibling differences')
        R_sdiff_indices = np.array([R_dict[x] for x in pheno_ids[in_fam_means]])
        R_sdiff_inv = np.linalg.inv(R[np.ix_(R_sdiff_indices, R_sdiff_indices)])
        XTX_sdiff = np.dot(G_sdiff[fam_means_indices, :].T, R_sdiff_inv.dot(G_sdiff[fam_means_indices, :]))
        XTY_sdiff = np.dot(G_sdiff[fam_means_indices, :].T, R_sdiff_inv.dot(y[in_fam_means]))
        alpha_sdiff = np.linalg.solve(XTX_sdiff, XTY_sdiff)
        alpha_sdiff_cov = np.linalg.inv(XTX_sdiff)

    if args.trios and args.sibdiff:
        print('Meta analysing trio and sib diff estimates')
        sdiff_side = alpha_sdiff_cov.dot(np.dot(G_sdiff[fam_means_indices, :].T,R_sdiff_inv))
        trio_side = alpha_bpg_cov.dot(np.dot(G[bpg_indices,:].T,R_trio_inv))
        sdiff_trio_cov = sdiff_side.dot(R[np.ix_(R_sdiff_indices,R_trio_indices)].dot(trio_side.T))
        ests = np.hstack((alpha_sdiff[1:alpha_sdiff.shape[0]],alpha_bpg[1:alpha_bpg.shape[0]])).reshape((5,1))
        vcov = np.zeros((5,5))
        vcov[0:2,0:2] = alpha_sdiff_cov[1:3,1:3]
        vcov[0:2,2:5] = sdiff_trio_cov[1:3,1:4]
        vcov[2:5,0:2] = vcov[0:2,2:5].T
        vcov[2:5,2:5] = alpha_bpg_cov[1:4,1:4]
        # Compute meta-analysis estimate of r
        r_meta = (G.shape[0]*r_trio+fam_sizes.shape[0]*r_sib)/float(G.shape[0]+fam_sizes.shape[0])
        print('meta analysis estimate of correlation between maternal and paternal PGS: '+str(r_meta))
        sib_c = (1+inv_fsize*(1-r_meta)/(1+r_meta))**(-1)
        A = np.array([[1,0,0],[1,sib_c,sib_c],[1,0,0],[0,1,0],[0,0,1]])
        vcov_inv = np.linalg.inv(vcov)
        AVA = np.dot(A.T,vcov_inv.dot(A))
        meta_est = np.linalg.solve(AVA,np.dot(A.T,vcov_inv.dot(ests)))
        meta_cov = np.linalg.inv(AVA)
        A_out = np.array([[1,0,0],[0,1,0],[0,0,1],[1,(1+r_meta)/2.0,(1+r_meta)/2.0],[0,(1+r_meta)/2.0,(1+r_meta)/2.0],[0,-1,1]])
        alpha_out = A_out.dot(meta_est)
        alpha_out_cov = A_out.dot(meta_cov.dot(A_out.T))
        alpha_out_ses = np.sqrt(np.diag(alpha_out_cov))
        np.savetxt(args.outprefix+'_effects.txt',np.hstack((alpha_out,alpha_out_ses.reshape((6,1)))))
        np.savetxt(args.outprefix+'_vcov.txt',alpha_out_cov[0:3,0:3])
        np.savetxt(args.outprefix+'_r.txt',r_meta.reshape(1,)) 
