from sibreg.sibreg import *
import argparse, code
from pysnptools.snpreader import Pheno

######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('pgs', type=str, help='Location of the PGS file')
    parser.add_argument('phenofile',type=str,help='Location of the phenotype file')
    parser.add_argument('ped',type=str,help='Location of pedigree file with FID giving sibships')
    parser.add_argument('outprefix',type=str,help='Location to output association statistic hdf5 file')
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
    print('Number of non-missing phenotype observations: ' + str(y.shape[0]))

    # Read pedigree
    print('Reading '+str(args.ped))
    ped = np.loadtxt(args.ped,dtype='U')
    controls = np.array([x[0]=='_' for x in ped[:,0]])
    ped = ped[np.logical_not(controls),:]
    ped_dict = make_id_dict(ped,1)

    # Estimate associative effect
    print('Estimating associative effect')
    in_pg = np.array([x in pg.id_dict for x in pheno_ids])
    pg_indices = np.array([pg.id_dict[x] for x in pheno_ids[in_pg]])
    pg_in_ped = np.array([x in ped_dict for x in pg.ids])
    pg.fams[pg_in_ped] = np.array([ped[ped_dict[x],0] for x in pg.ids[pg_in_ped]])
    pg_model = model(y[in_pg], pg.gts[pg_indices,0], pg.fams[pg_indices], add_intercept=True)
    sigma_2_init = np.var(y) * args.tau_init / (1 + args.tau_init)
    pg_optim = pg_model.optimize_model(np.array([sigma_2_init,args.tau_init]))
    pg_alpha = pg_model.alpha_mle(pg_optim['tau'],pg_optim['sigma2'],compute_cov = True)
    sibcor = 1/(1+pg_optim['tau'])
    print('Sibling correlation estimate: '+str(round(sibcor,4)))
    print('Associative effect: '+str(round(pg_alpha[0][1],6))+ ' (S.E. '+str(round(np.sqrt(pg_alpha[1][1,1]),7))+')')

    if args.trios:
        print('Analysing individuals with both parents genotyped')
        par_status, gt_indices, fam_labels = find_par_gts(pg.ids, ped, pg.id_dict)
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
        G = np.zeros((n, 3), dtype=np.float32)
        G[:, 0] = pg.gts[gt_indices[:, 0], 0]
        G[:, 1] = pg.gts[gt_indices[:, 1], 0]
        G[:, 2] = pg.gts[gt_indices[:, 2], 0]
        # Estimate
        bpg_id_dict = make_id_dict(bpg_ids)
        in_bpg = np.array([x in bpg_id_dict for x in pheno_ids])
        bpg_indices = np.array([bpg_id_dict[x] for x in pheno_ids[in_bpg]])
        print('Estimate model for individuals with both parents genotyped')
        bpg_model = model(y[in_bpg], G[bpg_indices, :], fam_labels[bpg_indices], add_intercept=True)
        alpha_bpg = bpg_model.alpha_mle(pg_optim['tau'],pg_optim['sigma2'],compute_cov = True)
        outcols = np.array(['proband', 'paternal', 'maternal']).reshape((3,1))
        # Save output
        alpha_bpg_out = np.zeros((3, 2))
        alpha_bpg_out[:, 0] = alpha_bpg[0][1:4]
        alpha_bpg_out[:, 1] = np.sqrt(np.diag(alpha_bpg[1])[1:4])
        np.savetxt(args.outprefix + '.bpg.pgs_effects.txt',
                   np.hstack((outcols, np.array(alpha_bpg_out, dtype='S20'))),
                   delimiter='\t', fmt='%s')
        np.savetxt(args.outprefix + '.bpg.pgs_vcov.txt', alpha_bpg[1][1:4, 1:4])

    if args.sibdiff:
        print('Analysing with sibling difference method')
        fam_means = get_fam_means(pg.ids, ped, pg.gts, pg.ids, remove_proband=False)
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
        G = np.zeros((fam_means.gts.shape[0],2),dtype = np.float32)
        pg_indices = np.array([pg.id_dict[x] for x in fam_means.ids])
        G[:,0] = pg.gts[pg_indices,0]
        G[:,1] = fam_means.gts[:,0]
        G[:,0] = G[:,0] - G[:,1]
        fam_labels = np.array([ped[ped_dict[x],0] for x in fam_means.ids])
        # Match with phenotype 
        in_fam_means = np.array([x in fam_means.id_dict for x in pheno_ids])
        fam_means_indices = np.array([fam_means.id_dict[x] for x in pheno_ids[in_fam_means]])
        print('Estimating model using sibling differences')
        sdiff_model = model(y[in_fam_means], G[fam_means_indices, :], fam_labels[fam_means_indices], add_intercept=True)
        alpha_sdiff = sdiff_model.alpha_mle(pg_optim['tau'],pg_optim['sigma2'],compute_cov = True)
        alpha_sdiff_out = np.zeros((2, 2))
        alpha_sdiff_out[:, 0] = alpha_sdiff[0][1:3]
        alpha_sdiff_out[:, 1] = np.sqrt(np.diag(alpha_sdiff[1])[1:3])
        np.savetxt(args.outprefix + '.sibdiff.pgs_effects.txt',
                   np.hstack((np.array(['direct', 'between-family']).reshape((2,1)), np.array(alpha_sdiff_out, dtype='S20'))),
                   delimiter='\t', fmt='%s')
