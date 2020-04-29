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
    args=parser.parse_args()

    if args.trios or args.sibdiff:
        pass
    else:
        raise ValueError('Must do at least one analysis from --trios or --sibdiff')

    # Read PGS
    print('Reading '+str(args.pgs))
    f = open(args.pgs, 'r')
    cols = f.readline()
    if cols.split('\t').shape[0] > cols.split(' ').shape[0]:
        cols = cols.split('\t')
    else:
        cols = cols.split(' ')
    if cols[0:2] == np.array(['FID', 'IID']):
        pass
    else:
        raise ValueError('First two columns of PGS must be FID, IID')
    f.close()
    ids = np.loadtxt(args.pgs, dtype='U', usecols=(1,), skiprows=1)
    pgs_vals = np.loadtxt(args.pgs, usecols=tuple([x for x in range(2, cols.shape[0] - 1)]))
    pg = gtarray(pgs_vals, ids[:, 1], sid=cols[2:cols.shape[0]], fams=ids[:, 0])

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
    ped_dict = make_id_dict(ped,1)

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
        G = G / np.std(G, axis=0)
        # Estimate
        bpg_id_dict = make_id_dict(bpg_ids)
        in_bpg = np.array([x in bpg_id_dict for x in pheno_ids])
        bpg_indices = np.array([bpg_id_dict[x] for x in pheno_ids[in_bpg]])
        print('Estimate model for individuals with both parents genotyped')
        alpha_bpg = get_alpha_mle(y[in_bpg], G[bpg_indices, :], fam_labels[bpg_indices], add_intercept=True)
        outcols = np.array(['proband', 'paternal', 'maternal'])
        # Save output
        alpha_bpg_out = np.zeros((3, 2))
        alpha_bpg_out[:, 0] = alpha_bpg[0][1:4]
        alpha_bpg_out[:, 1] = np.sqrt(np.diag(alpha_bpg[1])[1:4])
        np.savetxt(args.outprefix + '.bpg.pgs_effects.txt',
                   np.hstack((outcols, np.array(alpha_bpg_out, dtype='S20'))),
                   delimiter='\t', fmt='%s')
        np.savetxt(args.outprefix + '.bpg.pgs_vcov.txt', alpha_bpg_out[1][1:4, 1:4])

    if args.sibdiff:
        print('Analysing with sibling difference method')
        ped = np.zeros((pg.fams.shape[0], 2), dtype=pg.fams.dtype)
        ped[:, 0] = pg.fams
        ped[:, 1] = pg.ids
        G_sdiff = np.zeros((pg.gts.shape[0], 2))
        G_sdiff[:, 0] = pg.gts[:, 0]
        G_sdiff[:, 1] = get_fam_means(pg.ids, ped, pg.gts[:, 0], pg.ids, remove_proband=False)
        not_bpg = np.ones((G_sdiff.shape[0]), dtype=bool)
        if args.trios:
            print('Checking for overlap with individuals with both parents genotyped')
            not_bpg = np.array([x not in set(bpg_ids) for x in pg.ids])
            if np.sum(not_bpg) == 0:
                raise ValueError('No individuals')

        alpha_sdiff = get_alpha_mle(y, G_sdiff[gt_indices, :], pg.fams[gt_indices], add_intercept=True)
        alpha_sdiff_out = np.zeros((2, 2))
        alpha_sdiff_out[:, 0] = alpha_sdiff[0][1:3]
        alpha_sdiff_out[:, 1] = np.sqrt(np.diag(alpha_sdiff[1])[1:3])
        np.savetxt(args.outprefix + '.pgs_effects.txt',
                   np.hstack((np.array(['direct', 'between-family']), np.array(alpha_sdiff_out, dtype='S20'))),
                   delimiter='\t', fmt='%s')