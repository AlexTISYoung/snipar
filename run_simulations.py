from snipar.simulate import *
import snipar.pgs as pgs
from snipar.gtarray import gtarray
from numba import njit, prange

@njit(parallel=False)
def trio_gwas(y,gts_matrix):
    # Perform least squares regression
    ab = np.zeros((gts_matrix.shape[2],gts_matrix.shape[1]),dtype=np.float_)
    for i in range(gts_matrix.shape[2]):
        ab[i,:] = np.linalg.lstsq(gts_matrix[:,:,i],y)[0]
    return ab

def compute_pgis(a, new_haps, haps, father_indices, mother_indices, ped):
    causal = np.arange(0,a.shape[0])
    # Direct effect component of offspring
    G_males, G_females = compute_genetic_component(new_haps,causal,a)
    # Direct effect component of parents
    G_father, G_mother = compute_genetic_component(haps,causal,a)
    G_father, G_mother = G_father[father_indices], G_mother[mother_indices]
    # pgs array
    pgs_vals = np.zeros((G_males.shape[0]*2,3))
    pgs_vals[np.arange(0,pgs_vals.shape[0],step=2),0] = G_males
    pgs_vals[np.arange(1,pgs_vals.shape[0],step=2),0] = G_females
    pgs_vals[:,1] = np.repeat(G_father,2)
    pgs_vals[:,2] = np.repeat(G_mother,2)
    # pgarray
    pg = pgs.pgarray(pgs_vals,ped[:,1],sid=np.array(['proband','paternal','maternal']),
                    fams=ped[:,0],par_status = np.zeros((pgs_vals.shape[0],2)),
                    ped=ped)
    return pg

def estimate_direct_effect_pgi(pheno, ped, new_haps, haps, father_indices, mother_indices):
    ## Compute genotypes
    gts_matrix = np.zeros((new_haps[0].shape[0]*2,2,new_haps[0].shape[2]))
    proband_genotypes = np.sum(new_haps[0],axis=3)
    gts_matrix[:,0] = proband_genotypes.reshape((proband_genotypes.shape[0]*2,proband_genotypes.shape[2]),order='C')
    gts_matrix[np.arange(0,gts_matrix.shape[0],2),1] = np.sum(haps[0][father_indices,0,:,:],axis=2)
    gts_matrix[np.arange(0,gts_matrix.shape[0],2),1] += np.sum(haps[0][mother_indices,1,:,:],axis=2)
    gts_matrix[np.arange(1,gts_matrix.shape[0],2),1] = gts_matrix[np.arange(0,gts_matrix.shape[0],2),1]
    ## Estimate direct effects
    # Mean normalize
    y = pheno.gts[:,0]
    y = y-np.mean(y)
    for j in range(gts_matrix.shape[1]):
        gts_matrix[:,j,:]=gts_matrix[:,j,:]-np.mean(gts_matrix[:,j,:],axis=0)
    # Least squares regression
    ab = trio_gwas(y,gts_matrix)
    # Compute PGIs
    pg = compute_pgis(ab[:,0], new_haps, haps, father_indices, mother_indices, ped)
    return pg, ab

def pgi_analysis(a, ped, new_haps, haps, father_indices, mother_indices, h2f, h2f_se, outprefix):
    ## Extract phenotype from pedigree
    ped = ped[1:,:] # Remove header
    # Find last generation
    ngen = int(ped[ped.shape[0]-1,0].split('_')[0])
    ped = ped[[x.split('_')[0]==str(ngen) for x in ped[:,0]],:]
    ## Extract relevant pieces of information
    # pedigree
    pedigree = ped[:,0:4]
    # phenotype
    pheno = gtarray(ped[:,5].astype(float).reshape(ped.shape[0],1),ids=ped[:,1])
    pheno.scale()
    # Examine different noise levels
    for vratio in [0, 1, 10, 100]:
        a_est = a+np.random.normal(0,np.sqrt(vratio)*np.std(a),size=a.shape)
        pg = compute_pgis(a_est, new_haps, haps, father_indices, mother_indices, pedigree)
        # Estimate rk
        rk = np.corrcoef(pg.gts[:,2],pg.gts[:,1])[0,1]
        rk_se = (1-rk**2)/np.sqrt(pg.gts.shape[0]-2)
        ## Estimate direct effect and average NTC
        # scale
        pg.scale()
        # Estimate two-generation model
        estimate, estimate_cols = pgs.make_and_fit_model(pheno, pg, pg_cols=pg.sid)
        # Estimate parameters
        estimates, ses = pgs.am_adj_2gen(estimate, estimate_cols, h2f, h2f_se, rk=rk, rk_se=rk_se)
        pgs.write_2gen_adj_ests(estimates, ses, outprefix+'_pgi_v'+str(vratio))
    return pg

def extract_true_params(V):
    # Remove header
    V = np.array(V[1:V.shape[0],:],dtype=float)
    # Number of generations
    ngen = V.shape[0]
    params = {'h2_eq': V[ngen-1,0]/V[ngen-1,1],
              'r':V[ngen-1,2]}
    h2f = params['h2_eq']*(1-params['r'])
    params['h2f'] = h2f
    return params
    
# Fixed parameters
nfam = 30000
n_causal = 1000
maf = 0.5
h2 = 0.5
beta_vert = 0

###### No indirect effects #######
for r_y in [0, 0.25, 0.5, 0.75]:
    # Set random or assortative mating
    if r_y ==0:
        n_am = 0
        n_random = 20
    else:
        n_am = 20
        n_random = 0
    # Generate initial genotype data
    haps, maps, snp_ids, alleles, positions, chroms = simulate_first_gen(nfam, n_causal, maf=maf)
    # Run simulation
    new_haps, haps, father_indices, mother_indices, ibd, ped, a, V = forward_sim(haps, maps, n_random, n_am, True, n_causal, h2, 
                                                                        r_y=r_y)
    # Save output
    outprefix = 'simulations/r_y_{}'.format(r_y)
    print('Saving to '+outprefix)
    np.savetxt(outprefix+'_VCS.txt',V,fmt='%s')
    np.savetxt(outprefix+'.ped',ped,fmt='%s')
    np.savetxt(outprefix+'_causal.txt',a,fmt='%s')
    # Extract true parameters
    params = extract_true_params(V)
    # Perform direct effect PGI analysis
    print('Estimating parameters using direct effect PGIs with different noise levels')
    pg = pgi_analysis(a, ped, new_haps, haps, father_indices, mother_indices, params['h2f'], 0, outprefix)

######## Indirect effects #########
# Varying parameters
v_indir = 0.25
for r_direct_indirect in [0, 0.5, 1]:
    for r_y in [0, 0.25, 0.5, 0.75]:
        # Set random or assortative mating
        if r_y ==0:
            n_am = 0
            n_random = 20
        else:
            n_am = 20
            n_random = 0
        # Generate initial genotype data
        haps, maps, snp_ids, alleles, positions, chroms = simulate_first_gen(nfam, n_causal, maf=maf)
        # Run simulation
        new_haps, haps, father_indices, mother_indices, ibd, ped, a, V = forward_sim(haps, maps, n_random, n_am, True, n_causal, h2, 
                                                                            v_indirect=v_indir, r_direct_indirect=r_direct_indirect, 
                                                                            r_y=r_y, beta_vert=beta_vert)
        # Save output
        outprefix = 'simulations/v_indir_{}_r_dir_indir_{}_r_y_{}'.format(v_indir, r_direct_indirect, r_y)
        print('Saving to '+outprefix)
        np.savetxt(outprefix+'_VCS.txt',V,fmt='%s')
        np.savetxt(outprefix+'.ped',ped,fmt='%s')
        np.savetxt(outprefix+'_causal.txt',a,fmt='%s')
        # Extract true parameters
        params = extract_true_params(V)
        # Perform direct effect PGI analysis
        print('Estimating parameters using direct effect PGIs with different noise levels')
        pg = pgi_analysis(a[:,0], ped, new_haps, haps, father_indices, mother_indices, params['h2f'], 0, outprefix)
        # Population effect PGI
        pg = pgi_analysis(a[:,0]+a[:,1], ped, new_haps, haps, father_indices, mother_indices, params['h2f'], 0, outprefix+'_population')
        