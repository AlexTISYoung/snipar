import numpy as np
from numba import njit, prange
from snipar.utilities import *
from snipar.map import decode_map_from_pos
from bgen_reader import open_bgen

### Read haplotypes from bgen
def haps_from_bgen(bgenfiles,chr_range):
    bgenfiles, chroms = parse_obsfiles(bgenfiles, obsformat='bgen', chromosomes=chr_range)
    # Read genotypes
    haps = []
    maps = []
    snp_ids = []
    alleles = []
    positions = []
    for i in range(bgenfiles.shape[0]):
        print('Reading in chromosome '+str(chroms[i]))
        # Read bgen
        bgen = open_bgen(bgenfiles[i], verbose=True)
        # Read map
        map = decode_map_from_pos(chroms[i], bgen.positions)
        not_nan = np.logical_not(np.isnan(map))
        maps.append(map[not_nan])
        positions.append(bgen.positions[not_nan])
        # Snp
        snp_ids.append(bgen.rsids[not_nan])
        # Alleles
        alleles.append(np.array([x.split(',') for x in bgen.allele_ids[not_nan]]))
        # Read genotypes
        gts = bgen.read(([x for x in range(bgen.samples.shape[0])], [x for x in range(bgen.ids.shape[0])]), np.bool_)[:, :,
            np.array([0, 2])]
        gts = gts[:,not_nan]
        nfam = int(np.floor(gts.shape[0] / 2))
        ngen = np.zeros((nfam, 2, gts.shape[1], 2), dtype=np.bool_)
        ngen[:, 0, :, :] = gts[0:nfam, :, :]
        ngen[:, 1, :, :] = gts[nfam:(2 * nfam), :, :]
        del gts
        haps.append(ngen)
    return haps, maps, snp_ids, alleles, positions, chroms

### Simulate haplotypes ###
# Simulate frequencies: either all one value, or sampled from density proportional to 1/x
def simulate_freqs(nsnp, maf=None, min_maf=None):
    if maf is not None:
        if 0.5>=maf>0:
            freqs = maf*np.ones((nsnp))
        else:
            raise(ValueError('MAF must be between zero and 0.5'))
    else:
        if 0.5>=min_maf>0:
            C = -1/np.log(2*min_maf)
            u = np.random.uniform(0,1,nsnp)
            freqs = min_maf*np.exp(u/C)
        else:
            raise(ValueError('minimum MAF must be between zero and 0.5'))
    return freqs

@njit
def simulate_haplotypes(freqs,nfam):
    haps = np.zeros((nfam,2,freqs.shape[0],2),dtype=np.bool_)
    for i in prange(freqs.shape[0]):
        haps[:, :, i, :] = np.random.binomial(1,freqs[i],(nfam,2,2))
    return haps

def simulate_first_gen(nfam, ncausal, maf=None,min_maf=0.05):
        # Chromosomes
        chroms = [1]
        # Simulate allele frequencies
        freqs = simulate_freqs(ncausal, maf=maf, min_maf=min_maf)
        # Simulate haplotypes 
        haps = [simulate_haplotypes(freqs,nfam)]
        # SNP IDs
        snp_ids = [np.array(np.arange(1,ncausal+1),dtype=str)]
        # Alleles
        alleles = np.zeros((ncausal,2),dtype='S10')
        alleles[:, 0] = 'A1'
        alleles[:, 1] = 'A2'
        alleles = [alleles]
        # Positions
        positions = [np.arange(1,ncausal+1)]
        # Maps
        maps = [np.arange(1,ncausal+1)]
        return haps, maps, snp_ids, alleles, positions, chroms

# Simulate frequencies: either all one value, or sampled from density proportional to 1/x
def simulate_freqs(nsnp, maf=None, min_maf=None, Fst=0):
    if maf is not None:
        if 0.5>=maf>0:
            freqs = maf*np.ones((nsnp))
        else:
            raise(ValueError('MAF must be between zero and 0.5'))
    else:
        if 0.5>=min_maf>0:
            C = -1/np.log(2*min_maf)
            u = np.random.uniform(0,1,nsnp)
            freqs = min_maf*np.exp(u/C)
        else:
            raise(ValueError('minimum MAF must be between zero and 0.5'))
    return freqs

@njit
def simulate_haplotypes(freqs,nfam):
    haps = np.zeros((nfam,2,freqs.shape[0],2),dtype=np.bool_)
    for i in prange(freqs.shape[0]):
        haps[:, :, i, :] = np.random.binomial(1,freqs[i],(nfam,2,2))
    return haps

@njit
def simulate_recombinations(map):
    map_start = map[0]
    map_end = map[map.shape[0]-1]
    map_length = map_end-map_start
    n_recomb = np.random.poisson(map_length/100)
    recomb_points = map_start+np.sort(np.random.uniform(0,1,n_recomb))*map_length
    return n_recomb,recomb_points

######### Simulate meiosis with recombination #########
@njit
def meiosis(map,n=1):
    # Recomb vector
    recomb_vector = np.zeros((n,map.shape[0]), dtype=np.bool_)
    # Do recombinations
    for r in range(n):
        n_recomb, recomb_points = simulate_recombinations(map)
        # Initialize
        last_hap = np.bool_(np.random.binomial(1,0.5,1)[0])
        recomb_vector[r,:] = np.bool_(last_hap)
        # Recombine
        if n_recomb>0:
            for i in range(n_recomb):
                first_snp = np.argmax(map>recomb_points[i])
                recomb_vector[r,first_snp:recomb_vector.shape[1]] = ~recomb_vector[r,first_snp:recomb_vector.shape[1]]
    # Return
    return recomb_vector


##### Produce next generation for unlinked variants #####
@njit(parallel=True)
def produce_next_gen_unlinked(father_indices,mother_indices,males,females):
    ngen = np.zeros((father_indices.shape[0],2,males.shape[1],2),dtype=np.bool_)
    ibd = np.zeros((father_indices.shape[0],males.shape[1],2),dtype=np.bool_)
    for i in prange(ngen.shape[0]):
        meiosis = np.random.binomial(1,0.5,(4,ngen.shape[2]))
        for j in range(0,meiosis.shape[1]):
            # Paternal sib 1
            if meiosis[0,j] == 0:
                ngen[i, 0, j, 0] = males[father_indices[i], j, 0]
            else:
                ngen[i, 0, j, 0] = males[father_indices[i], j , 1]
            # Paternal sib 2
            if meiosis[1,j] == 0:
                ngen[i, 1, j, 0] = males[father_indices[i], j, 0]
            else:
                ngen[i, 1, j, 0] = males[father_indices[i], j , 1]
            # Maternal sib 1
            if meiosis[2,j] == 0:
                ngen[i, 0, j, 1] = females[mother_indices[i], j, 0]
            else:
                ngen[i, 0, j, 1] = females[mother_indices[i], j , 1]
            # Maternal sib 2
            if meiosis[3,j] == 0:
                ngen[i, 1, j, 1] = females[mother_indices[i], j, 0]
            else:
                ngen[i, 1, j, 1] = females[mother_indices[i], j , 1]
            # IBD
            ibd[i, j, 0] = meiosis[0,j]==meiosis[1,j]
            ibd[i, j, 1] = meiosis[2,j]==meiosis[3,j]
    return ngen, ibd

#### Produce next generation with linkage through recombination map #####
@njit(parallel=True)
def produce_next_gen(father_indices,mother_indices,males,females,map):
    ngen = np.zeros((father_indices.shape[0],2,males.shape[1],2),dtype=np.bool_)
    ibd = np.zeros((father_indices.shape[0],males.shape[1],2),dtype=np.bool_)
    for i in prange(father_indices.shape[0]):
        # recombinations
        recomb_i = meiosis(map,n=4)
        ## sib haplotypes
        # sib haplotypes and ibd
        for j in range(ibd.shape[1]):
            # paternal sib 1
            if recomb_i[0,j]:
                ngen[i, 0, j, 0] =  males[father_indices[i], j, 0]
            else:
                ngen[i, 0, j, 0] = males[father_indices[i], j, 1]
            # paternal sib 2
            if recomb_i[1,j]:
                ngen[i, 1, j, 0] =  males[father_indices[i], j, 0]
            else:
                ngen[i, 1, j, 0] = males[father_indices[i], j, 1]
            # maternal sib 1
            if recomb_i[2,j]:
                ngen[i, 0, j, 1] =  females[mother_indices[i], j, 0]
            else:
                ngen[i, 0, j, 1] = females[mother_indices[i], j, 1]
            # maternal sib 2
            if recomb_i[3,j]:
                ngen[i, 1, j, 1] =  females[mother_indices[i], j, 0]
            else:
                ngen[i, 1, j, 1] = females[mother_indices[i], j, 1]
            ibd[i,j,0] = recomb_i[0,j]==recomb_i[1,j]
            ibd[i,j,1] = recomb_i[2,j]==recomb_i[3,j]
    return ngen, ibd

@njit
def random_mating_indices(nfam):
    return np.random.choice(np.array([x for x in range(nfam)],dtype=np.int_),size=nfam,replace=False)

def am_indices(yp,ym,r):
    v = np.sqrt(np.var(yp)*np.var(ym))
    s2 = (1/r-1)*v
    zp = yp+np.sqrt(s2)*np.random.randn(yp.shape[0])
    zm = ym+np.sqrt(s2)*np.random.randn(ym.shape[0])
    rank_p = np.argsort(zp)
    rank_m = np.argsort(zm)
    return rank_p, rank_m

def compute_genetic_component(haps,causal,a):
    causal = set(causal)
    G_m = np.zeros((haps[0].shape[0]))
    G_p = np.zeros((haps[0].shape[0]))
    snp_count = 0
    for chr in range(len(haps)):
        snp_index = np.array([snp_count+x for x in range(haps[chr].shape[2])])
        in_causal = np.array([snp_index[x] in causal for x in range(haps[chr].shape[2])])
        causal_gts = np.sum(haps[chr][:,:,in_causal,:],axis=3)
        G_p += causal_gts[:, 0, :].dot(a[snp_index[in_causal]])
        G_m += causal_gts[:, 1, :].dot(a[snp_index[in_causal]])
        snp_count += haps[chr].shape[2]
    return G_p, G_m

def compute_phenotype(haps,causal,a,sigma2):
    G_p, G_m = compute_genetic_component(haps,causal,a)
    Y_p = G_p+np.sqrt(sigma2)*np.random.randn(G_p.shape[0])
    Y_m = G_m+np.sqrt(sigma2)*np.random.randn(G_m.shape[0])
    return G_p, G_m, Y_p, Y_m

def compute_phenotype_vert(haps,causal,a,sigma2,beta_vert,Y_p,Y_m):
    G_males, G_females = compute_genetic_component(haps, causal, a)
    Y_males = G_males + beta_vert*(Y_p+Y_m)+np.sqrt(sigma2) * np.random.randn(G_males.shape[0])
    Y_females = G_females + beta_vert*(Y_p+Y_m)+np.sqrt(sigma2) * np.random.randn(G_females.shape[0])
    return G_males, G_females, Y_males, Y_females

def compute_phenotype_indirect(haps,old_haps,father_indices,mother_indices,causal,a,b,sigma2):
    if father_indices is None or mother_indices is None:
        raise(ValueError('Must provide parental indices to compute indirect effect components'))
    # Direct effect component of offspring
    G_males, G_females = compute_genetic_component(haps,causal,a)
    # Direct effect component of parents
    G_father, G_mother = compute_genetic_component(old_haps,causal,a)
    # Indirect effect component of parents
    G_p, G_m = compute_genetic_component(old_haps, causal, b)
    # Male offspring phenotype
    Y_males = G_males+G_p[father_indices]+G_m[mother_indices]+np.sqrt(sigma2)*np.random.randn(G_males.shape[0])
    # Female offspring phenotype
    Y_females = G_females+G_p[father_indices]+G_m[mother_indices]+np.sqrt(sigma2)*np.random.randn(G_males.shape[0])
    # Direct males, direct females, phenotype males, phenotype females, indirect effect component father, 
    # indirect effect component mother, direct effect component father, direct effect component mother
    return G_males, G_females, Y_males, Y_females, G_p[father_indices], G_m[mother_indices], G_father[father_indices], G_mother[mother_indices]


def simulate_effects(haps, n_causal, h2, old_haps=None, v_indirect=0, r_direct_indirect=0, father_indices=None, mother_indices=None):
    # Select causal variants
    nsnp_chr = np.array([x.shape[2] for x in haps])
    nsnp = np.sum(nsnp_chr)
    if n_causal > nsnp:
        raise(ValueError('Not enough SNPs to simulate phenotype with '+str(n_causal)+' causal SNPs'))
    causal = np.sort(np.random.choice(np.arange(0, nsnp), n_causal, replace=False))
    # Check if indirect effect parameters are possible
    if v_indirect<0:
            raise(ValueError('Variance due to indirect effects must be non-negative'))
    elif v_indirect>0:
        h2_total = (1+v_indirect+r_direct_indirect*np.sqrt(2*v_indirect))*h2
        print('Direct and indirect variance components specified:')
        print('Random-mating heritability: '+str(round(h2,4)))
        print('Fraction of phenotypic variance explained by direct and indirect genetic effects in a random-mating population: '+str(round(h2_total,4)))
        if h2_total>1:
            raise(ValueError('Parameters imply total variance explained by direct and indirect effects is greater than phenotypic variance. Reduce variance due to indirect effects'))
        else:
            v_eg = v_indirect*h2
            c_ge = np.sqrt(2*v_indirect)*r_direct_indirect*h2
            print('Variance explained by indirect genetic effects: '+str(round(v_eg,4)))
            print('Variance explained by covariance between direct and indirect genetic effects: '+str(round(c_ge,4)))
    # Simulate direct effects
    if v_indirect==0:
        print('Simulating direct genetic effects')
        a = np.zeros((nsnp))
        causal = np.sort(np.random.choice(np.arange(0,nsnp),n_causal,replace=False))
        a[causal] = np.random.randn(n_causal)
        G_p, G_m = compute_genetic_component(haps,causal,a)
        scale_fctr = np.sqrt(h2/np.var(np.hstack((G_p,G_m))))
        a = a*scale_fctr
        h2_total = h2
    else:
        if old_haps is None:
            raise(ValueError('Must provide parental haplotypes to simulate indirect genetic effects'))
        if father_indices is None or mother_indices is None:
            raise(ValueError('Must provide parental indices to simulate indirect genetic effects'))
        print('Simulating indirect genetic effects')
        a = np.zeros((nsnp,2))
        a[causal,:] = np.random.multivariate_normal(np.zeros((2)),
                                                np.array([[1,r_direct_indirect],[r_direct_indirect,1]]),
                                                size=n_causal)
        G_males, G_females, Y_males, Y_females, eta_p, eta_m, delta_p, delta_m = compute_phenotype_indirect(haps,old_haps,father_indices,mother_indices,causal,a[:,0],a[:,1],0)
        # Rescale direct effects
        a[:,0] = a[:,0] * np.sqrt(h2 / np.var(np.hstack((G_males, G_females))))
        # Rescale indirect effects
        a[:,1] = a[:,1] * np.sqrt(v_eg / np.var(eta_p+eta_m))
    return a, causal, h2_total

def compute_vcomps(G_p ,G_m, Y_p, Y_m, delta_p=None, delta_m=None, v_indirect=0, eta_p=None, eta_m=None, verbose=True):
    if v_indirect==0:
        V = np.zeros((3))
    else:
        if delta_p is None or delta_m is None or eta_p is None or eta_m is None:
            raise(ValueError('Missing parental direct/indirect components. Cannot compute all variance components'))
        V = np.zeros((8))
    V[:] = np.nan
    # Print variances
    V[0:2] = np.array([np.var(np.hstack((G_p, G_m))), np.var(np.hstack((Y_p, Y_m)))])
    # Correlation between parents' direct effect components
    if delta_p is not None and delta_m is not None:
        V[2] = np.corrcoef(delta_p,delta_m)[0,1]
    if verbose:
        print('Genetic variance: ' + str(round(V[0], 4)))
        print('Phenotypic variance: ' + str(round(V[1], 4)))
        print('Heritability: ' + str(round(V[0] / V[1], 4)))
        if delta_p is not None and delta_m is not None:
            print('Parental genotypic correlation: '+str(round(V[2],4)))
    # Compute for indirect effect components
    if v_indirect>0:
        # Combined parental indirect effect component
        eta = eta_p+eta_m
        # Variance of parental indirect effect component
        V[3] = np.var(eta)
        # Variance due to covariance between direct and parental indirect effect components
        V[4] = np.cov(G_p,eta)[0,1]+np.cov(G_m,eta)[0,1]
        # Correlation between parents' indirect effect components
        V[5] = np.corrcoef(eta_p,eta_m)[0,1]
        # Correlation between direct and indirect effect components of parents
        V[6] = np.corrcoef(np.hstack((delta_p,delta_m)),np.hstack((eta_p,eta_m)))[0,1]
        # Correlation between direct component of father and indirect component of mother
        V[7] = np.corrcoef(np.hstack((delta_p,delta_m)),np.hstack((eta_m,eta_p)))[0,1]
        if verbose:
            print('Variance due to parental indirect genetic effects: ' + str(round(V[3], 4)))
            print('Variance due to covariance between direct effects and parental indirect genetic effects: ' + str(round(V[4], 4)))
            print('Correlation between maternal and paternal indirect effect components: '+str(round(V[5],4)))
            print('Correlation between direct and indirect effect components within parents: '+str(round(V[6],4)))
            print('Correlation between direct and indirect effect components across parents: '+str(round(V[7],4)))
    return V

def create_ped_output(G_males, G_females, Y_males, Y_females, Y_p, Y_m, delta_p, delta_m, 
                        eta_males=None, eta_females=None, eta_p=None, eta_m=None, father_indices=None, mother_indices=None,gen=None,
                        header=True):
    nfam = Y_males.shape[0]
    # Pedigree columns
    pedcols = ['FID','IID','FATHER_ID','MOTHER_ID','SEX',
                        'PHENO','FATHER_PHENO','MOTHER_PHENO',
                        'DIRECT','FATHER_DIRECT','MOTHER_DIRECT']
    if eta_males is not None and eta_females is not None:
        pedcols.append('INDIRECT')
    if eta_p is not None:
        pedcols.append('FATHER_INDIRECT')
    if eta_m is not None:
        pedcols.append('MOTHER_INDIRECT')
    pedcols = np.array(pedcols)
    ped = np.zeros((nfam*2,pedcols.shape[0]),dtype='U30')
    ## Fill in pedigree array
    # Record generation
    if gen is None:
        gen=''
        prev_gen = ''
    else:
        prev_gen = str(int(gen)-1)+'_'
        gen=str(gen)+'_'
    # FID
    ped[:, 0] = np.repeat(np.array([gen+str(x) for x in range(nfam)]),2)
    # IID
    id1 = np.array([gen+str(x)+'_0' for x in range(nfam)])
    id2 = np.array([gen+str(x)+'_1' for x in range(nfam)])
    ped[0::2, 1] = id1
    ped[1::2, 1] = id2
    # Father ID
    if father_indices is None:
        ped[:, 2] = np.repeat(np.array(['P_'+str(x) for x in range(nfam)]),2)
    else:
        ped[:, 2] = np.repeat(np.array([prev_gen+str(x)+'_0' for x in father_indices]),2)
    # Mother ID
    if mother_indices is None:
        ped[:, 3] = np.repeat(np.array(['M_'+str(x) for x in range(nfam)]),2)
    else:
        ped[:, 3] = np.repeat(np.array([prev_gen+str(x)+'_1' for x in mother_indices]),2)
    # Sex
    ped[0::2, 4] = '0'
    ped[1::2, 4] = '1'
    # Phenotype
    ped[0::2, 5] = Y_males
    ped[1::2, 5] = Y_females
    # Father phenotype
    ped[:, 6] = np.repeat(Y_p, 2)
    # Mother phenotype
    ped[:, 7] = np.repeat(Y_m, 2)
    # Genetic component
    ped[0::2, 8] = G_males
    ped[1::2, 8] = G_females
    # Paternal genetic component
    ped[:, 9] = np.repeat(delta_p, 2)
    # Maternal genetic component
    ped[:, 10] = np.repeat(delta_m, 2)
    # Indirect component
    if eta_males is not None and eta_females is not None:
        ped[0::2, 11] = eta_males
        ped[1::2, 11] = eta_females
    # Paternal indirect component
    if eta_p is not None:
        ped[:, 12] = np.repeat(eta_p, 2)
    if eta_m is not None:
        ped[:, 13] = np.repeat(eta_m, 2)
    # Add header
    if header:
        ped = np.vstack((pedcols,ped))
    return ped

def forward_sim(haps, maps, ngen_random, ngen_am, unlinked, n_causal, h2, v_indirect=0, r_direct_indirect=0, beta_vert=0, r_y=None,):
    nfam = haps[0].shape[0]
    if ngen_am>0 and r_y is None:
        raise(ValueError('Must provide parental phenotypic correlation for assortative mating'))
    ### Simulate population ###
    ### Produce first generation by random-mating
    print('Generating first generation by random-mating')
    father_indices = random_mating_indices(nfam)
    mother_indices = random_mating_indices(nfam)
    # Generate haplotypes of new generation
    new_haps = []
    ibd = []
    for chr in range(0,len(haps)):
        if unlinked:
            new_haps_chr, ibd_chr = produce_next_gen_unlinked(father_indices,mother_indices,haps[chr][:,0,:,:],haps[chr][:,1,:,:])
        else:
            new_haps_chr, ibd_chr = produce_next_gen(father_indices,mother_indices,haps[chr][:,0,:,:],haps[chr][:,1,:,:],maps[chr])
        new_haps.append(new_haps_chr)
        ibd.append(ibd_chr)
    # Pedigree columns
    pedcols = ['FID','IID','FATHER_ID','MOTHER_ID','SEX',
                        'PHENO','FATHER_PHENO','MOTHER_PHENO',
                        'DIRECT','FATHER_DIRECT','MOTHER_DIRECT']
    # Matings
    total_matings = ngen_random+ngen_am
    ## Record variance component evolution and pedigree
    if v_indirect == 0:
        V = np.zeros((total_matings+1,3))
        V_header = np.array(['v_g','v_y','r_delta']).reshape((1,3))
    else:
        V = np.zeros((total_matings+1,8))
        V_header = np.array(['v_g','v_y','r_delta','v_eg','c_ge','r_eta','r_delta_eta_c','r_delta_eta_tau']).reshape((1,8))
        pedcols += ['INDIRECT','FATHER_INDIRECT','MOTHER_INDIRECT']
    V[:] = np.nan
    ## Record full pedigree
    pedcols = np.array(pedcols)
    ped = np.zeros((nfam*2*(total_matings+1),pedcols.shape[0]),dtype='U30')
    # Count generations
    a_count = 0
    # Simulate ngen_random generations of random-mating then ngen_am generations of assortative mating
    for gen in range(0, total_matings):
        if v_indirect==0:
            if a_count==0:
                # Simulate direct effect component
                a, causal, h2 = simulate_effects(new_haps, n_causal, h2)
                # Compute parental phenotype
                delta_p, delta_m, Y_p, Y_m = compute_phenotype(haps, causal, a, 1 - h2)
                delta_p, delta_m, Y_p, Y_m = delta_p[father_indices], delta_m[mother_indices], Y_p[father_indices], Y_m[mother_indices]
            else:
                delta_p, delta_m, Y_p, Y_m = G_males[father_indices], G_females[mother_indices], Y_males[father_indices], Y_females[father_indices]
            # Compute phenotypes
            print('Computing phenotypes')
            if np.abs(beta_vert) > 0 and a_count>0:
                G_males, G_females, Y_males, Y_females = compute_phenotype_vert(new_haps, causal, a, 1-h2, beta_vert, Y_p, Y_m)
            else:
                G_males, G_females, Y_males, Y_females = compute_phenotype(new_haps, causal, a, 1 - h2)
        else:
            # Compute with indirect effects
            if a_count==0:
                # Compute indirect effect component
                a, causal, h2_total = simulate_effects(new_haps, n_causal, h2, old_haps=haps, v_indirect=v_indirect, 
                                                        r_direct_indirect=r_direct_indirect, father_indices=father_indices, mother_indices=mother_indices)
                # Compute parental phenotype
                delta_p, delta_m, Y_p, Y_m = compute_phenotype(haps, causal, a[:,0], 1 - h2)
                Y_p, Y_m = Y_p[father_indices], Y_m[mother_indices]
            # Compute phenotype
            print('Computing phenotypes')
            G_males, G_females, Y_males, Y_females, eta_p, eta_m, delta_p, delta_m = compute_phenotype_indirect(new_haps, haps, father_indices, mother_indices, 
                                                                                            causal, a[:,0], a[:,1], 1-h2_total)
            # Indirect effect component generation
            eta_males, eta_females = compute_genetic_component(new_haps, causal, a[:, 1])
        # Record variance components
        ## Final offspring generation variance components
        if v_indirect==0:
            V[a_count,:] = compute_vcomps(G_males, G_females, Y_males, Y_females, delta_p=delta_p, delta_m=delta_m)
        else:
            V[a_count,:] = compute_vcomps(G_males, G_females, Y_males, Y_females, delta_p=delta_p, delta_m=delta_m,
                                            v_indirect=v_indirect, eta_p=eta_p, eta_m=eta_m)
        # Record pedigree
        if v_indirect==0:
            ped[(a_count*2*nfam):((a_count+1)*2*nfam),:] = create_ped_output(G_males, G_females, Y_males, Y_females, Y_p, Y_m, delta_p, delta_m, 
                                                                                father_indices=father_indices, mother_indices=mother_indices, gen=a_count+1,header=False)
        else:
            ped[(a_count*2*nfam):((a_count+1)*2*nfam),:] = create_ped_output(G_males, G_females, Y_males, Y_females, Y_p, Y_m, delta_p, delta_m, 
                                    eta_males=eta_males, eta_females=eta_females, eta_p=eta_p, eta_m=eta_m, father_indices=father_indices, mother_indices=mother_indices, gen=a_count+1, header=False)
        a_count += 1
        ## Match current generation into parent-pairs for next generation
        print('Mating ' + str(gen + 2))
        # Random mating
        if gen<ngen_random:
            print('Matching randomly')
            father_indices = random_mating_indices(nfam)
            mother_indices = random_mating_indices(nfam)
        # Assortative mating
        if gen>=ngen_random:
            # Match assortatively
            print('Matching assortatively')
            father_indices, mother_indices = am_indices(Y_males, Y_females, r_y)
            print('Parental phenotypic correlation: '+str(round(np.corrcoef(Y_males[father_indices],Y_females[mother_indices])[0,1],4)))
        # Generate haplotypes of new generation
        haps = new_haps
        new_haps = []
        ibd = []
        for chr in range(0,len(haps)):
            if unlinked:
                new_haps_chr, ibd_chr = produce_next_gen_unlinked(father_indices,mother_indices,haps[chr][:,0,:,:],haps[chr][:,1,:,:])
            else:
                new_haps_chr, ibd_chr = produce_next_gen(father_indices,mother_indices,haps[chr][:,0,:,:],haps[chr][:,1,:,:],maps[chr])
            new_haps.append(new_haps_chr)
            ibd.append(ibd_chr)
    #### Compute final generation ####
    print('Computing final generation')
    if v_indirect==0:
        if a_count==0:
            # Simulate direct effect component
            a, causal, h2 = simulate_effects(new_haps, n_causal, h2)
            # Compute parental phenotype
            delta_p, delta_m, Y_p, Y_m = compute_phenotype(haps, causal, a, 1 - h2)
            delta_p, delta_m, Y_p, Y_m = delta_p[father_indices], delta_m[mother_indices], Y_p[father_indices], Y_m[mother_indices]
        else:
            delta_p, delta_m, Y_p, Y_m = G_males[father_indices], G_females[mother_indices], Y_males[father_indices], Y_females[father_indices]
        # Compute phenotypes
        print('Computing phenotypes')
        if np.abs(beta_vert) > 0 and a_count>0:
            G_males, G_females, Y_males, Y_females = compute_phenotype_vert(new_haps, causal, a, 1-h2, beta_vert, Y_p, Y_m)
        else:
            G_males, G_females, Y_males, Y_females = compute_phenotype(new_haps, causal, a, 1 - h2)
    else:
        # Compute with indirect effects
        if a_count==0:
            # Compute indirect effect component
            a, causal, h2_total = simulate_effects(new_haps, n_causal, h2, old_haps=haps, v_indirect=v_indirect, 
                                                    r_direct_indirect=r_direct_indirect, father_indices=father_indices, mother_indices=mother_indices)
            # Compute parental phenotype
            delta_p, delta_m, Y_p, Y_m = compute_phenotype(haps, causal, a[:,0], 1 - h2)
            Y_p, Y_m = Y_p[father_indices], Y_m[mother_indices]
        # Compute phenotype
        print('Computing phenotypes')
        G_males, G_females, Y_males, Y_females, eta_p, eta_m, delta_p, delta_m = compute_phenotype_indirect(new_haps, haps, father_indices, mother_indices, 
                                                                                        causal, a[:,0], a[:,1], 1-h2_total)
        # Indirect effect component generation
        eta_males, eta_females = compute_genetic_component(new_haps, causal, a[:, 1])
    print('Sibling genotypic correlation: ' + str(round(np.corrcoef(G_males, G_females)[0, 1], 4)))
    print('Sibling phenotypic correlation: ' + str(round(np.corrcoef(Y_males, Y_females)[0, 1], 4)))
    ## Final offspring generation variance components
    if v_indirect==0:
        V[a_count,:] = compute_vcomps(G_males, G_females, Y_males, Y_females, delta_p=delta_p, delta_m=delta_m)
    else:
        V[a_count,:] = compute_vcomps(G_males, G_females, Y_males, Y_females, delta_p=delta_p, delta_m=delta_m,
                                        v_indirect=v_indirect, eta_p=eta_p, eta_m=eta_m)
    V = np.vstack((V_header,V))
    # Record pedigree
    if v_indirect==0:
        ped[(a_count*2*nfam):((a_count+1)*2*nfam),:] = create_ped_output(G_males, G_females, Y_males, Y_females, Y_p, Y_m, delta_p, delta_m, 
                                                                            father_indices=father_indices, mother_indices=mother_indices, gen=a_count+1, header=False)
    else:
        ped[(a_count*2*nfam):((a_count+1)*2*nfam),:] = create_ped_output(G_males, G_females, Y_males, Y_females, Y_p, Y_m, delta_p, delta_m, 
                                eta_males=eta_males, eta_females=eta_females, eta_p=eta_p, eta_m=eta_m, father_indices=father_indices, mother_indices=mother_indices, gen=a_count+1, header=False)
    ped = np.vstack((pedcols,ped))
    ## Return output
    return new_haps, haps, father_indices, mother_indices, ibd, ped, a, V

##### Imputation functions #####
@njit
def impute_from_sibs(g1, g2, ibd, f):
    if ibd==2:
        return g1+2*f
    elif ibd==0:
        return g1+g2
    elif ibd==1:
        gsum = g1+g2
        if gsum==0:
            return f
        elif gsum==1:
            return 1+f
        elif gsum==2:
            return 1+2*f
        elif gsum==3:
            return 2+f
        elif gsum==4:
            return 3+f
@njit
def impute_from_sibs_phased(g1,g2,ibd,f):
    imp = 0.0
    if ibd[0]:
        imp += int(g1[0])+f
    else:
        imp += int(g1[0])+int(g2[0])
    if ibd[1]:
        imp += int(g1[1]) + f
    else:
        imp += int(g1[1]) + int(g2[1])
    return imp

@njit(parallel=True)
def impute_all_fams(gts,freqs,ibd):
    imp = np.zeros((gts.shape[0],gts.shape[2]),dtype=np.float_)
    for i in prange(gts.shape[0]):
        for j in range(gts.shape[2]):
            imp[i,j] = impute_from_sibs(gts[i,0,j],gts[i,1,j],ibd[i,j],freqs[j])
    return imp

@njit(parallel=True)
def impute_all_fams_phased(haps,freqs,ibd):
    imp = np.zeros((haps.shape[0],haps.shape[2]),dtype=np.float_)
    for i in prange(haps.shape[0]):
        for j in range(haps.shape[2]):
            imp[i,j] = impute_from_sibs_phased(haps[i,0,j,:],haps[i,1,j,:],ibd[i,j,:],freqs[j])
    return imp