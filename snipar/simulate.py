import numpy as np
from numba import njit, prange

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

@njit
def simulate_recombinations(map):
    map_start = map[0]
    map_end = map[map.shape[0]-1]
    map_length = map_end-map_start
    n_recomb = np.random.poisson(map_length/100)
    recomb_points = map_start+np.sort(np.random.uniform(0,1,n_recomb))*map_length
    return n_recomb,recomb_points

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


@njit(parallel=True)
def produce_next_gen(father_indices,mother_indices,males,females,map):
    ngen = np.zeros((father_indices.shape[0],2,males.shape[1],2),dtype=np.bool_)
    ibd = np.zeros((father_indices.shape[0],males.shape[1],2),dtype=np.bool_)
    for i in prange(father_indices.shape[0]):
        # recombinations
        recomb_i = meiosis(map,n=4)
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
    G_males, G_females = compute_genetic_component(haps,causal,a)
    G_p, G_m = compute_genetic_component(old_haps, causal, b)
    Y_males = G_males+G_p[father_indices]+G_m[mother_indices]+np.sqrt(sigma2)*np.random.randn(G_males.shape[0])
    Y_females = G_females+G_p[father_indices]+G_m[mother_indices]+np.sqrt(sigma2)*np.random.randn(G_males.shape[0])
    return G_males, G_females, Y_males, Y_females