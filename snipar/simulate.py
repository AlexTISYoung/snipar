import numpy as np
from numba import njit, prange
import gzip, h5py, os

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
def produce_next_gen(father_indices, mother_indices, males, females, map):
    ngen = np.zeros((father_indices.shape[0],2,males.shape[1],2), dtype=np.bool_)
    ibd = np.zeros((father_indices.shape[0],males.shape[1]), dtype=np.int_)
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
            ibd[i,j] = int(recomb_i[0,j]==recomb_i[1,j])+int(recomb_i[2,j]==recomb_i[3,j])
    return ngen, ibd

@njit
def random_mating_indices(nfam):
    return np.random.choice(np.array([x for x in range(nfam)],dtype=np.int_),size=nfam,replace=False)


class generation(object):
    def __init__(self, haps, phenotype=None):
        self.haps = haps
        self.phenotype = phenotype

    def produce_next_gen(self, father_indices, mother_indices):



def read_haps(bgenfiles):
    bgenfiles, chroms = preprocess.parse_obsfiles(bgenfiles, obsformat='bgen')
    # Find common IDs
    print('Finding common individuals between bgen files')
    bgen = open_bgen(bgenfiles[0])
    common_ids = set(bgen.samples)
    for i in range(1, bgenfiles.shape[0]):
        bgen = open_bgen(bgenfiles[i])
        common_ids = common_ids.intersection(bgen.samples)
    common_ids = list(common_ids)
    if len(common_ids) == 0:
        raise (ValueError('No IDs in common between bgen files'))
    else:
        print('Found ' + str(len(common_ids)) + ' IDs in common between bgen files')
    # Read genotypes
    haps = []
    for i in bgenfiles.shape[0]:
        print('Reading in chromosome ' + str(chroms[i]))
        # Read bgen
        bgen = open_bgen(bgenfiles[i], verbose=True)
        # Read map
        print('Finding genetic map positions')
        map_i = preprocess.decode_map_from_pos(chroms[i], bgen.positions)
        # Alleles
        alleles = np.array([x.split(',') for x in bgen.allele_ids])
        # Find indices of common IDs
        bgen_id_dict = preprocess.make_id_dict(bgen.samples)
        id_indices = [bgen_id_dict[x] for x in common_ids]
        # Read genotypes
        gts = bgen.read((id_indices,
                         [x for x in range(bgen.ids.shape[0])]), np.bool_)[:, :, np.array([0, 2])]
        nfam = int(np.floor(gts.shape[0] / 2))
        haps.append(gtarray(gts[0:(2 * nfam), :, :], ids=bgen.samples[id_indices], fams=np.arange(1, nfam),
                            sid=bgen.ids, map=map_i, alleles=alleles, chrom=chroms[i]))
        del gts
    return haps


def simulate_pop(haps, ngen_random, ngen_am, ncausal):
    total_matings = ngen_random + ngen_am
    V = np.zeros((total_matings + 1, 2))
    a_count = 0
    # Produce next generation
    old_haps = haps
    for gen in range(0, total_matings):
        # Simulate phenotype for AM
        if gen == 0 and h2_total == 0:
            print('Simulating phenotype')
            # Simulate phenotype
            nsnp_chr = np.array([x.shape[2] for x in old_haps])
            nsnp = np.sum(nsnp_chr)
            if ncausal > nsnp:
                raise (ValueError('Not enough SNPs to simulate phenotype with ' + str(ncausal) + ' causal SNPs'))
            a = np.zeros((nsnp))
            causal = np.sort(np.random.choice(np.arange(0, nsnp), ncausal, replace=False))
            a[causal] = np.random.randn(ncausal)
            G_p, G_m = compute_genetic_component(old_haps, causal, a)
            scale_fctr = np.sqrt(h2 / np.var(np.hstack((G_p, G_m))))
            a = a * scale_fctr
        # Compute parental phenotypes
        if np.abs(beta_vert) > 0 and gen > 0:
            print('Computing parental phenotypes')
            G_p, G_m, Y_p, Y_m = compute_phenotype_vert(old_haps, causal, a, 1 - h2, beta_vert, Y_p[father_indices],
                                                        Y_m[mother_indices])
        elif h2_total == 0:
            print('Computing parental phenotypes')
            G_p, G_m, Y_p, Y_m = compute_phenotype(old_haps, causal, a, 1 - h2)
        # Record variance components
        if gen > 0 or h2_total == 0:
            V[a_count, :] = np.array([np.var(np.hstack((G_p, G_m))), np.var(np.hstack((Y_p, Y_m)))])
        print('Genetic variance: ' + str(round(V[a_count, 0], 4)))
        print('Phenotypic variance: ' + str(round(V[a_count, 1], 4)))
        print('Heritability: ' + str(round(V[a_count, 0] / V[a_count, 1], 4)))
        a_count += 1
        ## Match parents
        print('Mating ' + str(gen + 1))
        # Random mating
        if gen < ngen_random:
            father_indices = random_mating_indices(nfam)
            mother_indices = random_mating_indices(nfam)
        # Assortative mating
        if gen >= ngen_random:
            # Compute parental phenotypes
            print('Computing parental phenotypes')
            # Match assortatively
            print('Matching assortatively')
            father_indices, mother_indices = am_indices(Y_p, Y_m, 0.5)
            # Print variances
            print('Parental phenotypic correlation: ' + str(
                round(np.corrcoef(Y_p[father_indices], Y_m[mother_indices])[0, 1], 4)))
            print('Parental genotypic correlation: ' + str(
                round(np.corrcoef(G_p[father_indices], G_m[mother_indices])[0, 1], 4)))
        # Generate haplotpyes of new generation
        new_haps = []
        ibd = []
        for chr in range(0, len(haps)):
            print('Chromosome ' + str(chr_start + chr))
            new_haps_chr, ibd_chr = produce_next_gen(father_indices, mother_indices, old_haps[chr][:, 0, :, :],
                                                     old_haps[chr][:, 1, :, :], maps[chr])
            new_haps.append(new_haps_chr)
            ibd.append(ibd_chr)
        # Compute indirect effect component
        if h2_total > 0:
            if gen == 0:
                print('Simulating indirect genetic effects')
                nsnp_chr = np.array([x.shape[2] for x in old_haps])
                nsnp = np.sum(nsnp_chr)
                if ncausal > nsnp:
                    raise (ValueError('Not enough SNPs to simulate phenotype with ' + str(ncausal) + ' causal SNPs'))
                ab = np.zeros((nsnp, 2))
                causal = np.sort(np.random.choice(np.arange(0, nsnp), ncausal, replace=False))
                ab[causal, :] = np.random.multivariate_normal(np.zeros((2)),
                                                              np.array([[1, r_dir_alpha], [r_dir_alpha, 1]]),
                                                              size=ncausal)
                G_males, G_females, Y_males, Y_females = compute_phenotype_indirect(new_haps, old_haps, father_indices,
                                                                                    mother_indices, causal, ab[:, 0],
                                                                                    ab[:, 1], 0)
                scale_fctr = np.sqrt(h2_total / np.var(np.hstack((Y_males, Y_females))))
                ab = ab * scale_fctr
            print('Computing parental phenotype')
            G_p, G_m, Y_p, Y_m = compute_phenotype_indirect(new_haps, old_haps, father_indices, mother_indices, causal,
                                                            ab[:, 0], ab[:, 1], 1 - h2_total)
        if gen < (total_matings - 1):
            old_haps = new_haps

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