import numpy as np
from numba import njit, prange
from snipar.utilities import *
from snipar.map import decode_map_from_pos

### Read haplotypes from bgen
def haps_from_bgen(bgenfiles, chr_range, nfam=None):
    bgenfiles, chroms = parse_obsfiles(bgenfiles, obsformat="bgen", chromosomes=chr_range)
    haps, maps, snp_ids, alleles, positions = [], [], [], [], []
    samples = None
    for i, file in enumerate(bgenfiles):
        bgen = open_bgen(file, verbose=True)
        if chroms[i] == 0:
            file_chroms = np.unique(bgen.chromosomes)
            if len(file_chroms) != 1 or not str(file_chroms[0]).isdigit():
                raise ValueError("Supply a single BGEN chromosome, or use a chromosome wildcard")
            chroms[i] = int(file_chroms[0])
        print(f"Reading in chromosome {chroms[i]}")
        if samples is None:
            n_pairs = len(bgen.samples) // 2 if nfam is None else nfam
            if not 1 <= n_pairs <= len(bgen.samples) // 2:
                raise ValueError("Number of requested families exceeds the available pairs")
            samples = bgen.samples[:2 * n_pairs]
            rows = np.arange(len(samples))
        else:
            sample_indices = make_id_dict(bgen.samples)
            try:
                rows = np.array([sample_indices[iid] for iid in samples])
            except KeyError as exc:
                raise ValueError(f"BGEN chromosomes disagree on sample {exc.args[0]}") from exc
        if len(set(samples)) != len(samples):
            raise ValueError("BGEN files must have unique sample IDs")
        genetic_map = decode_map_from_pos(chroms[i], bgen.positions)
        markers = np.flatnonzero(np.isfinite(genetic_map))
        if not len(markers):
            raise ValueError(f"No SNPs on chromosome {chroms[i]} have genetic map positions")
        order = np.argsort(bgen.positions[markers], kind="stable")
        markers = markers[order]
        probs = bgen.read((rows, markers), np.float32)
        if probs.shape[2] != 4 or not np.isfinite(probs).all():
            raise ValueError("Found missing or unphased genotypes; simulation requires complete phased biallelic BGEN input")
        calls = probs[:, :, [0, 2]] >= .5
        haps.append(calls.reshape((2, n_pairs, len(markers), 2)).transpose(1, 0, 2, 3).copy())
        maps.append(genetic_map[markers])
        positions.append(bgen.positions[markers])
        snp_ids.append(bgen.rsids[markers])
        alleles.append(np.array([x.split(",") for x in bgen.allele_ids[markers]]))
    return haps, maps, snp_ids, alleles, positions, chroms

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
            u = np.random.uniform(0, 1, nsnp)
            freqs = min_maf * np.exp(-u * np.log(2 * min_maf))
        else:
            raise(ValueError('minimum MAF must be between zero and 0.5'))
    return freqs

@njit(cache=True, parallel=True)
def simulate_haplotypes(freqs, nfam):
    haps = np.empty((nfam, 2, len(freqs), 2), dtype=np.bool_)
    for family in prange(nfam):
        for sex in range(2):
            for snp in range(len(freqs)):
                haps[family, sex, snp, 0] = np.random.random() < freqs[snp]
                haps[family, sex, snp, 1] = np.random.random() < freqs[snp]
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
        last_hap = np.random.randint(2) != 0
        recomb_vector[r,:] = np.bool_(last_hap)
        # Recombine
        if n_recomb>0:
            for i in range(n_recomb):
                first_snp = np.searchsorted(map, recomb_points[i], side="right")
                recomb_vector[r,first_snp:recomb_vector.shape[1]] = ~recomb_vector[r,first_snp:recomb_vector.shape[1]]
    # Return
    return recomb_vector


##### Produce next generation for unlinked variants #####
@njit(cache=True, parallel=True)
def produce_next_gen_unlinked(father_indices, mother_indices, males, females):
    n_families = len(father_indices)
    n_snps = males.shape[1]
    children = np.empty((n_families, 2, n_snps, 2), dtype=np.bool_)
    ibd = np.empty((n_families, n_snps, 2), dtype=np.bool_)
    for i in prange(n_families):
        father, mother = father_indices[i], mother_indices[i]
        # Four independent fair transmissions per SNP, eight SNPs per draw.
        for start in range(0, n_snps, 8):
            bits = np.random.randint(0, 2**32)
            for j in range(start, min(start + 8, n_snps)):
                father_a, father_b = bits & 1, (bits >> 1) & 1
                mother_a, mother_b = (bits >> 2) & 1, (bits >> 3) & 1
                children[i, 0, j, 0] = males[father, j, father_a]
                children[i, 1, j, 0] = males[father, j, father_b]
                children[i, 0, j, 1] = females[mother, j, mother_a]
                children[i, 1, j, 1] = females[mother, j, mother_b]
                ibd[i, j, 0] = father_a == father_b
                ibd[i, j, 1] = mother_a == mother_b
                bits >>= 4
    return children, ibd

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

def am_indices(yp, ym, r):
    if not -1 <= r <= 1:
        raise ValueError("Parental correlation must be between -1 and 1")
    if r == 0:
        return np.random.permutation(len(yp)), np.random.permutation(len(ym))
    # Adding independent noise attenuates same-rank mating to the requested correlation.
    noise_scale = np.sqrt(1 / abs(r) - 1)
    zp = yp + noise_scale * np.std(yp) * np.random.randn(len(yp))
    zm = ym + noise_scale * np.std(ym) * np.random.randn(len(ym))
    rank_p = np.argsort(zp)
    rank_m = np.argsort(zm)
    if r < 0:
        rank_m = rank_m[::-1]
    return rank_p, rank_m

@njit(cache=True, parallel=True)
def _add_chromosome_scores(chromosome, markers, effects, paternal, maternal):
    for row in prange(chromosome.shape[0]):
        p = m = 0.0
        for i in range(len(markers)):
            snp = markers[i]
            p += effects[i] * (int(chromosome[row, 0, snp, 0]) + int(chromosome[row, 0, snp, 1]))
            m += effects[i] * (int(chromosome[row, 1, snp, 0]) + int(chromosome[row, 1, snp, 1]))
        paternal[row] += p
        maternal[row] += m


def compute_genetic_component(haps, causal, a):
    causal = np.unique(np.fromiter(causal, dtype=np.intp))
    paternal = np.zeros(haps[0].shape[0])
    maternal = np.zeros(haps[0].shape[0])
    start = 0
    for chromosome in haps:
        end = start + chromosome.shape[2]
        markers = causal[np.searchsorted(causal, start):np.searchsorted(causal, end)]
        _add_chromosome_scores(chromosome, markers - start, a[markers], paternal, maternal)
        start = end
    return paternal, maternal

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
        a[causal] = np.random.randn(n_causal)
        G_p, G_m = compute_genetic_component(haps,causal,a)
        variance = np.var(np.hstack((G_p, G_m)))
        if h2 > 0 and variance <= 0:
            raise ValueError("Causal variants have no genetic variance in the simulated cohort")
        a *= np.sqrt(h2 / variance) if h2 > 0 else 0
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
        direct_variance = np.var(np.hstack((G_males, G_females)))
        indirect_variance = np.var(eta_p + eta_m)
        if h2 > 0 and (direct_variance <= 0 or indirect_variance <= 0):
            raise ValueError("Causal variants have no genetic variance in the simulated cohort")
        a[:, 0] *= np.sqrt(h2 / direct_variance) if h2 > 0 else 0
        a[:, 1] *= np.sqrt(v_eg / indirect_variance) if v_eg > 0 else 0
    return a, causal, h2_total


def stabilizing_selection_weights(phenotype, variance, optimum=0.0):
    """Gaussian viability fitness; variance is in fixed baseline phenotype units."""
    phenotype = np.asarray(phenotype, dtype=float)
    if phenotype.ndim != 1 or not np.isfinite(phenotype).all():
        raise ValueError('Selection phenotypes must be a finite one-dimensional array')
    if np.isnan(variance) or variance <= 0 or not np.isfinite(optimum):
        raise ValueError('Selection variance must be positive and optimum must be finite')
    return np.exp(-0.5 * (phenotype - optimum) ** 2 / variance)


def stabilizing_selection_diagnostics(haps, effects, genetic, phenotype, indices):
    """Partition genetic variance into diagonal and effect-signed LD components."""
    scores = np.concatenate([genetic[sex][indices[sex]] for sex in range(2)])
    values = np.concatenate([phenotype[sex][indices[sex]] for sex in range(2)])
    genic, diagonal, cursor = 0.0, 0.0, 0
    for chromosome in haps:
        dosage = np.concatenate([
            chromosome[indices[sex], sex].sum(axis=2, dtype=np.uint8) for sex in range(2)
        ])
        mean = dosage.mean(axis=0)
        variance = (dosage * dosage).mean(axis=0) - mean ** 2
        frequency = mean / 2
        squared_effect = effects[cursor:cursor + dosage.shape[1]] ** 2
        diagonal += np.dot(variance, squared_effect)
        genic += np.dot(2 * frequency * (1 - frequency), squared_effect)
        cursor += dosage.shape[1]
    genetic_variance = np.var(scores)
    return dict(n_individuals=len(scores), genetic_variance=float(genetic_variance),
                phenotypic_variance=float(np.var(values)), mean_phenotype=float(np.mean(values)),
                genic_variance=float(genic), genotypic_diagonal_variance=float(diagonal),
                signed_offdiagonal_ld_variance=float(genetic_variance - diagonal),
                residual_variance=float(np.var(values - scores)),
                genetic_residual_covariance=float(np.mean((scores - scores.mean()) *
                                                        (values - scores - np.mean(values - scores)))))


def forward_sim_stabilizing_selection(haps, maps, n_generations, n_causal, h2,
                                     selection_variance, nfam_output, unlinked=True,
                                     optimum=0.0):
    """Select breeders over generations, then sample offspring before selection.

    Effects, residual variance and the phenotype origin remain fixed. Neutral
    controls use selection_variance=inf. Final parents are unique within sex;
    earlier generations use survivor sampling with replacement at constant size.
    """
    if int(n_generations) != n_generations or n_generations < 1:
        raise ValueError('At least one integer generation of selection is required')
    if not 0 < h2 <= 1:
        raise ValueError('Selection simulations require 0 < h2 <= 1')
    n_candidates = haps[0].shape[0]
    if int(nfam_output) != nfam_output or not 1 <= nfam_output <= n_candidates:
        raise ValueError('Output family count must be a positive integer no larger than the breeding pool')
    stabilizing_selection_weights(np.zeros(1), selection_variance, optimum)
    effects, causal, _ = simulate_effects(haps, n_causal, h2)
    initial_genetic = compute_genetic_component(haps, causal, effects)
    origin = np.mean(np.concatenate(initial_genetic))
    genetic = tuple(x - origin for x in initial_genetic)
    phenotype = tuple(x + np.sqrt(1 - h2) * np.random.randn(len(x)) for x in genetic)
    history, variance_rows = [], []
    for generation in range(n_generations):
        all_indices = (np.arange(n_candidates), np.arange(n_candidates))
        history.append(dict(generation=generation, phase='before_selection',
                            **stabilizing_selection_diagnostics(haps, effects, genetic, phenotype, all_indices)))
        survivors = tuple(np.flatnonzero(np.random.random(n_candidates) <
                          stabilizing_selection_weights(phenotype[sex], selection_variance, optimum))
                          for sex in range(2))
        if any(len(x) == 0 for x in survivors):
            raise ValueError('No surviving breeders in one sex; increase the breeding pool or weaken selection')
        history.append(dict(generation=generation, phase='after_selection',
                            **stabilizing_selection_diagnostics(haps, effects, genetic, phenotype, survivors)))
        final = generation == n_generations - 1
        if final and any(len(x) < nfam_output for x in survivors):
            raise ValueError('Not enough unique surviving final parents; increase --selection_nfam')
        n_pairs = nfam_output if final else n_candidates
        father_indices, mother_indices = tuple(
            np.random.choice(x, n_pairs, replace=not final) for x in survivors
        )
        # From generation 1 onward each index labels a brother/sister pair.
        if generation > 0:
            same_family = father_indices == mother_indices
            attempts = 0
            while same_family.any():
                if final:
                    np.random.shuffle(mother_indices)
                else:
                    mother_indices[same_family] = np.random.choice(survivors[1], same_family.sum())
                same_family = father_indices == mother_indices
                attempts += 1
                if attempts > 1000:
                    raise ValueError('Unable to find matings without full-sibling pairs')
        new_haps, ibd = [], []
        for chromosome, recombination_map in zip(haps, maps):
            args = (father_indices, mother_indices, chromosome[:, 0], chromosome[:, 1])
            result = produce_next_gen_unlinked(*args) if unlinked else produce_next_gen(*args, recombination_map)
            new_haps.append(result[0])
            ibd.append(result[1])
        new_genetic = tuple(x - origin for x in compute_genetic_component(new_haps, causal, effects))
        new_phenotype = tuple(x + np.sqrt(1 - h2) * np.random.randn(n_pairs) for x in new_genetic)
        delta_p, delta_m = genetic[0][father_indices], genetic[1][mother_indices]
        variance_rows.append(compute_vcomps(*new_genetic, *new_phenotype,
                             delta_p=delta_p, delta_m=delta_m, verbose=False))
        print('Selection generation {}: {} male and {} female survivors; offspring genetic variance {:.6g}'.format(
            generation + 1, len(survivors[0]), len(survivors[1]), variance_rows[-1][0]))
        if final:
            missing = np.full(n_candidates, np.nan)
            parent_ped = create_ped_output(*genetic, *phenotype, missing, missing, missing, missing,
                                           gen=0, header=True)
            child_ped = create_ped_output(*new_genetic, *new_phenotype,
                                         phenotype[0][father_indices], phenotype[1][mother_indices],
                                         delta_p, delta_m, father_indices=father_indices,
                                         mother_indices=mother_indices, gen=1, header=False)
            ped = np.vstack((parent_ped, child_ped))
            history.append(dict(generation=n_generations, phase='offspring_before_selection',
                                **stabilizing_selection_diagnostics(new_haps, effects, new_genetic, new_phenotype,
                                                                  (np.arange(n_pairs), np.arange(n_pairs)))))
        else:
            haps, genetic, phenotype = new_haps, new_genetic, new_phenotype
    variance = np.vstack((np.array(['v_g', 'v_y', 'r_delta']), np.asarray(variance_rows)))
    metadata = dict(selection_variance=float(selection_variance), optimum=float(optimum),
                    n_generations=int(n_generations), baseline_h2=float(h2), residual_variance=1 - h2,
                    genetic_origin=float(origin), breeding_pairs=n_candidates, output_families=int(nfam_output),
                    offspring_sampling='before selection; no phenotype ascertainment',
                    effects_rescaled_after_selection=False, history=history)
    return new_haps, haps, father_indices, mother_indices, ibd, ped, effects, variance, metadata


def selection_family_relationships(ped, parent_siblings=True):
    """Sibling and first-cousin relationships when only some breeders are sampled."""
    children = ped[np.char.startswith(ped[:, 0], '1_')]
    rows = [['ID1', 'ID2', 'PropIBD', 'InfType']]
    rows += [[children[i, 1], children[i + 1, 1], '0.5', 'FS']
             for i in range(0, len(children), 2)]
    if parent_siblings:
        paternal, maternal = {}, {}
        for child in children:
            paternal.setdefault(child[2], []).append(child[1])
            maternal.setdefault(child[3], []).append(child[1])
        for father, paternal_children in paternal.items():
            sister = father.rsplit('_', 1)[0] + '_1'
            for first in paternal_children:
                for second in maternal.get(sister, []):
                    if first.rsplit('_', 1)[0] != second.rsplit('_', 1)[0]:
                        rows.append([first, second, '0.125', '3rd'])
    return np.asarray(rows)

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
                delta_p, delta_m = G_males[father_indices], G_females[mother_indices]
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
        # Save the phenotypes used for this mating before computing the offspring.
        Y_p, Y_m = Y_males[father_indices], Y_females[mother_indices]
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
            delta_p, delta_m = G_males[father_indices], G_females[mother_indices]
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
    imp = np.zeros((gts.shape[0],gts.shape[2]),dtype=np.float64)
    for i in prange(gts.shape[0]):
        for j in range(gts.shape[2]):
            imp[i,j] = impute_from_sibs(gts[i,0,j],gts[i,1,j],ibd[i,j],freqs[j])
    return imp

@njit(parallel=True)
def impute_all_fams_phased(haps,freqs,ibd):
    imp = np.zeros((haps.shape[0],haps.shape[2]),dtype=np.float64)
    for i in prange(haps.shape[0]):
        for j in range(haps.shape[2]):
            imp[i,j] = impute_from_sibs_phased(haps[i,0,j,:],haps[i,1,j,:],ibd[i,j,:],freqs[j])
    return imp

def find_sibs_and_cousins(ped):
    # Get last two generations
    n_last = ped[ped.shape[0]-1,0].split('_')[0]
    last_gen = ped[[x.split('_')[0]==n_last for x in ped[:,0]],:]
    parent_gen = ped[[x.split('_')[0]==str(int(n_last)-1) for x in ped[:,0]],:]
    # Fill in sibling relatedness
    sibs = np.zeros((int(last_gen.shape[0]/2),4),dtype=ped.dtype)
    sibs[:,0] = last_gen[0::2,1]
    sibs[:,1] = last_gen[1::2,1]
    sibs[:,2] = 0.5
    sibs[:,3] = 'FS'
    # Fill in cousin relatedness
    cousins = np.zeros((int(2*last_gen.shape[0]),4),dtype=ped.dtype)
    father_offspring_by_id, mother_offspring_by_id = {}, {}
    for _, iid, father, mother, *_ in last_gen:
        father_offspring_by_id.setdefault(father, []).append(iid)
        mother_offspring_by_id.setdefault(mother, []).append(iid)
    ccount = 0
    for i in range(0,parent_gen.shape[0],2):
        father_offspring = father_offspring_by_id[parent_gen[i, 1]]
        mother_offspring = mother_offspring_by_id[parent_gen[i + 1, 1]]
        for j in range(0,2):
            for k in range(0,2):
                cousins[ccount,:] = [father_offspring[j], mother_offspring[k], 0.125,
                                    '3rd']
                ccount += 1
    # Remove inbred 'cousins'
    cousins = cousins[cousins[:,0] != cousins[:,1],:]
    # Create relatedness file format output
    rel = np.vstack((np.array(['ID1', 'ID2', 'PropIBD', 'InfType']),sibs,cousins))
    return rel
