"""Make benchmark sibships from the expanded 1000 Genomes EUR haplotypes."""

import argparse
from collections import Counter
from contextlib import ExitStack
import gzip
from itertools import combinations
import json
from pathlib import Path

import numpy as np


PHASED_GT = np.array([b"0|0\t", b"0|1\t", b"1|0\t", b"1|1\t"])
UNPHASED_GT = np.char.replace(PHASED_GT, b"|", b"/")
BED_CODES = np.array([3, 2, 0], dtype=np.uint8)
BED_SHIFTS = np.array([0, 2, 4, 6], dtype=np.uint8)
BASES = {b"A", b"C", b"G", b"T"}


def vcf_header(source):
    for line in source:
        if line.startswith(b"#CHROM\t"):
            return [name.decode() for name in line.split()[9:]]
    raise ValueError("VCF sample header is missing")


def pedigrees(path, samples, count, unrelated, rng):
    with path.open() as file:
        columns = next(file).split()
        expected = "FamilyID SampleID FatherID MotherID Sex Population Superpopulation".split()
        if columns != expected:
            raise ValueError("Expected the IGSR 3,202-sample pedigree with population labels")
        eur = [row for line in file if (row := line.split())[-1] == "EUR"]
    eur_ids = {row[1] for row in eur if row[1] in samples}
    trios = [row for row in eur if row[1] in eur_ids and set(row[2:4]) <= eur_ids]
    trio_parent_pairs = Counter(tuple(row[2:4]) for row in trios)
    full_sibships = sum(children > 1 for children in trio_parent_pairs.values())
    # Keep one person per pedigree from families with no recorded EUR trio or duo.
    related_families = {row[0] for row in eur if any(parent != "0" for parent in row[2:4])}
    parent_ids = {parent for row in eur for parent in row[2:4]}
    candidates = {}
    for family, iid, _, _, _, pop, _ in eur:
        if family not in related_families and iid not in parent_ids and iid in eur_ids:
            candidates.setdefault(family, (iid, pop))
    pairs = []
    for pop in sorted({pop for _, pop in candidates.values()}):
        individuals = sorted(iid for iid, population in candidates.values() if population == pop)
        rng.shuffle(individuals)
        pairs.extend(zip(individuals[::2], individuals[1::2]))
    if len(pairs) < count:
        raise ValueError(f"Only {len(pairs)} population-matched EUR founder pairs are available")
    rng.shuffle(pairs)
    founders = np.array(pairs[:count])
    remainder = sorted({iid for iid, _ in candidates.values()} - set(founders.ravel()))
    if len(remainder) < unrelated:
        raise ValueError(f"Only {len(remainder)} independent EUR samples remain after pairing founders")
    rng.shuffle(remainder)
    return founders, remainder[:unrelated], eur_ids, len(trios), full_sibships


def write_bed_genotypes(output, haplotypes):
    packed = np.zeros((len(haplotypes) + 3) // 4 * 4, dtype=np.uint8)
    packed[:len(haplotypes)] = BED_CODES[haplotypes.sum(axis=1)]
    encoded = np.bitwise_or.reduce(packed.reshape(-1, 4) << BED_SHIFTS, axis=1)
    output.write(encoded.tobytes())


def write_vcf_header(output, sample_ids):
    output.write(b'##fileformat=VCFv4.2\n##contig=<ID=chr22,length=50818468>\n')
    output.write(b'##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count">\n')
    output.write(b'##INFO=<ID=AN,Number=1,Type=Integer,Description="Called allele count">\n')
    output.write(b'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    output.write(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
    output.write("\t".join(sample_ids).encode() + b"\n")


def write_vcf_genotypes(output, fields, haplotypes, codes, allele_count):
    genotypes = haplotypes[:, 0] * 2 + haplotypes[:, 1]
    counts = f"AC={allele_count};AN={haplotypes.size}".encode()
    output.write(b"\t".join(fields[:7] + [counts, b"GT", codes[genotypes].tobytes()[:-1]]) + b"\n")


def pedigree_files(destination, founders, visible, unrelated, rng):
    n, width = visible.shape
    families = np.array([f"F{index:04d}" for index in range(1, n + 1)])
    samples = np.array([[f"SIM22_{family}_{sibling}" for sibling in range(1, width - 1)]
                        for family in families])
    ids = np.concatenate((founders, samples), axis=1)
    with ((destination / "pedigree.txt").open("w") as ped,
          (destination / "chr_22.fam").open("w") as fam,
          (destination / "truth" / "chr_22.fam").open("w") as truth_fam):
        ped.write("FID IID FATHER_ID MOTHER_ID\n")
        for family, siblings, include in zip(families, ids, visible):
            for index, iid in enumerate(siblings):
                father, mother = siblings[:2] if index >= 2 else ("0", "0")
                sex = 1 if index == 0 else 2 if index == 1 else 0
                ped.write(f"{family} {iid} {father} {mother}\n")
                if include[index]:
                    fam.write(f"{family} {iid} {father} {mother} {sex} -9\n")
                if index < 2:
                    truth_fam.write(f"{family} {iid} 0 0 {sex} -9\n")
        for iid in unrelated:
            ped.write(f"U_{iid} {iid} 0 0\n")
            fam.write(f"U_{iid} {iid} 0 0 0 -9\n")
    # Keep whole families together in the GWAS discovery and PGS evaluation splits.
    discovery = rng.permutation(n)[:round(n * 0.7)]
    training = np.isin(np.arange(n), discovery)
    child_present = visible[:, 2:]
    discovery_ids = samples[training[:, None] & child_present]
    evaluation_ids = samples[~training[:, None] & child_present]
    np.savetxt(destination / "discovery.keep", np.concatenate((discovery_ids, unrelated)), fmt="%s")
    np.savetxt(destination / "target.keep", evaluation_ids, fmt="%s")
    return families, ids, training


def transmissions(path, families, ids, first_cm, last_cm, rng):
    shape = (len(families), ids.shape[1] - 2, 2)
    initial = rng.integers(0, 2, size=shape, dtype=np.uint8)
    switches = []
    events = {}
    with path.open("w") as file:
        file.write("FID\tIID\tFATHER_ID\tMOTHER_ID\tpaternal_start\tmaternal_start\t"
                   "paternal_crossover_cM\tmaternal_crossover_cM\n")
        for family_index, (family, siblings) in enumerate(zip(families, ids)):
            father, mother = siblings[:2]
            for child_index, child in enumerate(siblings[2:]):
                crossovers = []
                for parent in range(2):
                    n_crossovers = rng.poisson((last_cm - first_cm) / 100)
                    crossover = np.sort(rng.uniform(first_cm, last_cm, size=n_crossovers))
                    events[family_index, child_index, parent] = crossover
                    path_index = np.ravel_multi_index((family_index, child_index, parent), shape)
                    switches.extend((cm, path_index) for cm in crossover)
                    crossovers.append(",".join(format(cm, ".15g") for cm in crossover) or ".")
                paternal, maternal = initial[family_index, child_index]
                file.write(f"{family}\t{child}\t{father}\t{mother}\t{paternal}\t{maternal}\t"
                           f"{crossovers[0]}\t{crossovers[1]}\n")
    return initial, iter(sorted(switches)), events


def write_relationships(path, ids, visible, initial, events, first_cm, last_cm):
    with path.open("w") as file:
        file.write("ID1\tID2\tPropIBD\tInfType\n")
        for family_index, family in enumerate(ids):
            children = family[2:]
            genotyped_children = np.flatnonzero(visible[family_index, 2:])
            for parent, observed in zip(family[:2], visible[family_index, :2]):
                if observed:
                    for child in children[genotyped_children]:
                        file.write(f"{parent}\t{child}\t0.5\tPO\n")
            for left, right in combinations(genotyped_children, 2):
                relatedness = 0.
                for parent in range(2):
                    left_path = family_index, left, parent
                    right_path = family_index, right, parent
                    left_start, right_start = initial[left_path], initial[right_path]
                    if last_cm > first_cm:
                        left_cross, right_cross = events[left_path], events[right_path]
                        positions = np.sort(np.concatenate(([first_cm], left_cross, right_cross, [last_cm])))
                        lengths = np.diff(positions)
                        midpoints = positions[:-1] + lengths / 2
                        same = ((left_start + np.searchsorted(left_cross, midpoints)) % 2
                                == (right_start + np.searchsorted(right_cross, midpoints)) % 2)
                        relatedness += lengths[same].sum() / (last_cm - first_cm) / 2
                    else:
                        relatedness += (left_start == right_start) / 2
                # A short region can be fully IBD2 for two non-twin siblings.
                # Mix local sharing with pedigree expectation for the sparse GRM.
                grm = .25 + .5 * relatedness
                file.write(f"{children[left]}\t{children[right]}\t{grm:.12g}\tFS\n")


def read_diploid_gt(fields, sample_count, eur_index):
    calls = fields[9]
    full_size = len(calls) - (2 if calls.endswith(b"\r\n") else int(calls.endswith(b"\n")))
    if fields[8] != b"GT" or full_size != 4 * sample_count - 1:
        raise ValueError("Expected the NYGC panel's GT-only diploid phased calls")
    alleles = np.ndarray((sample_count, 3), dtype=np.uint8, buffer=calls, strides=(4, 1))
    return alleles[eur_index, ::2] - ord("0")


def simulate_variants(args, source, sample_count, eur_ids, eur_index, founder_index, unrelated_index,
                      cm_map, visible, present_ids, initial, switches, rng):
    destination = args.output
    parent_index = np.arange(2 * len(visible)).reshape(-1, 1, 2)
    transmission = initial.reshape(-1)
    switch_cm, switch_index = next(switches, (np.inf, -1))
    counts = dict(input_records=0, eur_biallelic_snps=0, family_snps=0,
                  family_monomorphic=0, family_maf_below_1_percent=0)
    genetic_score = np.zeros(visible.shape)
    unrelated_score = np.zeros(len(unrelated_index))
    common_variants = 0
    causal_variants = 0
    with ExitStack() as stack:
        def open_output(name, mode):
            return stack.enter_context((destination / name).open(mode))

        def open_gzip(name):
            return stack.enter_context(gzip.open(destination / name, "wb", compresslevel=1))

        vcf = open_gzip("chr_22.vcf.gz")
        unphased_vcf = open_gzip("unphased/chr_22.vcf.gz")
        source_eur = open_gzip("eur.chr_22.vcf.gz")
        write_vcf_header(vcf, present_ids)
        write_vcf_header(unphased_vcf, present_ids)
        write_vcf_header(source_eur, eur_ids)
        observed_bed = open_output("chr_22.bed", "wb")
        truth_bed = open_output("truth/chr_22.bed", "wb")
        for bed in (observed_bed, truth_bed):
            bed.write(b"\x6c\x1b\x01")
        bim = open_output("chr_22.bim", "w")
        weights = open_output("weights.txt", "w")
        weights.write("SNP A1 A2 b\n")

        for line in source:
            fields = line.split(b"\t", 9)
            position = int(fields[1])
            if position < args.start:
                continue
            if position > args.end:
                break
            counts["input_records"] += 1
            if fields[0] not in {b"chr22", b"22"} or fields[3] not in BASES or fields[4] not in BASES:
                continue
            eur = read_diploid_gt(fields, sample_count, eur_index)
            # The fully phased NYGC panel is expected to supply all EUR calls.
            if np.any(eur > 1):
                continue
            fields[0] = b"chr22"
            eur_count = eur.sum()
            write_vcf_genotypes(source_eur, fields, eur, PHASED_GT, eur_count)
            counts["eur_biallelic_snps"] += 1
            if args.min_family_maf and eur_count in (0, eur.size):
                counts["family_monomorphic"] += 1
                counts["family_maf_below_1_percent"] += 1
                continue
            cm = np.interp(position, cm_map[:, 0], cm_map[:, 1])
            while switch_cm <= cm:
                transmission[switch_index] ^= 1
                switch_cm, switch_index = next(switches, (np.inf, -1))
            founder = eur[founder_index]
            children = founder[parent_index, initial]
            all_calls = np.concatenate((founder.reshape(-1, 2, 2), children), axis=1)
            unrelated_calls = eur[unrelated_index]
            observed = np.concatenate((all_calls[visible], unrelated_calls))
            observed_count = observed.sum()
            observed_af = observed_count / observed.size
            maf = min(observed_af, 1 - observed_af)
            counts["family_monomorphic"] += int(maf == 0)
            counts["family_maf_below_1_percent"] += int(maf < 0.01)
            if maf < args.min_family_maf:
                continue
            counts["family_snps"] += 1
            write_vcf_genotypes(vcf, fields, observed, PHASED_GT, observed_count)
            write_vcf_genotypes(unphased_vcf, fields, observed, UNPHASED_GT, observed_count)
            write_bed_genotypes(observed_bed, observed)
            write_bed_genotypes(truth_bed, founder)
            sid, ref, alt = [value.decode() for value in fields[2:5]]
            bim.write(f"22\t{sid}\t{cm:.9f}\t{position}\t{alt}\t{ref}\n")
            if maf >= .05:
                common_variants += 1
                if (common_variants - 1) % args.causal_stride == 0:
                    beta = rng.normal()
                    genetic_score += beta * all_calls.sum(axis=2)
                    unrelated_score += beta * unrelated_calls.sum(axis=1)
                    weights.write(f"{sid} {alt} {ref} {beta:.17g}\n")
                    causal_variants += 1

    if causal_variants == 0:
        raise ValueError("No common SNPs in the selected chromosome interval")
    return counts, causal_variants, genetic_score, unrelated_score


def write_phenotype(destination, families, people, genetic_score, unrelated, unrelated_score,
                    training, genotyped_children, rng):
    direct = genetic_score[:, 2:]
    parent_mean = genetic_score[:, :2].mean(axis=1)
    # Standardize the simulated scores in the discovery cohort only.
    training_score = direct[training[:, None] & genotyped_children]
    training_mean, training_sd = training_score.mean(), training_score.std()
    direct = (direct - training_mean) / training_sd
    parent_mean = (parent_mean - parent_mean[training].mean()) / parent_mean[training].std()
    environment = rng.normal(scale=.25, size=(len(families), 1))
    residual = rng.normal(scale=.70, size=direct.shape)
    phenotype = direct + .30 * parent_mean[:, None] + environment + residual
    unrelated_y = (unrelated_score - training_mean) / training_sd + rng.normal(scale=.70, size=len(unrelated))
    with ((destination / "phenotype.txt").open("w") as full,
          (destination / "phenotype.target.txt").open("w") as target,
          (destination / "covariates.txt").open("w") as covariates):
        for out in (full, target):
            out.write("FID IID phenotype\n")
        covariates.write("FID IID covariate\n")
        for fid, siblings, values, train, present in zip(families, people[:, 2:], phenotype,
                                                        training, genotyped_children):
            for iid, value in zip(siblings[present], values[present]):
                row = f"{fid} {iid} {value:.17g}\n"
                full.write(row)
                if not train:
                    target.write(row)
                covariates.write(f"{fid} {iid} {rng.normal():.17g}\n")
        for iid, value in zip(unrelated, unrelated_y):
            full.write(f"U_{iid} {iid} {value:.17g}\n")
            covariates.write(f"U_{iid} {iid} {rng.normal():.17g}\n")


def write_marker_map(path, source_map, bim):
    positions = np.unique(np.loadtxt(bim, usecols=3, dtype=np.int64))
    cm = np.interp(positions, source_map[:, 0], source_map[:, 1])
    np.savetxt(path, np.column_stack((positions, cm)), fmt=("%d", "%.9f"),
               header="pposition gposition", comments="")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path,
                        help="NYGC 3,202-sample phased chr22 VCF.gz or derived biallelic EUR subset")
    parser.add_argument("pedigree", type=Path, help="IGSR 3,202-sample pedigree with population labels")
    parser.add_argument("map", type=Path, help="SHAPEIT5 GRCh38 chromosome 22 map (.gmap.gz)")
    parser.add_argument("output", type=Path)
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--end", type=int, default=50_818_468)
    parser.add_argument("--families", type=int, default=140)
    parser.add_argument("--children", type=int, default=3)
    parser.add_argument("--unrelated", type=int, default=20)
    parser.add_argument("--singleton-duos", action="store_true",
                        help="Mask all but one child in half of the paternal and maternal duos")
    parser.add_argument("--min-family-maf", type=float, default=.01,
                        help="Filter simulated families only. The EUR VCF retains all biallelic SNPs.")
    parser.add_argument("--causal-stride", type=int, default=200,
                        help="Select every kth SNP with sample MAF of at least .05 for the phenotype.")
    parser.add_argument("--seed", type=int, default=20261127)
    args = parser.parse_args()
    if not 1 <= args.start < args.end or args.families < 6 or args.children < 2 or args.causal_stride < 1:
        parser.error("Provide a valid interval, at least six families and two children, and a positive stride")
    if args.unrelated < 0:
        parser.error("The number of unrelated samples must be nonnegative")
    if not 0 <= args.min_family_maf <= .5:
        parser.error("The minimum MAF must be between zero and 0.5")

    for directory in (args.output, args.output / "unphased", args.output / "truth"):
        directory.mkdir(parents=True, exist_ok=True)
    cm_map = np.loadtxt(args.map, skiprows=1, usecols=(0, 2))
    rng = np.random.default_rng(args.seed)
    with gzip.open(args.source, "rb") as source:
        samples = vcf_header(source)
        founders, unrelated, eur, trios, sibships = pedigrees(args.pedigree, set(samples), args.families,
                                                              args.unrelated, rng)
        eur_ids = [sample for sample in samples if sample in eur]
        eur_index = np.flatnonzero(np.isin(samples, eur_ids))
        sample_index = {iid: index for index, iid in enumerate(eur_ids)}
        founder_index = np.array([sample_index[iid] for iid in founders.ravel()])
        unrelated_index = np.array([sample_index[iid] for iid in unrelated], dtype=int)
        visible = np.ones((args.families, args.children + 2), dtype=bool)
        visible[1::3, 0] = False
        visible[4::6, :2] = (True, False)
        visible[2::3, :2] = False
        parent_status = visible[:, :2]
        one_parent = parent_status.sum(axis=1) == 1
        if args.singleton_duos:
            for parent in parent_status.T:
                visible[np.flatnonzero(one_parent & parent)[::2], 3:] = False
        families, ids, training = pedigree_files(args.output, founders, visible, unrelated, rng)
        first_cm, last_cm = np.interp([args.start, args.end], cm_map[:, 0], cm_map[:, 1])
        initial, switches, events = transmissions(args.output / "truth" / "transmissions.tsv", families, ids,
                                                   first_cm, last_cm, rng)
        write_relationships(args.output / "relationships.seg", ids, visible, initial, events, first_cm, last_cm)
        present_ids = np.concatenate((ids[visible], unrelated))
        counts, causal, score, unrelated_score = simulate_variants(
            args, source, len(samples), eur_ids, eur_index, founder_index, unrelated_index,
            cm_map, visible, present_ids, initial, switches, rng,
        )
    bim = args.output / "chr_22.bim"
    (args.output / "truth" / "chr_22.bim").write_bytes(bim.read_bytes())
    write_marker_map(args.output / "map.txt", cm_map, bim)
    genotyped_children = visible[:, 2:]
    write_phenotype(args.output, families, ids, score, unrelated, unrelated_score,
                    training, genotyped_children, rng)
    manifest = dict(chromosome=22, start=args.start, end=args.end, seed=args.seed, eur_samples=len(eur),
                    eur_complete_trios=trios, eur_full_sibships=sibships, founder_samples=founders.size,
                    simulated_families=len(visible), observed_samples=len(present_ids),
                    unrelated_samples=len(unrelated), phenotyped_samples=int(genotyped_children.sum()) + len(unrelated),
                    causal_snps=causal,
                    discovery_families=int(training.sum()), evaluation_families=int((~training).sum()),
                    families_with_both_parents=int(visible[:, :2].all(axis=1).sum()),
                    families_with_one_parent=int((visible[:, :2].sum(axis=1) == 1).sum()),
                    singleton_duos=int((one_parent & (genotyped_children.sum(axis=1) == 1)).sum()),
                    families_without_parents=int((~visible[:, :2].any(axis=1)).sum()),
                    min_family_maf=args.min_family_maf, **counts)
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
