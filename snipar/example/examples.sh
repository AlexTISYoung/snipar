cd example_data

ibd.py --bgen chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 --ld_out --map genetic_map.txt
impute.py --ibd chr_@.ibd --bgen chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4

# robust
gwas.py phenotype.txt --bgen chr_@ --imp chr_@ --cpus 1 --out robust_@ --no_grm_var --robust

# unified
gwas.py phenotype.txt --bgen chr_trios_sibs_singletons_@ --imp chr_@ --cpus 1 --out unified_@ --no_grm_var --impute_unrel

# sib-difference
gwas.py phenotype.txt --bgen chr_trios_sibs_singletons_@ --imp chr_@ --cpus 1 --out sib_diff_@ --no_grm_var --sib_diff

# Young
gwas.py phenotype.txt --bgen chr_trios_sibs_singletons_@ --imp chr_@ --cpus 1 --out sib_diff_@ --no_grm_var
