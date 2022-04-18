# Simulate population genotypes
python ../snipar/simulate/simulate_pop_simple.py 1000 0.5 10000 0 10000 0 chr_1 --chrom 1 --blocksize 1
# Convert to bgen
qctool_v2.2.0-osx/qctool -g chr_1.haps -s chr_1.sample -og chr_1.bgen -os chr_1.sample -filetype shapeit_haplotypes -ofiletype bgen
# Convert to bed
./plink2 --bgen chr_1.bgen ref-last --sample chr_1.sample --make-bed --out chr_1 --oxford-single-chr 1
# Simulate phenotype
python ../snipar/simulate/simulate_trait_quad.py chr_1.bed chr_1.ped 0.8 phenotype --no_sib --dncor 0.5
# Remove some parents
qctool_v2.2.0-osx/qctool -g chr_1.bgen -s chr_1.sample -og chr_1_reduced.bgen -os chr_1_reduced.sample -filetype bgen -ofiletype bgen -excl-samples chr_1_remove.txt
# Convert reduced to bed
./plink2 --bgen chr_1_reduced.bgen ref-last --sample chr_1_reduced.sample --make-bed --out chr_1_reduced --oxford-single-chr 1
# Clean up
rm chr_1.haps
mv chr_1_reduced.bgen chr_1.bgen
mv chr_1_reduced.sample chr_1.sample
mv chr_1_reduced.bed chr_1.bed
mv chr_1_reduced.fam chr_1.fam
mv chr_1_reduced.bim chr_1.bim
rm chr_1_remove.txt
rm chr_1.log
rm chr_1_reduced.log
mv phenotype.direct_weights.txt direct_weights.txt
mv chr_1.genetic_map.txt genetic_map.txt
mv chr_1.king.kin0 king.kin0
mv chr_1.king.segments.gz king.segments.gz
mv chr_1.chr_1.ibd.segments.gz chr_1.ibd.segments.gz
mv chr_1.kingallsegs.txt kingallsegs.txt
mv chr_1.ped pedigree.txt
mv chr_1.agesex agesex.txt

## Analyse
impute.py --ibd chr_1.ibd --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 
#impute.py --ibd king --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 --ibd_is_king
gwas.py phenotype.txt --bed chr_@ --imp chr_@ --threads 4
python ../snipar/example/example_data/estimate_sim_effects.py chr_1.sumstats.hdf5 phenotype.effects.txt
pgs.py direct --bed chr_@ --imp chr_@ --weights direct_weights.txt
pgs.py direct --pgs direct.pgs.txt --phenofile phenotype.txt
pgs.py direct_sib --bed chr_1 --imp chr_1 --weights direct_weights.txt --phenofile phenotype.txt --fit_sib
## Clean
#rm *txt king* chr_1* 