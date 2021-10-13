# Simulate population genotypes
python example/simulate_pop_simple.py 1000 0.5 3000 1000 1000 0.5 example/sample --chrom 1
# Convert to bgen
/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7 -g example/sample.haps -s example/sample.sample -og example/sample.bgen -os example/sample.sample -filetype shapeit_haplotypes -ofiletype bgen
# Convert to bed
plink/plink2 --bgen example/sample.bgen ref-last --sample example/sample.sample --make-bed --out example/sample --oxford-single-chr 1
# Simulate phenotype
python example/simulate_trait_quad.py example/sample.bed example/sample.ped 0.8 example/h2_quad_0.8 --no_sib --dncor 0.5
# Remove some parents
/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7 -g example/sample.bgen -s example/sample.sample -og example/sample_reduced.bgen -os example/sample_reduced.sample -filetype bgen -ofiletype bgen -excl-samples example/sample_remove.txt
# Convert reduced to bed
plink/plink2 --bgen example/sample_reduced.bgen ref-last --sample example/sample_reduced.sample --make-bed --out example/sample_reduced --oxford-single-chr 1
# Clean up
rm example/sample.haps
mv example/sample_reduced.bgen example/sample.bgen
mv example/sample_reduced.sample example/sample.sample
mv example/sample_reduced.bed example/sample.bed
mv example/sample_reduced.fam example/sample.fam
mv example/sample_reduced.bim example/sample.bim
rm example/sample_remove.txt
rm example/sample.log
rm example/sample_reduced.log