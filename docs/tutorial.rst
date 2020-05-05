Tutorial
********
Tutorial on performing robust GWAS using family data

In the test_data/ directory, the file h2_quad_0.8.ped is a simulated trait with direct, paternal, and maternal effects, where 80% of the phenotypic
variance is explained by the combined direct, paternal and maternal effects of the SNPs; and the
pairwise correlations between the direct, paternal, and maternal effects is 0.5. The phenotype file is test_data/h2_quad_0.8.ped.

The genotype data has been simulated so that there are 1200 independent families, where 400 have two siblings but no parents genotyped,
400 have one parent genotyped and a 50% chance of having a genotyped sibling, and the final 400 have both parents genotyped and a 50%
chance of having a genotyped sibling.

To impute the missing parental genotypes, run the following command in the main sibreg directory:

python impute_runner.py 1 2 test_data/sample.segments.gz test_data/sample --king test_data/sample.king --agesex test_data/sample.agesex --out_prefix test_data/sample

python fGWAS.py test_data/sample1.bed test_data/sample1.hdf5 test_data/h2_quad_0.8.ped test_data/h2_quad_0.8

To impute the missing parental genotypes, type:

    ``python impute_runner.py 1 2 test_data/sample.segments.gz test_data/sample --king test_data/sample.king --agesex test_data/sample.agesex --out_prefix test_data/sample``

The script constructs a pedigree from the output of KING's relatedness inference (test_data/sample.king),
and age and sex information (test_data/sample.agesex). The pedigree along with the IBD segments shared between siblings recorded in test_data/sample.segments.gz are used to impute missing parental genotypes.
The imputed parental genotypes are in a HDF5 file test_data/sample1.hdf5.

To compute summary statistics for direct, paternal, and maternal effects for all SNPs in the .bed file, type:

    ``python fGWAS.py test_data/sample1.bed test_data/sample1.hdf5 test_data/h2_quad_0.8.ped test_data/h2_quad_0.8``

This takes the observed genotypes in test_data/sample1.bed and the imputed parental genotypes in test_data/sample1.hdf5 and uses
them to perform, for each SNP, a joint regression onto the proband's genotype, the father's (imputed) genotype, and the mother's
(imputed genotype). This is done using a random effects model that models phenotypic correlations between siblings,
where sibling relations are inferred from the pedigree stored in the output of the imputation script: test_data/sample1.hdf5.

Now we have estimated locus specific summary statistics. To estimate effects and compare to the true effects, run

    ``Rscript example/estimate_sim_effects.R test_data/h2_quad_0.8.hdf5 test_data/h2_quad_0.8.effects.txt test_data/h2_quad_0.8.estimates``

This should print estimates of the bias of the effect estimates, with output something like this:

    ``[1] "bias for direct effects: -0.021 (0.0626 S.E.)"``

    ``[1] "bias for paterntal effects: -0.012 (0.0826 S.E.)"``

    ``[1] "bias for maternal effects: -0.012 (0.0926 S.E.)"``

If everything has worked, the bias should not be statistically significantly different from zero (with high probability).

The meta-analysis estimates along with their standard errors are output in test_data/h2_quad_0.8.estimates.hdf5.