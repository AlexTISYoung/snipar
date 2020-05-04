Tutorial
********
Tutorial on performing robust GWAS using family data

To run this tutorial, you need plink installed and in your system path. (See http://zzz.bwh.harvard.edu/plink/download.shtml).

You also need R installed and in your system path. (See https://www.r-project.org/).

In the installation folder for sibreg, there is an example/ folder. Change to the example/ folder.

In this folder, there is a script 'simulate_pop.py'. At the command line, type this:

    ``python simulate_pop.py 1000 0.5 3000 1000 1000 0.5 chr_1``

This will produce a simulated population of 6,000 families genotyped at 2,000
independent SNPs with allele frequency 0.5. The genotypes will be in sim.ped (plink ped format)
and the pedigree will be in sim_fams.ped.

We need to convert this to plink .bed format (binary ped) for further use. To do this, type

    ``plink --file chr_1 --make-bed --out chr_1``

To simulate the test phenotype, at the command line, type

    ``python simulate_trait_quad.py chr_1.bed chr_1_fams.ped 0.8 h2_quad_0.8 --no_sib --dncor 0.5``

This simulates a trait with direct, paternal, and maternal effects, where 80% of the phenotypic
variance is explained by the combined direct, paternal and maternal effects of the SNPs; and the
pairwise correlations between the direct, paternal, and maternal effects is 0.5. The phenotype file is h2_quad_0.8.ped.

To simulate missing genotypes of parents/siblings, type the following command:

    ``plink --bfile chr_1 --remove chr_1_remove.txt --make-bed --out reduced_chr_1''

This removes the genotypes of both parents for 2,000 families, the genotype of one parent
for 2,000 other families, and randomly removes siblings with probability 0.5 from
the families with at least one genotyped parent. The resulting reduced bed file is 'sim_reduced.bed'

We note that the following commands can take some time to run depending on hardware.

To impute the missing parental genotypes, type:

    ``python ../impute_runner.py 1 2 chr_1.segments.gz reduced_chr_ --king chr_1.king.kin0 --agesex chr_1.agesex --out_prefix chr_``

The imputed parental genotypes are in a HDF5 file sib_impute.

To estimate effects, type:

    ``python ../fGWAS.py reduced_chr_1.bed chr_1.hdf5 h2_quad_0.8.ped h2_quad_0.8``

Now we have estimated locus specific summary statistics. To estimate effects and compare to the true effects, run

    ``Rscript estimate_and_meta_analyse.R pGWAS.hdf5 poGWAS_no_sib.hdf5 poGWAS_sib.hdf5 triGWAS_no_sib.hdf5 triGWAS_sib.hdf5 h2_quad_0.5.effects.txt h2_quad_0.5.estimates FALSE``

This should print estimates of the bias of the effect estimates, with output something like this:

    ``[1] "bias for direct effects: -0.021 (0.0426 S.E.)"``

    ``[1] "bias for sib effects: -0.012 (0.0426 S.E.)"``

If everything has worked, the bias should not be statistically significantly different from zero (with high probability).

The meta-analysis estimates along with their standard errors are output in h2_quad_0.5.estimates.hdf5.