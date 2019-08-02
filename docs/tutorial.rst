Tutorial
********
Tutorial on performing robust GWAS using family data

To run this tutorial, you need plink installed and in your system path. (See http://zzz.bwh.harvard.edu/plink/download.shtml).

You also need R installed and in your system path. (See https://www.r-project.org/).

In the installation folder for sibreg, there is an examples/ folder. Change to the examples/ folder.

In this folder, there is a script 'simulate_pop.py'. At the command line, type

    ``python simulate_pop.py 5000 0.5 12000 4000 4000 0.5 sim``

This will produce a simulated population of 12,000 families genotyped at 5,000
independent SNPs with allele frequency 0.5. The genotypes will be in sim.ped (plink ped format)
and the pedigree will be in sim_fam.ped.

We need to convert this to plink .bed format (binary ped) for further use. To do this, type

    ``plink --file sim --make-bed --out sim``

To simulate the test phenotype, at the command line, type

    ``python simulate_trait_quad.py sim.bed sim_fams.ped 0.5 h2_trio_0.5 --no_sib``

This simulates a trait with direct, paternal, and maternal effects, where 50% of the phenotypic
variance is explained by the combined direct, paternal and maternal effects of the 5,000 SNPs. In
the following tutorial, we do not fit indirect genetic effects from siblings to keep it simple.

To simulate missing genotypes of parents/siblings, type the following command:

    ``plink --bfile sim --remove sim_remove.txt --make-bed --out sim_reduced``

This removes the genotypes of both parents for 4,000 families, the genotype of one parent
for 4,000 other families, and randomly removes siblings with probability 0.5 from
the families with at least one genotyped parent. The resulting reduced bed file is 'sim_reduced.bed'

We note that the following commands can take some time to run (10s of minutes depending on hardware).

To impute the missing parental genotypes for the families with one parent genotyped, type:

    ``python ../bin/impute_po.py sim_reduced.bed sim_fams.ped po_impute``

To impute the sum of the missing parents' genotypes for the families without genotyped parents, type:

    ``python ../bin/impute_from_sibs_hdf5.py sim.hdf5 sim_reduced.bed sim_fams.ped sib_impute``

To estimate effects for the families with both parents genotyped, type:

    ``python ../bin/triGWAS.py sim_reduced.bed sim_fams.ped h2_trio_0.5.ped triGWAS``

To estimate effects for the families with only one parent genotyped, type:

    ``python ../bin/poGWAS.py sim_reduced.bed po_impute.hdf5 sim_fams.ped h2_trio_0.5.ped poGWAS``

To estimate effects for the families without genotyped parents, type:

    ``python ../bin/pGWAS.py sim_reduced.bed sib_impute.hdf5 sim_fams.ped h2_trio_0.5.ped pGWAS --no_sib``

The '--no_sib' option stops the script from fitting indirect effects from siblings (the default behaviour).

Now we have estimated effects from the three different subsets of data (no missing parents, one missing parent,
both parents missing). To meta-analyse the mutlidimensional parameter vector (direct, paternal, and maternal),
type:

    ``Rscript estimate_and_meta_analyse.R pGWAS.hdf5 poGWAS.hdf5 triGWAS.hdf5 h2_trio_0.5.effects.txt h2_trio_0.5.estimates TRUE``

This should print estimates of the bias of the effect estimates, with output something like this:

    ``[1] "bias for direct effects: -0.021 (0.0426 S.E.)"``

    ``[1] "bias for paternal effects: -0.008 (0.0509 S.E.)"``

    ``[1] "bias for maternal effects: -0.049 (0.0511 S.E.)"``

If everything has worked, the bias should not be statistically significant from zero (with high probability).

The meta-analysis estimates along with their sampling covariance matrices and standard errors are output in h2_trio_0.5.estimates.hdf5.