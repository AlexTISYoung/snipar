========
Tutorial
========
Tutorial on imputing missing parental genotypes and performing family based GWAS and polygenic score analyses. Before working through the tutorial, please first install the package and run the tests (see `github <https://github.com/AlexTISYoung/SNIPar>`_).

Generating test data
--------------------

To generate the test data, in the main SNIPar directory, run:

    ``bash tests/generate_test_population.sh``

In the test_data/ directory, the file h2_quad_0.8.ped is a simulated trait with direct, paternal, and maternal effects, where 80% of the phenotypic
variance is explained by the combined direct, paternal and maternal effects of the SNPs; and the
pairwise correlations between the direct, paternal, and maternal effects is 0.5. The phenotype file is test_data/h2_quad_0.8.ped.

The genotype data has been simulated so that there are 3000 independent families, where 1000 have two siblings but no parents genotyped,
1000 have one parent genotyped and a 50% chance of having a genotyped sibling, and the final 1000 have both parents genotyped and a 50%
chance of having a genotyped sibling.

Imputing missing parental genotypes
-----------------------------------

To impute the missing parental genotypes, type:

    ``python impute_runner.py test_data/sample.segments.gz --bed test_data/sample1 --king test_data/sample.king --agesex test_data/sample.agesex --output_address test_data/sample1 --threads 4``

The script constructs a pedigree from the output of KING's relatedness inference (test_data/sample.king),
and age and sex information (test_data/sample.agesex). The pedigree along with the IBD segments shared between siblings recorded in test_data/sample.segments.gz are used to impute missing parental genotypes
from the sibling and observed parental genotypes in test_data/sample1.bed. The imputed parental genotypes are in a HDF5 file test_data/sample1.hdf5. The --threads 4 argument
means the imputation will run on 4 threads.

If phased haplotypes are available in .bgen format, the imputation can use these as input, which improves the information gained by imputation
in certain situations. To perform imputation from the phased .bgen file in test_data/, use the following command:

    ``python impute_runner.py test_data/sample.segments.gz --bgen test_data/sample1 --king test_data/sample.king --agesex test_data/sample.agesex --output_address test_data/sample1 --threads 4 --from_chr 1 --to_chr 2``

It is necessary to provide the --from_chr and --to_chr arguments when imputing from .bgen files since they often do not contain information on which chromosome
the SNPs are located on, and we need to match up the IBD segments to the SNPs on the same chromosome.


Family based GWAS
-----------------

To compute summary statistics for direct, paternal, and maternal effects for all SNPs in the .bed file, type:

    ``python fGWAS.py test_data/sample1 test_data/h2_quad_0.8.ped --outprefix test_data/h2_quad_0.8 --bed test_data/sample1``

This takes the observed genotypes in test_data/sample1.bed and the imputed parental genotypes in test_data/sample1.hdf5 and uses
them to perform, for each SNP, a joint regression onto the proband's genotype, the father's (imputed) genotype, and the mother's
(imputed) genotype. This is done using a random effects model that models phenotypic correlations between siblings,
where sibling relations are inferred from the pedigree stored in the output of the imputation script: test_data/sample1.hdf5. The 'family variance estimate'
output is the  phenotypic variance explained by mean differences between sibships, and the residual variance is the remaining phenotypic variance.

To use the .bgen file instead, type:

    ``python fGWAS.py test_data/sample1 test_data/h2_quad_0.8.ped --outprefix test_data/h2_quad_0.8 --bgen test_data/sample1``

The script outputs summary statistics in a gzipped text file: h2_quad_0.8.sumstats.gz. This file gives the chromosome,
SNP id, position, alleles (A1, the allele that effects are given with respect to; and A2, the alternative allele),
the frequency of the A1 allele, then summary statistics for each type of effect. For each effect, we give the
effective N for each SNP; this differs from the actual N due to the fact that there are differing amounts of information
for each type of effect, and due to relatedness in the sample. We give the effect estimate in the first column for each effect, the column
'effect_Beta', where 'effect' can be direct, paternal, etc; this is followed by the standard error, the Z-score,
and the negative log10 P-value for a non-zero effect. In addition to effects directly estimated by the script,
we also output the average parental effect estimate (estimate of the average of maternal and paternal effects),
and the population effect estimate, which is equivalent to what is estimated by standard GWAS methods that
regress phenotype onto genotype without control for parental genotypes. The final columns give the sampling
correlations between the different effect estimates at that SNP.

In addition to the plain text output, the effects and their sampling variance-covariance matrices are output in test_data/h2_quad_0.8.sumstats.hdf5. The contents of HDF5 file can be read into Python (using `h5py <https://www.h5py.org>`_) and R (using `rhdf5 <https://www.bioconductor.org/packages/release/bioc/html/rhdf5.html>`_) easily.
The output contains different datasets:

1. *estimate*, the estimated SNP effect, where each row gives a SNP, and each column gives an effect
2. *bim*, equivalent to the bim file for plink, recording the information on each SNP
3. *estimate_cols*, gives the names of the effects estimate for each SNP: direct, paternal, maternal, etc.
4. *estimate_ses*, the standard errors for the effect estimates in *estimate*
5. *estimate_covariance*, 3 dimensional array with sampling variance-covariance matrices for each SNP's estimated effects, with SNPs indexed by the first axis
6. *freqs*, frequencies of the effect alleles
7. *sigma2*, maximum likelihood estimate of the residual variance in the null model
8. *tau*, maximum likelihood estimate of the ratio between the residual variance and family variance
9. *N*, the sample size
10. *NAs*, the number of missing values for each of SNPs, given for each relative in the regression (individual, father, mother, etc.)


Now we have estimated locus specific summary statistics. To compare to the true effects, run

    ``python example/estimate_sim_effects.py test_data/h2_quad_0.8.sumstats.hdf5 test_data/h2_quad_0.8.effects.txt``

This should print estimates of the bias of the effect estimates, with output something like this:

    ``Bias for direct genetic effects: -0.0582; S.E. 0.0375``

    ``Bias for paternal  effects: 0.0308; S.E. 0.0614``

    ``Bias for maternal effects: 0.058; S.E. 0.0597``

    ``Bias for average parental effects: 0.0177; S.E. 0.049``

    ``Bias for population effects as estimates of direct effects: 0.4561; S.E. 0.0413``

The bias estimates for direct, paternal, maternal, and average parental effects should not be statistically significantly different from zero (with high probability). Population effects (which are estimated by univariate regression of individuals' phenotypes onto their genotypes -- as in standard GWAS)
here are biased estimates of direct effects, since population effects include both direct and indirect parental effects.

If the imputation has been performed from siblings alone, then the regression onto proband (focal, phenotyped individual), imputed paternal, and imputed maternal becomes
co-linear. This is because the imputation is the same for paternal and maternal genotypes. In this case, the regression should be performed
onto proband and sum of imputed paternal and maternal genotypes. This can be achieved by providing the *--parsum* option to the script. 
The script can also estimate indirect sibling effects for each SNP by providing the *--fit_sib* option; however, this
will reduce power for estimating other effects.


Polygenic score analyses
------------------------

In addition to family based GWAS, SNIPar provides a script (fPGS.py) for computing polygenic scores (PGS) based on observed/imputed genotypes,
and for performing family based polygenic score analyses. Here, we give some examples of how to use this script. The script computes a PGS
from weights provided in `LD-pred <https://github.com/bvilhjal/ldpred>`_ format . The true direct genetic effects for the simulated trait are given as PGS weights in this format
in test_data/h2_quad_0.8.direct_weights.txt. This is a tab-delimited text file with a header and columns 'chrom' (chromosome), 'pos' (position), 'sid' (SNP ID), 'nt1' (allele 1),
'nt2' (allele 2), 'raw_beta' (raw effect estimates), 'ldpred_beta' (LD-pred adjusted weight). The script uses as weights the 'ldpred_beta' column.

To compute the PGS from the true direct effects, use the following command:

    ``python fPGS.py test_data/direct --bedfiles test_data/sample1 --impfiles test_data/sample1 --weights test_data/h2_quad_0.8.direct_weights.txt``
    
This uses the weights in the weights file to compute the polygenic scores for each genotyped individual for whom observed or imputed parental genotypes are available.
It outputs the PGS to test_data/direct.pgs.txt, which is a white-space delimited text file with columns FID (family ID, shared between siblings), IID (individual ID),
proband (PGS of individual with given IID), maternal (observed or imputed PGS of that individual's mother), paternal (observed or imputed PGS of that individual's father).
The script also supports bed files and imputed files split by chromosome. If you had bed files as chr_1.bed, chr_2.bed, ..., chr_22.bed; and imputed parental genotype files
as chr_1.hdf5, chr_2.hdf5, ..., chr_22.hdf5, then you can specify this in a command as:

    ``--bedfiles chr_~ --impfiles chr_~``

The script looks for all files that match the path given with '~' replaced by 1,2,...,22: chr_1.bed & chr_1.hdf5, chr_2.bed & chr_2.hdf5, etc.
To use .bgen input, replace the --bedfiles argument with --bgenfiles.

To estimate direct, paternal, and maternal effects of the PGS, use the following command:

    ``python fPGS.py test_data/direct --pgs test_data/direct.pgs.txt --phenofile test_data/h2_quad_0.8.ped``

This uses a linear mixed model that has a random effect for mean differences between families (defined as sibships here) and fixed effects for the direct,
paternal, and maternal effects of the PGS. It also estimates the 'population' effect of the PGS: the effect from regression of individuals' phenotypes onto their PGS values.
The estimated effects and their standard errors are output to test_data/direct.pgs_effects.txt, with the effect names (direct, paternal, maternal, population) in the first column,
their estimates in the second column, and their standard errors in the final column. The sampling variance-covariance matrix of direct, paternal, and maternal effects is output in test_data/direct.pgs_vcov.txt.

Estimates of the direct effect of the PGS should be equal to 1 in expectation since
we are using the true direct effects as the weights, so the PGS corresponds to the true direct effect component of the trait.
The parental effect estimates capture the correlation between the direct and indirect parental effects. The population effect estimate
should be greater than 1, since this captures both the direct effect of the PGS, and the correlation between direct and indirect parental effects.

If parental genotypes have been imputed from sibling data alone, then imputed paternal and maternal PGS are perfectly correlated, and the above regression on proband, paternal, and maternal
PGS becomes co-linear. To deal with this, add the --parsum option to the above command, which will estimate the average parental effect rather than separate maternal and paternal effects of the PGS.

It is also possible to estimate indirect effects from siblings. We can compute the PGS for genotyped individuals with genotyped siblings and estimate direct, indirect sibling, paternal and maternal effects in
one command with the addition of the --fit_sib option:

   ``python fPGS.py test_data/direct_sib --bedfiles test_data/sample1 --impfiles test_data/sample1 --weights test_data/h2_quad_0.8.direct_weights.txt --phenofile test_data/h2_quad_0.8.ped --fit_sib``

This outputs the PGS values for each individual along with the PGS value of their sibling, and imputed/observed paternal and maternal PGS to test_data/direct_sib.pgs.txt.
(If an individual has multiple genotyped siblings, the average of the siblings' PGS is used for the PGS of the sibling.)
It outputs estimates of direct, indirect sibling, paternal, and maternal effects of the PGS to test_data/direct_sib.pgs_effects.txt and their sampling variance-covariance matrix to test_data/direct_sib.pgs_vcov.txt.
Since indirect effects from siblings were zero in this simulation, the estimated sibling effect should be close to zero.

Note that the standard error for the direct effect estimate increases: this is due both to a drop in sample size since only those probands with genotyped siblings are included, and due to the fact that adding the sibling effect to the regression
decreases the independent information on the direct effect.
