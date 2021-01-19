Tutorial
********
Tutorial on performing robust GWAS using family data. Before working through the tutorial, please first install the package and run the tests.

To generate the test data, in the main SNIPar directory, run:

    ``bash tests/generate_test_population.sh``

In the test_data/ directory, the file h2_quad_0.8.ped is a simulated trait with direct, paternal, and maternal effects, where 80% of the phenotypic
variance is explained by the combined direct, paternal and maternal effects of the SNPs; and the
pairwise correlations between the direct, paternal, and maternal effects is 0.5. The phenotype file is test_data/h2_quad_0.8.ped.

The genotype data has been simulated so that there are 3000 independent families, where 1000 have two siblings but no parents genotyped,
1000 have one parent genotyped and a 50% chance of having a genotyped sibling, and the final 1000 have both parents genotyped and a 50%
chance of having a genotyped sibling.

To impute the missing parental genotypes, run the following command in the main SNIPar directory:

python impute_runner.py 1 2 test_data/sample.segments.gz test_data/sample --king test_data/sample.king --agesex test_data/sample.agesex --out_prefix test_data/sample

python fGWAS.py test_data/sample1.bed test_data/sample1.hdf5 test_data/h2_quad_0.8.ped test_data/h2_quad_0.8

To impute the missing parental genotypes, type:

    ``python impute_runner.py test_data/sample.segments.gz test_data/sample1 --king test_data/sample.king --agesex test_data/sample.agesex --output_address test_data/sample1 --threads 4``

The script constructs a pedigree from the output of KING's relatedness inference (test_data/sample.king),
and age and sex information (test_data/sample.agesex). The pedigree along with the IBD segments shared between siblings recorded in test_data/sample.segments.gz are used to impute missing parental genotypes
from the sibling and observed parental genotypes in test_data/sample1.bed. The imputed parental genotypes are in a HDF5 file test_data/sample1.hdf5. The --threads 4 argument
means the imputation will run on 4 threads. If imputing for more than 1 chromosome, the --processes n argument will use n different processes in parallel, one for
each chromosome, with the number of threads per process determined by the --threads argument.

To compute summary statistics for direct, paternal, and maternal effects for all SNPs in the .bed file, type:

    ``python fGWAS.py test_data/sample1 test_data/sample1.hdf5 test_data/h2_quad_0.8.ped test_data/h2_quad_0.8``

This takes the observed genotypes in test_data/sample1.bed and the imputed parental genotypes in test_data/sample1.hdf5 and uses
them to perform, for each SNP, a joint regression onto the proband's genotype, the father's (imputed) genotype, and the mother's
(imputed) genotype. This is done using a random effects model that models phenotypic correlations between siblings,
where sibling relations are inferred from the pedigree stored in the output of the imputation script: test_data/sample1.hdf5. The 'family variance estimate'
output is the  phenotypic variance explained by mean differences between sibships, and the residual variance is the remaining phenotypic variance.
The effects are output in test_data/h2_quad_0.8.hdf5 and are given with respect to the first allele in the .bim file.

Now we have estimated locus specific summary statistics. To estimate effects and compare to the true effects, run

    ``Rscript example/estimate_sim_effects.R test_data/h2_quad_0.8.hdf5 test_data/h2_quad_0.8.effects.txt test_data/h2_quad_0.8.estimates``

This should print estimates of the bias of the effect estimates, with output something like this:

    ``[1] "bias for direct effects: 0.0281 (0.0374 S.E.)"``
    ``[1] "bias for paternal effects: -0.0376 (0.0524 S.E.)"``
    ``[1] "bias for maternal effects: 0.0528 (0.05 S.E.)"``
    ``[1] "Chi-square test p-value: 0.4879"``

If everything has worked, the bias should not be statistically significantly different from zero (with high probability).

The Chi-Square test p-value should also not be statistically significant with high probability. (This tests that the sampling distribution
of the estimates is correct, as well as the estimates being unbiased.)

The estimates along with their sampling covariance matrices and standard errors are output in test_data/h2_quad_0.8.estimates.hdf5.

If the imputation has been performed from siblings alone, then the regression onto proband, imputed paternal, and imputed maternal becomes
co-linear. This is because the imputation is the same for paternal and maternal genotypes. In this case, the regression should be performed
onto proband and sum of imputed paternal and maternal genotypes. This can be achieved by providing the --parsum option to the script:

    ``python fGWAS.py test_data/sample1 test_data/sample1.hdf5 test_data/h2_quad_0.8.ped test_data/h2_quad_0.8_parsum --parsum``

This will output estimates of direct and average parental (average of maternal and paternal) effects, along with sampling covariance
matrices and standard errors.

In addition to family based GWAS, SNIPar provides a script (fPGS.py) for computing polygenic scores (PGS) based on observed/imputed genotypes,
and for performing family based polygenic score analyses. Here, we give an example of how to use this script. The script computes a PGS
from weights provided in LD-pred format. The true direct genetic effects for the simulated trait are given as PGS weights in this format
in test_data/h2_quad_0.8.direct_weights.txt. This is a tab-delimited text file with a header and columns 'chrom' (chromosome), 'pos' (position), 'sid' (SNP ID), 'nt1' (allele 1),
'nt2' (allele 2), 'raw_beta' (raw effect estimates), 'ldpred_beta' (LD-pred adjusted weight). The script uses as weights the 'ldpred_beta' column.

To compute the PGS from the true direct effect weights, use the following command:

    ``python fPGS.py test_data/direct --bedfiles test_data/sample1.bed --impfiles test_data/sample1.hdf5 --weights test_data/h2_quad_0.8.direct_weights.txt``

This uses the weights in the weights file to compute the polygenic scores for each genotyped individual for whom observed or imputed genotypes are available
for both parents. It outputs the PGS to test_data/direct.pgs.txt, which is a white-space delimited text file with columns FID (family ID, shared between siblings), IID (individual ID),
proband (PGS of individual with given IID), maternal (observed or imputed PGS of that individual's mother), paternal (observed or imputed PGS of that individual's father).

To estimate direct, paternal, and maternal effects of the PGS, use the following command:

    ``python fPGS.py test_data/direct --pgs test_data/direct.pgs.txt --phenofile test_data/h2_quad_0.8.ped``

This uses a linear mixed model that has a random effect for mean differences between families (defined as sibships here) and fixed effects for the direct,
paternal, and maternal effects of the PGS. It also estimates the 'population' effect of the PGS: the effect from regression of individual's phenotypes onto their PGS values.
The estimated effects and their standard errors are output to test_data/direct.pgs_effects.txt, with the effect names (direct, paternal, maternal, population) in the first column,
their estimates in the second column, and their standard errors in the final column. The sampling variance-covariance matrix of direct, paternal, and maternal effects is output in test_data/direct.pgs_vcov.txt.
Estimates of direct effect of the PGS should be equal to 1 in expectation since
we are using the true direct effects as the weights, so the PGS corresponds to the true direct effect component of the trait. The estimated direct effect here should be within 2 standard errors
of 1 approximately 95\% of the time. The parental effect estimates capture the correlation between the direct and indirect parental effects. The population effect estimate
should be greater than 1, since this captures both the direct effect of the PGS, and the correlation between direct and indirect parental effects.

If parental genotypes have been imputed from sibling data alone, then imputed paternal and maternal genotypes are perfectly correlated, and the above regression on proband, paternal, and maternal
PGS becomes co-linear. To deal with this, add the --parsum option to the above command, which will estimate the average parental effect rather than separate maternal and paternal effects of the PGS:

   ``python fPGS.py test_data/direct_avg_parental --pgs test_data/direct.pgs.txt --phenofile test_data/h2_quad_0.8.ped --parsum``

This outputs estimates of direct and average parental effects to test_data/direct_avg_parental.pgs_effects.txt, and their sampling variance-covariance matrix to test_data/direct_avg_parental.pgs_vcov.txt.

It is also possible to estimate indirect effects from siblings. We can compute the PGS for genotyped individuals with genotyped siblings and estimate direct, indirect sibling, paternal and maternal effects in
one command with the addition of the --fit_sib option:

   ``python fPGS.py test_data/direct_sib --bedfiles test_data/sample1.bed --impfiles test_data/sample1.hdf5 --weights test_data/h2_quad_0.8.direct_weights.txt --phenofile test_data/h2_quad_0.8.ped --fit_sib``

This outputs estimates of direct, indirect sibling, paternal, and maternal effects of the PGS to test_data/direct_sib.pgs_effects.txt and their sampling variance-covariance matrix to test_data/direct_sib.pgs_vcov.txt.
Since indirect effects from siblings were zero in this simulation, the estimated sibling effect should be within 2 standard errors of zero approximately 95% of the time. Note that the standard error for the direct
effect estimate increases: this is due both to a drop in sample size since only those probands with genotyped siblings are included, and due to the fact that adding the sibling effect to the regression
decreases the independent information on the direct effect.