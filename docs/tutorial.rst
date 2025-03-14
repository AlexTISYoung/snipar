.. _tutorial:
========
Tutorial
========

Tutorial on inferring IBD between siblings, imputing missing parental genotypes, and performing family-based GWAS and polygenic score analyses. 
Before working through the tutorial, please first install the package and read the :ref:`guide <guide>`. 

**family-based GWAS can be performed without imputed parental genotypes**

**snipar will meta-analyse siblings and trios by default when imputed parental genotypes are not provided**

**this tutorial guides through the entire workflow but gives examples of analyses without imputed genotypes**

Test data
--------------------

If *snipar* has been installed succesfully, the :ref:`command line scripts <scripts>` should be accessible as
executables in your terminal. A script that should be accessible loads the tutorial example data into a specified directory.
To create a directory called 'example_data/' in the current directory and load the example data into it, use the command:

    ``snipar_example_data.py --dest example_data``

You can create the example data directory elsewhere by changing the --dest argument. Please change your working directory to example_data/:

    ``cd example_data``

In this directory, there is some example data. 
The file phenotype.txt is a :ref:`phenotype file <phenotype>` containing a simulated phenotype with direct, paternal, and maternal effects, where 80% of the phenotypic
variance is explained by the combined direct, paternal and maternal effects of the SNPs; and the
pairwise correlations between the direct, paternal, and maternal effects are 0.5. 

The genotype data has been simulated so that there are 3000 independent families, where 1000 have two siblings but no parents genotyped,
1000 have one parent genotyped and a 50% chance of having a genotyped sibling, and the final 1000 have both parents genotyped and a 50%
chance of having a genotyped sibling. The example data includes :ref:`observed genotype data <observed genotypes>` formatted in both PLINK .bed format (chr_1.bed) and phased genotype
data in .bgen format (chr_1.bgen with associated sample file chr_1.sample). The folder also contains partial data—for example, chr_1_trios_sibs.bgen, which includes individuals with 
complete parental genotypes as well as individuals without parents but with siblings. This allows users to apply various estimators to different types of data.

Inferring IBD between siblings
------------------------------

The first step is to infer the identity-by-descent (IBD) segments shared between siblings.
*snipar* contains a script, :ref:`ibd.py <ibd.py>`, that employs a Hidden Markov Model (HMM) to infer the IBD segments for the sibling pairs.
The per-SNP genotyping error probability will be inferred from parent-offspring pairs when available;
alternatively, a genotyping error probability can be provided using the :code:`--p_error` option. By default, SNPs with
genotyping error rates greater than 0.01 will be filtered out, but this threshold can be changed with the :code:`--max_error` argument.
To infer the IBD segments from the genotype data in chr_1.bed,use the following command

    ``ibd.py --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 --ld_out``

This will output the IBD segments to a :ref:`gzipped text file <ibd_segments_file>` chr_1.ibd.segments.gz. 
Genotype files split over multiple chromosomes can be specified
using '@' as a numerical wildcard character: see :ref:`here <multichrom>`. 
In this example, :code:`--bed chr_@` instructs ibd.py to search for .bed files
chr_1.bed, chr_2.bed, ..., chr_22.bed, where each bed file contains SNPs from the numbered chromosome. 
In this case, only one bed file is in example_data/, chr_1.bed. 
If bed files for multiple chromosomes are found, IBD will be inferred separately for each chromosome, with one
output file per chromosome, with the chromosome number filling in the numerical wildcard in the --out argument. 
Alternatively, you can specify the chromosomes that you want to analyse: for example, you can include :code:`--chr_range 1-10`
for chromosome 1-10, or :code:`--chr_range 1 3 5-22` for chromosome 1, 3, and 5-22.

The :code:`--king` argument requires the address of the :ref:`KING kinship file <kinship>`, 
and the :code:`--agesex` argument requires the address of the :ref:`agesex file <agesex>`.
Age and sex information is needed to determine which is the mother/father in a parent-offspring relation.

The algorithm requires a genetic map to compute the probabilities of transitioning between different IBD states. 
If the genetic map positions (in cM) are provided in the .bim file, the script will use these. 
Alternatively, the ``--map`` argument allows the user to specify a genetic map in the same format as used by SHAPEIT 
(https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#formats) an example of which is 
provided in genetic_map.txt. 

If no genetic map is provided, then the deCODE sex-averaged map on GRCh38 coordinates (Halldorsson, Bjarni V., et al. "Characterizing mutagenic effects of recombination through a sequence-level genetic map." Science 363.6425 (2019).),
which is distributed as part of *snipar*, will be used. 

The algorithm computes LD scores of SNPs in order to account for correlations between SNPs. 
The ``--ld_out`` argument writes the LD scores to file in the same format as LDSC (https://github.com/bulik/ldsc). 

The user can also input a phased .bgen file. For example, to infer IBD from chr_1.bgen using the genetic map in genetic_map.txt, use this command:

    ``ibd.py --bgen chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 --ld_out --map genetic_map.txt``

If the user has a :ref:`pedigree file <pedigree>`, they can input that instead of the *--king* and *--agesex* arguments. 
Siblings are inferred as individuals in the pedigree that share both parents. 
Using the example pedigree in pedigree.txt, you can infer IBD using this command:

    ``ibd.py --bed chr_@ --pedigree pedigree.txt --map genetic_map.txt --out chr_@ --threads 4 --ld_out``

Imputing missing parental genotypes
-----------------------------------

This is performed using the :ref:`impute.py <impute.py>` script. 
To impute the missing parental genotypes without using phase information, use this command:

    ``impute.py --ibd chr_@.ibd --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4``

The script constructs a pedigree from the output of KING's relatedness inference (king.kin0),
and age and sex information (agesex.txt). 
The pedigree along with the IBD segments shared between siblings recorded in chr_1.ibd.segments.gz are used to impute missing parental genotypes
from the observed sibling and parental genotypes in chr_1.bed. 
The imputed parental genotypes are output to a :ref:`HDF5 file <imputed_file>`, chr_1.hdf5. 

If phased haplotypes are available in .bgen format, the imputation can use these as input, which improves the accuracy of the imputation. 
To perform imputation from the phased .bgen file in example_data/, use the following command:

    ``impute.py --ibd chr_@.ibd --bgen chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4``

As with the ibd.py script, the impute_runner.py script can use a user input :ref:`pedigree file <pedigree>` (with the *--pedigree* argument) rather than the *--king* and *--agesex* arguments.

Family-based GWAS with imputed parental genotypes
-------------------------------------------------

This is performed using the :ref:`gwas.py <gwas.py>` script, which implements different family-based GWAS designs.

The default design with imputed parental genotypes is described in Young et al. (2022) (https://www.nature.com/articles/s41588-022-01085-0).
This regresses the proband's phenotype jointly onto the proband's genotype, the father's imputed or observed genotype, and the mother's imputed or observed  genotype.
The regression produces variant level summary statistics on the direct genetic effect of the proband's genotype, the non-transmitted coefficient (NTC) for the father's genotype, and the NTC for the mother's genotype.
To use this design, supply the imputed parenta genotypes (above) to the gwas.py scrict. For example:

    ``gwas.py phenotype.txt --bed chr_@ --imp chr_@ --chr_range 1 --cpus 1 --out chr_@_young``

This regression is done using a linear mixed model that models phenotypic correlations between siblings,
where sibling relations are stored in the :ref:`output of the imputation script <imputed_file>`. 
The 'sibling variance component' output is the phenotypic variance explained by mean differences between sibships. 
:code:`--cpus` allows you to distribute computation across several processes to speed up analyses.

To use the .bgen file instead, use this command:

    ``gwas.py phenotype.txt --bgen chr_@ --imp chr_@ --cpus 1 --out chr_@_young``

The gwas.py script enables the use of different regression designs described in Guan et al. (https://www.nature.com/articles/s41588-025-02118-0).
For the homogeneous ancestry samples typically used in GWAS, you can increase statistical power for estimation of direct genetic effects by inclusion
of individuals without genotyped first-degree relatives (i.e., singletons) in the analysis. This also produces estimates of population effects (as estimated by standard GWAS)
that are almost identical to performing a standard GWAS in a linear mixed model framework. This is why we called this approach the 'unified estimator' in Guan et al. 
We give an example command here: 

    ``gwas.py phenotype.txt --bgen chr_@_trios_singletons --imp chr_@ --cpus 1 --impute_unrel --out chr_@_unified``

The ``--impute_unrel`` flag instructs *snipar* to linearly impute parental genotypes of singletons and include them into the analysis.

Approaches relying on imputed parental genotypes can be biased in strongly structured (Fst>0.01) and/or admixed samples. 
In Guan et al., we develop a 'robust' estimator that maximises power in such samples without introducing bias. 
Here's an example command invoking the robust estimator: 

    ``gwas.py phenotype.txt --bgen chr_@ --imp chr_@ --cpus 1 --robust --out chr_@_robust``

By default, the script outputs summary statistics in a :ref:`gzipped text file <sumstats_text>`: chr_1.sumstats.gz;
In addition to the text summary statistics, :ref:`HDF5 format summary statistics <sumstats_hdf5>` are also output to chr_1.sumstats.hdf5.
Alternatively, you can specify the output filename using the ``--out`` command: for example, with ``--out chr_@_X``, the script will
output the results to chr_1_X.sumstats.gz and chr_1_X.sumstats.hdf5; if '@' is not in the output suffix, e.g., ``--out gwas``, the results will
be stored in gwas_chr_1.sumstats.gz and gwas_chr_1.sumstats.hdf5.

Now we have estimated SNP effects. To compare the Young et al. (or robust or unified) estimates to the true effects, run
    
    ``python estimate_sim_effects.py chr_1_young.sumstats.hdf5 phenotype.effects.txt``

This should print estimates of the bias of the effect estimates.

The bias estimates for direct, paternal NTCs, maternal NTCs, and average NTCs should not be statistically significantly different from 
zero (with high probability). Population effects (as estimated by standard GWAS) are biased estimates of direct effects for this simulated 
phenotype because they also include indirect genetic effects and other confounding factors not modelled here. 

Family-based GWAS without imputed parental genotypes
----------------------------------------------------

Family-based GWAS can also be performed without imputed parental genotypes. In this case, only probands with genotypes for both parents and/or siblings available will be used.
In order to do this, one must provide a pedigree to gwas.py, as in:

    ``gwas.py phenotype.txt --out trios_sibs --bgen chr_@_trios_sibs --pedigree pedigree.txt --cpus 1``

With the above command, the script will default to meta-analysing samples with both parents genotyped (trios) and samples with siblings but without both parents genotyped. 
Alternatively, users can supply one of the following two options (``--robust`` is not applicable since it requires information derived from the imputation procedure):

- ``--sib_diff``: individuals with sibling genotypes will be used in a 'sib-GWAS' design using genetic differences between siblings, and those without will not be considered for the analysis;
- ``--impute_unrel``: individuals with both parents' genotypes and singletons will be used; individuals with sibling genotypes but incomplete parental genotypes will be ignored.

For example:

    ``gwas.py phenotype.txt --out sibs --bgen chr_@_trios_sibs --pedigree pedigree.txt --cpus 1 --sib_diff``

    ``gwas.py phenotype.txt --out trios_sibs_singletons --bgen chr_@_trios_sibs_singletons --pedigree pedigree.txt --cpus 1 --impute_unrel``

Correlations between effects
----------------------------

*snipar* provides a script (:ref:`correlate.py <correlate.py>`) to compute correlations between direct and population effects and between direct effects and average NTCs. 
To compute these correlations from the effects estimated in this tutorial (output by gwas.py to chr_1_young.sumstats.gz) 
using the LD scores computed by ibd.py (and output to chr_1.l2.ldscore.gz), use the following command: 

    ``correlate.py chr_@_young effect --ldscores chr_@``

This should give a correlation between direct effects and average NTCs of close to 0.5. The estimated correlations
and their standard errors, estimated by block-jacknife, are output to effect_corrs.txt. 

The method is similar to LDSC, but correlates the marginal effects (not joint-fit effects adjusted for population stratification, as LDSC attempts to use), 
adjusting for the known sampling variance-covariance matrix of the effects. The LD scores are used for weighting. LD scores output by LDSC can be input. If LD scores are not available, they can be
computed from .bed files by providing them through the --bed argument to :ref:`correlate.py <correlate.py>`. 

Polygenic score analyses
------------------------

For an exercise involving polygenic score analysis, please see the :ref:`Simulation Exercse <simulation>`.

.. In addition to family based GWAS, *snipar* provides a script (:ref:`pgs.py <pgs.py>`) for computing polygenic scores (PGS) based on observed/imputed genotypes,
.. and for performing family based polygenic score analyses. 
.. Here, we give some examples of how to use this script. The script computes a PGS
.. from a :ref:`weights file <weights>`. 
.. For the tutorial, we provide a weights file (direct_weights.txt) in `LD-pred <https://github.com/bvilhjal/ldpred>`_ format
.. where the weights are the true direct genetic effect of the SNP. 

.. To compute the PGS from the weights in direct_weights.txt, use the following command:

..     ``pgs.py direct --bed chr_@ --imp chr_@ --weights direct_weights.txt``
    
.. This uses the weights in the weights file to compute the PGS for each genotyped individual for whom observed or imputed parental genotypes are available.
.. It outputs the PGS to a :ref:`PGS file <pgs_file>`: direct.pgs.txt. 

.. To estimate direct, paternal, and maternal effects of the PGS, use the following command:

..     ``pgs.py direct --pgs direct.pgs.txt --phenofile phenotype.txt``

.. This uses a linear mixed model that has a random effect for mean differences between families (defined as sibships here) and fixed effects for the direct,
.. paternal, and maternal effects of the PGS. It also estimates the 'population' effect of the PGS: the effect from regression of individuals' phenotypes onto their PGS values.
.. The estimated effects and their standard errors are output to direct.effects.txt, described :ref:`here <pgs_effects>`. 
.. The sampling variance-covariance matrix of the direct effect and paternal and maternal NTCs is output to direct.vcov.txt, described :ref:`here <pgs_vcov>`.

.. Estimates of the direct effect of the PGS should be equal to 1 in expectation since
.. we are using the true direct effects as the weights, so the PGS corresponds to the true direct effect component of the trait.
.. The paternal/maternal NTC estimates capture the correlation between the direct and indirect parental effects. The population effect estimate
.. should be greater than 1, since this captures both the direct effect of the PGS, and the correlation between direct and indirect parental effects.

.. If parental genotypes have been imputed from sibling data alone, 
.. then imputed paternal and maternal PGS are perfectly correlated, 
.. and the above regression on proband, paternal, and maternal PGS becomes collinear. 
.. To deal with this, add the --parsum option to the above command, 
.. which will estimate the average NTC rather than separate maternal and paternal NTCs.