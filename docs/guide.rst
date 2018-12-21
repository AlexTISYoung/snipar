Guide
************

**Introduction**

sibreg is a python library for performing regression with correlated observations within-class.

A sibreg model (:class:`sibreg.model`) consists of a linear model for the mean along
with a vector of class labels that allows for correlations within-class.

The documentation for the sibreg module (:doc:`sibreg`) contains information on how to define a :class:`sibreg.model`,
how to optimise a :class:`sibreg.model`, how to predict from
a :class:`sibreg.model`, and how to simulate a :class:`sibreg.model`.

**Running tests**

To check that the code is working properly and computing likelihoods and gradients accurately, you can
run tests. In the sibreg/tests subdirectory, type

    ``python tests.py``

The output should say

    ``Ran 4 tests in...``

    ``OK``

**Running the code**

In the sibreg/bin subdirectory, there are two scripts: 'pGWAS.py' (:doc:`pGWAS`) and 'sGWAS.py' (:doc:`sGWAS`).

These perform different kinds of robust GWAS using genotype data from siblings (sGWAS)
or genotype data from siblings and imputed parental genotypes (pGWAS).



