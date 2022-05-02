.. _input files:
===========
Input files
===========

Observed genotypes
------------------
.. _observed genotypes:

Pedigree 
--------
.. _pedigree:
Pedigree file is a space(" ") seperated csv with columns 'FID', 'IID', 'FATHER_ID', 'MOTHER_ID'.
Default NaN value of Pedigree file is '0'. If your NaN value is something else be sure to specify it with --pedigree_nan option.

agesex file
-----------
.. _agesex: 
This is a space(" ") seperated CSV with columns "FID", "IID", "FATHER_ID", "MOTHER_ID", "sex", "age".
Each row contains the age and sex of one individual. Male and Female sex should be represented with 'M' and 'F'.
Age column is used for distinguishing between parent and child in a parent-offsring relationship inferred from the kinship file.
ID1 is a parent of ID2 if there is a 'PO' relationship between them and 'ID1' is at least 12 years older than ID2.

phenotype file
--------------
.. _phenotype:

weights file
------------
.. _weights: 