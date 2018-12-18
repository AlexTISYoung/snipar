Guide
************

**Introduction**

regrnd is a python library for fitting L2 regularised regression models with repeated observations from the same sample.
regrnd should be run with Anaconda built with Python 2.7.

A regrnd model (:class:`regrnd.model`) consists of a linear model for the mean along
with a vector of sample labels that allows for repeated observations from the same sample.

The documentation for the regrnd module (:doc:`regrnd`) contains information on how to define a :class:`regrnd.model`,
how to optimise a :class:`regrnd.model` given a L2 regularisation parameter, how to predict from
a :class:`regrnd.model`, and how to simulate a :class:`regrnd.model`.

**Running tests**

To check that the code is working properly and computing likelihoods and gradients accurately, you can
run tests. In the regrnd/ subdirectory, type

    ``python tests.py``

The output should say

    ``Ran 4 tests in...``

    ``OK``

**Running the code**

In the regrnd/ subdirectory, there is a script 'met_analysis.py'. The script takes one
command line argument that is the location of the training_data.csv file. Unless you wish to install
regrnd as a package, the 'met_analysis.py' script should be run from the regrnd/ subdirectory
so that it can import the code in regrnd/regrnd.py.

The cross-validation code is contained in this script, but it is commmented out as it
takes some time to run for all the investigated regularisation parameters.

Without modifying the script, it trains the model on the training data using the
regularisation parameters found by cross-validation. It then predicts
the 'mets' of the test set, and prints the mean squared error on the test set.

