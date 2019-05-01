.. _results:

FITS output
===========

General
^^^^^^^
FITS infers population genetics parameters using the Approximate Bayesian Computation (ABC) method. 
The output of this method is a distribution of values that explain the observed allele frequencies with the highest probabilities (also called *the posterior distribution*).
A common practice is to take the **median** of this distribution as the inferred value of the parameter under study.  

.. note :: If the p-value for Levene's test is not significant (>=0.05) or Nµ is small (<1), FITS will mark the relevant line in the results with an asterisk (*). This result is considered **unreliable**.

The results below are outputted for all inferences.

===================== ================================ 
Result header         Description
===================== ================================
median                  The median value of the posterior distribution. This is practically the inferred population genetics parameter.
--------------------- --------------------------------
MAD                      Median Absolute Deviation index (`MAD <https://en.wikipedia.org/wiki/Median_absolute_deviation>`_) of the posterior distribution. 
--------------------- --------------------------------
min                      The minimum value in the posterior distribution.
--------------------- --------------------------------
max                      the maximum value in the posterior distribution.
--------------------- --------------------------------
pval                  The result of a statistical test about the informativeness of the posterior distribution, with a null hypothesis that the posterior distribution is as informative as the prior distribution. 
===================== ================================ 

Fitness inference results
^^^^^^^^^^^^^^^^^^^^^^^^^
In addition to the `General`_ reported values, in fitness inference more data are available:
 
===================== ================================ 
Result header         Description
===================== ================================
allele                  The allele for which the results are reported
--------------------- --------------------------------
DEL(%)                  The proportion of the posterior distribution with values below 1. 
--------------------- --------------------------------
NEU(%)                  The proportion of the posterior distribution with values equal to 1.
--------------------- --------------------------------
ADV(%)                  The proportion of the posterior distribution with values above 1.
--------------------- --------------------------------
category              A possible classification of the allele into {LETHAL,DEL,NEU,ADV}, based on the inferred fitness value. 
===================== ================================ 

.. note :: Some `fitness priors <_static/priors.png>`_ rarely choose the exact value of 1 and therefore **NEU(%)** will approach zero, even for neutral alleles. 

Mutation rate inference results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FITS infers the mutation rates between all defined alleles. Accordingly, the output table contains the target allele in the first row and the source allele in the first column.

.. note :: In Evolve & Resequence (E&R) studies, when the population is homogeneous at first generation, 
           in the absence of more information the inference of the rates between the minor allele to the major will be insignificant, so the **pval** should be taken into account.   

Population size inference results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See the `General`_ inference results.  
