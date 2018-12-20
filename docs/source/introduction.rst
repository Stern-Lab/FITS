Introduction
============

FITS requires two types of input: `Data file`_ and `Parameters file`_.

.. _data_file:

Data file
^^^^^^^^^
This file is expected to hold observed allele information from the system under study. FITS expects a tab-delimited textual file, with three columns: 

* ``gen`` for the generation of the observation
* ``base`` for the observed state 
* ``freq`` for the measured frequency for that state

.. csv-table:: Example data file
  :file: examples/data_file_to_show.txt
  :delim: tab
  :header-rows: 1

You can also download an :download:`example <examples/data_file_example.txt>`. 

.. note:: For each generation, the sum of frequencies for the different alleles should be 1.
  
.. _parameters_file:
  
Parameters file
^^^^^^^^^^^^^^^
This file provides FITS with population genetics parameters information of the system under study. 

General parameters
******************
====================== ============== ================================ 
Parameter name         Type           Description
====================== ============== ================================
_N                     Integer        Size of the population
---------------------- -------------- --------------------------------
_sample_size           Integer        Size of observed population (e.g., sequenced genomes)
---------------------- -------------- --------------------------------
_bottleneck_size       Integer        Size of the population transferred on a bottleneck event
---------------------- -------------- --------------------------------
_bottleneck_interval\* Integer        Number of generations separating between bottleneck events (default: 0)
---------------------- -------------- --------------------------------
_num_alleles           Integer        Number of alleles observed in all loci
---------------------- -------------- --------------------------------
_mutation_rateX_Y      Float          Rate of mutation of allele X to allele Y
---------------------- -------------- --------------------------------
_single_mutation_rate  Float          Rate of mutation between alleles
---------------------- -------------- --------------------------------
_allele_fitnessX       Float          Fitness value assigned to allele X
---------------------- -------------- --------------------------------
_logistic_growth*      Float          1: assume the population growth in a logistic growth model throughout the generations (default: 0)
---------------------- -------------- --------------------------------
_logistic_growth_K     Float          Logistic model - upper bound.
---------------------- -------------- --------------------------------
_logistic_growth_r     Float          Logistic model - proportionality constant
====================== ============== ================================ 

\*parameter value of 0 means disabled/off; positive values mean enabled/on.

ABC parameters
**************
====================== ============== ================================ 
Parameter name         Type           Description
====================== ============== ================================
_num_repeats           Integer        Size of the population
---------------------- -------------- --------------------------------
_acceptance_rate       Float          Fraction of best simulations
====================== ============== ================================ 

Single simulation
*****************
====================== ============== ================================ 
Parameter name         Type           Description
====================== ============== ================================
_num_generations       Integer        Number of generations to simulate
---------------------- -------------- --------------------------------
_allele_init_freqX     Float          Initial frequency of allele X
====================== ============== ================================ 


Fitness inference parameters
****************************

====================== ============== ================================ 
Parameter name         Type           Description
====================== ============== ================================
_prior_distribution    Text           One of the following:
                                      
                                      uniform
                                      
                                      log_normal (based on `Bons et al. 2018 <https://doi.org/10.1093/ve/vey029>`_)
                                      
                                      fitness_composite

---------------------- -------------- --------------------------------
_allele_min_fitnessX   Float          The minimum fitness value (inclusive) that may be assigned to allele X
---------------------- -------------- --------------------------------
_allele_max_fitnessX   Float          The maximum fitness value (exclusive) that may be assigned to allele X
====================== ============== ================================ 

Mutation rate inference parameters
**********************************
X and Y are alleles defined in the data file (i.e., 0 and 1). 

============================= ============== ================================ 
Parameter name                Type           Description
============================= ============== ================================
_min_log_mutation_rateX_Y     Float          Minimum (inclusive) :math:`n` for mutation rate :math:`10^n` from alleleX to allele Y
----------------------------- -------------- --------------------------------
_max_log_mutation_rateX_Y     Float          Maximum (exclusive) :math:`n` for mutation rate :math:`10^n` from alleleX to allele Y
----------------------------- -------------- --------------------------------
_min_log_single_mutation_rate Float          Minimum (inclusive) :math:`n` for mutation rate :math:`10^n` from *any* alleleX to *any* allele Y
----------------------------- -------------- --------------------------------
_max_log_single_mutation_rate Float          Maximum (exclusive) :math:`n` for mutation rate :math:`10^n` from *any* alleleX to *any* allele Y
============================= ============== ================================

Population size inference parameters
************************************
====================== ============== ================================ 
Parameter name         Type           Description
====================== ============== ================================
_Nlog_min              Float          Minimum (inclusive) exponent :math:`n` for population size :math:`10^n` 
---------------------- -------------- --------------------------------
_Nlog_max              Float          Maximum (exclusive) exponent :math:`n` for population size :math:`10^n`
====================== ============== ================================
