FITS input
==========

FITS requires two types of input: `Data file`_ and `Parameters file`_.

.. _data_file:

Data file
^^^^^^^^^
This file is expected to hold observed allele information from the system under study. FITS expects a tab-delimited textual file, with following columns: 

#. ``gen`` for the generation of the observation
#. ``allele`` for the observed state 
#. ``freq`` for the measured frequency for that state
#. ``position`` for the position number for which the frequency data is given (optional)

.. note:: FITS assumes the columns to appear in the above order. 

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
Each line in this file represents a different parameter to set, where a space exists between the name of the parameter and its value: ``<parameter_name> <parameter value>``. 

.. note:: If you want to put comments within the parameters file, just add ``#`` at the beginning of the comments' lines. 

You can also download an :download:`example <examples/parameter_file_example.txt>`. 

General parameters
******************
===================== ============== ================================ 
Parameter name        Type           Description
===================== ============== ================================
N                     Integer        Size of population
--------------------- -------------- --------------------------------
sample_size           Integer        Size of observed population (e.g., sequenced genomes)
--------------------- -------------- --------------------------------
bottleneck_size       Integer        Size of the population transferred on a bottleneck event
--------------------- -------------- --------------------------------
bottleneck_interval\* Integer        Number of generations separating between bottleneck events (default: 0)
--------------------- -------------- --------------------------------
num_alleles           Integer        Number of alleles observed in all loci
--------------------- -------------- --------------------------------
mutation_rateX_Y      Float          Rate of mutation of allele X to allele Y. Not required if mutation rate is to be inferred
--------------------- -------------- --------------------------------
fitness_alleleX       Float          Fitness value assigned to allele X. Not required if fitness is to be inferred 
--------------------- -------------- --------------------------------
logistic_growth*      Float          1: model the population growth throughout the generations with a logistic growth model (default: 0)
--------------------- -------------- --------------------------------
logistic_growth_K     Float          Logistic model - upper bound
--------------------- -------------- --------------------------------
logistic_growth_r     Float          Logistic model - proportionality constant
===================== ============== ================================ 

\*parameter value of 0 means disabled/off; positive values mean enabled/on.

ABC parameters
**************
====================== ============== ================================ 
Parameter name         Type           Description
====================== ============== ================================
num_samples_from_prior Integer        How many simulations to perform
---------------------- -------------- --------------------------------
acceptance_rate        Float          Fraction of best simulations to utilize for the inference of the parameter. 
====================== ============== ================================ 

Single simulation
*****************
===================== ============== ================================ 
Parameter name        Type           Description
===================== ============== ================================
num_generations       Integer        Number of generations to simulate
--------------------- -------------- --------------------------------
init_freq_alleleX     Float          Initial frequency of allele X
===================== ============== ================================ 


Fitness inference parameters
****************************
===================== ============== ================================ 
Parameter name        Type           Description
===================== ============== ================================
fitness_prior         Text           | One of the following:
                                     | uniform (for Uniform distribution)
                                     | log_normal (based on `Bons et al. 2018 <https://doi.org/10.1093/ve/vey029>`_)
                                     | fitness_composite
                                     | smoothed_composite (default)
                                     | See the distribution of the above priors on a (0,2) fitness `here <_static/priors.png>`_  
--------------------- -------------- --------------------------------	
min_fitness_alleleX   Float          The minimum fitness value (inclusive) that may be assigned to allele X
--------------------- -------------- --------------------------------
max_fitness_alleleX   Float          The maximum fitness value (exclusive) that may be assigned to allele X
===================== ============== ================================ 

Mutation rate inference parameters
**********************************
X and Y are alleles defined in the data file (i.e., 0 and 1). 

============================ ============== ================================ 
Parameter name               Type           Description
============================ ============== ================================
min_log_mutation_rateX_Y     Float          Minimum (inclusive) :math:`n` for mutation rate :math:`10^n` from alleleX to allele Y
---------------------------- -------------- --------------------------------
max_log_mutation_rateX_Y     Float          Maximum (exclusive) :math:`n` for mutation rate :math:`10^n` from alleleX to allele Y
============================ ============== ================================

Population size inference parameters
************************************
===================== ============== ================================ 
Parameter name        Type           Description
===================== ============== ================================
Nlog_min              Float          Minimum (inclusive) exponent :math:`n` for population size :math:`10^n` 
--------------------- -------------- --------------------------------
Nlog_max              Float          Maximum (exclusive) exponent :math:`n` for population size :math:`10^n`
===================== ============== ================================
