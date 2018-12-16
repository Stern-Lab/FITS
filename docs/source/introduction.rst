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



 