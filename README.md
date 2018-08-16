# FITS
FITS (Flexible Inference from Time-Series) uses Approximate Bayesian Computation (ABC) to infer information from time-series data of allele frequencies. The frequency of an allele in a population depends on the following factors: (a) the replicative fitness of the allele as compared to the wild-type allele, (b) the population-wide mutation rate, and (c) the population size. If one has information on two of these three factors, FITS may be used to infer the missing one. A typical example may be an evolve-and-re-sequence experiment, where one obtains the frequency of a mutation over time. If the mutation rate and population size used in the experiment are known, one may use FITS to infer the fitness of the mutation. 

## Usage
### Graphical interface (GUI)
Download the software (GUI binaries). For Windows - extract the zip file and run the application. For Mac - you can run the app from within the dmg archive or extract it and then run it. A set of sample files are available to demonstrate the abilities of FITS (Sample_files.zip). After running FITS, you may load the parameters file (e.g. the sample file Params_infer_fitness_N5_u05_reps05.txt) and data file (e.g. the file sim_data_N5_u05_w0.5.txt). The parameters dictate which parameter FITS will infer (fitness in this example). Press Go!
FITS shows a progress bar and estimated time to finish running (typical running time is 30 seconds). A report is given in the bottom right side. Inferred values are the medians, given in the table. Version 1.0 also explicitly shows the value in the Quick Report.

#### Note: FITS expects the data file to be tab-delimited. If using Excel, you may need to save the file as Windows Formatted Tab Delimited file. Verify the content and format of the file if FITS fails to run.

#### Note: Parameter value may be changed by pressing the _Edit file_ button.

### Command line
Running fits with no parameters gives the help screen, listing the possible usage syntaxes. For fitness inference, as an example, the syntax is:
```
fits -fitness <param_file> <actual_data_file> <posterior_file> <summary_file> (optional: <prior_file>)
```
FITS needs to receive a parameter file, a data file containing a time-series of mutation frequencies, an output file for the posterior distribution and an output file for the summary (report) file. The prior distribution may also be outputted to a file. The summary file contains a summary of the inference. The posterior file is mostly helpful for downstream analysis. **All files except the prior file are mandatory for running FITS.**

## Compiling from source
### Command line (fitsX.XX_src\*)
In order to compile, FITS requires the [Boost library 1.61](https://sourceforge.net/projects/boost/files/boost/1.61.0/) or 1.60 and a C++11 supporting compiler. We used GCC5.3 on Linux (Centos), Clang supplied with Xcode9 on MacOS (High Sierra) and MinGW supplied with Qt 5.9 on Windows 10 & 7.
1. Download and install Boost (**Note:** For Windows - you need to install the [prebuilt binaries](https://sourceforge.net/projects/boost/files/boost-binaries/1.61.0/). For MinGW-32 we used **boost_1_61_0-msvc-14.0-32.exe**).
2. Compile all the *\*.cpp* files, referring the compiler to the Boost libraries, e.g.:
```
g++ -std=c++11 -O3 -o fits1.02 -I/boost_1_61/include -L/boost_1_61/libs *.cpp
```
(here we tell gcc to use c++11 standard, use optimization (-O3) and name the output file (-o) **fits1.02**.)

3. Run FITS with no command line arguments to get the help text
4. Run FITS with the proper syntax (see _Command line_) in order to generate data or infer the required parameter

### GUI (fits_guiX.XX_src\*)
In order to compile the GUI, you will need (in additional to Boost) the Qt framework, along with Qt Creator. Both are available in open source license from: https://www.qt.io/download.
1. After installing Qt creator, open the project file (**fits_gui.pro**)
2. Within the project file, replace placeholder text next to __INCLUDEPATH__ with the path to the boost library
3. Run qmake (Build\Run qmake)
4. Build (Build\Build All)
5. Run FITS

#### Note: If Qt is already installed on the system (e.g. through Anaconda) clashed may occur, leading to unexpected results. If possible, we recommend to remove the additional installations (e.g. by running _conda remove qt_)


# License
FITS is distributed under the [GPLv3 license](https://www.gnu.org/licenses/licenses.html#GPL), and uses the [Qt framework](https://www.qt.io/).
