# FITS
Flexible inference of population genetics parameters from time-series data

## Usage
### Graphical interface (GUI)
Download the software (GUI binaries). For Windows - extract the zip file and run the application. For Mac - you can run the app from within the dmg archive or extract it and then run it.
A set of sample files are available to demonstrate the abilities of FITS (**Sample_files.zip**).
After running FITS, you may load the parameters file (e.g. the sample file _Params_infer_fitness_N5_u05_reps05.txt_) and data file (e.g. the file _sim_data_N5_u05_w0.5.txt_). The available parameters dictate which parameters FITS would be able to infer (only fitness in this example). Select the desired parameter and press Go!. FITS shows a progress bar and estimated time to finish running (typical running time is 30 seconds).
A report is given in the bottom right side. Inferred values are the medians, given in the table. Version 1.0 also explicitly shows the value in the Quick Report.

#### Note: FITS expects the data file to be tab-delimited. If using Excel, you may need to save the file as _Windows Formatted Tab Delimited_ file. Verify the content of the file (e.g. with cat) if FITS fails to run.

#### Note: Parameter value may be changed by pressing the _Edit file_ button.

### Command line
Running fits with no parameters would give the help screen, listing the possible usage syntaxes. For fitness inference,as an example, the yntax is:
```
fits -fitness <param_file> <actual_data_file> <posterior_file> <summary_file> (optional: <prior_file>)
```
You need to give FITS the parameter file, data file, output file for posterior distribution and output file for summary (report) file. You may also store the prior distribution actually used in a file. The summary file would be useful for most users. The posterior file would be helpful for downstream analysis. **All but the prior file are mandatory.**

## Compiling from source
### Command line (fits1.0_src_20180814.\*)
In order to compile, FITS requires the [Boost library 1.61](https://sourceforge.net/projects/boost/files/boost/1.61.0/) or 1.60 and a C++11 supporting compiler (we used GCC5.3 and Clang supplied with Xcode9).
1. Download and install Boost (extractig from archive is all it takes)
2. Compile all the *\*.cpp* files, referring the compiler to the Boost libraries, e.g.:
```
g++ -std=c++11 -O3 -o fits100 -I/boost_1_61/include -L/boost_1_61/lib *.cpp
```
(here we tell gcc to use c++11 standard, use optimization (-O3) and name the output file (-o) **fits100**.)

3. Run FITS with no command line arguments to get the help text
4. run FITS with the proper syntax in order to generate data or infer the required parameter

### GUI (fits_gui1.00_src_20180814.\*)
In order to compile the GUI, you'll need the Qt framework, along with Qt Creator. It's available in open source license from: https://www.qt.io/download.
1. After installing Qt creator, open the project file (**fits_gui.pro**)
2. Within the project file, replace placeholder text next to __INCLUDEPATH__ with the path to the boost library
3. Run qmake (Build\Run qmake)
4. Build (Build\Build All)
5. Run FITS

# License
FITS is distibuted under the [GPLv3 license](https://www.gnu.org/licenses/licenses.html#GPL), and uses the [Qt framework](https://www.qt.io/).
