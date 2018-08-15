# FITS
Flexible inference of population genetics parameters from time-series data

## Usage
Download the software (binaries). For Windows - extract the zip file and run the application. For Mac - you can run the app from within the dmg archive or extract it and then run it.

### Sample files
A set of sample files are available to demonstrate the abilities of FITS (**Sample_files.zip**).
Different parameter files are given to infer fitness, mutation rate and population size.
The given datasets are for biallelic population of size (N) 10^5, mutation rate (u) of 10^-5 and fitness (w) of 0, 0.5, 1.0 and 1.3.


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
