# FITS
Flexible inference of population genetics parameters from time-series data

## Usage
Download the software (binaries). For Windows - extract the zip file and run the application. For Mac - you can run the app from within the dmg archive or extract it and then run it.

## Sample files
A set of sample files are available to demonstrate the abilities of FITS (**Sample_files.zip**).
Try loading the file 
## Motivation
With the advent of next generation sequencing techniques, it is now possible to observe evolution in action by tracking 
changes in allele frequencies in a population across time. 
Essentially, the frequency of an allele in a population depends on the following factors: 
(a) the replicative fitness of the allele as compared to the wild-type allele, 
(b) the population-wide mutation rate, and 
(c) the population size. In principle, if one has information on two of these three factors, it is possible to infer the missing one. 
Here we present a tool called FITS (Flexible Inference from Time-Series) that uses Approximate Bayesian Computation (ABC) to infer the 
value of a missing factor from time-series data of allele frequencies. 

## Compiling from source
### Command line
In order to compile, FITS requires the [Boost library 1.61](https://sourceforge.net/projects/boost/files/boost/1.61.0/) and a C++11 supporting compiler (we used GCC5.3 and Clang supplied with Xcode9).
1. Download and install Boost (extractig from archive is all it takes)
2. Compile all the *\*.cpp* files, referring the compiler to the Boost libraries, e.g.:
```
g++ -std=c++11 -O3 -o fits100 -I/boost_1_61/include -L/boost_1_61/lib *.cpp
```
(here we tell gcc to use c++11 standard, use optimization (-O3) and name the output file (-o) **fits100**.)

3. Run FITS with no command line arguments to get the help text
4. run FITS with the proper syntax in order to generate data or infer the required parameter

### GUI
In order to compile the GUI, you'll need the Qt framework, along with Qt Creator. It's available in open source license from: https://www.qt.io/download.
1. After installing Qt creator, open the project file (**fits_gui.pro**)
2. Within the project file, replace placeholder text next to __INCLUDEPATH__ with the path to the boost library
3. Run qmake (Build\Run qmake)
4. Build (Build\Build All)
5. Run FITS

## Download Binaries - Graphical (GUI)
### Note: version GUI v0.6 is based on fits0.94 (source available)
For __Windows__ - fits_gui0.6_windows.zip

For __macOS__ - fits_gui0.6_macos.dmg
