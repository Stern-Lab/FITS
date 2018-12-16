Compiling from source
=====================

In order to compile, FITS requires the `Boost library 1.61 <https://sourceforge.net/projects/boost/files/boost/1.61.0/>`_ or 1.60 and a C++11 supporting compiler. We used GCC5.3 on Linux (Centos), Clang provided with Xcode9 on MacOS (High Sierra) and MinGW provided with Qt 5.9 on Windows 10 & 7.

Compiling the command line interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Download and install Boost 

#. Compile all \*.cpp files, referring the compiler to the Boost libraries, e.g.: ``g++ -std=c++11 -O3 -o fits1.02 -I/boost_1_61/include -L/boost_1_61/libs *.cpp`` (here we tell gcc to use c++11 standard, use optimization (-O3) and name the output file (-o) **fits1.02**).

  .. note:: For Windows - you need to install the `prebuilt binaries <https://sourceforge.net/projects/boost/files/boost-binaries/1.61.0/>`_. For MinGW-32 we used `msvc-14.0-32 <https://sourceforge.net/projects/boost/files/boost-binaries/1.61.0/boost_1_61_0-msvc-14.0-32.exe/download>`_.

#. Run FITS with no command line arguments to get the help text

#. Run FITS with the proper syntax (see :ref:`cli`) in order to generate data or infer the required parameter

.. note:: In order to run g++, either run it from the installation folder (of MinGW, in Windows) or add it to ``PATH``.


Compiling the graphical user interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to compile the GUI, you will need (in addition to Boost) the Qt framework, along with Qt Creator. Both are available in open source license `here <https://www.qt.io/download>`_.

#) After installing Qt creator, open the project file (**fits_gui.pro**)

#) Within the project file, replace placeholder text next to ``__INCLUDEPATH__`` with the path to the boost library

#) Run qmake (Build/Run qmake)

#) Build (Build/Build All)

#) Run FITS

.. note:: If Qt is already installed on the system (e.g. through Anaconda) clashed may occur, leading to unexpected results. If possible, we recommend to remove the additional installations (e.g. by running _conda remove qt_)
