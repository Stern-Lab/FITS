Compiling from source
=====================

The most recent version of FITS is available `here <https://github.com/SternLabTAU/FITS/releases/latest>`_.
In order to compile, FITS requires the `Boost library 1.69 <https://www.boost.org/users/history/version_1_69_0.html>`_ and a C++11 supporting compiler. We used gcc 8.2 on Linux (Centos), Clang provided with Xcode9 on MacOS (High Sierra) and MinGW provided with Qt 5.12 on Windows 10 & 7.

.. note:: Before you start, make sure GCC is in the ``PATH``.

Compiling the command line interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Download and extract the `Boost library 1.69 <https://www.boost.org/users/history/version_1_69_0.html>`_ to a convenient location.

#. Compile all \*.cpp files, referring the compiler to the Boost libraries, e.g.: ``g++ -std=c++11 -O3 -o fits1.3 -I/path/to/boost/ -L/path/to/boost/libs *.cpp`` 

   (here we tell gcc to use c++11 standard, use optimization (-O3) and name the output file (-o) **fits1.1.16**).

#. Run FITS with no command line arguments to get the help text

#. Run FITS with the proper syntax (see :ref:`cli`) in order to generate data or infer the required parameter


Compiling the graphical user interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to compile the GUI, you will need (in addition to Boost) the Qt framework, along with Qt Creator. Both are available in open source license `here <https://www.qt.io/download>`_.

.. note:: When installing the Qt framework, make sure you also install the included MinGW compiler. 


#. Extract the source code zip folder in your favorite location.

#. Open **Qt Creator**. 

#. Click **File>Open file or project** and locate the project file (**fits_gui.pro**)

#. Within the project file, replace placeholder text next to ``__INCLUDEPATH__`` with the path to the boost library

#. Click "fits_gui" in the left toolbox and make sure the configuration is set to **Release**. 

#. Click **Build>Build All**. After build completion, FITS executable will be found in the Build directory (you can see where it is under **Projects>Build directory** that is available from the left toolbar.

#. Move fits_gui.exe to a new folder. Open the console and navigate to that folder.

#. - For **windows**: Locate windeployqt.exe under bin directory in the Qt installation folder. Run windeployqt.exe with fits_gui.exe as its sole argument: ``path/to/windeployqt.exe fits_gui.exe``. 
   - For **MacOS**: Locate macdeployqt.exe under bin directory in the Qt installation folder. Run macdeployqt with fits_gui as its sole argument: ``path/to/macdeployqt fits_gui``.
   
