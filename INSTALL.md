# INSTALLATION

This file contains installation instructions and examples for the Python and C++ modules. The discretization code works in Windows and Linux, however the C++ calibration functions are currently tested only on Linux.

## Python

Enter the laspated directory and run:

    pip install .

This command will download all Python dependencies and install laspated on the $PYTHONPATH, so you can use it in Python via:

    import laspated

See the demo directory for examples on how to run the discretization functions and how to run the C++ calibration functions from Python.


## C++

The C++ calibration functions have the following dependencies:

    - Boost 1.74 ([Home page](www.boost.org))
    - xtl ([Repository](https://github.com/xtensor-stack/xtl))
    - xtensor ([Repository](https://github.com/xtensor-stack/xtensor))
    - Gurobi 11.0 ([Home page](www.gurobi.com)) (Optional)

xtl and xtensor are provided in the Model_Calibration/Cpp directory. However, Boost and Gurobi must be installed. The libboost-all-dev package contains Boost and can be installed via the most famous package managers, or compiled from source, following instructions from its [Home page](www.boost.org).

The Gurobi solver can be downloaded from its [Home page](www.gurobi.com). The installation is done by extracting its content to some location on your system, and setting a $GUROBI_HOME environment variable to the location. Since LASPATED uses the C++ interface, the user must also compile the C++ classes following this [guide](https://support.gurobi.com/hc/en-us/articles/360039093112-How-do-I-resolve-undefined-reference-errors-while-linking-Gurobi-in-C).

### Compiling the laspated executable

Once the dependencies are installed, you must compile the laspated executable.

    cd Model_Calibration/Cpp
    make USE_GUROBI=1 GUROBI_VER=110

If you don't have Gurobi installed, you can still run the laspated functions for the model without covariates, by compiling with

    cd Model_Calibration/Cpp
    make USE_GUROBI=0


## Matlab/Octave

In the Model_Calibration directory we also provide a Matlab/Octave implementation of the calibration functions. It can be included as regular Matlab/Octave functions.





