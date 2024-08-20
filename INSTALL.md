# INSTALLATION

This file contains installation instructions and examples for the Python and C++ modules. The discretization code works in Windows and Linux, however the C++ calibration functions are currently tested only on Linux. We also provide instructions for creating a [Docker](#docker) container.

## Download

To download the C++ library and the python package, run:

    git clone https://github.com/vguigues/LASPATED.git

 
## Python

Enter the laspated directory and run:

    pip install .

This command will download all Python dependencies and install laspated on the PYTHONPATH, so you can use it in Python via:

    import laspated

See the demo directory for examples on how to run the discretization functions and how to run the C++ calibration functions from Python.

### PIP

The python package is also available via PIP, with the command

    pip install laspated





## C++

The C++ calibration functions have the following dependencies:

- Boost 1.74 ([Home page](www.boost.org))
- xtl ([Repository](https://github.com/xtensor-stack/xtl))
- xtensor ([Repository](https://github.com/xtensor-stack/xtensor))
- Gurobi 11.0 ([Home page](www.gurobi.com)) (Optional)

xtl and xtensor are provided in the Model_Calibration/Cpp directory. However, Boost must be installed. The libboost-all-dev package contains Boost and can be installed via the most famous package managers, for example:

    sudo apt install libboost-all-dev

Boost can also be compiled from source, following instructions from its [Home page](https://www.boost.org).

The Gurobi solver can be downloaded from its [Home page](https://www.gurobi.com). The installation is done by extracting its content to some location on your system, and setting an environment variable GUROBI_HOME to the location. Since LASPATED uses the C++ interface, the user must also compile the Gurobi C++ classes following this [guide](https://support.gurobi.com/hc/en-us/articles/360039093112-How-do-I-resolve-undefined-reference-errors-while-linking-Gurobi-in-C).

### Compiling the laspated executable

Once the dependencies are installed, you must compile the laspated executable. From the repository main directory, run:

    cd Model_Calibration/Cpp
    make USE_GUROBI=1 GUROBI_VER=110

If you don't have Gurobi installed, you can still run the laspated functions for the model without covariates, by compiling with

    cd Model_Calibration/Cpp
    make USE_GUROBI=0


## Docker

We also provide a Dockerfile. The docker container comes with all dependencies installed and can be built in the project root directory. To build the container, run:

```
docker build -t laspated-dock .
```

If you have a Gurobi Web License, you can build the container with Gurobi support by running:

```
docker build --build-arg USE_GUROBI=1 -t laspated-dock .
```

The above command will build the container with Gurobi 11.0.1 installed.

To run the container, use:
```
docker run -it laspated-dock
```

To run the container with Gurobi support, you can pass the license to the container with:
```
docker run --volume="/absolute/path/to/gurobi.lic:/opt/gurobi/gurobi.lic:ro" -it laspated-dock
```

This will open a shell environment with all dependencies installed. Once in the container environment you can run both the Python and C++ functions, as well as run the replication script.

## Matlab/Octave

In the Model_Calibration directory we also provide a Matlab/Octave implementation of the calibration functions. It can be included as regular Matlab/Octave functions.





