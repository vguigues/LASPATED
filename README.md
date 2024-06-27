# LASPATED: Library for Analysis of Spatio-Temporal Discrete Data

Arxiv paper: https://arxiv.org/abs/2401.04156

This repository contains the source code for the "LASPATED: A Library for the Analysis of
Spatio-Temporal Discrete Data" paper ([Arxiv link](https://arxiv.org/abs/2401.04156)).

## OS Support

LASPATED Python and Matlab code were tested under Windows 10/Linux systems. The C++ code were tested under Linux systems only.

## Data Directory
Directory Data contains Shapefiles and data files that are used in the code examples. It includes a Rio de Janeiro and New York border Shapefiles, a csv and Shapefile with location of ambulance bases from Rio, covariates Data from Rio, and a Shapefile with Rio's administrative regions. It also includes csv files with Ambulance emergency data from Rio, and a mock csv simulating emergency data from New York. 

## Discretization

Contains examples of how to generate Discretizations from data. File python_dependencies.txt describes the python modules that must be installed to run the discretization functions. sorted_events.csv is a csv example file with emergency data from Rio de Janeiro. laspated directory contains the discretization functions.


## laspated

Contains a python module that can be installed on your system via pip (currently under maintenance).

## Model_Calibration

Contains the source code of the Matlab, Octave and C++ functions.

## Replication

Contains a script and auxiliary files on how to replicate the results present in the paper. The README file inside contains a tutorial on how to compile and run the script


