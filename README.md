# Molecular Utilities

A suite of tools to assist with molecular simulations and analysis.

## Contents

* [Introduction](#Introduction)
* [Requirements](#Requirements)
* [Building the programs](#Building)

## <a name="Introduction"></a> Introduction

This software suite is intended to work with file formats from the Protein data bank ([PDB](https://www.rcsb.org/)) and the Large-scale Atomic/Molecular Massively Parallel Simulator ([LAMMPS](http://lammps.sandia.gov/), Sandia National Laboratories).

All programs and scripts can be run without command line options for user instructions and example usage.

Documentation for the compiled utilities can be found [here](https://github.com/JohnGrime/MolecularUtilities/tree/master/Programs), and documentation for the Python utilities can be found [here](https://github.com/JohnGrime/MolecularUtilities/tree/master/Scripts).

## <a name="Requirements"></a> Requirements

* Compiler supporting [c++11](https://en.wikipedia.org/wiki/C%2B%2B11)
* [Gnu make](https://www.gnu.org/software/make/)
* [Python](https://www.python.org/) (for using the PDB utility script)

## <a name="Building"></a> Building the programs

Simply call `make` in the project directory to build the programs; the binaries are placed into a `bin` subdirectory.

The fluctuation spectrum analysis program can be compiled to use OpenMP, which could offer a significant performance increase for large wave numbers. This is accomplished bby building the programs with `OMP=yes`:

`make  OMP=yes`

Note that you may need to change the `CC` variable if your default compiler does not support OpenMP (e.g. clang on macOS). If you have installed an OpenMP-compliant compiler, specify it as appropriate, e.g.:

`make  CC=/opt/local/bin/g++-mp-6  OMP=yes`

... and you should then have a fluctuation spectrum analysis program with OpenMP acceleration.
