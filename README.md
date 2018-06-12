# Molecular Utilities

A suite of tools to assist with molecular simulations and analysis.

## Requirements

* Compiler supporting [c++11](https://en.wikipedia.org/wiki/C%2B%2B11)
* [Python](https://www.python.org/) (for using the PDB utility script)

## Platform compatibility

In principle, anything with a c++11 compiler and a Python (2+) installation!

## Building the programs

Simply call `make` in the project directory to build the programs; the binaries are placed into a `bin` subdirectory.

The fluctuation spectrum analysis program can be compiled to use OpenMP, which could offer a significant performance increase for large wave numbers. This is accomplished bby building the programs with `OMP=yes`:

`make  OMP=yes`

Note that you may need to change the `CC` variable if your default compiler does not support OpenMP (e.g. clang on macOS). If you have installed an OpenMP-compliant compiler, specify it as appropriate, e.g.:

`make  CC=/opt/local/bin/g++-mp-6  OMP=yes`

... and you should then have a fluctuation spectrum analysis program with OpenMP acceleration.

## Running the programs

All program can be run without command line options for user instructions and example usage.

## Contents

* [AxisAlign][AxisAlign_link] : align PDB structure to specified Cartesian axes using filtered atom sets

**BestStructuralMatch** : find structure in a data set that matches most closesly (via RMSD) a reference structure

**Centroids** : extract centroid structure from a set of input structures

**Distances** : measure statistics on distances between PDB atoms in data sets

**FluctuationSpectrum** : generate bilayer membrane fluctuation spectra from LAMMPS trajectories (using LAMMPS trajectory frame class)

**FluctuationSpectrum2** : generate bilayer membrane fluctuation spectra from LAMMPS trajectories (using LAMMPS config frame class)

**Fuzzball** : generate fuzzball objects for Jesper Madsen's research.

**GenerateBilayers** : create bilayer systems for user specfied-lipids and geometries

**GenerateMembranes** : create both monolayer and bilayer systems for multiple user-specfied lipids and geometries. More powerful than GenerateBilayers, but slightly more complicated.

**LammpsCombiner** : combine LAMMPS config files, automatically renumbering bonds and angles etc.

**LammpsToXYZ** : convert combine LAMMPS config/trajectory files into xyz format

**SphereArbitrary** : generate a sphere using an arbitrary number of surface points

**SphereBySubdivision** : generate a sphere using poyhedral subdivision

**Superpose** : superpose arbitrary sets of PDB structures using filtered atom sets

**UnwrapTrajectory** : unwrap molecules in a LAMMPs trajectory so they are not broken across periodic boundaries

[AxisAlign_link]: ## AxisAlign

