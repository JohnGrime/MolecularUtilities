# Molecular Utilities

A suite of tools to assist with molecular simulations and analysis.

## Contents

* [Requirements](#Requirements)
* [Platform Compatibility](#Compatibility)
* [Building the programs](#Building)
* [Running the programs](#Running)
* The utilities:
  * [AxisAlign](#AxisAlign)
  * [BestStructuralMatch](#BestStructuralMatch)
  * [Centroids](#Centroids)
  * [Distances](#Distances)
  * [FluctuationSpectrum](#FluctuationSpectrum)
  * [FluctuationSpectrum2](#FluctuationSpectrum2)
  * [GenerateMembranes](#GenerateMembranes)
  * [LammpsCombiner](#LammpsCombiner)
  * [LammpsToXYZ](#LammpsToXYZ)
  * [SphereArbitrary](#SphereArbitrary)
  * [SphereBySubdivision](#SphereBySubdivision)
  * [Superpose](#Superpose)
  * [UnwrapTrajectory](#UnwrapTrajectory)

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

## <a name="AxisAlign"></a> AxisAlign

Align PDB structure to specified Cartesian axes using filtered atom sets

## <a name="BestStructuralMatch"></a> BestStructuralMatch

Find structure in a data set that matches most closesly (via RMSD) a reference structure

## <a name="Centroids"></a> Centroids

Extract centroid structure from a set of input structures

## <a name="Distances"></a> Distances

Measure statistics on distances between PDB atoms in data sets

## <a name="FluctuationSpectrum"></a> FluctuationSpectrum

Generate bilayer membrane fluctuation spectra from LAMMPS trajectories (using LAMMPS trajectory frame class)

## <a name="FluctuationSpectrum2"></a> FluctuationSpectrum2

Generate bilayer membrane fluctuation spectra from LAMMPS trajectories (using LAMMPS config frame class)

## <a name="GenerateMembranes"></a> GenerateMembranes

Create both monolayer and bilayer systems for multiple user-specfied lipids and geometries.

## <a name="LammpsCombiner"></a> LammpsCombiner

Combine LAMMPS config files, automatically renumbering bonds and angles etc.

## <a name="LammpsToXYZ"></a> LammpsToXYZ

## <a name="SphereArbitrary"></a> SphereArbitrary

Generate a sphere using an arbitrary number of surface points

## <a name="SphereBySubdivision"></a> SphereBySubdivision

Generate a sphere using poyhedral subdivision (cna provide more regular spacing that `SphereArbitrary`)

## <a name="Superpose"></a> Superpose

Superpose arbitrary sets of PDB structures using filtered atom sets

## <a name="UnwrapTrajectory"></a> UnwrapTrajectory

Unwrap molecules in a LAMMPs trajectory so they are not broken across periodic boundaries
