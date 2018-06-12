# Molecular Utilities

A suite of tools to assist with molecular simulations and analysis.

## Contents

* [Introduction](#Introduction)
* [Requirements](#Requirements)
* [Building the programs](#Building)
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

## <a name="Introduction"></a> Introduction

Protein data bank format ([PDB](https://www.rcsb.org/))
Large-scale Atomic/Molecular Massively Parallel Simulator ([LAMMPS](http://lammps.sandia.gov/))

All program can be run without command line options for user instructions and example usage.

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

## <a name="AxisAlign"></a> AxisAlign

Align major molecular axes of a PDB structure to the specified Cartesian axes. Supports alignment on a filtered subset of particles.

## <a name="BestStructuralMatch"></a> BestStructuralMatch

Find molecular structure in a data set that matches most closesly a refernece structure (via root-mean-squared-deviation, RMSD).

## <a name="Centroids"></a> Centroids

Extract the centroid "pseudoparticles" from a set of PDB input structures.

## <a name="Distances"></a> Distances

Measure statistics regarding the distances between PDB atoms in multiple data sets.

## <a name="FluctuationSpectrum"></a> FluctuationSpectrum

Generate the fluctuation spectra for a bilayer membrane from LAMMPS trajectories (using LAMMPS trajectory frame class)

## <a name="FluctuationSpectrum2"></a> FluctuationSpectrum2

Generate the fluctuation spectra for a bilayer membrane from LAMMPS trajectories (using LAMMPS config class)

## <a name="GenerateMembranes"></a> GenerateMembranes

Create monolayer and/or bilayer systems using user-defined lipid molecules and surface geometries.

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
