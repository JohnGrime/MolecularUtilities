# Contents

**AxisAlign** : align PDB structure to specified Cartesian axes using filtered atom sets

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

The `Misc/` directory also includes a simple implementation of a LAMMPS pair potential for the Brannigan-Philips-Brown lipid model.


# Building the programs

Simply call `make` in the directory to build the programs; the binaries are placed into the `bin` directory.


# Running the programs

Run the programs with no command line options for user instructions and example usage.


# Use OpenMP in the FluctuationSpectum program

As the fluctuation spectrum program can use OpenMP, the programs may be built with e.g. `OMP=yes` on the command line:

`make  OMP=yes`

Note that you may need to change the `CC` variable if your default compiler does not support OpenMP (e.g. clang on macOS). If you have installed an OpenMP-compliant compiler, specify it as appropriate, e.g.:

`make  CC=/opt/local/bin/g++-mp-6  OMP=yes`

... and you should then have a FluctuationSpectrum binary with OpenMP acceleration. This can be significantly faster when using large data/wave numbers!
