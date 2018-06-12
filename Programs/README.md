# Contents

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
