# Contents

This guide assumes basic familiarity with the [PDB file format](https://www.rcsb.org/), specifically the various data entries in `ATOM` fields (which are typically used when defining filters for atomic data in the utilities etc).

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

Running the program with no command line options revelas a brief user guide:

	MolecularUtilities $ bin/AxisAlign 

	Usage:

	bin/AxisAlign <input PDB> [-]x|y|z [-]x|y|z [-]x|y|z [ "key1=val1,val2-val3;key2=val4" ]

	Cartesian axis order overrides the defauly x,y,z. ALL axes must be
	specified, and a leading '-' indicates the axis is reflected.

	An optional set of filters can be included for the generation of the molecular axes.

	Output is 'mol_axes.pdb' (molecular axes) and 'aligned.pdb' (PDB aligned to specified axes).
	Note that 'mol_axes.pdb' refers to the original PDB file, not the aligned structure.

	Example:

	bin/AxisAlign blah.pdb y -z x "name=CA;resSeq=1,3,12-14,20-50"

	This aligns the major axis of the specified atoms in molecule in blah.pdb to the y axis,
	the secondary axis to the flipped z axis (i.e. {0,0,-1}) and the minor axis to the Cartesian
	x axis.

	MolecularUtilities $ 


**Example**: the alpha helical chain B of the [2ymk](https://www.rcsb.org/structure/2ymk) PDB structure (included in the `PDB_sources` directory) is initially aligned somewhat parallel to the Cartesian x-axis; we can instead align it such that the major molecular axis is parallel to the Cartesian y-axis, the secondary molecular axis is aligned to the Cartesian x-axis, and the minor axis is aligned to the Cartesian z-axis:

	 bin/AxisAlign PDB_sources/2ymk.B.pdb y x z "name=CA"

Here, we define a filter such that only the atoms with the name `CA` (i.e., carbon alpha atoms) are used to calculate the molecular axes. The results are shown below, where the initial structure is rendered in green and the aligned structure is in red. At the lower left of the image can be seen a set of axis indicators, where Cartesian x-, y- and z-axes are shown using red, green, and blue arrows respectively.

![aligned structure](../Images/AA_1.png)

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
