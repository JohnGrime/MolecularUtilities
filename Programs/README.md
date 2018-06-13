# Contents

This guide assumes basic familiarity with the [PDB file format](https://www.rcsb.org/), specifically the various data entries in [ATOM](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)) fields (which are typically used when defining filters for atomic data in the utilities etc).

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

This program will align the molecular axes of a PDB structure to the specified Cartesian axes, with the option to use only a filtred subset of the particles in the calculation of the molecular axes.

This functionality is useful for visualization, and aligned molecules are useful for e.g. building up complicated molecular systems with specific initial conditions.

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

Here, we define a filter such that only the atoms with the name `CA` (i.e., carbon alpha atoms) are used to calculate the molecular axes. The results are shown below, where the initial structure is rendered in red and the aligned structure is in green. At the lower left of the image can be seen a set of axis indicators, where Cartesian x-, y- and z-axes are shown using red, green, and blue arrows respectively.

![aligned structure](../Images/AA_1.png)


## <a name="BestStructuralMatch"></a> BestStructuralMatch

This program finds the molecular structure in a data set that has the lowest root-mean-squared deviation (RMSD) from a reference structure after superposing the two structures onto one another such that RMSD is minimized as much as possible.

This functionality is useful to identify one or more "typical" conformations of a molecule, given an ensemble of structures.

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/BestStructuralMatch 
	Usage: bin/BestStructuralMatch reference.pdb set_size data.pdb [filter_name=v1,v2,...] [filter_name=v1,v2,...]
	Where:
	  - reference.pdb : the reference structure(s) to compare to
	  - set_size : number of consecutive PDB molecules in the input file to group
	  - data.pdb : the potential candidate structures to compare to reference.pdb
	  - OPTIONAL filters : used in the superposition/RMSD comparison
	Notes:
	  Consecutive set_size entries from reference.pdb and data.pdb are grouped for superposition/RMSD calculation.
	  Closest match to reference in data set saved to 'best_match.pdb', with some info in REMARK lines.
	  ALL atom data is written in 'best_match.pdb', not just the atoms used in the superposition/RMSD calculation.
	MolecularUtilities $ 

*Example*: the [3P05](https://www.rcsb.org/structure/3P05) PDB file contains a ring of five proteins from the human immunodeficiency virus type 1 (HIV-1). Although the proteins in this ring are in principle identical, the experimental structures for each protein are slightly different. Let's see how different they are:

	MolecularUtilities $ bin/BestStructuralMatch PDB_sources/3P05.pdb 1 PDB_sources/3P05.pdb name=CA resSeq=1-145
	Set     1 : 136 common filtered atoms : RMSD        0.000 =>        0.000 <- best so far
	Set     2 : 126 common filtered atoms : RMSD       29.868 =>        0.209
	Set     3 : 126 common filtered atoms : RMSD       48.329 =>        0.159
	Set     4 : 126 common filtered atoms : RMSD       48.516 =>        0.192
	Set     5 : 125 common filtered atoms : RMSD       29.813 =>        0.164
	Unable to load 1 structures from 'PDB_sources/3P05.pdb'; only managed 0! Assuming EOF and stopping here.
	MolecularUtilities $

Here we specify to only use atoms with name `CA` for the comparison, and only use atoms from residues numbered 1 to 145 (these form the "head" of each protein).

Not surprisingly, the closest match to the reference structure (specified as the first entry in the reference PDB file, `3P05.pdb`) is ... the first entry in the `3P05.pdb` file! However, we can see that the structure of the "heads" of other proteins in the ring are very similar (< 1 Angstrom RMSD).

## <a name="Centroids"></a> Centroids

Extract the centroid "pseudoparticles" (i.e. the averaged position of each specified particle) from a set of PDB input structures.

This functionality is useful to generate the average location of particles whose locations fluctuate. It can be combined with the [BestStructuralMatch](#BestStructuralMatch) program to generate the "average" structure of a molecule given a set of similar natural conformations, and then extracting which of the input conformations is the closest match to this "average" structure.

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/Centroids 
	Usage: bin/Centroids input.pdb set_size
	Where:
	  - set_size : number of consecutive PDB molecules in the input file to group
	MolecularUtilities $ 

As shown in the output information above, the user may specify a `set_size` parameter to indicate how many consecutive entries in the PDB file are treated as a single "molecule" (with entries delineated by `TER` lines). This parameter is typically `1`.


## <a name="Distances"></a> Distances

Measure statistics regarding the distances between PDB atoms in multiple data sets.

This functionality can be useful when estimating effective particle volumes in a coarse-grained model using fine-grained structural data, and therefore for initial estimates for force field parameters etc.

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/Distances
	Usage: bin/Distances input=name:filepath.pdb[:min_samples[:set_size]] input=name:filepath.pdb ... rcut=X [filters="filter_string;filter_string;..."] [same="name:resSeq,name:resSeq;name:resSeq,name:resSeq;..."] [histogram_prefix=X] [histogram_res=X]
	Where:
	  input : define an input PDB:
	    -name : name for input set.
	    -min_samples : OPTIONAL minimum number of samples required to print distances (default = 1).
	    -set_size : OPTIONAL only measure over consecutive 'set_size' mol groups (useful for structures superposed onto common reference frame, default = all in file).
	  filters : PDB-style filtering.
	  same : OPTIONAL definition of atoms to consider the same, to generate/prints additional distance info.
	  histogram_prefix : OPTIONAL prefix for saved histograms of distances; if not specified, no histograms written.
	  histogram_res : OPTIONAL resolution (bins per unit distance) for histograms (ignored if histogram_prefix not defined).
	Examples:
	  bin/Distances input=test:blah.pdb:4 rcut=10.0 filters="name:CA;resSeq:3,6,12-45,112-116" same="CA:18,GCA:18;CA:45,GCA:45" 
	  bin/Distances input=test1:blah1.pdb:4 input=test2:blah2.pdb:4:2 rcut=10.0 filters="name:CA;resSeq:3,6,12-45,112-116" same="CA:18,GCA:18;CA:45,GCA:45" 
	MolecularUtilities $ 


## <a name="FluctuationSpectrum"></a> FluctuationSpectrum

Generate the fluctuation spectra for a bilayer membrane from LAMMPS trajectories (using LAMMPS trajectory frame class)

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/Distances
	MolecularUtilities $ 


## <a name="FluctuationSpectrum2"></a> FluctuationSpectrum2

Generate the fluctuation spectra for a bilayer membrane from LAMMPS trajectories (using LAMMPS config class)

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/Distances
	MolecularUtilities $ 


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
