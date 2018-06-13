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

*Example*: the [3P05](https://www.rcsb.org/structure/3P05) PDB file contains a ring of five proteins from the human immunodeficiency virus type 1 (HIV-1). Although the proteins in this ring are in principle identical, the experimental structures for each protein are slightly different due to e.g. thermal fluctuations. Let's see how different they are:

	MolecularUtilities $ bin/BestStructuralMatch PDB_sources/3P05.pdb 1 PDB_sources/3P05.pdb name=CA resSeq=1-145
	Set     1 : 136 common filtered atoms : RMSD        0.000 =>        0.000 <- best so far
	Set     2 : 126 common filtered atoms : RMSD       29.868 =>        0.209
	Set     3 : 126 common filtered atoms : RMSD       48.329 =>        0.159
	Set     4 : 126 common filtered atoms : RMSD       48.516 =>        0.192
	Set     5 : 125 common filtered atoms : RMSD       29.813 =>        0.164
	Unable to load 1 structures from 'PDB_sources/3P05.pdb'; only managed 0! Assuming EOF and stopping here.
	MolecularUtilities $

Here we specify to only use atoms with name `CA` for the comparison, and only use atoms from residues numbered 1 to 145 (these form the "head" of each protein).

The output tells us how many atoms pass the filter for each molecule in the data file, and the RMSD fro the reference structure before and after superposition. Not surprisingly, the closest match to the reference structure (specified as the first entry in the reference PDB file, `3P05.pdb`) is ... the first entry in the `3P05.pdb` file! However, we can see that the structure of the "heads" of other proteins in the ring are very similar (< 1 Angstrom RMSD).

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

*Example*: the [3P05](https://www.rcsb.org/structure/3P05) PDB file contains a ring of five proteins from the human immunodeficiency virus type 1 (HIV-1). The proteins in this ring are in principle identical, but the experimental structures for each protein are slightly different due to e.g. thermal fluctuations. This also means that the separations between pairs of particles in adjacent proteins "around" the ring will likewise differ slightly. Let's take a look at the carbon alpha atoms of the first 20 residues of each protein, ignoring any pairs further apart than 10 Angstrom (1 nanometer):

	MolecularUtilities $ bin/Distances input=test:PDB_sources/3P05.pdb rcut=10 filters="name:CA;resSeq:1-20"
	Inputs:
		 'input=test:PDB_sources/3P05.pdb'
	Parameters:
		 'filters' => 'name:CA;resSeq:1-20'
		 'rcut' => '10'

	*
	* 'test' : 'PDB_sources/3P05.pdb', min_samples 1, set_size <unspecified, using all molecules in file>
	*

	Read 5 molecules.
	79 atoms total passed filtering.
	Adjusting set_size from -1 to 5 (number of filtered molecules)

	  name_i name_j :      min     mean      max   stddev   stderr (   N) Ascending list of distances ...
	      CA:3 CA:6 :    8.804    8.804    8.804    0.000    0.000 (   1)    8.804 
	     CA:3 CA:13 :    8.431    8.431    8.431    0.000    0.000 (   1)    8.431 
	     CA:3 CA:14 :    8.310    8.310    8.310    0.000    0.000 (   1)    8.310 
	      CA:4 CA:6 :    8.362    8.362    8.362    0.000    0.000 (   1)    8.362 
	      CA:4 CA:7 :    9.797    9.797    9.797    0.000    0.000 (   1)    9.797 
	     CA:4 CA:11 :    9.215    9.215    9.215    0.000    0.000 (   1)    9.215 
	     CA:4 CA:12 :    8.097    8.097    8.097    0.000    0.000 (   1)    8.097 
	     CA:4 CA:13 :    7.238    8.555    9.872    1.862    1.317 (   2)    7.238    9.872 
	     CA:4 CA:14 :    8.726    8.726    8.726    0.000    0.000 (   1)    8.726 
	      CA:5 CA:5 :    9.639    9.639    9.639    0.000    0.000 (   1)    9.639 
	      CA:5 CA:6 :    5.844    5.844    5.844    0.000    0.000 (   1)    5.844 
	      CA:5 CA:7 :    6.943    6.943    6.943    0.000    0.000 (   1)    6.943 
	      CA:5 CA:8 :    9.221    9.221    9.221    0.000    0.000 (   1)    9.221 
	     CA:5 CA:11 :    8.944    8.944    8.944    0.000    0.000 (   1)    8.944 
	     CA:5 CA:12 :    9.050    9.050    9.050    0.000    0.000 (   1)    9.050 
	     CA:5 CA:13 :    8.007    8.007    8.007    0.000    0.000 (   1)    8.007 
	      CA:6 CA:9 :    9.417    9.417    9.417    0.000    0.000 (   1)    9.417 
	     CA:6 CA:10 :    9.228    9.228    9.228    0.000    0.000 (   1)    9.228 
	     CA:6 CA:11 :    6.671    6.671    6.671    0.000    0.000 (   1)    6.671 
	     CA:6 CA:12 :    8.042    8.042    8.042    0.000    0.000 (   1)    8.042 
	     CA:6 CA:13 :    7.742    7.742    7.742    0.000    0.000 (   1)    7.742 
	      CA:7 CA:9 :    8.772    8.772    8.772    0.000    0.000 (   1)    8.772 
	     CA:7 CA:10 :    9.572    9.572    9.572    0.000    0.000 (   1)    9.572 
	     CA:7 CA:11 :    7.564    7.564    7.564    0.000    0.000 (   1)    7.564 
	     CA:7 CA:12 :    9.784    9.784    9.784    0.000    0.000 (   1)    9.784 
	      CA:8 CA:9 :    9.932    9.932    9.932    0.000    0.000 (   1)    9.932 
	     CA:8 CA:10 :    9.488    9.488    9.488    0.000    0.000 (   1)    9.488 
	     CA:8 CA:11 :    6.587    6.587    6.587    0.000    0.000 (   1)    6.587 
	     CA:8 CA:12 :    7.798    7.798    7.798    0.000    0.000 (   1)    7.798 
	     CA:8 CA:13 :    9.048    9.048    9.048    0.000    0.000 (   1)    9.048 
	    CA:16 CA:17 :    9.317    9.520    9.820    0.207    0.092 (   5)    9.317    9.332    9.556    9.575    9.820 
	    CA:16 CA:18 :    9.289    9.502    9.852    0.269    0.134 (   4)    9.289    9.295    9.573    9.852 
	    CA:16 CA:19 :    8.430    8.657    8.890    0.203    0.091 (   5)    8.430    8.468    8.693    8.805    8.890 
	    CA:17 CA:17 :    9.829    9.913    9.997    0.119    0.084 (   2)    9.829    9.997 
	    CA:17 CA:18 :    6.530    6.871    7.129    0.260    0.116 (   5)    6.530    6.679    6.933    7.086    7.129 
	    CA:17 CA:19 :    5.464    5.618    5.739    0.113    0.051 (   5)    5.464    5.582    5.585    5.720    5.739 
	    CA:17 CA:20 :    9.035    9.175    9.352    0.137    0.061 (   5)    9.035    9.039    9.219    9.231    9.352 
	    CA:18 CA:18 :    6.592    6.747    7.021    0.173    0.077 (   5)    6.592    6.652    6.661    6.810    7.021 
	    CA:18 CA:19 :    7.220    7.527    7.767    0.245    0.109 (   5)    7.220    7.312    7.644    7.693    7.767 
	    CA:19 CA:20 :    9.745    9.822    9.989    0.114    0.057 (   4)    9.745    9.751    9.802    9.989 

Here we only examine distances from a single input set (`3P05`), but you can actually pass in multiple files and see a side-by-side comparison of the distance data to look for any interesting differences.


## <a name="FluctuationSpectrum"></a> FluctuationSpectrum

Generate the fluctuation spectra for a bilayer membrane from LAMMPS trajectories (using specialized LAMMPS file handling code).

This functionality can be useful in determining certain physical properties of a membrane system, for example the compressibility modulus or the bending modulus.

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/FluctuationSpectrum

	Usage: bin/FluctuationSpectum  traj=path  max_k=X max_l=X  head_type=X[,X,...] tail_type=X[,X,...]  gx=X gy=X  [delta_q=X] [scale=X] [which=X] [histogram=X] [out_prefix=X] [filter=X] [remap=X]

	Where:

	 - max_k, max_l : max integer wave numbers on x and y axes respectively.
	 - head_type, tail_type : LAMMPS atom types for head and terminal tail beads.
	 - gx, gy  : grid cell counts on x and y axes for midplane calculations.

	 - delta_q   : OPTIONAL resolution of output spectrum histogram (default: 0.05).
	 - scale     : OPTIONAL scaling for input->output length units (default: 1.0).
	 - which     : OPTIONAL setting for which monolayer: 'upper', 'lower', 'both', 'midplane' (default: 'both').
	 - histogram : OPTIONAL per-type histogram bin width in OUTPUT length units (ignored where <= 0.0).
	 - save_midplane : OPTIONAL flag to save the midplane coordinates as xyz (default: no midplane written).
	 - out_prefix: OPTIONAL output spectrum file prefix (default: 'spectrum').
	 - filter: OPTIONAL grid filter cell size, with contents of isolated cells ignored (default: no filtering).
	 - remap: OPTIONAL remap specifier, with most populated cell (via specified type) used recentre/wrap data (default: no remapping).
	 - start: OPTIONAL unit-based start frame in trajectory. Negative values ignored (default: -1).
	 - stop: OPTIONAL unit-based stop frame in trajectory. Negative values ignored (default: -1).
	 - save_raw: OPTIONAL flag to save raw (i.e. non-binned) spectral values. Negative values ignored (default: -1).

	Notes:

	Be careful if you use 'which=midplane', as the high frequency components of the spectrum (larger q values)
	will be limited by the resolution of the midplane grid (as specified by 'gx' and 'gy' parameters).

	The filtering assigns head group particles to a 3D grid with the specified size. Any cells with no neighbours
	have their contents ignored.
	MolecularUtilities $ 


## <a name="FluctuationSpectrum2"></a> FluctuationSpectrum2

This program is functionally identical to [FluctuationSpectrum](#FluctuationSpectrum), but instead of using the specialized LAMMPS file handling code it uses a more general representation of atoms and molecules etc. This in principle allows the straightforward addition of atomic filtering and support for other file formats etc, but it' also slightly slower and requires more memory. I've therefore included these two programs as separate entries for now.


## <a name="GenerateMembranes"></a> GenerateMembranes

Create monolayer and/or bilayer systems using user-defined lipid molecules and surface geometries.

This functionality can be useful when creating the initial system configuration for simulations of biological monolayer and bilayer interfaces.

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/GenerateMembranes

	Usage: bin/GenerateMembranes [bond_length=X] [leaflet_separation=X] [lipid=t1,t1,...[:b1,b2,...][:a1,a2,...], ...] [sphere=monolayer|bilayer:MOLSTRING:outer_r:APL, ...] [plane=bilayer|monolayer:MOLSTRING:N:APL, ...]

	Where:

	MOLSTRING : comma-separated lists of lipid types and relative proportions, separated by question mark.

	Examples:

	bin/GenerateMembranes lipid=1,2,3 sphere=bilayer:1:100:70
	bin/GenerateMembranes bond_length=7.5 leaflet_separation=7.5 lipid=1,2,3,3,3 lipid=1,2,3 sphere=bilayer:1,2?1,1:100:70
	bin/GenerateMembranes lipid=1,2,3,3,3 lipid=1,2,3 sphere=bilayer:1,2?10,12:100:70

	Notes:

	Molecule types are UNIT BASED and correspond to the order in which lipid definitions occurred.
	Molecule proportions are normalised internally, so they don't need to sum to 1 on the command line.
	If molecule proportions omitted, equal proportions used.

	MolecularUtilities $ 


## <a name="LammpsCombiner"></a> LammpsCombiner

Combine LAMMPS config files, automatically renumbering bonds and angles etc.

This functionality can be useful when you wish to combine separate LAMMPS config files to produce a cohesive configuration for simulation, ensuring all the atoms, molecules, bonds, angles, etc are correctly numbered and consistent.

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/LammpsCombiner

	Usage: bin/LAMMPSCombiner path:dx,dy,dz:adjust_topo_types path:dx,dy,dz:adjust_topo_types [check=X]

	Where:

	  path : LAMMPS config file
	  dx,dy,dz : offsets to translate system
	  adjust_topo_types : 1 where bond/angle types should be adjusted to follow previous data

	  check : OPTIONAL max bond length for sanity check of final data

	MolecularUtilities $ 


## <a name="LammpsToXYZ"></a> LammpsToXYZ

Read in a LAMMPS configuration or trajectory file, and write out a file in the simpler [XYZ](https://en.wikipedia.org/wiki/XYZ_file_format) format.

This functionality can be useful when you wish to visualize LAMMPS data, but the specific data file is either not commonly supported by visualization codes (e.g. LAMMPS configuration files) or you wish to reduce file sie and complexity for performance reasons (e.g. LAMMPS trajectory files).

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/FluctuationSpectrum2
	MolecularUtilities $ 


## <a name="SphereArbitrary"></a> SphereArbitrary

Generate a sphere using an arbitrary number of surface points

This functionality can be useful when ...

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/FluctuationSpectrum2
	MolecularUtilities $ 


## <a name="SphereBySubdivision"></a> SphereBySubdivision

Generate a sphere using poyhedral subdivision (cna provide more regular spacing that `SphereArbitrary`)

This functionality can be useful when ...

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/FluctuationSpectrum2
	MolecularUtilities $ 


## <a name="Superpose"></a> Superpose

Superpose arbitrary sets of PDB structures using filtered atom sets

This functionality can be useful when ...

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/FluctuationSpectrum2
	MolecularUtilities $ 


## <a name="UnwrapTrajectory"></a> UnwrapTrajectory

Unwrap molecules in a LAMMPs trajectory so they are not broken across periodic boundaries

This functionality can be useful when ...

Running the program with no command line options reveals a brief user guide:

	MolecularUtilities $ bin/FluctuationSpectrum2
	MolecularUtilities $ 
