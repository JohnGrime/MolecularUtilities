/*
	Author: John Grime, The University of Chicago.
*/

#include <stdio.h>
#include <string>
#include <vector>

#include "../Util/Util.h"

//
// Align PDB file to specified axes, via the calculated axes of a specified set of
// atoms in the PDB.
//

using namespace Util;

//
// Simple utility structure for eigen info.
//
struct EigenInfo
{
	double centroid[3], eigenvals[3], eigenvecs[3][3];
	
	EigenInfo()
	{
		Identity();
	}
	
	//
	// Cartesian axes in order x,y,z, centroid = origin.
	//
	void Identity()
	{
		for( int i=0; i<3; i++ )
		{
			centroid[i] = 0.0;
			for( int j=0; j<3; j++ ) eigenvecs[i][j] = 0.0;

			eigenvals[i] = 1.0; // eigenvals of identity matrix are 1.0
			eigenvecs[i][i] = 1.0; // identity matrix = 1.0 along diagonal
		}
	}
	
	void Populate( std::vector<double> &xyz, bool swap_orders = true )
	{
		EigenUtil<3>::Eigen( xyz, centroid, eigenvals, eigenvecs, swap_orders );
	}
	
	void Print( FILE *f ) const
	{
		fprintf( f, "Centroid:\n" );
		fprintf( f, "\t %12.3f %12.3f %12.3f\n", centroid[0], centroid[1], centroid[2] );
		fprintf( f, "Eigen values / vectors:\n" );
		fprintf( f, "\t %12.3f : %12.3f %12.3f %12.3f\n", eigenvals[0], eigenvecs[0][0], eigenvecs[0][1], eigenvecs[0][2] );
		fprintf( f, "\t %12.3f : %12.3f %12.3f %12.3f\n", eigenvals[1], eigenvecs[1][0], eigenvecs[1][1], eigenvecs[1][2] );
		fprintf( f, "\t %12.3f : %12.3f %12.3f %12.3f\n", eigenvals[2], eigenvecs[2][0], eigenvecs[2][1], eigenvecs[2][2] );
	}
	
	void PrintPDB( FILE *f, double scale[3] ) const
	{
		Molecule::Atom a;
		const char *names[4] = { "C", "O", "CLA", "N" }; // O,CLA,N = red,green,blue in VMD's CPK representation?
		const char *serials[4] = { "1", "2", "3", "4" };
	
		a.attr["recname"] = "ATOM";

		for( int i=0; i<4; i++ )
		{
			a.attr["name"] = names[i];
			a.attr["serial"] = serials[i];

			a.x = centroid[0];
			a.y = centroid[1];
			a.z = centroid[2];

			if( i > 0 )
			{
				a.x += eigenvecs[i-1][0] * scale[0];
				a.y += eigenvecs[i-1][1] * scale[1];
				a.z += eigenvecs[i-1][2] * scale[2];
			}
			
			Molecule::PDB::Print( f, a );
			fprintf( f, "\n" );
		}
		fprintf( f, "TER\n" );
		fprintf( f, "%6.6s%5d%5d\n", "CONECT", 1, 2 );
		fprintf( f, "%6.6s%5d%5d\n", "CONECT", 1, 3 );
		fprintf( f, "%6.6s%5d%5d\n", "CONECT", 1, 4 );
		fprintf( f, "END\n" );
	}
	
};


//
// Read desired "world" axis ordering from command line.
// This is the reference frame we're superposing the molecular axes onto.
//
void GetWorldReferenceFrame( int argc, char **argv, EigenInfo &ei )
{
	const char *names = "xyz";
	int defined[3] = { 0, 0, 0 }, which = 0;
	
	ei.Identity();
	
	//
	// Check parameters for axis ordering
	//
	for( int i=2; i<argc; i++ )
	{
		double val = 1.0;
		int index = 0;
		char c = argv[i][0];
		
		if( c == '-' )
		{
			val *= -1.0;
			c = argv[i][1];
		}
		else if( c == '+' )
		{
			val *= +1.0;
			c = argv[i][1];
		}
		
		index = -1;
		
		if( c == 'x' ) index = 0;
		else if( c == 'y' ) index = 1;
		else if( c == 'z' ) index = 2;
		
		if( index == -1 ) continue;
		
		if( defined[index] == 1 )
		{
			printf( "Error: axis %c multiply defined!\n", names[index] );
			exit( -1 );
		}
		defined[index] = 1;
		
		ei.eigenvecs[which][0] = 0.0;
		ei.eigenvecs[which][1] = 0.0;
		ei.eigenvecs[which][2] = 0.0;
		
		ei.eigenvecs[which][index] = val;
		
		which++;
	}
	
	//
	// Check all axes defined!
	//
	for( int i=0; i<3; i++ )
	{
		if( defined[i] == 1 ) continue;
		printf( "Error: order of axis %c is not described!\n", names[i] );
		exit( -1 );
	}
	
}

//
// Important note: we have a "chiral problem" if we treat the reference frames as 4 points to superpose. The final
// (least significant) axis could point "into" or "out of" the plane described by the origin and the first two axes.
//
// To avoid this issue in a simple way, we use only the first 3 points (i.e. the origin and two major axes) to align.
//
void AlignReferenceFrames( const EigenInfo &target, const EigenInfo &current, double M[3][3], double T[3] )
{
	std::vector<double> xyz, target_xyz;
	
	// We're only using 3 points!
	xyz.resize( 3*3 );
	target_xyz.resize( 3*3 );

	// This obviously ignores the least significant axis ...
	for( int i=0; i<3; i++ )
	{
		xyz[i*3 +0] = current.centroid[0];
		xyz[i*3 +1] = current.centroid[1];
		xyz[i*3 +2] = current.centroid[2];
		
		target_xyz[i*3 +0] = target.centroid[0];
		target_xyz[i*3 +1] = target.centroid[1];
		target_xyz[i*3 +2] = target.centroid[2];

		if( i > 0 )
		{
			xyz[i*3 +0] += current.eigenvecs[i-1][0];
			xyz[i*3 +1] += current.eigenvecs[i-1][1];
			xyz[i*3 +2] += current.eigenvecs[i-1][2];

			target_xyz[i*3 +0] += target.eigenvecs[i-1][0];
			target_xyz[i*3 +1] += target.eigenvecs[i-1][1];
			target_xyz[i*3 +2] += target.eigenvecs[i-1][2];
		}
	}

	Superposer::Calculate( target_xyz, xyz, M, T );
}


int main( int argc, char **argv )
{
	FILE *f;
	std::vector<double> xyz;
	
	EigenInfo world, molecular;

	Molecule::MolecularSystem ms;
	Molecule::AttributeFilter filter;
	Molecule::Molecules filtered_molecules;

	if( argc < 2 )
	{
		printf( "\n" );
		printf( "Usage:\n" );
		printf( "\n" );
		printf( "%s <input PDB> [-]x|y|z [-]x|y|z [-]x|y|z [ \"key1=val1,val2-val3;key2=val4\" ]\n", argv[0] );
		printf( "\n" );
		printf( "Cartesian axis order overrides the defauly x,y,z. ALL axes must be\n" );
		printf( "specified, and a leading '-' indicates the axis is reflected.\n" );
		printf( "\n" );
		printf( "An optional set of filters can be included for the generation of the molecular axes.\n" );
		printf( "\n" );
		printf( "Output is 'mol_axes.pdb' (molecular axes) and 'aligned.pdb' (PDB aligned to specified axes).\n" );
		printf( "Note that 'mol_axes.pdb' refers to the original PDB file, not the aligned structure.\n" );
		printf( "\n" );
		printf( "Example:\n" );
		printf( "\n" );
		printf( "%s blah.pdb y -z x \"name=CA;resSeq=1,3,12-14,20-50\"\n", argv[0] );
		printf( "\n" );
		printf( "This aligns the major axis of the specified atoms in molecule in blah.pdb to the y axis,\n" );
		printf( "the secondary axis to the flipped z axis (i.e. {0,0,-1}) and the minor axis to the Cartesian\n" );
		printf( "x axis.\n" );
		printf( "\n" );
		exit( -1 );
	}
	
	//
	// To what world reference frame are we aligning the molecular axes?
	//
	GetWorldReferenceFrame( argc, argv, world );

	//
	// Any filters specified?
	//
	if( argc > 5 )
	{
		std::vector<std::string> tokens;
		
		String::Tokenize( argv[5], tokens, ";" );
		for( size_t i=0; i<tokens.size(); i++ )
		{
			filter.AddFilter( tokens[i].c_str(), "=", ",", "-" );
		}
	}

	//
	// Load the input data; assumes PDB.
	//
	if( (f=fopen(argv[1],"r")) == NULL )
	{
		printf( "Unable to open file '%s'\n", argv[1] );
		exit( -1 );
	}
	Molecule::PDB::Load( f, ms );
	fclose( f );
	
	//
	// Calculate eigen info from molecule coordinates
	//
	int n_total_atoms = 0;
	for( auto& m : ms.molecules ) n_total_atoms += (int)m.size();
	int n_filtered_atoms = filter.Filter( ms.molecules, filtered_molecules );
	printf( "Atom filtering: %d remaining of %d.\n", n_filtered_atoms, n_total_atoms );
	if( n_filtered_atoms < 3 )
	{
		printf( "Too few atoms passed the filter.\n" );
		exit( -1 );
	}
	Molecule::Coords::Get( filtered_molecules, xyz );
	molecular.Populate( xyz );

	// Bounding box of initial system
	{
		double minx,maxx, miny,maxy, minz,maxz;
		Molecule::Coords::BoundingBox( ms.molecules, minx,maxx, miny,maxy, minz,maxz );
		double Lx = maxx-minx;
		double Ly = maxy-miny;
		double Lz = maxz-minz;
		double V = Lx*Ly*Lz;
		printf( "Bounding box of final system: %g %g %g -> %g %g %g ( %g %g %g , volume %g )\n",
			minx,miny,minz, maxx,maxy,maxz, Lx,Ly,Lz, V  );
	}

	molecular.Print( stdout );
	world.Print( stdout );

	//
	// Figure out appropriate lengths for the molecular axes, for when
	// we write axes structure to a PDB file.
	//
	double scale[3];
	for( size_t i=0; i<xyz.size(); i+=3 )
	{
		double x = xyz[i +0];
		double y = xyz[i +1];
		double z = xyz[i +2];
		
		if( i == 0 )
		{
			scale[0] = x;
			scale[1] = y;
			scale[2] = z;
		}
		else
		{
			if( x > scale[0] ) scale[0] = x;
			if( y > scale[1] ) scale[1] = y;
			if( z > scale[2] ) scale[2] = z;
		}
	}
	for( int i=0; i<3; i++ ) scale[i] = (scale[i]-molecular.centroid[i])*1.1;
	
	f = fopen( "mol_axes.pdb", "w" );
	molecular.PrintPDB( f, scale );
	fclose( f );

	//
	// Align molecular axes onto specified world axes.
	//
	{
		double M[3][3], T[3];
		AlignReferenceFrames( world, molecular, M, T );
		Molecule::Coords::Get( ms.molecules, xyz );
		Superposer::Apply( xyz, xyz, M, T );
		Molecule::Coords::Set( ms.molecules, xyz );
	}

	//
	// Comment this out if you want the molecular axes to lie on top of
	// the world axes, rather than just being parallel to them after
	// translating the molecule's centre of geometry to the origin.
	//
	Molecule::Coords::ResetCOG( ms.molecules );

	//
	// Calculate bounding box of final system, then translate for symmetry
	// on each axis.
	//
	{
		double minx,maxx, miny,maxy, minz,maxz;
		double Lx, Ly, Lz, V;

		Molecule::Coords::BoundingBox( ms.molecules, minx,maxx, miny,maxy, minz,maxz );
		Lx = maxx-minx;
		Ly = maxy-miny;
		Lz = maxz-minz;
		V = Lx*Ly*Lz;

		printf( "Bounding box of rotated system: %g %g %g -> %g %g %g ( %g %g %g , volume %g )\n",
			minx,miny,minz, maxx,maxy,maxz, Lx,Ly,Lz, V  );

		for( auto& mol : ms.molecules )
		{
			for( auto& a : mol )
			{
				a.x -= (minx + Lx/2);
				a.y -= (miny + Ly/2);
				a.z -= (minz + Lz/2);
			}
		}

		Molecule::Coords::BoundingBox( ms.molecules, minx,maxx, miny,maxy, minz,maxz );
		Lx = maxx-minx;
		Ly = maxy-miny;
		Lz = maxz-minz;
		V = Lx*Ly*Lz;

		printf( "Bounding box of recentred system: %g %g %g -> %g %g %g ( %g %g %g , volume %g )\n",
			minx,miny,minz, maxx,maxy,maxz, Lx,Ly,Lz, V  );
	}

	f = fopen( "aligned.pdb", "w" );
	Molecule::PDB::Print( f, ms );
	fclose( f );

}
