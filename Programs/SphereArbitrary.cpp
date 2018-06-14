/*
	Author: John Grime, The University of Chicago.
*/

#include "../Util/Util.h"

using namespace Util;

//
// Create a sphere composed of near-equally spaced points.
// Number of points may not actually be exactly the number
// specified on the command line!
//

template <typename T>
void write_xyz( FILE *f, std::vector<T> &xyz )
{
	int N = (int)(xyz.size()/3);
	
	fprintf( f, "%d\n", N );
	fprintf( f, "%s\n", "Test" );
	
	for( size_t i=0; i<xyz.size(); i+= 3 )
	{
		fprintf( f, "%s %f %f %f\n", "S", xyz[i+0], xyz[i+1], xyz[i+2] );
	}
}


void print_usage( const char *prog )
{
	printf( "\n" );
	printf( "Usage:\n" );
	printf( "\n" );
	printf( "%s  radius target_N aname rname [TER]\n", prog );
	printf( "\n" );
	printf( "Where:\n" );
	printf( "  - radius : radius of sphere, in Angstrom\n" );
	printf( "  - target_N : desired number of points on sphere surface\n" );
	printf( "  - aname : atom name for the surface beads\n" );
	printf( "  - rname : residue name for the surface beads\n" );
	printf( "  - TER   : OPTIONAL flag to insert TER lines after every ATOM\n" );
	printf( "\n" );
	printf( "Example:\n" );
	printf( "\n" );
	printf( " %s  100.0  350  PNT  SPH > my_sphere.pdb\n", prog );
	printf( "\n" );
	printf( "Note that the number of points may not actually be target_N; it should be close, though!\n" );
	printf( "\n" );
}


int main( int argc, char **argv )
{
	std::vector<double> xyz;
	Molecule::Atom a;
	char buffer[128];
	
	double radius;
	int target_N;
	bool TER_per_atom = false;


	//
	// Parse command line arguments
	//

	if( argc < 5 )
	{
		print_usage( argv[0] );
		exit( -1 );
	}
	
	if( (String::ToReal(argv[1],radius)!=String::ReturnValue::OK) || (radius<=0.0) )
	{
		printf( "Bad radius '%s'!\n", argv[1] );
		exit( -1 );
	}
	if( (String::ToInteger(argv[2],target_N)!=String::ReturnValue::OK) || (target_N<=0) )
	{
		printf( "Bad target_N '%s'!\n", argv[2] );
		exit( -1 );
	}
	
	const char *aname = argv[3];
	const char *rname = argv[4];

	TER_per_atom = (argc>5) ? (true) : (false);


	//
	// Generate sphere surface coordinates
	//

	Geometry::GetEquispacedSpherePoints( target_N, xyz, radius );

	//
	// Set up default PDB atom info
	//
	
	a.attr["recname"] = "ATOM";
	a.attr["name"] = aname;
	a.attr["resName"] = rname;
	a.attr["chainID"] = "A";
	a.attr["resSeq"] = "1";

	
	//
	// Print PDB to stdout, replacing atom indices/coordinates as we go
	//
	
	for( size_t i=0, max_i=xyz.size()/3; i<max_i; i++ )
	{
		sprintf( buffer, "%d", (int)i+1 );
		a.attr["serial"] = buffer;
		a.x = xyz[i*3 +0];
		a.y = xyz[i*3 +1];
		a.z = xyz[i*3 +2];
		Molecule::PDB::Print( stdout, a );
		printf( "\n" );
		if( TER_per_atom ) printf( "TER\n" );
	}
	printf( "TER\n" );
}

