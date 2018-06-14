/*
	Author: John Grime, The University of Chicago.
*/

#include "../Util/Util.h"

//
// Generate a spherical surface using polyhedral subdivision.
//

using namespace Util;

void print_usage( const char *prog )
{
	printf( "\n" );
	printf( "Usage: %s  radius n_subdiv aname rname\n", prog );
	printf( "\n" );
	printf( "Where:\n" );
	printf( "  - radius : radius of sphere, in Angstrom\n" );
	printf( "  - n_subdiv: number of polyhedral subdivisions to perform (larger number = more points)\n" );
	printf( "  - aname : atom name for the surface beads\n" );
	printf( "  - rname : residue name for the surface beads\n" );
	printf( "\n" );
	printf( "Notes:\n" );
	printf( "  - The number of points generated increases RAPIDLY with n_subdiv.\n" );
	printf( "  - n_subdiv between 2 and 5 is typically a good starting point.\n" );
	printf( "\n" );
}

int main( int argc, char **argv )
{
	const char* ATOM_format = "%-6.6s%5.5s %4.4s%1.1s%3.3s %1.1s%4d%1.1s   %8.2f%8.2f%8.2f%6.6s%6.6s          %2.2s%2.2s\n";

	double radius = 1.0;
	int n_subdiv = 3;
	
	const char *surface_aname = "?";
	const char *surface_rname = "?";
		
	Subdivider subdiv;
	
	if( argc < 5 )
	{
		print_usage( argv[0] );
		exit( -1 );
	}
	
	radius = atof( argv[1] );
	n_subdiv = atoi( argv[2] );
	
	surface_aname = argv[3];
	surface_rname = argv[4];
		
	subdiv.Init( Subdivider::InitType::Icosahedron );
	for( int i=0; i<n_subdiv; i++ ) subdiv.Subdivide();


	//
	// Outer beads, form surface of sphere
	//
	for( int i=0; i<subdiv.n_vertices; i++ )
	{
		const char* recname = "ATOM";
		char serial[16];
		const char* name = surface_aname;
		const char* altLoc = "";
		const char* resName = surface_rname;
		const char* chainID = "A";
		int resSeq = 1;
		const char* iCode = "";
		double x = subdiv.vertices[i*3+0] * radius;
		double y = subdiv.vertices[i*3+1] * radius;
		double z = subdiv.vertices[i*3+2] * radius;
		const char* occupancy = "";
		const char* tempFactor = "";
		const char* element = "";
		const char* charge = "";

		sprintf( serial, "%d", i+1 ); // serial

		printf(
			ATOM_format,
			recname,
			serial,
			name,
			altLoc,
			resName,
			chainID,
			resSeq,
			iCode,
			x,y,z,
			occupancy,
			tempFactor,
			element,
			charge );
	}

	printf( "TER\n" );
}
