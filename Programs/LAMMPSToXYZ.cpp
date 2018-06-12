/*
	Author: John Grime, The University of Chicago.
*/

#include <stdio.h>
#include <stdlib.h>

#include "../Util/Util.h"

//
// Convert a LAMMPS coordinate file into simple XYZ. This can be useful for visualizing large data sets
// in e.g. VMD.
//

using namespace Util;

int main( int argc, char** argv )
{
	std::vector< std::string > tokens;
	std::string type, in_path, out_path;
	FILE *in_f, *out_f;
	bool charges_present = false;

	LAMMPS::Config frame;

	if( argc < 3 )
	{
		printf( "\n" );
		printf( "Convert LAMMPS coordinate files into simple xyz format.\n" );
		printf( "\n" );
		printf( "Usage: %s  type  in_path  [out=path] [has_q]\n", argv[0] );
		printf( "\n" );
		printf( "Where:\n" );
		printf( "\n" );
		printf( "  type : either 'data' or 'traj'\n" );
		printf( "  has_q : optional switch to assume data files contain charges\n" );
		printf( "\n" );
		exit( -1 );
	}

	type = argv[1];
	in_path = argv[2];
	out_path = "output.xyz";

	for( int i=3; i<argc; i++ )
	{
		String::Tokenize( argv[i], tokens, "=" );		
		
		if( tokens.size() < 1 ) continue;
		if( tokens[0] == "has_q" ) { charges_present = true; continue; }
		
		if( tokens.size() < 2 ) continue;
		if( tokens[0] == "out" ) out_path = tokens[1];
	}

	if( (in_f=fopen(in_path.c_str(),"r")) == nullptr )
	{
		printf( "Unable to open input file '%s'.\n", in_path.c_str() );
		exit( -1 );
	}

	if( (out_f=fopen(out_path.c_str(),"w")) == nullptr )
	{
		printf( "Unable to open output file '%s'.\n", out_path.c_str() );
		exit( -1 );
	}

	int frame_no = 0;
	while( true )
	{
		frame_no++;

		//
		// Try to load some data
		//
		{
			int result;

			if( type == "data" ) result = LAMMPS::LoadData( in_f, frame, charges_present );
			else if( type == "traj" ) result = LAMMPS::LoadTrajectory( in_f, frame );
			else
			{
				printf( "Unknown input type '%s'\n", type.c_str() );
				return -1;
			}

			if( result == -1 || frame.molecules.size() < 1 )
			{
				printf( "Unable to load more data, so stopping here: frame %d\n", frame_no );
				break;
			}
			printf( "Frame %d, %d molecules, %d atoms.\n", frame_no, (int)frame.molecules.size(), result );
		}

		LAMMPS::SaveXYZ( out_f, frame );
	}

	fclose( in_f );
	fclose( out_f );

	return 0;
}

