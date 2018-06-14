/*
	Author: John Grime, The University of Chicago.
*/

#include <stdio.h>
#include <stdlib.h>

#include <map>

#include "../Util/Util.h"

using namespace Util;

//
// Unwraps molecules in a LAMMPS trajectory so that they do not
// "break" across periodic boundaries.
//

int main( int argc, char** argv )
{
	FILE *in_traj, *out_traj;

	LAMMPS::Config config;

	if( argc < 3 )
	{
		printf( "\n" );
		printf( "Usage: %s in_traj out_traj\n", argv[0] );
		printf( "\n" );
		printf( "Where:\n" );
		printf( "  - in_traj : input LAMMPS trajectory\n" );
		printf( "  - out_traj : output LAMMPS trajectory\n" );
		printf( "\n" );
		exit( -1 );
	}

	if( (in_traj=fopen(argv[1],"r") ) == nullptr )
	{
		printf( "Unable to open input trajectory '%s'.\n", argv[1] );
		exit( -1 );
	}
	if( (out_traj=fopen(argv[2],"w") ) == nullptr )
	{
		printf( "Unable to open output trajectory '%s'.\n", argv[2] );
		exit( -1 );
	}

	int frame_no = 0, timestep = 0;
	while( true )
	{
		frame_no++;

		if( frame_no%10 == 0 ) printf( "%d\n", frame_no );

		//
		// Try to load a trajectory frame
		//
		{
			int result = Util::LAMMPS::LoadTrajectory( in_traj, config, &timestep );
			if( (result==-1) || feof(in_traj) )
			{
				printf( "Unable to load trajectory frame. Stopping here, frame %d\n", frame_no );
				break;
			}
		}

		//
		// Wrap
		//
		{
			config.Unwrap( config.PBC );
		}

		//
		// Write
		//
		{
			Util::LAMMPS::SaveTrajectory( out_traj, config, timestep );
		}
	}

	fclose( in_traj );
	fclose( out_traj );

	return 0;
}

