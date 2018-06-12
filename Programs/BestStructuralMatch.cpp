/*
	Author: John Grime, The University of Chicago.
*/

#include "../Util/Util.h"

using namespace Util;

//
// Given an input reference set of PDB structures, a set size and a bunch of possible data, find the 
// best match in the data for the reference set on the basis of the RMSD (having applied appropriate
// filters)
//

int main( int argc, char **argv )
{
	const char *reference_path, *data_path;
	const char *best_path = "best_match.pdb";
	int set_size, N_loaded;
	
	FILE *f;

	Molecule::MolecularSystem reference_set, data_set;
	Molecule::AttributeFilter filter;
	
	Molecule::Molecules filtered_reference, filtered_data;
	std::vector<double> target_xyz, xyz;

	if( argc < 4 )
	{
		printf( "Usage: %s reference.pdb set_size data.pdb [filter_name=v1,v2,...] [filter_name=v1,v2,...]\n", argv[0] );
		printf( "Where:\n" );
		printf( "  - reference.pdb : the reference structure(s) to compare to\n" );
		printf( "  - set_size : number of consecutive PDB molecules in the input file to group\n" );
		printf( "  - data.pdb : the potential candidate structures to compare to reference.pdb\n" );
		printf( "  - OPTIONAL filters : used in the superposition/RMSD comparison\n" );
		printf( "Notes:\n" );
		printf( "  Consecutive set_size entries from reference.pdb and data.pdb are grouped for superposition/RMSD calculation.\n" );
		printf( "  Closest match to reference in data set saved to '%s', with some info in REMARK lines.\n", best_path );
		printf( "  ALL atom data is written in '%s', not just the atoms used in the superposition/RMSD calculation.\n", best_path );
		exit( -1 );
	}

	reference_path = argv[1];
	set_size = atoi( argv[2] );
	data_path = argv[3];
	
	//
	// Filters defined? Add them.
	// Assumes filters defined using one or more strings as: "name=val,val-val,...[:name=val,val,...]"
	//
	if( argc > 4 )
	{
		std::vector<std::string> tokens;
		for( int i=4; i<argc; i++ )
		{
			String::Tokenize( argv[i], tokens, ":" );
			for( size_t j=0; j<tokens.size(); j++ )
			{
				filter.AddFilter( tokens[j].c_str(), "=", ",", "-" ); // key and value separated by =, values separated, by comma, ranges use dash.
			}
		}
		//filter.Print();
	}
	
	//
	// Load input data.
	//
	if( (f=fopen(reference_path,"r")) == NULL )
	{
		printf( "Unable to open file '%s'\n", reference_path );
		exit( -1 );
	}
	N_loaded = Molecule::PDB::Load( f, reference_set, set_size );
	fclose( f );
	
	if( N_loaded < set_size )
	{
		printf( "Unable to load %d structures from '%s'; only managed %d!\n", set_size, reference_path, N_loaded );
		exit( -1 );
	}
	
	//
	// Walk through structures.
	//
	if( (f=fopen(data_path,"r")) == NULL )
	{
		printf( "Unable to open file '%s'\n", data_path );
		exit( -1 );
	}
	
	int set_number = 0;
	
	int best_set = 0, best_N_common = 0;
	double best_RMSD = 0.0;
	Molecule::Molecules best_molecules;
	
	while( true )
	{
		double rmsd[2], M[3][3], T[3];
		
		set_number++;
		
		N_loaded = Molecule::PDB::Load( f, data_set, set_size );
		if( N_loaded < set_size )
		{
			printf( "Unable to load %d structures from '%s'; only managed %d! Assuming EOF and stopping here.\n", set_size, data_path, N_loaded );
			break;
		}

		Molecule::Coords::Get( reference_set.molecules, target_xyz );
		Molecule::Coords::Get( data_set.molecules, xyz );
		
		int N_common = filter.Filter( reference_set.molecules, data_set.molecules, filtered_reference, filtered_data );

		Molecule::Coords::Get( filtered_reference, target_xyz );
		Molecule::Coords::Get( filtered_data, xyz );
		
		rmsd[0] = Superposer::RMSD( target_xyz, xyz ); // pre-superposition
		
		Superposer::Calculate( target_xyz, xyz, M, T );
		Superposer::Apply( xyz, xyz, M, T );

		rmsd[1] = Superposer::RMSD( target_xyz, xyz ); // post-superposition
				
		if( set_number == 1 || rmsd[1] < best_RMSD )
		{
			//
			// New best structure! Superpose entire atom set, and store.
			//
			Molecule::Coords::Get( data_set.molecules, xyz );
			Superposer::Apply( xyz, xyz, M, T );
			Molecule::Coords::Set( data_set.molecules, xyz );

			best_RMSD = rmsd[1];
			best_set = set_number;
			best_N_common = N_common;
			best_molecules = data_set.molecules;

			printf( "Set %5d : %d common filtered atoms : RMSD %12.3f => %12.3f <- best so far\n", set_number, N_common, rmsd[0], rmsd[1] );
		}
		else
		{
			printf( "Set %5d : %d common filtered atoms : RMSD %12.3f => %12.3f\n", set_number, N_common, rmsd[0], rmsd[1] );
		}
	}
	fclose( f );
	
	//
	// Save best set.
	//
	if( (f=fopen(best_path,"w")) == NULL )
	{
		printf( "Unable to open file '%s'\n", best_path );
		exit( -1 );
	}
	fprintf( f, "REMARK\n" );
	fprintf( f, "REMARK  Target:     '%s'\n", reference_path );
	fprintf( f, "REMARK  Structures: '%s'\n", data_path );
	if( argc > 4 )
	{
		for( int i=4; i<argc; i++ )
		{
			fprintf( f, "REMARK  %8.8s  %s\n", (i==4) ? ("Filters:") : (""), argv[i] );
		}
	}
	else
	{
		fprintf( f, "REMARK  No filtering was used\n" );
	}
	fprintf( f, "REMARK  Set size: %d\n", set_size );
	fprintf( f, "REMARK  Best set: %d (RMSD of %g using %d filtered common atoms)\n", best_set, best_RMSD, best_N_common );
	fprintf( f, "REMARK\n" );
	Molecule::PDB::Print( f, best_molecules );
	fclose( f );
	
	return 0;
}
