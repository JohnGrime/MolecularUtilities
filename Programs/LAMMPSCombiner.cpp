/*
	Author: John Grime, The University of Chicago.
*/

#include "../Util/Util.h"


//
// This program tests the LAMMPS handler by adding two LAMMPS data files together.
//

using namespace Util;

void PrintUsage( const char* prog_name )
{
	printf( "\n" );
	printf( "Usage: %s path:dx,dy,dz:adjust_topo_types path:dx,dy,dz:adjust_topo_types [check=X]\n", prog_name );
	printf( "\n" );
	printf( "Where:\n" );
	printf( "\n" );
	printf( "  path : LAMMPS config file\n" );
	printf( "  dx,dy,dz : offsets to translate system\n" );
	printf( "  adjust_topo_types : 1 where bond/angle types should be adjusted to follow previous data\n" );
	printf( "\n" );
	printf( "  check : OPTIONAL max bond length for sanity check of final data\n" );
	printf( "\n" );
	exit( -1 );
}

int ParsePath( const char* path_str, std::string& path, double& x, double& y, double& z, int& adjust_topo_types )
{
	double delta[3] = { 0.0, 0.0, 0.0 };

	std::vector<std::string> tokens, subtoks;

	if( String::Tokenize( path_str, tokens, ":" ) < 3 ) return -1;

	path = tokens[0];
	adjust_topo_types = 0;

	//
	// Translation offsets
	//
	if( String::Tokenize( tokens[1], subtoks, "," ) < 3 ) return -1;
	for( int j=0; j<3; j++ )
	{
		if( String::ToReal( subtoks[j], delta[j] ) != String::ReturnValue::OK ) return -1;
	}

	x = delta[0];
	y = delta[1];
	z = delta[2];

	//
	// Topo type adjust
	//
	if( String::ToInteger( tokens[2], adjust_topo_types ) != String::ReturnValue::OK ) return -1;
	
	return 1;
}

void PrintInfo( FILE* f, const LAMMPS::Config& config )
{
	const auto& b = config.bounds;
	size_t N_atom_types = config.mass.size();
	size_t N_molecules = config.molecules.size();
	size_t N_atoms = 0;
	size_t N_bonds = 0;
	size_t N_angles = 0;
	std::map< int, int > bond_types, angle_types;

	//
	// Count some system components.
	//
	for( const auto& m : config.molecules )
	{
		N_atoms += m.atoms.size();
		
		for( const auto& b : m.bonds ) bond_types[ b.type ] = 1;
		N_bonds += m.bonds.size();

		for( const auto& a : m.angles ) angle_types[ a.type ] = 1;
		N_angles += m.angles.size();
	}

	fprintf( f, "bounds: %g %g %g -> %g %g %g : %g %g %g\n",
		b.minx, b.miny, b.minz,
		b.maxx, b.maxy, b.maxz,
		b.maxx-b.minx, b.maxy-b.miny, b.maxz-b.minz );

	// There's a dummy zero entry in the mass array to allod direct lookup of mass from LAMMPs type, so -1 for count.
	fprintf( f, "%d atom types, %d atoms, %d molecules\n", (int)N_atom_types-1, (int)N_atoms, (int)N_molecules );
	fprintf( f, "%d bond types, %d bonds\n", (int)bond_types.size(), (int)N_bonds );
	fprintf( f, "%d angle types, %d angles\n", (int)angle_types.size(), (int)N_angles );
}

int main( int argc, char** argv )
{
	char buffer[1024];
	FILE* f;
	LAMMPS::Config combined;

	int max_bond_type = 0, max_angle_type = 0;

	double sanity_check_tolerance = -1.0; // optional bond length check for sanity check!

	if( argc < 2 ) PrintUsage( argv[0] );

	for( int i=1; i<argc; i++ )
	{
		std::string path;
		LAMMPS::Config config;
		double dx,dy,dz;
		int adjust_topo_types = 0;
		bool shift_bounds = false;

		//
		// Optional sanity check info?
		//
		{
			std::vector< std::string > toks;
			if( String::Tokenize( argv[i], toks, "=" ) == -1 )
			{
				printf( "Unable split string '%s'\n", argv[i] );
				exit( -1 );				
			}
			if( (toks.size()>1) && (toks[0]=="check") )
			{
				if( String::ToReal( toks[1], sanity_check_tolerance ) != String::ReturnValue::OK )
				{
					printf( "Unable to convert bond check value '%s' into real number\n", toks[1].c_str() );
					exit( -1 );				
				}
				continue; // this was not an input sequence, so move on to next command line param
			}
		}

		if( ParsePath( argv[i], path, dx,dy,dz, adjust_topo_types ) == -1 )
		{
			printf( "Unable to parse input string '%s'\n", argv[i] );
			exit( -1 );
		}

		if( (f=fopen( path.c_str(),"r")) == nullptr )
		{
			printf( "Unable to open file '%s'\n", path.c_str() );
			exit( -1 );
		}
		LAMMPS::LoadData( f, config );
		fclose( f );

		//
		// Unwrap molecules by default; we're not tiling, so maybe not needed?
		//
		//config.Unwrap();

		//
		// Adjust topological types if needed!
		//
		if( adjust_topo_types > 0 )
		{
			//
			// Get maximum bind/angle types in current combined system
			// 
			for( const auto& mol : combined.molecules )
			{
				for( const auto x : mol.bonds )
				{
					if( x.type > max_bond_type ) max_bond_type = x.type;
				}

				for( const auto x : mol.angles )
				{
					if( x.type > max_angle_type ) max_angle_type = x.type;
				}
			}

			printf( "Adjusting bond/angle type offsets for '%s':\n", path.c_str() );
			printf( "  using bond type offset = %d\n", max_bond_type );
			printf( "  using angle type offset = %d\n", max_angle_type );

			//
			// Modify bond/angle types in config we're adding to the combined system
			//
			for( auto& mol: config.molecules )
			{
				for( auto& x : mol.bonds )  x.type += max_bond_type;
				for( auto& x : mol.angles ) x.type += max_angle_type;
			}
		}

		//
		// Debug - write input data files as xyz
		//
		if( true )
		{
			sprintf( buffer, "config_input.%d.xyz", i );
			if( (f=fopen( buffer,"w")) == nullptr )
			{
				printf( "Unable to open file '%s'\n", buffer );
				exit( -1 );
			}
			LAMMPS::SaveXYZ( f, config );
			fclose( f );
		}

		config.Translate( dx, dy, dz, shift_bounds );

		printf( "\n" );
		printf( "Source: '%s'\n", path.c_str() );
		printf( "translate: %g %g %g\n", dx, dy, dz );
		PrintInfo( stdout, config );
		printf( "\n" );

		if( i == 1 ) combined = config;
		else combined += config;
	}

	printf( "\n" );
	printf( "Combined data:\n" );
	PrintInfo( stdout, combined );
	printf( "\n" );

	//
	// Sanity check topology?
	//
	if( sanity_check_tolerance > 0.0 )
	{
		int bad_mol_i = 0;
		bool sane = combined.SanityCheckTopology( sanity_check_tolerance, bad_mol_i );

		if( sane == false )
		{
			const auto& m = combined.molecules[bad_mol_i];

			//			
			// Print xyz of problem molecule!
			//
			FILE* f = fopen( "problem_molecule.xyz", "w" );
			fprintf( f, "%d\n", (int)m.atoms.size() );
			fprintf( f, "Problem molecule from topology sanity check, index %d\n", bad_mol_i );
			for( const auto& a : m.atoms ) fprintf( f, "%d %f %f %f\n", a.type, a.x, a.y, a.z );
			fclose( f );
		}
	}
	
	//
	// Save combined system, as LAMMPS data and xyz.
	//
	{
		sprintf( buffer, "combined.xyz" );
		if( (f=fopen( buffer,"w")) == nullptr )
		{
			printf( "Unable to open file '%s'\n", buffer );
			exit( -1 );
		}
		LAMMPS::SaveXYZ( f, combined );
		fclose( f );


		sprintf( buffer, "combined.lammps_config" );
		if( (f=fopen( buffer,"w")) == nullptr )
		{
			printf( "Unable to open file '%s'\n", buffer );
			exit( -1 );
		}
		LAMMPS::SaveData( f, combined );
		fclose( f );
	}

	return -1;
}
