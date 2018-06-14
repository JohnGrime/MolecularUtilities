/*
	Author: John Grime, The University of Chicago.
*/

#include <stdio.h>

#include <string>
#include <map>
#include <vector>

#include "../Util/Util.h"

using namespace Util;

/*
	This program allows the superposition of sets of molecules by defining target data and the data to superpose.
	
	There are two major use cases for something like this:
	
	1. You want to superpose a single example structure onto a larger group of structures.
	
	This is handy when you're converting a detailed system into a CG representation, or vice versa.
	
	2. You want to superpose multiple copies of a structure onto a single exemplar frame of reference.
	
	This is handy when you're showing the general conformational "spread" of a structure, and also for
	calculating average average deviations and positions etc from some set of data.
	
	This program therefore allows you to specify a "target" PDB file and a "structures" PDB file. The contents of
	the structures file is superposed onto the target, but you have a great deal of control over how this happens.
	
	Run the program without any parameters for a description of how to use this utility.
*/


//
// Description of a particular superposition entry:
//  - indices of source molecules to superpose ("superpose")
//  - indices of target molecules to superpose onto ("onto")
//  - molecules in source to apply the resultant superposition to ("apply_to")
//  - molecule filter for the superpose molecules ("mf")
//
struct Superposition
{
	std::vector< int > superpose, onto, apply_to; // mol indices
	Molecule::AttributeFilter mf;
	
	//
	// Expect str == "key=value[,value,value,...]:key=value[,value,value,...]:..."
	//
	int parse( const char *str )
	{
		std::vector< std::string > elements, tokens;
		
		superpose.clear();
		onto.clear();
		apply_to.clear();
		mf.filters.clear();

		if( String::Tokenize( str, elements, ":" ) < 3 )
		{
			fprintf( stderr, "Bad superpose definition '%s'\n", str );
			return -1;
		}

		//
		// Parse vectors of key=value[,value,value,...] elements
		//
		for( const auto& element : elements )
		{
			//
			// tokens[0] = key
			// tokens[1] = value[,value,value,...]
			//
			if( String::Tokenize( element.c_str(), tokens, "=" ) < 1 )
			{
				fprintf( stderr, "Bad superposition element '%s'!\n", element.c_str() );
				return -1;
			}
			
			const auto& key = tokens[0];
			
			if( key == "superpose" || key == "onto" || key == "apply_to" )
			{
				//
				// Defined key! Take appropriate action ...
				//
				
				if( tokens.size() < 2 )
				{
					fprintf( stderr, "Bad superposition element '%s'!\n", element.c_str() );
					return -1;
				}
				
				if( key == "superpose" )
				{
					if( String::ToIntegers( tokens[1].c_str(), superpose, ",", "-" ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "Bad superposition element '%s'!\n", element.c_str() );
						return -1;
					}
				}
				else if( key == "onto" )
				{
					if( String::ToIntegers( tokens[1].c_str(), onto, ",", "-" ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "Bad superposition element '%s'!\n", element.c_str() );
						return -1;
					}
				}
				else if( key == "apply_to" )
				{
					if( String::ToIntegers( tokens[1].c_str(), apply_to, ",", "-" ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "Bad superposition element '%s'!\n", element.c_str() );
						return -1;
					}
				}
			}
			else
			{
				std::vector< std::string > values;
				
				//
				// Assume this key=value[,value,...] elements describes a superposition filter,
				// with filter values separated by a comma.
				//
				
				// Only one token, so add "empty" filter
				if( tokens.size() < 2 )
				{
					mf.AddFilter( key.c_str(), NULL );
					continue;
				}
				
				//
				// Add filter values under specified key.
				// Special case: look for range separator, "-". If found, assume we're
				// trying to define a numeric range for the values.
				//
				String::Tokenize( tokens[1].c_str(), values, "," );
				for( const auto& value : values )
				{
					bool is_range = false;
					
					// Check for range character in the current value
					for( size_t i=0; i<value.size(); i++ )
					{
						if( value[i] == '-' )
						{
							is_range = true;
							break;
						}
					}
					
					if( is_range == true )
					{
						//
						// This is ugly. We try to get the list of integers in the range, and then
						// convert them back to strings to add to the filter.
						//
						char buffer[1024];
						std::vector<int> numbers;
						if( String::ToIntegers( value.c_str(), numbers, ",", "-" ) != String::ReturnValue::OK )
						{
							fprintf( stderr, "Bad range token in '%s' from superposition element '%s'!\n", value.c_str(), element.c_str() );
							return -1;
						}
						for( size_t i=0; i<numbers.size(); i++ )
						{
							sprintf( buffer, "%d", numbers[i] );
							mf.AddFilter( key.c_str(), buffer );
						}
					}
					else
					{
						// Value not assumed to describe a range, add directly.
						mf.AddFilter( key.c_str(), value.c_str() );
					}
				}
			}
		}
		return 1;
	}

	void print() const
	{
		printf( "Superposition:\n" );

		printf( "\t superpose: " );
		for( size_t i=0; i<superpose.size(); i++ ) printf( "%d ", superpose[i] );
		printf( "\n" );

		printf( "\t onto: " );
		for( size_t i=0; i<onto.size(); i++ ) printf( "%d ", onto[i] );
		printf( "\n" );

		printf( "\t apply_to: " );
		for( size_t i=0; i<apply_to.size(); i++ ) printf( "%d ", apply_to[i] );
		printf( "\n" );

		mf.Print();
	}
};

/*
	Parameters for superposition.
	
	Target and structures both have a path, a set size (default = 1) and
	a flag controlling whether we advance through the molecular data (default
	is to advance in both data sets.
	
	You're probably not going to want to advance through both data sets, so
	make sure you pay attention to the command-line parameters!
*/
struct Parameters
{
	struct Input
	{
		std::string path;
		int set_size, advance;
		
		Input()
		{
			path = "";
			set_size = 1;
			advance = 1;
		}

		//
		// Expects str == "path[:set_size][:noadvance]"
		// WILL NOT MODIFY THE "set_size" and "advance" PARAMETERS IF THEY ARE NOT
		// SPECIFIED! MAKE SURE YOU SET THESE TO DEFAULT VALUES BEFORE NOW!
		//
		int parse( const char *str )
		{
			std::vector< std::string > tokens;
			String::Tokenize( str, tokens, ":" );
			path = tokens[0];
			if( tokens.size() > 1 )
			{
				if( String::ToInteger( tokens[1], set_size ) != String::ReturnValue::OK ) return -1;
				if( tokens.size() > 2 )
				{
					advance = ( tokens[2] == "noadvance" ) ? (0) : (1);
				}
			}
			return 1;
		}
	};
	
	Input target, structures;
	std::vector< Superposition > superpositions;
	std::string output_path;
	bool print_rmsd;
	
	//
	// Check some basic parameters are reasonable
	//
	int check() const
	{
		int result = 1;
		
		if( target.path == "" )
		{
			fprintf( stderr, "No target path specified!\n" );
			result = -1;
		}
		if( target.set_size < 1 )
		{
			fprintf( stderr, "Bad target set size (%d)!\n", target.set_size );
			result = -1;
		}

		if( structures.path == "" )
		{
			fprintf( stderr, "No structures path specified!\n" );
			result = -1;
		}
		if( structures.set_size < 1 )
		{
			fprintf( stderr, "Bad structures set size (%d)!\n", structures.set_size );
			result = -1;
		}
		
		if( output_path == "" )
		{
			fprintf( stderr, "No output path specified!\n" );
			result = -1;
		}
		
		if( superpositions.size() < 1 )
		{
			fprintf( stderr, "No superpositions defined!\n" );
			result = -1;
		}
		
		//
		// Check the superposition molecule indices are valid for the set sizes
		//
		for( size_t i=0; i<superpositions.size(); i++ )
		{
			int index;
			const auto& s = superpositions[i];
			
			// "superpose" and "apply_to" describe UNIT BASED indices into "structures", so indices must be smaller than or equal
			// to "structures" set size!
			for( size_t j=0; j<s.superpose.size(); j++ )
			{
				index = s.superpose[j];
				if( index < 1 || index > structures.set_size )
				{
					fprintf( stderr, "Superposition %d has a bad superpose index! Must be 1 <= index <= %d\n", (int)i, structures.set_size );
					result = -1;
				}
			}
			for( size_t j=0; j<s.apply_to.size(); j++ )
			{
				index = s.apply_to[j];
				if( index < 1 || index > structures.set_size )
				{
					fprintf( stderr, "Superposition %d has a bad apply_to index! Must be 1 <= index <= %d\n", (int)i, structures.set_size );
					result = -1;
				}
			}

			// "onto" describes UNIT BASED indices into "target", so indices must be smaller than or equal
			// to "target" set size!
			for( size_t j=0; j<s.onto.size(); j++ )
			{
				int index = s.onto[j];
				if( index < 1 || index > target.set_size )
				{
					fprintf( stderr, "Superposition %d has a bad onto index! Must be 1 <= index <= %d\n", (int)i, target.set_size );
					result = -1;
				}
			}
		}
		
		return result;
	}
		
	
	//
	// Constructor!
	//
	Parameters( int argc, char **argv )
	{
		std::vector< std::string > tokens;
	
		//
		// By default the target advances and the structures do not.
		//
		target.path = "";
		target.set_size = 1;
		target.advance = 1;
		
		structures.path = "";
		structures.set_size = 1;
		structures.advance = 0;
		
		superpositions.clear();

		output_path = "";
		
		print_rmsd = false;
	
		for( int i=1; i<argc; i++ )
		{
			Superposition superpos;
			std::string key;
			std::vector< std::string > values;

			//
			// Specific parameter, starting with '-'?
			//
			if( argv[i][0] == '-' )
			{
				if( String::Tokenize( argv[i], tokens, "=" ) < 2 )
				{
					fprintf( stderr, "Bad parameter '%s'!\n", argv[i] );
					exit( -1 );
				}
				
				std::string &key = tokens[0];
				std::string &val = tokens[1];
				
				if( key == "-target" )
				{
					if( target.parse( val.c_str() ) == -1 )
					{
						fprintf( stderr, "Bad target parameter '%s'!\n", argv[i] );
						exit( -1 );
					}
				}
				else if( key == "-structures" )
				{
					if( structures.parse( val.c_str() ) == -1 )
					{
						fprintf( stderr, "Bad structures parameter '%s'!\n", argv[i] );
						exit( -1 );
					}
				}
				else if( key == "-output" )
				{
					output_path = val;
				}
				else if( key == "-print_rmsd" )
				{
					print_rmsd = (val == "yes") ? (true) : (false);
				}
			}
			else
			{
				if( superpos.parse( argv[i] ) == -1 )
				{
					exit( -1 );
				}
				superpositions.push_back( superpos );
			}
		}
		
		if( check() == -1 )
		{
			exit( -1 );
		}
	}
	
	void print() const
	{
		printf( "Parameters:\n" );
		
		printf( "\t target : %s set_size=%d advance=%d\n", target.path.c_str(), target.set_size, target.advance );
		printf( "\t structures : %s set_size=%d advance=%d\n", structures.path.c_str(), structures.set_size, structures.advance );
		printf( "\t output : %s\n", output_path.c_str() );
		
		for( size_t i=0; i<superpositions.size(); i++ )
		{
			superpositions[i].print();
		}
	}
	
};


void print_usage( const char *prog )
{
	printf( "\n" );

	printf( "Usage:\n" );
	printf( "\n" );

	printf( "\t %s -target=a.pdb[:set_size][:noadvance] -structures=b.pdb[:set_size][:noadvance] -output=c.pdb [ -print_rmsd=yes|no ] [superpose=i,j,...:onto=k,l,...:apply_to=m,n,...:attr1=val1,val2,...:attr2=val1,val2,...] ...\n", prog );
	printf( "\n" );
	
	printf( "Where:\n" );
	printf( "\n" );

	printf( "\t -target: input PDB with molecules to superpose ONTO\n" );
	printf( "\t -structures: input PDB with molecule to superpose ONTO those in \"target\"\n" );
	printf( "\t -output: results file\n" );
	printf( "\t -print_rmsd: print the before and after RMSD of superpositions\n" );
	printf( "\n" );

	printf( "And one or more superposition entries, with parameters (separated by colon, ':'):\n" );
	printf( "\n" );
	printf( "\t superpose: molecule indices into current \"structures\" set, UNIT BASED\n" );
	printf( "\t onto: molecule indices into current \"target\" set, UNIT BASED\n" );
	printf( "\t apply_to: molecule indices into current \"structures\" set, UNIT BASED\n" );
	printf( "\t key=values: restrictions on the atoms used in the superposition. For PDB files, these keys are the\n" );
	printf( "\t             PDB ATOM attributes such as \"name\", \"resName\", \"resSeq\" etc (see PDB file format).\n" );
	printf( "\n" );

	printf( "Notes:\n" );
	printf( "\n" );
	printf( "\t The optional 'set_size' and 'noadvance' parameters control how the contents of input files are considered.\n" );
	printf( "\t Groups of 'set_size' molecules are read in at a time, and it is on these sets that the superpositions act.\n" );
	printf( "\t The default 'set_size' is 1, and by default 'noadvance' is set for 'target' and 'structures' This results in\n" );
	printf( "\t consective sets from 'structures' being superposed onto consective sets from 'target'. Most of the time, you\n" );
	printf( "\t may wish to superpose only the first set of 'structures' onto all sets in 'target', which is specified as so:\n" );
	printf( "\n" );
	printf( "\t -target=a.pdb:1 -structures=b.pdb:1:noadvance\n" );
	printf( "\n" );
	printf( "\t Where ranges of values are provided in the values for superposition entries, using the dash character '-',\n" );
	printf( "\t the range is expanded into an INCLUSIVE set of indices. Checks are performed to ensure that the ranges\n" );
	printf( "\t are integers and that the start index is <= end index.\n" );
	printf( "\n" );
	printf( "\t Where values are provided for the keys in the superposition, the set of atoms to use for the superposition\n" );
	printf( "\t must have attributes in the list provided; where two atoms have attributes in the list specified, but the\n" );
	printf( "\t attributes are not identical for both atoms, those atoms are ignored.\n" );
	printf( "\n" );
	printf( "\t An empty value list for a particular key denotes that the named attribute must be the same for two atoms to be\n" );
	printf( "\t included in those atoms used for the superposition calculation. any attributes not specified in the superposition\n" );
	printf( "\t definition are simply ignored.\n" );
	printf( "\n" );

	printf( "Examples:\n" );
	printf( "\n" );
	printf( "\t 1. %s -target=tgt.pdb -structures=str.pdb:4:noadvance -output=out.pdb superpose=1,2:onto=1:apply_to=3,4:name=:resName=:resSeq=12-34\n", prog );
	printf( "\n" );
	printf( "\t Calculate superposition of molecules 1 and 2 (combined into single molecule) of str.pdb onto each sequential\n" );
	printf( "\t molecule of tgt.pdb where the atom names and residue names match, and the residue sequences are in the range of\n" );
	printf( "\t 12 to 34 (inclusive). This superposition transform is applied to molecules 3 and 4 of str.pdb and the results\n" );
	printf( "\t saved into out.pdb. We specify a set size of 4 for str.pdb, and prevent advancing beyond the first set.\n" );
	printf( "\n" );
	printf( "\t 2. %s -target=tgt.pdb:4 -structures=str.pdb:5:noadvance -output=out.pdb superpose=1,2:onto=1,2:apply_to=3,4,5:name=CA,CB:resName=:resSeq=\n", prog );
	printf( "\n" );
	printf( "\t Calculate superposition of molecules 1 and 2 (combined into single molecule) of the first set of 5 consective molecules in str.pdb\n" );
	printf( "\t onto the combined molecule formed from entries 1 and 2 in each set of tgt.pdb (with tgt.pdb processed using sets of 4 consecutive molecules).\n" );
	printf( "\t The atom names are one of either \"CA\" or \"CB\" but must be the same for the paired atoms used in the superposition. The residue\n" );
	printf( "\t names and sequence numbers must match, and the results are saved into out.pdb.\n" );
	
	exit( -1 );
}

int main( int argc, char **argv )
{
	std::string target_step;
	Molecule::MolecularSystem target_molset, structures_molset, output_molset;
	
	FILE *target_file, *structures_file, *output_file;
	
	Parameters *params;
	
	if( argc < 5 ) print_usage( argv[0] );
	
	params = new Parameters( argc, argv );
	
	printf( "\n" );
	params->print();
	printf( "\n" );

	if( (target_file = fopen( params->target.path.c_str(), "r" )) == NULL )
	{
		fprintf( stderr, "Unable to open target file '%s'!\n", params->target.path.c_str() );
		exit( -1 );
	}
	if( (structures_file = fopen( params->structures.path.c_str(), "r" )) == NULL )
	{
		fclose( target_file );

		fprintf( stderr, "Unable to open structures file '%s'!\n", params->structures.path.c_str() );
		exit( -1 );
	}
	if( (output_file = fopen( params->output_path.c_str(), "w" )) == NULL )
	{
		fclose( target_file );
		fclose( structures_file );

		fprintf( stderr, "Unable to open output file '%s'!\n", params->output_path.c_str() );
		exit( -1 );
	}

	bool have_targets = false;
	bool have_structures = false;

	//
	// Loop over target molecules
	//
	int target_set = 0;
	int structures_set = 0;
	int total_sets = 0;
	while( true )
	{
		int target_N, N;
		std::vector<double> xyz, target_xyz;
		double M[3][3], T[3], rmsd1, rmsd2;
		std::vector<Molecule::Molecule> superpose_mols, onto_mols, applyto_mols, filtered_superpose_mols, filtered_onto_mols;
		
		std::vector< int > N_vec;
		std::vector< double > rmsd1_vec, rmsd2_vec;
		
		//
		// Load next target molecule set, checking for problems. If we're not advancing through the target sets, we only do this once!
		//
		if( params->target.advance > 0 || have_targets == false )
		{
			N = Molecule::PDB::Load( target_file, target_molset, params->target.set_size );
			if( N == 0 ) break; // assume EOF, stop here.
			if( N != params->target.set_size )
			{
				fprintf( stderr, "Unable to read a set of %d molecules from target file! Could only read %d!\n", params->target.set_size, N );
				exit( -1 );
			}
			
			have_targets = true;
			target_set++;
		}
	
		//
		// Loop over structure molecules
		//
		while( true )
		{
			total_sets++;
			
			//
			// Load next structure molecule set, checking for problems. If we're not advancing through the structure sets, we only do this once!
			//
			if( params->structures.advance > 0 || have_structures == false )
			{
				N = Molecule::PDB::Load( structures_file, structures_molset, params->structures.set_size );
				if( N == 0 ) break; // assume EOF, stop here.
				if( N != params->structures.set_size )
				{
					fprintf( stderr, "Unable to read a set of %d molecules from structures file! Could only read %d!\n", params->structures.set_size, N );
					exit( -1 );
				}

				have_structures = true;
				structures_set++;
			}
			
			N_vec.clear();
			rmsd1_vec.clear();
			rmsd2_vec.clear();

			//
			// Loop over defined superpositions, using current target and structure data
			//
			for( const auto& s : params->superpositions )
			{
				// build sets we need.
				superpose_mols.clear();
				onto_mols.clear();
				applyto_mols.clear();
				filtered_superpose_mols.clear();
				filtered_onto_mols.clear();

				// get the molecules we're superposing from structures_molset
				for( size_t i=0; i<s.superpose.size(); i++ )
				{
					int index = s.superpose[i]-1;
					superpose_mols.push_back( structures_molset.molecules[index] );
				}

				// get the molecules we're superposing ONTO from target_molset, using current offset
				for( size_t i=0; i<s.onto.size(); i++ )
				{
					int index = s.onto[i]-1;
					onto_mols.push_back( target_molset.molecules[index] );
				}

				// get the molecules we're applying the superposition to, fron structures_molset
				for( size_t i=0; i<s.apply_to.size(); i++ )
				{
					int index = s.apply_to[i]-1;
					applyto_mols.push_back( structures_molset.molecules[index] );
				}

				// filter the subject and target atom sets so we only have the appropriate common atoms - error where no common atoms!
				if( s.mf.Filter( onto_mols, superpose_mols, filtered_onto_mols, filtered_superpose_mols ) < 1 )
				{
					printf( "%s(): unsuccessful filtering for target set %d, structures set %d\n", __func__, target_set, structures_set );
					break;
					//exit( -1 );
				}

				// get raw coord arrays for superposition target and subject atoms
				target_N = Molecule::Coords::Get( filtered_onto_mols, target_xyz );
				N = Molecule::Coords::Get( filtered_superpose_mols, xyz );

				// calculate superposition of filtered point sets, and RMSD of filtered superposition
				Superposer::Calculate( target_xyz, xyz, M, T );
				rmsd1 = Superposer::RMSD( xyz, target_xyz );
				Superposer::Apply( xyz, xyz, M, T );
				rmsd2 = Superposer::RMSD( xyz, target_xyz );
				
				printf( "RMSD: %g => %g\n", rmsd1, rmsd2 );

				// get raw coord array for atoms we're actually applying the superposition to, and transform those coords
				N = Molecule::Coords::Get( applyto_mols, xyz );
				Superposer::Apply( xyz, xyz, M, T );

				// copy newly transformed coords from raw arrays back into the molecule structures
				Molecule::Coords::Set( applyto_mols, xyz );

				// build our output data set as we go
				for( const auto& mol : applyto_mols ) output_molset.molecules.push_back( mol );

				N_vec.push_back( target_N );
				rmsd1_vec.push_back( rmsd1 );
				rmsd2_vec.push_back( rmsd2 );
			}

			//
			// If we're printing the RMSD, or we've processed 10 sets, print some RMSD info.
			//
			if( params->print_rmsd == true || total_sets % 10 == 0  )
			{
				printf( "Set %d : target set at %d (%d) structures at %d (%d)\n",
					total_sets,
					target_set, params->target.set_size,
					structures_set, params->structures.set_size );
				for( size_t i=0; i<N_vec.size(); i++ )
				{
					printf( "\t %d atoms, RMSD %f -> %f\n", N_vec[i], rmsd1_vec[i], rmsd2_vec[i] );
				}
			}

			//
			// Progressive generation of output file.
			//
			if( total_sets % 10 == 0 )
			{
				Molecule::PDB::Print( output_file, output_molset );
				output_molset.molecules.clear();
			}
			
			//
			// If we are only using a single structure set, stop this loop here.
			//
			if( params->structures.advance == 0 ) break;
		}


		//
		// If we are only using a single target set, stop this loop here.
		//
		if( params->target.advance == 0 ) break;
	}
	
	// save remaining data to output pdb
	if( output_molset.molecules.size() > 0 )
	{
		Molecule::PDB::Print( output_file, output_molset );
		output_molset.molecules.clear();
	}

	fclose( target_file );
	fclose( structures_file );
	fclose( output_file );
}
