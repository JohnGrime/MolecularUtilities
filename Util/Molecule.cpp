/*
	Author: John Grime, The University of Chicago.
*/

#include "Molecule.h"

//
// PDB entry information, from PDB spec.
//
static const char * pdb_keys[] = { "recname", "serial", "name", "altLoc", "resName", "chainID", "resSeq", "iCode", "x", "y", "z", "occupancy", "tempFactor", "element", "charge" };
static const int pdb_col_start[] = {  1,  7, 13, 17, 18, 22, 23, 27, 31, 39, 47, 55, 61, 77, 79 };
static const int pdb_col_end[]   = {  6, 11, 16, 17, 20, 22, 26, 27, 38, 46, 54, 60, 66, 78, 80 };
//static const char *chainIDs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz";
//static const char *pdb_format_string = "%-6.6s%5.5s %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6s%6s%6s%4s%2s%2s";
static const char *pdb_format_string = "%-6.6s%5.5s %-4s%1s%3s %1s%4s%1s   %8.2f%8.2f%8.2f%6s%6s%6s%4s%2s%2s"; // for additional leading digit!

namespace Util
{
namespace Molecule
{

void AttributeFilter::AddFilter( const char *key, const char *value )
{
	FilterValuesMap empty;
	
	assert( key != nullptr );

	// if we define no value, the map of acceptable values is empty
	// if we pass a value, ADD it to the acceptable values.
	if( value == nullptr )
	{
		filters[key] = empty;
	}
	else
	{
		filters[key][value] = 1;
	}
}
void AttributeFilter::AddFilter( const char *filter_str, const char *keyval_sep, const char *val_sep, const char *range_sep )
{
	//
	// Expect : key=v1,v2-v3,...
	// Where val_sep == "," and range_sep == "-" in the above.
	//
	std::vector< std::string > tokens, val_tokens, range_tokens;
	
	assert( filter_str != nullptr );
	assert( val_sep != nullptr );
	assert( range_sep != nullptr );
	
	String::Tokenize( filter_str, tokens, keyval_sep );
	
	//
	// Empty strings, or just a key?
	//
	if( tokens.size() == 0 ) return;
	if( tokens.size() == 1 )
	{
		AddFilter( tokens[0].c_str() );
		return;
	}
	
	const std::string &key = tokens[0];
	const std::string &values = tokens[1];

	String::Tokenize( values.c_str(), val_tokens, val_sep );
	for( const auto& value : val_tokens )
	{
		//
		// Could this value token denote an integer range?
		//
		String::Tokenize( value.c_str(), range_tokens, range_sep );
		if( range_tokens.size() > 1 )
		{
			int start, stop;
			char buffer[128];
			
			// check if we can convert the subtokens into integers. If not, assume string.
			const std::string &r1 = range_tokens[0];
			const std::string &r2 = range_tokens[1];
			
			if( (String::ToInteger(r1,start) != String::ReturnValue::OK) ||
				(String::ToInteger(r2,stop)  != String::ReturnValue::OK) )
			{
				// Assume not a range, add raw value token
				AddFilter( key.c_str(), value.c_str() );
			}
			else
			{
				// Add integers in range as strings to filter
				if( stop < start )
				{
					int temp = start;
					start = stop;
					stop = temp;
				}
				for( int i=start; i<=stop; i++ )
				{
					sprintf( buffer, "%d", i );
					AddFilter( key.c_str(), buffer );
				}
			}
		}
		else
		{
			// Assume not a range, add raw value token
			AddFilter( key.c_str(), value.c_str() );
		}
	}
	
}

int AttributeFilter::Filter( const Atom &in ) const
{
	FilterMap::const_iterator f_it;
	
	for( const auto& f_it : filters )
	{
		const std::string &key = f_it.first;
		const FilterValuesMap &values = f_it.second;
		
		//
		// No attribute matching this filter in the atom? automatic fail.
		//
		const auto& a_it1 = in.attr.find(key);
		if( a_it1 == in.attr.end() ) return 0;
		
		const auto& v1 = a_it1->second;

		//
		// If we have acceptable values defined for this attribute, test against them.
		//
		if( values.size() > 0 )
		{
			if( values.find(v1) == values.end() ) return 0;
		}
	}
	return 1;
}
int AttributeFilter::Filter( const Atom &in1, const Atom &in2 ) const
{
	FilterMap::const_iterator f_it;
	
	for( const auto& f_it : filters )
	{
		const std::string &key = f_it.first;
		const FilterValuesMap &values = f_it.second;
		
		//
		// No attribute matching this filter in the atom? automatic fail.
		//
		const auto& a_it1 = in1.attr.find(key);
		const auto& a_it2 = in2.attr.find(key);
		if( a_it1 == in1.attr.end() ) return 0;
		if( a_it2 == in2.attr.end() ) return 0;
		
		const std::string &v1 = a_it1->second;
		const std::string &v2 = a_it2->second;
		
		//
		// If we have acceptable values defined for this attribute, test against them.
		//
		if( values.size() > 0 )
		{
			if( values.find(v1) == values.end() ) return 0;
			if( values.find(v2) == values.end() ) return 0;
		}
		
		//
		// Finally, check the values for this attribute are the same in both atoms.
		//
		if( v1 != v2 ) return 0;
	}
	return 1;
}
int AttributeFilter::Filter( const Molecule &in, Molecule &out ) const
{
	Molecule temp; // allows in-place filtering if &in == &out
	
	for( const auto& a : in )
	{
		if( Filter(a) == 1 ) temp.push_back( a );
	}
	
	out = temp;
	
	return (int)out.size();
}
int AttributeFilter::Filter( const Molecules &in, Molecules &out ) const
{
	size_t N;
	Molecules temp; // allows in-place filtering if e.g. &in == &out
	
	N = 0;
	for( const auto in_m : in )
	{
		Molecule out_m;
		
		// Filter specific molecule
		Filter( in_m, out_m );
		
		N += out_m.size();
		
		// Build molecule vector
		temp.push_back( out_m );
	}
	
	out = temp;
	
	return (int)N;
}
int AttributeFilter::Filter( const Molecule &in1, const Molecule &in2, Molecule &out1, Molecule &out2 ) const
{
	std::string key;
	Molecule temp1, temp2, new1, new2; // allows in-place filtering if eg &in1 == &out1
	size_t start_j;
	
	//
	// First, reduce both molecules on the basis of the literal value filters.
	// this reduces the amount of work we need to do!
	//
	Filter( in1, temp1 );
	Filter( in2, temp2 );

	//
	// note: we could loop over smallest set, and simply use pointers to output etc
	// but this complicates the code even though it could be significantly faster.
	//
	start_j = 0;
	for( size_t i=0, max_i=temp1.size(); i<max_i; i++ )
	{
		const Atom &ai = temp1[i];

		int j = -1;
		for( size_t jj=start_j, max_jj=temp2.size(); jj<max_jj; jj++ )
		{
			const Atom &aj = temp2[jj];
			if( Filter( ai, aj ) == 1 )
			{
				j = (int)jj;
				break;
			}
		}
		if( j != -1 )
		{
			new1.push_back( temp1[i] );
			new2.push_back( temp2[j] );
			start_j = j+1;
		}
	}
	
	out1 = new1;
	out2 = new2;
	
	assert( out1.size() == out2.size() );
	
	return (int)out1.size();
}
int AttributeFilter::Filter( const Molecules &in1, const Molecules &in2, Molecules &out1, Molecules &out2 ) const
{
	size_t t1, t2;
	Molecules temp1, temp2; // allows in-place filtering if e.g. &in1 == &out1
	
	assert( in1.size() == in2.size() );
	
	t1 = t2 = 0;
	for( size_t i=0, max_i=in1.size(); i<max_i; i++ )
	{
		Molecule m1, m2;
		
		// Filter specific molecules
		Filter( in1[i], in2[i], m1, m2 );
		
		t1 += m1.size();
		t2 += m2.size();
		
		// Build molecule vectors
		temp1.push_back( m1 );
		temp2.push_back( m2 );
	}
	
	assert( t1 == t2 );
	assert( temp1.size() == temp2.size() );
	
	out1 = temp1;
	out2 = temp2;
	
	return t1;
}
void AttributeFilter::Print( FILE *f ) const
{
	if( f == nullptr ) f = stdout;
	
	fprintf( f, "MoleculeFilter::%s(): %d filters.\n", __func__, (int)filters.size() );
	for( const auto& f_it : filters )
	{
		const std::string &key = f_it.first;
		const FilterValuesMap &values = f_it.second;
		
		fprintf( f, "\t '%s' : ", key.c_str() );
		if( values.size() < 1 )
		{
			fprintf( f, "(unspecified)\n" );
			continue;
		}
		for( const auto& fv_it : values )
		{
			fprintf( f, "'%s' ", fv_it.first.c_str() );
		}
		fprintf( f, "\n" );
	}
}


void MolecularSystem::Print( FILE *f, int print_n_mols, int print_n_atoms ) const
{
	AttributeMap::const_iterator it;

	if( f == nullptr ) f = stdout;
	
	{
		size_t n_mols = molecules.size(), n_atoms = 0;
		for( const auto& m : molecules )
		{
			n_atoms += m.size();
		}
	
		fprintf( f, "MolSet::%s(): %d molecules, %d atoms\n", __func__, (int)n_mols, (int)n_atoms );
	}

	if( metadata.size() > 0 )
	{
		fprintf( f, "Metadata:\n" );
		for( it = metadata.begin(); it != metadata.end(); it++ )
		{
			fprintf( f, "\t '%s'='%s'\n", it->first.c_str(), it->second.c_str() );
		}
	}
	else
	{
		fprintf( f, "No metadata.\n" );
	}
	
	fprintf( f, "\n" );
	
	int max_mols = (int)molecules.size();
	int n_mols = (print_n_mols > max_mols) ? (max_mols) : (print_n_mols);

	for( int mol_i=0; mol_i<n_mols; mol_i++ )
	{
		const Molecule &mol = molecules[mol_i];

		int max_atoms = (int)mol.size();
		int n_atoms = (print_n_atoms > max_atoms) ? (max_atoms) : (print_n_atoms);
		
		fprintf( f, "\t %d (%d)\n", mol_i+1, max_atoms );
		
		for( int atom_i=0; atom_i<n_atoms; atom_i++ )
		{
			const Atom &a = mol[atom_i];
			
			fprintf( f, "\t\t %+.2e %+.2e %+.2e ", a.x, a.y, a.z );

			for( it = a.attr.begin(); it != a.attr.end(); it++ )
			{
				fprintf( f, "'%s'='%s' ", it->first.c_str(), it->second.c_str() );
			}
			fprintf( f, "\n" );
		}
	}
	
	fprintf( f, "\n" );
}


int PDB::Load( FILE *f, MolecularSystem &ms, int n_mols )
{
	bool include_HETATM = true;
	
	char line_buffer[1024], buffer[1024];
	const char *delimiters = " \t\n";
	std::vector< std::string > tokens;
	Molecule mol;
	Atom a;
	int line_length, n_mols_read;

	assert( f != nullptr );

	//
	// Generated directly from the PDB data arrays; this is convenient, as it keeps the code short and we can
	// automatically detect the number of entries in the arrays using sizeof().
	//
	std::vector< std::string > keys( pdb_keys, pdb_keys + sizeof(pdb_keys) / sizeof(pdb_keys[0]) );
	std::vector< int > starts( pdb_col_start, pdb_col_start + sizeof(pdb_col_start) / sizeof(pdb_col_start[0]) );
	std::vector< int > ends( pdb_col_end, pdb_col_end + sizeof(pdb_col_end) / sizeof(pdb_col_end[0]) );

	ms.molecules.clear();
	mol.clear();

	n_mols_read = 0;

	while( fgets( line_buffer, 1023, f ) != nullptr )
	{
		String::Tokenize( line_buffer, tokens, delimiters );

		if( tokens.size() < 1 ) continue;

		if( (tokens[0]=="ATOM") || (include_HETATM && (tokens[0]=="HETATM")) )
		{
			const char *tok;
			
			line_length = strlen( line_buffer );

			a.attr.clear();
			for( size_t i=0; i<keys.size(); i++ )
			{
				const std::string &key = keys[i];
				int start = starts[i];
				int stop = ends[i];
				
				a.attr[ key ] = "";

				if( start > line_length || stop > line_length ) break;

				//
				// Remember! PDB column start/stop values are UNIT BASED AND INCLUSIVE!
				//
				int k=0;
				for( int j=start-1; j<=stop-1; j++ )
				{
					buffer[k++] = line_buffer[j];
				}
				buffer[k] = '\0';
				String::Strip( buffer, delimiters );
				a.attr[ keys[i] ] = buffer;
			}
			
			tok = a.attr["x"].c_str();
			if( String::ToReal(tok,a.x) != String::ReturnValue::OK )
			{
				fprintf( stderr, "PDBLoader::Load() : unable to convert x token '%s' into a real number!\n", tok );
				return -1;
			}

			tok = a.attr["y"].c_str();
			if( String::ToReal(tok,a.y) != String::ReturnValue::OK )
			{
				fprintf( stderr, "PDBLoader::Load() : unable to convert y token '%s' into a real number!\n", tok );
				return -1;
			}

			tok = a.attr["z"].c_str();
			if( String::ToReal(tok,a.z) != String::ReturnValue::OK )
			{
				fprintf( stderr, "PDBLoader::Load() : unable to convert z token '%s' into a real number!\n", tok );
				return -1;
			}
			
			mol.push_back( a );
		}
		else if( tokens[0] == "TER" || tokens[0] == "ENDMDL" || tokens[0] == "END" )
		{
			if( mol.size() > 0 )
			{
				ms.molecules.push_back( mol );
				n_mols_read++;
			}
			mol.clear();
		}

		// stop if hit a specified limit to read.
		if( n_mols > 0 && n_mols_read >= n_mols ) break;
	}
	if( mol.size() > 0 ) ms.molecules.push_back( mol );

	return (int)ms.molecules.size();
}
int PDB::Print( FILE *f, const Atom &a_ )
{
	const int max_linebuffer = 128; // PDB ATOM/HETATM line is 80 columns.
	char linebuffer[max_linebuffer];

	std::string key;

	assert( f != nullptr );
	
	memset( linebuffer, ' ', max_linebuffer );
	
	Atom a = a_; // COPY, as we may be adding empty strings on the fly when we access the attr map below!

	// print PDB line. If atributes missing in the map, they'll be initialised to empty strings.
	sprintf( linebuffer, pdb_format_string,
		a.attr["recname"].c_str(),
		a.attr["serial"].c_str(),
		a.attr["name"].c_str(),
		a.attr["altLoc"].c_str(),
		a.attr["resName"].c_str(),
		a.attr["chainID"].c_str(),
		a.attr["resSeq"].c_str(),
		a.attr["iCode"].c_str(),
		a.x,
		a.y,
		a.z,
		a.attr["occupancy"].c_str(),
		a.attr["tempFactor"].c_str(),
		"",
		"",
		a.attr["element"].c_str(),
		a.attr["charge"].c_str() );

	linebuffer[80] = '\0';
	fprintf( f, "%s", linebuffer );

	return 1;
}
int PDB::Print( FILE *f, const Molecule &mol )
{
	assert( f != nullptr );

	for( const auto& a : mol )
	{
		Print( f, a );
		fprintf( f, "\n" );
	}
	return 1;
}
int PDB::Print( FILE *f, const Molecules &mols )
{
	assert( f != nullptr );

	for( const auto& m : mols )
	{
		Print( f, m );
		fprintf( f, "TER\n" );
	}
	return 1;
}
int PDB::Print( FILE *f, const MolecularSystem &ms )
{
	assert( f != nullptr );
	
	Print( f, ms.molecules );
	return 1;
}


/*
	Internal fast tokenize stuff, as LAMMPS trajectories can be slow to load via the default methods!
*/

//
//	Chop up "str", results in array "pointers". If too many substrings, this routine will just do as
//	many as possible.
//
static int lmp_isdelim( char test, const char *delim, int n_delim )
{
	for( int i=0; i<n_delim; i++ )
	{
		if( test == delim[i] ) return 1;
	}
	return 0;
}
static int lmp_tokenize( int str_len, const char *str, int n_delim, const char *delim, char **pointers, int maxptrs )
{
	int ntoks = 0, i, j;
	
	j = 0; // index into current token.
	for( i=0; i<str_len; i++ )
	{
		if( lmp_isdelim( str[i], delim, n_delim ) == 1 )
		{
			// If we've been building a token and hit whitespace, increase ntokens and reset token character index
			if( j > 0 ) ntoks++;
			j = 0;
			if( ntoks >= maxptrs ) break;
		}
		else
		{
			pointers[ntoks][j] = str[i];
			pointers[ntoks][j+1] = '\0';
			j++;
		}
	}
	
	// We may have a trailing token at the end of the string
	if( j > 0 ) ntoks++;
	
	return ntoks;
}

int LAMMPS::Load( FILE *f, MolecularSystem &ms )
{
	using lmp_atom = struct lmp_atom_
	{
		int index, type, mol;
		float x, y, z;
		
		// Sorts on mol, then particle index. Hopefully this is good enough.
		bool operator < (const lmp_atom_ &rhs ) const
		{
			if( mol != rhs.mol ) return (mol < rhs.mol);
//			return (index < rhs.index );
			return (type < rhs.type );
		}
		
		int set(
			const char *index_,
			const char *mol_,
			const char *type_,
			const char *x_,
			const char *y_,
			const char *z_
			)
		{
			if( String::ToInteger( index_, index ) != String::ReturnValue::OK ) return -1;
			if( String::ToInteger( mol_, mol )     != String::ReturnValue::OK ) return -1;
			if( String::ToInteger( type_, type )   != String::ReturnValue::OK ) return -1;

			if( String::ToReal( x_, x ) != String::ReturnValue::OK ) return -1;
			if( String::ToReal( y_, y ) != String::ReturnValue::OK ) return -1;
			if( String::ToReal( z_, z ) != String::ReturnValue::OK ) return -1;
			
			return 1;
		}
	};

	int ntoks, buffer_length, n_delim;
	int id_token, mol_token, type_token, x_token, y_token, z_token;
	char buffer[1024];
	fpos_t file_position;
	const char *delim = " \t\n";

	char tokens[100][128], *token_ptrs[100];

	// frame metadata - will be set as strings later on.
	int N = 0, timestep, pbc[3];
	double bounds[6];
	
	std::vector<lmp_atom> lmp_atoms;

	assert( f != nullptr );

	ms.molecules.clear();
	ms.metadata.clear();

	n_delim = strlen( delim );
	for( int i=0; i<100; i++ ) token_ptrs[i] = tokens[i];
	
	timestep = -1;
	
	bounds[0] = -1.0;
	bounds[1] = -1.0;
	bounds[2] = -1.0;
	bounds[3] = -1.0;
	bounds[4] = -1.0;
	bounds[5] = -1.0;
	
	pbc[0] = 0;
	pbc[1] = 0;
	pbc[2] = 0;
	
	lmp_atoms.clear();

	while( fgets( buffer, 1024, f ) != nullptr )
	{
		buffer_length = strlen( buffer );
		if( buffer_length > 0 && buffer[buffer_length-1] == '\n' ) buffer[ buffer_length-1 ] = '\0';
		
		if( (ntoks = lmp_tokenize( buffer_length, buffer, n_delim, delim, token_ptrs, 100 )) < 1 || strcmp( tokens[0], "ITEM:" ) != 0 ) continue;

		if( strcmp( tokens[1], "TIMESTEP" ) == 0 )
		{
			if( fgets( buffer, 1024, f ) == nullptr ) continue;
			buffer_length = strlen( buffer );
			if( (ntoks = lmp_tokenize( buffer_length, buffer, n_delim, delim, token_ptrs, 100 )) < 1 ) continue;

			if( String::ToInteger( tokens[0], timestep ) != String::ReturnValue::OK )
			{
				printf( "MolecularUtil::LAMMPS::%s(): unable to convert token '%s' into timestep!\n", __func__, tokens[0]  );
				printf( "'%s'\n", buffer );
				return -1;
			}
		}
		else if( strcmp( tokens[1], "BOX" ) == 0 )
		{
			// expects ITEM: BOX BOUNDS xx yy zz
			if( ntoks != 6 )
			{
				printf( "MolecularUtil::LAMMPS::%s(): bad box line\n", __func__ );
				printf( "'%s'\n", buffer );
				return -1;
			}
			
			if( tokens[3][0] == 'p' ) pbc[0] = 1;
			if( tokens[4][0] == 'p' ) pbc[1] = 1;
			if( tokens[5][0] == 'p' ) pbc[2] = 1;
			
			//
			// X bounds as min, max
			//
			if( fgets( buffer, 1024, f ) == nullptr ) continue;
			buffer_length = strlen( buffer );
			if( (ntoks = lmp_tokenize( buffer_length, buffer, n_delim, delim, token_ptrs, 100 )) < 1 ) continue;

			if( String::ToReal( tokens[0], bounds[0] ) != String::ReturnValue::OK )
			{
				printf( "MolecularUtil::LAMMPS::%s(): bad bounds value '%s'\n", __func__, tokens[0] );
				printf( "'%s'\n", buffer );
				return -1;
			}
			if( String::ToReal( tokens[1], bounds[1] ) != String::ReturnValue::OK )
			{
				printf( "MolecularUtil::LAMMPS::%s(): bad bounds value '%s'\n", __func__, tokens[1] );
				printf( "'%s'\n", buffer );
				return -1;
			}

			//
			// Y bounds as min, max
			//
			if( fgets( buffer, 1024, f ) == nullptr ) continue;
			buffer_length = strlen( buffer );
			if( (ntoks = lmp_tokenize( buffer_length, buffer, n_delim, delim, token_ptrs, 100 )) < 1 ) continue;

			if( String::ToReal( tokens[0], bounds[2] ) != String::ReturnValue::OK )
			{
				printf( "MolecularUtil::LAMMPS::%s(): bad bounds value '%s'\n", __func__, tokens[0] );
				printf( "'%s'\n", buffer );
				return -1;
			}
			if( String::ToReal( tokens[1], bounds[3] ) != String::ReturnValue::OK )
			{
				printf( "MolecularUtil::LAMMPS::%s(): bad bounds value '%s'\n", __func__, tokens[1] );
				printf( "'%s'\n", buffer );
				return -1;
			}

			//
			// Z bounds as min, max
			//
			if( fgets( buffer, 1024, f ) == nullptr ) continue;
			buffer_length = strlen( buffer );
			if( (ntoks = lmp_tokenize( buffer_length, buffer, n_delim, delim, token_ptrs, 100 )) < 1 ) continue;

			if( String::ToReal( tokens[0], bounds[4] ) != String::ReturnValue::OK )
			{
				printf( "MolecularUtil::LAMMPS::%s(): bad bounds value '%s'\n", __func__, tokens[0] );
				printf( "'%s'\n", buffer );
				return -1;
			}
			if( String::ToReal( tokens[1], bounds[5] ) != String::ReturnValue::OK )
			{
				printf( "MolecularUtil::LAMMPS::%s(): bad bounds value '%s'\n", __func__, tokens[1] );
				printf( "'%s'\n", buffer );
				return -1;
			}
		}
		else if( strcmp( tokens[1], "NUMBER" ) == 0 && strcmp( tokens[2], "OF" ) == 0 && strcmp( tokens[3], "ATOMS" ) == 0 )
		{
			if( fgets( buffer, 1024, f ) == nullptr ) continue;
			
			if( String::ToInteger( buffer, N ) != String::ReturnValue::OK )
			{
				printf( "MolecularUtil::LAMMPS::%s(): bad atom count '%s'\n", __func__, buffer );
				printf( "'%s'\n", buffer );
				return -1;
			}
			
			lmp_atoms.resize( N );
		}
		else if( strcmp( tokens[1], "ATOMS" ) == 0 )
		{
			id_token = -1;
			mol_token = -1;
			type_token = -1;
			x_token = -1;
			y_token = -1;
			z_token = -1;

			// determine where the appropriate tokens are!
			for( int i=0; i<ntoks; i++ )
			{
				if( strcmp( tokens[i], "id" ) == 0 ) id_token = i;
				else if( strcmp( tokens[i], "mol" ) == 0 ) mol_token = i;
				else if( strcmp( tokens[i], "type" ) == 0 ) type_token = i;
				else if( strcmp( tokens[i], "x" ) == 0 ) x_token = i;
				else if( strcmp( tokens[i], "xs" ) == 0 ) x_token = i;
				else if( strcmp( tokens[i], "y" ) == 0 ) y_token = i;
				else if( strcmp( tokens[i], "ys" ) == 0 ) y_token = i;
				else if( strcmp( tokens[i], "z" ) == 0 ) z_token = i;
				else if( strcmp( tokens[i], "zs" ) == 0 ) z_token = i;
			}
			
			if( id_token == -1 )
			{
				printf( "MolecularUtil::LAMMPS::%s(): unable to find id token in frame.\n", __func__ );
				return -1;
			}
			if( mol_token == -1 )
			{
				printf( "MolecularUtil::LAMMPS::%s(): unable to find mol token in frame.\n", __func__ );
				return -1;
			}
			if( type_token == -1 )
			{
				printf( "MolecularUtil::LAMMPS::%s(): unable to find type token in frame.\n", __func__ );
				return -1;
			}
			if( x_token == -1 )
			{
				printf( "MolecularUtil::LAMMPS::%s(): unable to find x token in frame.\n", __func__ );
				return -1;
			}
			if( y_token == -1 )
			{
				printf( "MolecularUtil::LAMMPS::%s(): unable to find y token in frame.\n", __func__ );
				return -1;
			}
			if( z_token == -1 )
			{
				printf( "MolecularUtil::LAMMPS::%s(): unable to find z token in frame.\n", __func__ );
				return -1;
			}
			
			
			
			int atom_upto = 0;
			while( 1 )
			{
				fgetpos( f, &file_position );
				if( fgets( buffer, 1024, f ) == nullptr ) break; // likely EOF
				
				buffer_length = strlen( buffer );
				if( buffer_length > 0 && buffer[buffer_length-1] == '\n' ) buffer[ buffer_length-1 ] = '\0';
				if( (ntoks = lmp_tokenize( buffer_length, buffer, n_delim, delim, token_ptrs, 100 )) >= 1 )
				{
					if( strcmp( tokens[0], "ITEM:" ) == 0 )
					{
						fsetpos( f, &file_position ); // move back to start of previous line, as this is no longer atom data.
						goto final_info; // I hate myself for this.
//						break;
					}
					
					lmp_atom &a = lmp_atoms[atom_upto];
					
					int result = a.set(
						tokens[id_token-2],
						tokens[mol_token-2],
						tokens[type_token-2],
						tokens[x_token-2],
						tokens[y_token-2],
						tokens[z_token-2] );
					
					if( result == -1 )
					{
						printf( "MolecularUtil::LAMMPS::%s(): bad atom\n", __func__ );
						printf( "'%s'\n", buffer );
						return -1;
					}
					
					atom_upto++;
					
					if( atom_upto > N )
					{
						printf( "LammpsTrajectoryFrame::%s(): too many atoms (%d > %d)\n", __func__, atom_upto, N );
						exit( -1 );
					}
				}
			}
		}
	}
	
	//
	// Ugly jump hack.
	//
	final_info:

	//
	// Hopefully, we have all the atom data present. Sort, and set up molecules etc.
	//
	std::sort( lmp_atoms.begin(), lmp_atoms.end() );
	
	Molecule mol;
	Atom a;

	mol.clear();
	int last_mol = -1;

	// defaults
	a.attr["recname"] = "ATOM";
	a.attr["resName"] = "LMP";

	for( size_t i=0, max_i=lmp_atoms.size(); i<max_i; i++ )
	{
		lmp_atom &lmp_a = lmp_atoms[i];
		
		if( i > 0 && lmp_a.mol != last_mol )
		{
			ms.molecules.push_back( mol );
			mol.clear();
		}
		
		a.x = lmp_a.x;
		a.y = lmp_a.y;
		a.z = lmp_a.z;
		
		sprintf( buffer, "%d", lmp_a.index );
		a.attr["serial"] = buffer;

		sprintf( buffer, "%d", lmp_a.mol );
		a.attr["resSeq"] = buffer;

		sprintf( buffer, "%d", lmp_a.type );
		a.attr["name"] = buffer;
		
		mol.push_back( a );
		
		last_mol = lmp_a.mol;
	}
	
	if( mol.size() > 0 )
	{
		ms.molecules.push_back( mol );
	}
	
	sprintf( buffer, "%d", timestep );
	ms.metadata["timestep"] = buffer;


	sprintf( buffer, "%d", pbc[0] );
	ms.metadata["pbcx"] = buffer;

	sprintf( buffer, "%d", pbc[1] );
	ms.metadata["pbcy"] = buffer;

	sprintf( buffer, "%d", pbc[2] );
	ms.metadata["pbcz"] = buffer;


	sprintf( buffer, "%g", bounds[0] );
	ms.metadata["minx"] = buffer;

	sprintf( buffer, "%g", bounds[1] );
	ms.metadata["maxx"] = buffer;


	sprintf( buffer, "%g", bounds[2] );
	ms.metadata["miny"] = buffer;

	sprintf( buffer, "%g", bounds[3] );
	ms.metadata["maxy"] = buffer;


	sprintf( buffer, "%g", bounds[4] );
	ms.metadata["minz"] = buffer;

	sprintf( buffer, "%g", bounds[5] );
	ms.metadata["maxz"] = buffer;


	return (int)ms.molecules.size();
}

int LAMMPS::Print( FILE *f, const MolecularSystem &ms )
{
	AttributeMap::const_iterator it;

	int timestep, N;
	double bounds[6], def_bound = 10.0;
	
	if( f == nullptr )
	{
		printf( "LammpsTrajectoryFrame::%s(): FILE pointer is nullptr.\n", __func__ );
		return -1;
	}
	
	timestep = -1;
	
	bounds[0] = -def_bound;
	bounds[1] = +def_bound;
	
	bounds[2] = -def_bound;
	bounds[3] = +def_bound;
	
	bounds[4] = -def_bound;
	bounds[5] = +def_bound;
	
	if( (it = ms.metadata.find("timestep") ) != ms.metadata.end() )
	{
		if( String::ToInteger(it->second.c_str(),timestep) != String::ReturnValue::OK )
		{
			printf( "MoleculeUtil::LAMMPS::%s(): unable to convert '%s' into timestep value\n", __func__, it->first.c_str() );
		}
	}

	if( (it = ms.metadata.find("minx") ) != ms.metadata.end() )
	{
		if( String::ToInteger(it->second.c_str(),bounds[0]) != String::ReturnValue::OK )
		{
			printf( "MoleculeUtil::LAMMPS::%s(): unable to convert '%s' into minx value\n", __func__, it->first.c_str() );
		}
	}
	if( (it = ms.metadata.find("maxx") ) != ms.metadata.end() )
	{
		if( String::ToInteger(it->second.c_str(),bounds[1]) != String::ReturnValue::OK )
		{
			printf( "MoleculeUtil::LAMMPS::%s(): unable to convert '%s' into maxx value\n", __func__, it->first.c_str() );
		}
	}

	if( (it = ms.metadata.find("miny") ) != ms.metadata.end() )
	{
		if( String::ToInteger(it->second.c_str(),bounds[2]) != String::ReturnValue::OK )
		{
			printf( "MoleculeUtil::LAMMPS::%s(): unable to convert '%s' into miny value\n", __func__, it->first.c_str() );
		}
	}
	if( (it = ms.metadata.find("maxy") ) != ms.metadata.end() )
	{
		if( String::ToInteger(it->second.c_str(),bounds[3]) != String::ReturnValue::OK )
		{
			printf( "MoleculeUtil::LAMMPS::%s(): unable to convert '%s' into maxy value\n", __func__, it->first.c_str() );
		}
	}

	if( (it = ms.metadata.find("minz") ) != ms.metadata.end() )
	{
		if( String::ToInteger(it->second.c_str(),bounds[4]) != String::ReturnValue::OK )
		{
			printf( "MoleculeUtil::LAMMPS::%s(): unable to convert '%s' into minz value\n", __func__, it->first.c_str() );
		}
	}
	if( (it = ms.metadata.find("maxz") ) != ms.metadata.end() )
	{
		if( String::ToInteger(it->second.c_str(),bounds[5]) != String::ReturnValue::OK )
		{
			printf( "MoleculeUtil::LAMMPS::%s(): unable to convert '%s' into maxz value\n", __func__, it->first.c_str() );
		}
	}
	

	N = 0;
	for( const auto& mol : ms.molecules )
	{
		N += (int)mol.size();
	}


	fprintf( f, "ITEM: TIMESTEP\n" );
	fprintf( f, "%d\n", timestep );

	fprintf( f, "ITEM: NUMBER OF ATOMS\n" );
	fprintf( f, "%d\n", N );

	fprintf( f, "ITEM: BOX BOUNDS ff ff ff\n" );
	fprintf( f, "%.3f %.3f\n", bounds[0], bounds[1] );
	fprintf( f, "%.3f %.3f\n", bounds[2], bounds[3] );
	fprintf( f, "%.3f %.3f\n", bounds[4], bounds[5] );

	/*
	if( write_velocities > 0 )
	{
		fprintf( f, "ITEM: ATOMS id mol type x y z vx vy vz \n" );
		for( int i=0; i<N; i++ ) fprintf( f, "%d %d %d %.6g %.6g %.6g %.6g %.6g %.6g \n", atoms[i].index, atoms[i].mol, atoms[i].type, atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].vx, atoms[i].vy, atoms[i].vz );
	}
	else
	*/
	
	{
		fprintf( f, "ITEM: ATOMS id mol type x y z \n" );
		
		int upto = 1;
		for( const auto& mol : ms.molecules )
		{
			for( const auto& a : mol )
			{
				int lmp_index = upto, lmp_mol = 1, lmp_type = 1;
				
				if( (it = a.attr.find("serial") ) != a.attr.end() )
				{
					if( String::ToInteger(it->second.c_str(),lmp_index) != String::ReturnValue::OK )
					{
						printf( "MoleculeUtil::LAMMPS::%s(): unable to convert '%s' into index value\n", __func__, it->first.c_str() );
					}
				}

				if( (it = a.attr.find("resSeq") ) != a.attr.end() )
				{
					if( String::ToInteger(it->second.c_str(),lmp_mol) != String::ReturnValue::OK )
					{
						printf( "MoleculeUtil::LAMMPS::%s(): unable to convert '%s' into mol value\n", __func__, it->first.c_str() );
					}
				}

				if( (it = a.attr.find("name") ) != a.attr.end() )
				{
					if( String::ToInteger(it->second.c_str(),lmp_type) != String::ReturnValue::OK )
					{
						printf( "MoleculeUtil::LAMMPS::%s(): unable to convert '%s' into type value\n", __func__, it->first.c_str() );
					}
				}
				
				fprintf( f, "%d %d %d %.6g %.6g %.6g \n", lmp_index, lmp_mol, lmp_type, a.x, a.y, a.z );

				upto++;
			}
		}
	}
	
	return 1;
	
}


}
}

