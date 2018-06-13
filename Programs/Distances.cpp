/*
	Author: John Grime, The University of Chicago.
*/

#include "../Util/Util.h"

using namespace Util;

//
// Measure specified distances between particles - designed for PDB files!
//

class Params
{
	public:
		
		std::vector< std::string > input_files;
		std::map< std::string, std::string > params;
		
		int Add( const char *s )
		{
			std::vector<std::string> tokens, subtoks;

			if( String::Tokenize( s, tokens, "=" ) < 2 ) return -1;

			//
			// Input data or other parameter? If former, sanity check number of params.
			//
			if( tokens[0] == "input" )
			{
				if( String::Tokenize( tokens[1].c_str(), subtoks, ":" ) < 2 ) return -1;
				input_files.push_back( std::string(s) );
			}
			else
			{
				params[ tokens[0] ] = tokens[1];
			}

			return 1;
		}
		const char *Get( const char *key ) const
		{
			const auto& it = params.find(key);
			if( it == params.end() ) return NULL;
			return it->second.c_str();
		}
};

class HistogramGenerator
{
	public:
	
	using Histogram = struct histogram_ {
		size_t N_samples;
		double min, max, delta;
		std::vector<int> counts;
		
		void Clear()
		{
			N_samples = 0;
			min = max = delta = 0.0;
			counts.clear();
		}
		void Resize( size_t N )
		{
			counts.resize( N );
			for( size_t i=0; i<N; i++ ) counts[i] = 0;
		}
		void AddSample( double r )
		{
			double u = (r-min)/(max-min);
			size_t bin = (size_t)floor( u*counts.size() );
			if( bin >= counts.size() )
			{
				printf( "HistogramGenerator::%s(): Bad bin index %d!\n", __func__, (int)bin );
				printf( "\t Generated from r = %g ( internal range: %g => %g ) => %g normalised\n", r, min,max, u );
				printf( "\t This corresponds to a bin index of %d ( 0 < index <= %d)\n", (int)bin, (int)counts.size() );
				exit( -1 );
			}
			counts[bin]++;
			N_samples++;
		}
		
	};
	
	static void Make( const std::vector<double>& vals, double bins_per_unit, Histogram& h )
	{
		h.Clear();
		
		for( size_t i=0; i<vals.size(); i++ )
		{
			double val = vals[i];
			
			if( i == 0 ) { h.min = h.max = val; }
			else
			{
				if( val < h.min ) h.min = val;
				if( val > h.max ) h.max = val;
			}
		}
		
		h.min = floor( h.min );
		h.max =  ceil( h.max );
		size_t N_bins = (size_t) floor( (h.max-h.min)*bins_per_unit );
		h.delta = (h.max-h.min)/N_bins;
		
		h.Resize( N_bins );
		for( size_t i=0; i<vals.size(); i++ ) h.AddSample( vals[i] );
	}
};

class Distances
{
	public:
		
		//
		// name:resSeq key for sorting atom info.
		//
		using Key = struct key_ {
			std::string name;
			int resSeq;
			
			bool operator < ( const key_ &rhs ) const
			{
				return (resSeq!=rhs.resSeq) ? (resSeq<rhs.resSeq) : (name<rhs.name);
			}
			bool operator == ( const key_ &k)
			{
				if( resSeq != k.resSeq ) return false; // fast, short-circuits need to check strings below
				if( name != k.name ) return false; // perform second, as likely slower than integer comparison
				return true;
			}
		};
		
		//
		// Pair of Keys, allows sorting on paired atom data (like distances between specified atoms)
		//
		using PairKey = struct pair_key_ {
			Key ki, kj;
			
			static void Generate( const Key &k1, const Key &k2, pair_key_ &pk )
			{
				bool k1_smaller = (k1 < k2);
				pk.ki = (k1_smaller) ? (k1) : (k2);
				pk.kj = (k1_smaller) ? (k2) : (k1);
			}
			bool operator < ( const pair_key_ &rhs ) const
			{
				if( ki.resSeq != rhs.ki.resSeq ) return (ki.resSeq < rhs.ki.resSeq);
				if( kj.resSeq != rhs.kj.resSeq ) return (kj.resSeq < rhs.kj.resSeq);
				
				if( ki.name != rhs.ki.name ) return (ki.name < rhs.ki.name);
				return (kj.name < rhs.kj.name);
			}
			bool operator == (const pair_key_ &pk )
			{
				if( (ki==pk.ki) && (kj==pk.kj) ) return true;
				return false;
			}
		};
		
		std::map< PairKey, std::vector<double> > dr_map;
				
		std::string name;
		int min_samples;
		bool skip_smallest;
		
		Distances()
		{
			Reset( "?", 1, false );
		}
		
		void Reset( const char* name_, int min_samples_, bool skip_smallest_ )
		{
			name = name_;
			min_samples = min_samples_;
			skip_smallest = skip_smallest_;
			dr_map.clear();
		}

		void AddDistances( Molecule::Molecule &mol1, Molecule::Molecule &mol2, double rcut );

		void AddDistances( const Distances& d );
		
		// Print raw distance info
		void Print();

		// For CG-style output
		void print_CG_nonbonded_line(
			const char *nb_name,
			const char *key1,
			const char *key2,
			double rcut,
			const Stats::Stats &stats,
			const std::vector<double> &values );

		// Print CG-style nonbonded distance info
		void PrintCG( const char *what_distances );

		// Prints CG-style nonbonded with equivalence replacements.
		void PrintCGEquivalances( const char *what_distances, const std::map< Key, Key > &m );
		
		// Debug!
		void PrintDistributions( const char* prefix, double bins_per_unit );
};
void Distances::AddDistances( Molecule::Molecule &mol1, Molecule::Molecule &mol2, double rcut )
{
	Molecule::AttributeMap::const_iterator attr_it;
	
	Key ki, kj;
	PairKey pk;
		
	double rcut2 = rcut*rcut;

	for( const auto& ai : mol1 )
	{
		//
		// Convert attribute strings from ai into a key, if found.
		//
		if( (attr_it = ai.attr.find("name")) == ai.attr.end() ) continue;
		ki.name = attr_it->second;

		if( (attr_it = ai.attr.find("resSeq")) == ai.attr.end() ) continue;
		if( (String::ToInteger(attr_it->second.c_str(),ki.resSeq)) != String::ReturnValue::OK ) continue;
		
		for( const auto& aj : mol2 )
		{
			double dx, dy, dz, dr2;

			dx = ai.x - aj.x;
			dy = ai.y - aj.y;
			dz = ai.z - aj.z;
			dr2 = dx*dx + dy*dy + dz*dz;
			
			if( dr2 < rcut2 )
			{
				//
				// Convert attribute strings from aj into a key, if found.
				//
				if( (attr_it = aj.attr.find("name")) == aj.attr.end() ) continue;
				kj.name = attr_it->second;

				if( (attr_it = aj.attr.find("resSeq")) == aj.attr.end() ) continue;
				if( (String::ToInteger(attr_it->second.c_str(),kj.resSeq)) != String::ReturnValue::OK ) continue;

				//
				// Generate paired key, and add this pair distance to the dr_map
				//
				PairKey::Generate( ki, kj, pk );
				dr_map[pk].push_back( sqrt(dr2) );
			}
		}
	}
}
void Distances::AddDistances( const Distances& d )
{
	std::vector<double> sorted_values;

	for( const auto& it : dr_map )
	{
		const auto& key = it.first;
		const auto& vals = it.second;

		// Use the min_samples / skip_smallest info from the other distance set!
		// This lets us combine a data set where we strip as appropriate in the sources
		// BEFORE the data is added!
		if( (int)vals.size() < d.min_samples ) continue;

		sorted_values = vals;
		sort( sorted_values.begin(), sorted_values.end() );

		size_t start = (d.skip_smallest==true) ? (1) : (0);
		for( size_t i = start; i<sorted_values.size(); i++ )
		{
			dr_map[ key ].push_back( sorted_values[i] );
		}
	}
}
void Distances::Print()
{
	std::vector<double> sorted_vals;
	char buffer[1024];

	Stats::Stats stats;
	
	printf( "%15.15s : %8.8s %8.8s %8.8s %8.8s %8.8s (%4s) Ascending list of distances ...\n",
		"name_i name_j",
		"min", "mean", "max",
		"stddev", "stderr",
		"N" );

	for( const auto& it : dr_map )
	{
		const PairKey &key = it.first;
		const std::vector<double> &vals = it.second;
		
		if( (int)vals.size() < min_samples ) continue;

		sorted_vals = vals;
		sort( sorted_vals.begin(), sorted_vals.end() );

		stats.Clear();
		size_t start = (skip_smallest==true) ? (1) : (0);
		for( size_t i = start; i<sorted_vals.size(); i++ ) stats.AddSample( sorted_vals[i] );

		sprintf( buffer, "%s:%d %s:%d", key.ki.name.c_str(), key.ki.resSeq, key.kj.name.c_str(), key.kj.resSeq );

		printf( "%15.15s : %8.3f %8.3f %8.3f %8.3f %8.3f (%4d) ",
			buffer,
			stats.min, stats.mean, stats.max,
			stats.StdDev(), stats.StdErr(), (int)stats.N );
		
		// Print ALL values, even if we're skipping the smallest!
		for( size_t j=0; j<sorted_vals.size(); j++ )
		{
			printf( "%8.3f ", sorted_vals[j] );
			if( j >= 4 ) break;
		}
		printf( "\n" );
	}
}
void Distances::print_CG_nonbonded_line(
	const char *nb_name,
	const char *key1,
	const char *key2,
	double rcut,
	const Stats::Stats &stats,
	const std::vector<double> &values )
{
	printf( "nonbonded  %10.10s %10.10s %10.10s  parameters %12.3f  # ", nb_name, key1, key2, rcut );
		
	// Print min, mean, max into output so we can easily swap between them in the file using column selects etc.
	printf( "%8d %8.3f %8.3f %8.3f : ", (int)stats.N, stats.min, stats.mean, stats.max );

	// Print first N values in the list for quick eyeball check (e.g. first very small vs rest).
	// This assumes values is sorted!
	for( size_t j=0; j<values.size(); j++ )
	{
		printf( "%8.3f ", values[j] );
		if( j >= 4 ) break;
	}
}
void Distances::PrintCG( const char *what_distances )
{
	//
	// Here we assume that the PDB keys [resSeq]:[name] => CA[resSeq] in CG naming terminology.
	//

	std::vector<double> sorted_vals;
	char ai_buf[128], aj_buf[128];

	Stats::Stats stats;
	
	std::string what = what_distances;
	
	printf( "#\n" );
	printf( "# Trailing data is # N_samples min mean max : val1 val2 ... \n" );
	printf( "#\n" );
	
	for( const auto& it : dr_map )
	{
		const auto& key = it.first;
		const auto& ki = key.ki;
		const auto& kj = key.kj;
		const auto& vals = it.second;
		
		if( (int)vals.size() < min_samples ) continue;

		sorted_vals = vals;
		sort( sorted_vals.begin(), sorted_vals.end() );

		stats.Clear();
		size_t start = (skip_smallest==true) ? (1) : (0);
		for( size_t i = start; i<sorted_vals.size(); i++ ) stats.AddSample( sorted_vals[i] );
		
		double r = (what == "min") ? (stats.min) : (stats.mean);
		
		sprintf( ai_buf, "${%s%d}", ki.name.c_str(), ki.resSeq );
		sprintf( aj_buf, "${%s%d}", kj.name.c_str(), kj.resSeq );
		
		print_CG_nonbonded_line( "${NB_NAME}", ai_buf, aj_buf, r, stats, sorted_vals );

		printf( "\n" );
	}
}
void Distances::PrintCGEquivalances( const char *what_distances, const std::map< Key, Key > &m )
{
	//
	// Here we assume that the PDB keys [resSeq]:[name] => CA[resSeq] in CG naming terminology.
	//

	std::vector<double> sorted_vals;
	char ai_buf[128], aj_buf[128];

	Stats::Stats stats;
	
	std::string what = what_distances;
	
	for( const auto& it : dr_map )
	{
		std::map< Key, Key >::const_iterator it2;
		const auto& key = it.first;
		const auto& ki = key.ki;
		const auto& kj = key.kj;
		const auto& vals = it.second;
		
		const Key *remapped_ki = NULL;
		const Key *remapped_kj = NULL;
		
		if( (it2 = m.find(ki)) != m.end() )
		{
			remapped_ki = &it2->second;
		}
		if( (it2 = m.find(kj)) != m.end() )
		{
			remapped_kj = &it2->second;
		}
		
		// If both NULL, skip.
		if( remapped_ki == NULL && remapped_kj == NULL )
		{
			continue;
		}
		
		if( (int)vals.size() < min_samples ) continue;

		sorted_vals = vals;
		sort( sorted_vals.begin(), sorted_vals.end() );

		stats.Clear();
		size_t start = (skip_smallest==true) ? (1) : (0);
		for( size_t i = start; i<sorted_vals.size(); i++ ) stats.AddSample( sorted_vals[i] );
		
		double r = (what == "min") ? (stats.min) : (stats.mean);
			
		if( remapped_kj != NULL )
		{
			sprintf( ai_buf, "${%s%d}", ki.name.c_str(), ki.resSeq );
			sprintf( aj_buf, "${%s%d}", remapped_kj->name.c_str(), remapped_kj->resSeq );
		
			print_CG_nonbonded_line( "${NB_NAME}", ai_buf, aj_buf, r, stats, sorted_vals );

			printf( " REMAP %s:%d => %s:%d\n", kj.name.c_str(), kj.resSeq, remapped_kj->name.c_str(), remapped_kj->resSeq );
		}
		
		if( remapped_ki != NULL )
		{
			sprintf( ai_buf, "${%s%d}", remapped_ki->name.c_str(), remapped_ki->resSeq );
			sprintf( aj_buf, "${%s%d}", kj.name.c_str(), kj.resSeq );
		
			print_CG_nonbonded_line( "${NB_NAME}", ai_buf, aj_buf, r, stats, sorted_vals );

			printf( " REMAP %s:%d => %s:%d\n", ki.name.c_str(), ki.resSeq, remapped_ki->name.c_str(), remapped_ki->resSeq );
		}

		if( remapped_ki != NULL && remapped_kj != NULL )
		{
			sprintf( ai_buf, "${%s%d}", remapped_ki->name.c_str(), remapped_ki->resSeq );
			sprintf( aj_buf, "${%s%d}", remapped_kj->name.c_str(), remapped_kj->resSeq );
		
			print_CG_nonbonded_line( "${NB_NAME}", ai_buf, aj_buf, r, stats, sorted_vals );

			printf( " REMAP %s:%d => %s:%d , %s:%d => %s:%d\n",
				ki.name.c_str(), ki.resSeq, remapped_ki->name.c_str(), remapped_ki->resSeq,
				kj.name.c_str(), kj.resSeq, remapped_kj->name.c_str(), remapped_kj->resSeq );
		}
	}
}
void Distances::PrintDistributions( const char* prefix, double bins_per_unit )
{
	char buffer[1024];
	HistogramGenerator::Histogram h;
	
	for( const auto& it : dr_map )
	{
		const auto& pk = it.first;
		const auto &vals = it.second;
		
		const auto& ki = pk.ki;
		const auto& kj = pk.kj;
		
		sprintf( buffer, "%s_%s%d_%s%d.txt",
			prefix,
			ki.name.c_str(), ki.resSeq,
			kj.name.c_str(), kj.resSeq );
		
			FILE *f = fopen( buffer, "w" );
			if( f == NULL )
			{
				fprintf( stderr, "Unable to open file '%s' for writing!\n", buffer );
				break;
			}
			
			HistogramGenerator::Make( vals, bins_per_unit, h );
			
			for( size_t i=0, N=h.counts.size(); i<N; i++ )
			{
				double pos = h.min + (0.5+i)*h.delta;
				int count = h.counts[i];
				double normalised = ((double)count) / h.N_samples;
				
				fprintf( f, "%e %d %e\n", pos, count, normalised );
			}
			
			fclose( f );
	}
}

void print_usage( const char *prog )
{
	printf( "Usage: %s input=name:filepath.pdb[:min_samples[:set_size]] input=name:filepath.pdb ... rcut=X [filters=\"filter_string;filter_string;...\"] [same=\"name:resSeq,name:resSeq;name:resSeq,name:resSeq;...\"] [histogram_prefix=X] [histogram_res=X]\n", prog );
	printf( "Where:\n" );
	printf( "  input : define an input PDB:\n" );
	printf( "    -name : name for input set.\n" );
	printf( "    -min_samples : OPTIONAL minimum number of samples required to print distances (default = 1).\n" );
	printf( "    -set_size : OPTIONAL only measure over consecutive 'set_size' mol groups (useful for structures superposed onto common reference frame, default = all in file).\n" );
	printf( "  filters : PDB-style filtering.\n" );
	printf( "  same : OPTIONAL definition of atoms to consider the same, to generate/prints additional distance info.\n" );
	printf( "  histogram_prefix : OPTIONAL prefix for saved histograms of distances; if not specified, no histograms written.\n" );
	printf( "  histogram_res : OPTIONAL resolution (bins per unit distance) for histograms (ignored if histogram_prefix not defined).\n" );
	printf( "Examples:\n" );
	printf( "  %s input=test:blah.pdb:4 rcut=10.0 filters=\"name:CA;resSeq:3,6,12-45,112-116\" same=\"CA:18,GCA:18;CA:45,GCA:45\" \n", prog );
	printf( "  %s input=test1:blah1.pdb:4 input=test2:blah2.pdb:4:2 rcut=10.0 filters=\"name:CA;resSeq:3,6,12-45,112-116\" same=\"CA:18,GCA:18;CA:45,GCA:45\" \n", prog );
	exit( -1 );
}

int main( int argc, char **argv )
{
	Molecule::AttributeFilter filter;
	double rcut;
	std::map< Distances::Key, Distances::Key > remap_info;
	bool ignore_smallest = false;

	Distances combined_distances;
	std::vector<Distances> distances_vec;
		
	FILE *f;
	
	size_t N_mols, N_filtered_atoms;
	
	Params params;
		
	if( argc < 3 ) print_usage( argv[0] );

	//
	// For combined data, no min_samples of filtering of smallest value
	// (we assume all that was done in the component data sets)
	//
	combined_distances.Reset( "Combined", 0, false );
	
	//
	// Add command line info
	//
	for( int i=1; i<argc; i++ )
	{
		params.Add( argv[i] );
	}

	
	if( params.input_files.size() < 1 ) print_usage( argv[0] );


	//
	// Show the user the info we think we have from the command line
	//
	printf( "Inputs:\n" );
	for( size_t i=0; i<params.input_files.size(); i++ )
	{
		printf( "\t '%s'\n", params.input_files[i].c_str() );
	}
	printf( "Parameters:\n" );
	for( auto it=params.params.begin(); it!=params.params.end(); it++ )
	{
		printf( "\t '%s' => '%s'\n", it->first.c_str(), it->second.c_str() );
	}
	

	//
	// Distance cutoff
	//
	if( params.Get( "rcut" ) == NULL ) print_usage( argv[0] );
	else
	{
		const char *val = params.Get( "rcut" );
		if( String::ToReal( val, rcut ) != String::ReturnValue::OK )
		{
			printf( "Bad rcut '%s'\n", val );
			exit( -1 );
		}
	}
	
	
	//
	// Process any specified filters.
	//
	if( params.Get( "ignore_smallest" ) != NULL ) ignore_smallest = true;


	//
	// Process any specified filters.
	//
	if( params.Get( "filters" ) != NULL )
	{
		std::vector<std::string> tokens;
		const char *val = params.Get( "filters" );
		
		String::Tokenize( val, tokens, ";" );
		for( size_t i=0; i<tokens.size(); i++ )
		{
			filter.AddFilter( tokens[i].c_str(), ":", ",", "-" );
		}
	}
	
	//
	// Equivalences defined?
	//
	if( params.Get( "same" ) != NULL )
	{
		Distances::Key k1, k2;
		std::vector<std::string> all_equiv_pairs, equiv_pair, name_resSeq;
		const char *val = params.Get( "same" );

		//
		// Add equivalence information
		//
		String::Tokenize( val, all_equiv_pairs, ";" );
		for( size_t i=0; i<all_equiv_pairs.size(); i++ )
		{
			// all_equiv_pairs[i] == "name:resSeq,name:resSeq"
			if( String::Tokenize( all_equiv_pairs[i].c_str(), equiv_pair, "," ) != 2 )
			{
				print_usage( argv[0] );
			}

			// equiv_pair[0] == "name:resSeq"
			if( String::Tokenize( equiv_pair[0].c_str(), name_resSeq, ":" ) != 2 )
			{
				print_usage( argv[0] );
			}
			k1.name = name_resSeq[0];
			if( String::ToInteger( name_resSeq[1].c_str(), k1.resSeq ) != String::ReturnValue::OK )
			{
				printf( "Bad equivalence pair '%s'\n", equiv_pair[0].c_str() );
				exit( -1 );
			}

			// equiv_pair[1] == "name:resSeq"
			if( String::Tokenize( equiv_pair[1].c_str(), name_resSeq, ":" ) != 2 )
			{
				print_usage( argv[0] );
			}
			k2.name = name_resSeq[0];
			if( String::ToInteger( name_resSeq[1].c_str(), k2.resSeq ) != String::ReturnValue::OK )
			{
				printf( "Bad equivalence pair '%s'\n", equiv_pair[1].c_str() );
				exit( -1 );
			}
			
			remap_info[ k1 ] = k2;
		}
	}

	//
	// Process input file stuff.
	//
	for( size_t ii=0; ii<params.input_files.size(); ii++ )
	{
		Molecule::MolecularSystem ms;
		Molecule::Molecules filtered_molecules;
		Distances distances;

		std::vector< std::string > tokens, subtoks;

		int set_size = -1;
		int min_samples = 1;

		String::Tokenize( params.input_files[ii].c_str(), tokens, "=" );
		String::Tokenize( tokens[1].c_str(), subtoks, ":" );

		const std::string& name = subtoks[0];
		const std::string& fpath = subtoks[1];
		
		//
		// Get min samples, if specified
		//
		if( subtoks.size() > 2 )
		{
			if( String::ToInteger( subtoks[2], min_samples ) != String::ReturnValue::OK )
			{
				printf( "Bad min_samples '%s' in '%s'\n", subtoks[2].c_str(), params.input_files[ii].c_str() );
				exit( -1 );
			}
		}

		//
		// Get set size, if specified
		//
		if( subtoks.size() > 3 )
		{
			if( String::ToInteger( subtoks[3], set_size ) != String::ReturnValue::OK )
			{
				printf( "Bad set_size '%s' in '%s'\n", subtoks[3].c_str(), params.input_files[ii].c_str() );
				exit( -1 );
			}
		}

		//
		// Reset data as appropriate to this input file
		//
		distances.Reset( name.c_str(), min_samples, ignore_smallest );

		//
		// Load the input data
		//
		if( (f=fopen(fpath.c_str(),"r")) == NULL )
		{
			printf( "Unable to open file '%s'\n", fpath.c_str() );
			exit( -1 );
		}
		Molecule::PDB::Load( f, ms );
		fclose( f );

		//
		// Filter the input data.
		//
		N_mols = ms.molecules.size();
		N_filtered_atoms = filter.Filter( ms.molecules, filtered_molecules );

		//
		// Print some info
		//
		printf( "\n" );
		printf( "*\n" );
		if( set_size < 1 )
		{
			printf( "* '%s' : '%s', min_samples %d, set_size <unspecified, using all molecules in file>\n", name.c_str(), fpath.c_str(), min_samples );
		}
		else
		{
			printf( "* '%s' : '%s', min_samples %d, set_size %d\n", name.c_str(), fpath.c_str(), min_samples, set_size );
		}
		if( ignore_smallest ) printf( "* IGNORING SMALLEST DISTANCE IN STATS\n" );
		printf( "*\n" );
		printf( "\n" );

		printf( "Read %d molecules.\n", (int)N_mols );
		printf( "%d atoms total passed filtering.\n", (int)N_filtered_atoms );
		
		//
		// Fix set size, if we need to.
		//
		if( set_size < 1 || set_size > (int)N_mols )
		{
			printf( "Adjusting set_size from %d to %d (number of filtered molecules)\n", set_size, (int)N_mols );
			set_size = (int)N_mols;
		}

		//
		// Iterate over input data, in chunks of set_size molecules. We only calculate distances between molecules
		// in a set. This is handy if we have e.g. sets of molecules superposed onto the same reference frame, as
		// we can specify the number of molecules in each set to prevent pollution from all the other superposed data.
		//
		for( size_t set_start=0; set_start<N_mols; set_start+=set_size )
		{
			size_t start = set_start;
			size_t stop = (start+set_size<N_mols) ? (start+set_size) : (N_mols);
			
			for( size_t i=start; i<stop; i++ )
			{
				if( (N_mols>=10) && i%(N_mols/10) == 0 ) printf( "Processing %d of %d ...\n", (int)i+1, (int)N_mols );
				for( size_t j=i+1; j<stop; j++ )
				{
					distances.AddDistances( filtered_molecules[i], filtered_molecules[j], rcut );
				}
			}
		}
			
		printf( "\n" );
		distances.Print();

		printf( "\n" );
		distances.PrintCG( "min" );
		
		//
		// Equivalences defined?
		//
		if( params.Get( "same" ) != NULL )
		{
			printf( "\n" );
			printf( "# \n" );
			printf( "# REMAPPED:\n" );
			for( auto it=remap_info.begin(); it!=remap_info.end(); it++ )
			{
				const auto& k1 = it->first;
				const auto& k2 = it->second;
				printf( "#   %s:%d => %s:%d\n", k1.name.c_str(),k1.resSeq, k2.name.c_str(),k2.resSeq );
			}
			printf( "#\n" );
			printf( "\n" );
			distances.PrintCGEquivalances( "min", remap_info );
		}
		
		//
		// Debug - print histograms
		//
		{
			double bins_per_unit = 2;
			const char* prefix = params.Get( "histogram_prefix" );
			const char* res = params.Get( "histogram_res" );
			
			if( prefix != NULL )
			{
				if( res != NULL )
				{
					if( String::ToReal( res, bins_per_unit ) != String::ReturnValue::OK )
					{
						printf( "Bad histogram_res value '%s'\n", res );
						exit( -1 );
					}
				}
				distances.PrintDistributions( prefix, bins_per_unit );
			}
		}

		combined_distances.AddDistances( distances );
		distances_vec.push_back( distances );
	}

	//
	// Combined distances should have been build from each input data set, only including
	// distances which pass the min_samples filter specified for each individual data set.
	// Therefore, use min_samples = 0 when printing, as we know the data has already passed
	// sample count filtering as appropriate to each input data set previously.
	//
	if( params.input_files.size() > 1 )
	{
		printf( "\n" );
		printf( "*\n" );
		printf( "* Combined data sets:\n" );
		printf( "*\n" );
		printf( "\n" );
		combined_distances.Print();
		printf( "\n" );
		combined_distances.PrintCG( "min" );
		
		//
		// Equivalences defined?
		//
		if( params.Get( "same" ) != NULL )
		{
			printf( "\n" );
			printf( "# \n" );
			printf( "# REMAPPED:\n" );
			for( auto it=remap_info.begin(); it!=remap_info.end(); it++ )
			{
				const auto& k1 = it->first;
				const auto& k2 = it->second;
				printf( "#   %s:%d => %s:%d\n", k1.name.c_str(),k1.resSeq, k2.name.c_str(),k2.resSeq );
			}
			printf( "#\n" );
			printf( "\n" );
			combined_distances.PrintCGEquivalances( "min", remap_info );
		}
		
		//
		// Debug - print histograms
		//
		char buffer[1024];
		double bins_per_unit = 2;
		const char* prefix = params.Get( "histogram_prefix" );
		const char* res = params.Get( "histogram_res" );
		
		if( prefix != NULL )
		{
			if( res != NULL )
			{
				if( String::ToReal( res, bins_per_unit ) != String::ReturnValue::OK )
				{
					printf( "Bad histogram_res value '%s'\n", res );
					exit( -1 );
				}
			}
			sprintf( buffer, "%s.combined", prefix );
			combined_distances.PrintDistributions( buffer, bins_per_unit );
		}
	}
	
	//
	// Test: attempt a combined table output?
	//
	{
		Stats::Stats stats;
		char buffer[512], buffer2[512];
		std::vector<double> sorted_vals;
		
		printf( "\n" );
		printf( "*\n" );
		printf( "* Side-by-side comparison\n" );
		if( ignore_smallest ) printf( "* ALL COLUMNS IGNORE SMALLEST DISTANCE IN DATA SETS\n" );
		printf( "*\n" );
		printf( "\n" );
		
		printf( "%17.17s ", "Pair" );
		printf( "%25.25s ", "Combined" );
		for( size_t i=0; i<distances_vec.size(); i++ )
		{
			printf( "%25.25s ", distances_vec[i].name.c_str() );
		}
		printf( "\n" );


		for( const auto& it : combined_distances.dr_map )
		{
			const auto& key = it.first;
			
			sprintf( buffer, "%s:%d", key.ki.name.c_str(), key.ki.resSeq );
			sprintf( buffer2, "%s:%d", key.kj.name.c_str(), key.kj.resSeq );
			printf( "%8.8s %8.8s ", buffer, buffer2 );
			
			//
			// Combined data ...
			//
			{
				const auto& values = it.second;
				stats.Clear();
				for( const auto v : values ) { stats.AddSample(v); }
				sprintf( buffer, "%8.3f %3.3d", stats.min, (int)stats.N );
				printf( "%25.25s ", buffer );
			}

			//
			// Iterate over input data...
			//
			for( size_t i=0; i<distances_vec.size(); i++ )
			{
				const auto& distances = distances_vec[i];
			
				sprintf( buffer, "%s", "..." );

				const auto& it2 = distances.dr_map.find( key );

				if( it2 != distances.dr_map.end() )
				{
					const auto& vals = it2->second;

					if( (int)vals.size() >= distances.min_samples )
					{
						sorted_vals = vals;
						sort( sorted_vals.begin(), sorted_vals.end() );
						
						stats.Clear();
						size_t start = (distances.skip_smallest==true) ? (1) : (0);
						for( size_t i = start; i<sorted_vals.size(); i++ ) stats.AddSample( sorted_vals[i] );

						sprintf( buffer, "%8.3f %3.3d", stats.min, (int)stats.N );
					}
				}
			
				printf( "%25.25s ", buffer );
			}
			printf( "\n" );
		}
		
	}

	return 0;
}
