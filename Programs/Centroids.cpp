/*
	Author: John Grime, The University of Chicago.
*/

#include "../Util/Util.h"

using namespace Util;


//
// Calculates centoids of point sets.
// Assumes inputs already superposed onto a common reference frame!
//

class AtomInfo
{
	public:
		
		//
		// Key for info map - sorts on resSeq, then resName, then atom name.
		//
		using Key = struct key_ {
			std::string name, resName;
			int resSeq;
			
			bool operator < ( const struct key_ &rhs ) const
			{
				if( resSeq != rhs.resSeq ) return (resSeq < rhs.resSeq);
				if( resName != rhs.resName ) return (resName < rhs.resName);
				return (name < rhs.name);
			}
		};
		
		//
		// Value for map - contains base atom info, and a list of coords from samples of that atom,
		// along with centroid and rmsf over the samples. Allows sorting on rmsf.
		//
		using Value = struct value_ {
			Molecule::Atom a; // atom this is based on

			std::vector<double> xyz; // set of xyz coords from all sampled data
			double centroid[3]; // centroid of above
			double rmsf; // rmsf of above
			
			bool operator < ( const struct value_ &rhs ) const
			{
				return (rmsf < rhs.rmsf);
			}
		};
				
		//
		// So we don't need to copy loads of data when sorting.
		//
		using Sort = struct sort_ {
			const void *ptr;
			double value;
			
			bool operator < ( const struct sort_ &rhs ) const
			{
				return (value < rhs.value);
			}
		};

		std::map< Key, Value > info;

		void AddAtom( const Molecule::Atom &a );
		void Print();
	
	protected:
						
		void get_info( const std::vector<double> &xyz, double *centroid, double &rmsf );
};
void AtomInfo::AddAtom( const Molecule::Atom &a_ )
{
	Key k;
	
	Molecule::Atom a = a_; // COPY, as we may modify the map below!
	
	k.name = a.attr["name"];
	k.resName = a.attr["resName"];
	String::ToInteger( a.attr["resSeq"], k.resSeq );
	
	auto it = info.find( k ); // used to check if we're adding a new atom, so we set the base atom below ...
	if( it == info.end() )
	{
		Value &v = info[k];

		v.a = a; // set base atom info for new entry.

		v.xyz.push_back( a.x );
		v.xyz.push_back( a.y );
		v.xyz.push_back( a.z );
		
	}
	else
	{
		Value &v = info[k];
		
		v.xyz.push_back( a.x );
		v.xyz.push_back( a.y );
		v.xyz.push_back( a.z );
	}
	
}
void AtomInfo::Print()
{
	std::vector<Sort> sorted_info;
	char buffer[128];
	
	//
	// Build sorted information vector from the info map
	//
	printf( "# Centroid atom info (n_samples): x y z RMSF:\n" );
	sorted_info.clear();
	for( auto& it : info )
	{
		Sort s;
		
		auto &k = it.first;
		auto &v = it.second;
		
		get_info( v.xyz, v.centroid, v.rmsf );
		
		
		s.ptr = (const void *)&k;
		s.value = v.rmsf;
		sorted_info.push_back( s );
		
		printf( "%4.4s %4.4s %5.5s (%3d): %8.3f %8.3f %8.3f  %8.3f\n",
			v.a.attr["name"].c_str(), v.a.attr["resName"].c_str(), v.a.attr["resSeq"].c_str(),
			(int)v.xyz.size()/3,
			v.centroid[0], v.centroid[1], v.centroid[2],
			v.rmsf );
	}
	sort( sorted_info.begin(), sorted_info.end() );
	printf( "\n" );

	printf( "# Sorted centroid atom info (n_samples): x y z RMSF:\n" );
	for( const auto& si : sorted_info )
	{
		const Key *k = (const Key *)si.ptr;
		Value &v = info[*k];

		printf( "%4.4s %4.4s %5.5s (%3d): %8.3f %8.3f %8.3f  %8.3f\n",
			v.a.attr["name"].c_str(), v.a.attr["resName"].c_str(), v.a.attr["resSeq"].c_str(),
			(int)v.xyz.size()/3,
			v.centroid[0], v.centroid[1], v.centroid[2],
			v.rmsf );
	}
	printf( "\n" );

	printf( "REMARK sorted PDB data: coords are centroids, occupancy is N_samples, tempFactor is RMSF (ie sort key)\n" );
	for( const auto& si : sorted_info )
	{
		const Key *k = (const Key *)si.ptr;
		Value &v = info[*k];

		v.a.x = v.centroid[0];
		v.a.y = v.centroid[1];
		v.a.z = v.centroid[2];

		sprintf( buffer, "%d", (int)v.xyz.size()/3 );
		v.a.attr["occupancy"] = buffer;
		
		sprintf( buffer, "%6.2f", v.rmsf );
		v.a.attr["tempFactor"] = buffer;
		
		Molecule::PDB::Print( stdout, v.a );
		fprintf( stdout, "\n" );
	}
	
}
void AtomInfo::get_info( const std::vector<double> &xyz, double *centroid, double &rmsf )
{
	size_t N;
	
	centroid[0] = 0.0;
	centroid[1] = 0.0;
	centroid[2] = 0.0;
	rmsf = 0.0;
	
	N = xyz.size()/3;
	
	if( N == 0 ) return;

	//
	// centroid (we also use this for RMSF)
	//
	for( size_t i=0; i<N; i++ )
	{
		centroid[0] += xyz[i*3 +0];
		centroid[1] += xyz[i*3 +1];
		centroid[2] += xyz[i*3 +2];
	}
	centroid[0] /= N;
	centroid[1] /= N;
	centroid[2] /= N;
	
	for( size_t i=0; i<N; i++ )
	{
		double dx, dy, dz;
		
		dx = centroid[0] - xyz[i*3 +0];
		dy = centroid[1] - xyz[i*3 +1];
		dz = centroid[2] - xyz[i*3 +2];
		rmsf += dx*dx + dy*dy + dz*dz;
	}
	rmsf = sqrt( rmsf/N );
}


int main( int argc, char **argv )
{
	const char *fpath;
	FILE *f;

	Molecule::MolecularSystem ms;

	int set_size;
	std::vector<AtomInfo> ai_vec;
	
	char buffer[256];

	if( argc < 2 )
	{
		printf( "Usage: %s input.pdb set_size\n", argv[0] );
		printf( "Where:\n" );
		printf( "  - set_size : number of consecutive PDB molecules in the input file to group\n" );
		exit( -1 );
	}

	fpath = argv[1];
	set_size = atoi( argv[2] );
	
	ai_vec.resize( set_size );
	
	//
	// Walk file in chunks of set_size, processing as we go.
	//
	if( (f=fopen(fpath,"r")) == NULL )
	{
		printf( "Unable to open file '%s'\n", fpath );
		exit( -1 );
	}
	int set_number = 0;
	while( true )
	{
		set_number++;
		Molecule::PDB::Load( f, ms, set_size );
		
		if( (int)ms.molecules.size() < set_size )
		{
			printf( "Unable to load %d molecules from file '%s' for set number %d (read %d): stopping here.\n",
				set_size,
				fpath,
				set_number,
				(int)ms.molecules.size() );
			break;
		}

		for( int i=0; i<set_size; i++ )
		{
			AtomInfo &ai = ai_vec[i];
			Molecule::Molecule &molecule = ms.molecules[i];
			for( const auto& a : molecule ) ai.AddAtom( a );
		}
	}
	fclose( f );
	
	//
	// Print some information to stdout ...
	//
	for( int i=0; i<set_size; i++ )
	{
		AtomInfo &ai = ai_vec[i];

		printf( "\n" );
		printf( "#\n" );
		printf( "# Molecule %d in set of %d ...\n", i+1, set_size );
		printf( "#\n" );
		printf( "\n" );
		ai.Print();
	}
	

	//
	// Save centroids
	//
	sprintf( buffer, "output.centroids.pdb" );
	fpath = buffer;
	if( (f=fopen(fpath,"w")) == NULL )
	{
		printf( "Unable to open file '%s'\n", fpath );
		exit( -1 );
	}
	fprintf( f, "REMARK sorted data: coords are centroids, occupancy is N_samples, tempFactor is RMSF (ie sort key)\n" );
	for( int i=0; i<set_size; i++ )
	{
		AtomInfo &ai = ai_vec[i];

		for( auto it=ai.info.begin(); it!=ai.info.end(); it++ )
		{
			Molecule::PDB::Print( f, it->second.a );
			fprintf( f, "\n" );
		}
		fprintf( f, "TER\n" );
	}
	fclose( f );
	
	return 0;
}
