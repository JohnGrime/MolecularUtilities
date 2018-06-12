/*
	Author: John Grime, The University of Chicago.
*/

#include <algorithm>
#include <numeric>
#include <random>

#include "../Util/Util.h"

using namespace Util;

const char* ATOM_format = "%-6.6s%5.5s %4.4s%1.1s%3.3s %1.1s%4d%1.1s   %8.2f%8.2f%8.2f%6.6s%6.6s          %2.2s%2.2s\n";		

//
// Generate a spherical surface using icosahedral subdivision, with specific bonds to maintain the structure.
// Then we graft on the specified polymer chains.
//

//
// Sphere metadata
//
struct SphereInfo
{
	double radius;
	int n_subdiv, atom_type, bond_type;
	double rcut;
	std::string atmName, resName;

	SphereInfo()
	{
		radius = -1.0;
		n_subdiv = 0;
		atom_type = 1;
		bond_type = 1; // start bond type ids at 1
		rcut = -1.0;
		atmName = resName = "SPH";
	}

	//
	// sphere=r;n_subdiv;type;rcut[;atmName][;resName]
	//
	// Prefix token ("sphere=") is stripped before passing in str.
	//
	int Parse( const char* str )
	{
		std::vector<std::string> toks, subtoks;

		if( String::Tokenize( str, toks, ":" ) < 4 ) return -1;

		//
		// Parse radius
		//
		if( String::ToReal( toks[0], radius ) != String::ReturnValue::OK ) return -1;
		if( radius <= 0.0 ) return -1;

		//
		// Parse n_subdiv
		//
		if( String::ToInteger( toks[1], n_subdiv ) != String::ReturnValue::OK ) return -1;
		if( n_subdiv < 0 ) return -1;

		//
		// Parse atom type offset
		//
		if( String::ToInteger( toks[2], atom_type ) != String::ReturnValue::OK ) return -1;
		if( atom_type <= 0 ) return -1;

		//
		// Parse bond type offset
		//
		if( String::ToInteger( toks[3], bond_type ) != String::ReturnValue::OK ) return -1;
		if( bond_type <= 0 ) return -1;

		//
		// Parse rcut
		//
		if( String::ToReal( toks[4], rcut ) != String::ReturnValue::OK ) return -1;
		if( rcut <= 0.0 ) return -1;

		//
		// Parse optional atom and res names
		//
		if( toks.size() > 5 ) atmName = toks[5];
		if( toks.size() > 6 ) resName = toks[6];

		return 1;
	}

	void Print() const
	{
		printf( "Sphere:\n" );
		printf( "  radius = %g\n", radius );
		printf( "  n_subdiv = %d\n", n_subdiv );
		printf( "  atom_type = %d\n", atom_type );
		printf( "  bond_type = %d\n", bond_type );
		printf( "  rcut = %g\n", rcut );
		printf( "  atom & res names = %s %s\n", atmName.c_str(), resName.c_str() );
	}
};

//
// Grafted polymer metadata
//
struct PolymerInfo
{
	int N;
	std::vector<int> atom_types;
	double dr;

	PolymerInfo()
	{
		atom_types.clear();
		dr = 1.0;
	}

	//
	// polymer=N;t1,t2,t3,...tN;dr
	//
	// Prefix token ("polymer=") is stripped before passing in str.
	//
	int Parse( const char* str )
	{
		std::vector<std::string> toks, subtoks;

		if( String::Tokenize( str, toks, ":" ) < 3 ) return -1;

		//
		// Parse count
		//
		if( String::ToInteger( toks[0], N ) != String::ReturnValue::OK ) return -1;
		if( N < 0 ) return -1;

		//
		// Parse atom types
		//
		atom_types.clear();
		String::Tokenize( toks[1], subtoks, "," );
		for( const auto& t : subtoks )
		{
			int what;
			if( String::ToInteger( t, what ) != String::ReturnValue::OK ) return -1;
			if( what < 0 ) return -1;
			atom_types.push_back( what );
		}

		//
		// Parse delta
		//
		if( String::ToReal( toks[2], dr ) != String::ReturnValue::OK ) return -1;
		if( dr <= 0.0 ) return -1;

		return 1;
	}

	void Print() const
	{
		printf( "Polymer:\n" );
		printf( "  N = %d\n", N );
		printf( "  atom types = " );
		for( const auto t : atom_types ) printf( "%d ", t );
		printf( "\n" );
		printf( "  dr = %g\n", dr );
	}
};

//
// Keep track of all the topological types (bonds, angles etc) we're using. Note: this does not
// track individual instance of the bonds etc - it tracks the TYPES!
//
struct TopologyInfo
{
	using Topo = struct topo_ {
		int type;
		double rest;
		std::string desc;
	};

	std::map<int,Topo> bonds, angles;

	void AddBond( int type, double rest, const char* desc )
	{
		if( bonds.find(type) != bonds.end() ) return;
		bonds[ type ] = { type, rest, desc };
	}
	int NextAvailableBondType() const
	{
		int t = 0;
		for( const auto& it : bonds )
		{
			if( it.first > t ) t = it.first;
		}
		return t+1;
	}

	void AddAngle( int type, double rest, const char* desc )
	{
		if( angles.find(type) != angles.end() ) return;
		angles[ type ] = { type, rest, desc };
	}
	int NextAvailableAngleType() const
	{
		int t = 0;
		for( const auto& it : angles )
		{
			if( it.first > t ) t = it.first;
		}
		return t+1;
	}

	void Print() const
	{
		if( bonds.size() > 0 )
		{
			printf( "#\n" );
			printf( "# %d bond types:\n", (int)bonds.size() );
			printf( "#\n" );
			for( const auto& it : bonds )
			{
				const auto& t = it.second;
				printf( "bond_coeff  %5d  K  %10g # %s\n", t.type, t.rest, t.desc.c_str() );
			}
			printf( "\n" );
		}

		if( angles.size() > 0 )
		{
			printf( "#\n" );
			printf( "# %d angle types:\n", (int)angles.size() );
			printf( "#\n" );
			for( const auto& it : angles )
			{
				const auto& t = it.second;
				printf( "angle_coeff  %5d  K  # rest: %g, %s\n", t.type, t.rest, t.desc.c_str() );
			}
			printf( "\n" );
		}
	}
};

//
// Save a LAMMPS molecule in PDB format.
//
void SaveMolecule( FILE* f, const LAMMPS::Molecule&mol )
{
	char buffer[16], serial[16];

	if( f == nullptr ) return;

	for( size_t i=0; i<mol.atoms.size(); i++ )
	{
		const auto& a = mol.atoms[i];

		sprintf( serial, "%d", (int)i+1 );
		sprintf( buffer, "%d", a.type );

		fprintf( f,
			ATOM_format,
			"ATOM",
			serial,
			buffer,
			"",
			buffer,
			"A",
			1,
			"",
			a.x, a.y, a.z,
			"", "", "", "" );
	}
	fprintf( f, "TER" ); // <- no EOL here; see below.

	int last_i = -1, count = 0; // "count" used to control 4 entries per CONECT line
	for( const auto& b : mol.bonds )
	{
		int i = b.i;
		int j = b.j;

		if( i != last_i || count >= 4 )
		{
			fprintf( f, "\nCONECT%5d", i+1 );
			count = 0;
		}
		last_i = i;

		fprintf( f, "%5d", j+1 );
		count++;
	}
	fprintf( f, "\n" );
}

//
// Find the closest point in specified hi_res data to the points in lo_res data, subject
// to the specified cutoff. Index array values of -1 denote that entry should be ignored.
//
struct ClosestPairs
{
	struct Pair
	{
		int i, j;
		double dr;
	};
	std::vector<Pair> pairs;

	template<typename T1, typename T2>
	int Generate(
		std::vector<int>& lo_res_idx,
		const std::vector<T1>& lo_res_xyz,
		std::vector<int>& hi_res_idx,
		const std::vector<T2>& hi_res_xyz,
		double rcut )
	{
		pairs.clear();

		for( const auto i : lo_res_idx )
		{
			if( i == -1 ) continue;

			auto xi = lo_res_xyz[i*3 +0];
			auto yi = lo_res_xyz[i*3 +1];
			auto zi = lo_res_xyz[i*3 +2];

			int cjj = -1;
			double check2 = 0.0;

			for( size_t jj=0; jj<hi_res_idx.size(); jj++ )
			{
				int j = hi_res_idx[jj];
				if( j == -1 ) continue;

				double dx = xi - hi_res_xyz[j*3 +0];
				double dy = yi - hi_res_xyz[j*3 +1];
				double dz = zi - hi_res_xyz[j*3 +2];
				double dr2 = dx*dx + dy*dy + dz*dz;

				if( (cjj==-1) || (dr2<check2) )
				{
					cjj = jj;
					check2 = dr2;
				}
			}

			int j = hi_res_idx[cjj];
			double check = sqrt(check2);
			if( check > rcut )
			{
				printf( "Unable to find partner for %d; closest was %d (%g)\n", i, j, check );
				continue;
			}

			pairs.push_back( {i,j,check} );
			hi_res_idx[cjj] = -1; // ignore this hi res point in future iterations.
		}
		return 1;
	}
};

//
// Generate subdivided sphere points, along with bonds to maintain structure. I've found that
// the bonds are typically slightly different lengths for a given subdivision stage, so I
// quantise the bond lenths to provide better structural reproduction while limiting the
// number of different bond types we need. Locking to the single average bond length seems
// to struggle to preserve an actual sphere - instead we end up with quasi-icosahedral thingies
// in simulation.
//
void GenerateSphere(
	const SphereInfo& sphere,
	const double quantum,
	LAMMPS::Molecule& mol,
	TopologyInfo& topo_info )
{
	Stats::Stats stats;

	//
	// Quantise bond lengths to 1% of the cutoff distance - seem to work okay!
	//
	std::map< int, int > quantised_to_type;

	//
	// PEr-stage bond info
	//
	std::vector< std::pair<int,int> > bonds;
	std::vector<int> types;
	std::vector<double> distances;

	Subdivider subdiv;
	auto& xyz = subdiv.vertices;

	char buffer[1024];

	double cut2 = sphere.rcut*sphere.rcut;
	double ratio = 1.25; // +25% of min rij between "new" beads seems to work well for per-pass bond cutoff

	subdiv.Init( Subdivider::InitType::Icosahedron );

	size_t start_point = 0; // start index for new beads added in current pass
	for( int pass=0; pass<sphere.n_subdiv; pass++ )
	{
		subdiv.Subdivide();
		size_t N_beads = xyz.size()/3;

		quantised_to_type.clear(); // new types for bonds in each pass, even if same quantised distance!

		bonds.clear();
		types.clear();
		distances.clear();

		stats.Clear();

		printf( "Pass %d ...\n", pass );

		//
		// Add new particles to molecule, and find min separation in new points.
		//
		double min_sep2 = (sphere.radius*2.0)*(sphere.radius*2.0);
		for( size_t i=start_point; i<N_beads; i++ )
		{
			double xi = xyz[i*3 +0] * sphere.radius;
			double yi = xyz[i*3 +1] * sphere.radius;
			double zi = xyz[i*3 +2] * sphere.radius;

			mol.AddAtom( i, sphere.atom_type, 0.0, xi, yi, zi );

			for( size_t j=i+1; j<N_beads; j++ )
			{
				double dx = xi - (xyz[j*3 +0] * sphere.radius);
				double dy = yi - (xyz[j*3 +1] * sphere.radius);
				double dz = zi - (xyz[j*3 +2] * sphere.radius);
				double dr2 = dx*dx + dy*dy + dz*dz;

				if( dr2 < min_sep2 )
				{
					min_sep2 = dr2;
				}
			}
		}

		//
		// Calculate bond cutoff for this pass
		//
		{
			min_sep2 = sqrt( min_sep2 );

			printf( "  bond cutoff is %g x %g = %g\n",
				min_sep2, ratio, ratio*min_sep2 );

			min_sep2 = ratio * min_sep2;
			min_sep2 = min_sep2 * min_sep2;
		}

		//
		// Find new bonds in this pass
		//
		for( size_t i=0; i<N_beads; i++ )
		{
			double xi = xyz[i*3 +0] * sphere.radius;
			double yi = xyz[i*3 +1] * sphere.radius;
			double zi = xyz[i*3 +2] * sphere.radius;

			size_t j0 = i+1;
			if( j0 < start_point ) j0 = start_point;

			for( size_t j=j0; j<N_beads; j++ )
			{
				double dx = xi - (xyz[j*3 +0] * sphere.radius);
				double dy = yi - (xyz[j*3 +1] * sphere.radius);
				double dz = zi - (xyz[j*3 +2] * sphere.radius);
				double dr2 = dx*dx + dy*dy + dz*dz;

				if( (dr2<cut2) && (dr2<min_sep2) )
				{
					double dr = sqrt( dr2 );
					int type;

					//
					// Figure out what bond type to use from the quantised length!
					//
					{
						int quantised = (int)round( dr / quantum );
						const auto& it = quantised_to_type.find( quantised );
						if( it == quantised_to_type.end() )
						{
							type = topo_info.NextAvailableBondType();
							if( type < sphere.bond_type ) type = sphere.bond_type;

							quantised_to_type[quantised] = type;

							double rest = (0.5+quantised)*quantum;
							sprintf( buffer, "Sphere pass %d, quantised bin %3d (quantum %g) => distance %g", pass, quantised, quantum, rest );
							topo_info.AddBond( type, rest, buffer );
							printf( "Registered new quantised bond: %s, type is %d\n", buffer, type );
						}
						else
						{
							type = it->second;
						}
					}

					bonds.push_back( {i,j} );
					types.push_back( type );

					distances.push_back( dr );
					stats.AddSample( dr );
				}
			}
		}

		double mean = stats.mean;
		double stddev = stats.StdDev();
		double stderr = stats.StdErr();

		if( bonds.size() > 0 )
		{
			//
			// Add new bonds to molecule
			//
			for( size_t i=0; i<bonds.size(); i++ )
			{
				LAMMPS::Bond bond;

				bond.type = types[i];
				bond.i = bonds[i].first;
				bond.j = bonds[i].second;

				mol.bonds.push_back( bond );
			}

			printf( "  %d beads, %d bonds (%d bond types).\n", (int)N_beads, (int)bonds.size(), (int)quantised_to_type.size() );
		}
		else
		{
			printf( "  %d beads.\n", (int)N_beads );
		}

		//
		// Save current pass information.
		//
		{
			FILE *f;

			sprintf( buffer, "pass.%d.pdb", pass );
			f = fopen( buffer, "w" );
			if( f == nullptr )
			{
				printf( "Unable to open file '%s' for writing\n", buffer );
				exit( -1 );
			}
			SaveMolecule( f, mol );
			fclose( f );

			sprintf( buffer, "pass.%d.rij", pass );
			f = fopen( buffer, "w" );
			if( f == nullptr )
			{
				printf( "Unable to open file '%s' for writing\n", buffer );
				exit( -1 );
			}
			fprintf( f, "# N %d mean %g stddev %g stderr %g\n", (int)stats.N, mean, stddev, stderr );
			for( size_t i=0; i<bonds.size(); i++ )
			{
				const auto b = bonds[i];
				fprintf( f, "%d %d %g\n", b.first, b.second, distances[i] );
			}
			fclose( f );
		}

		start_point = N_beads;
	}
}

//
// Add an instance of the specified polymer to a LAMMPS molecule data.
//
template<typename T>
void GeneratePolymer(
	T sx, T sy, T sz,
	T dx, T dy, T dz,
	const PolymerInfo& polymer,
	LAMMPS::Molecule& mol,
	int bond_type_offset,
	int angle_type_offset )
{
	//
	// Normalize deltas, scale by dr
	//
	double dr = sqrt( dx*dx + dy*dy + dz*dz );
	dx /= dr;
	dy /= dr;
	dz /= dr;
	dx *= polymer.dr;
	dy *= polymer.dr;
	dz *= polymer.dr;

	//
	// Relative offset into molecule & location of first bead
	//
	int bead_offset = (int)mol.atoms.size();
	double x = sx + dx;
	double y = sy + dy;
	double z = sz + dz;

	//
	// Add new polymer beads to the molecule
	//
	for( size_t ii=0; ii<polymer.atom_types.size(); ii++ )
	{
		mol.AddAtom( bead_offset+(int)ii, polymer.atom_types[ii], 0.0, x,y,z );

		x += dx;
		y += dy;
		z += dz;
	}

	//
	// Add new polymer bonds / angles to the molecule
	//

	for( size_t ii=0; ii<polymer.atom_types.size()-1; ii++ )
	{
		int i = bead_offset + (int)ii;
		int j = bead_offset + (int)ii+1;
		int type = bond_type_offset + (int)ii;
	
		mol.AddBond( type, i, j );
	}

	for( size_t ii=0; ii<polymer.atom_types.size()-2; ii++ )
	{
		int i = bead_offset + (int)ii;
		int j = bead_offset + (int)ii+1;
		int k = bead_offset + (int)ii+2;
		int type = angle_type_offset + (int)ii;
		
		mol.AddAngle( type, i, j, k );
	}
}

void print_usage( const char *prog )
{
	printf( "Usage: %s  sphere=r:n_subdiv:atom_type:bond_type_start:rcut[:atmName][resName]  [bond_quantum=X]  [polymer=N:t1,t2,t3,...tN:dr ...]\n", prog );
	printf( "\n" );
	printf( "Where:\n" );
	printf( "\n" );
	printf( "  sphere.r = radius\n" );
	printf( "  sphere.n_subdiv = number of icosahedral subdivisions in generation\n" );
	printf( "  sphere.atom_type = LAMMPS type of sphere beads\n" );
	printf( "  sphere.bond_type_start = bond type start\n" );
	printf( "  sphere.rcut = cutoff for bond detection in subdivision\n" );
	printf( "  sphere.atmName = OPTIONAL PDB name for sphere atoms (default is 'SPH')\n" );
	printf( "  sphere.resName = OPTIONAL PDB resName for sphere atoms (default is 'SPH')\n" );
	printf( "\n" );
	printf( "  bond_quantum = quatization unit for bond generation in sphere (default is 1%% of sphere.rcut)\n" );
	printf( "\n" );
	printf( "  polymer.N = number of instances of this specific polymer type on the sphere surface\n" );
	printf( "  polymer.t1,...tN = LAMMPS bead types for the polymer\n" );
	printf( "  polymer.dr = separation between adjacent beads in the polymer\n" );
	printf( "\n" );
	printf( "Notes:" );
	printf( "\n" );
	printf( "  Zero or more polymer types can be defined. They are randomly placed on the sphere.\n" );
	printf( "  Bond/angle types in polymers continue from sphere types.\n" );
	printf( "\n" );
	exit( -1 );
}

int main( int argc, char **argv )
{
	char buffer[1024];
	std::vector<std::string> toks;

	SphereInfo sphere_info;
	std::vector<PolymerInfo> polymer_info;
	TopologyInfo topo_info;

	double bond_quantum = -1;
	std::vector<double> hi_res_xyz, lo_res_xyz;

	ClosestPairs cp;

	LAMMPS::Config config;

	if( argc < 2 )
	{
		print_usage( argv[0] );
	}

	//
	// Set up an empty molecule in config to get us started.
	//
	{
		LAMMPS::Molecule mol;
		config.molecules.push_back( mol );
	}

	LAMMPS::Molecule& mol = config.molecules[0];

	//
	// Parse sphere & polymer definitions
	//
	{
		for( int i=0; i<argc; i++ )
		{
			String::Tokenize( argv[i], toks, "=" );
			
			if( toks.size() < 2 ) continue;

			if( toks[0] == "sphere" )
			{
				if( sphere_info.Parse( toks[1].c_str() ) == -1 )
				{
					printf( "Unable to parse sphere definition '%s'\n", argv[i] );
					exit( -1 );
				}
			}

			if( toks[0] == "polymer" )
			{
				PolymerInfo p;
				if( p.Parse( toks[1].c_str() ) == -1 )
				{
					printf( "Unable to parse polymer definition '%s'\n", argv[i] );
					exit( -1 );
				}
				polymer_info.push_back( p );
			}			

			if( toks[0] == "bond_quantum" )
			{
				if( String::ToReal( toks[1].c_str(), bond_quantum ) != String::ReturnValue::OK )
				{
					printf( "Unable to parse bond_quantum definition '%s'\n", argv[i] );
					exit( -1 );
				}
			}
		}

		if( sphere_info.radius <= 0.0 )
		{
			printf( "No sphere definition found\n" );
			exit( -1 );
		}
		
		printf( "\n" );
		printf( "*\n" );
		printf( "* System parameters ...\n" );
		printf( "*\n" );
		printf( "\n" );

		sphere_info.Print();
		for( const auto& p : polymer_info ) p.Print();
	}

	if( bond_quantum < 0.0 ) bond_quantum = 0.01 * sphere_info.rcut;

	//
	// Generate sphere part of molecule
	//
	{
		printf( "\n" );
		printf( "*\n" );
		printf( "* Subdividing & detecting bonds ...\n" );
		printf( "*\n" );
		printf( "\n" );

		GenerateSphere( sphere_info, bond_quantum, mol, topo_info );

		printf( "\n" );
		printf( "Sphere component has %d particles, %d bonds.\n", (int)mol.atoms.size(), (int)mol.bonds.size() );

		hi_res_xyz.clear();
		for( size_t i=0; i<mol.atoms.size(); i++ )
		{
			const auto& a = mol.atoms[i];
			hi_res_xyz.push_back( a.x );
			hi_res_xyz.push_back( a.y );
			hi_res_xyz.push_back( a.z );
		}
	}

	//
	// Figure out how many polymer graft location we need, and pair them up to the high-res sphere sites.
	//
	{
		int N_required = 0;

		//
		// Get randomly placed points on sphere surface for polymer grafting.
		//
		{
			printf( "\n" );
			printf( "*\n" );
			printf( "* Generating polymer graft locations & mapping onto sphere vertices ...\n" );
			printf( "*\n" );
			printf( "\n" );

			for( const auto& p : polymer_info ) N_required += p.N;
			
			printf( "Minimum of %d polymer graft locations needed.\n", N_required );

			if( (int)hi_res_xyz.size() < N_required*3 )
			{
				printf( "ERROR: Too few possible insertion points for polymers!\n" );
				exit( -1 );
			}

			int tmp = N_required;
			while( Geometry::GetEquispacedSpherePoints(
					tmp, lo_res_xyz, sphere_info.radius ) < N_required ) { tmp++; }

			printf( "Generated %d potential graft locations.\n", (int)lo_res_xyz.size()/3 );

			lo_res_xyz.resize( N_required*3 );
		}

		//
		// Get points in the high-res sphere closest to the specified low-res sphere points.
		// These are the beads onto which the polymers will be grafted.
		//
		{
			std::vector<int> lo_res_idx, hi_res_idx;
			auto rng = std::mt19937( 666 );
			double pairing_rcut = sphere_info.radius*2.0;

			for( size_t i=0; i<hi_res_xyz.size()/3; i++ ) hi_res_idx.push_back( i );
			for( size_t i=0; i<lo_res_xyz.size()/3; i++ ) lo_res_idx.push_back( i );

			std::shuffle( begin(lo_res_idx), end(lo_res_idx), rng );
			cp.Generate( lo_res_idx, lo_res_xyz, hi_res_idx, hi_res_xyz, pairing_rcut );

			if( (int)cp.pairs.size() < N_required )
			{
				printf( "ERROR: Not enough paired locations; need %d, got %d\n", N_required, (int)cp.pairs.size() );
				exit( -1 );
			}
		}
	}

	//
	// Generate grafted polymer data
	//
	{
		int graft_index = 0;

		printf( "\n" );
		printf( "*\n" );
		printf( "* Grafting polymers onto sphere ...\n" );
		printf( "*\n" );
		printf( "\n" );

		for( size_t poly_type=0; poly_type<polymer_info.size(); poly_type++ )
		{
			const auto& p = polymer_info[poly_type];

			//
			// Add required bond and angle types
			//

			//
			// Bond type grafting this polymer onto the sphere.
			//
			int graft_bond_type = topo_info.NextAvailableBondType();
			{
				sprintf( buffer, "Polymer type %d connection to sphere", (int)poly_type );
				topo_info.bonds[ graft_bond_type ] = { graft_bond_type, p.dr, buffer };
			}

			//
			// Figure out what the starting bond/angle types should be in the polymer
			// on the basis of what we've already defined in the system.
			//
			int bond_type_offset  = topo_info.NextAvailableBondType();
			int angle_type_offset = topo_info.NextAvailableAngleType();

			//
			// Bond types along the polymer
			//
			for( size_t ii=0; ii<p.atom_types.size()-1; ii++ )
			{
				int type = bond_type_offset + (int)ii;

				if( topo_info.bonds.find(type) == topo_info.bonds.end() )
				{
					sprintf( buffer, "Polymer type %d, bond %d", (int)poly_type, (int)ii );
					topo_info.bonds[type] = { type, p.dr, buffer };
				}
			}

			//
			// Angle types along the polymer
			//
			for( size_t ii=0; ii<p.atom_types.size()-2; ii++ )
			{
				int type = angle_type_offset + (int)ii;
				
				if( topo_info.angles.find(type) == topo_info.angles.end() )
				{
					sprintf( buffer, "Polymer type %d, angle %d", (int)poly_type, (int)ii );
					topo_info.angles[type] = { type, M_PI, buffer };
				}
			}

			//
			// Generate polymer instances
			//

			for( int poly_inst=0; poly_inst<p.N; poly_inst++ )
			{
				int sphere_idx = cp.pairs[graft_index].j; // site onto which we're grafting the polymer

				//
				// Add bond connecting polymer to sphere. VMD is stupid, so we
				// ensure lowest index listed first (i.e. sphere grat point)
				//
				{
					mol.AddBond( graft_bond_type, sphere_idx, mol.atoms.size() );
				}

				//
				// Generate polymer instance
				//
				double sx = hi_res_xyz[sphere_idx*3 +0];
				double sy = hi_res_xyz[sphere_idx*3 +1];
				double sz = hi_res_xyz[sphere_idx*3 +2];

				double abs = sqrt( sx*sx + sy*sy + sz*sz );
				double dx = sx / abs;
				double dy = sy / abs;
				double dz = sz / abs;

				GeneratePolymer( sx,sy,sz, dx,dy,dz, p, mol, bond_type_offset, angle_type_offset );

				graft_index++;
			}
		}
		printf( "done.\n" );

		printf( "\n" );
		printf( "*\n" );
		printf( "* Topological information ...\n" );
		printf( "*\n" );
		printf( "\n" );

		topo_info.Print();

		//
		// Final fix-up of LAMMPS info: expand the bounding box a little, and make sure
		// we have particle masses defined (use dummy masses for missing types).
		//
		{
			int max_type = 0;
			double default_mass = 10.0;
			double margin = 10.0;
			std::map< int, double > mass_map;

			config.RecalculateBounds();
			
			config.bounds.minx -= margin;
			config.bounds.maxx += margin;

			config.bounds.miny -= margin;
			config.bounds.maxy += margin;

			config.bounds.minz -= margin;
			config.bounds.maxz += margin;

			//
			// Which particle types are present?
			//
			for( const auto& mol : config.molecules )
			{
				for( const auto& a : mol.atoms )
				{
					int type = a.type+1; // zero -> unit based
					mass_map[ type ] = default_mass;
					if( type > max_type ) max_type = type;
				}
			}

			//
			// Add masses for present types, including dummy masses for omitted types.
			//
			for( int type=1; type<=max_type; type++ )
			{
				const auto& it = mass_map.find(type);
				double mass = (it!=mass_map.end()) ? (it->second) : (-1.0);
				config.mass.push_back( mass );
			}



		}

		//
		// Save system as both PDB and LAMMPS config.
		//
		{
			FILE *f;

			sprintf( buffer, "fuzzball.pdb" );
			f = fopen( buffer, "w" );
			if( f == nullptr )
			{
				printf( "Unable to open file '%s' for writing\n", buffer );
				exit( -1 );
			}
			SaveMolecule( f, mol );
			fclose( f );

			sprintf( buffer, "fuzzball.lammps_config" );
			f = fopen( buffer, "w" );
			if( f == nullptr )
			{
				printf( "Unable to open file '%s' for writing\n", buffer );
				exit( -1 );
			}
			LAMMPS::SaveData( f, config );
			fclose( f );
		}

		//
		// Save system as UCG topology description
		//
		{
			FILE *f;

			sprintf( buffer, "fuzzball.ucg_topo" );
			f = fopen( buffer, "w" );
			if( f == nullptr )
			{
				printf( "Unable to open file '%s' for writing\n", buffer );
				exit( -1 );
			}

			//
			// Set some variables.
			//
			{
				std::map< int, int > atm_types, bnd_types, ang_types;

				for( const auto& m : config.molecules )
				{
					for( const auto& x : m.atoms )  atm_types[ x.type ] = 1;
					for( const auto& x : m.bonds )  bnd_types[ x.type ] = 1;
					for( const auto& x : m.angles ) ang_types[ x.type ] = 1;
				}

				if( atm_types.size() > 0 )
				{
					fprintf( f, "\n" );
					fprintf( f, "#\n" );
					fprintf( f, "# Autodetected atom information\n" );
					fprintf( f, "#\n" );
					fprintf( f, "\n" );
					for( const auto& a : atm_types ) fprintf( f, "set  FB_ATM_TYPE%d  X\n", a.first );
				}

				if( bnd_types.size() > 0 )
				{
					fprintf( f, "\n" );
					fprintf( f, "#\n" );
					fprintf( f, "# Autodetected bond information\n" );
					fprintf( f, "#\n" );
					fprintf( f, "\n" );
					for( const auto& b : bnd_types ) fprintf( f, "set  FB_BND_TYPE%d  X\n", b.first );
					fprintf( f, "\n" );
					for( const auto& b : bnd_types ) fprintf( f, "set  FB_BND_K%d  X\n", b.first );
					fprintf( f, "\n" );
					for( const auto& b : bnd_types ) fprintf( f, "set  FB_BND_R%d  %g\n", b.first, topo_info.bonds[b.first].rest );
				}

				if( ang_types.size() > 0 )
				{
					fprintf( f, "\n" );
					fprintf( f, "#\n" );
					fprintf( f, "# Autodetected angle information\n" );
					fprintf( f, "#\n" );
					fprintf( f, "\n" );
					for( const auto& a : ang_types ) fprintf( f, "set  FB_ANG_TYPE%d  X\n", a.first );
					fprintf( f, "\n" );
					for( const auto& a : ang_types ) fprintf( f, "set  FB_ANG_K%d  X\n", a.first );
				}

			}

			fprintf( f, "\n" );
			fprintf( f, "subunit X Fuzzball\n" );

			//
			// Members.
			//
			{
				int serial = 0;

				fprintf( f, "\n" );
				fprintf( f, "\t #\n" );
				fprintf( f, "\t # Members\n" );
				fprintf( f, "\t #\n" );
				fprintf( f, "\n" );

				for( const auto& m : config.molecules )
				{
					for( const auto& a : m.atoms )
					{
						fprintf( f, "\t member ${FB_ATM_TYPE%d}  m%d\n", a.type, serial );
						serial++;
					}
				}
			}

			//
			// Bonds.
			//
			{
				for( const auto& m : config.molecules )
				{
					if( m.bonds.size() > 0 )
					{
						fprintf( f, "\n" );
						fprintf( f, "\t #\n" );
						fprintf( f, "\t # Bonds\n" );
						fprintf( f, "\t #\n" );
						fprintf( f, "\n" );
						for( const auto& b : m.bonds )
						{
							fprintf( f, "\t topology  ${FB_BND_TYPE%d}  m%d  m%d  parameters  ${FB_BND_K%d}  ${FB_BND_R%d}\n", b.type, b.i, b.j, b.type, b.type );
						}
					}

					if( m.angles.size() > 0 )
					{
						fprintf( f, "\n" );
						fprintf( f, "\t #\n" );
						fprintf( f, "\t # Angles\n" );
						fprintf( f, "\t #\n" );
						fprintf( f, "\n" );
						for( const auto& a : m.angles )
						{
							fprintf( f, "\t topology  ${FB_ANG_TYPE%d}  m%d  m%d  m%d  parameters  ${FB_ANG_K%d}\n", a.type, a.i, a.j, a.k, a.type );
						}
					}
				}
			}

			fprintf( f, "end\n" );

			fclose( f );
		}
	}

}
