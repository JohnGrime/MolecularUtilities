/*
	Author: John Grime, The University of Chicago.
*/

#include "stdio.h"

#include "../Util/Util.h"

using namespace Util;

using v3 = Geometry::Vec3<double>;

/*
	Generate a single-component bilayer. Consider using the GenerateMembranes program instead,
	as it's more flexible and powerful.
*/

//
// Surfaces defined by paired position and normal direction vector entries.
//
struct Surface
{
	std::vector< v3 > pos, dir;

	//
	// Origin is centre, normals "up" on z: transform later if needed.
	//
	static int ZPlane(
		Surface& s,
		int Nx, double minx, double maxx,
		int Ny, double miny, double maxy )
	{
		s.pos.resize( Nx*Ny );
		s.dir.resize( Nx*Ny );

		double dx = (maxx-minx)/Nx;
		double dy = (maxy-miny)/Ny;

		int i = 0;
		for( int x_=0; x_<Nx; x_++ )
		{
			for( int y_=0; y_<Ny; y_++ )
			{
				s.pos[i] = { minx+(0.5+x_)*dx, miny+(0.5+y_)*dy, 0.0 };
				s.dir[i] = { 0.0, 0.0, 1.0 };
				i++;
			}
		}
		return (int)s.pos.size();
	}

	//
	// Origin is centre, normals "out": transform later if needed.
	// Points on a unit sphere == normal vector at that location.
	//
	static int Sphere(
		Surface& s,
		int target_N, double radius )
	{
		Util::Geometry::GetEquispacedSpherePoints( target_N, s.dir, 1.0 );
		s.pos = s.dir;
		for( auto& p : s.pos ) p = p*radius;
		return (int)s.pos.size();
	}
};

//
// Some IO routines for PDB files.
//
struct PDB
{
	static constexpr const char* atm_fmt = "%-6.6s%5d %-4s%1s%3s %1s%4d%1s   %8.2f%8.2f%8.2f%6s%6s%6s%4s%2s%2s";
	static constexpr const char* bnd_fmt = "CONECT%5d%5d";

	static void WriteAtom( FILE *f, int serial, const char* name, int resSeq, const char* resName, const char* chainID, double x, double y, double z )
	{
		assert( f != nullptr );
		fprintf( f, atm_fmt, "ATOM",
			(int)serial, name,
			"",
			resName, chainID, resSeq,
			"",
			x, y, z,
			"", "", "", "", "", "" );
	}

	static void WriteBond( FILE *f, int i, int j )
	{
		assert( f != nullptr );
		fprintf( f, bnd_fmt, i, j );
	}

	static void WriteSurface( const Surface& s, FILE* f, double scale = 1.0 )
	{
		assert( f != nullptr );

		int serial = 0;
		for( size_t i=0, N=s.pos.size(); i<N; i++ )
		{
			const auto& r = s.pos[i];
			auto v = r + s.dir[i]*scale;

			serial++;
			WriteAtom( f, serial, "C", i, "SRF", "A", r.x, r.y, r.z );
			fprintf( f, "\n" );

			serial++;
			WriteAtom( f, serial, "H", i, "SRF", "A", v.x, v.y, v.z );
			fprintf( f, "\n" );
		}

		fprintf( f, "TER\n" );

		for( size_t i=0, N=s.pos.size(); i<N; i++ )
		{
			WriteBond( f, (int)(i*2+0)+1, (int)(i*2+1)+1 );
			fprintf( f, "\n" );
		}
	}

	static void WriteMolecules(
		FILE *f,
		const std::vector< Surface >& surfaces,
		const std::vector< int >& molecule_types,
		const std::vector< LAMMPS::Molecule >& molecules )
	{
		assert( surfaces.size() <= molecule_types.size() );

		std::vector<double> mol_axes = { 0,0,0, 1,0,0, 0,0,1 }; // don't need y: origin counts as one of the three points needed for 3D superposition
		std::vector<double> srf_axes = { 0,0,0, 1,0,0, 0,0,1 };
		std::vector<double> src_xyz, xyz;
		double M[3][3], T[3];

		char name[10], resName[10];

		//
		// Write atoms
		//
		int serial = 0;
		for( size_t si=0, sN=surfaces.size(); si<sN; si++ )
		{
			const auto surf = surfaces[si];
			const auto mol = molecules[ molecule_types[si] ];

			//
			// Source vertices that describe the molecule geometry.
			//
			src_xyz.clear();
			for( const auto& a : mol.atoms )
			{
				src_xyz.push_back( a.x );
				src_xyz.push_back( a.y );
				src_xyz.push_back( a.z );
			}
			xyz.resize( src_xyz.size() );

			sprintf( resName, "S%d", (int)si );

			for( size_t pi=0, pN=surf.pos.size(); pi<pN; pi++ )
			{
				//
				// Get local reference frame for direction vector and
				// translate origin onto position.
				//
				v3 z_ = surf.dir[pi];
				v3 x_ = Geometry::Orthogonal( z_ );

				z_ = z_ / Geometry::Abs( z_ );
				x_ = x_ / Geometry::Abs( x_ );

				srf_axes[ 0*3 +0 ] = surf.pos[pi].x;
				srf_axes[ 0*3 +1 ] = surf.pos[pi].y;
				srf_axes[ 0*3 +2 ] = surf.pos[pi].z;

				srf_axes[ 1*3 +0 ] = surf.pos[pi].x + x_.x;
				srf_axes[ 1*3 +1 ] = surf.pos[pi].y + x_.y;
				srf_axes[ 1*3 +2 ] = surf.pos[pi].z + x_.z;

				srf_axes[ 2*3 +0 ] = surf.pos[pi].x + z_.x;
				srf_axes[ 2*3 +1 ] = surf.pos[pi].y + z_.y;
				srf_axes[ 2*3 +2 ] = surf.pos[pi].z + z_.z;

				//
				// Superpose molecule onto reference frame
				//
				Superposer::Calculate( srf_axes, mol_axes, M, T );
				Superposer::Apply( src_xyz, xyz, M, T );

				//
				// Print molecule
				//
				for( size_t ai=0, aN = mol.atoms.size(); ai<aN; ai++ )
				{
					sprintf( name, "%d", mol.atoms[ai].type );
					WriteAtom( f, serial+1, name, (int)ai+1, resName, "A", xyz[ai*3+0], xyz[ai*3+1], xyz[ai*3+2] );
					fprintf( f, "\n" );
					serial++;
				}
				fprintf( f, "TER\n" );
			}
		}

		//
		// Write bonds
		//
		serial = 0;
		for( size_t si=0, sN=surfaces.size(); si<sN; si++ )
		{
			const auto surf = surfaces[si];
			const auto mol = molecules[ molecule_types[si] ];

			for( size_t pi=0, pN=surf.pos.size(); pi<pN; pi++ )
			{
				for( size_t bi=0, bN = mol.bonds.size(); bi<bN; bi++ )
				{
					auto i = serial + mol.bonds[bi].i;
					auto j = serial + mol.bonds[bi].j;
					WriteBond( f, i, j );
					fprintf( f, "\n" );
				}
				serial += mol.atoms.size();
			}
		}
	}

};


//
// Given a set of molecules and surfaces, write the resultant LAMMPS data.
//
void WriteLAMMPSMolecules(
	FILE *f,
	const std::vector< Surface >& surfaces,
	const std::vector< int >& molecule_types,
	const std::vector< LAMMPS::Molecule >& molecules,
	double minx, double maxx,
	double miny, double maxy,
	double minz, double maxz )
{
	assert( surfaces.size() <= molecule_types.size() );

	std::vector<double> mol_axes = { 0,0,0, 1,0,0, 0,0,1 }; // don't need y: origin counts as one of the three points needed for 3D superposition.
	std::vector<double> srf_axes = { 0,0,0, 1,0,0, 0,0,1 };
	std::vector<double> src_xyz, xyz;
	double M[3][3], T[3];

	int n_atm = 0, n_bnd = 0, n_ang = 0;
	std::map<int,int> atm_types, bnd_types, ang_types;

	//
	// Figure out number of atom/bond/angle types
	//
	{
		for( const auto& m: molecules )
		{
			for( const auto &x : m.atoms )
			{
				auto it = atm_types.find(x.type);
				if( it == end(atm_types) ) atm_types[x.type] = 1;
				else it->second++;
			}

			for( const auto &x : m.bonds )
			{
				auto it = bnd_types.find(x.type);
				if( it == end(bnd_types) ) bnd_types[x.type] = 1;
				else it->second++;
			}

			for( const auto &x : m.angles )
			{
				auto it = ang_types.find(x.type);
				if( it == end(ang_types) ) ang_types[x.type] = 1;
				else it->second++;
			}
		}

		for( size_t si=0, sN=surfaces.size(); si<sN; si++ )
		{
			auto t = molecule_types[si];
			const auto& mol = molecules[t];
			const auto N_mol = surfaces[si].pos.size();
			n_atm += (int)mol.atoms.size()  * N_mol;
			n_bnd += (int)mol.bonds.size()  * N_mol;
			n_ang += (int)mol.angles.size() * N_mol;
		}

		fprintf( f, "Autogenerated\n" );

		fprintf( f, "\n" );
		fprintf( f, "%d atoms\n", n_atm );
		fprintf( f, "%d atom types\n", (int)atm_types.size() );
		fprintf( f, "%d bonds\n", n_bnd );
		fprintf( f, "%d bond types\n", (int)bnd_types.size() );
		fprintf( f, "%d angles\n", n_ang );
		fprintf( f, "%d angle types\n", (int)ang_types.size() );

		fprintf( f, "\n" );
		fprintf( f, "%g %g xlo xhi\n", minx, maxx );
		fprintf( f, "%g %g ylo yhi\n", miny, maxy );
		fprintf( f, "%g %g zlo zhi\n", minz, maxz );

		fprintf( f, "\n" );
		fprintf( f, "Masses\n" );
		fprintf( f, "\n" );
		for( const auto& it : atm_types ) fprintf( f, "%d %g\n", it.first, 152.0 );
	}

	//
	// Print atoms
	//
	if( n_atm > 0 )
	{
		fprintf( f, "\n" );
		fprintf( f, "Atoms\n" );
		fprintf( f, "\n" );

		int serial = 0, mol_id = 0;
		for( size_t si=0, sN=surfaces.size(); si<sN; si++ )
		{
			const auto surf = surfaces[si];
			const auto mol = molecules[ molecule_types[si] ];

			src_xyz.clear();
			for( const auto& a : mol.atoms )
			{
				src_xyz.push_back( a.x );
				src_xyz.push_back( a.y );
				src_xyz.push_back( a.z );
			}
			xyz.resize( src_xyz.size() );

			for( size_t pi=0, pN=surf.pos.size(); pi<pN; pi++ )
			{
				//
				// Get local reference frame for direction vector and
				// translate origin onto position.
				//
				v3 z_ = surf.dir[pi];
				v3 x_ = Geometry::Orthogonal( z_ );

				srf_axes[ 0*3 +0 ] = surf.pos[pi].x;
				srf_axes[ 0*3 +1 ] = surf.pos[pi].y;
				srf_axes[ 0*3 +2 ] = surf.pos[pi].z;

				srf_axes[ 1*3 +0 ] = surf.pos[pi].x + x_.x;
				srf_axes[ 1*3 +1 ] = surf.pos[pi].y + x_.y;
				srf_axes[ 1*3 +2 ] = surf.pos[pi].z + x_.z;

				srf_axes[ 2*3 +0 ] = surf.pos[pi].x + z_.x;
				srf_axes[ 2*3 +1 ] = surf.pos[pi].y + z_.y;
				srf_axes[ 2*3 +2 ] = surf.pos[pi].z + z_.z;

				//
				// Superpose molecule onto reference frame
				//
				Superposer::Calculate( srf_axes, mol_axes, M, T );
				Superposer::Apply( src_xyz, xyz, M, T );

				//
				// Print molecule
				//
				for( size_t ai=0, aN = mol.atoms.size(); ai<aN; ai++ )
				{
					fprintf( f, "%d %d %d %.3f %.3f %.3f\n", serial+1, mol_id+1, mol.atoms[ai].type, xyz[ai*3+0], xyz[ai*3+1], xyz[ai*3+2] );
					serial++;
				}

				mol_id++;
			}
		}
	}

	//
	// Print bonds
	//
	if( n_bnd > 0 )
	{
		int serial = 0, offset = 0;

		fprintf( f, "\n" );
		fprintf( f, "Bonds\n" );
		fprintf( f, "\n" );

		for( size_t si=0, sN=surfaces.size(); si<sN; si++ )
		{
			const auto surf = surfaces[si];
			const auto mol = molecules[ molecule_types[si] ];

			for( size_t pi=0, pN=surf.pos.size(); pi<pN; pi++ )
			{
				for( size_t bi=0, bN = mol.bonds.size(); bi<bN; bi++ )
				{
					auto t = mol.bonds[bi].type;
					auto i = offset + mol.bonds[bi].i;
					auto j = offset + mol.bonds[bi].j;
					fprintf( f, "%d %d %d %d\n", serial+1, t, i, j );
					serial++;
				}
				offset += (int)mol.atoms.size();
			}
		}
	}

	//
	// Print angles
	//
	if( n_ang > 0 )
	{
		int serial = 0, offset = 0;

		fprintf( f, "\n" );
		fprintf( f, "Angles\n" );
		fprintf( f, "\n" );

		for( size_t si=0, sN=surfaces.size(); si<sN; si++ )
		{
			const auto surf = surfaces[si];
			const auto mol = molecules[ molecule_types[si] ];

			for( size_t pi=0, pN=surf.pos.size(); pi<pN; pi++ )
			{
				for( size_t ai=0, aN = mol.angles.size(); ai<aN; ai++ )
				{
					auto t = mol.angles[ai].type;
					auto i = offset + mol.angles[ai].i;
					auto j = offset + mol.angles[ai].j;
					auto k = offset + mol.angles[ai].k;
					fprintf( f, "%d %d %d %d %d\n", serial+1, t, i, j, k );
					serial++;
				}
				offset += (int)mol.atoms.size();
			}
		}
	}
}

//
// Volume enclosed between surfaces of two spheres
//
double GapVolume( double r0, double r1 )
{
	constexpr double pre = 4.0/3 * M_PI;
	double r0_3 = r0*r0*r0;
	double r1_3 = r1*r1*r1;
	return pre * ( r1_3 - r0_3 );
}

//
// Parse lipid description; assumes single-tailed for now.
// Format : t1,t2,...[:b1,b2,...][:a1,a2,...]
// Above specifies atom, bond, and angle types.
//
void ParseLipidDefinition( const char* desc, LAMMPS::Molecule& mol )
{
	int type;
	double q=0, x=0, y=0, z=0;

	std::vector<std::string> toks, subtoks;

	mol.atoms.clear();
	mol.bonds.clear();
	mol.angles.clear();

	String::Tokenize( desc, toks, ":" );

	//
	// Atom types
	//
	{
		String::Tokenize( toks[0], subtoks, "," );
		for( const auto& t : subtoks )
		{
			if( String::ToInteger( t, type ) != String::ReturnValue::OK )
			{
				printf( "Unable to convert atom type token '%s'\n", t.c_str() );
				exit( -1 );

			}
			mol.AddAtom( mol.atoms.size()+1, type, q, x,y,z );
		}
	}

	//
	// Bonds
	//
	if( mol.atoms.size() > 1 )
	{
		if( toks.size() > 1 )
		{
			String::Tokenize( toks[1], subtoks, "," );
			for( size_t i=0; i<subtoks.size(); i++ )
			{
				if( String::ToInteger( subtoks[i], type ) != String::ReturnValue::OK )
				{
					printf( "Unable to convert bond type token '%s'\n", subtoks[i].c_str() );
					exit( -1 );

				}
				mol.AddBond( type, i+1, i+2 );
			}		
		}
		else
		{
			for( size_t i=0; i<mol.atoms.size()-1; i++ ) mol.AddBond( 1, i+1, i+2 );
		}
	}

	//
	// Angles
	//	
	if( mol.atoms.size() > 2 )
	{
		if( toks.size() > 2 )
		{
			String::Tokenize( toks[2], subtoks, "," );
			for( size_t i=0; i<subtoks.size(); i++ )
			{
				if( String::ToInteger( subtoks[i], type ) != String::ReturnValue::OK )
				{
					printf( "Unable to convert angle type token '%s'\n", subtoks[i].c_str() );
					exit( -1 );

				}
				mol.AddAngle( type, i+1, i+2, i+3 );
			}		
		}
		else
		{
			for( size_t i=0; i<mol.atoms.size()-2; i++ ) mol.AddAngle( 1, i+1, i+2, i+3 );
		}
	}

	//
	// Sanity check
	//
	if( mol.bonds.size() > mol.atoms.size()-1 )
	{
		printf( "Warning: too many bonds defined, truncating.\n" );
		mol.bonds.resize( mol.atoms.size()-1 );
	}
	if( mol.angles.size() > mol.atoms.size()-2 )
	{
		printf( "Warning: too many angles defined, truncating.\n" );
		mol.angles.resize( mol.atoms.size()-2 );
	}

	printf( "Molecule:\n" );
	
	printf( "  atom types: " );
	for( const auto& a : mol.atoms ) printf( "%d ", a.type );
	printf( "\n" );

	if( mol.bonds.size() > 0 )
	{
		printf( "  bonds: " );
		for( const auto& b : mol.bonds ) printf( "[%d,%d,%d] ", b.type, b.i, b.j );
		printf( "\n" );
	}

	if( mol.angles.size() > 0 )
	{
		printf( "  angles: " );
		for( const auto& a : mol.angles ) printf( "[%d,%d,%d,%d] ", a.type, a.i, a.j, a.k );
		printf( "\n" );
	}
}


int main( int argc, char** argv )
{
	std::vector< Surface > surfaces;
	std::vector<int> molecule_types;
	std::vector< LAMMPS::Molecule > molecules;

	double bond_length = 7.5;
	const char* lipid_definition = "1,2,3,3,3"; // default: BPB-style 5-site lipid model.

	if( argc < 3 )
	{
		printf( "Usage: %s <outer radius> <APL> <bond_length> [t1,t1,...[:b1,b2,...][:a1,a2,...]]\n", argv[0] );
		exit( -1 );
	}
	double r, APL;

	if( String::ToReal(argv[1],r) != String::ReturnValue::OK )
	{
		printf( "Can't convert '%s' into radius\n", argv[1] );
		exit( -1 );
	}

	if( String::ToReal(argv[2],APL) != String::ReturnValue::OK )
	{
		printf( "Can't convert '%s' into APL\n", argv[2] );
		exit( -1 );
	}

	if( String::ToReal(argv[3],bond_length) != String::ReturnValue::OK )
	{
		printf( "Can't convert '%s' into bond_length\n", argv[2] );
		exit( -1 );
	}

	if( argc > 4 ) lipid_definition = argv[4];

	//
	// Generate template lipid molecule
	//
	{
		LAMMPS::Molecule mol;

		ParseLipidDefinition( lipid_definition, mol );

		for( size_t i=0; i<mol.atoms.size(); i++ )
		{
			auto& a = mol.atoms[i];
			a.q = 0.0;
			a.x = 0.0;
			a.y = 0.0;
			a.z = bond_length*i;
			printf( "%g,%g,%g\n", a.x, a.y, a.z );
		}

		molecules.push_back( mol );
	}

	double lipid_length = bond_length * (molecules[0].atoms.size()-1);
	if( lipid_length < bond_length ) lipid_length = bond_length;

	double outer_r = r;
	double inner_r = outer_r - (2.0*lipid_length);

	double outer_A = 4.0 * M_PI * (outer_r*outer_r);
	double inner_A = 4.0 * M_PI * (inner_r*inner_r);

	double outer_V = GapVolume( outer_r-lipid_length, outer_r );
	double inner_V = GapVolume( inner_r, inner_r+lipid_length );

	int outer_N = (int) ( outer_A / APL );
	int inner_N = (int) ( inner_A / APL );

	bool use_volume = true;
//	bool use_volume = false;
	if( use_volume == true )
	{

		//
		// Fixing lipid count via area of outer/inner radii
		// and area per lipid is a bad idea: breaks totally
		// for small vesicles.
		//
		//	Instead, preserve DENSITY, rho, in both monolayers?
		//	Assuming dz is the width of a single monolayer:
		//
		//	Volume per lipid, VPL = APL*Lz
		//
		//	Volume enclosed by spheres of radius r0 and r1 :
		//		= 4/3.pi.(r1^3) - 4/3.pi.(r0^3)
		//		= 4.3.pi.( r1^3 - r0^3 )
		//
		//	Therefore number of lipids in that volume should be:
		//		4.3.pi.( r1^3 - r0^3 ) / (APL.Lz)
		//
		outer_N = (int) ( outer_V / (APL*lipid_length) );
		inner_N = (int) ( inner_V / (APL*lipid_length) );
	}

	printf( "Outer leaflet : radius %g, %d lipids, VPL = %g, APL %g\n", outer_r, outer_N, outer_V/outer_N, outer_A/outer_N );
	printf( "Inner leaflet : radius %g, %d lipids, VPL = %g, APL %g\n", inner_r, inner_N, inner_V/inner_N, inner_A/inner_N );

	//
	// Generate monolayers
	//
	{
		Surface s;

		double leaflet_sep = bond_length; // safety gap between terminal tail beads?

		Surface::Sphere( s, outer_N, outer_r+leaflet_sep );
		for( auto& d : s.dir ) d = d * -1;
		surfaces.push_back( s );
		molecule_types.push_back( 0 );

		Surface::Sphere( s, inner_N, inner_r );
		for( auto& d : s.dir ) d = d * +1;
		surfaces.push_back( s );
		molecule_types.push_back( 0 );
	}

	//
	// Save individual surfaces as PDB
	//
	{
		char path[128];
		for( size_t i=0; i<surfaces.size(); i++ )
		{
			sprintf( path, "surface.%d.pdb", (int)i );
			FILE* f = fopen( path, "w" );
			PDB::WriteSurface( surfaces[i], f, 50.0 );
			fclose( f );
		}
	}


	//
	// Expand surface definition using molecules
	//
	{
		FILE* f = fopen( "molecules.pdb", "w" );
		PDB::WriteMolecules( f, surfaces, molecule_types, molecules );
		fclose( f );
	}

	//
	// Expand surface definition using molecules
	//
	{
		double margin = 100.0;

		double maxx = outer_r + margin;
		double minx = -maxx;

		double maxy = outer_r + margin;
		double miny = -maxy;

		double maxz = outer_r + margin;
		double minz = -maxz;

		FILE* f = fopen( "molecules.lammps_config", "w" );
		WriteLAMMPSMolecules( f, surfaces, molecule_types, molecules, minx,maxx, miny,maxy, minz,maxz );
		fclose( f );
	}

}
