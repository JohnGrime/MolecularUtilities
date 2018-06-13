/*
	Author: John Grime, The University of Chicago.
*/

#include "stdio.h"
#include <random>

#include "../Util/Util.h"

using namespace Util;

using v3 = Geometry::Vec3<double>;

//
// Surfaces defined by paired position and normal direction vector entries.
// We also store molecule types, so we don't worry about slicing.
//
struct MolecularSurface
{
	std::vector< v3 > pos, dir;
	std::vector<int> mol_types;

	void Clear()
	{
		pos.clear();
		dir.clear();
		mol_types.clear();
	}

	//
	// Generate by hand, rather than std::discrete_distribution from <random> in c++11
	//
	void RandomAssign( const std::vector<int>& types, const std::vector<double>& props )
	{
		std::vector<int> idx, target, actual;
		std::random_device rd;
		std::mt19937 mt( rd() );

		assert( types.size() == props.size() );

		int N_types = (int)types.size();
		int N_mols = (int)mol_types.size();

		double S = 0.0;
		for( auto& x : props ) S += x;
   	
		idx.resize( N_mols );
		std::iota( idx.begin(), idx.end(), 0 );
		std::shuffle( idx.begin(), idx.end(), mt );

		//
		// Assign mol_types
		//
		int upto = 0;
		target.resize( N_types );
		actual.assign( N_types, 0 );
		for( int ti=0; (ti<N_types) && (upto<N_mols); ti++ )
		{
			target[ti] = (int)ceil( (props[ti]/S) * N_mols );

			for( int i=0; (i<target[ti]) && (upto<N_mols); i++ )
			{
				mol_types[ idx[upto] ] = types[ti];
				actual[ti]++;
				upto++;
			}
		}

		/*
		int test = 0;
		for( int ti=0; ti<N_types; ti++ )
		{
			printf( "%d %d\n", target[ti], actual[ti] );
			test += actual[ti];
		}
		printf( "%d\n", test );
		*/
	}

	MolecularSurface() { Clear(); }
};

struct MolecularPlane : public MolecularSurface
{
	enum class Axis { X, Y, Z };

	int Generate( int mol_type, Axis axis, int N1, int N2, double min1, double max1, double min2, double max2 )
	{
		mol_types.resize( N1*N2 );
		pos.resize( N1*N2 );
		dir.resize( N1*N2 );

		double d1 = (max1-min1)/N1;
		double d2 = (max2-min2)/N2;

		int i = 0;
		for( int a_=0; a_<N1; a_++ )
		{
			for( int b_=0; b_<N2; b_++ )
			{
				double c1 = min1+(0.5+a_)*d1; 
				double c2 = min2+(0.5+b_)*d2; 

				switch( axis )
				{
					case Axis::X:
						pos[i] = { 0.0,  c1,  c2 };
						dir[i] = { 1.0, 0.0, 0.0 };
					break;

					case Axis::Y:
						pos[i] = {  c1, 0.0,  c2 };
						dir[i] = { 0.0, 1.0, 0.0 };
					break;

					case Axis::Z:
						pos[i] = {  c1,  c2, 0.0 };
						dir[i] = { 0.0, 0.0, 1.0 };
					break;
				}

				mol_types[i] = mol_type;

				i++;
			}
		}
		return (int)pos.size();
	};

	MolecularPlane( int mol_type, Axis axis, int N1, int N2, double min1, double max1, double min2, double max2 )
	{
		Generate( mol_type, axis, N1,N2, min1,max1, min2,max2 );
	}
	MolecularPlane()
	{
		Generate( 0, Axis::Z, 10,10, -5,+5, -5,+5 );
	}
};

struct MolecularSphere : public MolecularSurface
{
	int Generate( int mol_type, int target_N, double radius, const v3& origin )
	{
		Util::Geometry::GetEquispacedSpherePoints( target_N, dir, 1.0 );

		pos = dir;
		for( auto& p : pos ) p = (p*radius) + origin;

		mol_types.resize( pos.size() );
		for( auto& t : mol_types ) t = mol_type;

		return (int)pos.size();
	}

	MolecularSphere( int mol_type, int target_N, double radius, const v3& origin )
	{
		Generate( mol_type, target_N, radius, origin );
	}
	MolecularSphere()
	{
		Generate( 0, 10, 1.0, {0,0,0} );
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

	//
	// Write positions / normals
	//
	static void WriteSurface( const MolecularSurface& s, FILE* f, double scale = 1.0 )
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

	//
	// Write molecules
	//
	static void WriteMolecularSurface(
		FILE *f,
		const std::vector< MolecularSurface >& surfaces,
		const std::vector< LAMMPS::Molecule >& molecules )
	{
		static const char* chainIDs = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
		static const int N_chainIDs = strlen( chainIDs );

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
			const auto& surf = surfaces[si];
			const char chainID[] = { chainIDs[ si%N_chainIDs ], '\0' };

			for( size_t pi=0, pN=surf.pos.size(); pi<pN; pi++ )
			{
				auto mol_type = surf.mol_types[pi];
				assert( mol_type < (int)molecules.size() );
				const auto& mol = molecules[ mol_type ];

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

				sprintf( resName, "%d", mol_type );

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
					WriteAtom( f, serial+1, name, (int)ai+1, resName, chainID, xyz[ai*3+0], xyz[ai*3+1], xyz[ai*3+2] );
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
			const auto& surf = surfaces[si];

			for( size_t pi=0, pN=surf.pos.size(); pi<pN; pi++ )
			{
				auto mol_type = surf.mol_types[pi];
				assert( mol_type < (int)molecules.size() );
				const auto& mol = molecules[mol_type];

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
	const std::vector< MolecularSurface >& surfaces,
	const std::vector< LAMMPS::Molecule >& molecules,
	double minx, double maxx,
	double miny, double maxy,
	double minz, double maxz )
{
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
			const auto& surf = surfaces[si];
			for( auto t : surf.mol_types )
			{
				assert( t < (int)molecules.size() );
				const auto& mol = molecules[t];
				n_atm += (int)mol.atoms.size();
				n_bnd += (int)mol.bonds.size();
				n_ang += (int)mol.angles.size();
			}
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

			for( size_t pi=0, pN=surf.pos.size(); pi<pN; pi++ )
			{
				auto mol_type = surf.mol_types[pi];
				assert( mol_type < (int)molecules.size() );
				const auto mol = molecules[ mol_type ];

				src_xyz.clear();
				for( const auto& a : mol.atoms )
				{
					src_xyz.push_back( a.x );
					src_xyz.push_back( a.y );
					src_xyz.push_back( a.z );
				}
				xyz.resize( src_xyz.size() );

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

			for( size_t pi=0, pN=surf.pos.size(); pi<pN; pi++ )
			{
				auto mol_type = surf.mol_types[pi];
				assert( mol_type < (int)molecules.size() );
				const auto mol = molecules[ mol_type ];

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

			for( size_t pi=0, pN=surf.pos.size(); pi<pN; pi++ )
			{
				auto mol_type = surf.mol_types[pi];
				assert( mol_type < (int)molecules.size() );
				const auto mol = molecules[ mol_type ];

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
}

//
// type1,...[?prop1,...]
//
void ParseMOLSTRING( const char* molstring, std::vector<int>& mol_types, std::vector<double>& mol_props )
{
	std::vector<std::string> toks, types, props;

	String::Tokenize( molstring, toks, "?" );
	String::Tokenize( toks[0], types, "," );

	mol_types.clear();
	mol_props.clear();

	for( const auto& x : types )
	{
		int t;
		if( (String::ToInteger(x,t)!=String::ReturnValue::OK) || (t<1) )
		{
			printf( "Bad molecule type '%s' in '%s'\n", x.c_str(), molstring );
			exit( -1 );
		}
		mol_types.push_back( t-1 );
	}
	
	if( toks.size() > 1 )
	{
		String::Tokenize( toks[1], props, "," );
		if( props.size() != types.size() )
		{
			printf( "Too few molecular proportions defined: '%s' vs '%s' in '%s'\n", toks[0].c_str(), toks[1].c_str(), molstring );
			exit( -1 );
		}

		for( const auto& x : props )
		{
			double p;
			if( String::ToReal(x,p) != String::ReturnValue::OK )
			{
				printf( "Bad molecule proportion '%s' in '%s'\n", x.c_str(), molstring );
				exit( -1 );
			}
			mol_props.push_back( p );
		}
	}
	else
	{
		mol_props.assign( mol_types.size(), 1.0 );
	}

	//
	// Normalize the specified molecular proportions.
	//
	{
		double sum = 0.0;
		for( auto& p : mol_props ) sum += p;
		for( auto& p : mol_props ) p /= sum;
	}
}

//
// sphere=monolayer|bilayer:MOLSTRING:outer_r:APL
//
struct SphereDefinition
{
	bool is_monolayer;

	std::vector<int> mol_types;
	std::vector<double> mol_props;

	double outer_r, APL;
};
void ParseSphereDefinition( SphereDefinition& s, const char* def )
{
	std::vector<std::string> toks, subtoks, types, props;

	if( (String::Tokenize(def,toks,":")==-1) || (toks.size()<4) )
	{
		printf( "Bad sphere definition '%s'\n", def );
		exit( -1 );
	}

	if( toks[0] == "monolayer" ) s.is_monolayer = true;
	else if( toks[0] == "bilayer" ) s.is_monolayer = false;
	else
	{
		printf( "Unknown type '%s' in '%s'\n", toks[0].c_str(), def );
		exit( -1 );
	}

	ParseMOLSTRING( toks[1].c_str(), s.mol_types, s.mol_props );

	if( (String::ToReal(toks[2],s.outer_r)!=String::ReturnValue::OK) || (s.outer_r<=0.0) )
	{
		printf( "Bad sphere outer radius '%s' in '%s'\n", toks[2].c_str(), def );
		exit( -1 );
	}

	if( (String::ToReal(toks[3],s.APL)!=String::ReturnValue::OK) || (s.APL<=0.0) )
	{
		printf( "Bad sphere APL '%s' in '%s'\n", toks[3].c_str(), def );
		exit( -1 );
	}
}

//
// plane=monolayer|bilayer:MOLSTRING:N:APL
//
struct PlaneDefinition
{
	bool is_monolayer;

	std::vector<int> mol_types;
	std::vector<double> mol_props;

	int N;
	double APL;
};
void ParsePlaneDefinition( PlaneDefinition& p, const char* def )
{
	std::vector<std::string> toks, subtoks, types, props;

	if( (String::Tokenize(def,toks,":")==-1) || (toks.size()<4) )
	{
		printf( "Bad plane definition '%s'\n", def );
		exit( -1 );
	}

	if( toks[0] == "monolayer" ) p.is_monolayer = true;
	else if( toks[0] == "bilayer" ) p.is_monolayer = false;
	else
	{
		printf( "Unknown type '%s' in '%s'\n", toks[0].c_str(), def );
		exit( -1 );
	}

	ParseMOLSTRING( toks[1].c_str(), p.mol_types, p.mol_props );

	if( (String::ToInteger(toks[2],p.N)!=String::ReturnValue::OK) || (p.N<1) )
	{
		printf( "Bad plane N '%s' in '%s'\n", toks[2].c_str(), def );
		exit( -1 );
	}

	if( (String::ToReal(toks[3],p.APL)!=String::ReturnValue::OK) || (p.APL<=0.0) )
	{
		printf( "Bad plane APL '%s' in '%s'\n", toks[3].c_str(), def );
		exit( -1 );
	}
}



int main( int argc, char** argv )
{
	std::vector<std::string> toks;
	std::vector<SphereDefinition> sphere_defs;
	std::vector<PlaneDefinition> plane_defs;

	std::vector< MolecularSurface > surfaces;
	std::vector< LAMMPS::Molecule > molecules;

	double bond_length = 7.5;
	double leaflet_separation = -1.0;

	if( argc < 2 )
	{
		printf( "\n" );
		printf( "Usage: %s [bond_length=X] [leaflet_separation=X] [lipid=t1,t1,...[:b1,b2,...][:a1,a2,...], ...] [sphere=monolayer|bilayer:COMPOSITION:outer_r:APL, ...] [plane=bilayer|monolayer:COMPOSITION:N:APL, ...]\n", argv[0] );
		printf( "\n" );
		printf( "Where:\n" );
		printf( "\n" );
		printf( "bond_length: OPTIONAL default separation between bound particles (default: 7.5 Angstrom)\n" );
		printf( "leaflet_separation: OPTIONAL separation between inner and outer monolayere leaflets in a bilayer (default: 7.5 Angstrom)\n" );
		printf( "\n" );
		printf( "lipid: lists of particle types defining a lipid molecule, followed by optional lists of bond and angle types (default bond and angle type: 1)\n" );
		printf( "sphere: define a spherical mono- or bilayer of the specified lipid composition, radius, and area per lipid\n" );
		printf( "plane: define a planar mono- or bilayer of the specified lipid composition, NxN lipid grid, and area per lipid\n" );
		printf( "\n" );
		printf( "COMPOSITION : comma-separated lists of lipid types and their relative proportions, separated by question mark.\n" );
		printf( "\n" );
		printf( "Examples:\n" );
		printf( "\n" );
		printf( "1. %s lipid=1,2,3 sphere=bilayer:1:100:70\n", argv[0] );
		printf( "\n" );
		printf( "- Define a 3-site lipid (atom types 1,2,3), referred to in future as lipid type 1.\n" );
		printf( "- Create a spherical bilayer (using lipid type 1) of radius 100 Angstrom and area per lipid 70 Angstrom**2.\n" );
		printf( "\n" );
		printf( "2. %s bond_length=7.5 leaflet_separation=7.5 lipid=1,2,3,3,3 lipid=1,2,3 sphere=bilayer:1,2?1,1:100:70\n", argv[0] );
		printf( "\n" );
		printf( "- Specify a default bond length connecting particles of length 7.5 Angstrom.\n" );
		printf( "- Specify a default separation of 7.5 Angstroms between monolayer leaflets.\n" );
		printf( "- Define TWO lipid types:\n" );
		printf( "    - Lipid type 1 contains 5 particles (particle types: 1,2,3,3,3).\n" );
		printf( "    - Lipid type 2 contains 3 particles (particle types: 1,2,3).\n" );
		printf( "- Create a spherical bilayer with the following properties:\n" );
		printf( "    - Use lipid types 1 and 2 to create a 1:1 ratio of lipids in the bilayer.\n" );
		printf( "    - Bilayer radius is 100 Angstrom, area per lipid is 70 Angstrom**2.\n" );
		printf( "\n" );
		printf( "3. %s lipid=1,2,3,3,3 lipid=1,2,3 sphere=bilayer:1,2?10,12:100:70\n", argv[0] );
		printf( "\n" );
		printf( "- Define TWO lipid types:\n" );
		printf( "    - Lipid type 1 contains 5 particles (particle types: 1,2,3,3,3).\n" );
		printf( "    - Lipid type 2 contains 3 particles (particle types: 1,2,3).\n" );
		printf( "- Create a spherical bilayer with the following properties:\n" );
		printf( "    - Use lipid types 1 and 2 to create a 10:12 ratio of lipids in the bilayer.\n" );
		printf( "    - Bilayer radius is 100 Angstrom, area per lipid is 70 Angstrom**2.\n" );
		printf( "\n" );
		printf( "Notes:\n" );
		printf( "\n" );
		printf( "Lipid types are UNIT BASED and correspond to the order in which lipid definitions occurred.\n" );
		printf( "Molecule proportions are normalised internally, so they don't need to sum to 1 on the command line.\n" );
		printf( "If molecule proportions omitted, equal proportions used.\n" );
		printf( "\n" );
		exit( -1 );
	}

	//
	// sphere=monolayer|bilayer:outer_r:APL
	// plane=monolayer|bilayer:coord:N1:N2:L1:L2
	//

	for( int i=1; i<argc; i++ )
	{
		String::Tokenize( argv[i], toks, "=" );
		if( toks.size() < 2 ) continue;

		const auto& key = toks[0];
		const auto& val = toks[1];

		if( key == "sphere" )
		{
			SphereDefinition sd;
			ParseSphereDefinition( sd, val.c_str() );
			sphere_defs.push_back( sd );
		}
		else if( key == "plane" )
		{
			PlaneDefinition pd;
			ParsePlaneDefinition( pd, val.c_str() );
			plane_defs.push_back( pd );
		}
		else if( key == "lipid" )
		{
			LAMMPS::Molecule mol;
			ParseLipidDefinition( val.c_str(), mol );
			for( size_t i=0; i<mol.atoms.size(); i++ )
			{
				auto& a = mol.atoms[i];
				a.q = 0.0;
				a.x = 0.0;
				a.y = 0.0;
				a.z = bond_length*i;
			}
			molecules.push_back( mol );
		}
		else if( key == "bond_length" )
		{
			if( (String::ToReal(val,bond_length)!=String::ReturnValue::OK) || (bond_length<=0.0) )
			{
				printf( "Bad bond length '%s'\n", val.c_str() );
				exit( -1 );
			}
		}
		else if( key == "leaflet_separation" )
		{
			if( (String::ToReal(val,leaflet_separation)!=String::ReturnValue::OK) || (leaflet_separation<=0.0) )
			{
				printf( "Bad leaflet separation '%s'\n", val.c_str() );
				exit( -1 );
			}
		}
	}

	if( leaflet_separation < 0.0 ) leaflet_separation = bond_length;

	//
	// Print some information
	//
	printf( "\n" );
	printf( "*\n" );
	printf( "* Input parameters\n" );
	printf( "*\n" );
	printf( "\n" );

	for( const auto& x : molecules )
	{
		printf( "Molecule:\n" );

		printf( "  atom types: " );
		for( const auto& a : x.atoms ) printf( "%d ", a.type );
		printf( "\n" );

		printf( "  bonds: "  );
		for( const auto& b : x.bonds ) printf( "[%d,%d,%d] ", b.type, b.i, b.j );
		printf( "\n" );

		printf( "  angles: "  );
		for( const auto& a : x.angles ) printf( "[%d,%d,%d,%d] ", a.type, a.i, a.j, a.k );
		printf( "\n" );

		printf( "  coords:\n"  );
		for( const auto& a : x.atoms ) printf( "    %g %g %g\n", a.x, a.y, a.z );
	}
	printf( "\n" );

	for( const auto& x : sphere_defs )
	{
		printf( "Sphere defined:\n" );
		printf( "  type: %s\n", (x.is_monolayer) ? ("monolayer") : ("bilayer") );
		printf( "  mol types/props:\n" );
		for( size_t i=0, N=x.mol_types.size(); i<N; i++ ) printf( "    %d %g\n", x.mol_types[i]+1, x.mol_props[i] );
		printf( "  outer_r: %g\n", x.outer_r );
		printf( "  target APL: %g\n", x.APL );
		printf( "\n" );
	}

	for( const auto& x : plane_defs )
	{
		printf( "Plane defined:\n" );
		printf( "  type: %s\n", (x.is_monolayer) ? ("monolayer") : ("bilayer") );
		printf( "  mol types/props:\n" );
		for( size_t i=0, N=x.mol_types.size(); i<N; i++ ) printf( "    %d %g\n", x.mol_types[i]+1, x.mol_props[i] );
		printf( "  N: %d\n", x.N );
		printf( "  APL: %g\n", x.APL );
		printf( "\n" );
	}

	printf( "bond_length = %g\n", bond_length );
	printf( "leaflet_separation = %g\n", leaflet_separation );

	printf( "\n" );
	printf( "*\n" );
	printf( "* Generating surfaces\n" );
	printf( "*\n" );
	printf( "\n" );

	//
	// Add spheres to surfaces
	// Fixing lipid count via area of outer/inner radii
	// and area per lipid is a bad idea: breaks totally
	// for small vesicles.
	// 
	// Instead, preserve DENSITY, rho, in both monolayers?
	// Assuming dz is the width of a single monolayer:
	// 
	// Volume per lipid, VPL = APL*Lz
	// 
	// Volume enclosed by spheres of radius r0 and r1 :
	// 	= 4/3.pi.(r1^3) - 4/3.pi.(r0^3)
	// 	= 4.3.pi.( r1^3 - r0^3 )
	// 
	// Therefore number of lipids in that volume should be:
	// 	4.3.pi.( r1^3 - r0^3 ) / (APL.Lz)
	//
	{
		for( size_t def_i=0, def_N=sphere_defs.size(); def_i<def_N; def_i++ )
		{
			printf( "Sphere.\n" );

			const auto& def = sphere_defs[def_i];

			//
			// Determine longest molecule used.
			//
			double lipid_length = 0.0;
			{
				for( auto& t : def.mol_types )
				{
					if( t > (int)molecules.size() )
					{
						printf( "Bad molecule type %d in sphere %d\n", t, (int)def_i );
						exit( -1 );
					}
					const auto& mol = molecules[t];
					double len = bond_length * (mol.atoms.size()-1);
					lipid_length = std::max( lipid_length, std::max(len,bond_length) );
				}
			}

			double r = def.outer_r;
			double A = 4.0 * M_PI * (r*r);
			double V = GapVolume( r-lipid_length, r );
			int N = (int) ( V / (def.APL*lipid_length) );

			printf( "  Outer leaflet : radius %g, %d lipids, VPL = %g, APL %g\n", r, N, V/N, A/N );

			MolecularSphere sph;

			sph.Generate( 0, N, r, {0,0,0} );
			sph.RandomAssign( def.mol_types, def.mol_props );
			for( auto& d : sph.dir ) d = d * -1;
			surfaces.push_back( sph );

			if( def.is_monolayer == false )
			{
				r = def.outer_r - (2.0*lipid_length) - leaflet_separation;
				A = 4.0 * M_PI * (r*r);
				V = GapVolume( r, r+lipid_length );
				N = (int) ( V / (def.APL*lipid_length) );

				if( r < 0.0 )
				{
					printf( "Bad inner radius for bilayer sphere: %g in sphere %d\n", r, (int)def_i );
					exit( -1 );
				}

				printf( "  Inner leaflet : radius %g, %d lipids, VPL = %g, APL %g\n", r, N, V/N, A/N );

				sph.Generate( 0, N, r, {0,0,0} );
				sph.RandomAssign( def.mol_types, def.mol_props );
				for( auto& d : sph.dir ) d = d * +1;
				surfaces.push_back( sph );
			}

			printf( "\n" );
		}
	}

	//
	// Add planes to surfaces
	//
	{
		for( size_t def_i=0, def_N=plane_defs.size(); def_i<def_N; def_i++ )
		{
			printf( "Plane.\n" );

			const auto& def = plane_defs[def_i];

			//
			// Determine longest molecule used.
			//
			double lipid_length = 0.0;
			{
				for( auto& t : def.mol_types )
				{
					if( t > (int)molecules.size() )
					{
						printf( "Bad molecule type %d in plane %d\n", t, (int)def_i );
						exit( -1 );
					}
					const auto& mol = molecules[t];
					double len = bond_length * (mol.atoms.size()-1);
					lipid_length = std::max( lipid_length, std::max(len,bond_length) );
				}
			}

			int N = def.N;
			double APL = def.APL;

			double delta = sqrt(APL);
			double L = delta*N;


			MolecularPlane pln;
			pln.Generate( 0, MolecularPlane::Axis::Z, N,N, -L/2,+L/2, -L/2,+L/2 );
			pln.RandomAssign( def.mol_types, def.mol_props );
			for( auto& d : pln.dir ) d = d * -1;
			if( def.is_monolayer == false )
			{
				for( auto& p : pln.pos )
				{
					p.z += lipid_length + leaflet_separation/2;
				}
			}
			surfaces.push_back( pln );
			printf( "  Upper leaflet : %d x %d on [%g,%g] x [%g,%g], APL %g\n", N,N, -L/2,+L/2, -L/2,+L/2, (L*L)/(N*N) );

			if( def.is_monolayer == false )
			{
				MolecularPlane pln;
				pln.Generate( 0, MolecularPlane::Axis::Z, N,N, -L/2,+L/2, -L/2,+L/2 );
				pln.RandomAssign( def.mol_types, def.mol_props );
				for( auto& d : pln.dir ) d = d * +1;
				for( auto& p : pln.pos )
				{
					p.z -= lipid_length + leaflet_separation/2;
				}
				surfaces.push_back( pln );
				printf( "  Lower leaflet : %d x %d on [%g,%g] x [%g,%g], APL %g\n", N,N, -L/2,+L/2, -L/2,+L/2, (L*L)/(N*N) );
			}

			printf( "\n" );
		}
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
		PDB::WriteMolecularSurface( f, surfaces, molecules );
		fclose( f );
	}

	//
	// Figure out margin?
	//
	double minx = 0, maxx = 0, miny = 0, maxy = 0, minz = 0, maxz = 0;
	{
		bool first = true;

		for( const auto& s : sphere_defs )
		{
			double minx_ = -s.outer_r;
			double maxx_ = +s.outer_r;

			double miny_ = -s.outer_r;
			double maxy_ = +s.outer_r;

			double minz_ = -s.outer_r;
			double maxz_ = +s.outer_r;

			if( first == true )
			{
				minx = minx_;
				maxx = maxx_;

				miny = miny_;
				maxy = maxy_;

				minz = minz_;
				maxz = maxz_;

				first = false;
				continue;
			}

			if( minx_ < minx ) minx = minx_;
			if( maxx_ > maxx ) maxx = maxx_;

			if( miny_ < miny ) miny = miny_;
			if( maxy_ > maxy ) maxy = maxy_;

			if( minz_ < minz ) minz = minz_;
			if( maxz_ > maxz ) maxz = maxz_;
		}

		for( const auto& p : plane_defs )
		{
			double L = sqrt(p.APL)*p.N;

			double minx_ = -L/2;
			double maxx_ = +L/2;

			double miny_ = -L/2;
			double maxy_ = +L/2;

			if( first )
			{
				minx = minx_;
				maxx = maxx_;

				miny = miny_;
				maxy = maxy_;

				first = false;
				continue;
			}

			if( minx_ < minx ) minx = minx_;
			if( maxx_ > maxx ) maxx = maxx_;

			if( miny_ < miny ) miny = miny_;
			if( maxy_ > maxy ) maxy = maxy_;
		}
	}

	//
	// Expand surface definition using molecules
	//
	{
		FILE* f = fopen( "molecules.lammps_config", "w" );
		WriteLAMMPSMolecules( f, surfaces, molecules, minx,maxx, miny,maxy, minz,maxz );
		fclose( f );
	}

}
