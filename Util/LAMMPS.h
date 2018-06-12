/*
	Author: John Grime, The University of Chicago.
*/

//
// Routines specialized for IO and manipulation of LAMMPS data.
//

#if !defined(UTIL_LAMMPS_DEFINED)
#define UTIL_LAMMPS_DEFINED

#include <math.h>
#include <vector>
#include <map>

#include "String.h"

namespace Util
{
namespace LAMMPS
{

//
// Used throughout
//
using Real = double;


//
// Atom
//
struct Atom
{
	int serial;
	int type;
	Real x, y, z;
	Real q;
};


//
// Bond
//
struct Bond
{
	int type, i, j;
};


//
// Angle
//
struct Angle
{
	int type, i, j, k;
};


//
// Cuboid, with utility methods.
//
struct Box
{
	Real minx, maxx;
	Real miny, maxy;
	Real minz, maxz;

	void Translate( Real& dx, Real& dy, Real& dz );

	//
	// Wrap coordinate and separation vector into box.
	// If PBC == nullptr, assume periodicity on all axes.
	//
	void WrapPosition(   Real& x, Real& y, Real& z, const int* PBC = nullptr ) const;
	void WrapSeparation( Real& x, Real& y, Real& z, const int* PBC = nullptr ) const;

	//
	// Check if a point lies inside the box. We assume range is [min,max) on axes.
	//
	bool IsInside( Real x, Real y, Real z ) const;
};


//
// A molecule contains atoms and topological information. By storing the
// minimum atom serial in the molecule, we can trivially express the bonds,
// angles etc as being relative to that minimum to removes global indices and
// create a modular system: we can add and remove topological info, molecules etc
// without needing to manually adjust all the indices to retain consistency.
//
struct Molecule
{
	int id;         // original molecule id, if available.
	int min_serial; // for trivial "localization" of topology info
	Box bounds;     // molecular bounding box

	std::map<int,int> type_contents; // unique particle types in molecule

	std::vector< Atom > atoms;
	std::vector< Bond > bonds;
	std::vector< Angle > angles;

	Molecule();

	//
	// Add components to the molecule. AddAtom() updates "min_serial", and
	// "type_contents" as we go, to retain consistency.
	//
	void AddAtom( int serial, int type, Real q, Real x, Real y, Real z );
	void AddBond( int type, int i, int j );
	void AddAngle( int type, int i, int j, int k );

	//
	// Translate all atoms by dx,dy,dz. Also adjusts "bounds" by default.
	//
	void Translate( Real dx, Real dy, Real dz, bool shift_bounds = true );

	//
	// Wrap or unwrap molecular coordinates in the specified box. Where
	// PBC == nullptr, assume periodicity on all axes.
	//
	// WARNING:
	//  - These do not modify "bounds"
	//  - Unwrap() only works for small molecules. Recursive version needed
	//    for e.g. long polymers that wrap across boundaries multiple times.
	//
	void Wrap( const Box& cell, const int* PBC = nullptr );
	void Unwrap( size_t pivot_i, const Box& cell, const int* PBC = nullptr );

	//
	// Walk the atoms in the molecule, to recalculate the "bounds" variable.
	//
	void RecalculateBounds();
};


//
// A LAMMPS configuration.
//
class Config
{
	public:

		int PBC[3];
		Box bounds; // simulation cell
		std::vector< Real > mass; // mass of atom types; mass.size() == number of atom types present
		std::vector< Molecule > molecules;
		std::vector< std::string > bond_coeffs, angle_coeffs;

		Config();

		void ClearAll();
		void ClearMolecules();

		//
		// Extact all molecules that contain one of the interesting atom types.
		//
		void Extract( size_t N, const int* interesting_types, std::vector<Molecule>& extracted_molecules );
		void Extract( const std::vector<int>& interesting_types, std::vector<Molecule>& extracted_molecules );

		//
		// Extract molecules inside the box. If "strict" == true, only return molecules where ALL
		// particles lie inside the box. NOTE: we use the actual particle coords rather than the
		// bounding boxes (which can be predictably weird for wrapped molecules).
		//
		void Extract( const Box& b, std::vector<Molecule>& extracted_molecules, bool strict = true );

		//
		// Translate whole system by specified deltas. Also updates "bounds" by default.
		//
		void Translate( Real dx, Real dy, Real dz, bool shift_bounds = true );

		//
		// Get system centre of geometry, mass
		//
		void GetCOG( Real& x, Real& y, Real& z ) const;
		void GetCOM( Real& x, Real& y, Real& z ) const;

		//
		// Wrap/unwrap all molecules to the simulation cell. Assumes periodicity on all axes where
		// PBC == nullptr. These do not modify config "bounds", or bounds of any component molecule.
		//
		void Wrap( const int* PBC = nullptr );
		void Unwrap( const int* PBC = nullptr );

		//
		// Recalculates bounding boxes of system and any component molecules.
		//
		void RecalculateBounds();

		//
		// Add another configuration to this one.
		//
		Config& operator += ( const Config& rhs );

		//
		// These methods should not be used in the += operator, as e.g. will merge atoms
		// into the same molecules where same mol ids exist etc!
		//
		void AddAtom( int serial, int mol, int type, Real q, Real x, Real y, Real z );
		void AddBond( int type, int i, int j );
		void AddAngle( int type, int i, int j, int k );

		//
		// Check the separations between particles in bonds and angles to ensure molecules are
		// contiguous across boundaries. This is useful where you're tiling systems.
		//
		bool SanityCheckTopology( double tolerance, int& problem_molecule_index ) const;

	protected:

		std::map<int,int> mol_to_index;    // LAMMPS mol tag     => index into molecules[]
		std::map<int,int> serial_to_index; // LAMMPS atom serial => index into molecules[]

		//
		// Return index into molecules[] for the molecule with the specified LAMMPS molecule
		// id, or the molecule containing the atom with the specified LAMMPS serial number.
		//
		int mol_index_from_mol_id( int mol );
		int mol_index_from_atom_serial( int serial );

		//
		// Get centre via summation over atom positions. If use_mass == true, this
		// generates the centre of mass, otherwise it's the centre of geometry.
		//
		void get_centre( Real& x, Real& y, Real& z, bool use_mass ) const;
};


//
// File IO
//
int LoadData( FILE* f, Config& config, bool has_q = false, bool noisy = false );
int LoadTrajectory( FILE *f, Config& config, int* timestep = nullptr );
int SaveTrajectory( FILE *f, Config& config, int timestep = 0 );

int SaveData( FILE* f, const Config& config, bool has_q = false );
int SaveXYZ(  FILE* f, const Config& config, const char *heading = "" );

//
// Simplified LAMMPS trajectory frame class.
// Can be faster and lighter in memory vs Config class.
//
class TrajectoryFrame
{

	public:
	
		int timestep;
		
		double mins[3], maxs[3];
		int PBC[3];
		
		std::vector<double> x_vec,y_vec,z_vec;
		std::vector<int> atom_ids, atom_types; // per-atom
		
		void Clear()
		{
			timestep = -1;
			for( int i=0; i<3; i++ )
			{
				mins[i] = -0.5;
				maxs[i] = +0.5;
				PBC[i] = 1;
			}

			x_vec.clear();
			y_vec.clear();
			z_vec.clear();
			atom_ids.clear();
			atom_types.clear();
		}
		
		void Wrap( double &x, double &y, double &z ) const
		{
			double Lx = maxs[0]-mins[0];
			double Ly = maxs[1]-mins[1];
			double Lz = maxs[2]-mins[2];
			
			// Wrap deltas
			if( PBC[0] == 1 ) x += Lx*( (x<mins[0]) - (x>=maxs[0]) );
			if( PBC[1] == 1 ) y += Ly*( (y<mins[1]) - (y>=maxs[1]) );
			if( PBC[2] == 1 ) z += Lz*( (z<mins[2]) - (z>=maxs[2]) );
		}
		
		// Basic loading of trajectory frame
		int Load( FILE *f )
		{
			const char *delimiters = " \t\n\r";

			int maxbuf = 1023;
			char buffer[1024];
			
			int ntoks;
			std::vector< std::string > tokens;
			
			fpos_t file_position;

			std::map< std::string, int > token_indices;
			int number_of_atoms = 0; // internal check
			
			token_indices[ "id" ] = -1;
			token_indices[ "type" ] = -1;
			token_indices[ "x" ] = -1;
			token_indices[ "y" ] = -1;
			token_indices[ "z" ] = -1;

			if( f == NULL )
			{
				fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): file pointer is NULL!\n", __func__ );
				return -1;
			}
			
			Clear();

			while( fgets( buffer, maxbuf, f ) != NULL )
			{
				if( (ntoks = String::Tokenize( buffer, tokens, delimiters )) < 1 || tokens[0] != "ITEM:" ) continue;

				//
				// Trajectory frame timestep number
				//
				if( tokens[1] == "TIMESTEP" )
				{
					if( fgets( buffer, maxbuf, f ) == NULL ) continue;
					if( String::ToInteger( buffer, timestep ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert timestep '%s' into an integer.\n", __func__, buffer );
						return -1;
					}
				}
				//
				// Simulation cell
				//
				else if( tokens[1] == "BOX" )
				{
					// expects ITEM: BOX BOUNDS xx yy zz
					if( tokens.size() != 6 )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): bad box line: '%s'\n", __func__, buffer  );
						return -1;
					}
					
					for( int i=0; i<3; i++ ) PBC[i] = ( tokens[3+i][0] == 'p' ) ? 1 : 0;

					if( fgets( buffer, maxbuf, f ) == NULL ) continue;
					if( (ntoks = String::Tokenize( buffer, tokens, delimiters)) < 2 ) continue;

					//
					// Get x bounds
					//
					if( String::ToReal( tokens[0], mins[0] ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert min x value '%s' into a number.\n", __func__, tokens[0].c_str() );
						return -1;
					}
					if( String::ToReal( tokens[1], maxs[0] ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert max x value '%s' into a number.\n", __func__, tokens[1].c_str() );
						return -1;
					}

					//
					// Get y bounds
					//
					if( fgets( buffer, maxbuf, f ) == NULL ) continue;
					if( (ntoks = String::Tokenize( buffer, tokens, delimiters)) < 2 ) continue;

					if( String::ToReal( tokens[0], mins[1] ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert min y value '%s' into a number.\n", __func__, tokens[0].c_str() );
						return -1;
					}
					if( String::ToReal( tokens[1], maxs[1] ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert max y value '%s' into a number.\n", __func__, tokens[1].c_str() );
						return -1;
					}

					//
					// Get z bounds
					//
					if( fgets( buffer, maxbuf, f ) == NULL ) continue;
					if( (ntoks = String::Tokenize( buffer, tokens, delimiters)) < 2 ) continue;

					if( String::ToReal( tokens[0], mins[2] ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert min z value '%s' into a number.\n", __func__, tokens[0].c_str() );
						return -1;
					}
					if( String::ToReal( tokens[1], maxs[2] ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert max z value '%s' into a number.\n", __func__, tokens[1].c_str() );
						return -1;
					}
				}
				//
				// Number of atoms in frame
				//
				else if( tokens[1] == "NUMBER" && tokens[2] == "OF" && tokens[3] == "ATOMS" )
				{
					if( fgets( buffer, maxbuf, f ) == NULL ) continue;
					if( String::ToInteger( buffer, number_of_atoms ) != String::ReturnValue::OK )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert number of atoms '%s' into a number.\n", __func__, buffer );
						return -1;
					}
				}
				//
				// Atom coordinates etc
				//
				else if( tokens[1] == "ATOMS" )
				{
					bool xscale = false, yscale = false, zscale = false;

					std::map< std::string, int >::iterator it;
				
					//
					// Clear token indices
					//
					for( it=token_indices.begin(); it!=token_indices.end(); it++ )
					{
						it->second = -1;
					}

					//
					// Determine column indices for the appropriate atomic data
					//
					for( int i=0; i<ntoks; i++ )
					{
						std::string &tok = tokens[i];
						
						if( tok == "xu" ) { tok = "x"; }
						if( tok == "yu" ) { tok = "y"; }
						if( tok == "zu" ) { tok = "z"; }

						if( tok == "xs" ) { tok = "x"; xscale = true; }
						if( tok == "ys" ) { tok = "y"; yscale = true; }
						if( tok == "zs" ) { tok = "z"; zscale = true; }
						
						it = token_indices.find(tok);
						if( it == token_indices.end() ) continue;
						
						token_indices[tok] = i-2; // tokens start with "ITEM:", "ATOMS"
					}

					//
					// Check we've actually defined all the columns for the appropriate atom data
					//
					for( it=token_indices.begin(); it!=token_indices.end(); it++ )
					{
						if( it->second >= 0 ) continue;
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to find '%s' info in frame.\n", __func__, it->first.c_str() );
						return -1;
					}
					
					int id_i = token_indices["id"];
					int type_i = token_indices["type"];
					int x_i = token_indices["x"];
					int y_i = token_indices["y"];
					int z_i = token_indices["z"];

					//
					// Store file position, so we can rewind if we hit an ITEM line.
					//
					fgetpos( f, &file_position );

					//
					// Read atom data, using the columns we determined previously
					//
					while( fgets( buffer, maxbuf, f ) != NULL )
					{				
						ntoks = String::Tokenize( buffer, tokens, delimiters );
						if( ntoks < 1 ) continue;

						//
						// ITEM indicates we've hit the next frame; rewind th start of line and return.
						//
						if( tokens[0] == "ITEM:" )
						{
							fsetpos( f, &file_position );

							if( number_of_atoms != (int)atom_ids.size() )
							{
								fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): WARNING: odd number of atoms (read %d, expected %d)\n", __func__, (int)atom_ids.size(), number_of_atoms );
							}

							return (int)atom_ids.size();
						}
						
						int atom_id, atom_type;
						double x, y, z;

						if( String::ToInteger( tokens[id_i], atom_id ) != String::ReturnValue::OK )
						{
							fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert id '%s' into a number.\n", __func__, tokens[id_i].c_str() );
							return -1;
						}
						if( String::ToInteger( tokens[type_i], atom_type ) != String::ReturnValue::OK )
						{
							fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert type '%s' into a number.\n", __func__, tokens[type_i].c_str() );
							return -1;
						}

						if( String::ToReal( tokens[x_i], x ) != String::ReturnValue::OK )
						{
							fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert x '%s' into a number.\n", __func__, tokens[x_i].c_str() );
							return -1;
						}
						if( String::ToReal( tokens[y_i], y ) != String::ReturnValue::OK )
						{
							fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert y '%s' into a number.\n", __func__, tokens[y_i].c_str() );
							return -1;
						}
						if( String::ToReal( tokens[z_i], z ) != String::ReturnValue::OK )
						{
							fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert z '%s' into a number.\n", __func__, tokens[z_i].c_str() );
							return -1;
						}

						atom_ids.push_back( atom_id );
						atom_types.push_back( atom_type );

						if( xscale == true ) x = mins[0] + x*(maxs[0]-mins[0]);
						if( yscale == true ) y = mins[1] + y*(maxs[1]-mins[1]);
						if( zscale == true ) z = mins[2] + z*(maxs[2]-mins[2]);

						x_vec.push_back( x );
						y_vec.push_back( y );
						z_vec.push_back( z );
					}
				}
			}

			return (int)atom_ids.size();
		}		

		//
		// Write xyz file fo the specific particles in the trajectory frame as indicated by "indices".
		//
		void WriteXYZ( FILE *f, const std::vector<int>& indices ) const
		{
			if( f == nullptr ) return;

			fprintf( f, "%d\n", (int)indices.size() );
			fprintf( f, "\n" );
			for( size_t i=0, max_i=indices.size(); i<max_i; i++ )
			{
				int index = indices[i];

				if( index >= (int)atom_ids.size() ) continue; // should probably be an error!

				double x = x_vec[index];
				double y = y_vec[index];
				double z = z_vec[index];
				fprintf( f, "%s %8.3f %8.3f %8.3f\n", "?", x, y, z );
			}	
		}
};


}
}

#endif
