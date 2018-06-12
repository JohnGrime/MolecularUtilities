/*
	Author: John Grime, The University of Chicago.
*/

#include "LAMMPS.h"

namespace Util
{
namespace LAMMPS
{

void Box::Translate( Real& dx, Real& dy, Real& dz )
{
	minx += dx;
	maxx += dx;

	miny += dy;
	maxy += dy;

	minz += dz;
	maxz += dz;
}

//
// PBC == nullptr = assume periodic on all boundaries.
// All internal coords are double, for precision.
//
void Box::WrapPosition( Real& x, Real& y, Real& z, const int* PBC ) const
{
	double L[] = { maxx-minx, maxy-miny, maxz-minz };
	double n[] = { (x-minx)/L[0], (y-miny)/L[1], (z-minz)/L[2] };

	for( int i=0; i<3; i++ )
	{
		if( (PBC!=nullptr) && (PBC[i]==0) ) continue;
		if( n[i] < 0.0 ) n[i] += 1.0;
		else if( n[i] >= 1.0 ) n[i] -= 1.0;
	}

	x = minx + (n[0]*L[0]);
	y = miny + (n[1]*L[1]);
	z = minz + (n[2]*L[2]);
}

//
// PBC == nullptr = assume periodic on all boundaries.
// All internal coords are double, for precision.
//
void Box::WrapSeparation( Real& dx, Real& dy, Real& dz, const int* PBC ) const
{
	double L[] = { maxx-minx, maxy-miny, maxz-minz };
	double n[] = { dx/L[0], dy/L[1], dz/L[2] };

	for( int i=0; i<3; i++ )
	{
		if( (PBC!=nullptr) && (PBC[i]==0) ) continue;
		if( n[i] < -0.5 ) n[i] += 1.0;
		else if( n[i] >= 0.5 ) n[i] -= 1.0;
	}

	dx = n[0]*L[0];
	dy = n[1]*L[1];
	dz = n[2]*L[2];
}

bool Box::IsInside( Real x, Real y, Real z ) const
{
	if( (x<minx) || (x>=maxx) ) return false;
	if( (y<miny) || (y>=maxy) ) return false;
	if( (z<minz) || (z>=maxz) ) return false;
	return true;
}


Molecule::Molecule()
{
	id = -1;
	min_serial = -1;
	bounds = { 0,0, 0,0, 0,0 };
}

void Molecule::AddAtom( int serial, int type, Real q, Real x, Real y, Real z )
{
	Atom a;

	a.serial = serial;
	a.type = type;
	a.q = q;
	a.x = x;
	a.y = y;
	a.z = z;

	if( min_serial == -1 )
	{
		min_serial = serial;
		bounds = { x,x, y,y, z,z };
	}
	else
	{
		if( serial < min_serial ) min_serial = serial;

		if( x < bounds.minx ) bounds.minx = x;
		if( x > bounds.maxx ) bounds.maxx = x;

		if( y < bounds.miny ) bounds.miny = y;
		if( y > bounds.maxy ) bounds.maxy = y;

		if( z < bounds.minz ) bounds.minz = z;
		if( z > bounds.maxz ) bounds.maxz = z;
	}

	atoms.push_back( a );

	type_contents[type] = 1;
}
//
// Note: we retain the "global" indices in the angle, to allow
// progressive generation of the molecules. Subtract min_serial
// from the indices to convert into molecule-local indices.
//
void Molecule::AddBond( int type, int i, int j )
{
	Bond b;
	b.type = type;
	b.i = i;
	b.j = j;
	bonds.push_back( b );		
}
void Molecule::AddAngle( int type, int i, int j, int k )
{
	Angle a;
	a.type = type;
	a.i = i;
	a.j = j;
	a.k = k;
	angles.push_back( a );		
}

void Molecule::Translate( Real dx, Real dy, Real dz, bool shift_bounds )
{
	for( size_t i=0, N=atoms.size(); i<N; i++ )
	{
		atoms[i].x += dx;
		atoms[i].y += dy;
		atoms[i].z += dz;
	}
	if( shift_bounds ) bounds.Translate( dx, dy, dz );
}

void Molecule::RecalculateBounds()
{
	for( size_t i=0, N=atoms.size(); i<N; i++ )
	{
		const auto& a = atoms[i];
		Real x = a.x;
		Real y = a.y;
		Real z = a.z;

		if( i == 0 )
		{
			bounds = { x,x, y,y, z,z };
			continue;
		}

		if( x < bounds.minx ) bounds.minx = x;
		if( x > bounds.maxx ) bounds.maxx = x;

		if( y < bounds.miny ) bounds.miny = y;
		if( y > bounds.maxy ) bounds.maxy = y;

		if( z < bounds.minz ) bounds.minz = z;
		if( z > bounds.maxz ) bounds.maxz = z;
	}
}

//
// PBC == nullptr = all boundaries periodic. All internal coords double prec.
// BE CAREFUL: this can obviously change the bounds of the molecule, but we
// do not modify the information in the Molecule "bounds" structure!
//
void Molecule::Wrap( const Box& cell, const int* PBC )
{
	for( auto& a : atoms ) cell.WrapPosition( a.x, a.y, a.z, PBC );
}

//
// PBC == nullptr = all boundaries periodic. All internal coords double prec.
// BE CAREFUL: this can obviously change the bounds of the molecule, but we
// do not modify the information in the Molecule "bounds" structure!
//
void Molecule::Unwrap( size_t pivot_i, const Box& cell, const int* PBC )
{
	if( pivot_i >= atoms.size() )
	{
		printf( "Molecule::%s(): bad pivot %d, require 0 < pivot <= %d\n", __func__, (int)pivot_i, (int)atoms.size() );
		exit( -1 );
	}

	double x0 = atoms[pivot_i].x;
	double y0 = atoms[pivot_i].y;
	double z0 = atoms[pivot_i].z;

	//
	// Convert coords to be relative to pivot atom
	//
	for( size_t i=0, N=atoms.size(); i<N; i++ )
	{
		if( i == pivot_i ) continue;

		auto& a = atoms[i];

		double dx = a.x - x0;
		double dy = a.y - y0;
		double dz = a.z - z0;

		cell.WrapSeparation( dx, dy, dz, PBC );

		a.x = x0+dx;
		a.y = y0+dy;
		a.z = z0+dz;
	}
}




Config::Config()
{
	ClearAll();
}

void Config::ClearAll()
{
	mass.clear();
	bond_coeffs.clear();
	angle_coeffs.clear();

	ClearMolecules();

	bounds = { -0.5,+0.5, -0.5,+0.5, -0.5,+0.5 };
	PBC[0] = PBC[1] = PBC[2] = 1;
}
void Config::ClearMolecules()
{
	molecules.clear();
	mol_to_index.clear();
	serial_to_index.clear();
}

//
// Extacts all molecules that contain one of the interesting particle types.
//
void Config::Extract( size_t N, const int* interesting_types, std::vector<Molecule>& extracted_molecules )
{
	std::vector<Molecule> temp_mols;

	temp_mols.clear();

	if( N>0 && interesting_types==nullptr ) return;

	for( const auto& m : molecules )
	{
		for( size_t j=0; j<N; j++ )
		{
			int type = interesting_types[j];
			if( m.type_contents.find(type) != m.type_contents.end() )
			{
				temp_mols.push_back( m );
				break;
			}
		}
	}

	extracted_molecules = temp_mols;
}
void Config::Extract( const std::vector<int>& interesting_types, std::vector<Molecule>& extracted_molecules )
{
	Extract( interesting_types.size(), &interesting_types[0], extracted_molecules );
}

//
// Extracts molecules according to the geometric box region
// If "strict" == true, only extract molecules where ALL particles inside the region.
// NOTE: we use the particle coords, rather than the bounding boxes, as the latter
// can go weird for wrapped molecules (which provide a bounding box spanning the entire
// simulation domain on the axis they wrap across).
//
void Config::Extract(
	const Box& b,
	std::vector<Molecule>& extracted_molecules,
	bool strict
	)
{
	std::vector<Molecule> temp_mols;

	temp_mols.clear();

	for( const auto& m : molecules )
	{
		size_t N_inside = 0;
		size_t N_outside = 0;

		for( const auto& a : m.atoms )
		{
			int OK = 0;

			if( (a.x>=b.minx) && (a.x<=b.maxx) ) OK++;
			if( (a.y>=b.miny) && (a.y<=b.maxy) ) OK++;
			if( (a.z>=b.minz) && (a.z<=b.maxz) ) OK++;

			if( OK == 3 ) N_inside++;
			else N_outside++;
		}

		bool add = (strict) ? (N_inside==m.atoms.size()) : (N_inside>0);
		if( add ) temp_mols.push_back(m);
	}
	extracted_molecules = temp_mols;
}

void Config::Translate( Real dx, Real dy, Real dz, bool shift_bounds )
{
	for( auto& m : molecules ) m.Translate( dx, dy, dz );
	if( shift_bounds ) bounds.Translate( dx, dy, dz );
}

void Config::GetCOG( Real& x, Real& y, Real& z ) const { get_centre( x,y,z, false ); }
void Config::GetCOM( Real& x, Real& y, Real& z ) const { get_centre( x,y,z, true ); }

void Config::RecalculateBounds()
{
	for( size_t i=0, max_i=molecules.size(); i<max_i; i++ )
	{
		molecules[i].RecalculateBounds();

		if( i == 0 )
		{
			bounds = molecules[i].bounds;
			continue;
		}

		const auto& b = molecules[i].bounds;

		if( b.minx < bounds.minx ) bounds.minx = b.minx;
		if( b.maxx > bounds.maxx ) bounds.maxx = b.maxx;

		if( b.miny < bounds.miny ) bounds.miny = b.miny;
		if( b.maxy > bounds.maxy ) bounds.maxy = b.maxy;

		if( b.minz < bounds.minz ) bounds.minz = b.minz;
		if( b.maxz > bounds.maxz ) bounds.maxz = b.maxz;
	}
}

void Config::Wrap( const int* PBC )
{
	for( auto& m : molecules ) m.Wrap( bounds, PBC );
}

//
// PBC == nullptr = wrap on all bounds.
// BE CAREFUL: this can change the effective system bounds, but will not
// modify the values in the Config "bounds" structure!
//
void Config::Unwrap( const int* PBC )
{
	for( auto& m : molecules ) m.Unwrap( 0, bounds, PBC );
}

Config& Config::operator += ( const Config& rhs )
{
	//
	// Check whether we need to adjust the box dims.
	//
	{
		const auto& b = rhs.bounds;

		if( b.minx < bounds.minx ) bounds.minx = b.minx;
		if( b.maxx > bounds.maxx ) bounds.maxx = b.maxx;

		if( b.miny < bounds.miny ) bounds.miny = b.miny;
		if( b.maxy > bounds.maxy ) bounds.maxy = b.maxy;

		if( b.minz < bounds.minz ) bounds.minz = b.minz;
		if( b.maxz > bounds.maxz ) bounds.maxz = b.maxz;
	}

	//
	// Copy molecules. As we declared the bonds etc to be local, we should not need to
	// modify any topological indices.
	//
	for( const auto& m : rhs.molecules ) molecules.push_back( m );

	return *this;
}

bool Config::SanityCheckTopology( double tolerance, int& problem_molecule_index ) const
{
	size_t bonds_checked = 0;
	size_t angles_checked = 0;

	double tol2 = tolerance*tolerance;

	for( size_t mol_i=0; mol_i<molecules.size(); mol_i++ )
	{
		const auto& m = molecules[mol_i];

		for( size_t bond_i=0; bond_i<m.bonds.size(); bond_i++ )
		{
			int local_i = m.bonds[bond_i].i - m.min_serial; // local offsets of bond atoms into current molecule
			int local_j = m.bonds[bond_i].j - m.min_serial;
			const auto& a1 = m.atoms[local_i];
			const auto& a2 = m.atoms[local_j];

			double dx = a1.x - a2.x;
			double dy = a1.y - a2.y;
			double dz = a1.z - a2.z;

			bounds.WrapSeparation( dx, dy, dz );
			double dr2 = dx*dx + dy*dy + dz*dz;

			if( dr2 > tol2 )
			{
				printf( "\n" );
				printf( "Bad bond %d in molecule %d: local i=%d j=%d, type %d, dr = %g\n", (int)bond_i, (int)mol_i, local_i, local_j, m.bonds[bond_i].type, sqrt(dr2) );
				printf( "dx, dy, dz: %g, %g, %g\n", dx, dy, dz );
				printf( "Involving:\n" );
				printf( "%d %d %f %f %f\n", (int)(m.min_serial+local_i), a1.type, a1.x, a1.y, a1.z );
				printf( "%d %d %f %f %f\n", (int)(m.min_serial+local_j), a2.type, a2.x, a2.y, a2.z );
				printf( "All atoms in molecule:\n" );
				for( size_t i=0; i<m.atoms.size(); i++ )
				{
					const auto&a = m.atoms[i];
					printf( "%d %d %f %f %f\n", (int)(m.min_serial+i), a.type, a.x, a.y, a.z );
				}

				printf( "\n" );
			
				problem_molecule_index = (int)mol_i;
				return false;
			}
			bonds_checked++;
		}


		for( size_t angle_i=0; angle_i<m.angles.size(); angle_i++ )
		{
			int local_i = m.angles[angle_i].i - m.min_serial; // local offsets of angle atoms into current molecule
			int local_j = m.angles[angle_i].j - m.min_serial;
			int local_k = m.angles[angle_i].k - m.min_serial;
			const auto& a1 = m.atoms[local_i];
			const auto& a2 = m.atoms[local_j];
			const auto& a3 = m.atoms[local_k];

			double dx = a1.x - a2.x;
			double dy = a1.y - a2.y;
			double dz = a1.z - a2.z;

			bounds.WrapSeparation( dx, dy, dz );
			double dr2 = dx*dx + dy*dy + dz*dz;

			if( dr2 > tol2 )
			{
				printf( "\n" );
				printf( "Bad angle %d in molecule %d: type %d, dr = %g\n", (int)angle_i, (int)mol_i, m.angles[angle_i].type, sqrt(dr2) );
				for( size_t i=0; i<m.atoms.size(); i++ )
				{
					const auto&a = m.atoms[i];
					printf( "%d %d %f %f %f\n", (int)(m.min_serial+i), a.type, a.x, a.y, a.z );
				}
				printf( "\n" );
				exit( -1 );
			}

			dx = a3.x - a2.x;
			dy = a3.y - a2.y;
			dz = a3.z - a2.z;

			bounds.WrapSeparation( dx, dy, dz );
			dr2 = dx*dx + dy*dy + dz*dz;

			if( dr2 > tol2 )
			{
				printf( "\n" );
				printf( "Bad angle %d in molecule %d: type %d, dr = %g\n", (int)angle_i, (int)mol_i, m.angles[angle_i].type, sqrt(dr2) );
				for( size_t i=0; i<m.atoms.size(); i++ )
				{
					const auto&a = m.atoms[i];
					printf( "%d %d %f %f %f\n", (int)(m.min_serial+i), a.type, a.x, a.y, a.z );
				}
				printf( "\n" );

				problem_molecule_index = (int)mol_i;
				return false;
			}
			angles_checked++;
		}
	}
	printf( "checked %d bonds, %d angles.\n", (int)bonds_checked, (int)angles_checked );
	return true;
}



void Config::AddAtom( int serial, int mol, int type, Real q, Real x, Real y, Real z )
{
	//
	// Which molecule do we use?
	//
	int index = mol_index_from_mol_id( mol );

	if( index == -1 )
	{
		Molecule m;
		m.id = mol;
		index = (int)molecules.size();
		mol_to_index[mol] = index;
		molecules.push_back( m );
	}

	auto& m = molecules[index];
	m.AddAtom( serial, type, q, x, y, z );
	serial_to_index[serial] = index;
}
void Config::AddBond( int type, int i, int j )
{
	//
	// Figure out which molecule this bond should be placed in.
	//
	int i_mol = mol_index_from_atom_serial( i );
	int j_mol = mol_index_from_atom_serial( j );

	if( i_mol == -1 || j_mol == -1 )
	{
		printf( "Unable to locate molecule for bond : %d %d %d\n", type, i, j );
		exit( -1 );
	}
	if( i_mol != j_mol )
	{
		printf( "Mismatch between location of molecules: %d vs %d, from %d and %d\n", i_mol, j_mol, i, j );
		exit( -1 );
	}

	auto& m = molecules[i_mol];
	m.AddBond( type, i, j );
}
void Config::AddAngle( int type, int i, int j, int k )
{
	//
	// Figure out which molecule this angle should be placed in.
	//
	int i_mol = mol_index_from_atom_serial( i );
	int j_mol = mol_index_from_atom_serial( j );
	int k_mol = mol_index_from_atom_serial( k );

	if( i_mol == -1 || j_mol == -1 || k_mol == -1 )
	{
		printf( "Unable to locate molecule for angle : %d %d %d %d\n", type, i, j, k );
		exit( -1 );
	}
	if( i_mol != j_mol || i_mol != k_mol )
	{
		printf( "Mismatch between location of molecules: %d vs %d vs %d\n", i_mol, j_mol, k_mol );
		exit( -1 );
	}

	auto& m = molecules[i_mol];
	m.AddAngle( type, i, j, k );
}

int Config::mol_index_from_mol_id( int mol )
{
	auto it = mol_to_index.find( mol );
	if( it == mol_to_index.end() ) return -1;
	return it->second;
}
int Config::mol_index_from_atom_serial( int serial )
{
	auto it = serial_to_index.find( serial );
	if( it == serial_to_index.end() ) return -1;
	return it->second;
}

void Config::get_centre( Real& x, Real& y, Real& z, bool use_mass ) const
{
	double M = 0.0;
	double x_ = 0.0, y_ = 0.0, z_ = 0.0;
	size_t N = 0;

	x = y = z = 0.0;

	for( const auto& m : molecules )
	{
		for( const auto& a : m.atoms )
		{
			double factor = 1.0; // assign equal weighting to all atoms

			if( use_mass )
			{
				if( a.type > (int)mass.size() )
				{
					printf( "Config::%s(): bad type %d, must be 1 <= type < %d\n", __func__, a.type, (int)mass.size() );
					exit( -1 );
				}
				factor = mass[a.type];
			}

			x_ += a.x * factor;
			y_ += a.y * factor;
			z_ += a.z * factor;

			M += factor;
			N++;
		}
	}

	if( N < 1 ) return;

	// Try to avoid inaccuracies via accumulation etc
	x = x_ / ( (use_mass) ? (M) : (N) );
	y = y_ / ( (use_mass) ? (M) : (N) );
	z = z_ / ( (use_mass) ? (M) : (N) );
}







int LoadData( FILE* f, Config& config, bool has_q, bool noisy )
{
	if( f == nullptr ) return -1;

	int line_no = 0;

	const int max_line = 1024;
	char line[max_line+1];

	std::string state = "";
	std::vector<std::string> tokens;

	config.ClearAll();

	size_t min_atom_tokens = 6;
	if( has_q == true ) min_atom_tokens = 7;

	while( fgets( line, max_line, f ) != nullptr )
	{
		line_no++;

		if( line_no == 1 ) continue; // skip header

		//
		// Remove any leading/trailing whitespace. This also
		// removes any pesky linefeed formatting, hopefully.
		//
		if( String::Strip( line, " \t\n\r" ) != String::ReturnValue::OK )
		{
			printf( "Error on line %d: '%s'\n", line_no, line );
			printf( "Unable to strip string\n" );
			return -1;
		}
		if( line[0] == '\0' ) continue;

		//
		// Get tokens
		//
		if( String::Tokenize( line, tokens, " \t\n\r" ) < 1 )
		{
			printf( "Error on line %d: '%s'\n", line_no, line );
			printf( "Unable to tokenize string\n" );
			return -1;
		}
		if( tokens.size() < 1 ) continue;

		//
		// Check for state transitions
		//
		{
			bool change = false;
			if( tokens[0] == "Masses" ) change = true;
			else if( tokens[0] == "Atoms" ) change = true;
			else if( tokens[0] == "Velocities" ) change = true;
			else if( tokens[0] == "Bonds" ) change = true;
			else if( tokens[0] == "Angles" ) change = true;
			
			if( change )
			{
				//printf( "State change: '%s' => '%s'\n", state.c_str(), tokens[0].c_str() );
				state = tokens[0];
				continue;
			}

			if( tokens.size() >= 2 && tokens[1] == "Coeffs" )
			{
				std::string new_state;

				new_state = tokens[0];
				new_state += " ";
				new_state += tokens[1];

				if( noisy ) printf( "State change: '%s' => '%s', via '%s'\n", state.c_str(), new_state.c_str(), line );

				state = new_state;
				continue;
			}
		}


		//
		// Not in a defined state?
		//
		if( state == "" )
		{
			//
			// Check for freeform info - we only need "atom types" and box dims.
			//
			if( tokens.size() > 2 )
			{
				if( tokens[1] == "atom" && tokens[2] == "types" )
				{
					int N = -1;
					if( String::ToInteger( tokens[0], N ) != String::ReturnValue::OK )
					{
						printf( "Error on line %d: '%s'\n", line_no, line );
						printf( "Unable to convert atom types string\n" );
						return -1;
					}
					config.mass.resize( N+1 ); // +1, as we're keeping the unit-based indexing
					for( size_t i=0, max_i=config.mass.size(); i<max_i; i++ ) config.mass[i] = 0.0;
					continue;
				}

				//
				// Box bounds
				//
				if( tokens[2] == "xlo" || tokens[2] == "ylo" || tokens[2] == "zlo" )
				{
					Real lo = -1.0;
					Real hi = -1.0;
					if( String::ToReal( tokens[0], lo ) != String::ReturnValue::OK )
					{
						printf( "Error on line %d: '%s'\n", line_no, line );
						printf( "Unable to convert low value\n" );
						return -1;
					}
					if( String::ToReal( tokens[1], hi ) != String::ReturnValue::OK )
					{
						printf( "Error on line %d: '%s'\n", line_no, line );
						printf( "Unable to convert high value\n" );
						return -1;
					}

					if( tokens[2] == "xlo" )
					{
						config.bounds.minx = lo;
						config.bounds.maxx = hi;
					}
					else if( tokens[2] == "ylo" )
					{
						config.bounds.miny = lo;
						config.bounds.maxy = hi;
					}
					else
					{
						config.bounds.minz = lo;
						config.bounds.maxz = hi;
					}

					if( noisy ) printf( "%s : %f %f\n", tokens[2].c_str(), lo, hi );

					continue;
				}
			}
		}

		//
		// We're in a defined state, so process info.
		//

		if( state == "Masses" )
		{
			int index = -1;
			Real m = 0.0;

			if( String::ToInteger( tokens[0], index ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert index\n" );
				return -1;
			}

			if( String::ToReal( tokens[1], m ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert mass\n" );
				return -1;
			}

			if( index >= (int)config.mass.size() )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Bad index %d : 1 <= index < %d\n", index, (int)config.mass.size() );
				return -1;
			}

			if( noisy ) printf( "Mass %d : %f\n", index, m );

			config.mass[index] = m;

			continue;
		}

		if( state == "Atoms" )
		{
			if( tokens.size() < min_atom_tokens )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Too few tokens: got %d, wanted at least %d\n", (int)tokens.size(), (int)min_atom_tokens );
				return -1;
			}

			int token_no;
			int serial, mol, type;
			Real q, x, y, z;

			token_no = 0;
			serial = mol = type = 1;
			q = x = y = z = 0.0;

			if( String::ToInteger( tokens[token_no], serial ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert serial token\n" );
				return -1;
			}
			token_no++;

			if( String::ToInteger( tokens[token_no], mol ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert mol token\n" );
				return -1;
			}
			token_no++;

			if( String::ToInteger( tokens[token_no], type ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert type token\n" );
				return -1;
			}
			token_no++;

			if( has_q )
			{
				if( String::ToReal( tokens[token_no], q ) != String::ReturnValue::OK )
				{
					printf( "Error on line %d: '%s'\n", line_no, line );
					printf( "Unable to convert q token\n" );
					return -1;
				}
				token_no++;
			}

			if( String::ToReal( tokens[token_no], x ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert x token\n" );
				return -1;
			}
			token_no++;

			if( String::ToReal( tokens[token_no], y ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert y token\n" );
				return -1;
			}
			token_no++;

			if( String::ToReal( tokens[token_no], z ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert z token\n" );
				return -1;
			}
			token_no++;

			/*
			if( tokens.size() > min_atom_tokens )
			{
				double L[] = {
					config.bounds.maxx - config.bounds.minx,
					config.bounds.maxy - config.bounds.miny,
					config.bounds.maxz - config.bounds.minz };
				int ix, iy, iz;

				if( String::ToInteger( tokens[token_no], ix ) != String::ReturnValue::OK )
				{
					printf( "Error on line %d: '%s'\n", line_no, line );
					printf( "Unable to convert x image flag\n" );
					return -1;
				}
				token_no++;

				if( String::ToInteger( tokens[token_no], iy ) != String::ReturnValue::OK )
				{
					printf( "Error on line %d: '%s'\n", line_no, line );
					printf( "Unable to convert y image flag\n" );
					return -1;
				}
				token_no++;

				if( String::ToInteger( tokens[token_no], iz ) != String::ReturnValue::OK )
				{
					printf( "Error on line %d: '%s'\n", line_no, line );
					printf( "Unable to convert z image flag\n" );
					return -1;
				}
				token_no++;

				config.AddAtom( serial, mol, type, q, x+L[0]*ix, y+L[1]*iy, z+L[2]*iz );
			}
			else*/
			{
				config.AddAtom( serial, mol, type, q, x, y, z );
			}

			continue;
		}
		
		if( state == "Bonds" )
		{
			if( tokens.size() < 4 )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Too few tokens\n" );
				return -1;
			}

			int serial, type, i, j;

			if( String::ToInteger( tokens[0], serial ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert serial token\n" );
				return -1;
			}

			if( String::ToInteger( tokens[1], type ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert type token\n" );
				return -1;
			}

			if( String::ToInteger( tokens[2], i ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert i token\n" );
				return -1;
			}

			if( String::ToInteger( tokens[3], j ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert j token\n" );
				return -1;
			}

			config.AddBond( type, i, j );

			continue;
		}

		if( state == "Angles" )
		{
			if( tokens.size() < 5 )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Too few tokens\n" );
				return -1;
			}

			int serial, type, i, j, k;

			if( String::ToInteger( tokens[0], serial ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert serial token\n" );
				return -1;
			}

			if( String::ToInteger( tokens[1], type ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert type token\n" );
				return -1;
			}

			if( String::ToInteger( tokens[2], i ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert i token\n" );
				return -1;
			}

			if( String::ToInteger( tokens[3], j ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert j token\n" );
				return -1;
			}

			if( String::ToInteger( tokens[4], k ) != String::ReturnValue::OK )
			{
				printf( "Error on line %d: '%s'\n", line_no, line );
				printf( "Unable to convert k token\n" );
				return -1;
			}

			config.AddAngle( type, i, j, k );

			continue;
		}

		if( state == "Bond Coeffs" )
		{
			config.bond_coeffs.push_back( line );
			continue;
		}

		if( state == "Angle Coeffs" )
		{
			config.angle_coeffs.push_back( line );
			continue;
		}
	}

	return 0;			
}

int SaveData( FILE* f, const Config& config, bool has_q )
{
	if( f == nullptr ) return -1;

	const auto& bounds = config.bounds;
	const auto& mass = config.mass;
	const auto& molecules = config.molecules;
	const auto& bond_coeffs = config.bond_coeffs;
	const auto& angle_coeffs = config.angle_coeffs;

	size_t N_atoms, N_bonds, N_angles, serial;
	std::map< int, int > bond_types, angle_types;

	fprintf( f, "Auotogenerated LAMMPS data file.\n" );

	fprintf( f, "\n" );

	//
	// Count some system components.
	//
	{
		N_atoms = N_bonds = N_angles = 0;
		bond_types.clear();
		angle_types.clear();
		for( const auto& m : molecules )
		{
			N_atoms += m.atoms.size();
			
			for( const auto& b : m.bonds ) bond_types[ b.type ] = 1;
			N_bonds += m.bonds.size();

			for( const auto& a : m.angles ) angle_types[ a.type ] = 1;
			N_angles += m.angles.size();
		}
	}

	fprintf( f, "%12d atoms\n", (int)N_atoms );
	fprintf( f, "%12d bonds\n", (int)N_bonds );
	fprintf( f, "%12d angles\n", (int)N_angles );

	fprintf( f, "\n" );

	fprintf( f, "%12d atom types\n", (int)config.mass.size()-1 ); // as we used +1 in size for unit-based indexing
	fprintf( f, "%12d bond types\n", (int)bond_types.size() );
	fprintf( f, "%12d angle types\n", (int)angle_types.size() );

	fprintf( f, "\n" );

	fprintf( f, "  %g   %g   xlo xhi\n", config.bounds.minx, bounds.maxx );
	fprintf( f, "  %g   %g   ylo yhi\n", config.bounds.miny, bounds.maxy );
	fprintf( f, "  %g   %g   zlo zhi\n", config.bounds.minz, bounds.maxz );

	fprintf( f, "\n" );

	fprintf( f, "Masses\n" );

	fprintf( f, "\n" );

	// Note - we start at index 1!
	for( size_t i=1; i<mass.size(); i++ ) fprintf( f, "%12d  %12.3f\n", (int)i, mass[i] );

	fprintf( f, "\n" );

	fprintf( f, "Atoms\n" );

	fprintf( f, "\n" );

	serial = 0;
	for( size_t i=0, max_i=molecules.size(); i<max_i; i++ )
	{
		const auto& m = molecules[i];
		for( size_t j=0, max_j=m.atoms.size(); j<max_j; j++ )
		{
			const auto& a = m.atoms[j];

			if( has_q )
			{
				fprintf( f, "%12d  %12d  %12d  %12.3f  %12.3f  %12.3f  %12.3f\n",
					(int)serial+1,
					(int)i+1,
					a.type,
					a.q,
					a.x, a.y, a.z );
			}
			else
			{
				fprintf( f, "%12d  %12d  %12d  %12.3f  %12.3f  %12.3f\n",
					(int)serial+1,
					(int)i+1,
					a.type,
					a.x, a.y, a.z );
			}
			serial++;
		}
	}

	//
	// Okay, careful now: we may have added multiple configs together,
	// so we can't rely on the raw "min_serial" variable in each molecule to be
	// valid. We should generate a new one as we go, but modify the topology
	// indices to reflect this!
	//
	serial = 0;
	std::vector<size_t> mol_offsets;
	for( auto& m : molecules )
	{
		mol_offsets.push_back( serial );
		serial += m.atoms.size();
	}

	fprintf( f, "\n" );

	fprintf( f, "Bonds\n" );

	fprintf( f, "\n" );

	serial = 0;
	for( size_t i=0, N=molecules.size(); i<N; i++ )
	{
		const auto& m = molecules[i];
		size_t mol_start = mol_offsets[i];
		for( const auto& b : m.bonds )
		{
			fprintf( f, "%12d  %12d  %12d  %12d\n",
				(int)serial+1,
				b.type,
				(int)mol_start + (b.i-m.min_serial) + 1,
				(int)mol_start + (b.j-m.min_serial) + 1 );
			serial++;
		}
	}


	fprintf( f, "\n" );

	fprintf( f, "Angles\n" );

	fprintf( f, "\n" );

	serial = 0;
	for( size_t i=0, N=molecules.size(); i<N; i++ )
	{
		const auto& m = molecules[i];
		size_t mol_start = mol_offsets[i];
		for( const auto& a : m.angles )
		{
			fprintf( f, "%12d  %12d  %12d  %12d  %12d\n",
				(int)serial+1,
				a.type,
				(int)mol_start + (a.i-m.min_serial) + 1,
				(int)mol_start + (a.j-m.min_serial) + 1,
				(int)mol_start + (a.k-m.min_serial) + 1 );
			serial++;
		}
	}

	if( bond_coeffs.size() > 0 )
	{
		fprintf( f, "\n" );
		fprintf( f, "Bond Coeffs\n" );
		fprintf( f, "\n" );
		for( size_t i=0, max_i=bond_coeffs.size(); i<max_i; i++ )
		{
			fprintf( f, "%s\n", bond_coeffs[i].c_str() );
		}
	}

	if( angle_coeffs.size() > 0 )
	{
		fprintf( f, "\n" );
		fprintf( f, "Angle Coeffs\n" );
		fprintf( f, "\n" );
		for( size_t i=0, max_i=angle_coeffs.size(); i<max_i; i++ )
		{
			fprintf( f, "%s\n", angle_coeffs[i].c_str() );
		}
	}

	return 1;
}
//
// Does pretty much what you'd expect!
//
int SaveXYZ( FILE* f, const Config& config, const char *heading )
{
	static const char* abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	static size_t N_abc = strlen( abc );

	if( f == nullptr ) return -1;

	size_t N = 0;
	for( const auto& m : config.molecules ) N += m.atoms.size();

	fprintf( f, "%d\n", (int)N );
	fprintf( f, "%s\n", heading );

	for( const auto& m : config.molecules )
	{
		for( const auto& a : m.atoms ) fprintf( f, "%c %g %g %g\n", abc[(a.type-1)%N_abc], a.x, a.y, a.z );
	}
	return 1;
}



// Basic loading of trajectory frame
int LoadTrajectory( FILE *f, Config& config, int* timestep )
{
	const char *delimiters = " \t\n\r";

	int maxbuf = 1023;
	char buffer[1024];
	
	int ntoks;
	std::vector< std::string > tokens;
	
	fpos_t file_position;

	std::map< std::string, int > token_indices {
		{ "id", -1 }, { "type", -1 }, { "mol", -1 },
		{ "x",  -1 }, { "y",    -1 }, { "z",   -1 } };

	int number_of_atoms = 0; // internal check
	
	if( f == NULL )
	{
		fprintf( stderr, "LAMMPS::%s(): file pointer is NULL!\n", __func__ );
		return -1;
	}
	
	//
	// Only clear the molecules - leave any metadata intact!
	//
	config.ClearMolecules();

	while( fgets( buffer, maxbuf, f ) != NULL )
	{
		if( (ntoks = String::Tokenize( buffer, tokens, delimiters )) < 1 || tokens[0] != "ITEM:" ) continue;

		//
		// Trajectory frame timestep number
		//
		if( tokens[1] == "TIMESTEP" )
		{
			int ts;

			if( fgets( buffer, maxbuf, f ) == NULL ) continue;
			if( String::ToInteger( buffer, ts ) != String::ReturnValue::OK )
			{
				fprintf( stderr, "LAMMPS::%s(): unable to convert timestep '%s' into an integer.\n", __func__, buffer );
				return -1;
			}
			if( timestep != nullptr ) *timestep = ts;
		}
		//
		// Simulation cell
		//
		else if( tokens[1] == "BOX" )
		{
			const char* xyz = "xyz";
			double* dims[] = {
				&config.bounds.minx, &config.bounds.maxx,
				&config.bounds.miny, &config.bounds.maxy,
				&config.bounds.minz, &config.bounds.maxz };

			// expects ITEM: BOX BOUNDS xx yy zz
			if( tokens.size() != 6 )
			{
				fprintf( stderr, "LAMMPS::%s(): bad box line: '%s'\n", __func__, buffer  );
				return -1;
			}
			
			for( int axis=0; axis<3; axis++ ) config.PBC[axis] = ( tokens[3+axis][0] == 'p' ) ? (1) : (0);

			for( int axis=0; axis<3; axis++ )
			{
				if( fgets( buffer, maxbuf, f ) == NULL ) continue;
				if( (ntoks = String::Tokenize( buffer, tokens, delimiters)) < 2 ) continue;

				if( String::ToReal( tokens[0], *dims[axis*2 +0] ) != String::ReturnValue::OK )
				{
					fprintf( stderr, "LAMMPS::%s(): unable to convert min %c value '%s' into a number.\n", __func__, xyz[axis], tokens[0].c_str() );
					return -1;
				}
				if( String::ToReal( tokens[1], *dims[axis*2 +1] ) != String::ReturnValue::OK )
				{
					fprintf( stderr, "LAMMPS::%s(): unable to convert max %c value '%s' into a number.\n", __func__, xyz[axis], tokens[1].c_str() );
					return -1;
				}
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
				fprintf( stderr, "LAMMPS::%s(): unable to convert number of atoms '%s' into a number.\n", __func__, buffer );
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
			/*
			for( it=token_indices.begin(); it!=token_indices.end(); it++ )
			{
				if( it->second >= 0 ) continue;
				fprintf( stderr, "LAMMPS::%s(): unable to find '%s' info in frame.\n", __func__, it->first.c_str() );
				return -1;
			}
			*/
			
			int id_i   = token_indices["id"];
			int type_i = token_indices["type"];
			int mol_i  = token_indices["mol"];
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

					int check = 0;
					for( const auto& m : config.molecules ) check += (int)m.atoms.size();

					if( check != number_of_atoms )
					{
						fprintf( stderr, "LAMMPS::%s(): WARNING: odd number of atoms (read %d, expected %d)\n", __func__, check, number_of_atoms );
					}

					return check;
				}
				
				int atom_id=1, mol_id=1, atom_type=1;
				double x=0, y=0, z=0;

				if( (id_i>=0) && String::ToInteger( tokens[id_i], atom_id ) != String::ReturnValue::OK )
				{
					fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert atom id '%s' into a number.\n", __func__, tokens[id_i].c_str() );
					return -1;
				}
				if( (type_i>=0) && String::ToInteger( tokens[type_i], atom_type ) != String::ReturnValue::OK )
				{
					fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert type '%s' into a number.\n", __func__, tokens[type_i].c_str() );
					return -1;
				}
				if( (mol_i>=0) && String::ToInteger( tokens[mol_i], mol_id ) != String::ReturnValue::OK )
				{
					fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert mol id '%s' into a number.\n", __func__, tokens[mol_i].c_str() );
					return -1;
				}

				if( (x_i>=0) && String::ToReal( tokens[x_i], x ) != String::ReturnValue::OK )
				{
					fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert x '%s' into a number.\n", __func__, tokens[x_i].c_str() );
					return -1;
				}
				if( (y_i>=0) && String::ToReal( tokens[y_i], y ) != String::ReturnValue::OK )
				{
					fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert y '%s' into a number.\n", __func__, tokens[y_i].c_str() );
					return -1;
				}
				if( (z_i>=0) && String::ToReal( tokens[z_i], z ) != String::ReturnValue::OK )
				{
					fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert z '%s' into a number.\n", __func__, tokens[z_i].c_str() );
					return -1;
				}

				double Lx = config.bounds.maxx - config.bounds.minx;
				double Ly = config.bounds.maxy - config.bounds.miny;
				double Lz = config.bounds.maxz - config.bounds.minz;

				if( xscale == true ) x = config.bounds.minx + x*Lx;
				if( yscale == true ) y = config.bounds.miny + y*Ly;
				if( zscale == true ) z = config.bounds.minz + z*Lz;

				config.AddAtom( atom_id, mol_id, atom_type, 0.0, x, y, z );
			}
		}
	}

	int check = 0;
	for( const auto& m : config.molecules ) check += (int)m.atoms.size();
	return check;
}		

// Basic saving of trajectory frame
int SaveTrajectory( FILE *f, Config& config, int timestep )
{
	fprintf( f, "ITEM: TIMESTEP\n" );
	fprintf( f, "%d\n", timestep );

	{
		int N = 0;
		for( const auto& m : config.molecules ) N += (int)m.atoms.size();
		fprintf( f, "ITEM: NUMBER OF ATOMS\n" );
		fprintf( f, "%d\n", N );
	}

	fprintf( f, "ITEM: BOX BOUNDS pp pp pp\n" );
	fprintf( f, "%g %g\n", config.bounds.minx, config.bounds.maxx );
	fprintf( f, "%g %g\n", config.bounds.miny, config.bounds.maxy );
	fprintf( f, "%g %g\n", config.bounds.minz, config.bounds.maxz );

	fprintf( f, "ITEM: ATOMS id type mol x y z\n" );
	for( const auto& m : config.molecules )
	{
		for( size_t i=0, N=m.atoms.size(); i<N; i++ )
		{
			const auto&a = m.atoms[i];
			fprintf( f, "%d %d %d %g %g %g\n", a.serial, a.type, m.id, a.x, a.y, a.z );
		}
	}
	return 1;
}		


}
}
