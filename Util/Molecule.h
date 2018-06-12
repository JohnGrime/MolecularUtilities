/*
	Author: John Grime, The University of Chicago.
*/

#if !defined( UTIL_MOLECULE_DEFINED )
#define UTIL_MOLECULE_DEFINED

#include <assert.h>
#include <stdio.h>

#include <string>
#include <vector>
#include <map>

#include <algorithm>

#include "String.h"
#include "Geometry.h"

//
// Utility routines for handling arbitrary molecular data, including filtering and IO.
// Can be slower and require more memory vs code for specific data (e.g. LAMMPS.h/cpp).
//

namespace Util
{
namespace Molecule
{

//
// Per-atom and metadata map type
//
using AttributeMap = std::map< std::string, std::string >;

//
// We're using xyz values frequently, and we assume they're guaranteed to exist,
// so store as raw numbers. Other per-atom data stored in attribute map, for
// arbitrary information. Slow, but simple.
//
struct Atom
{
		double x, y, z;
		AttributeMap attr;
		
		Atom()
		{
			x = y = z = 0.0;
			attr.clear();
		}
		Atom( double x_, double y_, double z_, const AttributeMap *def_attr = nullptr )
		{
			x = x_;
			y = y_;
			z = z_;
			if( def_attr != nullptr ) attr = *def_attr;
		}
};
using Molecule = std::vector<Atom>;
using Molecules = std::vector<Molecule>;

class Coords
{
	public:
	
	//
	// Flatten one or more Molecule(s) into xyz array.
	//
	// For single molecules, we specify whether to clear the xyz vector (allows progressive generation).
	//
	template<typename T>
	static int Get( const Molecule &mol, std::vector<T> &xyz, bool clear_xyz = true )
	{
		int N = 0;
		
		if( clear_xyz == true ) xyz.clear();
		
		for( const auto& a : mol )
		{
			xyz.push_back( a.x );
			xyz.push_back( a.y );
			xyz.push_back( a.z );		
			N++;
		}
		return N;
	}

	template<typename T>
	static int Get( const Molecules &mols, std::vector<T> &xyz )
	{
		int N = 0;
		
		xyz.clear(); // Assume clearing xyz vector to start
		for( const auto& m : mols )
		{
			N += Get( m, xyz, false ); // Note - no clearing of the xyz vector for each call!
		}
		return N;
	}

	template<typename T>
	static int Get( const Molecules &mols, const std::vector<int> &mol_indices, std::vector<T> &xyz )
	{
		int N = 0;
		
		xyz.clear(); // Assume clearing xyz vector to start
		for( const auto index : mol_indices )
		{
			N += Get( mols[index], xyz, false ); // Note - no clearing of xyz vector for each call!
		}
		return N;
	}

	//
	// Expand xyz arrays back into Atom coords
	//
	template<typename T>
	static int Set( Molecule &mol, const std::vector<T> &xyz )
	{
		assert( xyz.size()%3 == 0 );
		assert( xyz.size()/3 == mol.size() );
	
		for( size_t ai=0, max_ai=mol.size(); ai<max_ai; ai++ )
		{
			Atom &a = mol[ai];

			a.x = xyz[ai*3 +0];
			a.y = xyz[ai*3 +1];
			a.z = xyz[ai*3 +2];
		}
		return (int)xyz.size()/3;
	}

	template<typename T>
	static int Set( Molecules &mols, const std::vector<T> &xyz )
	{
		size_t offset = 0;
	
		assert( xyz.size()%3 == 0 );
		for( const auto& m : mols )
		{
			offset += m.size();
		}
		assert( xyz.size()/3 == offset );
	
		offset = 0;
		for( auto& mol : mols )
		{
			for( auto& a : mol )
			{
				a.x = xyz[offset*3 +0];
				a.y = xyz[offset*3 +1];
				a.z = xyz[offset*3 +2];
				offset++;
			}
		}
		return (int)xyz.size()/3;
	}

	template<typename T>
	static int Set( Molecules &mols, const std::vector<int> &mol_indices, const std::vector<T> &xyz )
	{
		size_t offset = 0;
	
		assert( xyz.size()%3 == 0 );
		for( const auto index : mol_indices )
		{
			offset += mols[index].size();
		}
		assert( xyz.size()/3 == offset );
	
		offset = 0;
		for( const auto index : mol_indices )
		{
			auto& mol = mols[index];
			for( auto& a : mol )
			{
				a.x = xyz[offset*3 +0];
				a.y = xyz[offset*3 +1];
				a.z = xyz[offset*3 +2];
				offset++;
			}
		}
		return (int)xyz.size()/3;
	}
	
	//
	// Translate
	//
	template<typename T>
	static void Translate( Molecule &mol, T tx, T ty, T tz )
	{
		size_t N = mol.size();
	
		if( N == 0 ) return;

		for( auto& a : mol )
		{
			a.x += tx;
			a.y += ty;
			a.z += tz;
		}
	}
	template<typename T>
	static void Translate( Molecules &mols, T tx, T ty, T tz )
	{
		size_t N_mols = mols.size();
	
		if( N_mols == 0 ) return;
		
		for( auto& mol : mols )
		{
			Translate( mol, tx, ty, tz );
		}
	}
	
	//
	// Rotate around point as per GeometryUtil : if point == nullptr, rotation axis passes through origin
	//
	template<typename T>
	static void Rotate( const Molecule &mol, const T *point, const T *axis, T theta )
	{
		size_t N = mol.size();
		
		if( N == 0 ) return;
		
		Geometry::Rotate( N, &mol[0], point, axis, theta );
	}

	template<typename T>
	static void Rotate( const Molecules &mols, const T *point, const T *axis, T theta )
	{
		size_t N = mols.size();
	
		if( N == 0 ) return;
		
		for( size_t i=0; i<N; i++ )
		{
			Rotate( mols[i], point, axis, theta );
		}
	}

	//
	// Get centre of geometry.
	//
	template<typename T>
	static void GetCOG( const Molecule &mol, T &cx, T &cy, T &cz )
	{
		size_t N = mol.size();
		
		cx = cy = cz = 0.0;
	
		if( N == 0 ) return;
		
		for( auto& a : mol )
		{
			cx += a.x;
			cy += a.y;
			cz += a.z;
		}
		
		cx /= N;
		cy /= N;
		cz /= N;
	}

	template<typename T>
	static void GetCOG( const Molecules &mols, T &cx, T &cy, T &cz )
	{
		size_t N = 0;
		
		cx = cy = cz = 0.0;
	
		if( mols.size() == 0 ) return;
		
		for( const auto& mol : mols )
		{
			for( const auto& a : mol )
			{
				cx += a.x;
				cy += a.y;
				cz += a.z;
				N++;
			}
		}
		
		if( N < 1 ) return; // avoid divide by zero.
		
		cx /= N;
		cy /= N;
		cz /= N;
	}

	//
	// Reset centre of geometry. Returns maximum distance from COG of any particle.
	//
	static double ResetCOG( Molecule &mol )
	{
		double cx, cy, cz;
		double max_dr2 = -1.0;
		
		if( mol.size() == 0 ) return 0.0;
		
		GetCOG( mol, cx, cy, cz );
		
		for( auto& a : mol )
		{
			a.x -= cx;
			a.y -= cy;
			a.z -= cz;

			double dr2 = (a.x*a.x) + (a.y*a.y) + (a.z*a.z);
			if( dr2 > max_dr2 ) max_dr2 = dr2;
		}

		return sqrt( max_dr2 );
	}
	static double ResetCOG( Molecules &mols )
	{
		double cx, cy, cz;
		double max_dr2 = -1.0;
		
		if( mols.size() == 0 ) return 0.0;
		
		GetCOG( mols, cx, cy, cz );
		
		for( auto& mol : mols )
		{
			for( auto& a : mol )			
			{
				a.x -= cx;
				a.y -= cy;
				a.z -= cz;

				double dr2 = (a.x*a.x) + (a.y*a.y) + (a.z*a.z);
				if( dr2 > max_dr2 ) max_dr2 = dr2;
			}
		}

		return sqrt( max_dr2 );
	}

	//
	// Get bounding box of one or more Molecule(s).
	//
	// _bounding_box() assumes you have set the initial min/max values
	// in advance as appropriate.
	//
	template<typename T>
	static int _bounding_box( const Molecule &mol,
		T& minx, T &maxx,
		T& miny, T &maxy,
		T& minz, T &maxz )
	{
		int N = 0;

		for( const auto& a : mol )
		{
			minx = std::min( minx, a.x );
			maxx = std::max( maxx, a.x );

			miny = std::min( miny, a.y );
			maxy = std::max( maxy, a.y );

			minz = std::min( minz, a.z );
			maxz = std::max( maxz, a.z );

			N++;
		}

		return N;
	}

	template<typename T>
	static int BoundingBox( const Molecule &mol,
		T& minx, T &maxx,
		T& miny, T &maxy,
		T& minz, T &maxz )
	{
		minx = maxx = miny = maxy = minz = maxz = 0;

		if( mol.size() < 1 ) return 0;

		minx = maxx = mol[0].x;
		miny = maxy = mol[0].y;
		minz = maxz = mol[0].z;

		return _bounding_box( mol, minx,maxx, miny,maxy, minz,maxz );
	}

	template<typename T>
	static int BoundingBox( const Molecules &mols,
		T& minx, T &maxx,
		T& miny, T &maxy,
		T& minz, T &maxz )
	{
		int N = 0;

		minx = maxx = miny = maxy = minz = maxz = 0;

		for( const auto& mol : mols )
		{
			if( mol.size() < 1 ) continue;

			if( N == 0 )
			{
				minx = maxx = mol[0].x;
				miny = maxy = mol[0].y;
				minz = maxz = mol[0].z;
			}

			N += _bounding_box( mol, minx,maxx, miny,maxy, minz,maxz );
		}

		return N;
	}

	template<typename T>
	static int BoundingBox( const Molecules &mols, const std::vector<int> &mol_indices,
		T& minx, T &maxx,
		T& miny, T &maxy,
		T& minz, T &maxz )
	{
		int N = 0;

		minx = maxx = miny = maxy = minz = maxz = 0;

		for( const auto index : mol_indices )
		{
			if( (index<0) || (index>=(int)mols.size()) ) continue;

			const auto& mol = mols[index];
			if( mol.size() < 1 ) continue;

			if( N == 0 )
			{
				minx = maxx = mol[0].x;
				miny = maxy = mol[0].y;
				minz = maxz = mol[0].z;
			}

			N += _bounding_box( mol, minx,maxx, miny,maxy, minz,maxz );
		}

		return N;
	}

};


/*
	Molecule filters allow to to generate reduced sets of atoms on the basis of attribute keys.
	
	There are two ways to do this:
	
	1. Filter on literal values (extract atoms that match explicitly specified attributes)
	2. Filter on attribute equality (extract atoms whose attributes match)
	
	AddFilter() lets you set up the filters; if you specify a key alone (with no explicit values), any atom will pass that
	filter provided it has an attribute of that name (i.e. the value of that attribute is irrelevant). However, when
	we filter 2 atoms against one another, not only must both atoms have the same attribute name present but the values
	of this attribute must also match across the two atoms.
	
	If you provide a set of explicit values to the filter, atoms must have not only attributes of that name, but also the value
	of those atributes must be present in the explicitly specified values. For filtering 2 atoms against one another, both atoms must
	have an attribute of that name, with equal values that are present in the explicitly specified filter values.
	
	Eg:
	
	Attribute key = "blah", empty value:
	
		single-atom filter: atom must have an attribute of name "blah", but the value can be anything.
		dual-atom filter: both atoms must have attributes of name "blah", and their values must be equal.

	Attribute key = "blah", values = [ "wibble1", "wibble2" ]:

		single-atom filter: atom must have an attribute of name "blah", with value in [ "wibble1", "wibble2" ].
		dual-atom filter: both atoms must have an attribute of name "blah", with identical values in [ "wibble1", "wibble2" ].
	
	The Filter() methods use temporary storage before copying to the output molecule, so "in" and "out" can be the same.
*/
class AttributeFilter
{
	public:
				
		//
		// filter_values : acceptable values for a given filter key (empty map == any value, provided matches across atoms)
		// filters : associates filter keys with acceptable values (see above)
		//
		using FilterValuesMap = std::map< std::string, int >;
		using FilterMap = std::map< std::string, FilterValuesMap >;

		FilterMap filters;

		void AddFilter( const char *key, const char *value = nullptr );
		void AddFilter( const char *filter_str, const char *keyval_sep, const char *val_sep, const char *range_sep );

		//
		// Return number of atom(s) that passed the filter. For single atom filters,
		// value of 1 therefore indicates the specified atom(s) passed the filter.
		//
		int Filter( const Atom &in ) const;
		int Filter( const Atom &in1, const Atom &in2 ) const;

		int Filter( const Molecule &in, Molecule &out ) const;
		int Filter( const Molecule &in1, const Molecule &in2, Molecule &out1, Molecule &out2 ) const;
		
		int Filter( const Molecules &in, Molecules &out ) const;
		int Filter( const Molecules &in1, const Molecules &in2, Molecules &out1, Molecules &out2 ) const;

		void Print( FILE *f = nullptr ) const;
};


/*
	Collection of Molecules, and some metadata about the system; this allows us to store not only
	arbitrary atom data (assuming we can parse it into attributes etc), but other data depending on
	the atom data source (eg sim cell dimensions, timestep etc if from a trajectory).
*/
class MolecularSystem
{
	public:
		
		AttributeMap metadata;
		Molecules molecules;
		
		// if n_mols and n_atoms control printing of first n_atoms of first n_mols for info.
		void Print( FILE *f, int n_mols = 3, int n_atoms = 3 ) const;
};

/*
	Load data from a PDB file into a MoleculeSet, and print from a MoleculeSet into PDB format.
	
	To do: add default PDB column values. As it stands, the Print() routine will put empty columns
	in for any missing attributes - it would be nice to specify defaults for missing attributes.
*/
class PDB
{
	public:
		
		static int Load( FILE *f, MolecularSystem &ms, int n_mols = -1 );

		static int Print( FILE *f, const Atom &a ); // USES RAW x|y|z, NOT DATA IN ATTRIBUTE MAPS!
		static int Print( FILE *f, const Molecule &mol ); // USES RAW x|y|z, NOT DATA IN ATTRIBUTE MAPS!
		static int Print( FILE *f, const Molecules &mols ); // USES RAW x|y|z, NOT DATA IN ATTRIBUTE MAPS!
		static int Print( FILE *f, const MolecularSystem &ms ); // USES RAW x|y|z, NOT DATA IN ATTRIBUTE MAPS!
};

/*
	Load data from a LAMMPS trajectory file into a MoleculeSet, and print MoleculeSet into LAMMPS trajectory.
*/
class LAMMPS
{
	public:
		
		static int Load( FILE *f, MolecularSystem &ms );

		static int Print( FILE *f, const MolecularSystem &ms ); // USES RAW x|y|z, NOT DATA IN ATTRIBUTE MAPS!
};


}
}

#endif
