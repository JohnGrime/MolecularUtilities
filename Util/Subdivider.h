/*
	Author: John Grime, The University of Chicago.
	Based on public domain code from AT&T
*/

#if !defined( UTIL_SUBDIVIDER_DEFINED )
#define UTIL_SUBDIVIDER_DEFINED

#include <stdio.h>
#include <stdlib.h>

#include <assert.h>
#include <math.h>
#include <string.h> // for memcpy

#include <vector>

namespace Util
{

class Subdivider
{
	public:

		enum class InitType { Tetrahedron, Octahedron, Icosahedron };

		using Real = double;
	
		int n_vertices, n_faces, n_edges;

		std::vector<Real> vertices;
		std::vector<int> faces;
	
		std::vector<int> start, end, midpoint;

		//
		// Refinement of initial structures
		//
		void Init( InitType init_as )
		{
			n_vertices = n_faces = n_edges = 0;
			
			switch( init_as )
			{
				case InitType::Tetrahedron:
					init_tetrahedron();
				break;

				case InitType::Octahedron:
					init_octahedron();
				break;

				case InitType::Icosahedron:
					init_icosahedron();
				break;

				default:
					printf( "Subdivider::%s(): unknown init type.\n", __func__ );
					exit( -1 );
				break;
			}
		}

		void Subdivide()
		{ 
			int n_vertices_new = n_vertices+2*n_edges; 
			int n_faces_new = 4*n_faces; 
			int edge_walk = 0; 
			
			n_edges = 2*n_vertices + 3*n_faces;

			start.resize( n_edges );
			end.resize( n_edges );
			midpoint.resize( n_edges );

			auto faces_old = faces; // "hidden" allocation & copy.
			
			vertices.resize( 3*n_vertices_new );
			faces.resize( 3*n_faces_new );
			
			n_faces_new = 0; 

			for( int i=0; i<n_faces; i++ )
			{
				int a = faces_old[3*i+0];
				int b = faces_old[3*i+1];
				int c = faces_old[3*i+2];
		        
				int ab_midpoint = search_midpoint( b, a, edge_walk ); 
				int bc_midpoint = search_midpoint( c, b, edge_walk ); 
				int ca_midpoint = search_midpoint( a, c, edge_walk ); 

				faces[3*n_faces_new+0] = a; 
				faces[3*n_faces_new+1] = ab_midpoint; 
				faces[3*n_faces_new+2] = ca_midpoint; 
				n_faces_new++; 

				faces[3*n_faces_new+0] = ca_midpoint; 
				faces[3*n_faces_new+1] = ab_midpoint; 
				faces[3*n_faces_new+2] = bc_midpoint; 
				n_faces_new++; 

				faces[3*n_faces_new+0] = ca_midpoint; 
				faces[3*n_faces_new+1] = bc_midpoint; 
				faces[3*n_faces_new+2] = c; 
				n_faces_new++; 

				faces[3*n_faces_new+0] = ab_midpoint; 
				faces[3*n_faces_new+1] = b; 
				faces[3*n_faces_new+2] = bc_midpoint; 
				n_faces_new++; 
			} 
			n_faces = n_faces_new; 

			//
			// Adjust size of internal storage, as algorithm over-allocates! This is not a 
			// problem if you're using n_vertices to determine length of vertex array, but
			// you can't rely on vertices.size() for the right number of coordinates.
			//
			// Is this due to a potential bug in the original code?
			//
			vertices.resize( n_vertices*3 );
		} 
		

		Subdivider()
		{
			n_vertices = n_faces = n_edges = 0;
			Init( InitType::Tetrahedron );
		}

	protected:

		//
		// Initialise from specific structure types
		//
		void init_tetrahedron()
		{
			const Real sqrt3 = 1.0 / sqrt(3.0);
			
			vertices = {
				 sqrt3,  sqrt3,  sqrt3,
				-sqrt3, -sqrt3,  sqrt3,
				-sqrt3,  sqrt3, -sqrt3,
				 sqrt3, -sqrt3, -sqrt3 }; 
			
			faces = { 0,2,1, 0,1,3, 2,3,1, 3,2,0 };

			n_vertices = vertices.size()/3;
			n_faces    = faces.size()/3;
			n_edges    = 6;
		}

		void init_octahedron()
		{
			vertices = {
				 0.0,  0.0, -1.0,
				 1.0,  0.0,  0.0,
				 0.0, -1.0,  0.0,
				-1.0,  0.0,  0.0,
				 0.0,  1.0,  0.0,
				 0.0,  0.0,  1.0 };

			faces = { 0,1,2, 0,2,3, 0,3,4, 0,4,1, 5,2,1, 5,3,2, 5,4,3, 5,1,4 };

			n_vertices = vertices.size()/3;
			n_faces    = faces.size()/3;
			n_edges    = 12;
		}

		void init_icosahedron()
		{
			const Real t   = (1.0+sqrt(5))/2;
			const Real tau = t/sqrt(1+t*t);
			const Real one = 1.0/sqrt(1+t*t);

			vertices = {
				 tau,   one,   0.0,
				-tau,   one,   0.0,
				-tau,  -one,   0.0,
				 tau,  -one,   0.0,
				 one,   0.0 ,  tau,
				 one,   0.0 , -tau,
				-one,   0.0 , -tau,
				-one,   0.0 ,  tau,
				 0.0 ,  tau,   one,
				 0.0 , -tau,   one,
				 0.0 , -tau,  -one,
				 0.0 ,  tau,  -one };
				
			faces = {
				4,  8,  7,
				4,  7,  9,
				5,  6,  11,
				5,  10, 6,
				0,  4,  3,
				0,  3,  5,
				2,  7,  1,
				2,  1,  6,
				8,  0,  11,
				8,  11, 1,
				9,  10, 3,
				9,  2,  10,
				8,  4,  0,
				11, 0,  5,
				4,  9,  3,
				5,  3,  10,
				7,  8,  1,
				6,  1,  11,
				7,  2,  9,
				6,  10, 2 };

			n_vertices = vertices.size()/3;
			n_faces    = faces.size()/3;
			n_edges    = 30;
		}

		int search_midpoint( int index_start, int index_end, int& edge_walk )
		{ 
			int res;
			float length;
			
			for( int i=0; i<edge_walk; i++ )
			{
				if( (start[i] == index_start && end[i] == index_end) || (start[i] == index_end && end[i] == index_start) )
				{
					res = midpoint[i];

					// update the arrays
					start[i]    = start[edge_walk-1];
					end[i]      = end[edge_walk-1];
					midpoint[i] = midpoint[edge_walk-1];
					edge_walk--;

					return res; 
				}
			}

			// vertex not in the list, so we add it
			start[edge_walk] = index_start;
			end[edge_walk] = index_end; 
			midpoint[edge_walk] = n_vertices; 

			// create new vertex
			vertices[3*n_vertices+0] = (vertices[3*index_start+0] + vertices[3*index_end+0]) / 2.0;
			vertices[3*n_vertices+1] = (vertices[3*index_start+1] + vertices[3*index_end+1]) / 2.0;
			vertices[3*n_vertices+2] = (vertices[3*index_start+2] + vertices[3*index_end+2]) / 2.0;

			// normalize the new vertex
			length = sqrt(
				vertices[3*n_vertices+0] * vertices[3*n_vertices+0] +
				vertices[3*n_vertices+1] * vertices[3*n_vertices+1] +
				vertices[3*n_vertices+2] * vertices[3*n_vertices+2] );
			length = 1.0/length;
			vertices[3*n_vertices+0] *= length;
			vertices[3*n_vertices+1] *= length;
			vertices[3*n_vertices+2] *= length;

			n_vertices++;
			edge_walk++;
			return midpoint[edge_walk-1];
		} 
};

}

#endif
