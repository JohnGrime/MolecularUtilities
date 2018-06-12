/*
	Author: John Grime, The University of Chicago.
*/

#if !defined( UTIL_SUPERPOSE_DEFINED )
#define UTIL_SUPERPOSE_DEFINED

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <vector>

namespace Util
{

//
// Stand-alone implementation of the Kabsch algorithm.
//
class Superposer
{
	protected:
		
		typedef double T_internal;
		
		//
		// Converted from the gnarly original fortran u3b routine. I have no idea who the original author was. Gandalf perhaps? Merlin?
		//
		// N : number of atoms
		// w : per-atom weights
		// x, y : contiguous xyz values for point sets x and y (i.e. lengths are 3*N)
		// mode: 0 == calculate RMS alone, 1 == calculate RMS and rotation matrix/translation vector
		// rmsd : RMS deviation between x and y after superposition
		// M, T : rotation matrix & translation vector
		// ier : error flag
		//
		template< typename T_external >
		static int u3b( int N, const T_external *weights, const T_external *x, const T_external *y, const int mode, T_external *rmsd, T_external M[3][3], T_external T[3] )
		{
			int i, j, k, l, m1, m;
			int ip[9] = { 0, 1, 3, 1, 2, 4, 3, 4, 5 }; // NOTE - changed from unit to zero-based indices!
			int ip1201[4] = { 1, 2, 0, 1 }; // again, switched to zero based. Name changed to reflect this.
			T_internal r[3][3], xc[3], yc[3], a[3][3], b[3][3], e[3], rr[6], ss[6];
			T_internal sum_weights, e0, d, spur, det, cof, h, g, cth, sth, sqrth, p, sigma;

			const T_internal sqrt3 = 1.73205080756888e+00;
			const T_internal tol = 1.0e-2;

			//
			// Clear initial data ( M = identity )
			//
			sum_weights = 0.0;
			*rmsd = 0.0;
			e0 = 0.0;
			for( i=0; i<3; i++ )
			{
				xc[i] = 0.0;
				yc[i] = 0.0;
				T[i] = 0.0;
				for( j=0; j<3; j++ )
				{
					if( i == j )
					{
						M[i][j] = 1.0;
						a[i][j] = 1.0;
					}
					else
					{
						M[i][j] = 0.0;
						a[i][j] = 0.0;
					}
					r[i][j] = 0.0;

					a[i][j] = 0.0;
					b[i][j] = 0.0;
				}
			}

			if( N < 1 ) return -1;

			//
			// Get centroid of the point sets
			//
			for( m=0; m<N; m++ )
			{
				T_internal w = (weights == NULL) ? 1.0 : weights[m];
				
				if( w < 0.0 ) return -2;
				sum_weights += w;
				
				for( i=0; i<3; i++ )
				{
					xc[i] += w*x[m*3+i];
					yc[i] += w*y[m*3+i];
				}
			}
			for( i=0; i<3; i++ )
			{
				xc[i] /= sum_weights;
				yc[i] /= sum_weights;
			}

			//
			// Correlation matrix of two point sets
			//
			for( m=0; m<N; m++ )
			{
				T_internal w = (weights == NULL) ? 1.0 : weights[m];
				
				for( i=0; i<3; i++ )
				{
					e0 += w * ( (x[m*3+i]-xc[i])*(x[m*3+i]-xc[i]) + (y[m*3+i]-yc[i])*(y[m*3+i]-yc[i]) );
					d = w * ( y[m*3+i] - yc[i] );
					for( j=0; j<3; j++ )
					{
						r[i][j] += d*( x[m*3+j] - xc[j] );
					}
				}
			}

			//
			// Determinant of correlation matrix
			//
			det = r[0][0] * ( (r[1][1]*r[2][2]) - (r[1][2]*r[2][1]) )
			    - r[0][1] * ( (r[1][0]*r[2][2]) - (r[1][2]*r[2][0]) )
			    + r[0][2] * ( (r[1][0]*r[2][1]) - (r[1][1]*r[2][0]) );
			sigma = det;

			//
			// rr is correlation matrix multiplied by itself; as any such matrix multiplied by itself is guaranteed symmetrical,
			// we only store the upper triangular.
			//
			m = 0;
			for( j=0; j<3; j++ )
			{
				for( i=0; i<=j; i++ )
				{
					rr[m] = r[0][i]*r[0][j] + r[1][i]*r[1][j] + r[2][i]*r[2][j];
					m++;
				}
			}

			//
			// Calculate some useful info, and get the eigenvalues by solving the cubic equation:
			//
			// Y**3 - 3HY + 2G = 0
			//
			// via setting X=Y+SPUR
			//
			spur = (rr[0]+rr[2]+rr[5]) / 3.0;
			cof = (((((rr[2]*rr[5] - rr[4]*rr[4]) + rr[0]*rr[5]) - rr[3]*rr[3]) + rr[0]*rr[2]) - rr[1]*rr[1]) / 3.0;
			det = det*det;

			for( i=0; i<3; i++ )
			{
				e[i] = spur;
			}
			if( spur <= 0.0 ) goto label_40; // JUMP POINT HERE => build translation vector and get rms
			d = spur*spur;
			h = d - cof;
			g = (spur*cof - det)/2.0 - spur*h;
			if( h <= 0.0 )
			{
				if( mode == 0 ) goto label_50; // JUMP POINT HERE => get rms
				else goto label_30; // JUMP POINT HERE => skip part of the generation of matrix a
			}
			sqrth = sqrt( h );
			d = h*h*h - g*g;
			if( d < 0.0 ) d = 0.0;
			d = atan2( sqrt(d), -g ) / 3.0;
			cth = sqrth * cos( d );
			sth = sqrth * sqrt3 * sin(d);
			e[0] = (spur + cth) + cth;
			e[1] = (spur - cth) + sth;
			e[2] = (spur - cth) - sth;

			//
			// Do we need to calculate the rotation and translation vector etc, or can we
			// make do with only the eigenvalues? If so, skip the time consuming eigenvector
			// generations, and jump right to the final rms code.
			//
			if( mode == 0 ) goto label_50; // JUMP POINT HERE => get rms


			//
			// Eigenvector time!
			//
			
			
			for( l=0; l<=3; l+=2 )
			{
				d = e[l];
				ss[0] = (d-rr[2]) * (d-rr[5])  - rr[4]*rr[4];
				ss[1] = (d-rr[5]) * rr[1]      + rr[3]*rr[4];
				ss[2] = (d-rr[0]) * (d-rr[5])  - rr[3]*rr[3];
				ss[3] = (d-rr[2]) * rr[3]      + rr[1]*rr[4];
				ss[4] = (d-rr[0]) * rr[4]      + rr[1]*rr[3];
				ss[5] = (d-rr[0]) * (d-rr[2])  - rr[1]*rr[1];

				if( fabs(ss[0]) >= fabs(ss[2]) )
				{
					j = 1;
					if( fabs(ss[0]) < fabs(ss[5]) ) j = 3;
				}
				else if( fabs(ss[2]) >= fabs(ss[5]) )
				{
					j = 2;
				}
				else
				{
					j = 3;
				}
				d = 0.0;
				j = 3*(j-1);
				for( i=0; i<3; i++ )
				{
					k = ip[i+j];
					a[i][l] = ss[k];
					d += ss[k]*ss[k];
				}
				if( d > 0.0 ) d = 1.0 / sqrt( d );
				for( i=0; i<3; i++ ) a[i][l] *= d;
			}

			//
			// Construct matrix a
			//
			d = a[0][0]*a[0][2] + a[1][0]*a[1][2] + a[2][0]*a[2][2];
			if( (e[0]-e[1]) > (e[1]-e[2]) )
			{
				m1 = 2;
				m = 0;
			}
			else
			{
				m1 = 0;
				m = 2;
			}

			p = 0.0;
			for( i=0; i<3; i++ )
			{
				a[i][m1] -= d*a[i][m];
				p += a[i][m1]*a[i][m1];
			}
			if( p <= tol )	
			{
				p = 1.0;
				for( i=0; i<3; i++ )
				{
					if( p < fabs( a[i][m] ) ) continue; // was "cycle" in fortran code
					p = fabs( a[i][m] );
					j = i;
				}
				k = ip1201[j];
				l = ip1201[j+1];
				p = sqrt( a[k][m]*a[k][m] + a[l][m]*a[l][m] );
				if( p <= tol ) goto label_40; // JUMP POINT HERE => build translation vector and get rms.
				a[j][m1] = 0.0;
				a[k][m1] = -a[l][m]/p;
				a[l][m1] =  a[k][m]/p;
			}
			else
			{
				p = 1.0 / sqrt( p );
				for( i=0; i<3; i++ ) a[i][m1] *= p;
			}


			a[0][1] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
			a[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
			a[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];

			label_30:
			for( l=0; l<2; l++ )
			{
				d = 0.0;
				for( i=0; i<3; i++ )
				{
					b[i][l] = r[i][0]*a[0][l] + r[i][1]*a[1][l] + r[i][2]*a[2][l];
					d += b[i][l]*b[i][l];
				}
				if( d > 0.0 ) d = 1.0 / sqrt( d );
				for( i=0; i<3; i++ )
				{
					b[i][l] *= d;
				}
			}

			//
			// Construct matrix b.
			// Note that although this appears similar to the construction of a,
			// it is not the same - different indices etc are used!
			//
			d = b[0][0]*b[0][1] + b[1][0]*b[1][1] + b[2][0]*b[2][1];
			p = 0.0;

			for( i=0; i<3; i++ )
			{
				b[i][1] -= d*b[i][0];
				p += b[i][1]*b[i][1];
			}
			if( p <= tol )
			{
				p = 1.0;
				for( i=0; i<3; i++ )
				{
					if( p < fabs( b[i][0] ) ) continue; // was "cycle" in fortran, but "continue" does same in c.
					p = fabs( b[i][0] );
					j = i;
				}
				k = ip1201[j];
				l = ip1201[j+1];
				p = sqrt( b[k][0]*b[k][0] + b[l][0]*b[l][0] );
				if( p <= tol ) goto label_40; // JUMP POINT HERE => build translation vector and get rms.
				b[j][1] = 0.0;
				b[k][1] = -b[l][0]/p;
				b[l][1] =  b[k][0]/p;
			}
			else
			{
				p = 1.0 / sqrt(p);
				for( i=0; i<3; i++ )
				{
					b[i][1] *= p;
				}
			}

			b[0][2] = b[1][0]*b[2][1] - b[1][1]*b[2][0];
			b[1][2] = b[2][0]*b[0][1] - b[2][1]*b[0][0];
			b[2][2] = b[0][0]*b[1][1] - b[0][1]*b[1][0];

			//
			// Construct rotation matrix M
			//
			for( i=0; i<3; i++ )
			{
				for( j=0; j<3; j++ )
				{
					M[i][j] = b[i][0]*a[j][0] + b[i][1]*a[j][1] + b[i][2]*a[j][2];
				}
			}

			//
			// Construct translation vector T
			//
			label_40:
			for( i=0; i<3; i++ ) T[i] = ((yc[i] - M[i][0]*xc[0]) - M[i][1]*xc[1] ) - M[i][2]*xc[2];

			//
			// Get RMSD, having cheked the eigen info. I'm not convinced by this.
			//
			label_50:
			for( i=0; i<3; i++ )
			{
				if( e[i] < 0.0 ) e[i] = 0.0;
				e[i] = sqrt( e[i] );
			}

			if( e[1] <= (e[0]*1.0e-5) ) return -1; // FIX THIS! LITERAL TOLERANCE VALUE

			d = e[2];
			if( sigma < 0.0 )
			{
				d = -d;
				if( (e[1]-e[2]) < (e[1]*1.0e-5) ) return -1; // FIX THIS! LITERAL TOLERANCE VALUE
			}
			d = (d+e[1]) + e[0];

			*rmsd = (e0-d) -d;
			if( *rmsd < 0.0 ) *rmsd = 0.0;
			
			return 0;
		}

		
	public:
		
		//
		// Calculate optimal superposition of point sets.
		//
		template< typename T_external >
		static T_external Calculate( const std::vector<T_external> &target_xyz, const std::vector<T_external> &current_xyz, T_external M[3][3], T_external T[3] )
		{
			T_external rmsd;
			
			assert( target_xyz.size() == current_xyz.size() );
			assert( target_xyz.size()%3 == 0 );
			// Insert check for at least 3 points?

			size_t N = target_xyz.size()/3;
			
			u3b( N, (T_external *)NULL, &target_xyz[0], &current_xyz[0], 1, &rmsd, M, T ); // mode fixed at 1 (calculate transform)!
			return rmsd;
		}

		//
		// Apply superposition described by M, T. Safe to pass same source_xyz and dest_xyz, as temporary x,y,z used in transform
		//
		template< typename T_external >
		static void Apply( const std::vector<T_external> &source_xyz, std::vector<T_external> &dest_xyz, const T_external M[3][3], const T_external T[3] )
		{
			T_internal x, y, z;
			
			assert( source_xyz.size() == dest_xyz.size() );
			assert( source_xyz.size()%3 == 0 );
			
			size_t N = source_xyz.size()/3;

			for( size_t i=0; i<N; i++ )
			{
				x = source_xyz[i*3 +0];
				y = source_xyz[i*3 +1];
				z = source_xyz[i*3 +2];

				x -= T[0];
				y -= T[1];
				z -= T[2];

				dest_xyz[i*3 +0] = x*M[0][0] + y*M[1][0] + z*M[2][0];
				dest_xyz[i*3 +1] = x*M[0][1] + y*M[1][1] + z*M[2][1];
				dest_xyz[i*3 +2] = x*M[0][2] + y*M[1][2] + z*M[2][2];
			}
		}
		
		//
		// Utility RMSD routines
		//
		template< typename T_external >
		static T_external RMSD( const std::vector<T_external> &xyz1, const std::vector<T_external> &xyz2, std::vector<T_external> *per_site_dr2 = NULL )
		{
			T_internal dx, dy, dz, dr2, acc = 0.0;
			
			assert( xyz1.size() == xyz2.size() );
			assert( xyz1.size()%3 == 0 );
			
			size_t N = xyz1.size()/3;
			
			if( per_site_dr2 != NULL ) (*per_site_dr2).resize( N );

			for( size_t i=0; i<N; i++ )
			{
				dx = xyz1[i*3 +0] - xyz2[i*3 +0];
				dy = xyz1[i*3 +1] - xyz2[i*3 +1];
				dz = xyz1[i*3 +2] - xyz2[i*3 +2];

				dr2 = dx*dx + dy*dy + dz*dz;

				acc += dr2;
				if( per_site_dr2 != NULL ) (*per_site_dr2)[i] = dr2;
			}
			
			return sqrt(acc/N);
		}
};

}

#endif
