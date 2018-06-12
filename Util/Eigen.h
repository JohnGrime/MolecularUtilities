/*
	Author: John Grime, The University of Chicago.
*/

#if !defined( UTIL_EIGEN_DEFINED )
#define UTIL_EIGEN_DEFINED

//
// Based on Connelly Barnes' public domain conversion of the
// public domain JAMA library code (originally in Java).
//

#include <math.h>
#include <vector>

namespace Util
{

template <int dim>
class EigenUtil
{
	protected:
		
		static double hypot2( double x, double y )
		{
			return sqrt( x*x + y*y );
		}

		//
		// Symmetric Householder reduction to tridiagonal form.
		//
		static void tred2( double V[dim][dim], double d[dim], double e[dim] )
		{
			//
			//  This is derived from the Algol procedures tred2 by
			//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
			//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
			//  Fortran subroutine in EISPACK.
			//

			for( int j = 0; j < dim; j++ )
			{
				d[j] = V[dim-1][j];
			}

			//
			// Householder reduction to tridiagonal form.
			//

			for( int i = dim-1; i > 0; i-- )
			{
				//
				// Scale to avoid under/overflow.
				//

				double scale = 0.0;
				double h = 0.0;
				
				for( int k = 0; k < i; k++ )
				{
					scale = scale + fabs( d[k] );
				}
			
				if( scale == 0.0 )
				{
					e[i] = d[i-1];
					for( int j = 0; j < i; j++ )
					{
						d[j] = V[i-1][j];
						V[i][j] = 0.0;
						V[j][i] = 0.0;
					}
				}
				else
				{
					//
					// Generate Householder vector.
					//

					for( int k = 0; k < i; k++ )
					{
						d[k] /= scale;
						h += d[k] * d[k];
					}
					
					double f = d[i-1];
					double g = sqrt( h );
					
					if( f > 0 )
					{
						g = -g;
					}
					e[i] = scale * g;
					h = h - f * g;
					d[i-1] = f - g;
					for( int j = 0; j < i; j++ )
					{
						e[j] = 0.0;
					}

					//
					// Apply similarity transformation to remaining columns.
					//

					for( int j = 0; j < i; j++ )
					{
						f = d[j];
						V[j][i] = f;
						g = e[j] + V[j][j] * f;
						for( int k = j+1; k <= i-1; k++ )
						{
							g += V[k][j] * d[k];
							e[k] += V[k][j] * f;
						}
						e[j] = g;
					}
					f = 0.0;
					for( int j = 0; j < i; j++ )
					{
						e[j] /= h;
						f += e[j] * d[j];
					}
					double hh = f / (h + h);
					for( int j = 0; j < i; j++ )
					{
						e[j] -= hh * d[j];
					}
					for( int j = 0; j < i; j++ )
					{
						f = d[j];
						g = e[j];
						for( int k = j; k <= i-1; k++ )
						{
							V[k][j] -= (f * e[k] + g * d[k]);
						}
						d[j] = V[i-1][j];
						V[i][j] = 0.0;
					}
				}
				d[i] = h;
			}

			//
			// Accumulate transformations.
			//

			for( int i = 0; i < dim-1; i++ )
			{
				V[dim-1][i] = V[i][i];
				V[i][i] = 1.0;
				
				double h = d[i+1];
				
				if( h != 0.0 )
				{
					for( int k = 0; k <= i; k++ )
					{
						d[k] = V[k][i+1] / h;
					}
					for( int j = 0; j <= i; j++ )
					{
						double g = 0.0;
						for( int k = 0; k <= i; k++ )
						{
							g += V[k][i+1] * V[k][j];
						}
						for( int k = 0; k <= i; k++ )
						{
							V[k][j] -= g * d[k];
						}
					}
				}
				for( int k = 0; k <= i; k++ )
				{
					V[k][i+1] = 0.0;
				}
			}
			for( int j = 0; j < dim; j++ )
			{
				d[j] = V[dim-1][j];
				V[dim-1][j] = 0.0;
			}
			V[dim-1][dim-1] = 1.0;
			e[0] = 0.0;
		} 
		
		//
		// Symmetric tridiagonal QL algorithm.
		//
		static void tql2( double V[dim][dim], double d[dim], double e[dim] )
		{
			//
			//  This is derived from the Algol procedures tql2, by
			//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
			//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
			//  Fortran subroutine in EISPACK.
			//

			for( int i = 1; i < dim; i++ )
			{
				e[i-1] = e[i];
			}
			e[dim-1] = 0.0;

			double f = 0.0;
			double tst1 = 0.0;
			double eps = pow( 2.0, -52.0 );
			for( int l = 0; l < dim; l++ )
			{
				//
				// Find small subdiagonal element
				//

				tst1 = ( tst1 > (fabs(d[l])+fabs(e[l])) ) ? (tst1) : ( fabs(d[l])+fabs(e[l]) );
				int m = l;
				while( m < dim )
				{
					if( fabs(e[m]) <= eps*tst1 )
					{
						break;
					}
					m++;
				}

				//
				// If m == l, d[l] is an eigenvalue,
				// otherwise, iterate.
				//

				if( m > l )
				{
					int iter = 0;
					do
					{
						iter = iter + 1;  // (Could check iteration count here.)

						//
						// Compute implicit shift
						//

						double g = d[l];
						double p = (d[l+1] - g) / (2.0 * e[l]);
						double r = hypot2(p,1.0);
						if( p < 0 )
						{
							r = -r;
						}
						d[l] = e[l] / (p + r);
						d[l+1] = e[l] * (p + r);
						double dl1 = d[l+1];
						double h = g - d[l];
						for( int i = l+2; i < dim; i++ )
						{
							d[i] -= h;
						}
						f = f + h;

						//
						// Implicit QL transformation.
						//

						p = d[m];
						double c = 1.0;
						double c2 = c;
						double c3 = c;
						double el1 = e[l+1];
						double s = 0.0;
						double s2 = 0.0;
						for( int i = m-1; i >= l; i-- )
						{
							c3 = c2;
							c2 = c;
							s2 = s;
							g = c * e[i];
							h = c * p;
							r = hypot2(p,e[i]);
							e[i+1] = s * r;
							s = e[i] / r;
							c = p / r;
							p = c * d[i] - s * g;
							d[i+1] = h + s * (c * g + s * d[i]);

							//
							// Accumulate transformation.
							//

							for( int k = 0; k < dim; k++ )
							{
								h = V[k][i+1];
								V[k][i+1] = s * V[k][i] + c * h;
								V[k][i] = c * V[k][i] - s * h;
							}
						}
						p = -s * s2 * c3 * el1 * e[l] / dl1;
						e[l] = s * p;
						d[l] = c * p;

						//
						// Check for convergence.
						//
					} while ( fabs(e[l]) > eps*tst1 );
				}
				d[l] = d[l] + f;
				e[l] = 0.0;
			}

			//
			// Sort eigenvalues and corresponding vectors.
			//

			for( int i = 0; i < dim-1; i++ )
			{
				int k = i;
				double p = d[i];
				for( int j = i+1; j < dim; j++ )
				{
					if( d[j] < p )
					{
						k = j;
						p = d[j];
					}
				}
				if( k != i )
				{
					d[k] = d[i];
					d[i] = p;
					for( int j = 0; j < dim; j++ )
					{
						p = V[j][i];
						V[j][i] = V[j][k];
						V[j][k] = p;
					}
				}
			}
		}
		
	public:

		//
		// Generate centroid / covariance matrix from point set.
		// In:
		//    N   : number of points in point set xyz
		//    xyz : flat contiguous array of xyz coords for point set
		// Out:
		//    centroid : centroid of point set xyz
		//    covar    : covariance matrix of point set xyz
		//
		static int Covariance( int N, const double *coords, double centroid[dim], double covar[dim][dim] )
		{
			//
			// Sanity checks.
			//
			if( N < 2 || coords == nullptr )
			{
				return -1;
			}
	
			//
			// Calculate centroid
			//
			for( int i=0; i<dim; i++ ) { centroid[i] = 0.0; }
	
			for( int i=0; i<N; i++ )
			{
				for( int j=0; j<dim; j++ )
				{
					centroid[j] += coords[i*dim +j];
				}
			}
	
			for( int i=0; i<dim; i++ ) { centroid[i] /= N; }

			//
			// Covariance matrix
			//
			for( int i=0; i<dim; i++ )
			{
				for( int j=0; j<dim; j++ )
				{
					covar[i][j] = 0.0;
				}
			}

			for( int i=0; i<N; i++ )
			{
				for( int j=0; j<dim; j++ )
				{
					for( int k=0; k<dim; k++ )
					{
						covar[j][k] += (centroid[j]-coords[i*dim +j]) * (centroid[k]-coords[i*dim +k]);
					}
				}
			}
	
			for( int i=0; i<dim; i++ )
			{
				for( int j=0; j<dim; j++ )
				{
					covar[i][j] /= N-1; // should avoid DBZ, as checked if N < 2 earlier.
				}
			}
	
			return 1;
		}
		static int Covariance( const std::vector<double> &coords, double centroid[dim], double covar[dim][dim] )
		{
			return Covariance( (int)coords.size()/dim, &coords[0], centroid, covar );
		}
		

		//
		// Assumes A[][] is symmetrical, as are guaranteed in e.g. a covariance matrix.
		// Returns paired eigenvalues and eigenvectors.
		//
		// IMPORTANT: by default, this routine returns the eigenvals in ASCENDING order
		// and eigenvecs are COLUMN MAJOR. That is, the least significant eigenvalue is the
		// first entry in eigenvals with corresponding eigenvector in first column of eigenvecs.
		//
		// More useful might be the eigenvals in DESCENDING order, with corresponding vectors as
		// ROWS in eigenvals. This is enabled by the "swap_orders" parameter, so we have the most
		// significant eigenvalue as the first entry, with corresponding vector in the first row.
		//
		static void Eigen( double A[dim][dim], double eigenvals[dim], double eigenvecs[dim][dim], bool swap_orders = false )
		{
			double e[dim];
			
			for( int i = 0; i < dim; i++)
			{
				for( int j = 0; j < dim; j++)
				{
					eigenvecs[i][j] = A[i][j];
				}
			}
			tred2( eigenvecs, eigenvals, e );
			tql2(  eigenvecs, eigenvals, e );
			
			//
			// At this point, eigenvals are ascending order and eigenvecs are column-major.
			// Swap for more convenient descending order, and row major.
			//
			if( swap_orders == true )
			{
				for( int i=0; i<dim; i++ )
				{
					int X = (dim-1)-i;
					e[i] = eigenvals[X]; // <- e used as temporary storage
					for( int j=0; j<dim; j++ )
					{
						A[i][j] = eigenvecs[j][X]; // <- A used as temporary storage
					}
				}

				// copy back.
				for( int i=0; i<dim; i++ )
				{
					eigenvals[i] = e[i];
					for( int j=0; j<dim; j++ )
					{
						eigenvecs[i][j] = A[i][j];
					}
				}
			}
		}
		//
		// From a flat vector of coords
		//
		static void Eigen( const std::vector<double> &coords, double centroid[dim], double eigenvals[dim], double eigenvecs[dim][dim], bool swap_orders = false )
		{
			double covar[dim][dim];
			if( Covariance( coords, centroid, covar ) == -1 )
			{
				return;
			}
			Eigen( covar, eigenvals, eigenvecs, swap_orders );
		}
		
};

}

#endif
