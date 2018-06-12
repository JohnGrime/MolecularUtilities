/*
	Author: John Grime, The University of Chicago.
*/

#if !defined( UTIL_GEOMETRY_DEFINED )
#define UTIL_GEOMETRY_DEFINED

#include <math.h>
#include <tuple>

namespace Util
{
namespace Geometry
{

//
// Simple 3 element vector.
//
template <typename T>
struct Vec3
{
	T x, y, z;

	//
	// Define some standard operators, so we can use Vec3 structures as std::map keys etc.
	// All comparisons performed in order: z, y, x (ie z assumed to be "most significant").
	//
	template<typename T2> bool operator == (const Vec3<T2> &rhs) const { return std::tie(z,y,x) == std::tie(rhs.z,rhs.y,rhs.x); }
	template<typename T2> bool operator != (const Vec3<T2> &rhs) const { return std::tie(z,y,x) != std::tie(rhs.z,rhs.y,rhs.x); }
	template<typename T2> bool operator <  (const Vec3<T2> &rhs) const { return std::tie(z,y,x)  < std::tie(rhs.z,rhs.y,rhs.x); }
	template<typename T2> bool operator >  (const Vec3<T2> &rhs) const { return std::tie(z,y,x)  > std::tie(rhs.z,rhs.y,rhs.x); }
	template<typename T2> bool operator <= (const Vec3<T2> &rhs) const { return std::tie(z,y,x) <= std::tie(rhs.z,rhs.y,rhs.x); }
	template<typename T2> bool operator >= (const Vec3<T2> &rhs) const { return std::tie(z,y,x) >= std::tie(rhs.z,rhs.y,rhs.x); }
	template<typename T2> struct Vec3<T> operator + (const Vec3<T2> &rhs) const { return Vec3<T> { (T)(x+rhs.x), (T)(y+rhs.y), (T)(z+rhs.z) }; }
	template<typename T2> struct Vec3<T> operator - (const Vec3<T2> &rhs) const { return Vec3<T> { (T)(x-rhs.x), (T)(y-rhs.y), (T)(z-rhs.z) }; }
	template<typename T2> struct Vec3<T> operator / (const T2 &rhs      ) const { return Vec3<T> { (T)(x/rhs),   (T)(y/rhs),   (T)(z/rhs)   }; }
	template<typename T2> struct Vec3<T> operator * (const T2 &rhs      ) const { return Vec3<T> { (T)(x*rhs),   (T)(y*rhs),   (T)(z*rhs)   }; }
};


//
// Some handy routines, e.g. Orthogonal() followed by Cross() gives
// orthonormal coordinate system if single vector provided.
//

template<typename T>
T Dot( const Vec3<T>& u, const Vec3<T>& v )
{
	return ( u.x*v.x + u.y*v.y + u.z*v.z );
}

template<typename T>
Vec3<T> Cross( const Vec3<T>& u, const Vec3<T>& v )
{
	return { (u.y*v.z)-(u.z*v.y), (u.z*v.x)-(u.x*v.z), (u.x*v.y)-(u.y*v.x) };
}

template<typename T>
T Abs( const Vec3<T>& v )
{
	return sqrt( Dot(v,v) );
}

template<typename T>
Vec3<T> Orthogonal( const Vec3<T>& v_ )
{
	constexpr Vec3<T> zero = {0,0,0}, x = {1,0,0}, y = {0,1,0};
	
	if( v_ == zero ) return v_;

	Vec3<T> v = v_ / Abs(v_);

	//
	// Cross of v and x is orthogonal to both, unless
	// x is parallel to v (indicated by result of cross
	// being zero). If so, use cross of v and y: v cannot
	// be parallel to both as x and y are orthogonal.
	//
	Vec3<T> u = Cross(v,x);
	return (Dot(u,u) < 1.0e-3) ? (Cross(v,y)) : (u);
}


//
// Methods assuming a flat vector of coords etc
//

//
// Pass point = nullptr for rotation about system origin, else rotation axis will pass through that point.
//
template <typename T>
void Rotate( int N, T *xyz, const T *point, const T *axis, T theta )
{
	double a, b, c, u, v, w, u2, v2, w2, K, L, ct, st;

	if( point != nullptr )
	{
		a = point[0];
		b = point[1];
		c = point[2];
	}
	else
	{
		a = 0.0;
		b = 0.0;
		c = 0.0;
	}
	
	ct = cos( theta );
	st = sin( theta );

	u = axis[0];
	v = axis[1];
	w = axis[2];
	u2 = u*u;
	v2 = v*v;
	w2 = w*w;
	K = u2 + v2 + w2;
	L = sqrt( K );
	
	for( int i=0; i<N; i++ )
	{
		double x = xyz[i*3 +0] - a;
		double y = xyz[i*3 +1] - b;
		double z = xyz[i*3 +2] - c;

		double rx = a*(v2+w2) + u*(-b*v -c*w + u*x + v*y + w*z) + ( (x-a)*(v2+w2) + u*(b*v + c*w - v*y - w*z) ) * ct + L*(b*w - c*v - w*y + v*z)*st;
		rx = rx / K;

		double ry = b*(u2+w2) + v*(-a*u -c*w + u*x + v*y + w*z) + ( (y-b)*(u2+w2) + v*(a*u + c*w - u*x - w*z) ) * ct + L*(-a*w + c*u + w*x - u*z)*st;
		ry = ry / K;

		double rz = c*(u2+v2) + w*(-a*u -b*v + u*x + v*y + w*z) + ( (z-c)*(u2+v2) + w*(a*u + b*v - u*x - v*y) ) * ct + L*(a*v - b*u - v*x + u*y)*st;
		rz = rz / K;

		xyz[i*3 +0] = rx + a;
		xyz[i*3 +1] = ry + b;
		xyz[i*3 +2] = rz + c;
	}
}

//
// Hack-y trick; we assume the template data type T1 has members "x", "y" and "z"!
//

template <typename T1, typename T2>
void Rotate( int N, T1 *xyz, const T2 *point, const T2 *axis, T2 theta )
{
	double a, b, c, u, v, w, u2, v2, w2, K, L, ct, st;

	if( point != nullptr )
	{
		a = point[0];
		b = point[1];
		c = point[2];
	}
	else
	{
		a = 0.0;
		b = 0.0;
		c = 0.0;
	}
	
	ct = cos( theta );
	st = sin( theta );

	u = axis[0];
	v = axis[1];
	w = axis[2];
	u2 = u*u;
	v2 = v*v;
	w2 = w*w;
	K = u2 + v2 + w2;
	L = sqrt( K );
	
	for( int i=0; i<N; i++ )
	{
		double x = xyz[i].x - a;
		double y = xyz[i].y - b;
		double z = xyz[i].z - c;

		double rx = a*(v2+w2) + u*(-b*v -c*w + u*x + v*y + w*z) + ( (x-a)*(v2+w2) + u*(b*v + c*w - v*y - w*z) ) * ct + L*(b*w - c*v - w*y + v*z)*st;
		rx = rx / K;

		double ry = b*(u2+w2) + v*(-a*u -c*w + u*x + v*y + w*z) + ( (y-b)*(u2+w2) + v*(a*u + c*w - u*x - w*z) ) * ct + L*(-a*w + c*u + w*x - u*z)*st;
		ry = ry / K;

		double rz = c*(u2+v2) + w*(-a*u -b*v + u*x + v*y + w*z) + ( (z-c)*(u2+v2) + w*(a*u + b*v - u*x - v*y) ) * ct + L*(a*v - b*u - v*x + u*y)*st;
		rz = rz / K;

		xyz[i].x = rx + a;
		xyz[i].y = ry + b;
		xyz[i].z = rz + c;
	}
}

//
// This algorithm is via:
// http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere
// Default output is points on the unit sphere (radius = 1).
//
template<typename T>
int GetEquispacedSpherePoints( int target_N, std::vector<T>& xyz, double radius = 1.0 )
{
	static const double dl = M_PI * ( 3.0 - sqrt(5.0) );

	xyz.clear();

	double dz = 2.0 / target_N;
	double l = 0.0;
	double z = 1.0 - dz/2;

	for( int i=0; i<target_N; i++ )
	{
		double r = sqrt( 1.0 - z*z );
		double x = cos(l)*r;
		double y = sin(l)*r;

		xyz.push_back( x*radius );
		xyz.push_back( y*radius );
		xyz.push_back( z*radius );

		z -= dz;
		l += dl;
	}

	return (int)xyz.size()/3;
}

template<typename T>
int GetEquispacedSpherePoints( int target_N, std::vector< Vec3<T> >& xyz, double radius = 1.0 )
{
	static const double dl = M_PI * ( 3.0 - sqrt(5.0) );

	xyz.clear();

	double dz = 2.0 / target_N;
	double l = 0.0;
	double z = 1.0 - dz/2;

	for( int i=0; i<target_N; i++ )
	{
		double r = sqrt( 1.0 - z*z );
		double x = cos(l)*r;
		double y = sin(l)*r;

		xyz.push_back( {x*radius, y*radius, z*radius} );

		z -= dz;
		l += dl;
	}

	return (int)xyz.size();
}

}
}

#endif
