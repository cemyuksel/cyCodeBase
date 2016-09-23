// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
/// \file		cyPoint.h 
/// \author		Cem Yuksel
/// \brief		2D, 3D and 4D point classes.
//-------------------------------------------------------------------------------
//
// Copyright (c) 2016, Cem Yuksel <cem@cemyuksel.com>
// All rights reserved.
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy 
// of this software and associated documentation files (the "Software"), to deal 
// in the Software without restriction, including without limitation the rights 
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
// copies of the Software, and to permit persons to whom the Software is 
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
// SOFTWARE.
// 
//-------------------------------------------------------------------------------

#ifndef _CY_POINT_H_INCLUDED_
#define _CY_POINT_H_INCLUDED_

//-------------------------------------------------------------------------------

#include <math.h>

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

/// 2D point class
template <typename TYPE>
class Point2
{
	friend Point2 operator + ( const TYPE v, const Point2 &pt ) { return pt+v; }	///< Addition with a constant
	friend Point2 operator - ( const TYPE v, const Point2 &pt ) { return -(pt-v); }	///< Subtraction from a constant
	friend Point2 operator * ( const TYPE v, const Point2 &pt ) { return pt*v; }	///< Multiplication with a constant

public:

	TYPE x, y;

	///@name Constructors
	Point2() {}
	Point2( const Point2<TYPE> &pt ) { x=pt.x; y=pt.y; }
	explicit Point2( TYPE _x, TYPE _y ) { x=_x; y=_y; }
	explicit Point2( const TYPE *pt ) { x=pt[0]; y=pt[1]; }
	explicit Point2( TYPE _xy ) { x=_xy; y=_xy; }
	template <typename T> Point2( const Point2<T> &pt ) { x=(TYPE)pt.x; y=(TYPE)pt.y; }

	///@name Set & Get value functions
	Point2& Zero() { x=0; y=0; return *this; }						///< Sets x and y coordinates as zero
	Point2& Set( TYPE _x, TYPE _y ) { x=_x; y=_y; return *this; }	///< Sets x and y coordinates as given
	Point2& Set( const TYPE *v ) { x=v[0]; y=v[1]; return *this; }	///< Sets x and y coordinates using the values in the given array
	void    GetValue( TYPE *v ) const { v[0]=x; v[1]=y; }			///< Puts x and y values into the array

	///@name Length and Normalize functions
	TYPE   LengthSquared () const { return x*x + y*y; }
	TYPE   Length        () const { return (TYPE) sqrt(LengthSquared()); }
	void   Normalize     ()       { TYPE s = Length(); *this /= s; }
	Point2 GetNormalized () const { TYPE s = Length(); return *this / s; }
	TYPE   MaxComponent  () const { return x>y ? x : y; }
	int    MaxComponentID() const { return x>y ? 0 : 1; }

	///@name Limit functions
	void ClampMinMax( TYPE minValue, TYPE maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( TYPE n ) { if(x<n)x=n; if(y<n)y=n; }
	void ClampMax( TYPE n ) { if(x>n)x=n; if(y>n)y=n; }

	///@name Unary operators
	Point2 operator - () const { return Point2(-x,-y); } 
	Point2 operator + () const { return *this; }

	///@name Binary operators
	Point2 operator + ( const Point2 &pt ) const { return Point2(x+pt.x, y+pt.y); }
	Point2 operator - ( const Point2 &pt ) const { return Point2(x-pt.x, y-pt.y); }
	Point2 operator * ( const Point2 &pt ) const { return Point2(x*pt.x, y*pt.y); }
	Point2 operator / ( const Point2 &pt ) const { return Point2(x/pt.x, y/pt.y); }
	Point2 operator + ( TYPE n ) const { return Point2(x+n, y+n); }
	Point2 operator - ( TYPE n ) const { return Point2(x-n, y-n); }
	Point2 operator * ( TYPE n ) const { return Point2(x*n, y*n); }
	Point2 operator / ( TYPE n ) const { return Point2(x/n, y/n); }

	///@name Assignment operators
	Point2& operator += ( const Point2 &pt ) { x+=pt.x; y+=pt.y; return *this; }
	Point2& operator -= ( const Point2 &pt ) { x-=pt.x; y-=pt.y; return *this; }
	Point2& operator *= ( const Point2 &pt ) { x*=pt.x; y*=pt.y; return *this; }
	Point2& operator /= ( const Point2 &pt ) { x/=pt.x; y/=pt.y; return *this; }
	Point2& operator += ( TYPE n ) { x+=n; y+=n; return *this; }
	Point2& operator -= ( TYPE n ) { x-=n; y-=n; return *this; }
	Point2& operator *= ( TYPE n ) { x*=n; y*=n; return *this; }
	Point2& operator /= ( TYPE n ) { x/=n; y/=n; return *this; }

	///@name Test operators
	int operator == ( const Point2& pt ) const { return ( (pt.x==x) && (pt.y==y) ); }
	int operator != ( const Point2& pt ) const { return ( (pt.x!=x) || (pt.y!=y) ); }

	///@name Access operators
	TYPE& operator [] ( int i )       { return (&x)[i]; }
	TYPE  operator [] ( int i ) const { return (&x)[i]; }

	///@name Cross product and dot product
	TYPE Cross      ( const Point2 &pt ) const { return x*pt.y-y*pt.x; }	///< Cross product
	TYPE operator ^ ( const Point2 &pt ) const { return Cross(pt); }		///< Cross product operator
	TYPE Dot        ( const Point2 &pt ) const { return x*pt.x + y*pt.y; }	///< Dot product
	TYPE operator % ( const Point2 &pt ) const { return Dot(pt); }			///< Dot product operator

};

float Point2<float>::Length() const { return sqrtf(LengthSquared()); }

//-------------------------------------------------------------------------------

/// 3D point class
template <typename TYPE>
class Point3
{
	friend Point3 operator+( const TYPE v, const Point3 &pt ) { return pt+v; }		///< Addition with a constant
	friend Point3 operator-( const TYPE v, const Point3 &pt ) { return -(pt-v); }	///< Subtraction from a constant
	friend Point3 operator*( const TYPE v, const Point3 &pt ) { return pt*v; }		///< Multiplication with a constant

public:

	TYPE x, y, z;

	///@name Constructors
	Point3() { }
	Point3( const Point3<TYPE> &pt ) { x=pt.x; y=pt.y; z=pt.z; }
	Point3( TYPE _x, TYPE _y, TYPE _z ) { x=_x; y=_y; z=_z; }
	explicit Point3( const TYPE *pt ) { x=pt[0]; y=pt[1]; z=pt[2]; }
	explicit Point3( TYPE _xyz ) { x=_xyz; y=_xyz; z=_xyz; }
	explicit Point3( const Point2<TYPE> &pt ) { x=pt.x; y=pt.y; z=0; }
	explicit Point3( const Point2<TYPE> &pt, TYPE _z ) { x=pt.x; y=pt.y; z=_z; }
	template <typename T> Point3( const Point3<T> &pt ) { x=(TYPE)pt.x; y=(TYPE)pt.y; z=(TYPE)pt.z; }

	///@name Set & Get value functions
	Point3& Zero() { x=0; y=0; z=0; return *this; }									///< Sets x, y and z coordinates as zero
	Point3& Set( TYPE _x, TYPE _y, TYPE _z ) { x=_x; y=_y; z=_z; return *this; }	///< Sets x, y and z coordinates as given
	Point3& Set( const TYPE *v ) { x=v[0]; y=v[1]; z=v[2]; return *this; }			///< Sets x, y and z coordinates using the values in the given array
	void    GetValue( TYPE *v ) const { v[0]=x; v[1]=y; v[2]=z; }							///< Puts x, y and z values into the array

	///@name Length and Normalize functions
	TYPE   LengthSquared () const { return x*x + y*y + z*z; }
	TYPE   Length        () const { return (TYPE) sqrt(LengthSquared()); }
	void   Normalize     ()       { TYPE s = Length(); *this /= s; }
	Point3 GetNormalized () const { TYPE s = Length(); return *this / s; }
	TYPE   MaxComponent  () const { return x>y ? (x>z ? x : z) : (y>z ? y : z); }
	int    MaxComponentID() const { return x>y ? (x>z ? 0 : 2) : (y>z ? 1 : 2); }

	///@name Limit functions
	void ClampMinMax( TYPE min, TYPE max ) { ClampMin(min); ClampMax(max); }
	void ClampMin( TYPE n ) { if(x<n)x=n; if(y<n)y=n; if(z<n)z=n; }
	void ClampMax( TYPE n ) { if(x>n)x=n; if(y>n)y=n; if(z>n)z=n; }

	///@name Unary operators
	Point3 operator - () const { return Point3(-x,-y,-z); } 
	Point3 operator + () const { return *this; }

	///@name Binary operators
	Point3 operator + ( const Point3 &pt ) const { return Point3(x+pt.x, y+pt.y, z+pt.z); }
	Point3 operator - ( const Point3 &pt ) const { return Point3(x-pt.x, y-pt.y, z-pt.z); }
	Point3 operator * ( const Point3 &pt ) const { return Point3(x*pt.x, y*pt.y, z*pt.z); }
	Point3 operator / ( const Point3 &pt ) const { return Point3(x/pt.x, y/pt.y, z/pt.z); }
	Point3 operator + ( TYPE n ) const { return Point3(x+n, y+n, z+n); }
	Point3 operator - ( TYPE n ) const { return Point3(x-n, y-n, z-n); }
	Point3 operator * ( TYPE n ) const { return Point3(x*n, y*n, z*n); }
	Point3 operator / ( TYPE n ) const { return Point3(x/n, y/n, z/n); }

	///@name Assignment operators
	Point3& operator += ( const Point3 &pt ) { x+=pt.x; y+=pt.y; z+=pt.z; return *this; }
	Point3& operator -= ( const Point3 &pt ) { x-=pt.x; y-=pt.y; z-=pt.z; return *this; }
	Point3& operator *= ( const Point3 &pt ) { x*=pt.x; y*=pt.y; z*=pt.z; return *this; }
	Point3& operator /= ( const Point3 &pt ) { x/=pt.x; y/=pt.y; z/=pt.z; return *this; }
	Point3& operator += ( TYPE n ) { x+=n; y+=n; z+=n; return *this; }
	Point3& operator -= ( TYPE n ) { x-=n; y-=n; z-=n; return *this; }
	Point3& operator *= ( TYPE n ) { x*=n; y*=n; z*=n; return *this; }
	Point3& operator /= ( TYPE n ) { x/=n; y/=n; z/=n; return *this; }

	///@name Test operators
	int operator == ( const Point3& pt ) const { return ( (pt.x==x) && (pt.y==y) && (pt.z==z) ); }
	int operator != ( const Point3& pt ) const { return ( (pt.x!=x) || (pt.y!=y) || (pt.z!=z) ); }

	///@name Access operators
	TYPE& operator [] ( int i )       { return (&x)[i]; }
	TYPE  operator [] ( int i ) const { return (&x)[i]; }

	///@name Cross product and dot product
	Point3 Cross      ( const Point3 &pt ) const { return Point3(y*pt.z-z*pt.y, z*pt.x-x*pt.z, x*pt.y-y*pt.x); }	///< Cross product
	Point3 operator ^ ( const Point3 &pt ) const { return Cross(pt); }					///< Cross product
	TYPE   Dot        ( const Point3 &pt ) const { return x*pt.x + y*pt.y + z*pt.z; }	///< Dot product
	TYPE   operator % ( const Point3 &pt ) const { return Dot(pt); }					///< Dot product

	///@name Conversion Methods
	Point2<TYPE> XY() const { return Point2<TYPE>(x,y); }
};

float Point3<float>::Length() const { return sqrtf(LengthSquared()); }

//-------------------------------------------------------------------------------

/// 4D point class
template <typename TYPE>
class Point4
{
	friend Point4 operator + ( const TYPE v, const Point4 &pt ) { return pt+v; }		///< Addition with a constant
	friend Point4 operator - ( const TYPE v, const Point4 &pt ) { return -(pt-v); }	///< Subtraction from a constant
	friend Point4 operator * ( const TYPE v, const Point4 &pt ) { return pt*v; }		///< Multiplication with a constant

public:

	TYPE x, y, z, w;

	///@name Constructors
	Point4() { }
	Point4( const Point4<TYPE> &pt ) { x=pt.x; y=pt.y; z=pt.z; w=pt.w; }
	Point4( TYPE _x, TYPE _y, TYPE _z, TYPE _w ) { x=_x; y=_y; z=_z; w=_w; }
	explicit Point4( const TYPE *pt ) { x=pt[0]; y=pt[1]; z=pt[2]; w=pt[3]; }
	explicit Point4( TYPE _xyzw ) { x=_xyzw; y=_xyzw; z=_xyzw; w=_xyzw; }
	explicit Point4( const Point3<TYPE> &pt ) { x=pt.x; y=pt.y; z=pt.z; w=1; }
	explicit Point4( const Point3<TYPE> &pt, TYPE _w ) { x=pt.x; y=pt.y; z=pt.z; w=_w; }
	template <typename T> Point4( const Point4<T> &pt ) { x=(TYPE)pt.x; y=(TYPE)pt.y; z=(TYPE)pt.z; w=(TYPE)pt.w; }

	///@name Set & Get value functions
	Point4& Zero() { x=0; y=0; z=0; w=0; return *this; }												///< Sets x, y, z and w coordinates as zero
	Point4& Set( TYPE _x, TYPE _y, TYPE _z, TYPE _w ) { x=_x; y=_y; z=_z; w=_w; return *this; }	///< Sets x, y, z and w coordinates as given
	Point4& Set( const TYPE *v ) { x=v[0]; y=v[1]; z=v[2]; w=v[3]; return *this; }					///< Sets x, y, z and w coordinates using the values in the given array
	void    GetValue( TYPE *v ) const { v[0]=x; v[1]=y; v[2]=z; v[3]=w; }									///< Puts x, y, z and w values into the array

	///@name Length and Normalize functions
	TYPE   LengthSquared () const { return x*x + y*y + z*z + w*w; }
	TYPE   Length        () const { return (TYPE) sqrt(LengthSquared()); }
	void   Normalize     ()       { TYPE s = Length(); *this /= s; }
	Point4 GetNormalized () const { TYPE s = Length(); return *this / s; }
	TYPE   MaxComponent  () const { TYPE mxy = x>y ? x : y; TYPE mzw = z>w ? z : w; return mxy>mzw ? mxy : mzw; }
	int    MaxComponentID() const { int  ixy = x>y ? 0 : 1; int  izw = z>w ? 2 : 3; return (&x)[ixy]>(&x)[izw] ? ixy : izw; }

	///@name Limit functions
	void ClampMinMax( TYPE min, TYPE max ) { ClampMin(min); ClampMax(max); }
	void ClampMin( TYPE n ) { if(x<n)x=n; if(y<n)y=n; if(z<n)z=n; if(w<n)w=n; }
	void ClampMax( TYPE n ) { if(x>n)x=n; if(y>n)y=n; if(z>n)z=n; if(w>n)w=n; }

	///@name Unary operators
	Point4 operator - () const { return Point4(-x, -y, -z, -w); }
	Point4 operator + () const { return *this; }

	///@name Binary operators
	Point4 operator + ( const Point4 &pt ) const { return Point4(x+pt.x, y+pt.y, z+pt.z, w+pt.w); }
	Point4 operator - ( const Point4 &pt ) const { return Point4(x-pt.x, y-pt.y, z-pt.z, w-pt.w); }
	Point4 operator * ( const Point4 &pt ) const { return Point4(x*pt.x, y*pt.y, z*pt.z, w*pt.w); }
	Point4 operator / ( const Point4 &pt ) const { return Point4(x/pt.x, y/pt.y, z/pt.z, w/pt.w); }
	Point4 operator + ( TYPE n ) const { return Point4(x+n, y+n, z+n, w+n); }
	Point4 operator - ( TYPE n ) const { return Point4(x-n, y-n, z-n, w-n); }
	Point4 operator * ( TYPE n ) const { return Point4(x*n, y*n, z*n, w*n); }
	Point4 operator / ( TYPE n ) const { return Point4(x/n, y/n, z/n, w/n); }

	///@name Assignment operators
	Point4& operator += ( const Point4 &pt ) { x+=pt.x; y+=pt.y; z+=pt.z; w+=pt.w; return *this; }
	Point4& operator -= ( const Point4 &pt ) { x-=pt.x; y-=pt.y; z-=pt.z; w-=pt.w; return *this; }
	Point4& operator *= ( const Point4 &pt ) { x*=pt.x; y*=pt.y; z*=pt.z; w*=pt.w; return *this; }
	Point4& operator /= ( const Point4 &pt ) { x/=pt.x; y/=pt.y; z/=pt.z; w/=pt.w; return *this; }
	Point4& operator += ( TYPE n ) { x+=n; y+=n; z+=n; w+=n; return *this; }
	Point4& operator -= ( TYPE n ) { x-=n; y-=n; z-=n; w-=n; return *this; }
	Point4& operator *= ( TYPE n ) { x*=n; y*=n; z*=n; w*=n; return *this; }
	Point4& operator /= ( TYPE n ) { x/=n; y/=n; z/=n; w/=n; return *this; }

	///@name Test operators
	int operator == ( const Point4& pt ) const { return ( (pt.x==x) && (pt.y==y) && (pt.z==z) && (pt.w==w) ); }
	int operator != ( const Point4& pt ) const { return ( (pt.x!=x) || (pt.y!=y) || (pt.z!=z) || (pt.w!=w) ); }

	///@name Access operators
	TYPE& operator [] (int i)       { return ( &x )[i]; }
	TYPE  operator [] (int i) const { return ( &x )[i]; }

	///@ Dot product
	TYPE Dot		( const Point4 &pt ) const { return x*pt.x + y*pt.y + z*pt.z + w*pt.w; }	///< Dot product
	TYPE operator % ( const Point4 &pt ) const { return Dot(pt); }								///< Dot product

	///@name Conversion Methods
	Point2<TYPE> XY () const { return Point2<TYPE>(x,y); }
	Point3<TYPE> XYZ() const { return Point3<TYPE>(x,y,z); }
	Point3<TYPE> GetNonHomogeneous() const { return Point3<TYPE>(x,y,z)/w; }
};

float Point4<float>::Length() const { return sqrtf(LengthSquared()); }

//-------------------------------------------------------------------------------

typedef Point2<float>  Point2f;
typedef Point3<float>  Point3f;
typedef Point4<float>  Point4f;

typedef Point2<double> Point2d;
typedef Point3<double> Point3d;
typedef Point4<double> Point4d;

typedef Point2<int>    Point2i;
typedef Point3<int>    Point3i;
typedef Point4<int>    Point4i;

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Point2f cyPoint2f;
typedef cy::Point3f cyPoint3f;
typedef cy::Point4f cyPoint4f;

typedef cy::Point2d cyPoint2d;
typedef cy::Point3d cyPoint3d;
typedef cy::Point4d cyPoint4d;

typedef cy::Point2i cyPoint2i;
typedef cy::Point3i cyPoint3i;
typedef cy::Point4i cyPoint4i;

//-------------------------------------------------------------------------------

#endif

