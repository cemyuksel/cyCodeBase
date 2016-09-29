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

#include "cyCore.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

// Forward declarations
template <typename TYPE> class Point2;
template <typename TYPE> class Point3;
template <typename TYPE> class Point4;

//-------------------------------------------------------------------------------

/// A general class for N-dimensional points (vectors).
template <typename TYPE, int N>
class Point
{
	friend Point operator + ( const TYPE v, const Point &pt ) { return pt+v; }	///< Addition with a constant
	friend Point operator - ( const TYPE v, const Point &pt ) { return -(pt-v); }	///< Subtraction from a constant
	friend Point operator * ( const TYPE v, const Point &pt ) { return pt*v; }	///< Multiplication with a constant

public:

	TYPE data[N];

	///@name Constructors
	Point() {}
	Point( const Point &pt )        { CY_MEMCOPY(TYPE,data,pt.data,N); }
	explicit Point( const TYPE *p ) { CY_MEMCOPY(TYPE,data,pt,N); }
	explicit Point( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i]=v; }
	template <typename T> explicit Point( const Point<T,N> &pt ) { CY_MEMCONVERT(TYPE,data,pt.data,N); }
	template <int M> explicit Point( const Point<TYPE,M> &pt )
	{
		if ( N <= M ) CY_MEMCOPY(TYPE,data,pt.data,N);
		else {        CY_MEMCOPY(TYPE,data,pt.data,M); CY_MEMCLEAR(TYPE,data,N-M); }
	}
	template <typename T, int M> explicit Point( const Point<T,M> &pt )
	{
		if ( N <= M ) CY_MEMCONVERT(TYPE,data,pt.data,N);
		else {        CY_MEMCONVERT(TYPE,data,pt.data,M); CY_MEMCLEAR(TYPE,data,N-M); }
	}
	explicit Point( const Point2<TYPE> &pt );
	explicit Point( const Point3<TYPE> &pt );
	explicit Point( const Point4<TYPE> &pt );
	template <typename T> explicit Point( const Point2<T> &pt );
	template <typename T> explicit Point( const Point3<T> &pt );
	template <typename T> explicit Point( const Point4<T> &pt );

	///@name Set & Get value functions
	void Zero()               { CY_MEMCLEAR(TYPE,data,N); }		///< Sets the coordinates as zero
	void Get( TYPE *p ) const { CY_MEMCOPY(TYPE,p,data,N); }	///< Puts the coordinate values into the array
	void Set( const TYPE *p ) { CY_MEMCOPY(TYPE,data,p,N); }	///< Sets the coordinates using the values in the given array

	///@name Length and Normalize functions
	TYPE  LengthSquared () const { Point p=operator*(*this); return p.Sum(); }	///< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	TYPE  Length        () const { return (TYPE) cySqrt(LengthSquared()); }		///< Returns the length of the vector.
	void  Normalize     ()       { *this /= Length(); }							///< Normalizes the vector, such that its length becomes 1.
	Point GetNormalized () const { return *this / Length(); }					///< Returns a normalized copy of the vector.
	TYPE  Sum           () const { TYPE v=data[0]; for ( int i=1; i<N; ++i ) v+=data[i]; return v; }		///< Returns the sum of its components
	bool  IsZero        () const { Point zero(TYPE(0)); return memcmp(data,zero.data,sizeof(data))==0; }	///< Returns true if all components are exactly zero

	///@name Limit functions
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i] = (data[i]<v) ? v : data[i]; }
	void ClampMax( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i] = (data[i]>v) ? v : data[i]; }
	void Abs() { for ( int i=0; i<N; i++ ) data[i] = cyAbs(data[i]); }	///< Converts all negative components to positive values

	///@name Unary operators
	Point operator - () const { Point pt; for ( int i=0; i<N; ++i ) pt.data[i]=-data[i]; return pt; } 
	Point operator + () const { return *this; }

	///@name Binary operators
	Point operator + ( const Point &p ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] + p.data[i]; return r; }
	Point operator - ( const Point &p ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] - p.data[i]; return r; }
	Point operator * ( const Point &p ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] * p.data[i]; return r; }
	Point operator / ( const Point &p ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] / p.data[i]; return r; }
	Point operator + ( const TYPE  &v ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] + v; return r; }
	Point operator - ( const TYPE  &v ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] - v; return r; }
	Point operator * ( const TYPE  &v ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] * v; return r; }
	Point operator / ( const TYPE  &v ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] / v; return r; }

	///@name Assignment operators
	const Point& operator  = ( const Point &p ) { CY_MEMCOPY(TYPE,data,p.data,N); return *this; }	
	const Point& operator += ( const Point &p ) { for ( int i=0; i<N; ++i ) data[i] += p.data[i]; return *this; }
	const Point& operator -= ( const Point &p ) { for ( int i=0; i<N; ++i ) data[i] -= p.data[i]; return *this; }
	const Point& operator *= ( const Point &p ) { for ( int i=0; i<N; ++i ) data[i] *= p.data[i]; return *this; }
	const Point& operator /= ( const Point &p ) { for ( int i=0; i<N; ++i ) data[i] /= p.data[i]; return *this; }
	const Point& operator += ( const TYPE   v ) { for ( int i=0; i<N; ++i ) data[i] += v; return *this; }
	const Point& operator -= ( const TYPE   v ) { for ( int i=0; i<N; ++i ) data[i] -= v; return *this; }
	const Point& operator *= ( const TYPE   v ) { for ( int i=0; i<N; ++i ) data[i] *= v; return *this; }
	const Point& operator /= ( const TYPE   v ) { for ( int i=0; i<N; ++i ) data[i] /= v; return *this; }

	///@name Test operators
	bool operator == ( const Point& pt ) const { Point p=operator-(pt); return  p.IsZero(); }
	bool operator != ( const Point& pt ) const { Point p=operator-(pt); return !p.IsZero(); }

	///@name Access operators
	TYPE& operator [] ( int i )       { return data[i]; }
	TYPE  operator [] ( int i ) const { return data[i]; }

	///@name Dot product
	TYPE Dot        ( const Point &pt ) const { Point p=operator*(pt); return p.Sum(); }	///< Dot product
	TYPE operator % ( const Point &pt ) const { return Dot(pt); }							///< Dot product operator
};

//-------------------------------------------------------------------------------

/// 2D point (vector) class
template <typename TYPE>
class Point2
{
	friend Point2 operator + ( const TYPE v, const Point2 &pt ) { return   pt+v;  }	///< Addition with a constant
	friend Point2 operator - ( const TYPE v, const Point2 &pt ) { return -(pt-v); }	///< Subtraction from a constant
	friend Point2 operator * ( const TYPE v, const Point2 &pt ) { return   pt*v;  }	///< Multiplication with a constant

public:

	union {
		struct {
			TYPE x, y;	///< Components of the vector
		};
		TYPE data[2];	///< Components of the vector in array form
	};

	///@name Constructors
	Point2() {}
	Point2( const TYPE &_x, const TYPE &_y ) : x(_x), y(_y) {}
	Point2( const Point2 &pt )        { CY_MEMCOPY(TYPE,data,pt.data,2); }
	explicit Point2( const TYPE *pt ) { CY_MEMCOPY(TYPE,data,pt,2); }
	explicit Point2( const TYPE &v  ) { for ( int i=0; i<2; ++i ) data[i]=v; }
	explicit Point2( const Point3<TYPE> &pt );
	explicit Point2( const Point4<TYPE> &pt );
	template <typename T> explicit Point2( const Point2<T> &pt ) { CY_MEMCONVERT(TYPE,data,pt.data,2); }
	template <typename T> explicit Point2( const Point3<T> &pt );
	template <typename T> explicit Point2( const Point4<T> &pt );
	template <int M> explicit Point2( const Point<TYPE,M> &pt )
	{
		if ( 2 <= M ) CY_MEMCOPY(TYPE,data,pt.data,2);
		else {        CY_MEMCOPY(TYPE,data,pt.data,M); CY_MEMCLEAR(TYPE,data,2-M); }
	}
	template <typename T, int M> explicit Point2( const Point<T,M> &pt )
	{
		if ( 2 <= M ) CY_MEMCONVERT(TYPE,data,pt.data,2);
		else {        CY_MEMCONVERT(TYPE,data,pt.data,M); CY_MEMCLEAR(TYPE,data,2-M); }
	}

	///@name Set & Get value functions
	void Zero()               { CY_MEMCLEAR(TYPE,data,2); }		///< Sets the coordinates as zero.
	void Get( TYPE *p ) const { CY_MEMCOPY(TYPE,p,data,2); }	///< Puts the coordinate values into the array.
	void Set( const TYPE *p ) { CY_MEMCOPY(TYPE,data,p,2); }	///< Sets the coordinates using the values in the given array.
	void Set( const TYPE &_x, const TYPE &_y ) { x=_x; y=_y; }	///< Sets the coordinates using the given values

	///@name Length and Normalize functions
	TYPE   LengthSquared () const { Point2 p=operator*(*this); return p.Sum(); }	///< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	TYPE   Length        () const { return (TYPE) cySqrt(LengthSquared()); }			///< Returns the length of the vector.
	void   Normalize     ()       { *this /= Length(); }							///< Normalizes the vector, such that its length becomes 1.
	Point2 GetNormalized () const { return *this / Length(); }						///< Returns a normalized copy of the vector.
	TYPE   Sum           () const { return x+y; }									///< Returns the sum of its components
	bool   IsZero        () const { TYPE zero[2]={0,0}; return memcmp(data,zero,sizeof(data))==0; }	///< Returns true if all components are exactly zero
	TYPE   MaxComponent  () const { return x>y ? x : y; }
	int    MaxComponentID() const { return x>y ? 0 : 1; }

	///@name Limit functions
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { for ( int i=0; i<2; ++i ) data[i] = (data[i]<v) ? v : data[i]; }
	void ClampMax( const TYPE &v ) { for ( int i=0; i<2; ++i ) data[i] = (data[i]>v) ? v : data[i]; }
	void Abs() { for ( int i=0; i<2; i++ ) data[i] = cyAbs(data[i]); }	///< Converts all negative components to positive values

	///@name Unary operators
	Point2 operator - () const { Point2 pt; for ( int i=0; i<2; ++i ) pt.data[i]=-data[i]; return pt; } 
	Point2 operator + () const { return *this; }

	///@name Binary operators
	Point2 operator + ( const Point2 &p ) const { Point2 r; for ( int i=0; i<2; ++i ) r.data[i] = data[i] + p.data[i]; return r; }
	Point2 operator - ( const Point2 &p ) const { Point2 r; for ( int i=0; i<2; ++i ) r.data[i] = data[i] - p.data[i]; return r; }
	Point2 operator * ( const Point2 &p ) const { Point2 r; for ( int i=0; i<2; ++i ) r.data[i] = data[i] * p.data[i]; return r; }
	Point2 operator / ( const Point2 &p ) const { Point2 r; for ( int i=0; i<2; ++i ) r.data[i] = data[i] / p.data[i]; return r; }
	Point2 operator + ( const TYPE   &v ) const { Point2 r; for ( int i=0; i<2; ++i ) r.data[i] = data[i] + v;         return r; }
	Point2 operator - ( const TYPE   &v ) const { Point2 r; for ( int i=0; i<2; ++i ) r.data[i] = data[i] - v;         return r; }
	Point2 operator * ( const TYPE   &v ) const { Point2 r; for ( int i=0; i<2; ++i ) r.data[i] = data[i] * v;         return r; }
	Point2 operator / ( const TYPE   &v ) const { Point2 r; for ( int i=0; i<2; ++i ) r.data[i] = data[i] / v;         return r; }

	///@name Assignment operators
	const Point2& operator  = ( const Point2 &p ) { CY_MEMCOPY(TYPE,data,p.data,2); return *this; }	
	const Point2& operator += ( const Point2 &p ) { for ( int i=0; i<2; ++i ) data[i] += p.data[i]; return *this; }
	const Point2& operator -= ( const Point2 &p ) { for ( int i=0; i<2; ++i ) data[i] -= p.data[i]; return *this; }
	const Point2& operator *= ( const Point2 &p ) { for ( int i=0; i<2; ++i ) data[i] *= p.data[i]; return *this; }
	const Point2& operator /= ( const Point2 &p ) { for ( int i=0; i<2; ++i ) data[i] /= p.data[i]; return *this; }
	const Point2& operator += ( const TYPE   &v ) { for ( int i=0; i<2; ++i ) data[i] += v;         return *this; }
	const Point2& operator -= ( const TYPE   &v ) { for ( int i=0; i<2; ++i ) data[i] -= v;         return *this; }
	const Point2& operator *= ( const TYPE   &v ) { for ( int i=0; i<2; ++i ) data[i] *= v;         return *this; }
	const Point2& operator /= ( const TYPE   &v ) { for ( int i=0; i<2; ++i ) data[i] /= v;         return *this; }

	///@name Test operators
	bool operator == ( const Point2& pt ) const { Point2 p=operator-(pt); return  p.IsZero(); }
	bool operator != ( const Point2& pt ) const { Point2 p=operator-(pt); return !p.IsZero(); }

	///@name Access operators
	TYPE& operator [] ( int i )       { return data[i]; }
	TYPE  operator [] ( int i ) const { return data[i]; }

	///@name Cross product and dot product
	TYPE Cross      ( const Point2 &pt ) const { Point2 p(-y,x); return p.Dot(pt); }		///< Cross product
	TYPE operator ^ ( const Point2 &pt ) const { return Cross(pt); }						///< Cross product operator
	TYPE Dot        ( const Point2 &pt ) const { Point2 p=operator*(pt); return p.Sum(); }	///< Dot product
	TYPE operator % ( const Point2 &pt ) const { return Dot(pt); }							///< Dot product operator
};

//-------------------------------------------------------------------------------

/// 3D point class
template <typename TYPE>
class Point3
{
	friend Point3 operator + ( const TYPE v, const Point3 &pt ) { return   pt+v;  }	///< Addition with a constant
	friend Point3 operator - ( const TYPE v, const Point3 &pt ) { return -(pt-v); }	///< Subtraction from a constant
	friend Point3 operator * ( const TYPE v, const Point3 &pt ) { return   pt*v;  }	///< Multiplication with a constant

public:

	union {
		struct {
			TYPE x, y, z;	///< Components of the vector
		};
		TYPE data[3];		///< Components of the vector in array form
	};

	///@name Constructors
	Point3() { }
	Point3( const TYPE &_x, const TYPE &_y, const TYPE &_z ) : x(_x), y(_y), z(_z) {}
	Point3( const Point3 &pt )        { CY_MEMCOPY(TYPE,data,pt.data,3); }
	explicit Point3( const TYPE *pt ) { CY_MEMCOPY(TYPE,data,pt,3); }
	explicit Point3( const TYPE &v  ) { for ( int i=0; i<3; ++i ) data[i]=v; }
	explicit Point3( const Point2<TYPE> &pt, TYPE _z=0 ) { CY_MEMCOPY(TYPE,data,pt.data,2); z=_z; }
	explicit Point3( const Point4<TYPE> &pt );
	template <typename T> explicit Point3( const Point3<T> &pt )            { CY_MEMCONVERT(TYPE,data,pt.data,3); }
	template <typename T> explicit Point3( const Point2<T> &pt, TYPE _z=0 ) { CY_MEMCONVERT(TYPE,data,pt.data,2); z=_z; }
	template <typename T> explicit Point3( const Point4<T> &pt );
	template <int M> explicit Point3( const Point<TYPE,M> &pt )
	{
		if ( 3 <= M ) CY_MEMCOPY(TYPE,data,pt.data,3);
		else {        CY_MEMCOPY(TYPE,data,pt.data,M); CY_MEMCLEAR(TYPE,data,3-M); }
	}
	template <typename T, int M> explicit Point3( const Point<T,M> &pt )
	{
		if ( 3 <= M ) CY_MEMCONVERT(TYPE,data,pt.data,3);
		else {        CY_MEMCONVERT(TYPE,data,pt.data,M); CY_MEMCLEAR(TYPE,data,3-M); }
	}

	///@name Set & Get value functions
	void Zero()               { CY_MEMCLEAR(TYPE,data,3); }		///< Sets the coordinates as zero
	void Get( TYPE *p ) const { CY_MEMCOPY(TYPE,p,data,3); }	///< Puts the coordinate values into the array
	void Set( const TYPE *p ) { CY_MEMCOPY(TYPE,data,p,3); }	///< Sets the coordinates using the values in the given array
	void Set( const TYPE &_x, const TYPE &_y, const TYPE &_z ) { x=_x; y=_y; z=_z; }	///< Sets the coordinates using the given values

	///@name Length and Normalize functions
	TYPE   LengthSquared () const { Point3 p=operator*(*this); return p.Sum(); }	///< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	TYPE   Length        () const { return (TYPE) cySqrt(LengthSquared()); }			///< Returns the length of the vector.
	void   Normalize     ()       { *this /= Length(); }							///< Normalizes the vector, such that its length becomes 1.
	Point3 GetNormalized () const { return *this / Length(); }						///< Returns a normalized copy of the vector.
	TYPE   Sum           () const { return x+y+z; }									///< Returns the sum of its components
	bool   IsZero        () const { TYPE zero[3]={0,0,0}; return memcmp(data,zero,sizeof(data))==0; }	///< Returns true if all components are exactly zero
	TYPE   MaxComponent  () const { return x>y ? (x>z ? x : z) : (y>z ? y : z); }
	int    MaxComponentID() const { return x>y ? (x>z ? 0 : 2) : (y>z ? 1 : 2); }

	///@name Limit functions
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { for ( int i=0; i<3; ++i ) data[i] = (data[i]<v) ? v : data[i]; }
	void ClampMax( const TYPE &v ) { for ( int i=0; i<3; ++i ) data[i] = (data[i]>v) ? v : data[i]; }
	void Abs() { for ( int i=0; i<3; i++ ) data[i] = cyAbs(data[i]); }	///< Converts all negative components to positive values

	///@name Unary operators
	Point3 operator - () const { Point3 pt; for ( int i=0; i<3; ++i ) pt.data[i]=-data[i]; return pt; } 
	Point3 operator + () const { return *this; }

	///@name Binary operators
	Point3 operator + ( const Point3 &p ) const { Point3 r; for ( int i=0; i<3; ++i ) r.data[i] = data[i] + p.data[i]; return r; }
	Point3 operator - ( const Point3 &p ) const { Point3 r; for ( int i=0; i<3; ++i ) r.data[i] = data[i] - p.data[i]; return r; }
	Point3 operator * ( const Point3 &p ) const { Point3 r; for ( int i=0; i<3; ++i ) r.data[i] = data[i] * p.data[i]; return r; }
	Point3 operator / ( const Point3 &p ) const { Point3 r; for ( int i=0; i<3; ++i ) r.data[i] = data[i] / p.data[i]; return r; }
	Point3 operator + ( const TYPE   &v ) const { Point3 r; for ( int i=0; i<3; ++i ) r.data[i] = data[i] + v;         return r; }
	Point3 operator - ( const TYPE   &v ) const { Point3 r; for ( int i=0; i<3; ++i ) r.data[i] = data[i] - v;         return r; }
	Point3 operator * ( const TYPE   &v ) const { Point3 r; for ( int i=0; i<3; ++i ) r.data[i] = data[i] * v;         return r; }
	Point3 operator / ( const TYPE   &v ) const { Point3 r; for ( int i=0; i<3; ++i ) r.data[i] = data[i] / v;         return r; }

	///@name Assignment operators
	const Point3& operator  = ( const Point3 &p ) { CY_MEMCOPY(TYPE,data,p.data,3); return *this; }	
	const Point3& operator += ( const Point3 &p ) { for ( int i=0; i<3; ++i ) data[i] += p.data[i]; return *this; }
	const Point3& operator -= ( const Point3 &p ) { for ( int i=0; i<3; ++i ) data[i] -= p.data[i]; return *this; }
	const Point3& operator *= ( const Point3 &p ) { for ( int i=0; i<3; ++i ) data[i] *= p.data[i]; return *this; }
	const Point3& operator /= ( const Point3 &p ) { for ( int i=0; i<3; ++i ) data[i] /= p.data[i]; return *this; }
	const Point3& operator += ( const TYPE   &v ) { for ( int i=0; i<3; ++i ) data[i] += v;         return *this; }
	const Point3& operator -= ( const TYPE   &v ) { for ( int i=0; i<3; ++i ) data[i] -= v;         return *this; }
	const Point3& operator *= ( const TYPE   &v ) { for ( int i=0; i<3; ++i ) data[i] *= v;         return *this; }
	const Point3& operator /= ( const TYPE   &v ) { for ( int i=0; i<3; ++i ) data[i] /= v;         return *this; }

	///@name Test operators
	bool operator == ( const Point3& pt ) const { Point3 p=operator-(pt); return  p.IsZero(); }
	bool operator != ( const Point3& pt ) const { Point3 p=operator-(pt); return !p.IsZero(); }

	///@name Access operators
	TYPE& operator [] ( int i )       { return data[i]; }
	TYPE  operator [] ( int i ) const { return data[i]; }

	///@name Cross product and dot product
	Point3 Cross      ( const Point3 &pt ) const { return Point3(y*pt.z-z*pt.y, z*pt.x-x*pt.z, x*pt.y-y*pt.x); }	///< Cross product
	Point3 operator ^ ( const Point3 &pt ) const { return Cross(pt); }							///< Cross product
	TYPE   Dot        ( const Point3 &pt ) const { Point3 p=operator*(pt); return p.Sum(); }	///< Dot product
	TYPE   operator % ( const Point3 &pt ) const { return Dot(pt); }							///< Dot product

	///@name Conversion Methods
	Point2<TYPE> XY() const { return Point2<TYPE>(x,y); }
};

//-------------------------------------------------------------------------------

/// 4D point class
template <typename TYPE>
class Point4
{
	friend Point4 operator + ( const TYPE v, const Point4 &pt ) { return   pt+v;  }	///< Addition with a constant
	friend Point4 operator - ( const TYPE v, const Point4 &pt ) { return -(pt-v); }	///< Subtraction from a constant
	friend Point4 operator * ( const TYPE v, const Point4 &pt ) { return   pt*v;  }	///< Multiplication with a constant

public:

	union {
		struct {
			TYPE x, y, z, w;	///< Components of the vector
		};
		TYPE data[4];			///< Components of the vector in array form
	};

	///@name Constructors
	Point4() { }
	Point4( const TYPE &_x, const TYPE &_y, const TYPE &_z, const TYPE &_w ) : x(_x), y(_y), z(_z), w(_w) {}
	Point4( const Point4 &pt )        { CY_MEMCOPY(TYPE,data,pt.data,4); }
	explicit Point4( const TYPE *pt ) { CY_MEMCOPY(TYPE,data,pt,4); }
	explicit Point4( const TYPE &v  ) { for ( int i=0; i<4; ++i ) data[i]=v; }
	explicit Point4( const Point3<TYPE> &pt,            TYPE _w=1 ) { CY_MEMCOPY(TYPE,data,pt.data,3);       w=_w; }
	explicit Point4( const Point2<TYPE> &pt, TYPE _z=0, TYPE _w=1 ) { CY_MEMCOPY(TYPE,data,pt.data,2); z=_z; w=_w; }
	template <typename T> explicit Point4( const Point4<T> &pt )                       { CY_MEMCONVERT(TYPE,data,pt.data,4); }
	template <typename T> explicit Point4( const Point3<T> &pt,            TYPE _w=1 ) { CY_MEMCONVERT(TYPE,data,pt.data,3);       w=_w; }
	template <typename T> explicit Point4( const Point2<T> &pt, TYPE _z=0, TYPE _w=1 ) { CY_MEMCONVERT(TYPE,data,pt.data,2); z=_z; w=_w; }
	template <int M> explicit Point4( const Point<TYPE,M> &pt )
	{
		if ( 4 <= M ) CY_MEMCOPY(TYPE,data,pt.data,4);
		else {        CY_MEMCOPY(TYPE,data,pt.data,M); CY_MEMCLEAR(TYPE,data,4-M); }
	}
	template <typename T, int M> explicit Point4( const Point<T,M> &pt )
	{
		if ( 4 <= M ) CY_MEMCONVERT(TYPE,data,pt.data,4);
		else {        CY_MEMCONVERT(TYPE,data,pt.data,M); CY_MEMCLEAR(TYPE,data,4-M); }
	}

	///@name Set & Get value functions
	void Zero()               { CY_MEMCLEAR(TYPE,data,4); }		///< Sets the coordinates as zero
	void Get( TYPE *p ) const { CY_MEMCOPY(TYPE,p,data,4); }	///< Puts the coordinate values into the array
	void Set( const TYPE *p ) { CY_MEMCOPY(TYPE,data,p,4); }	///< Sets the coordinates using the values in the given array
	void Set( const TYPE &_x, const TYPE &_y, const TYPE &_z, const TYPE &_w=1 ) { x=_x; y=_y; z=_z; w=_w; }	///< Sets the coordinates using the given values

	///@name Length and Normalize functions
	TYPE   LengthSquared () const { Point4 p=operator*(*this); return p.Sum(); }	///< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	TYPE   Length        () const { return (TYPE) cySqrt(LengthSquared()); }			///< Returns the length of the vector.
	void   Normalize     ()       { *this /= Length(); }							///< Normalizes the vector, such that its length becomes 1.
	Point4 GetNormalized () const { return *this / Length(); }						///< Returns a normalized copy of the vector.
	TYPE   Sum           () const { return x+y+z+w; }								///< Returns the sum of its components
	bool   IsZero        () const { TYPE zero[4]={0,0,0,0}; return memcmp(data,zero,sizeof(data))==0; }	///< Returns true if all components are exactly zero
	TYPE   MaxComponent  () const { TYPE mxy = x>y ? x : y; TYPE mzw = z>w ? z : w; return mxy>mzw ? mxy : mzw; }
	int    MaxComponentID() const { int  ixy = x>y ? 0 : 1; int  izw = z>w ? 2 : 3; return data[ixy]>data[izw] ? ixy : izw; }

	///@name Limit functions
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { for ( int i=0; i<4; ++i ) data[i] = (data[i]<v) ? v : data[i]; }
	void ClampMax( const TYPE &v ) { for ( int i=0; i<4; ++i ) data[i] = (data[i]>v) ? v : data[i]; }
	void Abs() { for ( int i=0; i<4; i++ ) data[i] = cyAbs(data[i]); }	///< Converts all negative components to positive values

	///@name Unary operators
	Point4 operator - () const { Point4 pt; for ( int i=0; i<4; ++i ) pt.data[i]=-data[i]; return pt; } 
	Point4 operator + () const { return *this; }

	///@name Binary operators
	Point4 operator + ( const Point4 &p ) const { Point4 r; for ( int i=0; i<4; ++i ) r.data[i] = data[i] + p.data[i]; return r; }
	Point4 operator - ( const Point4 &p ) const { Point4 r; for ( int i=0; i<4; ++i ) r.data[i] = data[i] - p.data[i]; return r; }
	Point4 operator * ( const Point4 &p ) const { Point4 r; for ( int i=0; i<4; ++i ) r.data[i] = data[i] * p.data[i]; return r; }
	Point4 operator / ( const Point4 &p ) const { Point4 r; for ( int i=0; i<4; ++i ) r.data[i] = data[i] / p.data[i]; return r; }
	Point4 operator + ( const TYPE   &v ) const { Point4 r; for ( int i=0; i<4; ++i ) r.data[i] = data[i] + v;         return r; }
	Point4 operator - ( const TYPE   &v ) const { Point4 r; for ( int i=0; i<4; ++i ) r.data[i] = data[i] - v;         return r; }
	Point4 operator * ( const TYPE   &v ) const { Point4 r; for ( int i=0; i<4; ++i ) r.data[i] = data[i] * v;         return r; }
	Point4 operator / ( const TYPE   &v ) const { Point4 r; for ( int i=0; i<4; ++i ) r.data[i] = data[i] / v;         return r; }

	///@name Assignment operators
	const Point4& operator  = ( const Point4 &p ) { CY_MEMCOPY(TYPE,data,p.data,4); return *this; }	
	const Point4& operator += ( const Point4 &p ) { for ( int i=0; i<4; ++i ) data[i] += p.data[i]; return *this; }
	const Point4& operator -= ( const Point4 &p ) { for ( int i=0; i<4; ++i ) data[i] -= p.data[i]; return *this; }
	const Point4& operator *= ( const Point4 &p ) { for ( int i=0; i<4; ++i ) data[i] *= p.data[i]; return *this; }
	const Point4& operator /= ( const Point4 &p ) { for ( int i=0; i<4; ++i ) data[i] /= p.data[i]; return *this; }
	const Point4& operator += ( const TYPE   &v ) { for ( int i=0; i<4; ++i ) data[i] += v;         return *this; }
	const Point4& operator -= ( const TYPE   &v ) { for ( int i=0; i<4; ++i ) data[i] -= v;         return *this; }
	const Point4& operator *= ( const TYPE   &v ) { for ( int i=0; i<4; ++i ) data[i] *= v;         return *this; }
	const Point4& operator /= ( const TYPE   &v ) { for ( int i=0; i<4; ++i ) data[i] /= v;         return *this; }

	///@name Test operators
	bool operator == ( const Point4& pt ) const { Point4 p=operator-(pt); return  p.IsZero(); }
	bool operator != ( const Point4& pt ) const { Point4 p=operator-(pt); return !p.IsZero(); }

	///@name Access operators
	TYPE& operator [] ( int i )       { return data[i]; }
	TYPE  operator [] ( int i ) const { return data[i]; }

	///@ Dot product
	TYPE Dot		( const Point4 &pt ) const { Point4 p=operator*(pt); return p.Sum(); }	///< Dot product
	TYPE operator % ( const Point4 &pt ) const { return Dot(pt); }							///< Dot product

	///@name Conversion Methods
	Point2<TYPE> XY () const { return Point2<TYPE>(data); }
	Point3<TYPE> XYZ() const { return Point3<TYPE>(data); }
	Point3<TYPE> GetNonHomogeneous() const { return Point3<TYPE>(data)/w; }
};

//-------------------------------------------------------------------------------

// Definitions of the conversion constructors
template <typename TYPE, int N> Point<TYPE,N>::Point( const Point2<TYPE> &pt ) { if ( N <= 2 ) CY_MEMCOPY(TYPE,data,pt.data,N); else { CY_MEMCOPY(TYPE,data,pt.data,2);CY_MEMCLEAR(TYPE,data,N-2); } }
template <typename TYPE, int N> Point<TYPE,N>::Point( const Point3<TYPE> &pt ) { if ( N <= 3 ) CY_MEMCOPY(TYPE,data,pt.data,N); else { CY_MEMCOPY(TYPE,data,pt.data,3);CY_MEMCLEAR(TYPE,data,N-3); } }
template <typename TYPE, int N> Point<TYPE,N>::Point( const Point4<TYPE> &pt ) { if ( N <= 4 ) CY_MEMCOPY(TYPE,data,pt.data,N); else { CY_MEMCOPY(TYPE,data,pt.data,4);CY_MEMCLEAR(TYPE,data,N-4); } }
template <typename TYPE, int N> template <typename T> Point<TYPE,N>::Point( const Point2<T> &pt ) { if ( N <= 2 ) CY_MEMCONVERT(TYPE,data,pt.data,N); else { CY_MEMCONVERT(TYPE,data,pt.data,2); CY_MEMCLEAR(TYPE,data,N-2); } }
template <typename TYPE, int N> template <typename T> Point<TYPE,N>::Point( const Point3<T> &pt ) { if ( N <= 3 ) CY_MEMCONVERT(TYPE,data,pt.data,N); else { CY_MEMCONVERT(TYPE,data,pt.data,3); CY_MEMCLEAR(TYPE,data,N-3); } }
template <typename TYPE, int N> template <typename T> Point<TYPE,N>::Point( const Point4<T> &pt ) { if ( N <= 4 ) CY_MEMCONVERT(TYPE,data,pt.data,N); else { CY_MEMCONVERT(TYPE,data,pt.data,4); CY_MEMCLEAR(TYPE,data,N-4); } }
template <typename TYPE> Point2<TYPE>::Point2( const Point3<TYPE> &pt ) { CY_MEMCOPY(TYPE,data,pt.data,2); }
template <typename TYPE> Point2<TYPE>::Point2( const Point4<TYPE> &pt ) { CY_MEMCOPY(TYPE,data,pt.data,2); }
template <typename TYPE> Point3<TYPE>::Point3( const Point4<TYPE> &pt ) { CY_MEMCOPY(TYPE,data,pt.data,3); }
template <typename TYPE> template <typename T> Point2<TYPE>::Point2( const Point3<T> &pt ) { CY_MEMCONVERT(TYPE,data,pt.data,2); }
template <typename TYPE> template <typename T> Point2<TYPE>::Point2( const Point4<T> &pt ) { CY_MEMCONVERT(TYPE,data,pt.data,2); }
template <typename TYPE> template <typename T> Point3<TYPE>::Point3( const Point4<T> &pt ) { CY_MEMCONVERT(TYPE,data,pt.data,3); }

//-------------------------------------------------------------------------------

typedef Point2<float>    Point2f;
typedef Point3<float>    Point3f;
typedef Point4<float>    Point4f;

typedef Point2<double>   Point2d;
typedef Point3<double>   Point3d;
typedef Point4<double>   Point4d;

typedef Point2<int32_t>  Point2i;
typedef Point3<int32_t>  Point3i;
typedef Point4<int32_t>  Point4i;

typedef Point2<uint32_t> Point2ui;
typedef Point3<uint32_t> Point3ui;
typedef Point4<uint32_t> Point4ui;

typedef Point2<int64_t>  Point2l;
typedef Point3<int64_t>  Point3l;
typedef Point4<int64_t>  Point4l;

typedef Point2<uint64_t> Point2ul;
typedef Point3<uint64_t> Point3ul;
typedef Point4<uint64_t> Point4ul;

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Point2f  cyPoint2f;
typedef cy::Point3f  cyPoint3f;
typedef cy::Point4f  cyPoint4f;

typedef cy::Point2d  cyPoint2d;
typedef cy::Point3d  cyPoint3d;
typedef cy::Point4d  cyPoint4d;

typedef cy::Point2i  cyPoint2i;
typedef cy::Point3i  cyPoint3i;
typedef cy::Point4i  cyPoint4i;

typedef cy::Point2ui cyPoint2ui;
typedef cy::Point3ui cyPoint3ui;
typedef cy::Point4ui cyPoint4ui;

typedef cy::Point2l  cyPoint2l;
typedef cy::Point3l  cyPoint3l;
typedef cy::Point4l  cyPoint4l;

typedef cy::Point2ul cyPoint2ul;
typedef cy::Point3ul cyPoint3ul;
typedef cy::Point4ul cyPoint4ul;

//-------------------------------------------------------------------------------

#endif

