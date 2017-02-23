// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyPoint.h 
//! \author Cem Yuksel
//! 
//! \brief  2D, 3D, 4D, and ND point classes.
//! 
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
//!	\cond HIDDEN_SYMBOLS
template <typename TYPE> class Point2;
template <typename TYPE> class Point3;
template <typename TYPE> class Point4;
//! \endcond

//-------------------------------------------------------------------------------

//! A general class for N-dimensional points (vectors).

template <typename TYPE, int N>
class Point
{
	friend Point operator + ( const TYPE v, const Point &p ) { return   p+v; }	//!< Addition with a constant
	friend Point operator - ( const TYPE v, const Point &p ) { return -(p-v); }	//!< Subtraction from a constant
	friend Point operator * ( const TYPE v, const Point &p ) { return   p*v; }	//!< Multiplication with a constant

public:

	//!@name Components of the point/vector
	TYPE data[N];

	//!@name Constructors
	Point() {}
	Point( const Point &p )         { CY_MEMCOPY(TYPE,data,p.data,N); }
	explicit Point( const TYPE *p ) { CY_MEMCOPY(TYPE,data,p,N); }
	explicit Point( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i]=v; }
	template <typename T> explicit Point( const Point<T,N> &p ) { CY_MEMCONVERT(TYPE,data,p.data,N); }
	template <int M> explicit Point( const Point<TYPE,M> &p )
	{
		if ( N <= M ) { CY_MEMCOPY(TYPE,data,p.data,N); }
		else          { CY_MEMCOPY(TYPE,data,p.data,M); CY_MEMCLEAR(TYPE,data,N-M); }
	}
	template <typename T, int M> explicit Point( const Point<T,M> &p )
	{
		if ( N <= M ) { CY_MEMCONVERT(TYPE,data,p.data,N); }
		else          { CY_MEMCONVERT(TYPE,data,p.data,M); CY_MEMCLEAR(TYPE,data,N-M); }
	}
	explicit Point( const Point2<TYPE> &p );
	explicit Point( const Point3<TYPE> &p );
	explicit Point( const Point4<TYPE> &p );
	template <typename T> explicit Point( const Point2<T> &p );
	template <typename T> explicit Point( const Point3<T> &p );
	template <typename T> explicit Point( const Point4<T> &p );
	template <typename P> explicit Point( const P &p ) { for ( int i=0; i<N; ++i ) data[i]=(TYPE)p[i]; }

	//!@name Set & Get value methods
	void Zero()               { CY_MEMCLEAR(TYPE,data,N); }					//!< Sets the coordinates as zero
	void Get( TYPE *p ) const { CY_MEMCOPY(TYPE,p,data,N); }				//!< Puts the coordinate values into the array
	void Set( const TYPE *p ) { CY_MEMCOPY(TYPE,data,p,N); }				//!< Sets the coordinates using the values in the given array
	void Set( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i] = v; }	//!< Sets all coordinates using the given value
	template <int M> void CopyData( TYPE *p ) { if ( M <= N ) { CY_MEMCOPY(TYPE,p,data,M); } else { CY_MEMCOPY(TYPE,p,data,N); CY_MEMCLEAR(TYPE,p+N,M-N); }	}
	template <typename T, int M> void ConvertData( T *p ) { if ( M <= N ) { CY_MEMCONVERT(T,p,data,M); } else { CY_MEMCONVERT(T,p,data,N); CY_MEMCLEAR(T,p+N,M-N); }	}

	//!@name General methods
	TYPE  LengthSquared() const { Point p=operator*(*this); return p.Sum(); }	//!< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	TYPE  Length       () const { return (TYPE) cySqrt(LengthSquared()); }		//!< Returns the length of the vector.
	void  Normalize    ()       { *this /= Length(); }							//!< Normalizes the vector, such that its length becomes 1.
	Point GetNormalized() const { return *this / Length(); }					//!< Returns a normalized copy of the vector.
	TYPE  Sum          () const { TYPE v=data[0]; for ( int i=1; i<N; ++i ) v+=data[i]; return v; }		//!< Returns the sum of its components
	bool  IsZero       () const { for ( int i=0; i<N; ++i ) if ( data[i] != TYPE(0) ) return false; return true; }	//!< Returns true if all components are exactly zero
	TYPE  Min          () const { TYPE m = data[0]; for ( int i=1; i<N; ++i ) if ( m > data[i] ) m = data[i]; return m; }
	TYPE  Max          () const { TYPE m = data[0]; for ( int i=1; i<N; ++i ) if ( m < data[i] ) m = data[i]; return m; }
	int   MinID        () const { TYPE m = data[0]; int ix=0; for ( int i=1; i<N; ++i ) if ( m > data[i] ) { m = data[i]; ix = i; } return m; }
	int   MaxID        () const { TYPE m = data[0]; int ix=0; for ( int i=1; i<N; ++i ) if ( m < data[i] ) { m = data[i]; ix = i; } return m; }

	//!@name Limit methods
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i] = (data[i]<v) ? v : data[i]; }
	void ClampMax( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i] = (data[i]>v) ? v : data[i]; }
	void Abs() { for ( int i=0; i<N; i++ ) data[i] = cyAbs(data[i]); }	//!< Converts all negative components to positive values

	//!@name Unary operators
	Point operator - () const { Point r; for ( int i=0; i<N; ++i ) r.data[i]=-data[i]; return r; } 

	//!@name Binary operators
	Point operator + ( const Point &p ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] + p.data[i]; return r; }
	Point operator - ( const Point &p ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] - p.data[i]; return r; }
	Point operator * ( const Point &p ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] * p.data[i]; return r; }
	Point operator / ( const Point &p ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] / p.data[i]; return r; }
	Point operator + ( const TYPE  &v ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] + v; return r; }
	Point operator - ( const TYPE  &v ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] - v; return r; }
	Point operator * ( const TYPE  &v ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] * v; return r; }
	Point operator / ( const TYPE  &v ) const { Point r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] / v; return r; }

	//!@name Assignment operators
	const Point& operator  = ( const Point &p ) { CY_MEMCOPY(TYPE,data,p.data,N); return *this; }	
	const Point& operator += ( const Point &p ) { for ( int i=0; i<N; ++i ) data[i] += p.data[i]; return *this; }
	const Point& operator -= ( const Point &p ) { for ( int i=0; i<N; ++i ) data[i] -= p.data[i]; return *this; }
	const Point& operator *= ( const Point &p ) { for ( int i=0; i<N; ++i ) data[i] *= p.data[i]; return *this; }
	const Point& operator /= ( const Point &p ) { for ( int i=0; i<N; ++i ) data[i] /= p.data[i]; return *this; }
	const Point& operator += ( const TYPE  &v ) { for ( int i=0; i<N; ++i ) data[i] += v; return *this; }
	const Point& operator -= ( const TYPE  &v ) { for ( int i=0; i<N; ++i ) data[i] -= v; return *this; }
	const Point& operator *= ( const TYPE  &v ) { for ( int i=0; i<N; ++i ) data[i] *= v; return *this; }
	const Point& operator /= ( const TYPE  &v ) { for ( int i=0; i<N; ++i ) data[i] /= v; return *this; }

	//!@name Test operators
	bool operator == ( const Point& p ) const { for ( int i=0; i<N; ++i ) if ( data[i] != p.data[i] ) return false; return true; }
	bool operator != ( const Point& p ) const { for ( int i=0; i<N; ++i ) if ( data[i] != p.data[i] ) return true; return false; }

	//!@name Access operators
	TYPE& operator [] ( int i )       { return data[i]; }
	TYPE  operator [] ( int i ) const { return data[i]; }

	//!@name Dot product
	TYPE Dot        ( const Point &p ) const { Point r=operator*(p); return r.Sum(); }	//!< Dot product
	TYPE operator % ( const Point &p ) const { return Dot(p); }						//!< Dot product operator
};

//-------------------------------------------------------------------------------

//! 2D point (vector) class

template <typename TYPE>
class Point2
{
	friend Point2 operator + ( const TYPE v, const Point2 &p ) { return   p+v;  }	//!< Addition with a constant
	friend Point2 operator - ( const TYPE v, const Point2 &p ) { return -(p-v); }	//!< Subtraction from a constant
	friend Point2 operator * ( const TYPE v, const Point2 &p ) { return   p*v;  }	//!< Multiplication with a constant

public:

	//!@name Components of the point/vector
	TYPE x, y;

	//!@name Constructors
	Point2() {}

	Point2( const TYPE &_x, const TYPE &_y ) : x( _x), y( _y) {}
	Point2( const Point2 &p )                : x(p.x), y(p.y) {}
	explicit Point2( const TYPE &v )         : x(v  ), y(v  ) {}
	explicit Point2( const Point3<TYPE> &p );
	explicit Point2( const Point4<TYPE> &p );
	template <typename T> explicit Point2( const Point2<T> &p ) : x(TYPE(p.x)), y(TYPE(p.y)) {}
	template <typename T> explicit Point2( const Point3<T> &p );
	template <typename T> explicit Point2( const Point4<T> &p );
	template <int M> explicit Point2( const Point<TYPE,M> &p ) { p.CopyData<2>(&x); }
	template <typename T, int M> explicit Point2( const Point<T,M> &p ) { p.template ConvertData<TYPE,2>(&x); }
	template <typename P> explicit Point2( const P &p ) : x((TYPE)p[0]), y((TYPE)p[1]) {}

	//!@name Set & Get value methods
	void Zero()               { CY_MEMCLEAR(TYPE,Data(),2); }		//!< Sets the coordinates as zero.
	void Get( TYPE *p ) const { ((Point2*)p)->operator=(*this); }	//!< Puts the coordinate values into the array.
	void Set( const TYPE *p ) { operator=(*((Point2*)p)); }			//!< Sets the coordinates using the values in the given array.
	void Set( const TYPE &v ) { x=v; y=v; }							//!< Sets all coordinates using the given value
	void Set( const TYPE &_x, const TYPE &_y ) { x=_x; y=_y; }		//!< Sets the coordinates using the given values

	//!@name General methods
	TYPE   LengthSquared() const { Point2 p=operator*(*this); return p.Sum(); }	//!< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	TYPE   Length       () const { return (TYPE) cySqrt(LengthSquared()); }		//!< Returns the length of the vector.
	void   Normalize    ()       { *this /= Length(); }							//!< Normalizes the vector, such that its length becomes 1.
	Point2 GetNormalized() const { return *this / Length(); }					//!< Returns a normalized copy of the vector.
	TYPE   Sum          () const { return x+y; }								//!< Returns the sum of its components
	bool   IsZero       () const { return x==TYPE(0) && y==TYPE(0); }			//!< Returns true if all components are exactly zero
	TYPE   Min          () const { return x<y ? x : y; }
	TYPE   Max          () const { return x>y ? x : y; }
	int    MinID        () const { return x<y ? 0 : 1; }
	int    MaxID        () const { return x>y ? 0 : 1; }

	//!@name Limit methods
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { x=(x<v)?v:x; y=(y<v)?v:y; }
	void ClampMax( const TYPE &v ) { x=(x>v)?v:x; y=(y>v)?v:y; }
	void Abs() { x=cyAbs(x); y=cyAbs(y); }	//!< Converts all negative components to positive values

	//!@name Unary operators
	Point2 operator - () const { Point2 r; r.x=-x; r.y=-y; return r; } 

	//!@name Binary operators
	Point2 operator + ( const Point2 &p ) const { Point2 r; r.x=x+p.x; r.y=y+p.y; return r; }
	Point2 operator - ( const Point2 &p ) const { Point2 r; r.x=x-p.x; r.y=y-p.y; return r; }
	Point2 operator * ( const Point2 &p ) const { Point2 r; r.x=x*p.x; r.y=y*p.y; return r; }
	Point2 operator / ( const Point2 &p ) const { Point2 r; r.x=x/p.x; r.y=y/p.y; return r; }
	Point2 operator + ( const TYPE   &v ) const { Point2 r; r.x=x+v;   r.y=y+v;   return r; }
	Point2 operator - ( const TYPE   &v ) const { Point2 r; r.x=x-v;   r.y=y-v;   return r; }
	Point2 operator * ( const TYPE   &v ) const { Point2 r; r.x=x*v;   r.y=y*v;   return r; }
	Point2 operator / ( const TYPE   &v ) const { Point2 r; r.x=x/v;   r.y=y/v;   return r; }

	//!@name Assignment operators
	const Point2& operator  = ( const Point2 &p ) { x =p.x; y =p.y; return *this; }	
	const Point2& operator += ( const Point2 &p ) { x+=p.x; y+=p.y; return *this; }
	const Point2& operator -= ( const Point2 &p ) { x-=p.x; y-=p.y; return *this; }
	const Point2& operator *= ( const Point2 &p ) { x*=p.x; y*=p.y; return *this; }
	const Point2& operator /= ( const Point2 &p ) { x/=p.x; y/=p.y; return *this; }
	const Point2& operator += ( const TYPE   &v ) { x+=v;   y+=v;   return *this; }
	const Point2& operator -= ( const TYPE   &v ) { x-=v;   y-=v;   return *this; }
	const Point2& operator *= ( const TYPE   &v ) { x*=v;   y*=v;   return *this; }
	const Point2& operator /= ( const TYPE   &v ) { x/=v;   y/=v;   return *this; }

	//!@name Test operators
	bool operator == ( const Point2& p ) const { return x==p.x && y==p.y; }
	bool operator != ( const Point2& p ) const { return x!=p.x && y!=p.y; }

	//!@name Access operators
	TYPE&       operator [] ( int i )       { return Element(i); }
	const TYPE& operator [] ( int i ) const { return Element(i); }
	TYPE&       Element     ( int i )       { return (&x)[i]; }
	const TYPE& Element     ( int i ) const { return (&x)[i]; }
	TYPE*       Data        ()              { return &x; }
	const TYPE* Data        ()        const { return &x; }

	//!@name Cross product and dot product
	TYPE Cross      ( const Point2 &p ) const { Point2 r(-y,x); return r.Dot(p); }			//!< Cross product
	TYPE operator ^ ( const Point2 &p ) const { return Cross(p); }							//!< Cross product operator
	TYPE Dot        ( const Point2 &p ) const { Point2 r=operator*(p); return r.Sum(); }	//!< Dot product
	TYPE operator % ( const Point2 &p ) const { return Dot(p); }							//!< Dot product operator
};

//-------------------------------------------------------------------------------

//! 3D point (vector) class

template <typename TYPE>
class Point3
{
	friend Point3 operator + ( const TYPE v, const Point3 &p ) { return   p+v;  }	//!< Addition with a constant
	friend Point3 operator - ( const TYPE v, const Point3 &p ) { return -(p-v); }	//!< Subtraction from a constant
	friend Point3 operator * ( const TYPE v, const Point3 &p ) { return   p*v;  }	//!< Multiplication with a constant

public:

	//!@name Components of the point/vector
	TYPE x, y, z;

	//!@name Constructors
	Point3() { }
	Point3( const TYPE &_x, const TYPE &_y, const TYPE &_z ) : x( _x), y( _y), z( _z) {}
	Point3( const Point3 &p )                                : x(p.x), y(p.y), z(p.z) {}
	explicit Point3( const TYPE &v )                         : x(v  ), y(v  ), z(v  ) {}
	explicit Point3( const Point2<TYPE> &p, TYPE _z=0 )      : x(p.x), y(p.y), z( _z) {}
	explicit Point3( const Point4<TYPE> &p );
	template <typename T> explicit Point3( const Point3<T> &p )            : x(TYPE(p.x)), y(TYPE(p.y)), z(TYPE(p.z)) {}
	template <typename T> explicit Point3( const Point2<T> &p, TYPE _z=0 ) : x(TYPE(p.x)), y(TYPE(p.y)), z(      _z  ) {}
	template <typename T> explicit Point3( const Point4<T> &p );
	template <int M> explicit Point3( const Point<TYPE,M> &p ) { p.CopyData<3>(&x); }
	template <typename T, int M> explicit Point3( const Point<T,M> &p ) { p.template ConvertData<TYPE,3>(&x); }
	template <typename P> explicit Point3( const P &p ) : x((TYPE)p[0]), y((TYPE)p[1]), z((TYPE)p[2]) {}

	//!@name Set & Get value methods
	void Zero()               { CY_MEMCLEAR(TYPE,Data(),3); }		//!< Sets the coordinates as zero
	void Get( TYPE *p ) const { ((Point3*)p)->operator=(*this); }	//!< Puts the coordinate values into the array
	void Set( const TYPE *p ) { operator=(*((Point3*)p)); }			//!< Sets the coordinates using the values in the given array
	void Set( const TYPE &v ) { x=v; y=v; z=v; }					//!< Sets all coordinates using the given value
	void Set( const TYPE &_x, const TYPE &_y, const TYPE &_z ) { x=_x; y=_y; z=_z; }	//!< Sets the coordinates using the given values

	//!@name General methods
	TYPE   LengthSquared() const { Point3 p=operator*(*this); return p.Sum(); }		//!< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	TYPE   Length       () const { return (TYPE) cySqrt(LengthSquared()); }			//!< Returns the length of the vector.
	void   Normalize    ()       { *this /= Length(); }								//!< Normalizes the vector, such that its length becomes 1.
	Point3 GetNormalized() const { return *this / Length(); }						//!< Returns a normalized copy of the vector.
	TYPE   Sum          () const { return x+y+z; }									//!< Returns the sum of its components
	bool   IsZero       () const { return x==TYPE(0) && y==TYPE(0) && z==TYPE(0); }	//!< Returns true if all components are exactly zero
	TYPE   Min          () const { return x<y ? (x<z ? x : z) : (y<z ? y : z); }
	TYPE   Max          () const { return x>y ? (x>z ? x : z) : (y>z ? y : z); }
	int    MinID        () const { return x<y ? (x<z ? 0 : 2) : (y<z ? 1 : 2); }
	int    MaxID        () const { return x>y ? (x>z ? 0 : 2) : (y>z ? 1 : 2); }

	//!@name Limit methods
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { x=(x<v)?v:x; y=(y<v)?v:y; z=(z<v)?v:z; }
	void ClampMax( const TYPE &v ) { x=(x>v)?v:x; y=(y>v)?v:y; z=(z>v)?v:z; }
	void Abs() { x=cyAbs(x); y=cyAbs(y); z=cyAbs(z); }	//!< Converts all negative components to positive values

	//!@name Unary operators
	Point3 operator - () const { Point3 r; r.x=-x; r.y=-y; r.z=-z; return r; } 

	//!@name Binary operators
	Point3 operator + ( const Point3 &p ) const { Point3 r; r.x=x+p.x; r.y=y+p.y; r.z=z+p.z; return r; }
	Point3 operator - ( const Point3 &p ) const { Point3 r; r.x=x-p.x; r.y=y-p.y; r.z=z-p.z; return r; }
	Point3 operator * ( const Point3 &p ) const { Point3 r; r.x=x*p.x; r.y=y*p.y; r.z=z*p.z; return r; }
	Point3 operator / ( const Point3 &p ) const { Point3 r; r.x=x/p.x; r.y=y/p.y; r.z=z/p.z; return r; }
	Point3 operator + ( const TYPE   &v ) const { Point3 r; r.x=x+v;   r.y=y+v;   r.z=z+v;   return r; }
	Point3 operator - ( const TYPE   &v ) const { Point3 r; r.x=x-v;   r.y=y-v;   r.z=z-v;   return r; }
	Point3 operator * ( const TYPE   &v ) const { Point3 r; r.x=x*v;   r.y=y*v;   r.z=z*v;   return r; }
	Point3 operator / ( const TYPE   &v ) const { Point3 r; r.x=x/v;   r.y=y/v;   r.z=z/v;   return r; }

	//!@name Assignment operators
	const Point3& operator  = ( const Point3 &p ) { x =p.x; y =p.y; z =p.z; return *this; }	
	const Point3& operator += ( const Point3 &p ) { x+=p.x; y+=p.y; z+=p.z; return *this; }
	const Point3& operator -= ( const Point3 &p ) { x-=p.x; y-=p.y; z-=p.z; return *this; }
	const Point3& operator *= ( const Point3 &p ) { x*=p.x; y*=p.y; z*=p.z; return *this; }
	const Point3& operator /= ( const Point3 &p ) { x/=p.x; y/=p.y; z/=p.z; return *this; }
	const Point3& operator += ( const TYPE   &v ) { x+=v;   y+=v;   z+=v;   return *this; }
	const Point3& operator -= ( const TYPE   &v ) { x-=v;   y-=v;   z-=v;   return *this; }
	const Point3& operator *= ( const TYPE   &v ) { x*=v;   y*=v;   z*=v;   return *this; }
	const Point3& operator /= ( const TYPE   &v ) { x/=v;   y/=v;   z/=v;   return *this; }

	//!@name Test operators
	bool operator == ( const Point3& p ) const { return x==p.x && y==p.y && z==p.z; }
	bool operator != ( const Point3& p ) const { return x!=p.x && y!=p.y && z!=p.z; }

	//!@name Access operators
	TYPE&       operator [] ( int i )       { return Element(i); }
	const TYPE& operator [] ( int i ) const { return Element(i); }
	TYPE&       Element     ( int i )       { return (&x)[i]; }
	const TYPE& Element     ( int i ) const { return (&x)[i]; }
	TYPE*       Data        ()              { return &x; }
	const TYPE* Data        ()        const { return &x; }

	//!@name Cross product and dot product
	Point3 Cross      ( const Point3 &p ) const { return Point3(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x); }	//!< Cross product
	Point3 operator ^ ( const Point3 &p ) const { return Cross(p); }										//!< Cross product
	TYPE   Dot        ( const Point3 &p ) const { Point3 r=operator*(p); return r.Sum(); }					//!< Dot product
	TYPE   operator % ( const Point3 &p ) const { return Dot(p); }											//!< Dot product

	//!@name Conversion Methods
	Point2<TYPE> XY() const { return Point2<TYPE>(*this); }
};

//-------------------------------------------------------------------------------

//! 4D point (vector) class

template <typename TYPE>
class Point4
{
	friend Point4 operator + ( const TYPE v, const Point4 &p ) { return   p+v;  }	//!< Addition with a constant
	friend Point4 operator - ( const TYPE v, const Point4 &p ) { return -(p-v); }	//!< Subtraction from a constant
	friend Point4 operator * ( const TYPE v, const Point4 &p ) { return   p*v;  }	//!< Multiplication with a constant

public:

	//!@name Components of the point/vector
	TYPE x, y, z, w;

	//!@name Constructors
	Point4() { }
	Point4( const TYPE &_x, const TYPE &_y, const TYPE &_z, const TYPE &_w ) : x( _x), y( _y), z( _z), w( _w) {}
	Point4( const Point4 &p )                                                : x(p.x), y(p.y), z(p.z), w(p.w) {}
	explicit Point4( const TYPE &v )                                         : x(v  ), y(v  ), z(v  ), w(v  ) {}
	explicit Point4( const Point3<TYPE> &p,            TYPE _w=1 )           : x(p.x), y(p.y), z(p.z), w( _w) {}
	explicit Point4( const Point2<TYPE> &p, TYPE _z=0, TYPE _w=1 )           : x(p.x), y(p.y), z( _z), w( _w) {}
	template <typename T> explicit Point4( const Point4<T> &p )                       : x(TYPE(p.x)), y(TYPE(p.y)), z(TYPE(p.z)), w(TYPE(p.w)) {}
	template <typename T> explicit Point4( const Point3<T> &p,            TYPE _w=1 ) : x(TYPE(p.x)), y(TYPE(p.y)), z(TYPE(p.z)), w(      _w ) {}
	template <typename T> explicit Point4( const Point2<T> &p, TYPE _z=0, TYPE _w=1 ) : x(TYPE(p.x)), y(TYPE(p.y)), z(      _z ), w(      _w ) {}
	template <int M> explicit Point4( const Point<TYPE,M> &p ) { p.CopyData<4>(&x); }
	template <typename T, int M> explicit Point4( const Point<T,M> &p ) { p.template ConvertData<TYPE,4>(&x); }
	template <typename P> explicit Point4( const P &p ) : x((TYPE)p[0]), y((TYPE)p[1]), z((TYPE)p[2]), w((TYPE)p[3]) {}

	//!@name Set & Get value methods
	void Zero()               { CY_MEMCLEAR(TYPE,Data(),4); }		//!< Sets the coordinates as zero
	void Get( TYPE *p ) const { ((Point4*)p)->operator=(*this); }	//!< Puts the coordinate values into the array
	void Set( const TYPE *p ) { operator=(*((Point4*)p)); }			//!< Sets the coordinates using the values in the given array
	void Set( const TYPE &v ) { x=v; y=v; z=v; w=v; }				//!< Sets all coordinates using the given value
	void Set( const TYPE &_x, const TYPE &_y, const TYPE &_z, const TYPE &_w=1 ) { x=_x; y=_y; z=_z; w=_w; }	//!< Sets the coordinates using the given values

	//!@name General methods
	TYPE   LengthSquared() const { Point4 p=operator*(*this); return p.Sum(); }	//!< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	TYPE   Length       () const { return (TYPE) cySqrt(LengthSquared()); }		//!< Returns the length of the vector.
	void   Normalize    ()       { *this /= Length(); }							//!< Normalizes the vector, such that its length becomes 1.
	Point4 GetNormalized() const { return *this / Length(); }					//!< Returns a normalized copy of the vector.
	TYPE   Sum          () const { return x+y+z+w; }							//!< Returns the sum of its components
	bool   IsZero       () const { return x==TYPE(0) && y==TYPE(0) && z==TYPE(0) && w==TYPE(0); }	//!< Returns true if all components are exactly zero
	TYPE   Min          () const { TYPE mxy = x<y ? x : y; TYPE mzw = z<w ? z : w; return mxy<mzw ? mxy : mzw; }
	TYPE   Max          () const { TYPE mxy = x>y ? x : y; TYPE mzw = z>w ? z : w; return mxy>mzw ? mxy : mzw; }
	int    MinID        () const { int  ixy = x<y ? 0 : 1; int  izw = z<w ? 2 : 3; return (&x)[ixy]<(&x)[izw] ? ixy : izw; }
	int    MaxID        () const { int  ixy = x>y ? 0 : 1; int  izw = z>w ? 2 : 3; return (&x)[ixy]>(&x)[izw] ? ixy : izw; }

	//!@name Limit methods
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { x=(x<v)?v:x; y=(y<v)?v:y; z=(z<v)?v:z; w=(w<v)?v:w; }
	void ClampMax( const TYPE &v ) { x=(x>v)?v:x; y=(y>v)?v:y; z=(z>v)?v:z; w=(w>v)?v:w; }
	void Abs() { x=cyAbs(x); y=cyAbs(y); z=cyAbs(z); w=cyAbs(w); }	//!< Converts all negative components to positive values

	//!@name Unary operators
	Point4 operator - () const { Point4 r; r.x=-x; r.y=-y; r.z=-z; r.w=-w; return r; } 

	//!@name Binary operators
	Point4 operator + ( const Point4 &p ) const { Point4 r; r.x=x+p.x; r.y=y+p.y; r.z=z+p.z; r.w=w+p.w; return r; }
	Point4 operator - ( const Point4 &p ) const { Point4 r; r.x=x-p.x; r.y=y-p.y; r.z=z-p.z; r.w=w-p.w; return r; }
	Point4 operator * ( const Point4 &p ) const { Point4 r; r.x=x*p.x; r.y=y*p.y; r.z=z*p.z; r.w=w*p.w; return r; }
	Point4 operator / ( const Point4 &p ) const { Point4 r; r.x=x/p.x; r.y=y/p.y; r.z=z/p.z; r.w=w/p.w; return r; }
	Point4 operator + ( const TYPE   &v ) const { Point4 r; r.x=x+v;   r.y=y+v;   r.z=z+v;   r.w=w+v;   return r; }
	Point4 operator - ( const TYPE   &v ) const { Point4 r; r.x=x-v;   r.y=y-v;   r.z=z-v;   r.w=w-v;   return r; }
	Point4 operator * ( const TYPE   &v ) const { Point4 r; r.x=x*v;   r.y=y*v;   r.z=z*v;   r.w=w*v;   return r; }
	Point4 operator / ( const TYPE   &v ) const { Point4 r; r.x=x/v;   r.y=y/v;   r.z=z/v;   r.w=w/v;   return r; }

	//!@name Assignment operators
	const Point4& operator  = ( const Point4 &p ) { x =p.x; y =p.y; z =p.z; w =p.w; return *this; }	
	const Point4& operator += ( const Point4 &p ) { x+=p.x; y+=p.y; z+=p.z; w+=p.w; return *this; }
	const Point4& operator -= ( const Point4 &p ) { x-=p.x; y-=p.y; z-=p.z; w-=p.w; return *this; }
	const Point4& operator *= ( const Point4 &p ) { x*=p.x; y*=p.y; z*=p.z; w*=p.w; return *this; }
	const Point4& operator /= ( const Point4 &p ) { x/=p.x; y/=p.y; z/=p.z; w/=p.w; return *this; }
	const Point4& operator += ( const TYPE   &v ) { x+=v;   y+=v;   z+=v;   w+=v;   return *this; }
	const Point4& operator -= ( const TYPE   &v ) { x-=v;   y-=v;   z-=v;   w-=v;   return *this; }
	const Point4& operator *= ( const TYPE   &v ) { x*=v;   y*=v;   z*=v;   w*=v;   return *this; }
	const Point4& operator /= ( const TYPE   &v ) { x/=v;   y/=v;   z/=v;   w/=v;   return *this; }

	//!@name Test operators
	bool operator == ( const Point4& p ) const { return x==p.x && y==p.y && z==p.z && w==p.w; }
	bool operator != ( const Point4& p ) const { return x!=p.x && y!=p.y && z!=p.z && w!=p.w; }

	//!@name Access operators
	TYPE&       operator [] ( int i )       { return Element(i); }
	const TYPE& operator [] ( int i ) const { return Element(i); }
	TYPE&       Element     ( int i )       { return (&x)[i]; }
	const TYPE& Element     ( int i ) const { return (&x)[i]; }
	TYPE*       Data        ()              { return &x; }
	const TYPE* Data        ()        const { return &x; }

	//!@name Dot product
	TYPE Dot		( const Point4 &p ) const { Point4 r=operator*(p); return r.Sum(); }	//!< Dot product
	TYPE operator % ( const Point4 &p ) const { return Dot(p); }							//!< Dot product

	//!@name Conversion Methods
	Point2<TYPE> XY () const { return Point2<TYPE>(*this); }
	Point3<TYPE> XYZ() const { return Point3<TYPE>(*this); }
	Point3<TYPE> GetNonHomogeneous() const { return Point3<TYPE>(*this)/w; }
};

//-------------------------------------------------------------------------------

// Definitions of the conversion constructors
template <typename TYPE, int N> Point<TYPE,N>::Point( const Point2<TYPE> &p ) { if ( N <= 2 ) { CY_MEMCOPY(TYPE,data,&p.x,N); } else { CY_MEMCOPY(TYPE,data,&p.x,2); CY_MEMCLEAR(TYPE,data,N-2); } }
template <typename TYPE, int N> Point<TYPE,N>::Point( const Point3<TYPE> &p ) { if ( N <= 3 ) { CY_MEMCOPY(TYPE,data,&p.x,N); } else { CY_MEMCOPY(TYPE,data,&p.x,3); CY_MEMCLEAR(TYPE,data,N-3); } }
template <typename TYPE, int N> Point<TYPE,N>::Point( const Point4<TYPE> &p ) { if ( N <= 4 ) { CY_MEMCOPY(TYPE,data,&p.x,N); } else { CY_MEMCOPY(TYPE,data,&p.x,4); CY_MEMCLEAR(TYPE,data,N-4); } }
template <typename TYPE, int N> template <typename T> Point<TYPE,N>::Point( const Point2<T> &p ) { if ( N <= 2 ) { CY_MEMCONVERT(TYPE,data,&p.x,N); } else { CY_MEMCONVERT(TYPE,data,&p.x,2); CY_MEMCLEAR(TYPE,data,N-2); } }
template <typename TYPE, int N> template <typename T> Point<TYPE,N>::Point( const Point3<T> &p ) { if ( N <= 3 ) { CY_MEMCONVERT(TYPE,data,&p.x,N); } else { CY_MEMCONVERT(TYPE,data,&p.x,3); CY_MEMCLEAR(TYPE,data,N-3); } }
template <typename TYPE, int N> template <typename T> Point<TYPE,N>::Point( const Point4<T> &p ) { if ( N <= 4 ) { CY_MEMCONVERT(TYPE,data,&p.x,N); } else { CY_MEMCONVERT(TYPE,data,&p.x,4); CY_MEMCLEAR(TYPE,data,N-4); } }
template <typename TYPE> Point2<TYPE>::Point2( const Point3<TYPE> &p ) : x(p.x), y(p.y)         {}
template <typename TYPE> Point2<TYPE>::Point2( const Point4<TYPE> &p ) : x(p.x), y(p.y)         {}
template <typename TYPE> Point3<TYPE>::Point3( const Point4<TYPE> &p ) : x(p.x), y(p.y), z(p.z) {}
template <typename TYPE> template <typename T> Point2<TYPE>::Point2( const Point3<T> &p ) : x(TYPE(p.x)), y(TYPE(p.y))               {}
template <typename TYPE> template <typename T> Point2<TYPE>::Point2( const Point4<T> &p ) : x(TYPE(p.x)), y(TYPE(p.y))               {}
template <typename TYPE> template <typename T> Point3<TYPE>::Point3( const Point4<T> &p ) : x(TYPE(p.x)), y(TYPE(p.y)), z(TYPE(p.z)) {}

//-------------------------------------------------------------------------------

typedef Point2<float>    Point2f;	//!< 2D point (vector) class with float type elements
typedef Point3<float>    Point3f;	//!< 3D point (vector) class with float type elements
typedef Point4<float>    Point4f;	//!< 4D point (vector) class with float type elements

typedef Point2<double>   Point2d;	//!< 2D point (vector) class with double type elements
typedef Point3<double>   Point3d;	//!< 3D point (vector) class with double type elements
typedef Point4<double>   Point4d;	//!< 4D point (vector) class with double type elements

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Point2f  cyPoint2f;		//!< 2D point (vector) class with float type elements
typedef cy::Point3f  cyPoint3f;		//!< 3D point (vector) class with float type elements
typedef cy::Point4f  cyPoint4f;		//!< 4D point (vector) class with float type elements

typedef cy::Point2d  cyPoint2d;		//!< 2D point (vector) class with double type elements
typedef cy::Point3d  cyPoint3d;		//!< 3D point (vector) class with double type elements
typedef cy::Point4d  cyPoint4d;		//!< 4D point (vector) class with double type elements

//-------------------------------------------------------------------------------

#endif

