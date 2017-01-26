// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyIPoint.h 
//! \author Cem Yuksel
//! 
//! \brief  2D, 3D, 4D, and ND integer point classes.
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

#ifndef _CY_IPOINT_H_INCLUDED_
#define _CY_IPOINT_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyCore.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

// Forward declarations
//!	\cond HIDDEN_SYMBOLS
template <typename TYPE> class IPoint2;
template <typename TYPE> class IPoint3;
template <typename TYPE> class IPoint4;
//! \endcond

//-------------------------------------------------------------------------------

//! A general class for N-dimensional integer points.

template <typename TYPE, int N>
class IPoint
{
	friend IPoint operator + ( const TYPE v, const IPoint &p ) { return   p+v; }	//!< Addition with a constant
	friend IPoint operator - ( const TYPE v, const IPoint &p ) { return -(p-v); }	//!< Subtraction from a constant
	friend IPoint operator * ( const TYPE v, const IPoint &p ) { return   p*v; }	//!< Multiplication with a constant

public:

	//!@name Components of the point
	TYPE data[N];

	//!@name Constructors
	IPoint() {}
	IPoint( const IPoint &p )        { CY_MEMCOPY(TYPE,data,p.data,N); }
	explicit IPoint( const TYPE *p ) { CY_MEMCOPY(TYPE,data,p,N); }
	explicit IPoint( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i]=v; }
	template <typename T> explicit IPoint( const IPoint<T,N> &p ) { CY_MEMCONVERT(TYPE,data,p.data,N); }
	template <int M> explicit IPoint( const IPoint<TYPE,M> &p )
	{
		if ( N <= M ) { CY_MEMCOPY(TYPE,data,p.data,N); }
		else          { CY_MEMCOPY(TYPE,data,p.data,M); CY_MEMCLEAR(TYPE,data,N-M); }
	}
	template <typename T, int M> explicit IPoint( const IPoint<T,M> &p )
	{
		if ( N <= M ) { CY_MEMCONVERT(TYPE,data,p.data,N); }
		else          { CY_MEMCONVERT(TYPE,data,p.data,M); CY_MEMCLEAR(TYPE,data,N-M); }
	}
	explicit IPoint( const IPoint2<TYPE> &p );
	explicit IPoint( const IPoint3<TYPE> &p );
	explicit IPoint( const IPoint4<TYPE> &p );
	template <typename T> explicit IPoint( const IPoint2<T> &p );
	template <typename T> explicit IPoint( const IPoint3<T> &p );
	template <typename T> explicit IPoint( const IPoint4<T> &p );
	template <typename P> explicit IPoint( const P &p ) { for ( int i=0; i<N; ++i ) data[i]=(TYPE)p[i]; }

	//!@name Set & Get value methods
	void Zero()               { CY_MEMCLEAR(TYPE,data,N); }					//!< Sets the coordinates as zero
	void Get( TYPE *p ) const { CY_MEMCOPY(TYPE,p,data,N); }				//!< Puts the coordinate values into the array
	void Set( const TYPE *p ) { CY_MEMCOPY(TYPE,data,p,N); }				//!< Sets the coordinates using the values in the given array
	void Set( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i] = v; }	//!< Sets all coordinates using the given value
	template <int M> void CopyData( TYPE *p ) { if ( M <= N ) { CY_MEMCOPY(TYPE,p,data,M); } else { CY_MEMCOPY(TYPE,p,data,N); CY_MEMCLEAR(TYPE,p+N,M-N); }	}
	template <typename T, int M> void ConvertData( T *p ) { if ( M <= N ) { CY_MEMCONVERT(T,p,data,M); } else { CY_MEMCONVERT(T,p,data,N); CY_MEMCLEAR(T,p+N,M-N); }	}

	//!@name General methods
	TYPE  Sum   () const { TYPE v=data[0]; for ( int i=1; i<N; ++i ) v+=data[i]; return v; }		//!< Returns the sum of its components
	bool  IsZero() const { for ( int i=0; i<N; ++i ) if ( data[i] != TYPE(0) ) return false; return true; }	//!< Returns true if all components are exactly zero
	TYPE  Min   () const { TYPE m = data[0]; for ( int i=1; i<N; ++i ) if ( m > data[i] ) m = data[i]; return m; }
	TYPE  Max   () const { TYPE m = data[0]; for ( int i=1; i<N; ++i ) if ( m < data[i] ) m = data[i]; return m; }
	int   MinID () const { TYPE m = data[0]; int ix=0; for ( int i=1; i<N; ++i ) if ( m > data[i] ) { m = data[i]; ix = i; } return m; }
	int   MaxID () const { TYPE m = data[0]; int ix=0; for ( int i=1; i<N; ++i ) if ( m < data[i] ) { m = data[i]; ix = i; } return m; }

	//!@name Limit methods
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i] = (data[i]<v) ? v : data[i]; }
	void ClampMax( const TYPE &v ) { for ( int i=0; i<N; ++i ) data[i] = (data[i]>v) ? v : data[i]; }
	void Abs() { for ( int i=0; i<N; i++ ) data[i] = cyAbs(data[i]); }	//!< Converts all negative components to positive values

	//!@name Unary operators
	IPoint operator - () const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i]=-data[i]; return r; } 

	//!@name Binary operators
	IPoint operator + ( const IPoint &p ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] + p.data[i]; return r; }
	IPoint operator - ( const IPoint &p ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] - p.data[i]; return r; }
	IPoint operator * ( const IPoint &p ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] * p.data[i]; return r; }
	IPoint operator / ( const IPoint &p ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] / p.data[i]; return r; }
	IPoint operator + ( const TYPE   &v ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] + v; return r; }
	IPoint operator - ( const TYPE   &v ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] - v; return r; }
	IPoint operator * ( const TYPE   &v ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] * v; return r; }
	IPoint operator / ( const TYPE   &v ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] / v; return r; }

	//!@name Assignment operators
	const IPoint& operator  = ( const IPoint &p ) { CY_MEMCOPY(TYPE,data,p.data,N); return *this; }	
	const IPoint& operator += ( const IPoint &p ) { for ( int i=0; i<N; ++i ) data[i] += p.data[i]; return *this; }
	const IPoint& operator -= ( const IPoint &p ) { for ( int i=0; i<N; ++i ) data[i] -= p.data[i]; return *this; }
	const IPoint& operator *= ( const IPoint &p ) { for ( int i=0; i<N; ++i ) data[i] *= p.data[i]; return *this; }
	const IPoint& operator /= ( const IPoint &p ) { for ( int i=0; i<N; ++i ) data[i] /= p.data[i]; return *this; }
	const IPoint& operator += ( const TYPE    v ) { for ( int i=0; i<N; ++i ) data[i] += v; return *this; }
	const IPoint& operator -= ( const TYPE    v ) { for ( int i=0; i<N; ++i ) data[i] -= v; return *this; }
	const IPoint& operator *= ( const TYPE    v ) { for ( int i=0; i<N; ++i ) data[i] *= v; return *this; }
	const IPoint& operator /= ( const TYPE    v ) { for ( int i=0; i<N; ++i ) data[i] /= v; return *this; }

	//!@name Bitwise operators
	IPoint operator << ( const IPoint &p ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] << p.data[i]; return r; }
	IPoint operator >> ( const IPoint &p ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] >> p.data[i]; return r; }
	IPoint operator  & ( const IPoint &p ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i]  & p.data[i]; return r; }
	IPoint operator  | ( const IPoint &p ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i]  | p.data[i]; return r; }
	IPoint operator  ^ ( const IPoint &p ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i]  ^ p.data[i]; return r; }
	IPoint operator << ( const TYPE   &v ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] << v; return r; }
	IPoint operator >> ( const TYPE   &v ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i] >> v; return r; }
	IPoint operator  & ( const TYPE   &v ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i]  & v; return r; }
	IPoint operator  | ( const TYPE   &v ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i]  | v; return r; }
	IPoint operator  ^ ( const TYPE   &v ) const { IPoint r; for ( int i=0; i<N; ++i ) r.data[i] = data[i]  ^ v; return r; }

	//!@name Bitwise Assignment operators
	const IPoint& operator <<= ( const IPoint &p ) { for ( int i=0; i<N; ++i ) data[i] <<= p.data[i]; return *this; }
	const IPoint& operator >>= ( const IPoint &p ) { for ( int i=0; i<N; ++i ) data[i] >>= p.data[i]; return *this; }
	const IPoint& operator  &= ( const IPoint &p ) { for ( int i=0; i<N; ++i ) data[i]  &= p.data[i]; return *this; }
	const IPoint& operator  |= ( const IPoint &p ) { for ( int i=0; i<N; ++i ) data[i]  |= p.data[i]; return *this; }
	const IPoint& operator  ^= ( const IPoint &p ) { for ( int i=0; i<N; ++i ) data[i]  ^= p.data[i]; return *this; }
	const IPoint& operator <<= ( const TYPE   &v ) { for ( int i=0; i<N; ++i ) data[i] <<= v; return *this; }
	const IPoint& operator >>= ( const TYPE   &v ) { for ( int i=0; i<N; ++i ) data[i] >>= v; return *this; }
	const IPoint& operator  &= ( const TYPE   &v ) { for ( int i=0; i<N; ++i ) data[i]  &= v; return *this; }
	const IPoint& operator  |= ( const TYPE   &v ) { for ( int i=0; i<N; ++i ) data[i]  |= v; return *this; }
	const IPoint& operator  ^= ( const TYPE   &v ) { for ( int i=0; i<N; ++i ) data[i]  ^= v; return *this; }

	//!@name Test operators
	bool operator == ( const IPoint& p ) const { for ( int i=0; i<N; ++i ) if ( data[i] != p.data[i] ) return false; return true; }
	bool operator != ( const IPoint& p ) const { for ( int i=0; i<N; ++i ) if ( data[i] != p.data[i] ) return true; return false; }

	//!@name Access operators
	TYPE& operator [] ( int i )       { return data[i]; }
	TYPE  operator [] ( int i ) const { return data[i]; }

	//!@name Dot product
	TYPE Dot        ( const IPoint &p ) const { IPoint r=operator*(p); return r.Sum(); }	//!< Dot product
	TYPE operator % ( const IPoint &p ) const { return Dot(p); }							//!< Dot product operator
};

//-------------------------------------------------------------------------------

//! 2D integer point class

template <typename TYPE>
class IPoint2
{
	friend IPoint2 operator + ( const TYPE v, const IPoint2 &p ) { return   p+v;  }	//!< Addition with a constant
	friend IPoint2 operator - ( const TYPE v, const IPoint2 &p ) { return -(p-v); }	//!< Subtraction from a constant
	friend IPoint2 operator * ( const TYPE v, const IPoint2 &p ) { return   p*v;  }	//!< Multiplication with a constant

public:

	//!@name Components of the point
	TYPE x, y;

	//!@name Constructors
	IPoint2() {}

	IPoint2( const TYPE &_x, const TYPE &_y ) : x( _x), y( _y) {}
	IPoint2( const IPoint2 &p )               : x(p.x), y(p.y) {}
	explicit IPoint2( const TYPE &v )         : x(v  ), y(v  ) {}
	explicit IPoint2( const IPoint3<TYPE> &p );
	explicit IPoint2( const IPoint4<TYPE> &p );
	template <typename T> explicit IPoint2( const IPoint2<T> &p ) : x(TYPE(p.x)), y(TYPE(p.y)) {}
	template <typename T> explicit IPoint2( const IPoint3<T> &p );
	template <typename T> explicit IPoint2( const IPoint4<T> &p );
	template <int M> explicit IPoint2( const IPoint<TYPE,M> &p ) { p.CopyData<2>(&x); }
	template <typename T, int M> explicit IPoint2( const IPoint<T,M> &p ) { p.template ConvertData<TYPE,2>(&x); }
	template <typename P> explicit IPoint2( const P &p ) : x((TYPE)p[0]), y((TYPE)p[1]) {}

	//!@name Set & Get value methods
	void Zero()               { CY_MEMCLEAR(TYPE,Data(),2); }		//!< Sets the coordinates as zero.
	void Get( TYPE *p ) const { ((IPoint2*)p)->operator=(*this); }	//!< Puts the coordinate values into the array.
	void Set( const TYPE *p ) { operator=(*((IPoint2*)p)); }		//!< Sets the coordinates using the values in the given array.
	void Set( const TYPE &v ) { x=v; y=v; }							//!< Sets all coordinates using the given value
	void Set( const TYPE &_x, const TYPE &_y ) { x=_x; y=_y; }		//!< Sets the coordinates using the given values

	//!@name General methods
	TYPE   Sum   () const { return x+y; }							//!< Returns the sum of its components
	bool   IsZero() const { return x==TYPE(0) && y==TYPE(0); }		//!< Returns true if all components are exactly zero
	TYPE   Min   () const { return x<y ? x : y; }
	TYPE   Max   () const { return x>y ? x : y; }
	int    MinID () const { return x<y ? 0 : 1; }
	int    MaxID () const { return x>y ? 0 : 1; }

	//!@name Limit methods
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { x=(x<v)?v:x; y=(y<v)?v:y; }
	void ClampMax( const TYPE &v ) { x=(x>v)?v:x; y=(y>v)?v:y; }
	void Abs() { x=cyAbs(x); y=cyAbs(y); }	//!< Converts all negative components to positive values

	//!@name Unary operators
	IPoint2 operator - () const { IPoint2 r; r.x=-x; r.y=-y; return r; } 

	//!@name Binary operators
	IPoint2 operator + ( const IPoint2 &p ) const { IPoint2 r; r.x=x+p.x; r.y=y+p.y; return r; }
	IPoint2 operator - ( const IPoint2 &p ) const { IPoint2 r; r.x=x-p.x; r.y=y-p.y; return r; }
	IPoint2 operator * ( const IPoint2 &p ) const { IPoint2 r; r.x=x*p.x; r.y=y*p.y; return r; }
	IPoint2 operator / ( const IPoint2 &p ) const { IPoint2 r; r.x=x/p.x; r.y=y/p.y; return r; }
	IPoint2 operator + ( const TYPE    &v ) const { IPoint2 r; r.x=x+v;   r.y=y+v;   return r; }
	IPoint2 operator - ( const TYPE    &v ) const { IPoint2 r; r.x=x-v;   r.y=y-v;   return r; }
	IPoint2 operator * ( const TYPE    &v ) const { IPoint2 r; r.x=x*v;   r.y=y*v;   return r; }
	IPoint2 operator / ( const TYPE    &v ) const { IPoint2 r; r.x=x/v;   r.y=y/v;   return r; }

	//!@name Assignment operators
	const IPoint2& operator  = ( const IPoint2 &p ) { x =p.x; y =p.y; return *this; }	
	const IPoint2& operator += ( const IPoint2 &p ) { x+=p.x; y+=p.y; return *this; }
	const IPoint2& operator -= ( const IPoint2 &p ) { x-=p.x; y-=p.y; return *this; }
	const IPoint2& operator *= ( const IPoint2 &p ) { x*=p.x; y*=p.y; return *this; }
	const IPoint2& operator /= ( const IPoint2 &p ) { x/=p.x; y/=p.y; return *this; }
	const IPoint2& operator += ( const TYPE    &v ) { x+=v;   y+=v;   return *this; }
	const IPoint2& operator -= ( const TYPE    &v ) { x-=v;   y-=v;   return *this; }
	const IPoint2& operator *= ( const TYPE    &v ) { x*=v;   y*=v;   return *this; }
	const IPoint2& operator /= ( const TYPE    &v ) { x/=v;   y/=v;   return *this; }

	//!@name Bitwise operators
	IPoint2 operator << ( const IPoint2 &p ) const { IPoint2 r; r.x = x << p.x; r.y = y << p.y; return r; }
	IPoint2 operator >> ( const IPoint2 &p ) const { IPoint2 r; r.x = x >> p.x; r.y = y >> p.y; return r; }
	IPoint2 operator  & ( const IPoint2 &p ) const { IPoint2 r; r.x = x  & p.x; r.y = y  & p.y; return r; }
	IPoint2 operator  | ( const IPoint2 &p ) const { IPoint2 r; r.x = x  | p.x; r.y = y  | p.y; return r; }
	IPoint2 operator  ^ ( const IPoint2 &p ) const { IPoint2 r; r.x = x  ^ p.x; r.y = y  ^ p.y; return r; }
	IPoint2 operator << ( const TYPE    &v ) const { IPoint2 r; r.x = x << v;   r.y = y << v;   return r; }
	IPoint2 operator >> ( const TYPE    &v ) const { IPoint2 r; r.x = x >> v;   r.y = y >> v;   return r; }
	IPoint2 operator  & ( const TYPE    &v ) const { IPoint2 r; r.x = x  & v;   r.y = y  & v;   return r; }
	IPoint2 operator  | ( const TYPE    &v ) const { IPoint2 r; r.x = x  | v;   r.y = y  | v;   return r; }
	IPoint2 operator  ^ ( const TYPE    &v ) const { IPoint2 r; r.x = x  ^ v;   r.y = y  ^ v;   return r; }

	//!@name Bitwise Assignment operators
	const IPoint2& operator <<= ( const IPoint2 &p ) { x<<=p.x; y<<=p.y; return *this; }
	const IPoint2& operator >>= ( const IPoint2 &p ) { x>>=p.x; y>>=p.y; return *this; }
	const IPoint2& operator  &= ( const IPoint2 &p ) { x &=p.x; y &=p.y; return *this; }
	const IPoint2& operator  |= ( const IPoint2 &p ) { x |=p.x; y |=p.y; return *this; }
	const IPoint2& operator  ^= ( const IPoint2 &p ) { x ^=p.x; y ^=p.y; return *this; }
	const IPoint2& operator <<= ( const TYPE    &v ) { x<<=v;   y<<=v;   return *this; }
	const IPoint2& operator >>= ( const TYPE    &v ) { x>>=v;   y>>=v;   return *this; }
	const IPoint2& operator  &= ( const TYPE    &v ) { x &=v;   y &=v;   return *this; }
	const IPoint2& operator  |= ( const TYPE    &v ) { x |=v;   y |=v;   return *this; }
	const IPoint2& operator  ^= ( const TYPE    &v ) { x ^=v;   y ^=v;   return *this; }

	//!@name Test operators
	bool operator == ( const IPoint2& p ) const { return x==p.x && y==p.y; }
	bool operator != ( const IPoint2& p ) const { return x!=p.x && y!=p.y; }

	//!@name Access operators
	TYPE&       operator [] ( int i )       { return Element(i); }
	const TYPE& operator [] ( int i ) const { return Element(i); }
	TYPE&       Element     ( int i )       { return (&x)[i]; }
	const TYPE& Element     ( int i ) const { return (&x)[i]; }
	TYPE*       Data        ()              { return &x; }
	const TYPE* Data        ()        const { return &x; }

	//!@name Cross product and dot product
	TYPE Dot        ( const IPoint2 &p ) const { IPoint2 r=operator*(p); return r.Sum(); }	//!< Dot product
	TYPE operator % ( const IPoint2 &p ) const { return Dot(p); }							//!< Dot product operator
};

//-------------------------------------------------------------------------------

//! 3D integer point class

template <typename TYPE>
class IPoint3
{
	friend IPoint3 operator + ( const TYPE v, const IPoint3 &p ) { return   p+v;  }	//!< Addition with a constant
	friend IPoint3 operator - ( const TYPE v, const IPoint3 &p ) { return -(p-v); }	//!< Subtraction from a constant
	friend IPoint3 operator * ( const TYPE v, const IPoint3 &p ) { return   p*v;  }	//!< Multiplication with a constant

public:

	//!@name Components of the point
	TYPE x, y, z;

	//!@name Constructors
	IPoint3() { }
	IPoint3( const TYPE &_x, const TYPE &_y, const TYPE &_z ) : x( _x), y( _y), z( _z) {}
	IPoint3( const IPoint3 &p )                               : x(p.x), y(p.y), z(p.z) {}
	explicit IPoint3( const TYPE &v )                         : x(v  ), y(v  ), z(v  ) {}
	explicit IPoint3( const IPoint2<TYPE> &p, TYPE _z=0 )     : x(p.x), y(p.y), z( _z) {}
	explicit IPoint3( const IPoint4<TYPE> &p );
	template <typename T> explicit IPoint3( const IPoint3<T> &p )            : x(TYPE(p.x)), y(TYPE(p.y)), z(TYPE(p.z)) {}
	template <typename T> explicit IPoint3( const IPoint2<T> &p, TYPE _z=0 ) : x(TYPE(p.x)), y(TYPE(p.y)), z(      _z) {}
	template <typename T> explicit IPoint3( const IPoint4<T> &p );
	template <int M> explicit IPoint3( const IPoint<TYPE,M> &p ) { p.CopyData<3>(&x); }
	template <typename T, int M> explicit IPoint3( const IPoint<T,M> &p ) { p.template ConvertData<TYPE,3>(&x); }
	template <typename P> explicit IPoint3( const P &p ) : x((TYPE)p[0]), y((TYPE)p[1]), z((TYPE)p[2]) {}

	//!@name Set & Get value methods
	void Zero()               { CY_MEMCLEAR(TYPE,Data(),3); }		//!< Sets the coordinates as zero
	void Get( TYPE *p ) const { ((IPoint3*)p)->operator=(*this); }	//!< Puts the coordinate values into the array
	void Set( const TYPE *p ) { operator=(*((IPoint3*)p)); }		//!< Sets the coordinates using the values in the given array
	void Set( const TYPE &v ) { x=v; y=v; z=v; }					//!< Sets all coordinates using the given value
	void Set( const TYPE &_x, const TYPE &_y, const TYPE &_z ) { x=_x; y=_y; z=_z; }	//!< Sets the coordinates using the given values

	//!@name Length and Normalize methods
	TYPE   Sum   () const { return x+y+z; }										//!< Returns the sum of its components
	bool   IsZero() const { return x==TYPE(0) && y==TYPE(0) && z==TYPE(0); }	//!< Returns true if all components are exactly zero
	TYPE   Min   () const { return x<y ? (x<z ? x : z) : (y<z ? y : z); }
	TYPE   Max   () const { return x>y ? (x>z ? x : z) : (y>z ? y : z); }
	int    MinID () const { return x<y ? (x<z ? 0 : 2) : (y<z ? 1 : 2); }
	int    MaxID () const { return x>y ? (x>z ? 0 : 2) : (y>z ? 1 : 2); }

	//!@name Limit methods
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { x=(x<v)?v:x; y=(y<v)?v:y; z=(z<v)?v:z; }
	void ClampMax( const TYPE &v ) { x=(x>v)?v:x; y=(y>v)?v:y; z=(z>v)?v:z; }
	void Abs() { x=cyAbs(x); y=cyAbs(y); z=cyAbs(z); }	//!< Converts all negative components to positive values

	//!@name Unary operators
	IPoint3 operator - () const { IPoint3 r; r.x=-x; r.y=-y; r.z=-z; return r; } 

	//!@name Binary operators
	IPoint3 operator + ( const IPoint3 &p ) const { IPoint3 r; r.x=x+p.x; r.y=y+p.y; r.z=z+p.z; return r; }
	IPoint3 operator - ( const IPoint3 &p ) const { IPoint3 r; r.x=x-p.x; r.y=y-p.y; r.z=z-p.z; return r; }
	IPoint3 operator * ( const IPoint3 &p ) const { IPoint3 r; r.x=x*p.x; r.y=y*p.y; r.z=z*p.z; return r; }
	IPoint3 operator / ( const IPoint3 &p ) const { IPoint3 r; r.x=x/p.x; r.y=y/p.y; r.z=z/p.z; return r; }
	IPoint3 operator + ( const TYPE    &v ) const { IPoint3 r; r.x=x+v;   r.y=y+v;   r.z=z+v;   return r; }
	IPoint3 operator - ( const TYPE    &v ) const { IPoint3 r; r.x=x-v;   r.y=y-v;   r.z=z-v;   return r; }
	IPoint3 operator * ( const TYPE    &v ) const { IPoint3 r; r.x=x*v;   r.y=y*v;   r.z=z*v;   return r; }
	IPoint3 operator / ( const TYPE    &v ) const { IPoint3 r; r.x=x/v;   r.y=y/v;   r.z=z/v;   return r; }

	//!@name Assignment operators
	const IPoint3& operator  = ( const IPoint3 &p ) { x =p.x; y =p.y; z =p.z; return *this; }	
	const IPoint3& operator += ( const IPoint3 &p ) { x+=p.x; y+=p.y; z+=p.z; return *this; }
	const IPoint3& operator -= ( const IPoint3 &p ) { x-=p.x; y-=p.y; z-=p.z; return *this; }
	const IPoint3& operator *= ( const IPoint3 &p ) { x*=p.x; y*=p.y; z*=p.z; return *this; }
	const IPoint3& operator /= ( const IPoint3 &p ) { x/=p.x; y/=p.y; z/=p.z; return *this; }
	const IPoint3& operator += ( const TYPE    &v ) { x+=v;   y+=v;   z+=v;   return *this; }
	const IPoint3& operator -= ( const TYPE    &v ) { x-=v;   y-=v;   z-=v;   return *this; }
	const IPoint3& operator *= ( const TYPE    &v ) { x*=v;   y*=v;   z*=v;   return *this; }
	const IPoint3& operator /= ( const TYPE    &v ) { x/=v;   y/=v;   z/=v;   return *this; }

	//!@name Bitwise operators
	IPoint3 operator << ( const IPoint3 &p ) const { IPoint3 r; r.x = x << p.x; r.y = y << p.y; r.z = z << p.z; return r; }
	IPoint3 operator >> ( const IPoint3 &p ) const { IPoint3 r; r.x = x >> p.x; r.y = y >> p.y; r.z = z >> p.z; return r; }
	IPoint3 operator  & ( const IPoint3 &p ) const { IPoint3 r; r.x = x  & p.x; r.y = y  & p.y; r.z = z  & p.z; return r; }
	IPoint3 operator  | ( const IPoint3 &p ) const { IPoint3 r; r.x = x  | p.x; r.y = y  | p.y; r.z = z  | p.z; return r; }
	IPoint3 operator  ^ ( const IPoint3 &p ) const { IPoint3 r; r.x = x  ^ p.x; r.y = y  ^ p.y; r.z = z  ^ p.z; return r; }
	IPoint3 operator << ( const TYPE    &v ) const { IPoint3 r; r.x = x << v;   r.y = y << v;   r.z = z << v;   return r; }
	IPoint3 operator >> ( const TYPE    &v ) const { IPoint3 r; r.x = x >> v;   r.y = y >> v;   r.z = z >> v;   return r; }
	IPoint3 operator  & ( const TYPE    &v ) const { IPoint3 r; r.x = x  & v;   r.y = y  & v;   r.z = z  & v;   return r; }
	IPoint3 operator  | ( const TYPE    &v ) const { IPoint3 r; r.x = x  | v;   r.y = y  | v;   r.z = z  | v;   return r; }
	IPoint3 operator  ^ ( const TYPE    &v ) const { IPoint3 r; r.x = x  ^ v;   r.y = y  ^ v;   r.z = z  ^ v;   return r; }

	//!@name Bitwise Assignment operators
	const IPoint3& operator <<= ( const IPoint3 &p ) { x<<=p.x; y<<=p.y; z<<=p.z; return *this; }
	const IPoint3& operator >>= ( const IPoint3 &p ) { x>>=p.x; y>>=p.y; z>>=p.z; return *this; }
	const IPoint3& operator  &= ( const IPoint3 &p ) { x &=p.x; y &=p.y; z &=p.z; return *this; }
	const IPoint3& operator  |= ( const IPoint3 &p ) { x |=p.x; y |=p.y; z |=p.z; return *this; }
	const IPoint3& operator  ^= ( const IPoint3 &p ) { x ^=p.x; y ^=p.y; z ^=p.z; return *this; }
	const IPoint3& operator <<= ( const TYPE    &v ) { x<<=v;   y<<=v;   z<<=v;   return *this; }
	const IPoint3& operator >>= ( const TYPE    &v ) { x>>=v;   y>>=v;   z>>=v;   return *this; }
	const IPoint3& operator  &= ( const TYPE    &v ) { x &=v;   y &=v;   z &=v;   return *this; }
	const IPoint3& operator  |= ( const TYPE    &v ) { x |=v;   y |=v;   z |=v;   return *this; }
	const IPoint3& operator  ^= ( const TYPE    &v ) { x ^=v;   y ^=v;   z ^=v;   return *this; }

	//!@name Test operators
	bool operator == ( const IPoint3& p ) const { return x==p.x && y==p.y && z==p.z; }
	bool operator != ( const IPoint3& p ) const { return x!=p.x && y!=p.y && z!=p.z; }

	//!@name Access operators
	TYPE&       operator [] ( int i )       { return Element(i); }
	const TYPE& operator [] ( int i ) const { return Element(i); }
	TYPE&       Element     ( int i )       { return (&x)[i]; }
	const TYPE& Element     ( int i ) const { return (&x)[i]; }
	TYPE*       Data        ()              { return &x; }
	const TYPE* Data        ()        const { return &x; }

	//!@name Cross product and dot product
	TYPE   Dot        ( const IPoint3 &p ) const { IPoint3 r=operator*(p); return r.Sum(); }	//!< Dot product
	TYPE   operator % ( const IPoint3 &p ) const { return Dot(p); }								//!< Dot product

	//!@name Conversion Methods
	IPoint2<TYPE> XY() const { return IPoint2<TYPE>(x,y); }
};

//-------------------------------------------------------------------------------

//! 4D integer point class

template <typename TYPE>
class IPoint4
{
	friend IPoint4 operator + ( const TYPE v, const IPoint4 &p ) { return   p+v;  }	//!< Addition with a constant
	friend IPoint4 operator - ( const TYPE v, const IPoint4 &p ) { return -(p-v); }	//!< Subtraction from a constant
	friend IPoint4 operator * ( const TYPE v, const IPoint4 &p ) { return   p*v;  }	//!< Multiplication with a constant

public:

	//!@name Components of the point
	TYPE x, y, z, w;

	//!@name Constructors
	IPoint4() { }
	IPoint4( const TYPE &_x, const TYPE &_y, const TYPE &_z, const TYPE &_w ) : x( _x), y( _y), z( _z), w( _w) {}
	IPoint4( const IPoint4 &p )                                               : x(p.x), y(p.y), z(p.z), w(p.w) {}
	explicit IPoint4( const TYPE &v )                                         : x(v  ), y(v  ), z(v  ), w(v  ) {}
	explicit IPoint4( const IPoint3<TYPE> &p,            TYPE _w=0 )          : x(p.x), y(p.y), z(p.z), w( _w) {}
	explicit IPoint4( const IPoint2<TYPE> &p, TYPE _z=0, TYPE _w=0 )          : x(p.x), y(p.y), z( _z), w( _w) {}
	template <typename T> explicit IPoint4( const IPoint4<T> &p )                        : x(TYPE(p.x)), y(TYPE(p.y)), z(TYPE(p.z)), z(TYPE(p.w)) {}
	template <typename T> explicit IPoint4( const IPoint3<T> &p,            TYPE _w=0 )  : x(TYPE(p.x)), y(TYPE(p.y)), z(TYPE(p.z)), z(      _w ) {}
	template <typename T> explicit IPoint4( const IPoint2<T> &p, TYPE _z=0, TYPE _w=0 )  : x(TYPE(p.x)), y(TYPE(p.y)), z(      _z ), z(      _w ) {}
	template <int M> explicit IPoint4( const IPoint<TYPE,M> &p ) { p.CopyData<3>(&x); }
	template <typename T, int M> explicit IPoint4( const IPoint<T,M> &p ) { p.template ConvertData<TYPE,3>(&x); }
	template <typename P> explicit IPoint4( const P &p ) : x((TYPE)p[0]), y((TYPE)p[1]), z((TYPE)p[2]) {}

	//!@name Set & Get value methods
	void Zero()               { CY_MEMCLEAR(TYPE,Data(),4); }		//!< Sets the coordinates as zero
	void Get( TYPE *p ) const { ((IPoint4*)p)->operator=(*this); }	//!< Puts the coordinate values into the array
	void Set( const TYPE *p ) { operator=(*((IPoint4*)p)); }		//!< Sets the coordinates using the values in the given array
	void Set( const TYPE &v ) { x=v; y=v; z=v; w=v; }				//!< Sets all coordinates using the given value
	void Set( const TYPE &_x, const TYPE &_y, const TYPE &_z, const TYPE &_w=0 ) { x=_x; y=_y; z=_z; w=_w; }	//!< Sets the coordinates using the given values

	//!@name Length and Normalize methods
	TYPE   Sum   () const { return x+y+z+w; }						//!< Returns the sum of its components
	bool   IsZero() const { return x==TYPE(0) && y==TYPE(0) && z==TYPE(0); && w==TYPE(0); }	//!< Returns true if all components are exactly zero
	TYPE   Min   () const { TYPE mxy = x<y ? x : y; TYPE mzw = z<w ? z : w; return mxy<mzw ? mxy : mzw; }
	TYPE   Max   () const { TYPE mxy = x>y ? x : y; TYPE mzw = z>w ? z : w; return mxy>mzw ? mxy : mzw; }
	int    MinID () const { int  ixy = x<y ? 0 : 1; int  izw = z<w ? 2 : 3; return (&x)[ixy]<(&x)[izw] ? ixy : izw; }
	int    MaxID () const { int  ixy = x>y ? 0 : 1; int  izw = z>w ? 2 : 3; return (&x)[ixy]>(&x)[izw] ? ixy : izw; }

	//!@name Limit methods
	void Clamp( const TYPE &minValue, const TYPE &maxValue ) { ClampMin(minValue); ClampMax(maxValue); }
	void ClampMin( const TYPE &v ) { x=(x<v)?v:x; y=(y<v)?v:y; z=(z<v)?v:z; w=(w<v)?v:w; }
	void ClampMax( const TYPE &v ) { x=(x>v)?v:x; y=(y>v)?v:y; z=(z>v)?v:z; w=(w>v)?v:w; }
	void Abs() { x=cyAbs(x); y=cyAbs(y); z=cyAbs(z); w=cyAbs(w); }	//!< Converts all negative components to positive values

	//!@name Unary operators
	IPoint4 operator - () const { IPoint4 r; r.x=-x; r.y=-y; r.z=-z; r.w=-w; return r; } 

	//!@name Binary operators
	IPoint4 operator + ( const IPoint4 &p ) const { IPoint4 r; r.x=x+p.x; r.y=y+p.y; r.z=z+p.z; r.w=w+p.w; return r; }
	IPoint4 operator - ( const IPoint4 &p ) const { IPoint4 r; r.x=x-p.x; r.y=y-p.y; r.z=z-p.z; r.w=w-p.w; return r; }
	IPoint4 operator * ( const IPoint4 &p ) const { IPoint4 r; r.x=x*p.x; r.y=y*p.y; r.z=z*p.z; r.w=w*p.w; return r; }
	IPoint4 operator / ( const IPoint4 &p ) const { IPoint4 r; r.x=x/p.x; r.y=y/p.y; r.z=z/p.z; r.w=w/p.w; return r; }
	IPoint4 operator + ( const TYPE    &v ) const { IPoint4 r; r.x=x+v;   r.y=y+v;   r.z=z+v;   r.w=w+v;   return r; }
	IPoint4 operator - ( const TYPE    &v ) const { IPoint4 r; r.x=x-v;   r.y=y-v;   r.z=z-v;   r.w=w-v;   return r; }
	IPoint4 operator * ( const TYPE    &v ) const { IPoint4 r; r.x=x*v;   r.y=y*v;   r.z=z*v;   r.w=w*v;   return r; }
	IPoint4 operator / ( const TYPE    &v ) const { IPoint4 r; r.x=x/v;   r.y=y/v;   r.z=z/v;   r.w=w/v;   return r; }

	//!@name Assignment operators
	const IPoint4& operator  = ( const IPoint4 &p ) { x =p.x; y =p.y; z =p.z; w =p.w; return *this; }	
	const IPoint4& operator += ( const IPoint4 &p ) { x+=p.x; y+=p.y; z+=p.z; w+=p.w; return *this; }
	const IPoint4& operator -= ( const IPoint4 &p ) { x-=p.x; y-=p.y; z-=p.z; w-=p.w; return *this; }
	const IPoint4& operator *= ( const IPoint4 &p ) { x*=p.x; y*=p.y; z*=p.z; w*=p.w; return *this; }
	const IPoint4& operator /= ( const IPoint4 &p ) { x/=p.x; y/=p.y; z/=p.z; w/=p.w; return *this; }
	const IPoint4& operator += ( const TYPE    &v ) { x+=v;   y+=v;   z+=v;   w+=v;   return *this; }
	const IPoint4& operator -= ( const TYPE    &v ) { x-=v;   y-=v;   z-=v;   w-=v;   return *this; }
	const IPoint4& operator *= ( const TYPE    &v ) { x*=v;   y*=v;   z*=v;   w*=v;   return *this; }
	const IPoint4& operator /= ( const TYPE    &v ) { x/=v;   y/=v;   z/=v;   w/=v;   return *this; }

	//!@name Bitwise operators
	IPoint4 operator << ( const IPoint4 &p ) const { IPoint4 r; r.x = x << p.x; r.y = y << p.y; r.z = z << p.z; r.w = w << p.w; return r; }
	IPoint4 operator >> ( const IPoint4 &p ) const { IPoint4 r; r.x = x >> p.x; r.y = y >> p.y; r.z = z >> p.z; r.w = w >> p.w; return r; }
	IPoint4 operator  & ( const IPoint4 &p ) const { IPoint4 r; r.x = x  & p.x; r.y = y  & p.y; r.z = z  & p.z; r.w = w  & p.w; return r; }
	IPoint4 operator  | ( const IPoint4 &p ) const { IPoint4 r; r.x = x  | p.x; r.y = y  | p.y; r.z = z  | p.z; r.w = w  | p.w; return r; }
	IPoint4 operator  ^ ( const IPoint4 &p ) const { IPoint4 r; r.x = x  ^ p.x; r.y = y  ^ p.y; r.z = z  ^ p.z; r.w = w  ^ p.w; return r; }
	IPoint4 operator << ( const TYPE    &v ) const { IPoint4 r; r.x = x << v;   r.y = y << v;   r.z = z << v;   r.w = w << v;   return r; }
	IPoint4 operator >> ( const TYPE    &v ) const { IPoint4 r; r.x = x >> v;   r.y = y >> v;   r.z = z >> v;   r.w = w >> v;   return r; }
	IPoint4 operator  & ( const TYPE    &v ) const { IPoint4 r; r.x = x  & v;   r.y = y  & v;   r.z = z  & v;   r.w = w  & v;   return r; }
	IPoint4 operator  | ( const TYPE    &v ) const { IPoint4 r; r.x = x  | v;   r.y = y  | v;   r.z = z  | v;   r.w = w  | v;   return r; }
	IPoint4 operator  ^ ( const TYPE    &v ) const { IPoint4 r; r.x = x  ^ v;   r.y = y  ^ v;   r.z = z  ^ v;   r.w = w  ^ v;   return r; }

	//!@name Bitwise Assignment operators
	const IPoint4& operator <<= ( const IPoint4 &p ) { x<<=p.x; y<<=p.y; z<<=p.z; w<<=p.w; return *this; }
	const IPoint4& operator >>= ( const IPoint4 &p ) { x>>=p.x; y>>=p.y; z>>=p.z; w>>=p.w; return *this; }
	const IPoint4& operator  &= ( const IPoint4 &p ) { x &=p.x; y &=p.y; z &=p.z; w &=p.w; return *this; }
	const IPoint4& operator  |= ( const IPoint4 &p ) { x |=p.x; y |=p.y; z |=p.z; w |=p.w; return *this; }
	const IPoint4& operator  ^= ( const IPoint4 &p ) { x ^=p.x; y ^=p.y; z ^=p.z; w ^=p.w; return *this; }
	const IPoint4& operator <<= ( const TYPE    &v ) { x<<=v;   y<<=v;   z<<=v;   w<<=v;   return *this; }
	const IPoint4& operator >>= ( const TYPE    &v ) { x>>=v;   y>>=v;   z>>=v;   w>>=v;   return *this; }
	const IPoint4& operator  &= ( const TYPE    &v ) { x &=v;   y &=v;   z &=v;   w &=v;   return *this; }
	const IPoint4& operator  |= ( const TYPE    &v ) { x |=v;   y |=v;   z |=v;   w |=v;   return *this; }
	const IPoint4& operator  ^= ( const TYPE    &v ) { x ^=v;   y ^=v;   z ^=v;   w ^=v;   return *this; }

	//!@name Test operators
	bool operator == ( const IPoint4& p ) const { return x==p.x && y==p.y && z==p.z; && w==p.w; }
	bool operator != ( const IPoint4& p ) const { return x!=p.x && y!=p.y && z!=p.z; && w!=p.w; }

	//!@name Access operators
	TYPE&       operator [] ( int i )       { return Element(i); }
	const TYPE& operator [] ( int i ) const { return Element(i); }
	TYPE&       Element     ( int i )       { return (&x)[i]; }
	const TYPE& Element     ( int i ) const { return (&x)[i]; }
	TYPE*       Data        ()              { return &x; }
	const TYPE* Data        ()        const { return &x; }

	//!@name Cross product and dot product
	TYPE   Dot        ( const IPoint4 &p ) const { IPoint4 r=operator*(p); return r.Sum(); }	//!< Dot product
	TYPE   operator % ( const IPoint4 &p ) const { return Dot(p); }								//!< Dot product

	//!@name Conversion Methods
	IPoint2<TYPE> XY()  const { return IPoint2<TYPE>(x,y); }
	IPoint3<TYPE> XYZ() const { return IPoint3<TYPE>(*this); }
};

//-------------------------------------------------------------------------------

// Definitions of the conversion constructors
template <typename TYPE, int N> IPoint<TYPE,N>::IPoint( const IPoint2<TYPE> &p ) { if ( N <= 2 ) { CY_MEMCOPY(TYPE,data,&p.x,N); } else { CY_MEMCOPY(TYPE,data,&p.x,2); CY_MEMCLEAR(TYPE,data,N-2); } }
template <typename TYPE, int N> IPoint<TYPE,N>::IPoint( const IPoint3<TYPE> &p ) { if ( N <= 3 ) { CY_MEMCOPY(TYPE,data,&p.x,N); } else { CY_MEMCOPY(TYPE,data,&p.x,3); CY_MEMCLEAR(TYPE,data,N-3); } }
template <typename TYPE, int N> IPoint<TYPE,N>::IPoint( const IPoint4<TYPE> &p ) { if ( N <= 4 ) { CY_MEMCOPY(TYPE,data,&p.x,N); } else { CY_MEMCOPY(TYPE,data,&p.x,4); CY_MEMCLEAR(TYPE,data,N-4); } }
template <typename TYPE, int N> template <typename T> IPoint<TYPE,N>::IPoint( const IPoint2<T> &p ) { if ( N <= 2 ) { CY_MEMCONVERT(TYPE,data,&p.x,N); } else { CY_MEMCONVERT(TYPE,data,&p.x,2); CY_MEMCLEAR(TYPE,data,N-2); } }
template <typename TYPE, int N> template <typename T> IPoint<TYPE,N>::IPoint( const IPoint3<T> &p ) { if ( N <= 3 ) { CY_MEMCONVERT(TYPE,data,&p.x,N); } else { CY_MEMCONVERT(TYPE,data,&p.x,3); CY_MEMCLEAR(TYPE,data,N-3); } }
template <typename TYPE, int N> template <typename T> IPoint<TYPE,N>::IPoint( const IPoint4<T> &p ) { if ( N <= 4 ) { CY_MEMCONVERT(TYPE,data,&p.x,N); } else { CY_MEMCONVERT(TYPE,data,&p.x,4); CY_MEMCLEAR(TYPE,data,N-4); } }
template <typename TYPE> IPoint2<TYPE>::IPoint2( const IPoint3<TYPE> &p ) : x(p.x), y(p.y)         {}
template <typename TYPE> IPoint2<TYPE>::IPoint2( const IPoint4<TYPE> &p ) : x(p.x), y(p.y)         {}
template <typename TYPE> IPoint3<TYPE>::IPoint3( const IPoint4<TYPE> &p ) : x(p.x), y(p.y), z(p.z) {}
template <typename TYPE> template <typename T> IPoint2<TYPE>::IPoint2( const IPoint3<T> &p ) : x(TYPE(p.x)), y(TYPE(p.y))               {}
template <typename TYPE> template <typename T> IPoint2<TYPE>::IPoint2( const IPoint4<T> &p ) : x(TYPE(p.x)), y(TYPE(p.y))               {}
template <typename TYPE> template <typename T> IPoint3<TYPE>::IPoint3( const IPoint4<T> &p ) : x(TYPE(p.x)), y(TYPE(p.y)), z(TYPE(p.z)) {}

//-------------------------------------------------------------------------------

typedef IPoint2<int8_t>   IPoint2b;		//!< 8-bit  signed   integer (int8_t)   2D integer point class
typedef IPoint3<int8_t>   IPoint3b;		//!< 8-bit  signed   integer (int8_t)   3D integer point class
typedef IPoint4<int8_t>   IPoint4b;		//!< 8-bit  signed   integer (int8_t)   4D integer point class
typedef IPoint2<uint8_t>  IPoint2ub;	//!< 8-bit  unsigned integer (uint8_t)  2D integer point class
typedef IPoint3<uint8_t>  IPoint3ub;	//!< 8-bit  unsigned integer (uint8_t)  3D integer point class
typedef IPoint4<uint8_t>  IPoint4ub;	//!< 8-bit  unsigned integer (uint8_t)  4D integer point class
typedef IPoint2<int16_t>  IPoint2s;		//!< 16-bit signed   integer (int16_t)  2D integer point class
typedef IPoint3<int16_t>  IPoint3s;		//!< 16-bit signed   integer (int16_t)  3D integer point class
typedef IPoint4<int16_t>  IPoint4s;		//!< 16-bit signed   integer (int16_t)  4D integer point class
typedef IPoint2<uint16_t> IPoint2us;	//!< 16-bit unsigned integer (uint16_t) 2D integer point class
typedef IPoint3<uint16_t> IPoint3us;	//!< 16-bit unsigned integer (uint16_t) 3D integer point class
typedef IPoint4<uint16_t> IPoint4us;	//!< 16-bit unsigned integer (uint16_t) 4D integer point class
typedef IPoint2<int32_t>  IPoint2i;		//!< 32-bit signed   integer (int32_t)  2D integer point class
typedef IPoint3<int32_t>  IPoint3i;		//!< 32-bit signed   integer (int32_t)  3D integer point class
typedef IPoint4<int32_t>  IPoint4i;		//!< 32-bit signed   integer (int32_t)  4D integer point class
typedef IPoint2<uint32_t> IPoint2ui;	//!< 32-bit unsigned integer (uint32_t) 2D integer point class
typedef IPoint3<uint32_t> IPoint3ui;	//!< 32-bit unsigned integer (uint32_t) 3D integer point class
typedef IPoint4<uint32_t> IPoint4ui;	//!< 32-bit unsigned integer (uint32_t) 4D integer point class
typedef IPoint2<int64_t>  IPoint2l;		//!< 64-bit signed   integer (int64_t)  2D integer point class
typedef IPoint3<int64_t>  IPoint3l;		//!< 64-bit signed   integer (int64_t)  3D integer point class
typedef IPoint4<int64_t>  IPoint4l;		//!< 64-bit signed   integer (int64_t)  4D integer point class
typedef IPoint2<uint64_t> IPoint2ul;	//!< 64-bit unsigned integer (uint64_t) 2D integer point class
typedef IPoint3<uint64_t> IPoint3ul;	//!< 64-bit unsigned integer (uint64_t) 3D integer point class
typedef IPoint4<uint64_t> IPoint4ul;	//!< 64-bit unsigned integer (uint64_t) 4D integer point class

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::IPoint2b  cyIPoint2b;		//!< 8-bit  signed   integer (int8_t)   2D integer point class
typedef cy::IPoint3b  cyIPoint3b;		//!< 8-bit  signed   integer (int8_t)   3D integer point class
typedef cy::IPoint4b  cyIPoint4b;		//!< 8-bit  signed   integer (int8_t)   4D integer point class
typedef cy::IPoint2ub cyIPoint2ub;		//!< 8-bit  unsigned integer (uint8_t)  2D integer point class
typedef cy::IPoint3ub cyIPoint3ub;		//!< 8-bit  unsigned integer (uint8_t)  3D integer point class
typedef cy::IPoint4ub cyIPoint4ub;		//!< 8-bit  unsigned integer (uint8_t)  4D integer point class
typedef cy::IPoint2s  cyIPoint2s;		//!< 16-bit signed   integer (int16_t)  2D integer point class
typedef cy::IPoint3s  cyIPoint3s;		//!< 16-bit signed   integer (int16_t)  3D integer point class
typedef cy::IPoint4s  cyIPoint4s;		//!< 16-bit signed   integer (int16_t)  4D integer point class
typedef cy::IPoint2us cyIPoint2us;		//!< 16-bit unsigned integer (uint16_t) 2D integer point class
typedef cy::IPoint3us cyIPoint3us;		//!< 16-bit unsigned integer (uint16_t) 3D integer point class
typedef cy::IPoint4us cyIPoint4us;		//!< 16-bit unsigned integer (uint16_t) 4D integer point class
typedef cy::IPoint2i  cyIPoint2i;		//!< 32-bit signed   integer (int32_t)  2D integer point class
typedef cy::IPoint3i  cyIPoint3i;		//!< 32-bit signed   integer (int32_t)  3D integer point class
typedef cy::IPoint4i  cyIPoint4i;		//!< 32-bit signed   integer (int32_t)  4D integer point class
typedef cy::IPoint2ui cyIPoint2ui;		//!< 32-bit unsigned integer (uint32_t) 2D integer point class
typedef cy::IPoint3ui cyIPoint3ui;		//!< 32-bit unsigned integer (uint32_t) 3D integer point class
typedef cy::IPoint4ui cyIPoint4ui;		//!< 32-bit unsigned integer (uint32_t) 4D integer point class
typedef cy::IPoint2l  cyIPoint2l;		//!< 64-bit signed   integer (int64_t)  2D integer point class
typedef cy::IPoint3l  cyIPoint3l;		//!< 64-bit signed   integer (int64_t)  3D integer point class
typedef cy::IPoint4l  cyIPoint4l;		//!< 64-bit signed   integer (int64_t)  4D integer point class
typedef cy::IPoint2ul cyIPoint2ul;		//!< 64-bit unsigned integer (uint64_t) 2D integer point class
typedef cy::IPoint3ul cyIPoint3ul;		//!< 64-bit unsigned integer (uint64_t) 3D integer point class
typedef cy::IPoint4ul cyIPoint4ul;		//!< 64-bit unsigned integer (uint64_t) 4D integer point class

//-------------------------------------------------------------------------------

#endif

