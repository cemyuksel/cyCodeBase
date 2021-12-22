// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyIVector.h 
//! \author Cem Yuksel
//! 
//! \brief  2D, 3D, 4D, and ND integer vector classes.
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

#ifndef _CY_IVECTOR_H_INCLUDED_
#define _CY_IVECTOR_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyCore.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

// Forward declarations
//!	\cond HIDDEN_SYMBOLS
template <typename T> class IVec2;
template <typename T> class IVec3;
template <typename T> class IVec4;
//! \endcond

//-------------------------------------------------------------------------------

//! A general class for N-dimensional integer vectors.

template <typename T, int N>
class IVec
{
	CY_NODISCARD friend IVec operator - ( T v, IVec const &p ) { return IVec<T,N>(v)-p; }	//!< Subtraction from a constant
	CY_NODISCARD friend IVec operator + ( T v, IVec const &p ) { return p+v; }				//!< Addition with a constant
	CY_NODISCARD friend IVec operator * ( T v, IVec const &p ) { return p*v; }				//!< Multiplication with a constant

public:

	//!@name Components of the vector
	T elem[N];

	//!@name Constructors
	IVec() CY_CLASS_FUNCTION_DEFAULT
	explicit IVec( T const * restrict p ) { MemCopy(elem,p,N); }
	explicit IVec( T v ) { for ( int i=0; i<N; ++i ) elem[i]=v; }
	template <typename S> explicit IVec( IVec<S,N> const &p ) { MemConvert(elem,p.elem,N); }
	template <int M> explicit IVec( IVec<T,M> const &p )
	{
		if ( N <= M ) { MemCopy(elem,p.elem,N); }
		else          { MemCopy(elem,p.elem,M); MemClear(elem,N-M); }
	}
	template <typename S, int M> explicit IVec( IVec<S,M> const &p )
	{
		if ( N <= M ) { MemConvert(elem,p.elem,N); }
		else          { MemConvert(elem,p.elem,M); MemClear(elem,N-M); }
	}
	explicit IVec( IVec2<T> const &p );
	explicit IVec( IVec3<T> const &p );
	explicit IVec( IVec4<T> const &p );
	template <typename S> explicit IVec( IVec2<S> const &p );
	template <typename S> explicit IVec( IVec3<S> const &p );
	template <typename S> explicit IVec( IVec4<S> const &p );
	template <typename P> explicit IVec( P const &p ) { for ( int i=0; i<N; ++i ) elem[i]=(T)p[i]; }

	//!@name Set & Get value methods
	void Zero()                      { MemClear(elem,N); }						//!< Sets the coordinates as zero
	void Get( T * restrict p ) const { MemCopy(p,elem,N); }						//!< Puts the coordinate values into the array
	void Set( T const * restrict p ) { MemCopy(elem,p,N); }						//!< Sets the coordinates using the values in the given array
	void Set( T v )                  { for ( int i=0; i<N; ++i ) elem[i] = v; }	//!< Sets all coordinates using the given value
	template <int M> void CopyData( T * restrict p ) { if ( M <= N ) { MemCopy(p,elem,M); } else { MemCopy(p,elem,N); MemClear(p+N,M-N); }	}
	template <typename S, int M> void ConvertData( S * restrict p ) { if ( M <= N ) { MemConvert(p,elem,M); } else { MemConvert(p,elem,N); MemClear(p+N,M-N); }	}

	//!@name General methods
	CY_NODISCARD T    Sum    () const { T v=elem[0]; for ( int i=1; i<N; ++i ) v+=elem[i]; return v; }		//!< Returns the sum of its components
	CY_NODISCARD bool IsZero () const { for ( int i=0; i<N; ++i ) if ( elem[i] != T(0) ) return false; return true; }	//!< Returns true if all components are exactly zero
	CY_NODISCARD T    Min    () const { T m = elem[0]; for ( int i=1; i<N; ++i ) if ( m > elem[i] ) m = elem[i]; return m; }
	CY_NODISCARD T    Max    () const { T m = elem[0]; for ( int i=1; i<N; ++i ) if ( m < elem[i] ) m = elem[i]; return m; }
	CY_NODISCARD int  MinComp() const { T m = elem[0]; int ix=0; for ( int i=1; i<N; ++i ) if ( m > elem[i] ) { m = elem[i]; ix = i; } return ix; }
	CY_NODISCARD int  MaxComp() const { T m = elem[0]; int ix=0; for ( int i=1; i<N; ++i ) if ( m < elem[i] ) { m = elem[i]; ix = i; } return ix; }

	//!@name Limit methods
	void Clamp   ( T minValue, T maxValue ) { ClampMin(minValue); ClampMax(maxValue); }		//!< Ensures that all components of the vector are within the given limits.
	void ClampMin( T v ) { for ( int i=0; i<N; ++i ) elem[i] = (elem[i]<v) ? v : elem[i]; }	//!< Ensures that all components of the vector are greater than or equal to the given limit.
	void ClampMax( T v ) { for ( int i=0; i<N; ++i ) elem[i] = (elem[i]>v) ? v : elem[i]; }	//!< Ensures that all components of the vector are smaller than or equal to the given limit.
	void SetAbs  ()      { for ( int i=0; i<N; i++ ) elem[i] = std::abs(elem[i]); }			//!< Converts all negative components to positive values

	//!@name Unary operators
	CY_NODISCARD IVec operator - () const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i]=-elem[i]; return r; } 

	//!@name Binary operators
	CY_NODISCARD IVec operator + ( IVec const &p ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] + p.elem[i]; return r; }
	CY_NODISCARD IVec operator - ( IVec const &p ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] - p.elem[i]; return r; }
	CY_NODISCARD IVec operator * ( IVec const &p ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] * p.elem[i]; return r; }
	CY_NODISCARD IVec operator / ( IVec const &p ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] / p.elem[i]; return r; }
	CY_NODISCARD IVec operator + ( T    const  v ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] + v; return r; }
	CY_NODISCARD IVec operator - ( T    const  v ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] - v; return r; }
	CY_NODISCARD IVec operator * ( T    const  v ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] * v; return r; }
	CY_NODISCARD IVec operator / ( T    const  v ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] / v; return r; }

	//!@name Assignment operators
	const IVec& operator += ( IVec const &p ) { for ( int i=0; i<N; ++i ) elem[i] += p.elem[i]; return *this; }
	const IVec& operator -= ( IVec const &p ) { for ( int i=0; i<N; ++i ) elem[i] -= p.elem[i]; return *this; }
	const IVec& operator *= ( IVec const &p ) { for ( int i=0; i<N; ++i ) elem[i] *= p.elem[i]; return *this; }
	const IVec& operator /= ( IVec const &p ) { for ( int i=0; i<N; ++i ) elem[i] /= p.elem[i]; return *this; }
	const IVec& operator += ( T    const  v ) { for ( int i=0; i<N; ++i ) elem[i] += v; return *this; }
	const IVec& operator -= ( T    const  v ) { for ( int i=0; i<N; ++i ) elem[i] -= v; return *this; }
	const IVec& operator *= ( T    const  v ) { for ( int i=0; i<N; ++i ) elem[i] *= v; return *this; }
	const IVec& operator /= ( T    const  v ) { for ( int i=0; i<N; ++i ) elem[i] /= v; return *this; }

	//!@name Bitwise operators
	CY_NODISCARD IVec operator << ( IVec const &p ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] << p.elem[i]; return r; }
	CY_NODISCARD IVec operator >> ( IVec const &p ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] >> p.elem[i]; return r; }
	CY_NODISCARD IVec operator  & ( IVec const &p ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i]  & p.elem[i]; return r; }
	CY_NODISCARD IVec operator  | ( IVec const &p ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i]  | p.elem[i]; return r; }
	CY_NODISCARD IVec operator  ^ ( IVec const &p ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i]  ^ p.elem[i]; return r; }
	CY_NODISCARD IVec operator << ( T    const  v ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] << v; return r; }
	CY_NODISCARD IVec operator >> ( T    const  v ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] >> v; return r; }
	CY_NODISCARD IVec operator  & ( T    const  v ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i]  & v; return r; }
	CY_NODISCARD IVec operator  | ( T    const  v ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i]  | v; return r; }
	CY_NODISCARD IVec operator  ^ ( T    const  v ) const { IVec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i]  ^ v; return r; }

	//!@name Bitwise Assignment operators
	const IVec& operator <<= ( IVec const &p ) { for ( int i=0; i<N; ++i ) elem[i] <<= p.elem[i]; return *this; }
	const IVec& operator >>= ( IVec const &p ) { for ( int i=0; i<N; ++i ) elem[i] >>= p.elem[i]; return *this; }
	const IVec& operator  &= ( IVec const &p ) { for ( int i=0; i<N; ++i ) elem[i]  &= p.elem[i]; return *this; }
	const IVec& operator  |= ( IVec const &p ) { for ( int i=0; i<N; ++i ) elem[i]  |= p.elem[i]; return *this; }
	const IVec& operator  ^= ( IVec const &p ) { for ( int i=0; i<N; ++i ) elem[i]  ^= p.elem[i]; return *this; }
	const IVec& operator <<= ( T    const  v ) { for ( int i=0; i<N; ++i ) elem[i] <<= v; return *this; }
	const IVec& operator >>= ( T    const  v ) { for ( int i=0; i<N; ++i ) elem[i] >>= v; return *this; }
	const IVec& operator  &= ( T    const  v ) { for ( int i=0; i<N; ++i ) elem[i]  &= v; return *this; }
	const IVec& operator  |= ( T    const  v ) { for ( int i=0; i<N; ++i ) elem[i]  |= v; return *this; }
	const IVec& operator  ^= ( T    const  v ) { for ( int i=0; i<N; ++i ) elem[i]  ^= v; return *this; }

	//!@name Test operators
	CY_NODISCARD bool operator == ( IVec const &p ) const { for ( int i=0; i<N; ++i ) if ( elem[i] != p.elem[i] ) return false; return true; }
	CY_NODISCARD bool operator != ( IVec const &p ) const { for ( int i=0; i<N; ++i ) if ( elem[i] != p.elem[i] ) return true; return false; }

	//!@name Access operators
	CY_NODISCARD T& operator [] ( int i )       { return elem[i]; }
	CY_NODISCARD T  operator [] ( int i ) const { return elem[i]; }

	//!@name Dot product
	CY_NODISCARD T Dot        ( IVec const &p ) const { IVec r=operator*(p); return r.Sum(); }	//!< Dot product
	CY_NODISCARD T operator % ( IVec const &p ) const { return Dot(p); }						//!< Dot product operator
};

//-------------------------------------------------------------------------------

//! 2D integer vector class

template <typename T>
class IVec2
{
	CY_NODISCARD friend IVec2 operator - ( T v, IVec2 const &p ) { return IVec2<T>(v)-p; }	//!< Subtraction from a constant
	CY_NODISCARD friend IVec2 operator + ( T v, IVec2 const &p ) { return p+v; }			//!< Addition with a constant
	CY_NODISCARD friend IVec2 operator * ( T v, IVec2 const &p ) { return p*v; }			//!< Multiplication with a constant

public:

	//!@name Components of the vector
	union {
		struct { T x, y; };
		T elem[2];
	};

	//!@name Constructors
	IVec2() CY_CLASS_FUNCTION_DEFAULT

	IVec2( T _x, T _y )   : x( _x), y( _y) {}
	explicit IVec2( T v ) : x(v  ), y(v  ) {}
	explicit IVec2( IVec3<T> const &p );
	explicit IVec2( IVec4<T> const &p );
	explicit IVec2( T const * restrict v ) { Set( v ); }
	template <typename S> explicit IVec2( IVec2<S> const &p ) : x(T(p.x)), y(T(p.y)) {}
	template <typename S> explicit IVec2( IVec3<S> const &p );
	template <typename S> explicit IVec2( IVec4<S> const &p );
	template <            int M> explicit IVec2( IVec<T,M> const &p ) { p.CopyData<2>(elem); }
	template <typename S, int M> explicit IVec2( IVec<S,M> const &p ) { p.ConvertData<T,2>(elem); }
	template <typename P> explicit IVec2( P const &p ) : x(T(p[0])), y(T(p[1])) {}

	//!@name Conversion
#ifdef _CY_VECTOR_H_INCLUDED_
	template <typename T> explicit operator Vec2<T> () const { return Vec2<T>(T(x),T(y)); }
#endif

	//!@name Set & Get value methods
	void Zero()                      { MemClear(elem,2); }				//!< Sets the coordinates as zero.
	void Get( T * restrict p ) const { ((IVec2*)p)->operator=(*this); }	//!< Puts the coordinate values into the array.
	void Set( T const * restrict p ) { operator=(*((IVec2*)p)); }		//!< Sets the coordinates using the values in the given array.
	void Set( T v )                  { x=v; y=v; }						//!< Sets all coordinates using the given value
	void Set( T _x, T _y )           { x=_x; y=_y; }					//!< Sets the coordinates using the given values

	//!@name General methods
	CY_NODISCARD T     Sum     () const { return x+y; }						//!< Returns the sum of its components
	CY_NODISCARD bool  IsZero  () const { return x==T(0) && y==T(0); }		//!< Returns true if all components are exactly zero
	CY_NODISCARD T     Min     () const { return x<y ? x : y; }				//!< Returns the minimum component of the vector.
	CY_NODISCARD T     Max     () const { return x>y ? x : y; }				//!< Returns the maximum component of the vector.
	CY_NODISCARD int   MinComp () const { return x<y ? 0 : 1; }				//!< Returns the index of the minimum component of the vector.
	CY_NODISCARD int   MaxComp () const { return x>y ? 0 : 1; }				//!< Returns the index of the maximum component of the vector.
	CY_NODISCARD IVec2 Abs     () const { return IVec2(std::abs(x),std::abs(y)); }	//!< Returns a vector containing the absolute values of all components.
	CY_NODISCARD IVec2 SortAsc () const { IVec2 v; Sort2<true >( v.elem, elem ); return v; }	//!< Returns a vector with components sorted in ascending order.
	CY_NODISCARD IVec2 SortDesc() const { IVec2 v; Sort2<false>( v.elem, elem ); return v; }	//!< Returns a vector with components sorted in descending order.

	//!@name Limit methods
	void Clamp   ( T minValue, T maxValue ) { ClampMin(minValue); ClampMax(maxValue); }	//!< Ensures that all components of the vector are within the given limits.
	void ClampMin( T v ) { x=(x<v)?v:x; y=(y<v)?v:y; }									//!< Ensures that all components of the vector are greater than or equal to the given limit.
	void ClampMax( T v ) { x=(x>v)?v:x; y=(y>v)?v:y; }									//!< Ensures that all components of the vector are smaller than or equal to the given limit.
	void SetAbs  ()      { x=std::abs(x); y=std::abs(y); }								//!< Converts all negative components to positive values

	//!@name Unary operators
	CY_NODISCARD IVec2 operator - () const { IVec2 r; r.x=-x; r.y=-y; return r; } 

	//!@name Binary operators
	CY_NODISCARD IVec2 operator + ( IVec2 const &p ) const { IVec2 r; r.x=x+p.x; r.y=y+p.y; return r; }
	CY_NODISCARD IVec2 operator - ( IVec2 const &p ) const { IVec2 r; r.x=x-p.x; r.y=y-p.y; return r; }
	CY_NODISCARD IVec2 operator * ( IVec2 const &p ) const { IVec2 r; r.x=x*p.x; r.y=y*p.y; return r; }
	CY_NODISCARD IVec2 operator / ( IVec2 const &p ) const { IVec2 r; r.x=x/p.x; r.y=y/p.y; return r; }
	CY_NODISCARD IVec2 operator + ( T     const  v ) const { IVec2 r; r.x=x+v;   r.y=y+v;   return r; }
	CY_NODISCARD IVec2 operator - ( T     const  v ) const { IVec2 r; r.x=x-v;   r.y=y-v;   return r; }
	CY_NODISCARD IVec2 operator * ( T     const  v ) const { IVec2 r; r.x=x*v;   r.y=y*v;   return r; }
	CY_NODISCARD IVec2 operator / ( T     const  v ) const { IVec2 r; r.x=x/v;   r.y=y/v;   return r; }

	//!@name Assignment operators
	IVec2 const& operator += ( IVec2 const &p ) { x+=p.x; y+=p.y; return *this; }
	IVec2 const& operator -= ( IVec2 const &p ) { x-=p.x; y-=p.y; return *this; }
	IVec2 const& operator *= ( IVec2 const &p ) { x*=p.x; y*=p.y; return *this; }
	IVec2 const& operator /= ( IVec2 const &p ) { x/=p.x; y/=p.y; return *this; }
	IVec2 const& operator += ( T     const  v ) { x+=v;   y+=v;   return *this; }
	IVec2 const& operator -= ( T     const  v ) { x-=v;   y-=v;   return *this; }
	IVec2 const& operator *= ( T     const  v ) { x*=v;   y*=v;   return *this; }
	IVec2 const& operator /= ( T     const  v ) { x/=v;   y/=v;   return *this; }

	//!@name Bitwise operators
	CY_NODISCARD IVec2 operator << ( IVec2 const &p ) const { IVec2 r; r.x = x << p.x; r.y = y << p.y; return r; }
	CY_NODISCARD IVec2 operator >> ( IVec2 const &p ) const { IVec2 r; r.x = x >> p.x; r.y = y >> p.y; return r; }
	CY_NODISCARD IVec2 operator  & ( IVec2 const &p ) const { IVec2 r; r.x = x  & p.x; r.y = y  & p.y; return r; }
	CY_NODISCARD IVec2 operator  | ( IVec2 const &p ) const { IVec2 r; r.x = x  | p.x; r.y = y  | p.y; return r; }
	CY_NODISCARD IVec2 operator  ^ ( IVec2 const &p ) const { IVec2 r; r.x = x  ^ p.x; r.y = y  ^ p.y; return r; }
	CY_NODISCARD IVec2 operator << ( T     const  v ) const { IVec2 r; r.x = x << v;   r.y = y << v;   return r; }
	CY_NODISCARD IVec2 operator >> ( T     const  v ) const { IVec2 r; r.x = x >> v;   r.y = y >> v;   return r; }
	CY_NODISCARD IVec2 operator  & ( T     const  v ) const { IVec2 r; r.x = x  & v;   r.y = y  & v;   return r; }
	CY_NODISCARD IVec2 operator  | ( T     const  v ) const { IVec2 r; r.x = x  | v;   r.y = y  | v;   return r; }
	CY_NODISCARD IVec2 operator  ^ ( T     const  v ) const { IVec2 r; r.x = x  ^ v;   r.y = y  ^ v;   return r; }

	//!@name Bitwise Assignment operators
	IVec2 const& operator <<= ( IVec2 const &p ) { x<<=p.x; y<<=p.y; return *this; }
	IVec2 const& operator >>= ( IVec2 const &p ) { x>>=p.x; y>>=p.y; return *this; }
	IVec2 const& operator  &= ( IVec2 const &p ) { x &=p.x; y &=p.y; return *this; }
	IVec2 const& operator  |= ( IVec2 const &p ) { x |=p.x; y |=p.y; return *this; }
	IVec2 const& operator  ^= ( IVec2 const &p ) { x ^=p.x; y ^=p.y; return *this; }
	IVec2 const& operator <<= ( T     const  v ) { x<<=v;   y<<=v;   return *this; }
	IVec2 const& operator >>= ( T     const  v ) { x>>=v;   y>>=v;   return *this; }
	IVec2 const& operator  &= ( T     const  v ) { x &=v;   y &=v;   return *this; }
	IVec2 const& operator  |= ( T     const  v ) { x |=v;   y |=v;   return *this; }
	IVec2 const& operator  ^= ( T     const  v ) { x ^=v;   y ^=v;   return *this; }

	//!@name Test operators
	CY_NODISCARD bool operator == ( IVec2 const &p ) const { return x==p.x && y==p.y; }
	CY_NODISCARD bool operator != ( IVec2 const &p ) const { return x!=p.x && y!=p.y; }

	//!@name Access operators
	CY_NODISCARD T&       operator [] ( int i )       { return Element(i); }
	CY_NODISCARD T const& operator [] ( int i ) const { return Element(i); }
	CY_NODISCARD T&       Element     ( int i )       { return elem[i]; }
	CY_NODISCARD T const& Element     ( int i ) const { return elem[i]; }
	CY_NODISCARD T*       Data        ()              { return elem; }
	CY_NODISCARD T const* Data        ()        const { return elem; }

	//!@name Cross product and dot product
	CY_NODISCARD T Dot        ( IVec2 const &p ) const { IVec2 r=operator*(p); return r.Sum(); }	//!< Dot product
	CY_NODISCARD T operator % ( IVec2 const &p ) const { return Dot(p); }							//!< Dot product operator
};

//-------------------------------------------------------------------------------

//! 3D integer vector class

template <typename T>
class IVec3
{
	CY_NODISCARD friend IVec3 operator - ( T v, IVec3 const &p ) { return IVec3<T>(v)-p; }	//!< Subtraction from a constant
	CY_NODISCARD friend IVec3 operator + ( T v, IVec3 const &p ) { return p+v; }			//!< Addition with a constant
	CY_NODISCARD friend IVec3 operator * ( T v, IVec3 const &p ) { return p*v; }			//!< Multiplication with a constant

public:

	//!@name Components of the vector
	union {
		struct { T x, y, z; };
		T elem[3];
	};

	//!@name Constructors
	IVec3() CY_CLASS_FUNCTION_DEFAULT
	IVec3( T _x, T _y, T _z )                   : x( _x), y( _y), z( _z) {}
	explicit IVec3( T v )                       : x(v  ), y(v  ), z(v  ) {}
	explicit IVec3( IVec2<T> const &p, T _z=0 ) : x(p.x), y(p.y), z( _z) {}
	explicit IVec3( IVec4<T> const &p );
	explicit IVec3( T const * restrict v ) { Set( v ); }
	template <typename S> explicit IVec3( IVec3<S> const &p )         : x(T(p.x)), y(T(p.y)), z(T(p.z)) {}
	template <typename S> explicit IVec3( IVec2<S> const &p, T _z=0 ) : x(T(p.x)), y(T(p.y)), z(   _z ) {}
	template <typename S> explicit IVec3( IVec4<S> const &p );
	template <            int M> explicit IVec3( IVec<T,M> const &p ) { p.CopyData<3>(elem); }
	template <typename S, int M> explicit IVec3( IVec<S,M> const &p ) { p.ConvertData<T,3>(elem); }
	template <typename P> explicit IVec3( P const &p ) : x((T)p[0]), y((T)p[1]), z((T)p[2]) {}

	//!@name Conversion
#ifdef _CY_VECTOR_H_INCLUDED_
	template <typename T> explicit operator Vec3<T> () const { return Vec3<T>(T(x),T(y),T(z)); }
#endif

	//!@name Set & Get value methods
	void Zero()                      { MemClear(elem,3); }				//!< Sets the coordinates as zero
	void Get( T * restrict p ) const { ((IVec3*)p)->operator=(*this); }	//!< Puts the coordinate values into the array
	void Set( T const * restrict p ) { operator=(*((IVec3*)p)); }		//!< Sets the coordinates using the values in the given array
	void Set( T v )                  { x=v; y=v; z=v; }					//!< Sets all coordinates using the given value
	void Set( T _x, T _y, T _z )     { x=_x; y=_y; z=_z; }				//!< Sets the coordinates using the given values

	//!@name Length and Normalize methods
	CY_NODISCARD T     Sum     () const { return x+y+z; }										//!< Returns the sum of its components
	CY_NODISCARD bool  IsZero  () const { return x==T(0) && y==T(0) && z==T(0); }				//!< Returns true if all components are exactly zero
	CY_NODISCARD T     Min     () const { return elem[MinComp()]; }								//!< Returns the minimum component of the vector.
	CY_NODISCARD T     Max     () const { return elem[MaxComp()]; }								//!< Returns the maximum component of the vector.
	CY_NODISCARD int   MinComp () const { int yx=y<x; int zx=z<x; int zy=z<y; return (yx|zx)+(zx&zy); }	//!< Returns the index of the minimum component of the vector.
	CY_NODISCARD int   MaxComp () const { int xy=x<y; int xz=x<z; int yz=y<z; return (xy|xz)+(xz&yz); }	//!< Returns the index of the maximum component of the vector.
	CY_NODISCARD IVec3 Abs     () const { return IVec3(std::abs(x),std::abs(y),std::abs(z)); }	//!< Returns a vector containing the absolute values of all components.
	CY_NODISCARD IVec3 SortAsc () const { IVec3 v; Sort3<true >( v.elem, elem ); return v; }	//!< Returns a vector with components sorted in ascending order.
	CY_NODISCARD IVec3 SortDesc() const { IVec3 v; Sort3<false>( v.elem, elem ); return v; }	//!< Returns a vector with components sorted in descending order.

	//!@name Limit methods
	void Clamp   ( T minValue, T maxValue ) { ClampMin(minValue); ClampMax(maxValue); }	//!< Ensures that all components of the vector are within the given limits.
	void ClampMin( T v ) { x=(x<v)?v:x; y=(y<v)?v:y; z=(z<v)?v:z; }						//!< Ensures that all components of the vector are greater than or equal to the given limit.
	void ClampMax( T v ) { x=(x>v)?v:x; y=(y>v)?v:y; z=(z>v)?v:z; }						//!< Ensures that all components of the vector are smaller than or equal to the given limit.
	void SetAbs  ()      { x=std::abs(x); y=std::abs(y); z=std::abs(z); }				//!< Converts all negative components to positive values

	//!@name Unary operators
	CY_NODISCARD IVec3 operator - () const { IVec3 r; r.x=-x; r.y=-y; r.z=-z; return r; } 

	//!@name Binary operators
	CY_NODISCARD IVec3 operator + ( IVec3 const &p ) const { IVec3 r; r.x=x+p.x; r.y=y+p.y; r.z=z+p.z; return r; }
	CY_NODISCARD IVec3 operator - ( IVec3 const &p ) const { IVec3 r; r.x=x-p.x; r.y=y-p.y; r.z=z-p.z; return r; }
	CY_NODISCARD IVec3 operator * ( IVec3 const &p ) const { IVec3 r; r.x=x*p.x; r.y=y*p.y; r.z=z*p.z; return r; }
	CY_NODISCARD IVec3 operator / ( IVec3 const &p ) const { IVec3 r; r.x=x/p.x; r.y=y/p.y; r.z=z/p.z; return r; }
	CY_NODISCARD IVec3 operator + ( T     const  v ) const { IVec3 r; r.x=x+v;   r.y=y+v;   r.z=z+v;   return r; }
	CY_NODISCARD IVec3 operator - ( T     const  v ) const { IVec3 r; r.x=x-v;   r.y=y-v;   r.z=z-v;   return r; }
	CY_NODISCARD IVec3 operator * ( T     const  v ) const { IVec3 r; r.x=x*v;   r.y=y*v;   r.z=z*v;   return r; }
	CY_NODISCARD IVec3 operator / ( T     const  v ) const { IVec3 r; r.x=x/v;   r.y=y/v;   r.z=z/v;   return r; }

	//!@name Assignment operators
	IVec3 const& operator += ( IVec3 const &p ) { x+=p.x; y+=p.y; z+=p.z; return *this; }
	IVec3 const& operator -= ( IVec3 const &p ) { x-=p.x; y-=p.y; z-=p.z; return *this; }
	IVec3 const& operator *= ( IVec3 const &p ) { x*=p.x; y*=p.y; z*=p.z; return *this; }
	IVec3 const& operator /= ( IVec3 const &p ) { x/=p.x; y/=p.y; z/=p.z; return *this; }
	IVec3 const& operator += ( T     const  v ) { x+=v;   y+=v;   z+=v;   return *this; }
	IVec3 const& operator -= ( T     const  v ) { x-=v;   y-=v;   z-=v;   return *this; }
	IVec3 const& operator *= ( T     const  v ) { x*=v;   y*=v;   z*=v;   return *this; }
	IVec3 const& operator /= ( T     const  v ) { x/=v;   y/=v;   z/=v;   return *this; }

	//!@name Bitwise operators
	CY_NODISCARD IVec3 operator << ( IVec3 const &p ) const { IVec3 r; r.x = x << p.x; r.y = y << p.y; r.z = z << p.z; return r; }
	CY_NODISCARD IVec3 operator >> ( IVec3 const &p ) const { IVec3 r; r.x = x >> p.x; r.y = y >> p.y; r.z = z >> p.z; return r; }
	CY_NODISCARD IVec3 operator  & ( IVec3 const &p ) const { IVec3 r; r.x = x  & p.x; r.y = y  & p.y; r.z = z  & p.z; return r; }
	CY_NODISCARD IVec3 operator  | ( IVec3 const &p ) const { IVec3 r; r.x = x  | p.x; r.y = y  | p.y; r.z = z  | p.z; return r; }
	CY_NODISCARD IVec3 operator  ^ ( IVec3 const &p ) const { IVec3 r; r.x = x  ^ p.x; r.y = y  ^ p.y; r.z = z  ^ p.z; return r; }
	CY_NODISCARD IVec3 operator << ( T     const  v ) const { IVec3 r; r.x = x << v;   r.y = y << v;   r.z = z << v;   return r; }
	CY_NODISCARD IVec3 operator >> ( T     const  v ) const { IVec3 r; r.x = x >> v;   r.y = y >> v;   r.z = z >> v;   return r; }
	CY_NODISCARD IVec3 operator  & ( T     const  v ) const { IVec3 r; r.x = x  & v;   r.y = y  & v;   r.z = z  & v;   return r; }
	CY_NODISCARD IVec3 operator  | ( T     const  v ) const { IVec3 r; r.x = x  | v;   r.y = y  | v;   r.z = z  | v;   return r; }
	CY_NODISCARD IVec3 operator  ^ ( T     const  v ) const { IVec3 r; r.x = x  ^ v;   r.y = y  ^ v;   r.z = z  ^ v;   return r; }

	//!@name Bitwise Assignment operators
	IVec3 const& operator <<= ( IVec3 const &p ) { x<<=p.x; y<<=p.y; z<<=p.z; return *this; }
	IVec3 const& operator >>= ( IVec3 const &p ) { x>>=p.x; y>>=p.y; z>>=p.z; return *this; }
	IVec3 const& operator  &= ( IVec3 const &p ) { x &=p.x; y &=p.y; z &=p.z; return *this; }
	IVec3 const& operator  |= ( IVec3 const &p ) { x |=p.x; y |=p.y; z |=p.z; return *this; }
	IVec3 const& operator  ^= ( IVec3 const &p ) { x ^=p.x; y ^=p.y; z ^=p.z; return *this; }
	IVec3 const& operator <<= ( T     const  v ) { x<<=v;   y<<=v;   z<<=v;   return *this; }
	IVec3 const& operator >>= ( T     const  v ) { x>>=v;   y>>=v;   z>>=v;   return *this; }
	IVec3 const& operator  &= ( T     const  v ) { x &=v;   y &=v;   z &=v;   return *this; }
	IVec3 const& operator  |= ( T     const  v ) { x |=v;   y |=v;   z |=v;   return *this; }
	IVec3 const& operator  ^= ( T     const  v ) { x ^=v;   y ^=v;   z ^=v;   return *this; }

	//!@name Test operators
	CY_NODISCARD bool operator == ( IVec3 const& p ) const { return x==p.x && y==p.y && z==p.z; }
	CY_NODISCARD bool operator != ( IVec3 const& p ) const { return x!=p.x && y!=p.y && z!=p.z; }

	//!@name Access operators
	CY_NODISCARD T&       operator [] ( int i )       { return Element(i); }
	CY_NODISCARD T const& operator [] ( int i ) const { return Element(i); }
	CY_NODISCARD T&       Element     ( int i )       { return elem[i]; }
	CY_NODISCARD T const& Element     ( int i ) const { return elem[i]; }
	CY_NODISCARD T*       Data        ()              { return elem; }
	CY_NODISCARD T const* Data        ()        const { return elem; }

	//!@name Cross product and dot product
	CY_NODISCARD T   Dot        ( IVec3 const &p ) const { IVec3 r=operator*(p); return r.Sum(); }	//!< Dot product
	CY_NODISCARD T   operator % ( IVec3 const &p ) const { return Dot(p); }							//!< Dot product

	//!@name Conversion Methods
	CY_NODISCARD IVec2<T> XY() const { return IVec2<T>(x,y); }
};

//-------------------------------------------------------------------------------

//! 4D integer vector class

template <typename T>
class IVec4
{
	CY_NODISCARD friend IVec4 operator - ( T v, IVec4 const &p ) { return IVec4(v)-p; }	//!< Subtraction from a constant
	CY_NODISCARD friend IVec4 operator + ( T v, IVec4 const &p ) { return p+v; }		//!< Addition with a constant
	CY_NODISCARD friend IVec4 operator * ( T v, IVec4 const &p ) { return p*v; }		//!< Multiplication with a constant

public:

	//!@name Components of the vector
	union {
		struct { T x, y, z, w; };
		T elem[4];
	};

	//!@name Constructors
	IVec4() CY_CLASS_FUNCTION_DEFAULT
	IVec4( T _x, T _y, T _z,   T _w )                   : x( _x), y( _y), z( _z), w( _w) {}
	explicit IVec4( T v )                               : x(v  ), y(v  ), z(v  ), w(v  ) {}
	explicit IVec4( IVec3<T> const &p,         T _w=0 ) : x(p.x), y(p.y), z(p.z), w( _w) {}
	explicit IVec4( IVec2<T> const &p, T _z=0, T _w=0 ) : x(p.x), y(p.y), z( _z), w( _w) {}
	explicit IVec4( T const * restrict v ) { Set( v ); }
	template <typename S> explicit IVec4( IVec4<S> const &p )                  : x(T(p.x)), y(T(p.y)), z(T(p.z)), z(T(p.w)) {}
	template <typename S> explicit IVec4( IVec3<S> const &p,         T _w=0 )  : x(T(p.x)), y(T(p.y)), z(T(p.z)), z(   _w ) {}
	template <typename S> explicit IVec4( IVec2<S> const &p, T _z=0, T _w=0 )  : x(T(p.x)), y(T(p.y)), z(   _z ), z(   _w ) {}
	template <            int M> explicit IVec4( IVec<T,M> const &p ) { p.CopyData<3>(elem); }
	template <typename S, int M> explicit IVec4( IVec<S,M> const &p ) { p.ConvertData<T,3>(elem); }
	template <typename P> explicit IVec4( P const &p ) : x((T)p[0]), y((T)p[1]), z((T)p[2]) {}

	//!@name Conversion
#ifdef _CY_VECTOR_H_INCLUDED_
	template <typename T> explicit operator Vec4<T> () const { return Vec4<T>(T(x),T(y),T(z),T(w)); }
#endif

	//!@name Set & Get value methods
	void Zero()                            { MemClear(elem,4); }				//!< Sets the coordinates as zero
	void Get( T       * restrict p ) const { ((IVec4*)p)->operator=(*this); }	//!< Puts the coordinate values into the array
	void Set( T const * restrict p )       { operator=(*((IVec4*)p)); }			//!< Sets the coordinates using the values in the given array
	void Set( T v )                        { x=v; y=v; z=v; w=v; }				//!< Sets all coordinates using the given value
	void Set( T _x, T _y, T _z, T _w=0 )   { x=_x; y=_y; z=_z; w=_w; }			//!< Sets the coordinates using the given values

	//!@name Length and Normalize methods
	CY_NODISCARD T     Sum     () const { return x+y+z+w; }										//!< Returns the sum of its components
	CY_NODISCARD bool  IsZero  () const { return x==T(0) && y==T(0) && z==T(0) && w==T(0); }	//!< Returns true if all components are exactly zero
	CY_NODISCARD T     Min     () const { T   mxy = x<y ? x : y; T   mzw = z<w ? z : w; return mxy<mzw ? mxy : mzw; }				//!< Returns the minimum component of the vector.
	CY_NODISCARD T     Max     () const { T   mxy = x>y ? x : y; T   mzw = z>w ? z : w; return mxy>mzw ? mxy : mzw; }				//!< Returns the maximum component of the vector.
	CY_NODISCARD int   MinComp () const { int ixy = x<y ? 0 : 1; int izw = z<w ? 2 : 3; return elem[ixy]<elem[izw] ? ixy : izw; }	//!< Returns the index of the minimum component of the vector.
	CY_NODISCARD int   MaxComp () const { int ixy = x>y ? 0 : 1; int izw = z>w ? 2 : 3; return elem[ixy]>elem[izw] ? ixy : izw; }	//!< Returns the index of the maximum component of the vector.
	CY_NODISCARD IVec4 Abs     () const { return IVec4(std::abs(x),std::abs(y),std::abs(z),std::abs(w)); }	//!< Returns a vector containing the absolute values of all components.
	CY_NODISCARD IVec4 SortAsc () const { IVec4 v; Sort4<true >( v.elem, elem ); return v; }	//!< Returns a vector with components sorted in ascending order.
	CY_NODISCARD IVec4 SortDesc() const { IVec4 v; Sort4<false>( v.elem, elem ); return v; }	//!< Returns a vector with components sorted in descending order.

	//!@name Limit methods
	void Clamp   ( T minValue, T maxValue ) { ClampMin(minValue); ClampMax(maxValue); }		//!< Ensures that all components of the vector are within the given limits.
	void ClampMin( T v ) { x=(x<v)?v:x; y=(y<v)?v:y; z=(z<v)?v:z; w=(w<v)?v:w; }			//!< Ensures that all components of the vector are greater than or equal to the given limit.
	void ClampMax( T v ) { x=(x>v)?v:x; y=(y>v)?v:y; z=(z>v)?v:z; w=(w>v)?v:w; }			//!< Ensures that all components of the vector are smaller than or equal to the given limit.
	void SetAbs  ()      { x=std::abs(x); y=std::abs(y); z=std::abs(z); w=std::abs(w); }	//!< Converts all negative components to positive values

	//!@name Unary operators
	CY_NODISCARD IVec4 operator - () const { IVec4 r; r.x=-x; r.y=-y; r.z=-z; r.w=-w; return r; } 

	//!@name Binary operators
	CY_NODISCARD IVec4 operator + ( IVec4 const &p ) const { IVec4 r; r.x=x+p.x; r.y=y+p.y; r.z=z+p.z; r.w=w+p.w; return r; }
	CY_NODISCARD IVec4 operator - ( IVec4 const &p ) const { IVec4 r; r.x=x-p.x; r.y=y-p.y; r.z=z-p.z; r.w=w-p.w; return r; }
	CY_NODISCARD IVec4 operator * ( IVec4 const &p ) const { IVec4 r; r.x=x*p.x; r.y=y*p.y; r.z=z*p.z; r.w=w*p.w; return r; }
	CY_NODISCARD IVec4 operator / ( IVec4 const &p ) const { IVec4 r; r.x=x/p.x; r.y=y/p.y; r.z=z/p.z; r.w=w/p.w; return r; }
	CY_NODISCARD IVec4 operator + ( T     const  v ) const { IVec4 r; r.x=x+v;   r.y=y+v;   r.z=z+v;   r.w=w+v;   return r; }
	CY_NODISCARD IVec4 operator - ( T     const  v ) const { IVec4 r; r.x=x-v;   r.y=y-v;   r.z=z-v;   r.w=w-v;   return r; }
	CY_NODISCARD IVec4 operator * ( T     const  v ) const { IVec4 r; r.x=x*v;   r.y=y*v;   r.z=z*v;   r.w=w*v;   return r; }
	CY_NODISCARD IVec4 operator / ( T     const  v ) const { IVec4 r; r.x=x/v;   r.y=y/v;   r.z=z/v;   r.w=w/v;   return r; }

	//!@name Assignment operators
	IVec4 const& operator += ( IVec4 const &p ) { x+=p.x; y+=p.y; z+=p.z; w+=p.w; return *this; }
	IVec4 const& operator -= ( IVec4 const &p ) { x-=p.x; y-=p.y; z-=p.z; w-=p.w; return *this; }
	IVec4 const& operator *= ( IVec4 const &p ) { x*=p.x; y*=p.y; z*=p.z; w*=p.w; return *this; }
	IVec4 const& operator /= ( IVec4 const &p ) { x/=p.x; y/=p.y; z/=p.z; w/=p.w; return *this; }
	IVec4 const& operator += ( T     const  v ) { x+=v;   y+=v;   z+=v;   w+=v;   return *this; }
	IVec4 const& operator -= ( T     const  v ) { x-=v;   y-=v;   z-=v;   w-=v;   return *this; }
	IVec4 const& operator *= ( T     const  v ) { x*=v;   y*=v;   z*=v;   w*=v;   return *this; }
	IVec4 const& operator /= ( T     const  v ) { x/=v;   y/=v;   z/=v;   w/=v;   return *this; }

	//!@name Bitwise operators
	CY_NODISCARD IVec4 operator << ( IVec4 const &p ) const { IVec4 r; r.x = x << p.x; r.y = y << p.y; r.z = z << p.z; r.w = w << p.w; return r; }
	CY_NODISCARD IVec4 operator >> ( IVec4 const &p ) const { IVec4 r; r.x = x >> p.x; r.y = y >> p.y; r.z = z >> p.z; r.w = w >> p.w; return r; }
	CY_NODISCARD IVec4 operator  & ( IVec4 const &p ) const { IVec4 r; r.x = x  & p.x; r.y = y  & p.y; r.z = z  & p.z; r.w = w  & p.w; return r; }
	CY_NODISCARD IVec4 operator  | ( IVec4 const &p ) const { IVec4 r; r.x = x  | p.x; r.y = y  | p.y; r.z = z  | p.z; r.w = w  | p.w; return r; }
	CY_NODISCARD IVec4 operator  ^ ( IVec4 const &p ) const { IVec4 r; r.x = x  ^ p.x; r.y = y  ^ p.y; r.z = z  ^ p.z; r.w = w  ^ p.w; return r; }
	CY_NODISCARD IVec4 operator << ( T     const  v ) const { IVec4 r; r.x = x << v;   r.y = y << v;   r.z = z << v;   r.w = w << v;   return r; }
	CY_NODISCARD IVec4 operator >> ( T     const  v ) const { IVec4 r; r.x = x >> v;   r.y = y >> v;   r.z = z >> v;   r.w = w >> v;   return r; }
	CY_NODISCARD IVec4 operator  & ( T     const  v ) const { IVec4 r; r.x = x  & v;   r.y = y  & v;   r.z = z  & v;   r.w = w  & v;   return r; }
	CY_NODISCARD IVec4 operator  | ( T     const  v ) const { IVec4 r; r.x = x  | v;   r.y = y  | v;   r.z = z  | v;   r.w = w  | v;   return r; }
	CY_NODISCARD IVec4 operator  ^ ( T     const  v ) const { IVec4 r; r.x = x  ^ v;   r.y = y  ^ v;   r.z = z  ^ v;   r.w = w  ^ v;   return r; }

	//!@name Bitwise Assignment operators
	IVec4 const& operator <<= ( IVec4 const &p ) { x<<=p.x; y<<=p.y; z<<=p.z; w<<=p.w; return *this; }
	IVec4 const& operator >>= ( IVec4 const &p ) { x>>=p.x; y>>=p.y; z>>=p.z; w>>=p.w; return *this; }
	IVec4 const& operator  &= ( IVec4 const &p ) { x &=p.x; y &=p.y; z &=p.z; w &=p.w; return *this; }
	IVec4 const& operator  |= ( IVec4 const &p ) { x |=p.x; y |=p.y; z |=p.z; w |=p.w; return *this; }
	IVec4 const& operator  ^= ( IVec4 const &p ) { x ^=p.x; y ^=p.y; z ^=p.z; w ^=p.w; return *this; }
	IVec4 const& operator <<= ( T     const  v ) { x<<=v;   y<<=v;   z<<=v;   w<<=v;   return *this; }
	IVec4 const& operator >>= ( T     const  v ) { x>>=v;   y>>=v;   z>>=v;   w>>=v;   return *this; }
	IVec4 const& operator  &= ( T     const  v ) { x &=v;   y &=v;   z &=v;   w &=v;   return *this; }
	IVec4 const& operator  |= ( T     const  v ) { x |=v;   y |=v;   z |=v;   w |=v;   return *this; }
	IVec4 const& operator  ^= ( T     const  v ) { x ^=v;   y ^=v;   z ^=v;   w ^=v;   return *this; }

	//!@name Test operators
	CY_NODISCARD bool operator == ( IVec4 const &p ) const { return x==p.x && y==p.y && z==p.z && w==p.w; }
	CY_NODISCARD bool operator != ( IVec4 const &p ) const { return x!=p.x && y!=p.y && z!=p.z && w!=p.w; }

	//!@name Access operators
	CY_NODISCARD T&       operator [] ( int i )       { return Element(i); }
	CY_NODISCARD T const& operator [] ( int i ) const { return Element(i); }
	CY_NODISCARD T&       Element     ( int i )       { return elem[i]; }
	CY_NODISCARD T const& Element     ( int i ) const { return elem[i]; }
	CY_NODISCARD T*       Data        ()              { return elem; }
	CY_NODISCARD T const* Data        ()        const { return elem; }

	//!@name Cross product and dot product
	CY_NODISCARD T   Dot        ( IVec4 const &p ) const { IVec4 r=operator*(p); return r.Sum(); }	//!< Dot product
	CY_NODISCARD T   operator % ( IVec4 const &p ) const { return Dot(p); }							//!< Dot product

	//!@name Conversion Methods
	CY_NODISCARD IVec2<T> XY () const { return IVec2<T>(x,y); }
	CY_NODISCARD IVec3<T> XYZ() const { return IVec3<T>(*this); }
};

//-------------------------------------------------------------------------------

// Definitions of the conversion constructors
template <typename T, int N>                       IVec<T,N>::IVec( IVec2<T> const &p ) { if ( N <= 2 ) { MemCopy   (elem,&p.x,N); } else { MemCopy   (elem,&p.x,2); MemClear(elem,N-2); } }
template <typename T, int N>                       IVec<T,N>::IVec( IVec3<T> const &p ) { if ( N <= 3 ) { MemCopy   (elem,&p.x,N); } else { MemCopy   (elem,&p.x,3); MemClear(elem,N-3); } }
template <typename T, int N>                       IVec<T,N>::IVec( IVec4<T> const &p ) { if ( N <= 4 ) { MemCopy   (elem,&p.x,N); } else { MemCopy   (elem,&p.x,4); MemClear(elem,N-4); } }
template <typename T, int N> template <typename S> IVec<T,N>::IVec( IVec2<S> const &p ) { if ( N <= 2 ) { MemConvert(elem,&p.x,N); } else { MemConvert(elem,&p.x,2); MemClear(elem,N-2); } }
template <typename T, int N> template <typename S> IVec<T,N>::IVec( IVec3<S> const &p ) { if ( N <= 3 ) { MemConvert(elem,&p.x,N); } else { MemConvert(elem,&p.x,3); MemClear(elem,N-3); } }
template <typename T, int N> template <typename S> IVec<T,N>::IVec( IVec4<S> const &p ) { if ( N <= 4 ) { MemConvert(elem,&p.x,N); } else { MemConvert(elem,&p.x,4); MemClear(elem,N-4); } }
template <typename T>                              IVec2<T>::IVec2( IVec3<T> const &p ) : x(  p.x ), y(  p.y )            {}
template <typename T>                              IVec2<T>::IVec2( IVec4<T> const &p ) : x(  p.x ), y(  p.y )            {}
template <typename T>                              IVec3<T>::IVec3( IVec4<T> const &p ) : x(  p.x ), y(  p.y ), z(  p.z ) {}
template <typename T> template <typename S>        IVec2<T>::IVec2( IVec3<S> const &p ) : x(T(p.x)), y(T(p.y))            {}
template <typename T> template <typename S>        IVec2<T>::IVec2( IVec4<S> const &p ) : x(T(p.x)), y(T(p.y))            {}
template <typename T> template <typename S>        IVec3<T>::IVec3( IVec4<S> const &p ) : x(T(p.x)), y(T(p.y)), z(T(p.z)) {}

//-------------------------------------------------------------------------------

typedef IVec2< int8_t>  IVec2b;		//!< 8-bit  signed   byte    (int8_t)   2D integer vector class
typedef IVec3< int8_t>  IVec3b;		//!< 8-bit  signed   byte    (int8_t)   3D integer vector class
typedef IVec4< int8_t>  IVec4b;		//!< 8-bit  signed   byte    (int8_t)   4D integer vector class
typedef IVec2<uint8_t>  IVec2ub;	//!< 8-bit  unsigned byte    (uint8_t)  2D integer vector class
typedef IVec3<uint8_t>  IVec3ub;	//!< 8-bit  unsigned byte    (uint8_t)  3D integer vector class
typedef IVec4<uint8_t>  IVec4ub;	//!< 8-bit  unsigned byte    (uint8_t)  4D integer vector class

typedef IVec2< int16_t> IVec2s;		//!< 16-bit signed   short   (int16_t)  2D integer vector class
typedef IVec3< int16_t> IVec3s;		//!< 16-bit signed   short   (int16_t)  3D integer vector class
typedef IVec4< int16_t> IVec4s;		//!< 16-bit signed   short   (int16_t)  4D integer vector class
typedef IVec2<uint16_t> IVec2us;	//!< 16-bit unsigned short   (uint16_t) 2D integer vector class
typedef IVec3<uint16_t> IVec3us;	//!< 16-bit unsigned short   (uint16_t) 3D integer vector class
typedef IVec4<uint16_t> IVec4us;	//!< 16-bit unsigned short   (uint16_t) 4D integer vector class

typedef IVec2< int32_t> IVec2i;		//!< 32-bit signed   integer (int32_t)  2D integer vector class
typedef IVec3< int32_t> IVec3i;		//!< 32-bit signed   integer (int32_t)  3D integer vector class
typedef IVec4< int32_t> IVec4i;		//!< 32-bit signed   integer (int32_t)  4D integer vector class
typedef IVec2<uint32_t> IVec2ui;	//!< 32-bit unsigned integer (uint32_t) 2D integer vector class
typedef IVec3<uint32_t> IVec3ui;	//!< 32-bit unsigned integer (uint32_t) 3D integer vector class
typedef IVec4<uint32_t> IVec4ui;	//!< 32-bit unsigned integer (uint32_t) 4D integer vector class

typedef IVec2< int64_t> IVec2l;		//!< 64-bit signed   long    (int64_t)  2D integer vector class
typedef IVec3< int64_t> IVec3l;		//!< 64-bit signed   long    (int64_t)  3D integer vector class
typedef IVec4< int64_t> IVec4l;		//!< 64-bit signed   long    (int64_t)  4D integer vector class
typedef IVec2<uint64_t> IVec2ul;	//!< 64-bit unsigned long    (uint64_t) 2D integer vector class
typedef IVec3<uint64_t> IVec3ul;	//!< 64-bit unsigned long    (uint64_t) 3D integer vector class
typedef IVec4<uint64_t> IVec4ul;	//!< 64-bit unsigned long    (uint64_t) 4D integer vector class

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::IVec2b  cyIVec2b;		//!< 8-bit  signed   integer (int8_t)   2D integer vector class
typedef cy::IVec3b  cyIVec3b;		//!< 8-bit  signed   integer (int8_t)   3D integer vector class
typedef cy::IVec4b  cyIVec4b;		//!< 8-bit  signed   integer (int8_t)   4D integer vector class
typedef cy::IVec2ub cyIVec2ub;		//!< 8-bit  unsigned integer (uint8_t)  2D integer vector class
typedef cy::IVec3ub cyIVec3ub;		//!< 8-bit  unsigned integer (uint8_t)  3D integer vector class
typedef cy::IVec4ub cyIVec4ub;		//!< 8-bit  unsigned integer (uint8_t)  4D integer vector class

typedef cy::IVec2s  cyIVec2s;		//!< 16-bit signed   integer (int16_t)  2D integer vector class
typedef cy::IVec3s  cyIVec3s;		//!< 16-bit signed   integer (int16_t)  3D integer vector class
typedef cy::IVec4s  cyIVec4s;		//!< 16-bit signed   integer (int16_t)  4D integer vector class
typedef cy::IVec2us cyIVec2us;		//!< 16-bit unsigned integer (uint16_t) 2D integer vector class
typedef cy::IVec3us cyIVec3us;		//!< 16-bit unsigned integer (uint16_t) 3D integer vector class
typedef cy::IVec4us cyIVec4us;		//!< 16-bit unsigned integer (uint16_t) 4D integer vector class

typedef cy::IVec2i  cyIVec2i;		//!< 32-bit signed   integer (int32_t)  2D integer vector class
typedef cy::IVec3i  cyIVec3i;		//!< 32-bit signed   integer (int32_t)  3D integer vector class
typedef cy::IVec4i  cyIVec4i;		//!< 32-bit signed   integer (int32_t)  4D integer vector class
typedef cy::IVec2ui cyIVec2ui;		//!< 32-bit unsigned integer (uint32_t) 2D integer vector class
typedef cy::IVec3ui cyIVec3ui;		//!< 32-bit unsigned integer (uint32_t) 3D integer vector class
typedef cy::IVec4ui cyIVec4ui;		//!< 32-bit unsigned integer (uint32_t) 4D integer vector class

typedef cy::IVec2l  cyIVec2l;		//!< 64-bit signed   integer (int64_t)  2D integer vector class
typedef cy::IVec3l  cyIVec3l;		//!< 64-bit signed   integer (int64_t)  3D integer vector class
typedef cy::IVec4l  cyIVec4l;		//!< 64-bit signed   integer (int64_t)  4D integer vector class
typedef cy::IVec2ul cyIVec2ul;		//!< 64-bit unsigned integer (uint64_t) 2D integer vector class
typedef cy::IVec3ul cyIVec3ul;		//!< 64-bit unsigned integer (uint64_t) 3D integer vector class
typedef cy::IVec4ul cyIVec4ul;		//!< 64-bit unsigned integer (uint64_t) 4D integer vector class

//-------------------------------------------------------------------------------

#endif

