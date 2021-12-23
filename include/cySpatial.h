// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cySpatial.h 
//! \author Cem Yuksel
//! 
//! \brief  Spatial vector algebra classes
//! 
//! This file includes spatial vector algebra classes intended for the
//! implementation of Featherstone's articulated rigid body dynamics method.
//! SpatialVector6 class is both for spatial motion vectors and spatial
//! force vectors, SpatialTrans6 is a spatial matrix class for coordinate
//! transformations only, and SpatialMatrix6 is the general spatial
//! matrix class.
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

#ifndef _CY_SPATIAL_H_INCLUDED_
#define _CY_SPATIAL_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyVector.h"
#include "cyMatrix.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

//! 6D spatial vector (for 3D).
//!
//! This class is both for spatial motion vectors and spatial force vectors.

template <typename T>
class SpatialVector6
{
public:

	//!@name Components
	Vec3<T> a, b;

	//!@name Constructors
	SpatialVector6() CY_CLASS_FUNCTION_DEFAULT
	SpatialVector6( Vec3<T> const &p1, Vec3<T> const &p2 ) { Set(p1,p2); }
	SpatialVector6( T a1, T a2, T a3, T b1, T b2, T b3 ) { Set(a1,a2,a3,b1,b2,b3); }
	SpatialVector6( SpatialVector6 const &v ) { a=v.a; b=v.b; }

	//!@name Initialization methods
	void Set( Vec3<T> const &p1, Vec3<T> const &p2 ) { a = p1; b = p2; }
	void Set( T a1, T a2, T a3, T b1, T b2, T b3 ) { a.Set(a1,a2,a3); b.Set(b1,b2,b3); }
	void Zero() { a.Zero(); b.Zero(); }

	//!@name Transpose methods
	void SetTranspose() { Vec3<T> p=a; a=b; b=p; }
	SpatialVector6 Transpose() const { return SpatialVector6( b, a ); }

	//!@name Unary operators
	SpatialVector6 operator-() const { return SpatialVector6(-a,-b); } 

	//!@name Binary operators
	SpatialVector6	operator + ( SpatialVector6 const &s ) const { return SpatialVector6(a+s.a, b+s.b); }
	SpatialVector6	operator - ( SpatialVector6 const &s ) const { return SpatialVector6(a-s.a, b-s.b); }
	SpatialVector6	operator * ( T              const &t ) const { return SpatialVector6( a*t, b*t ); }

	//! Scalar product of two vectors.
	//! Note that one of the vectors should be motion vector ant the other should be a force vector.
	//! Otherwise, scalar product is not defined in spatial vector algebra.
	T operator * ( SpatialVector6 const &s ) const { return a.Dot(s.a) + b.Dot(s.b); }

	//!@name Assignment operators
	void operator =  ( SpatialVector6 const &v ) { a=v.a; b=v.b; }
	void operator += ( SpatialVector6 const &s ) { a+=s.a; b+=s.b; }
	void operator -= ( SpatialVector6 const &s ) { a-=s.a; b-=s.b; }
	void operator *= ( T              const &t ) { a*=t; b*=t; }

};


//-------------------------------------------------------------------------------

//! 6D spatial matrix for coordinate transforms.
//!
//! This is a special case for SpatialMatrix6 class,
//! where the matrix represents a coordinate transformation.
//! In this case, instead of keeping a full 6x6 matrix values,
//! we can keep a 3x3 matrix for rotation, and a 3D point
//! for translation. This compact representation simplifies
//! some computations, therefore you should use this class
//! instead of SpatialMatrix6 whenever you represent
//! a coordinate transformation. However, for general matrix
//! operations, you have to use SpatialMatrix6.
//!

template <typename T>
class SpatialTrans6
{
public:

	// |    R     0 |
	// | -r x R   R |

	Matrix3<T> R;	//!< Rotation matrix
	Vec3<T>    r;	//!< Transformation

	//!@name Constructors
	SpatialTrans6() CY_CLASS_FUNCTION_DEFAULT
	SpatialTrans6( SpatialTrans6 const &mat ) { R=mat.R; r=mat.r; }
	SpatialTrans6( Matrix3<T> const &_R, Vec3<T> const &_r ) { Set(_R,_r); }

	//!@name Initialization methods
	void Set( Matrix3<T> const &_R, Vec3<T> const &_r ) { R=_R; r=_r; }
	void SetIdentity() { R.SetIdentity(); r.Zero(); }

	//!@name Unary operators
	SpatialTrans6 operator - () const { return SpatialTrans6( -R, -r ); }


	//!@name Binary operators
	SpatialVector6<T> operator * ( SpatialVector6<T> const &p ) const { Vec3<T> Ra = R*p.a; return SpatialVector6<T>( Ra, Matrix3<T>::Scale(-r)*Ra + R*p.b ); }

	SpatialTrans6 operator * ( SpatialTrans6 const &mat ) const { return SpatialTrans6( R*mat.R, r + R*mat.r ); }
	SpatialTrans6 operator + ( SpatialTrans6 const &mat ) const { return SpatialTrans6( R + mat.R, r + mat.r ); }
	SpatialTrans6 operator - ( SpatialTrans6 const &mat ) const { return SpatialTrans6( R - mat.R, r - mat.r ); }

	SpatialTrans6 operator * ( T t ) const { return SpatialTrans6(R*t,r*t); }
	SpatialTrans6 operator / ( T t ) const { T d=1.0f/t; return *this * d; }


	//!@name Assignment operators
	void operator *= ( SpatialTrans6 const &mat ) { *this = *this * mat; }
	void operator += ( SpatialTrans6 const &mat ) { R+=mat.R; r+=mat.r; }
	void operator -= ( SpatialTrans6 const &mat ) { R-=mat.R; r-=mat.r; }
	void operator *= ( T             const &t   ) { *this = *this * t; }

};

//-------------------------------------------------------------------------------

//! 6D spatial matrix.
//!
//! This is the general class for 6D spatial matrices.
//! For representing coordinate transformation matrices
//! use SpatialTrans6 instead, since it is more efficient.
//! However, SpatialTrans6 cannot be used for general
//! matrix operations that do not correspond to a
//! coordinate transformation.

template <typename T>
class SpatialMatrix6
{
public:

	// | m[0]  m[2] |
	// | m[1]  m[3] |

	Matrix3<T> m[4];	//!< Matrix data in column major order


	//!@name Constructors
	SpatialMatrix6() CY_CLASS_FUNCTION_DEFAULT
	SpatialMatrix6( SpatialMatrix6 const &mat ) { m[0]=mat.m[0]; m[1]=mat.m[1]; m[2]=mat.m[2]; m[3]=mat.m[3]; }
	explicit SpatialMatrix6( Matrix3<T> const &_R, Vec3<T> const &_r ) { Set(_R,_r); }
	explicit SpatialMatrix6( Matrix3<T> const &m11, Matrix3<T> const &m21, Matrix3<T> const &m12, Matrix3<T> const &m22 ) { Set(m11,m21,m12,m22); }
	explicit SpatialMatrix6( SpatialTrans6<T> const &tm ) { Set(tm); }


	//!@name Initialization methods
	void Set( Matrix3<T> const &_R, Vec3<T> const &_r ) { m[0]=m[3]=_R; m[1]=Matrix3<T>::Scale(-_r)*_R; m[2].Zero(); }
	void Set( Matrix3<T> const &m11, Matrix3<T> const &m21, Matrix3<T> const &m12, Matrix3<T> const &m22 ) { m[0]=m11; m[1]=m21; m[2]=m12; m[3]=m22; }
	void Set( SpatialTrans6<T> const &tm ) { m[0]=m[3]=tm.R; m[1]=Matrix3<T>::Scale(-tm.r)*tm.R; m[2].Zero(); }

	//! Sets the matrix as the outer product of two vectors.
	void SetTensorProduct( SpatialVector6<T> const &p1, SpatialVector6<T> const &p2 )
	{
			SetMatrix( m[0], p1.a, p2.a );
			SetMatrix( m[1], p1.b, p2.a );
			SetMatrix( m[2], p1.a, p2.b );
			SetMatrix( m[3], p1.b, p2.b );
	}

	void SetIdentity() { m[0].SetIdentity(); m[1].Zero(); m[2].Zero(); m[3].SetIdentity(); }
	void Zero() { m[0].Zero(); m[1].Zero(); m[2].Zero(); m[3].Zero(); }


	//!@name Unary operators
	SpatialMatrix6 operator - () const { return SpatialMatrix6( -m[0], -m[1], -m[2], -m[3] ); }


	//!@name Unary operators

	SpatialVector6<T> operator * ( SpatialVector6<T> const &p ) const { return SpatialVector6<T>( m[0]*p.a + m[2]*p.b, m[1]*p.a + m[3]*p.b ); }

	SpatialMatrix6 operator * ( SpatialMatrix6 const &mat ) const { return SpatialMatrix6( m[0]*mat.m[0]+m[2]*mat.m[1], m[1]*mat.m[0]+m[3]*mat.m[1], m[0]*mat.m[2]+m[2]*mat.m[3], m[1]*mat.m[2]+m[3]*mat.m[3] ); }
	SpatialMatrix6 operator + ( SpatialMatrix6 const &mat ) const { return SpatialMatrix6( m[0]+mat.m[0], m[1]+mat.m[1], m[2]+mat.m[2], m[3]+mat.m[3] ); }
	SpatialMatrix6 operator - ( SpatialMatrix6 const &mat ) const { return SpatialMatrix6( m[0]-mat.m[0], m[1]-mat.m[1], m[2]-mat.m[2], m[3]-mat.m[3] ); }

	SpatialMatrix6 operator * ( T t ) const { return SpatialMatrix6(m[0]*t,m[1]*t,m[2]*t,m[3]*t); }
	SpatialMatrix6 operator / ( T t ) const { T d=1.0f/t; return *this * d; }


	//!@name Assignment operators
	void	operator *= ( SpatialMatrix6 const &mat ) { *this = *this * mat; }
	void	operator += ( SpatialMatrix6 const &mat ) { m[0]+=mat.m[0]; m[1]+=mat.m[1]; m[2]+=mat.m[2]; m[3]+=mat.m[3]; }
	void	operator -= ( SpatialMatrix6 const &mat ) { m[0]-=mat.m[0]; m[1]-=mat.m[1]; m[2]-=mat.m[2]; m[3]-=mat.m[3]; }
	void	operator *= ( T t ) { *this = *this * t; }


protected:

	//! \internal

	// Sets the given matrix as the outer product of the given two vectors.
	void SetMatrix( Matrix3<T> &m, Vec3<T> const &p1, Vec3<T> const &p2 )
	{ 
		T val[] = {p1.x * p2.x, p1.y * p2.x, p1.z * p2.x,   p1.x * p2.y, p1.y * p2.y, p1.z * p2.y,   p1.x * p2.z, p1.y * p2.z, p1.z * p2.z};
		m.Set( val ); 
	}

};

//-------------------------------------------------------------------------------

template<typename T> inline SpatialMatrix6<T> operator & ( SpatialVector6<T> const &v0, SpatialVector6<T> const &v1 ) { Matrix2<T> buffer; buffer.SetTensorProduct(v0,v1); return buffer; }		//!< tensor product (outer product) of two vectors

//-------------------------------------------------------------------------------

typedef SpatialVector6<float>  SpatialVector6f;	//!< 6D spatial vector (for 3D) with float type elements
typedef SpatialTrans6 <float>  SpatialTrans6f;	//!< 6D spatial matrix for coordinate transforms with float type elements
typedef SpatialMatrix6<float>  SpatialMatrix6f;	//!< 6D spatial matrix with float type elements

typedef SpatialVector6<double> SpatialVector6d;	//!< 6D spatial vector (for 3D) with double type elements
typedef SpatialTrans6 <double> SpatialTrans6d;	//!< 6D spatial matrix for coordinate transforms with double type elements
typedef SpatialMatrix6<double> SpatialMatrix6d;	//!< 6D spatial matrix with double type elements

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::SpatialVector6f cySpatialVector6f;	//!< 6D spatial vector (for 3D) with float type elements
typedef cy::SpatialTrans6f  cySpatialTrans6f;	//!< 6D spatial matrix for coordinate transforms with float type elements
typedef cy::SpatialMatrix6f cySpatialMatrix6f;	//!< 6D spatial matrix with float type elements

typedef cy::SpatialVector6d cySpatialVector6d;	//!< 6D spatial vector (for 3D) with double type elements
typedef cy::SpatialTrans6d  cySpatialTrans6d;	//!< 6D spatial matrix for coordinate transforms with double type elements
typedef cy::SpatialMatrix6d cySpatialMatrix6d;	//!< 6D spatial matrix with double type elements

//-------------------------------------------------------------------------------

#endif
