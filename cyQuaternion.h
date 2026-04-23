// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyQuaternion.h 
//! \author Cem Yuksel
//! 
//! \brief  Quaternion class
//! 
//! This file includes a quaternion class that can be used to rotate 3D vectors.
//! It works with Vec3, Matrix3, and Matrix4 classes.
//! Special thanks to Can Yuksel for his contributions in writing this class.
//!
//-------------------------------------------------------------------------------
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

#ifndef _CY_QUATERNION_H_INCLUDED_
#define _CY_QUATERNION_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyVector.h"
#include "cyMatrix.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

template <typename T> class UnitQuaternion;

//-------------------------------------------------------------------------------

//! Quaternion class
template <typename T>
class alignas(sizeof(T)*4) Quaternion
{
	CY_NODISCARD friend Quaternion operator * ( T          const &f, const Quaternion &q ) { return q * f; }
	CY_NODISCARD friend Matrix3<T> operator * ( Matrix3<T> const &m, const Quaternion &q ) { return m * q.GetMatrix3(); }
	CY_NODISCARD friend Quaternion Inverse( Quaternion const &q ) { return q.GetInverse(); }	//!< return the inverse of the quaternion

public:
	Vec3<T> v;	//!< vector part
	T       s;	//!< scalar part

	//!@name Constructors
	Quaternion() CY_CLASS_FUNCTION_DEFAULT
	Quaternion( T const &x, T const &y, T const &z, T const &_s ) : v(x,y,z),       s(_s)  {}
	Quaternion( Vec3<T> const &_v, T _s=0 )                       : v(_v),          s(_s)  {}
	explicit Quaternion( Vec4<T> const &q )                       : v(q.x,q.y,q.z), s(q.w) {}
	explicit Quaternion( UnitQuaternion<T> const &q );
	Quaternion( Vec3<T> const &from, Vec3<T> const &to ) { SetRotation(from,to); }

	//!@name Set & Get value functions
	void Zero       () { v.Zero(); s=T(1); }										//!< Sets the scalar part to one and vector part to zero, meaning no rotation.
	void Set        ( Vec3<T> const &_v, T _s ) { v=_v; s=_s; }						//!< Sets the scalar and vector components of the quaternion.
	void Set        ( T x, T y, T z, T _s ) { v.Set(x,y,z); s=_s; }					//!< Sets the scalar and vector components of the quaternion.
	bool IsFinite   () const { return v.IsFinite() && cy::IsFinite(s); }			//!< Returns true if all components are non-infinite numbers.
	bool IsUnit     () const { return std::abs(LengthSquared()-T(1)) < T(0.001); }	//!< Returns true if the length of the quaternion is close to 1.
	void SetRotation( Vec3<T> const &axis, T angle );								//!< Sets a rotation around the given axis.
	void SetRotation( Matrix3<T> const &rotation );									//!< Sets the quaternion using a rotation matrix. This method only works for rotation matrices.
	void SetRotation( Vec3<T> const &from, Vec3<T> const &to );						//!< Sets the quaternion as a rotation from one direction to another. The resulting quaternion is not normalized.
	void SetHalfRotation ()       { s += Length(); }								//!< Reduces the rotation of the quaternion by half

	CY_NODISCARD T       GetRotationAngle() const { return T(2)*std::acos(s/Length()); }	//!< Returns rotation angle in radians
	CY_NODISCARD Vec3<T> GetRotationAxis () const { return v.GetNormalized(); }				//!< Returns the normalized rotation axis

	//!@name Normalize, Length, Conjugate, and Inverse functions
	void                    Normalize    ()       { *this /= Length(); }				//!< Normalizes the quaternion, such that its length becomes 1.
	CY_NODISCARD Quaternion GetNormalized() const { return *this / Length(); }			//!< Returns a normalized copy of the quaternion.
	CY_NODISCARD T          LengthSquared() const { return s*s + v.LengthSquared(); }	//!< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	CY_NODISCARD T          Length       () const { return Sqrt(LengthSquared()); }		//!< Returns the length of the quaternion.
	void                    Conjugate    ()       { v = -v; }							//!< Converts the quaternion to its conjugate
	CY_NODISCARD Quaternion GetConjugate () const { return Quaternion(-v,s); }			//!< Returns the conjugate of the quaternion
	void                    Invert       ()       { T l2=LengthSquared(); v/=-l2; s/=l2; }					//!< Converts to quaternion to its inverse
	CY_NODISCARD Quaternion GetInverse   () const { T l2=LengthSquared(); return Quaternion(-v/l2,s/l2); }	//!< Returns the inverse of the quaternion

	//!@name Conversion functions
	CY_NODISCARD Matrix3<T> GetRotationMatrix() const;								//!< Returns a matrix representing the rotation of the quaternion.
	CY_NODISCARD Matrix3<T> GetMatrix3() const;										//!< Returns a matrix representing the rotation and scale of the quaternion.
	CY_NODISCARD Matrix4<T> GetMatrix4() const { return Matrix4<T>(GetMatrix3()); }	//!< Returns a matrix representing the rotation and scale of the quaternion.
	CY_NODISCARD explicit operator Matrix3<T>() const { return GetMatrix3(); }		//!< Returns a matrix representing the rotation ans scale of the quaternion.
	CY_NODISCARD explicit operator Matrix4<T>() const { return GetMatrix4(); }		//!< Returns a matrix representing the rotation ans scale of the quaternion.
	CY_NODISCARD explicit operator Vec4<T>   () const { return Vec4(v,s); }			//!< Returns a vector containing the values of the quaternion with the scalar in the w coordinate

	//!@name Unary operators
	CY_NODISCARD Quaternion operator - () const { return Quaternion(-v,-s); }

	//!@name Binary operators
	CY_NODISCARD Quaternion operator * ( Quaternion const &q ) const { return Quaternion( s*q.v + v*q.s + v.Cross(q.v), s*q.s - v.Dot(q.v) ); }
	CY_NODISCARD Quaternion operator + ( Quaternion const &q ) const { return Quaternion( v + q.v, s + q.s ); }
	CY_NODISCARD Quaternion operator - ( Quaternion const &q ) const { return Quaternion( v - q.v, s - q.s ); }
	CY_NODISCARD Quaternion operator * ( T          const &f ) const { return Quaternion( v*f, s*f ); }
	CY_NODISCARD Quaternion operator / ( T          const &f ) const { return Quaternion( v/f, s*f ); }
	CY_NODISCARD Vec3<T>    operator * ( Vec3<T>    const &p ) const { Vec3<T> sp = s*p; return v.Dot(p)*v + s*sp + 2*v.Cross(sp) + v.Cross(v.Cross(p)); }
	CY_NODISCARD Matrix3<T> operator * ( Matrix3<T> const &m ) const { return GetMatrix3()*m; }

	//!@name Assignment operators
	Quaternion& operator *= ( Quaternion const &q ) { Set( s*q.v + q.s*v + v.Cross(q.v), s*q.s - v.Dot(q.v) ); return *this; }
	Quaternion& operator += ( Quaternion const &q ) { v+=q.v; s+=q.s; return *this; }
	Quaternion& operator -= ( Quaternion const &q ) { v-=q.v; s-=q.s; return *this; }
	Quaternion& operator *= ( T          const &f ) { v*=f;   s*=f;   return *this; }
	Quaternion& operator /= ( T          const &f ) { v/=f;   s/=f;   return *this; }

	//!@name Test operators
	CY_NODISCARD int operator == ( Quaternion const &q ) const { return q.v == v && q.s == s; }
	CY_NODISCARD int operator != ( Quaternion const &q ) const { return q.v != v || q.s != s; }

	//!@name Dot product
	CY_NODISCARD T Dot        ( Quaternion const &p ) const { return v.Dot(p.v) + s*p.s; }	//!< Dot product
	CY_NODISCARD T operator % ( Quaternion const &p ) const { return Dot(p); }				//!< Dot product
};

//-------------------------------------------------------------------------------

//! A unit-length quaternion class.
//! The methods of this class guarantee that the quaternion remains unit length.
template <typename T>
class alignas(sizeof(T)*4) UnitQuaternion : private Quaternion<T>
{
	CY_NODISCARD friend Quaternion<T> operator * ( T          const &f, const UnitQuaternion &q ) { return q * f; }
	CY_NODISCARD friend Matrix3<T>    operator * ( Matrix3<T> const &m, const UnitQuaternion &q ) { return m * q.GetMatrix3(); }
	CY_NODISCARD friend UnitQuaternion Inverse( UnitQuaternion const &q ) { return q.GetInverse(); }	//!< return the inverse of the quaternion

public:

	//!@name Constructors
	UnitQuaternion() CY_CLASS_FUNCTION_DEFAULT
	UnitQuaternion( T const &x, T const &y, T const &z, T const &_s ) : Quaternion<T>(x,y,z,_s) { Quaternion<T>::Normalize(); }
	UnitQuaternion( Vec3<T> const &_v, T _s=0 )                       : Quaternion<T>(_v,   _s) { Quaternion<T>::Normalize(); }
	explicit UnitQuaternion( Vec4<T> const &q )                       : Quaternion<T>(q)        { Quaternion<T>::Normalize(); }
	explicit UnitQuaternion( Quaternion<T> const &q )                 : Quaternion<T>(q)        { Quaternion<T>::Normalize(); }
	UnitQuaternion( Vec3<T> const &from, Vec3<T> const &to ) { SetRotation(from,to); }

	//!@name Set & Get value functions
	using Quaternion<T>::Zero;
	using Quaternion<T>::IsFinite;
	using Quaternion<T>::GetRotationAngle;
	using Quaternion<T>::GetRotationAxis;
	void SetRotation( Vec3<T> const &axis, T angle );											//!< Sets a rotation around the given axis.
	void SetRotation( Matrix3<T> const &rotation ) { Quaternion<T>::SetRotation(rotation); }	//!< Sets the quaternion using a rotation matrix. This method only works for rotation matrices.
	void SetRotation( Vec3<T> const &from, Vec3<T> const &to );									//!< Sets the quaternion as a rotation from one direction to another. The input must be unit vectors.
	void SetHalfRotation() { Quaternion<T>::s += T(1); Quaternion<T>::Normalize(); }			//!< Reduces the rotation of the quaternion by half
	Vec3<T> const & GetVector() const { return Quaternion<T>::v; }
	T       const & GetScalar() const { return Quaternion<T>::s; }

	CY_NODISCARD T       GetRotationAngle() const { return T(2)*std::acos(Quaternion<T>::s); }	//!< Returns rotation angle in radians
	CY_NODISCARD Vec3<T> GetRotationAxis () const { return Quaternion<T>::v.GetNormalized(); }	//!< Returns the normalized rotation axis

	//!@name Normalize, Length, and Inverse functions
	void                        Normalize    ()       {}							//!< Does nothing, since the unit quaternion is already normalized.
	CY_NODISCARD UnitQuaternion GetNormalized() const { return *this; }				//!< Returns the quaternion itself, since it is already normalized.
	CY_NODISCARD T              LengthSquared() const { return T(1); }				//!< Returns 1, since the unit quaternion always has length 1.
	CY_NODISCARD T              Length       () const { return T(1); }				//!< Returns 1, since the unit quaternion always has length 1.
	void                        Conjugate    ()       { Quaternion::Conjugate(); }	//!< Converts the quaternion to its conjugate
	CY_NODISCARD UnitQuaternion GetConjugate () const { return UnitQuaternion(-Quaternion<T>::v,Quaternion<T>::s); }	//!< Returns the conjugate of the quaternion
	void                        Invert       ()       { Conjugate(); }				//!< Converts to quaternion to its inverse
	CY_NODISCARD UnitQuaternion GetInverse   () const { return GetConjugate(); }	//!< Returns the inverse of the quaternion

	//!@name Conversion functions
	CY_NODISCARD Matrix3<T> GetRotationMatrix() const { return Quaternion<T>::GetMatrix3(); }	//!< Returns a matrix representing the rotation of the quaternion.
	CY_NODISCARD Matrix3<T> GetMatrix3()        const { return Quaternion<T>::GetMatrix3(); }	//!< Returns a matrix representing the rotation of the quaternion.
	CY_NODISCARD Matrix4<T> GetMatrix4()        const { return Quaternion<T>::GetMatrix4(); }	//!< Returns a matrix representing the rotation of the quaternion.
	CY_NODISCARD explicit operator Matrix3<T>() const { return GetMatrix3(); }		//!< Returns a matrix representing the rotation of the quaternion.
	CY_NODISCARD explicit operator Matrix4<T>() const { return GetMatrix4(); }		//!< Returns a matrix representing the rotation of the quaternion.
	CY_NODISCARD operator Quaternion<T>() const { return Quaternion<T>(*this); }	//!< Returns a generic quaternion representation of the unit quaternion.
	using Quaternion<T>::operator Vec4<T>;

	//!@name Unary and Binary operators
	CY_NODISCARD UnitQuaternion operator * ( UnitQuaternion const &q ) const { Quaternion<T> m=Quaternion<T>::operator*(q); UnitQuaternion r; r.v=m.v; r.s=m.s; return r; }
	using Quaternion<T>::operator +;
	using Quaternion<T>::operator -;
	using Quaternion<T>::operator *;
	using Quaternion<T>::operator /;

	//! Rotates the given vector using the quaternion. If multiple vectors will be rotated, converting the quaternion to a matrix would be more efficient.
	CY_NODISCARD Vec3<T>    operator * ( Vec3<T>    const &p ) const { Vec3<T> vxp = Quaternion<T>::v.Cross(p); vxp += vxp; return p + vxp*Quaternion<T>::s + Quaternion<T>::v.Cross(vxp); }
	CY_NODISCARD Matrix3<T> operator * ( Matrix3<T> const &m ) const { return GetMatrix3()*m; }

	//!@name Assignment operators
	UnitQuaternion& operator *= ( UnitQuaternion const &q ) { *this = *this * q; return *this; }

	//!@name Test operators
	using Quaternion<T>::operator ==;
	using Quaternion<T>::operator !=;

	//!@name Dot product
	using Quaternion<T>::Dot;
	using Quaternion<T>::operator %;

	//!@name Slerp (Spherical Linear Interpolation)
	CY_NODISCARD static UnitQuaternion<T> Slerp( UnitQuaternion<T> const &q0, UnitQuaternion<T> const &q1, T t );
	CY_NODISCARD static UnitQuaternion<T> Slerp( UnitQuaternion<T> const &q0, UnitQuaternion<T> const &q1, UnitQuaternion<T> const &q2, Vec3<T> const &bary );
};

//-------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////
// Method Definitions
/////////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------------

template <typename T> Quaternion<T>::Quaternion( UnitQuaternion<T> const &q ) : v(q.GetVector()), s(q.GetScalar()) {}

//-------------------------------------------------------------------------------

template <typename T>
inline void Quaternion<T>::SetRotation( Matrix3<T> const &r )
{
	assert( std::abs(r.GetDeterminant()) - T(1) < T(0.0001) );	// The given matrix must be a rotation matrix.
	s   = Sqrt( Max<T>( 0, 1 + r[0] + r[4] + r[8] ) ) / T(2);
	v.x = Sqrt( Max<T>( 0, 1 + r[0] - r[4] - r[8] ) ) / T(2);
	v.y = Sqrt( Max<T>( 0, 1 - r[0] + r[4] - r[8] ) ) / T(2);
	v.z = Sqrt( Max<T>( 0, 1 - r[0] - r[4] + r[8] ) ) / T(2);
	v.x = CopySign( v.x, r(2,1) - r(1,2) );
	v.y = CopySign( v.y, r(0,2) - r(2,0) );
	v.z = CopySign( v.z, r(1,0) - r(0,1) );
}

//-------------------------------------------------------------------------------

template <typename T>
inline void Quaternion<T>::SetRotation( Vec3<T> const &axis, T angle )
{
	s = std::cos( angle*T(0.5) );
	v = std::sin( angle*T(0.5) )*axis.GetNormalized();
}

//-------------------------------------------------------------------------------

template <typename T>
inline void UnitQuaternion<T>::SetRotation( Vec3<T> const &axis, T angle )
{
	assert( axis.IsUnit() );
	Quaternion<T>::s = std::cos( angle*T(0.5) );
	Quaternion<T>::v = std::sin( angle*T(0.5) )*axis;
}

//-------------------------------------------------------------------------------

template <typename T>
inline void Quaternion<T>::SetRotation( Vec3<T> const &from, Vec3<T> const &to )
{
	T len = Sqrt( from.LengthSquared() * to.LengthSquared() );
	s = len + (from % to);
	v = from ^ to;
}

//-------------------------------------------------------------------------------

template <typename T>
inline void UnitQuaternion<T>::SetRotation( Vec3<T> const &from, Vec3<T> const &to )
{
	assert( from.IsUnit() );
	assert( to  .IsUnit() );
	Quaternion<T>::s = from % to;
	Quaternion<T>::v = from ^ to;
	SetHalfRotation();
}

//-------------------------------------------------------------------------------

template <typename T>
inline Matrix3<T> Quaternion<T>::GetRotationMatrix() const
{
	T const l2 = LengthSquared();
	T const x_l2 = v.x / l2;
	T const y_l2 = v.y / l2;
	T const z_l2 = v.z / l2;
	T const xx = v.x * x_l2;
	T const yy = v.y * y_l2;
	T const zz = v.z * z_l2;
	T const xz = v.x * z_l2;
	T const xy = v.x * y_l2;
	T const yz = v.y * z_l2;
	T const sx = s   * x_l2;
	T const sy = s   * y_l2;
	T const sz = s   * z_l2;

	Matrix3<T> m;
	m[0] = 1 - 2*(yy + zz);  m[1] =     2*(xy + sz);  m[2] =     2*(xz - sy);
	m[3] =     2*(xy - sz);  m[4] = 1 - 2*(xx + zz);  m[5] =     2*(yz + sx);
	m[6] =     2*(xz + sy);  m[7] =     2*(yz - sx);  m[8] = 1 - 2*(xx + yy);
	return m;
}

//-------------------------------------------------------------------------------

template <typename T>
inline Matrix3<T> Quaternion<T>::GetMatrix3() const
{
	T const xx = v.x * v.x;
	T const yy = v.y * v.y;
	T const zz = v.z * v.z;
	T const xz = v.x * v.z;
	T const xy = v.x * v.y;
	T const yz = v.y * v.z;
	T const sx = s   * v.x;
	T const sy = s   * v.y;
	T const sz = s   * v.z;

	Matrix3<T> m;
	m[0] = 1 - 2*(yy + zz);  m[1] =     2*(xy + sz);  m[2] =     2*(xz - sy);
	m[3] =     2*(xy - sz);  m[4] = 1 - 2*(xx + zz);  m[5] =     2*(yz + sx);
	m[6] =     2*(xz + sy);  m[7] =     2*(yz - sx);  m[8] = 1 - 2*(xx + yy);
	return m;
}

//-------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////
// Support Functions
/////////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------------

///! Interpolates two quaternions using parameter t in [0,1].
template <typename T>
inline UnitQuaternion<T> UnitQuaternion<T>::Slerp( UnitQuaternion<T> const &q0, UnitQuaternion<T> const &q1, T t )
{
	T cosTheta = q0 % q1;
	T flip = cosTheta < T(0) ? T(-1) : T(1);	// Slerp is not going to take the shortest path, we should effectively reverse one of the quaternions.
	T theta = std::acos( Min<T>(1,cosTheta*flip) );
	Quaternion<T> q;
	if ( std::abs(theta) < T(0.0001) ) {
		q = (q0 * (1-t)) + (q1 * t);
		q.Normalize();
	} else {
		T theta0 = theta * (1-t);
		T theta1 = theta * t;
		T sinTheta  = std::sin(theta);
		T sinTheta0 = std::sin(theta0);
		T s0 = sinTheta0 / sinTheta;
		T s1 = std::cos(theta0)*flip - cosTheta*sinTheta0 / sinTheta;
		q = (q0 * s0) + (q1 * s1);
	}
	UnitQuaternion<T> uq;
	uq.s = q.s;
	uq.v = q.v;
	return uq;
}

//-------------------------------------------------------------------------------

///! Interpolates three quaternions using the given barycentric coordinates
template <typename T>
inline UnitQuaternion<T> UnitQuaternion<T>::Slerp( UnitQuaternion<T> const &q0, UnitQuaternion<T> const &q1, UnitQuaternion<T> const &q2, Vec3<T> const &bc )
{
	T ab = bc.x + bc.y;
	if ( std::abs(ab) < T(0.000001) ) return q2;
	T t = bc.y / ab;
	UnitQuaternion<T> q01 = Slerp( q0, q1, t );
	return Slerp( q01, q2, T(1)-ab );
}

//-------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////
// Defined types
/////////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------------

typedef Quaternion    <float>  Quatf;	//!< Quaternion class with float  elements
typedef Quaternion    <double> Quatd;	//!< Quaternion class with double elements
typedef UnitQuaternion<float>  UQuatf;	//!< Unit quaternion class with float  elements
typedef UnitQuaternion<double> UQuatd;	//!< Unit quaternion class with double elements

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Quatf  cyQuatf;		//!< Quaternion class with float  elements
typedef cy::Quatd  cyQuatd;		//!< Quaternion class with double elements
typedef cy::UQuatf cyUQuatf;	//!< Unit quaternion class with float  elements
typedef cy::UQuatd cyUQuatd;	//!< Unit quaternion class with double elements

//-------------------------------------------------------------------------------

#endif
