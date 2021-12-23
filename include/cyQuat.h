// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyQuat.h 
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

#ifndef _CY_QUAT_H_INCLUDED_
#define _CY_QUAT_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyVector.h"
#include "cyMatrix.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

//! Quaternion class

template <typename T>
class Quat
{
	friend Quat operator*( T f, Quat const &q ) { return Quat( q.s*f, q.v*f); }
public:
	T       s;	//!< scaler part
	Vec3<T> v;	//!< vector part

	//!@name Constructors
	Quat() CY_CLASS_FUNCTION_DEFAULT
	Quat( T _s, Vec3<T> const &_v ) : s(_s), v(_v) {}
	Quat( T _s, T _x, T _y, T _z ) : s(_s), v(_x,_y,_z) {}
	Quat( Quat const &q ) : s(q.s), v(q.v) {}

	//!@name Set & Get value functions
	void    Zero() { s=0; v.Zero(); }
	void    Reset() { s=1; v.Zero(); }	//!< Sets the scalar part to 1 and vector part to zero vector.
	void    Set( T _s, Vec3<T> const &_v ) { s=_s; v=_v; }
	void    Set( T _s, T _x, T _y, T _z ) { s=_s; v.Set(_x,_y,_z); }
	void    Set( T *array ) { s=array[0]; v.Set(&array[1]); }
	void    SetRotation( T angle, Vec3<T> const &axis ) { s=cos(angle*0.5f); v=(T)sin(angle*0.5f)*(axis.GetNormalized()); }
	void    SetRotation( T angle, T axisX, T axisY, T axisZ ) { SetRotation( angle, Vec3<T>(axisX,axisY,axisZ) );	}
	void    Get( T *array ) { array[0]=s; v.Get(&array[1]); }
	T       GetRotationAngle() { return 2.0f*(T)acos(s); }	//!< Returns rotation angle in radiants
	Vec3<T> GetRotationAxis () { return v.GetNormalized(); }

	//!@name Length and Normalize functions
	void Normalize    ()       { T f=1.0f/Length(); *this *= f; }
	Quat GetNormalized() const { T f=1.0f/Length(); return *this*f; }
	T    LengthSquared() const { return s*s + v.LengthSquared(); }
	T    Length       () const { return (T) sqrt(LengthSquared()); }

	//!@name Matrix conversion functions
	void       ToMatrix3( Matrix3<T> &m ) const { FillMatrix( m.cell, &m.cell[3], &m.cell[6] ); }
	void       ToMatrix4( Matrix4<T> &m ) const { FillMatrix( m.cell, &m.cell[4], &m.cell[8] ); m.cell[3]=m.cell[7]=m.cell[11]=m.cell[12]=m.cell[13]=m.cell[14]=0; m.cell[15]=1; }
	Matrix3<T> ToMatrix3() const { Matrix3<T> m; ToMatrix3(m); return m; }
	Matrix4<T> ToMatrix4() const { Matrix4<T> m; ToMatrix4(m); return m; }

	//!@name Unary operators
	Quat operator - () const { return Quat(-s,-v); }

	//!@name Binary operators
	Quat operator * ( Quat const &q ) const { return Quat( s*q.s - v.Dot(q.v), s*q.v + q.s*v + v.Cross(q.v)); }
	Quat operator + ( Quat const &q ) const { return Quat( s + q.s, v + q.v ); }
	Quat operator - ( Quat const &q ) const { return Quat( s - q.s, v - q.v ); }
	Quat operator * ( T    const &f ) const { return Quat( s*f, v*f); }

	//!@name Assignment operators
	Quat& operator *= ( Quat const &q ) { Set( s*q.s - v.Dot(q.v), s*q.v + q.s*v + v.Cross(q.v)); return *this; }
	Quat& operator += ( Quat const &q ) { s+=q.s; v+=q.v; return *this; }
	Quat& operator -= ( Quat const &q ) { s-=q.s; v-=q.v; return *this; }
	Quat& operator *= ( T    const &f ) { s*=f; v*=f; return *this; }

	//!@name Test operators
	int operator == ( Quat const &q ) const { return ( (q.s==s) && (q.v==v) ); }
	int operator != ( Quat const &q ) const { return ( (q.s!=s) || (q.v!=v) ); }

	//!@name Vector rotations
	//! Rotates the given vector using the quaternion.
	//! Note that the quaternion must be a unit quaternion.
	Vec3<T> GetRotatedVector( Vec3<T> const &p )
	{
		T t2  =  s*v.x;
		T t3  =  s*v.y;
		T t4  =  s*v.z;
		T t5  = -v.x*v.x;
		T t6  =  v.x*v.y;
		T t7  =  v.x*v.z;
		T t8  = -v.y*v.y;
		T t9  =  v.y*v.z;
		T t10 = -v.z*v.z;
		Vec3<T> pnew;
		pnew.x = 2*( (t8 + t10)*p.x + (t6 -  t4)*p.y + (t3 + t7)*p.z ) + p.x;
		pnew.y = 2*( (t4 +  t6)*p.x + (t5 + t10)*p.y + (t9 - t2)*p.z ) + p.y;
		pnew.z = 2*( (t7 -  t3)*p.x + (t2 +  t9)*p.y + (t5 + t8)*p.z ) + p.z;
		return pnew;
	}
	void RotateVector( Vec3<T> &p ) { p = GetRotatedVector(p); }

private:
	//! \internal
	void FillMatrix( T *col1, T *col2, T *col3 ) const
	{
		col1[0] = s*s + v.x*v.x - v.y*v.y - v.z*v.z;
		col1[1] = 2.0f*v.x*v.y + 2.0f*s*v.z;
		col1[2] = 2.0f*v.x*v.z - 2.0f*s*v.y;

		col2[0] = 2.0f*v.x*v.y - 2.0f*s*v.z;
		col2[1] = s*s - v.x*v.x + v.y*v.y - v.z*v.z;
		col2[2] = 2.0f*v.y*v.z + 2.0f*s*v.x;

		col3[0] = 2.0f*v.x*v.z + 2.0f*s*v.y;
		col3[1] = 2.0f*v.y*v.z - 2.0f*s*v.x;
		col3[2] = s*s - v.x*v.x - v.y*v.y + v.z*v.z;
	}
};

//-------------------------------------------------------------------------------

typedef Quat<float>  Quatf;	//!< Quaternion class with float  elements
typedef Quat<double> Quatd;	//!< Quaternion class with double elements

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Quatf cyQuatf;	//!< Quaternion class with float  elements
typedef cy::Quatd cyQuatd;	//!< Quaternion class with double elements

//-------------------------------------------------------------------------------

#endif
