// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyQuat.h 
//! \author Cem Yuksel
//! 
//! \brief  Quaternion class
//! 
//! This file includes a quaternion class that can be used to rotate 3D vectors.
//! It works with Point3, Matrix3, and Matrix4 classes.
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

#include "cyPoint.h"
#include "cyMatrix.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

//! Quaternion class

template <typename TYPE>
class Quat
{
	friend Quat operator*( TYPE f, const Quat &q ) { return Quat( q.s*f, q.v*f); }
public:
	TYPE         s;	//!< scaler part
	Point3<TYPE> v;	//!< vector part

	//!@name Constructors
	Quat() {}
	Quat( TYPE _s, const Point3<TYPE> &_v ) : s(_s), v(_v) {}
	Quat( TYPE _s, TYPE _x, TYPE _y, TYPE _z ) : s(_s), v(_x,_y,_z) {}
	Quat( const Quat &q ) : s(q.s), v(q.v) {}

	//!@name Set & Get value functions
	void         Zero() { s=0; v.Zero(); }
	void         Reset() { s=1; v.Zero(); }	//!< Sets the scalar part to 1 and vector part to zero vector.
	void         Set( TYPE _s, const Point3<TYPE> &_v ) { s=_s; v=_v; }
	void         Set( TYPE _s, TYPE _x, TYPE _y, TYPE _z ) { s=_s; v.Set(_x,_y,_z); }
	void         Set( TYPE *array ) { s=array[0]; v.Set(&array[1]); }
	void         SetRotation( TYPE angle, const Point3<TYPE> &axis ) { s=cos(angle*0.5f); v=(TYPE)sin(angle*0.5f)*(axis.GetNormalized()); }
	void         SetRotation( TYPE angle, TYPE axisX, TYPE axisY, TYPE axisZ ) { SetRotation( angle, Point3<TYPE>(axisX,axisY,axisZ) );	}
	void         GetValue( TYPE *array ) { array[0]=s; v.GetValue(&array[1]); }
	TYPE         GetRotationAngle() { return 2.0f*(TYPE)acos(s); }	//!< Returns rotation angle in radiants
	Point3<TYPE> GetRotationAxis () { return v.GetNormalized(); }

	//!@name Length and Normalize functions
	void         Normalize    ()       { TYPE f=1.0f/Length(); *this *= f; }
	Quat         GetNormalized() const { TYPE f=1.0f/Length(); return *this*f; }
	TYPE         LengthSquared() const { return s*s + v.LengthSquared(); }
	TYPE         Length       () const { return (TYPE) sqrt(LengthSquared()); }

	//!@name Matrix conversion functions
	void          ToMatrix3( Matrix3<TYPE> &m ) const { FillMatrix( m.data, &m.data[3], &m.data[6] ); }
	void          ToMatrix4( Matrix4<TYPE> &m ) const { FillMatrix( m.data, &m.data[4], &m.data[8] ); m.data[3]=m.data[7]=m.data[11]=m.data[12]=m.data[13]=m.data[14]=0; m.data[15]=1; }
	Matrix3<TYPE> ToMatrix3() const { Matrix3<TYPE> m; ToMatrix3(m); return m; }
	Matrix4<TYPE> ToMatrix4() const { Matrix4<TYPE> m; ToMatrix4(m); return m; }

	//!@name Unary operators
	Quat operator - () const { return Quat(-s,-v); }

	//!@name Binary operators
	Quat operator * ( const Quat &q ) const { return Quat( s*q.s - v.Dot(q.v), s*q.v + q.s*v + v.Cross(q.v)); }
	Quat operator + ( const Quat &q ) const { return Quat( s + q.s, v + q.v ); }
	Quat operator - ( const Quat &q ) const { return Quat( s - q.s, v - q.v ); }
	Quat operator * ( TYPE f ) const { return Quat( s*f, v*f); }

	//!@name Assignment operators
	Quat& operator *= ( const Quat &q ) { Set( s*q.s - v.Dot(q.v), s*q.v + q.s*v + v.Cross(q.v)); return *this; }
	Quat& operator += ( const Quat &q ) { s+=q.s; v+=q.v; return *this; }
	Quat& operator -= ( const Quat &q ) { s-=q.s; v-=q.v; return *this; }
	Quat& operator *= ( TYPE f ) { s*=f; v*=f; return *this; }

	//!@name Test operators
	int operator == ( const Quat& q ) const { return ( (q.s==s) && (q.v==v) ); }
	int operator != ( const Quat& q ) const { return ( (q.s!=s) || (q.v!=v) ); }

	//!@name Vector rotations
	//! Rotates the given vector using the quaternion.
	//! Note that the quaternion must be a unit quaternion.
	Point3<TYPE> GetRotatedVector( const Point3<TYPE> &p )
	{
		TYPE t2  =  s*v.x;
		TYPE t3  =  s*v.y;
		TYPE t4  =  s*v.z;
		TYPE t5  = -v.x*v.x;
		TYPE t6  =  v.x*v.y;
		TYPE t7  =  v.x*v.z;
		TYPE t8  = -v.y*v.y;
		TYPE t9  =  v.y*v.z;
		TYPE t10 = -v.z*v.z;
		Point3<TYPE> pnew;
		pnew.x = 2*( (t8 + t10)*p.x + (t6 -  t4)*p.y + (t3 + t7)*p.z ) + p.x;
		pnew.y = 2*( (t4 +  t6)*p.x + (t5 + t10)*p.y + (t9 - t2)*p.z ) + p.y;
		pnew.z = 2*( (t7 -  t3)*p.x + (t2 +  t9)*p.y + (t5 + t8)*p.z ) + p.z;
		return pnew;
	}
	void RotateVector( Point3<TYPE> &p ) { p = GetRotatedVector(p); }

private:
	//! \internal
	void FillMatrix( TYPE *col1, TYPE *col2, TYPE *col3 ) const
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
