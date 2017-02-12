// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyMatrix.h 
//! \author Cem Yuksel
//! 
//! \brief  2x2, 3x3, 3x4, and 4x4 matrix classes
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

#ifndef _CY_MATRIX_H_INCLUDED_
#define _CY_MATRIX_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyPoint.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

// Forward declarations
//!	\cond HIDDEN_SYMBOLS
template <typename TYPE> class Matrix3;
template <typename TYPE> class Matrix34;
template <typename TYPE> class Matrix4;
//! \endcond

//-------------------------------------------------------------------------------

//! 2x2 matrix class.
//!
//! Its data stores 4-value array of column-major matrix elements.
//! You can use Matrix2 with Point2<TYPE> to transform 2D points.

template <typename TYPE>
class Matrix2
{
	friend Matrix2 operator * ( const TYPE value, const Matrix2 &right ) { Matrix2 r; for (int i=0; i<4; i++) r.data[i] = value * right.data[i]; return r; }	//!< multiply matrix by a value
	friend Matrix2 Inverse( const Matrix2 &m ) { return m.GetInverse(); }	//!< return the inverse of the matrix

public:

	//! Elements of the matrix are column-major: \n
	//! | 0  2 | \n
	//! | 1  3 | \n
	TYPE data[4];

	//////////////////////////////////////////////////////////////////////////
	//!@name Constructors

	Matrix2() {}																//!< Default constructor
	Matrix2( const Matrix2 &matrix ) { CY_MEMCOPY(TYPE,data,matrix.data,4); }	//!< Copy constructor
	template <typename T> explicit Matrix2<TYPE>( const Matrix2<T> &matrix ) { CY_MEMCONVERT(TYPE,data,matrix.data,4); }	//!< Copy constructor for different types
	explicit Matrix2( const TYPE *values ) { Set(values); }									//!< Initialize the matrix using an array of 4 values
	explicit Matrix2( const TYPE &v ) { SetScaledIdentity(v); }								//!< Initialize the matrix as identity scaled by v
	explicit Matrix2( const Point2<TYPE> &x, const Point2<TYPE> &y ) { Set(x,y); }			//!< Initialize the matrix using two vectors as columns
	explicit Matrix2( const Matrix3<TYPE>  &m );
	explicit Matrix2( const Matrix34<TYPE> &m );
	explicit Matrix2( const Matrix4<TYPE>  &m );


	//////////////////////////////////////////////////////////////////////////
	//!@name Set & Get Methods

	//! Set all the values as zero
	void Zero() { CY_MEMCLEAR(TYPE,data,4); }
	//! Returns true if the matrix is exactly zero
	bool IsZero() const { for ( int i=0; i<4; i++ ) if ( data[i] != 0 ) return false; return true; }
	//! Copies the matrix data to the given values array of size 4
	void Get( TYPE *values ) { CY_MEMCOPY(TYPE,values,data,4); } 
	//! Set Matrix using an array of 4 values
	void Set( const TYPE *values ) { CY_MEMCOPY(TYPE,data,values,4); } 
	//! Set Matrix using two vectors as columns
	void Set( const Point2<TYPE> &x, const Point2<TYPE> &y ) { x.Get(data); y.Get(data+2); }
	//! Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	//! Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity(TYPE v) { SetScale(v); }
	//! Sets the matrix as the tensor product (outer product) of two vectors
	void SetTensorProduct( const Point2<TYPE> &v0, const Point2<TYPE> &v1 )
	{
		for ( int i=0; i<2; i++ ) data[  i] = v0[i] * v1.x;
		for ( int i=0; i<2; i++ ) data[2+i] = v0[i] * v1.y;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Affine transformations

	//! Sets a uniform scale matrix
	void SetScale( const TYPE &uniformScale ) { SetScale(uniformScale,uniformScale); }
	//! Sets a scale matrix
	void SetScale( const TYPE &scaleX, const TYPE &scaleY ) { data[0]=scaleX; data[1]=0; data[2]=0; data[3]=scaleY;}
	//! Sets a scale matrix
	void SetScale( const Point2<TYPE> &scale ) { SetScale(scale.x,scale.y); }
	//! Removes the scale component of the matrix
	void SetNoScale() { Point2<TYPE> *p = (Point2<TYPE>*)data; p[0].Normalize(); p[1].Normalize(); }
	//! Set a rotation matrix by angle
	void SetRotation( TYPE angle ) { SetRotation( cySin(angle), cyCos(angle) ); }
	//! Set a rotation matrix by cos and sin of angle
	void SetRotation( TYPE sinAngle, TYPE cosAngle ) { data[0]=cosAngle; data[1]=-sinAngle; data[2]=sinAngle; data[3]=cosAngle; }


	//////////////////////////////////////////////////////////////////////////
	//!@name Set Row, Column, or Diagonal

	void SetRow( int row, TYPE x, TYPE y ) { data[row]=x; data[row+2]=y; }						//!< Sets a row of the matrix
	void SetColumn( int column, TYPE x, TYPE y ) { data[2*column]=x; data[2*column+1]=y; }		//!< Sets a column of the matrix
	void SetDiagonal( const TYPE &xx, const TYPE &yy ) { data[0]=xx; data[3]=yy; }				//!< Sets the diagonal values of the matrix
	void SetDiagonal( const Point2<TYPE> &p ) { SetDiagonal( p.x, p.y ); }						//!< Sets the diagonal values of the matrix
	void SetDiagonal( const TYPE *values ) { SetDiagonal(values[0],values[1]); }				//!< Sets the diagonal values of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Row, Column, or Diagonal

	Point2<TYPE>  GetRow   ( int row )                  const { return Point2<TYPE>( data[row], data[row+2] ); }	//!< Returns a row of the matrix
	void          GetRow   ( int row, Point2<TYPE> &p ) const { p.Set( data[row], data[row+1] ); }					//!< Returns a row of the matrix
	void          GetRow   ( int row, TYPE *values )    const { values[0]=data[row]; values[1]=data[row+2]; }		//!< Returns a row of the matrix
	Point2<TYPE>  GetColumn( int col )                  const { return Point2<TYPE>( &data[col*2] ); }				//!< Returns a column of the matrix
	void          GetColumn( int col, Point2<TYPE> &p ) const { p.Set( &data[col*2] ); }							//!< Returns a column of the matrix
	void          GetColumn( int col, TYPE *values )    const { values[0]=data[col*2]; values[1]=data[col*2+1]; }	//!< Returns a column of the matrix
	Point2<TYPE>  GetDiagonal()                         const { return Point2<TYPE>( data[0], data[3] ); }			//!< Returns the diagonal of the matrix
	void          GetDiagonal( Point2<TYPE> &p )        const { p.Set( data[0], data[3] ); }						//!< Returns the diagonal of the matrix
	void	      GetDiagonal( TYPE *values )           const { values[0]=data[0]; values[1]=data[3]; }				//!< Returns the diagonal of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Comparison Operators

	bool operator == ( const Matrix2 &right ) const { for ( int i=0; i<4; i++ ) if ( data[i] != right.data[i] ) return false; return true; } //!< compare equal
	bool operator != ( const Matrix2 &right ) const { for ( int i=0; i<4; i++ ) if ( data[i] != right.data[i] ) return true; return false; } //!< compare not equal


	//////////////////////////////////////////////////////////////////////////
	//!@name Access Operators

	TYPE&       operator () ( int row, int column )       { return data[ column * 2 + row ]; }	//!< subscript operator
	const TYPE& operator () ( int row, int column ) const { return data[ column * 2 + row ]; }	//!< constant subscript operator
	TYPE&       operator [] ( int i )       { return data[i]; }	//!< subscript operator
	const TYPE& operator [] ( int i ) const { return data[i]; }	//!< constant subscript operator
	
	//////////////////////////////////////////////////////////////////////////
	//!@name Unary and Binary Operators

	// Unary operators
	Matrix2 operator - () const { Matrix2 r; for (int i=0; i<4; i++) r.data[i]=- data[i]; return r; }	//!< negative matrix

	// Binary operators
	Matrix2 operator * ( const TYPE    &value ) const { Matrix2 r; for (int i=0; i<4; i++) r.data[i] = data[i] * value;         return r; }	//!< multiply matrix by a value
	Matrix2 operator / ( const TYPE    &value ) const { Matrix2 r; for (int i=0; i<4; i++) r.data[i] = data[i] / value;         return r; }	//!< divide matrix by a value;
	Matrix2 operator + ( const Matrix2 &right ) const { Matrix2 r; for (int i=0; i<4; i++) r.data[i] = data[i] + right.data[i]; return r; }	//!< add two Matrices
	Matrix2 operator - ( const Matrix2 &right ) const { Matrix2 r; for (int i=0; i<4; i++) r.data[i] = data[i] - right.data[i]; return r; }	//!< subtract one Matrix2 from another
	Matrix2 operator * ( const Matrix2 &right ) const	//!< multiply a matrix with another
	{
		Matrix2 r;
		r[0] = data[0] * right.data[0] + data[2] * right.data[1];
		r[1] = data[1] * right.data[0] + data[3] * right.data[1];
		r[2] = data[0] * right.data[2] + data[2] * right.data[3];
		r[3] = data[1] * right.data[2] + data[3] * right.data[3];
		return r;
	}
	Point2<TYPE> operator * ( const Point2<TYPE> &p ) const { return Point2<TYPE>( p.x*data[0] + p.y*data[2], p.x*data[1] + p.y*data[3] ); }


	//////////////////////////////////////////////////////////////////////////
	//!@name Assignment Operators

	const Matrix2& operator  = ( const Matrix2 &right ) { CY_MEMCOPY(TYPE,data,right.data,4); return *this; }	
	const Matrix2& operator += ( const Matrix2 &right ) { for (int i=0; i<4; i++) data[i] += right.data[i]; return *this; }	//!< add two Matrices modify this
	const Matrix2& operator -= ( const Matrix2 &right ) { for (int i=0; i<4; i++) data[i] -= right.data[i]; return *this; }	//!< subtract one Matrix2 from another matrix and modify this matrix
	const Matrix2& operator *= ( const Matrix2 &right ) { *this = operator*(right); return *this; }							//!< multiply a matrix with another matrix and modify this matrix
	const Matrix2& operator *= ( const TYPE value )     { for (int i=0; i<4; i++) data[i] *= value;         return *this; }	//!< multiply a matrix with a value modify this matrix
	const Matrix2& operator /= ( const TYPE value )     { for (int i=0; i<4; i++) data[i] /= value;         return *this; }	//!< divide the matrix by a value modify the this matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Other Methods

	void Transpose() { TYPE tmp=data[0]; data[0]=data[3]; data[3]=tmp; }	//!< Transpose this matrix
	void GetTranspose( Matrix2 &m ) const									//!< return Transpose of this matrix
	{
		m.data[0] = data[0];   m.data[1] = data[2];
		m.data[2] = data[1];   m.data[3] = data[3];
	}
	Matrix2 GetTranspose() const { Matrix2 t; GetTranspose(t); return t; }	//!< return Transpose of this matrix

	//! Multiply the give vector with the transpose of the matrix
	Point2<TYPE> TransposeMult( const Point2<TYPE> &p ) const { return Point2<TYPE>( p.x*data[0] + p.y*data[1], p.x*data[2] + p.y*data[3] ); }

	TYPE GetDeterminant() const { return data[0]*data[3]-data[2]*data[1]; }	//!< Get the determinant of this matrix

	void Invert()					//!< Invert this matrix
	{
		TYPE det = GetDeterminant();
		TYPE d0 =  data[0] / det;
		data[0] =  data[3] / det;
		data[1] = -data[1] / det;
		data[2] = -data[2] / det;
		data[3] =  d0;
	}
	void GetInverse( Matrix2 &inverse ) const { inverse=*this; inverse.Invert(); }	//!< Get the inverse of this matrix
	Matrix2 GetInverse() const { Matrix2 inv(*this); inv.Invert(); return inv; }	//!< Get the inverse of this matrix

	//! Orthogonalizes the matrix and removes the scale component, preserving the x direction
	void OrthogonalizeX()
	{
		Point2<TYPE> *p = (Point2<TYPE>*)data;
		p[0].Normalize();
		p[1] -= p[0] * (p[1]%p[0]);
		p[1].Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the y direction
	void OrthogonalizeY()
	{
		Point2<TYPE> *p = (Point2<TYPE>*)data;
		p[1].Normalize();
		p[0] -= p[1] * (p[0]%p[1]);
		p[0].Normalize();
	}

	//! Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const { return cyAbs(data[0] - TYPE(1)) < tollerance && cyAbs(data[1]) < tollerance && cyAbs(data[2]) < tollerance && cyAbs(data[3] - TYPE(1)) < tollerance; }

	//! Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( TYPE tollerance=TYPE(0.001) ) const { return cyAbs(data[0] - data[2]) < tollerance; }


	//////////////////////////////////////////////////////////////////////////
	//!@name Static Methods

	//! Returns an identity matrix
	static Matrix2 MatrixIdentity() { Matrix2 m; m.SetIdentity(); return m; }
	//! Returns a rotation matrix about the given axis by angle in radians
	static Matrix2 MatrixRotation( TYPE angle ) { Matrix2 m; m.SetRotation(angle); return m; }
	//! Returns a uniform scale matrix
	static Matrix2 MatrixScale( const TYPE &uniformScale ) { Matrix2 m; m.SetScale(uniformScale); return m; }
	//! Returns a scale matrix
	static Matrix2 MatrixScale( const TYPE &scaleX, const TYPE &scaleY ) { Matrix2 m; m.SetScale(scaleX,scaleY); return m; }
	//! Returns a scale matrix
	static Matrix2 MatrixScale( const Point2<TYPE> &scale ) { Matrix2 m; m.SetScale(scale); return m; }


	//////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

//! 3x3 matrix class.
//!
//! Its data stores 9-value array of column-major matrix elements.
//! You can use Matrix3 with Point3<TYPE> to transform 3D points.

template <typename TYPE>
class Matrix3
{
#ifdef CY_NONVECTORIZED_MATRIX3
	friend Matrix3 operator * ( const TYPE value, const Matrix3 &right ) { Matrix3 r; for (int i=0; i<9; i++) r.data[i] = value * right.data[i]; return r; }	//!< multiply matrix by a value
#else
	friend Matrix3 operator * ( const TYPE value, const Matrix3 &right ) { Matrix3 r; _CY_IVDEP_FOR (int i=0; i<8; i++) r.data[i] = value * right.data[i]; r.data[8] = value * right.data[8]; return r; }	//!< multiply matrix by a value
#endif
	friend Matrix3 Inverse( const Matrix3 &m ) { return m.GetInverse(); }	//!< return the inverse of the matrix

public:

	//! Elements of the matrix are column-major: \n
	//! | 0  3  6 | \n
	//! | 1  4  7 | \n
	//! | 2  5  8 | \n
	TYPE data[9];

	//////////////////////////////////////////////////////////////////////////
	//!@name Constructors

	Matrix3() {}																							//!< Default constructor
	Matrix3( const Matrix3 &matrix ) { CY_MEMCOPY(TYPE,data,matrix.data,9); }								//!< Copy constructor
	template <typename T> explicit Matrix3<TYPE>( const Matrix3<T> &matrix ) { CY_MEMCONVERT(TYPE,data,matrix.data,9); }		//!< Copy constructor for different types
	explicit Matrix3( const TYPE *values ) { Set(values); }													//!< Initialize the matrix using an array of 9 values
	explicit Matrix3( const TYPE &v ) { SetScaledIdentity(v); }												//!< Initialize the matrix as identity scaled by v
	explicit Matrix3( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z ) { Set(x,y,z); }	//!< Initialize the matrix using x,y,z vectors as columns
	explicit Matrix3( const Matrix2<TYPE> &m ) { 
		data[0] = m.data[0]; data[1] = m.data[1]; data[2] = TYPE(0);
		data[3] = m.data[2]; data[4] = m.data[3]; data[5] = TYPE(0);
		data[6] = TYPE(0);   data[7] = TYPE(0);   data[8] = TYPE(1);
	}
	explicit Matrix3( const Matrix34<TYPE> &m );
	explicit Matrix3( const Matrix4<TYPE>  &m );


	//////////////////////////////////////////////////////////////////////////
	//!@name Set & Get Methods

	//! Set all the values as zero
	void Zero() { CY_MEMCLEAR(TYPE,data,9); }
	//! Returns true if the matrix is exactly zero
	bool IsZero() const { for ( int i=0; i<9; i++ ) if ( data[i] != 0 ) return false; return true; }
	//! Copies the matrix data to the given values array of size 9
	void Get( TYPE *values ) { CY_MEMCOPY(TYPE,values,data,9); } 
	//! Set matrix using an array of 9 values
	void Set( const TYPE *values ) { CY_MEMCOPY(TYPE,data,values,9); } 
	//! Set matrix using x,y,z vectors as columns
	void Set( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z ) { x.Get(&data[0]); y.Get(&data[3]); z.Get(&data[6]); }
	//! Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	//! Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity(TYPE v) { SetScale(v); }
	//! Sets the matrix as the tensor product (outer product) of two vectors
	void SetTensorProduct( const Point3<TYPE> &v0, const Point3<TYPE> &v1 )
	{
		for ( int i=0; i<3; i++ ) data[  i] = v0[i] * v1.x;
		for ( int i=0; i<3; i++ ) data[3+i] = v0[i] * v1.y;
		for ( int i=0; i<3; i++ ) data[6+i] = v0[i] * v1.z;
	}
	//! Matrix representation of the cross product ( a x b)
	void SetCrossProd( const Point3<TYPE> &p ) { data[0]=TYPE(0); data[1]=p.z; data[2]=-p.y; data[3]=-p.z; data[4]=TYPE(0); data[5]=p.x; data[6]=p.y; data[7]=-p.x; data[8]=TYPE(0); }


	//////////////////////////////////////////////////////////////////////////
	//!@name Affine transformations

	//! Sets a uniform scale matrix
	void SetScale( const TYPE &uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	//! Sets a scale matrix
	void SetScale( const TYPE &scaleX, const TYPE &scaleY, const TYPE &scaleZ )
	{
		data[0] = scaleX; data[1] = 0;      data[2]=0;     
		data[3] = 0;      data[4] = scaleY; data[5]=0;     
		data[6] = 0;      data[7] = 0;      data[8]=scaleZ;
	}
	//! Sets a scale matrix
	void SetScale( const Point3<TYPE> &scale ) { SetScale(scale.x,scale.y,scale.z); }
	//! Removes the scale component of the matrix
	void SetNoScale() { Point3<TYPE> *p = (Point3<TYPE>*)data; p[0].Normalize(); p[1].Normalize(); p[2].Normalize(); }
	//! Set as rotation matrix around x axis
	void SetRotationX( TYPE angle ) { SetRotationX( cySin(angle), cyCos(angle) ); }
	//! Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = TYPE(1);   data[1] =  TYPE(0);    data[2] = TYPE(0); 
		data[3] = TYPE(0);   data[4] =  cosAngle;   data[5] = sinAngle;
		data[6] = TYPE(0);   data[7] = -sinAngle;   data[8] = cosAngle;
	}
	//! Set as rotation matrix around y axis
	void SetRotationY( TYPE angle ) { SetRotationY( cySin(angle), cyCos(angle) ); }
	//! Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = cosAngle;   data[1] = TYPE(0);   data[2] = -sinAngle;
		data[3] = TYPE(0);    data[4] = TYPE(1);   data[5] =  TYPE(0); 
		data[6] = sinAngle;   data[7] = TYPE(0);   data[8] =  cosAngle;
	}
	//! Set as rotation matrix around z axis
	void SetRotationZ( TYPE angle ) { SetRotationZ( cySin(angle), cyCos(angle) ); }
	//! Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] =  cosAngle;   data[1] = sinAngle;   data[2] = TYPE(0);
		data[3] = -sinAngle;   data[4] = cosAngle;   data[5] = TYPE(0);
		data[6] =  TYPE(0);	   data[7] = TYPE(0);    data[8] = TYPE(1);
	}
	//! Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
	void SetRotationXYZ( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = cySin(angleX);
		const TYPE cx = cyCos(angleX);
		const TYPE sy = cySin(angleY);
		const TYPE cy = cyCos(angleY);
		const TYPE sz = cySin(angleZ);
		const TYPE cz = cyCos(angleZ);
		data[0] = cy*cz; 		      data[1] = cy*sz; 			    data[2] =-sy;   
		data[3] = cz*sx*sy - cx*sz;   data[4] = cx*cz + sx*sy*sz;   data[5] = cy*sx;
		data[6] = cx*cz*sy + sx*sz;   data[7] =-cz*sx + cx*sy*sz;   data[8] = cx*cy;
	}
	//! Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
	void SetRotationZYX( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = cySin(angleX);
		const TYPE cx = cyCos(angleX);
		const TYPE sy = cySin(angleY);
		const TYPE cy = cyCos(angleY);
		const TYPE sz = cySin(angleZ);
		const TYPE cz = cyCos(angleZ);
		data[0] = cy*cz;    data[1] = cx*sz + sx*sy*cz;   data[2] = sx*sz - cx*sy*cz;
		data[3] = -cy*sz;   data[4] = cx*cz - sx*sy*sz;   data[5] = sx*cz + cx*sy*sz;
		data[6] = sy;	    data[7] = -sx*cy;		      data[8] = cx*cy;
	}
	//! Set a rotation matrix about the given axis by angle
	void SetRotation( const Point3<TYPE> &axis, TYPE angle ) { SetRotation(axis,cySin(angle),cyCos(angle)); }
	//! Set a rotation matrix about the given axis by cos and sin of angle
	void SetRotation( const Point3<TYPE> &axis, TYPE sinAngle, TYPE cosAngle )
	{
		const TYPE t = TYPE(1) - cosAngle;
		const TYPE tx = t * axis.x;
		const TYPE ty = t * axis.y;
		const TYPE tz = t * axis.z;
		const TYPE txy = tx * axis.y;
		const TYPE txz = tx * axis.z;
		const TYPE tyz = ty * axis.z;
		const TYPE sx = sinAngle * axis.x;
		const TYPE sy = sinAngle * axis.y;
		const TYPE sz = sinAngle * axis.z;

		data[0] = tx * axis.x + cosAngle;   data[1] = txy + sz;                 data[2] = txz - sy;
		data[3] = txy - sz;                 data[4] = ty * axis.y + cosAngle;   data[5] = tyz + sx;
		data[6] = txz + sy;                 data[7] = tyz - sx;                 data[8] = tz * axis.z + cosAngle;
	}
	//! Set a rotation matrix that sets [from] unit vector to [to] unit vector
	void SetRotation( const Point3<TYPE> &from, const Point3<TYPE> &to )
	{
		TYPE c = from.Dot(to);
		if ( c > TYPE(0.9999999) ) SetIdentity();
		else {
			TYPE s = cySqrt(TYPE(1) - c*c);
			Point3<TYPE> axis = from.Cross(to).GetNormalized();
			SetRotation(axis, s, c);
		}
	}
	//! Set view matrix using position, target and approximate up vector
	void SetView( const Point3<TYPE> &target, const Point3<TYPE> &up )
	{
		Point3<TYPE> f = target;
		f.Normalize();
		Point3<TYPE> s = f.Cross(up);
		s.Normalize();
		Point3<TYPE> u = s.Cross(f);
		data[0] = s.x; data[1] = u.x; data[2] = -f.x;
		data[3] = s.y; data[4] = u.y; data[5] = -f.y;
		data[6] = s.z; data[7] = u.z; data[8] = -f.z;
	}
	//! Set matrix using normal, and approximate x direction
	void SetNormal(const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y = normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal); }


	//////////////////////////////////////////////////////////////////////////
	//!@name Set Row, Column, or Diagonal

	void SetRow   ( int row,    TYPE x, TYPE y, TYPE z ) { data[row]=x;      data[row+3]=y;      data[row+6]=z; }		//!< Sets a row of the matrix
	void SetColumn( int column, TYPE x, TYPE y, TYPE z ) { data[3*column]=x; data[3*column+1]=y; data[3*column+2]=z; }	//!< Sets a column of the matrix
	void SetDiagonal( const TYPE &xx, const TYPE &yy, const TYPE &zz ) { data[0]=xx; data[4]=yy; data[8]=zz; }			//!< Sets the diagonal values of the matrix
	void SetDiagonal( const Point3<TYPE> &p ) { SetDiagonal( p.x, p.y, p.z ); }											//!< Sets the diagonal values of the matrix
	void SetDiagonal( const TYPE *values )    { SetDiagonal(values[0],values[1],values[2]); }							//!< Sets the diagonal values of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Row, Column, or Diagonal
	
	Point3<TYPE>  GetRow   ( int row )                  const { return Point3<TYPE>( data[row], data[row+3], data[row+6] ); }				//!< Returns a row of the matrix
	void          GetRow   ( int row, Point3<TYPE> &p ) const { p.Set( data[row], data[row+3], data[row+6] ); }								//!< Returns a row of the matrix
	void          GetRow   ( int row, TYPE *values )    const { values[0]=data[row]; values[1]=data[row+3]; values[2]=data[row+6]; }		//!< Returns a row of the matrix
	Point3<TYPE>  GetColumn( int col )                  const { return Point3<TYPE>( &data[col*3] ); }										//!< Returns a column of the matrix
	void          GetColumn( int col, Point3<TYPE> &p ) const { p.Set( &data[col*3] ); }													//!< Returns a column of the matrix
	void          GetColumn( int col, TYPE *values )    const { values[0]=data[col*3]; values[1]=data[col*3+1]; values[2]=data[col*3+2]; }	//!< Returns a column of the matrix
	Point3<TYPE>  GetDiagonal()                         const { Point3<TYPE> r; GetDiagonal(r); return r; }									//!< Returns the diagonal of the matrix
	void          GetDiagonal( Point3<TYPE> &p )        const { GetDiagonal(&p.x); }														//!< Returns the diagonal of the matrix
	void	      GetDiagonal( TYPE *values )           const { values[0]=data[0]; values[1]=data[4]; values[2]=data[8]; }					//!< Returns the diagonal of the matrix

	
	//////////////////////////////////////////////////////////////////////////
	//!@name Get Sub-matrix data
	
	void          GetSubMatrix ( Matrix2<TYPE> &m )     const { GetSubMatrix2(m); }														//!< Returns the 2x2 portion of the matrix
	Matrix2<TYPE> GetSubMatrix2()                       const { Matrix2<TYPE> m; GetSubMatrix2(m.data); return m; }						//!< Returns the 2x2 portion of the matrix
	void          GetSubMatrix2( Matrix2<TYPE> &m )     const { GetSubMatrix2(m.data); }												//!< Returns the 2x2 portion of the matrix
	void          GetSubMatrix2( TYPE *mdata )          const { CY_MEMCOPY(TYPE,mdata,data,2); CY_MEMCOPY(TYPE,mdata+2,data+3,2); }		//!< Returns the 2x2 portion of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Comparison Operators

	bool operator == ( const Matrix3 &right ) const { for ( int i=0; i<9; i++ ) if ( data[i] != right.data[i] ) return false; return true; } //!< compare equal
	bool operator != ( const Matrix3 &right ) const { for ( int i=0; i<9; i++ ) if ( data[i] != right.data[i] ) return true; return false; } //!< compare not equal


	//////////////////////////////////////////////////////////////////////////
	//!@name Access Operators

	TYPE&       operator () ( int row, int column )       { return data[ column * 3 + row ]; }	//!< subscript operator
	const TYPE& operator () ( int row, int column ) const { return data[ column * 3 + row ]; }	//!< constant subscript operator
	TYPE&       operator [] ( int i )       { return data[i]; }	//!< subscript operator
	const TYPE& operator [] ( int i ) const { return data[i]; }	//!< constant subscript operator
	

	//////////////////////////////////////////////////////////////////////////
	//!@name Unary and Binary Operators

	// Unary operators
#ifdef CY_NONVECTORIZED_MATRIX3
	Matrix3 operator - () const { Matrix3 r; for (int i=0; i<9; i++) r.data[i] = -data[i]; return r; }	//!< negative matrix
#else
	Matrix3 operator - () const	{ Matrix3 r; _CY_IVDEP_FOR (int i=0; i<8; i++) r.data[i] = -data[i]; r.data[8] = -data[8]; return r; }	//!< negative matrix
#endif

	// Binary operators
#ifdef CY_NONVECTORIZED_MATRIX3
	Matrix3 operator * ( const TYPE    &value ) const { Matrix3 r; for (int i=0; i<9; i++) r.data[i] = data[i] * value;         return r; }	//!< multiply matrix by a value
	Matrix3 operator / ( const TYPE    &value ) const { Matrix3 r; for (int i=0; i<9; i++) r.data[i] = data[i] / value;         return r; }	//!< divide matrix by a value;
	Matrix3 operator + ( const Matrix3 &right ) const { Matrix3 r; for (int i=0; i<9; i++) r.data[i] = data[i] + right.data[i]; return r; }	//!< add two Matrices
	Matrix3 operator - ( const Matrix3 &right ) const { Matrix3 r; for (int i=0; i<9; i++) r.data[i] = data[i] - right.data[i]; return r; }	//!< subtract one Matrix3 from another
#else
	Matrix3 operator * ( const TYPE    &value ) const { Matrix3 r; _CY_IVDEP_FOR (int i=0; i<8; i++) r.data[i] = data[i] * value;         r.data[8] = data[8] * value;         return r; }	//!< multiply matrix by a value
	Matrix3 operator / ( const TYPE    &value ) const { Matrix3 r; _CY_IVDEP_FOR (int i=0; i<8; i++) r.data[i] = data[i] / value;         r.data[8] = data[8] / value;         return r; }	//!< divide matrix by a value;
	Matrix3 operator + ( const Matrix3 &right ) const { Matrix3 r; _CY_IVDEP_FOR (int i=0; i<8; i++) r.data[i] = data[i] + right.data[i]; r.data[8] = data[8] + right.data[8]; return r; }	//!< add two Matrices
	Matrix3 operator - ( const Matrix3 &right ) const { Matrix3 r; _CY_IVDEP_FOR (int i=0; i<8; i++) r.data[i] = data[i] - right.data[i]; r.data[8] = data[8] - right.data[8]; return r; }	//!< subtract one Matrix3 from another
#endif

	Matrix3 operator * ( const Matrix3 &right ) const	//!< multiply a matrix with another
	{
		Matrix3 r;
		TYPE *rd = r.data;
		for ( int i=0; i<9; i+=3, rd+=3 ) {
			TYPE a[3], b[3], c[3];
			for ( int j=0; j<3; ++j ) a [j] = data[  j] * right.data[i  ];
			for ( int j=0; j<3; ++j ) b [j] = data[3+j] * right.data[i+1];
			for ( int j=0; j<3; ++j ) c [j] = data[6+j] * right.data[i+2];
			for ( int j=0; j<3; ++j ) rd[j] = a[j] + b[j] + c[j];
		}
		return r;
	}
	Point3<TYPE> operator * ( const Point3<TYPE> &p ) const
	{
		return Point3<TYPE>( p.x*data[0] + p.y*data[3] + p.z*data[6], 
							 p.x*data[1] + p.y*data[4] + p.z*data[7],
							 p.x*data[2] + p.y*data[5] + p.z*data[8] );
		//TYPE a[3], b[3], c[3];
		//Point3<TYPE> rr;
		//for ( int i=0; i<3; ++i ) a [i] = p[0] * data[  i];
		//for ( int i=0; i<3; ++i ) b [i] = p[1] * data[3+i];
		//for ( int i=0; i<3; ++i ) c [i] = p[2] * data[6+i];
		//for ( int i=0; i<3; ++i ) rr[i] = a[i] + b[i] + c[i];	
		//return rr;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Assignment Operators

	const Matrix3& operator  = ( const Matrix3 &right ) { CY_MEMCOPY(TYPE,data,right.data,9); return *this; }				//!< assignment operator
	const Matrix3& operator *= ( const Matrix3 &right ) { *this = operator*(right);           return *this; }				//!< multiply a matrix with another matrix and modify this matrix
#ifdef CY_NONVECTORIZED_MATRIX3
	const Matrix3& operator += ( const Matrix3 &right ) { for (int i=0; i<9; i++) data[i] += right.data[i]; return *this; }	//!< add two Matrices modify this
	const Matrix3& operator -= ( const Matrix3 &right ) { for (int i=0; i<9; i++) data[i] -= right.data[i]; return *this; }	//!< subtract one Matrix3 from another matrix and modify this matrix
	const Matrix3& operator *= ( const TYPE    &value ) { for (int i=0; i<9; i++) data[i] *= value;         return *this; }	//!< multiply a matrix with a value modify this matrix
	const Matrix3& operator /= ( const TYPE    &value ) { for (int i=0; i<9; i++) data[i] /= value;         return *this; }	//!< divide the matrix by a value modify the this matrix
#else
	const Matrix3& operator += ( const Matrix3 &right ) { _CY_IVDEP_FOR (int i=0; i<8; i++) data[i] += right.data[i]; data[8] += right.data[8]; return *this; }	//!< add two Matrices modify this
	const Matrix3& operator -= ( const Matrix3 &right ) { _CY_IVDEP_FOR (int i=0; i<8; i++) data[i] -= right.data[i]; data[8] -= right.data[8]; return *this; }	//!< subtract one Matrix3 from another matrix and modify this matrix
	const Matrix3& operator *= ( const TYPE    &value ) { _CY_IVDEP_FOR (int i=0; i<8; i++) data[i] *= value;         data[8] *= value;         return *this; }	//!< multiply a matrix with a value modify this matrix
	const Matrix3& operator /= ( const TYPE    &value ) { _CY_IVDEP_FOR (int i=0; i<8; i++) data[i] /= value;         data[8] /= value;         return *this; }	//!< divide the matrix by a value modify the this matrix
#endif

	//////////////////////////////////////////////////////////////////////////
	//!@name Other Methods

	void Transpose()														//!< Transpose this matrix
	{
		for (int i = 1; i < 3; i++) {
			for (int j = 0; j < i; j++) {
				TYPE temp = data[i * 3 + j];
				data[i * 3 + j] = data[j * 3 + i];
				data[j * 3 + i] = temp;
			}
		}
	}
	void GetTranspose( Matrix3 &m ) const									//!< return Transpose of this matrix
	{
		m.data[0] = data[0];   m.data[1] = data[3];   m.data[2] = data[6];
		m.data[3] = data[1];   m.data[4] = data[4];   m.data[5] = data[7];
		m.data[6] = data[2];   m.data[7] = data[5];   m.data[8] = data[8];
	}
	Matrix3 GetTranspose() const { Matrix3 t; GetTranspose(t); return t; }	//!< return Transpose of this matrix

	//! Multiply the give vector with the transpose of the matrix
	Point3<TYPE> TransposeMult( const Point3<TYPE> &p ) const
	{
		return Point3<TYPE>( p.x*data[0] + p.y*data[1] + p.z*data[2], 
							 p.x*data[3] + p.y*data[4] + p.z*data[5],
							 p.x*data[6] + p.y*data[7] + p.z*data[8] );
	}

	TYPE GetDeterminant() const {	//!< Get the determinant of this matrix
		// 0 (4 8 - 5 7) + 1 (5 6 - 3 8) + 2 (3 7 - 4 6)
		return data[0] * ( data[4] * data[8] - data[5] * data[7] ) + 
		       data[1] * ( data[5] * data[6] - data[3] * data[8] ) + 
		       data[2] * ( data[3] * data[7] - data[4] * data[6] );
	}

	void Invert() { Matrix3 inv; GetInverse(inv); *this=inv; }				//!< Invert this matrix
	void GetInverse( Matrix3 &inverse ) const								//!< Get the inverse of this matrix
	{
		//  ( 4 8 - 5 7    5 6 - 3 8    3 7 - 4 6 ) 
		//  ( 2 7 - 1 8    0 8 - 2 6    1 6 - 0 7 )  / det
		//  ( 1 5 - 2 4    2 3 - 0 5    0 4 - 1 3 ) 

		inverse.data[0] = (data[4]*data[8] - data[5]*data[7]);
		inverse.data[1] = (data[2]*data[7] - data[1]*data[8]);
		inverse.data[2] = (data[1]*data[5] - data[2]*data[4]);

		inverse.data[3] = (data[5]*data[6] - data[3]*data[8]);
		inverse.data[4] = (data[0]*data[8] - data[2]*data[6]);
		inverse.data[5] = (data[2]*data[3] - data[0]*data[5]);

		inverse.data[6] = (data[3]*data[7] - data[4]*data[6]);
		inverse.data[7] = (data[1]*data[6] - data[0]*data[7]);
		inverse.data[8] = (data[0]*data[4] - data[1]*data[3]);

		TYPE det = data[0] * inverse.data[0] + data[1] * inverse.data[3] + data[2] * inverse.data[6];
		inverse /= det;

	}
	Matrix3 GetInverse() const { Matrix3 inv; GetInverse(inv); return inv; }	//!< Get the inverse of this matrix

	//! Orthogonalizes the matrix and removes the scale component, preserving the x direction
	void OrthogonalizeX()
	{
		Point3<TYPE> *p = (Point3<TYPE>*)data;
		p[0].Normalize();
		p[1] -= p[0] * (p[1]%p[0]);
		p[1].Normalize();
		p[2] -= p[0] * (p[2]%p[0]);
		p[2] -= p[1] * (p[2]%p[1]);
		p[2].Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the y direction
	void OrthogonalizeY()
	{
		Point3<TYPE> *p = (Point3<TYPE>*)data;
		p[1].Normalize();
		p[0] -= p[1] * (p[0]%p[1]);
		p[0].Normalize();
		p[2] -= p[1] * (p[2]%p[1]);
		p[2] -= p[0] * (p[2]%p[0]);
		p[2].Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the z direction
	void OrthogonalizeZ()
	{
		Point3<TYPE> *p = (Point3<TYPE>*)data;
		p[2].Normalize();
		p[0] -= p[2] * (p[0]%p[2]);
		p[0].Normalize();
		p[1] -= p[2] * (p[1]%p[2]);
		p[1] -= p[0] * (p[1]%p[0]);
		p[1].Normalize();
	}


	//! Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const
	{
		return cyAbs(data[0]-TYPE(1)) < tollerance && cyAbs(data[1])         < tollerance && cyAbs(data[2])         < tollerance && 
			   cyAbs(data[3])         < tollerance && cyAbs(data[4]-TYPE(1)) < tollerance && cyAbs(data[5])         < tollerance &&
			   cyAbs(data[6])         < tollerance && cyAbs(data[7])         < tollerance && cyAbs(data[8]-TYPE(1)) < tollerance;
	}

	//! Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( TYPE tollerance=TYPE(0.001) ) const { return cyAbs(data[1] - data[3]) < tollerance && cyAbs(data[2] - data[6]) < tollerance && cyAbs(data[5] - data[7]) < tollerance; }

	
	//////////////////////////////////////////////////////////////////////////
	//!@name Static Methods

	//! Returns an identity matrix
	static Matrix3 MatrixIdentity() { Matrix3 m; m.SetIdentity(); return m; }
	//! Returns a view matrix using position, target and approximate up vector
	static Matrix3 MatrixView( const Point3<TYPE> &target, Point3<TYPE> &up ) { Matrix3 m; m.SetView(target,up); return m; }
	//! Returns a matrix using normal, and approximate x direction
	static Matrix3 MatrixNormal( const Point3<TYPE> &normal, Point3<TYPE> &dir ) { Matrix3 m; m.SetNormal(normal,dir); return m; }
	//! Returns a rotation matrix around x axis by angle in radians
	static Matrix3 MatrixRotationX( TYPE angle ) { Matrix3 m; m.SetRotationX(angle); return m; }
	//! Returns a rotation matrix around y axis by angle in radians
	static Matrix3 MatrixRotationY( TYPE angle ) { Matrix3 m; m.SetRotationY(angle); return m; }
	//! Returns a rotation matrix around z axis by angle in radians
	static Matrix3 MatrixRotationZ( TYPE angle ) { Matrix3 m; m.SetRotationZ(angle); return m; }
	//! Returns a rotation matrix around x, y, and then z axes by angle in radians (Rz * Ry * Rx)
	static Matrix3 MatrixRotationXYZ( TYPE angleX, TYPE angleY, TYPE angleZ ) { Matrix3 m; m.SetRotationXYZ(angleX,angleY,angleZ); return m; }
	//! Returns a rotation matrix around z, y, and then x axes by angle in radians (Rx * Ry * Rz)
	static Matrix3 MatrixRotationZYX( TYPE angleX, TYPE angleY, TYPE angleZ ) { Matrix3 m; m.SetRotationZYX(angleX,angleY,angleZ); return m; }
	//! Returns a rotation matrix about the given axis by angle in radians
	static Matrix3 MatrixRotation( const Point3<TYPE> &axis, TYPE angle ) { Matrix3 m; m.SetRotation(axis,angle); return m; }
	//! Returns a rotation matrix that sets [from] unit vector to [to] unit vector
	static Matrix3 MatrixRotation( const Point3<TYPE> &from, const Point3<TYPE> &to ) { Matrix3 m; m.SetRotation(from,to); return m; }
	//! Returns a uniform scale matrix
	static Matrix3 MatrixScale( const TYPE &uniformScale ) { Matrix3 m; m.SetScale(uniformScale); return m; }
	//! Returns a scale matrix
	static Matrix3 MatrixScale( const TYPE &scaleX, const TYPE &scaleY, const TYPE &scaleZ ) { Matrix3 m; m.SetScale(scaleX,scaleY,scaleZ); return m; }
	//! Returns a scale matrix
	static Matrix3 MatrixScale( const Point3<TYPE> &scale ) { Matrix3 m; m.SetScale(scale); return m; }
	//! Returns the matrix representation of cross product ( a x b )
	static Matrix3 MatrixCrossProd( const Point3<TYPE> &a ) { Matrix3 m; m.SetCrossProd(a); return m; }

	//////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

//! 3x4 matrix class.
//!
//! Its data stores 12-value array of column-major matrix elements.
//! I chose column-major format to be compatible with OpenGL
//! You can use Matrix34 with Point3<TYPE> and Point4<TYPE>
//! to transform 3D and 4D points.

template <typename TYPE>
class Matrix34
{
	friend Matrix34 operator * ( const TYPE value, const Matrix34 &right ) { Matrix34 r; for (int i=0; i<12; i++) r.data[i] = value * right.data[i]; return r; }	//!< multiply matrix by a value
	friend Matrix34 Inverse( const Matrix34 &m ) { return m.GetInverse(); }	//!< return the inverse of the matrix

public:

	//! Elements of the matrix are column-major: \n
	//! | 0   3   6   9 | \n
	//! | 1   4   7  10 | \n
	//! | 2   5   8  11 | \n
	TYPE data[12];


	//////////////////////////////////////////////////////////////////////////
	//!@name Constructors

	Matrix34() {}																				//!< Default constructor
	Matrix34( const Matrix34 &matrix ) { CY_MEMCOPY(TYPE,data,matrix.data,12); }				//!< Copy constructor
	template <typename T> explicit Matrix34<TYPE>( const Matrix34<T> &matrix ) { CY_MEMCONVERT(TYPE,data,matrix.data,12); }		//!< Copy constructor for different types
	explicit Matrix34( const TYPE *values ) { Set(values); }									//!< Initialize the matrix using an array of 9 values
	explicit Matrix34( const TYPE &v ) { SetScaledIdentity(v); }								//!< Initialize the matrix as identity scaled by v
	explicit Matrix34( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z, const Point3<TYPE> &pos ) { Set(x,y,z,pos); }	//!< Initialize the matrix using x,y,z vectors and coordinate center
	explicit Matrix34( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Set(pos,normal,dir); }				//!< Initialize the matrix using position, normal, and approximate x direction
	explicit Matrix34( const Matrix3<TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,9); CY_MEMCLEAR(TYPE,data+9,3); }
	explicit Matrix34( const Matrix3<TYPE> &m, const Point3<TYPE> &pos ) { CY_MEMCOPY(TYPE,data,m.data,9); data[9]=pos.x; data[10]=pos.y; data[11]=pos.z; }
	explicit Matrix34( const Matrix2<TYPE> &m ) { 
		data[ 0] = m.data[ 0]; data[ 1] = m.data[ 1]; data[ 2] = TYPE(0);
		data[ 3] = m.data[ 2]; data[ 4] = m.data[ 3]; data[ 5] = TYPE(0);
		data[ 6] = TYPE(0);    data[ 7] = TYPE(0);    data[ 8] = TYPE(1);
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	explicit Matrix34( const Matrix4<TYPE> &m );


	//////////////////////////////////////////////////////////////////////////
	//!@name Set & Get Methods

	//! Set all the values as zero
	void Zero() { CY_MEMCLEAR(TYPE,data,12); }
	//! Returns true if the matrix is exactly zero
	bool IsZero() const { for ( int i=0; i<12; i++ ) if ( data[i] != 0 ) return false; return true; }
	//! Copies the matrix data to the given values array of size 12
	void Get( TYPE *values ) { CY_MEMCOPY(TYPE,values,data,12); } 
	//! Set Matrix using an array of 12 values
	void Set( const TYPE *values ) { CY_MEMCOPY(TYPE,data,values,12); } 
	//! Set matrix using x,y,z vectors and coordinate center
	void Set( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z, const Point3<TYPE> &pos ) { x.Get(data); y.Get(data+3); z.Get(data+6); pos.Get(data+9); }
	//! Set matrix using position, normal, and approximate x direction
	void Set( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,pos); }
	//! Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	//! Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity( TYPE v ) { SetScale(v); }


	//////////////////////////////////////////////////////////////////////////
	//!@name Affine transformations

	//! Sets a uniform scale matrix
	void SetScale( const TYPE &uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	//! Sets a scale matrix
	void SetScale( const TYPE &scaleX, const TYPE &scaleY, const TYPE &scaleZ )
	{
		data[ 0] = scaleX; data[ 1] = 0;      data[ 2]=0;     
		data[ 3] = 0;      data[ 4] = scaleY; data[ 5]=0;     
		data[ 6] = 0;      data[ 7] = 0;      data[ 8]=scaleZ;
		data[ 9] = 0;      data[10] = 0;      data[11]=0;
	}
	//! Sets a scale matrix
	void SetScale( const Point3<TYPE> &scale ) { SetScale(scale.x,scale.y,scale.z); }
	//! Removes the scale component of the matrix
	void SetNoScale() { Point3<TYPE> *p = (Point3<TYPE>*)data; p[0].Normalize(); p[1].Normalize(); p[2].Normalize(); }
	//! Set as rotation matrix around x axis
	void SetRotationX( TYPE angle ) { SetRotationX( cySin(angle), cyCos(angle) ); }
	//! Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = TYPE(1);   data[1] = TYPE(0);     data[2] = TYPE(0); 
		data[3] = TYPE(0);   data[4] = cosAngle;    data[5] = sinAngle;
		data[6] = TYPE(0);   data[7] = -sinAngle;   data[8] = cosAngle;
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	//! Set as rotation matrix around y axis
	void SetRotationY( TYPE angle ) { SetRotationY( cySin(angle), cyCos(angle) ); }
	//! Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = cosAngle;   data[1] = TYPE(0);   data[2] = -sinAngle;  
		data[3] = TYPE(0);    data[4] = TYPE(1);   data[5] = TYPE(0);  
		data[6] = sinAngle;   data[7] = TYPE(0);   data[8] = cosAngle; 
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	//! Set as rotation matrix around z axis
	void SetRotationZ( TYPE angle ) { SetRotationZ( cySin(angle), cyCos(angle) ); }
	//! Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] =  cosAngle;   data[1] = sinAngle;   data[2] = TYPE(0);
		data[3] = -sinAngle;   data[4] = cosAngle;   data[5] = TYPE(0); 
		data[6] =  TYPE(0);	   data[7] = TYPE(0);    data[8] = TYPE(1);
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	//! Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
	void SetRotationXYZ( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = cySin(angleX);
		const TYPE cx = cyCos(angleX);
		const TYPE sy = cySin(angleY);
		const TYPE cy = cyCos(angleY);
		const TYPE sz = cySin(angleZ);
		const TYPE cz = cyCos(angleZ);
		data[0] = cy*cz; 		      data[1] = cy*sz; 			    data[2] =-sy;   
		data[3] = cz*sx*sy - cx*sz;   data[4] = cx*cz + sx*sy*sz;   data[5] = cy*sx; 
		data[6] = cx*cz*sy + sx*sz;   data[7] =-cz*sx + cx*sy*sz;   data[8] = cx*cy;
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	//! Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
	void SetRotationZYX( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = cySin(angleX);
		const TYPE cx = cyCos(angleX);
		const TYPE sy = cySin(angleY);
		const TYPE cy = cyCos(angleY);
		const TYPE sz = cySin(angleZ);
		const TYPE cz = cyCos(angleZ);
		data[0] = cy*cz;   data[1] = cx*sz + sx*sy*cz;   data[2] = sx*sz - cx*sy*cz;            
		data[3] =-cy*sz;   data[4] = cx*cz - sx*sy*sz;   data[5] = sx*cz + cx*sy*sz;            
		data[6] = sy;      data[7] =-sx*cy;			     data[8] = cx*cy;
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	//! Set a rotation matrix about the given axis by angle
	void SetRotation( const Point3<TYPE> &axis, TYPE angle ) { SetRotation(axis,cySin(angle),cyCos(angle)); }
	//! Set a rotation matrix about the given axis by cos and sin of angle
	void SetRotation( const Point3<TYPE> &axis, TYPE sinAngle, TYPE cosAngle )
	{
		const TYPE t = TYPE(1) - cosAngle;
		const TYPE tx = t * axis.x;
		const TYPE ty = t * axis.y;
		const TYPE tz = t * axis.z;
		const TYPE txy = tx * axis.y;
		const TYPE txz = tx * axis.z;
		const TYPE tyz = ty * axis.z;
		const TYPE sx = sinAngle * axis.x;
		const TYPE sy = sinAngle * axis.y;
		const TYPE sz = sinAngle * axis.z;
		data[ 0] = tx * axis.x + cosAngle;   data[ 1] = txy + sz;                 data[ 2] = txz - sy;
		data[ 3] = txy - sz;                 data[ 4] = ty * axis.y + cosAngle;   data[ 5] = tyz + sx;
		data[ 6] = txz + sy;                 data[ 7] = tyz - sx;                 data[ 8] = tz * axis.z + cosAngle;
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	//! Set a rotation matrix that sets [from] unit vector to [to] unit vector
	void SetRotation( const Point3<TYPE> &from, const Point3<TYPE> &to )
	{
		TYPE c = from.Dot(to);
		if ( c > TYPE(0.9999999) ) SetIdentity();
		else {
			TYPE s = cySqrt(TYPE(1) - c*c);
			Point3<TYPE> axis = from.Cross(to).GetNormalized();
			SetRotation(axis, s, c);
		}
	}
	//! Sets a translation matrix with no rotation or scale
	void SetTrans( const Point3<TYPE> &move ) { TYPE d[12]={1,0,0, 0,1,0, 0,0,1 }; CY_MEMCOPY(TYPE,data,d,9); data[9]=move.x; data[10]=move.y; data[11]=move.z; }
	//! Adds a translation to the matrix
	void AddTrans( const Point3<TYPE> &move ) { data[9]+=move.x; data[10]+=move.y; data[11]+=move.z; }
	//! Sets the translation component of the matrix
	void SetTransComponent( const Point3<TYPE> &move ) { data[9]=move.x; data[10]=move.y; data[11]=move.z; }
	//! Set view matrix using position, target and approximate up vector
	void SetView( const Point3<TYPE> &pos, const Point3<TYPE> &target, const Point3<TYPE> &up )
	{
		Point3<TYPE> f = target - pos;
		f.Normalize();
		Point3<TYPE> s = f.Cross(up);
		s.Normalize();
		Point3<TYPE> u = s.Cross(f);
		data[ 0]=s.x; data[ 1]=u.x; data[ 2]=-f.x;
		data[ 3]=s.y; data[ 4]=u.y; data[ 5]=-f.y;
		data[ 6]=s.z; data[ 7]=u.z; data[ 8]=-f.z;
		data[ 9]= -s % pos;
		data[10]= -u % pos;
		data[11]=  f % pos;
	}
	//! Set matrix using normal and approximate x direction
	void SetNormal(const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,Point3<TYPE>(TYPE(0),TYPE(0),TYPE(0))); }


	//////////////////////////////////////////////////////////////////////////
	//!@name Set Row, Column, or Diagonal

	void SetRow( int row, TYPE x, TYPE y, TYPE z, TYPE w ) { data[row]=x; data[row+3]=y; data[row+6]=z; data[row+9]=w; }	//!< Sets a row of the matrix
	void SetColumn( int column, TYPE x, TYPE y, TYPE z ) { data[3*column]=x; data[3*column+1]=y; data[3*column+2]=z; }		//!< Sets a column of the matrix
	void SetDiagonal( const TYPE &xx, const TYPE &yy, const TYPE &zz ) { data[0]=xx; data[4]=yy; data[8]=zz; }				//!< Sets the diagonal values of the matrix
	void SetDiagonal( const Point3<TYPE> &p ) { SetDiagonal( p.x, p.y, p.z ); }												//!< Sets the diagonal values of the matrix
	void SetDiagonal( const TYPE *values ) { SetDiagonal(values[0],values[1],values[2]); }									//!< Sets the diagonal values of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Row, Column, or Diagonal

	Point4<TYPE>  GetRow   ( int row )                  const { return Point4<TYPE>( data[row], data[row+3], data[row+6], data[row+9] ); }					//!< Returns a row of the matrix
	void          GetRow   ( int row, Point4<TYPE> &p ) const { p.Set( data[row], data[row+3], data[row+6], data[row+9] ); }								//!< Returns a row of the matrix
	void          GetRow   ( int row, TYPE *values )    const { values[0]=data[row]; values[1]=data[row+3]; values[2]=data[row+6]; values[3]=data[row+9]; }	//!< Returns a row of the matrix
	Point3<TYPE>  GetColumn( int col )                  const { return Point3<TYPE>( &data[col*3] ); }														//!< Returns a column of the matrix
	void          GetColumn( int col, Point3<TYPE> &p ) const { p.Set( &data[col*3] ); }																	//!< Returns a column of the matrix
	void          GetColumn( int col, TYPE *values )    const { values[0]=data[col*3]; values[1]=data[col*3+1]; values[2]=data[col*3+2]; }					//!< Returns a column of the matrix
	Point3<TYPE>  GetDiagonal()                         const { Point3<TYPE> r; GetDiagonal(r); return r; }													//!< Returns the diagonal of the matrix
	void          GetDiagonal( Point3<TYPE> &p )        const { GetDiagonal(&p.x); }																		//!< Returns the diagonal of the matrix
	void	      GetDiagonal( TYPE *values )           const { values[0]=data[0]; values[1]=data[4]; values[2]=data[8]; }									//!< Returns the diagonal of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Sub-matrix data

	void          GetSubMatrix ( Matrix3<TYPE> &m )     const { GetSubMatrix3(m); }														//!< Returns the 3x3 portion of the matrix
	void          GetSubMatrix ( Matrix2<TYPE> &m )     const { GetSubMatrix2(m); }														//!< Returns the 2x2 portion of the matrix
	Matrix3<TYPE> GetSubMatrix3()                       const { Matrix3<TYPE> m; GetSubMatrix3(m.data); return m; }						//!< Returns the 3x3 portion of the matrix
	void          GetSubMatrix3( Matrix3<TYPE> &m )     const { GetSubMatrix3(m.data); }												//!< Returns the 3x3 portion of the matrix
	void          GetSubMatrix3( TYPE *mdata )          const { CY_MEMCOPY(TYPE,mdata,data,9); }										//!< Returns the 3x3 portion of the matrix
	Matrix2<TYPE> GetSubMatrix2()                       const { Matrix2<TYPE> m; GetSubMatrix2(m.data); return m; }						//!< Returns the 2x2 portion of the matrix
	void          GetSubMatrix2( Matrix2<TYPE> &m )     const { GetSubMatrix2(m.data); }												//!< Returns the 2x2 portion of the matrix
	void          GetSubMatrix2( TYPE *mdata )          const { CY_MEMCOPY(TYPE,mdata,data,2); CY_MEMCOPY(TYPE,mdata+2,data+3,2); }		//!< Returns the 2x2 portion of the matrix
	Point3<TYPE>  GetTrans()                            const { Point3<TYPE> p; GetTrans(p); return p; }								//! Returns the translation component of the matrix
	void          GetTrans( Point3<TYPE> &p )           const { p.x=data[9]; p.y=data[10]; p.z=data[11]; }								//! Returns the translation component of the matrix
	void          GetTrans( TYPE *trans )               const { CY_MEMCOPY(TYPE,trans,data+9,3); }										//! Returns the translation component of the matrix

	//////////////////////////////////////////////////////////////////////////
	//!@name Comparison Operators

	bool operator == ( const Matrix34 &right ) const { for ( int i=0; i<12; i++ ) if ( data[i] != right.data[i] ) return false; return true; } //!< compare equal
	bool operator != ( const Matrix34 &right ) const { for ( int i=0; i<12; i++ ) if ( data[i] != right.data[i] ) return true; return false; } //!< compare not equal


	//////////////////////////////////////////////////////////////////////////
	//!@name Access Operators

	TYPE&       operator () ( int row, int column )       { return data[ column * 3 + row ]; }	//!< subscript operator
	const TYPE& operator () ( int row, int column ) const { return data[ column * 3 + row ]; }	//!< constant subscript operator
	TYPE&       operator [] ( int i )       { return data[i]; }	//!< subscript operator
	const TYPE& operator [] ( int i ) const { return data[i]; }	//!< constant subscript operator


	//////////////////////////////////////////////////////////////////////////
	//!@name Unary and Binary Operators

	// Unary operators
	Matrix34 operator - () const { Matrix34 r; for (int i=0; i<12; i++) r.data[i]=-data[i]; return r; }	//!< negative matrix

	// Binary operators
	Matrix34 operator * ( const TYPE     &value ) const { Matrix34 r; for (int i=0; i<12; i++) r.data[i] = data[i] * value;         return r; }	//!< multiply matrix by a value
	Matrix34 operator / ( const TYPE     &value ) const { Matrix34 r; for (int i=0; i<12; i++) r.data[i] = data[i] / value;         return r; }	//!< divide matrix by a value;
	Matrix34 operator + ( const Matrix34 &right ) const { Matrix34 r; for (int i=0; i<12; i++) r.data[i] = data[i] + right.data[i]; return r; }	//!< add two Matrices
	Matrix34 operator - ( const Matrix34 &right ) const { Matrix34 r; for (int i=0; i<12; i++) r.data[i] = data[i] - right.data[i]; return r; }	//!< subtract one Matrix4 from another
	Matrix34 operator * ( const Matrix34 &right ) const	//!< multiply a matrix with another
	{
		Matrix34 r;
		TYPE *rd = r.data;
		for ( int i=0; i<12; i+=3, rd+=3 ) {
			TYPE a[3], b[3], c[3];
			for ( int k=0; k<3; ++k ) a[k] = data[  k] * right.data[i  ];
			for ( int k=0; k<3; ++k ) b[k] = data[3+k] * right.data[i+1];
			for ( int k=0; k<3; ++k ) c[k] = data[6+k] * right.data[i+2];
			for ( int j=0; j<3; ++j ) rd[j] = a[j] + b[j] + c[j];
		}
		for ( int j=0; j<3; ++j ) r.data[9+j] += data[9+j];
		return r;
	}
	Matrix34 operator * ( const Matrix3<TYPE> &right ) const	//!< multiply a matrix with another
	{
		Matrix34 r;
		TYPE *rd = r.data;
		for ( int i=0; i<9; i+=3, rd+=3 ) {
			TYPE a[3], b[3], c[3];
			for ( int k=0; k<3; ++k ) a[k] = data[  k] * right.data[i  ];
			for ( int k=0; k<3; ++k ) b[k] = data[3+k] * right.data[i+1];
			for ( int k=0; k<3; ++k ) c[k] = data[6+k] * right.data[i+2];
			for ( int j=0; j<3; ++j ) rd[j] = a[j] + b[j] + c[j];
		}
		CY_MEMCOPY(TYPE,r.data+9,data+9,3);
		return r;
	}
	Point3<TYPE> operator * ( const Point3<TYPE> &p ) const
	{
		//return Point3<TYPE>( p.x*data[0] + p.y*data[3] + p.z*data[6] + data[ 9], 
		//					 p.x*data[1] + p.y*data[4] + p.z*data[7] + data[10],
		//					 p.x*data[2] + p.y*data[5] + p.z*data[8] + data[11] );
		TYPE a[4], b[4], c[4];
		for ( int i=0; i<3; ++i ) a[i] = p[0] * data[  i];
		for ( int i=0; i<3; ++i ) b[i] = p[1] * data[3+i];
		for ( int i=0; i<3; ++i ) c[i] = p[2] * data[6+i];
		Point3<TYPE> rr;
		for ( int i=0; i<3; ++i ) rr[i] = a[i] + b[i] + c[i] + data[12+i];	
		return rr;
	}
	Point4<TYPE> operator * ( const Point4<TYPE> &p ) const
	{
		//return Point4<TYPE>( p.x*data[0] + p.y*data[3] + p.z*data[6] + p.w*data[ 9],
		//					 p.x*data[1] + p.y*data[4] + p.z*data[7] + p.w*data[10],
		//					 p.x*data[2] + p.y*data[5] + p.z*data[8] + p.w*data[11],
		//					 0           + 0           + 0           + p.w          );
		TYPE a[6], b[6];
		for ( int i=0; i<3; ++i ) a[  i] = p[0] * data[  i];
		for ( int i=0; i<3; ++i ) a[3+i] = p[1] * data[3+i];
		for ( int i=0; i<3; ++i ) b[  i] = p[2] * data[6+i];
		for ( int i=0; i<3; ++i ) b[3+i] = p[3] * data[9+i];
		for ( int i=0; i<6; ++i ) a[i] += b[i];
		Point4<TYPE> rr;
		for ( int i=0; i<3; ++i ) rr[i] = a[i] + a[3+i];
		rr.w = p.w;
		return rr;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Assignment Operators

	const Matrix34& operator  = ( const Matrix34 &right ) { CY_MEMCOPY(TYPE,data,right.data,12); return *this; }	
	const Matrix34& operator += ( const Matrix34 &right ) { for (int i=0; i<12; i++) data[i] += right.data[i]; return *this; }	//!< add two Matrices modify this
	const Matrix34& operator -= ( const Matrix34 &right ) { for (int i=0; i<12; i++) data[i] -= right.data[i]; return *this; }	//!< subtract one Matrix4 from another matrix and modify this matrix
	const Matrix34& operator *= ( const Matrix34 &right )      { *this = operator*(right); return *this; }						//!< multiply a matrix with another matrix and modify this matrix
	const Matrix34& operator *= ( const Matrix3<TYPE> &right ) { *this = operator*(right); return *this; }						//!< multiply a matrix with another matrix and modify this matrix
	const Matrix34& operator *= ( const TYPE     &value ) { for (int i=0; i<12; i++) data[i] *= value;         return *this; }	//!< multiply a matrix with a value modify this matrix
	const Matrix34& operator /= ( const TYPE     &value ) { for (int i=0; i<12; i++) data[i] /= value;         return *this; }	//!< divide the matrix by a value modify the this matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Other Methods

	//! Transpose this matrix
	void Transpose()
	{
		for (int i = 1; i < 3; i++) {
			for (int j = 0; j < i; j++) {
				TYPE temp = data[i * 3 + j];
				data[i * 3 + j] = data[j * 3 + i];
				data[j * 3 + i] = temp;
			}
		}
		CY_MEMCLEAR(TYPE,data+9,3);
	}

	//! return Transpose of this matrix
	void GetTranspose( Matrix34 &m ) const
	{
		m.data[0] = data[0];  m.data[1] = data[3];  m.data[2] = data[6];
		m.data[3] = data[1];  m.data[4] = data[4];  m.data[5] = data[7];
		m.data[6] = data[2];  m.data[7] = data[5];  m.data[8] = data[8];
		CY_MEMCLEAR(TYPE,m.data+9,3);
	}
	Matrix34 GetTranspose() const { Matrix34 t; GetTranspose(t); return t; }	//!< return Transpose of this matrix

	//! Multiply the give vector with the transpose of the matrix
	Point4<TYPE> TransposeMult( const Point3<TYPE> &p ) const
	{
		return Point4<TYPE>( p.x*data[ 0] + p.y*data[ 1] + p.z*data[ 2], 
							 p.x*data[ 3] + p.y*data[ 4] + p.z*data[ 5],
							 p.x*data[ 6] + p.y*data[ 7] + p.z*data[ 8],
							 p.x*data[ 9] + p.y*data[10] + p.z*data[11] + 1 );
	}

	//! Multiply the give vector with the transpose of the matrix
	Point4<TYPE> TransposeMult( const Point4<TYPE> &p ) const
	{
		return Point4<TYPE>( p.x*data[ 0] + p.y*data[ 1] + p.z*data[ 2], 
							 p.x*data[ 3] + p.y*data[ 4] + p.z*data[ 5],
							 p.x*data[ 6] + p.y*data[ 7] + p.z*data[ 8],
							 p.x*data[ 9] + p.y*data[10] + p.z*data[11] + p.w );
	}

	TYPE GetDeterminant() const	//!< Get the determinant of this matrix
	{
		// 0 (4 8 - 5 7) + 1 (5 6 - 3 8) + 2 (3 7 - 4 6)
		return data[0] * ( data[4] * data[8] - data[5] * data[7] ) + 
		       data[1] * ( data[5] * data[6] - data[3] * data[8] ) + 
		       data[2] * ( data[3] * data[7] - data[4] * data[6] );
	}
	void Invert() { Matrix34 inv; GetInverse(inv); *this=inv; }	//!< Invert this matrix
	void GetInverse( Matrix34 &inverse ) const					//!< Get the inverse of this matrix
	{
		//  (4 8 - 5 7)    (5 6 - 3 8)    (3 7 - 4 6)    (3 (8 10 - 7 11) + 4 (6 11 - 8  9) + 5 (7  9 - 6 10))
		//  (2 7 - 1 8)    (0 8 - 2 6)    (1 6 - 0 7)    (0 (7 11 - 8 10) + 1 (8  9 - 6 11) + 2 (6 10 -  7 9))    / det
		//  (1 5 - 2 4)    (2 3 - 0 5)    (0 4 - 1 3)    (0 (5 10 - 4 11) + 1 (3 11 - 5  9) + 2 (4  9 - 3 10))

		const TYPE data_8_10__7_11 = data[8] * data[10] - data[7] * data[11];
		const TYPE data_6_11__8__9 = data[6] * data[11] - data[8] * data[ 9];
		const TYPE data_7__9__6_10 = data[7] * data[ 9] - data[6] * data[10];

		inverse.data[ 0] = (data[4]*data[8] - data[5]*data[7]);
		inverse.data[ 1] = (data[2]*data[7] - data[1]*data[8]);
		inverse.data[ 2] = (data[1]*data[5] - data[2]*data[4]);

		inverse.data[ 3] = (data[5]*data[6] - data[3]*data[8]);
		inverse.data[ 4] = (data[0]*data[8] - data[2]*data[6]);
		inverse.data[ 5] = (data[2]*data[3] - data[0]*data[5]);

		inverse.data[ 6] = (data[3]*data[7] - data[4]*data[6]);
		inverse.data[ 7] = (data[1]*data[6] - data[0]*data[7]);
		inverse.data[ 8] = (data[0]*data[4] - data[1]*data[3]);

		inverse.data[ 9] = data[3] * data_8_10__7_11 + data[4] * data_6_11__8__9 + data[5] * data_7__9__6_10;
		inverse.data[10] = data[0] *-data_8_10__7_11 + data[1] *-data_6_11__8__9 + data[2] *-data_7__9__6_10;
		inverse.data[11] = data[0] * (data[5] * data[10] - data[4] * data[11]) + 
		                   data[1] * (data[3] * data[11] - data[5] * data[ 9]) +
		                   data[2] * (data[4] * data[ 9] - data[3] * data[10]);

		const TYPE det = data[0] * inverse.data[0] + data[1] * inverse.data[3] + data[2] * inverse.data[6];
		inverse /= det;
	}
	Matrix34 GetInverse() const { Matrix34 inv; GetInverse(inv); return inv; }	//!< Get the inverse of this matrix

	//! Orthogonalizes the matrix and removes the scale component, preserving the x direction
	void OrthogonalizeX()
	{
		Point3<TYPE> *p = (Point3<TYPE>*)data;
		p[0].Normalize();
		p[1] -= p[0] * (p[1]%p[0]);
		p[1].Normalize();
		p[2] -= p[0] * (p[2]%p[0]);
		p[2] -= p[1] * (p[2]%p[1]);
		p[2].Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the y direction
	void OrthogonalizeY()
	{
		Point3<TYPE> *p = (Point3<TYPE>*)data;
		p[1].Normalize();
		p[0] -= p[1] * (p[0]%p[1]);
		p[0].Normalize();
		p[2] -= p[1] * (p[2]%p[1]);
		p[2] -= p[0] * (p[2]%p[0]);
		p[2].Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the z direction
	void OrthogonalizeZ()
	{
		Point3<TYPE> *p = (Point3<TYPE>*)data;
		p[2].Normalize();
		p[0] -= p[2] * (p[0]%p[2]);
		p[0].Normalize();
		p[1] -= p[2] * (p[1]%p[2]);
		p[1] -= p[0] * (p[1]%p[0]);
		p[1].Normalize();
	}

	//! Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const
	{
		return cyAbs(data[0]-TYPE(1)) < tollerance && cyAbs(data[ 1])         < tollerance && cyAbs(data[ 2])         < tollerance && 
			   cyAbs(data[3])         < tollerance && cyAbs(data[ 4]-TYPE(1)) < tollerance && cyAbs(data[ 5])         < tollerance &&
			   cyAbs(data[6])         < tollerance && cyAbs(data[ 7])         < tollerance && cyAbs(data[ 8]-TYPE(1)) < tollerance &&
			   cyAbs(data[9])         < tollerance && cyAbs(data[10])         < tollerance && cyAbs(data[11])         < tollerance;
	}

	//! Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( TYPE tollerance=TYPE(0.001) ) const
	{
		return cyAbs(data[ 1] - data[3]) < tollerance && 
			   cyAbs(data[ 2] - data[6]) < tollerance &&
			   cyAbs(data[ 5] - data[7]) < tollerance &&
			   cyAbs(data[ 9])           < tollerance &&
			   cyAbs(data[10])           < tollerance &&
			   cyAbs(data[11])           < tollerance;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Static Methods

	//! Returns an identity matrix
	static Matrix34 MatrixIdentity() { Matrix34 m; m.SetIdentity(); return m; }
	//! Returns a matrix using normal, and approximate x direction
	static Matrix34 MatrixNormal( const Point3<TYPE> &normal, Point3<TYPE> &dir ) { Matrix34 m; m.SetNormal(normal,dir); return m; }
	//! Returns a rotation matrix around x axis by angle in radians
	static Matrix34 MatrixRotationX( TYPE angle ) { Matrix34 m; m.SetRotationX(angle); return m; }
	//! Returns a rotation matrix around y axis by angle in radians
	static Matrix34 MatrixRotationY( TYPE angle ) { Matrix34 m; m.SetRotationY(angle); return m; }
	//! Returns a rotation matrix around z axis by angle in radians
	static Matrix34 MatrixRotationZ( TYPE angle ) { Matrix34 m; m.SetRotationZ(angle); return m; }
	//! Returns a rotation matrix around x, y, and then z axes by angle in radians (Rz * Ry * Rx)
	static Matrix34 MatrixRotationXYZ( TYPE angleX, TYPE angleY, TYPE angleZ ) { Matrix34 m; m.SetRotationXYZ(angleX,angleY,angleZ); return m; }
	//! Returns a rotation matrix around z, y, and then x axes by angle in radians (Rx * Ry * Rz)
	static Matrix34 MatrixRotationZYX( TYPE angleX, TYPE angleY, TYPE angleZ ) { Matrix34 m; m.SetRotationZYX(angleX,angleY,angleZ); return m; }
	//! Returns a rotation matrix about the given axis by angle in radians
	static Matrix34 MatrixRotation( const Point3<TYPE> &axis, TYPE angle ) { Matrix34 m; m.SetRotation(axis,angle); return m; }
	//! Returns a rotation matrix that sets [from] unit vector to [to] unit vector
	static Matrix34 MatrixRotation( const Point3<TYPE> &from, const Point3<TYPE> &to ) { Matrix34 m; m.SetRotation(from,to); return m; }
	//! Returns a uniform scale matrix
	static Matrix34 MatrixScale( const TYPE &uniformScale ) { Matrix34 m; m.SetScale(uniformScale); return m; }
	//! Returns a scale matrix
	static Matrix34 MatrixScale( const TYPE &scaleX, const TYPE &scaleY, const TYPE &scaleZ ) { Matrix34 m; m.SetScale(scaleX,scaleY,scaleZ); return m; }
	//! Returns a scale matrix
	static Matrix34 MatrixScale( const Point3<TYPE> &scale ) { Matrix34 m; m.SetScale(scale); return m; }
	//! Returns a translation matrix with no rotation or scale
	static Matrix34 MatrixTrans( const Point3<TYPE> &move ) { Matrix34 m; m.SetTrans(move); return m; }


	//////////////////////////////////////////////////////////////////////////
};


//-------------------------------------------------------------------------------

//! 4x4 matrix class.
//!
//! Its data stores 16-value array of column-major matrix elements.
//! I chose column-major format to be compatible with OpenGL
//! You can use Matrix4 with Point3<TYPE> and Point4<TYPE>
//! to transform 3D and 4D points.

template <typename TYPE>
class Matrix4
{
	friend Matrix4 operator * ( const TYPE value, const Matrix4 &right ) { Matrix4 r; for (int i=0; i<16; i++) r.data[i] = value * right.data[i]; return r; }	//!< multiply matrix by a value
	friend Matrix4 Inverse( const Matrix4 &m ) { return m.GetInverse(); }	//!< return the inverse of the matrix

	//! multiply a 4x4 matrix with a 3x4 matrix, treating it as a 4x4 matrix with last row 0,0,0,1
	friend Matrix4 operator * ( const Matrix34<TYPE> &left, const Matrix4 &right )
	{
		Matrix4 r;
		TYPE *rd = r.data;
		for ( int i=0; i<16; i+=4, rd+=4 ) {
			TYPE a[3], b[3], c[3], d[3];
			for ( int j=0; j<3; ++j ) a[j] = left.data[  j] * right.data[i  ];
			for ( int j=0; j<3; ++j ) b[j] = left.data[3+j] * right.data[i+1];
			for ( int j=0; j<3; ++j ) c[j] = left.data[6+j] * right.data[i+2];
			for ( int j=0; j<3; ++j ) d[j] = left.data[9+j] * right.data[i+3];
			_CY_IVDEP_FOR ( int j=0; j<3; ++j ) rd[j] = a[j] + b[j] + c[j] + d[j];
		}
		r.data[ 3] = right.data[ 3];
		r.data[ 7] = right.data[ 7];
		r.data[11] = right.data[11];
		r.data[15] = right.data[15];
		return r;
	}

public:

	//! Elements of the matrix are column-major: \n
	//! | 0   4   8  12 | \n
	//! | 1   5   9  13 | \n
	//! | 2   6  10  14 | \n
	//! | 3   7  11  15 | \n
	TYPE data[16];


	//////////////////////////////////////////////////////////////////////////
	//!@name Constructors

	Matrix4() {}																	//!< Default constructor
	Matrix4( const Matrix4 &matrix ) { CY_MEMCOPY(TYPE,data,matrix.data,16); }		//!< Copy constructor
	template <typename T> explicit Matrix4<TYPE>( const Matrix4<T> &matrix ) { CY_MEMCONVERT(TYPE,data,matrix.data,16); }	//!< Copy constructor for different types
	explicit Matrix4( const TYPE *values ) { Set(values); }							//!< Initialize the matrix using an array of 9 values
	explicit Matrix4( const TYPE &v ) { SetScaledIdentity(v); }						//!< Initialize the matrix as identity scaled by v
	explicit Matrix4( const Point3<TYPE> &x,   const Point3<TYPE> &y,      const Point3<TYPE> &z, const Point3<TYPE> &pos ) { Set(x,y,z,pos); }	//!< Initialize the matrix using x,y,z vectors and coordinate center
	explicit Matrix4( const Point4<TYPE> &x,   const Point4<TYPE> &y,      const Point4<TYPE> &z, const Point4<TYPE> &w   ) { Set(x,y,z,w);   }	//!< Initialize the matrix using x,y,z vectors as columns
	explicit Matrix4( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Set(pos,normal,dir); }					//!< Initialize the matrix using position, normal, and approximate x direction
	explicit Matrix4( const Matrix34<TYPE> &m ) { 
		data[ 0]=m.data[ 0]; data[ 1]=m.data[ 1]; data[ 2]=m.data[ 2]; data[ 3]=TYPE(0); 
		data[ 4]=m.data[ 3]; data[ 5]=m.data[ 4]; data[ 6]=m.data[ 5]; data[ 7]=TYPE(0); 
		data[ 8]=m.data[ 6]; data[ 9]=m.data[ 7]; data[10]=m.data[ 8]; data[11]=TYPE(0); 
		data[12]=m.data[ 9]; data[13]=m.data[10]; data[14]=m.data[11]; data[15]=TYPE(1);
	}
	explicit Matrix4( const Matrix3<TYPE> &m ) { 
		data[ 0]=m.data[ 0]; data[ 1]=m.data[ 1]; data[ 2]=m.data[ 2]; data[ 3]=TYPE(0); 
		data[ 4]=m.data[ 3]; data[ 5]=m.data[ 4]; data[ 6]=m.data[ 5]; data[ 7]=TYPE(0); 
		data[ 8]=m.data[ 6]; data[ 9]=m.data[ 7]; data[10]=m.data[ 8]; 
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15]=TYPE(1);
	}
	explicit Matrix4( const Matrix2<TYPE> &m ) { 
		data[ 0]=m.data[ 0]; data[ 1]=m.data[ 1]; data[ 2]=TYPE(0); data[ 3]=TYPE(0); 
		data[ 4]=m.data[ 2]; data[ 5]=m.data[ 3]; 
		CY_MEMCLEAR(TYPE,data+6,4);
		data[10]=TYPE(1);
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15]=TYPE(1);
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Set & Get Methods

	//! Set all the values as zero
	void Zero() { CY_MEMCLEAR(TYPE,data,16); }
	//! Returns true if the matrix is exactly zero
	bool IsZero() const { for ( int i=0; i<16; i++ ) if ( data[i] != 0 ) return false; return true; }
	//! Copies the matrix data to the given values array of size 16
	void Get( TYPE *values ) { CY_MEMCOPY(TYPE,values,data,16); } 
	//! Set Matrix using an array of 16 values
	void Set( const TYPE *values ) { CY_MEMCOPY(TYPE,data,values,16); } 
	//! Set matrix using x,y,z vectors and coordinate center
	void Set( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z, const Point3<TYPE> &pos ) { x.Get(data); data[3]=TYPE(0); y.Get(data+4); data[7]=TYPE(0); z.Get(data+8); data[11]=TYPE(0); pos.Get(data+12); data[15]=TYPE(1); }
	//! Set matrix using x,y,z,w vectors
	void Set( const Point4<TYPE> &x, const Point4<TYPE> &y, const Point4<TYPE> &z, const Point4<TYPE> &w ) { x.Get(data); y.Get(data+4); z.Get(data+8); w.Get(data+12); }
	//! Set matrix using position, normal, and approximate x direction
	void Set( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,pos); }
	//! Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	//! Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity(TYPE v) { SetScale(v); }
	//! Sets the matrix as the tensor product (outer product) of two vectors
	void SetTensorProduct( const Point4<TYPE> &v0, const Point4<TYPE> &v1 )
	{
		for ( int i=0; i<4; ++i ) data[   i] = v0[i] * v1.x;	 // data[0]=v0.x*v1.x;  data[4]=v0.x*v1.y;  data[ 8]=v0.x*v1.z;  data[12]=v0.x*v1.w;
		for ( int i=0; i<4; ++i ) data[ 4+i] = v0[i] * v1.y;	 // data[1]=v0.y*v1.x;  data[5]=v0.y*v1.y;  data[ 9]=v0.y*v1.z;  data[13]=v0.y*v1.w;
		for ( int i=0; i<4; ++i ) data[ 8+i] = v0[i] * v1.z;	 // data[2]=v0.z*v1.x;  data[6]=v0.z*v1.y;  data[10]=v0.z*v1.z;  data[14]=v0.z*v1.w;
		for ( int i=0; i<4; ++i ) data[12+i] = v0[i] * v1.w;	 // data[3]=v0.w*v1.x;  data[7]=v0.w*v1.y;  data[11]=v0.w*v1.z;  data[15]=v0.w*v1.w;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Affine transformations

	//! Sets a uniform scale matrix
	void SetScale( const TYPE &uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	//! Sets a scale matrix
	void SetScale( const TYPE &scaleX, const TYPE &scaleY, const TYPE &scaleZ, TYPE scaleW=1 )
	{
		data[ 0] = scaleX; data[ 1] = 0;      data[ 2]=0;      data[ 3]=0; 
		data[ 4] = 0;      data[ 5] = scaleY; data[ 6]=0;      data[ 7]=0; 
		data[ 8] = 0;      data[ 9] = 0;      data[10]=scaleZ; data[11]=0;
		data[12] = 0;      data[13] = 0;      data[14]=0;      data[15]=scaleW;
	}
	//! Sets a scale matrix
	void SetScale( const Point3<TYPE> &scale ) { SetScale(scale.x,scale.y,scale.z); }
	//! Removes the scale component of the matrix
	void SetNoScale() { ((Point3<TYPE>*)&data[0])->Normalize(); ((Point3<TYPE>*)&data[4])->Normalize(); ((Point3<TYPE>*)&data[8])->Normalize(); }
	//! Set as rotation matrix around x axis
	void SetRotationX( TYPE angle ) { SetRotationX( cySin(angle), cyCos(angle) ); }
	//! Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( TYPE sinAngle, TYPE cosAngle )
	{
		data[ 0] = TYPE(1);  data[ 1] =  TYPE(0);    data[ 2] = TYPE(0);   data[ 3] = TYPE(0);
		data[ 4] = TYPE(0);  data[ 5] =  cosAngle;   data[ 6] = sinAngle;  data[ 7] = TYPE(0);
		data[ 8] = TYPE(0);  data[ 9] = -sinAngle;   data[10] = cosAngle;
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15] = TYPE(1);
	}
	//! Set as rotation matrix around y axis
	void SetRotationY( TYPE angle ) { SetRotationY( cySin(angle), cyCos(angle) ); }
	//! Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( TYPE sinAngle, TYPE cosAngle )
	{
		data[ 0] = cosAngle;  data[ 1] = TYPE(0);  data[ 2] = -sinAngle;  data[ 3] = TYPE(0);
		data[ 4] = TYPE(0);   data[ 5] = TYPE(1);  data[ 6] =  TYPE(0);   data[ 7] = TYPE(0);
		data[ 8] = sinAngle;  data[ 9] = TYPE(0);  data[10] =  cosAngle;
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15] = TYPE(1);
	}
	//! Set as rotation matrix around z axis
	void SetRotationZ( TYPE angle ) { SetRotationZ( cySin(angle), cyCos(angle) ); }
	//! Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( TYPE sinAngle, TYPE cosAngle )
	{
		data[ 0] =  cosAngle;  data[ 1] = sinAngle;  data[ 2] = TYPE(0);  data[ 3] = TYPE(0);
		data[ 4] = -sinAngle;  data[ 5] = cosAngle;  data[ 6] = TYPE(0);  data[ 7] = TYPE(0); 
		data[ 8] =  TYPE(0);   data[ 9] = TYPE(0);   data[10] = TYPE(1);
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15] = TYPE(1);
	}
	//! Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
	void SetRotationXYZ( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = cySin(angleX);
		const TYPE cx = cyCos(angleX);
		const TYPE sy = cySin(angleY);
		const TYPE cy = cyCos(angleY);
		const TYPE sz = cySin(angleZ);
		const TYPE cz = cyCos(angleZ);
		data[ 0] = cy*cz;             data[ 1] = cy*sz;             data[ 2] =-sy;     data[ 3] = TYPE(0);
		data[ 4] = cz*sx*sy - cx*sz;  data[ 5] = cx*cz + sx*sy*sz;  data[ 6] = cy*sx;  data[ 7] = TYPE(0);
		data[ 8] = cx*cz*sy + sx*sz;  data[ 9] =-cz*sx + cx*sy*sz;  data[10] = cx*cy;
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15] = TYPE(1);
	}
	//! Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
	void SetRotationZYX( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = cySin(angleX);
		const TYPE cx = cyCos(angleX);
		const TYPE sy = cySin(angleY);
		const TYPE cy = cyCos(angleY);
		const TYPE sz = cySin(angleZ);
		const TYPE cz = cyCos(angleZ);
		data[ 0] = cy*cz;  data[ 1] = cx*sz + sx*sy*cz;  data[ 2] = sx*sz - cx*sy*cz;  data[ 3] = TYPE(0);           
		data[ 4] =-cy*sz;  data[ 5] = cx*cz - sx*sy*sz;  data[ 6] = sx*cz + cx*sy*sz;  data[ 7] = TYPE(0);           
		data[ 8] = sy;     data[ 9] =-sx*cy;             data[10] = cx*cy;
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15] = TYPE(1);
	}
	//! Set a rotation matrix about the given axis by angle
	void SetRotation( const Point3<TYPE> &axis, TYPE angle ) { SetRotation(axis,cySin(angle),cyCos(angle)); }
	//! Set a rotation matrix about the given axis by cos and sin of angle
	void SetRotation( const Point3<TYPE> &axis, TYPE sinAngle, TYPE cosAngle )
	{
		const TYPE t = TYPE(1) - cosAngle;
		const TYPE tx = t * axis.x;
		const TYPE ty = t * axis.y;
		const TYPE tz = t * axis.z;
		const TYPE txy = tx * axis.y;
		const TYPE txz = tx * axis.z;
		const TYPE tyz = ty * axis.z;
		const TYPE sx = sinAngle * axis.x;
		const TYPE sy = sinAngle * axis.y;
		const TYPE sz = sinAngle * axis.z;
		data[ 0] = tx * axis.x + cosAngle;  data[ 1] = txy + sz;                data[ 2] = txz - sy;                data[ 3] = TYPE(0);
		data[ 4] = txy - sz;                data[ 5] = ty * axis.y + cosAngle;  data[ 6] = tyz + sx;                data[ 7] = TYPE(0);
		data[ 8] = txz + sy;                data[ 9] = tyz - sx;                data[10] = tz * axis.z + cosAngle;
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15] = TYPE(1);
	}
	//! Set a rotation matrix that sets [from] unit vector to [to] unit vector
	void SetRotation( const Point3<TYPE> &from, const Point3<TYPE> &to )
	{
		TYPE c = from.Dot(to);
		if ( c > TYPE(0.9999999) ) SetIdentity();
		else {
			TYPE s = cySqrt(TYPE(1) - c*c);
			Point3<TYPE> axis = from.Cross(to).GetNormalized();
			SetRotation(axis, s, c);
		}
	}
	//! Sets a translation matrix with no rotation or scale
	void SetTrans( const Point3<TYPE> &move ) { TYPE d[12]={1,0,0,0, 0,1,0,0, 0,0,1,0}; CY_MEMCOPY(TYPE,data,d,12); data[12]=move.x; data[13]=move.y; data[14]=move.z; data[15]=TYPE(1); }
	//! Adds a translation to the matrix
	void AddTrans( const Point3<TYPE> &move ) { data[12]+=move.x; data[13]+=move.y; data[14]+=move.z; }
	//! Sets the translation component of the matrix
	void SetTransComponent( const Point3<TYPE> &move ) { data[12]=move.x; data[13]=move.y; data[14]=move.z; }
	//! Set view matrix using position, target and approximate up vector
	void SetView( const Point3<TYPE> &pos, const Point3<TYPE> &target, const Point3<TYPE> &up )
	{
		Point3<TYPE> f = target - pos;
		f.Normalize();
		Point3<TYPE> s = f.Cross(up);
		s.Normalize();
		Point3<TYPE> u = s.Cross(f);
		data[ 0]=s.x; data[ 1]=u.x; data[ 2]=-f.x; data[ 3]=TYPE(0);
		data[ 4]=s.y; data[ 5]=u.y; data[ 6]=-f.y; data[ 7]=TYPE(0);
		data[ 8]=s.z; data[ 9]=u.z; data[10]=-f.z; data[11]=TYPE(0);
		data[12]= -s % pos;
		data[13]= -u % pos;
		data[14]=  f % pos;
		data[15]=TYPE(1);
	}
	//! Set matrix using normal and approximate x direction
	void SetNormal(const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,Point3<TYPE>(TYPE(0),TYPE(0),TYPE(0))); }
	//! Set a project matrix with field of view in radians
	void SetPerspective( TYPE fov, TYPE aspect, TYPE znear, TYPE zfar ) { SetPerspectiveTan(cyTan(fov*TYPE(0.5)),aspect,znear,zfar); }
	//! Set a project matrix with the tangent of the half field of view (tan_fov_2)
	void SetPerspectiveTan( TYPE tan_fov_2, TYPE aspect, TYPE znear, TYPE zfar )
	{
		const TYPE yScale = TYPE(1) / tan_fov_2;
		const TYPE xScale = yScale / aspect;
		const TYPE zdif = znear - zfar;
		TYPE d[16] = { xScale,0,0,0,  0,yScale,0,0,  0,0,(zfar+znear)/zdif,-1,  0,0,(2*zfar*znear)/zdif,0 };
		CY_MEMCOPY(TYPE,data,d,16);
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Set Row, Column, or Diagonal

	void SetRow( int row, TYPE x, TYPE y, TYPE z, TYPE w ) { data[row]=x; data[row+4]=y; data[row+8]=z; data[row+12]=w; }							//!< Sets a row of the matrix
	void SetColumn( int column, TYPE x, TYPE y, TYPE z, TYPE w ) { data[4*column]=x; data[4*column+1]=y; data[4*column+2]=z; data[4*column+3]=w; }	//!< Sets a column of the matrix
	void SetDiagonal( const TYPE &xx, const TYPE &yy, const TYPE &zz, const TYPE &ww=1 ) { data[0]=xx; data[5]=yy; data[10]=zz; data[15]=ww; }		//!< Sets the diagonal values of the matrix
	void SetDiagonal( const Point4<TYPE> &p ) { SetDiagonal( p.x, p.y, p.z, p.w ); }																//!< Sets the diagonal values of the matrix
	void SetDiagonal( const Point3<TYPE> &p ) { SetDiagonal( p.x, p.y, p.z, TYPE(1) ); }															//!< Sets the diagonal values of the matrix
	void SetDiagonal( const TYPE *values ) { SetDiagonal(values[0],values[1],values[2],values[3]); }												//!< Sets the 4 diagonal values of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Row, Column, or Diagonal

	Point4<TYPE>   GetRow   ( int row )                  const { return Point4<TYPE>( data[row], data[row+4], data[row+8], data[row+12] ); }							//!< Returns a row of the matrix
	void           GetRow   ( int row, Point4<TYPE> &p ) const { p.Set( data[row], data[row+4], data[row+8], data[row+12] ); }											//!< Returns a row of the matrix
	void           GetRow   ( int row, TYPE *values )    const { values[0]=data[row]; values[1]=data[row+4]; values[2]=data[row+8]; values[3]=data[row+12]; }			//!< Returns a row of the matrix
	Point4<TYPE>   GetColumn( int col )                  const { return Point4<TYPE>( &data[col*4] ); }																	//!< Returns a column of the matrix
	void           GetColumn( int col, Point4<TYPE> &p ) const { p.Set( &data[col*4] ); }																				//!< Returns a column of the matrix
	void           GetColumn( int col, TYPE *values )    const { values[0]=data[col*4]; values[1]=data[col*4+1]; values[2]=data[col*4+2]; values[3]=data[col*4+3]; }	//!< Returns a column of the matrix
	Point4<TYPE>   GetDiagonal()                         const { Point4<TYPE> r; GetDiagonal(r); return r; }															//!< Returns the diagonal of the matrix
	void           GetDiagonal( Point4<TYPE> &p )        const { GetDiagonal(&p.x); }																					//!< Returns the diagonal of the matrix
	void	       GetDiagonal( TYPE *values )           const { values[0]=data[0]; values[1]=data[5]; values[2]=data[10]; values[3]=data[15]; }						//!< Returns the diagonal of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Sub-matrix data

	void           GetSubMatrix ( Matrix34<TYPE> &m )    const { GetSubMatrix34(m); }																					//!< Returns the 3x4 portion of the matrix
	void           GetSubMatrix ( Matrix3 <TYPE> &m )    const { GetSubMatrix3 (m); }																					//!< Returns the 3x3 portion of the matrix
	void           GetSubMatrix ( Matrix2 <TYPE> &m )    const { GetSubMatrix2 (m); }																					//!< Returns the 2x2 portion of the matrix
	Matrix34<TYPE> GetSubMatrix34()                      const { Matrix34<TYPE> m; GetSubMatrix34(m.data); return m; }													//!< Returns the 3x4 portion of the matrix
	void           GetSubMatrix34( Matrix34<TYPE> &m )   const { GetSubMatrix34(m.data); }																				//!< Returns the 3x4 portion of the matrix
	void           GetSubMatrix34( TYPE *mdata )         const { CY_MEMCOPY(TYPE,mdata,data,3); CY_MEMCOPY(TYPE,mdata+3,data+4,3); CY_MEMCOPY(TYPE,mdata+6,data+8,3); CY_MEMCOPY(TYPE,mdata+9,data+12,3); }	//!< Returns the 3x4 portion of the matrix
	Matrix3<TYPE>  GetSubMatrix3 ()                      const { Matrix3<TYPE> m; GetSubMatrix3(m.data); return m; }													//!< Returns the 3x3 portion of the matrix
	void           GetSubMatrix3 ( Matrix3<TYPE> &m )    const { GetSubMatrix3(m.data); }																				//!< Returns the 3x3 portion of the matrix
	void           GetSubMatrix3 ( TYPE *mdata )         const { CY_MEMCOPY(TYPE,mdata,data,3); CY_MEMCOPY(TYPE,mdata+3,data+4,3); CY_MEMCOPY(TYPE,mdata+6,data+8,3); }	//!< Returns the 3x3 portion of the matrix
	Matrix2<TYPE>  GetSubMatrix2 ()                      const { Matrix2<TYPE> m; GetSubMatrix2(m.data); return m; }													//!< Returns the 2x2 portion of the matrix
	void           GetSubMatrix2 ( Matrix2<TYPE> &m )    const { GetSubMatrix2(m.data); }																				//!< Returns the 2x2 portion of the matrix
	void           GetSubMatrix2 ( TYPE *mdata )         const { CY_MEMCOPY(TYPE,mdata,data,2); CY_MEMCOPY(TYPE,mdata+2,data+4,2); }									//!< Returns the 2x2 portion of the matrix
	Point3<TYPE>   GetTrans()                            const { Point3<TYPE> p; GetTrans(p); return p; }																//! Returns the translation component of the matrix
	void           GetTrans( Point3<TYPE> &p )           const { GetTrans(&p.x); }																						//! Returns the translation component of the matrix
	void           GetTrans( TYPE *trans )               const { CY_MEMCOPY(TYPE,trans,data+12,3); }																	//! Returns the translation component of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Comparison Operators

	bool operator == ( const Matrix4 &right ) const { for ( int i=0; i<16; i++ ) if ( data[i] != right.data[i] ) return false; return true; } //!< compare equal
	bool operator != ( const Matrix4 &right ) const { for ( int i=0; i<16; i++ ) if ( data[i] != right.data[i] ) return true; return false; } //!< compare not equal


	//////////////////////////////////////////////////////////////////////////
	//!@name Access Operators

	TYPE&       operator () ( int row, int column )       { return data[ column * 4 + row ]; }	//!< subscript operator
	const TYPE& operator () ( int row, int column ) const { return data[ column * 4 + row ]; }	//!< constant subscript operator
	TYPE&       operator [] ( int i )       { return data[i]; }	//!< subscript operator
	const TYPE& operator [] ( int i ) const { return data[i]; }	//!< constant subscript operator


	//////////////////////////////////////////////////////////////////////////
	//!@name Unary and Binary Operators

	// Unary operators
	Matrix4 operator - () const { Matrix4 r; for (int i=0; i<16; i++) r.data[i]=-data[i]; return r; }	//!< negative matrix

	// Binary operators
	Matrix4 operator * ( const TYPE    &value ) const { Matrix4 r; for (int i=0; i<16; ++i) r.data[i] = data[i] * value;         return r; }	//!< multiply matrix by a value
	Matrix4 operator / ( const TYPE    &value ) const { Matrix4 r; for (int i=0; i<16; ++i) r.data[i] = data[i] / value;         return r; }	//!< divide matrix by a value;
	Matrix4 operator + ( const Matrix4 &right ) const { Matrix4 r; for (int i=0; i<16; i++) r.data[i] = data[i] + right.data[i]; return r; }	//!< add two Matrices
	Matrix4 operator - ( const Matrix4 &right ) const { Matrix4 r; for (int i=0; i<16; i++) r.data[i] = data[i] - right.data[i]; return r; }	//!< subtract one Matrix4 from another
	Matrix4 operator * ( const Matrix4 &right ) const	//!< multiply a matrix with another
	{
		Matrix4 r;
		TYPE *rd = r.data;
		for ( int i=0; i<16; i+=4, rd+=4 ) {
			TYPE a[4], b[4], c[4], d[4];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) a[j] = data[   j] * right.data[i  ];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) b[j] = data[ 4+j] * right.data[i+1];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) c[j] = data[ 8+j] * right.data[i+2];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) d[j] = data[12+j] * right.data[i+3];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) rd[j] = a[j] + b[j] + c[j] + d[j];
		}
		return r;
	}
	Matrix4 operator * ( const Matrix34<TYPE> &right ) const	//!< multiply a matrix with another
	{
		TYPE a[4], b[4], c[4];
		Matrix4 r;
		TYPE *rd = r.data;
		for ( int i=0; i<9; i+=3, rd+=4 ) {
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) a[k] = data[  k] * right.data[i  ];
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) b[k] = data[4+k] * right.data[i+1];
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) c[k] = data[8+k] * right.data[i+2];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) rd[j] = a[j] + b[j] + c[j];
		}
		_CY_IVDEP_FOR ( int k=0; k<4; ++k ) a[k] = data[  k] * right.data[ 9];
		_CY_IVDEP_FOR ( int k=0; k<4; ++k ) b[k] = data[4+k] * right.data[10];
		_CY_IVDEP_FOR ( int k=0; k<4; ++k ) c[k] = data[8+k] * right.data[11];
		_CY_IVDEP_FOR ( int j=0; j<4; ++j ) rd[j] = (a[j] + b[j]) + (c[j] + data[12+j]);
		return r;
	}
	Matrix4 operator * ( const Matrix3<TYPE> &right ) const	//!< multiply a matrix with another
	{
		TYPE a[4], b[4], c[4];
		Matrix4 r;
		TYPE *rd = r.data;
		for ( int i=0; i<9; i+=3, rd+=4 ) {
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) a[k] = data[  k] * right.data[i  ];
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) b[k] = data[4+k] * right.data[i+1];
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) c[k] = data[8+k] * right.data[i+2];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) rd[j] = a[j] + b[j] + c[j];
		}
		CY_MEMCOPY(TYPE,r.data+12,data+12,4);
		return r;
	}
	Point4<TYPE> operator * ( const Point3<TYPE>& p ) const 
	{
		//return Point4<TYPE>(	p.x*data[0] + p.y*data[4] + p.z*data[ 8] + data[12], 
		//						p.x*data[1] + p.y*data[5] + p.z*data[ 9] + data[13],
		//						p.x*data[2] + p.y*data[6] + p.z*data[10] + data[14],
		//						p.x*data[3] + p.y*data[7] + p.z*data[11] + data[15] );
		TYPE a[4], b[4], c[4];
		Point4<TYPE> rr;
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) a[i] = p[0] * data[   i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) b[i] = p[1] * data[ 4+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) c[i] = p[2] * data[ 8+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) rr[i] = a[i] + b[i] + c[i] + data[12+i];	
		return rr;
	}
	Point4<TYPE> operator * ( const Point4<TYPE>& p ) const 
	{
		//return Point4<TYPE>(	p.x*data[0] + p.y*data[4] + p.z*data[ 8] + p.w*data[12],
		//						p.x*data[1] + p.y*data[5] + p.z*data[ 9] + p.w*data[13],
		//						p.x*data[2] + p.y*data[6] + p.z*data[10] + p.w*data[14],
		//						p.x*data[3] + p.y*data[7] + p.z*data[11] + p.w*data[15] );
		TYPE a[8], b[8];
		const TYPE *pd = p.Data();
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) a[  i] = pd[0] * data[   i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) a[4+i] = pd[1] * data[ 4+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) b[  i] = pd[2] * data[ 8+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) b[4+i] = pd[3] * data[12+i];
		_CY_IVDEP_FOR ( int i=0; i<8; ++i ) a[i] += b[i];
		Point4<TYPE> rr;
		TYPE *rd = rr.Data();
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) rd[i] = a[i] + a[4+i];
		return rr;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Assignment Operators

	const Matrix4& operator  = ( const Matrix4 &right ) { CY_MEMCOPY(TYPE,data,right.data,16); return *this; }	
	const Matrix4& operator += ( const Matrix4 &right ) { for (int i=0; i<16; i++) data[i] += right.data[i]; return *this; }	//!< add two Matrices modify this
	const Matrix4& operator -= ( const Matrix4 &right ) { for (int i=0; i<16; i++) data[i] -= right.data[i]; return *this; }	//!< subtract one Matrix4 from another matrix and modify this matrix
	const Matrix4& operator *= ( const Matrix4 &right )        { *this = operator*(right); return *this; }						//!< multiply a matrix with another matrix and modify this matrix
	const Matrix4& operator *= ( const Matrix34<TYPE> &right ) { *this = operator*(right); return *this; }						//!< multiply a matrix with another matrix and modify this matrix
	const Matrix4& operator *= ( const Matrix3<TYPE>  &right ) { *this = operator*(right); return *this; }						//!< multiply a matrix with another matrix and modify this matrix
	const Matrix4& operator *= ( const TYPE    &value ) { for (int i=0; i<16; i++) data[i] *= value;         return *this; }	//!< multiply a matrix with a value modify this matrix
	const Matrix4& operator /= ( const TYPE    &value ) { for (int i=0; i<16; i++) data[i] /= value;         return *this; }	//!< divide the matrix by a value modify the this matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Other Methods

	void Transpose()															//!< Transpose this matrix
	{
		for (int i = 1; i < 4; i++) {
			for (int j = 0; j < i; j++) {
				TYPE temp = data[i * 4 + j];
				data[i * 4 + j] = data[j * 4 + i];
				data[j * 4 + i] = temp;
			}
		}
	}
	void GetTranspose( Matrix4 &m ) const										//!< return Transpose of this matrix
	{
		m.data[ 0] = data[0];   m.data[ 1] = data[4];   m.data[ 2] = data[ 8];  m.data[ 3] = data[12];
		m.data[ 4] = data[1];   m.data[ 5] = data[5];   m.data[ 6] = data[ 9];  m.data[ 7] = data[13];
		m.data[ 8] = data[2];   m.data[ 9] = data[6];   m.data[10] = data[10];  m.data[11] = data[14];
		m.data[12] = data[3];   m.data[13] = data[7];   m.data[14] = data[11];  m.data[15] = data[15];
	}
	Matrix4 GetTranspose() const { Matrix4 t; GetTranspose(t); return t; }	//!< return Transpose of this matrix

	//! Multiply the give vector with the transpose of the matrix
	Point4<TYPE> TransposeMult( const Point3<TYPE>& p ) const { return TransposeMult( Point4<TYPE>(p.x, p.y, p.z, TYPE(1)) ); }

	//! Multiply the give vector with the transpose of the matrix
	Point4<TYPE> TransposeMult( const Point4<TYPE>& p ) const 
	{
		//return Point4<TYPE>(	p.x*data[ 0] + p.y*data[ 1] + p.z*data[ 2] + p.w*data[ 3],
		//						p.x*data[ 4] + p.y*data[ 5] + p.z*data[ 6] + p.w*data[ 7],
		//						p.x*data[ 8] + p.y*data[ 9] + p.z*data[10] + p.w*data[11],
		//						p.x*data[12] + p.y*data[13] + p.z*data[14] + p.w*data[15] );
		TYPE a[4], b[4], c[4], d[4];
		const TYPE *pd = p.Data();
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) a[i] = pd[i] * data[   i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) b[i] = pd[i] * data[ 4+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) c[i] = pd[i] * data[ 8+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) d[i] = pd[i] * data[12+i];
		Point4<TYPE> rr;
		rr.x = a[0] + a[1] + a[2] + a[3];
		rr.y = b[0] + b[1] + b[2] + b[3];
		rr.z = c[0] + c[1] + c[2] + c[3];
		rr.w = d[0] + d[1] + d[2] + d[3];
		return rr;
	}

	TYPE GetDeterminant() const	//!< Get the determinant of this matrix
	{
		// 0 (                      5 ( 10 15 - 11 14) + 6 ( 11 13 -  9 15) + 7 (  9 14 - 10 13)) +  
		// 1 ( 4 ( 11 14 - 10 15) +                      6 (  8 15 - 11 12) + 7 ( 10 12 -  8 14)) +  
		// 2 ( 4 (  9 15 - 11 13) + 5 ( 11 12 -  8 15) +                      7 (  8 13 -  9 12)) + 
		// 3 ( 4 ( 10 13 -  9 14) + 5 (  8 14 - 10 12) + 6 (  9 12 -  8 13)) 

		const TYPE data_11_14__10_15 = data[11] * data[14] - data[10] * data[15];
		const TYPE data__9_15__11_13 = data[ 9] * data[15] - data[11] * data[13];
		const TYPE data_10_13___9_14 = data[10] * data[13] - data[ 9] * data[14];
		const TYPE data_11_12___8_15 = data[11] * data[12] - data[ 8] * data[15];
		const TYPE data__8_14__10_12 = data[ 8] * data[14] - data[10] * data[12];
		const TYPE data__9_12___8_13 = data[ 9] * data[12] - data[ 8] * data[13];
		return data[0] * ( data[5] * (-data_11_14__10_15 ) + data[6] * (-data__9_15__11_13 ) + data[7] * (-data_10_13___9_14 ) ) +  
		       data[1] * ( data[4] * ( data_11_14__10_15 ) + data[6] * (-data_11_12___8_15 ) + data[7] * (-data__8_14__10_12 ) ) +  
		       data[2] * ( data[4] * ( data__9_15__11_13 ) + data[5] * ( data_11_12___8_15 ) + data[7] * (-data__9_12___8_13 ) ) + 
		       data[3] * ( data[4] * ( data_10_13___9_14 ) + data[5] * ( data__8_14__10_12 ) + data[6] * ( data__9_12___8_13 ) );
	}

	void Invert() { Matrix4 inv; GetInverse(inv); *this=inv; }					//!< Invert this matrix
	void GetInverse( Matrix4 &inverse ) const									//!< Get the inverse of this matrix
	{
		//                       5 ( 10 15 - 11 14 ) + 6 ( 11 13 -  9 15 ) + 7 (  9 14 - 10 13 )
		//                       1 ( 11 14 - 10 15 ) + 2 (  9 15 - 11 13 ) + 3 ( 10 13 -  9 14 ) 
		//                       1 (  6 15 -  7 14 ) + 2 (  7 13 -  5 15 ) + 3 (  5 14 -  6 13 )
		//                       1 (  7 10 -  6 11 ) + 2 (  5 11 -  7  9 ) + 3 (  6  9 -  5 10 )
		//					 					   						 					   
		// 4 ( 11 14 - 10 15 ) +                       6 (  8 15 - 11 12 ) + 7 ( 10 12 -  8 14 )
		// 0 ( 10 15 - 11 14 ) +                       2 ( 11 12 -  8 15 ) + 3 (  8 14 - 10 12 )
		// 0 (  7 14 -  6 15 ) +                       2 (  4 15 -  7 12 ) + 3 (  6 12 -  4 14 )      / det
		// 0 (  6 11 -  7 10 ) +                       2 (  8  7 -  4 11 ) + 3 (  4 10 -  6  8 )
		//					 					   						 					   
		// 4 (  9 15 - 11 13 ) + 5 ( 11 12 -  8 15 ) +                       7 (  8 13 -  9 12 ) 
		// 0 ( 11 13 -  9 15 ) + 1 (  8 15 - 11 12 ) +                       3 (  9 12 -  8 13 )
		// 0 (  5 15 -  7 13 ) + 1 (  7 12 -  4 15 ) +                       3 (  4 13 -  5 12 )
		// 0 (  7  9 -  5 11 ) + 1 (  4 11 -  7  8 ) +                       3 (  5  8 -  4  9 )
		//					 					   						 
		// 4 ( 10 13 -  9 14 ) + 5 (  8 14 - 10 12 ) + 6 (  9 12 -  8 13 )
		// 0 (  9 14 - 10 13 ) + 1 ( 10 12 -  8 14 ) + 2 (  8 13 -  9 12 )
		// 0 (  6 13 -  5 14 ) + 1 (  4 14 -  6 12 ) + 2 (  5 12 -  4 13 )
		// 0 (  5 10 -  6  9 ) + 1 (  6  8 -  4 10 ) + 2 (  4  9 -  5  8 )

		const TYPE data_11_14__10_15 = data[11] * data[14] - data[10] * data[15];
		const TYPE data_10_15__11_14 = data[10] * data[15] - data[11] * data[14];
		const TYPE data__7_14___6_15 = data[ 7] * data[14] - data[ 6] * data[15];
		const TYPE data__6_11___7_10 = data[ 6] * data[11] - data[ 7] * data[10];

		const TYPE data__9_15__11_13 = data[ 9] * data[15] - data[11] * data[13];
		const TYPE data_11_13___9_15 = data[11] * data[13] - data[ 9] * data[15];
		const TYPE data__5_15___7_13 = data[ 5] * data[15] - data[ 7] * data[13];
		const TYPE data__7__9___5_11 = data[ 7] * data[ 9] - data[ 5] * data[11];
		
		const TYPE data_10_13___9_14 = data[10] * data[13] - data[ 9] * data[14];
		const TYPE data__9_14__10_13 = data[ 9] * data[14] - data[10] * data[13];
		const TYPE data__6_13___5_14 = data[ 6] * data[13] - data[ 5] * data[14];
		const TYPE data__5_10___6__9 = data[ 5] * data[10] - data[ 6] * data[ 9];
		
		const TYPE data_11_12___8_15 = data[11] * data[12] - data[ 8] * data[15];
		const TYPE data__8_15__11_12 = data[ 8] * data[15] - data[11] * data[12];
		const TYPE data__7_12___4_15 = data[ 7] * data[12] - data[ 4] * data[15];
		const TYPE data__4_11___7__8 = data[ 4] * data[11] - data[ 7] * data[ 8];
		
		const TYPE data__8_14__10_12 = data[ 8] * data[14] - data[10] * data[12];
		const TYPE data_10_12___8_14 = data[10] * data[12] - data[ 8] * data[14];
		const TYPE data__4_14___6_12 = data[ 4] * data[14] - data[ 6] * data[12];
		const TYPE data__6__8___4_10 = data[ 6] * data[ 8] - data[ 4] * data[10];
		
		const TYPE data__9_12___8_13 = data[ 9] * data[12] - data[ 8] * data[13];
		const TYPE data__8_13___9_12 = data[ 8] * data[13] - data[ 9] * data[12];
		const TYPE data__5_12___4_13 = data[ 5] * data[12] - data[ 4] * data[13];
		const TYPE data__4__9___5__8 = data[ 4] * data[ 9] - data[ 5] * data[ 8];

		inverse.data[ 0] = data[5] * (-data_11_14__10_15) + data[6] * (-data__9_15__11_13) + data[7] * (-data_10_13___9_14);
		inverse.data[ 1] = data[1] * (-data_10_15__11_14) + data[2] * (-data_11_13___9_15) + data[3] * (-data__9_14__10_13);
		inverse.data[ 2] = data[1] * (-data__7_14___6_15) + data[2] * (-data__5_15___7_13) + data[3] * (-data__6_13___5_14);
		inverse.data[ 3] = data[1] * (-data__6_11___7_10) + data[2] * (-data__7__9___5_11) + data[3] * (-data__5_10___6__9);
		
		inverse.data[ 4] = data[4] * ( data_11_14__10_15) + data[6] * (-data_11_12___8_15) + data[7] * (-data__8_14__10_12);
		inverse.data[ 5] = data[0] * ( data_10_15__11_14) + data[2] * (-data__8_15__11_12) + data[3] * (-data_10_12___8_14);
		inverse.data[ 6] = data[0] * ( data__7_14___6_15) + data[2] * (-data__7_12___4_15) + data[3] * (-data__4_14___6_12);
		inverse.data[ 7] = data[0] * ( data__6_11___7_10) + data[2] * (-data__4_11___7__8) + data[3] * (-data__6__8___4_10);
		
		inverse.data[ 8] = data[4] * ( data__9_15__11_13) + data[5] * ( data_11_12___8_15) + data[7] * (-data__9_12___8_13);
		inverse.data[ 9] = data[0] * ( data_11_13___9_15) + data[1] * ( data__8_15__11_12) + data[3] * (-data__8_13___9_12);
		inverse.data[10] = data[0] * ( data__5_15___7_13) + data[1] * ( data__7_12___4_15) + data[3] * (-data__5_12___4_13);
		inverse.data[11] = data[0] * ( data__7__9___5_11) + data[1] * ( data__4_11___7__8) + data[3] * (-data__4__9___5__8);

		inverse.data[12] = data[4] * ( data_10_13___9_14) + data[5] * ( data__8_14__10_12) + data[6] * ( data__9_12___8_13);
		inverse.data[13] = data[0] * ( data__9_14__10_13) + data[1] * ( data_10_12___8_14) + data[2] * ( data__8_13___9_12);
		inverse.data[14] = data[0] * ( data__6_13___5_14) + data[1] * ( data__4_14___6_12) + data[2] * ( data__5_12___4_13);
		inverse.data[15] = data[0] * ( data__5_10___6__9) + data[1] * ( data__6__8___4_10) + data[2] * ( data__4__9___5__8);

		const TYPE det = data[0] * inverse.data[0] + data[1] * inverse.data[4] + data[2] * inverse.data[8] + data[3] * inverse.data[12];
		inverse /= det;
	}
	Matrix4 GetInverse() const { Matrix4 inv; GetInverse(inv); return inv; }	//!< Get the inverse of this matrix

	//! Orthogonalizes the matrix and removes the scale component, preserving the x direction
	void OrthogonalizeX()
	{
		Point3<TYPE> &px = *((Point3<TYPE>*)&data[0]);
		Point3<TYPE> &py = *((Point3<TYPE>*)&data[4]);
		Point3<TYPE> &pz = *((Point3<TYPE>*)&data[8]);
		px.Normalize();
		py -= px * (py%px);
		py.Normalize();
		pz -= px * (pz%px);
		pz -= py * (pz%py);
		pz.Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the y direction
	void OrthogonalizeY()
	{
		Point3<TYPE> &px = *((Point3<TYPE>*)&data[0]);
		Point3<TYPE> &py = *((Point3<TYPE>*)&data[4]);
		Point3<TYPE> &pz = *((Point3<TYPE>*)&data[8]);
		py.Normalize();
		px -= py * (px%py);
		px.Normalize();
		pz -= py * (pz%py);
		pz -= px * (pz%px);
		pz.Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the z direction
	void OrthogonalizeZ()
	{
		Point3<TYPE> &px = *((Point3<TYPE>*)&data[0]);
		Point3<TYPE> &py = *((Point3<TYPE>*)&data[4]);
		Point3<TYPE> &pz = *((Point3<TYPE>*)&data[8]);
		pz.Normalize();
		px -= pz * (px%pz);
		px.Normalize();
		py -= pz * (py%pz);
		py -= px * (py%px);
		py.Normalize();
	}

	//! Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const
	{
		return cyAbs(data[ 0]-TYPE(1)) < tollerance && cyAbs(data[ 1])         < tollerance && cyAbs(data[ 2])         < tollerance && cyAbs(data[ 3])         < tollerance && 
			   cyAbs(data[ 4])         < tollerance && cyAbs(data[ 5]-TYPE(1)) < tollerance && cyAbs(data[ 6])         < tollerance && cyAbs(data[ 7])         < tollerance &&
			   cyAbs(data[ 8])         < tollerance && cyAbs(data[ 9])         < tollerance && cyAbs(data[10]-TYPE(1)) < tollerance && cyAbs(data[11])         < tollerance &&
			   cyAbs(data[12])         < tollerance && cyAbs(data[13])         < tollerance && cyAbs(data[14])         < tollerance && cyAbs(data[15]-TYPE(1)) < tollerance;
	}

	//! Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( TYPE tollerance=TYPE(0.001) ) const
	{
		return cyAbs(data[ 1] - data[ 4]) < tollerance && 
			   cyAbs(data[ 2] - data[ 8]) < tollerance &&
			   cyAbs(data[ 3] - data[12]) < tollerance &&
			   cyAbs(data[ 6] - data[ 9]) < tollerance &&
			   cyAbs(data[ 7] - data[13]) < tollerance &&
			   cyAbs(data[11] - data[14]) < tollerance;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Static Methods

	//! Returns an identity matrix
	static Matrix4 MatrixIdentity() { Matrix4 m; m.SetIdentity(); return m; }
	//! Returns a view matrix using position, target and approximate up vector
	static Matrix4 MatrixView( const Point3<TYPE> &pos, const Point3<TYPE> &target, Point3<TYPE> &up ) { Matrix4 m; m.SetView(pos,target,up); return m; }
	//! Returns a matrix using normal, and approximate x direction
	static Matrix4 MatrixNormal( const Point3<TYPE> &normal, Point3<TYPE> &dir ) { Matrix4 m; m.SetNormal(normal,dir); return m; }
	//! Returns a rotation matrix around x axis by angle in radians
	static Matrix4 MatrixRotationX( TYPE angle ) { Matrix4 m; m.SetRotationX(angle); return m; }
	//! Returns a rotation matrix around y axis by angle in radians
	static Matrix4 MatrixRotationY( TYPE angle ) { Matrix4 m; m.SetRotationY(angle); return m; }
	//! Returns a rotation matrix around z axis by angle in radians
	static Matrix4 MatrixRotationZ( TYPE angle ) { Matrix4 m; m.SetRotationZ(angle); return m; }
	//! Returns a rotation matrix about the given axis by angle in radians
	static Matrix4 MatrixRotation( const Point3<TYPE> &axis, TYPE angle ) { Matrix4 m; m.SetRotation(axis,angle); return m; }
	//! Returns a rotation matrix that sets [from] unit vector to [to] unit vector
	static Matrix4 MatrixRotation( const Point3<TYPE> &from, const Point3<TYPE> &to ) { Matrix4 m; m.SetRotation(from,to); return m; }
	//! Returns a uniform scale matrix
	static Matrix4 MatrixScale( const TYPE &uniformScale ) { Matrix4 m; m.SetScale(uniformScale); return m; }
	//! Returns a scale matrix
	static Matrix4 MatrixScale( const TYPE &scaleX, const TYPE &scaleY, const TYPE &scaleZ ) { Matrix4 m; m.SetScale(scaleX,scaleY,scaleZ); return m; }
	//! Returns a scale matrix
	static Matrix4 MatrixScale( const Point3<TYPE> &scale ) { Matrix4 m; m.SetScale(scale); return m; }
	//! Returns a translation matrix with no rotation or scale
	static Matrix4 MatrixTrans( const Point3<TYPE> &move ) { Matrix4 m; m.SetTrans(move); return m; }
	//! Returns a project matrix with field of view in radians
	static Matrix4 MatrixPerspective( TYPE fov, TYPE aspect, TYPE znear, TYPE zfar ) { Matrix4 m; m.SetPerspective(fov,aspect,znear,zfar); return m; }
	//! Returns a project matrix with the tangent of the half field of view (tan_fov_2)
	static Matrix4 MatrixPerspectiveTan( TYPE tan_fov_2, TYPE aspect, TYPE znear, TYPE zfar ) { Matrix4 m; m.SetPerspectiveTan(tan_fov_2,aspect,znear,zfar); return m; }


	//////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

template<typename TYPE> inline Matrix2<TYPE> operator & ( const Point2<TYPE> &v0, const Point2<TYPE> &v1 ) { Matrix2<TYPE> r; r.SetTensorProduct(v0,v1); return r; }	//!< tensor product (outer product) of two vectors
template<typename TYPE> inline Matrix3<TYPE> operator & ( const Point3<TYPE> &v0, const Point3<TYPE> &v1 ) { Matrix3<TYPE> r; r.SetTensorProduct(v0,v1); return r; }	//!< tensor product (outer product) of two vectors
template<typename TYPE> inline Matrix4<TYPE> operator & ( const Point4<TYPE> &v0, const Point4<TYPE> &v1 ) { Matrix4<TYPE> r; r.SetTensorProduct(v0,v1); return r; }	//!< tensor product (outer product) of two vectors

//-------------------------------------------------------------------------------

// Definitions of the conversion constructors
template <typename TYPE>  Matrix2 <TYPE>::Matrix2 ( const Matrix3 <TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,2); CY_MEMCOPY(TYPE,data+2,m.data+3,2); }
template <typename TYPE>  Matrix2 <TYPE>::Matrix2 ( const Matrix34<TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,2); CY_MEMCOPY(TYPE,data+2,m.data+3,2); }
template <typename TYPE>  Matrix2 <TYPE>::Matrix2 ( const Matrix4 <TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,2); CY_MEMCOPY(TYPE,data+2,m.data+4,2); }
template <typename TYPE>  Matrix3 <TYPE>::Matrix3 ( const Matrix34<TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,9); }
template <typename TYPE>  Matrix3 <TYPE>::Matrix3 ( const Matrix4 <TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,3); CY_MEMCOPY(TYPE,data+3,m.data+4,3); CY_MEMCOPY(TYPE,data+6,m.data+8,3); }
template <typename TYPE>  Matrix34<TYPE>::Matrix34( const Matrix4 <TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,3); CY_MEMCOPY(TYPE,data+3,m.data+4,3); CY_MEMCOPY(TYPE,data+6,m.data+8,3); CY_MEMCOPY(TYPE,data+9,m.data+12,3); }

//-------------------------------------------------------------------------------

typedef Matrix2 <float>  Matrix2f;	//!< Single precision (float) 2x2 Matrix class
typedef Matrix3 <float>  Matrix3f;	//!< Single precision (float) 3x3 Matrix class
typedef Matrix34<float>  Matrix34f;	//!< Single precision (float) 3x4 Matrix class
typedef Matrix4 <float>  Matrix4f;	//!< Single precision (float) 4x4 Matrix class

typedef Matrix2 <double> Matrix2d;	//!< Double precision (double) 2x2 Matrix class
typedef Matrix3 <double> Matrix3d;	//!< Double precision (double) 3x3 Matrix class
typedef Matrix34<double> Matrix34d;	//!< Double precision (double) 3x4 Matrix class
typedef Matrix4 <double> Matrix4d;	//!< Double precision (double) 4x4 Matrix class

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Matrix2f  cyMatrix2f;	//!< Single precision (float) 2x2 Matrix class
typedef cy::Matrix3f  cyMatrix3f;	//!< Single precision (float) 3x3 Matrix class
typedef cy::Matrix34f cyMatrix34f;	//!< Single precision (float) 3x4 Matrix class
typedef cy::Matrix4f  cyMatrix4f;	//!< Single precision (float) 4x4 Matrix class

typedef cy::Matrix2d  cyMatrix2d;	//!< Double precision (double) 2x2 Matrix class
typedef cy::Matrix3d  cyMatrix3d;	//!< Double precision (double) 3x3 Matrix class
typedef cy::Matrix34d cyMatrix34d;	//!< Double precision (double) 3x4 Matrix class
typedef cy::Matrix4d  cyMatrix4d;	//!< Double precision (double) 4x4 Matrix class

//-------------------------------------------------------------------------------

#endif
