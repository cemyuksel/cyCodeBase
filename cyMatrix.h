// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
/// \file		cyMatrix.h 
/// \author		Cem Yuksel
/// \brief		2x2, 3x3, 3x4, and 4x4 matrix classes
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

#ifndef _CY_MATRIX2_H_INCLUDED_
#define _CY_MATRIX2_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyPoint.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

/// 2x2 matrix class.
/// Its data stores 4-value array of column-major matrix elements.
/// You can use Matrix2 with Point2<TYPE> to transform 2D points.
/// Both post-multiplication and pre-multiplication are supported.

template <typename TYPE>
class Matrix2
{
	
	friend Matrix2 operator * ( const TYPE value, const Matrix2 &right ) { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = value * right.data[i]; return buffer; }	///< multiple matrix by a value
	friend Matrix2 operator / ( const TYPE value, const Matrix2 &right ) { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = value / right.data[i]; return buffer; }	///< divide a value by a matrix
	friend Matrix2 Inverse( Matrix2 &m ) { return m.GetInverse(); }	///< return the inverse of the matrix

public:

	/// Elements of the matrix are column-major:
	/// | 0  2 |
	/// | 1  3 |
	TYPE data[4];

	//////////////////////////////////////////////////////////////////////////
	///@name Constructors

	Matrix2() {}																			///< Default constructor
	Matrix2( const Matrix2 &matrix ) { for ( int i=0; i<4; i++ ) data[i]=matrix.data[i]; }	///< Copy constructor
	template <typename T> explicit Matrix2<TYPE>( const Matrix2<T> &matrix ) { for ( int i=0; i<4; i++ ) data[i]=(TYPE)matrix.data[i]; }		///< Copy constructor for different types
	explicit Matrix2( const TYPE *array ) { Set(array); }									///< Initialize the matrix using an array of 4 values
	explicit Matrix2( TYPE v ) { SetScaledIdentity(v); }									///< Initialize the matrix as identity scaled by v
	explicit Matrix2( const Point2<TYPE> &x, const Point2<TYPE> &y ) { Set(x,y); }			///< Initialize the matrix using two vectors as columns


	//////////////////////////////////////////////////////////////////////////
	///@name Set & Get Methods

	/// Set all the values as zero
	void Zero() { for ( int i=0; i<4; i++ ) data[i] = TYPE(0); }
	/// Set Matrix using an array of 4 values
	void Set( const TYPE *array ) { for ( int i=0; i<4; i++ ) data[i] = array[i]; } 
	/// Set Matrix using two vectors as columns
	void Set( const Point2<TYPE> &x, const Point2<TYPE> &y ) { x.GetValue( &data[0] ); y.GetValue( &data[2] ); }
	/// Set a row of the matrix
	void SetRow( int row, TYPE x, TYPE y ) { data[row]=x; data[row+2]=y; }
	/// Set a column of the matrix
	void SetColumn( int column, TYPE x, TYPE y ) { data[2*column]=x; data[2*column+1]=y; }
	/// Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	/// Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity(TYPE v) { data[0]=v; data[1]=TYPE(0); data[2]=TYPE(0); data[3]=v; }
	/// Set a rotation matrix by angle
	void SetRotation( TYPE angle ) { SetRotation( sin(angle), cos(angle) ); }
	/// Set a rotation matrix by cos and sin of angle
	void SetRotation( TYPE sinAngle, TYPE cosAngle ) { data[0]=cosAngle; data[1]=-sinAngle; data[2]=sinAngle; data[3]=cosAngle; }
	/// Sets a uniform scale matrix
	void SetScale( TYPE uniformScale ) { SetScale(uniformScale,uniformScale); }
	/// Sets a scale matrix
	void SetScale( TYPE scaleX, TYPE scaleY ) { data[0]=scaleX; data[1]=TYPE(0); data[2]=scaleY; data[3]=TYPE(0); }
	/// Sets a scale matrix
	void SetScale( const Point2<TYPE> &scale ) { SetScale(scale.x,scale.y); }
	/// Removes the scale component of the matrix
	void SetNoScale() { Point2<TYPE> *p = (Point2<TYPE>*)data; p[0].Normalize(); p[1].Normalize(); }
	/// Sets the matrix as the tensor product (outer product) of two vectors
	void SetTensorProduct( const Point2<TYPE> &v0, const Point2<TYPE> &v1 )
	{
		data[0]=v0.x*v1.x;  data[2]=v0.x*v1.y;
		data[1]=v0.y*v1.x;  data[3]=v0.y*v1.y;
	}

	// Get Row and Column
	Point2<TYPE> GetRow   ( int row )                  const { return Point2<TYPE>( data[row], data[row+2] ); }
	void         GetRow   ( int row, Point2<TYPE> &p ) const { p.Set( data[row], data[row+1] ); }
	void         GetRow   ( int row, TYPE *array )     const { array[0]=data[row]; array[1]=data[row+2]; }
	Point2<TYPE> GetColumn( int col )                  const { return Point2<TYPE>( &data[col*2] ); }
	void         GetColumn( int col, Point2<TYPE> &p ) const { p.Set( &data[col*2] ); }
	void         GetColumn( int col, TYPE *array )     const { array[0]=data[col*2]; array[1]=data[col*2+1]; }
	// These methods return the diagonal component of the matrix
	Point2<TYPE>  GetDiagonal()                        const { return Point2<TYPE>( data[0], data[3] ); }
	void          GetDiagonal( Point2<TYPE> &p )       const { p.Set( data[0], data[3] ); }
	void	      GetDiagonal( TYPE *array )           const { array[0]=data[0]; array[1]=data[3]; }

	//////////////////////////////////////////////////////////////////////////
	///@name Overloaded Operators

	/// assign matrix
	const Matrix2& operator = ( const Matrix2 &right ) { for (int i = 0; i < 4; i++) data[i]=right.data[i]; return *this; }	

	// Overloaded comparison operators 
	bool operator == ( const Matrix2 &right ) const { for ( int i=0; i<4; i++ ) if ( data[i] != right.data[i] ) return false; return true; } ///< compare equal
	bool operator != ( const Matrix2 &right ) const { return ! ( *this == right ); } ///< compare not equal

	// Overloaded subscript operators
	TYPE&       operator () ( int row, int column )       { return data[ column * 2 + row ]; }	///< subscript operator
	const TYPE& operator () ( int row, int column ) const { return data[ column * 2 + row ]; }	///< constant subscript operator
	TYPE&       operator [] ( int i )       { return data[i]; }	///< subscript operator
	const TYPE& operator [] ( int i ) const { return data[i]; }	///< constant subscript operator
	
	// Unary operators
	Matrix2 operator - () const { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i]=- data[i]; return buffer; }	///< negative matrix

	// Binary operators
	Matrix2      operator * ( const TYPE value )      const { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = data[i] * value; return buffer; }	///< multiple matrix by a value
	Matrix2      operator / ( const TYPE value )      const { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = data[i] / value; return buffer; }	///< divide matrix by a value;
	Point2<TYPE> operator * ( const Point2<TYPE> &p ) const { return Point2<TYPE>( p.x * data[0] + p.y * data[2], p.x * data[1] + p.y * data[3] ); }
	Matrix2      operator + ( const Matrix2 &right  ) const { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = data[i] + right.data[i]; return buffer; }	///< add two Matrices
	Matrix2      operator - ( const Matrix2 &right  ) const { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = data[i] - right.data[i]; return buffer; }	///< subtract one Matrix2 from an other
	Matrix2      operator * ( const Matrix2 &right  ) const	///< multiply a matrix with an other
	{
		Matrix2 buffer;
		buffer[0] = data[0] * right.data[0] + data[2] * right.data[1];
		buffer[1] = data[1] * right.data[0] + data[3] * right.data[1];
		buffer[2] = data[0] * right.data[2] + data[2] * right.data[3];
		buffer[3] = data[1] * right.data[2] + data[3] * right.data[3];
		return buffer;
	}

	// Assignment operators
	void operator *= ( const TYPE value )     { for (int i=0; i<4; i++) data[i] *= value; }			///< multiply a matrix with a value modify this matrix
	void operator /= ( const TYPE value )     { for (int i=0; i<4; i++) data[i] /= value; }			///< divide the matrix by a value modify the this matrix
	void operator += ( const Matrix2 &right ) { for (int i=0; i<4; i++) data[i] += right.data[i]; }	///< add two Matrices modify this
	void operator -= ( const Matrix2 &right ) { for (int i=0; i<4; i++) data[i] -= right.data[i]; }	///< subtract one Matrix2 from an other modify this matrix
	void operator *= ( const Matrix2 &right )	///< multiply a matrix with an other modify this matrix
	{
		Matrix2 buffer;
		buffer[0] = data[0] * right.data[0] + data[2] * right.data[1];
		buffer[1] = data[1] * right.data[0] + data[3] * right.data[1];
		buffer[2] = data[0] * right.data[2] + data[2] * right.data[3];
		buffer[3] = data[1] * right.data[2] + data[3] * right.data[3];
		*this = buffer;
	}

	//////////////////////////////////////////////////////////////////////////
	///@name Other Public Methods

	void Transpose() { TYPE tmp=data[0]; data[0]=data[3]; data[3]=tmp; }	///< Transpose this matrix
	void GetTranspose( Matrix2 &m ) const									///< return Transpose of this matrix
	{
		m.data[0] = data[0];
		m.data[2] = data[1];
		m.data[1] = data[2];
		m.data[3] = data[3];
	}
	Matrix2 GetTranspose() const { Matrix2 t; GetTranspose(t); return t; }	///< return Transpose of this matrix

	TYPE GetDeterminant() const { return data[0]*data[3]-data[2]*data[1]; }		///< Get the determinant of this matrix

	void Invert()					///< Invert this matrix
	{
		TYPE det = GetDeterminant();
		TYPE data_0 = data[0];
		data[0] =  data[3] / det;
		data[1] = -data[1] / det;
		data[2] = -data[2] / det;
		data[3] =  data_0  / det;
	}
	void GetInverse( Matrix2 &inverse ) const { inverse=*this; inverse.Invert(); }	///< Get the inverse of this matrix
	Matrix2 GetInverse() const { Matrix2 inv(*this); inv.Invert(); return inv; }	///< Get the inverse of this matrix

	/// Orthogonalizes the matrix and removes the scale component, preserving the x direction
	void OrthogonalizeX()
	{
		Point2<TYPE> *p = (Point2<TYPE>*)data;
		p[0].Normalize();
		p[1] -= p[0] * (p[1]%p[0]);
		p[1].Normalize();
	}
	/// Orthogonalizes the matrix and removes the scale component, preserving the y direction
	void OrthogonalizeY()
	{
		Point2<TYPE> *p = (Point2<TYPE>*)data;
		p[1].Normalize();
		p[0] -= p[1] * (p[0]%p[1]);
		p[0].Normalize();
	}

	/// Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const { return abs(data[0] - TYPE(1)) < tollerance && abs(data[1]) < tollerance && abs(data[2]) < tollerance && abs(data[3] - TYPE(1)) < tollerance; }

	/// Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( TYPE tollerance=TYPE(0.001) ) const { return abs(data[0] - data[2]) < tollerance; }

	//////////////////////////////////////////////////////////////////////////
	///@name Static Methods

	/// Returns an identity matrix
	static Matrix2 MatrixIdentity() { Matrix2 m; m.SetIdentity(); return m; }
	/// Returns a rotation matrix about the given axis by angle in radians
	static Matrix2 MatrixRotation( TYPE angle ) { Matrix2 m; m.SetRotation(angle); return m; }
	/// Returns a uniform scale matrix
	static Matrix2 MatrixScale( TYPE uniformScale ) { Matrix2 m; m.SetScale(uniformScale); return m; }
	/// Returns a scale matrix
	static Matrix2 MatrixScale( TYPE scaleX, TYPE scaleY ) { Matrix2 m; m.SetScale(scaleX,scaleY); return m; }
	/// Returns a scale matrix
	static Matrix2 MatrixScale( const Point2<TYPE> &scale ) { Matrix2 m; m.SetScale(scale); return m; }


	/////////////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

/// 3x3 matrix class.
/// Its data stores 9-value array of column-major matrix elements.
/// You can use Matrix3 with Point3<TYPE> to transform 3D points.
/// Both post-multiplication and pre-multiplication are supported.

template <typename TYPE>
class Matrix3
{
	
	friend Matrix3 operator * ( const TYPE value, const Matrix3 &right ) { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = value * right.data[i]; return buffer; }	///< multiple matrix by a value
	friend Matrix3 operator / ( const TYPE value, const Matrix3 &right ) { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = value / right.data[i]; return buffer; }	///< divide a value by a matrix
	friend Matrix3 Inverse( Matrix3 &m ) { return m.GetInverse(); }	///< return the inverse of the matrix

public:

	/// Elements of the matrix are column-major:
	/// | 0  3  6 |
	/// | 1  4  7 |
	/// | 2  5  8 |
	TYPE data[9];

	//////////////////////////////////////////////////////////////////////////
	///@name Constructors

	Matrix3() {}																							///< Default constructor
	Matrix3( const Matrix3 &matrix ) { for ( int i=0; i<9; i++ ) data[i]=matrix.data[i]; }					///< Copy constructor
	template <typename T> explicit Matrix3<TYPE>( const Matrix3<T> &matrix ) { for ( int i=0; i<9; i++ ) data[i]=(TYPE)matrix.data[i]; }		///< Copy constructor for different types
	explicit Matrix3( const TYPE *array ) { Set(array); }													///< Initialize the matrix using an array of 9 values
	explicit Matrix3( TYPE v ) { SetScaledIdentity(v); }													///< Initialize the matrix as identity scaled by v
	explicit Matrix3( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z ) { Set(x,y,z); }	///< Initialize the matrix using x,y,z vectors as columns
	explicit Matrix3( const Matrix2<TYPE> &m ) { 
		data[0] = m.data[0]; data[1] = m.data[1]; data[2] = TYPE(0);
		data[3] = m.data[2]; data[4] = m.data[3]; data[5] = TYPE(0);
		data[6] = TYPE(0);   data[7] = TYPE(0);   data[8] = TYPE(1);
	}


	//////////////////////////////////////////////////////////////////////////
	///@name Set & Get Methods

	/// Set all the values as zero
	void Zero() { for ( int i=0; i<9; i++ ) data[i] = TYPE(0); }
	/// Set matrix using an array of 9 values
	void Set( const TYPE *array ) { for ( int i=0; i<9; i++ ) data[i] = array[i]; } 
	/// Matrix formulation of the cross product
	void Set( const Point3<TYPE> &p ) { data[0]=TYPE(0); data[1]=p.z; data[2]=-p.y; data[3]=-p.z; data[4]=TYPE(0); data[5]=p.x; data[6]=p.y; data[7]=-p.x; data[8]=TYPE(0); }
	/// Set matrix using x,y,z vectors as columns
	void Set( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z ) { x.GetValue(&data[0]); y.GetValue(&data[3]); z.GetValue(&data[6]); }
	/// Set a row of the matrix
	void SetRow( int row, TYPE x, TYPE y, TYPE z ) { data[row]=x; data[row+3]=y; data[row+6]=z; }
	/// Set a column of the matrix
	void SetColumn( int column, TYPE x, TYPE y, TYPE z ) { data[3*column]=x; data[3*column+1]=y; data[3*column+2]=z; }
	/// Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	/// Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity(TYPE v) { for(int i=0; i<9; i++ ) data[i]=(i%4==0) ? v : TYPE(0); }
	/// Set view matrix using position, target and approximate up vector
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
	/// Set matrix using normal, and approximate x direction
	void SetNormal(const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y = normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal); }
	/// Set as rotation matrix around x axis
	void SetRotationX( TYPE angle ) { SetRotationX( sin(angle), cos(angle) ); }
	/// Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = TYPE(1); data[3] = TYPE(0);  data[6] = TYPE(0);
		data[1] = TYPE(0); data[4] = cosAngle; data[7] = -sinAngle;
		data[2] = TYPE(0); data[5] = sinAngle; data[8] = cosAngle;
	}
	/// Set as rotation matrix around y axis
	void SetRotationY( TYPE angle ) { SetRotationY( sin(angle), cos(angle) ); }
	/// Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = cosAngle;  data[3] = TYPE(0); data[6] = sinAngle;
		data[1] = TYPE(0);   data[4] = TYPE(1); data[7] = TYPE(0);
		data[2] = -sinAngle; data[5] = TYPE(0); data[8] = cosAngle;
	}
	/// Set as rotation matrix around z axis
	void SetRotationZ( TYPE angle ) { SetRotationZ( sin(angle), cos(angle) ); }
	/// Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = cosAngle; data[3] = -sinAngle; data[6] = TYPE(0);
		data[1] = sinAngle; data[4] = cosAngle;  data[7] = TYPE(0);
		data[2] = TYPE(0);  data[5] = TYPE(0);   data[8] = TYPE(1);
	}
	/// Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
	void SetRotationXYZ( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = sin(angleX);
		const TYPE cx = cos(angleX);
		const TYPE sy = sin(angleY);
		const TYPE cy = cos(angleY);
		const TYPE sz = sin(angleZ);
		const TYPE cz = cos(angleZ);
		data[0] = cy*cz; data[3] = cz*sx*sy - cx*sz; data[6] = cx*cz*sy + sx*sz;
		data[1] = cy*sz; data[4] = cx*cz + sx*sy*sz; data[7] =-cz*sx + cx*sy*sz;
		data[2] =-sy;    data[5] = cy*sx;            data[8] = cx*cy;
	}
	/// Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
	void SetRotationZYX( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = sin(angleX);
		const TYPE cx = cos(angleX);
		const TYPE sy = sin(angleY);
		const TYPE cy = cos(angleY);
		const TYPE sz = sin(angleZ);
		const TYPE cz = cos(angleZ);
		data[0] = cy*cz;            data[3] = -cy*sz;           data[6] = sy;
		data[1] = cx*sz + sx*sy*cz; data[4] = cx*cz - sx*sy*sz; data[7] = -sx*cy;
		data[2] = sx*sz - cx*sy*cz; data[5] = sx*cz + cx*sy*sz; data[8] = cx*cy;
	}
	/// Set a rotation matrix about the given axis by angle
	void SetRotation( const Point3<TYPE> &axis, TYPE angle ) { SetRotation(axis,sin(angle),cos(angle)); }
	/// Set a rotation matrix about the given axis by cos and sin of angle
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

		data[0] = tx * axis.x + cosAngle;
		data[1] = txy + sz;
		data[2] = txz - sy;

		data[3] = txy - sz;
		data[4] = ty * axis.y + cosAngle;
		data[5] = tyz + sx;

		data[6] = txz + sy;
		data[7] = tyz - sx;
		data[8] = tz * axis.z + cosAngle;
	}
	/// Set a rotation matrix that sets [from] unit vector to [to] unit vector
	void SetRotation( const Point3<TYPE> &from, const Point3<TYPE> &to )
	{
		TYPE c = from.Dot(to);
		if ( c > TYPE(0.9999999) ) SetIdentity();
		else {
			TYPE s = sqrt(TYPE(1) - c*c);
			Point3<TYPE> axis = from.Cross(to).GetNormalized();
			SetRotation(axis, s, c);
		}
	}
	/// Sets a uniform scale matrix
	void SetScale( TYPE uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	/// Sets a scale matrix
	void SetScale( TYPE scaleX, TYPE scaleY, TYPE scaleZ ) { data[0]=scaleX; data[1]=TYPE(0); data[2]=TYPE(0); data[3]=TYPE(0); data[4]=scaleY; data[5]=TYPE(0); data[6]=TYPE(0); data[7]=TYPE(0); data[8]=scaleZ; }
	/// Sets a scale matrix
	void SetScale( const Point3<TYPE> &scale ) { SetScale(scale.x,scale.y,scale.z); }
	/// Removes the scale component of the matrix
	void SetNoScale() { Point3<TYPE> *p = (Point3<TYPE>*)data; p[0].Normalize(); p[1].Normalize(); p[2].Normalize(); }
	/// Sets the matrix as the tensor product (outer product) of two vectors
	void SetTensorProduct( const Point3<TYPE> &v0, const Point3<TYPE> &v1 )
	{
		data[0]=v0.x*v1.x;  data[3]=v0.x*v1.y;  data[6]=v0.x*v1.z;
		data[1]=v0.y*v1.x;  data[4]=v0.y*v1.y;  data[7]=v0.y*v1.z;
		data[2]=v0.z*v1.x;  data[5]=v0.z*v1.y;  data[8]=v0.z*v1.z;
	}

	// Get Row and Column
	Point3<TYPE>  GetRow   ( int row )                  const { return Point3<TYPE>( data[row], data[row+3], data[row+6] ); }
	void          GetRow   ( int row, Point3<TYPE> &p ) const { p.Set( data[row], data[row+3], data[row+6] ); }
	void          GetRow   ( int row, TYPE *array )     const { array[0]=data[row]; array[1]=data[row+3]; array[2]=data[row+6]; }
	Point3<TYPE>  GetColumn( int col )                  const { return Point3<TYPE>( &data[col*3] ); }
	void          GetColumn( int col, Point3<TYPE> &p ) const { p.Set( &data[col*3] ); }
	void          GetColumn( int col, TYPE *array )     const { array[0]=data[col*3]; array[1]=data[col*3+1]; array[2]=data[col*3+2]; }
	// These methods return the diagonal component of the matrix
	Point3<TYPE>  GetDiagonal()                         const { return Point3<TYPE>( data[0], data[4], data[8] ); }
	void          GetDiagonal( Point3<TYPE> &p )        const { p.Set( data[0], data[4], data[8] ); }
	void	      GetDiagonal( TYPE *array )            const { array[0]=data[0]; array[1]=data[4]; array[2]=data[8]; }
	// These methods can be used for converting the 2x2 portion of the matrix into a Matrix
	void          GetSubMatrix2( TYPE mdata[4] )        const { mdata[0]=data[0]; mdata[1]=data[1]; mdata[2]=data[3]; mdata[3]=data[4]; }
	void          GetSubMatrix2( Matrix2<TYPE> &m )     const { GetSubMatrix2(m.data); }
	Matrix2<TYPE> GetSubMatrix2()                       const { Matrix2<TYPE> m; GetSubMatrix2(m.data); return m; }


	//////////////////////////////////////////////////////////////////////////
	///@name Overloaded Operators

	const Matrix3& operator = ( const Matrix3 &right ) { for (int i=0; i<9; i++) data[i] = right.data[i]; return *this; }

	// Overloaded comparison operators 
	bool operator == ( const Matrix3 &right ) const { for (int i=0; i<9; i++) if ( data[i] != right.data[i] ) return false; return true; }		///< compare equal
	bool operator != ( const Matrix3 &right ) const { return ! ( *this == right ); } ///< compare not equal

	// Overloaded subscript operators
	TYPE&       operator () ( int row, int column )       { return data[ column * 3 + row ]; }	///< subscript operator
	const TYPE& operator () ( int row, int column ) const { return data[ column * 3 + row ]; }	///< constant subscript operator
	TYPE&       operator [] ( int i )       { return data[i]; }	///< subscript operator
	const TYPE& operator [] ( int i ) const { return data[i]; }	///< constant subscript operator
	
	// Unary operators
	Matrix3 operator - () const { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = -data[i]; return buffer; }	///< negative matrix

	// Binary operators
	Matrix3      operator * ( const TYPE value )      const { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = data[i] * value; return buffer; }	///< multiple matrix by a value
	Matrix3      operator / ( const TYPE value )      const { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = data[i] / value; return buffer; }	///< divide matrix by a value;
	Point3<TYPE> operator * ( const Point3<TYPE> &p ) const { return Point3<TYPE>( p.x*data[0] + p.y*data[3] + p.z*data[6], p.x*data[1] + p.y*data[4] + p.z*data[7], p.x*data[2] + p.y*data[5] + p.z*data[8] ); }
	Matrix3      operator + ( const Matrix3 &right  ) const { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = data[i] + right.data[i]; return buffer; }	///< add two Matrices
	Matrix3      operator - ( const Matrix3 &right  ) const { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = data[i] - right.data[i]; return buffer; }	///< subtract one Matrix3 from an other
	Matrix3      operator * ( const Matrix3 &right  ) const	///< multiply a matrix with an other
	{
		Matrix3 buffer;
		for (int i = 0; i < 3; i++) {
			for (int k = 0; k < 3; k++) {
				TYPE v = TYPE(0);
				for (int j = 0; j < 3; j++) {
					v += data[i + 3*j] * right.data[j + 3*k];
				}
				buffer.data[i + 3*k] = v;
			}
		}
		return buffer;
	}

	// Assignment operators
	void operator *= ( const TYPE value )     { for (int i=0; i<9; i++) data[i] *= value; }			///< multiply a matrix with a value modify this matrix
	void operator /= ( const TYPE value )     { for (int i=0; i<9; i++) data[i] /= value; }			///< divide the matrix by a value modify the this matrix
	void operator += ( const Matrix3 &right ) { for (int i=0; i<9; i++) data[i] += right.data[i]; }	///< add two Matrices modify this
	void operator -= ( const Matrix3 &right ) { for (int i=0; i<9; i++) data[i] -= right.data[i]; }	///< subtract one Matrix3 from an other modify this matrix
	void operator *= ( const Matrix3 &right ) { *this = *this * right; }	///< multiply a matrix with an other modify this matrix

	//////////////////////////////////////////////////////////////////////////
	///@name Other Public Methods

	void Transpose()															///< Transpose this matrix
	{
		for (int i = 1; i < 3; i++) {
			for (int j = 0; j < i; j++) {
				TYPE temp = data[i * 3 + j];
				data[i * 3 + j] = data[j * 3 + i];
				data[j * 3 + i] = temp;
			}
		}
	}
	void GetTranspose( Matrix3 &m ) const										///< return Transpose of this matrix
	{
		m.data[0] = data[0];
		m.data[3] = data[1];
		m.data[6] = data[2];
		m.data[1] = data[3];
		m.data[4] = data[4];
		m.data[7] = data[5];
		m.data[2] = data[6];
		m.data[5] = data[7];
		m.data[8] = data[8];
	}
	Matrix3 GetTranspose() const { Matrix3 t; GetTranspose(t); return t; }	///< return Transpose of this matrix

	TYPE GetDeterminant() const { return data[0]*(data[4]*data[8]-data[7]*data[5]) + data[3]*(data[7]*data[2]-data[1]*data[8]) + data[6]*(data[1]*data[5]-data[4]*data[2]); }	///< Get the determinant of this matrix

	void Invert() { Matrix3 inv; GetInverse(inv); *this=inv; }					///< Invert this matrix
	void GetInverse( Matrix3 &inverse ) const									///< Get the inverse of this matrix
	{
		TYPE det = GetDeterminant();
		inverse.data[0] = (data[4]*data[8] - data[7]*data[5]) / det;
		inverse.data[1] = (data[7]*data[2] - data[1]*data[8]) / det;
		inverse.data[2] = (data[1]*data[5] - data[4]*data[2]) / det;
		inverse.data[3] = (data[6]*data[5] - data[3]*data[8]) / det;
		inverse.data[4] = (data[0]*data[8] - data[6]*data[2]) / det;
		inverse.data[5] = (data[3]*data[2] - data[0]*data[5]) / det;
		inverse.data[6] = (data[3]*data[7] - data[6]*data[4]) / det;
		inverse.data[7] = (data[6]*data[1] - data[0]*data[7]) / det;
		inverse.data[8] = (data[0]*data[4] - data[3]*data[1]) / det;
	}
	Matrix3 GetInverse() const { Matrix3 inv; GetInverse(inv); return inv; }	///< Get the inverse of this matrix

	/// Orthogonalizes the matrix and removes the scale component, preserving the x direction
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
	/// Orthogonalizes the matrix and removes the scale component, preserving the y direction
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
	/// Orthogonalizes the matrix and removes the scale component, preserving the z direction
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


	/// Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const
	{
		return abs(data[0]-TYPE(1)) < tollerance && abs(data[1])         < tollerance && abs(data[2])         < tollerance && 
			   abs(data[3])         < tollerance && abs(data[4]-TYPE(1)) < tollerance && abs(data[5])         < tollerance &&
			   abs(data[6])         < tollerance && abs(data[7])         < tollerance && abs(data[8]-TYPE(1)) < tollerance;
	}

	/// Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( TYPE tollerance=TYPE(0.001) ) const { return abs(data[1] - data[3]) < tollerance && abs(data[2] - data[6]) < tollerance && abs(data[5] - data[7]) < tollerance; }

	
	//////////////////////////////////////////////////////////////////////////
	///@name Static Methods

	/// Returns an identity matrix
	static Matrix3 MatrixIdentity() { Matrix3 m; m.SetIdentity(); return m; }
	/// Returns a view matrix using position, target and approximate up vector
	static Matrix3 MatrixView( const Point3<TYPE> &target, Point3<TYPE> &up ) { Matrix3 m; m.SetView(target,up); return m; }
	/// Returns a matrix using normal, and approximate x direction
	static Matrix3 MatrixNormal( const Point3<TYPE> &normal, Point3<TYPE> &dir ) { Matrix3 m; m.SetNormal(normal,dir); return m; }
	/// Returns a rotation matrix around x axis by angle in radians
	static Matrix3 MatrixRotationX( TYPE angle ) { Matrix3 m; m.SetRotationX(angle); return m; }
	/// Returns a rotation matrix around y axis by angle in radians
	static Matrix3 MatrixRotationY( TYPE angle ) { Matrix3 m; m.SetRotationY(angle); return m; }
	/// Returns a rotation matrix around z axis by angle in radians
	static Matrix3 MatrixRotationZ( TYPE angle ) { Matrix3 m; m.SetRotationZ(angle); return m; }
	/// Returns a rotation matrix around x, y, and then z axes by angle in radians (Rz * Ry * Rx)
	static Matrix3 MatrixRotationXYZ( TYPE angleX, TYPE angleY, TYPE angleZ ) { Matrix3 m; m.SetRotationXYZ(angleX,angleY,angleZ); return m; }
	/// Returns a rotation matrix around z, y, and then x axes by angle in radians (Rx * Ry * Rz)
	static Matrix3 MatrixRotationZYX( TYPE angleX, TYPE angleY, TYPE angleZ ) { Matrix3 m; m.SetRotationZYX(angleX,angleY,angleZ); return m; }
	/// Returns a rotation matrix about the given axis by angle in radians
	static Matrix3 MatrixRotation( const Point3<TYPE> &axis, TYPE angle ) { Matrix3 m; m.SetRotation(axis,angle); return m; }
	/// Returns a rotation matrix that sets [from] unit vector to [to] unit vector
	static Matrix3 MatrixRotation( const Point3<TYPE> &from, const Point3<TYPE> &to ) { Matrix3 m; m.SetRotation(from,to); return m; }
	/// Returns a uniform scale matrix
	static Matrix3 MatrixScale( TYPE uniformScale ) { Matrix3 m; m.SetScale(uniformScale); return m; }
	/// Returns a scale matrix
	static Matrix3 MatrixScale( TYPE scaleX, TYPE scaleY, TYPE scaleZ ) { Matrix3 m; m.SetScale(scaleX,scaleY,scaleZ); return m; }
	/// Returns a scale matrix
	static Matrix3 MatrixScale( const Point3<TYPE> &scale ) { Matrix3 m; m.SetScale(scale); return m; }


	/////////////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

/// 3x4 matrix class.
/// Its data stores 12-value array of column-major matrix elements.
/// I chose column-major format to be compatible with OpenGL
/// You can use Matrix34 with Point3<TYPE> and Point4<TYPE>
/// to transform 3D and 4D points.
/// Both post-multiplication and pre-multiplication are supported.

template <typename TYPE>
class Matrix34
{
	
	friend Matrix34 operator * ( const TYPE value, const Matrix34 &right ) { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = value * right.data[i]; return buffer; }	///< multiple matrix by a value
	friend Matrix34 operator / ( const TYPE value, const Matrix34 &right ) { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = value / right.data[i]; return buffer; }	///< multiple matrix by a value
	friend Matrix34 Inverse( Matrix34 &m ) { return m.GetInverse(); }	///< return the inverse of the matrix

public:

	/// Elements of the matrix are column-major:
	/// | 0   3   6   9 |
	/// | 1   4   7  10 |
	/// | 2   5   8  11 |
	TYPE data[12];


	//////////////////////////////////////////////////////////////////////////
	///@name Constructors

	Matrix34() {}																				///< Default constructor
	Matrix34( const Matrix34 &matrix ) { for ( int i=0; i<12; i++ ) data[i]=matrix.data[i]; }	///< Copy constructor
	template <typename T> explicit Matrix34<TYPE>( const Matrix34<T> &matrix ) { for ( int i=0; i<12; i++ ) data[i]=(TYPE)matrix.data[i]; }		///< Copy constructor for different types
	explicit Matrix34( const TYPE *array ) { Set(array); }										///< Initialize the matrix using an array of 9 values
	explicit Matrix34( TYPE v ) { SetScaledIdentity(v); }										///< Initialize the matrix as identity scaled by v
	explicit Matrix34( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z, const Point3<TYPE> &pos ) { Set(x,y,z,pos); }	///< Initialize the matrix using x,y,z vectors and coordinate center
	explicit Matrix34( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Set(pos,normal,dir); }				///< Initialize the matrix using position, normal, and approximate x direction
	explicit Matrix34( const Matrix3<TYPE> &m ) { for (int i=0; i<9; i++) data[i]=m.data[i]; data[9]=TYPE(0); data[10]=TYPE(0); data[11]=TYPE(0); }
	explicit Matrix34( const Matrix3<TYPE> &m, const Point3<TYPE> &pos ) { for (int i=0; i<9; i++) data[i]=m.data[i]; data[9]=pos.x; data[10]=pos.y; data[11]=pos.z; }
	explicit Matrix34( const Matrix2<TYPE> &m ) { 
		data[ 0] = m.data[ 0]; data[ 1] = m.data[ 1]; data[ 2] = TYPE(0);
		data[ 3] = m.data[ 2]; data[ 4] = m.data[ 3]; data[ 5] = TYPE(0);
		data[ 6] = TYPE(0);    data[ 7] = TYPE(0);    data[ 8] = TYPE(1);
		data[ 9] = TYPE(0);    data[10] = TYPE(0);    data[11] = TYPE(0);
	}


	//////////////////////////////////////////////////////////////////////////
	///@name Set & Get Methods

	/// Set all the values as zero
	void Zero() { for ( int i=0; i<12; i++ ) data[i] = TYPE(0); }
	/// Set Matrix using an array of 12 values
	void Set( const TYPE *array ) { for ( int i=0; i<12; i++ ) data[i] = array[i]; } 
	/// Set matrix using x,y,z vectors and coordinate center
	void Set( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z, const Point3<TYPE> &pos ) { x.GetValue(&data[0]); y.GetValue(&data[3]); z.GetValue(&data[6]); pos.GetValue(&data[9]); }
	/// Set matrix using position, normal, and approximate x direction
	void Set( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,pos); }
	/// Set a row of the matrix
	void SetRow( int row, TYPE x, TYPE y, TYPE z, TYPE w ) { data[row]=x; data[row+3]=y; data[row+6]=z; data[row+9]=w; }
	/// Set a column of the matrix
	void SetColumn( int column, TYPE x, TYPE y, TYPE z ) { data[3*column]=x; data[3*column+1]=y; data[3*column+2]=z; }
	/// Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	/// Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity(TYPE v) { for(int i=0; i<12; i++) data[i]=(i%4==0) ? v : TYPE(0); }
	/// Set view matrix using position, target and approximate up vector
	void SetView( const Point3<TYPE> &pos, const Point3<TYPE> &target, const Point3<TYPE> &up )
	{
		Point3<TYPE> f = target - pos;
		f.Normalize();
		Point3<TYPE> s = f.Cross(up);
		s.Normalize();
		Point3<TYPE> u = s.Cross(f);
		Matrix34 m;
		m.data[0] = s.x; m.data[1] = u.x; m.data[2] = -f.x;
		m.data[3] = s.y; m.data[4] = u.y; m.data[5] = -f.y;
		m.data[6] = s.z; m.data[7] = u.z; m.data[8] = -f.z;
		Matrix34 t;
		t.data[ 9] = -pos.x;
		t.data[10] = -pos.y;
		t.data[11] = -pos.z;
		*this = m * t;
	}
	/// Set matrix using normal and approximate x direction
	void SetNormal(const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,Point3<TYPE>(TYPE(0),TYPE(0),TYPE(0))); }
	/// Set as rotation matrix around x axis
	void SetRotationX( TYPE angle ) { SetRotationX( sin(angle), cos(angle) ); }
	/// Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = TYPE(1); data[3] = TYPE(0);  data[6] = TYPE(0);   data[ 9] = TYPE(0);
		data[1] = TYPE(0); data[4] = cosAngle; data[7] = -sinAngle; data[10] = TYPE(0);
		data[2] = TYPE(0); data[5] = sinAngle; data[8] = cosAngle;  data[11] = TYPE(0);
	}
	/// Set as rotation matrix around y axis
	void SetRotationY( TYPE angle ) { SetRotationY( sin(angle), cos(angle) ); }
	/// Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = cosAngle;  data[3] = TYPE(0); data[6] = sinAngle; data[ 9] = TYPE(0);
		data[1] = TYPE(0);   data[4] = TYPE(1); data[7] = TYPE(0);  data[10] = TYPE(0);
		data[2] = -sinAngle; data[5] = TYPE(0); data[8] = cosAngle; data[11] = TYPE(0);
	}
	/// Set as rotation matrix around z axis
	void SetRotationZ( TYPE angle ) { SetRotationZ( sin(angle), cos(angle) ); }
	/// Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = cosAngle; data[3] = -sinAngle; data[6] = TYPE(0); data[ 9] = TYPE(0);
		data[1] = sinAngle; data[4] = cosAngle;  data[7] = TYPE(0); data[10] = TYPE(0);
		data[2] = TYPE(0);  data[5] = TYPE(0);   data[8] = TYPE(1); data[11] = TYPE(0);
	}
	/// Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
	void SetRotationXYZ( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = sin(angleX);
		const TYPE cx = cos(angleX);
		const TYPE sy = sin(angleY);
		const TYPE cy = cos(angleY);
		const TYPE sz = sin(angleZ);
		const TYPE cz = cos(angleZ);
		data[0] = cy*cz; data[3] = cz*sx*sy - cx*sz; data[6] = cx*cz*sy + sx*sz; data[ 9] = TYPE(0);
		data[1] = cy*sz; data[4] = cx*cz + sx*sy*sz; data[7] =-cz*sx + cx*sy*sz; data[10] = TYPE(0);
		data[2] =-sy;    data[5] = cy*sx;            data[8] = cx*cy;            data[11] = TYPE(0);
	}
	/// Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
	void SetRotationZYX( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = sin(angleX);
		const TYPE cx = cos(angleX);
		const TYPE sy = sin(angleY);
		const TYPE cy = cos(angleY);
		const TYPE sz = sin(angleZ);
		const TYPE cz = cos(angleZ);
		data[0] = cy*cz;            data[3] =-cy*sz;            data[6] = sy;    data[ 9] = TYPE(0);
		data[1] = cx*sz + sx*sy*cz; data[4] = cx*cz - sx*sy*sz; data[7] =-sx*cy; data[10] = TYPE(0);
		data[2] = sx*sz - cx*sy*cz; data[5] = sx*cz + cx*sy*sz; data[8] = cx*cy; data[11] = TYPE(0);
	}
	/// Set a rotation matrix about the given axis by angle
	void SetRotation( const Point3<TYPE> &axis, TYPE angle ) { SetRotation(axis,sin(angle),cos(angle)); }
	/// Set a rotation matrix about the given axis by cos and sin of angle
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
		data[ 0] = tx * axis.x + cosAngle;
		data[ 1] = txy + sz;
		data[ 2] = txz - sy;
		data[ 3] = txy - sz;
		data[ 4] = ty * axis.y + cosAngle;
		data[ 5] = tyz + sx;
		data[ 6] = txz + sy;
		data[ 7] = tyz - sx;
		data[ 8] = tz * axis.z + cosAngle;
		data[ 9] = TYPE(0);
		data[10] = TYPE(0);
		data[11] = TYPE(0);
	}
	/// Set a rotation matrix that sets [from] unit vector to [to] unit vector
	void SetRotation( const Point3<TYPE> &from, const Point3<TYPE> &to )
	{
		TYPE c = from.Dot(to);
		if ( c > TYPE(0.9999999) ) SetIdentity();
		else {
			TYPE s = sqrt(TYPE(1) - c*c);
			Point3<TYPE> axis = from.Cross(to).GetNormalized();
			SetRotation(axis, s, c);
		}
	}
	/// Sets a uniform scale matrix
	void SetScale( TYPE uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	/// Sets a scale matrix
	void SetScale( TYPE scaleX, TYPE scaleY, TYPE scaleZ ) { data[0]=scaleX; data[1]=TYPE(0); data[2]=TYPE(0); data[3]=TYPE(0); data[4]=scaleY; data[5]=TYPE(0); data[6]=TYPE(0); data[7]=TYPE(0); data[8]=scaleZ; data[9]=TYPE(0); data[10]=TYPE(0); data[11]=TYPE(0); }
	/// Sets a scale matrix
	void SetScale( const Point3<TYPE> &scale ) { SetScale(scale.x,scale.y,scale.z); }
	/// Removes the scale component of the matrix
	void SetNoScale() { Point3<TYPE> *p = (Point3<TYPE>*)data; p[0].Normalize(); p[1].Normalize(); p[2].Normalize(); }
	/// Sets a translation matrix with no rotation or scale
	void SetTrans( const Point3<TYPE> &move ) { for(int i=0; i<9; i++) data[i]=(i%4==0)?TYPE(1):TYPE(0); data[9]=move.x; data[10]=move.y; data[11]=move.z; }
	/// Adds a translation to the matrix
	void AddTrans( const Point3<TYPE> &move ) { data[9]+=move.x; data[10]+=move.y; data[11]+=move.z; }
	/// Sets the translation component of the matrix
	void SetTransComponent( const Point3<TYPE> &move ) { data[9]=move.x; data[10]=move.y; data[11]=move.z; }

	// Get Row and Column
	Point4<TYPE>  GetRow   ( int row )                  const { return Point4<TYPE>( data[row], data[row+3], data[row+6], data[row+9] ); }
	void          GetRow   ( int row, Point4<TYPE> &p ) const { p.Set( data[row], data[row+3], data[row+6], data[row+9] ); }
	void          GetRow   ( int row, TYPE *array )     const { array[0]=data[row]; array[1]=data[row+3]; array[2]=data[row+6]; array[3]=data[row+9]; }
	Point3<TYPE>  GetColumn( int col )                  const { return Point3<TYPE>( &data[col*3] ); }
	void          GetColumn( int col, Point3<TYPE> &p ) const { p.Set( &data[col*3] ); }
	void          GetColumn( int col, TYPE *array )     const { array[0]=data[col*3]; array[1]=data[col*3+1]; array[2]=data[col*3+2]; }
	// These methods return the diagonal component of the matrix
	Point4<TYPE>  GetDiagonal()                         const { return Point4<TYPE>( data[0], data[4], data[8], 1 ); }
	void          GetDiagonal( Point4<TYPE> &p )        const { p.Set( data[0], data[4], data[8], 1 ); }
	void	      GetDiagonal( TYPE *array )            const { array[0]=data[0]; array[1]=data[4]; array[2]=data[8]; array[3]=1; }
	// These methods can be used for converting the 3x3 portion of the matrix into a Matrix3
	void          GetSubMatrix3( TYPE mdata[9] )        const { for ( int i=0; i<9; i++ ) mdata[i] = data[i]; }
	void          GetSubMatrix3( Matrix3<TYPE> &m )     const { GetSubMatrix3(m.data); }
	Matrix3<TYPE> GetSubMatrix3()                       const { Matrix3<TYPE> m; GetSubMatrix3(m.data); return m; }
	// These methods can be used for converting the 2x2 portion of the matrix into a Matrix2
	void          GetSubMatrix2( TYPE mdata[4] )        const { mdata[0]=data[0]; mdata[1]=data[1]; mdata[2]=data[3]; mdata[3]=data[4]; }
	void          GetSubMatrix2( Matrix2<TYPE> &m )     const { GetSubMatrix2(m.data); }
	Matrix2<TYPE> GetSubMatrix2()                       const { Matrix2<TYPE> m; GetSubMatrix2(m.data); return m; }
	// These methods return the translation component of the matrix
	void          GetTrans( TYPE trans[3] )             const { trans[0]=data[9]; trans[1]=data[10]; trans[2]=data[11]; }
	void          GetTrans( Point3<TYPE> &p )           const { p.x=data[9]; p.y=data[10]; p.z=data[11]; }
	Point3<TYPE>  GetTrans()                            const { Point3<TYPE> p; GetTrans(p); return p; }

	//////////////////////////////////////////////////////////////////////////
	///@name Overloaded Operators

	const Matrix34& operator = ( const Matrix34 &right ) { for (int i=0; i<12; i++) data[i] = right.data[i]; return *this; }	///< assign matrix

	// Overloaded comparison operators 
	bool operator == ( const Matrix34 &right ) const { for (int i=0; i<12; i++) if ( data[i] != right.data[i] ) return false; return true; }	///< compare equal
	bool operator != ( const Matrix34 &right ) const { return ! ( *this == right ); } ///< compare not equal

	// Overloaded subscript operators
	TYPE&       operator () ( int row, int column )       { return data[ column * 3 + row ]; }	///< subscript operator
	const TYPE& operator () ( int row, int column ) const { return data[ column * 3 + row ]; }	///< constant subscript operator
	TYPE&       operator [] ( int i )       { return data[i]; }	///< subscript operator
	const TYPE& operator [] ( int i ) const { return data[i]; }	///< constant subscript operator

	// Unary operators
	Matrix34 operator - () const { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i]=-data[i]; return buffer; }	///< negative matrix

	// Binary operators
	Matrix34     operator * ( const TYPE value )      const { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = data[i] * value; return buffer; }	///< multiple matrix by a value
	Matrix34     operator / ( const TYPE value )      const { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = data[i] / value; return buffer; }	///< divide matrix by a value;
	Point3<TYPE> operator * ( const Point3<TYPE> &p ) const { return Point3<TYPE>( p.x*data[0] + p.y*data[3] + p.z*data[6] +     data[9], p.x*data[1] + p.y*data[4] + p.z*data[7] +     data[10], p.x*data[2] + p.y*data[5] + p.z*data[8] +     data[11] ); }
	Point4<TYPE> operator * ( const Point4<TYPE> &p ) const { return Point4<TYPE>( p.x*data[0] + p.y*data[3] + p.z*data[6] + p.w*data[9], p.x*data[1] + p.y*data[4] + p.z*data[7] + p.w*data[10], p.x*data[2] + p.y*data[5] + p.z*data[8] + p.w*data[11], p.w ); }
	Matrix34     operator + ( const Matrix34 &right ) const { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = data[i] + right.data[i]; return buffer; }	///< add two Matrices
	Matrix34     operator - ( const Matrix34 &right ) const { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = data[i] - right.data[i]; return buffer; }	///< subtract one Matrix4 from an other
	Matrix34     operator * ( const Matrix34 &right ) const	///< multiply a matrix with an other
	{
		Matrix34 buffer;
		for (int myRow = 0; myRow < 3; myRow++) {
			for (int rtCol = 0; rtCol < 4; rtCol++) {
				TYPE v = TYPE(0);
				for (int j = 0; j < 3; j++) {
					v += data[myRow + 3*j] * right.data[j + 3*rtCol];
				}
				buffer.data[myRow + 3*rtCol] = v;
			}
			buffer.data[myRow + 9] += data[myRow + 9];
		}
		return buffer;
	}

	// Assignment operators
	void operator *= ( const TYPE value )      { for (int i=0; i<12; i++) data[i] *= value; }			///< multiply a matrix with a value modify this matrix
	void operator /= ( const TYPE value )      { for (int i=0; i<12; i++) data[i] /= value; }			///< divide the matrix by a value modify the this matrix
	void operator += ( const Matrix34 &right ) { for (int i=0; i<12; i++) data[i] += right.data[i]; }	///< add two Matrices modify this
	void operator -= ( const Matrix34 &right ) { for (int i=0; i<12; i++) data[i] -= right.data[i]; }	///< subtract one Matrix4 from an other modify this matrix
	void operator *= ( const Matrix34 &right ) { *this = *this * right; }	///< multiply a matrix with an other modify this matrix

	//////////////////////////////////////////////////////////////////////////
	///@name Other Public Methods

	/// Transpose this matrix
	void Transpose()
	{
		for (int i = 1; i < 3; i++) {
			for (int j = 0; j < i; j++) {
				TYPE temp = data[i * 3 + j];
				data[i * 3 + j] = data[j * 3 + i];
				data[j * 3 + i] = temp;
			}
		}
		data[ 9] = TYPE(0);
		data[10] = TYPE(0);
		data[11] = TYPE(0);
	}

	/// return Transpose of this matrix
	void GetTranspose( Matrix34 &m ) const
	{
		m.data[ 0] = data[0];
		m.data[ 1] = data[3];
		m.data[ 2] = data[6];
		m.data[ 3] = data[1];
		m.data[ 4] = data[4];
		m.data[ 5] = data[7];
		m.data[ 6] = data[2];
		m.data[ 7] = data[5];
		m.data[ 8] = data[8];
		m.data[ 9] = TYPE(0);
		m.data[10] = TYPE(0);
		m.data[11] = TYPE(0);
	}
	Matrix34 GetTranspose() const { Matrix34 t; GetTranspose(t); return t; }	///< return Transpose of this matrix

	TYPE GetDeterminant() const	///< Get the determinant of this matrix
	{
		return	data[0]*( data[4]*data[8] - data[7]*data[5] ) + 
				data[3]*( data[7]*data[2] - data[1]*data[8] ) + 
				data[6]*( data[1]*data[5] - data[4]*data[2] );
	}
	void Invert() { Matrix34 inv; GetInverse(inv); *this=inv; }	///< Invert this matrix
	void GetInverse( Matrix34 &inverse ) const						///< Get the inverse of this matrix
	{
		TYPE det = GetDeterminant();
		inverse.data[ 0] = ( data[4]*data[8] - data[7]*data[5] ) / det;
		inverse.data[ 1] = ( data[7]*data[2] - data[1]*data[8] ) / det;
		inverse.data[ 2] = ( data[1]*data[5] - data[4]*data[2] ) / det;
		inverse.data[ 3] = ( data[6]*data[5] - data[3]*data[8] ) / det;
		inverse.data[ 4] = ( data[0]*data[8] - data[6]*data[2] ) / det;
		inverse.data[ 5] = ( data[3]*data[2] - data[0]*data[5] ) / det;
		inverse.data[ 6] = ( data[3]*data[7] - data[6]*data[4] ) / det;
		inverse.data[ 7] = ( data[6]*data[1] - data[0]*data[7] ) / det;
		inverse.data[ 8] = ( data[0]*data[4] - data[3]*data[1] ) / det;
		inverse.data[ 9] = - ( data[9]*inverse.data[0] + data[10]*inverse.data[3] + data[11]*inverse.data[6] );
		inverse.data[10] = - ( data[9]*inverse.data[1] + data[10]*inverse.data[4] + data[11]*inverse.data[7] );
		inverse.data[11] = - ( data[9]*inverse.data[2] + data[10]*inverse.data[5] + data[11]*inverse.data[8] );
	}
	Matrix34 GetInverse() const { Matrix34 inv; GetInverse(inv); return inv; }	///< Get the inverse of this matrix

	/// Orthogonalizes the matrix and removes the scale component, preserving the x direction
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
	/// Orthogonalizes the matrix and removes the scale component, preserving the y direction
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
	/// Orthogonalizes the matrix and removes the scale component, preserving the z direction
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

	bool IsCloseToIdentity( TYPE tollerance=TYPE(0.001) ) const		///< Returns if the matrix is close to identity. Closeness is determined by the tollerance parameter.
	{
		for (int i = 0; i < 12; i++) {
			TYPE v = (i % 5 == 0) ? TYPE(1) : TYPE(0);
			if (fabsf(data[i] - v) > tollerance) return false;
		}
		return true;
	}

	/// Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const
	{
		return abs(data[0]-TYPE(1)) < tollerance && abs(data[ 1])         < tollerance && abs(data[ 2])         < tollerance && 
			   abs(data[3])         < tollerance && abs(data[ 4]-TYPE(1)) < tollerance && abs(data[ 5])         < tollerance &&
			   abs(data[6])         < tollerance && abs(data[ 7])         < tollerance && abs(data[ 8]-TYPE(1)) < tollerance &&
			   abs(data[9])         < tollerance && abs(data[10])         < tollerance && abs(data[11])         < tollerance;
	}

	/// Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( TYPE tollerance=TYPE(0.001) ) const
	{
		return abs(data[ 1] - data[3]) < tollerance && 
			   abs(data[ 2] - data[6]) < tollerance &&
			   abs(data[ 5] - data[7]) < tollerance &&
			   abs(data[ 9])           < tollerance &&
			   abs(data[10])           < tollerance &&
			   abs(data[11])           < tollerance;
	}

	//////////////////////////////////////////////////////////////////////////
	///@name Static Methods

	/// Returns an identity matrix
	static Matrix34 MatrixIdentity() { Matrix34 m; m.SetIdentity(); return m; }
	/// Returns a matrix using normal, and approximate x direction
	static Matrix34 MatrixNormal( const Point3<TYPE> &normal, Point3<TYPE> &dir ) { Matrix34 m; m.SetNormal(normal,dir); return m; }
	/// Returns a rotation matrix around x axis by angle in radians
	static Matrix34 MatrixRotationX( TYPE angle ) { Matrix34 m; m.SetRotationX(angle); return m; }
	/// Returns a rotation matrix around y axis by angle in radians
	static Matrix34 MatrixRotationY( TYPE angle ) { Matrix34 m; m.SetRotationY(angle); return m; }
	/// Returns a rotation matrix around z axis by angle in radians
	static Matrix34 MatrixRotationZ( TYPE angle ) { Matrix34 m; m.SetRotationZ(angle); return m; }
	/// Returns a rotation matrix around x, y, and then z axes by angle in radians (Rz * Ry * Rx)
	static Matrix34 MatrixRotationXYZ( TYPE angleX, TYPE angleY, TYPE angleZ ) { Matrix34 m; m.SetRotationXYZ(angleX,angleY,angleZ); return m; }
	/// Returns a rotation matrix around z, y, and then x axes by angle in radians (Rx * Ry * Rz)
	static Matrix34 MatrixRotationZYX( TYPE angleX, TYPE angleY, TYPE angleZ ) { Matrix34 m; m.SetRotationZYX(angleX,angleY,angleZ); return m; }
	/// Returns a rotation matrix about the given axis by angle in radians
	static Matrix34 MatrixRotation( const Point3<TYPE> &axis, TYPE angle ) { Matrix34 m; m.SetRotation(axis,angle); return m; }
	/// Returns a rotation matrix that sets [from] unit vector to [to] unit vector
	static Matrix34 MatrixRotation( const Point3<TYPE> &from, const Point3<TYPE> &to ) { Matrix34 m; m.SetRotation(from,to); return m; }
	/// Returns a uniform scale matrix
	static Matrix34 MatrixScale( TYPE uniformScale ) { Matrix34 m; m.SetScale(uniformScale); return m; }
	/// Returns a scale matrix
	static Matrix34 MatrixScale( TYPE scaleX, TYPE scaleY, TYPE scaleZ ) { Matrix34 m; m.SetScale(scaleX,scaleY,scaleZ); return m; }
	/// Returns a scale matrix
	static Matrix34 MatrixScale( const Point3<TYPE> &scale ) { Matrix34 m; m.SetScale(scale); return m; }
	/// Returns a translation matrix with no rotation or scale
	static Matrix34 MatrixTrans( const Point3<TYPE> &move ) { Matrix34 m; m.SetTrans(move); return m; }


	/////////////////////////////////////////////////////////////////////////////////
};


//-------------------------------------------------------------------------------

/// 4x4 matrix class.
/// Its data stores 16-value array of column-major matrix elements.
/// I chose column-major format to be compatible with OpenGL
/// You can use Matrix4 with Point3<TYPE> and Point4<TYPE>
/// to transform 3D and 4D points.
/// Both post-multiplication and pre-multiplication are supported.

template <typename TYPE>
class Matrix4
{
	
	friend Matrix4 operator + ( const TYPE value, const Matrix4 &right ) { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i] = value + right.data[i]; return buffer; }	///< add a value to a matrix
	friend Matrix4 operator - ( const TYPE value, const Matrix4 &right ) { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i] = value - right.data[i]; return buffer; }	///< subtract the matrix from a value
	friend Matrix4 operator * ( const TYPE value, const Matrix4 &right ) { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i] = value * right.data[i]; return buffer; }	///< multiple matrix by a value
	friend Matrix4 operator / ( const TYPE value, const Matrix4 &right ) { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i] = value / right.data[i]; return buffer; }	///< multiple matrix by a value
	friend Matrix4 Inverse( Matrix4 &m ) { return m.GetInverse(); }	///< return the inverse of the matrix

public:

	/// Elements of the matrix are column-major:
	/// | 0   4   8  12 |
	/// | 1   5   9  13 |
	/// | 2   6  10  14 |
	/// | 3   7  11  15 |
	TYPE data[16];


	//////////////////////////////////////////////////////////////////////////
	///@name Constructors

	Matrix4() {}																				///< Default constructor
	Matrix4( const Matrix4 &matrix ) { for ( int i=0; i<16; i++ ) data[i]=matrix.data[i]; }		///< Copy constructor
	template <typename T> explicit Matrix4<TYPE>( const Matrix4<T> &matrix ) { for ( int i=0; i<16; i++ ) data[i]=(TYPE)matrix.data[i]; }		///< Copy constructor for different types
	explicit Matrix4( const TYPE *array ) { Set(array); }										///< Initialize the matrix using an array of 9 values
	explicit Matrix4( TYPE v ) { SetScaledIdentity(v); }										///< Initialize the matrix as identity scaled by v
	explicit Matrix4( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z, const Point3<TYPE> &pos ) { Set(x,y,z,pos); }	///< Initialize the matrix using x,y,z vectors and coordinate center
	explicit Matrix4( const Point4<TYPE> &x, const Point4<TYPE> &y, const Point4<TYPE> &z, const Point4<TYPE> &w   ) { Set(x,y,z,w);   }	///< Initialize the matrix using x,y,z vectors as columns
	explicit Matrix4( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Set(pos,normal,dir); }				///< Initialize the matrix using position, normal, and approximate x direction
	explicit Matrix4( const Matrix34<TYPE> &m ) { 
		data[ 0]=m.data[ 0]; data[ 1]=m.data[ 1]; data[ 2]=m.data[ 2]; data[ 3]=TYPE(0); 
		data[ 4]=m.data[ 3]; data[ 5]=m.data[ 4]; data[ 6]=m.data[ 5]; data[ 7]=TYPE(0); 
		data[ 8]=m.data[ 6]; data[ 9]=m.data[ 7]; data[10]=m.data[ 8]; data[11]=TYPE(0); 
		data[12]=m.data[ 9]; data[13]=m.data[10]; data[14]=m.data[11]; data[15]=TYPE(1);
	}
	explicit Matrix4( const Matrix3<TYPE> &m ) { 
		data[ 0]=m.data[ 0]; data[ 1]=m.data[ 1]; data[ 2]=m.data[ 2]; data[ 3]=TYPE(0); 
		data[ 4]=m.data[ 3]; data[ 5]=m.data[ 4]; data[ 6]=m.data[ 5]; data[ 7]=TYPE(0); 
		data[ 8]=m.data[ 6]; data[ 9]=m.data[ 7]; data[10]=m.data[ 8]; data[11]=TYPE(0); 
		data[12]=TYPE(0);    data[13]=TYPE(0);    data[14]=TYPE(0);    data[15]=TYPE(1);
	}
	explicit Matrix4( const Matrix2<TYPE> &m ) { 
		data[ 0]=m.data[ 0]; data[ 1]=m.data[ 1]; data[ 2]=TYPE(0); data[ 3]=TYPE(0); 
		data[ 4]=m.data[ 2]; data[ 5]=m.data[ 3]; data[ 6]=TYPE(0); data[ 7]=TYPE(0); 
		data[ 8]=TYPE(0);    data[ 9]=TYPE(0);    data[10]=TYPE(1); data[11]=TYPE(0); 
		data[12]=TYPE(0);    data[13]=TYPE(0);    data[14]=TYPE(0); data[15]=TYPE(1);
	}


	//////////////////////////////////////////////////////////////////////////
	///@name Set & Get Methods

	/// Set all the values as zero
	void Zero() { for ( int i=0; i<16; i++ ) data[i] = TYPE(0); }
	/// Set Matrix using an array of 16 values
	void Set( const TYPE *array ) { for ( int i=0; i<16; i++ ) data[i] = array[i]; } 
	/// Set matrix using x,y,z vectors and coordinate center
	void Set( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z, const Point3<TYPE> &pos ) { x.GetValue(&data[0]); data[3]=TYPE(0); y.GetValue(&data[4]); data[7]=TYPE(0); z.GetValue(&data[8]); data[11]=TYPE(0); pos.GetValue(&data[12]); data[15]=TYPE(1); }
	/// Set matrix using x,y,z,w vectors
	void Set( const Point4<TYPE> &x, const Point4<TYPE> &y, const Point4<TYPE> &z, const Point4<TYPE> &w ) { x.GetValue(&data[0]); y.GetValue(&data[4]); z.GetValue(&data[8]); w.GetValue(&data[12]); }
	/// Set matrix using position, normal, and approximate x direction
	void Set( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,pos); }
	/// Set a row of the matrix
	void SetRow( int row, TYPE x, TYPE y, TYPE z, TYPE w ) { data[row]=x; data[row+4]=y; data[row+8]=z; data[row+12]=w; }
	/// Set a column of the matrix
	void SetColumn( int column, TYPE x, TYPE y, TYPE z, TYPE w ) { data[4*column]=x; data[4*column+1]=y; data[4*column+2]=z; data[4*column+3]=w; }
	/// Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	/// Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity(TYPE v) { for(int i=0; i<16; i++) data[i]=(i%5==0) ? v : TYPE(0); }
	/// Set view matrix using position, target and approximate up vector
	void SetView( const Point3<TYPE> &pos, const Point3<TYPE> &target, const Point3<TYPE> &up )
	{
		Point3<TYPE> f = target - pos;
		f.Normalize();
		Point3<TYPE> s = f.Cross(up);
		s.Normalize();
		Point3<TYPE> u = s.Cross(f);
		Matrix4 m;
		m.data[ 0]=s.x; m.data[ 1]=u.x; m.data[ 2]=-f.x; m.data[ 3]=TYPE(0);
		m.data[ 4]=s.y; m.data[ 5]=u.y; m.data[ 6]=-f.y; m.data[ 7]=TYPE(0);
		m.data[ 8]=s.z; m.data[ 9]=u.z; m.data[10]=-f.z; m.data[11]=TYPE(0);
		Matrix4 t;
		t.data[12]=-pos.x;
		t.data[13]=-pos.y;
		t.data[14]=-pos.z;
		t.data[15]=TYPE(1);
		*this = m * t;
	}
	/// Set matrix using normal and approximate x direction
	void SetNormal(const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,Point3<TYPE>(TYPE(0),TYPE(0),TYPE(0))); }
	/// Set as rotation matrix around x axis
	void SetRotationX( TYPE angle ) { SetRotationX( sin(angle), cos(angle) ); }
	/// Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = TYPE(1); data[4] = TYPE(0);  data[ 8] = TYPE(0);   data[12] = TYPE(0);
		data[1] = TYPE(0); data[5] = cosAngle; data[ 9] = -sinAngle; data[13] = TYPE(0);
		data[2] = TYPE(0); data[6] = sinAngle; data[10] = cosAngle;  data[14] = TYPE(0);
		data[3] = TYPE(0); data[7] = TYPE(0);  data[11] = TYPE(0);   data[15] = TYPE(1);
	}
	/// Set as rotation matrix around y axis
	void SetRotationY( TYPE angle ) { SetRotationY( sin(angle), cos(angle) ); }
	/// Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = cosAngle;  data[4] = TYPE(0); data[ 8] = sinAngle; data[12] = TYPE(0);
		data[1] = TYPE(0);   data[5] = TYPE(1); data[ 9] = TYPE(0);  data[13] = TYPE(0);
		data[2] = -sinAngle; data[6] = TYPE(0); data[10] = cosAngle; data[14] = TYPE(0);
		data[3] = TYPE(0);   data[7] = TYPE(0); data[11] = TYPE(0);  data[15] = TYPE(1);
	}
	/// Set as rotation matrix around z axis
	void SetRotationZ( TYPE angle ) { SetRotationZ( sin(angle), cos(angle) ); }
	/// Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = cosAngle; data[4] = -sinAngle; data[ 8] = TYPE(0); data[12] = TYPE(0);
		data[1] = sinAngle; data[5] = cosAngle;  data[ 9] = TYPE(0); data[13] = TYPE(0);
		data[2] = TYPE(0);  data[6] = TYPE(0);   data[10] = TYPE(1); data[14] = TYPE(0);
		data[3] = TYPE(0);  data[7] = TYPE(0);   data[11] = TYPE(0); data[15] = TYPE(1);
	}
	/// Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
	void SetRotationXYZ( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = sin(angleX);
		const TYPE cx = cos(angleX);
		const TYPE sy = sin(angleY);
		const TYPE cy = cos(angleY);
		const TYPE sz = sin(angleZ);
		const TYPE cz = cos(angleZ);
		data[0] = cy*cz;   data[4] = cz*sx*sy - cx*sz; data[ 8] = cx*cz*sy + sx*sz; data[12] = TYPE(0);
		data[1] = cy*sz;   data[5] = cx*cz + sx*sy*sz; data[ 9] =-cz*sx + cx*sy*sz; data[13] = TYPE(0);
		data[2] =-sy;      data[6] = cy*sx;            data[10] = cx*cy;            data[14] = TYPE(0);
		data[3] = TYPE(0); data[7] = TYPE(0);          data[11] = TYPE(0);          data[15] = TYPE(1);
	}
	/// Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
	void SetRotationZYX( TYPE angleX, TYPE angleY, TYPE angleZ )
	{
		const TYPE sx = sin(angleX);
		const TYPE cx = cos(angleX);
		const TYPE sy = sin(angleY);
		const TYPE cy = cos(angleY);
		const TYPE sz = sin(angleZ);
		const TYPE cz = cos(angleZ);
		data[0] = cy*cz;            data[4] =-cy*sz;            data[ 8] = sy;      data[12] = TYPE(0);
		data[1] = cx*sz + sx*sy*cz; data[5] = cx*cz - sx*sy*sz; data[ 9] =-sx*cy;   data[13] = TYPE(0);
		data[2] = sx*sz - cx*sy*cz; data[6] = sx*cz + cx*sy*sz; data[10] = cx*cy;   data[14] = TYPE(0);
		data[3] = TYPE(0);          data[7] = TYPE(0);          data[11] = TYPE(0); data[15] = TYPE(1);
	}
	/// Set a rotation matrix about the given axis by angle
	void SetRotation( const Point3<TYPE> &axis, TYPE angle ) { SetRotation(axis,sin(angle),cos(angle)); }
	/// Set a rotation matrix about the given axis by cos and sin of angle
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
		data[ 0] = tx * axis.x + cosAngle;
		data[ 1] = txy + sz;
		data[ 2] = txz - sy;
		data[ 3] = TYPE(0);
		data[ 4] = txy - sz;
		data[ 5] = ty * axis.y + cosAngle;
		data[ 6] = tyz + sx;
		data[ 7] = TYPE(0);
		data[ 8] = txz + sy;
		data[ 9] = tyz - sx;
		data[10] = tz * axis.z + cosAngle;
		data[11] = TYPE(0);
		data[12] = TYPE(0);
		data[13] = TYPE(0);
		data[14] = TYPE(0);
		data[15] = TYPE(1);
	}
	/// Set a rotation matrix that sets [from] unit vector to [to] unit vector
	void SetRotation( const Point3<TYPE> &from, const Point3<TYPE> &to )
	{
		TYPE c = from.Dot(to);
		if ( c > TYPE(0.9999999) ) SetIdentity();
		else {
			TYPE s = sqrt(TYPE(1) - c*c);
			Point3<TYPE> axis = from.Cross(to).GetNormalized();
			SetRotation(axis, s, c);
		}
	}
	/// Sets a uniform scale matrix
	void SetScale( TYPE uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	/// Sets a scale matrix
	void SetScale( TYPE scaleX, TYPE scaleY, TYPE scaleZ ) { data[0]=scaleX; data[1]=TYPE(0); data[2]=TYPE(0); data[3]=TYPE(0); data[4]=TYPE(0); data[5]=scaleY; data[6]=TYPE(0); data[7]=TYPE(0); data[8]=TYPE(0); data[9]=TYPE(0); data[10]=scaleZ; data[11]=TYPE(0); data[12]=TYPE(0); data[13]=TYPE(0); data[14]=TYPE(0); data[15]=TYPE(1); }
	/// Sets a scale matrix
	void SetScale( const Point3<TYPE> &scale ) { SetScale(scale.x,scale.y,scale.z); }
	/// Removes the scale component of the matrix
	void SetNoScale() { ((Point3<TYPE>*)&data[0])->Normalize(); ((Point3<TYPE>*)&data[4])->Normalize(); ((Point3<TYPE>*)&data[8])->Normalize(); }
	/// Sets a translation matrix with no rotation or scale
	void SetTrans( const Point3<TYPE> &move ) { for(int i=0; i<12; i++) data[i]=(i%5==0)?TYPE(1):TYPE(0); data[12]=move.x; data[13]=move.y; data[14]=move.z; data[15]=TYPE(1); }
	/// Adds a translation to the matrix
	void AddTrans( const Point3<TYPE> &move ) { data[12]+=move.x; data[13]+=move.y; data[14]+=move.z; }
	/// Sets the translation component of the matrix
	void SetTransComponent( const Point3<TYPE> &move ) { data[12]=move.x; data[13]=move.y; data[14]=move.z; }
	/// Set a project matrix with field of view in radians
	void SetPerspective( TYPE fov, TYPE aspect, TYPE znear, TYPE zfar ) { SetPerspectiveTan(tan(fov*TYPE(0.5)),aspect,znear,zfar); }
	/// Set a project matrix with the tangent of the half field of view (tan_fov_2)
	void SetPerspectiveTan( TYPE tan_fov_2, TYPE aspect, TYPE znear, TYPE zfar )
	{
		const TYPE yScale = TYPE(1) / tan_fov_2;
		const TYPE xScale = yScale / aspect;
		const TYPE zdif = znear - zfar;
		data[ 0] = xScale;
		data[ 1] = TYPE(0);
		data[ 2] = TYPE(0);
		data[ 3] = TYPE(0);
		data[ 4] = TYPE(0);
		data[ 5] = yScale;
		data[ 6] = TYPE(0);
		data[ 7] = TYPE(0);
		data[ 8] = TYPE(0);
		data[ 9] = TYPE(0);
		data[10] = (zfar + znear) / zdif;
		data[11] = -TYPE(1);
		data[12] = TYPE(0);
		data[13] = TYPE(0);
		data[14] = (2 * zfar*znear) / zdif;
		data[15] = TYPE(0);
	}
	/// Sets the matrix as the tensor product (outer product) of two vectors
	void SetTensorProduct( const Point4<TYPE> &v0, const Point4<TYPE> &v1 )
	{
		data[0]=v0.x*v1.x;  data[4]=v0.x*v1.y;  data[ 8]=v0.x*v1.z;  data[12]=v0.x*v1.w;
		data[1]=v0.y*v1.x;  data[5]=v0.y*v1.y;  data[ 9]=v0.y*v1.z;  data[13]=v0.y*v1.w;
		data[2]=v0.z*v1.x;  data[6]=v0.z*v1.y;  data[10]=v0.z*v1.z;  data[14]=v0.z*v1.w;
		data[3]=v0.w*v1.x;  data[7]=v0.w*v1.y;  data[11]=v0.w*v1.z;  data[15]=v0.w*v1.w;
	}

	// Get Row and Column
	Point4<TYPE>   GetRow   ( int row )                  const { return Point4<TYPE>( data[row], data[row+4], data[row+8], data[row+12] ); }
	void           GetRow   ( int row, Point4<TYPE> &p ) const { p.Set( data[row], data[row+4], data[row+8], data[row+12] ); }
	void           GetRow   ( int row, TYPE *array )     const { array[0]=data[row]; array[1]=data[row+4]; array[2]=data[row+8]; array[3]=data[row+12]; }
	Point4<TYPE>   GetColumn( int col )                  const { return Point4<TYPE>( &data[col*4] ); }
	void           GetColumn( int col, Point4<TYPE> &p ) const { p.Set( &data[col*4] ); }
	void           GetColumn( int col, TYPE *array )     const { array[0]=data[col*4]; array[1]=data[col*4+1]; array[2]=data[col*4+2]; array[3]=data[col*4+3]; }
	// These methods return the diagonal component of the matrix
	Point4<TYPE>   GetDiagonal()                         const { return Point4<TYPE>( data[0], data[5], data[10], data[15] ); }
	void           GetDiagonal( Point4<TYPE> &p )        const { p.Set( data[0], data[5], data[10], data[15] ); }
	void	       GetDiagonal( TYPE *array )            const { array[0]=data[0]; array[1]=data[5]; array[2]=data[10]; array[3]=data[15]; }
	// These methods can be used for converting the 3x4 portion of the matrix into a Matrix34
	void           GetSubMatrix34( TYPE mdata[12] )      const { mdata[0]=data[0]; mdata[1]=data[1]; mdata[2]=data[2]; mdata[3]=data[4]; mdata[4]=data[5]; mdata[5]=data[6]; mdata[6]=data[8]; mdata[7]=data[9]; mdata[8]=data[10]; mdata[9]=data[12]; mdata[10]=data[13]; mdata[11]=data[14]; }
	void           GetSubMatrix34( Matrix34<TYPE> &m )   const { GetSubMatrix34(m.data); }
	Matrix34<TYPE> GetSubMatrix34()                      const { Matrix34<TYPE> m; GetSubMatrix34(m.data); return m; }
	// These methods can be used for converting the 3x3 portion of the matrix into a Matrix3
	void           GetSubMatrix3 ( TYPE mdata[9] )       const { mdata[0]=data[0]; mdata[1]=data[1]; mdata[2]=data[2]; mdata[3]=data[4]; mdata[4]=data[5]; mdata[5]=data[6]; mdata[6]=data[8]; mdata[7]=data[9]; mdata[8]=data[10]; }
	void           GetSubMatrix3 ( Matrix3<TYPE> &m )    const { GetSubMatrix3(m.data); }
	Matrix3<TYPE>  GetSubMatrix3 ()                      const { Matrix3<TYPE> m; GetSubMatrix3(m.data); return m; }
	// These methods can be used for converting the 2x2 portion of the matrix into a Matrix2
	void           GetSubMatrix2 ( TYPE mdata[4] )       const { mdata[0]=data[0]; mdata[1]=data[1]; mdata[2]=data[4]; mdata[3]=data[5]; }
	void           GetSubMatrix2 ( Matrix2<TYPE> &m )    const { GetSubMatrix2(m.data); }
	Matrix2<TYPE>  GetSubMatrix2 ()                      const { Matrix2<TYPE> m; GetSubMatrix2(m.data); return m; }
	// These methods return the translation component of the matrix
	void           GetTrans( TYPE trans[3] )             const { trans[0]=data[12]; trans[1]=data[13]; trans[2]=data[14]; }
	void           GetTrans( Point3<TYPE> &p )           const { p.x=data[12]; p.y=data[13]; p.z=data[14]; }
	Point3<TYPE>   GetTrans()                            const { Point3<TYPE> p; GetTrans(p); return p; }


	//////////////////////////////////////////////////////////////////////////
	///@name Overloaded Operators

	const Matrix4& operator = ( const Matrix4 &right ) { for (int i=0; i<16; i++) data[i] = right.data[i]; return *this; }	///< assign matrix

	// Overloaded comparison operators 
	bool operator == ( const Matrix4 &right ) const { for (int i=0; i<16; i++) if ( data[i] != right.data[i] ) return false; return true; }	///< compare equal
	bool operator != ( const Matrix4 &right ) const { return ! ( *this == right ); } ///< compare not equal

	// Overloaded subscript operators
	TYPE&       operator () ( int row, int column )       { return data[ column * 4 + row ]; }	///< subscript operator
	const TYPE& operator () ( int row, int column ) const { return data[ column * 4 + row ]; }	///< constant subscript operator
	TYPE&       operator [] ( int i )       { return data[i]; }	///< subscript operator
	const TYPE& operator [] ( int i ) const { return data[i]; }	///< constant subscript operator

	// Unary operators
	Matrix4 operator - () const { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i]=-data[i]; return buffer; }	///< negative matrix

	// Binary operators
	Matrix4      operator * ( const TYPE value )      const { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i] = data[i] * value; return buffer; }	///< multiple matrix by a value
	Matrix4      operator / ( const TYPE value )      const { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i] = data[i] / value; return buffer; }	///< divide matrix by a value;
	Point4<TYPE> operator * ( const Point3<TYPE>& p ) const { return Point4<TYPE>( p.x*data[0] + p.y*data[4] + p.z*data[8] +     data[12], p.x*data[1] + p.y*data[5] + p.z*data[9] +     data[13], p.x*data[2] + p.y*data[6] + p.z*data[10] +     data[14], p.x*data[3] + p.y*data[7] + p.z*data[11] +     data[15] ); }
	Point4<TYPE> operator * ( const Point4<TYPE>& p ) const { return Point4<TYPE>( p.x*data[0] + p.y*data[4] + p.z*data[8] + p.w*data[12], p.x*data[1] + p.y*data[5] + p.z*data[9] + p.w*data[13], p.x*data[2] + p.y*data[6] + p.z*data[10] + p.w*data[14], p.x*data[3] + p.y*data[7] + p.z*data[11] + p.w*data[15] ); }
	Matrix4      operator + ( const Matrix4 &right  ) const { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i] = data[i] + right.data[i]; return buffer; }	///< add two Matrices
	Matrix4      operator - ( const Matrix4 &right  ) const { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i] = data[i] - right.data[i]; return buffer; }	///< subtract one Matrix4 from an other
	Matrix4      operator * ( const Matrix4 &right  ) const	///< multiply a matrix with an other
	{
		Matrix4 buffer;
		for (int i = 0; i < 4; i++) {
			for (int k = 0; k < 4; k++) {
				TYPE v = TYPE(0);
				for (int j = 0; j < 4; j++) {
					v += data[i + 4*j] * right.data[j + 4*k];
				}
				buffer.data[i + 4*k] = v;
			}
		}
		return buffer;
	}
	Matrix4 operator * ( const Matrix34<TYPE> &right ) const	///< multiply a matrix with an other
	{
		Matrix4 buffer;
		for (int i = 0; i < 4; i++) {
			for (int k = 0; k < 4; k++) {
				TYPE v = TYPE(0);
				for (int j = 0; j < 3; j++) {
					v += data[i + 4*j] * right.data[j + 3*k];
				}
				buffer.data[i + 4*k] = v;
			}
			buffer.data[i + 12] += data[i + 12];
		}
		return buffer;
	}

	// Assignment operators
	void operator *= ( const TYPE value )     { for (int i=0; i<16; i++) data[i] *= value; }	///< multiply a matrix with a value modify this matrix
	void operator /= ( const TYPE value )     { for (int i=0; i<16; i++) data[i] /= value; }	///< divide the matrix by a value modify the this matrix
	void operator += ( const Matrix4 &right ) { for (int i=0; i<16; i++) data[i] += right.data[i]; }	///< add two Matrices modify this
	void operator -= ( const Matrix4 &right ) { for (int i=0; i<16; i++) data[i] -= right.data[i]; }	///< subtract one Matrix4 from an other modify this matrix
	void operator *= ( const Matrix4 &right ) { *this = *this * right; }			///< multiply a matrix with an other modify this matrix
	void operator *= ( const Matrix34<TYPE> &right ) { *this = *this * right; }	///< multiply a matrix with an other modify this matrix

	//////////////////////////////////////////////////////////////////////////
	///@name Other Public Methods

	void Transpose()															///< Transpose this matrix
	{
		for (int i = 1; i < 4; i++) {
			for (int j = 0; j < i; j++) {
				TYPE temp = data[i * 4 + j];
				data[i * 4 + j] = data[j * 4 + i];
				data[j * 4 + i] = temp;
			}
		}
	}
	void GetTranspose( Matrix4 &m ) const										///< return Transpose of this matrix
	{
		m.data[ 0] = data[0];
		m.data[ 1] = data[4];
		m.data[ 2] = data[8];
		m.data[ 3] = data[12];
		m.data[ 4] = data[1];
		m.data[ 5] = data[5];
		m.data[ 6] = data[9];
		m.data[ 7] = data[13];
		m.data[ 8] = data[2];
		m.data[ 9] = data[6];
		m.data[10] = data[10];
		m.data[11] = data[14];
		m.data[12] = data[3];
		m.data[13] = data[7];
		m.data[14] = data[11];
		m.data[15] = data[15];
	}
	Matrix4 GetTranspose() const { Matrix4 t; GetTranspose(t); return t; }	///< return Transpose of this matrix

	TYPE GetDeterminant() const	///< Get the determinant of this matrix
	{
		const TYPE data0914_1013 = data[ 9]*data[14] - data[10]*data[13];
		const TYPE data0613_0514 = data[ 6]*data[13] - data[ 5]*data[14];
		const TYPE data0510_0609 = data[ 5]*data[10] - data[ 6]*data[ 9];
		const TYPE data0114_0213 = data[ 1]*data[14] - data[ 2]*data[13];
		const TYPE data0209_0110 = data[ 2]*data[ 9] - data[ 1]*data[10];
		const TYPE data0106_0205 = data[ 1]*data[ 6] - data[ 2]*data[ 5];
		return	data[ 0] * ( data[ 7]*( data0914_1013) + data[11]*( data0613_0514) + data[15]*( data0510_0609) ) + 
				data[ 4] * ( data[ 3]*(-data0914_1013) + data[11]*( data0114_0213) + data[15]*( data0209_0110) ) + 
				data[ 8] * ( data[ 3]*(-data0613_0514) + data[ 7]*(-data0114_0213) + data[15]*( data0106_0205) ) + 
				data[12] * ( data[ 3]*(-data0510_0609) + data[ 7]*(-data0209_0110) + data[11]*(-data0106_0205) );
	}
	void Invert() { Matrix4 inv; GetInverse(inv); *this=inv; }					///< Invert this matrix
	void GetInverse( Matrix4 &inverse ) const									///< Get the inverse of this matrix
	{
		const TYPE data0914_1013 = data[ 9]*data[14] - data[10]*data[13];
		const TYPE data0613_0514 = data[ 6]*data[13] - data[ 5]*data[14];
		const TYPE data0510_0609 = data[ 5]*data[10] - data[ 6]*data[ 9];
		const TYPE data0114_0213 = data[ 1]*data[14] - data[ 2]*data[13];
		const TYPE data0209_0110 = data[ 2]*data[ 9] - data[ 1]*data[10];
		const TYPE data0106_0205 = data[ 1]*data[ 6] - data[ 2]*data[ 5];
		const TYPE data0815_1112 = data[ 8]*data[15] - data[11]*data[12];
		const TYPE data0712_0415 = data[ 7]*data[12] - data[ 4]*data[15];
		const TYPE data0411_0708 = data[ 4]*data[11] - data[ 7]*data[ 8];
		const TYPE data0015_0312 = data[ 0]*data[15] - data[ 3]*data[12];
		const TYPE data0308_0011 = data[ 3]*data[ 8] - data[ 0]*data[11];
		const TYPE data0007_0304 = data[ 0]*data[ 7] - data[ 3]*data[ 4];
		const TYPE det = data[ 0] * ( data[ 7]*( data0914_1013) + data[11]*( data0613_0514) + data[15]*( data0510_0609) ) + 
		                 data[ 4] * ( data[ 3]*(-data0914_1013) + data[11]*( data0114_0213) + data[15]*( data0209_0110) ) + 
		                 data[ 8] * ( data[ 3]*(-data0613_0514) + data[ 7]*(-data0114_0213) + data[15]*( data0106_0205) ) + 
		                 data[12] * ( data[ 3]*(-data0510_0609) + data[ 7]*(-data0209_0110) + data[11]*(-data0106_0205) );
		inverse.data[ 0] = ( data[ 7]*( data0914_1013) + data[11]*( data0613_0514) + data[15]*( data0510_0609) ) / det;
		inverse.data[ 1] = ( data[ 3]*(-data0914_1013) + data[11]*( data0114_0213) + data[15]*( data0209_0110) ) / det;
		inverse.data[ 2] = ( data[ 3]*(-data0613_0514) + data[ 7]*(-data0114_0213) + data[15]*( data0106_0205) ) / det;
		inverse.data[ 3] = ( data[ 3]*(-data0510_0609) + data[ 7]*(-data0209_0110) + data[11]*(-data0106_0205) ) / det;
		inverse.data[ 4] = ( data[ 6]*( data0815_1112) + data[10]*( data0712_0415) + data[14]*( data0411_0708) ) / det;
		inverse.data[ 5] = ( data[ 2]*(-data0815_1112) + data[10]*( data0015_0312) + data[14]*( data0308_0011) ) / det;
		inverse.data[ 6] = ( data[ 2]*(-data0712_0415) + data[ 6]*(-data0015_0312) + data[14]*( data0007_0304) ) / det;
		inverse.data[ 7] = ( data[ 2]*(-data0411_0708) + data[ 6]*(-data0308_0011) + data[10]*(-data0007_0304) ) / det;
		inverse.data[ 8] = ( data[ 5]*(-data0815_1112) + data[ 9]*(-data0712_0415) + data[13]*(-data0411_0708) ) / det;
		inverse.data[ 9] = ( data[ 1]*( data0815_1112) + data[ 9]*(-data0015_0312) + data[13]*(-data0308_0011) ) / det;
		inverse.data[10] = ( data[ 1]*( data0712_0415) + data[ 5]*( data0015_0312) + data[13]*(-data0007_0304) ) / det;
		inverse.data[11] = ( data[ 1]*( data0411_0708) + data[ 5]*( data0308_0011) + data[ 9]*( data0007_0304) ) / det;
		inverse.data[12] = ( data[ 4]*(-data0914_1013) + data[ 8]*(-data0613_0514) + data[12]*(-data0510_0609) ) / det;
		inverse.data[13] = ( data[ 0]*( data0914_1013) + data[ 8]*(-data0114_0213) + data[12]*(-data0209_0110) ) / det;
		inverse.data[14] = ( data[ 0]*( data0613_0514) + data[ 4]*( data0114_0213) + data[12]*(-data0106_0205) ) / det;
		inverse.data[15] = ( data[ 0]*( data0510_0609) + data[ 4]*( data0209_0110) + data[ 8]*( data0106_0205) ) / det;
	}
	Matrix4 GetInverse() const { Matrix4 inv; GetInverse(inv); return inv; }	///< Get the inverse of this matrix

	/// Orthogonalizes the matrix and removes the scale component, preserving the x direction
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
	/// Orthogonalizes the matrix and removes the scale component, preserving the y direction
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
	/// Orthogonalizes the matrix and removes the scale component, preserving the z direction
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

	bool IsCloseToIdentity( TYPE tollerance=TYPE(0.001) ) const		///< Returns if the matrix is close to identity. Closeness is determined by the tollerance parameter.
	{
		for (int i = 0; i < 16; i++) {
			TYPE v = (i % 5 == 0) ? TYPE(1) : TYPE(0);
			if (fabsf(data[i] - v) > tollerance) return false;
		}
		return true;
	}

	/// Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const
	{
		return abs(data[ 0]-TYPE(1)) < tollerance && abs(data[ 1])         < tollerance && abs(data[ 2])         < tollerance && abs(data[ 3])         < tollerance && 
			   abs(data[ 4])         < tollerance && abs(data[ 5]-TYPE(1)) < tollerance && abs(data[ 6])         < tollerance && abs(data[ 7])         < tollerance &&
			   abs(data[ 8])         < tollerance && abs(data[ 9])         < tollerance && abs(data[10]-TYPE(1)) < tollerance && abs(data[11])         < tollerance &&
			   abs(data[12])         < tollerance && abs(data[13])         < tollerance && abs(data[14])         < tollerance && abs(data[15]-TYPE(1)) < tollerance;
	}

	/// Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( TYPE tollerance=TYPE(0.001) ) const
	{
		return abs(data[ 1] - data[ 4]) < tollerance && 
			   abs(data[ 2] - data[ 8]) < tollerance &&
			   abs(data[ 3] - data[12]) < tollerance &&
			   abs(data[ 6] - data[ 9]) < tollerance &&
			   abs(data[ 7] - data[13]) < tollerance &&
			   abs(data[11] - data[14]) < tollerance;
	}

	//////////////////////////////////////////////////////////////////////////
	///@name Static Methods

	/// Returns an identity matrix
	static Matrix4 MatrixIdentity() { Matrix4 m; m.SetIdentity(); return m; }
	/// Returns a view matrix using position, target and approximate up vector
	static Matrix4 MatrixView( const Point3<TYPE> &pos, const Point3<TYPE> &target, Point3<TYPE> &up ) { Matrix4 m; m.SetView(pos,target,up); return m; }
	/// Returns a matrix using normal, and approximate x direction
	static Matrix4 MatrixNormal( const Point3<TYPE> &normal, Point3<TYPE> &dir ) { Matrix4 m; m.SetNormal(normal,dir); return m; }
	/// Returns a rotation matrix around x axis by angle in radians
	static Matrix4 MatrixRotationX( TYPE angle ) { Matrix4 m; m.SetRotationX(angle); return m; }
	/// Returns a rotation matrix around y axis by angle in radians
	static Matrix4 MatrixRotationY( TYPE angle ) { Matrix4 m; m.SetRotationY(angle); return m; }
	/// Returns a rotation matrix around z axis by angle in radians
	static Matrix4 MatrixRotationZ( TYPE angle ) { Matrix4 m; m.SetRotationZ(angle); return m; }
	/// Returns a rotation matrix about the given axis by angle in radians
	static Matrix4 MatrixRotation( const Point3<TYPE> &axis, TYPE angle ) { Matrix4 m; m.SetRotation(axis,angle); return m; }
	/// Returns a rotation matrix that sets [from] unit vector to [to] unit vector
	static Matrix4 MatrixRotation( const Point3<TYPE> &from, const Point3<TYPE> &to ) { Matrix4 m; m.SetRotation(from,to); return m; }
	/// Returns a uniform scale matrix
	static Matrix4 MatrixScale( TYPE uniformScale ) { Matrix4 m; m.SetScale(uniformScale); return m; }
	/// Returns a scale matrix
	static Matrix4 MatrixScale( TYPE scaleX, TYPE scaleY, TYPE scaleZ ) { Matrix4 m; m.SetScale(scaleX,scaleY,scaleZ); return m; }
	/// Returns a scale matrix
	static Matrix4 MatrixScale( const Point3<TYPE> &scale ) { Matrix4 m; m.SetScale(scale); return m; }
	/// Returns a translation matrix with no rotation or scale
	static Matrix4 MatrixTrans( const Point3<TYPE> &move ) { Matrix4 m; m.SetTrans(move); return m; }
	/// Returns a project matrix with field of view in radians
	static Matrix4 MatrixPerspective( TYPE fov, TYPE aspect, TYPE znear, TYPE zfar ) { Matrix4 m; m.SetPerspective(fov,aspect,znear,zfar); return m; }
	/// Returns a project matrix with the tangent of the half field of view (tan_fov_2)
	static Matrix4 MatrixPerspectiveTan( TYPE tan_fov_2, TYPE aspect, TYPE znear, TYPE zfar ) { Matrix4 m; m.SetPerspectiveTan(tan_fov_2,aspect,znear,zfar); return m; }


	/////////////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

template<typename TYPE> inline Matrix2<TYPE> operator & ( const Point2<TYPE> &v0, const Point2<TYPE> &v1 ) { Matrix2<TYPE> buffer; buffer.SetTensorProduct(v0,v1); return buffer; }		///< tensor product (outer product) of two vectors
template<typename TYPE> inline Matrix3<TYPE> operator & ( const Point3<TYPE> &v0, const Point3<TYPE> &v1 ) { Matrix3<TYPE> buffer; buffer.SetTensorProduct(v0,v1); return buffer; }		///< tensor product (outer product) of two vectors
template<typename TYPE> inline Matrix4<TYPE> operator & ( const Point4<TYPE> &v0, const Point4<TYPE> &v1 ) { Matrix4<TYPE> buffer; buffer.SetTensorProduct(v0,v1); return buffer; }		///< tensor product (outer product) of two vectors

//-------------------------------------------------------------------------------

typedef Matrix2 <float>  Matrix2f;
typedef Matrix3 <float>  Matrix3f;
typedef Matrix34<float>  Matrix34f;
typedef Matrix4 <float>  Matrix4f;

typedef Matrix2 <double> Matrix2d;
typedef Matrix3 <double> Matrix3d;
typedef Matrix34<double> Matrix34d;
typedef Matrix4 <double> Matrix4d;

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Matrix2f  cyMatrix2f;
typedef cy::Matrix3f  cyMatrix3f;
typedef cy::Matrix34f cyMatrix34f;
typedef cy::Matrix4f  cyMatrix4f;

typedef cy::Matrix2d  cyMatrix2d;
typedef cy::Matrix3d  cyMatrix3d;
typedef cy::Matrix34d cyMatrix34d;
typedef cy::Matrix4d  cyMatrix4d;

//-------------------------------------------------------------------------------

#endif
