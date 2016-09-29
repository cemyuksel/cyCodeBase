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

#ifndef _CY_MATRIX_H_INCLUDED_
#define _CY_MATRIX_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyPoint.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

// Forward declarations
template <typename TYPE> class Matrix3;
template <typename TYPE> class Matrix34;
template <typename TYPE> class Matrix4;

//-------------------------------------------------------------------------------

/// 2x2 matrix class.
/// Its data stores 4-value array of column-major matrix elements.
/// You can use Matrix2 with Point2<TYPE> to transform 2D points.
/// Both post-multiplication and pre-multiplication are supported.

template <typename TYPE>
class Matrix2// : protected Math<TYPE>
{
	
	friend Matrix2 operator + ( const TYPE value, const Matrix2 &right ) { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = value + right.data[i]; return buffer; }	///< add a value to a matrix
	friend Matrix2 operator - ( const TYPE value, const Matrix2 &right ) { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = value - right.data[i]; return buffer; }	///< subtract the matrix from a value
	friend Matrix2 operator * ( const TYPE value, const Matrix2 &right ) { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = value * right.data[i]; return buffer; }	///< multiple matrix by a value
	friend Matrix2 Inverse( Matrix2 &m ) { return m.GetInverse(); }	///< return the inverse of the matrix

public:

	/// Elements of the matrix are column-major:
	/// | 0  2 |
	/// | 1  3 |
	TYPE data[4];

	//////////////////////////////////////////////////////////////////////////
	///@name Constructors

	Matrix2() {}																			///< Default constructor
	Matrix2( const Matrix2 &matrix ) { CY_MEMCOPY(TYPE,data,matrix.data,4); }	///< Copy constructor
	template <typename T> explicit Matrix2<TYPE>( const Matrix2<T> &matrix ) { CY_MEMCONVERT(TYPE,data,matrix.data,4); }		///< Copy constructor for different types
	explicit Matrix2( const TYPE *values ) { Set(values); }									///< Initialize the matrix using an values of 4 values
	explicit Matrix2( const TYPE &v ) { SetScaledIdentity(v); }								///< Initialize the matrix as identity scaled by v
	explicit Matrix2( const Point2<TYPE> &x, const Point2<TYPE> &y ) { Set(x,y); }			///< Initialize the matrix using two vectors as columns
	explicit Matrix2( const Matrix3<TYPE>  &m );
	explicit Matrix2( const Matrix34<TYPE> &m );
	explicit Matrix2( const Matrix4<TYPE>  &m );


	//////////////////////////////////////////////////////////////////////////
	///@name Set & Get Methods

	/// Set all the values as zero
	void Zero() { CY_MEMCLEAR(TYPE,data,4); }
	/// Returns true if the matrix is exactly zero
	bool IsZero() const { for ( int i=0; i<4; i++ ) if ( data[i] != 0 ) return false; return true; }
	/// Copies the matrix data to the given values array of size 4
	void Get( TYPE *values ) { CY_MEMCOPY(TYPE,values,data,4); } 
	/// Set Matrix using an values of 4 values
	void Set( const TYPE *values ) { CY_MEMCOPY(TYPE,data,values,4); } 
	/// Set Matrix using two vectors as columns
	void Set( const Point2<TYPE> &x, const Point2<TYPE> &y ) { x.Get(data); y.Get(data+2); }
	/// Set a row of the matrix
	void SetRow( int row, TYPE x, TYPE y ) { data[row]=x; data[row+2]=y; }
	/// Set a column of the matrix
	void SetColumn( int column, TYPE x, TYPE y ) { data[2*column]=x; data[2*column+1]=y; }
	/// Set the diagonal values of the matrix
	void SetDiagonal( const TYPE &xx, const TYPE &yy ) { data[0]=xx; data[3]=yy; }
	void SetDiagonal( const Point2<TYPE> &p ) { SetDiagonal( p.x, p.y ); }
	void SetDiagonal( const TYPE *values ) { SetDiagonal(values[0],values[1]); }
	/// Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	/// Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity(TYPE v) { SetScale(v); }
	/// Sets a uniform scale matrix
	void SetScale( TYPE uniformScale ) { SetScale(uniformScale,uniformScale); }
	/// Sets a scale matrix
	void SetScale( TYPE scaleX, TYPE scaleY ) { data[0]=scaleX; data[1]=0; data[2]=0; data[3]=scaleY;}
	/// Sets a scale matrix
	void SetScale( const Point2<TYPE> &scale ) { SetScale(scale.x,scale.y); }
	/// Removes the scale component of the matrix
	void SetNoScale() { Point2<TYPE> *p = (Point2<TYPE>*)data; p[0].Normalize(); p[1].Normalize(); }
	/// Set a rotation matrix by angle
	void SetRotation( TYPE angle ) { SetRotation( cySin(angle), cyCos(angle) ); }
	/// Set a rotation matrix by cos and sin of angle
	void SetRotation( TYPE sinAngle, TYPE cosAngle ) { data[0]=cosAngle; data[1]=-sinAngle; data[2]=sinAngle; data[3]=cosAngle; }
	/// Sets the matrix as the tensor product (outer product) of two vectors
	void SetTensorProduct( const Point2<TYPE> &v0, const Point2<TYPE> &v1 )
	{
		for ( int i=0; i<2; i++ ) data[  i] = v0.data[i] * v1.x;
		for ( int i=0; i<2; i++ ) data[2+i] = v0.data[i] * v1.y;
	}

	// Get Row and Column
	Point2<TYPE>  GetRow   ( int row )                  const { return Point2<TYPE>( data[row], data[row+2] ); }
	void          GetRow   ( int row, Point2<TYPE> &p ) const { p.Set( data[row], data[row+1] ); }
	void          GetRow   ( int row, TYPE *values )    const { values[0]=data[row]; values[1]=data[row+2]; }
	Point2<TYPE>  GetColumn( int col )                  const { return Point2<TYPE>( &data[col*2] ); }
	void          GetColumn( int col, Point2<TYPE> &p ) const { p.Set( &data[col*2] ); }
	void          GetColumn( int col, TYPE *values )    const { values[0]=data[col*2]; values[1]=data[col*2+1]; }
	/// Returns the diagonal component of the matrix
	Point2<TYPE>  GetDiagonal()                         const { return Point2<TYPE>( data[0], data[3] ); }			///< Returns the diagonal component of the matrix
	void          GetDiagonal( Point2<TYPE> &p )        const { p.Set( data[0], data[3] ); }						///< Returns the diagonal component of the matrix
	void	      GetDiagonal( TYPE *values )           const { values[0]=data[0]; values[1]=data[3]; }				///< Returns the diagonal component of the matrix

	//////////////////////////////////////////////////////////////////////////
	///@name Overloaded Operators

	// Overloaded comparison operators 
	bool operator == ( const Matrix2 &right ) const { for ( int i=0; i<4; i++ ) if ( data[i] != right.data[i] ) return false; return true; } ///< compare equal
	bool operator != ( const Matrix2 &right ) const { for ( int i=0; i<4; i++ ) if ( data[i] != right.data[i] ) return true; return false; } ///< compare not equal

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
	Point2<TYPE> operator * ( const Point2<TYPE> &p ) const
	{
		/*
		return Point2<TYPE>(	p.x*data[0] + p.y*data[2],
								p.x*data[1] + p.y*data[3] );
		//*/
		/*
		Point2<TYPE> r;
		TYPE pt[4] = { p.x, p.x, p.y, p.y };
		TYPE a[4];
		for ( int i=0; i<4; ++i ) a[i] = pt[i] * data[i];		// return Point2<TYPE>(	p.x*data[0] + p.y*data[2],
		#pragma _CY_IVDEP										// 						p.x*data[1] + p.y*data[3] );
		for ( int i=0; i<2; ++i ) r.data[i] = a[i] + a[i+2];
		return r;
		//*/
		//*
		Point2<TYPE> r;
		TYPE a[2], b[2];
		#pragma _CY_IVDEP
		for ( int i=0; i<2; ++i ) a[i] = p.data[0] * data[i];		// return Point2<TYPE>(	p.x*data[0] + p.y*data[2],
		#pragma _CY_IVDEP
		for ( int i=0; i<2; ++i ) b[i] = p.data[1] * data[i];		// 						p.x*data[1] + p.y*data[3] );
		#pragma _CY_IVDEP
		for ( int i=0; i<2; ++i ) r.data[i] = a[i] + b[i];
		return r;
		//*/
	}
	Matrix2 operator + ( const Matrix2 &right  ) const { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = data[i] + right.data[i]; return buffer; }	///< add two Matrices
	Matrix2 operator - ( const Matrix2 &right  ) const { Matrix2 buffer; for (int i=0; i<4; i++) buffer.data[i] = data[i] - right.data[i]; return buffer; }	///< subtract one Matrix2 from an other
	Matrix2 operator * ( const Matrix2 &right  ) const	///< multiply a matrix with an other
	{
		Matrix2 r;
		r[0] = data[0] * right.data[0] + data[2] * right.data[1];
		r[1] = data[1] * right.data[0] + data[3] * right.data[1];
		r[2] = data[0] * right.data[2] + data[2] * right.data[3];
		r[3] = data[1] * right.data[2] + data[3] * right.data[3];
		return r;
	}

	// Assignment operators
	const Matrix2& operator  = ( const Matrix2 &right ) { CY_MEMCOPY(TYPE,data,right.data,4); return *this; }	
	const Matrix2& operator += ( const Matrix2 &right ) { for (int i=0; i<4; i++) data[i] += right.data[i]; return *this; }	///< add two Matrices modify this
	const Matrix2& operator -= ( const Matrix2 &right ) { for (int i=0; i<4; i++) data[i] -= right.data[i]; return *this; }	///< subtract one Matrix2 from another matrix and modify this matrix
	const Matrix2& operator *= ( const Matrix2 &right ) { *this = operator*(right); return *this; }							///< multiply a matrix with another matrix and modify this matrix
	const Matrix2& operator *= ( const TYPE value )     { for (int i=0; i<4; i++) data[i] *= value;         return *this; }	///< multiply a matrix with a value modify this matrix
	const Matrix2& operator /= ( const TYPE value )     { for (int i=0; i<4; i++) data[i] /= value;         return *this; }	///< divide the matrix by a value modify the this matrix

	//////////////////////////////////////////////////////////////////////////
	///@name Other Public Methods

	void Transpose() { TYPE tmp=data[0]; data[0]=data[3]; data[3]=tmp; }	///< Transpose this matrix
	void GetTranspose( Matrix2 &m ) const									///< return Transpose of this matrix
	{
		m.data[0] = data[0];   m.data[1] = data[2];
		m.data[2] = data[1];   m.data[3] = data[3];
	}
	Matrix2 GetTranspose() const { Matrix2 t; GetTranspose(t); return t; }	///< return Transpose of this matrix

	TYPE GetDeterminant() const { return data[0]*data[3]-data[2]*data[1]; }		///< Get the determinant of this matrix

	void Invert()					///< Invert this matrix
	{
		TYPE det = GetDeterminant();
		TYPE d0 =  data[0] / det;
		data[0] =  data[3] / det;
		data[1] = -data[1] / det;
		data[2] = -data[2] / det;
		data[3] =  d0;
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
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const { return cyAbs(data[0] - TYPE(1)) < tollerance && cyAbs(data[1]) < tollerance && cyAbs(data[2]) < tollerance && cyAbs(data[3] - TYPE(1)) < tollerance; }

	/// Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( TYPE tollerance=TYPE(0.001) ) const { return cyAbs(data[0] - data[2]) < tollerance; }

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
	
#ifdef CY_NONVECTORIZED_POINT3
	friend Matrix3 operator + ( const TYPE value, const Matrix3 &right ) { Matrix3 r; for (int i=0; i<9; i++) r.data[i] = value + right.data[i]; return r; }	///< add a value to a matrix
	friend Matrix3 operator - ( const TYPE value, const Matrix3 &right ) { Matrix3 r; for (int i=0; i<9; i++) r.data[i] = value - right.data[i]; return r; }	///< subtract the matrix from a value
	friend Matrix3 operator * ( const TYPE value, const Matrix3 &right ) { Matrix3 r; for (int i=0; i<9; i++) r.data[i] = value * right.data[i]; return r; }	///< multiple matrix by a value
#else
	friend Matrix3 operator + ( const TYPE value, const Matrix3 &right )	///< add a value to a matrix
	{
		Matrix3 r;
		#pragma _CY_IVDEP
		for (int i=0; i<8; i++) r.data[i] = value + right.data[i];
		r.data[8] = value + right.data[8];
		return r;
	}
	friend Matrix3 operator - ( const TYPE value, const Matrix3 &right )	///< subtract the matrix from a value
	{
		Matrix3 r;
		#pragma _CY_IVDEP
		for (int i=0; i<8; i++) r.data[i] = value - right.data[i];
		r.data[8] = value - right.data[8];
		return r;
	}
	friend Matrix3 operator * ( const TYPE value, const Matrix3 &right )	///< multiple matrix by a value
	{
		Matrix3 r;
		#pragma _CY_IVDEP
		for (int i=0; i<8; i++) r.data[i] = value * right.data[i];
		r.data[8] = value * right.data[8];
		return r;
	}
#endif
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
	Matrix3( const Matrix3 &matrix ) { CY_MEMCOPY(TYPE,data,matrix.data,9); }								///< Copy constructor
	template <typename T> explicit Matrix3<TYPE>( const Matrix3<T> &matrix ) { CY_MEMCONVERT(TYPE,data,matrix.data,9); }		///< Copy constructor for different types
	explicit Matrix3( const TYPE *values ) { Set(values); }													///< Initialize the matrix using an values of 9 values
	explicit Matrix3( const TYPE &v ) { SetScaledIdentity(v); }												///< Initialize the matrix as identity scaled by v
	explicit Matrix3( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z ) { Set(x,y,z); }	///< Initialize the matrix using x,y,z vectors as columns
	explicit Matrix3( const Matrix2<TYPE> &m ) { 
		data[0] = m.data[0]; data[1] = m.data[1]; data[2] = TYPE(0);
		data[3] = m.data[2]; data[4] = m.data[3]; data[5] = TYPE(0);
		data[6] = TYPE(0);   data[7] = TYPE(0);   data[8] = TYPE(1);
	}
	explicit Matrix3( const Matrix34<TYPE> &m );
	explicit Matrix3( const Matrix4<TYPE>  &m );


	//////////////////////////////////////////////////////////////////////////
	///@name Set & Get Methods

	/// Set all the values as zero
	void Zero() { CY_MEMCLEAR(TYPE,data,9); }
	/// Returns true if the matrix is exactly zero
	bool IsZero() const { for ( int i=0; i<9; i++ ) if ( data[i] != 0 ) return false; return true; }
	/// Copies the matrix data to the given values array of size 9
	void Get( TYPE *values ) { CY_MEMCOPY(TYPE,values,data,9); } 
	/// Set matrix using an values of 9 values
	void Set( const TYPE *values ) { CY_MEMCOPY(TYPE,data,values,9); } 
	/// Set matrix using x,y,z vectors as columns
	void Set( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z ) { x.Get(&data[0]); y.Get(&data[3]); z.Get(&data[6]); }
	/// Set a row of the matrix
	void SetRow( int row, TYPE x, TYPE y, TYPE z ) { data[row]=x; data[row+3]=y; data[row+6]=z; }
	/// Set a column of the matrix
	void SetColumn( int column, TYPE x, TYPE y, TYPE z ) { data[3*column]=x; data[3*column+1]=y; data[3*column+2]=z; }
	/// Set the diagonal values of the matrix
	void SetDiagonal( const TYPE &xx, const TYPE &yy, const TYPE &zz ) { data[0]=xx; data[4]=yy; data[8]=zz; }
	void SetDiagonal( const Point3<TYPE> &p ) { SetDiagonal( p.x, p.y, p.z ); }
	void SetDiagonal( const TYPE *values ) { SetDiagonal(values[0],values[1],values[2]); }
	/// Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	/// Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity(TYPE v) { SetScale(v); }
	/// Sets a uniform scale matrix
	void SetScale( TYPE uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	/// Sets a scale matrix
	void SetScale( TYPE scaleX, TYPE scaleY, TYPE scaleZ )
	{
		data[0] = scaleX; data[1] = 0;      data[2]=0;     
		data[3] = 0;      data[4] = scaleY; data[5]=0;     
		data[6] = 0;      data[7] = 0;      data[8]=scaleZ;
	}
	/// Sets a scale matrix
	void SetScale( const Point3<TYPE> &scale ) { SetScale(scale.x,scale.y,scale.z); }
	/// Removes the scale component of the matrix
	void SetNoScale() { Point3<TYPE> *p = (Point3<TYPE>*)data; p[0].Normalize(); p[1].Normalize(); p[2].Normalize(); }
	/// Set as rotation matrix around x axis
	void SetRotationX( TYPE angle ) { SetRotationX( cySin(angle), cyCos(angle) ); }
	/// Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = TYPE(1);   data[1] =  TYPE(0);    data[2] = TYPE(0); 
		data[3] = TYPE(0);   data[4] =  cosAngle;   data[5] = sinAngle;
		data[6] = TYPE(0);   data[7] = -sinAngle;   data[8] = cosAngle;
	}
	/// Set as rotation matrix around y axis
	void SetRotationY( TYPE angle ) { SetRotationY( cySin(angle), cyCos(angle) ); }
	/// Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = cosAngle;   data[1] = TYPE(0);   data[2] = -sinAngle;
		data[3] = TYPE(0);    data[4] = TYPE(1);   data[5] =  TYPE(0); 
		data[6] = sinAngle;   data[7] = TYPE(0);   data[8] =  cosAngle;
	}
	/// Set as rotation matrix around z axis
	void SetRotationZ( TYPE angle ) { SetRotationZ( cySin(angle), cyCos(angle) ); }
	/// Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] =  cosAngle;   data[1] = sinAngle;   data[2] = TYPE(0);
		data[3] = -sinAngle;   data[4] = cosAngle;   data[5] = TYPE(0);
		data[6] =  TYPE(0);	   data[7] = TYPE(0);    data[8] = TYPE(1);
	}
	/// Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
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
	/// Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
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
	/// Set a rotation matrix about the given axis by angle
	void SetRotation( const Point3<TYPE> &axis, TYPE angle ) { SetRotation(axis,cySin(angle),cyCos(angle)); }
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

		data[0] = tx * axis.x + cosAngle;   data[1] = txy + sz;                 data[2] = txz - sy;
		data[3] = txy - sz;                 data[4] = ty * axis.y + cosAngle;   data[5] = tyz + sx;
		data[6] = txz + sy;                 data[7] = tyz - sx;                 data[8] = tz * axis.z + cosAngle;
	}
	/// Set a rotation matrix that sets [from] unit vector to [to] unit vector
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
	/// Sets the matrix as the tensor product (outer product) of two vectors
	void SetTensorProduct( const Point3<TYPE> &v0, const Point3<TYPE> &v1 )
	{
#ifdef CY_NONVECTORIZED_POINT3
		for ( int i=0; i<3; i++ ) data[  i] = v0.data[i] * v1.x;
		for ( int i=0; i<3; i++ ) data[3+i] = v0.data[i] * v1.y;
		for ( int i=0; i<3; i++ ) data[6+i] = v0.data[i] * v1.z;
#else
		TYPE d[10], v[4];
		v0.Get(v);
		for ( int i=0; i<4; i++ ) d[  i] = v[i] * v1.x;
		for ( int i=0; i<4; i++ ) d[3+i] = v[i] * v1.y;
		for ( int i=0; i<4; i++ ) d[6+i] = v[i] * v1.z;
		CY_MEMCOPY(TYPE,data,d,9);
#endif
	}
	/// Matrix representation of the cross product ( a x b)
	void SetCrossProd( const Point3<TYPE> &p ) { data[0]=TYPE(0); data[1]=p.z; data[2]=-p.y; data[3]=-p.z; data[4]=TYPE(0); data[5]=p.x; data[6]=p.y; data[7]=-p.x; data[8]=TYPE(0); }

	// Get Row and Column
	Point3<TYPE>  GetRow   ( int row )                  const { return Point3<TYPE>( data[row], data[row+3], data[row+6] ); }
	void          GetRow   ( int row, Point3<TYPE> &p ) const { p.Set( data[row], data[row+3], data[row+6] ); }
	void          GetRow   ( int row, TYPE *values )    const { values[0]=data[row]; values[1]=data[row+3]; values[2]=data[row+6]; }
	Point3<TYPE>  GetColumn( int col )                  const { return Point3<TYPE>( &data[col*3] ); }
	void          GetColumn( int col, Point3<TYPE> &p ) const { p.Set( &data[col*3] ); }
	void          GetColumn( int col, TYPE *values )    const { values[0]=data[col*3]; values[1]=data[col*3+1]; values[2]=data[col*3+2]; }
	/// Returns the diagonal component of the matrix
	Point3<TYPE>  GetDiagonal()                         const { Point3<TYPE> p; GetDiagonal(p); return p; }
	void          GetDiagonal( Point3<TYPE> &p )        const { GetDiagonal(p.data); }
	void	      GetDiagonal( TYPE *values )           const { values[0]=data[0]; values[1]=data[4]; values[2]=data[8]; }
	/// Converts the 2x2 portion of the matrix into a Matrix2
	Matrix2<TYPE> GetSubMatrix2()                       const { Matrix2<TYPE> m; GetSubMatrix2(m.data); return m; }
	void          GetSubMatrix2( Matrix2<TYPE> &m )     const { GetSubMatrix2(m.data); }
	void          GetSubMatrix2( TYPE *mdata )          const { CY_MEMCOPY(TYPE,mdata,data,2); CY_MEMCOPY(TYPE,mdata+2,data+3,2); }


	//////////////////////////////////////////////////////////////////////////
	///@name Overloaded Operators

	// Overloaded comparison operators 
	bool operator == ( const Matrix3 &right ) const { for ( int i=0; i<9; i++ ) if ( data[i] != right.data[i] ) return false; return true; } ///< compare equal
	bool operator != ( const Matrix3 &right ) const { for ( int i=0; i<9; i++ ) if ( data[i] != right.data[i] ) return true; return false; } ///< compare not equal

	// Overloaded subscript operators
	TYPE&       operator () ( int row, int column )       { return data[ column * 3 + row ]; }	///< subscript operator
	const TYPE& operator () ( int row, int column ) const { return data[ column * 3 + row ]; }	///< constant subscript operator
	TYPE&       operator [] ( int i )       { return data[i]; }	///< subscript operator
	const TYPE& operator [] ( int i ) const { return data[i]; }	///< constant subscript operator
	
#ifdef CY_NONVECTORIZED_POINT3
	// Unary operators
	Matrix3 operator - () const { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = -data[i]; return buffer; }	///< negative matrix

	// Binary operators
	Matrix3      operator * ( const TYPE value )      const { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = data[i] * value; return buffer; }	///< multiple matrix by a value
	Matrix3      operator / ( const TYPE value )      const { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = data[i] / value; return buffer; }	///< divide matrix by a value;
	Point3<TYPE> operator * ( const Point3<TYPE> &p ) const
	{
		TYPE a[3], b[3], c[3];
		Point3<TYPE> rr;
		for ( int i=0; i<3; ++i ) a[i] = p.data[0] * data[  i];		// return Point3<TYPE>(	p.x*data[0] + p.y*data[3] + p.z*data[6], 
		for ( int i=0; i<3; ++i ) b[i] = p.data[1] * data[3+i];		// 						p.x*data[1] + p.y*data[4] + p.z*data[7],
		for ( int i=0; i<3; ++i ) c[i] = p.data[2] * data[6+i];		// 						p.x*data[2] + p.y*data[5] + p.z*data[8] );
		#pragma _CY_IVDEP
		for ( int i=0; i<3; ++i ) rr.data[i] = a[i] + b[i] + c[i];	
		return rr;
	}
	Matrix3 operator + ( const Matrix3 &right  ) const { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = data[i] + right.data[i]; return buffer; }	///< add two Matrices
	Matrix3 operator - ( const Matrix3 &right  ) const { Matrix3 buffer; for (int i=0; i<9; i++) buffer.data[i] = data[i] - right.data[i]; return buffer; }	///< subtract one Matrix3 from an other
#else
	// Unary operators
	Matrix3 operator - () const		///< negative matrix
	{
		Matrix3 buffer;
		#pragma _CY_IVDEP
		for (int i=0; i<8; i++) buffer.data[i] = -data[i];
		buffer.data[8] = -data[8];
		return buffer;
	}

	// Binary operators
	Matrix3 operator * ( const TYPE value ) const	///< multiple matrix by a value
	{
		Matrix3 buffer;
		#pragma _CY_IVDEP
		for (int i=0; i<8; i++) buffer.data[i] = data[i] * value;
		buffer.data[8] = data[8] * value;
		return buffer;
	}
	Matrix3 operator / ( const TYPE value ) const	///< divide matrix by a value;
	{
		Matrix3 buffer;
		#pragma _CY_IVDEP
		for (int i=0; i<8; i++) buffer.data[i] = data[i] / value;
		buffer.data[8] = data[8] / value;
		return buffer;
	}
	Point3<TYPE> operator * ( const Point3<TYPE> &p ) const
	{
		TYPE a[4], b[4], c[5], r[4];
		#pragma _CY_IVDEP
		for ( int i=0; i<4; ++i ) a[i] = p.data[0] * data[  i];		// return Point3<TYPE>(	p.x*data[0] + p.y*data[3] + p.z*data[6], 
		#pragma _CY_IVDEP
		for ( int i=0; i<4; ++i ) b[i] = p.data[1] * data[3+i];		// 						p.x*data[1] + p.y*data[4] + p.z*data[7],
		#pragma _CY_IVDEP
		for ( int i=0; i<4; ++i ) c[i] = p.data[2] * data[5+i];		// 						p.x*data[2] + p.y*data[5] + p.z*data[8] );
		#pragma _CY_IVDEP
		for ( int i=0; i<4; ++i ) r[i] = a[i] + b[i] + c[i+1];	
		return Point3<TYPE>(r);
	}
	Matrix3 operator + ( const Matrix3 &right  ) const		///< add two Matrices
	{
		Matrix3 buffer; 
		#pragma _CY_IVDEP
		for (int i=0; i<8; i++) buffer.data[i] = data[i] + right.data[i];
		buffer.data[8] = data[8] + right.data[8];
		return buffer;
	}
	Matrix3 operator - ( const Matrix3 &right  ) const		///< subtract one Matrix3 from an other
	{
		Matrix3 buffer;
		#pragma _CY_IVDEP
		for (int i=0; i<8; i++) buffer.data[i] = data[i] - right.data[i];
		buffer.data[8] = data[8] - right.data[8];
		return buffer;
	}
#endif
	Matrix3 operator * ( const Matrix3 &right  ) const	///< multiply a matrix with an other
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

	// Assignment operators
	const Matrix3& operator  = ( const Matrix3 &right ) { CY_MEMCOPY(TYPE,data,right.data,9); return *this; }	
	const Matrix3& operator += ( const Matrix3 &right ) { for (int i=0; i<9; i++) data[i] += right.data[i]; return *this; }	///< add two Matrices modify this
	const Matrix3& operator -= ( const Matrix3 &right ) { for (int i=0; i<9; i++) data[i] -= right.data[i]; return *this; }	///< subtract one Matrix3 from another matrix and modify this matrix
	const Matrix3& operator *= ( const Matrix3 &right ) { *this = operator*(right); return *this; }							///< multiply a matrix with another matrix and modify this matrix
	const Matrix3& operator *= ( const TYPE value )     { for (int i=0; i<9; i++) data[i] *= value;         return *this; }	///< multiply a matrix with a value modify this matrix
	const Matrix3& operator /= ( const TYPE value )     { for (int i=0; i<9; i++) data[i] /= value;         return *this; }	///< divide the matrix by a value modify the this matrix

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
		m.data[0] = data[0];   m.data[1] = data[3];   m.data[2] = data[6];
		m.data[3] = data[1];   m.data[4] = data[4];   m.data[5] = data[7];
		m.data[6] = data[2];   m.data[7] = data[5];   m.data[8] = data[8];
	}
	Matrix3 GetTranspose() const { Matrix3 t; GetTranspose(t); return t; }	///< return Transpose of this matrix

	TYPE GetDeterminant() const {	///< Get the determinant of this matrix
		// 0 (4 8 - 5 7) + 1 (5 6 - 3 8) + 2 (3 7 - 4 6)
		return data[0] * ( data[4] * data[8] - data[5] * data[7] ) + 
		       data[1] * ( data[5] * data[6] - data[3] * data[8] ) + 
		       data[2] * ( data[3] * data[7] - data[4] * data[6] );
	}

	void Invert() { Matrix3 inv; GetInverse(inv); *this=inv; }					///< Invert this matrix
	void GetInverse( Matrix3 &inverse ) const									///< Get the inverse of this matrix
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
		return cyAbs(data[0]-TYPE(1)) < tollerance && cyAbs(data[1])         < tollerance && cyAbs(data[2])         < tollerance && 
			   cyAbs(data[3])         < tollerance && cyAbs(data[4]-TYPE(1)) < tollerance && cyAbs(data[5])         < tollerance &&
			   cyAbs(data[6])         < tollerance && cyAbs(data[7])         < tollerance && cyAbs(data[8]-TYPE(1)) < tollerance;
	}

	/// Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( TYPE tollerance=TYPE(0.001) ) const { return cyAbs(data[1] - data[3]) < tollerance && cyAbs(data[2] - data[6]) < tollerance && cyAbs(data[5] - data[7]) < tollerance; }

	
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
	/// Returns the matrix representation of cross product ( a x b )
	static Matrix3 MatrixCrossProd( const Point3<TYPE> &a ) { Matrix3 m; m.SetCrossProd(a); return m; }

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
	
	friend Matrix34 operator + ( const TYPE value, const Matrix34 &right ) { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = value + right.data[i]; return buffer; }	///< add a value to a matrix
	friend Matrix34 operator - ( const TYPE value, const Matrix34 &right ) { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = value - right.data[i]; return buffer; }	///< subtract the matrix from a value
	friend Matrix34 operator * ( const TYPE value, const Matrix34 &right ) { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = value * right.data[i]; return buffer; }	///< multiple matrix by a value
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
	Matrix34( const Matrix34 &matrix ) { CY_MEMCOPY(TYPE,data,matrix.data,12); }				///< Copy constructor
	template <typename T> explicit Matrix34<TYPE>( const Matrix34<T> &matrix ) { CY_MEMCONVERT(TYPE,data,matrix.data,12); }		///< Copy constructor for different types
	explicit Matrix34( const TYPE *values ) { Set(values); }									///< Initialize the matrix using an values of 9 values
	explicit Matrix34( const TYPE &v ) { SetScaledIdentity(v); }										///< Initialize the matrix as identity scaled by v
	explicit Matrix34( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z, const Point3<TYPE> &pos ) { Set(x,y,z,pos); }	///< Initialize the matrix using x,y,z vectors and coordinate center
	explicit Matrix34( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Set(pos,normal,dir); }				///< Initialize the matrix using position, normal, and approximate x direction
	explicit Matrix34( const Matrix3<TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,9); CY_MEMCLEAR(TYPE,data+9,3); }
	explicit Matrix34( const Matrix3<TYPE> &m, const Point3<TYPE> &pos ) { CY_MEMCOPY(TYPE,data,m.data,9); CY_MEMCOPY(TYPE,data+9,pos.data,3); }
	explicit Matrix34( const Matrix2<TYPE> &m ) { 
		data[ 0] = m.data[ 0]; data[ 1] = m.data[ 1]; data[ 2] = TYPE(0);
		data[ 3] = m.data[ 2]; data[ 4] = m.data[ 3]; data[ 5] = TYPE(0);
		data[ 6] = TYPE(0);    data[ 7] = TYPE(0);    data[ 8] = TYPE(1);
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	explicit Matrix34( const Matrix4<TYPE> &m );


	//////////////////////////////////////////////////////////////////////////
	///@name Set & Get Methods

	/// Set all the values as zero
	void Zero() { CY_MEMCLEAR(TYPE,data,12); }
	/// Returns true if the matrix is exactly zero
	bool IsZero() const { for ( int i=0; i<12; i++ ) if ( data[i] != 0 ) return false; return true; }
	/// Copies the matrix data to the given values array of size 12
	void Get( TYPE *values ) { CY_MEMCOPY(TYPE,values,data,12); } 
	/// Set Matrix using an values of 12 values
	void Set( const TYPE *values ) { CY_MEMCOPY(TYPE,data,values,12); } 
	/// Set matrix using x,y,z vectors and coordinate center
	void Set( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z, const Point3<TYPE> &pos ) { x.Get(data); y.Get(data+3); z.Get(data+6); pos.Get(data+9); }
	/// Set matrix using position, normal, and approximate x direction
	void Set( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,pos); }
	/// Set a row of the matrix
	void SetRow( int row, TYPE x, TYPE y, TYPE z, TYPE w ) { data[row]=x; data[row+3]=y; data[row+6]=z; data[row+9]=w; }
	/// Set a column of the matrix
	void SetColumn( int column, TYPE x, TYPE y, TYPE z ) { data[3*column]=x; data[3*column+1]=y; data[3*column+2]=z; }
	/// Set the diagonal values of the matrix
	void SetDiagonal( const TYPE &xx, const TYPE &yy, const TYPE &zz ) { data[0]=xx; data[4]=yy; data[8]=zz; }
	void SetDiagonal( const Point3<TYPE> &p ) { SetDiagonal( p.x, p.y, p.z ); }
	void SetDiagonal( const TYPE *values ) { SetDiagonal(values[0],values[1],values[2]); }
	/// Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	/// Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity( TYPE v ) { SetScale(v); }
	/// Sets a uniform scale matrix
	void SetScale( TYPE uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	/// Sets a scale matrix
	void SetScale( TYPE scaleX, TYPE scaleY, TYPE scaleZ )
	{
		data[ 0] = scaleX; data[ 1] = 0;      data[ 2]=0;     
		data[ 3] = 0;      data[ 4] = scaleY; data[ 5]=0;     
		data[ 6] = 0;      data[ 7] = 0;      data[ 8]=scaleZ;
		data[ 9] = 0;      data[10] = 0;      data[11]=0;
	}
	/// Sets a scale matrix
	void SetScale( const Point3<TYPE> &scale ) { SetScale(scale.x,scale.y,scale.z); }
	/// Removes the scale component of the matrix
	void SetNoScale() { Point3<TYPE> *p = (Point3<TYPE>*)data; p[0].Normalize(); p[1].Normalize(); p[2].Normalize(); }
	/// Set as rotation matrix around x axis
	void SetRotationX( TYPE angle ) { SetRotationX( cySin(angle), cyCos(angle) ); }
	/// Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = TYPE(1);   data[1] = TYPE(0);     data[2] = TYPE(0); 
		data[3] = TYPE(0);   data[4] = cosAngle;    data[5] = sinAngle;
		data[6] = TYPE(0);   data[7] = -sinAngle;   data[8] = cosAngle;
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	/// Set as rotation matrix around y axis
	void SetRotationY( TYPE angle ) { SetRotationY( cySin(angle), cyCos(angle) ); }
	/// Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] = cosAngle;   data[1] = TYPE(0);   data[2] = -sinAngle;  
		data[3] = TYPE(0);    data[4] = TYPE(1);   data[5] = TYPE(0);  
		data[6] = sinAngle;   data[7] = TYPE(0);   data[8] = cosAngle; 
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	/// Set as rotation matrix around z axis
	void SetRotationZ( TYPE angle ) { SetRotationZ( cySin(angle), cyCos(angle) ); }
	/// Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( TYPE sinAngle, TYPE cosAngle )
	{
		data[0] =  cosAngle;   data[1] = sinAngle;   data[2] = TYPE(0);
		data[3] = -sinAngle;   data[4] = cosAngle;   data[5] = TYPE(0); 
		data[6] =  TYPE(0);	   data[7] = TYPE(0);    data[8] = TYPE(1);
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	/// Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
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
	/// Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
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
	/// Set a rotation matrix about the given axis by angle
	void SetRotation( const Point3<TYPE> &axis, TYPE angle ) { SetRotation(axis,cySin(angle),cyCos(angle)); }
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
		data[ 0] = tx * axis.x + cosAngle;   data[ 1] = txy + sz;                 data[ 2] = txz - sy;
		data[ 3] = txy - sz;                 data[ 4] = ty * axis.y + cosAngle;   data[ 5] = tyz + sx;
		data[ 6] = txz + sy;                 data[ 7] = tyz - sx;                 data[ 8] = tz * axis.z + cosAngle;
		CY_MEMCLEAR(TYPE,data+9,3);
	}
	/// Set a rotation matrix that sets [from] unit vector to [to] unit vector
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
	/// Sets a translation matrix with no rotation or scale
	void SetTrans( const Point3<TYPE> &move ) { TYPE d[12]={1,0,0, 0,1,0, 0,0,1 }; CY_MEMCOPY(TYPE,data,d,9); CY_MEMCOPY(TYPE,data+9,move.data,3); }
	/// Adds a translation to the matrix
	void AddTrans( const Point3<TYPE> &move ) { for ( int i=0; i<3; ++i ) data[9+i] += move.data[i]; }
	/// Sets the translation component of the matrix
	void SetTransComponent( const Point3<TYPE> &move ) { CY_MEMCOPY(TYPE,data+9,move.data,3); }
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

	// Get Row and Column
	Point4<TYPE>  GetRow   ( int row )                  const { return Point4<TYPE>( data[row], data[row+3], data[row+6], data[row+9] ); }
	void          GetRow   ( int row, Point4<TYPE> &p ) const { p.Set( data[row], data[row+3], data[row+6], data[row+9] ); }
	void          GetRow   ( int row, TYPE *values )    const { values[0]=data[row]; values[1]=data[row+3]; values[2]=data[row+6]; values[3]=data[row+9]; }
	Point3<TYPE>  GetColumn( int col )                  const { return Point3<TYPE>( &data[col*3] ); }
	void          GetColumn( int col, Point3<TYPE> &p ) const { p.Set( &data[col*3] ); }
	void          GetColumn( int col, TYPE *values )    const { values[0]=data[col*3]; values[1]=data[col*3+1]; values[2]=data[col*3+2]; }
	/// Returns the diagonal component of the matrix
	Point3<TYPE>  GetDiagonal()                         const { Point3<TYPE> p; GetDiagonal(p); return p; }
	void          GetDiagonal( Point3<TYPE> &p )        const { GetDiagonal(p.data); }
	void	      GetDiagonal( TYPE *values )           const { values[0]=data[0]; values[1]=data[4]; values[2]=data[8]; }
	/// Converts the 3x3 portion of the matrix into a Matrix3
	Matrix3<TYPE> GetSubMatrix3()                       const { Matrix3<TYPE> m; GetSubMatrix3(m.data); return m; }
	void          GetSubMatrix3( Matrix3<TYPE> &m )     const { GetSubMatrix3(m.data); }
	void          GetSubMatrix3( TYPE *mdata )          const { CY_MEMCOPY(TYPE,mdata,data,9); }
	/// Converts the 2x2 portion of the matrix into a Matrix2
	Matrix2<TYPE> GetSubMatrix2()                       const { Matrix2<TYPE> m; GetSubMatrix2(m.data); return m; }
	void          GetSubMatrix2( Matrix2<TYPE> &m )     const { GetSubMatrix2(m.data); }
	void          GetSubMatrix2( TYPE *mdata )          const { CY_MEMCOPY(TYPE,mdata,data,2); CY_MEMCOPY(TYPE,mdata+2,data+3,2); }
	/// Returns the translation component of the matrix
	Point3<TYPE>  GetTrans()                            const { Point3<TYPE> p; GetTrans(p); return p; }
	void          GetTrans( Point3<TYPE> &p )           const { p.x=data[9]; p.y=data[10]; p.z=data[11]; }
	void          GetTrans( TYPE *trans )               const { CY_MEMCOPY(TYPE,trans,data+9,3); }

	//////////////////////////////////////////////////////////////////////////
	///@name Overloaded Operators

	// Overloaded comparison operators 
	bool operator == ( const Matrix34 &right ) const { for ( int i=0; i<12; i++ ) if ( data[i] != right.data[i] ) return false; return true; } ///< compare equal
	bool operator != ( const Matrix34 &right ) const { for ( int i=0; i<12; i++ ) if ( data[i] != right.data[i] ) return true; return false; } ///< compare not equal

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
	Point3<TYPE> operator * ( const Point3<TYPE> &p ) const
	{
		TYPE a[4], b[4], c[4];
		for ( int i=0; i<3; ++i ) a[i] = p.data[0] * data[  i];		// return Point3<TYPE>(	p.x*data[0] + p.y*data[3] + p.z*data[6] + data[ 9], 
		for ( int i=0; i<3; ++i ) b[i] = p.data[1] * data[3+i];		// 						p.x*data[1] + p.y*data[4] + p.z*data[7] + data[10],
		for ( int i=0; i<3; ++i ) c[i] = p.data[2] * data[6+i];		// 						p.x*data[2] + p.y*data[5] + p.z*data[8] + data[11] );
		Point3<TYPE> rr;
		for ( int i=0; i<3; ++i ) rr.data[i] = a[i] + b[i] + c[i] + data[12+i];	
		return rr;
	}
	Point4<TYPE> operator * ( const Point4<TYPE> &p ) const
	{
		TYPE a[6], b[6];
		for ( int i=0; i<3; ++i ) a[  i] = p.data[0] * data[  i];	// return Point4<TYPE>(	p.x*data[0] + p.y*data[3] + p.z*data[6] + p.w*data[ 9],
		for ( int i=0; i<3; ++i ) a[3+i] = p.data[1] * data[3+i];	// 						p.x*data[1] + p.y*data[4] + p.z*data[7] + p.w*data[10],
		for ( int i=0; i<3; ++i ) b[  i] = p.data[2] * data[6+i];	// 						p.x*data[2] + p.y*data[5] + p.z*data[8] + p.w*data[11],
		for ( int i=0; i<3; ++i ) b[3+i] = p.data[3] * data[9+i];	// 						0           + 0           + 0           + p.w          );
		for ( int i=0; i<6; ++i ) a[i] += b[i];
		Point4<TYPE> rr;
		for ( int i=0; i<3; ++i ) rr.data[i] = a[i] + a[3+i];
		rr.w = p.w;
		return rr;
	}

	Matrix34 operator + ( const Matrix34 &right ) const { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = data[i] + right.data[i]; return buffer; }	///< add two Matrices
	Matrix34 operator - ( const Matrix34 &right ) const { Matrix34 buffer; for (int i=0; i<12; i++) buffer.data[i] = data[i] - right.data[i]; return buffer; }	///< subtract one Matrix4 from an other
	Matrix34 operator * ( const Matrix34 &right ) const	///< multiply a matrix with an other
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
		for ( int j=0; j<3; ++j ) r.data[9+j] += data[9+j];
		return r;
	}
	Matrix34 operator * ( const Matrix3<TYPE> &right ) const	///< multiply a matrix with an other
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

	// Assignment operators
	const Matrix34& operator  = ( const Matrix34 &right ) { CY_MEMCOPY(TYPE,data,right.data,12); return *this; }	
	const Matrix34& operator += ( const Matrix34 &right ) { for (int i=0; i<12; i++) data[i] += right.data[i]; return *this; }	///< add two Matrices modify this
	const Matrix34& operator -= ( const Matrix34 &right ) { for (int i=0; i<12; i++) data[i] -= right.data[i]; return *this; }	///< subtract one Matrix4 from another matrix and modify this matrix
	const Matrix34& operator *= ( const Matrix34 &right ) { *this = operator*(right); return *this; }							///< multiply a matrix with another matrix and modify this matrix
	const Matrix34& operator *= ( const Matrix3<TYPE> &right ) { *this = operator*(right); return *this; }						///< multiply a matrix with another matrix and modify this matrix
	const Matrix34& operator *= ( const TYPE value )      { for (int i=0; i<12; i++) data[i] *= value;         return *this; }	///< multiply a matrix with a value modify this matrix
	const Matrix34& operator /= ( const TYPE value )      { for (int i=0; i<12; i++) data[i] /= value;         return *this; }	///< divide the matrix by a value modify the this matrix

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
		CY_MEMCLEAR(TYPE,data+9,3);
	}

	/// return Transpose of this matrix
	void GetTranspose( Matrix34 &m ) const
	{
		m.data[0] = data[0];  m.data[1] = data[3];  m.data[2] = data[6];
		m.data[3] = data[1];  m.data[4] = data[4];  m.data[5] = data[7];
		m.data[6] = data[2];  m.data[7] = data[5];  m.data[8] = data[8];
		CY_MEMCLEAR(TYPE,m.data+9,3);
	}
	Matrix34 GetTranspose() const { Matrix34 t; GetTranspose(t); return t; }	///< return Transpose of this matrix

	TYPE GetDeterminant() const	///< Get the determinant of this matrix
	{
		// 0 (4 8 - 5 7) + 1 (5 6 - 3 8) + 2 (3 7 - 4 6)
		return data[0] * ( data[4] * data[8] - data[5] * data[7] ) + 
		       data[1] * ( data[5] * data[6] - data[3] * data[8] ) + 
		       data[2] * ( data[3] * data[7] - data[4] * data[6] );
	}
	void Invert() { Matrix34 inv; GetInverse(inv); *this=inv; }	///< Invert this matrix
	void GetInverse( Matrix34 &inverse ) const						///< Get the inverse of this matrix
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

	/// Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const
	{
		return cyAbs(data[0]-TYPE(1)) < tollerance && cyAbs(data[ 1])         < tollerance && cyAbs(data[ 2])         < tollerance && 
			   cyAbs(data[3])         < tollerance && cyAbs(data[ 4]-TYPE(1)) < tollerance && cyAbs(data[ 5])         < tollerance &&
			   cyAbs(data[6])         < tollerance && cyAbs(data[ 7])         < tollerance && cyAbs(data[ 8]-TYPE(1)) < tollerance &&
			   cyAbs(data[9])         < tollerance && cyAbs(data[10])         < tollerance && cyAbs(data[11])         < tollerance;
	}

	/// Returns if the matrix is symmetric within the given error tollerance.
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
	friend Matrix4 Inverse( Matrix4 &m ) { return m.GetInverse(); }	///< return the inverse of the matrix

	friend Matrix4 operator * ( const Matrix34<TYPE> &left, const Matrix4 &right  )	///< multiply a matrix with an other
	{
		Matrix4 r;
		TYPE *rd = r.data;
		for ( int i=0; i<16; i+=4, rd+=4 ) {
			TYPE a[3], b[3], c[3], d[3];
			for ( int j=0; j<3; ++j ) a[j] = left.data[  j] * right.data[i  ];
			for ( int j=0; j<3; ++j ) b[j] = left.data[3+j] * right.data[i+1];
			for ( int j=0; j<3; ++j ) c[j] = left.data[6+j] * right.data[i+2];
			for ( int j=0; j<3; ++j ) d[j] = left.data[9+j] * right.data[i+3];
			#pragma _CY_IVDEP
			for ( int j=0; j<3; ++j ) rd[j] = a[j] + b[j] + c[j] + d[j];
		}
		r.data[ 3] = right.data[ 3];
		r.data[ 7] = right.data[ 7];
		r.data[11] = right.data[11];
		r.data[15] = right.data[15];
		return r;
	}

public:

	/// Elements of the matrix are column-major:
	/// | 0   4   8  12 |
	/// | 1   5   9  13 |
	/// | 2   6  10  14 |
	/// | 3   7  11  15 |
	TYPE data[16];


	//////////////////////////////////////////////////////////////////////////
	///@name Constructors

	Matrix4() {}																	///< Default constructor
	Matrix4( const Matrix4 &matrix ) { CY_MEMCOPY(TYPE,data,matrix.data,16); }		///< Copy constructor
	template <typename T> explicit Matrix4<TYPE>( const Matrix4<T> &matrix ) { CY_MEMCONVERT(TYPE,data,matrix.data,16); }		///< Copy constructor for different types
	explicit Matrix4( const TYPE *values ) { Set(values); }							///< Initialize the matrix using an values of 9 values
	explicit Matrix4( const TYPE &v ) { SetScaledIdentity(v); }						///< Initialize the matrix as identity scaled by v
	explicit Matrix4( const Point3<TYPE> &x,   const Point3<TYPE> &y,      const Point3<TYPE> &z, const Point3<TYPE> &pos ) { Set(x,y,z,pos); }	///< Initialize the matrix using x,y,z vectors and coordinate center
	explicit Matrix4( const Point4<TYPE> &x,   const Point4<TYPE> &y,      const Point4<TYPE> &z, const Point4<TYPE> &w   ) { Set(x,y,z,w);   }	///< Initialize the matrix using x,y,z vectors as columns
	explicit Matrix4( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Set(pos,normal,dir); }					///< Initialize the matrix using position, normal, and approximate x direction
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
	///@name Set & Get Methods

	/// Set all the values as zero
	void Zero() { CY_MEMCLEAR(TYPE,data,16); }
	/// Returns true if the matrix is exactly zero
	bool IsZero() const { for ( int i=0; i<16; i++ ) if ( data[i] != 0 ) return false; return true; }
	/// Copies the matrix data to the given values array of size 16
	void Get( TYPE *values ) { CY_MEMCOPY(TYPE,values,data,16); } 
	/// Set Matrix using an values of 16 values
	void Set( const TYPE *values ) { CY_MEMCOPY(TYPE,data,values,16); } 
	/// Set matrix using x,y,z vectors and coordinate center
	void Set( const Point3<TYPE> &x, const Point3<TYPE> &y, const Point3<TYPE> &z, const Point3<TYPE> &pos ) { x.Get(data); data[3]=TYPE(0); y.Get(data+4); data[7]=TYPE(0); z.Get(data+8); data[11]=TYPE(0); pos.Get(data+12); data[15]=TYPE(1); }
	/// Set matrix using x,y,z,w vectors
	void Set( const Point4<TYPE> &x, const Point4<TYPE> &y, const Point4<TYPE> &z, const Point4<TYPE> &w ) { x.Get(data); y.Get(data+4); z.Get(data+8); w.Get(data+12); }
	/// Set matrix using position, normal, and approximate x direction
	void Set( const Point3<TYPE> &pos, const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,pos); }
	/// Set a row of the matrix
	void SetRow( int row, TYPE x, TYPE y, TYPE z, TYPE w ) { data[row]=x; data[row+4]=y; data[row+8]=z; data[row+12]=w; }
	/// Set a column of the matrix
	void SetColumn( int column, TYPE x, TYPE y, TYPE z, TYPE w ) { data[4*column]=x; data[4*column+1]=y; data[4*column+2]=z; data[4*column+3]=w; }
	/// Set the diagonal values of the matrix
	void SetDiagonal( const TYPE &xx, const TYPE &yy, const TYPE &zz, const TYPE &ww=1 ) { data[0]=xx; data[5]=yy; data[10]=zz; data[15]=ww; }
	void SetDiagonal( const Point4<TYPE> &p ) { SetDiagonal( p.x, p.y, p.z, p.w ); }
	void SetDiagonal( const Point3<TYPE> &p ) { SetDiagonal( p.x, p.y, p.z, TYPE(1) ); }
	void SetDiagonal( const TYPE *values ) { SetDiagonal(values[0],values[1],values[2],values[3]); }
	/// Converts the matrix to an identity matrix
	void SetIdentity() { SetScaledIdentity(TYPE(1)); }
	/// Converts the matrix to an identity matrix scaled by a scalar
	void SetScaledIdentity(TYPE v) { SetScale(v); }
	/// Sets a uniform scale matrix
	void SetScale( TYPE uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	/// Sets a scale matrix
	void SetScale( TYPE scaleX, TYPE scaleY, TYPE scaleZ, TYPE scaleW=1 )
	{
		data[ 0] = scaleX; data[ 1] = 0;      data[ 2]=0;      data[ 3]=0; 
		data[ 4] = 0;      data[ 5] = scaleY; data[ 6]=0;      data[ 7]=0; 
		data[ 8] = 0;      data[ 9] = 0;      data[10]=scaleZ; data[11]=0;
		data[12] = 0;      data[13] = 0;      data[14]=0;      data[15]=scaleW;
	}
	/// Sets a scale matrix
	void SetScale( const Point3<TYPE> &scale ) { SetScale(scale.x,scale.y,scale.z); }
	/// Removes the scale component of the matrix
	void SetNoScale() { ((Point3<TYPE>*)&data[0])->Normalize(); ((Point3<TYPE>*)&data[4])->Normalize(); ((Point3<TYPE>*)&data[8])->Normalize(); }
	/// Set as rotation matrix around x axis
	void SetRotationX( TYPE angle ) { SetRotationX( cySin(angle), cyCos(angle) ); }
	/// Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( TYPE sinAngle, TYPE cosAngle )
	{
		data[ 0] = TYPE(1);  data[ 1] =  TYPE(0);    data[ 2] = TYPE(0);   data[ 3] = TYPE(0);
		data[ 4] = TYPE(0);  data[ 5] =  cosAngle;   data[ 6] = sinAngle;  data[ 7] = TYPE(0);
		data[ 8] = TYPE(0);  data[ 9] = -sinAngle;   data[10] = cosAngle;
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15] = TYPE(1);
	}
	/// Set as rotation matrix around y axis
	void SetRotationY( TYPE angle ) { SetRotationY( cySin(angle), cyCos(angle) ); }
	/// Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( TYPE sinAngle, TYPE cosAngle )
	{
		data[ 0] = cosAngle;  data[ 1] = TYPE(0);  data[ 2] = -sinAngle;  data[ 3] = TYPE(0);
		data[ 4] = TYPE(0);   data[ 5] = TYPE(1);  data[ 6] =  TYPE(0);   data[ 7] = TYPE(0);
		data[ 8] = sinAngle;  data[ 9] = TYPE(0);  data[10] =  cosAngle;
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15] = TYPE(1);
	}
	/// Set as rotation matrix around z axis
	void SetRotationZ( TYPE angle ) { SetRotationZ( cySin(angle), cyCos(angle) ); }
	/// Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( TYPE sinAngle, TYPE cosAngle )
	{
		data[ 0] =  cosAngle;  data[ 1] = sinAngle;  data[ 2] = TYPE(0);  data[ 3] = TYPE(0);
		data[ 4] = -sinAngle;  data[ 5] = cosAngle;  data[ 6] = TYPE(0);  data[ 7] = TYPE(0); 
		data[ 8] =  TYPE(0);   data[ 9] = TYPE(0);   data[10] = TYPE(1);
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15] = TYPE(1);
	}
	/// Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
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
	/// Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
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
	/// Set a rotation matrix about the given axis by angle
	void SetRotation( const Point3<TYPE> &axis, TYPE angle ) { SetRotation(axis,cySin(angle),cyCos(angle)); }
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
		data[ 0] = tx * axis.x + cosAngle;  data[ 1] = txy + sz;                data[ 2] = txz - sy;                data[ 3] = TYPE(0);
		data[ 4] = txy - sz;                data[ 5] = ty * axis.y + cosAngle;  data[ 6] = tyz + sx;                data[ 7] = TYPE(0);
		data[ 8] = txz + sy;                data[ 9] = tyz - sx;                data[10] = tz * axis.z + cosAngle;
		CY_MEMCLEAR(TYPE,data+11,4);
		data[15] = TYPE(1);
	}
	/// Set a rotation matrix that sets [from] unit vector to [to] unit vector
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
	/// Sets a translation matrix with no rotation or scale
	void SetTrans( const Point3<TYPE> &move ) { TYPE d[12]={1,0,0,0, 0,1,0,0, 0,0,1,0}; CY_MEMCOPY(TYPE,data,d,12); CY_MEMCOPY(TYPE,data+12,move.data,3); data[15]=TYPE(1); }
	/// Adds a translation to the matrix
	void AddTrans( const Point3<TYPE> &move ) { for ( int i=0; i<3; i++ ) data[12+i] += move.data[i]; }
	/// Sets the translation component of the matrix
	void SetTransComponent( const Point3<TYPE> &move ) { CY_MEMCOPY(TYPE,data+12,move.data,3); }
	/// Set view matrix using position, target and approximate up vector
	void SetView( const Point3<TYPE> &pos, const Point3<TYPE> &target, const Point3<TYPE> &up )
	{
		Point3<TYPE> f = target - pos;
		f.Normalize();
		Point3<TYPE> s = f.Cross(up);
		s.Normalize();
		Point3<TYPE> u = s.Cross(f);
		Matrix34<TYPE> m;
		m.data[ 0]=s.x; m.data[ 1]=u.x; m.data[ 2]=-f.x; m.data[ 3]=TYPE(0);
		m.data[ 4]=s.y; m.data[ 5]=u.y; m.data[ 6]=-f.y; m.data[ 7]=TYPE(0);
		m.data[ 8]=s.z; m.data[ 9]=u.z; m.data[10]=-f.z; m.data[11]=TYPE(0);
		Matrix4 t;
		t.SetTrans(-pos);
		*this = m * t;
	}
	/// Set matrix using normal and approximate x direction
	void SetNormal(const Point3<TYPE> &normal, const Point3<TYPE> &dir ) { Point3<TYPE> y=normal.Cross(dir); y.Normalize(); Point3<TYPE> newdir=y.Cross(normal); Set(newdir,y,normal,Point3<TYPE>(TYPE(0),TYPE(0),TYPE(0))); }
	/// Set a project matrix with field of view in radians
	void SetPerspective( TYPE fov, TYPE aspect, TYPE znear, TYPE zfar ) { SetPerspectiveTan(cyTan(fov*TYPE(0.5)),aspect,znear,zfar); }
	/// Set a project matrix with the tangent of the half field of view (tan_fov_2)
	void SetPerspectiveTan( TYPE tan_fov_2, TYPE aspect, TYPE znear, TYPE zfar )
	{
		const TYPE yScale = TYPE(1) / tan_fov_2;
		const TYPE xScale = yScale / aspect;
		const TYPE zdif = znear - zfar;
		TYPE d[16] = { xScale,0,0,0,  0,yScale,0,0,  0,0,(zfar+znear)/zdif,-1,  0,0,(2*zfar*znear)/zdif,0 };
		CY_MEMCOPY(TYPE,data,d,16);
	}
	/// Sets the matrix as the tensor product (outer product) of two vectors
	void SetTensorProduct( const Point4<TYPE> &v0, const Point4<TYPE> &v1 )
	{
		for ( int i=0; i<4; ++i ) data[   i] = v0.data[i] * v1.x;	 // data[0]=v0.x*v1.x;  data[4]=v0.x*v1.y;  data[ 8]=v0.x*v1.z;  data[12]=v0.x*v1.w;
		for ( int i=0; i<4; ++i ) data[ 4+i] = v0.data[i] * v1.y;	 // data[1]=v0.y*v1.x;  data[5]=v0.y*v1.y;  data[ 9]=v0.y*v1.z;  data[13]=v0.y*v1.w;
		for ( int i=0; i<4; ++i ) data[ 8+i] = v0.data[i] * v1.z;	 // data[2]=v0.z*v1.x;  data[6]=v0.z*v1.y;  data[10]=v0.z*v1.z;  data[14]=v0.z*v1.w;
		for ( int i=0; i<4; ++i ) data[12+i] = v0.data[i] * v1.w;	 // data[3]=v0.w*v1.x;  data[7]=v0.w*v1.y;  data[11]=v0.w*v1.z;  data[15]=v0.w*v1.w;
	}

	// Get Row and Column
	Point4<TYPE>   GetRow   ( int row )                  const { return Point4<TYPE>( data[row], data[row+4], data[row+8], data[row+12] ); }
	void           GetRow   ( int row, Point4<TYPE> &p ) const { p.Set( data[row], data[row+4], data[row+8], data[row+12] ); }
	void           GetRow   ( int row, TYPE *values )    const { values[0]=data[row]; values[1]=data[row+4]; values[2]=data[row+8]; values[3]=data[row+12]; }
	Point4<TYPE>   GetColumn( int col )                  const { return Point4<TYPE>( &data[col*4] ); }
	void           GetColumn( int col, Point4<TYPE> &p ) const { p.Set( &data[col*4] ); }
	void           GetColumn( int col, TYPE *values )    const { values[0]=data[col*4]; values[1]=data[col*4+1]; values[2]=data[col*4+2]; values[3]=data[col*4+3]; }
	/// Returns the diagonal component of the matrix
	Point4<TYPE>   GetDiagonal()                         const { Point4<TYPE> r; GetDiagonal(r.data); return r; }
	void           GetDiagonal( Point4<TYPE> &p )        const { GetDiagonal(p.data); }
	void	       GetDiagonal( TYPE *values )           const { values[0]=data[0]; values[1]=data[5]; values[2]=data[10]; values[3]=data[15]; }
	/// Converts the 3x4 portion of the matrix into a Matrix34
	Matrix34<TYPE> GetSubMatrix34()                      const { Matrix34<TYPE> m; GetSubMatrix34(m.data); return m; }
	void           GetSubMatrix34( Matrix34<TYPE> &m )   const { GetSubMatrix34(m.data); }
	void           GetSubMatrix34( TYPE *mdata )         const { CY_MEMCOPY(TYPE,mdata,data,3); CY_MEMCOPY(TYPE,mdata+3,data+4,3); CY_MEMCOPY(TYPE,mdata+6,data+8,3); CY_MEMCOPY(TYPE,mdata+9,data+12,3); }
	/// Converts the 3x3 portion of the matrix into a Matrix3
	Matrix3<TYPE>  GetSubMatrix3 ()                      const { Matrix3<TYPE> m; GetSubMatrix3(m.data); return m; }
	void           GetSubMatrix3 ( Matrix3<TYPE> &m )    const { GetSubMatrix3(m.data); }
	void           GetSubMatrix3 ( TYPE *mdata )         const { CY_MEMCOPY(TYPE,mdata,data,3); CY_MEMCOPY(TYPE,mdata+3,data+4,3); CY_MEMCOPY(TYPE,mdata+6,data+8,3); }
	/// Converts the 2x2 portion of the matrix into a Matrix2
	Matrix2<TYPE>  GetSubMatrix2 ()                      const { Matrix2<TYPE> m; GetSubMatrix2(m.data); return m; }
	void           GetSubMatrix2 ( Matrix2<TYPE> &m )    const { GetSubMatrix2(m.data); }
	void           GetSubMatrix2 ( TYPE *mdata )         const { CY_MEMCOPY(TYPE,mdata,data,2); CY_MEMCOPY(TYPE,mdata+2,data+4,2); }
	/// Returns the translation component of the matrix
	Point3<TYPE>   GetTrans()                            const { Point3<TYPE> p; GetTrans(p); return p; }
	void           GetTrans( Point3<TYPE> &p )           const { GetTrans(p.data); }
	void           GetTrans( TYPE *trans )               const { CY_MEMCOPY(TYPE,trans,data+12,3); }


	//////////////////////////////////////////////////////////////////////////
	///@name Overloaded Operators

	// Overloaded comparison operators 
	bool operator == ( const Matrix4 &right ) const { for ( int i=0; i<16; i++ ) if ( data[i] != right.data[i] ) return false; return true; } ///< compare equal
	bool operator != ( const Matrix4 &right ) const { for ( int i=0; i<16; i++ ) if ( data[i] != right.data[i] ) return true; return false; } ///< compare not equal

	// Overloaded subscript operators
	TYPE&       operator () ( int row, int column )       { return data[ column * 4 + row ]; }	///< subscript operator
	const TYPE& operator () ( int row, int column ) const { return data[ column * 4 + row ]; }	///< constant subscript operator
	TYPE&       operator [] ( int i )       { return data[i]; }	///< subscript operator
	const TYPE& operator [] ( int i ) const { return data[i]; }	///< constant subscript operator

	// Unary operators
	Matrix4 operator - () const { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i]=-data[i]; return buffer; }	///< negative matrix

	// Binary operators
	Matrix4      operator * ( const TYPE value )      const { Matrix4 buffer; for (int i=0; i<16; ++i) buffer.data[i] = data[i] * value; return buffer; }	///< multiple matrix by a value
	Matrix4      operator / ( const TYPE value )      const { Matrix4 buffer; for (int i=0; i<16; ++i) buffer.data[i] = data[i] / value; return buffer; }	///< divide matrix by a value;
	Point4<TYPE> operator * ( const Point3<TYPE>& p ) const 
	{
		TYPE a[4], b[4], c[4];
		Point4<TYPE> rr;
		#pragma _CY_IVDEP
		for ( int i=0; i<4; ++i ) a[i] = p.data[0] * data[   i];	// return Point4<TYPE>(	p.x*data[0] + p.y*data[4] + p.z*data[ 8] + data[12], 
		#pragma _CY_IVDEP											// 						p.x*data[1] + p.y*data[5] + p.z*data[ 9] + data[13],
		for ( int i=0; i<4; ++i ) b[i] = p.data[1] * data[ 4+i];	// 						p.x*data[2] + p.y*data[6] + p.z*data[10] + data[14],
		#pragma _CY_IVDEP											// 						p.x*data[3] + p.y*data[7] + p.z*data[11] + data[15] );
		for ( int i=0; i<4; ++i ) c[i] = p.data[2] * data[ 8+i];	
		#pragma _CY_IVDEP											
		for ( int i=0; i<4; ++i ) rr.data[i] = a[i] + b[i] + c[i] + data[12+i];	
		return rr;
	}
	Point4<TYPE> operator * ( const Point4<TYPE>& p ) const 
	{
		TYPE a[8], b[8];
		#pragma _CY_IVDEP
		for ( int i=0; i<4; ++i ) a[  i] = p.data[0] * data[   i];		// return Point4<TYPE>(	p.x*data[0] + p.y*data[4] + p.z*data[ 8] + p.w*data[12],
		#pragma _CY_IVDEP												// 						p.x*data[1] + p.y*data[5] + p.z*data[ 9] + p.w*data[13],
		for ( int i=0; i<4; ++i ) a[4+i] = p.data[1] * data[ 4+i];		// 						p.x*data[2] + p.y*data[6] + p.z*data[10] + p.w*data[14],
		#pragma _CY_IVDEP												// 						p.x*data[3] + p.y*data[7] + p.z*data[11] + p.w*data[15] );
		for ( int i=0; i<4; ++i ) b[  i] = p.data[2] * data[ 8+i];		
		#pragma _CY_IVDEP
		for ( int i=0; i<4; ++i ) b[4+i] = p.data[3] * data[12+i];		
		#pragma _CY_IVDEP
		for ( int i=0; i<8; ++i ) a[i] += b[i];
		Point4<TYPE> rr;
		#pragma _CY_IVDEP
		for ( int i=0; i<4; ++i ) rr.data[i] = a[i] + a[4+i];
		return rr;
	}
	Matrix4 operator + ( const Matrix4 &right  ) const { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i] = data[i] + right.data[i]; return buffer; }	///< add two Matrices
	Matrix4 operator - ( const Matrix4 &right  ) const { Matrix4 buffer; for (int i=0; i<16; i++) buffer.data[i] = data[i] - right.data[i]; return buffer; }	///< subtract one Matrix4 from an other
	Matrix4 operator * ( const Matrix4 &right  ) const	///< multiply a matrix with an other
	{
		Matrix4 r;
		TYPE *rd = r.data;
		for ( int i=0; i<16; i+=4, rd+=4 ) {
			TYPE a[4], b[4], c[4], d[4];
			#pragma _CY_IVDEP
			for ( int j=0; j<4; ++j ) a[j] = data[   j] * right.data[i  ];
			#pragma _CY_IVDEP
			for ( int j=0; j<4; ++j ) b[j] = data[ 4+j] * right.data[i+1];
			#pragma _CY_IVDEP
			for ( int j=0; j<4; ++j ) c[j] = data[ 8+j] * right.data[i+2];
			#pragma _CY_IVDEP
			for ( int j=0; j<4; ++j ) d[j] = data[12+j] * right.data[i+3];
			#pragma _CY_IVDEP
			for ( int j=0; j<4; ++j ) rd[j] = a[j] + b[j] + c[j] + d[j];
		}
		return r;
	}
	Matrix4 operator * ( const Matrix34<TYPE> &right ) const	///< multiply a matrix with an other
	{
		TYPE a[4], b[4], c[4];
		Matrix4 r;
		TYPE *rd = r.data;
		for ( int i=0; i<9; i+=3, rd+=4 ) {
			#pragma _CY_IVDEP
			for ( int k=0; k<4; ++k ) a[k] = data[  k] * right.data[i  ];
			#pragma _CY_IVDEP
			for ( int k=0; k<4; ++k ) b[k] = data[4+k] * right.data[i+1];
			#pragma _CY_IVDEP
			for ( int k=0; k<4; ++k ) c[k] = data[8+k] * right.data[i+2];
			#pragma _CY_IVDEP
			for ( int j=0; j<4; ++j ) rd[j] = a[j] + b[j] + c[j];
		}
		#pragma _CY_IVDEP
		for ( int k=0; k<4; ++k ) a[k] = data[  k] * right.data[ 9];
		#pragma _CY_IVDEP
		for ( int k=0; k<4; ++k ) b[k] = data[4+k] * right.data[10];
		#pragma _CY_IVDEP
		for ( int k=0; k<4; ++k ) c[k] = data[8+k] * right.data[11];
		#pragma _CY_IVDEP
		for ( int j=0; j<4; ++j ) rd[j] = (a[j] + b[j]) + (c[j] + data[12+j]);
		return r;
	}
	Matrix4 operator * ( const Matrix3<TYPE> &right ) const	///< multiply a matrix with an other
	{
		TYPE a[4], b[4], c[4];
		Matrix4 r;
		TYPE *rd = r.data;
		for ( int i=0; i<9; i+=3, rd+=4 ) {
			#pragma _CY_IVDEP
			for ( int k=0; k<4; ++k ) a[k] = data[  k] * right.data[i  ];
			#pragma _CY_IVDEP
			for ( int k=0; k<4; ++k ) b[k] = data[4+k] * right.data[i+1];
			#pragma _CY_IVDEP
			for ( int k=0; k<4; ++k ) c[k] = data[8+k] * right.data[i+2];
			#pragma _CY_IVDEP
			for ( int j=0; j<4; ++j ) rd[j] = a[j] + b[j] + c[j];
		}
		CY_MEMCOPY(TYPE,r.data+12,data+12,4);
		return r;
	}

	// Assignment operators
	const Matrix4& operator  = ( const Matrix4 &right ) { CY_MEMCOPY(TYPE,data,right.data,16); return *this; }	
	const Matrix4& operator += ( const Matrix4 &right ) { for (int i=0; i<16; i++) data[i] += right.data[i]; return *this; }	///< add two Matrices modify this
	const Matrix4& operator -= ( const Matrix4 &right ) { for (int i=0; i<16; i++) data[i] -= right.data[i]; return *this; }	///< subtract one Matrix4 from another matrix and modify this matrix
	const Matrix4& operator *= ( const Matrix4 &right )        { *this = operator*(right); return *this; }						///< multiply a matrix with another matrix and modify this matrix
	const Matrix4& operator *= ( const Matrix34<TYPE> &right ) { *this = operator*(right); return *this; }						///< multiply a matrix with another matrix and modify this matrix
	const Matrix4& operator *= ( const Matrix3<TYPE>  &right ) { *this = operator*(right); return *this; }						///< multiply a matrix with another matrix and modify this matrix
	const Matrix4& operator *= ( const TYPE value )     { for (int i=0; i<16; i++) data[i] *= value;         return *this; }	///< multiply a matrix with a value modify this matrix
	const Matrix4& operator /= ( const TYPE value )     { for (int i=0; i<16; i++) data[i] /= value;         return *this; }	///< divide the matrix by a value modify the this matrix

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
		m.data[ 0] = data[0];   m.data[ 1] = data[4];   m.data[ 2] = data[ 8];  m.data[ 3] = data[12];
		m.data[ 4] = data[1];   m.data[ 5] = data[5];   m.data[ 6] = data[ 9];  m.data[ 7] = data[13];
		m.data[ 8] = data[2];   m.data[ 9] = data[6];   m.data[10] = data[10];  m.data[11] = data[14];
		m.data[12] = data[3];   m.data[13] = data[7];   m.data[14] = data[11];  m.data[15] = data[15];
	}
	Matrix4 GetTranspose() const { Matrix4 t; GetTranspose(t); return t; }	///< return Transpose of this matrix

	TYPE GetDeterminant() const	///< Get the determinant of this matrix
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

	void Invert() { Matrix4 inv; GetInverse(inv); *this=inv; }					///< Invert this matrix
	void GetInverse( Matrix4 &inverse ) const									///< Get the inverse of this matrix
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

	/// Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( TYPE tollerance=TYPE(0.001) ) const
	{
		return cyAbs(data[ 0]-TYPE(1)) < tollerance && cyAbs(data[ 1])         < tollerance && cyAbs(data[ 2])         < tollerance && cyAbs(data[ 3])         < tollerance && 
			   cyAbs(data[ 4])         < tollerance && cyAbs(data[ 5]-TYPE(1)) < tollerance && cyAbs(data[ 6])         < tollerance && cyAbs(data[ 7])         < tollerance &&
			   cyAbs(data[ 8])         < tollerance && cyAbs(data[ 9])         < tollerance && cyAbs(data[10]-TYPE(1)) < tollerance && cyAbs(data[11])         < tollerance &&
			   cyAbs(data[12])         < tollerance && cyAbs(data[13])         < tollerance && cyAbs(data[14])         < tollerance && cyAbs(data[15]-TYPE(1)) < tollerance;
	}

	/// Returns if the matrix is symmetric within the given error tollerance.
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

// Definitions of the conversion constructors
template <typename TYPE>  Matrix2 <TYPE>::Matrix2 ( const Matrix3 <TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,2); CY_MEMCOPY(TYPE,data+2,m.data+3,2); }
template <typename TYPE>  Matrix2 <TYPE>::Matrix2 ( const Matrix34<TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,2); CY_MEMCOPY(TYPE,data+2,m.data+3,2); }
template <typename TYPE>  Matrix2 <TYPE>::Matrix2 ( const Matrix4 <TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,2); CY_MEMCOPY(TYPE,data+2,m.data+4,2); }
template <typename TYPE>  Matrix3 <TYPE>::Matrix3 ( const Matrix34<TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,9); }
template <typename TYPE>  Matrix3 <TYPE>::Matrix3 ( const Matrix4 <TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,3); CY_MEMCOPY(TYPE,data+3,m.data+4,3); CY_MEMCOPY(TYPE,data+6,m.data+8,3); }
template <typename TYPE>  Matrix34<TYPE>::Matrix34( const Matrix4 <TYPE> &m ) { CY_MEMCOPY(TYPE,data,m.data,3); CY_MEMCOPY(TYPE,data+3,m.data+4,3); CY_MEMCOPY(TYPE,data+6,m.data+8,3); CY_MEMCOPY(TYPE,data+9,m.data+12,3); }

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
