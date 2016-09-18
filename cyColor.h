// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
/// \file		cyColor.h 
/// \author		Cem Yuksel
/// \brief		Color classes.
///
/// @copydoc Color
///
/// Color is color class that holds floating point color values
///
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

#ifndef _CY_COLOR_H_INCLUDED_
#define _CY_COLOR_H_INCLUDED_

//-------------------------------------------------------------------------------

#include <math.h>

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

class ColorA;
class Color24;
class Color32;

//-------------------------------------------------------------------------------

/// Color class

class Color
{
	friend Color operator+( const float v, const Color &c ) { return c+v; }		///< Addition with a constant
	friend Color operator-( const float v, const Color &c ) { return -(c-v); }	///< Subtraction from a constant
	friend Color operator*( const float v, const Color &c ) { return c*v; }		///< Multiplication with a constant

public:

	float r, g, b;

	///@name Constructors
	Color() { }
	Color( float rgb ) : r(rgb), g(rgb), b(rgb) {}
	Color( float _r, float _g, float _b ) : r(_r), g(_g), b(_b) {}
	Color( const Color &c ) : r(c.r), g(c.g), b(c.b) {}
	explicit Color( const float *c ) : r(c[0]), g(c[1]), b(c[2]) {}
	explicit Color( const ColorA  &c );
	explicit Color( const Color24 &c );
	explicit Color( const Color32 &c );

	///@name Set & Get value functions
	Color& Black () { r=0; g=0; b=0; return *this; }									///< Sets r, g and b components as zero
	Color& White () { r=1; g=1; b=1; return *this; }									///< Sets r, g and b components as one
	Color& Set   ( float _r, float _g, float _b ) { r=_r; g=_g; b=_b; return *this; }	///< Sets r, g and b components as given
	Color& Set   ( const float *v ) { r=v[0]; g=v[1]; b=v[2]; return *this; }			///< Sets r, g and b components using the values in the given array
	void GetValue( float *v ) const { v[0]=r; v[1]=g; v[2]=b; }							///< Puts r, g and b values into the array

	///@name Gray-scale functions
	float Sum() const   { return r + g + b; }
	float Gray() const  { return Sum() / 3.0f; }
	float Luma1() const { return 0.299f *r + 0.587f *g + 0.114f *b; }
	float Luma2() const { return 0.2126f*r + 0.7152f*g + 0.0722f*b; }

	///@name Limit functions
	void ClampMinMax( float min=0, float max=1 ) { ClampMin(min); ClampMax(max); }
	void ClampMin   ( float n=0 ) { if(r<n)r=n; if(g<n)g=n; if(b<n)b=n; }
	void ClampMax   ( float n=1 ) { if(r>n)r=n; if(g>n)g=n; if(b>n)b=n; }
	void Abs() { r = fabsf(r); g = fabsf(g); b = fabsf(b); }

	///@name Unary operators
	Color  operator - () const { return Color(-r,-g,-b); } 
	Color  operator + () const { return *this; }

	///@name Binary operators
	Color  operator + ( const Color &c ) const { return Color(r+c.r, g+c.g, b+c.b); }
	Color  operator - ( const Color &c ) const { return Color(r-c.r, g-c.g, b-c.b); }
	Color  operator * ( const Color &c ) const { return Color(r*c.r, g*c.g, b*c.b); }
	Color  operator / ( const Color &c ) const { return Color(r/c.r, g/c.g, b/c.b); }
	Color  operator + ( float n ) const { return Color(r+n, g+n, b+n); }
	Color  operator - ( float n ) const { return Color(r-n, g-n, b-n); }
	Color  operator * ( float n ) const { return Color(r*n, g*n, b*n); }
	Color  operator / ( float n ) const { return Color(r/n, g/n, b/n); }

	///@name Assignment operators
	Color& operator += ( const Color &c ) { r+=c.r; g+=c.g; b+=c.b; return *this; }
	Color& operator -= ( const Color &c ) { r-=c.r; g-=c.g; b-=c.b; return *this; }
	Color& operator *= ( const Color &c ) { r*=c.r; g*=c.g; b*=c.b; return *this; }
	Color& operator /= ( const Color &c ) { r/=c.r; g/=c.g; b/=c.b; return *this; }
	Color& operator += ( float n ) { r+=n; g+=n; b+=n; return *this; }
	Color& operator -= ( float n ) { r-=n; g-=n; b-=n; return *this; }
	Color& operator *= ( float n ) { r*=n; g*=n; b*=n; return *this; }
	Color& operator /= ( float n ) { r/=n; g/=n; b/=n; return *this; }

	///@name Test operators
	int operator == ( const Color& c ) const { return ( (c.r==r) && (c.g==g) && (c.b==b) ); }
	int operator != ( const Color& c ) const { return ( (c.r!=r) || (c.g!=g) || (c.b!=b) ); }

	///@name Access operators
	float& operator [] ( int i )       { return (&r)[i]; }
	float  operator [] ( int i ) const { return (&r)[i]; }
};

//-------------------------------------------------------------------------------


/// Color class with alpha
class ColorA
{
	friend ColorA operator+( const float v, const ColorA &c ) { return c+v; }		///< Addition with a constant
	friend ColorA operator-( const float v, const ColorA &c ) { return -(c-v); }	///< Subtraction from a constant
	friend ColorA operator*( const float v, const ColorA &c ) { return c*v; }		///< Multiplication with a constant

public:

	float r, g, b, a;

	///@name Constructors
	ColorA() { }
	ColorA( float rgb, float _a=1 ) : r(rgb), g(rgb), b(rgb), a(_a) {}
	ColorA( float _r, float _g, float _b, float _a=1 ) : r(_r), g(_g), b(_b), a(_a) {}
	ColorA( const ColorA  &c ) : r(c.r), g(c.g), b(c.b), a(c.a) {}
	explicit ColorA( const float *c ) : r(c[0]), g(c[1]), b(c[2]), a(c[3]) {}
	explicit ColorA( const Color   &c, float _a=1 ) : r(c.r), g(c.g), b(c.b), a(_a) {}
	explicit ColorA( const Color24 &c, float _a=1 );
	explicit ColorA( const Color32 &c );

	///@name Set & Get value functions
	ColorA& Black( float alpha=1 ) { r=0; g=0; b=0; a=alpha; return *this; }							///< Sets r, g, and b components as zero and a component as given
	ColorA& White( float alpha=1 ) { r=0; g=0; b=0; a=alpha; return *this; }							///< Sets r, g, and b components as one and a component as given
	ColorA& Set  ( float _r, float _g, float _b, float _a=1 ) { r=_r; g=_g; b=_b; a=_a; return *this; }	///< Sets r, g, b and a components as given
	ColorA& Set  ( const float *v ) { r=v[0]; g=v[1]; b=v[2]; a=v[3]; return *this; }					///< Sets r, g, b and a components using the values in the given array
	void GetValue( float *v ) const { v[0]=r; v[1]=g; v[2]=b; v[3]=a; }									///< Puts r, g, b and a values into the array

	///@name Gray-scale functions
	float Sum  () const { return r + g + b; }
	float Gray () const { return Sum() / 3.0f; }
	float Luma1() const { return 0.299f *r + 0.587f *g + 0.114f *b; }
	float Luma2() const { return 0.2126f*r + 0.7152f*g + 0.0722f*b; }

	///@name Limit functions
	void ClampMinMax( float min, float max ) { ClampMin(min); ClampMax(max); }
	void ClampMin   ( float n ) { if(r<n)r=n; if(g<n)g=n; if(b<n)b=n; if(a<n)a=n; }
	void ClampMax   ( float n ) { if(r>n)r=n; if(g>n)g=n; if(b>n)b=n; if(a>n)a=n; }
	void Abs() { r = fabsf(r); g = fabsf(g); b = fabsf(b); a = fabsf(a); }

	///@name Unary operators
	ColorA  operator - () const { return ColorA(-r,-g,-b,-a); } 
	ColorA  operator + () const { return *this; }

	///@name Binary operators
	ColorA  operator + ( const ColorA &c ) const { return ColorA(r+c.r, g+c.g, b+c.b, a+c.a); }
	ColorA  operator - ( const ColorA &c ) const { return ColorA(r-c.r, g-c.g, b-c.b, a-c.a); }
	ColorA  operator * ( const ColorA &c ) const { return ColorA(r*c.r, g*c.g, b*c.b, a*c.a); }
	ColorA  operator / ( const ColorA &c ) const { return ColorA(r/c.r, g/c.g, b/c.b, a/c.a); }
	ColorA  operator + ( float n ) const { return ColorA(r+n, g+n, b+n, a+n); }
	ColorA  operator - ( float n ) const { return ColorA(r-n, g-n, b-n, a-n); }
	ColorA  operator * ( float n ) const { return ColorA(r*n, g*n, b*n, a*n); }
	ColorA  operator / ( float n ) const { return ColorA(r/n, g/n, b/n, a/n); }

	///@name Assignment operators
	ColorA& operator += ( const ColorA &c ) { r+=c.r; g+=c.g; b+=c.b; a+=c.a; return *this; }
	ColorA& operator -= ( const ColorA &c ) { r-=c.r; g-=c.g; b-=c.b; a-=c.a; return *this; }
	ColorA& operator *= ( const ColorA &c ) { r*=c.r; g*=c.g; b*=c.b; a*=c.a; return *this; }
	ColorA& operator /= ( const ColorA &c ) { r/=c.r; g/=c.g; b/=c.b; a/=c.a; return *this; }
	ColorA& operator += ( float n ) { r+=n; g+=n; b+=n; a+=n; return *this; }
	ColorA& operator -= ( float n ) { r-=n; g-=n; b-=n; a-=n; return *this; }
	ColorA& operator *= ( float n ) { r*=n; g*=n; b*=n; a*=n; return *this; }
	ColorA& operator /= ( float n ) { r/=n; g/=n; b/=n; a/=n; return *this; }

	///@name Test operators
	int operator == ( const ColorA& c ) const { return ( (c.r==r) && (c.g==g) && (c.b==b) && (c.a==a) ); }
	int operator != ( const ColorA& c ) const { return ( (c.r!=r) || (c.g!=g) || (c.b!=b) || (c.a!=a) ); }

	///@name Access operators
	float& operator [] ( int i )       { return (&r)[i]; }
	float  operator [] ( int i ) const { return (&r)[i]; }
};

//-------------------------------------------------------------------------------

/// 24-bit color class
class Color24
{
public:
	unsigned char r, g, b;

	///@name Constructors
	Color24() {}
	Color24( unsigned char _r, unsigned char _g, unsigned char _b ) : r(_r), g(_g), b(_b) {}
	Color24( const Color24 &c ) : r(c.r), g(c.g), b(c.b) {}
	explicit Color24( const Color   &c ) { r=FloatToByte(c.r); g=FloatToByte(c.g); b=FloatToByte(c.b); }
	explicit Color24( const ColorA  &c ) { r=FloatToByte(c.r); g=FloatToByte(c.g); b=FloatToByte(c.b); }
	explicit Color24( const Color32 &c );

	///@name Conversion Methods
	Color ToColor() const { return Color(r/255.0f,g/255.0f,b/255.0f); }

	///@name Set & Get value functions
	Color24& Black() { r=0; g=0; b=0; return *this; }															///< Sets r, g, and b components as zero
	Color24& White() { r=255; g=255; b=255; return *this; }													///< Sets r, g, and b components as 255
	Color24& Set  ( unsigned char _r, unsigned char _g, unsigned char _b ) { r=_r; g=_g; b=_b; return *this; }	///< Sets r, g, and b components as given
	Color24& Set  ( const unsigned char *v ) { r=v[0]; g=v[1]; b=v[2]; return *this; }							///< Sets r, g, and b components using the values in the given array
	void GetValue ( unsigned char *v ) const { v[0]=r; v[1]=g; v[2]=b; }											///< Puts r, g, and b values into the array

protected:
	static unsigned char FloatToByte(float r) { return Clamp(int(r*255)); }
	static unsigned char Clamp(int v) { return v<0 ? 0 : (v>255 ? 255 : v); }
};

//-------------------------------------------------------------------------------

/// 32-bit color class
class Color32
{
public:
	unsigned char r, g, b, a;

	///@name Constructors
	Color32() {}
	Color32( unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a=255 ) : r(_r), g(_g), b(_b), a(_a) {}
	Color32( const Color32 &c ) : r(c.r), g(c.g), b(c.b), a(c.a) {}
	explicit Color32( const Color  &c, float _a=1.0f ) { r=FloatToByte(c.r); g=FloatToByte(c.g); b=FloatToByte(c.b); a=FloatToByte(_a); }
	explicit Color32( const ColorA &c ) { r=FloatToByte(c.r); g=FloatToByte(c.g); b=FloatToByte(c.b); a=FloatToByte(c.a); }
	explicit Color32( const Color24 &c, unsigned char _a=255 ) : r(c.r), g(c.g), b(c.b), a(_a) {}

	///@name Conversion Methods
	Color  ToColor () const { return Color(r/255.0f,g/255.0f,b/255.0f); }
	ColorA ToColorA() const { return ColorA(r/255.0f,g/255.0f,b/255.0f,a/255.0f); }

	///@name Set & Get value functions
	Color32& Black( unsigned char _a=255 ) { r=0; g=0; b=0; a=_a; return *this; }														///< Sets r, g, and b components as zero and a component as given
	Color32& White( unsigned char _a=255 ) { r=255; g=255; b=255; a=_a; return *this; }													///< Sets r, g, and b components as one and a component as given
	Color32& Set  ( unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a ) { r=_r; g=_g; b=_b; a=_a; return *this; }	///< Sets r, g, b and a components as given
	Color32& Set  ( const unsigned char *v ) { r=v[0]; g=v[1]; b=v[2]; a=v[3]; return *this; }											///< Sets r, g, b and a components using the values in the given array
	void GetValue ( unsigned char *v ) const { v[0]=r; v[1]=g; v[2]=b; v[3]=a; }															///< Puts r, g, b and a values into the array

protected:
	static unsigned char FloatToByte(float r) { return Clamp(int(r*255)); }
	static unsigned char Clamp(int v) { return v<0 ? 0 : (v>255 ? 255 : v); }
};

//-------------------------------------------------------------------------------

inline Color  ::Color  ( const ColorA  &c ) : r(c.r), g(c.g), b(c.b) {}
inline Color  ::Color  ( const Color24 &c ) { *this = c.ToColor(); }
inline Color  ::Color  ( const Color32 &c ) { *this = c.ToColor(); }
inline ColorA ::ColorA ( const Color24 &c, float _a ) { Color rgb = c.ToColor(); r = rgb.r; g = rgb.g; b = rgb.b; a = _a; }
inline ColorA ::ColorA ( const Color32 &c ) { *this = c.ToColorA(); }
inline Color24::Color24( const Color32 &c ) : r(c.r), g(c.g), b(c.b) {}

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Color   cyColor;
typedef cy::ColorA  cyColorA;
typedef cy::Color24 cyColor24;
typedef cy::Color32 cyColor32;

//-------------------------------------------------------------------------------

#endif

