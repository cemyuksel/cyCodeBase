// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyColor.h 
//! \author Cem Yuksel
//! 
//! \brief  Color classes.
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

#ifndef _CY_COLOR_H_INCLUDED_
#define _CY_COLOR_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyCore.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

class ColorA;
class Color24;
class Color32;

//-------------------------------------------------------------------------------

//! RGB color class with 3 float components

class Color
{
	friend Color operator+( float const& v, Color const &c ) { return   c+v;  }	//!< Addition with a constant
	friend Color operator-( float const& v, Color const &c ) { return -(c-v); }	//!< Subtraction from a constant
	friend Color operator*( float const& v, Color const &c ) { return   c*v;  }	//!< Multiplication with a constant

public:

	//!@name Color components
	float r, g, b;

	//!@name Constructors
	Color() CY_CLASS_FUNCTION_DEFAULT
	Color( Color const &c ) : r(c.r), g(c.g), b(c.b) {}
	explicit Color( float _r, float _g, float _b ) : r(_r), g(_g), b(_b) {}
	explicit Color( float const *c ) : r(c[0]), g(c[1]), b(c[2]) {}
	explicit Color( float rgb ) : r(rgb), g(rgb), b(rgb) {}
	explicit Color( ColorA  const &c );
	explicit Color( Color24 const &c );
	explicit Color( Color32 const &c );

	//!@name Set & Get value functions
	void SetBlack() { r=0.0f; g=0.0f; b=0.0f; }							//!< Sets r, g and b components as zero
	void SetWhite() { r=1.0f; g=1.0f; b=1.0f; }							//!< Sets r, g and b components as one
	void Set     ( float _r, float _g, float _b ) { r=_r; g=_g; b=_b; }	//!< Sets r, g and b components as given
	void Set     ( float const *v ) { r=v[0]; g=v[1]; b=v[2]; }			//!< Sets r, g and b components using the values in the given array
	void GetValue( float *v ) const { v[0]=r; v[1]=g; v[2]=b; }			//!< Puts r, g and b values into the array

	//!@name Gray-scale functions
	float Sum    () const { return r + g + b; }
	float Gray   () const { return Sum() / 3.0f; }
	float Luma1  () const { return 0.299f *r + 0.587f *g + 0.114f *b; }
	float Luma2  () const { return 0.2126f*r + 0.7152f*g + 0.0722f*b; }
	bool  IsBlack() const { return r==0.0f && g==0.0f && b==0.0f; }			//!< Returns true if all components are exactly zero
	float Min    () const { return r<g ? (r<b ? r : b) : (g<b ? g : b); }
	float Max    () const { return r>g ? (r>b ? r : b) : (g>b ? g : b); }

	//!@name Limit functions
	void Clamp   ( float limitMin=0.0f, float limitMax=1.0f ) { ClampMin(limitMin); ClampMax(limitMax); }
	void ClampMin( float limitMin=0.0f ) { r=cy::Max(r,limitMin); g=cy::Max(g,limitMin); b=cy::Max(b,limitMin); }
	void ClampMax( float limitMax=1.0f ) { r=cy::Min(r,limitMax); g=cy::Min(g,limitMax); b=cy::Min(b,limitMax); }
	void Abs() { r = std::abs(r); g = std::abs(g); b = std::abs(b); }

	//!@name Unary operators
	Color  operator - () const { return Color(-r,-g,-b); } 

	//!@name Binary operators
	Color  operator + ( Color const &c ) const { return Color(r+c.r, g+c.g, b+c.b); }
	Color  operator - ( Color const &c ) const { return Color(r-c.r, g-c.g, b-c.b); }
	Color  operator * ( Color const &c ) const { return Color(r*c.r, g*c.g, b*c.b); }
	Color  operator / ( Color const &c ) const { return Color(r/c.r, g/c.g, b/c.b); }
	Color  operator + ( float const &n ) const { return Color(r+n, g+n, b+n); }
	Color  operator - ( float const &n ) const { return Color(r-n, g-n, b-n); }
	Color  operator * ( float const &n ) const { return Color(r*n, g*n, b*n); }
	Color  operator / ( float const &n ) const { return Color(r/n, g/n, b/n); }

	//!@name Assignment operators
	Color& operator += ( Color const &c ) { r+=c.r; g+=c.g; b+=c.b; return *this; }
	Color& operator -= ( Color const &c ) { r-=c.r; g-=c.g; b-=c.b; return *this; }
	Color& operator *= ( Color const &c ) { r*=c.r; g*=c.g; b*=c.b; return *this; }
	Color& operator /= ( Color const &c ) { r/=c.r; g/=c.g; b/=c.b; return *this; }
	Color& operator += ( float const &n ) { r+=n; g+=n; b+=n; return *this; }
	Color& operator -= ( float const &n ) { r-=n; g-=n; b-=n; return *this; }
	Color& operator *= ( float const &n ) { r*=n; g*=n; b*=n; return *this; }
	Color& operator /= ( float const &n ) { r/=n; g/=n; b/=n; return *this; }

	//!@name Test operators
	bool operator == ( Color const &c ) const { return ( (c.r==r) && (c.g==g) && (c.b==b) ); }
	bool operator != ( Color const &c ) const { return ( (c.r!=r) || (c.g!=g) || (c.b!=b) ); }

	//!@name Access operators
	float& operator [] ( int i )       { return (&r)[i]; }
	float  operator [] ( int i ) const { return (&r)[i]; }

	//!@name Static Methods
	static Color Black() { return Color(0.0f,0.0f,0.0f); }	//!< Returns a black color
	static Color White() { return Color(1.0f,1.0f,1.0f); }	//!< Returns a white color
};

//-------------------------------------------------------------------------------

//! RGBA color class with 4 float components

class ColorA
{
	friend ColorA operator + ( float const &v, ColorA const &c ) { return   c+v;  }	//!< Addition with a constant
	friend ColorA operator - ( float const &v, ColorA const &c ) { return -(c-v); }	//!< Subtraction from a constant
	friend ColorA operator * ( float const &v, ColorA const &c ) { return   c*v;  }	//!< Multiplication with a constant

public:

	//!@name Color components
	float r, g, b, a;

	//!@name Constructors
	ColorA() CY_CLASS_FUNCTION_DEFAULT
	ColorA( ColorA const &c ) : r(c.r), g(c.g), b(c.b), a(c.a) {}
	explicit ColorA( float _r, float _g, float _b, float _a=1 ) : r(_r), g(_g), b(_b), a(_a) {}
	explicit ColorA( float const *c ) : r(c[0]), g(c[1]), b(c[2]), a(c[3]) {}
	explicit ColorA( float rgb, float _a=1 ) : r(rgb), g(rgb), b(rgb), a(_a) {}
	explicit ColorA( Color   const &c, float _a=1 ) : r(c.r), g(c.g), b(c.b), a(_a) {}
	explicit ColorA( Color24 const &c, float _a=1 );
	explicit ColorA( Color32 const &c );

	//!@name Set & Get value functions
	void SetBlack( float alpha=1.0f ) { r=0.0f; g=0.0f; b=0.0f; a=alpha; }					//!< Sets r, g, and b components as zero and a component as given
	void SetWhite( float alpha=1.0f ) { r=0.0f; g=0.0f; b=0.0f; a=alpha; }					//!< Sets r, g, and b components as one and a component as given
	void Set     ( float _r, float _g, float _b, float _a=1 ) { r=_r; g=_g; b=_b; a=_a; }	//!< Sets r, g, b and a components as given
	void Set     ( float const *v ) { r=v[0]; g=v[1]; b=v[2]; a=v[3]; }						//!< Sets r, g, b and a components using the values in the given array
	void GetValue( float *v ) const { v[0]=r; v[1]=g; v[2]=b; v[3]=a; }						//!< Puts r, g, b and a values into the array

	//!@name Gray-scale functions
	float Sum    () const { return r + g + b; }
	float Gray   () const { return Sum() / 3.0f; }
	float Luma1  () const { return 0.299f *r + 0.587f *g + 0.114f *b; }
	float Luma2  () const { return 0.2126f*r + 0.7152f*g + 0.0722f*b; }
	bool  IsBlack() const { return r==0.0f && g==0.0f && b==0.0f; }			//!< Returns true if the r, g, and b components are exactly zero
	float Min    () const { float mrg = r<g ? r : g; float mba = b<a ? b : a; return mrg<mba ? mrg : mba; }
	float Max    () const { float mrg = r>g ? r : g; float mba = b>a ? b : a; return mrg>mba ? mrg : mba; }

	//!@name Limit functions
	void Clamp   ( float limitMin=0.0f, float limitMax=1.0f ) { ClampMin(limitMin); ClampMax(limitMax); }
	void ClampMin( float limitMin=0.0f ) { r=cy::Max(r,limitMin); g=cy::Max(g,limitMin); b=cy::Max(b,limitMin); a=cy::Max(a,limitMin); }
	void ClampMax( float limitMax=1.0f ) { r=cy::Min(r,limitMax); g=cy::Min(g,limitMax); b=cy::Min(b,limitMax); a=cy::Min(a,limitMax); }
	void Abs() { r = std::abs(r); g = std::abs(g); b = std::abs(b); a = std::abs(a); }

	//!@name Unary operators
	ColorA  operator - () const { return ColorA(-r,-g,-b,-a); } 

	//!@name Binary operators
	ColorA  operator + ( ColorA const &c ) const { return ColorA(r+c.r, g+c.g, b+c.b, a+c.a); }
	ColorA  operator - ( ColorA const &c ) const { return ColorA(r-c.r, g-c.g, b-c.b, a-c.a); }
	ColorA  operator * ( ColorA const &c ) const { return ColorA(r*c.r, g*c.g, b*c.b, a*c.a); }
	ColorA  operator / ( ColorA const &c ) const { return ColorA(r/c.r, g/c.g, b/c.b, a/c.a); }
	ColorA  operator + ( float  const &n ) const { return ColorA(r+n, g+n, b+n, a+n); }
	ColorA  operator - ( float  const &n ) const { return ColorA(r-n, g-n, b-n, a-n); }
	ColorA  operator * ( float  const &n ) const { return ColorA(r*n, g*n, b*n, a*n); }
	ColorA  operator / ( float  const &n ) const { return ColorA(r/n, g/n, b/n, a/n); }

	//!@name Assignment operators
	ColorA& operator += ( ColorA const &c ) { r+=c.r; g+=c.g; b+=c.b; a+=c.a; return *this; }
	ColorA& operator -= ( ColorA const &c ) { r-=c.r; g-=c.g; b-=c.b; a-=c.a; return *this; }
	ColorA& operator *= ( ColorA const &c ) { r*=c.r; g*=c.g; b*=c.b; a*=c.a; return *this; }
	ColorA& operator /= ( ColorA const &c ) { r/=c.r; g/=c.g; b/=c.b; a/=c.a; return *this; }
	ColorA& operator += ( float  const &n ) { r+=n; g+=n; b+=n; a+=n; return *this; }
	ColorA& operator -= ( float  const &n ) { r-=n; g-=n; b-=n; a-=n; return *this; }
	ColorA& operator *= ( float  const &n ) { r*=n; g*=n; b*=n; a*=n; return *this; }
	ColorA& operator /= ( float  const &n ) { r/=n; g/=n; b/=n; a/=n; return *this; }

	//!@name Test operators
	bool operator == ( ColorA const &c ) const { return ( (c.r==r) && (c.g==g) && (c.b==b) && (c.a==a) ); }
	bool operator != ( ColorA const &c ) const { return ( (c.r!=r) || (c.g!=g) || (c.b!=b) || (c.a!=a) ); }

	//!@name Access operators
	float& operator [] ( int i )       { return (&r)[i]; }
	float  operator [] ( int i ) const { return (&r)[i]; }

	//!@name Static Methods
	static ColorA Black( float alpha=1.0f ) { return ColorA(0.0f,0.0f,0.0f,alpha); }	//!< Returns a black color
	static ColorA White( float alpha=1.0f ) { return ColorA(1.0f,1.0f,1.0f,alpha); }	//!< Returns a white color
};

//-------------------------------------------------------------------------------

//! 24-bit RGB color class with 3 unsigned byte components

class Color24
{
public:

	//!@name Color components
	uint8_t r, g, b;

	//!@name Constructors
	Color24() CY_CLASS_FUNCTION_DEFAULT
	Color24( Color24 const &c ) : r(c.r), g(c.g), b(c.b) {}
	explicit Color24( uint8_t _r, uint8_t _g, uint8_t _b ) : r(_r), g(_g), b(_b) {}
	explicit Color24( Color   const &c ) { r=FloatToByte(c.r); g=FloatToByte(c.g); b=FloatToByte(c.b); }
	explicit Color24( ColorA  const &c ) { r=FloatToByte(c.r); g=FloatToByte(c.g); b=FloatToByte(c.b); }
	explicit Color24( Color32 const &c );

	//!@name Conversion Methods
	Color  ToColor () const { return Color (r/255.0f,g/255.0f,b/255.0f); }
	ColorA ToColorA() const { return ColorA(r/255.0f,g/255.0f,b/255.0f,1.0f); }

	//!@name Set & Get value functions
	void SetBlack() { r=  0; g=  0; b=  0; }									//!< Sets r, g, and b components as zero
	void SetWhite() { r=255; g=255; b=255; }									//!< Sets r, g, and b components as 255
	void Set     ( uint8_t _r, uint8_t _g, uint8_t _b ) { r=_r; g=_g; b=_b; }	//!< Sets r, g, and b components as given
	void Set     ( uint8_t const *v ) { r=v[0]; g=v[1]; b=v[2]; }				//!< Sets r, g, and b components using the values in the given array
	void GetValue( uint8_t *v ) const { v[0]=r; v[1]=g; v[2]=b; }				//!< Puts r, g, and b values into the array

	//!@name Gray-scale functions
	int     Sum    () const { return int(r) + int(g) + int(b); }
	uint8_t Gray   () const { return uint8_t( (Sum()+1) / 3 ); }
	bool    IsBlack() const { return r==0 && g==0 && b==0; }			//!< Returns true if all components are exactly zero
	uint8_t Min    () const { return r<g ? (r<b ? r : b) : (g<b ? g : b); }
	uint8_t Max    () const { return r>g ? (r>b ? r : b) : (g>b ? g : b); }

	//!@name Limit functions
	void Clamp   ( uint8_t limitMin=  0, uint8_t limitMax=255 ) { ClampMin(limitMin); ClampMax(limitMax); }
	void ClampMin( uint8_t limitMin=  0 ) { r=cy::Max(r,limitMin); g=cy::Max(g,limitMin); b=cy::Max(b,limitMin); }
	void ClampMax( uint8_t limitMax=255 ) { r=cy::Min(r,limitMax); g=cy::Min(g,limitMax); b=cy::Min(b,limitMax); }

	//!@name Test operators
	bool operator == ( Color24 const &c ) const { return ( (c.r==r) && (c.g==g) && (c.b==b) ); }
	bool operator != ( Color24 const &c ) const { return ( (c.r!=r) || (c.g!=g) || (c.b!=b) ); }

	//!@name Access operators
	uint8_t& operator [] ( int i )       { return (&r)[i]; }
	uint8_t  operator [] ( int i ) const { return (&r)[i]; }

	//!@name Static Methods
	static Color24 Black() { return Color24(  0,  0,  0); }	//!< Returns a black color
	static Color24 White() { return Color24(255,255,255); }	//!< Returns a white color

protected:
	static uint8_t FloatToByte(float r) { return ClampInt(int(r*255+0.5f)); }
	static uint8_t ClampInt(int v) { return v<0 ? 0 : (v>255 ? 255 : v); }
};

//-------------------------------------------------------------------------------

//! 32-bit RGBA color class with 4 unsigned byte components

class Color32
{
public:

	//!@name Color components
	uint8_t r, g, b, a;

	//!@name Constructors
	Color32() CY_CLASS_FUNCTION_DEFAULT
	Color32( Color32 const &c ) : r(c.r), g(c.g), b(c.b), a(c.a) {}
	explicit Color32( uint8_t _r, uint8_t _g, uint8_t _b, uint8_t _a=255 ) : r(_r), g(_g), b(_b), a(_a) {}
	explicit Color32( Color   const &c, float _a=1.0f ) { r=FloatToByte(c.r); g=FloatToByte(c.g); b=FloatToByte(c.b); a=FloatToByte( _a); }
	explicit Color32( ColorA  const &c )                { r=FloatToByte(c.r); g=FloatToByte(c.g); b=FloatToByte(c.b); a=FloatToByte(c.a); }
	explicit Color32( Color24 const &c, uint8_t _a=255 ) : r(c.r), g(c.g), b(c.b), a(_a) {}

	//!@name Conversion Methods
	Color  ToColor () const { return Color (r/255.0f,g/255.0f,b/255.0f); }
	ColorA ToColorA() const { return ColorA(r/255.0f,g/255.0f,b/255.0f,a/255.0f); }

	//!@name Set & Get value functions
	void SetBlack( uint8_t _a=255 ) { r=  0; g=  0; b=  0; a=_a; }								//!< Sets r, g, and b components as zero and a component as given
	void SetWhite( uint8_t _a=255 ) { r=255; g=255; b=255; a=_a; }								//!< Sets r, g, and b components as 255 and a component as given
	void Set     ( uint8_t _r, uint8_t _g, uint8_t _b, uint8_t _a ) { r=_r; g=_g; b=_b; a=_a; }	//!< Sets r, g, b and a components as given
	void Set     ( uint8_t const *v ) { r=v[0]; g=v[1]; b=v[2]; a=v[3]; }						//!< Sets r, g, b and a components using the values in the given array
	void GetValue( uint8_t *v ) const { v[0]=r; v[1]=g; v[2]=b; v[3]=a; }						//!< Puts r, g, b and a values into the array

	//!@name Gray-scale functions
	int     Sum    () const { return int(r) + int(g) + int(b); }
	uint8_t Gray   () const { return uint8_t( (Sum()+1) / 3 ); }
	bool    IsBlack() const { return r==0 && g==0 && b==0; }			//!< Returns true if the r, g, and b components are exactly zero
	uint8_t Min    () const { uint8_t mrg = r<g ? r : g; uint8_t mba = b<a ? b : a; return mrg<mba ? mrg : mba; }
	uint8_t Max    () const { uint8_t mrg = r>g ? r : g; uint8_t mba = b>a ? b : a; return mrg>mba ? mrg : mba; }

	//!@name Limit functions
	void Clamp   ( uint8_t limitMin=  0, uint8_t limitMax=255 ) { ClampMin(limitMin); ClampMax(limitMax); }
	void ClampMin( uint8_t limitMin=  0 ) { r=cy::Max(r,limitMin); g=cy::Max(g,limitMin); b=cy::Max(b,limitMin); a=cy::Max(a,limitMin); }
	void ClampMax( uint8_t limitMax=255 ) { r=cy::Min(r,limitMax); g=cy::Min(g,limitMax); b=cy::Min(b,limitMax); a=cy::Min(a,limitMax); }

	//!@name Test operators
	bool operator == ( Color32 const &c ) const { return ( (c.r==r) && (c.g==g) && (c.b==b) && (c.a==a) ); }
	bool operator != ( Color32 const &c ) const { return ( (c.r!=r) || (c.g!=g) || (c.b!=b) || (c.a!=a) ); }

	//!@name Access operators
	uint8_t& operator [] ( int i )       { return (&r)[i]; }
	uint8_t  operator [] ( int i ) const { return (&r)[i]; }

	//!@name Static Methods
	static Color32 Black( uint8_t alpha=255 ) { return Color32(  0,  0,  0,alpha); }	//!< Returns a black color
	static Color32 White( uint8_t alpha=255 ) { return Color32(255,255,255,alpha); }	//!< Returns a white color

protected:
	static uint8_t FloatToByte(float r) { return ClampInt(int(r*255+0.5f)); }
	static uint8_t ClampInt(int v) { return v<0 ? 0 : (v>255 ? 255 : v); }
};

//-------------------------------------------------------------------------------

inline Color  ::Color  ( ColorA  const &c ) : r(c.r), g(c.g), b(c.b) {}
inline Color  ::Color  ( Color24 const &c ) { *this = c.ToColor(); }
inline Color  ::Color  ( Color32 const &c ) { *this = c.ToColor(); }
inline ColorA ::ColorA ( Color24 const &c, float alpha ) { Color rgb = c.ToColor(); r = rgb.r; g = rgb.g; b = rgb.b; a = alpha; }
inline ColorA ::ColorA ( Color32 const &c ) { *this = c.ToColorA(); }
inline Color24::Color24( Color32 const &c ) : r(c.r), g(c.g), b(c.b) {}

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Color   cyColor;	//!< RGB color class with 3 float components
typedef cy::ColorA  cyColorA;	//!< RGBA color class with 4 float components
typedef cy::Color24 cyColor24;	//!< 24-bit RGB color class with 3 unsigned byte components
typedef cy::Color32 cyColor32;	//!< 32-bit RGBA color class with 4 unsigned byte components

//-------------------------------------------------------------------------------

#endif

