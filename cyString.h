// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyString.h 
//! \author Cem Yuksel
//! 
//! \brief  String class
//! 
//! This file includes a general purpose string class for char arrays.
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

#ifndef _CY_STRING_H_INCLUDED_
#define _CY_STRING_H_INCLUDED_

//-------------------------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

//! A general-purpose string class for char arrays.
//!
//! This is a general purpose string class, which supports many useful string operations.
//! It has methods for numeric conversions, character and string search,
//! formating, editing, and file operations.

class String
{

#ifdef _OSTREAM_
	//! Overloaded stream operator. It only works if you include iostream before this file.
	friend std::ostream& operator << ( std::ostream &os, const String &str )
	{
		os<<str.GetString();
		return os;
	}
#endif

public:

	//////////////////////////////////////////////////////////////////////////
	//!@name Constructors and Destructor

	//! Default constructor
	String() : string(nullptr), length(0) { EmptyString(); }

	//! Copy constructor
	String( const String &src ) : string(nullptr), length(0) { Set(src); }

	//! Sets the string using the given array.
	String( const char *src ) : string(nullptr), length(0) { Set(src); }

	//! Sets the string using the given number.
	String( char   src ) : string(nullptr), length(0) { Set(src); }
	String( int    src ) : string(nullptr), length(0) { Set(src); }
	String( long   src ) : string(nullptr), length(0) { Set(src); }
	String( float  src ) : string(nullptr), length(0) { Set(src); }
	String( double src ) : string(nullptr), length(0) { Set(src); }

	//! Sets the string using the first 'count' characters of the given array.
	String( const char *src, unsigned int count ) : string(nullptr), length(0) { Set(src,count); }

	//! Sets the string as the given double number 'src' using given number of digits and precision.
	String( double src, int digits, int precisition ) : string(nullptr), length(0) { Set(src,digits,precisition); }

	//! Sets the string using the given format and predicted string size.
	//! The predicted size should be greater than or equal to the final size.
	String( int size, const char *format, ... ) : string(nullptr), length(0)
	{
		va_list args;
		va_start(args,format);
		FormatString(size,format,args);
	}

	//! Destructor
	virtual ~String()
	{
		DeleteString( string );
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Set methods (Assign value)

	//! Sets the string using the given string. Returns the string.
	String& Set( const String &str )
	{
		SetStrLen( str.length );
		strncpy( string, str.string, length );
		return *this;
	}

	//! Sets the string using the given array. Returns the string.
	String& Set( const char *src )
	{
		int n = (int) strlen(src);
		SetStrLen( n );
		strncpy( string, src, n );
		return *this;
	}

	//! Sets the string as the given character. Returns the string.
	String& Set( char src )
	{
		SetStrLen(1);
		string[0] = src;
		return *this;
	}

	//! Sets the string using the given number. Returns the string.
	String& Set( int src )
	{
		SetStrLen(30);
		sprintf( string, "%d", src );
		Shrink();
		return *this;
	}

	//! Sets the string using the given number. Returns the string.
	String& Set( long src )
	{
		SetStrLen(30);
		sprintf( string, "%d", src );
		Shrink();
		return *this;
	}

	//! Sets the string using the given number. Returns the string.
	String& Set( float src )
	{
		SetStrLen(30);
		sprintf( string, "%f", src );
		Shrink();
		return *this;
	}

	//! Sets the string using the given number. Returns the string.
	String& Set( double src )
	{
		SetStrLen(30);
		sprintf( string, "%f", src );
		Shrink();
		return *this;
	}

	//! Sets the string using the first 'count' characters of the given array. Return this string.
	String& Set( const char *src, unsigned int count )
	{
		SetStrLen(count);
		strncpy( string, src, length );
		return *this;
	}

	//! Sets the string as the given double number 'src'
	//! using given number of digits and precision. Returns this string.
	String& Set( double src, int digits, int precision )
	{
		SetStrLen(30);
		char format[30];
		sprintf( format, "%%%d.%df", digits, precision );
		sprintf( string, format, src );
		Shrink();
		return *this;
	}

	//! Sets the string using the given format and predicted string size.
	//! The predicted size should be greater than or equal to the final size.
	//! Returns this string.
	String& Set( int size, const char *format, ... )
	{
		va_list args;
		va_start( args, format );
		FormatString( size, format, args );
		return *this; 
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Assignment operators

	String& operator = ( const String &src ){ return Set(src); }
	String& operator = ( const char *src )	{ return Set(src); }
	String& operator = ( char src )			{ return Set(src); }
	String& operator = ( int src )			{ return Set(src); }
	String& operator = ( long src )			{ return Set(src); }
	String& operator = ( float src )			{ return Set(src); }
	String& operator = ( double src )			{ return Set(src); }


	//////////////////////////////////////////////////////////////////////////
	//!@name Editing methods

	//! Converts to empty string.
	void EmptyString() { SetStrLen(0); }

	//! Changes the length of the array by keeping the values in the string.
	//! If the new size is smaller, a part of the string is lost.
	//! If the new size is greater, space characters are added at the end of the string.
	int SetLength( int newLength )
	{
		if ( newLength < 0 || newLength == length ) return length;
		if ( newLength > 0 ) {
			char *str = NewString(newLength);
			int l = (length < newLength) ? length : newLength;
			strncpy(str,string,l);
			memset(str+l,' ',newLength-l);
			ReplaceString( str, newLength );
		} else {
			EmptyString();
		}
		return length;
	}

	//! Creates a string of length 'count' with the given character 'chr'.
	void SetCharString( char chr, int count )
	{
		SetStrLen( count );
		memset( string, chr, length );
	}

	//! Deletes 'count' number of characters from the array, starting from 'start'.
	//! If 'count' is too big, deletes till the end of the string.
	void Delete( int start, int count )
	{
		if ( start > length || count <= 0 ) return;
		int end = start + count;
		if ( end > length ) { count -= end - length; end = length; }
		length -= count;
		char *string2 = NewString(length);
		strncpy ( string2, string, start );
		strncpy ( string2+start, string+end, length - start );
		ReplaceString( string2 );
	}

	//! Inserts the given string to the given index.
	void Insert( const String &str, int index )
	{
		int count = length - index;
		if ( count < 0 ) { index = length; count = 0; }
		length += str.length;
		char *string2 = NewString(length);
		strncpy( string2, string, index );
		strncpy( string2+index, str.string, str.length );
		strncpy( string2+index+str.length, string+index, count );
		ReplaceString( string2 );
	}

	//! Deletes all carriage return characters.
	//! Returns the number of characters deleted.
	int DeleteCarriageReturns() { return DeleteChar('\r'); }

	//! Insert a carriage return character after every new line character
	//! if it doesn't already exits. Returns the number of carriage return characters inserted.
	int InsertCarriageReturns()
	{
		int nn = CountChar('\n');
		if ( nn == 0 ) return 0;
		int nr = CountChar('\r');
		char *string2 = NewString( length + nn - nr );
		int i = length - 1;
		int j = length + nn - nr - 1;
		while ( j >= 0 ) {
			bool needCR = false;
			if ( string[i] == '\n' ) needCR = true;
			string[j--] = string[i--];
			if ( needCR ) {
				if ( i > 0 && string[i] == '\r' ) continue;
				string[j--] = '\r';
			}
		}
		ReplaceString( string2, length + nn - nr );
		Shrink();
		return nn - nr;

	}

	//! Deletes all the characters in the string that are equal to the given character 'c'.
	//! Returns the number of characters deleted.
	int DeleteChar( char c )
	{
		int n = 0;
		char *str = strchr( string, c );
		while ( str != nullptr ) {
			n++;
			int count;
			char *str2 = strchr( ++str, c );
			if ( str2 != nullptr ) count = int(str2 - str);
			else count = int(string - str) + length + 1;
			memmove( str-n, str, count );
			str = str2;
		}
		Shrink();
		return n;
	}

	//! Deletes all the characters in the string that exits in the given set.
	//! Returns the number of characters deleted.
	int DeleteChars( const String &set )
	{
		int n = 0;
		char *str = strpbrk( string, set.string );
		while ( str != nullptr ) {
			n++;
			int count;
			char *str2 = strpbrk( ++str, set.string );
			if ( str2 != nullptr ) count = int(str2 - str);
			else count = int(string - str) + length + 1;
			memmove( str-n, str, count );
			str = str2;
		}
		Shrink();
		return n;
	}

	//! Reverses the string.
	void Reverse()
	{
		for ( int i=0; i<length/2; i++ ) {
			SwapChars( i, length-i-1 );
		}
	}

	//! Swap 'count' number of characters starting from 'index1' and 'index2'.
	void SwapChars( int index1, int index2, int count=1 )
	{
		_swab( string + index1, string + index2, count );
	}

	//! Sets the string using the given format and predicted string size.
	//! The predicted size should be greater than or equal to the final size.
	//! Returns this string.
	String& Format( int size, const char *format, ... )
	{
		va_list args;
		va_start( args, format );
		FormatString( size, format, args );
		return *this; 
	}

	//! Returns the character of the string at the specified index.
	char& operator [] ( const int index ) { return string[index]; }


	//////////////////////////////////////////////////////////////////////////
	//!@name Get value methods

	//! Returns the length of the string.
	int Length() const { return length; }

	//! Returns the char array that keep the string data.
	const char* GetString() const { return string; }

	//! Checks if the string is empty (length is zero).
	int IsEmpty() const { return (length==0); }

	//! Returns a pointer to the last character before the first null character.
	//! If null character is not found, or the string is empty returns nullptr pointer.
	char* LastChar() const
	{
		char * lc = strchr(string,'\0') - 1;
		if ( lc < string ) return nullptr;
		else return lc;
	}

	//! If the character at the given index exists in the given set,
	//! returns the index of the character.
	//! If the character is not found in the set, returns -1.
	int SearchCharInSet( int index, const String &set ) const
	{
		return set.GetPosition( string[index] );
	}

	//! Searches the string starting from the beginning for a character that exists in the given set.
	//! Returns the first position of a character from the set.
	//! If not found, returns -1.
	int SearchSet( const String &set ) const
	{
		int pos = (int) strcspn( string, set.string );
		if ( pos < length ) return pos;
		else return -1;
	}

	//! Searches the string starting from the end for a character that exists in the given set.
	//! Returns the last position of a character from the set.
	//! If not found, returns -1.
	int ReverseSearchSet( const String &set ) const
	{
		for ( int i=length-1; i >= 0; i-- ) {
			int pos = set.GetPosition( string[i] );
			if ( pos >= 0 ) return i;
		}
		return -1;
	}

	//! Searches the string starting from the beginning for a character that exists in the given set.
	//! Returns the first position of a character from the set.
	//! If not found, returns -1.
	int SearchNextSet( int start, const String &set ) const
	{
		if ( start > length ) return -1;
		int pos = (int) strcspn( string + start, set.string );
		if ( pos < length ) return pos + start;
		else return -1;
	}

	//! Returns the first position of the character.
	//! If the character is not found, returns -1.
	int GetPosition( char c ) const
	{
		char *cp = strchr( string, c );
		if ( cp != nullptr ) return int(cp - string);
		else return -1;
	}

	//! Returns the first position of the given string.
	//! If the string is not found, returns -1.
	int GetPosition( const String &str ) const
	{
		char *sub = strstr( string, str.GetString() );
		if ( sub != nullptr ) return int(sub - string);
		else return -1;
	}

	//! Returns the last position of the character.
	//! If the character is not found, returns -1.
	int GetLastPosition( char c ) const
	{
		char *cp = strrchr( string, c );
		if ( cp != nullptr ) return int(cp - string);
		else return -1;
	}

	//! Returns the last position of the given string.
	//! If the string is not found, returns -1.
	int GetLastPosition( const String &str ) const
	{
		String s(*this);
		String str2(str);
		s.Reverse();
		str2.Reverse();
		int pos = s.GetPosition(str2);
		if ( pos >= 0 )	return length - pos - str.length;
		else return -1;
	}

	//! Returns the next position of the character.
	//! If the character is not found, returns -1.
	int GetNextPosition( int start, char c ) const
	{
		if ( start > length ) return -1;
		char *cp = strchr( string+start, c );
		if ( cp != nullptr ) return int(cp - string);
		else return -1;
	}

	//! Returns the next position of the given string.
	//! If the string is not found, returns -1.
	int GetNextPosition( int start, const String &str ) const
	{
		if ( start > length ) return -1;
		char *sub = strstr( string+start, str.GetString() );
		if ( sub != nullptr ) return int(sub - string);
		else return -1;
	}

	//! Returns the line number of the given position in the string
	int GetLineNumber( int pos ) const { return CountChar('\n',0,pos) + 1; }

	//! Returns the number of lines in the string.
	int GetNumLines() const { return CountChar('\n') + 1; }

	//! Returns the number of characters in the string that are equal to the given character 'c'.
	int CountChar( char c ) const
	{
		int n = 0;
		char *str = strchr( string, c );
		while ( str != nullptr ) {
			n++;
			str = strchr( str+1, c );
		}
		return n;
	}

	//! Returns the number of characters in the string that are equal to the given character 'c'.between start and end
	int CountChar( char c, int start, int end ) const
	{
		int n = 0;

		if ( start > length ) return 0;
		if ( end > length ) end = length;
		if ( end < start ) return 0;

		char *str = string + start;

		str = strchr( str, c );
		while ( str != nullptr &&  str <= (string + end) ) {
			n++;
			str = strchr( str+1, c );
		}
		return n;
	}

	//! Returns the number of characters in the string that exists in the given character set.
	int CountChars( const String &set ) const
	{
		return (int) strspn( string, set.string );
	}

	//! Returns the character of the string at the specified index.
	char operator [] ( const int index ) const { return string[index]; }


	//////////////////////////////////////////////////////////////////////////
	//!@name Comparison methods

	//! Checks if this string is the same as 'str'.
	int operator == ( const String &str ) const { return ( Compare(str) == 0 ); }

	//! Checks if this string is not the same as 'str'.
	int operator != ( const String &str ) const { return ( Compare(str) != 0 ); }

	//! Checks if this string is smaller than 'str'.
	int operator <  ( const String &str ) const { return ( Compare(str) <  0 ); }

	//! Checks if this string is smaller than or equal to 'str'.
	int operator <= ( const String &str ) const { return ( Compare(str) <= 0 ); }

	//! Checks if this string is greater than 'str'.
	int operator >  ( const String &str ) const { return ( Compare(str) >  0 ); }

	//! Checks if this string is greater than or equal to 'str'.
	int operator >= ( const String &str ) const { return ( Compare(str) >= 0 ); }

	//! Compares this string to 'str' (case sensitive).
	//! Returns zero if the strings are the same,
	//! smaller than zero if this string is smaller than 'str',
	//! greater than zero if this string is greater than 'str'.
	int Compare( const String &str ) const
	{
		if ( length == str.length ) {
			if ( length > 0 ) return strcmp( string, str.string );
			else if ( str.length == 0 ) return 0;
			else return -1;
		} else return length - str.length;
	}

	//! Compares this string to 'str' (case sensitive) limited to the first 'count' characters.
	//! Returns zero if the strings are the same,
	//! smaller than zero if this string is smaller than 'str',
	//! greater than zero if this string is greater than 'str'.
	int Compare( const String &str, int count ) const
	{
		return strncmp( string, str.string, count );
	}

	//! Compares this string to 'str' (case insensitive).
	//! Returns zero if the strings are the same,
	//! smaller than zero if this string is smaller than 'str',
	//! greater than zero if this string is greater than 'str'.
	int CompareIC( const String &str ) const
	{
		String str1(*this);
		String str2(str);
		str1.LowerCase();
		str2.LowerCase();
		return str1.Compare( str2 );
	}

	//! Compares this string to 'str' (case insensitive) limited to the first 'count' characters.
	//! Returns zero if the strings are the same,
	//! smaller than zero if this string is smaller than 'str',
	//! greater than zero if this string is greater than 'str'.
	int CompareIC( const String &str, int count ) const
	{
		String str1(*this);
		String str2(str);
		str1.LowerCase();
		str2.LowerCase();
		return str1.Compare( str2, count );
	}

	//! Returns if the char in the given position is a capital letter
	int IsCapitalChar( int pos ) const
	{
		if ( string[pos] < 65 ) return false;
		if ( string[pos] > 90 ) return false;
		return true;
	}

	//! Returns if the char in the given position is a lower case letter
	int IsLowerCaseChar( int pos ) const
	{
		if ( string[pos] < 97 ) return false;
		if ( string[pos] > 122 ) return false;
		return true;
	}

	//! Returns if the char in the given position is a numeric character
	int IsNumericChar( int pos ) const
	{
		if ( string[pos] < 48 ) return false;
		if ( string[pos] > 57 ) return false;
		return true;
	}

	//! Returns if the char in the given position is a lower case or capital letter
	int IsAlphaChar( int pos ) const
	{
		if ( IsLowerCaseChar(pos) ) return true;
		if ( IsCapitalChar(pos) ) return true;
		return false;
	}

	//! Returns if the char in the given position is a lower case or capital letter or a numeric character
	int IsAlphaNumericChar( int pos ) const
	{
		if ( IsNumericChar(pos) ) return true;
		if ( IsAlphaChar(pos) ) return true;
		return false;
	}

	//! Returns if the char in the given position is a lower case or capital letter, a numeric character or '_'
	int IsNameChar( int pos ) const
	{
		if ( string[pos] == '-') return true;
		else return IsAlphaNumericChar( pos );
	}

	//! Returns if the char in the given position is a space character
	int IsSpaceChar( int pos ) const
	{
		switch ( string[pos] ) {
		case ' ':
		case '\t':
		case '\n':
		case '\r':
			return true;
		default:
			return false;
		}
	}

	//! Returns if the char in the given position is a control char (including new line etc.)
	int IsControlChar( int pos ) const
	{
		if ( string[pos] < 32 ) return true;
		if ( string[pos] > 126 ) return true;
		return false;
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name Append and substring

	//! Appends to the end of the string.
	void Append( const String &str )
	{
		int length2 = length + str.length;
		char *string2 = NewString(length2);
		strncpy(string2,string,length);
		strncpy(string2+length,str.string,str.length);
		ReplaceString( string2, length2 );
	}

	//! Appends to the end of the string and returns this string.
	String& operator += ( const String &str ) { Append(str); return *this; }

	//! Appends two strings and returns the result
	String operator + ( const String &right ) {
		String r;
		r.SetStrLen( length + right.length );
		strncpy( r.string, string, length );
		strncpy( r.string+length, right.string, right.length );
		return r;
	}

	//! Converts the string to its sub-string starting from 'start' with length of 'count'.
	//! If 'count' is too large, returns the sub-string till the end of the string.
	void SubString( int start, int count )
	{
		if ( start > length ) start = length;
		int end = start + count;
		if ( end > length ) { count -= end - length; end = length; }
		length = count;
		char *string2 = NewString(length);
		strncpy( string2, string+start, length );
		ReplaceString( string2 );
	}

	//! Returns the sub-string from 'start' with the length of 'count'.
	//! If 'count' is too big, returns the sub-string till the end of the string.
	String GetSubString( int start, int count ) const
	{
		if ( start > length ) start = length;
		int end = start + count;
		if ( end > length ) { count -= end - length; end = length; }
		String str(string+start,count);
		return str;
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name Conversion to numbers and number arrays

	//! Converts the string to int.
	int ToInt() const { return atoi(string); }

	//! Converts the string to long.
	long ToLong() const { return atol(string); }

	//! Converts the string to float.
	float ToFloat() const { return (float) atof(string); }

	//! Converts the string to double.
	double ToDouble() const { return atof(string); }

	//! Converts the string (ex: "23,43,14") to int array.
	void ToIntArray(int n, int *v) const
	{
		int p1 = 0;
		for ( int i=0; i<n; i++ ) {
			int p2 = GetNextPosition( p1, ',' );
			String sub;
			if ( p2 > 0 ) {
				sub = GetSubString(p1,p2-p1);
				p1 = p2+1;
			} else {
				sub = GetSubString(p1,length);
			}
			v[i] = sub.ToInt();
		}
	}

	//! Converts the string to long.
	void ToLongArray(int n, long *v) const
	{
		int p1 = 0;
		for ( int i=0; i<n; i++ ) {
			int p2 = GetNextPosition( p1, ',' );
			String sub;
			if ( p2 > 0 ) {
				sub = GetSubString(p1,p2-p1);
				p1 = p2+1;
			} else {
				sub = GetSubString(p1,length);
			}
			v[i] = sub.ToLong();
		}
	}

	//! Converts the string (ex: "23,43,14") to float array.
	void ToFloatArray(int n, float *v) const
	{
		int p1 = 0;
		for ( int i=0; i<n; i++ ) {
			int p2 = GetNextPosition( p1, ',' );
			String sub;
			if ( p2 > 0 ) {
				sub = GetSubString(p1,p2-p1);
				p1 = p2+1;
			} else {
				sub = GetSubString(p1,length);
			}
			v[i] = sub.ToFloat();
		}
	}

	//! Converts the string (ex: "23,43,14") to double array.
	void ToDoubleArray(int n, double *v) const
	{
		int p1 = 0;
		for ( int i=0; i<n; i++ ) {
			int p2 = GetNextPosition( p1, ',' );
			String sub;
			if ( p2 > 0 ) {
				sub = GetSubString(p1,p2-p1);
				p1 = p2+1;
			} else {
				sub = GetSubString(p1,length);
			}
			v[i] = sub.ToDouble();
		}
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name Trim and shrink
	
	//! Deletes empty characters from the beginning and end of the string.
	void Trim() { TrimRight(); TrimLeft(); }

	//! Deletes empty characters from the beginning of the string.
	void TrimLeft()
	{
		int i;
		unsigned char *uchr = (unsigned char *) string; // unsigned to support extended ASCII characters
		for ( i=0; i<length; i++ ) {
			if ( *uchr > 32 && *uchr != 127 && *uchr != 255 ) break;
			uchr++;
		}
		Delete(0,i);
	}

	//! Deletes empty characters from the end of the string.
	void TrimRight()
	{
		unsigned char *uchr = (unsigned char *) LastChar(); // unsigned to support extended ASCII characters
		if ( uchr ) {
			unsigned char *str = (unsigned char *) string;
			for ( ; uchr >= str; uchr--) {
				if ( *uchr > 32 && *uchr != 127 && *uchr != 255 ) break;
			}
			Delete( int(uchr-str)+1, length );
		}
	}

	//! Deletes the rest of the string after the first null character.
	void Shrink()
	{
		char *end = strchr(string,'\0');
		if ( end != nullptr ) {
			int count = int(end - string);
			if ( count < length ) SubString(0,count);
		}
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name Case change

	//! Coverts the string to lower case.
	void LowerCase()
	{
		for ( int i=0; i<length; i++ )
			if ( string[i] >= 65 && string[i] <= 90 ) string[i] += 32;
	}

	//! Coverts the string to upper case.
	void UpperCase()
	{
		for ( int i=0; i<length; i++ )
			if ( string[i] >= 97 && string[i] <= 122 ) string[i] -= 32;
	}

	//! Coverts the string to swapped case.
	void SwapCase()
	{
		for ( int i=0; i<length; i++ ) {
			if ( string[i] >= 65 && string[i] <= 90 ) string[i] += 32;
			else if ( string[i] >= 97 && string[i] <= 122 ) string[i] -= 32;
		}
	}

	//! Coverts the string to title case.
	void TitleCase()
	{
		int wordstart = true;
		for ( int i=0; i<length; i++ ) {
			if ( wordstart ) {
				if ( string[i] >= 97 && string[i] <= 122 ) {
					string[i] -= 32;
					wordstart = false;
				}
				else if( string[i] >= 65 && string[i] <= 90 ) wordstart = false;
			} else {
				if ( string[i] >= 65 && string[i] <= 90 ) string[i] += 32;
				else if( string[i] < 97 || string[i] > 122 ) wordstart = true;
			}
		}
	}
	
	//! Returns the lower case version of the string.
	String GetLowerCase() const { String str(*this); str.LowerCase(); return str; }

	//! Returns the upper case version of the string.
	String GetUpperCase() const { String str(*this); str.UpperCase(); return str; }

	//! Returns the swapped case version of the string.
	String GetSwapCase()  const { String str(*this); str.SwapCase();  return str; }

	//! Returns the title case version of the string.
	String GetTitleCase() const { String str(*this); str.TitleCase(); return str; }


	//////////////////////////////////////////////////////////////////////////
	//!@name Load/save and file name methods

	//! Parses the string as file name string and returns path, filename.
	void DecodeFileName( String &path, String &filename )
	{
		int pos = ReverseSearchSet("\\/:");
		if ( pos >= 0 ) {
			path = GetSubString(0,pos+1);
			filename = GetSubString(pos+1,length);
		} else {
			path.EmptyString();
			filename.Set(*this);
		}
	}

	//! Parses the string as file name string and returns path, filename and extension.
	void DecodeFileName( String &path, String &filename, String &extention )
	{
		String fname;
		DecodeFileName( path, fname );
		int pos = fname.GetLastPosition( '.' );
		if ( pos > 0 && pos < length - 1 ) {
			filename = fname.GetSubString( 0, pos );
			extention = fname.GetSubString( pos+1, fname.length );
		} else {
			filename.Set(fname);
			extention.EmptyString();
		}
	}

	//! Reads 'count' number of characters from the given file stream onto the string.
	//! Returns final length of the string.
	int LoadFromStream( FILE *stream, unsigned int count )
	{
		SetStrLen(count);
		if ( count ) {
			count = (unsigned int) fread( string, 1, count, stream );
			string[count] = '\0';
			Shrink();
		}
		return length;
	}

	//! Loads the string from the file. Returns final length of the string.
	int LoadFromFile( const String &filename )
	{
		FILE *fp = fopen( filename.GetString(), "r" );
		if ( fp != nullptr ) {
			int seek = fseek( fp, 0, SEEK_END );
			if ( seek == 0 ) {
				int count = ftell(fp);
				fseek( fp, 0, SEEK_SET );
				LoadFromStream( fp, count );
				fclose ( fp );
				return length;
			}
		}
		EmptyString();
		return length;
	}

	//! Writes string data to the given file stream. Returns the number of bytes written.
	int SaveToStream( FILE *stream ) const
	{
		return (int) fwrite( string, 1, length, stream );
	}

	//! Saves the string to the text file. Returns the number of bytes written.
	int SaveToFile( const String &filename ) const
	{
		FILE *fp = fopen( filename.GetString(), "w" );
		if ( fp != nullptr ) {
			int count = SaveToStream( fp );
			fclose ( fp );
			return count;
		}
		return 0;
	}

private:

	//! \internal
	//////////////////////////////////////////////////////////////////////////
	//!@name Internal methods

	char *string;	// keeps the string data
	int  length;	// keeps the length of the string (the actual length of the array is length +1)

	// Deletes string and allocates memory for new string with the given length.
	// Returns the length of the string.
	int SetStrLen( int newLength )
	{
		if ( newLength < 0 || newLength == length ) return length;
		ReplaceString( NewString(newLength), newLength );
		return length;
	}

	// Deletes string data and sets as the given string.
	void ReplaceString( char *str )
	{
		DeleteString( string );
		string = str;
	}

	// Deletes string data and sets as the given string and changes the length.
	void ReplaceString( char *str, int len )
	{
		DeleteString( string );
		string = str;
		length = len;
	}

	// Allocates memory for a string of given size and returns the new string.
	// Sets the last character at [size] as null.
	char* NewString( int size )
	{
		char *str = new char[size+1];
		str[size] = '\0';
		return str;
	}

	// Deletes allocated memory for a string.
	void DeleteString( char *str )
	{
		if ( str ) delete [] str;
	}

	// Given a predicted size (should be greater than or equal to the final size)
	// format and arguments list (args), formats the string.
	void FormatString( int size, const char *format, va_list args )
	{
		SetStrLen(size);
		vsprintf( string, format, args );
		Shrink();
	}

};

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::String cyString;	//!< String class for char arrays.

//-------------------------------------------------------------------------------

#endif
