// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyIPoint.h 
//! \author Cem Yuksel
//! 
//! \brief  2D, 3D, 4D, and ND integer point classes.
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

#ifndef _CY_IPOINT_H_INCLUDED_
#define _CY_IPOINT_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyIVector.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

#define IPoint  IVec
#define IPoint2 IVec2
#define IPoint3 IVec3
#define IPoint4 IVec4

//-------------------------------------------------------------------------------

typedef IPoint2<int8_t>   IPoint2b;		//!< 8-bit  signed   integer (int8_t)   2D integer point class
typedef IPoint3<int8_t>   IPoint3b;		//!< 8-bit  signed   integer (int8_t)   3D integer point class
typedef IPoint4<int8_t>   IPoint4b;		//!< 8-bit  signed   integer (int8_t)   4D integer point class
typedef IPoint2<uint8_t>  IPoint2ub;	//!< 8-bit  unsigned integer (uint8_t)  2D integer point class
typedef IPoint3<uint8_t>  IPoint3ub;	//!< 8-bit  unsigned integer (uint8_t)  3D integer point class
typedef IPoint4<uint8_t>  IPoint4ub;	//!< 8-bit  unsigned integer (uint8_t)  4D integer point class
typedef IPoint2<int16_t>  IPoint2s;		//!< 16-bit signed   integer (int16_t)  2D integer point class
typedef IPoint3<int16_t>  IPoint3s;		//!< 16-bit signed   integer (int16_t)  3D integer point class
typedef IPoint4<int16_t>  IPoint4s;		//!< 16-bit signed   integer (int16_t)  4D integer point class
typedef IPoint2<uint16_t> IPoint2us;	//!< 16-bit unsigned integer (uint16_t) 2D integer point class
typedef IPoint3<uint16_t> IPoint3us;	//!< 16-bit unsigned integer (uint16_t) 3D integer point class
typedef IPoint4<uint16_t> IPoint4us;	//!< 16-bit unsigned integer (uint16_t) 4D integer point class
typedef IPoint2<int32_t>  IPoint2i;		//!< 32-bit signed   integer (int32_t)  2D integer point class
typedef IPoint3<int32_t>  IPoint3i;		//!< 32-bit signed   integer (int32_t)  3D integer point class
typedef IPoint4<int32_t>  IPoint4i;		//!< 32-bit signed   integer (int32_t)  4D integer point class
typedef IPoint2<uint32_t> IPoint2ui;	//!< 32-bit unsigned integer (uint32_t) 2D integer point class
typedef IPoint3<uint32_t> IPoint3ui;	//!< 32-bit unsigned integer (uint32_t) 3D integer point class
typedef IPoint4<uint32_t> IPoint4ui;	//!< 32-bit unsigned integer (uint32_t) 4D integer point class
typedef IPoint2<int64_t>  IPoint2l;		//!< 64-bit signed   integer (int64_t)  2D integer point class
typedef IPoint3<int64_t>  IPoint3l;		//!< 64-bit signed   integer (int64_t)  3D integer point class
typedef IPoint4<int64_t>  IPoint4l;		//!< 64-bit signed   integer (int64_t)  4D integer point class
typedef IPoint2<uint64_t> IPoint2ul;	//!< 64-bit unsigned integer (uint64_t) 2D integer point class
typedef IPoint3<uint64_t> IPoint3ul;	//!< 64-bit unsigned integer (uint64_t) 3D integer point class
typedef IPoint4<uint64_t> IPoint4ul;	//!< 64-bit unsigned integer (uint64_t) 4D integer point class

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::IPoint2b  cyIPoint2b;		//!< 8-bit  signed   integer (int8_t)   2D integer point class
typedef cy::IPoint3b  cyIPoint3b;		//!< 8-bit  signed   integer (int8_t)   3D integer point class
typedef cy::IPoint4b  cyIPoint4b;		//!< 8-bit  signed   integer (int8_t)   4D integer point class
typedef cy::IPoint2ub cyIPoint2ub;		//!< 8-bit  unsigned integer (uint8_t)  2D integer point class
typedef cy::IPoint3ub cyIPoint3ub;		//!< 8-bit  unsigned integer (uint8_t)  3D integer point class
typedef cy::IPoint4ub cyIPoint4ub;		//!< 8-bit  unsigned integer (uint8_t)  4D integer point class
typedef cy::IPoint2s  cyIPoint2s;		//!< 16-bit signed   integer (int16_t)  2D integer point class
typedef cy::IPoint3s  cyIPoint3s;		//!< 16-bit signed   integer (int16_t)  3D integer point class
typedef cy::IPoint4s  cyIPoint4s;		//!< 16-bit signed   integer (int16_t)  4D integer point class
typedef cy::IPoint2us cyIPoint2us;		//!< 16-bit unsigned integer (uint16_t) 2D integer point class
typedef cy::IPoint3us cyIPoint3us;		//!< 16-bit unsigned integer (uint16_t) 3D integer point class
typedef cy::IPoint4us cyIPoint4us;		//!< 16-bit unsigned integer (uint16_t) 4D integer point class
typedef cy::IPoint2i  cyIPoint2i;		//!< 32-bit signed   integer (int32_t)  2D integer point class
typedef cy::IPoint3i  cyIPoint3i;		//!< 32-bit signed   integer (int32_t)  3D integer point class
typedef cy::IPoint4i  cyIPoint4i;		//!< 32-bit signed   integer (int32_t)  4D integer point class
typedef cy::IPoint2ui cyIPoint2ui;		//!< 32-bit unsigned integer (uint32_t) 2D integer point class
typedef cy::IPoint3ui cyIPoint3ui;		//!< 32-bit unsigned integer (uint32_t) 3D integer point class
typedef cy::IPoint4ui cyIPoint4ui;		//!< 32-bit unsigned integer (uint32_t) 4D integer point class
typedef cy::IPoint2l  cyIPoint2l;		//!< 64-bit signed   integer (int64_t)  2D integer point class
typedef cy::IPoint3l  cyIPoint3l;		//!< 64-bit signed   integer (int64_t)  3D integer point class
typedef cy::IPoint4l  cyIPoint4l;		//!< 64-bit signed   integer (int64_t)  4D integer point class
typedef cy::IPoint2ul cyIPoint2ul;		//!< 64-bit unsigned integer (uint64_t) 2D integer point class
typedef cy::IPoint3ul cyIPoint3ul;		//!< 64-bit unsigned integer (uint64_t) 3D integer point class
typedef cy::IPoint4ul cyIPoint4ul;		//!< 64-bit unsigned integer (uint64_t) 4D integer point class

//-------------------------------------------------------------------------------

#endif

