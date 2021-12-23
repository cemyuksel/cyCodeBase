// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyPoint.h 
//! \author Cem Yuksel
//! 
//! \brief  2D, 3D, 4D, and ND point (vector) classes. This file is deprecated.
//! 
//! This file provides backwards compatibility for vector classes.
//! It is deprecated. Use cyVector.h and classes in that header file instead.
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

#ifndef _CY_POINT_H_INCLUDED_
#define _CY_POINT_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyVector.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

#define Point  Vec	//!< A general class for N-dimensional points (vectors).
#define Point2 Vec2	//!< 2D point (vector) class
#define Point3 Vec3	//!< 3D point (vector) class
#define Point4 Vec4	//!< 4D point (vector) class

//-------------------------------------------------------------------------------

typedef Point2<float>  Point2f;	//!< 2D point (vector) class with float type elements
typedef Point3<float>  Point3f;	//!< 3D point (vector) class with float type elements
typedef Point4<float>  Point4f;	//!< 4D point (vector) class with float type elements

typedef Point2<double> Point2d;	//!< 2D point (vector) class with double type elements
typedef Point3<double> Point3d;	//!< 3D point (vector) class with double type elements
typedef Point4<double> Point4d;	//!< 4D point (vector) class with double type elements

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Point2f cyPoint2f;	//!< 2D point (vector) class with float type elements
typedef cy::Point3f cyPoint3f;	//!< 3D point (vector) class with float type elements
typedef cy::Point4f cyPoint4f;	//!< 4D point (vector) class with float type elements

typedef cy::Point2d cyPoint2d;	//!< 2D point (vector) class with double type elements
typedef cy::Point3d cyPoint3d;	//!< 3D point (vector) class with double type elements
typedef cy::Point4d cyPoint4d;	//!< 4D point (vector) class with double type elements

//-------------------------------------------------------------------------------

#endif

