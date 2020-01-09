// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyPointCloud.h 
//! \author Cem Yuksel
//!
//! \brief  Point cloud using a k-d tree
//! 
//! This file includes a class that keeps a point cloud as a k-d tree
//! for quickly finding n-nearest points to a given location.
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

#ifndef _CY_POINT_CLOUD_H_INCLUDED_
#define _CY_POINT_CLOUD_H_INCLUDED_

//-------------------------------------------------------------------------------

#ifdef max
# define _CY_POP_MACRO_max
# pragma push_macro("max")
# undef max
#endif

#ifdef min
# define _CY_POP_MACRO_min
# pragma push_macro("min")
# undef min
#endif

//-------------------------------------------------------------------------------

#ifndef _CY_PARALLEL_LIB
# ifdef __TBB_tbb_H
#  define _CY_PARALLEL_LIB tbb
# elif defined(_PPL_H)
#  define _CY_PARALLEL_LIB concurrency
# endif
#endif

//-------------------------------------------------------------------------------

#include <assert.h>
#include <algorithm>
#include <stdint.h>

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

//! A point cloud class that uses a k-d tree for storing points.
//!
//! The GetPoints and GetClosest methods return the neighboring points to a given location.

template <typename PointType, typename FType, uint32_t DIMENSIONS, typename SIZE_TYPE=uint32_t>
class PointCloud
{
public:
	/////////////////////////////////////////////////////////////////////////////////
	//!@name Constructors and Destructor

	PointCloud() : points(nullptr), pointCount(0) {}
	PointCloud( SIZE_TYPE numPts, PointType const *pts, SIZE_TYPE const *customIndices=nullptr ) : points(nullptr), pointCount(0) { Build(numPts,pts,customIndices); }
	~PointCloud() { delete [] points; }

	/////////////////////////////////////////////////////////////////////////////////
	//!@ Access to internal data

	SIZE_TYPE GetPointCount() const { return pointCount-1; }					//!< Returns the point count
	PointType const & GetPoint(SIZE_TYPE i) const { return points[i+1].Pos(); }	//!< Returns the point at position i
	SIZE_TYPE GetPointIndex(SIZE_TYPE i) const { return points[i+1].Index(); }	//!< Returns the index of the point at position i

	/////////////////////////////////////////////////////////////////////////////////
	//!@ Initialization

	//! Builds a k-d tree for the given points.
	//! The positions are stored internally.
	//! The build is parallelized using Intel's Thread Building Library (TBB) or Microsoft's Parallel Patterns Library (PPL),
	//! if ttb.h or ppl.h is included prior to including cyPointCloud.h.
	void Build( SIZE_TYPE numPts, PointType const *pts ) { BuildWithFunc( numPts, [&pts](SIZE_TYPE i){ return pts[i]; } ); }

	//! Builds a k-d tree for the given points.
	//! The positions are stored internally, along with the indices to the given array.
	//! The build is parallelized using Intel's Thread Building Library (TBB) or Microsoft's Parallel Patterns Library (PPL),
	//! if ttb.h or ppl.h is included prior to including cyPointCloud.h.
	void Build( SIZE_TYPE numPts, PointType const *pts, SIZE_TYPE const *customIndices ) { BuildWithFunc( numPts, [&pts](SIZE_TYPE i){ return pts[i]; }, [&customIndices](SIZE_TYPE i){ return customIndices[i]; } ); }

	//! Builds a k-d tree for the given points.
	//! The positions are stored internally, retrieved from the given function.
	//! The build is parallelized using Intel's Thread Building Library (TBB) or Microsoft's Parallel Patterns Library (PPL),
	//! if ttb.h or ppl.h is included prior to including cyPointCloud.h.
	template <typename PointPosFunc>
	void BuildWithFunc( SIZE_TYPE numPts, PointPosFunc ptPosFunc ) { BuildWithFunc(numPts, ptPosFunc, [](SIZE_TYPE i){ return i; }); }

	//! Builds a k-d tree for the given points.
	//! The positions are stored internally, along with the indices to the given array.
	//! The positions and custom indices are retrieved from the given functions.
	//! The build is parallelized using Intel's Thread Building Library (TBB) or Microsoft's Parallel Patterns Library (PPL),
	//! if ttb.h or ppl.h is included prior to including cyPointCloud.h.
	template <typename PointPosFunc, typename CustomIndexFunc>
	void BuildWithFunc( SIZE_TYPE numPts, PointPosFunc ptPosFunc, CustomIndexFunc custIndexFunc )
	{
		if ( points ) delete [] points;
		pointCount = numPts;
		if ( pointCount == 0 ) { points = nullptr; return; }
		points = new PointData[(pointCount|1)+1];
		PointData *orig = new PointData[pointCount];
		PointType boundMin( (std::numeric_limits<FType>::max)() ), boundMax( (std::numeric_limits<FType>::min)() );
		for ( SIZE_TYPE i=0; i<pointCount; i++ ) {
			PointType p = ptPosFunc(i);
			orig[i].Set( p, custIndexFunc(i) );
			for ( int j=0; j<DIMENSIONS; j++ ) {
				if ( boundMin[j] > p[j] ) boundMin[j] = p[j];
				if ( boundMax[j] < p[j] ) boundMax[j] = p[j];
			}
		}
		BuildKDTree( orig, boundMin, boundMax, 1, 0, pointCount );
		delete [] orig;
		if ( (pointCount & 1) == 0 ) {
			// if the point count is even, we should add a bogus point
			points[ pointCount+1 ].Set( PointType( std::numeric_limits<FType>::infinity() ), 0, 0 );
		}
		numInternal = pointCount / 2;
	}

	//! Returns true if the Build or BuildWithFunc methods would perform the build in parallel using multi-threading.
	//! The build is parallelized using Intel's Thread Building Library (TBB) or Microsoft's Parallel Patterns Library (PPL),
	//! if ttb.h or ppl.h are included prior to including cyPointCloud.h.
	static bool IsBuildParallel()
	{
#ifdef _CY_PARALLEL_LIB
		return true;
#else
		return false;
#endif
	}

	/////////////////////////////////////////////////////////////////////////////////
	//!@ General search methods

	//! Returns all points to the given position within the given radius.
	//! Calls the given pointFound function for each point found.
	//!
	//! The given pointFound function can reduce the radiusSquared value.
	//! However, increasing the radiusSquared value can have unpredictable results.
	//! The callback function must be in the following form:
	//!
	//! void _CALLBACK(SIZE_TYPE index, PointType const &p, FType distanceSquared, FType &radiusSquared)
	template <typename _CALLBACK>
	void GetPoints( PointType const &position, FType radius, _CALLBACK pointFound ) const
	{
		FType r2 = radius*radius;
		GetPoints( position, r2, pointFound, 1 );
	}

	//! Used by one of the PointCloud::GetPoints() methods.
	//!
	//! Keeps the point index, position, and distance squared to a given search position.
	//! Used by one of the GetPoints methods.
	struct PointInfo {
		SIZE_TYPE index;			//!< The index of the point
		PointType pos;				//!< The position of the point
		FType     distanceSquared;	//!< Squared distance from the search position
		bool operator < ( PointInfo const &b ) const { return distanceSquared < b.distanceSquared; }	//!< Comparison operator
	};

	//! Returns the closest points to the given position within the given radius.
	//! It returns the number of points found.
	int GetPoints( PointType const &position, FType radius, SIZE_TYPE maxCount, PointInfo *closestPoints ) const
	{
		int pointsFound = 0;
		GetPoints( position, radius, [&](SIZE_TYPE i, PointType const &p, FType d2, FType &r2) {
			if ( pointsFound == maxCount ) {
				std::pop_heap( closestPoints, closestPoints+maxCount );
				closestPoints[maxCount-1].index = i;
				closestPoints[maxCount-1].pos = p;
				closestPoints[maxCount-1].distanceSquared = d2;
				std::push_heap( closestPoints, closestPoints+maxCount );
				r2 = closestPoints[0].distanceSquared;
			} else {
				closestPoints[pointsFound].index = i;
				closestPoints[pointsFound].pos = p;
				closestPoints[pointsFound].distanceSquared = d2;
				pointsFound++;
				if ( pointsFound == maxCount ) {
					std::make_heap( closestPoints, closestPoints+maxCount );
					r2 = closestPoints[0].distanceSquared;
				}
			}
		} );
		return pointsFound;
	}

	//! Returns the closest points to the given position.
	//! It returns the number of points found.
	int GetPoints( PointType const &position, SIZE_TYPE maxCount, PointInfo *closestPoints ) const
	{
		return GetPoints( position, (std::numeric_limits<FType>::max)(), maxCount, closestPoints );
	}

	/////////////////////////////////////////////////////////////////////////////////
	//!@name Closest point methods

	//! Returns the closest point to the given position within the given radius.
	//! It returns true, if a point is found.
	bool GetClosest( PointType const &position, FType radius, SIZE_TYPE &closestIndex, PointType &closestPosition, FType &closestDistanceSquared ) const
	{
		bool found = false;
		FType dist2 = radius * radius;
		GetPoints( position, dist2, [&](SIZE_TYPE i, PointType const &p, FType d2, FType &r2){ found=true; closestIndex=i; closestPosition=p; closestDistanceSquared=d2; r2=d2; }, 1 );
		return found;
	}

	//! Returns the closest point to the given position.
	//! It returns true, if a point is found.
	bool GetClosest( PointType const &position, SIZE_TYPE &closestIndex, PointType &closestPosition, FType &closestDistanceSquared ) const
	{
		return GetClosest( position, (std::numeric_limits<FType>::max)(), closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point index and position to the given position within the given index.
	//! It returns true, if a point is found.
	bool GetClosest( PointType const &position, FType radius, SIZE_TYPE &closestIndex, PointType &closestPosition ) const
	{
		FType closestDistanceSquared;
		return GetClosest( position, radius, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point index and position to the given position.
	//! It returns true, if a point is found.
	bool GetClosest( PointType const &position, SIZE_TYPE &closestIndex, PointType &closestPosition ) const
	{
		FType closestDistanceSquared;
		return GetClosest( position, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point index to the given position within the given radius.
	//! It returns true, if a point is found.
	bool GetClosestIndex( PointType const &position, FType radius, SIZE_TYPE &closestIndex ) const
	{
		FType closestDistanceSquared;
		PointType closestPosition;
		return GetClosest( position, radius, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point index to the given position.
	//! It returns true, if a point is found.
	bool GetClosestIndex( PointType const &position, SIZE_TYPE &closestIndex ) const
	{
		FType closestDistanceSquared;
		PointType closestPosition;
		return GetClosest( position, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point position to the given position within the given radius.
	//! It returns true, if a point is found.
	bool GetClosestPosition( PointType const &position, FType radius, PointType &closestPosition ) const
	{
		SIZE_TYPE closestIndex;
		FType closestDistanceSquared;
		return GetClosest( position, radius, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point position to the given position.
	//! It returns true, if a point is found.
	bool GetClosestPosition( PointType const &position, PointType &closestPosition ) const
	{
		SIZE_TYPE closestIndex;
		FType closestDistanceSquared;
		return GetClosest( position, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point distance squared to the given position within the given radius.
	//! It returns true, if a point is found.
	bool GetClosestDistanceSquared( PointType const &position, FType radius, FType &closestDistanceSquared ) const
	{
		SIZE_TYPE closestIndex;
		PointType closestPosition;
		return GetClosest( position, radius, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point distance squared to the given position.
	//! It returns true, if a point is found.
	bool GetClosestDistanceSquared( PointType const &position, FType &closestDistanceSquared ) const
	{
		SIZE_TYPE closestIndex;
		PointType closestPosition;
		return GetClosest( position, closestIndex, closestPosition, closestDistanceSquared );
	}

	/////////////////////////////////////////////////////////////////////////////////

private:

	/////////////////////////////////////////////////////////////////////////////////
	//!@name Internal Structures and Methods

	class PointData
	{
	private:
		SIZE_TYPE indexAndSplitPlane;	// first NBits bits indicates the splitting plane, the rest of the bits store the point index.
		PointType p;					// point position
	public:
		void Set( PointType const &pt, SIZE_TYPE index, uint32_t plane=0 ) { p=pt; indexAndSplitPlane = (index<<NBits()) | (plane&((1<<NBits())-1)); }
		void SetPlane( uint32_t plane ) { indexAndSplitPlane = (indexAndSplitPlane & (~((SIZE_TYPE(1)<<NBits())-1))) | plane; }
		int       Plane() const { return indexAndSplitPlane & ((1<<NBits())-1); }
		SIZE_TYPE Index() const { return indexAndSplitPlane >> NBits(); }
		PointType const & Pos() const { return p; }
	private:
#if defined(__cpp_constexpr) || (defined(_MSC_VER) && _MSC_VER >= 1900)
		constexpr int NBits(uint32_t v=DIMENSIONS) const { return v < 2 ? v : 1+NBits(v>>1); }
#else
		int NBits() const { int v = DIMENSIONS-1, r, s; r=(v>0xF)<<2; v>>=r; s=(v>0x3)<<1; v>>=s; r|=s|(v>>1); return r+1; }	// Supports up to 256 dimensions
#endif
	};

	PointData *points;		// Keeps the points as a k-d tree.
	SIZE_TYPE  pointCount;	// Keeps the point count.
	SIZE_TYPE  numInternal;	// Keeps the number of internal k-d tree nodes.

	// The main method for recursively building the k-d tree.
	void BuildKDTree( PointData *orig, PointType boundMin, PointType boundMax, SIZE_TYPE kdIndex, SIZE_TYPE ixStart, SIZE_TYPE ixEnd )
	{
		SIZE_TYPE n = ixEnd - ixStart;
		if ( n > 1 ) {
			int axis = SplitAxis( boundMin, boundMax );
			SIZE_TYPE leftSize = LeftSize(n);
			SIZE_TYPE ixMid = ixStart+leftSize;
			std::nth_element( orig+ixStart, orig+ixMid, orig+ixEnd, [axis](PointData const &a, PointData const &b){ return a.Pos()[axis] < b.Pos()[axis]; } );
			points[kdIndex] = orig[ixMid];
			points[kdIndex].SetPlane( axis );
			PointType bMax = boundMax;
			bMax[axis] = orig[ixMid].Pos()[axis];
			PointType bMin = boundMin;
			bMin[axis] = orig[ixMid].Pos()[axis];
#ifdef _CY_PARALLEL_LIB
			SIZE_TYPE const parallel_invoke_threshold = 256;
			if ( ixMid-ixStart > parallel_invoke_threshold && ixEnd - ixMid+1 > parallel_invoke_threshold ) {
				_CY_PARALLEL_LIB::parallel_invoke(
					[&]{ BuildKDTree( orig, boundMin, bMax, kdIndex*2,   ixStart, ixMid ); },
					[&]{ BuildKDTree( orig, bMin, boundMax, kdIndex*2+1, ixMid+1, ixEnd ); }
				);
			} else 
#endif
			{
				BuildKDTree( orig, boundMin, bMax, kdIndex*2,   ixStart, ixMid );
				BuildKDTree( orig, bMin, boundMax, kdIndex*2+1, ixMid+1, ixEnd );
			}
		} else if ( n > 0 ) {
			points[kdIndex] = orig[ixStart];
		}
	}

	// Returns the total number of nodes on the left sub-tree of a complete k-d tree of size n.
	static SIZE_TYPE LeftSize( SIZE_TYPE n )
	{
		SIZE_TYPE f = n; // Size of the full tree
		for ( SIZE_TYPE s=1; s<8*sizeof(SIZE_TYPE); s*=2 ) f |= f >> s;
		SIZE_TYPE l = f >> 1; // Size of the full left child
		SIZE_TYPE r = l >> 1; // Size of the full right child without leaf nodes
		return (l+r+1 <= n) ? l : n-r-1;
	}

	// Returns axis with the largest span, used as the splitting axis for building the k-d tree
	static int SplitAxis( PointType const &boundMin, PointType const &boundMax )
	{
		PointType d = boundMax - boundMin;
		int axis = 0;
		FType dmax = d[0];
		for ( int j=1; j<DIMENSIONS; j++ ) {
			if ( dmax < d[j] ) {
				axis = j;
				dmax = d[j];
			}
		}
		return axis;
	}

	template <typename _CALLBACK>
	void GetPoints( PointType const &position, FType &dist2, _CALLBACK pointFound, SIZE_TYPE nodeID ) const
	{
		SIZE_TYPE stack[sizeof(SIZE_TYPE)*8];
		SIZE_TYPE stackPos = 0;

		TraverseCloser( position, dist2, pointFound, nodeID, stack, stackPos );

		// empty the stack
		while ( stackPos > 0 ) {
			SIZE_TYPE nodeID = stack[ --stackPos ];
			// check the internal node point
			PointData const &p = points[nodeID];
			PointType const pos = p.Pos();
			int axis = p.Plane();
			float dist1 = position[axis] - pos[axis];
			if ( dist1*dist1 < dist2 ) {
				// check its point
				FType d2 = (position - pos).LengthSquared();
				if ( d2 < dist2 ) pointFound( p.Index(), pos, d2, dist2 );
				// traverse down the other child node
				SIZE_TYPE child = 2*nodeID;
				nodeID = dist1 < 0 ? child+1 : child;
				TraverseCloser( position, dist2, pointFound, nodeID, stack, stackPos );
			}
		}
	}

	template <typename _CALLBACK>
	void TraverseCloser( PointType const &position, FType &dist2, _CALLBACK pointFound, SIZE_TYPE nodeID, SIZE_TYPE *stack, SIZE_TYPE &stackPos ) const
	{
		// Traverse down to a leaf node along the closer branch
		while ( nodeID <= numInternal ) {
			stack[stackPos++] = nodeID;
			PointData const &p = points[nodeID];
			PointType const pos = p.Pos();
			int axis = p.Plane();
			float dist1 = position[axis] - pos[axis];
			SIZE_TYPE child = 2*nodeID;
			nodeID = dist1 < 0 ? child : child + 1;
		}
		// Now we are at a leaf node, do the test
		PointData const &p = points[nodeID];
		PointType const pos = p.Pos();
		FType d2 = (position - pos).LengthSquared();
		if ( d2 < dist2 ) pointFound( p.Index(), pos, d2, dist2 );
	}

	/////////////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

#ifdef _CY_VECTOR_H_INCLUDED_
template <typename T> _CY_TEMPLATE_ALIAS( PointCloud2, (PointCloud<Vec2<T>,T,2>) );	//!< A 2D point cloud using a k-d tree
template <typename T> _CY_TEMPLATE_ALIAS( PointCloud3, (PointCloud<Vec3<T>,T,3>) );	//!< A 3D point cloud using a k-d tree
template <typename T> _CY_TEMPLATE_ALIAS( PointCloud4, (PointCloud<Vec4<T>,T,4>) );	//!< A 4D point cloud using a k-d tree

typedef PointCloud<Vec2f,float,2>  PointCloud2f;	//!< A 2D point cloud using a k-d tree with float  type elements
typedef PointCloud<Vec3f,float,3>  PointCloud3f;	//!< A 3D point cloud using a k-d tree with float  type elements
typedef PointCloud<Vec4f,float,4>  PointCloud4f;	//!< A 4D point cloud using a k-d tree with float  type elements

typedef PointCloud<Vec2d,double,2> PointCloud2d;	//!< A 2D point cloud using a k-d tree with double type elements
typedef PointCloud<Vec3d,double,3> PointCloud3d;	//!< A 3D point cloud using a k-d tree with double type elements
typedef PointCloud<Vec4d,double,4> PointCloud4d;	//!< A 4D point cloud using a k-d tree with double type elements

template <typename T, uint32_t DIMENSIONS> _CY_TEMPLATE_ALIAS( PointCloudN, (PointCloud<Vec<T,DIMENSIONS>,T,DIMENSIONS>) );	//!< A multi-dimensional point cloud using a k-d tree
template <uint32_t DIMENSIONS> _CY_TEMPLATE_ALIAS( PointCloudNf , (PointCloudN<float,   DIMENSIONS>) );	//!< A multi-dimensional point cloud using a k-d tree with single precision (float)
template <uint32_t DIMENSIONS> _CY_TEMPLATE_ALIAS( PointCloudNd , (PointCloudN<double,  DIMENSIONS>) );	//!< A multi-dimensional point cloud using a k-d tree with double precision (double)
#endif

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

#ifdef _CY_VECTOR_H_INCLUDED_
template <typename T> _CY_TEMPLATE_ALIAS( cyPointCloud2, (cy::PointCloud<cy::Vec2<T>,T,2>) );	//!< A 2D point cloud using a k-d tree
template <typename T> _CY_TEMPLATE_ALIAS( cyPointCloud3, (cy::PointCloud<cy::Vec3<T>,T,3>) );	//!< A 3D point cloud using a k-d tree
template <typename T> _CY_TEMPLATE_ALIAS( cyPointCloud4, (cy::PointCloud<cy::Vec4<T>,T,4>) );	//!< A 4D point cloud using a k-d tree

typedef cy::PointCloud<cy::Vec2f,float,2>  cyPointCloud2f;	//!< A 2D point cloud using a k-d tree with float  type elements
typedef cy::PointCloud<cy::Vec3f,float,3>  cyPointCloud3f;	//!< A 3D point cloud using a k-d tree with float  type elements
typedef cy::PointCloud<cy::Vec4f,float,4>  cyPointCloud4f;	//!< A 4D point cloud using a k-d tree with float  type elements

typedef cy::PointCloud<cy::Vec2d,double,2> cyPointCloud2d;	//!< A 2D point cloud using a k-d tree with double type elements
typedef cy::PointCloud<cy::Vec3d,double,3> cyPointCloud3d;	//!< A 3D point cloud using a k-d tree with double type elements
typedef cy::PointCloud<cy::Vec4d,double,4> cyPointCloud4d;	//!< A 4D point cloud using a k-d tree with double type elements

template <typename T, uint32_t DIMENSIONS> _CY_TEMPLATE_ALIAS( cyPointCloudN, (cy::PointCloud<cy::Vec<T,DIMENSIONS>,T,DIMENSIONS>) );	//!< A multi-dimensional point cloud using a k-d tree
template <uint32_t DIMENSIONS> _CY_TEMPLATE_ALIAS( cyPointCloudNf , (cyPointCloudN<float,   DIMENSIONS>) );	//!< A multi-dimensional point cloud using a k-d tree with float  type elements
template <uint32_t DIMENSIONS> _CY_TEMPLATE_ALIAS( cyPointCloudNd , (cyPointCloudN<double,  DIMENSIONS>) );	//!< A multi-dimensional point cloud using a k-d tree with double type elements
#endif

//-------------------------------------------------------------------------------

#ifdef _CY_POP_MACRO_max
# pragma pop_macro("max")
# undef _CY_POP_MACRO_max
#endif

#ifdef _CY_POP_MACRO_min
# pragma pop_macro("min")
# undef _CY_POP_MACRO_min
#endif

//-------------------------------------------------------------------------------

#endif
