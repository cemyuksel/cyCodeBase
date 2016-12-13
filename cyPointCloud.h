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
# undef max
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
	//////////////////////////////////////////////////////////////////////////!//!//!
	//!@name Constructors and Destructor

	PointCloud() : points(nullptr), pointCount(0) {}
	PointCloud( SIZE_TYPE numPts, const PointType *pts, const SIZE_TYPE *customIndices=nullptr ) : points(nullptr), pointCount(0) { Build(numPts,pts,customIndices); }
	~PointCloud() { delete [] points; }

	//////////////////////////////////////////////////////////////////////////!//!//!
	//!@ Initialization

	//! Builds a k-d tree for the given points.
	//! The point locations are stored internally, along with the indices to the given array.
	void Build( SIZE_TYPE numPts, const PointType *pts, const SIZE_TYPE *customIndices=nullptr )
	{
		if ( points ) delete [] points;
		pointCount = numPts;
		if ( pointCount == 0 ) { points = nullptr; return; }
		points = new PointData[pointCount+1];
		SIZE_TYPE *order = new SIZE_TYPE[pointCount];
		for ( SIZE_TYPE i=0; i<pointCount; i++ ) order[i] = i;
		BuildKDTree( pts, customIndices, order, 1, 0, pointCount );
		delete [] order;
	}

	//////////////////////////////////////////////////////////////////////////!//!//!
	//!@ General search methods

	//! Returns all points to the given position within the given radius.
	//! Calls the given pointFound function for each point found.
	//!
	//! The given pointFound function can reduce the radiusSquared value.
	//! However, increasing the radiusSquared value can have unpredictable results.
	//! The callback function must be in the following form:
	//!
	//! void _CALLBACK(SIZE_TYPE index, const PointType &p, FType distanceSquared, FType &radiusSquared)
	template <typename _CALLBACK>
	void GetPoints( const PointType &position, FType radius, _CALLBACK pointFound )
	{
		SIZE_TYPE internalNodes = (pointCount+1) >> 1;
		SIZE_TYPE stack[60];	// deep enough for 2^30 points
		int stackPos = 0;
		stack[0] = 1;	// root node
		FType dist2 = radius * radius;
		while ( stackPos >= 0 ) {
			SIZE_TYPE ix = stack[ stackPos-- ];
			const PointData *p = &points[ix];
			const PointType pos = p->Pos();
			if ( ix < internalNodes ) {
				int axis = p->Plane();
				FType d = position[axis] - pos[axis];
				if( d > 0 ) {	// if dist1 is positive search right child first
					stack[++stackPos] = 2*ix + 1;
					if ( d*d < dist2 ) stack[++stackPos] = 2*ix;
				} else {	// dist1 is negative, search left child first
					stack[++stackPos] = 2*ix;
					if ( d*d < dist2 ) stack[++stackPos] = 2*ix + 1;
				}
			}
			FType d2 = (pos - position).LengthSquared();
			if ( d2 < dist2 ) pointFound( p->Index(), pos, d2, dist2 );
		}
		if ( (pointCount & SIZE_TYPE(1)) == 0 ) {
			const PointData *p = &points[pointCount];
			FType d2 = (p->Pos() - position).LengthSquared();
			if ( d2 < dist2 ) pointFound( p->Index(), p->Pos(), d2, dist2 );
		}
	}

	//! Used by one of the PointCloud::GetPoints() methods.
	//!
	//! Keeps the point index, position, and distance squared to a given search position.
	//! Used by one of the GetPoints methods.
	struct PointInfo {
		SIZE_TYPE index;			//!< The index of the point
		PointType pos;				//!< The position of the point
		FType     distanceSquared;	//!< Squared distance from the search position
		bool operator < (const PointInfo &b) const { return distanceSquared < b.distanceSquared; }	//!< Comparison operator
	};

	//! Returns the closest points to the given position within the given radius.
	//! The returned value is the number of points found.
	int GetPoints( const PointType &position, FType radius, SIZE_TYPE maxCount, PointInfo *closestPoints )
	{
		bool tooManyPoints = false;
		int pointsFound = 0;
		GetPoints( position, radius, [&](SIZE_TYPE i, const PointType &p, FType d2, FType &r2) {
			if ( pointsFound == maxCount ) {
				if ( !tooManyPoints ) {
					std::make_heap( closestPoints, closestPoints+maxCount );
				}
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
			}
		} );
		return pointsFound;
	}

	//! Returns the closest points to the given position.
	//! The returned value is the number of points found.
	int GetPoints( const PointType &position, SIZE_TYPE maxCount, PointInfo *closestPoints )
	{
		return GetPoints( position, std::numeric_limits<FType>::max(), maxCount, closestPoints );
	}

	//////////////////////////////////////////////////////////////////////////!//!//!
	//!@name Closest point methods

	//! Returns the closest point to the given position within the given radius.
	//! The returned value is true, if a point is found.
	bool GetClosest( const PointType &position, FType radius, SIZE_TYPE &closestIndex, PointType &closestPosition, FType &closestDistanceSquared )
	{
		bool found = false;
		GetPoints( position, radius, [&](SIZE_TYPE i, const PointType &p, FType d2, FType &r2){ found=true; closestIndex=i; closestPosition=p; closestDistanceSquared=d2; r2=d2; } );
		return found;
	}

	//! Returns the closest point to the given position.
	//! The returned value is true, if a point is found.
	bool GetClosest( const PointType &position, SIZE_TYPE &closestIndex, PointType &closestPosition, FType &closestDistanceSquared )
	{
		return GetClosest( position, std::numeric_limits<FType>::max(), closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point index and position to the given position within the given index.
	//! The returned value is true, if a point is found.
	bool GetClosest( const PointType &position, FType radius, SIZE_TYPE &closestIndex, PointType &closestPosition )
	{
		FType closestDistanceSquared;
		return GetClosest( position, radius, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point index and position to the given position.
	//! The returned value is true, if a point is found.
	bool GetClosest( const PointType &position, SIZE_TYPE &closestIndex, PointType &closestPosition )
	{
		FType closestDistanceSquared;
		return GetClosest( position, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point index to the given position within the given radius.
	//! The returned value is true, if a point is found.
	bool GetClosestIndex( const PointType &position, FType radius, SIZE_TYPE &closestIndex )
	{
		FType closestDistanceSquared;
		PointType closestPosition;
		return GetClosest( position, radius, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point index to the given position.
	//! The returned value is true, if a point is found.
	bool GetClosestIndex( const PointType &position, SIZE_TYPE &closestIndex )
	{
		FType closestDistanceSquared;
		PointType closestPosition;
		return GetClosest( position, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point position to the given position within the given radius.
	//! The returned value is true, if a point is found.
	bool GetClosestPosition( const PointType &position, FType radius, PointType &closestPosition )
	{
		SIZE_TYPE closestIndex;
		FType closestDistanceSquared;
		return GetClosest( position, radius, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point position to the given position.
	//! The returned value is true, if a point is found.
	bool GetClosestPosition( const PointType &position, PointType &closestPosition )
	{
		SIZE_TYPE closestIndex;
		FType closestDistanceSquared;
		return GetClosest( position, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point distance squared to the given position within the given radius.
	//! The returned value is true, if a point is found.
	bool GetClosestDistanceSquared( const PointType &position, FType radius, FType &closestDistanceSquared )
	{
		SIZE_TYPE closestIndex;
		PointType closestPosition;
		return GetClosest( position, radius, closestIndex, closestPosition, closestDistanceSquared );
	}

	//! Returns the closest point distance squared to the given position.
	//! The returned value is true, if a point is found.
	bool GetClosestDistanceSquared( const PointType &position, FType &closestDistanceSquared )
	{
		SIZE_TYPE closestIndex;
		PointType closestPosition;
		return GetClosest( position, closestIndex, closestPosition, closestDistanceSquared );
	}

	//////////////////////////////////////////////////////////////////////////!//!//!

private:

	//////////////////////////////////////////////////////////////////////////!//!//!
	//!@name Internal Structures and Methods

	class PointData
	{
	private:
		SIZE_TYPE indexAndSplitPlane;	// first NBits bits indicates the splitting plane, the rest of the bits store the point index.
		PointType p;					// point position
	public:
		void Set( const PointType &pt, SIZE_TYPE index, uint32_t plane=0 ) { p=pt; indexAndSplitPlane = (index<<NBits()) | (plane&((1<<NBits())-1)); }
		int       Plane() const { return indexAndSplitPlane & ((1<<NBits())-1); }
		SIZE_TYPE Index() const { return indexAndSplitPlane >> NBits(); }
		const PointType& Pos() const { return p; }
	private:
		constexpr int NBits(uint32_t v=DIMENSIONS) const { return v < 2 ? v : 1+NBits(v>>1); }
	};

	PointData *points;		// Keeps the points as a k-d tree.
	SIZE_TYPE  pointCount;	// Keeps the point count.

	// The main method for recursively building the k-d tree.
	void BuildKDTree( const PointType *pts, const SIZE_TYPE *indices, SIZE_TYPE *order, SIZE_TYPE kdIndex, SIZE_TYPE ixStart, SIZE_TYPE ixEnd )
	{
		SIZE_TYPE n = ixEnd - ixStart;
		if ( n <= 1 ) {
			if ( n > 0 ) {
				SIZE_TYPE ix = order[ixStart];
				if ( indices ) points[kdIndex].Set( pts[ix], indices[ix] );
				else           
					points[kdIndex].Set( pts[ix], ix );
			}
		} else {
			int axis = SplitAxis( pts, order, ixStart, ixEnd );
			SIZE_TYPE leftSize = LeftSize(n);
			SIZE_TYPE ixMid = ixStart+leftSize;
			std::nth_element( order+ixStart, order+ixMid, order+ixEnd, [&pts,axis](const int &a, const int &b){ return pts[a][axis] < pts[b][axis]; } );
			SIZE_TYPE ix = order[ixMid];
			if ( indices ) points[kdIndex].Set( pts[ix], indices[ix], axis );
			else           
				points[kdIndex].Set( pts[ix], ix, axis );
			BuildKDTree( pts, indices, order, kdIndex*2,   ixStart, ixMid );
			BuildKDTree( pts, indices, order, kdIndex*2+1, ixMid+1, ixEnd );
		}
	}

	// Returns the total number of nodes on the left sub-tree of a complete k-d tree of size n.
	SIZE_TYPE LeftSize( SIZE_TYPE n )
	{
		SIZE_TYPE f = n; // Size of the full tree
		for ( SIZE_TYPE s=1; s<8*sizeof(SIZE_TYPE); s*=2 ) f |= f >> s;
		SIZE_TYPE l = f >> 1; // Size of the full left child
		SIZE_TYPE r = l >> 1; // Size of the full right child without leaf nodes
		return (l+r+1 <= n) ? l : n-r-1;
	}

	// Returns axis with the largest span, used as the splitting axis for building the k-d tree
	int SplitAxis( const PointType *pts, SIZE_TYPE *indices, SIZE_TYPE ixStart, SIZE_TYPE ixEnd )
	{
		PointType box_min = pts[ indices[ixStart] ];
		PointType box_max = box_min;
		for ( SIZE_TYPE i=ixStart+1; i<ixEnd; i++ ) {
			PointType p = pts[ indices[i] ];
			for ( SIZE_TYPE d=0; d<DIMENSIONS; d++ ) {
				if ( box_min[d] > p[d] ) box_min[d] = p[d];
				if ( box_max[d] < p[d] ) box_max[d] = p[d];
			}
		}
		int axis = 0;
		{
			FType axisSize = box_max[0] - box_min[0];
			for ( SIZE_TYPE d=1; d<DIMENSIONS; d++ ) {
				FType s = box_max[d] - box_min[d];
				if ( axisSize < s ) {
					axis = d;
					axisSize = s;
				}
			}
		}
		return axis;
	}

	//////////////////////////////////////////////////////////////////////////!//!//!
};

//-------------------------------------------------------------------------------

#ifdef _CY_POINT_H_INCLUDED_
template <typename TYPE> using PointCloud2 = PointCloud<Point2<TYPE>,TYPE,2>;	//!< A 2D point cloud using a k-d tree
template <typename TYPE> using PointCloud3 = PointCloud<Point2<TYPE>,TYPE,3>;	//!< A 3D point cloud using a k-d tree
template <typename TYPE> using PointCloud4 = PointCloud<Point2<TYPE>,TYPE,4>;	//!< A 4D point cloud using a k-d tree

typedef PointCloud<Point2f,float,2>  PointCloud2f;	//!< A 2D point cloud using a k-d tree with single precision (float)
typedef PointCloud<Point3f,float,3>  PointCloud3f;	//!< A 3D point cloud using a k-d tree with single precision (float)
typedef PointCloud<Point4f,float,4>  PointCloud4f;	//!< A 4D point cloud using a k-d tree with single precision (float)

typedef PointCloud<Point2f,double,2> PointCloud2d;	//!< A 2D point cloud using a k-d tree with double precision (double)
typedef PointCloud<Point3f,double,3> PointCloud3d;	//!< A 3D point cloud using a k-d tree with double precision (double)
typedef PointCloud<Point4f,double,4> PointCloud4d;	//!< A 4D point cloud using a k-d tree with double precision (double)

typedef PointCloud<Point2f,int32_t,2>  PointCloud2i;	//!< A 2D point cloud using a k-d tree with 32-bit signed integer (int32_t)
typedef PointCloud<Point3f,int32_t,3>  PointCloud3i;	//!< A 3D point cloud using a k-d tree with 32-bit signed integer (int32_t)
typedef PointCloud<Point4f,int32_t,4>  PointCloud4i;	//!< A 4D point cloud using a k-d tree with 32-bit signed integer (int32_t)

typedef PointCloud<Point2f,uint32_t,2> PointCloud2ui;	//!< A 2D point cloud using a k-d tree with 32-bit unsigned integer (uint32_t)
typedef PointCloud<Point3f,uint32_t,3> PointCloud3ui;	//!< A 3D point cloud using a k-d tree with 32-bit unsigned integer (uint32_t)
typedef PointCloud<Point4f,uint32_t,4> PointCloud4ui;	//!< A 4D point cloud using a k-d tree with 32-bit unsigned integer (uint32_t)

typedef PointCloud<Point2f,int64_t,2>  PointCloud2l;	//!< A 2D point cloud using a k-d tree with 64-bit signed integer (int64_t)
typedef PointCloud<Point3f,int64_t,3>  PointCloud3l;	//!< A 3D point cloud using a k-d tree with 64-bit signed integer (int64_t)
typedef PointCloud<Point4f,int64_t,4>  PointCloud4l;	//!< A 4D point cloud using a k-d tree with 64-bit signed integer (int64_t)

typedef PointCloud<Point2f,uint64_t,2> PointCloud2ul;	//!< A 2D point cloud using a k-d tree with 64-bit unsigned integer (uint64_t)
typedef PointCloud<Point3f,uint64_t,3> PointCloud3ul;	//!< A 3D point cloud using a k-d tree with 64-bit unsigned integer (uint64_t)
typedef PointCloud<Point4f,uint64_t,4> PointCloud4ul;	//!< A 4D point cloud using a k-d tree with 64-bit unsigned integer (uint64_t)

template <typename TYPE, uint32_t DIMENSIONS> using PointCloudN = PointCloud<Point<TYPE,DIMENSIONS>,TYPE,DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree
template <uint32_t DIMENSIONS> using PointCloudNf  = PointCloudN<float,   DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with single precision (float)
template <uint32_t DIMENSIONS> using PointCloudNd  = PointCloudN<double,  DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with double precision (double)
template <uint32_t DIMENSIONS> using PointCloudNi  = PointCloudN<int32_t, DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with 32-bit signed integer (int32_t)
template <uint32_t DIMENSIONS> using PointCloudNui = PointCloudN<uint32_t,DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with 32-bit unsigned integer (uint32_t)
template <uint64_t DIMENSIONS> using PointCloudNl  = PointCloudN<int64_t, DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with 64-bit signed integer (int64_t)
template <uint64_t DIMENSIONS> using PointCloudNul = PointCloudN<uint64_t,DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with 64-bit unsigned integer (uint64_t)
#endif

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

#ifdef _CY_POINT_H_INCLUDED_
template <typename TYPE> using cyPointCloud2 = cy::PointCloud<cy::Point2<TYPE>,TYPE,2>;	//!< A 2D point cloud using a k-d tree
template <typename TYPE> using cyPointCloud3 = cy::PointCloud<cy::Point2<TYPE>,TYPE,3>;	//!< A 3D point cloud using a k-d tree
template <typename TYPE> using cyPointCloud4 = cy::PointCloud<cy::Point2<TYPE>,TYPE,4>;	//!< A 4D point cloud using a k-d tree

typedef cy::PointCloud<cy::Point2f,float,2>  cyPointCloud2f;	//!< A 2D point cloud using a k-d tree with single precision (float)
typedef cy::PointCloud<cy::Point3f,float,3>  cyPointCloud3f;	//!< A 3D point cloud using a k-d tree with single precision (float)
typedef cy::PointCloud<cy::Point4f,float,4>  cyPointCloud4f;	//!< A 4D point cloud using a k-d tree with single precision (float)

typedef cy::PointCloud<cy::Point2f,double,2> cyPointCloud2d;	//!< A 2D point cloud using a k-d tree with double precision (double)
typedef cy::PointCloud<cy::Point3f,double,3> cyPointCloud3d;	//!< A 3D point cloud using a k-d tree with double precision (double)
typedef cy::PointCloud<cy::Point4f,double,4> cyPointCloud4d;	//!< A 4D point cloud using a k-d tree with double precision (double)

typedef cy::PointCloud<cy::Point2f,int32_t,2>  cyPointCloud2i;	//!< A 2D point cloud using a k-d tree with 32-bit signed integer (int32_t)
typedef cy::PointCloud<cy::Point3f,int32_t,3>  cyPointCloud3i;	//!< A 3D point cloud using a k-d tree with 32-bit signed integer (int32_t)
typedef cy::PointCloud<cy::Point4f,int32_t,4>  cyPointCloud4i;	//!< A 4D point cloud using a k-d tree with 32-bit signed integer (int32_t)

typedef cy::PointCloud<cy::Point2f,uint32_t,2> cyPointCloud2ui;	//!< A 2D point cloud using a k-d tree with 32-bit unsigned integer (uint32_t)
typedef cy::PointCloud<cy::Point3f,uint32_t,3> cyPointCloud3ui;	//!< A 3D point cloud using a k-d tree with 32-bit unsigned integer (uint32_t)
typedef cy::PointCloud<cy::Point4f,uint32_t,4> cyPointCloud4ui;	//!< A 4D point cloud using a k-d tree with 32-bit unsigned integer (uint32_t)

typedef cy::PointCloud<cy::Point2f,int64_t,2>  cyPointCloud2l;	//!< A 2D point cloud using a k-d tree with 64-bit signed integer (int64_t)
typedef cy::PointCloud<cy::Point3f,int64_t,3>  cyPointCloud3l;	//!< A 3D point cloud using a k-d tree with 64-bit signed integer (int64_t)
typedef cy::PointCloud<cy::Point4f,int64_t,4>  cyPointCloud4l;	//!< A 4D point cloud using a k-d tree with 64-bit signed integer (int64_t)

typedef cy::PointCloud<cy::Point2f,uint64_t,2> cyPointCloud2ul;	//!< A 2D point cloud using a k-d tree with 64-bit unsigned integer (uint64_t)
typedef cy::PointCloud<cy::Point3f,uint64_t,3> cyPointCloud3ul;	//!< A 3D point cloud using a k-d tree with 64-bit unsigned integer (uint64_t)
typedef cy::PointCloud<cy::Point4f,uint64_t,4> cyPointCloud4ul;	//!< A 4D point cloud using a k-d tree with 64-bit unsigned integer (uint64_t)

template <typename TYPE, uint32_t DIMENSIONS> using cyPointCloudN = cy::PointCloud<cy::Point<TYPE,DIMENSIONS>,TYPE,DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree
template <uint32_t DIMENSIONS> using cyPointCloudNf  = cyPointCloudN<float,   DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with single precision (float)
template <uint32_t DIMENSIONS> using cyPointCloudNd  = cyPointCloudN<double,  DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with double precision (double)
template <uint32_t DIMENSIONS> using cyPointCloudNi  = cyPointCloudN<int32_t, DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with 32-bit signed integer (int32_t)
template <uint32_t DIMENSIONS> using cyPointCloudNui = cyPointCloudN<uint32_t,DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with 32-bit unsigned integer (uint32_t)
template <uint64_t DIMENSIONS> using cyPointCloudNl  = cyPointCloudN<int64_t, DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with 64-bit signed integer (int64_t)
template <uint64_t DIMENSIONS> using cyPointCloudNul = cyPointCloudN<uint64_t,DIMENSIONS>;	//!< A multi-dimensional point cloud using a k-d tree with 64-bit unsigned integer (uint64_t)
#endif

//-------------------------------------------------------------------------------

#endif
