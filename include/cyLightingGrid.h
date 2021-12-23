// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyLightingGrid.h 
//! \author Cem Yuksel
//! 
//! \brief  Implementation of the Lighting Grid Hierarchy method.
//! 
//! This file includes an implementation of the Lighting Grid Hierarchy 
//! method for efficiently handling many lights.
//!
//! The Lighting Grid Hierarchy method was originally designed for explosion
//! rendering, but it can be used for general purpose lighting computations
//! involving a large number of (point) lights. It is particularly effective
//! when shadows can be precomputed.
//! 
//! More details can be found in the original publication:
//!
//! Can Yuksel and Cem Yuksel. 2017. Lighting Grid Hierarchy for
//! Self-illuminating Explosions. ACM Transactions on Graphics
//! (Proceedings of SIGGRAPH 2017) 36, 4, Article 110 (July 2017).
//! http://www.cemyuksel.com/research/lgh/
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

#ifndef _CY_LIGHTING_GRID_H_INCLUDED_
#define _CY_LIGHTING_GRID_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyPoint.h"
#include "cyIPoint.h"
#include "cyColor.h"
#include "cyPointCloud.h"
#include <random>

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

//#define CY_LIGHTING_GRID_ORIG_POS

//-------------------------------------------------------------------------------

//! An implementation of the Lighting Grid Hierarchy method.
//!
//! Can Yuksel and Cem Yuksel. 2017. Lighting Grid Hierarchy for
//! Self-illuminating Explosions. ACM Transactions on Graphics
//! (Proceedings of SIGGRAPH 2017) 36, 4, Article 110 (July 2017).
//! http://www.cemyuksel.com/research/lgh/
//!
//! This class builds the Lighting Grid Hierarchy and uses it for lighting.

class LightingGridHierarchy
{
public:
	LightingGridHierarchy() : levels(nullptr), numLevels(0) {}	//!< Constructor
	virtual ~LightingGridHierarchy() { Clear(); }				//!< Destructor

	int   GetNumLevels() const { return numLevels; }	//!< Returns the number of levels in the hierarchy.
	float GetCellSize () const { return cellSize; }		//!< Returns the size of a cell in the lowest (finest) level of the hierarchy.
	const Point3f& GetLightPos   ( int level, int i  ) const { return levels[level].pc.GetPoint(i); }							//!< Returns the i^th light position at the given level. Note that this is not the position of the light with index i.
	const int      GetLightIndex ( int level, int i  ) const { return levels[level].pc.GetPointIndex(i); }						//!< Returns the i^th light index at the given level.
	const Color&   GetLightIntens( int level, int ix ) const { return levels[level].colors[ix]; }								//!< Returns the intensity of the light with index ix at the given level.
	const Point3f& GetLightPosDev( int level, int ix ) const { return level > 0 ? levels[level].pDev[ix] : Point3f(0,0,0); }	//!< Returns the position variation of the light with index ix at the given level.

	void Clear() { if ( levels ) delete [] levels; levels = nullptr; numLevels = 0; }	//!< Deletes all data.

	//! Builds the Lighting Grid Hierarchy for the given point light positions and intensities using the given parameters.
	//! This method builds the hierarchy using the given cellSize as the size of the lowest (finest) level grid cells.
	bool Build( const Point3f *lightPos,			//!< Light positions.
	            const Color   *lightIntensities, 	//!< Light intensities.
	            int            numLights, 			//!< Number of lights.
	            int            minLevelLights,		//!< The minimum number of lights permitted for the highest (coarsest) level of the hierarchy. The build stops if a higher (coarser) level would have fewer lights.
	            float          cellSize,			//!< The size of a grid cell in the lowest (finest) level of the hierarchy.
	            int            highestLevel			//!< The highest level permitted, where level 0 contains the original lights.
	          )
	{
		return DoBuild( lightPos, lightIntensities, numLights, 0, minLevelLights, cellSize, highestLevel );
	}

	//! Builds the Lighting Grid Hierarchy for the given point light positions and intensities using the given parameters.
	//! This method automatically determines the grid cell size based on the bounding box of the light positions.
	bool Build( const Point3f *lightPos,			//!< Light positions.
	            const Color   *lightIntensities, 	//!< Light intensities.
	            int            numLights, 			//!< Number of lights.
	            int            minLevelLights, 		//!< The minimum number of lights permitted for the highest (coarsest) level of the hierarchy. The build stops if a higher (coarser) level would have fewer lights.
	            float          autoFitScale = 1.01f	//!< Extends the bounding box of the light positions using the given scale. This value must be 1 or greater. A value slightly greater than 1 often provides a good fit for the grid.
	          )
	{
		return DoBuild( lightPos, lightIntensities, numLights, autoFitScale, minLevelLights );
	}

	//! Computes the illumination at the given position using the given accuracy parameter alpha.
	template <typename LightingFunction>
	void Light( const Point3f    &pos,						//!< The position where the lighting will be evaluated.
	            float            alpha,						//!< The accuracy parameter. It should be 1 or greater. Larger values produce more accurate results with substantially more computation.
	            LightingFunction lightingFunction			//!< This function is called for each light used for lighting computation. It should be in the form void LightingFunction(int level, int light_id, const Point3f &light_position, const Color &light_intensity).
	          )
	{
		Light( pos, alpha, 0, lightingFunction );
	}

	//! Computes the illumination at the given position using the given accuracy parameter alpha.
	//! This method provides stochastic sampling by randomly changing the given light position when calling the lighting function.
	template <typename LightingFunction>
	void Light( const Point3f    &pos,						//!< The position where the lighting will be evaluated.
	            float            alpha,						//!< The accuracy parameter. It should be 1 or greater. Larger values produce more accurate results with substantially more computation.
	            int              stochasticShadowSamples,   //!< When this parameter is zero, the given lightingFunction is called once per light, using the position of the light. Otherwise, it is called as many times as this parameter specifies, using random positions around each light position.
	            LightingFunction lightingFunction			//!< This function is called for each light used for lighting computation. It should be in the form void LightingFunction(int level, int light_id, const Point3f &light_position, const Color &light_intensity).
	          )
	{
		if ( numLevels > 1 ) {

			// First level
			float r = alpha * cellSize;
			float rr = r * r;
			levels[0].pc.GetPoints( pos, r*2, [&](int i, const Point3f &p, float dist2, float &radius2) {
				Color c = levels[0].colors[i];
				if ( dist2 > rr ) c *= 1 - (sqrtf(dist2)-r)/r;
				lightingFunction( 0, i, p, c );
			} );

			auto callLightingFunc = [&]( int level, int i, const Point3f &p, const Color &c )
			{
				if ( stochasticShadowSamples > 0 ) {
					Color cc = c / (float) stochasticShadowSamples;
					for ( int j=0; j<stochasticShadowSamples; j++ ) {
						Point3f pj = p + RandomPos() * levels[level].pDev[i];
						lightingFunction( level, i, pj, cc );
					}
				} else {
					lightingFunction( level, i, p, c );
				}
			};

			// Middle levels
			for ( int level=1; level<numLevels-1; level++ ) {
				float r_min = r;
				float rr_min = r * r;
				r *= 2;
				levels[level].pc.GetPoints( pos, r*2, [&](int i, const Point3f &p, float dist2, float &radius2) {
					if ( dist2 <= rr_min ) return;
					Color c = levels[level].colors[i];
					float d = sqrtf(dist2);
					if ( d > r ) c *= 1 - (d-r)/r;
					else c *= (d-r_min)/r_min;
					callLightingFunc( level, i, p, c );
				} );
			}

			// Last level
			float r_min = r;
			float rr_min = r * r;
			r *= 2;
			rr = r * r;
			int n = levels[numLevels-1].pc.GetPointCount();
			for ( int i=0; i<n; i++ ) {
				const Point3f &p = levels[numLevels-1].pc.GetPoint(i);
				float dist2 = (pos - p).LengthSquared();
				if ( dist2 <= rr_min ) continue;
				int id = levels[numLevels-1].pc.GetPointIndex(i);
				Color c = levels[numLevels-1].colors[id];
				if ( dist2 < rr ) c *= (sqrtf(dist2)-r_min)/r_min;
				callLightingFunc( numLevels-1, i, p, c );
			}

		} else {
			// Single-level (a.k.a. brute-force)
			int n = levels[0].pc.GetPointCount();
			for ( int i=0; i<n; i++ ) {
				const Point3f &p = levels[0].pc.GetPoint(i);
				int id = levels[0].pc.GetPointIndex(i);
				Color c = levels[0].colors[id];
				lightingFunction( 0, i, p, c );
			}
		}
	}

private:
	struct Level {
		Level() : colors(nullptr), pDev(nullptr) {}
		~Level() { delete [] colors; delete [] pDev; }
		cy::PointCloud<Point3f,float,3,int> pc;
		Color   *colors;
		Point3f *pDev; // position deviation for random shadow sampling
	};
	Level *levels;
	int    numLevels;
	float  cellSize;

	float RandomX()
	{
		static thread_local std::mt19937 generator;
		std::uniform_real_distribution<float> distribution;
		float x = distribution(generator);
		float y = distribution(generator);
		if ( y > (cosf(x*Pi<float>())+1)*0.5f ) x -= 1;
		return x;
	}

	Point3f RandomPos()
	{
		Point3f p;
		p.x = RandomX();
		p.y = RandomX();
		p.z = RandomX();
		return p;
	}

	bool DoBuild( const Point3f *lightPos, const Color *lightColor, int numLights, float autoFitScale, int minLevelLights, float cellSize=0, int highestLevel=10 )
	{
		Clear();
		if ( numLights <= 0 || highestLevel <= 0 ) return false;

		// Compute the bounding box for the lighPoss
		Point3f boundMin = Point3f(lightPos[0]);
		Point3f boundMax = Point3f(lightPos[0]);
		for ( int i=1; i<numLights; i++ ) {
			for ( int d=0; d<3; d++ ) {
				if ( boundMin[d] > lightPos[i][d] ) boundMin[d] = lightPos[i][d];
				if ( boundMax[d] < lightPos[i][d] ) boundMax[d] = lightPos[i][d];
			}
		}
		Point3f boundDif = boundMax - boundMin;
		float boundDifMin = boundDif.Min();

		// Determine the actual highest level
		float highestCellSize;
		cy::IPoint3i highestGridRes;
		if ( autoFitScale > 0 ) {
			highestCellSize = boundDif.Max() * autoFitScale;
			int s = int(1.0f/autoFitScale) + 2;
			if ( s < 2 ) s = 2;
			highestGridRes.Set(s,s,s);
		} else {
			int highestLevelMult = 1 << (highestLevel-1);
			highestCellSize = cellSize * highestLevelMult;
			while ( highestLevel>1 && highestCellSize > boundDifMin*2 ) {
				highestLevel--;
				highestLevelMult = 1 << (highestLevel-1);
				highestCellSize = cellSize * highestLevelMult;
			}
			highestGridRes = cy::IPoint3i(boundDif / highestCellSize) + 2;
		}

		struct Node {
			Node() : position(0,0,0), color(0,0,0), weight(0), firstChild(-1) {}
			Point3f position;
			Color color;
#ifdef CY_LIGHTING_GRID_ORIG_POS
			Point3f origPos;
#endif // CY_LIGHTING_GRID_ORIG_POS
			Point3f stdev;
			float weight;
			int firstChild;
			void AddLight( float w, const Point3f &p, const Color &c )
			{
				weight   += w;
				position += w * p;
				color    += w * c;
				stdev    += w * (p*p);
			}
			void Normalize()
			{ 
				if ( weight > 0 ) {
					position /= weight;
					stdev = stdev/weight - position*position;
				}
			}
		};

		// Allocate the temporary nodes
		numLevels = highestLevel+1;
		std::vector< std::vector<Node> > nodes(numLevels);

		auto gridIndex = []( cy::IPoint3i &index, const Point3f &pos, float cellSize )
		{
			Point3f normP = pos / cellSize;
			index = cy::IPoint3i(normP);
			return normP - Point3f(index);
		};

		auto addLightToNodes = []( std::vector<Node> &nds, const int nodeIDs[8], const Point3f &interp, const Point3f &light_pos, const Color &light_color )
		{
			for ( int j=0; j<8; j++ ) {
				float w = ((j&1) ? interp.x : (1-interp.x)) * ((j&2) ? interp.y : (1-interp.y)) * ((j&4) ? interp.z : (1-interp.z));
				nds[ nodeIDs[j] ].AddLight( w, light_pos, light_color );
			}
		};

		// Generate the grid for the highest level
		Point3f highestGridSize = Point3f(highestGridRes-1) * highestCellSize;
		Point3f center = (boundMax + boundMin) / 2;
		Point3f corner = center - highestGridSize/2;
		nodes[highestLevel].resize( highestGridRes.x * highestGridRes.y * highestGridRes.z );
#ifdef CY_LIGHTING_GRID_ORIG_POS
		for ( int z=0, j=0; z<highestGridRes.z; z++ ) {
			for ( int y=0; y<highestGridRes.y; y++ ) {
				for ( int x=0; x<highestGridRes.x; x++, j++ ) {
					nodes[highestLevel][j].origPos = corner + Point3f(x,y,z)*highestCellSize;
				}
			}
		}
#endif // CY_LIGHTING_GRID_ORIG_POS
		for ( int i=0; i<numLights; i++ ) {
			cy::IPoint3i index;
			Point3f interp = gridIndex( index, Point3f(lightPos[i])-corner, highestCellSize );
			int is = index.z*highestGridRes.y*highestGridRes.x + index.y*highestGridRes.x + index.x;
			int nodeIDs[8] = {
				is,
				is + 1,
				is + highestGridRes.x,
				is + highestGridRes.x + 1,
				is + highestGridRes.x*highestGridRes.y,
				is + highestGridRes.x*highestGridRes.y + 1,
				is + highestGridRes.x*highestGridRes.y + highestGridRes.x,
				is + highestGridRes.x*highestGridRes.y + highestGridRes.x + 1,
			};
			for ( int j=0; j<8; j++ ) assert( nodeIDs[j] >= 0 && nodeIDs[j] < (int)nodes[highestLevel].size() );
			addLightToNodes( nodes[highestLevel], nodeIDs, interp, lightPos[i], lightColor[i] );
		}
		for ( int i=0; i<(int)nodes[highestLevel].size(); i++ ) nodes[highestLevel][i].Normalize();

		// Generate the lower levels
		float nodeCellSize = highestCellSize;
		cy::IPoint3i gridRes = highestGridRes;
		int levelSkip = 0;
		for ( int level=highestLevel-1; level>0; level-- ) {
			// Find the number of nodes for this level
			int nodeCount = 0;
			for ( int i=0; i<(int)nodes[level+1].size(); i++ ) {
				if ( nodes[level+1][i].weight > 0 ) {
					nodes[level+1][i].firstChild = nodeCount;
					nodeCount += 8;
				}
			}

			if ( nodeCount > numLights/4 ) {
				levelSkip = level;
				break;
			}

			nodes[level].resize( nodeCount );
			// Add the lights to the nodes
			nodeCellSize /= 2;
			gridRes *= 2;
#ifdef CY_LIGHTING_GRID_ORIG_POS
			for ( int i=0; i<(int)nodes[level+1].size(); i++ ) {
				int fc = nodes[level+1][i].firstChild;
				if ( fc < 0 ) continue;
				for ( int z=0, j=0; z<2; z++ ) {
					for ( int y=0; y<2; y++ ) {
						for ( int x=0; x<2; x++, j++ ) {
							nodes[level][fc+j].origPos = nodes[level+1][i].origPos + Point3f(x,y,z)*nodeCellSize;
						}
					}
				}
			}
#endif // CY_LIGHTING_GRID_ORIG_POS
			for ( int i=0; i<numLights; i++ ) {
				cy::IPoint3i index;
				Point3f interp = gridIndex( index, Point3f(lightPos[i])-corner, nodeCellSize );
				// find the node IDs
				int nodeIDs[8];
				index <<= level+2;
				for ( int z=0, j=0; z<2; z++ ) {
					int iz = index.z + z;
					for ( int y=0; y<2; y++ ) {
						int iy = index.y + y;
						for ( int x=0; x<2; x++, j++ ) {
							int ix = index.x + x;
							int hix = ix >> (highestLevel+2);
							int hiy = iy >> (highestLevel+2);
							int hiz = iz >> (highestLevel+2);
							int nid = hiz*highestGridRes.y*highestGridRes.x + hiy*highestGridRes.x + hix;
							for ( int l=highestLevel-1; l>=level; l-- ) {
								int ii = ((index.z >> l)&4) | ((index.y >> (l+1))&2) |  ((index.x >> (l+2))&1);
								assert( nodes[l+1][nid].firstChild >= 0 );
								nid = nodes[l+1][nid].firstChild + ii;
								assert( nid >= 0 && nid < (int)nodes[l].size() );
							}
							nodeIDs[j] = nid;
						}
					}
				}
				addLightToNodes( nodes[level], nodeIDs, interp, lightPos[i], lightColor[i] );
			}
			for ( int i=0; i<(int)nodes[level].size(); i++ ) nodes[level][i].Normalize();
		}

		// Copy light data
		numLevels = highestLevel + 1 - levelSkip;
		int levelBaseSkip = 0;
		// Skip levels that have two few lights (based on minLevelLights).
		for ( int level=1; level<numLevels; level++ ) {
			std::vector<Node> &levelNodes = nodes[level+levelSkip];
			int count = 0;
			for ( int i=0; i<(int)levelNodes.size(); i++ ) {
				if ( levelNodes[i].weight > 0 ) {
					count++;
				}
			}
			if ( count < minLevelLights ) {
				numLevels = level;
				break;
			}
		}

		levels = new Level[ numLevels ];
		for ( int level=1; level<numLevels; level++ ) {
			std::vector<Node> &levelNodes = nodes[level+levelSkip];
			Level &thisLevel = levels[level];
			std::vector<Point3f> pos( levelNodes.size() );
			int lightCount = 0;
			for ( int i=0; i<(int)levelNodes.size(); i++ ) {
				if ( levelNodes[i].weight > 0 ) {
					pos[lightCount++] = levelNodes[i].position;
				}
			}
			thisLevel.pc.Build( lightCount, pos.data() );
			thisLevel.colors = new Color[ lightCount ];
			thisLevel.pDev = new Point3f[ lightCount ];
			for ( int i=0, j=0; i<(int)levelNodes.size(); i++ ) {
				if ( levelNodes[i].weight > 0 ) {
					assert( j < lightCount );
					thisLevel.colors[j] = levelNodes[i].color;
					thisLevel.pDev[j].x = sqrtf( levelNodes[i].stdev.x ) * Pi<float>();
					thisLevel.pDev[j].y = sqrtf( levelNodes[i].stdev.y ) * Pi<float>();
					thisLevel.pDev[j].z = sqrtf( levelNodes[i].stdev.z ) * Pi<float>();
					j++;
				}
			}
			levelNodes.resize(0);
			levelNodes.shrink_to_fit();
		}
		std::vector<Point3f> pos( numLights );
		levels[0].colors = new Color[ numLights ];
		for ( int i=0; i<numLights; i++ ) {
			pos[i] = lightPos[i];
			levels[0].colors[i] = lightColor[i];
		}
		levels[0].pc.Build( numLights, pos.data() );
		this->cellSize = nodeCellSize;

		return true;
	}
};

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

#endif
