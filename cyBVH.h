// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyBVH.h 
//! \author Cem Yuksel
//! 
//! \brief  Bounding Volume Hierarchy class.
//!
//! BVH is a storage class for Bounding Volume Hierarchies.
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

#ifndef _CY_BVH_H_INCLUDED_
#define _CY_BVH_H_INCLUDED_

//-------------------------------------------------------------------------------

#include <vector>

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

#ifndef CY_BVH_ELEMENT_COUNT_BITS
#define CY_BVH_ELEMENT_COUNT_BITS   3	//!< Determines the number of bits needed to represent the maximum number of elements in a node (8)
#endif

#ifndef CY_BVH_MAX_ELEMENT_COUNT
#define CY_BVH_MAX_ELEMENT_COUNT    (1<<CY_BVH_ELEMENT_COUNT_BITS)	//!< Determines the maximum number of elements in a node (8)
#endif

#define _CY_BVH_NODE_DATA_BITS      (sizeof(uint32_t)*8)
#define _CY_BVH_ELEMENT_COUNT_MASK  ((1<<CY_BVH_ELEMENT_COUNT_BITS)-1)
#define _CY_BVH_LEAF_BIT_MASK       ((uint32_t)1<<(_CY_BVH_NODE_DATA_BITS-1))
#define _CY_BVH_CHILD_ID_BITS       (_CY_BVH_NODE_DATA_BITS-1)
#define _CY_BVH_CHILD_ID_MASK       (_CY_BVH_LEAF_BIT_MASK-1)
#define _CY_BVH_ELEMENT_OFFSET_BITS (_CY_BVH_NODE_DATA_BITS-1-CY_BVH_ELEMENT_COUNT_BITS)
#define _CY_BVH_ELEMENT_OFFSET_MASK ((1<<_CY_BVH_ELEMENT_OFFSET_BITS)-1)

//-------------------------------------------------------------------------------

//! Bounding Volume Hierarchy class

class BVH
{
public:

	/////////////////////////////////////////////////////////////////////////////////
	//@ Node Access Methods
	/////////////////////////////////////////////////////////////////////////////////

	//! Returns the id of the root node.
	uint32_t GetRootNodeID() const { return 1; }

	//! Returns the bounding box of the node as 6 float values.
	//! The first 3 values are the minimum x, y, and z coordinates and
	//! the last 3 values are the maximum x, y, and z coordinates of the box.
	float const * GetNodeBounds( uint32_t nodeID ) const { return nodes[nodeID].GetBounds(); }

	//! Returns true if the node is a leaf node.
	bool IsLeafNode( uint32_t nodeID ) const { return nodes[nodeID].IsLeafNode(); }

	//! Returns the id of the child node with the given index (parent must be an internal node).
	//! The child index can only be 0 or 1, since this is a binary BVH.
	uint32_t GetFirstChildNode( uint32_t parentNodeID, int childIndex ) const { return nodes[parentNodeID].ChildID()+childIndex; }

	//! Returns the id of the first child node (parent must be an internal node).
	uint32_t GetFirstChildNode( uint32_t parentNodeID ) const { return nodes[parentNodeID].ChildID(); }

	//! Returns the id of the second child node (parent must be an internal node).
	uint32_t GetSecondChildNode( uint32_t parentNodeID ) const { return nodes[parentNodeID].ChildID()+1; }

	//! Given the first child node id, returns the id of the second child node.
	uint32_t GetSiblingNode( uint32_t firstChildNodeID ) const { return firstChildNodeID+1; }

	//! Returns the child nodes of the given node (parent must be an internal node).
	void GetChildNodes( uint32_t parent, uint32_t &child1, uint32_t &child2 ) const
	{
		child1 = GetFirstChildNode(parent);
		child2 = GetSiblingNode(child1);
	}

	//! Returns the number of elements inside the given node (must be a leaf node).
	uint32_t GetNodeElementCount( uint32_t nodeID ) const  { return nodes[nodeID].ElementCount(); }

	//! Returns the list of element inside the given node (must be a leaf node).
	uint32_t const * GetNodeElements( uint32_t nodeID ) const { return &elements[nodes[nodeID].ElementOffset()]; }

	/////////////////////////////////////////////////////////////////////////////////
	//@ Clear and Build Methods
	/////////////////////////////////////////////////////////////////////////////////

	//! Clears the tree structure
	void Clear()
	{
		nodes.clear();
		elements.clear();
	}

	//! Builds the tree structure by recursively splitting the nodes. maxElementsPerNode cannot be larger than 8.
	//! The given functions should be in the following forms:
	//! void getElementBounds( uint32_t i, float box[6] )
	//! uint32_t findSplit( uint32_t elementCount, uint32_t *_elements, float const *box, int maxElementsPerNode )
	template <typename TGetElementBounds, typename TFindSplit>
	void Build( TGetElementBounds getElementBounds, TFindSplit findSplit, uint32_t numElements, int maxElementsPerNode=CY_BVH_MAX_ELEMENT_COUNT )
	{
		if ( numElements <= 0 ) return Clear();
		SetNumElements( numElements );
		for ( uint32_t i=0; i<numElements; i++ ) SetElement( i, i );
		BuildElements<TGetElementBounds,TFindSplit>( getElementBounds, findSplit, maxElementsPerNode );
	}

	//! Builds the tree structure by recursively splitting the nodes at the centers of bounding boxes. maxElementsPerNode cannot be larger than 8.
	//! The given functions should be in the following forms:
	//! void getElementBounds( uint32_t i, float box[6] )
	//! float getElementCenter( uint32_t i, int dimension )
	template <typename TGetElementBounds, typename TGetElementCenter>
	void BuildCenterSplits( TGetElementBounds getElementBounds, TGetElementCenter getElementCenter, uint32_t numElements, int maxElementsPerNode=CY_BVH_MAX_ELEMENT_COUNT )
	{
		Build( getElementBounds, 
			[&]( uint32_t elementCount, uint32_t *_elements, float const *box, int maxElementsPerNode ) {
				return CenterSplit<TGetElementCenter>(getElementCenter,elementCount,_elements,box,maxElementsPerNode);
			}, numElements, maxElementsPerNode );
	}

	//! Sets the number of elements prior to calling BuildElements.
	void SetNumElements( uint32_t numElements ) { elements.resize(numElements); }

	//! Set the element data prior to calling BuildElements.
	//! Calling this function after the BVH is built can break the BVH data.
	void SetElement( uint32_t id, uint32_t data ) { elements[id] = data; }

	//! Builds the tree structure using the previously set elements by recursively splitting the nodes. maxElementsPerNode cannot be larger than 8.
	//! The given functions should be in the following forms:
	//! void getElementBounds( uint32_t i, float box[6] )
	//! uint32_t findSplit( uint32_t elementCount, uint32_t *_elements, float const *box, int maxElementsPerNode )
	template <typename TGetElementBounds, typename TFindSplit>
	void BuildElements( TGetElementBounds getElementBounds, TFindSplit findSplit, int maxElementsPerNode=CY_BVH_MAX_ELEMENT_COUNT )
	{
		if ( maxElementsPerNode > CY_BVH_MAX_ELEMENT_COUNT ) maxElementsPerNode = CY_BVH_MAX_ELEMENT_COUNT;
		uint32_t numElements = (uint32_t) elements.size();
		Box box;
		box.Init();
		for ( uint32_t i=0; i<numElements; i++ ) {
			Box b;
			getElementBounds(i,b.b);
			box += b;
		}
		TempNode *tempRoot = new TempNode( numElements, 0, box );
		SplitTempNode( getElementBounds, findSplit, tempRoot, maxElementsPerNode );
		uint32_t numNodes = tempRoot->GetNumNodes();
		nodes.resize( numNodes + 1 );
		ConvertTempData( 1, tempRoot, 2 );
		delete tempRoot;
	}

	//! Builds the tree structure using the previously set elements by recursively splitting the nodes at the centers of bounding boxes. maxElementsPerNode cannot be larger than 8.
	//! The given functions should be in the following forms:
	//! void getElementBounds( uint32_t i, float box[6] )
	//! float getElementCenter( uint32_t i, int dimension )
	template <typename TGetElementBounds, typename TGetElementCenter>
	void BuildElementsCenterSplits( TGetElementBounds getElementBounds, TGetElementCenter getElementCenter, int maxElementsPerNode=CY_BVH_MAX_ELEMENT_COUNT )
	{
		BuildElements( getElementBounds, 
			[&]( uint32_t elementCount, uint32_t *_elements, float const *box, int maxElementsPerNode ) {
				return CenterSplit<TGetElementCenter>(getElementCenter,elementCount,_elements,box,maxElementsPerNode);
			}, maxElementsPerNode );
	}

	//! Recomputes the bounding boxes of all BVH nodes without modifying the tree structure.
	//! BVH must already be built. The given function should be in the following form:
	//! void getElementBounds( uint32_t i, float box[6] )
	template <typename TGetElementBounds>
	void Refit( TGetElementBounds getElementBounds )
	{
		for ( uint32_t i=(uint32_t)nodes.size()-1; i>0; --i ) {
			Box box;
			if ( IsLeafNode(i) ) {
				uint32_t n = GetNodeElementCount(i);
				uint32_t const *elem = GetNodeElements(i);
				box.Init();
				for ( uint32_t j=0; j<n; ++j ) {
					Box b;
					getElementBounds( elem[j], b.b );
					box += b;
				}
			} else {
				box  = GetNodeBounds( GetFirstChildNode (i) );
				box += GetNodeBounds( GetSecondChildNode(i) );
			}
			nodes[i].SetBounds(box);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////

private:

	/////////////////////////////////////////////////////////////////////////////////
	//@ Internal storage
	/////////////////////////////////////////////////////////////////////////////////

	struct Box
	{
		float b[6];
		Box() { Init(); }
		Box( Box const &box ) { for(int i=0; i<6; i++) b[i]=box.b[i]; }
		void Init() { b[0]=b[1]=b[2]=1e30f; b[3]=b[4]=b[5]=-1e30f; }
		void operator  = ( float const *x ) { for(int i=0; i<6; i++) { b[i]=x[i]; } }
		void operator += ( float const *x ) { for(int i=0; i<3; i++) { b[i]=Min(b[i],    x[i]); b[i+3]=Max(b[i+3],    x[i+3]); } }
		void operator += ( Box const &box ) { for(int i=0; i<3; i++) { b[i]=Min(b[i],box.b[i]); b[i+3]=Max(b[i+3],box.b[i+3]); } }
	};

	class Node
	{
	public:
		void SetLeafNode( Box const &bound, uint32_t elemCount, uint32_t elemOffset ) { box=bound; data=(elemOffset&_CY_BVH_ELEMENT_OFFSET_MASK)|((elemCount-1)<<_CY_BVH_ELEMENT_OFFSET_BITS)|_CY_BVH_LEAF_BIT_MASK; }
		void SetInternalNode( Box const &bound, uint32_t chilID ) { box=bound; data=(chilID&_CY_BVH_CHILD_ID_MASK); }
		uint32_t      ChildID      () const { return (data&_CY_BVH_CHILD_ID_MASK); }										//!< returns the id to the first child (must be internal node)
		uint32_t      ElementOffset() const { return (data&_CY_BVH_ELEMENT_OFFSET_MASK); }									//!< returns the offset to the first element (must be leaf node)
		uint32_t      ElementCount () const { return ((data>>_CY_BVH_ELEMENT_OFFSET_BITS)&_CY_BVH_ELEMENT_COUNT_MASK)+1; }	//!< returns the number of elements in this node (must be leaf node)
		bool          IsLeafNode   () const { return (data&_CY_BVH_LEAF_BIT_MASK)>0; }										//!< returns true if this is a leaf node
		float const * GetBounds    () const { return box.b; }																//!< returns the bounding box of the node
		void          SetBounds    ( Box const &bound ) { box = bound; }													//!< updates the bounds of the node
	private:
		Box      box;	//!< bounding box of the node
		uint32_t data;	//!< node data bits that keep the leaf node flag and the child node id or element count and element offset.
	};

	std::vector<Node>     nodes;	//!< the tree structure that keeps all the node data (nodeData[0] is not used for cache coherency)
	std::vector<uint32_t> elements;	//!< indices of all elements in all nodes

	/////////////////////////////////////////////////////////////////////////////////
	//@ Internal methods for building the BVH tree
	/////////////////////////////////////////////////////////////////////////////////

	//! Temporary node class used for building the hierarchy and then converted to NodeData.
	class TempNode
	{
	public:
		TempNode( uint32_t count, uint32_t offset, Box const &boundBox) : child1(0), child2(0), elementCount(count), elementOffset(offset), box(boundBox) {}
		~TempNode() { if ( child1 ) delete child1; if ( child2 ) delete child2; }

		void Split( uint32_t child1ElementCount, Box const &child1Box, Box const &child2Box )
		{
			child1 = new TempNode(child1ElementCount,elementOffset,child1Box);
			child2 = new TempNode(ElementCount()-child1ElementCount,elementOffset+child1ElementCount,child2Box);
		}
		uint32_t GetNumNodes() const
		{
			uint32_t n = 1;
			if ( child1 ) n += child1->GetNumNodes();
			if ( child2 ) n += child2->GetNumNodes();
			return n;
		}
		bool      IsLeafNode   () const { return child1==0; }
		uint32_t  ElementCount () const { return elementCount; }
		uint32_t  ElementOffset() const { return elementOffset; }
		TempNode* GetChild1    ()       { return child1; }
		TempNode* GetChild2    ()       { return child2; }
		Box const & GetBounds  () const { return box; }
	private:
		TempNode *child1, *child2;
		Box      box;
		uint32_t elementCount;
		uint32_t elementOffset;
	};

	//! Recursively splits the given temporary node.
	template <typename TGetElementBounds, typename TFindSplit>
	void SplitTempNode( TGetElementBounds getElementBounds, TFindSplit findSplit, TempNode *tNode, int maxElementsPerNode )
	{
		float const *box = tNode->GetBounds().b;
		uint32_t *nodeElements = &elements[tNode->ElementOffset()];
		uint32_t child1ElemCount = findSplit(tNode->ElementCount(),nodeElements,box,maxElementsPerNode);

		// If the FindSplit call does not return a valid split position
		if ( child1ElemCount == 0 || child1ElemCount >= tNode->ElementCount() ) {
			// if we must split anyway
			if ( tNode->ElementCount() > CY_BVH_MAX_ELEMENT_COUNT ) {
				// we split in half arbitrarily.
				child1ElemCount = tNode->ElementCount() / 2;
			} else {
				// otherwise, we reached a leaf node and no more split is necessary.
				return;
			}
		}

		// Compute child bounding boxes
		Box child1Box;
		Box child2Box;
		for ( uint32_t i=0; i<child1ElemCount; i++ ) {
			Box eBox;
			getElementBounds( nodeElements[i], eBox.b );
			child1Box += eBox;
		}
		for ( uint32_t i=child1ElemCount; i<tNode->ElementCount(); i++ ) {
			Box eBox;
			getElementBounds( nodeElements[i], eBox.b );
			child2Box += eBox;
		}

		// Split recursively
		tNode->Split( child1ElemCount, child1Box, child2Box );
		SplitTempNode<TGetElementBounds,TFindSplit>(getElementBounds,findSplit,tNode->GetChild1(),maxElementsPerNode);
		SplitTempNode<TGetElementBounds,TFindSplit>(getElementBounds,findSplit,tNode->GetChild2(),maxElementsPerNode);
	}

	//! Recursively converts the temporary node data to NodeData.
	uint32_t ConvertTempData( uint32_t nodeID, TempNode *tNode, uint32_t childID )
	{
		if ( tNode->IsLeafNode() ) {
			nodes[nodeID].SetLeafNode( tNode->GetBounds(), tNode->ElementCount(), tNode->ElementOffset() );
			return childID;
		} else {
			nodes[nodeID].SetInternalNode( tNode->GetBounds(), childID );
			uint32_t newChildID = ConvertTempData( childID, tNode->GetChild1(), childID+2 );
			return ConvertTempData( childID+1, tNode->GetChild2(), newChildID );
		}
	}

	//! Called by the default implementation of FindSplit.
	//! Splits the elements using the widest axis of the given bounding box.
	template <typename TGetElementCenter>
	static uint32_t CenterSplit( TGetElementCenter getElementCenter, uint32_t elementCount, uint32_t *nodeElements, float const *box, int maxElementsPerNode )
	{
		if ( elementCount <= (uint32_t)maxElementsPerNode ) return 0;
		float d[3] = { box[3]-box[0], box[4]-box[1], box[5]-box[2] };
		int sd[3]; // split dimensions
		sd[0] = d[0] >= d[1] ? ( d[0] >= d[2] ? 0 : 2 ) : ( d[1] >= d[2] ? 1 : 2 );
		sd[1] = (sd[0]+1) % 3;
		sd[2] = (sd[0]+2) % 3;
		if ( d[sd[1]] < d[sd[2]] ) { int t=sd[1]; sd[1]=sd[2]; sd[2]=t; }

		uint32_t child1ElemCount = 0;
		for ( int s=0; s<3; s++ ) {
			int splitDim = sd[s];
			float splitPos = 0.5f * ( box[splitDim] + box[splitDim+3] );
			uint32_t i=0, j=elementCount;
			while ( i<j ) {
				float center = getElementCenter( nodeElements[i], splitDim );
				if ( center <= splitPos ) {
					i++;
				} else {
					j--;
					uint32_t t = nodeElements[i];
					nodeElements[i] = nodeElements[j];
					nodeElements[j] = t;
				}
			}
			if ( i < elementCount && i > 0 ) {
				child1ElemCount = i;
				break;
			}
		}

		return child1ElemCount;
	}

	/////////////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

#ifdef _CY_TRIMESH_H_INCLUDED_

//! Bounding Volume Hierarchy for triangular meshes (TriMesh)

class BVHTriMesh : public BVH
{
public:
	//!@name Constructors
	BVHTriMesh() : mesh(0) {}
	BVHTriMesh( TriMesh const *m ) { SetMesh(m); }

	//! Sets the mesh pointer and builds the BVH structure.
	void SetMesh( TriMesh const *m, int maxElementsPerNode=CY_BVH_MAX_ELEMENT_COUNT )
	{
		mesh = m;
		Clear();
		auto getElementBounds = [&]( uint32_t i, float box[6] )
		{
			TriMesh::TriFace const &f = mesh->F(i);
			cyVec3f p = mesh->V( f.v[0] );
			box[0]=box[3]=p.x; box[1]=box[4]=p.y; box[2]=box[5]=p.z;
			for ( int j=1; j<3; j++ ) { // for each triangle
				cyVec3f q = mesh->V( f.v[j] );
				for ( int k=0; k<3; k++ ) { // for each dimension
					if ( box[k] > q[k] ) box[k] = q[k];
					if ( box[k+3] < q[k] ) box[k+3] = q[k];
				}
			}
		};
		auto getElementCenter = [&]( uint32_t i, int dim )
		{
			TriMesh::TriFace const &f = mesh->F(i);
			return ( mesh->V(f.v[0])[dim] + mesh->V(f.v[1])[dim] + mesh->V(f.v[2])[dim] ) / 3.0f;
		};
		BuildCenterSplits( getElementBounds, getElementCenter, mesh->NF(), maxElementsPerNode );
	}

private:
	TriMesh const *mesh;
};

#endif

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::BVH cyBVH;	//!< Bounding Volume Hierarchy class

#ifdef _CY_TRIMESH_H_INCLUDED_
typedef cy::BVHTriMesh cyBVHTriMesh;	//!< BVH hierarchy for triangular meshes (TriMesh)
#endif

//-------------------------------------------------------------------------------

#endif

