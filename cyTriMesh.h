// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyTriMesh.h 
//! \author Cem Yuksel
//! 
//! \brief  Triangular Mesh class.
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

#ifndef _CY_TRIMESH_H_INCLUDED_
#define _CY_TRIMESH_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyPoint.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <vector>

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

//! Triangular Mesh Class

class TriMesh
{
public:
	//! Triangular Mesh Face
	struct TriFace
	{
		unsigned int v[3];	//!< vertex indices
	};

	//! Material definition
	struct Mtl
	{
		//! Texture map information
		struct Map
		{
			char name[256];	//!< filename of the texture map
			Map() { name[0] = '\0'; }
		};
		char name[256];	//!< Material name
		float Ka[3];	//!< Ambient color
		float Kd[3];	//!< Diffuse color
		float Ks[3];	//!< Specular color
		float Tf[3];	//!< Transmission color
		float Ns;		//!< Specular exponent
		float Ni;		//!< Index of refraction
		int illum;		//!< Illumination model
		Map map_Ka;		//!< Ambient texture map
		Map map_Kd;		//!< Diffuse texture map
		Map map_Ks;		//!< Specular texture map

		Mtl()
		{
			name[0] = '\0';
			Ka[0]=Ka[1]=Ka[2]=0;
			Kd[0]=Kd[1]=Kd[2]=1;
			Ks[0]=Ks[1]=Ks[2]=0;
			Tf[0]=Tf[1]=Tf[2]=0;
			Ns=0;
			Ni=1;
			illum=2;
		}
	};

protected:
	Point3f *v;		//!< vertices
	TriFace *f;		//!< faces
	Point3f *vn;	//!< vertex normal
	TriFace *fn;	//!< normal faces
	Point3f *vt;	//!< texture vertices
	TriFace *ft;	//!< texture faces
	Mtl     *m;		//!< materials
	int     *mcfc;	//!< material cumulative face count

	unsigned int nv;	//!< number of vertices
	unsigned int nf;	//!< number of faces
	unsigned int nvn;	//!< number of vertex normals
	unsigned int nvt;	//!< number of texture vertices
	unsigned int nm;	//!< number of materials

	//!@{
	//! Bounding box
	Point3f boundMin, boundMax;
	//!@}

public:

	//!@name Constructor and destructor
	TriMesh() : v(NULL), f(NULL), vn(NULL), fn(NULL), vt(NULL), ft(NULL), m(NULL), mcfc(NULL)
				, nv(0), nf(0), nvn(0), nvt(0), nm(0),boundMin(0,0,0), boundMax(0,0,0) {}
	virtual ~TriMesh() { Clear(); }

	//!@name Component Access Methods
	const Point3f& V (int i) const { return v[i]; }		//!< returns the i^th vertex
	Point3f&       V (int i)       { return v[i]; }		//!< returns the i^th vertex
	const TriFace& F (int i) const { return f[i]; }		//!< returns the i^th face
	TriFace&       F (int i)       { return f[i]; }		//!< returns the i^th face
	const Point3f& VN(int i) const { return vn[i]; }	//!< returns the i^th vertex normal
	Point3f&       VN(int i)       { return vn[i]; }	//!< returns the i^th vertex normal
	const TriFace& FN(int i) const { return fn[i]; }	//!< returns the i^th normal face
	TriFace&       FN(int i)       { return fn[i]; }	//!< returns the i^th normal face
	const Point3f& VT(int i) const { return vt[i]; }	//!< returns the i^th vertex texture
	Point3f&       VT(int i)       { return vt[i]; }	//!< returns the i^th vertex texture
	const TriFace& FT(int i) const { return ft[i]; }	//!< returns the i^th texture face
	TriFace&       FT(int i)       { return ft[i]; }	//!< returns the i^th texture face
	const Mtl&     M (int i) const { return m[i]; }		//!< returns the i^th material
	Mtl&           M (int i)       { return m[i]; }		//!< returns the i^th material

	unsigned int NV () const { return nv; }		//!< returns the number of vertices
	unsigned int NF () const { return nf; }		//!< returns the number of faces
	unsigned int NVN() const { return nvn; }	//!< returns the number of vertex normals
	unsigned int NVT() const { return nvt; }	//!< returns the number of texture vertices
	unsigned int NM () const { return nm; }		//!< returns the number of materials

	bool HasNormals() const { return NVN() > 0; }			//!< returns true if the mesh has vertex normals
	bool HasTextureVertices() const { return NVT() > 0; }	//!< returns true if the mesh has texture vertices

	//!@name Set Component Count
	void Clear() { SetNumVertex(0); SetNumFaces(0); SetNumNormals(0); SetNumTexVerts(0); SetNumMtls(0); boundMin.Zero(); boundMax.Zero(); }
	void SetNumVertex  (unsigned int n) { Allocate(n,v,nv); }
	void SetNumFaces   (unsigned int n) { Allocate(n,f,nf); if (fn||vn) Allocate(n,fn); if (ft||vt) Allocate(n,ft); }
	void SetNumNormals (unsigned int n) { Allocate(n,vn,nvn); if (!fn) Allocate(nf,fn); }
	void SetNumTexVerts(unsigned int n) { Allocate(n,vt,nvt); if (!ft) Allocate(nf,ft); }
	void SetNumMtls    (unsigned int n) { Allocate(n,m,nm); Allocate(n,mcfc); }

	//!@name Get Property Methods
	bool    IsBoundBoxReady() const { return boundMin.x!=0 && boundMin.y!=0 && boundMin.z!=0 && boundMax.x!=0 && boundMax.y!=0 && boundMax.z!=0; }
	Point3f GetBoundMin() const { return boundMin; }		//!< Returns the minimum values of the bounding box
	Point3f GetBoundMax() const { return boundMax; }		//!< Returns the maximum values of the bounding box
	Point3f GetPoint   (int faceID, const Point3f &bc) const { return Interpolate(faceID,v,f,bc); }	//!< Returns the point on the given face with the given barycentric coordinates (bc).
	Point3f GetNormal  (int faceID, const Point3f &bc) const { return Interpolate(faceID,vn,fn,bc); }	//!< Returns the the surface normal on the given face at the given barycentric coordinates (bc). The returned vector is not normalized.
	Point3f GetTexCoord(int faceID, const Point3f &bc) const { return Interpolate(faceID,vt,ft,bc); }	//!< Returns the texture coordinate on the given face at the given barycentric coordinates (bc).
	int     GetMaterialIndex(int faceID) const;				//!< Returns the material index of the face. This method goes through material counts of all materials to find the material index of the face. Returns a negaive number if the face as no material
	int     GetMaterialFaceCount(int mtlID) const { return mtlID>0 ? mcfc[mtlID]-mcfc[mtlID-1] : mcfc[0]; }	//!< Returns the number of faces associated with the given material ID.
	int     GetMaterialFirstFace(int mtlID) const { return mtlID>0 ? mcfc[mtlID-1] : 0; }	//!< Returns the first face index associated with the given material ID. Other faces associated with the same material are placed are placed consecutively.

	//!@name Compute Methods
	void ComputeBoundingBox();						//!< Computes the bounding box
	void ComputeNormals(bool clockwise=false);		//!< Computes and stores vertex normals

	//!@name Load and Save methods
	bool LoadFromFileObj( const char *filename, bool loadMtl=true );	//!< Loads the mesh from an OBJ file. Automatically converts all faces to triangles.
	bool SaveToFileObj( const char *filename );

private:
	template <class T> void Allocate(unsigned int n, T* &t) { if (t) delete [] t; if (n>0) t = new T[n]; else t=NULL; }
	template <class T> bool Allocate(unsigned int n, T* &t, unsigned int &nt) { if (n==nt) return false; nt=n; Allocate(n,t); return true; }
	static Point3f Interpolate( int i, const Point3f *v, const TriFace *f, const Point3f &bc ) { return v[f[i].v[0]]*bc.x + v[f[i].v[1]]*bc.y + v[f[i].v[2]]*bc.z; }
};

//-------------------------------------------------------------------------------

inline int TriMesh::GetMaterialIndex(int faceID) const
{
	for ( unsigned int i=0; i<nm; i++ ) {
		if ( faceID < mcfc[i] ) return (int) i;
	}
	return -1;
}

inline void TriMesh::ComputeBoundingBox()
{
	boundMin=v[0];
	boundMax=v[0];
	for ( unsigned int i=1; i<nv; i++ ) {
		if ( boundMin.x > v[i].x ) boundMin.x = v[i].x;
		if ( boundMin.y > v[i].y ) boundMin.y = v[i].y;
		if ( boundMin.z > v[i].z ) boundMin.z = v[i].z;
		if ( boundMax.x < v[i].x ) boundMax.x = v[i].x;
		if ( boundMax.y < v[i].y ) boundMax.y = v[i].y;
		if ( boundMax.z < v[i].z ) boundMax.z = v[i].z;
	}
}

inline void TriMesh::ComputeNormals(bool clockwise)
{
	SetNumNormals(nv);
	for ( unsigned int i=0; i<nvn; i++ ) vn[i].Set(0,0,0);	// initialize all normals to zero
	for ( unsigned int i=0; i<nf; i++ ) {
		Point3f N = (v[f[i].v[1]]-v[f[i].v[0]]) ^ (v[f[i].v[2]]-v[f[i].v[0]]);	// face normal (not normalized)
		if ( clockwise ) N = -N;
		vn[f[i].v[0]] += N;
		vn[f[i].v[1]] += N;
		vn[f[i].v[2]] += N;
		fn[i] = f[i];
	}
	for ( unsigned int i=0; i<nvn; i++ ) vn[i].Normalize();
}

struct MtlData
{
	char mtlName[256];
	unsigned int firstFace;
	unsigned int faceCount;
	MtlData() { mtlName[0]='\0'; faceCount=0; firstFace=0; }
};

struct MtlLibName { char filename[1024]; };

inline bool TriMesh::LoadFromFileObj( const char *filename, bool loadMtl )
{
	FILE *fp = fopen(filename,"r");
	if ( !fp ) return false;

	Clear();

	class Buffer
	{
		char data[1024];
		int readLine;
	public:
		int ReadLine(FILE *fp)
		{
			char c = fgetc(fp);
			while ( !feof(fp) ) {
				while ( isspace(c) && ( !feof(fp) || c!='\0' ) ) c = fgetc(fp);	// skip empty space
				if ( c == '#' ) while ( !feof(fp) && c!='\n' && c!='\r' && c!='\0' ) c = fgetc(fp);	// skip comment line
				else break;
			}
			int i=0;
			bool inspace = false;
			while ( i<1024-1 ) {
				if ( feof(fp) || c=='\n' || c=='\r' || c=='\0' ) break;
				if ( isspace(c) ) {	// only use a single space as the space character
					inspace = true;
				} else {
					if ( inspace ) data[i++] = ' ';
					inspace = false;
					data[i++] = c;
				}
				c = fgetc(fp);
			}
			data[i] = '\0';
			readLine = i;
			return i;
		}
		char& operator[](int i) { return data[i]; }
		void ReadVertex( Point3f &v ) const { sscanf( data+2, "%f %f %f", &v.x, &v.y, &v.z ); }
		void ReadFloat3( float f[3] ) const { sscanf( data+2, "%f %f %f", &f[0], &f[1], &f[2] ); }
		void ReadFloat( float *f ) const { sscanf( data+2, "%f", f ); }
		void ReadInt( int *i, int start ) const { sscanf( data+start, "%d", i ); }
		bool IsCommand( const char *cmd ) const {
			int i=0;
			while ( cmd[i]!='\0' ) {
				if ( cmd[i] != data[i] ) return false;
				i++;
			}
			return (data[i]=='\0' || data[i]==' ');
		}
		void Copy( char *a, int count, int start=0 ) const {
			strncpy( a, data+start, count-1 );
			a[count-1] = '\0';
		}
	};
	Buffer buffer;


	struct MtlList {
		std::vector<MtlData> mtlData;
		int GetMtlIndex(const char *mtlName)
		{
			for ( unsigned int i=0; i<mtlData.size(); i++ ) {
				if ( strcmp(mtlName,mtlData[i].mtlName) == 0 ) return (int)i;
			}
			return -1;
		}
		int CreateMtl(const char *mtlName, unsigned int firstFace)
		{
			if ( mtlName[0] == '\0' ) return 0;
			int i = GetMtlIndex(mtlName);
			if ( i >= 0 ) return i;
			MtlData m;
			strncpy(m.mtlName,mtlName,256);
			m.firstFace = firstFace;
			mtlData.push_back(m);
			return (int)mtlData.size()-1;
		}
	};
	MtlList mtlList;
	MtlData *currentMtlData = NULL;
	std::vector<Point3f>	_v;		// vertices
	std::vector<TriFace>	_f;		// faces
	std::vector<Point3f>	_vn;	// vertex normal
	std::vector<TriFace>	_fn;	// normal faces
	std::vector<Point3f>	_vt;	// texture vertices
	std::vector<TriFace>	_ft;	// texture faces
	std::vector<MtlLibName> mtlFiles;
	std::vector<int> faceMtlIndex;

	int currentMtlIndex = -1;
	bool hasTextures=false, hasNormals=false;

	while ( int rb = buffer.ReadLine(fp) ) {
		if ( buffer.IsCommand("v") ) {
			Point3f vertex;
			buffer.ReadVertex(vertex);
			_v.push_back(vertex);
		}
		else if ( buffer.IsCommand("vt") ) {
			Point3f texVert;
			buffer.ReadVertex(texVert);
			_vt.push_back(texVert);
			hasTextures = true;
		}
		else if ( buffer.IsCommand("vn") ) {
			Point3f normal;
			buffer.ReadVertex(normal);
			_vn.push_back(normal);
			hasNormals = true;
		}
		else if ( buffer.IsCommand("f") ) {
			int facevert = -1;
			bool inspace = true;
			int type = 0;
			unsigned int index;
			TriFace face, textureFace, normalFace;
			unsigned int nFacesBefore = _f.size();
			for ( int i=2; i<rb; i++ ) {
				if ( buffer[i] == ' ' ) inspace = true;
				else {
					if ( inspace ) {
						inspace=false;
						type=0;
						index=0;
						switch ( facevert ) {
							case -1:
								// initialize face
								face.v[0] = face.v[1] = face.v[2] = 0;
								textureFace.v[0] = textureFace.v[1] = textureFace.v[2] = 0;
								normalFace.v[0] = normalFace.v[1] = normalFace.v[2] = 0;
							case 0:
							case 1:
								facevert++;
								break;
							case 2:
								// copy the first two vertices from the previous face
								_f.push_back(face);
								face.v[1] = face.v[2];
								if ( hasTextures ) {
									_ft.push_back(textureFace);
									textureFace.v[1] = textureFace.v[2];
								}
								if ( hasNormals ) {
									_fn.push_back(normalFace);
									normalFace.v[1] = normalFace.v[2];
								}
								faceMtlIndex.push_back(currentMtlIndex);
								break;
						}
					}
					if ( buffer[i] == '/' ) { type++; index=0; }
					if ( buffer[i] >= '0' && buffer[i] <= '9' ) {
						index = index*10 + (buffer[i]-'0');
						switch ( type ) {
							case 0: face.v[facevert] = index-1; break;
							case 1: textureFace.v[facevert] = index-1; hasTextures=true; break;
							case 2: normalFace.v[facevert] = index-1; hasNormals=true; break;
						}
					}
				}
			}
			_f.push_back(face);
			if ( hasTextures ) _ft.push_back(textureFace);
			if ( hasNormals  ) _fn.push_back(normalFace);
			faceMtlIndex.push_back(currentMtlIndex);
			if ( currentMtlIndex>=0 ) mtlList.mtlData[currentMtlIndex].faceCount += _f.size() - nFacesBefore;
		}
		else if ( loadMtl ) {
			if ( buffer.IsCommand("usemtl") ) {
				char mtlName[256];
				buffer.Copy(mtlName,256,7);
				currentMtlIndex = mtlList.CreateMtl(mtlName, _f.size());
			}
			if ( buffer.IsCommand("mtllib") ) {
				MtlLibName libName;
				buffer.Copy(libName.filename,1024,7);
				mtlFiles.push_back(libName);
			}
		}
		if ( feof(fp) ) break;
	}

	fclose(fp);


	if ( _f.size() == 0 ) return true; // No faces found
	SetNumVertex(_v.size());
	SetNumFaces(_f.size());
	SetNumTexVerts(_vt.size());
	SetNumNormals(_vn.size());
	if ( loadMtl ) SetNumMtls(mtlList.mtlData.size());

	// Copy data
	memcpy(v, _v.data(), sizeof(Point3f)*_v.size());
	if ( _vt.size() > 0 ) memcpy(vt, _vt.data(), sizeof(Point3f)*_vt.size());
	if ( _vn.size() > 0 ) memcpy(vn, _vn.data(), sizeof(Point3f)*_vn.size());

	if ( mtlList.mtlData.size() > 0 ) {
		unsigned int fid = 0;
		for ( int m=0; m<(int)mtlList.mtlData.size(); m++ ) {
			for ( unsigned int i=mtlList.mtlData[m].firstFace, j=0; j<mtlList.mtlData[m].faceCount && i<_f.size(); i++ ) {
				if ( faceMtlIndex[i] == m ) {
					f[fid] = _f[i];
					if ( fn ) fn[fid] = _fn[i];
					if ( ft ) ft[fid] = _ft[i];
					fid++;
				}
			}
		}
		if ( fid <_f.size() ) {
			for ( unsigned int i=0; i<_f.size(); i++ ) {
				if ( faceMtlIndex[i] < 0 ) {
					f[fid] = _f[i];
					if ( fn ) fn[fid] = _fn[i];
					if ( ft ) ft[fid] = _ft[i];
					fid++;
				}
			}
		}
	} else {
		memcpy(f, _f.data(), sizeof(TriFace)*_f.size());
		if ( ft ) memcpy(ft, _ft.data(), sizeof(TriFace)*_ft.size());
		if ( fn ) memcpy(fn, _fn.data(), sizeof(TriFace)*_fn.size());
	}


	// Load the .mtl files
	if ( loadMtl ) {
		// get the path from filename
		char *mtlFullFilename = NULL;
		char *mtlFilename = NULL;
		const char* pathEnd = strrchr(filename,'\\');
		if ( !pathEnd ) pathEnd = strrchr(filename,'/');
		if ( pathEnd ) {
			int n = pathEnd-filename + 1;
			mtlFullFilename = new char[n+1024];
			strncpy(mtlFullFilename,filename,n);
			mtlFilename = &mtlFullFilename[n];
		} else {
			mtlFullFilename = new char[1024];
			mtlFilename = mtlFullFilename;
		}
		for ( unsigned int mi=0; mi<mtlFiles.size(); mi++ ) {
			strncpy( mtlFilename, mtlFiles[mi].filename, 1024 );
			FILE *fp = fopen(mtlFullFilename,"r");
			if ( !fp ) continue;
			int mtlID = -1;
			while ( int rb = buffer.ReadLine(fp) ) {
				if ( buffer.IsCommand("newmtl") ) {
					char mtlName[256];
					buffer.Copy(mtlName,256,7);
					mtlID = mtlList.GetMtlIndex(mtlName);
					if ( mtlID >= 0 ) strncpy(m[mtlID].name,mtlName,256);
				} else if ( mtlID >= 0 ) {
					if ( buffer.IsCommand("Ka") ) buffer.ReadFloat3( m[mtlID].Ka );
					else if ( buffer.IsCommand("Kd") ) buffer.ReadFloat3( m[mtlID].Kd );
					else if ( buffer.IsCommand("Ks") ) buffer.ReadFloat3( m[mtlID].Ks );
					else if ( buffer.IsCommand("Tf") ) buffer.ReadFloat3( m[mtlID].Tf );
					else if ( buffer.IsCommand("Ns") ) buffer.ReadFloat( &m[mtlID].Ns );
					else if ( buffer.IsCommand("Ni") ) buffer.ReadFloat( &m[mtlID].Ni );
					else if ( buffer.IsCommand("illum") ) buffer.ReadInt( &m[mtlID].illum, 5 );
					else if ( buffer.IsCommand("map_Ka") ) buffer.Copy( m[mtlID].map_Ka.name, 256, 7 );
					else if ( buffer.IsCommand("map_Kd") ) buffer.Copy( m[mtlID].map_Kd.name, 256, 7 );
					else if ( buffer.IsCommand("map_Ks") ) buffer.Copy( m[mtlID].map_Ks.name, 256, 7 );
				}
			}
			fclose(fp);
		}
		delete [] mtlFullFilename;
	}

	return true;
}

//-------------------------------------------------------------------------------

inline bool TriMesh::SaveToFileObj( const char *filename )
{
	FILE *fp = fopen(filename,"w");
	if ( !fp ) return false;

	for ( unsigned int i=0; i<nv; i++ ) {
		fprintf(fp,"v %f %f %f\n",v[i].x, v[i].y, v[i].z);
	}
	for ( unsigned int i=0; i<nvt; i++ ) {
		fprintf(fp,"vt %f %f %f\n",vt[i].x, vt[i].y, vt[i].z);
	}
	for ( unsigned int i=0; i<nvn; i++ ) {
		fprintf(fp,"vn %f %f %f\n",vn[i].x, vn[i].y, vn[i].z);
	}
	int faceFormat = ((nvn>0)<<1) | (nvt>0);
	switch ( faceFormat ) {
	case 0:
		for ( unsigned int i=0; i<nf; i++ ) {
			fprintf(fp,"f %d %d %d\n", f[i].v[0]+1, f[i].v[1]+1, f[i].v[2]+1);
		}
		break;
	case 1:
		for ( unsigned int i=0; i<nf; i++ ) {
			fprintf(fp,"f %d/%d %d/%d %d/%d\n", f[i].v[0]+1, ft[i].v[0]+1, f[i].v[1]+1, ft[i].v[1]+1, f[i].v[2]+1, ft[i].v[2]+1);
		}
		break;
	case 2:
		for ( unsigned int i=0; i<nf; i++ ) {
			fprintf(fp,"f %d//%d %d//%d %d//%d\n", f[i].v[0]+1, fn[i].v[0]+1, f[i].v[1]+1, fn[i].v[1]+1, f[i].v[2]+1, fn[i].v[2]+1);
		}
		break;
	case 3:
		for ( unsigned int i=0; i<nf; i++ ) {
			fprintf(fp,"f %d/%d/%d %d/%d/%d %d/%d/%d\n", f[i].v[0]+1, ft[i].v[0]+1, fn[i].v[0]+1, f[i].v[1]+1, ft[i].v[1]+1, fn[i].v[1]+1, f[i].v[2]+1, ft[i].v[2]+1, fn[i].v[2]+1);
		}
		break;
	}

	fclose(fp);

	return true;
}

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::TriMesh cyTriMesh;	//!< Triangular Mesh Class

//-------------------------------------------------------------------------------

#endif

