// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyGL.h 
//! \author Cem Yuksel
//! 
//! \brief  OpenGL helper classes
//! 
//! The classes in this file are designed to provide convenient interfaces for
//! some OpenGL features. They are not intended to provide the full flexibility
//! of the underlying OpenGL functions, but they greatly simplify the 
//! implementation of some general OpenGL operations.
//!
//-------------------------------------------------------------------------------
//
// Copyright (c) 2017, Cem Yuksel <cem@cemyuksel.com>
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

#ifndef _CY_GL_H_INCLUDED_
#define _CY_GL_H_INCLUDED_

//-------------------------------------------------------------------------------
#if !defined(__gl_h_) && !defined(__GL_H__) && !defined(_GL_H) && !defined(__X_GL_H)
#error gl.h not included before cyGL.h
#endif
#ifndef GL_VERSION_2_0
#error OpenGL 2.0 extensions are required for cyGL.h
#endif
//-------------------------------------------------------------------------------
#ifndef GL_GEOMETRY_SHADER
#define GL_GEOMETRY_SHADER 0x8DD9
#endif
#ifndef GL_TESS_EVALUATION_SHADER
#define GL_TESS_EVALUATION_SHADER 0x8E87
#endif
#ifndef GL_TESS_CONTROL_SHADER
#define GL_TESS_CONTROL_SHADER 0x8E88
#endif
//-------------------------------------------------------------------------------

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

#define CY_GL_INVALID_ID 0xFFFFFFFF

//-------------------------------------------------------------------------------

//! GLSL shader class.
//!
//! This class provides basic functionality for compiling GLSL shaders
//! either from a given source string or a given file.
//! It only stores the shader ID and it can be safely deleted after
//! the shader is used for building (linking) a GLSL program.

class GLSLShader
{
private:
	GLuint shaderID;	//!< The shader ID

public:
	GLSLShader() : shaderID(CY_GL_INVALID_ID) {}	//!< Constructor.
	virtual ~GLSLShader() { Delete(); }				//!< Destructor that deletes the shader.

	//!@name General Methods
	void   Delete();								//!< Deletes the shader.
	GLuint GetID () const { return shaderID; }		//!< Returns the shader ID.

	//!@name Compilation Methods

	//! Compiles the shader using the given file.
	//! If the shader was previously compiled, it is deleted.
	bool CompileFile( const char *filename, GLenum shaderType, std::ostream *outStream=&std::cout ) { return CompileFile(filename,shaderType,0,NULL,outStream); }

	//! Compiles the shader using the given file.
	//! If the shader was previously compiled, it is deleted.
	//! The prependSource string is added to the beginning of the shader code, so it must begin with the "#version" statement.
	bool CompileFile( const char *filename, GLenum shaderType, const char *prependSource, std::ostream *outStream=&std::cout ) { return CompileFile(filename,shaderType,1,&prependSource,outStream); }

	//! Compiles the shader using the given file.
	//! If the shader was previously compiled, it is deleted.
	//! The prependSources strings are added to the beginning of the shader code, so the first string must begin with "#version" statement.
	bool CompileFile( const char *filename, GLenum shaderType, int prependSourceCount, const char **prependSources, std::ostream *outStream=&std::cout );

	//! Compiles the shader using the given source code.
	//! If the shader was previously compiled, it is deleted.
	bool Compile( const char *shaderSourceCode, GLenum shaderType, std::ostream *outStream=&std::cout ) { return Compile(shaderSourceCode,shaderType,0,NULL,outStream); }

	//! Compiles the shader using the given source code.
	//! If the shader was previously compiled, it is deleted.
	//! The prependSource string is added to the beginning of the shader code, so it must begin with the "#version" statement.
	bool Compile( const char *shaderSourceCode, GLenum shaderType, const char *prependSource, std::ostream *outStream=&std::cout ) { return Compile(shaderSourceCode,shaderType,1,&prependSource,outStream); }

	//! Compiles the shader using the given source code.
	//! If the shader was previously compiled, it is deleted.
	//! The prependSources strings are added to the beginning of the shader code, so the first string must begin with "#version" statement.
	bool Compile( const char *shaderSourceCode, GLenum shaderType, int prependSourceCount, const char **prependSources, std::ostream *outStream=&std::cout );
};

//-------------------------------------------------------------------------------

//! GLSL program class.
//!
//! This class provides basic functionality for building GLSL programs
//! using vertex and fragment shaders, along with optionally geometry and tessellation shaders.
//! The shader sources can be provides as GLSLShader class objects, source strings, or file names.
//! This class also stores a vector of registered uniform parameter IDs.

class GLSLProgram
{
private:
	GLuint programID;			//!< The program ID
	std::vector<GLint> params;	//!< A list of registered uniform parameter IDs

public:
	GLSLProgram() : programID(CY_GL_INVALID_ID) {}		//!< Constructor
	virtual ~GLSLProgram() { Delete(); }				//!< Destructor that deletes the program

	//!@name General Methods
	void   Delete();									//!< Deletes the program.
	GLuint GetID () const { return programID; }			//!< Returns the program ID
	void   Bind  () const { glUseProgram(programID); }	//!< Binds the program for rendering

	//! Attaches the given shader to the program.
	//! This function must be called before calling Link.
	void CreateProgram() { Delete(); programID = glCreateProgram(); }

	//! Attaches the given shader to the program.
	//! This function must be called before calling Link.
	void AttachShader( const GLSLShader &shader ) { AttachShader(shader.GetID()); }

	//! Attaches the given shader to the program.
	//! This function must be called before calling Link.
	void AttachShader( GLuint shaderID ) { glAttachShader(programID,shaderID); }

	//! Links the program.
	//! The shaders must be attached before calling this function.
	//! Returns true if the link operation is successful.
	//! Writes any error or warning messages to the given output stream.
	bool Link( std::ostream *outStream=&std::cout );

	//!@name Build Methods

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	bool Build( const GLSLShader *vertexShader, 
                const GLSLShader *fragmentShader,
	            const GLSLShader *geometryShader=NULL,
	            const GLSLShader *tessControlShader=NULL,
	            const GLSLShader *tessEvaluationShader=NULL,
	            std::ostream *outStream=&std::cout );

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	bool BuildFiles( const char *vertexShaderFile, 
                     const char *fragmentShaderFile,
	                 const char *geometryShaderFile=NULL,
	                 const char *tessControlShaderFile=NULL,
	                 const char *tessEvaluationShaderFile=NULL,
	                 std::ostream *outStream=&std::cout )
	{ return BuildFiles(vertexShaderFile,fragmentShaderFile,geometryShaderFile,tessControlShaderFile,tessEvaluationShaderFile,0,NULL,outStream); }

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	//! The prependSource string is added to the beginning of each shader code, so it must begin with the "#version" statement.
	bool BuildFiles( const char *vertexShaderFile, 
                     const char *fragmentShaderFile,
	                 const char *geometryShaderFile,
	                 const char *tessControlShaderFile,
	                 const char *tessEvaluationShaderFile,
	                 const char *prependSource,
	                 std::ostream *outStream=&std::cout )
	{ return BuildFiles(vertexShaderFile,fragmentShaderFile,geometryShaderFile,tessControlShaderFile,tessEvaluationShaderFile,1,&prependSource,outStream); }

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	//! The prependSources strings are added to the beginning of each shader code, so the first string must begin with "#version" statement.
	bool BuildFiles( const char *vertexShaderFile, 
                     const char *fragmentShaderFile,
	                 const char *geometryShaderFile,
	                 const char *tessControlShaderFile,
	                 const char *tessEvaluationShaderFile,
	                 int         prependSourceCount,
	                 const char **prependSource,
	                 std::ostream *outStream=&std::cout );

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	bool BuildSources( const char *vertexShaderSourceCode, 
                       const char *fragmentShaderSourceCode,
	                   const char *geometryShaderSourceCode=NULL,
	                   const char *tessControlShaderSourceCode=NULL,
	                   const char *tessEvaluationShaderSourceCode=NULL,
	                   std::ostream *outStream=&std::cout )
	{ return BuildSources(vertexShaderSourceCode,fragmentShaderSourceCode,geometryShaderSourceCode,tessControlShaderSourceCode,tessEvaluationShaderSourceCode,0,NULL,outStream); }

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	//! The prependSource string is added to the beginning of each shader code, so it must begin with the "#version" statement.
	bool BuildSources( const char *vertexShaderSourceCode, 
                       const char *fragmentShaderSourceCode,
	                   const char *geometryShaderSourceCode,
	                   const char *tessControlShaderSourceCode,
	                   const char *tessEvaluationShaderSourceCode,
	                   const char *prependSource,
	                   std::ostream *outStream=&std::cout )
	{ return BuildSources(vertexShaderSourceCode,fragmentShaderSourceCode,geometryShaderSourceCode,tessControlShaderSourceCode,tessEvaluationShaderSourceCode,1,&prependSource,outStream); }

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	//! The prependSources strings are added to the beginning of each shader code, so the first string must begin with "#version" statement.
	bool BuildSources( const char *vertexShaderSourceCode, 
                       const char *fragmentShaderSourceCode,
	                   const char *geometryShaderSourceCode,
	                   const char *tessControlShaderSourceCode,
	                   const char *tessEvaluationShaderSourceCode,
	                   int         prependSourceCount,
	                   const char **prependSource,
	                   std::ostream *outStream=&std::cout );

	//!@name Uniform Parameter Methods

	//! Registers a single uniform parameter.
	//! The index must be unique and the name should match a uniform parameter name in one of the shaders.
	//! The index values for different parameters don't have to be consecutive, but unused index values waste memory.
	void RegisterUniform( unsigned int index, const char *name, std::ostream *outStream=&std::cout );

	//! Registers multiple parameters.
	//! The names should be separated by a space character.
	void RegisterUniforms( const char *names, unsigned int startingIndex=0, std::ostream *outStream=&std::cout );

	void SetUniform(int param, float x)                                { glUniform1f  (params[param],x); }				//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, float x, float y)                       { glUniform2f  (params[param],x,y); }			//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, float x, float y, float z)              { glUniform3f  (params[param],x,y,z); }			//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, float x, float y, float z, float w)     { glUniform4f  (params[param],x,y,z,w); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, int count, const float  *data)          { glUniform2fv (params[param],count,data); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, int x)                                  { glUniform1i  (params[param],x); }				//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, int x, int y)                           { glUniform2i  (params[param],x,y); }			//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, int x, int y, int z)                    { glUniform3i  (params[param],x,y,z); }			//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, int x, int y, int z, int w)             { glUniform4i  (params[param],x,y,z,w); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, int count, const int    *data)          { glUniform2iv (params[param],count,data); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
#ifdef GL_VERSION_3_0
	void SetUniform(int param, GLuint x)                               { glUniform1ui (params[param],x); }				//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, GLuint x, GLuint y)                     { glUniform2ui (params[param],x,y); }			//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, GLuint x, GLuint y, GLuint z)           { glUniform3ui (params[param],x,y,z); }			//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, GLuint x, GLuint y, GLuint z, GLuint w) { glUniform4ui (params[param],x,y,z,w); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, int count, const GLuint *data)          { glUniform2uiv(params[param],count,data); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
#endif
#ifdef GL_VERSION_4_0
	void SetUniform(int param, double x)                               { glUniform1d  (params[param],x); }				//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, double x, double y)                     { glUniform2d  (params[param],x,y); }			//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, double x, double y, double z)           { glUniform3d  (params[param],x,y,z); }			//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, double x, double y, double z, double w) { glUniform4d  (params[param],x,y,z,w); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, int count, const double *data)          { glUniform2dv (params[param],count,data); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
#endif

	void SetUniformMatrix2  (int param, const float  *m, bool transpose=false) { glUniformMatrix2fv  (params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix3  (int param, const float  *m, bool transpose=false) { glUniformMatrix3fv  (params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix4  (int param, const float  *m, bool transpose=false) { glUniformMatrix4fv  (params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
#ifdef GL_VERSION_2_1
	void SetUniformMatrix2x3(int param, const float  *m, bool transpose=false) { glUniformMatrix2x3fv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix2x4(int param, const float  *m, bool transpose=false) { glUniformMatrix2x4fv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix3x2(int param, const float  *m, bool transpose=false) { glUniformMatrix3x2fv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix3x4(int param, const float  *m, bool transpose=false) { glUniformMatrix3x4fv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix4x2(int param, const float  *m, bool transpose=false) { glUniformMatrix4x2fv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix4x3(int param, const float  *m, bool transpose=false) { glUniformMatrix4x3fv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
#endif
#ifdef GL_VERSION_4_0
	void SetUniformMatrix2  (int param, const double *m, bool transpose=false) { glUniformMatrix2dv  (params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix3  (int param, const double *m, bool transpose=false) { glUniformMatrix3dv  (params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix4  (int param, const double *m, bool transpose=false) { glUniformMatrix4dv  (params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix2x3(int param, const double *m, bool transpose=false) { glUniformMatrix2x3dv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix2x4(int param, const double *m, bool transpose=false) { glUniformMatrix2x4dv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix3x2(int param, const double *m, bool transpose=false) { glUniformMatrix3x2dv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix3x4(int param, const double *m, bool transpose=false) { glUniformMatrix3x4dv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix4x2(int param, const double *m, bool transpose=false) { glUniformMatrix4x2dv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniformMatrix4x3(int param, const double *m, bool transpose=false) { glUniformMatrix4x3dv(params[param],1,transpose,m); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
#endif


#ifdef _CY_POINT_H_INCLUDED_
	void SetUniform(int param, const Point2<float>   &p) { glUniform2fv (params[param],2,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Point3<float>   &p) { glUniform3fv (params[param],3,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Point4<float>   &p) { glUniform4fv (params[param],4,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	template <int N> 
	void SetUniform(int param, const Point<float, N> &p) { glUniform4fv (params[param],N,p.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Point2<int>     &p) { glUniform2iv (params[param],2,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Point3<int>     &p) { glUniform3iv (params[param],3,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Point4<int>     &p) { glUniform4iv (params[param],4,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	template <int N> 
	void SetUniform(int param, const Point<int,   N> &p) { glUniform4iv (params[param],N,p.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
# ifdef GL_VERSION_3_0
	void SetUniform(int param, const Point2<GLuint>  &p) { glUniform2uiv(params[param],2,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Point3<GLuint>  &p) { glUniform3uiv(params[param],3,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Point4<GLuint>  &p) { glUniform4uiv(params[param],4,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	template <int N> 
	void SetUniform(int param, const Point<GLuint,N> &p) { glUniform4uiv(params[param],N,p.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
# endif
# ifdef GL_VERSION_4_0
	void SetUniform(int param, const Point2<double>  &p) { glUniform2dv (params[param],2,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Point3<double>  &p) { glUniform3dv (params[param],3,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Point4<double>  &p) { glUniform4dv (params[param],4,&p.x); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	template <int N> 
	void SetUniform(int param, const Point<double,N> &p) { glUniform4dv (params[param],N,p.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
# endif
#endif

#ifdef _CY_IPOINT_H_INCLUDED_
	void SetUniform(int param, const IPoint2<int>     &p) { glUniform2iv (params[param],2,&p.x); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const IPoint3<int>     &p) { glUniform3iv (params[param],3,&p.x); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	template <int N> 
	void SetUniform(int param, const IPoint<int,   N> &p) { glUniform4iv (params[param],N,p.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
# ifdef GL_VERSION_3_0
	void SetUniform(int param, const IPoint2<GLuint>  &p) { glUniform2uiv(params[param],2,&p.x); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const IPoint3<GLuint>  &p) { glUniform3uiv(params[param],3,&p.x); }		//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	template <int N> 
	void SetUniform(int param, const IPoint<GLuint,N> &p) { glUniform4uiv(params[param],N,p.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
# endif
#endif

#ifdef _CY_MATRIX_H_INCLUDED_
	void SetUniform(int param, const Matrix2 <float>  &m) { glUniformMatrix2fv  (params[param],1,GL_FALSE,m.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Matrix3 <float>  &m) { glUniformMatrix3fv  (params[param],1,GL_FALSE,m.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Matrix4 <float>  &m) { glUniformMatrix4fv  (params[param],1,GL_FALSE,m.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
# ifdef GL_VERSION_2_1
	void SetUniform(int param, const Matrix34<float>  &m) { glUniformMatrix3x4fv(params[param],1,GL_FALSE,m.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
# endif
# ifdef GL_VERSION_4_0
	void SetUniform(int param, const Matrix2 <double> &m) { glUniformMatrix2dv  (params[param],1,GL_FALSE,m.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Matrix3 <double> &m) { glUniformMatrix3dv  (params[param],1,GL_FALSE,m.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Matrix4 <double> &m) { glUniformMatrix4dv  (params[param],1,GL_FALSE,m.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
	void SetUniform(int param, const Matrix34<double> &m) { glUniformMatrix3x4dv(params[param],1,GL_FALSE,m.data); }	//!< Sets the value of the uniform parameter with the given index. The parameter must be registered before.
# endif
#endif

};

//-------------------------------------------------------------------------------
// GLSLShader Implementation
//-------------------------------------------------------------------------------

inline void GLSLShader::Delete()
{
	if ( shaderID != CY_GL_INVALID_ID ) {
		try { 
			// If the OpenGL context is deleted before calling this function,
			// it can cause crash.
			glDeleteShader(shaderID);
		} catch (...) {}
		shaderID = CY_GL_INVALID_ID;
	}
}

inline bool GLSLShader::CompileFile( const char *filename, GLenum shaderType, int prependSourceCount, const char **prependSources, std::ostream *outStream )
{
	std::ifstream shaderStream(filename, std::ios::in);
	if(! shaderStream.is_open()) {
		if ( outStream ) *outStream << "ERROR: Cannot open file." << std::endl;
		return false;
	}

	std::string shaderSourceCode((std::istreambuf_iterator<char>(shaderStream)), std::istreambuf_iterator<char>());
	shaderStream.close();

	return Compile( shaderSourceCode.data(), shaderType, prependSourceCount, prependSources, outStream );
}

inline bool GLSLShader::Compile( const char *shaderSourceCode, GLenum shaderType, int prependSourceCount, const char **prependSources, std::ostream *outStream )
{
	Delete();

	shaderID = glCreateShader( shaderType );
	if ( prependSourceCount > 0 ) {
		std::vector<const char*> sources(prependSourceCount+1);
		for ( int i=0; i<prependSourceCount; i++ ) sources[i] = prependSources[i];
		sources[prependSourceCount] = shaderSourceCode;
		glShaderSource(shaderID, prependSourceCount+1, sources.data(), NULL);
	} else {
		glShaderSource(shaderID, 1, &shaderSourceCode, NULL);
	}
	glCompileShader(shaderID);

	GLint result = GL_FALSE;
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &result);

	int infoLogLength;
	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
	if ( infoLogLength > 1 ) {
		std::vector<char> compilerMessage(infoLogLength);
		glGetShaderInfoLog( shaderID, infoLogLength, NULL, compilerMessage.data() );
		if ( outStream ) {
			if ( !result ) *outStream << "ERROR: Cannot compile shader." << std::endl;
			int versionMajor, versionMinor;
#ifdef GL_VERSION_3_0
			glGetIntegerv(GL_MAJOR_VERSION,&versionMajor);
			glGetIntegerv(GL_MAJOR_VERSION,&versionMinor);
			*outStream << "OpenGL version " << versionMajor << "." << versionMinor << std::endl;
#endif
			*outStream << compilerMessage.data() << std::endl;
		}
	}

	if ( result ) {
		GLint stype;
		glGetShaderiv(shaderID, GL_SHADER_TYPE, &stype);
		if ( stype != shaderType ) {
			if ( outStream ) *outStream << "ERROR: Incorrect shader type." << std::endl;
			return false;
		}
	}

	return result == GL_TRUE;
}

//-------------------------------------------------------------------------------
// GLSLProgram Implementation
//-------------------------------------------------------------------------------

inline void GLSLProgram::Delete()
{
	if ( programID != CY_GL_INVALID_ID ) {
		try { 
			// If the OpenGL context is deleted before calling this function,
			// it can cause crash.
			glDeleteProgram(programID);
		} catch (...) {}
		programID = CY_GL_INVALID_ID;
	}
}

inline bool GLSLProgram::Link( std::ostream *outStream )
{
	glLinkProgram(programID);

	GLint result = GL_FALSE;
	glGetProgramiv(programID, GL_LINK_STATUS, &result);

	int infoLogLength;
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &infoLogLength);
	if ( infoLogLength > 1 ) {
		std::vector<char> compilerMessage(infoLogLength);
		glGetProgramInfoLog( programID, infoLogLength, NULL, compilerMessage.data() );
		if ( outStream ) *outStream << "ERROR: " << compilerMessage.data() << std::endl;
	}

	return result == GL_TRUE;
}

inline bool GLSLProgram::Build( const GLSLShader *vertexShader, 
                                const GLSLShader *fragmentShader,
	                            const GLSLShader *geometryShader,
	                            const GLSLShader *tessControlShader,
	                            const GLSLShader *tessEvaluationShader,
	                            std::ostream *outStream )
{
	CreateProgram();
	std::stringstream output;
	AttachShader(*vertexShader);
	AttachShader(*fragmentShader);
	if ( geometryShader ) AttachShader(*geometryShader);
	if ( tessControlShader ) AttachShader(*tessControlShader);
	if ( tessEvaluationShader ) AttachShader(*tessEvaluationShader);
	return Link(outStream);
}

inline bool GLSLProgram::BuildFiles( const char *vertexShaderFile, 
                                     const char *fragmentShaderFile,
	                                 const char *geometryShaderFile,
	                                 const char *tessControlShaderFile,
	                                 const char *tessEvaluationShaderFile,
	                                 int         prependSourceCount,
	                                 const char **prependSource,
	                                 std::ostream *outStream )
{
	CreateProgram();
	GLSLShader vs, fs, gs, tcs, tes;
	std::stringstream shaderOutput;
	if ( ! vs.CompileFile(vertexShaderFile, GL_VERTEX_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
		if ( outStream ) *outStream << "ERROR: Failed compiling vertex shader \"" << vertexShaderFile << ".\"" << std::endl << shaderOutput.str();
		return false;
	}
	AttachShader(vs);
	if ( ! fs.CompileFile(fragmentShaderFile, GL_FRAGMENT_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
		if ( outStream ) *outStream << "ERROR: Failed compiling fragment shader \"" << fragmentShaderFile << ".\"" <<  std::endl << shaderOutput.str();
		return false;
	}
	AttachShader(fs);
	if ( geometryShaderFile ) {
		if ( ! gs.CompileFile(geometryShaderFile, GL_GEOMETRY_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling geometry shader \"" << geometryShaderFile << ".\"" <<  std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(gs);
	}
	if ( tessControlShaderFile ) {
		if ( ! tcs.CompileFile(tessControlShaderFile, GL_TESS_CONTROL_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling tessellation control shader \"" << tessControlShaderFile << ".\"" <<  std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(tcs);
	}
	if ( tessEvaluationShaderFile ) {
		if ( ! tes.CompileFile(tessEvaluationShaderFile, GL_TESS_EVALUATION_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling tessellation evaluation shader \"" << tessEvaluationShaderFile << ".\"" <<  std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(tes);
	}
	return Link(outStream);
}


inline bool GLSLProgram::BuildSources( const char *vertexShaderSourceCode, 
                                       const char *fragmentShaderSourceCode,
	                                   const char *geometryShaderSourceCode,
	                                   const char *tessControlShaderSourceCode,
	                                   const char *tessEvaluationShaderSourceCode,
	                                   int         prependSourceCount,
	                                   const char **prependSource,
	                                   std::ostream *outStream )
{
	CreateProgram();
	GLSLShader vs, fs, gs, tcs, tes;
	std::stringstream shaderOutput;
	if ( ! vs.Compile(vertexShaderSourceCode, GL_VERTEX_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
		if ( outStream ) *outStream << "ERROR: Failed compiling vertex shader." << std::endl << shaderOutput.str();
		return false;
	}
	AttachShader(vs);
	if ( ! fs.Compile(fragmentShaderSourceCode, GL_FRAGMENT_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
		if ( outStream ) *outStream << "ERROR: Failed compiling fragment shader." << std::endl << shaderOutput.str();
		return false;
	}
	AttachShader(fs);
	if ( geometryShaderSourceCode ) {
		if ( ! gs.Compile(geometryShaderSourceCode, GL_GEOMETRY_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling geometry shader." << std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(gs);
	}
	if ( tessControlShaderSourceCode ) {
		if ( ! tcs.Compile(tessControlShaderSourceCode, GL_TESS_CONTROL_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling tessellation control shader." << std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(tcs);
	}
	if ( tessEvaluationShaderSourceCode ) {
		if ( ! tes.Compile(tessEvaluationShaderSourceCode, GL_TESS_EVALUATION_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling tessellation evaluation shader." << std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(tes);
	}
	return Link(outStream);
}

inline void GLSLProgram::RegisterUniform( unsigned int index, const char *name, std::ostream *outStream )
{
	if ( params.size() <= index ) params.resize( index+1, -1 );
	params[index] = glGetUniformLocation( programID, name );
	if ( params[index] < 0 ) {
		GLenum error = glGetError();
		GLenum newError;
		while ( (newError = glGetError()) != GL_NO_ERROR ) error = newError; // get the latest error.
		if ( outStream ) {
			*outStream << "ERROR: ";
			switch (error) {
				case GL_NO_ERROR:					   *outStream << "GL_NO_ERROR.";					  break;
				case GL_INVALID_ENUM:				   *outStream << "GL_INVALID_ENUM.";				  break;
				case GL_INVALID_VALUE:				   *outStream << "GL_INVALID_VALUE.";				  break;
				case GL_INVALID_OPERATION:			   *outStream << "GL_INVALID_OPERATION.";			  break;
				case GL_STACK_OVERFLOW:				   *outStream << "GL_STACK_OVERFLOW.";				  break;
				case GL_STACK_UNDERFLOW:			   *outStream << "GL_STACK_UNDERFLOW.";				  break;
				case GL_OUT_OF_MEMORY:				   *outStream << "GL_OUT_OF_MEMORY.";				  break;
#ifdef GL_VERSION_3_0
				case GL_INVALID_FRAMEBUFFER_OPERATION: *outStream << "GL_INVALID_FRAMEBUFFER_OPERATION."; break;
#endif
			}
			*outStream << " Parameter \"" << name << "\" could not be registered." << std::endl;
		}
	}
}

inline void GLSLProgram::RegisterUniforms( const char *names, unsigned int startingIndex, std::ostream *outStream )
{
	std::stringstream ss(names);
	unsigned int index = startingIndex;
	while ( ss.good() ) {
		std::string name;
		ss >> name;
		RegisterUniform( index++, name.c_str(), outStream );
	}
}

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::GLSLShader  cyGLSLShader;		//!< GLSL shader class
typedef cy::GLSLProgram cyGLSLProgram;		//!< GLSL program class

//-------------------------------------------------------------------------------
#endif