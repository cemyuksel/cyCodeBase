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
#error OpenGL 2.0 extensions are required for cyGL.h. You must include an OpenGL extensions header before including cyGL.h.
#endif
#ifndef GL_VERSION_3_0
# define _CY_GL_VERSION_3_0_WARNING "OpenGL version 3 is required for using geometry and tessellation shaders, but the OpenGL extensions header included before cyGL.h does not include OpenGL version 3.0 definitions."
# if _MSC_VER
#  pragma message ("Warning: " _CY_GL_VERSION_3_0_WARNING)
# elif   __GNUC__
#  warning (_CY_GL_VERSION_3_0_WARNING)
# endif
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
#ifdef APIENTRY
# define _CY_APIENTRY APIENTRY
#else
# if defined(__MINGW32__) || defined(__CYGWIN__) || (_MSC_VER >= 800) || defined(_STDCALL_SUPPORTED) || defined(__BORLANDC__)
#  define _CY_APIENTRY __stdcall
# else
#  define _CY_APIENTRY
# endif
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

#define CY_GL_INVALID_ID 0xFFFFFFFF	//!< Invalid ID value

//-------------------------------------------------------------------------------

//! This global function prints the OpenGL version to the given stream.
inline void GLPrintVersion(std::ostream *outStream=&std::cout)
{
	const GLubyte* version = glGetString(GL_VERSION);
	if ( version ) *outStream << version;
	else {
		int versionMajor=0, versionMinor=0;
#ifdef GL_VERSION_3_0
		glGetIntegerv(GL_MAJOR_VERSION,&versionMajor);
		glGetIntegerv(GL_MINOR_VERSION,&versionMinor);
#endif
		if ( versionMajor > 0 && versionMajor < 100 ) {
			*outStream << versionMajor << "." << versionMinor;
		} else *outStream << "Unknown";
	}
}

//-------------------------------------------------------------------------------

//! Checks all previously triggered OpenGL errors and prints them to the given output stream.
static void GLCheckError(const char *sourcefile, int line, const char *call=nullptr, std::ostream *outStream=&std::cout)
{
	GLenum error;
    while ( (error = glGetError()) != GL_NO_ERROR) {
        *outStream << "OpenGL ERROR: " << sourcefile << " (line " << line << "): ";
		if ( call ) *outStream << call << " triggered ";
		*outStream << gluErrorString(error) << std::endl;
    }
}

//! Checks and prints OpenGL error messages to the default output stream.
#define CY_GL_ERROR _CY_GL_ERROR
#define _CY_GL_ERROR cy::GLCheckError(__FILE__,__LINE__)

//! Checks OpenGL errors before calling the given function,
//! calls the function, and then checks OpenGL errors again.
//! If an error is found, it is printed to the default output stream.
#define CY_GL_ERR(gl_function_call) _CY_GL_ERR(gl_function_call)
#define _CY_GL_ERR(f) cy::GLCheckError(__FILE__,__LINE__,"a prior call"); f; cy::GLCheckError(__FILE__,__LINE__,#f)

#ifdef _DEBUG
# define CY_GL_ERROR_D CY_GL_ERROR //!< Checks and prints OpenGL error messages using CY_GL_ERROR in code compiled in debug mode only (with _DEBUG defined).
# define CY_GL_ERR_D   CY_GL_ERR   //!< Checks and prints OpenGL error messages using CY_GL_ERR in code compiled in debug mode only (with _DEBUG defined).
#else
# define CY_GL_ERROR_D			   //!< Checks and prints OpenGL error messages using CY_GL_ERROR in code compiled in debug mode only (with _DEBUG defined).
# define CY_GL_ERR_D			   //!< Checks and prints OpenGL error messages using CY_GL_ERR in code compiled in debug mode only (with _DEBUG defined).
#endif

//-------------------------------------------------------------------------------

#ifdef GL_KHR_debug

//! OpenGL debug callback class.
//!
//! For this class to work, you may need to initialize the OpenGL context in debug mode.
//! This class registers an OpenGL debug callback function, which is called when
//! there is an OpenGL generates a debug message. 
//! The class has no local storage, so it can be safely deleted.
//! Deleting an object of this class, however, does not automatically disable the debug callbacks.
class GLDebugCallback
{
public:
	//! Constructor can register the callback, but only if the OpenGL context is created
	//! before the constructor is called. If the regsiterCallback argument is false,
	//! the other arguments are ignored.
	GLDebugCallback(bool registerCallback=false, bool ignoreNotifications=false, std::ostream *outStream=&std::cout)
		{ if ( registerCallback ) { Register(outStream); IgnoreNotifications(ignoreNotifications); } }

	//! Registers the debug callback function.
	//! The callback function outputs the debug data to the given stream.
	//! Note that if the OpenGL context is not created in debug mode, 
	//! OpenGL may not call the callback function.
	//! If there is a previously registered callback function,
	//! calling this function overwrites the previous callback registration.
	void Register(std::ostream *outStream=&std::cout) { glEnable(GL_DEBUG_OUTPUT); glDebugMessageCallback(Callback,outStream); }

	//! Unregisters the OpenGL debug callback function.
	void Unregister() { glDisable(GL_DEBUG_OUTPUT); glDebugMessageCallback(0,0); }

	//! Sets which type of non-critical debug messages should be ignored.
	//! By default, no debug message type is ignored.
	void SetIgnoredTypes( bool deprecated_behavior, bool portability, bool performance, bool other )
	{
		glDebugMessageControl( GL_DONT_CARE, GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR, GL_DONT_CARE, 0, 0, deprecated_behavior );
		glDebugMessageControl( GL_DONT_CARE, GL_DEBUG_TYPE_PORTABILITY,         GL_DONT_CARE, 0, 0, portability );
		glDebugMessageControl( GL_DONT_CARE, GL_DEBUG_TYPE_PERFORMANCE,         GL_DONT_CARE, 0, 0, performance );
		glDebugMessageControl( GL_DONT_CARE, GL_DEBUG_TYPE_OTHER,               GL_DONT_CARE, 0, 0, other );
	}

	//! Sets weather notification messages should be ignored.
	void IgnoreNotifications(bool ignore=true) { glDebugMessageControl( GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_NOTIFICATION, 0, 0, !ignore ); }

protected:
	//! This callback function that is called by OpenGL whenever there is a debug message.
	//! See the OpenGL documentation for glDebugMessageCallback for details.
	//! Placing the break point in this function allows easily identifying the
	//! OpenGL call that triggered the debug message (using the call stack).
	static void _CY_APIENTRY Callback( GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam );
	static void _CY_APIENTRY Callback( GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, void* userParam ) { Callback(source,type,id,severity,length,message,(const void*)userParam); }
};

//! Registers the OpenGL callback by ignoring notifications.
//! After this macro is called, the debug messages get printed to the default output stream.
//! The OpenGL context must be created before using this macro.
//! Note that OpenGL may not generate debug messages if it is not created in debug mode.
#define CY_GL_REGISTER_DEBUG_CALLBACK _CY_GL_REGISTER_DEBUG_CALLBACK
#define _CY_GL_REGISTER_DEBUG_CALLBACK { cy::GLDebugCallback callback(true,true); }

#endif // GL_KHR_debug

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
	bool CompileFile( const char *filename, GLenum shaderType, std::ostream *outStream=&std::cout ) { return CompileFile(filename,shaderType,0,nullptr,outStream); }

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
	bool Compile( const char *shaderSourceCode, GLenum shaderType, std::ostream *outStream=&std::cout ) { return Compile(shaderSourceCode,shaderType,0,nullptr,outStream); }

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
	            const GLSLShader *geometryShader=nullptr,
	            const GLSLShader *tessControlShader=nullptr,
	            const GLSLShader *tessEvaluationShader=nullptr,
	            std::ostream *outStream=&std::cout );

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	bool BuildFiles( const char *vertexShaderFile, 
                     const char *fragmentShaderFile,
	                 const char *geometryShaderFile=nullptr,
	                 const char *tessControlShaderFile=nullptr,
	                 const char *tessEvaluationShaderFile=nullptr,
	                 std::ostream *outStream=&std::cout )
	{ return BuildFiles(vertexShaderFile,fragmentShaderFile,geometryShaderFile,tessControlShaderFile,tessEvaluationShaderFile,0,nullptr,outStream); }

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
	                   const char *geometryShaderSourceCode=nullptr,
	                   const char *tessControlShaderSourceCode=nullptr,
	                   const char *tessEvaluationShaderSourceCode=nullptr,
	                   std::ostream *outStream=&std::cout )
	{ return BuildSources(vertexShaderSourceCode,fragmentShaderSourceCode,geometryShaderSourceCode,tessControlShaderSourceCode,tessEvaluationShaderSourceCode,0,nullptr,outStream); }

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

	//!@{
	//! Sets the value of the uniform parameter with the given index. 
	//! The uniform parameter must be registered before using RegisterUniform() or RegisterUniforms().
	//! The program must be bind by calling Bind() before this method.
	void SetUniform (int index, float x)                                { glUniform1f  (params[index],x); }
	void SetUniform (int index, float x, float y)                       { glUniform2f  (params[index],x,y); }
	void SetUniform (int index, float x, float y, float z)              { glUniform3f  (params[index],x,y,z); }
	void SetUniform (int index, float x, float y, float z, float w)     { glUniform4f  (params[index],x,y,z,w); }
	void SetUniform1(int index, int count, const float  *data)          { glUniform1fv (params[index],count,data); }
	void SetUniform2(int index, int count, const float  *data)          { glUniform2fv (params[index],count,data); }
	void SetUniform3(int index, int count, const float  *data)          { glUniform3fv (params[index],count,data); }
	void SetUniform4(int index, int count, const float  *data)          { glUniform4fv (params[index],count,data); }
	void SetUniform (int index, int x)                                  { glUniform1i  (params[index],x); }
	void SetUniform (int index, int x, int y)                           { glUniform2i  (params[index],x,y); }
	void SetUniform (int index, int x, int y, int z)                    { glUniform3i  (params[index],x,y,z); }
	void SetUniform (int index, int x, int y, int z, int w)             { glUniform4i  (params[index],x,y,z,w); }
	void SetUniform1(int index, int count, const int    *data)          { glUniform1iv (params[index],count,data); }
	void SetUniform2(int index, int count, const int    *data)          { glUniform2iv (params[index],count,data); }
	void SetUniform3(int index, int count, const int    *data)          { glUniform3iv (params[index],count,data); }
	void SetUniform4(int index, int count, const int    *data)          { glUniform4iv (params[index],count,data); }
#ifdef GL_VERSION_3_0
	void SetUniform (int index, GLuint x)                               { glUniform1ui (params[index],x); }
	void SetUniform (int index, GLuint x, GLuint y)                     { glUniform2ui (params[index],x,y); }
	void SetUniform (int index, GLuint x, GLuint y, GLuint z)           { glUniform3ui (params[index],x,y,z); }
	void SetUniform (int index, GLuint x, GLuint y, GLuint z, GLuint w) { glUniform4ui (params[index],x,y,z,w); }
	void SetUniform1(int index, int count, const GLuint *data)          { glUniform1uiv(params[index],count,data); }
	void SetUniform2(int index, int count, const GLuint *data)          { glUniform2uiv(params[index],count,data); }
	void SetUniform3(int index, int count, const GLuint *data)          { glUniform3uiv(params[index],count,data); }
	void SetUniform4(int index, int count, const GLuint *data)          { glUniform4uiv(params[index],count,data); }
#endif
#ifdef GL_VERSION_4_0
	void SetUniform (int index, double x)                               { glUniform1d  (params[index],x); }
	void SetUniform (int index, double x, double y)                     { glUniform2d  (params[index],x,y); }
	void SetUniform (int index, double x, double y, double z)           { glUniform3d  (params[index],x,y,z); }
	void SetUniform (int index, double x, double y, double z, double w) { glUniform4d  (params[index],x,y,z,w); }
	void SetUniform1(int index, int count, const double *data)          { glUniform1dv (params[index],count,data); }
	void SetUniform2(int index, int count, const double *data)          { glUniform2dv (params[index],count,data); }
	void SetUniform3(int index, int count, const double *data)          { glUniform3dv (params[index],count,data); }
	void SetUniform4(int index, int count, const double *data)          { glUniform4dv (params[index],count,data); }
#endif

	void SetUniformMatrix2  (int index, const float  *m, int count=1, bool transpose=false) { glUniformMatrix2fv  (params[index],count,transpose,m); }
	void SetUniformMatrix3  (int index, const float  *m, int count=1, bool transpose=false) { glUniformMatrix3fv  (params[index],count,transpose,m); }
	void SetUniformMatrix4  (int index, const float  *m, int count=1, bool transpose=false) { glUniformMatrix4fv  (params[index],count,transpose,m); }
#ifdef GL_VERSION_2_1
	void SetUniformMatrix2x3(int index, const float  *m, int count=1, bool transpose=false) { glUniformMatrix2x3fv(params[index],count,transpose,m); }
	void SetUniformMatrix2x4(int index, const float  *m, int count=1, bool transpose=false) { glUniformMatrix2x4fv(params[index],count,transpose,m); }
	void SetUniformMatrix3x2(int index, const float  *m, int count=1, bool transpose=false) { glUniformMatrix3x2fv(params[index],count,transpose,m); }
	void SetUniformMatrix3x4(int index, const float  *m, int count=1, bool transpose=false) { glUniformMatrix3x4fv(params[index],count,transpose,m); }
	void SetUniformMatrix4x2(int index, const float  *m, int count=1, bool transpose=false) { glUniformMatrix4x2fv(params[index],count,transpose,m); }
	void SetUniformMatrix4x3(int index, const float  *m, int count=1, bool transpose=false) { glUniformMatrix4x3fv(params[index],count,transpose,m); }
#endif
#ifdef GL_VERSION_4_0
	void SetUniformMatrix2  (int index, const double *m, int count=1, bool transpose=false) { glUniformMatrix2dv  (params[index],count,transpose,m); }
	void SetUniformMatrix3  (int index, const double *m, int count=1, bool transpose=false) { glUniformMatrix3dv  (params[index],count,transpose,m); }
	void SetUniformMatrix4  (int index, const double *m, int count=1, bool transpose=false) { glUniformMatrix4dv  (params[index],count,transpose,m); }
	void SetUniformMatrix2x3(int index, const double *m, int count=1, bool transpose=false) { glUniformMatrix2x3dv(params[index],count,transpose,m); }
	void SetUniformMatrix2x4(int index, const double *m, int count=1, bool transpose=false) { glUniformMatrix2x4dv(params[index],count,transpose,m); }
	void SetUniformMatrix3x2(int index, const double *m, int count=1, bool transpose=false) { glUniformMatrix3x2dv(params[index],count,transpose,m); }	
	void SetUniformMatrix3x4(int index, const double *m, int count=1, bool transpose=false) { glUniformMatrix3x4dv(params[index],count,transpose,m); }	
	void SetUniformMatrix4x2(int index, const double *m, int count=1, bool transpose=false) { glUniformMatrix4x2dv(params[index],count,transpose,m); }	
	void SetUniformMatrix4x3(int index, const double *m, int count=1, bool transpose=false) { glUniformMatrix4x3dv(params[index],count,transpose,m); }	
#endif

#ifdef _CY_POINT_H_INCLUDED_
	void SetUniform(int index, const Point2<float>  &p)              { glUniform2fv (params[index],1,    &p.x ); }
	void SetUniform(int index, const Point3<float>  &p)              { glUniform3fv (params[index],1,    &p.x ); }
	void SetUniform(int index, const Point4<float>  &p)              { glUniform4fv (params[index],1,    &p.x ); }
	void SetUniform(int index, const Point2<int>    &p)              { glUniform2iv (params[index],1,    &p.x ); }
	void SetUniform(int index, const Point3<int>    &p)              { glUniform3iv (params[index],1,    &p.x ); }
	void SetUniform(int index, const Point4<int>    &p)              { glUniform4iv (params[index],1,    &p.x ); }
	void SetUniform(int index, const Point2<float>  *p, int count=1) { glUniform2fv (params[index],count,&p->x); }
	void SetUniform(int index, const Point3<float>  *p, int count=1) { glUniform3fv (params[index],count,&p->x); }
	void SetUniform(int index, const Point4<float>  *p, int count=1) { glUniform4fv (params[index],count,&p->x); }
	void SetUniform(int index, const Point2<int>    *p, int count=1) { glUniform2iv (params[index],count,&p->x); }
	void SetUniform(int index, const Point3<int>    *p, int count=1) { glUniform3iv (params[index],count,&p->x); }
	void SetUniform(int index, const Point4<int>    *p, int count=1) { glUniform4iv (params[index],count,&p->x); }
# ifdef GL_VERSION_3_0
	void SetUniform(int index, const Point2<GLuint> &p)              { glUniform2uiv(params[index],1,    &p.x ); }
	void SetUniform(int index, const Point3<GLuint> &p)              { glUniform3uiv(params[index],1,    &p.x ); }
	void SetUniform(int index, const Point4<GLuint> &p)              { glUniform4uiv(params[index],1,    &p.x ); }
	void SetUniform(int index, const Point2<GLuint> *p, int count=1) { glUniform2uiv(params[index],count,&p->x); }
	void SetUniform(int index, const Point3<GLuint> *p, int count=1) { glUniform3uiv(params[index],count,&p->x); }
	void SetUniform(int index, const Point4<GLuint> *p, int count=1) { glUniform4uiv(params[index],count,&p->x); }
# endif
# ifdef GL_VERSION_4_0
	void SetUniform(int index, const Point2<double> &p)              { glUniform2dv (params[index],1,    &p.x ); }
	void SetUniform(int index, const Point3<double> &p)              { glUniform3dv (params[index],1,    &p.x ); }
	void SetUniform(int index, const Point4<double> &p)              { glUniform4dv (params[index],1,    &p.x ); }
	void SetUniform(int index, const Point2<double> *p, int count=1) { glUniform2dv (params[index],count,&p->x); }
	void SetUniform(int index, const Point3<double> *p, int count=1) { glUniform3dv (params[index],count,&p->x); }
	void SetUniform(int index, const Point4<double> *p, int count=1) { glUniform4dv (params[index],count,&p->x); }
# endif
#endif

#ifdef _CY_IPOINT_H_INCLUDED_
	void SetUniform(int index, const IPoint2<int>    &p)              { glUniform2iv (params[index],1,    &p.x );
	void SetUniform(int index, const IPoint3<int>    &p)              { glUniform3iv (params[index],1,    &p.x );
	void SetUniform(int index, const IPoint4<int>    &p)              { glUniform4iv (params[index],1,    &p.x );
	void SetUniform(int index, const IPoint2<int>    *p, int count=1) { glUniform2iv (params[index],count,&p->x);
	void SetUniform(int index, const IPoint3<int>    *p, int count=1) { glUniform3iv (params[index],count,&p->x);
	void SetUniform(int index, const IPoint4<int>    *p, int count=1) { glUniform4iv (params[index],count,&p->x);
# ifdef GL_VERSION_3_0
	void SetUniform(int index, const IPoint2<GLuint> &p)              { glUniform2uiv(params[index],1,    &p.x );
	void SetUniform(int index, const IPoint3<GLuint> &p)              { glUniform3uiv(params[index],1,    &p.x );
	void SetUniform(int index, const IPoint4<GLuint> &p)              { glUniform4uiv(params[index],1,    &p.x );
	void SetUniform(int index, const IPoint2<GLuint> *p, int count=1) { glUniform2uiv(params[index],count,&p->x);
	void SetUniform(int index, const IPoint3<GLuint> *p, int count=1) { glUniform3uiv(params[index],count,&p->x);
	void SetUniform(int index, const IPoint4<GLuint> *p, int count=1) { glUniform4uiv(params[index],count,&p->x);
# endif
#endif

#ifdef _CY_MATRIX_H_INCLUDED_
	void SetUniform(int index, const Matrix2 <float>  &m)              { glUniformMatrix2fv  (params[index],1,    GL_FALSE,m.data ); }
	void SetUniform(int index, const Matrix3 <float>  &m)              { glUniformMatrix3fv  (params[index],1,    GL_FALSE,m.data ); }
	void SetUniform(int index, const Matrix4 <float>  &m)              { glUniformMatrix4fv  (params[index],1,    GL_FALSE,m.data ); }
	void SetUniform(int index, const Matrix2 <float>  *m, int count=1) { glUniformMatrix2fv  (params[index],count,GL_FALSE,m->data); }
	void SetUniform(int index, const Matrix3 <float>  *m, int count=1) { glUniformMatrix3fv  (params[index],count,GL_FALSE,m->data); }
	void SetUniform(int index, const Matrix4 <float>  *m, int count=1) { glUniformMatrix4fv  (params[index],count,GL_FALSE,m->data); }
# ifdef GL_VERSION_2_1
	void SetUniform(int index, const Matrix34<float>  &m)              { glUniformMatrix3x4fv(params[index],1,    GL_FALSE,m.data ); }
	void SetUniform(int index, const Matrix34<float>  *m, int count=1) { glUniformMatrix3x4fv(params[index],count,GL_FALSE,m->data); }
# endif
# ifdef GL_VERSION_4_0
	void SetUniform(int index, const Matrix2 <double> &m)              { glUniformMatrix2dv  (params[index],1,    GL_FALSE,m.data ); }
	void SetUniform(int index, const Matrix3 <double> &m)              { glUniformMatrix3dv  (params[index],1,    GL_FALSE,m.data ); }
	void SetUniform(int index, const Matrix4 <double> &m)              { glUniformMatrix4dv  (params[index],1,    GL_FALSE,m.data ); }
	void SetUniform(int index, const Matrix34<double> &m)              { glUniformMatrix3x4dv(params[index],1,    GL_FALSE,m.data ); }
	void SetUniform(int index, const Matrix2 <double> *m, int count=1) { glUniformMatrix2dv  (params[index],count,GL_FALSE,m->data); }
	void SetUniform(int index, const Matrix3 <double> *m, int count=1) { glUniformMatrix3dv  (params[index],count,GL_FALSE,m->data); }
	void SetUniform(int index, const Matrix4 <double> *m, int count=1) { glUniformMatrix4dv  (params[index],count,GL_FALSE,m->data); }
	void SetUniform(int index, const Matrix34<double> *m, int count=1) { glUniformMatrix3x4dv(params[index],count,GL_FALSE,m->data); }
# endif
#endif
	//!@}
};

//-------------------------------------------------------------------------------
// Implementation of GLDebugCallback
//-------------------------------------------------------------------------------

inline void _CY_APIENTRY GLDebugCallback::Callback( GLenum source,
                                                    GLenum type,
                                                    GLuint id,
                                                    GLenum severity,
                                                    GLsizei length,
                                                    const GLchar* message,
                                                    const void* userParam )
{
	std::ostream *outStream = (std::ostream*) userParam;

	*outStream << std::endl;
	*outStream << "OpenGL Debug Output:" << std::endl;
	*outStream << "VERSION:  ";
	GLPrintVersion(outStream);
	*outStream << std::endl;

	*outStream << "SOURCE:   ";
	switch (source) {
		case GL_DEBUG_SOURCE_API:             *outStream << "API";             break;
		case GL_DEBUG_SOURCE_WINDOW_SYSTEM:   *outStream << "Window System";   break;
		case GL_DEBUG_SOURCE_SHADER_COMPILER: *outStream << "Shader Compiler"; break;
		case GL_DEBUG_SOURCE_THIRD_PARTY:     *outStream << "Third Party";     break;
		case GL_DEBUG_SOURCE_APPLICATION:     *outStream << "Application";     break;
		case GL_DEBUG_SOURCE_OTHER:           *outStream << "Other";           break;
		default:                              *outStream << "Unknown";         break;
	}
	*outStream << std::endl;

	*outStream << "TYPE:     ";
	switch (type) {
		case GL_DEBUG_TYPE_ERROR:               *outStream << "Error";               break;
		case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: *outStream << "Deprecated Behavior"; break;
		case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:  *outStream << "Undefined Behavior";  break;
		case GL_DEBUG_TYPE_PORTABILITY:         *outStream << "Portability";         break;
		case GL_DEBUG_TYPE_PERFORMANCE:         *outStream << "Performance";         break;
		case GL_DEBUG_TYPE_OTHER:               *outStream << "Other";               break;
		default:                                *outStream << "Unknown";             break;
	}
	*outStream << std::endl;

	*outStream << "ID:       " << id << std::endl;

	*outStream << "SEVERITY: ";
	switch (severity) {
		case GL_DEBUG_SEVERITY_HIGH:         *outStream << "High";         break;
		case GL_DEBUG_SEVERITY_MEDIUM:       *outStream << "Medium";       break;
		case GL_DEBUG_SEVERITY_LOW:          *outStream << "Low";          break;
		case GL_DEBUG_SEVERITY_NOTIFICATION: *outStream << "Notification"; break;
		default:                             *outStream << "Unknown";      break;
	}
	*outStream << std::endl;

	*outStream << "MESSAGE:  ";
	*outStream << message << std::endl;

	// You can set a breakpoint at the following line. Your debugger will stop the execution,
	// and the call stack will show the OpenGL call causing this callback.
	*outStream << std::endl;
}

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
		glShaderSource(shaderID, prependSourceCount+1, sources.data(), nullptr);
	} else {
		glShaderSource(shaderID, 1, &shaderSourceCode, nullptr);
	}
	glCompileShader(shaderID);

	GLint result = GL_FALSE;
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &result);

	int infoLogLength;
	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
	if ( infoLogLength > 1 ) {
		std::vector<char> compilerMessage(infoLogLength);
		glGetShaderInfoLog( shaderID, infoLogLength, nullptr, compilerMessage.data() );
		if ( outStream ) {
			if ( !result ) *outStream << "ERROR: Cannot compile shader." << std::endl;
			*outStream << "OpenGL Version: ";
			GLPrintVersion(outStream);
			*outStream << std::endl;
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
		glGetProgramInfoLog( programID, infoLogLength, nullptr, compilerMessage.data() );
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
				case GL_INVALID_VALUE:     *outStream << "GL_INVALID_VALUE.";     break;
				case GL_INVALID_OPERATION: *outStream << "GL_INVALID_OPERATION."; break;
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

#ifdef GL_KHR_debug
typedef cy::GLDebugCallback cyGLDebugCallback;	//!< OpenGL debug callback class.
#endif

typedef cy::GLSLShader  cyGLSLShader;		//!< GLSL shader class
typedef cy::GLSLProgram cyGLSLProgram;		//!< GLSL program class

//-------------------------------------------------------------------------------
#endif
