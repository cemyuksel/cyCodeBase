// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyAlphaDistribution.h 
//! \author Cem Yuksel
//!
//! \brief  Implementation of the alpha distribution methods.
//!
//! This file includes an implementation of the alpha distribution methods using 
//! alpha pyramid and error diffusion used for pre-computing textures to be used
//! with alpha testing or sampleMask-to-alpha.
//!
//! Alpha testing using the original alpha values of a texture cannot handle
//! semi-transparent regions and often leads to problems with mipmapping. Alpha
//! distribution is a pre-processing approach that modifies the alpha values of
//! mipmap levels, such that they can properly handle semi-transparent regions.
//! 
//! The AlphaDistribution class provided in this file implements two methods
//! that can be used for alpha distribution: error diffusion and alpha pyramid.
//! Both methods produce similar results with minor differences.
//!
//! More details can be found in the original publication:
//!
//! Cem Yuksel. 2017. Alpha Distribution for Alpha Testing. PACM on CGIT (I3D 2018).
//! http://www.cemyuksel.com/research/alphadistribution/
//!
//-------------------------------------------------------------------------------
//
// Copyright (c) 2018, Cem Yuksel <cem@cemyuksel.com>
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

#ifndef _CY_ALPHA_DISTRIBUTION_H_INCLUDED_
#define _CY_ALPHA_DISTRIBUTION_H_INCLUDED_

//-------------------------------------------------------------------------------

#include <vector>
#include <random>
#include <assert.h>

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

//! An implementation of alpha distribution methods.
//! This implementation only works for textures with 8-bit channels.
//!
//! Cem Yuksel. 2017. Alpha Distribution for Alpha Testing. PACM on CGIT (I3D 2018).
//! http://www.cemyuksel.com/research/alphadistribution/

class AlphaDistribution
{
public:

	enum Method {
		METHOD_ERROR_DIF,	//!< Error diffusion using Floyd–Steinberg dithering
		METHOD_PYRAMID,		//!< Alpha Pyramid
	};

	//! Fixes the alpha values of an image with the given width and height using the specified method.
	//! If the image will be used with alpha-to-coverage, the spp parameter
	//! should indicate the number of alpha samples; otherwise, it should be 1.
	static void FixAlpha( Method method, unsigned char *alpha, int width, int height, int spp=1 )
	{
		assert(spp>=1);
		if ( width*height == 1 ) {
			alpha[0] = 255;  // fix for single pixel
		} else {
			switch( method ) {
				case METHOD_ERROR_DIF: ErrorDiffusion<1>(alpha,width,height,spp); break;
				case METHOD_PYRAMID:   AlphaPyramid  <1>(alpha,width,height,spp); break;
			}
		}
	}

	//! Fixes the alpha values of an RGBA image with the given width and height using the specified method.
	//! Only the alpha channel of the image is modified, the RGB values are not altered.
	//! If the image will be used with alpha-to-coverage, the spp parameter
	//! should indicate the number of alpha samples; otherwise, it should be 1.
	static void FixAlphaRGBA( Method method, unsigned char *image, int width, int height, int spp=1 )
	{
		assert(spp>=1);
		if ( width*height == 1 ) {
			image[3] = 255;  // fix for single pixel
		} else {
			switch( method ) {
				case METHOD_ERROR_DIF: ErrorDiffusion<4>(image,width,height,spp); break;
				case METHOD_PYRAMID:   AlphaPyramid  <4>(image,width,height,spp); break;
			}
		}
	}

	//! Generates a sample mask texture from the alpha values of an image with the given width and height.
	//! This sample mask texture can be used with the original texture for handling alpha-to-coverage.
	//! The spp parameter is the number of alpha samples that will be used.
	template <typename SAMPLE_MASK_TYPE=unsigned char>
	static void GenerateSampleMaskTexture( SAMPLE_MASK_TYPE *sampleMask, unsigned char const *alpha, int width, int height, int spp )
	{
		GenSampleMaskTexture<1>(sampleMask,alpha,width,height,spp);
	}

	//! Generates a sample mask texture from the alpha values of an RGBA image with the given width and height.
	//! This sample mask texture can be used with the original texture for handling alpha-to-coverage.
	//! The spp parameter is the number of alpha samples that will be used.
	template <typename SAMPLE_MASK_TYPE=unsigned char>
	static void GenerateSampleMaskTextureRGBA( SAMPLE_MASK_TYPE *sampleMask, unsigned char const *image, int width, int height, int spp )
	{
		GenSampleMaskTexture<4>(sampleMask,alpha,width,height,spp);
	}

#if defined(__gl_h_) || defined(__GL_H__) || defined(_GL_H) || defined(__X_GL_H)

	//!@name OpenGL Methods

	//! Fixes the alpha values of a mipmap level of an OpenGL texture.
	//! If the image will be used with alpha-to-coverage, the spp parameter
	//! should indicate the number of alpha samples; otherwise, it should be 1.
	static void FixTextureLevelAlpha( Method method, GLuint textureID, int level, int spp=1 )
	{
		int width=0, height=0;
		glBindTexture(GL_TEXTURE_2D, textureID);
		glGetTexLevelParameteriv(GL_TEXTURE_2D,level,GL_TEXTURE_WIDTH,&width);
		glGetTexLevelParameteriv(GL_TEXTURE_2D,level,GL_TEXTURE_HEIGHT,&height);
		if ( width * height == 0 ) return;

		std::vector<unsigned char> image( width * height *4 );
		glGetTexImage(GL_TEXTURE_2D,level,GL_RGBA,GL_UNSIGNED_BYTE,image.data());
		FixAlphaRGBA( method, image.data(), width, height, spp );
		glTexImage2D(GL_TEXTURE_2D,level,GL_RGBA,width,height,0,GL_RGBA,GL_UNSIGNED_BYTE,image.data());
	}

	//! Fixes the alpha values of all mipmap levels of an OpenGL texture, starting with the given level.
	//! If the texture does not contain semi-transparent regions, modifying the first level (level zero)
	//! is not advisable, since the original values might work better with magnification filtering.
	//! If the image will be used with alpha-to-coverage, the spp parameter
	//! should indicate the number of alpha samples; otherwise, it should be 1.
	static void FixTextureAlpha( Method method, GLuint textureID, int startingLevel=0, int spp=1 )
	{
		int level = startingLevel;
		int width, height;
		glGetTexLevelParameteriv(GL_TEXTURE_2D,level,GL_TEXTURE_WIDTH, &width );
		glGetTexLevelParameteriv(GL_TEXTURE_2D,level,GL_TEXTURE_HEIGHT,&height);

		while ( width >= 1 && height >= 1 ) {
			FixTextureLevelAlpha( method, textureID, level, spp );
			level++;
			glGetTexLevelParameteriv(GL_TEXTURE_2D,level,GL_TEXTURE_WIDTH, &width );
			glGetTexLevelParameteriv(GL_TEXTURE_2D,level,GL_TEXTURE_HEIGHT,&height);
		}
	}

	//! Fixes the alpha values of all mipmap levels of an OpenGL texture, starting with the given level.
	//! If the texture does not contain semi-transparent regions, modifying the first level (level zero)
	//! is not advisable, since the original values might work better with magnification filtering.
	//! The number of alpha samples for alpha-to-coverage is obtained from the current OpenGL context.
	static void FixTextureAlphaToSampleMask( Method method, GLuint textureID, int startingLevel=0 )
	{
		GLint spp;
        glGetIntegerv(GL_SAMPLES, &spp);
		if ( spp < 1 ) spp = 1;
		FixTextureAlpha( method, textureID, startingLevel, spp );
	}

#endif

private:

	template <int NUM_CHANNELS>
	static void ErrorDiffusion( unsigned char *image, int width, int height, int spp )
	{
		auto addError = [&]( int ix, int err ) {
			int a = image[ix] + err;
			if ( a < 0 ) a = 0;
			if ( a > 255 ) a = 255;
			image[ix] = a;
		};

		for ( int i=0, ih=0; ih<height; ih++ ) {
			for ( int iw=0; iw<width; iw++, i++ ) {
				int a0 = image[i*NUM_CHANNELS+(NUM_CHANNELS-1)];	// current value
				int a1 = a0 >= 128 ? 255 : 0;
				if ( spp > 1 ) {
					for ( int j=1; j<=spp; j++ ) {
						int cutoff = (256*(j*2-1)) / (spp*2);
						if ( a0 < cutoff ) break;
						a1 = (256*j) / spp;
					}
					if ( a1 > 255 ) a1 = 255;
				}
				image[i*NUM_CHANNELS+(NUM_CHANNELS-1)] = a1;
				int err = a0 - a1;
				int e[4] = { 7*err/16, 3*err/16, 5*err/16, 1*err/16 };
				int de = err - (e[0]+e[1]+e[2]+e[3]);
				e[0] += de;
				if ( iw < width-1 ) addError( (i+1)*4+3, e[0] );
				if ( ih < height-1 ) {
					if ( iw > 0 ) addError( (width+i-1)*4+3, e[1] );
					addError( (width+i)*4+3, e[2] );
					if ( iw < width-1 ) addError( (width+i+1)*4+3, e[3] );
				}
			}
		}
	}

	struct AlphaPyramidLevel
	{
		int width, height;
		std::vector<uint32_t> alpha;
		uint32_t total_alpha;

		uint32_t GetAlpha(int x, int y) const { return alpha[ y*width + x ]; }
		template <typename ARG_SCALE> void SetData( int w, int h, int prev_width, int prev_height, ARG_SCALE accessor )
		{
			total_alpha = 0;
			width  = w;
			height = h;
			alpha.resize(width*height);
			for ( int ih=0; ih<height; ih++ ) {
				for ( int iw=0; iw<width; iw++ ) {
					uint32_t a0 = accessor( (ih*2    )*prev_width + iw*2      );
					uint32_t a1 = accessor( (ih*2    )*prev_width + iw*2 + 1  );
					uint32_t a2 = accessor( (ih*2 + 1)*prev_width + iw*2      );
					uint32_t a3 = accessor( (ih*2 + 1)*prev_width + iw*2 + 1  );
					alpha[ih*width+iw] = a0 + a1 + a2 + a3;
					total_alpha += a0 + a1 + a2 + a3;
				}
				if ( width*2 < prev_width ) {
					uint32_t a0 = accessor( (ih*2    )*prev_width + width*2 );
					uint32_t a1 = accessor( (ih*2 + 1)*prev_width + width*2 );
					alpha[(ih+1)*width-1] += a0 + a1;
					total_alpha += a0 + a1;
				}
			}
			if ( height*2 < prev_height ) {
				int ii = (height-1)*width;
				for ( int iw=0; iw<width; iw++ ) {
					uint32_t a0 = accessor( (height*2)*prev_width + iw*2     );
					uint32_t a1 = accessor( (height*2)*prev_width + iw*2 + 1 );
					alpha[ii+iw] += a0 + a1;
					total_alpha += a0 + a1;
				}
				if ( width*2 < prev_width ) {
					uint32_t a0 = accessor( (height*2)*prev_width + width*2 );
					alpha[ii+width-1] += a0;
					total_alpha += a0;
				}
			}
		}
		void Alpha2CountBlock( int const *ix, int n, uint32_t count, int spp )
		{
			uint32_t sum = 0, remSum = 0;
			for ( int j=0; j<n; j++ ) { 
				uint32_t v = alpha[ ix[j] ];
				uint32_t c = v / 255;
				uint32_t r = v % 255;
				if ( spp > 1 ) {
					c = (v*spp) / 255;
					r = (v*spp) % 255;
				}
				sum += c;
				remSum += r;
				alpha[ ix[j] ] = (c << 8) | r;
			}
			assert( sum <= count );
			assert( sum + remSum/255 + 1 >= count);

			for ( uint32_t rem = count - sum; rem > 0; rem-- ) {
				int max_i = 0;
				uint32_t max_v;
				int eqCount = 1;
				max_v = alpha[ ix[0] ] & 255;
				for ( int j=1; j<n; j++ ) {
					uint32_t v = alpha[ ix[j] ] & 255;
					if ( max_v < v ) {
						max_v = v;
						max_i = j;
						eqCount = 1;
					}
				}
				assert(max_v > 0);
				alpha[ ix[max_i] ] += 256 - max_v;
			}
			for ( int j=0; j<n; j++ ) alpha[ ix[j] ] >>= 8;
		}
		void Alpha2Count( AlphaPyramidLevel const *parent, int spp )
		{
			int hLim = (height&1) ? height-3 : height;
			int wLim = (width &1) ? width -3 : width;
			for ( int ih=0; ih<hLim; ih+=2 ) {
				for ( int iw=0; iw<wLim; iw+=2 ) {
					uint32_t count = parent->GetAlpha(iw/2,ih/2);
					int i = ih*width + iw;
					int ix[] = { i, i+width+1, i+1, i+width };
					Alpha2CountBlock( ix, 4, count, spp );
				}
				if ( wLim < width ) {
					uint32_t count = parent->GetAlpha(wLim/2,ih/2);
					int i = ih*width + wLim;
					int ix[] = { i, i+width+1, i+2, i+width, i+1, i+width+2 };
					Alpha2CountBlock( ix, 6, count, spp );
				}
			}
			if ( hLim < height ) {
				for ( int iw=0; iw<wLim; iw+=2 ) {
					uint32_t count = parent->GetAlpha(iw/2,hLim/2);
					int i = hLim*width + iw;
					int ix[] = { i, i+width+1, i+width+width, i+1, i+width, i+width+width+1 };
					Alpha2CountBlock( ix, 6, count, spp );
				}
				if ( wLim < width ) {
					uint32_t count = parent->GetAlpha(wLim/2,hLim/2);
					int i = hLim*width + wLim;
					int ix[] = { i, i+width+width+2, i+2, i+width+width, i+width+1, i+1, i+width+width+1, i+width, i+width+2 };
					Alpha2CountBlock( ix, 9, count, spp );
				}
			}
		}
		void Alpha2CountSimple( uint32_t count, int spp )
		{
			int n = width*height;
			std::vector<int> ix(n);
			for ( int i=0; i<n; i++ ) ix[i] = i;
			Alpha2CountBlock( ix.data(), n, count, spp );
		}
	};

	template <int NUM_CHANNELS>
	static void AlphaPyramid( unsigned char *image, int width, int height, int spp )
	{
		// Step 1: Compute the alpha pyramid
		std::vector<AlphaPyramidLevel*> pyramid;
		int pw = width / 2;
		int ph = height / 2;
		if ( pw>0 && ph>0 ) {
			// First level
			uint32_t total_alpha = 0;
			AlphaPyramidLevel *lev = new AlphaPyramidLevel;
			lev->SetData( pw, ph, width, height, [&](int i){ return image[i*NUM_CHANNELS+(NUM_CHANNELS-1)]; } );
			pyramid.push_back(lev);
			// Higher levels
			AlphaPyramidLevel *pLev = lev;
			while ( pw>1 && ph>1 ) {
				int ww = pw;
				int hh = ph;
				pw /= 2;
				ph /= 2;
				AlphaPyramidLevel *lev = new AlphaPyramidLevel;
				lev->SetData( pw, ph, ww, hh, [&](int i){ return pLev->alpha[i]; } );
				pyramid.push_back(lev);
				pLev = lev;
			}
		}
		// Step 2: Compute the number of texels that should pass the alpha test
		uint32_t total_alpha = 0;
		if ( width > 1 && height > 1 ) {
			AlphaPyramidLevel *level = pyramid.back();
			for ( int i=0; i<(int)level->alpha.size(); i++ ) {
				total_alpha += level->alpha[i];
			}
		} else {
			for ( int i=0; i<height*width; i++ ) {
				uint32_t a0 = image[i*NUM_CHANNELS+(NUM_CHANNELS-1)];	// current value
				total_alpha += a0;
			}
		}
		int on_texels = (total_alpha*spp + 254) / 255;

		// Step 3: Alpha to Count
		int lvl = (int) pyramid.size() - 1;
		if ( lvl >= 0 ) {
			pyramid[lvl]->Alpha2CountSimple( on_texels, spp );
			for ( lvl--; lvl>=0; lvl-- ) {
				pyramid[lvl]->Alpha2Count( pyramid[lvl+1], spp );
			}
		}
		// Step 4: Update texture
		{
			auto setImgAlpha = [&]( int const *ix, int n, uint32_t count )
			{
				if ( spp > 1 ) {
					count *= 256/spp;
				}

				unsigned char a[9];
				assert(n<=9);
				for ( int j=0; j<n; j++ ) a[j] = 0;
				uint32_t remSum = 0;
				uint32_t rem = count;
				while ( rem > 0 ) {
					int max_i = 0;
					int eqCount = 1;
					int max_v = image[ ix[0]*NUM_CHANNELS+(NUM_CHANNELS-1) ] - a[0];
					if ( max_v < 0 ) max_v = 0;
					for ( int j=1; j<n; j++ ) {
						int v = image[ ix[j]*NUM_CHANNELS+(NUM_CHANNELS-1) ] - a[j];
						if ( v < 0 ) v = 0;
						if ( max_v < v ) {
							max_v = v;
							max_i = j;
							eqCount = 1;
						}
					}
					assert(max_v > 0);
												
					if ( spp > 1 ) {
						unsigned char r = image[ ix[max_i]*NUM_CHANNELS+(NUM_CHANNELS-1) ] - a[max_i];
						int inc = (r > 256/spp) ? (256/spp) : r;
						if ( inc == 0 ) break;
						int new_a = a[max_i] + (256/spp);
						if ( new_a > 255 ) new_a = 255;
						a[max_i] = new_a;
						remSum -= inc;
						if ( rem < uint32_t(256/spp) ) rem = 0;
						else rem -= uint32_t(256/spp);
					} else {
						a[max_i] = image[ ix[max_i]*NUM_CHANNELS+(NUM_CHANNELS-1) ];
						remSum -= max_v;
						rem--;
					}
				}

				if ( spp > 1 ) {
					for ( int j=0; j<n; j++ ) {
						int av = ((a[j] + 128 / spp) / (256/spp)) * (256/spp) + 1;
						if ( av > 255 ) av = 255;
						image[ ix[j]*NUM_CHANNELS+(NUM_CHANNELS-1) ] = av;
					}
				} else {
					for ( int j=0; j<n; j++ ) image[ ix[j]*NUM_CHANNELS+(NUM_CHANNELS-1) ] = a[j] ? 255 : 0;
				}
			};

			if ( pyramid.size() > 0 ) {
				int hLim = (height&1) ? height-3 : height;
				int wLim = (width&1) ? width -3 : width;
				for ( int ih=0; ih<hLim; ih+=2 ) {
					for ( int iw=0; iw<wLim; iw+=2 ) {
						uint32_t count = pyramid[0]->GetAlpha(iw/2,ih/2);
						int i = ih*width + iw;
						int ix[] = { i, i+width+1, i+1, i+width };
						setImgAlpha( ix, 4, count );
					}
					if ( wLim < width ) {
						uint32_t count = pyramid[0]->GetAlpha(wLim/2,ih/2);
						int i = ih*width + wLim;
						int ix[] = { i, i+width+1, i+2, i+width, i+1, i+width+2 };
						setImgAlpha( ix, 6, count );
					}
				}
				if ( hLim < height ) {
					for ( int iw=0; iw<wLim; iw+=2 ) {
						uint32_t count = pyramid[0]->GetAlpha(iw/2,hLim/2);
						int i = hLim*width + iw;
						int ix[] = { i, i+width+1, i+width+width, i+1, i+width, i+width+width+1 };
						setImgAlpha( ix, 6, count );
					}
					if ( wLim < width ) {
						uint32_t count = pyramid[0]->GetAlpha(wLim/2,hLim/2);
						int i = hLim*width + wLim;
						int ix[] = { i, i+width+width+2, i+2, i+width+width, i+width+1, i+1, i+width+width+1, i+width, i+width+2 };
						setImgAlpha( ix, 9, count );
					}
				}
			} else {
				// no pyramid
				int ix[9];
				for ( int i=0; i<9; i++ ) ix[i] = i;
				setImgAlpha( ix, width*height, on_texels );
			}
		}
		// Step 5: Clean up
		for ( int i=0; i<(int)pyramid.size(); i++ ) delete pyramid[i];
		pyramid.clear();
	}


	template <int NUM_CHANNELS, typename SAMPLE_MASK_TYPE=unsigned char>
	static void GenSampleMaskTexture( SAMPLE_MASK_TYPE *sampleMask, unsigned char const *image, int width, int height, int spp )
	{
		std::random_device rd;
		std::mt19937 gen( rd() );
		std::uniform_int_distribution<int> rnd(0,spp-1);

		for ( int i=0; i<width * height; i++ ) {
			int a = image[i*NUM_CHANNELS+(NUM_CHANNELS-1)];
			int a1 = 0;
			for ( int j=1; j<=spp; j++ ) {
				int cutoff = (256*(j*2-1)) / (spp*2);
				if ( a < cutoff ) break;
				a1 = (256*j) / spp;
			}
			int n = a1 / (256 / spp);
			int sc = ( n > spp/2 ) ? spp - n : n;
			int bits = 0;
			for ( int j=0; j<sc; j++ ) {
				while ( true ) {
					int r = rnd(gen) % spp;
					if ( (bits & (1<<r)) == 0 ) {
						bits |= (1<<r);
						break;
					}
				}
			}
			if ( n > spp/2 ) {
				bits = ( (1 << spp) - 1 ) & (~bits);
			}
			sampleMask[i] = SAMPLE_MASK_TYPE(bits);
		}
	}
};

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::AlphaDistribution cyAlphaDistribution;		//!< An implementation of alpha distribution methods

//-------------------------------------------------------------------------------

#endif
