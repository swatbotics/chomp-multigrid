/*
* Copyright (c) 2008-2014, Matt Zucker
*
* This file is provided under the following "BSD-style" License:
*
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the following
* conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and/or other materials provided
* with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#ifndef _QUICKPNG_H_
#define _QUICKPNG_H_

#include <stdio.h>
#include <png.h>
#include "vec2u.h"
#include "Gradient.h"
#include <vector>

template<class Tval> inline bool quickPNG(const std::string& filename,
                                          const Tval* data,
                                          size_t ncols,
                                          size_t nrows,
                                          size_t rowsz,
                                          bool yflip,
                                          const Gradient& g,
                                          Tval minVal=0, Tval maxVal=0) {


  if (minVal >= maxVal) {
    minVal = maxVal = data[0];
    const Tval* rowptr = data;
    for (size_t y=0; y<nrows; ++y) {
      for (size_t x=0; x<ncols; ++x) {
        const Tval* val = rowptr+x;
        minVal = std::min(minVal, *val);
        maxVal = std::max(maxVal, *val);
      }
      rowptr += rowsz;
    }
  }

  FILE* fp = fopen(filename.c_str(), "wb");
  if (!fp) { 
    std::cerr << "couldn't open " << filename << " for output!\n";
    return false;
  }
  
  png_structp png_ptr = png_create_write_struct
    (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr) {
    std::cerr << "error creating png write struct\n";
    return false;
  }
  
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    std::cerr << "error creating png info struct\n";
    png_destroy_write_struct(&png_ptr,
			     (png_infopp)NULL);
    fclose(fp);
    return false;
  }  

  if (setjmp(png_jmpbuf(png_ptr))) {
    std::cerr << "error in png processing\n";
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    return false;
  }

  png_init_io(png_ptr, fp);

  png_set_IHDR(png_ptr, info_ptr, 
	       ncols, nrows,
	       8, 
	       PNG_COLOR_TYPE_RGB_ALPHA,
	       PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT,
	       PNG_FILTER_TYPE_DEFAULT);

  png_write_info(png_ptr, info_ptr);

  std::vector<unsigned char> buf(ncols*4);
  float rng = maxVal-minVal;
  if (!rng) { rng = 1; }

  const Tval* rowptr = data + (yflip ? rowsz*(nrows-1) : 0);

  for (size_t y=0; y<nrows; ++y) {
    buf.clear(); 
    for (size_t x=0; x<ncols; ++x) {
      const Tval* pval = rowptr+x;
      float u = (*pval-minVal)/rng;
      vec3f c = g.lookup(u);
      for (int i=0; i<3; ++i) {
        int v = std::min(255, (int)ceil(c[i]*255));
        buf.push_back(v);
      }
      buf.push_back(255);
    }
    png_write_row(png_ptr, (png_bytep)&(buf[0]));
    if (yflip) {
      rowptr -= rowsz;
    } else {
      rowptr += rowsz;
    }
  }

  png_write_end(png_ptr, info_ptr);

  png_destroy_write_struct(&png_ptr, &info_ptr);

  fclose(fp);

  return true;


}

inline bool quickPNG(const std::string& filename,
                     const vec3f* data,
                     size_t ncols,
                     size_t nrows,
                     size_t rowsz,
                     bool yflip) {


  FILE* fp = fopen(filename.c_str(), "wb");
  if (!fp) { 
    std::cerr << "couldn't open " << filename << " for output!\n";
    return false;
  }
  
  png_structp png_ptr = png_create_write_struct
    (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr) {
    std::cerr << "error creating png write struct\n";
    return false;
  }
  
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    std::cerr << "error creating png info struct\n";
    png_destroy_write_struct(&png_ptr,
			     (png_infopp)NULL);
    fclose(fp);
    return false;
  }  

  if (setjmp(png_jmpbuf(png_ptr))) {
    std::cerr << "error in png processing\n";
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    return false;
  }

  png_init_io(png_ptr, fp);

  png_set_IHDR(png_ptr, info_ptr, 
	       ncols, nrows,
	       8, 
	       PNG_COLOR_TYPE_RGB_ALPHA,
	       PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT,
	       PNG_FILTER_TYPE_DEFAULT);

  png_write_info(png_ptr, info_ptr);

  std::vector<unsigned char> buf(ncols*4);

  const vec3f* rowptr = data + (yflip ? rowsz*(nrows-1) : 0);

  for (size_t y=0; y<nrows; ++y) {
    buf.clear(); 
    for (size_t x=0; x<ncols; ++x) {
      const vec3f& c = rowptr[x];
      for (int i=0; i<3; ++i) {
        int v = std::min(255, (int)ceil(c[i]*255));
        buf.push_back(v);
      }
      buf.push_back(255);
    }
    png_write_row(png_ptr, (png_bytep)&(buf[0]));
    if (yflip) {
      rowptr -= rowsz;
    } else {
      rowptr += rowsz;
    }
  }

  png_write_end(png_ptr, info_ptr);

  png_destroy_write_struct(&png_ptr, &info_ptr);

  fclose(fp);

  return true;


}


#endif
