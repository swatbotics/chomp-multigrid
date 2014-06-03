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

#include "GlFrameGrabber.h"

#include "glstuff.h"
#include <png.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <iostream>

class GlFrameGrabberPrivate {
public:
  std::string format;
  unsigned int w;
  unsigned int h;
  unsigned int stride;
  unsigned int count;
  std::vector<unsigned char> buffer;
};

GlFrameGrabber::GlFrameGrabber(unsigned int width,
			       unsigned int height,
			       const char* format) {

  _pimpl = new GlFrameGrabberPrivate;
  _pimpl->format = format;
  _pimpl->w = width;
  _pimpl->h = height;
  _pimpl->count = 0;
  _pimpl->stride = width * 4;
  _pimpl->buffer.resize(width * height * 4);

  for (int i=0; i<999; ) {
    while (1) {
      char buf[1024];
      snprintf(buf, 1024, _pimpl->format.c_str(), i++);
      struct stat sb;
      if (stat(buf, &sb) || !S_ISREG(sb.st_mode)) { break; }
      _pimpl->count = i;
    }
  }

}

GlFrameGrabber::~GlFrameGrabber() {
  delete _pimpl; 
  _pimpl = 0;
}

void GlFrameGrabber::deleteAll() {
  for (int i=0; i<999; ) {
    while (1) {
      char buf[1024];
      snprintf(buf, 1024, _pimpl->format.c_str(), i++);
      struct stat sb;
      if (stat(buf, &sb) || !S_ISREG(sb.st_mode)) { break; }
      unlink(buf);
    }
  }
}

void GlFrameGrabber::resetCount() {
  _pimpl->count = 0;
}

unsigned int GlFrameGrabber::count() const {
  return _pimpl->count;
}

void GlFrameGrabber::writeFrame() {

  glFlush();

  glReadPixels(0, 0,
	       _pimpl->w, _pimpl->h, 
	       GL_RGBA, GL_UNSIGNED_BYTE, 
	       &(_pimpl->buffer[0]));


  // now write the png!
  char buf[1024];
  snprintf(buf, 1024, _pimpl->format.c_str(), _pimpl->count++);

  FILE *fp = fopen(buf, "wb");
  if (!fp) {
    std::cerr << "error opening PNG for output\n";
    return;
  }

  png_structp png_ptr = png_create_write_struct
    (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr) {
    std::cerr << "error creating png write struct\n";
    return;
  }
  
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    std::cerr << "error creating png info struct\n";
    png_destroy_write_struct(&png_ptr,
			     (png_infopp)NULL);
    fclose(fp);
    return;
  }  

  if (setjmp(png_jmpbuf(png_ptr))) {
    std::cerr << "error in png processing\n";
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    return;
  }

  png_init_io(png_ptr, fp);

  png_set_IHDR(png_ptr, info_ptr, 
	       _pimpl->w,
	       _pimpl->h,
	       8, 
	       PNG_COLOR_TYPE_RGB_ALPHA,
	       PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT,
	       PNG_FILTER_TYPE_DEFAULT);

  png_write_info(png_ptr, info_ptr);
  
  unsigned char* rowptr = &(_pimpl->buffer[_pimpl->stride * (_pimpl->h-1)]);
  //unsigned char* rowptr = &(_pimpl->buffer[0]);

  for (unsigned int y=0; y<_pimpl->h; ++y) {
    png_write_row(png_ptr, (png_bytep)rowptr);
    rowptr -= _pimpl->stride;
    //rowptr += _pimpl->stride;
  }
  
  png_write_end(png_ptr, info_ptr);

  png_destroy_write_struct(&png_ptr, &info_ptr);

  fclose(fp);

}

unsigned int GlFrameGrabber::width() const {
  return _pimpl->w;
}

unsigned int GlFrameGrabber::height() const {
  return _pimpl->h;
}


  
