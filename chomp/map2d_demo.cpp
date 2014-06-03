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

#include "Map2D.h"
#include <png.h>

#ifdef MZ_HAVE_CAIRO
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>
//#include <cairo/cairo-image.h>
#endif

bool savePNG_RGB24(const std::string& filename,
                   size_t ncols, size_t nrows, 
                   size_t rowsz,
                   const unsigned char* src,
                   bool yflip=false) {

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
	       PNG_COLOR_TYPE_RGB,
	       PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT,
	       PNG_FILTER_TYPE_DEFAULT);

  png_write_info(png_ptr, info_ptr);

  std::vector<unsigned char> buf(ncols*4);

  const unsigned char* rowptr;
  if (yflip) {
    rowptr = src + rowsz * (nrows-1);
  } else {
    rowptr = src;
  }

  for (size_t y=0; y<nrows; ++y) {
    const unsigned char* pxptr = rowptr;
    buf.clear(); 
    for (size_t x=0; x<ncols; ++x) {
      buf.push_back(*pxptr++);
      buf.push_back(*pxptr++);
      buf.push_back(*pxptr++);
      pxptr++;
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

int main(int argc, char** argv) {

  Map2D map;

  map.load(argv[1]);

  std::vector<unsigned char> buf;

  map.rasterize(Map2D::RASTER_DISTANCE, buf, 0);
  savePNG_RGB24("dist.png", map.grid.nx(), map.grid.ny(), 
                map.grid.nx()*4, &buf[0]);

  map.rasterize(Map2D::RASTER_COST, buf, 0);
  savePNG_RGB24("cost.png", map.grid.nx(), map.grid.ny(), 
                map.grid.nx()*4, &buf[0]);

  map.rasterize(Map2D::RASTER_OCCUPANCY, buf, 0);
  savePNG_RGB24("occupancy.png", map.grid.nx(), map.grid.ny(), 
                map.grid.nx()*4, &buf[0]);

#ifdef MZ_HAVE_CAIRO

  Box3f bbox = map.grid.bbox();
  vec3f dims = bbox.p1 - bbox.p0;

  // mscl = DISPLAY / MAP
  float mscl = (400) / std::max(dims.x(), dims.y());
  float cs = map.grid.cellSize();
  float mpt = 1.0/mscl;

  int width = int(mscl * dims.x());
  int height = int(mscl * dims.y());
  std::cout << "image will be " << width << "x" << height << "\n";

  cairo_surface_t* surface = cairo_pdf_surface_create("map2d.pdf",
                                                      width, height);

  cairo_t* cr = cairo_create(surface);

  size_t stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, map.grid.nx());
  map.rasterize(Map2D::RASTER_OCCUPANCY, buf, stride);

  cairo_surface_t* image = 
    cairo_image_surface_create_for_data(&(buf[0]), 
                                        CAIRO_FORMAT_RGB24,
                                        map.grid.nx(), map.grid.ny(),
                                        stride);

  cairo_identity_matrix(cr);
  cairo_translate(cr, 0, height);
  cairo_scale(cr, mscl, -mscl);
  cairo_translate(cr, -bbox.p0.x(), -bbox.p0.y());

  cairo_save(cr);
  cairo_scale(cr, cs, cs);
  cairo_set_source_surface(cr, image, bbox.p0.x()/cs, bbox.p0.y()/cs);
  cairo_paint(cr);
  cairo_restore(cr);
  
  cairo_set_line_width(cr, 0.5*cs);

  float x0 = bbox.p0.x();
  float y0 = bbox.p0.y();

  float x1 = 0.5 * (bbox.p0.x() + bbox.p1.x());
  float y1 = 0.5 * (bbox.p0.y() + bbox.p1.y());

  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_arc(cr, x0, y0, 10.0*mpt, 0.0, 2*M_PI);
  cairo_stroke(cr);

  cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
  cairo_arc(cr, x1, y0, 10.0*mpt, 0.0, 2*M_PI);
  cairo_stroke(cr);

  cairo_set_source_rgb(cr, 0.0, 1.0, 0.0);
  cairo_arc(cr, x0, y1, 10.0*mpt, 0.0, 2*M_PI);
  cairo_stroke(cr);

  
  for (int i=0; i<500; ++i) {
    vec3f p(0.0);
    for (int j=0; j<2; ++j) {
      p[j] = (double(rand())/RAND_MAX)*dims[j] + bbox.p0[j];
    } 
    float d = map.grid.sample(p);
    if (d < 0) {
      cairo_set_source_rgb(cr, 0.6, 0.0, 0.0);
    } else {
      cairo_set_source_rgb(cr, 0.0, 0.0, 0.6);
    }
    cairo_arc(cr, p.x(), p.y(), cs, 0.0, 2*M_PI);
    cairo_stroke(cr);
  }

  cairo_surface_destroy(image);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);

#endif

  return 0;

}
