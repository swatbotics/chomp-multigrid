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
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <png.h>
#include <mzcommon/QuickPNG.h>
#include <mzcommon/strutils.h>

struct Image {
public:

  size_t nx, ny;
  float ox, oy, sx, sy;
  std::vector<unsigned char> buf;

  bool contains(const vec2f& v) {
    return (v.x() >= ox &&
            v.x() <= ox + sx &&
            v.y() >= oy && 
            v.y() <= oy + sy);
  }

  vec2f p2v(size_t x, size_t y) {
    return vec2f(ox + (x+0.5)*sx/nx, oy + sy - (y+0.5)*sy/ny);
  }

  vec2u v2p(const vec2f& v) {
    return vec2u( (v.x() - ox)*nx/sx, (sy - v.y() + oy)*ny/sy );
  }

  vec2u ind2sub(size_t i) {
    return vec2u(i % nx, i / nx);
  }

  size_t sub2ind(const vec2u& s) {
    return s.x() + s.y()*nx;
  }

  void open(const char* filename) {

    FILE* fp = fopen(filename, "rb");
    if (!fp) {
      std::cerr << "couldn't open " << filename << " for input!\n";
      return;
    }

    png_byte header[8];
    fread(header, 1, 8, fp);

    if (png_sig_cmp(header, 0, 8) != 0) {
      std::cerr << "not a PNG file!\n";
      fclose(fp);
      return;
    }

    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, 
                                                 NULL, NULL, NULL);

    if (!png_ptr) {
      std::cerr << "error creating png_ptr!\n";
      fclose(fp);
      return;
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
      std::cerr << "error creating info_ptr!\n";
    
      png_destroy_read_struct(&png_ptr,
                              (png_infopp)NULL, (png_infopp)NULL);
      fclose(fp);
      return;
    }

    png_infop end_info = png_create_info_struct(png_ptr);
    if (!end_info) {
      std::cerr << "error creating end_info!\n";
      png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
      fclose(fp);
      return;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
      std::cerr << "error reading PNG!\n";
      png_destroy_read_struct(&png_ptr, &info_ptr,
                              &end_info);
      fclose(fp);
      return;
    }

    png_init_io(png_ptr, fp);

    png_set_sig_bytes(png_ptr, 8);
  
    png_read_info(png_ptr, info_ptr);

    png_uint_32 width, height;
    int bit_depth, color_type, interlace_type, compression_type, filter_method;

    png_get_IHDR(png_ptr, info_ptr, &width, &height,
                 &bit_depth, &color_type, &interlace_type,
                 &compression_type, &filter_method);

    if (bit_depth != 8 ||(color_type != PNG_COLOR_TYPE_GRAY)) {
      std::cerr << "not an 8-bit uninterlaced image!\n";
      std::cerr << "width = " << width << "\n";
      std::cerr << "height = " << height << "\n";
      std::cerr << "bit depth = " << bit_depth << "\n";
      std::cerr << "color type = " << color_type << "\n";
      png_destroy_read_struct(&png_ptr, &info_ptr,
                              &end_info);
      fclose(fp);
      return;
    }

    buf.resize(width*height);

    int number_passes = png_set_interlace_handling(png_ptr);
  
    for (int pass = 0; pass<number_passes; ++pass) {
      unsigned char* dst = &(buf[0]);
      for (size_t y=0; y<height; ++y) {
        png_bytep row_pointers = (png_bytep)dst;
        png_read_rows(png_ptr, &row_pointers, NULL, 1);
        dst += width;
      }
    }
  
    png_destroy_read_struct(&png_ptr, &info_ptr,
                            &end_info);
    fclose(fp);

    nx = width;
    ny = height;

  }

};

void Map2D::load(const char* filename) {

  grid.clear();
  eps = 1;

  std::string dir = directoryOf(filename);

  std::ifstream istr(filename);

  if (!istr.is_open()) { 
    std::cerr << "error opening " << filename << "\n";
    exit(1);
  }

  float xmin, ymin, xsz, ysz, csz;
  std::string stag;
  
  if (!(istr >> stag) || stag != "scene") {
    std::cerr << "expected scene!\n";
    exit(1);
  }

  if (!(istr >> xmin >> ymin >> xsz >> ysz >> csz)) {
    std::cerr << "expected scene dims!\n";
    exit(1);
  }

  std::cout << "got scene " << xmin << " " << ymin << " " << xsz << " " << ysz << " " << csz << "\n";

  vec3f min(xmin, ymin, -0.5*csz);
  vec3f max = min + vec3f(xsz, ysz, 0.5*csz);
  

  grid.resize(min, max, DtGridf::AXIS_Z, csz);

  std::cout << "grid origin " << grid.origin() << " " 
            << "with size " << grid.nx() << "x" << grid.ny() << "x" << grid.nz() 
            << "\n";

  std::string otag;

  while (istr >> otag) {
    if (otag == "obs2d") {
      std::string ofile;
      Image img;
      if (!(istr >> ofile >> img.ox >> img.oy >> img.sx >> img.sy)) {
        std::cout << "couldn't read obs2d stuff!\n";
        exit(1);
      } else {
        std::cout << "got obs2d " << ofile << " " 
                  << img.ox << " " 
                  << img.oy << " " 
                  << img.sx << " " 
                  << img.sy << "\n";

        if (!dir.empty()) {
          ofile = dir + "/" + ofile;
        }
        img.open(ofile.c_str());
        

        for (size_t y=0; y<grid.ny(); ++y) {
          for (size_t x=0; x<grid.nx(); ++x) {
            vec3u sg(x,y,0);
            vec2f vg = grid.cellCenter(sg).trunc();
            if (img.contains(vg)) {
              vec2u si = img.v2p(vg);
              size_t ii = img.sub2ind(si);
              if (img.buf[ii] >= 128) { grid(sg) = -1; }
            }
          }

        }

      }
    } else if (otag == "eps") {
      if (!(istr >> eps)) {
        std::cerr << "error reading eps\n";
      }
      std::cout << "read eps " << eps << "\n";
    } else {
      std::cerr << "bad tag " << otag << "\n";
      exit(1);
    }
  }

  if (!istr.eof()) {
    std::cerr << "unexpected eof!\n";
    exit(1);
  }

  grid.computeDistsFromBinary();

  /*
  quickPNG("foo.png", &(grid(0,0,0)),
           grid.nx(), grid.ny(), grid.nx(), true,
           Gradient::jet());
  */

}

float Map2D::sampleCost(const vec3f& p) const {
  
  float d = grid.sample(p);
  return distToCost(d);

}

float Map2D::distToCost(float d) const {

  if (d < 0) {
    return -d + 0.5*eps;
  } else if (d <= eps) {
    float f = d-eps;
    return f*f*0.5/eps;
  } else {
    return 0;
  }

}

float Map2D::distToCost(float d, float& g) const {

  if (d < 0) {

    g = -1;
    return -d + 0.5*eps;

  } else if (d <= eps) {

    float f = d-eps;
    g = f*0.5/eps;
    return f*f*0.5/eps;

  } else {

    g = 0;
    return 0;

  }

}

float Map2D::sampleCost(const vec3f& p, vec3f& grad) const {

  float d = grid.sample(p, grad);
  float g;

  float c = distToCost(d, g);

  grad *= g;
  return c;

}

void Map2D::rasterize(RasterType type,
                      std::vector<unsigned char>& map,
                      size_t stride) const {

  // get limits
  float vmin = 1e10;
  float vmax = -1e10;

  if (type == RASTER_DISTANCE || type == RASTER_COST) {

    for (size_t y=0; y<grid.ny(); ++y) {
      for (size_t x=0; x<grid.nx(); ++x) {
        float v = grid(x,y,0);
        if (type == RASTER_COST) {
          v = distToCost(v);
        }
        vmin = std::min(vmin, v);
        vmax = std::max(vmax, v);
      }
    }

  }

  if (stride < 4*grid.nx()) {
    stride = 4*grid.nx();
  }

  map.resize( grid.ny() * stride, 0xff );

  unsigned char* rowptr = &(map[0]);

  Gradient grad;

  if (type == RASTER_DISTANCE) {
    grad = Gradient::jet();
  } else if (type == RASTER_COST) {
    grad = Gradient::hot();
  } else {
    grad.stops[0.0f] = vec3f(1.0f);
    grad.stops[1.0f] = vec3f(0.5f);
  }

  for (size_t y=0; y<grid.ny(); ++y) {

    unsigned char* pxptr = rowptr;

    for (size_t x=0; x<grid.nx(); ++x) {

      float v = grid(x,y,0);

      if (type == RASTER_DISTANCE || type == RASTER_COST) {
        if (type == RASTER_COST) {
          v = distToCost(v);
        } 
        v = (v-vmin)/(vmax-vmin);
      } else {
        v = (v < 0.0);
      }

      vec3f c = grad.lookup(v);

      for (int i=0; i<3; ++i) {
        *pxptr++ = std::max(0.0f, std::min(c[2-i], 1.0f))*255;
      }

      *pxptr++ = 0xff;

    }

    rowptr += stride;

  }

}
