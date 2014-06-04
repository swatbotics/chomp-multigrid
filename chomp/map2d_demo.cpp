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
#include <getopt.h>
#include "Chomp.h"

using namespace chomp;

#ifdef MZ_HAVE_CAIRO
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>
#endif

// Utility function to visualize maps
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
      buf.push_back(pxptr[2]);
      buf.push_back(pxptr[1]);
      buf.push_back(pxptr[0]);
      pxptr += 4;
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

  std::cout << "wrote " << filename << "\n";

  return true;

}

//////////////////////////////////////////////////////////////////////
// class to help evaluate collisons for gradients

class Map2DCHelper: public ChompCollisionHelper {
public: 

  enum {
    NUM_CSPACE = 2,
    NUM_WKSPACE = 3, // actually just 2 but this way I can test matrix dims better
    NUM_BODIES = 1,
  };

  const Map2D& map;

  Map2DCHelper(const Map2D& m): 
    ChompCollisionHelper(NUM_CSPACE, NUM_WKSPACE, NUM_BODIES), 
    map(m) {}

  virtual ~Map2DCHelper() {}

  virtual double getCost(const MatX& q, 
                         size_t body_index,
                         MatX& dx_dq,
                         MatX& cgrad) {

    assert( (q.rows() == 2 && q.cols() == 1) ||
            (q.rows() == 1 && q.cols() == 2) );

    dx_dq.conservativeResize(3, 2);
    dx_dq.setZero();

    dx_dq << 1, 0, 0, 1, 0, 0;

    cgrad.conservativeResize(3, 1);

    vec3f g;
    float c = map.sampleCost(vec3f(q(0), q(1), 0.0), g);

    cgrad << g[0], g[1], 0.0;

    return c;

  }


};

//////////////////////////////////////////////////////////////////////
// help generate an an initial trajectory

void generateInitialTraj(int N, 
                         const Map2D& map, 
                         const vec2f& p0, 
                         const vec2f& p1,
                         MatX& xi,
                         MatX& q0,
                         MatX& q1) {
  
  xi.resize(N, 2);
  q0.resize(1, 2);
  q1.resize(1, 2);

  q0 << p0.x(), p0.y();
  q1 << p1.x(), p1.y();

  for (int i=0; i<N; ++i) {
    float u = float(i+1)/(N+1);
    vec2f pi = p0 + u*(p1-p0);
    xi(i,0) = pi.x();
    xi(i,1) = pi.y();
  }

}

//////////////////////////////////////////////////////////////////////

void usage(int status) {
  std::ostream& ostr = status ? std::cerr : std::cout;
  ostr <<
    "usage: map2d_demo OPTIONS map.txt\n"
    "Also, checkout the map2d_tests.sh script!\n"
    "\n"
    "OPTIONS:\n"
    "\n"
    "  -c, --coords             Set start, goal (x0,y0,x1,y1)\n"
    "  -n, --num                Number of steps for trajectory\n"
    "  -a, --alpha              Overall step size for CHOMP\n"
    "  -g, --gamma              Step size for collisions\n"
    "  -m, --max-iter           Set maximum iterations\n"
    "  -e, --error-tol          Relative error tolerance\n"
    "  -o, --objective          Quantity to minimize (vel|accel)\n"
    "  -p, --pdf                Output PDF's every I iterations (0=only init/final)\n"
    "      --help               See this message.\n";
  exit(status);
}

//////////////////////////////////////////////////////////////////////
// helper class to visualize stuff

#ifdef MZ_HAVE_CAIRO

class PdfEmitter: public DebugChompObserver {
public:

  const Map2D& map;
  const MatX& xi_init;
  
  int dump_every;
  int count;

  const char* filename;
  
  cairo_surface_t* surface;
  cairo_surface_t* image;
  cairo_t* cr;
  int width, height;
  float mscl;

  PdfEmitter(const Map2D& m, const MatX& x, int de, const char* f):
    map(m), xi_init(x), dump_every(de), count(0), filename(f)
  {
    
    Box3f bbox = map.grid.bbox();
    vec3f dims = bbox.p1 - bbox.p0;

    mscl = (400) / std::max(dims.x(), dims.y());
    
    width = int(mscl * dims.x());
    height = int(mscl * dims.y());
    std::cout << "image will be " << width << "x" << height << "\n";
    
    surface = cairo_pdf_surface_create(filename, width, height);
    
    cr = cairo_create(surface);
    
    size_t stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, map.grid.nx());

    std::vector<unsigned char> buf;
    map.rasterize(Map2D::RASTER_OCCUPANCY, buf, stride);
    
    image = cairo_image_surface_create_for_data(&(buf[0]), 
                                                CAIRO_FORMAT_RGB24,
                                                map.grid.nx(), map.grid.ny(),
                                                stride);
  }

  virtual ~PdfEmitter() {

    cairo_surface_destroy(surface);
    cairo_surface_destroy(image);
    cairo_destroy(cr);

    std::cout << "wrote " << filename << "\n\n";

  }

  virtual int notify(const Chomp& chomper, 
                     ChompEventType event,
                     size_t iter,
                     double curObjective,
                     double lastObjective,
                     double hmag) {
 
    DebugChompObserver::notify(chomper, event, iter, 
                               curObjective, lastObjective, hmag);

    bool dump_pdf = (event == CHOMP_FINISH ||
                     dump_every == 0 || 
                     (iter % dump_every == 0));

    if (!dump_pdf) {
      return 0;
    }

    if (count++) {
      cairo_show_page(cr);
    }

    float cs = map.grid.cellSize();
    Box3f bbox = map.grid.bbox();

    cairo_identity_matrix(cr);
    cairo_translate(cr, 0, height);
    cairo_scale(cr, mscl, -mscl);
    cairo_translate(cr, -bbox.p0.x(), -bbox.p0.y());

    cairo_save(cr);
    cairo_scale(cr, cs, cs);
    cairo_set_source_surface(cr, image, bbox.p0.x()/cs, bbox.p0.y()/cs);
    cairo_paint(cr);
    cairo_restore(cr);
  
    cairo_set_line_width(cr, 1.0*cs);
    cairo_set_source_rgb(cr, 0.0, 0.1, 0.5);

    int N = chomper.N;

    for (int i=0; i<N; ++i) {
      vec3f pi(xi_init(i,0), xi_init(i,1), 0.0);
      cairo_arc(cr, pi.x(), pi.y(), 1.0*cs, 0.0, 2*M_PI);
      cairo_fill(cr);
    }

    const MatX& q0 = chomper.q0;
    const MatX& q1 = chomper.q1;

    cairo_set_source_rgb(cr, 0.5, 0.0, 1.0);
    cairo_arc(cr, q0(0), q0(1), 4*cs, 0.0, 2*M_PI);
    cairo_fill(cr);
    cairo_arc(cr, q1(0), q1(1), 4*cs, 0.0, 2*M_PI);
    cairo_fill(cr);

    for (int i=0; i<N; ++i) {
      vec3f pi(chomper.xi(i,0), chomper.xi(i,1), 0.0);
      cairo_arc(cr, pi.x(), pi.y(), 2*cs, 0.0, 2*M_PI);
      cairo_fill(cr);
    }

    return 0;

 }

  
};

#endif

//////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

  const struct option long_options[] = {
    { "coords",            required_argument, 0, 'c' },
    { "num",               required_argument, 0, 'n' },
    { "alpha",             required_argument, 0, 'a' },
    { "gamma",             required_argument, 0, 'g' },
    { "error-tol",         required_argument, 0, 'e' },
    { "max-iter",          required_argument, 0, 'm' },
    { "objective",         required_argument, 0, 'o' },
    { "pdf",               required_argument, 0, 'p' },
    { "help",              no_argument,       0, 'h' },
    { 0,                   0,                 0,  0  }
  };

  const char* short_options = "c:n:a:g:e:m:o:p:h";
  int opt, option_index;

  int N = 127;
  double gamma = 0.5;
  double alpha = 0.02;
  double errorTol = 1e-6;
  size_t max_iter = 500;
  ChompObjectiveType otype = MINIMIZE_VELOCITY;
  int pdf = -1;
  float x0=0, y0=0, x1=0, y1=0;

  while ( (opt = getopt_long(argc, argv, short_options, 
                             long_options, &option_index) ) != -1 ) {

    switch (opt) {
    case 'c': 
      if (sscanf(optarg, "%f,%f,%f,%f", &x0, &y0, &x1, &y1) != 4) {
        std::cerr << "error parsing coords!\n\n";
        usage(1);
      } 
      break;
    case 'n':
      N = atoi(optarg);
      break;
    case 'a':
      alpha = atof(optarg);
      break;
    case 'g':
      gamma = atof(optarg);
      break;
    case 'e':
      errorTol = atof(optarg);
      break;
    case 'm':
      max_iter = atoi(optarg);
      break;
    case 'o':
      if (!strcasecmp(optarg, "vel")) {
        otype = MINIMIZE_VELOCITY;
      } else if (!strcasecmp(optarg, "accel")) {
        otype = MINIMIZE_ACCELERATION;
      } else {
        std::cerr << "error parsing opt. type: " << optarg << "\n\n";
        usage(1);
      }
      break;
    case 'p':
      pdf = atoi(optarg);
      break;
    case 'h':
      usage(0);
      break;
    default:
      usage(1);
      break;
    }

  }

  if (argc < 2) {
    usage(1);
  }

  Map2D map;

  map.load(argv[argc-1]);

  std::vector<unsigned char> buf;

  map.rasterize(Map2D::RASTER_DISTANCE, buf, 0);
  savePNG_RGB24("dist.png", map.grid.nx(), map.grid.ny(), 
                map.grid.nx()*4, &buf[0], true);

  map.rasterize(Map2D::RASTER_COST, buf, 0);
  savePNG_RGB24("cost.png", map.grid.nx(), map.grid.ny(), 
                map.grid.nx()*4, &buf[0], true);

  map.rasterize(Map2D::RASTER_OCCUPANCY, buf, 0);
  savePNG_RGB24("occupancy.png", map.grid.nx(), map.grid.ny(), 
                map.grid.nx()*4, &buf[0], true);

  MatX q0, q1, xi;

  Map2DCHelper mhelper(map);
  ChompCollGradHelper cghelper(&mhelper, gamma);

  vec2f p0(x0, y0), p1(x1, y1);
  if (p0.x() == p0.y() && p0 == p1) {
    p0 = map.grid.bbox().p0.trunc();
    p1 = map.grid.bbox().p1.trunc();
  }

  generateInitialTraj(N, map, p0, p1, xi, q0, q1);

  Chomp chomper(NULL, xi, q0, q1, N, alpha, errorTol, max_iter);
  chomper.objective_type = otype;
  chomper.ghelper = &cghelper;

  DebugChompObserver dobs;
  chomper.observer = &dobs;

#ifdef MZ_HAVE_CAIRO

  PdfEmitter* pe = NULL;

  if (pdf >= 0) {
    char buf[1024];
    sprintf(buf, "map2d_n%d_g%f_a%f_e%f_m%d_o%s_%f,%f,%f,%f.pdf",
            N, gamma, alpha, errorTol, (int)max_iter,
            otype == MINIMIZE_VELOCITY ? "vel" : "accel",
            p0.x(), p0.y(), p1.x(), p1.y());
    pe = new PdfEmitter(map, xi, pdf, buf);
    chomper.observer = pe;
  }

#endif
  
  chomper.solve(true, false);

#ifdef MZ_HAVE_CAIRO
  delete pe;
#endif

  return 0;

}
