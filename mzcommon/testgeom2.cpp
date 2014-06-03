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

#include <mzcommon/mersenne.h>
#include <mzcommon/TimeUtil.h>
#include <mzcommon/Affine2.h>
#include <mzcommon/Geom2.h>
#include <algorithm>
#include <assert.h>

#ifdef MZ_HAVE_CAIRO
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>
#include <stdio.h>
#endif


typedef double real;
typedef vec2_t<real> vec2;
typedef vec3_t<real> vec3;
typedef Box2_t<real> Box2;
typedef Affine2_t<real> Affine2;
typedef Geom2_t<real> Geom2;
typedef vec2_t<real> Point;
typedef Polygon2_t<Point> Polygon2;
typedef Polygon2::PointArray PointArray;

#ifdef MZ_HAVE_CAIRO

class CairoWrapper {
public:

  cairo_t* cr;
  Affine2 tx;

  CairoWrapper(cairo_t* c, const Affine2& a): cr(c), tx(a) {}

  void moveTo(vec2 p) {
    p = tx*p;
    cairo_move_to(cr, p.x(), p.y());
  }

  void lineTo(vec2 p) {
    p = tx*p;
    cairo_line_to(cr, p.x(), p.y());
  }

  void circle(vec2 p, real r) {
    p = tx*p;
    cairo_arc(cr, p.x(), p.y(), r, 0, 2*M_PI);
  }

  void makePath(const PointArray& p) {
    if (p.empty()) { return; }
    // draw convex hull
    moveTo(p[0]);
    for (size_t i=1; i<p.size(); ++i) { lineTo(p[i]); }
    cairo_close_path(cr);
  }

  void fillPreserve(const vec3& c) {
    cairo_set_source_rgb(cr, c[0], c[1], c[2]);
    cairo_fill_preserve(cr);
  }

  void fill(const vec3& c) {
    cairo_set_source_rgb(cr, c[0], c[1], c[2]);
    cairo_fill(cr);
  }

  void stroke(const vec3& c) {
    cairo_set_source_rgb(cr, c[0], c[1], c[2]);
    cairo_stroke(cr);
  }


  void drawPolygon(const PointArray& p,
                   bool doFill, const vec3& fillColor,
                   bool doStroke, const vec3& strokeColor) {
    if (p.size() < 2) { return; }
    cairo_new_path(cr);
    makePath(p);
    if (doFill) { fillPreserve(fillColor); }
    if (doStroke) { stroke(strokeColor); }
  }

  void drawDots(const PointArray& p, real rad, 
                const vec3& color) {
    cairo_new_path(cr);
    for (size_t i=0; i<p.size(); ++i) {
      cairo_new_sub_path(cr);
      circle(p[i], rad);
    }
    fill(color);
  }

};

void drawConvex(const PointArray& orig, const PointArray& convex, size_t iter) {
  
  // get the bounding box of the polygon
  Box2 b = Polygon2::computeBBox(orig);
  vec2 dims = b.p1 - b.p0;
  vec2 center = b.p0 + 0.5*dims;

  real psz = 1.2 * ( std::max(dims[0], dims[1]) + 0.05 );
  real offs = psz * 0.05;
  int dsz = 288;
  real scl = dsz/psz;

  cairo_surface_t *surface;
  cairo_t *cr;

  Affine2 a = 
    Affine2::translation(vec2(dsz/2)) * 
    Affine2::scale(scl, scl) *
    Affine2::translation(-center);

  char filename[1024];
  snprintf(filename, 1024, "convex_%02d_%06d.pdf", int(orig.size()), int(iter));

  surface = cairo_pdf_surface_create(filename, dsz, dsz);
  cr = cairo_create(surface);
  cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);

  CairoWrapper w(cr, a);

  //Polygon offset;
  //offsetConvexPolygon(convex, offs, offset);

  PointArray cpts, opts, ipts;
  cpts.insert(cpts.end(), convex.begin(), convex.end());

  Polygon2::offsetConvexPolygon(cpts, offs, opts);
  Polygon2::offsetConvexPolygon(cpts, -offs, ipts);

  vec3 red(1.0, 0.7, 0.7);
  w.drawPolygon(opts, true, red, true, 0.5*red);
  
  vec3 blue(0.7, 0.7, 1.0);
  w.drawPolygon(convex, true, blue, true, 0.5*blue);

  vec3 green(0.7, 1.0, 0.7);
  w.drawPolygon(ipts, true, green, true, 0.5*green);

  vec2 ctr;
  Polygon2::area(convex, &ctr);

  cairo_new_path(cr);
  cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
  w.circle(ctr, 6);
  cairo_fill(cr);

  w.drawDots(orig, 2, vec3(0));
  
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
  
  

}

#endif

void generateRandomPoints(size_t size, PointArray& p) {

  p.clear();

  for (size_t i=0; i<size; ++i) {
    p.push_back(vec2( mt_genrand_real1(),
                      mt_genrand_real1() ));
  }

}

void generateDegenerateSet(size_t size, PointArray& p) {

  p.clear();

  size_t nbasis = mt_genrand_int32() % 4 + 1; // between 1 and 4 basis vectors
  PointArray basis;

  for (size_t j=0; j<nbasis; ++j) {
    basis.push_back(vec2(mt_genrand_int32() % 3,
                         mt_genrand_int32() % 3));
  }

  for (size_t i=0; i<size; ++i) {
    p.push_back(vec2(0));
    for (size_t j=0; j<nbasis; ++j) {
      int c = mt_genrand_int32() % 4;
      p.back() += c * basis[j];
    }
  }

}

bool isInside(const vec2& p, const PointArray& convex) {

  for (size_t i=0; i<convex.size(); ++i) {

    size_t j=(i+1)%convex.size();

    vec2 pi = convex[i] - p;
    vec2 pj = convex[j] - p;
    real c = vec2::cross(pi, pj);

    if (c < 0) { 
      return false; 
    }

  }

  return true;

}

std::ostream& operator<<(std::ostream& ostr, const PointArray& p) {
  ostr << "[" << p.size() << "] = {\n";
  for (size_t i=0; i<p.size(); ++i) {
    ostr << "  vec2" << p[i];
    if (i+1 == p.size()) {
      ostr << "\n";
    } else {
      ostr << ",\n";
    }
  }
  ostr << "}";
  return ostr;
}

void testConvex(const PointArray& orig, 
                const PointArray& convex, 
                bool draw, size_t iter) {

#ifdef MZ_HAVE_CAIRO
  if (draw && orig.size() > 2) {
    drawConvex(orig, convex, iter);
  }
#endif

  assert( convex.size() <= orig.size() );

  assert( Polygon2::isConvexCCW( convex ) );
  
  for (size_t i=0; i<orig.size(); ++i) {
    if (!isInside(orig[i], convex)) {
      std::cout << "FAILURE!\n";
      std::cout << "orig" << orig << "\n";
      std::cout << "convex" << convex << "\n";
      std::cout << "bad point is " << orig[i];
      exit(1);
    }
  }

}

int main(int argc, char** argv) {

  PointArray orig, convex;
  
  int niter = 10000;

#ifndef MZ_HAVE_CAIRO
  const bool draw = false;
#else
  bool draw = (argc > 1 && std::string(argv[1]) == "--draw");
  if (draw) { 
    niter = 10; 
  } else {
    std::cout << "you can run this program with the --draw "
              << "argument to emit PDF examples\n";
  }
#endif

  mt_init_genrand(0xDEADBEEF);
    
  for (size_t size=1; size<=32; ++size) {

    TimeStamp start = TimeStamp::now();

    for (int iter=0; iter<niter; ++iter) {
      if (iter < niter/2) {
        generateRandomPoints(size, orig);
      } else {
        generateDegenerateSet(size, orig);
      }
      assert(orig.size() == size);
      Polygon2::convexHull(orig, convex);
      testConvex(orig, convex, draw, iter);
    }
    
    TimeStamp end = TimeStamp::now();
    double elapsed = (end - start).toDouble();
    std::cout << "generated and tested " << niter << " convex hulls "
              << "of size " << size << " in " << elapsed << "s ("
              << (elapsed/niter) << "s per)\n";
      

  }

  return 0;

}
