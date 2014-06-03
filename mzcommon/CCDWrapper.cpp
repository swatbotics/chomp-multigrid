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

#include "CCDWrapper.h"
#include <assert.h>
#include <string.h>
#include <sstream>

namespace ccdw {

  static inline vec3& ccd2mz(ccd_vec3_t* c) { 
    return *((vec3*)c); 
  }

  static inline const vec3& ccd2mz(const ccd_vec3_t* c) { 
    return *((const vec3*)c); 
  }

  static inline ccd_vec3_t* mz2ccd(vec3& v) { 
    return (ccd_vec3_t*)(&v); 
  }

  /*
  static inline const ccd_vec3_t* mz2ccd(const vec3& v) { 
    return (const ccd_vec3_t*)(&v); 
  }

  static inline ccd_vec3_t* mz2ccd(vec3* v) { 
    return (ccd_vec3_t*)(v); 
  }

  static inline const ccd_vec3_t* mz2ccd(const vec3* v) { 
    return (const ccd_vec3_t*)(v); 
  }
  */

  DrawHelper::DrawHelper(): 
    quadric(0),
    slices(32),
    sstacks(16),
    cstacks(1)
  {}

  DrawHelper::~DrawHelper() {
    if (quadric) { 
      gluDeleteQuadric(quadric);
    }
  }

  GLUquadric* DrawHelper::getQuadric() {
    if (!quadric) {
      quadric = gluNewQuadric();
    }
    return quadric;
  }

  //////////////////////////////////////////////////////////////////////

  void support(const void* obj,
               const ccd_vec3_t* dir,
               ccd_vec3_t* vec) {

    const Convex* c = (const Convex*)obj;
    assert( c );
    c->support( ccd2mz(dir), ccd2mz(vec) );

  }
  
  void center(const void* obj,
              ccd_vec3_t* center) {

    const Convex* c = (const Convex*)obj;
    assert( c );
    
    c->center( ccd2mz(center) );

  }

  //////////////////////////////////////////////////////////////////////

  Convex::~Convex() {}

  std::string Convex::description() const {
    std::ostringstream ostr;
    describe(ostr);
    return ostr.str();
  }
    

  void Convex::center(vec3& c) const {
    c = vec3(0);
  }

  bool Convex::isDilated() const {
    return false;
  }

  //////////////////////////////////////////////////////////////////////

  Box::Box(): extents(1) {}

  Box::Box(const vec3& e): extents(e) {}

  Box::~Box() {}

  void Box::describe(std::ostream& ostr) const {
    ostr << "Box(" << extents << ")";
  }

  /*
  void Box::closest(vec3& v) const {

    bool inside = true;

    vec3 absdists(0), vout(0);
    int minaxis = 0;
    
    for (int axis=0; axis<3; ++axis) {
      if (v[axis] < -0.5*extents[axis]) {
        inside = false;
        v[axis] = -0.5*extents[axis];
      } else if (v[axis] > 0.5*extents[axis]) {
        inside = false;
        v[axis] = 0.5*extents[axis];
      } else {
        if (v[axis] < 0) {
          absdists[axis] = 0.5*extents[axis] + v[axis];
          vout[axis] = -0.5*extents[axis];
        } else {
          absdists[axis] = 0.5*extents[axis] - v[axis];
          vout[axis] = 0.5*extents[axis];
        }
        if (absdists[axis] < absdists[minaxis]) {
          minaxis = axis;
        }
      }
    }

    if (inside) { 
      v[minaxis] = vout[minaxis];
    }

  };
  */

  ccd_real_t Box::maxDist() const {
    return 0.5*extents.norm();
  }
  
  void Box::support(const vec3& dir, vec3& s) const {

    s = vec3( ccdSign(dir[0]) * extents[0] * 0.5,
              ccdSign(dir[1]) * extents[1] * 0.5,
              ccdSign(dir[2]) * extents[2] * 0.5 );

  }

  /*
  static inline vec3 cwiseproduct(const vec3& a, const vec3& b) {
    return vec3(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
  }
  */


  void Box::render(DrawHelper& h, ccd_real_t dilation) const {


    Box3_t<ccd_real_t> box(-0.5*extents, 0.5*extents);

    if (dilation) {

      glstuff::draw_round_box(box, dilation, h.slices, h.sstacks);
      
    } else {

      glstuff::draw_box( box );

    }

  }
  
  //////////////////////////////////////////////////////////////////////


  Point::Point() {}

  Point::~Point() {}

  void Point::describe(std::ostream& ostr) const {
    ostr << "Point()";
  }

  /*
  void Point::closest(vec3& v) const {
    v = vec3(0);
  }
  */

  ccd_real_t Point::maxDist() const {
    return 0;
  }

  void Point::support(const vec3& dir, vec3& s) const {

    s = vec3(0);

  }

  void Point::render(DrawHelper& h, ccd_real_t dilation) const {

    if (dilation) {
      gluSphere(h.getQuadric(), dilation, h.slices, h.sstacks);
    } else {
      glPushAttrib(GL_ENABLE_BIT);
      glDisable(GL_LIGHTING);
      glBegin(GL_POINTS);
      glVertex3f(0,0,0);
      glEnd();
      glPopAttrib();
    }

  }

  //////////////////////////////////////////////////////////////////////

  Line::Line(): length(1) {}

  Line::Line(ccd_real_t l): length(l) {}
  
  Line::~Line() {}

  void Line::describe(std::ostream& ostr) const {
    ostr << "Line(" << length << ")";
  }

  /*
  void Line::closest(vec3& v) const {

    if (v[2] < -0.5*length) {
      v[2] = -0.5*length;
    } else if (v[2] > 0.5*length) {
      v[2] = 0.5*length;
    }
    v[0] = v[1] = 0;

  }
  */

  ccd_real_t Line::maxDist() const {
    return 0.5*length;
  }

  void Line::support(const vec3& dir, vec3& s) const {
    s = vec3( 0, 0, ccdSign(dir[2]) * length * 0.5 );
  }

  void Line::render(DrawHelper& h, ccd_real_t dilation) const {

    vec3 p0(0, 0, -0.5*length);
    vec3 p1(0, 0,  0.5*length);

    if (dilation) {
      glstuff::draw_capsule(h.getQuadric(), p0, p1, dilation, 
                            h.slices, h.sstacks, h.cstacks);
    } else {
      glPushAttrib(GL_ENABLE_BIT);
      glDisable(GL_LIGHTING);
      glBegin(GL_LINES);
      glstuff::vertex(p0);
      glstuff::vertex(p1);
      glEnd();
      glPopAttrib();
    }

  }
    
  //////////////////////////////////////////////////////////////////////


  Cylinder::Cylinder(): length(1), radius(0.5) {}

  Cylinder::Cylinder(ccd_real_t l, ccd_real_t r): length(l), radius(r) {}
  
  Cylinder::~Cylinder() {}

  void Cylinder::describe(std::ostream& ostr) const {
    ostr << "Cylinder(" << length << ", " << radius << ")";
  }

  ccd_real_t Cylinder::maxDist() const {
    return sqrt(0.25*length*length + radius*radius);
  }

  void Cylinder::support(const vec3& dir, vec3& s) const {
    vec2_t<ccd_real_t> p = dir.trunc();
    s = vec3( p * (radius / p.norm()), ccdSign(dir[2]) * length * 0.5 );
  }

  void Cylinder::render(DrawHelper& h, ccd_real_t dilation) const {

    vec3 p0(0, 0, -0.5*length);
    vec3 p1(0, 0,  0.5*length);

    if (dilation) {
      // TODO: deal
    } else {
      glstuff::draw_cylinder(h.getQuadric(), p0, p1, radius, 
                             h.slices, h.cstacks);
    }

  }

  //////////////////////////////////////////////////////////////////////

  DilatedConvex::DilatedConvex(): child(0), dilation(0) {}

  DilatedConvex::DilatedConvex(const Convex* c, ccd_real_t r): 
    child(c), dilation(r) 
  {
    assert(child);
    assert(dilation > 0);
  }

  DilatedConvex::~DilatedConvex() {}

  void DilatedConvex::describe(std::ostream& ostr) const {
    ostr << "DilatedConvex(";
    child->describe(ostr);
    ostr << ", " << dilation << ")";
  }

  /*
  void DilatedConvex::closest(vec3& v) const {
    assert(child);
    assert(dilation > 0);
    vec3 vc;
    child->closest(vc);
    vec3 dir = v - vc;
    v = vc + dir * (dilation/dir.norm());
  }
  */


  ccd_real_t DilatedConvex::maxDist() const {
    assert(child);
    assert(dilation > 0);
    return child->maxDist() + dilation;
  }

  void DilatedConvex::support(const vec3& dir, vec3& s) const {
    assert(child);
    assert(dilation > 0);
    child->support(dir, s);
    s += dir * (dilation / dir.norm());
  }

  void DilatedConvex::render(DrawHelper& h, float extra) const {
    assert(child);
    child->render(h, dilation+extra);
  }

  bool DilatedConvex::isDilated() const {
    return true;
  }

  //////////////////////////////////////////////////////////////////////

  TransformedConvex::TransformedConvex(): child(0) {}
    
  TransformedConvex::TransformedConvex(const Convex* c, const Transform3& x):
    child(c), xform(x)
  {
    assert(child);
  }

  TransformedConvex::~TransformedConvex() {}

  void TransformedConvex::describe(std::ostream& ostr) const {
    ostr << "TransformedConvex(";
    child->describe(ostr);
    ostr << ", " << xform << ")";
  }

  /*
  void TransformedConvex::closest(vec3& v) const {
    assert(child);
    v = xform.transformInv(v);
    child->closest(v);
    v = xform.transformFwd(v);
  }
  */

  ccd_real_t TransformedConvex::maxDist() const {
    assert(child);
    return child->maxDist();
  }

  void TransformedConvex::support(const vec3& dir, vec3& s) const {
    assert(child);
    child->support( xform.rotInv()*dir, s );
    s = xform * s;
  }

  void TransformedConvex::center(vec3& c) const {
    assert(child);
    child->center(c);
    c = xform * c;
  }

  void TransformedConvex::render(DrawHelper& h, float radius) const {
    assert(child);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glstuff::mult_transform(xform);
    child->render(h, radius);
    glPopMatrix();
  }

  bool TransformedConvex::isDilated() const {
    return child->isDilated();
  }

  //////////////////////////////////////////////////////////////////////

  DilatedConvex* sphere(ccd_real_t radius) {
    return new DilatedConvex(new Point(), radius);
  }

  DilatedConvex* capsule(ccd_real_t length, ccd_real_t radius) {
    return new DilatedConvex(new Line(length), radius);
  }

  TransformedConvex* capsule(const vec3& p0,
                             const vec3& p1, 
                             ccd_real_t radius) {

    vec3 z = p1-p0;
    ccd_real_t length = z.norm();

    vec3 mid = 0.5*(p1+p0);

    return transform(capsule(length, radius),
                     Transform3(quat::fromOneVector(z, 2), mid));
    
  }

  TransformedConvex* cylinder(const vec3& p0,
                              const vec3& p1, 
                              ccd_real_t radius) {
    
    vec3 z = p1-p0;
    ccd_real_t length = z.norm();

    vec3 mid = 0.5*(p1+p0);
    
    return transform(new Cylinder(length, radius),
                     Transform3(quat::fromOneVector(z, 2), mid));
    
  }

  TransformedConvex* transform(const Convex* c, const Transform3& xform) {

    const TransformedConvex* t = dynamic_cast<const TransformedConvex*>(c);

    if (t) {
      return transform(t->child, xform * t->xform);
    } else {
      return new TransformedConvex(c, xform);
    }

  }

  Convex* dilate(const Convex* c, ccd_real_t radius) {

    const TransformedConvex* t = dynamic_cast<const TransformedConvex*>(c);
    
    if (t) {
      return transform(dilate(t->child, radius), t->xform);
    } else {
      const DilatedConvex* d = dynamic_cast<const DilatedConvex*>(c);
      if (d) {
        return dilate(d->child, radius + d->dilation);
      } else {
        return new DilatedConvex(c, radius);
      }
    }

  }

  //////////////////////////////////////////////////////////////////////

  Report::Report(const Convex* a, const Convex* b): 
    c1(a), c2(b), flags(0), distance(0), algorithm(NUM_ALGORITHMS) {}

  Checker::Checker() {
    CCD_INIT(&ccd);
    ccd.support1 = ccdw::support;
    ccd.support2 = ccdw::support;
    ccd.center1 = ccdw::center;
    ccd.center2 = ccdw::center;
    ccd.max_iterations = 1000;
    algorithm = ALGORITHM_MPR;
  }

  bool Checker::intersect(const Convex* c1, const Convex* c2, 
                           ccd_real_t dmin) const {

    return query(QUERY_INTERSECT, 0, c1, c2, dmin);

  }
  
  bool Checker::separate(Report& report,
                          const Convex* c1, const Convex* c2, 
                          ccd_real_t dmin) const {

    return query(QUERY_SEPARATION, &report, c1, c2, dmin);

  }
  
  bool Checker::penetration(Report& report,
                             const Convex* c1, const Convex* c2,
                             ccd_real_t dmin) const {
    
    return query(QUERY_PENETRATION, &report, c1, c2, dmin);

  }
  
  bool Checker::query(QueryType qtype, Report* report,
                         const Convex* orig_c1, const Convex* orig_c2,
                         ccd_real_t dmin) const {

    assert( orig_c1 );
    assert( orig_c2 );
    assert( dmin >= 0 );

    if (report) {
      *report = Report(orig_c1, orig_c2);
    }

    // do bounding sphere test
    vec3 ctr1, ctr2;
    orig_c1->center(ctr1);
    orig_c2->center(ctr2);

    ccd_real_t r1 = orig_c1->maxDist();
    ccd_real_t r2 = orig_c2->maxDist();

    ccd_real_t d = dmin + r1 + r2;

    if ((ctr2-ctr1).norm2()  > d*d) {
      return false;
    }

    const Convex* c1 = orig_c1;
    const Convex* c2 = orig_c2;

    // do relative transform if necessary
    TransformedConvex* transformed = 0;

    const TransformedConvex* t1 = dynamic_cast<const TransformedConvex*>(orig_c1);
    const TransformedConvex* t2 = dynamic_cast<const TransformedConvex*>(orig_c2);

    if (t1 && t2) {
      transformed = transform(t2->child, t1->xform.inverse() * t2->xform);
      c1 = t1->child;
      c2 = transformed;
    }



    AlgorithmType actual_algorithm = algorithm;
    if (qtype == QUERY_SEPARATION) {
      actual_algorithm = ALGORITHM_GJK;
    }

    Convex* dilated = 0;

    if (dmin && qtype != QUERY_INTERSECT) {
      ccd_real_t actual = ccdGJKDist(c1, c2, &ccd);
      if (actual > 0) {
        if (actual < dmin + ccd.dist_tolerance) {
          if (qtype == QUERY_PENETRATION) {
            actual_algorithm = ALGORITHM_MPR;
          }
          ccd_real_t alg_tolerance;
          if (actual_algorithm == ALGORITHM_GJK) {
            alg_tolerance = 100*ccd.epa_tolerance;
          } else {
            alg_tolerance = 10*ccd.mpr_tolerance;
          }
          dmin = actual + alg_tolerance;
        } else {
          // separated by more than dmin.
          // disable query, it'll return false anyways
          qtype = NUM_QUERY_TYPES;
        }
      } else {
        // they are intersecting
        dmin = 0;
      }
    }

    if (dmin) {
      if (c1->isDilated()) {
        c1 = dilated = dilate(c1, dmin);
      } else {
        c2 = dilated = dilate(c2, dmin);
      }
    }


    bool intersect = false;

    switch (qtype) {
    case QUERY_INTERSECT:
      if (actual_algorithm == ALGORITHM_GJK) {
        intersect = ccdGJKIntersect(c1, c2, &ccd);
      } else {
        intersect = ccdMPRIntersect(c1, c2, &ccd);
      }
      if (report && intersect) { 
        report->flags = INTERSECT;
      }
      break;
    case QUERY_SEPARATION: {
      vec3 sep;
      assert(actual_algorithm == ALGORITHM_GJK);
      intersect = (ccdGJKSeparate(c1, c2, &ccd, mz2ccd(sep)) == 0);
      if (report && intersect) {

        report->flags = INTERSECT | HAVE_SEPARATION;

        vec3 s1, s2;
        c1->support( sep, s1);
        c2->support(-sep, s2);

        if (c2->isDilated()) {
          if (c1->isDilated()) {
            report->pos1 = 0.5*(s1+s2+sep);
            report->pos2 = 0.5*(s1+s2-sep);
          } else {
            report->pos1 = s2 + sep;
            report->pos2 = s2;
          }
        } else {
          report->pos1 = s1;
          report->pos2 = s1 - sep;
        }


        report->distance = sep.norm();
        report->direction = sep / report->distance;

      }
      break;
    }
    case QUERY_PENETRATION: {
      ccd_real_t depth;
      vec3 dir, pos;
      if (actual_algorithm == ALGORITHM_GJK) {
        intersect = (ccdGJKPenetration(c1, c2, &ccd, &depth, 
                                       mz2ccd(dir), mz2ccd(pos)) == 0);
      } else {
        intersect = (ccdMPRPenetration(c1, c2, &ccd, &depth, 
                                       mz2ccd(dir), mz2ccd(pos)) == 0);
      }
      if (report && intersect) {
        report->flags = INTERSECT | HAVE_SEPARATION | HAVE_POSITION;
        report->distance = depth;
        report->direction = dir;
        report->pos1 = pos + (0.5*depth)*dir;
        report->pos2 = pos - (0.5*depth)*dir;
      }
      break;
    }
    default:
      break;
    };

    if (transformed && report) {
      report->pos1 = t1->xform * report->pos1;
      report->pos2 = t1->xform * report->pos2;
      report->direction = t1->xform.rotFwd() * report->direction;
    }

    if (dilated && report) {
      if (dilated == c1) {
        report->pos1 -= dmin*report->direction;
      } else {
        report->pos2 += dmin*report->direction;
      }
      report->distance -= dmin;
    }

    if (report) {
      report->algorithm = actual_algorithm;
    }

    delete dilated;
    delete transformed;
    
    return intersect;

  }

  //////////////////////////////////////////////////////////////////////

  static inline ccd_real_t frac(size_t i, size_t n) {
    return ccd_real_t(i) / ccd_real_t(n-1);
  }

  void cubePoints(size_t n, std::vector<vec3>& points) {

    points.clear();

    
    for (int side=0; side<6; ++side) {

      int ax0 = (side + 0) % 3;
      int ax1 = (side + 1) % 3;
      int ax2 = (side + 2) % 3;

      int sign = (side / 3) ? -1 : 1;
      
      vec3 p;
      
      p[ax0] = sign;

      for (size_t u=0; u<n; ++u) {
        p[ax1] = sign*(frac(u, n)*2-1);
        for (size_t v=0; v<n; ++v) {
          p[ax2] = sign*(frac(v, n)*2-1);
          points.push_back(p);
        }
      }

    }

  }


}
