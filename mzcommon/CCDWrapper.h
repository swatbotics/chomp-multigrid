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

#ifndef _CCDWRAPPER_H_
#define _CCDWRAPPER_H_

#include <ccd/ccd.h>
#include "Transform3.h"
#include <vector>
#include "glstuff.h"
#include <string>

namespace ccdw {

  typedef Transform3_t<ccd_real_t> Transform3;
  typedef vec3_t<ccd_real_t> vec3;
  typedef quat_t<ccd_real_t> quat;

  void support(const void* obj,
               const ccd_vec3_t* dir,
               ccd_vec3_t* vec);
    
  void center(const void* obj,
              ccd_vec3_t* center);

  class DrawHelper {
  public:
    GLUquadric* quadric;
    int slices;
    int sstacks;
    int cstacks;
    DrawHelper();
    ~DrawHelper();
    GLUquadric* getQuadric();
  };
  
  //////////////////////////////////////////////////////////////////////

  class Convex {
  public:

    virtual ~Convex();

    std::string description() const;

    virtual void describe(std::ostream& ostr) const =0;

    virtual ccd_real_t maxDist() const =0;

    virtual void support(const vec3& dir, vec3& s) const=0;

    virtual void center(vec3& c) const;

    virtual bool isDilated() const;

    virtual void render(DrawHelper& h, ccd_real_t radius=0) const =0;


  };

  static inline std::ostream& operator<<(std::ostream& ostr, 
                                         const Convex& c) {
    c.describe(ostr);
    return ostr;
  }

  //////////////////////////////////////////////////////////////////////

  class Box: public Convex {
  public:

    vec3 extents;

    Box();
    Box(const vec3& extents);

    virtual ~Box();

    virtual void describe(std::ostream& ostr) const;
    virtual ccd_real_t maxDist() const;
    virtual void support(const vec3& dir, vec3& s) const;
    virtual void render(DrawHelper& h, ccd_real_t radius=0) const;

  };

  //////////////////////////////////////////////////////////////////////

  class Point: public Convex {
  public:

    Point();

    virtual ~Point();

    virtual void describe(std::ostream& ostr) const;
    virtual ccd_real_t maxDist() const;
    virtual void support(const vec3& dir, vec3& s) const;
    virtual void render(DrawHelper& h, ccd_real_t radius=0) const;

  };

  //////////////////////////////////////////////////////////////////////

  class Line: public Convex {
  public:

    ccd_real_t length;

    Line();
    Line(ccd_real_t length);

    virtual ~Line();

    virtual void describe(std::ostream& ostr) const;
    virtual ccd_real_t maxDist() const;
    virtual void support(const vec3& dir, vec3& s) const;
    virtual void render(DrawHelper& h, ccd_real_t radius=0) const;

  };

  //////////////////////////////////////////////////////////////////////

  class Cylinder: public Convex {
  public:

    ccd_real_t length;
    ccd_real_t radius;

    Cylinder();
    Cylinder(ccd_real_t length, ccd_real_t radius);

    virtual ~Cylinder();

    virtual void describe(std::ostream& ostr) const;
    virtual ccd_real_t maxDist() const;
    virtual void support(const vec3& dir, vec3& s) const;
    virtual void render(DrawHelper& h, ccd_real_t radius=0) const;

  };

  //////////////////////////////////////////////////////////////////////

  class DilatedConvex: public Convex {
  public:

    const Convex* child;
    ccd_real_t dilation;

    DilatedConvex();
    DilatedConvex(const Convex* child, ccd_real_t dilation);

    virtual ~DilatedConvex();

    virtual void describe(std::ostream& ostr) const;
    virtual ccd_real_t maxDist() const;
    virtual void support(const vec3& dir, vec3& s) const;
    virtual void render(DrawHelper& h, ccd_real_t radius=0) const;
    virtual bool isDilated() const;

  };

  //////////////////////////////////////////////////////////////////////

  class TransformedConvex: public Convex {
  public:

    const Convex* child;
    Transform3 xform;

    TransformedConvex();
    TransformedConvex(const Convex* child, const Transform3& xform);

    virtual ~TransformedConvex();

    virtual void describe(std::ostream& ostr) const;
    virtual ccd_real_t maxDist() const;
    virtual void support(const vec3& dir, vec3& s) const;
    virtual void center(vec3& c) const;
    virtual void render(DrawHelper& h, ccd_real_t radius=0) const;
    virtual bool isDilated() const;

  };


  //////////////////////////////////////////////////////////////////////

  DilatedConvex* sphere(ccd_real_t radius);
  DilatedConvex* capsule(ccd_real_t length, ccd_real_t radius);

  TransformedConvex* capsule(const vec3& p0,
                             const vec3& p1, 
                             ccd_real_t radius);

  TransformedConvex* cylinder(const vec3& p0,
                              const vec3& p1, 
                              ccd_real_t radius);

  TransformedConvex* transform(const Convex* c, const Transform3& xform);

  Convex* dilate(const Convex* c, ccd_real_t radius);

  //////////////////////////////////////////////////////////////////////

  enum AlgorithmType {
    ALGORITHM_GJK=0,
    ALGORITHM_MPR,
    NUM_ALGORITHMS
  };

  enum QueryType {
    QUERY_INTERSECT=0,
    QUERY_SEPARATION,
    QUERY_PENETRATION,    
    NUM_QUERY_TYPES,
  };

  enum ReportFlags {
    INTERSECT       = 0x01,
    HAVE_SEPARATION = 0x02,
    HAVE_POSITION   = 0x04
  };

  class Report {
  public:

    const Convex* c1;
    const Convex* c2;

    int flags;
    ccd_real_t distance;
    vec3 direction;
    vec3 pos1;
    vec3 pos2;
    AlgorithmType algorithm;

    Report(const Convex* c1=0, const Convex* c2=0);
    
  };

  class Checker {
  public:

    ccd_t ccd;
    AlgorithmType algorithm;

    Checker();

    bool intersect(const Convex* c1, const Convex* c2, 
                   ccd_real_t dmin=0) const;

    bool separate(Report& report,
                  const Convex* c1, const Convex* c2, 
                  ccd_real_t dmin=0) const;

    bool penetration(Report& report,
                     const Convex* c1, const Convex* c2,
                     ccd_real_t dmin=0) const;

    bool query(QueryType qtype, Report* report,
               const Convex* c1, const Convex* c2,
               ccd_real_t dmin=0) const;

  };

  //////////////////////////////////////////////////////////////////////

  void cubePoints(size_t n, std::vector<vec3>& points);



};

#endif
