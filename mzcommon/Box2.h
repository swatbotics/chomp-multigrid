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

#ifndef _BBOX2_H_
#define _BBOX2_H_

#include "vec2.h"
#include "bittwiddle.h"

template <class real> class Box2_t {
public:

  typedef vec2_t<real> vec2;

  vec2 p0;
  vec2 p1;

  Box2_t(): p0(1), p1(-1) {}

  Box2_t(const vec2& pmin, 
            const vec2& pmax): p0(pmin), p1(pmax) {}

  Box2_t(real x0, real y0, real x1, real y1):
    p0(x0, y0), p1(x1, y1) {}

  template <class Tval2>
  explicit Box2_t(const Box2_t<Tval2>& b): p0(b.p0), p1(b.p1) {}

  vec2 lerp(real x, real y) const {
    return vec2( p0.x() + x*(p1.x()-p0.x()),
                 p0.y() + y*(p1.y()-p0.y()) );
  }
  
  vec2 lerp(const vec2& p) const {
    return lerp(p.x(), p.y());
  }

  vec2 corner(unsigned int i) const {
    unsigned int gray = bin2gray(i);
    return lerp( gray & 0x01, (gray & 0x02)>>1 );
  }

  bool empty() const {
    return 
      (p0.x() > p1.x()) || 
      (p0.y() > p1.y());
  }
  
  vec2 center() const {
    return 0.5*(p0+p1);
  }

  bool contains(const vec2& v) const {
    return 
      (v.x() >= p0.x() && v.x() <= p1.x()) &&
      (v.y() >= p0.y() && v.y() <= p1.y());
  }

  void addPoint(const vec2& v) {
    if (empty()) {
      p0 = p1 = v;
    } else {
      for (int i=0; i<2; ++i) {
	if (v[i] < p0[i]) { p0[i] = v[i]; }
	if (v[i] > p1[i]) { p1[i] = v[i]; }
      }
    }
  }

  void dilate(double d) {
    p0 -= vec2(d);
    p1 += vec2(d);
  }

  void clear() {
    p0 = vec2(1);
    p1 = vec2(-1);
  }

  static bool intersects(const Box2_t& b1, const Box2_t& b2) {
    if (b1.empty() || b2.empty()) {
      return false;
    }
    for (int i=0; i<2; ++i) {
      if (b1.p0[i] > b2.p1[i] || b1.p1[i] < b2.p0[i]) { return false; }
    }
    return true;
  }

  static Box2_t unite(const Box2_t& b1, const Box2_t& b2) {
    if (b1.empty()) {
      return b2;
    } else if (b2.empty()) {
      return b1;
    } else {
      Box2_t rval;
      for (unsigned int i=0; i<2; ++i) {
	rval.p0[i] = std::min(b1.p0[i], b2.p0[i]);
	rval.p1[i] = std::max(b1.p1[i], b2.p1[i]);
      }
      return rval;
    }
  }

  static Box2_t intersect(const Box2_t& b1, const Box2_t& b2) {
    if (b1.empty()) {
      return b1;
    } else if (b2.empty()) {
      return b2;
    } else {
      Box2_t rval;
      for (unsigned int i=0; i<2; ++i) {
	rval.p0[i] = std::max(b1.p0[i], b2.p0[i]);
	rval.p1[i] = std::min(b1.p1[i], b2.p1[i]);
      }
      return rval;
    }
  }

  vec2 closest(const vec2& v) const {
    vec2 rval;
    for (int i=0; i<2; ++i) {
      rval[i] = std::max(v[i], p0[i]);
      rval[i] = std::min(rval[i], p1[i]);
    }
    return rval;
  }

  bool clipLine(const vec2& v0, const vec2& v1, real& u0, real& u1) const {
    vec2 delta = v1-v0;
    return _clipLine(v0, v1, delta, u0, u1);
  }

  bool clipLine(vec2& v0, vec2& v1) const {
    vec2 delta = v1-v0;
    real u0, u1;
    if (!clipLine(v0, v1, u0, u1)) { return false; }
    v1 = v0 + delta*u1;
    v0 = v0 + delta*u0;
    return true;
  }

private:

  static bool _clipTest(real p, real q, real& u0, real& u1) {
    if (p == 0 && q < 0) { 
      return false;
    } else if (p < 0) {
      real u = q / p;
      if (u > u1) {
        return false;
      } else if (u > u0) {
        u0 = u;
      }
    } else if (p > 0) {
      real u = q / p;
      if (u < u0) {
        return false;
      } else if (u < u1) {
        u1 = u;
      }
    }
    return true;
  }

  bool _clipLine(const vec2& v0, const vec2& v1, 
                 const vec2& delta,
                 real& u0, real& u1) const {

    u0 = 0;
    u1 = 1;

    for (int i=0; i<2; ++i) {
      if (!_clipTest(-delta[i], v0[i]-p0[i], u0, u1) ||
          !_clipTest( delta[i], p1[i]-v0[i], u0, u1)) {
        return false;
      }
    }

    return true;
    
  }


};

typedef Box2_t<double> Box2d;
typedef Box2_t<float>  Box2f;

#endif
