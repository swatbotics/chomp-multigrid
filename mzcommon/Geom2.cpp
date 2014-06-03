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

#include "Geom2.h"
#include <algorithm>
#include <assert.h>
#include <float.h>

namespace {
  template <class real>
  class max_traits {
  public:
    static real max();
  };
  template <> float max_traits<float>::max() { return FLT_MAX; }
  template <> double max_traits<double>::max() { return DBL_MAX; }
}

template <class real>
real Geom2_t<real>::turn(const vec2& p, const vec2& q, const vec2& r) {
  return ( (q[0] - p[0])*(r[1] - p[1]) - (r[0] - p[0])*(q[1] - p[1]) );
}


template <class real>
vec3_t<real> Geom2_t<real>::lineFromPoints(const vec2& p0, 
                                           const vec2& p1,
                                           bool normalize) {

  vec3 l = vec3::cross( vec3(p0, 1), vec3(p1, 1) );
  if (normalize) { l /= l.trunc().norm(); }
  return l;

}

template <class real>
vec3_t<real> Geom2_t<real>::lineIntersect(const vec3& l0,
                                          const vec3& l1) {

  return vec3::cross(l0, l1);

}

template <class real>
real Geom2_t<real>::pointLineDist(const vec2& p,
                                  const vec3& l,
                                  bool normalize) {

  real dp = vec3::dot(vec3(p, 1), l);
  if (normalize) { dp /= l.trunc().norm(); }
  return dp;

}


template <class real>
vec2_t<real> Geom2_t<real>::pointSegmentClosest(const vec2& p,
                                                const vec2& l0,
                                                const vec2& l1) {

  // dl points from l0 to l1
  vec2 dl = l1 - l0;
  
  // dp points from l0 to p
  vec2 dp = p - l0;
  
  // numerator is dot(dp, dl)
  real num = vec2::dot(dp, dl);

  // denominator is dot(dl, dl)
  real denom = vec2::dot(dl, dl);

  // get the closest point
  real u = num / denom;
  u = std::max( real(0), std::min( u, real(1) ) );

  // closest point is l0 + u * dl
  return  l0 + u * dl;

}



template <class Point>
static bool pointLessThan(const Point& a, const Point& b) {
  return ( (a.x() < b.x()) ||
           (a.x() == b.x() && a.y() < b.y()) );
}

template <class Point>
static bool pointEqual(const Point& a, const Point& b) {
  return ( (a.x() == b.x()) && (a.y() == b.y()) );
}

template <class real>
static inline const vec2_t<real>& p2d(const vec2_t<real>& p) { 
  return p;
}
template <class real>
static inline const vec2_t<real>& p2d(const vec3_t<real>& p) {
  return *(reinterpret_cast<const vec2_t<real>*>(&p));
}

template <class Point>
const vec2_t<typename Point::value_type>& Polygon2_t<Point>::p2d(const Point& p) {
  return ::p2d(p);
}


template <class Point>
Point Polygon2_t<Point>::makePoint(value_type x, value_type y) {
  Point rval;
  rval.x() = x;
  rval.y() = y;
  return rval;
}

template <class Point>
Point Polygon2_t<Point>::makePoint(const vec2& v) {
  return makePoint(v.x(), v.y());
}


template <class Point>
bool Polygon2_t<Point>::isConvexCCW(const PointArray& v) {

  if (v.size() <= 2) { return true; }

  for (size_t i=0; i<v.size(); ++i) {
    size_t j = (i+1)%v.size();
    size_t k = (i+2)%v.size();
    if (Geom2::turn(p2d(v[i]), p2d(v[j]), p2d(v[k])) <= 0) { 
      return false;
    }
  }

  return true;

}

template <class Point>
void Polygon2_t<Point>::convexHull(const PointArray& src, PointArray& dst) {

  assert( &src != &dst );

  dst.clear();
  if (src.empty()) { return; }

  PointArray sorted = src;
  std::sort(sorted.begin(), sorted.end(), pointLessThan<Point>);

  typename PointArray::const_iterator send = 
    std::unique(sorted.begin(), sorted.end(), pointEqual<Point>);

  sorted.resize( send - sorted.begin() );

  if (sorted.size() < 3) {
    dst.insert(dst.end(), sorted.begin(), sorted.end());
    return;
  }

  // build the convex hull for the bottom
  for (typename PointArray::const_iterator i=sorted.begin(); 
       i!=sorted.end(); ++i) {
    while (dst.size() >= 2 && Geom2::turn( p2d(dst[dst.size()-2]),
						   p2d(dst[dst.size()-1]),
						   p2d(*i) ) <= 0) {
      dst.pop_back();
    }
    dst.push_back(*i);
  }

  // build the convex hull for the top - this will add two duplicate points
  size_t tsize = dst.size();

  for (typename PointArray::const_reverse_iterator i=sorted.rbegin(); 
       i!=sorted.rend(); ++i) {
    while (dst.size() >= tsize+2 && Geom2::turn( p2d(dst[dst.size()-2]),
							 p2d(dst[dst.size()-1]),
							 p2d(*i) ) <= 0) {
      dst.pop_back();
    }
    dst.push_back(*i);
  }

  // kill the last thing
  dst.pop_back();

  // move everyone back one space
  for (size_t i=tsize; i<dst.size()-1; ++i) {
    dst[i] = dst[i+1];
  }

  // kill the last thing
  dst.pop_back();

}

template <class Point>
void Polygon2_t<Point>::boxPolygon(const Box2& b, PointArray& dst) {
  dst.clear();
  if (!b.empty()) {
    for (int i=0; i<4; ++i) {
      dst.push_back(makePoint(b.corner(i)));
    }
  }
}

  // compute the centroid of a set of points
template <class Point>
typename Point::value_type 
Polygon2_t<Point>::area(const PointArray& src, Point* centroid) {

  if (centroid) { *centroid = Point(0); }

  if (src.empty()) {
    return 0;
  } else if (src.size() < 3) {
    if (centroid) { 
      for (size_t i=0; i<src.size(); ++i) { *centroid += src[i]; }
      *centroid /= src.size();
    }
    return 0;
  }

  value_type A = 0;

  for (size_t i0=0; i0<src.size(); ++i0) {

    size_t i1 = (i0+1)%src.size();

    const Point& p0 = src[i0];
    const Point& p1 = src[i1];

    value_type t = p0.x()*p1.y() - p1.x()*p0.y();
    A += t;

    if (centroid) {
      *centroid += (p0 + p1) * t;
    }

  }

  A *= 0.5;

  if (centroid) {
    *centroid *= 1.0 / (6*A);
  }

  return A;

}
  
template <class Point>
Box2_t<typename Point::value_type> 
Polygon2_t<Point>::computeBBox(const PointArray& src) {
  Box2 rval;
  for (typename PointArray::const_iterator i=src.begin();
       i!=src.end(); ++i) {
    rval.addPoint(p2d(*i));
  }
  return rval;
}

template <class real>
vec2_t<real> handleIntersect(const vec2_t<real>& l1, 
			     const vec2_t<real>& l2,
			     const vec2_t<real>& p) {
  return p;
}

template <class real>
vec3_t<real> handleIntersect(const vec3_t<real>& l1, 
			     const vec3_t<real>& l2,
			     const vec2_t<real>& p) {
  vec2_t<real> dl = l2.trunc() - l1.trunc();
  vec2_t<real> dp = p - l1.trunc();
  int idx = ( fabs(dl.x()) > fabs(dl.y()) ? 0 : 1 );
  real u = dp[idx] / dl[idx];
  return vec3_t<real>(p, l1.z() + u*(l2.z() - l1.z()));
}

template <class Point>
void Polygon2_t<Point>::offsetConvexPolygon(const PointArray& src,
						  value_type offset,
						  PointArray& dst) {

  assert( &src != &dst );

  dst.clear();

  if (src.size() < 3) {

    return;

  } else if (offset == 0) { 

    dst = src;
    return;

  } else if (offset < 0) {

    if (src.size() <= 4 || Point::size == 3) {
      dst = src;
    } else {
      Box2 b = computeBBox(src);
      b.dilate(-2*offset);
      boxPolygon(b, dst);
    }

    PointArray tmp;

    for (size_t i0=0; i0<src.size() && !dst.empty(); ++i0) {
      size_t i1 = (i0 + 1) % src.size();
      vec3 line = Geom2::lineFromPoints(p2d(src[i0]), p2d(src[i1]), true);
      line[2] += offset;
      clipPolygon(dst, line, tmp);
      dst.swap(tmp);
    }

  } else {

    for (size_t i0=0; i0<src.size(); ++i0) {

      size_t i1 = (i0 + 1) % src.size();
      size_t i2 = (i0 + 2) % src.size();

      vec3 l0 = Geom2::lineFromPoints(p2d(src[i0]), p2d(src[i1]), true);
      vec3 l1 = Geom2::lineFromPoints(p2d(src[i1]), p2d(src[i2]), true);

      l0[2] += offset;
      l1[2] += offset;

      vec3 h = Geom2::lineIntersect(l0, l1);
      assert( fabs(h[2]) > 1e-9 );
      dst.push_back( handleIntersect(src[i0], src[i1], h.proj()) );
      
    }

  }



}

template <class Point>
void Polygon2_t<Point>::clipPolygon(const PointArray& src,
					  const vec3& line,
					  PointArray& dst) {

  assert( &src != &dst );

  // TODO: writeme
  dst.clear();

  if (src.empty()) { return; }

  bool was_inside = false;

  for (size_t i=0; i<=src.size(); ++i) {

    const Point& p1 = src[i % src.size()];

    // see if current point is inside the clip region
    bool is_inside = (Geom2::pointLineDist(p2d(p1), line, true) > 1e-6);

    bool add_intersect = (i && (is_inside != was_inside));

    if (add_intersect) {
      const Point& p0 = src[i-1];
      vec3 edge = Geom2::lineFromPoints(p2d(p0), p2d(p1));
      vec3 h = Geom2::lineIntersect(line, edge);
      if (fabs(h[2]) <= 1e-9) {
	std::cerr << "WARNING: CAN'T HANDLE PARALLEL EDGES IN CLIP YET!\n";
	dst.clear();
	return;
      }
      dst.push_back(handleIntersect(p0, p1, h.proj()));
    }

    if (is_inside && i < src.size()) {
      dst.push_back(p1);
    }

    was_inside = is_inside;

  }

}


// normal points into the polygon
template <class Point>
typename Point::value_type 
Polygon2_t<Point>::convexPolygonDist(const Point& p,
				     const PointArray& points,
				     Point* normal) {

  if (points.size() == 1) {

    vec2 diff = p2d(p) - p2d(points.front());
    value_type d = diff.norm();
    if (normal) { *normal = makePoint(diff/d); }
    return d;

  } else if (points.size() == 2) {

    vec2 pc = Geom2::pointSegmentClosest(p2d(p), p2d(points[0]), p2d(points[1]));
    vec2 diff = p2d(p) - pc;
    value_type d = diff.norm();
    if (normal) { *normal = makePoint(diff/d); }
    return d;

  } else {
    
    bool inside = true;
    value_type dmin = max_traits<value_type>::max();
    
    for (size_t j0=0; j0<points.size(); ++j0) {

      size_t j1 = (j0+1) % points.size();

      const Point& l0 = points[j0];
      const Point& l1 = points[j1];
      
      // positive distance is INSIDE
      vec3 l = Geom2::lineFromPoints(p2d(l0), p2d(l1));

      value_type sd = Geom2::pointLineDist(p2d(p), l);
      if (sd <= 0) { inside = false; }

      vec2 pc = Geom2::pointSegmentClosest(p2d(p), p2d(l0), p2d(l1));

      vec2 diff = p2d(p)-pc;
      value_type d = diff.norm();

      if (d < dmin) {
        dmin = d;
        if (normal) { *normal = makePoint(diff / ((sd > 0) ? -d : d)); }
      }

    }

    if (inside) { dmin = -dmin; }
    
    return dmin;

  }

}


template class Geom2_t<float>;
template class Geom2_t<double>;

template class Polygon2_t<vec2f>;
template class Polygon2_t<vec3f>;
template class Polygon2_t<vec2d>;
template class Polygon2_t<vec3d>;

