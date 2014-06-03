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

#ifndef _GEOM22_H_
#define _GEOM22_H_

#include "vec3.h"
#include "Box2.h"
#include <vector>

template <class real>
class Geom2_t {
public:

  typedef real value_type;
  typedef vec2_t<real> vec2;
  typedef vec3_t<real> vec3;

  // return the homogeneous coordinates (a,b,c) of the line connecting
  // p0 and p1. if normalize is true, then a^2 + b^2 = 1 and the function
  // pointLineDist will return true signed distances below, and the value
  // of c is the signed distance from the line to (x,y) = (0,0)
  static vec3 lineFromPoints(const vec2& p0, 
                             const vec2& p1,
                             bool normalize=false);

  // return the intersection of two lines as a homogeneous coordinate
  // set.  use proj() to project it back to R^2
  static vec3 lineIntersect(const vec3& l0,
                            const vec3& l1);

  // return the (possibly normalized) signed distance from point p to
  // the line l. if normalized, then 
  static real pointLineDist(const vec2& p,
                            const vec3& l,
                            bool normalize=false);

  // return the closest point on the segment l0-l1 to point p
  static vec2 pointSegmentClosest(const vec2& p,
                                  const vec2& l0,
                                  const vec2& l1);

  // return positive if angle pqr is CCW, negative if CW
  static real turn(const vec2& p, const vec2& q, const vec2& r);

private:

  Geom2_t();
  
};

template <class Point>
class Polygon2_t {
public:
  
  typedef typename Point::value_type value_type;
  typedef vec2_t<value_type> vec2;
  typedef vec3_t<value_type> vec3;
  typedef Box2_t<value_type> Box2;
  typedef std::vector<Point> PointArray;
  typedef Geom2_t<value_type> Geom2;

  static const vec2& p2d(const Point& p);
  static Point makePoint(const vec2& v);
  static Point makePoint(value_type x, value_type y);

  // return true if the input is a convex polygon, wound CCW
  static bool isConvexCCW(const PointArray& v);

  // set dst to the polygon whose outline is the given box
  static void boxPolygon(const Box2& b, PointArray& dst);
  
  // compute the bounding box of a set of points
  static Box2 computeBBox(const PointArray& src);

  // compute the signed area (CW = negative, CCW = positive) for a
  // simple (non-intersecting) polygon.
  static value_type area(const PointArray& src, Point* centroid=0);

  // compute a convex hull of a set of points in O(n lg n) time.
  static void convexHull(const PointArray& src, PointArray& dst);

  // compute the offset (positive outside, negative inside) of a
  // convex polygon
  static void offsetConvexPolygon(const PointArray& src,
                                  value_type offset,
                                  PointArray& dst);

  // clip the given polygon by the given line, keeping any points
  // which have positive signed distance to the line
  static void clipPolygon(const PointArray& src,
                          const vec3& line,
                          PointArray& dst);

  // return the signed distance (positive outside, negative inside)
  // from the convex polygon to the point p. points must represent a
  // convex polygon as produced by convexHull above or verified by
  // isConvexCCW above.
  // 
  // if desired, normal is a vector pointing from p in the direction
  // of the exterior of the polygon (i.e. to increase the distance,
  // travel in the direction of normal).
  static value_type convexPolygonDist(const Point& p, 
				      const PointArray& points,
				      Point* normal=0);

private:

  Polygon2_t();

};

typedef Geom2_t<float> Geom2f;
typedef Geom2_t<double> Geom2d;

typedef Polygon2_t<vec2f> Polygon22f;
typedef Polygon2_t<vec3f> Polygon23f;
typedef Polygon2_t<vec2f> Polygon22d;
typedef Polygon2_t<vec3f> Polygon23d;

#endif
