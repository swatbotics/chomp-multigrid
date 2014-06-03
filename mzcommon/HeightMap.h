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

#ifndef _HEIGHTMAP_H_
#define _HEIGHTMAP_H_

#include <vector>
#include "Grid2.h"
#include "vec3.h"
#include "vec3u.h"
#include "Box2.h"
#include "TriMesh3.h"

template <class real>
class HeightMap_t: public Grid2_t<real> {
public:

  typedef Grid2_t<real> Grid2;
  typedef TriMesh3_t<real> TriMesh3;
  typedef vec2_t<real> vec2;
  typedef vec3_t<real> vec3;
  typedef Box2_t<real> Box2;
  typedef Box3_t<real> Box3;
  typedef std::vector<real> RealArray;
  typedef std::vector<RealArray*> BinStorage;
  static const real INVALID_HEIGHT;

  struct Split {
    vec2u  s0;              // lower bound (inclusive)
    vec2u  s1;              // upper bound (exclusive)
    Box3   box;             // bounding box of points in this split
    size_t parent_index;    // index of parent split (or -1 if root)
    size_t child_index[2];  // indices of child splits (or -1 if no children)
  };

  // array of splits
  typedef std::vector<Split> SplitArray;

  HeightMap_t();


  static Box2 boundPoints(const std::vector<vec3>& points,
                          const vec3u& ax=vec3u(0,1,2));

  Box3 computeBBox() const;

  static void clearBins(BinStorage& bins);

  void binPoints(BinStorage& bins,
                 const std::vector<vec3>& points,
                 const vec3u& ax=vec3u(0,1,2)) const;

  void medianMap(BinStorage& bins, size_t minCount=1);
  
  void computeExtents(const vec2u& s0, 
		      const vec2u& s1,
		      real& hmin,
		      real& hmax) const;

  void clear();

  void resize(size_t nx, size_t ny, 
              real cellSize, 
              const vec2& origin);

  void resize(const vec2& min,
              const vec2& max,
              real cellSize);

  const real& operator[](size_t idx) const;
  real& operator[](size_t idx);

  const real& operator()(const vec2u& s) const;
  real& operator()(const vec2u& s);

  const real& operator()(size_t x, size_t y) const;
  real& operator()(size_t x, size_t y);

  static real correctHeight(real h, real defaultValue);
  real correctHeight(size_t x, size_t y, real defaultValue) const;
  real correctHeight(const vec2u& s, real defaultValue) const;

  vec2 slope(size_t x, size_t y) const;
  vec3 normal(size_t x, size_t y) const;

  vec2 slope(const vec2u& s) const;
  vec3 normal(const vec2u& s) const;

  void recomputeExtents();

  bool load(const char* filename);
  void save(const char* filename) const;

  real minHeight() const;
  real maxHeight() const;

  void generateMesh(TriMesh3& mesh,
                    real dropEdges=0,
                    const vec3u& ax=vec3u(0,1,2)) const;

  void subdivide(SplitArray& splits) const;
  

private:

  void _medianMap(const std::vector<vec3>& points,
                  const vec3u& ax=vec3u(0,1,2));

  void _subdivide(SplitArray& splits, 
		  size_t parent_index,
		  size_t which_child,
		  const vec2u& s0,
		  const vec2u& s1) const;

  RealArray _hdata;
  real _minHeight;
  real _maxHeight;
  
};

typedef HeightMap_t<float> HeightMapf;
typedef HeightMap_t<double> HeightMapd;

#endif
