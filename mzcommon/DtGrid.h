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

#ifndef _DTGRID_H_
#define _DTGRID_H_

#include "TriMesh3.h"
#include "Grid3.h"
#include "vec2u.h"
#include "vec2.h"
#include "HeightMap.h"

template <class real>
class DtGrid_t: public Grid3_t<real> {
public:

  typedef TriMesh3_t<real>   TriMesh3;
  typedef Transform3_t<real> Transform3;
  typedef HeightMap_t<real>  HeightMap;
  typedef vec3_t<real>       vec3;
  typedef vec2_t<real>       vec2;
  typedef std::vector<real>  RealArray;
  typedef std::vector<vec3>  Vec3Array;
  typedef std::vector<int>   IntArray;
  typedef std::vector<bool>  BoolArray;

  static const real DT_INF;

  enum Axis {
    AXIS_X = 0,
    AXIS_Y = 1,
    AXIS_Z = 2,
  };

  //////////////////////////////////////////////////////////////////////

  DtGrid_t();

  ~DtGrid_t();

  //////////////////////////////////////////////////////////////////////

  const vec3u& referenceAxes() const;

  Axis referenceAxis() const;

  void clear();

  void resize(size_t nx, size_t ny, size_t nz, 
              Axis referenceAxis,
              real cellSize, const vec3& origin);

  void resize(const vec3& min,
              const vec3& max,
              Axis referenceAxis,
              real cellSize);
  //////////////////////////////////////////////////////////////////////

  const real& operator[](size_t idx) const;

  real& operator[](size_t idx);

  const real& operator()(const vec3u& s) const;

  real& operator()(const vec3u& s);

  const real& operator()(size_t x, size_t y, size_t z) const;

  real& operator()(size_t x, size_t y, size_t z);

  vec3 gradient(size_t x, size_t y, size_t z) const;

  vec3 gradient(vec3u s) const;

  vec3 gradient(size_t idx) const;

  vec3 normal(size_t x, size_t y, size_t z) const;

  vec3 normal(const vec3u& s) const;

  vec3 normal(size_t idx) const;

  real sample(const vec3& v) const;

  real sample(const vec3& v, vec3& gradient) const;

  // returns the position of minimum value along the line connecting
  // s1 and s2
  real lineMin(const vec3u& s1, const vec3u& s2, vec3u& smin) const;

  // returns the minimum value along the line connecting v1 and v2
  real lineMin(const vec3& v1, const vec3& v2, 
                vec3& vmin, vec3& gmin) const;

  //////////////////////////////////////////////////////////////////////

  const real& height(size_t u, size_t v) const;
  real& height(size_t u, size_t v);

  const real& height(const vec2u& s) const;
  real& height(const vec2u& s);

  real height(const vec3& v) const;
  
  vec3 heightVec(size_t u, size_t v) const;
  vec3 heightVec(const vec2u& s) const;
  vec3 heightVec(const vec3& v) const;

  real minHeight() const;
  real maxHeight() const;

  const HeightMap& heightMap() const;
  HeightMap& heightMap();
  
  //////////////////////////////////////////////////////////////////////

  real minDist() const;
  real maxDist() const;

  void scanConvert(const TriMesh3& model, 
                   bool asHeightmap=false);

  void scanConvert(const TriMesh3& model,
                   const Transform3& transform, 
                   bool asHeightmap=false);

  void processPointCloudAsHeightmap(const std::vector<vec3>& points,
                                    size_t minBinSize=5);


  void computeDists(bool storeGradients=true);

  void computeDistsFromHeightMap(bool storeGradients=true);

  void computeDistsFromBinary(bool storeGradients=true);

  void recomputeExtents();

  static void dt(const RealArray& f,
                 size_t n,
                 RealArray& ft,
                 IntArray& v,
                 RealArray& z);

  bool load(const char* filename, bool storeGradients=true);
  void save(const char* filename) const;
   
private:

  void _create();
  void _createGradients();
  void _computeEDT();
  
  real _sample(const vec3& v, vec3* gradient) const;

  void _scanConvert(const TriMesh3& model, 
                    const Transform3* transform,
                    bool asHeightmap);

  vec3u  _ax;
  
  RealArray _data;

  // gradient data
  Vec3Array _gdata;

  real _minDist;
  real _maxDist;

  // heightmap data
  HeightMap _hmap;

};

typedef DtGrid_t<float> DtGridf;
typedef DtGrid_t<double> DtGridd;

#endif
