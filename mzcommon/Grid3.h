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

#ifndef _GRID3_H_
#define _GRID3_H_

#include "Box3.h"
#include "vec3u.h"

template <class real>
class Grid3_t {
public:

  typedef vec3_t<real> vec3;
  typedef Box3_t<real> Box3;

  Grid3_t() { _clear(); }

  const vec3u& dims() const { return _dims; }
  size_t nx() const { return _dims.x(); }
  size_t ny() const { return _dims.y(); }
  size_t nz() const { return _dims.z(); }

  size_t size() const { return _size; }

  const vec3& origin() const { return _origin; }
  real cellSize() const { return _cellSize; }

  bool empty() const { return !_size; }


  Box3 bbox() const {
    Box3 rval;
    if (empty()) { return rval; }
    rval.p0 = _origin;
    rval.p1 = _origin + vec3(_cellSize*_dims[0], 
                             _cellSize*_dims[1],
                             _cellSize*_dims[2]);
    return rval;
  }

  vec3 center() const {
    real c = real(0.5)*_cellSize;
    return _origin + vec3(c*_dims[0], c*_dims[1], c*_dims[2]);
  }

  vec3u floorCell(const vec3& pos) const { 
    return _vec2sub(pos, real(-0.5));
  }

  vec3u ceilCell(const vec3& pos) const {
    return _vec2sub(pos, real(0.5));
  }

  vec3u nearestCell(const vec3& pos) const {
    return _vec2sub(pos, 0);
  }

  vec3 cellCenter(size_t x, size_t y, size_t z) const { 
    return cellCenter(vec3u(x,y,z));
  }

  vec3 cellCenter(const vec3u& s) const {
    return _origin + vec3(s.x()+real(0.5), 
                          s.y()+real(0.5), 
                          s.z()+real(0.5))*_cellSize;
  }

  vec3 cellCenter(size_t idx) const {
    return cellCenter(ind2sub(idx));
  }

  size_t sub2ind(const vec3u& s) const {
    return vec3u::sub2ind(_dims, s);
  }

  size_t sub2ind(size_t x, size_t y, size_t z) const {
    return vec3u::sub2ind(_dims, x, y, z);
  }
  
  vec3u ind2sub(size_t idx) const {
    return vec3u::ind2sub(_dims, idx); 
  }

  void sampleCoeffs(const vec3& pos,
                    vec3u& fs,
                    vec3& alpha) const {

    fs = floorCell(pos);
    vec3 fv = cellCenter(fs);
    real invCS = 1/cellSize();
    
    for (int j=0; j<3; ++j) {
      real diff = pos[j] - fv[j];
      if (diff < 0 || diff >= _cellSize || fs[j] + 1 >= _dims[j]) {
        alpha[j] = 1;
      } else {
        real u = diff * invCS;
        alpha[j] = 1-u;
      }
    }
    
  }
  
  template <class Tval>
  Tval sample(const vec3u& fs, const vec3& alpha, const Tval* data) const {

    static const vec3u disp[8] = {
      vec3u(0,0,0),
      vec3u(0,0,1),
      vec3u(0,1,0),
      vec3u(0,1,1),
      vec3u(1,0,0),
      vec3u(1,0,1),
      vec3u(1,1,0),
      vec3u(1,1,1)
    };

    Tval f(0);

    for (int i=0; i<8; ++i) {
      const vec3u d=disp[i];
      vec3u s = fs + d;
      real coeff = 1;
      for (int j=0; j<3; ++j) {
        coeff *= d[j] ? 1-alpha[j] : alpha[j];
      }
      if (!coeff) { continue; }
      f = f + coeff * data[sub2ind(s)];
    }

    return f;

  }

  template <class Tval>
  Tval sample(const vec3& pos, const Tval* raster) const {
    vec3u fs;
    vec3  alpha;
    sampleCoeffs(pos, fs, alpha);
    return sample(fs, alpha, raster);
  }

protected:

  void _clear() {
    _dims = vec3u(0);
    _size = 0;
    _origin = vec3(0);
    _cellSize = 0;
  }

  vec3u _vec2sub(const vec3& pos, real delta) const {
    vec3 v = (pos - _origin) * (1/cellSize());
    vec3u s;
    for (int i=0; i<3; ++i) {
      v[i] = std::max(v[i]+delta, real(0));
      s[i] = std::min(size_t(v[i]), _dims[i]-1);
    }
    return s;
  }

  // straightforward
  void _resize(const vec3u& dims, 
               real cellSize,
               const vec3& origin) {
    _clear();
    _size = dims.prod();
    if (!_size) { return; }
    _dims = dims;
    _cellSize = cellSize;
    _origin = origin;
  }

  void _resize(size_t nx, size_t ny, size_t nz,
               real cellSize, const vec3& origin) {
    _resize(vec3u(nx, ny, nz), cellSize, origin);
  }
    
  // automatically computes origin and dims so that the bottom left
  // cell contains min and the top right cell contains max and that
  // the center of this is the center of min and max
  void _resize(const vec3& min,
               const vec3& max,
               real cellSize) {
    _clear();
    vec3u dims;
    vec3  origin;
    vec3  center = real(0.5)*(max+min);
    for (int i=0; i<3; ++i) {
      real f = (max[i]-min[i]) / cellSize;
      if (f < 0) { return; }
      dims[i] = size_t(ceil(f));
      origin[i] = center[i] - real(0.5)*cellSize*dims[i];
    }
    _resize(dims, cellSize, origin);
  }

  void _resize(const Box3& bbox, real cellSize) {
    _resize(bbox.p0, bbox.p1, cellSize);
  }
  
  vec3u  _dims;
  size_t _size;
  vec3   _origin;
  real   _cellSize;

};

typedef Grid3_t<float> Grid3f;
typedef Grid3_t<double> Grid3d;

#endif
