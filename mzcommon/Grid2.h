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

#ifndef _GRID2_H_
#define _GRID2_H_

#include "Box2.h"
#include "vec2u.h"

template <class real>
class Grid2_t {
public:

  typedef vec2_t<real> vec2;
  typedef Box2_t<real> Box2;

  Grid2_t() { _clear(); }

  const vec2u& dims() const { return _dims; }
  size_t nx() const { return _dims.x(); }
  size_t ny() const { return _dims.y(); }

  size_t size() const { return _size; }

  const vec2& origin() const { return _origin; }
  real cellSize() const { return _cellSize; }

  bool empty() const { return !_size; }

  Box2 bbox() const {
    Box2 rval;
    if (empty()) { return rval; }
    rval.p0 = _origin;
    rval.p1 = _origin + vec2(_cellSize*_dims[0], _cellSize*_dims[1]);
    return rval;
  }

  vec2 center() const {
    real c = real(0.5)*_cellSize;
    return _origin + vec2(c*_dims[0], c*_dims[1]);
  }

  vec2u floorCell(const vec2& pos) const { 
    return _vec2sub(pos, real(-0.5));
  }

  vec2u ceilCell(const vec2& pos) const {
    return _vec2sub(pos, real(0.5));
  }

  vec2u nearestCell(const vec2& pos) const {
    return _vec2sub(pos, 0);
  }

  vec2 cellCenter(size_t x, size_t y) const { 
    return cellCenter(vec2u(x,y));
  }

  vec2 cellCenter(const vec2u& s) const {
    return _origin + vec2(s.x()+real(0.5), s.y()+real(0.5))*_cellSize;
  }

  vec2 cellCenter(size_t idx) const {
    return cellCenter(ind2sub(idx));
  }

  size_t sub2ind(const vec2u& s) const {
    return vec2u::sub2ind(_dims, s);
  }

  size_t sub2ind(size_t x, size_t y) const {
    return vec2u::sub2ind(_dims, x, y);
  }
  
  vec2u ind2sub(size_t idx) const {
    return vec2u::ind2sub(_dims, idx); 
  }

  // alpha holds the blending coefficients for the floor cell
  void sampleCoeffs(const vec2& pos,
                    vec2u& fs,
                    vec2& alpha) const {

    fs = floorCell(pos);
    vec2 fv = cellCenter(fs);
    real invCS = 1/cellSize();
    
    for (int j=0; j<2; ++j) {
      real diff = pos[j] - fv[j];
      if (diff < 0 || diff >= _cellSize || fs[j] + 1 >= _dims[j]) {
        alpha[j] = 1;
      } else {
        real u = diff * invCS;
        alpha[j] = 1-u;
      }
    }
    
  }
  
  // alpha holds the blending coefficients for the floor cell
  template <class Tval>
  Tval sample(const vec2u& fs, const vec2& alpha, const Tval* data) const {

    static const vec2u disp[4] = {
      vec2u(0,0),
      vec2u(0,1),
      vec2u(1,0),
      vec2u(1,1)
    };

    Tval f(0);

    for (int i=0; i<4; ++i) {
      const vec2u d=disp[i];
      vec2u s = fs + d;
      real coeff = 1;
      for (int j=0; j<2; ++j) {
        coeff *= d[j] ? 1-alpha[j] : alpha[j];
      }
      if (!coeff) { continue; }
      f += coeff * data[sub2ind(s)];
    }

    return f;

  }

  template <class Tval>
  Tval sample(const vec2& pos, const Tval* raster) const {

    vec2u fs;
    vec2  alpha;
    sampleCoeffs(pos, fs, alpha);
    return sample(fs, alpha, raster);

  }

  template <class Tfunc>
  bool raycast(const vec2& p0, const vec2& p1, Tfunc& func) const {

    const vec2 dp = p1-p0;

    vec2 c0 = p0;
    vec2 c1 = p1;

    const vec2& o = _origin;
    const float& cs = _cellSize;

    const size_t nx = this->nx();
    const size_t ny = this->ny();

    Box2 brect(o, o+cs*vec2(nx,ny));

    if (!brect.clipLine(c0, c1)) {
      return false;
    }

    vec2u s = nearestCell(c0);
    vec2u sgoal = nearestCell(c1);

    vec2  edge = o + vec2(s[0],s[1])*cs;
    
    vec2 coord;
    real sdir[2], mdir[2], odir[2];

    real mx = fabs(dp[1]);
    real my = fabs(dp[0]);

    for (int i=0; i<2; ++i) {
      coord[i] = (c0[i]-edge[i])/cs;
      if (coord[i] < 0) {
        coord[i] = 0;
      } else if (coord[i] > 1) {
        coord[i] = 1;
      }
      sdir[i] = (dp[i] < 0) ? -1 : (dp[i] > 0) ? 1 : 0;
      odir[i] = sdir[i] > 0 ? 1 : 0;
      mdir[i] = sdir[i] > 0 ? -1 : 1;
    }

    if (func(s, coord)) { return true; }

    while (s != sgoal) {

      // find the distance to horizontal edge
      real dx = odir[0] + mdir[0]*coord[0];
      real dy = odir[1] + mdir[1]*coord[1];

      if (mx*dx < my*dy) {
        // move horizontally
        s[0] += sdir[0];
        coord[0] = (sdir[0] > 0) ? 0 : 1;
        coord[1] += sdir[1]*mx*dx/my;
        //assert(coord[1] >= 0 && coord[1] <= 1);
      } else {
        // move vertically
        s[1] += sdir[1];
        coord[1] = (sdir[1] > 0) ? 0 : 1;
        coord[0] += sdir[0]*my*dy/mx;
        //assert(coord[0] >= 0 && coord[0] <= 1);
      }
    
      if (s[0] >= nx || s[1] >= ny) { // this should never happen but i'm paranoid
        return false;
      } 

      if (func(s, coord)) {
        return true;
      }

    }

    return false;


  }

protected:

  void _clear() {
    _dims = vec2u(0);
    _size = 0;
    _origin = vec2(0);
    _cellSize = 0;
  }

  vec2u _vec2sub(const vec2& pos, real delta) const {
    vec2 v = (pos - _origin) * (1/cellSize());
    vec2u s;
    for (int i=0; i<2; ++i) {
      v[i] = std::max(v[i]+delta, real(0));
      s[i] = std::min(size_t(v[i]), _dims[i]-1);
    }
    return s;
  }

  // straightforward
  void _resize(const vec2u& dims, 
               real cellSize,
               const vec2& origin) {
    _clear();
    _size = dims.prod();
    if (!_size) { return; }
    _dims = dims;
    _cellSize = cellSize;
    _origin = origin;
  }

  void _resize(size_t nx, size_t ny, 
               real cellSize, const vec2& origin) {
    _resize(vec2u(nx, ny), cellSize, origin);
  }
    
  // automatically computes origin and dims so that the bottom left
  // cell contains min and the top right cell contains max and that
  // the center of this is the center of min and max
  void _resize(const vec2& min,
               const vec2& max,
               real cellSize) {
    _clear();
    vec2u dims;
    vec2  origin;
    vec2  center = real(0.5)*(max+min);
    for (int i=0; i<2; ++i) {
      real f = (max[i]-min[i]) / cellSize;
      if (f < 0) { return; }
      dims[i] = size_t(ceil(f));
      origin[i] = center[i] - real(0.5)*cellSize*dims[i];
    }
    _resize(dims, cellSize, origin);
  }

  void _resize(const Box2& bbox, real cellSize) {
    _resize(bbox.p0, bbox.p1, cellSize);
  }
  
  vec2u  _dims;
  size_t _size;
  vec2   _origin;
  real   _cellSize;

};

typedef Grid2_t<float> Grid2f;
typedef Grid2_t<double> Grid2d;

#endif
