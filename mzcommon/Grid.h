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

#ifndef _GRID_H_
#define _GRID_H_

#include "vec.h"

template <size_t ndims, class real>
class Grid_t {
public:

  typedef vecu_t<ndims> vecu;
  typedef vec_t<ndims, real> vec;

  Grid_t() { _clear(); }

  const vecu& dims() const { return _dims; }

  const vecu& rprod() const { return _rprod; }

  size_t size() const { return _rprod[ndims-1]; }

  const vec& origin() const { return _origin; }

  const vec& cellSizes() const { return _cellSizes; }

  bool empty() const { return !size(); }

  vecu floorCell(const vec& pos) const {
    return _vec2sub(pos, real(-0.5));
  }

  vecu ceilCell(const vec& pos) const {
    return _vec2sub(pos, real(0.5));
  }

  vecu nearestCell(const vec& pos) const {
    return _vec2sub(pos, real(0));
  }

  vec center() const {
    vec rval;
    for (size_t i=0; i<ndims; ++i) {
      rval[i] = _origin[i] + _dims[i]*0.5*_cellSizes[i];
    }
    return rval;
  }
  
  vec max() const {
    vec rval;
    for (size_t i=0; i<ndims; ++i) {
      rval[i] = _origin[i] + _dims[i] * _cellSizes[i];
    }
    return rval;
  }

  vec cellCenter(size_t i) { return cellCenter(ind2sub(i)); }

  vec cellCenter(const vecu& s) const {
    vec rval;
    for (size_t i=0; i<ndims; ++i) {
      rval[i] = _origin[i] + (s[i] + real(0.5)) * _cellSizes[i];
    }
    return rval;
  }

  size_t sub2ind(const vecu& s) const {
    return vecu::sub2ind(_rprod, s);
  }

  vecu ind2sub(size_t idx) const {
    return vecu::ind2sub(_rprod, idx);
  }

  void fracCell(const vec& p, vecu& s, vec& u) const {
    s = floorCell(p);
    vec c = cellCenter(s);
    for (size_t i=0; i<ndims; ++i) {
      real di = p[i] - c[i];
      if (di < 0 || di >= _cellSizes[i] || s[i] + 1 >= _dims[i]) {
        u[i] = 0;
      } else {
        u[i] = di / _cellSizes[i];
      }
    }
  }

  template <class Tval>
  Tval sample(const vecu& s, const vec& u, const Tval* data) const {
    Tval f(0);
    vecu sn;
    size_t expval = 1<<ndims;
    for (size_t n=0; n<expval; ++n) {
      real w = 1;
      for (size_t i=0; i<ndims; ++i) {
        bool di = (n & (1<<i));
        sn[i] = s[i] + (di ? 1 : 0);
        w *= (di ? u[i] : (real(1) - u[i]));
        if (!w) { break; }
      }
      if (w) { f = f + w * data[sub2ind(sn)]; }
    }
    return f;
  }

  template <class Tval>
  Tval sample(const vec& pos, const Tval* data) const {
    vecu s;
    vec  u;
    fracCell(pos, s, u);
    return sample(s, u, data);
  }
  
  void simplexCoeffs(const vec& pos, 
                     vecu points[ndims+1], 
                     real coeffs[ndims+1]) const {

    vec u;

    fracCell(pos, points[0], u);

    bool flip[ndims];

    for (size_t i=0; i<ndims; ++i) {
      if (points[0][i] % 2) { 
        flip[i] = true;
        u[i] = 1-u[i];
        ++points[0][i];
      } else {
        flip[i] = false;
      }
    }

    vecu sorted = u.rsort();

    coeffs[0] = 1;
    
    for (size_t i=0; i<ndims; ++i) {
      points[i+1] = points[i];
      if (flip[sorted[i]]) {
        --points[i+1][sorted[i]];
      } else {
        ++points[i+1][sorted[i]];
      }
      if (i+1 < ndims) {
        coeffs[i+1] = u[sorted[i]] - u[sorted[i+1]];
      } else {
        coeffs[i+1] = u[sorted[i]];
      }
      coeffs[0] -= coeffs[i+1];
    }

  }


  template <class Tval>
  Tval sampleSimplex(const vecu points[ndims+1],
                     const real coeffs[ndims+1],
                     const Tval* data) const {
    Tval f(0);
    for (size_t i=0; i<ndims+1; ++i) {
      if (coeffs[i]) {  f = f + coeffs[i] * data[sub2ind(points[i])]; }
    }
    return f;
  }

  template <class Tval>
  Tval sampleSimplex(const vec& pos, const Tval* data) const {
    vecu points[ndims+1];
    real coeffs[ndims+1];
    simplexCoeffs(pos, points, coeffs);
    return sampleSimplex(points, coeffs, data);
  }
  
protected:

  void _clear() {
    _dims = vecu(0);
    _rprod = vecu(0);
    _origin = vec(0);
    _cellSizes = vec(0);
  }

  vecu _vec2sub(const vec& pos, real delta) const {
    vec v = (pos - _origin);
    for (size_t i=0; i<ndims; ++i) { v[i] /= _cellSizes[i]; }
    vecu s;
    for (int i=0; i<2; ++i) {
      v[i] = std::max(v[i]+delta, real(0));
      s[i] = std::min(size_t(v[i]), _dims[i]-1);
    }
    return s;
  }

  void _resize(const vecu& dims, 
               const vec& cellSizes,
               const vec& origin) {
    _clear();
    _rprod = dims.rprod();
    if (!size()) { return; }
    _dims = dims;
    _cellSizes = cellSizes;
    _origin = origin;
  }

  void _resize(const vec& min,
               const vec& max,
               const vec& cellSizes) {
    _clear();
    vecu dims;
    vec  origin;
    vec  center = real(0.5)*(max+min);
    for (int i=0; i<2; ++i) {
      real f = (max[i]-min[i]) / cellSizes[i];
      if (f < 0) { return; }
      dims[i] = size_t(ceil(f));
      origin[i] = center[i] - real(0.5)*cellSizes[i]*dims[i];
    }
    _resize(dims, cellSizes, origin);
  }

  vecu  _dims;
  vecu  _rprod;
  vec   _origin;
  vec   _cellSizes;

};

#endif
