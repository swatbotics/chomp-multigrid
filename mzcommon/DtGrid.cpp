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

#include "DtGrid.h"
#include <float.h>
#include <iostream>
#include <stdio.h>
#include <stdexcept>
#include <assert.h>
#include <fstream>
#include "Bresenham.h"


template <> const float DtGrid_t<float>::DT_INF = FLT_MAX;
template <> const double DtGrid_t<double>::DT_INF = DBL_MAX;

#define debug if (0) std::cerr

#include <map>

template <class real>
DtGrid_t<real>::DtGrid_t() {
  clear();
}

template <class real>
DtGrid_t<real>::~DtGrid_t() {}

template <class real>
real DtGrid_t<real>::minDist() const { return _minDist; }

template <class real>
real DtGrid_t<real>::maxDist() const { return _maxDist; }


template <class real>
void DtGrid_t<real>::clear() { 
  this->_clear();
  _data.clear();
  _gdata.clear();
  _hmap.clear();
  _minDist = DT_INF;
  _maxDist = -DT_INF;
}

template <class real>
real& DtGrid_t<real>::operator()(size_t x, size_t y, size_t z) {
  return _data[this->sub2ind(x,y,z)];
}

template <class real>
const real& DtGrid_t<real>::operator()(size_t x, size_t y, size_t z) const {
  return _data[this->sub2ind(x,y,z)];
}

template <class real>
const real& DtGrid_t<real>::operator()(const vec3u& s) const {
  return _data[this->sub2ind(s)];
}

template <class real>
real& DtGrid_t<real>::operator()(const vec3u& s) {
  return _data[this->sub2ind(s)];
}

template <class real>
const real& DtGrid_t<real>::operator[](size_t idx) const {
  return _data[idx];
}

template <class real>
real& DtGrid_t<real>::operator[](size_t idx) {
  return _data[idx];
}

template <class real>
real DtGrid_t<real>::_sample(const vec3& v, vec3* grad) const {

  vec3u fs = this->floorCell(v);
  vec3 fv = this->cellCenter(fs);

  vec3 alpha[2];

  real invCS = 1/this->_cellSize;

  for (int j=0; j<3; ++j) {
    if (v[j] < fv[j]) {
      alpha[0][j] = 1.0;
      alpha[1][j] = 0.0;
    } else {
      real diff = v[j] - fv[j];
      if (diff >= this->_cellSize || fs[j] + 1 >= this->_dims[j]) {
        alpha[0][j] = 1.0;
        alpha[1][j] = 0.0;
      } else {
        real u = diff * invCS;
        alpha[0][j] = 1.0f-u;
        alpha[1][j] = u;
      }
    }
  }

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

  real f=0;
  vec3 g(0);

  for (int i=0; i<8; ++i) {
    const vec3u d=disp[i];
    vec3u s = fs + d;
    real coeff = 1.0f;
    for (int j=0; j<3; ++j) { 
      coeff *= alpha[d[j]][j];
    }
    if (!coeff) { continue; }
    f += coeff * _data[this->sub2ind(s)];
    if (grad) {
      g += coeff * gradient(s);
    }
  }

  if (grad) { *grad = g; }

  return f;

}

template <class real>
real DtGrid_t<real>::sample(const vec3& v) const {
  return _sample(v, 0);
}

template <class real>
real DtGrid_t<real>::sample(const vec3& v, vec3& gradient) const {
  return _sample(v, &gradient);
}
  

template <class real>
vec3_t<real> DtGrid_t<real>::gradient(size_t x, size_t y, size_t z) const {
  return gradient(vec3u(x,y,z));
}

template <class real>
vec3_t<real> DtGrid_t<real>::gradient(vec3u s) const {

  vec3 g(0);
  if (!this->_size) { return g; }

  if (!_gdata.empty()) {
    return _gdata[this->sub2ind(s)];
  }

  real vcur = _data[this->sub2ind(s)];
  real invCS = 1/this->_cellSize;

  for (int d=0; d<3; ++d) {

    if (this->_dims[d] == 1) { continue; }

    if (s[d] > 0) { 

      --s[d];
      real vprev = _data[this->sub2ind(s)];
      ++s[d];

      if (s[d]+1 < this->_dims[d]) {

        ++s[d];
        real vnext = _data[this->sub2ind(s)];
        --s[d];

        g[d] = (vnext-vprev) * invCS * 0.5f;

      } else {

        g[d] = (vcur-vprev) * invCS;

      }

    } else {

      assert( s[d]+1 < this->_dims[d] );

      ++s[d];
      real vnext = _data[this->sub2ind(s)];
      --s[d];

      g[d] = (vnext-vcur) * invCS;

    }
      
  }

  return g;

}

template <class real>
vec3_t<real> DtGrid_t<real>::gradient(size_t idx) const {
  return gradient(this->ind2sub(idx));
}

template <class real>
static vec3_t<real> safeNormalize(vec3_t<real> v) {
  real vn = v.norm();
  if (vn > 1e-12) {
    return v / vn;
  } else {
    return vec3_t<real>(0);
  }
}

template <class real>
vec3_t<real> DtGrid_t<real>::normal(size_t x, size_t y, size_t z) const {
  return normal(vec3u(x,y,z));
}

template <class real>
vec3_t<real> DtGrid_t<real>::normal(const vec3u& s) const {
  return safeNormalize(gradient(s));
}

template <class real>
vec3_t<real> DtGrid_t<real>::normal(size_t idx) const {
  return safeNormalize(gradient(this->ind2sub(idx)));
}

template <class real>
void DtGrid_t<real>::resize(size_t nx, size_t ny, size_t nz,
                    Axis referenceAxis,
                    real cellSize, const vec3& origin) {
  clear();
  this->_resize(nx,ny,nz,cellSize,origin);
  if (this->empty()) { return; }
  _ax[2] = referenceAxis;
  _create();
}

template <class real>
void DtGrid_t<real>::resize(const vec3& min,
                    const vec3& max,
                    Axis referenceAxis,
                    real cellSize) {
  clear();
  this->_resize(min,max,cellSize);
  if (this->empty()) { return; }
  _ax[2] = referenceAxis;
  
  _create();

}


template <class real>
void DtGrid_t<real>::_create() {

  _ax[0] = (_ax[2] + 1) % 3;
  _ax[1] = (_ax[2] + 2) % 3;
  _ax[2] = (_ax[2] + 0) % 3;
  _data.clear();
  _gdata.clear();
  
  this->_size = this->_dims.prod();

  _data.resize(this->_size, DT_INF);

  _hmap.resize(this->_dims[_ax[0]], this->_dims[_ax[1]], 
               this->_cellSize,
               vec2(this->_origin[_ax[0]], this->_origin[_ax[1]]));

  real mb = real(this->_size*sizeof(real)) / (1024.0*1024.0);
  mb = ceil(mb*100)/100;
  
  //std::cerr << "dims = " << this->_dims << ", size = " << _size << " (" << mb << "MB)\n";
  //std::cerr << "origin = " << _origin << "\n";
  //std::cerr << "cellsize = " << _cellSize << "\n";

}

template <class real>
static inline real sqr(real x) { return x*x; }

template <class real>
void DtGrid_t<real>::dt(const RealArray& f,
                        size_t nn,
                        RealArray& ft,
                        IntArray& v,
                        RealArray& z) {

  assert(f.size() >= nn);
  assert(ft.size() >= nn);
  assert(v.size() >= nn);
  assert(z.size() >= nn+1);

  int n = nn;

  int k = 0;

  v[0] = 0;

  z[0] = -DT_INF;
  z[1] =  DT_INF;

  for (int q=1; q<n; ++q) {
    real s = ((f[q]+sqr(q))-(f[v[k]]+sqr(v[k])))/(2*q-2*v[k]);
    while (s <= z[k]) {
      --k;
      s = ((f[q]+sqr(q))-(f[v[k]]+sqr(v[k])))/(2*q-2*v[k]);
    }
    ++k;
    v[k] = q;
    z[k] = s;
    z[k+1] = DT_INF;
  }

  k = 0;
  for (int q=0; q<n; ++q) {
    while (z[k+1] < q) { ++k; }
    ft[q] = sqr(q-v[k]) + f[v[k]];
  }
  
}

template <class real>
static inline real gaussian(real x, real sigma2) {
  return exp(-(x*x) / (2*sigma2));
}


template <class real>
void DtGrid_t<real>::computeDistsFromBinary(bool storeGradients) {

  if (this->empty()) { return; }

  // set bintmp to true wherever data <= 0 and false wherever data > 0
  // also set data to 0 wherever bintmp is true and inf otherwise
  std::vector<bool> bintmp(this->_size);

  for (size_t i=0; i<this->_size; ++i) {
    bintmp[i] = (_data[i] <= 0);
    _data[i] = (bintmp[i] ? 0 : DT_INF);
  }

  // compute the EDT of the data so far
  _computeEDT();

  // copy data into a temp vector
  std::vector<real> tmp = _data;

  // now set data to 0 wherever bintmp is false and inf otherwise
  for (size_t i=0; i<this->_size; ++i) {
    _data[i] = (bintmp[i] ? DT_INF : 0);
  }

  // compute the EDT of the data again
  _computeEDT();

  _minDist = DT_INF;
  _maxDist = -DT_INF;

  for (size_t i=0; i<this->_size; ++i) {
    real d1 = std::max(tmp[i]-0.5f*this->_cellSize, real(0));
    real d2 = std::max(_data[i]-0.5f*this->_cellSize, real(0));
    _data[i] = d1-d2;
    _minDist = std::min(_minDist, _data[i]);
    _maxDist = std::max(_maxDist, _data[i]);
  }

  _gdata.clear();

  if (storeGradients) {
    _createGradients();
  }

}

template <class real>
void DtGrid_t<real>::computeDists(bool storeGradients) {

  if (this->empty()) { return; }

  // clear min/max dist
  _minDist = DT_INF;
  _maxDist = -DT_INF;

  // store things near the boundary
  typedef std::map<size_t, real> NearLookup;

  // store whether we are inside
  std::vector<bool> inside(this->_size, false);
  NearLookup near;

  for (size_t i=0; i<this->_size; ++i) {
    if (_data[i] < 0) { 
      inside[i] = true; 
      _data[i] = -_data[i];
    } 
    // now positive
    if (_data[i] != DT_INF) {
      if (_data[i] < 0.87*this->_cellSize) {
        near[i] = _data[i];
        //_data[i] /= this->_cellSize;
        //_data[i] *= _data[i];
        _data[i] = 0;
      } else {
        _data[i] = DT_INF;
      }
    }
  }
  
  _computeEDT();

  for (size_t z=0; z<this->nz(); ++z) {
    for (size_t y=0; y<this->ny(); ++y) {
      for (size_t x=0; x<this->nx(); ++x) {
        size_t idx = this->sub2ind(x,y,z);
        typename NearLookup::const_iterator i = near.find(idx);
        if (i != near.end()) {
          _data[idx] = i->second;
        }
        if (inside[idx]) { _data[idx] = -_data[idx]; }
        _minDist = std::min(_minDist, _data[idx]);
        _maxDist = std::max(_maxDist, _data[idx]);
      }
    }
  }

  _gdata.clear();

  if (storeGradients) {
    _createGradients();
  }

}

template <class real>
void DtGrid_t<real>::_computeEDT() {

  size_t maxdim = this->_dims.max();
  RealArray f(maxdim), ft(maxdim), zz(maxdim+1);
  std::vector<int> v(maxdim);

  // along z
  for (size_t y=0; y<this->ny(); ++y) {
    for (size_t x=0; x<this->nx(); ++x) {
      for (size_t z=0; z<this->nz(); ++z) {
        f[z] = _data[this->sub2ind(x,y,z)];
        assert(f[z] >= 0);
      }
      dt(f, this->nz(), ft, v, zz);
      for (size_t z=0; z<this->nz(); ++z) {
        _data[this->sub2ind(x,y,z)] = ft[z];
      }
    }
  }

  // along y
  for (size_t z=0; z<this->nz(); ++z) {
    for (size_t x=0; x<this->nx(); ++x) {
      for (size_t y=0; y<this->ny(); ++y) {
        f[y] = _data[this->sub2ind(x,y,z)];
        assert(f[y] >= 0);
      }
      dt(f, this->ny(), ft, v, zz);
      for (size_t y=0; y<this->ny(); ++y) {
        _data[this->sub2ind(x,y,z)] = ft[y];
      }
    }
  }

  // along x
  for (size_t z=0; z<this->nz(); ++z) {
    for (size_t y=0; y<this->ny(); ++y) {
      for (size_t x=0; x<this->nx(); ++x) {
        f[x] = _data[this->sub2ind(x,y,z)];
        assert(f[x] >= 0);
      }
      dt(f, this->nx(), ft, v, zz);
      for (size_t x=0; x<this->nx(); ++x) {
        _data[this->sub2ind(x,y,z)] = sqrt(ft[x])*this->_cellSize;
      }
    }
  }

}

template <class real>
void DtGrid_t<real>::_createGradients() {

  _gdata.clear();

  Vec3Array tmp(this->_size);
    
  for (size_t z=0; z<this->nz(); ++z) {
    for (size_t y=0; y<this->ny(); ++y) {
      for (size_t x=0; x<this->nx(); ++x) {
        tmp[this->sub2ind(x,y,z)] = gradient(x,y,z);
      }
    }
  }
  
  tmp.swap(_gdata);

}


template <class real>
static inline vec4_t<real> calcPlane(const vec3_t<real> v[3]) {
  vec3_t<real> normal = vec3_t<real>::cross((v[1]-v[0]), (v[2] - v[0]));
  normal = safeNormalize(normal);
  if (!normal.norm2()) { return vec4f(0); }
  real d = -vec3_t<real>::dot(normal, v[0]);
  return vec4_t<real>(normal[0], normal[1], normal[2], d);
}

template <class real>
static inline real planeDist(const vec3_t<real>& v, const vec4_t<real>& p) {
  return v[0]*p[0] + v[1]*p[1] + v[2]*p[2] + p[3];
}
 
template <class real>
struct HitStruct_t {
  real height;
  int   orient;
};

template <class real>
inline bool operator<(const HitStruct_t<real>& a, const HitStruct_t<real>& b) {
  return a.height > b.height;
}


template <class real>
void DtGrid_t<real>::scanConvert(const TriMesh3& g, bool asHeightmap) {

  _scanConvert(g, 0, asHeightmap);

}

template <class real>
void DtGrid_t<real>::scanConvert(const TriMesh3& g, 
                                 const Transform3& tx, 
                                 bool asHeightmap) {

  _scanConvert(g, &tx, asHeightmap);

}
                     

template <class real>
void DtGrid_t<real>::_scanConvert(const TriMesh3& g, 
                                  const Transform3* ptx, 
                                  bool asHeightmap) {


  throw std::runtime_error("_scanConvert not implemented yet");

}

template <class real>
const real& DtGrid_t<real>::height(const vec2u& s) const {
  return _hmap(s);
}

template <class real>
real& DtGrid_t<real>::height(const vec2u& s) {
  return _hmap(s);
}

template <class real>
const real& DtGrid_t<real>::height(size_t u, size_t v) const {
  return _hmap(u,v);
}

template <class real>
real& DtGrid_t<real>::height(size_t u, size_t v) {
  return _hmap(u,v);
}

template <class real>
vec3_t<real> DtGrid_t<real>::heightVec(const vec3& v) const {
  vec3 rval = v;
  rval[_ax[2]] = height(v);
  return rval;
}

template <class real>
real DtGrid_t<real>::height(const vec3& v) const {
  vec3u s = this->nearestCell(v);
  return height(s[_ax[0]], s[_ax[1]]);
}

template <class real>
vec3_t<real> DtGrid_t<real>::heightVec(size_t u, size_t v) const {
  vec3u s;
  s[_ax[0]] = u;
  s[_ax[1]] = v;
  s[_ax[2]] = 0;
  vec3 rval = this->cellCenter(s);
  rval[_ax[2]] = height(u,v);
  return rval;
}
  

template <class real>
vec3_t<real> DtGrid_t<real>::heightVec(const vec2u& s) const {
  return heightVec(s[0], s[1]);
}

template <class real>
void DtGrid_t<real>::recomputeExtents() {

  for (size_t z=0; z<this->nz(); ++z) {
    for (size_t y=0; y<this->ny(); ++y) {
      for (size_t x=0; x<this->nx(); ++x) {
        vec3u s(x,y,z);
        size_t i = this->sub2ind(s);
        _minDist = std::min(_data[i], _minDist);
        _maxDist = std::max(_data[i], _maxDist);
      }
    }
  }

  _hmap.recomputeExtents();

}

template <class real>
real DtGrid_t<real>::minHeight() const { return _hmap.minHeight(); }

template <class real>
real DtGrid_t<real>::maxHeight() const { return _hmap.maxHeight(); }

template <class real>
const vec3u& DtGrid_t<real>::referenceAxes() const {
  return _ax;
}

template <class real>
typename DtGrid_t<real>::Axis DtGrid_t<real>::referenceAxis() const {
  return Axis(_ax[2]);
}


template <class real>
class LineMin3_t {
public:
  const vec3u& dims;
  const std::vector<real>& data;
  real& minVal;
  vec3u& minPixel;
  
  LineMin3_t(const vec3u& sz,
             const std::vector<real>& d,
             real& mv, 
             vec3u& mp):
    
    dims(sz), data(d), minVal(mv), minPixel(mp) {}
  
  void operator()(const vec3u& pixel) {
    real d = data[vec3u::sub2ind(dims,pixel)];            
    if (d < minVal) {                           
      minVal = d;                               
      minPixel = pixel;                         
    }                                           
  }
           
  
};

template <class real>
real DtGrid_t<real>::lineMin(const vec3& v0, 
                             const vec3& v1, 
                             vec3& vmin, 
                             vec3& gmin) const {

  vec3u s1 = this->nearestCell(v0);
  vec3u s2 = this->nearestCell(v1);
  vec3u smin;

  real minVal = lineMin(s1, s2, smin);

  // TODO: project onto line
  vmin = this->cellCenter(smin);
  gmin = this->gradient(smin);

  return minVal;
  


}

template <class real>
real DtGrid_t<real>::lineMin(const vec3u& s1, const vec3u& s2, vec3u& minPixel) const {

  minPixel = s1;
  if (this->empty()) { return 0; }
  
  real minVal = _data[this->sub2ind(s1)];

  LineMin3_t<real> lm(this->_dims, _data, minVal, minPixel);
  
  bresenham3D(s1, s2, lm);

  return minVal;

}

template <class Tval>
void get(std::istream& istr, Tval& value) {
  if (!istr.read((char*)&value, sizeof(value))) {
    throw std::runtime_error("error reading!\n");
  }
}

template <class Tval>
void put(std::ostream& ostr, const Tval& value) {
  if (!ostr.write((const char*)&value, sizeof(value))) {
    throw std::runtime_error("error writing!\n");
  }
}

template <class real>
bool DtGrid_t<real>::load(const char* filename, bool storeGradients) {

  std::ifstream istr(filename);
  if (!istr.is_open()) { return false; }

  clear();

  try {
    
    for (int i=0; i<3; ++i) { get(istr, this->_dims[i]); }
    get(istr, this->_size);
    
    size_t hsize;

    for (int i=0; i<3; ++i) { get(istr, _ax[i]); }
    get(istr, hsize);
    
    for (int i=0; i<3; ++i) { get(istr, this->_origin[i]); }
    get(istr, this->_cellSize);
    
    if (this->_size != this->_dims.prod()) { 
      clear();
      return false; 
    }
    if (hsize != this->_dims[_ax[0]] * this->_dims[_ax[1]]) {
      clear();
      return false;
    }
    
    _data.resize(this->_size);

    _hmap.resize(this->_dims[_ax[0]], this->_dims[_ax[1]], 
                 this->_cellSize,
                 vec2(this->_origin[_ax[0]], this->_origin[_ax[1]]));
    
    for (size_t i=0; i<this->_size; ++i) {
      get(istr, _data[i]);
    }
    
    for (size_t i=0; i<_hmap.size(); ++i) {
      get(istr, _hmap[i]);
    }
    
    get(istr, _minDist);
    get(istr, _maxDist);
    _hmap.recomputeExtents();

  } catch (...) {
    
    clear();
    return false;

  }

  if (storeGradients) {
    _createGradients();
  }

  return true;

}

template <class real>
void DtGrid_t<real>::save(const char* filename) const { 

  std::ofstream ostr(filename);
  if (!ostr.is_open()) { return; }


  for (int i=0; i<3; ++i) { put(ostr, this->_dims[i]); }
  put(ostr, this->_size);

  for (int i=0; i<3; ++i) { put(ostr, _ax[i]); }

  put(ostr, _hmap.size());

  for (int i=0; i<3; ++i) { put(ostr, this->_origin[i]); }
  put(ostr, this->_cellSize);

  for (size_t i=0; i<this->_size; ++i) {
    put(ostr, _data[i]);
  }

  for (size_t i=0; i<_hmap.size(); ++i) {
    put(ostr, _hmap[i]);
  }

  put(ostr, _minDist);
  put(ostr, _maxDist);

}

template <class real>
const HeightMap_t<real>& DtGrid_t<real>::heightMap() const { return _hmap; }

template <class real>
HeightMap_t<real>& DtGrid_t<real>::heightMap() { return _hmap; }

// TODO: check this!
/*
template <class real>
void DtGrid_t<real>::processPointCloudAsHeightmap(const std::vector<vec3>& points,
                                          size_t minBinSize) {
  
  HeightMap mmap;
  typename HeightMap::BinStorage bins;

  size_t nu = _hmap.nx();
  size_t nv = _hmap.ny();
  size_t nc = this->_dims[_ax[2]];
  
  assert( nu == this->_dims[_ax[0]] );
  assert( nv == this->_dims[_ax[1]] );

  mmap.resize(nu, nv, _hmap.cellSize(), _hmap.origin());
  mmap.binPoints(bins, points, _ax);
  mmap.medianMap(bins, minBinSize);

  // for each u,v in new heightmap
  for (size_t v=0; v<nv; ++v) {
    for (size_t u=0; u<nu; ++u) {

      // get height of this u,v location
      real hm = mmap(u,v);
      _hmap(u,v) = std::max(_hmap(u,v), hm);

      // get min and max of neighbors (including this)
      real nmin = hm;
      real nmax = hm;

      const int du[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
      const int dv[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };

      for (int n=0; n<8; ++n) {
        real hn = mmap.correctHeight(u+du[n], v+dv[n], hm);
        nmin = std::min(nmin, hn);
        nmax = std::max(nmax, hn);
      }

      // for each column in the distance field
      for (size_t c=0; c<nc; ++c) {

        vec3u s;

        s[_ax[0]] = u;
        s[_ax[1]] = v;
        s[_ax[2]] = c;

        vec3 cc = this->cellCenter(s);

        real hc = cc[_ax[2]];

        real hdiff = (hc - hm) / this->_cellSize;
        
        real d = DT_INF;

        if (fabs(hdiff) < 0.5) {
          d = hdiff;
        } else if (hc > hm && hc < nmax) {
          d = -0.5;
        } else if (hc < hm && hc > nmin) {
          d = 0.5;
        } else if (hc < hm) {
          d = -DT_INF;
        }

        // choose magnitude of smaller and negate if either negative
        size_t idx = this->sub2ind(s);
        real& cur = _data[idx];

        bool neg = (d < 0 || cur < 0);
        cur = std::min(fabs(d), fabs(cur)) * (neg ? -1 : 1);
        
      } // foreach cell in u,v column
    } // foreach v
  } // foreach u

  _hmap.recomputeExtents();

}
*/

template <class real>
void DtGrid_t<real>::computeDistsFromHeightMap(bool storeGradients) {


  size_t nu = _hmap.nx();
  size_t nv = _hmap.ny();
  size_t nc = this->_dims[_ax[2]];
  
  assert( nu == this->_dims[_ax[0]] );
  assert( nv == this->_dims[_ax[1]] );

  // for each u,v in heightmap
  for (size_t v=0; v<nv; ++v) {
    for (size_t u=0; u<nu; ++u) {

      // get the height
      real h = _hmap(u,v);

      // get min and max of neighbors (including this)
      real nmin =  DT_INF;
      real nmax = -DT_INF;

      const int du[8] = { -1,  0,  1, -1, 1, -1, 0, 1 };
      const int dv[8] = { -1, -1, -1,  0, 0,  1, 1, 1 };

      for (int n=0; n<8; ++n) {
        size_t uu = u + du[n];
        size_t vv = v + dv[n];
        if (uu < nu && vv < nv) {
          real hn = _hmap(u+du[n], v+dv[n]);
          nmin = std::min(nmin, hn);
          nmax = std::max(nmax, hn);
        }
      }

      // for each column in the distance field
      for (size_t c=0; c<nc; ++c) {

        vec3u s;

        s[_ax[0]] = u;
        s[_ax[1]] = v;
        s[_ax[2]] = c;

        // get the 3d position of u,v,c
        vec3 cc = this->cellCenter(s);

        // get its height
        real hc = cc[_ax[2]];

        real hdiff = (hc - h);

        real d = DT_INF;

        if (fabs(hdiff) < 0.5*this->_cellSize) {
          // if within cell of height of this column
          d = hdiff;
        } else if (hdiff > 0 && hc < nmax) { 
          // above terrain, below neighbor
          d = 0.5 * this->_cellSize;
        } else if (hdiff < 0 && hc > nmin) {
          // below terrain, above neighbor
          d = -0.5 * this->_cellSize;
        } else if (hdiff < 0) {
          // below terrain
          d = -DT_INF;
        }

        _data[this->sub2ind(s)] = d;

      }


    }
  }

  computeDists(storeGradients);
  
}

template class DtGrid_t<float>;
template class DtGrid_t<double>;


