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

#include "HeightMap.h"
#include <float.h>
#include "Box2.h"
#include <assert.h>
#include <fstream>
#include <algorithm>

template <> const float HeightMap_t<float>::INVALID_HEIGHT = -FLT_MAX;
template <> const double HeightMap_t<double>::INVALID_HEIGHT = -DBL_MAX;

template <class real>
Box3_t<real> HeightMap_t<real>::computeBBox() const {
  Box2 box2 = this->bbox();
  return Box3(vec3(box2.p0, _minHeight),
	      vec3(box2.p1, _maxHeight));
}

template <class real>
HeightMap_t<real>::HeightMap_t() {
  clear();
}

template <class real>
void HeightMap_t<real>::clear() {
  this->_clear();
  _hdata.clear();
  _minHeight = INVALID_HEIGHT;
  _maxHeight = INVALID_HEIGHT;
}

template <class real>
void HeightMap_t<real>::resize(size_t nx, size_t ny, 
                               real cellSize, 
                               const vec2& origin) {
  clear();
  this->_resize(nx, ny, cellSize, origin);
  _hdata.resize(this->size(), INVALID_HEIGHT);
}

template <class real>
void HeightMap_t<real>::resize(const vec2& min,
            const vec2& max,
            real cellSize) {
  clear();
  this->_resize(min, max, cellSize);
  _hdata.resize(this->size(), INVALID_HEIGHT);
}

template <class real>
real HeightMap_t<real>::correctHeight(real h, real defaultValue) {
  return (h == INVALID_HEIGHT) ? defaultValue : h;
}

template <class real>
real HeightMap_t<real>::correctHeight(size_t x, size_t y, real defaultValue) const {
  if (x >= this->nx() || y >= this->ny()) { return defaultValue; }
  return correctHeight(_hdata[this->sub2ind(x,y)], defaultValue);
}

template <class real>
real HeightMap_t<real>::correctHeight(const vec2u& s, real defaultValue) const {
  if (s.x() >= this->nx() || s.y() >= this->ny()) { return defaultValue; }
  return correctHeight(_hdata[this->sub2ind(s)], defaultValue);
}

template <class real>
const real& HeightMap_t<real>::operator[](size_t idx) const { return _hdata[idx]; }

template <class real>
real& HeightMap_t<real>::operator[](size_t idx) { return _hdata[idx]; }

template <class real>
const real& HeightMap_t<real>::operator()(const vec2u& s) const { return _hdata[this->sub2ind(s)]; }

template <class real>
real& HeightMap_t<real>::operator()(const vec2u& s) { return _hdata[this->sub2ind(s)]; }

template <class real>
const real& HeightMap_t<real>::operator()(size_t x, size_t y) const {
  return _hdata[this->sub2ind(x,y)];
}

template <class real>
real& HeightMap_t<real>::operator()(size_t x, size_t y) {
  return _hdata[this->sub2ind(x,y)];
}

template <class real>
void HeightMap_t<real>::recomputeExtents() {
  computeExtents(vec2u(0,0), vec2u(this->nx(), this->ny()), _minHeight, _maxHeight);
}

template <class real>
void HeightMap_t<real>::computeExtents(const vec2u& s0, 
				       const vec2u& s1,
				       real& hmin,
				       real& hmax) const {

  hmin=INVALID_HEIGHT;
  hmax=INVALID_HEIGHT;

  for (size_t y=s0.y(); y<s1.y(); ++y) {
    for (size_t x=s0.x(); x<s1.x(); ++x) {
      real h = (*this)(x,y);
      if (h != INVALID_HEIGHT) {
	if (hmin == INVALID_HEIGHT || h < hmin) { hmin = h; }
	if (hmax == INVALID_HEIGHT || h > hmax) { hmax = h; }
      }
    }
  }


}


template <class real>
bool HeightMap_t<real>::load(const char* filename) {

  std::ifstream istr(filename);
  if (!istr) { return false; }

  vec2u d;
  vec2 o;
  real cs;
  size_t s;
  RealArray htmp;

  istr.read((char*)&(d[0]), sizeof(d[0]));
  istr.read((char*)&(d[1]), sizeof(d[1]));
  istr.read((char*)&(s), sizeof(s));

  if (d.prod() != s) { return false; }
  
  istr.read((char*)&(o[0]), sizeof(o[0]));
  istr.read((char*)&(o[1]), sizeof(o[1]));
  istr.read((char*)&cs, sizeof(cs));

  htmp.resize(s);
  
  for (size_t i=0; i<s; ++i) {
    istr.read((char*)&(htmp[i]), sizeof(real));
  }

  this->_dims = d;
  this->_size = s;
  this->_origin = o;
  this->_cellSize = cs;

  _hdata.swap(htmp);

  recomputeExtents();
  
  return true;

}

template <class real>
void HeightMap_t<real>::save(const char* filename) const {
  std::ofstream ostr(filename);
  if (!ostr) { return; }
  ostr.write((const char*)&(this->_dims[0]), sizeof(this->_dims[0]));
  ostr.write((const char*)&(this->_dims[1]), sizeof(this->_dims[1]));
  ostr.write((const char*)&(this->_size), sizeof(this->_size));
  ostr.write((const char*)&(this->_origin[0]), sizeof(this->_origin[0]));
  ostr.write((const char*)&(this->_origin[1]), sizeof(this->_origin[1]));
  ostr.write((const char*)&(this->_cellSize), sizeof(this->_cellSize));
  for (size_t i=0; i<this->_size; ++i) {
    ostr.write((const char*)&(_hdata[i]), sizeof(real));
  }

}

template <class real>
real HeightMap_t<real>::minHeight() const { return _minHeight; }

template <class real>
real HeightMap_t<real>::maxHeight() const { return _maxHeight; }

template <class real>
void HeightMap_t<real>::clearBins(BinStorage& bins) {
  for (size_t i=0; i<bins.size(); ++i) { delete bins[i]; }
  bins.clear();
}

template <class real>
Box2_t<real>
HeightMap_t<real>::boundPoints(const std::vector<vec3>& points,
                               const vec3u& ax) {
  

  Box2 box;

  // first, compute the bounds
  for (size_t j=0; j<points.size(); ++j) {
    const vec3& pj = points[j];
    vec2 pj2(pj[ax[0]], pj[ax[1]]);
    box.addPoint(pj2);
  }

  return box;

}

template <class real>
void HeightMap_t<real>::binPoints(BinStorage& bins,
                          const std::vector<vec3>& points,
                          const vec3u& ax) const {

  if (bins.size() != this->size()) {
    clearBins(bins);
    bins.resize(this->size(),0);
  }
  
  for (size_t j=0; j<points.size(); ++j) {
    const vec3& pj = points[j];
    vec2 pj2(pj[ax[0]], pj[ax[1]]);
    real pjz = pj[ax[2]];
    vec2u s = this->nearestCell(pj2);
    size_t idx = this->sub2ind(s);
    if (!bins[idx]) {
      bins[idx] = new RealArray(1, pjz);
    } else {
      bins[idx]->push_back(pjz);
    }
  }

}

template <class real>
void HeightMap_t<real>::medianMap(BinStorage& bins, size_t minCount) {

  assert(bins.size() == this->size());

  // get median for each bin
  for (size_t i=0; i<this->size(); ++i) {

    if (!bins[i]) { continue; }

    size_t binsize = bins[i]->size();
    if (binsize < minCount) { continue; }

    typename RealArray::iterator begin = bins[i]->begin();
    typename RealArray::iterator nth = begin + binsize/2;
    typename RealArray::iterator end = bins[i]->end();
    std::nth_element(begin, nth, end);

    _hdata[i] = std::max(_hdata[i], *nth);

  }

  // cleanup
  clearBins(bins);
  
  recomputeExtents();

}

template <class real>
vec2_t<real> HeightMap_t<real>::slope(size_t x, size_t y) const {

  const vec2u& d = this->_dims;

  if (x >= d[0] || y >= d[1]) { return vec2(0); }
  
  real h11 = (*this)(x,y);

  vec2 slope(0);

  {

    real h01 = x ? (*this)(x-1, y) : INVALID_HEIGHT;
    real h21 = x+1 < d[0] ? (*this)(x+1, y) : INVALID_HEIGHT;
    
    if (h21 != INVALID_HEIGHT) {
      if (h01 != INVALID_HEIGHT) {
        slope[0] = (h21 - h01) / (2*this->_cellSize);
      } else if (h11 != INVALID_HEIGHT) {
        slope[0] = (h21 - h11) / (this->_cellSize);
      }
    } else if (h01 != INVALID_HEIGHT && h11 != INVALID_HEIGHT) {
      slope[0] = (h11 - h01) / (this->_cellSize);
    }

  }

  {

    real h10 = y ? (*this)(x, y-1) : INVALID_HEIGHT;
    real h12 = y+1 < d[1] ? (*this)(x, y+1) : INVALID_HEIGHT;

    if (h12 != INVALID_HEIGHT) {
      if (h10 != INVALID_HEIGHT) {
        slope[1] = (h12 - h10) / (2*this->_cellSize);
      } else if (h11 != INVALID_HEIGHT) {
        slope[1] = (h12 - h11) / (this->_cellSize);
      }
    } else if (h10 != INVALID_HEIGHT && h11 != INVALID_HEIGHT) {
      slope[1] = (h11 - 10) / (this->_cellSize);
    }
  }
      
  return slope;

}

template <class real>
vec3_t<real> HeightMap_t<real>::normal(size_t x, size_t y) const {
  vec2 s = slope(x,y);
  vec3 tx(1, 0, s[0]);
  vec3 ty(0, 1, s[1]);
  vec3 n = vec3::cross(tx,ty);
  return n / n.norm();
}

template <class real>
vec2_t<real> HeightMap_t<real>::slope(const vec2u& s) const {
  return slope(s[0], s[1]);
}

template <class real>
vec3_t<real> HeightMap_t<real>::normal(const vec2u& s) const {
  return normal(s[0], s[1]);
}

template <class real>
void HeightMap_t<real>::generateMesh(TriMesh3& mesh, real dropEdges, const vec3u& ax) const {
  
  mesh.clear();
  mesh.verts.resize(this->size());

  for (size_t i=0; i<this->size(); ++i) {
    vec3 vi;
    vec2 cc = this->cellCenter(i);
    real h = _hdata[i];
    if (dropEdges > 0) { h = correctHeight(h, _minHeight-dropEdges); }
    vi[ax[0]] = cc[0];
    vi[ax[1]] = cc[1];
    vi[ax[2]] = h;
    mesh.verts[i] = vi;
  }

  //////////////////////////////////////////////////
  //
  // 01 *-------* 11
  //    |\     /|
  //    | \   / |
  //    |  \ /  |
  //    |   .   |
  //    |  / \  |
  //    | /   \ |
  //    |/     \|
  // 00 *-------* 10       
  //    
  //////////////////////////////////////////////////


  for (size_t v=1; v<this->ny(); ++v) {
    for (size_t u=1; u<this->nx(); ++u) {

      size_t i00 = this->sub2ind(u-1, v-1);
      size_t i10 = this->sub2ind(u,   v-1);
      size_t i01 = this->sub2ind(u-1, v  );
      size_t i11 = this->sub2ind(u,   v  );

      real h00 = _hdata[i00];
      real h10 = _hdata[i10];
      real h01 = _hdata[i01];
      real h11 = _hdata[i11];

      if (dropEdges > 0) {
        if (h00 != INVALID_HEIGHT || h10 != INVALID_HEIGHT ||
            h01 != INVALID_HEIGHT || h11 != INVALID_HEIGHT) {
          h00 = correctHeight(h00, _minHeight-dropEdges);
          h10 = correctHeight(h10, _minHeight-dropEdges);
          h01 = correctHeight(h01, _minHeight-dropEdges);
          h11 = correctHeight(h11, _minHeight-dropEdges);
        }
      }

      if (h00 != INVALID_HEIGHT && h10 != INVALID_HEIGHT &&
          h01 != INVALID_HEIGHT && h11 != INVALID_HEIGHT) {
        if (fabs(h00 - h11) < fabs(h01 - h10)) {
          mesh.addTriangle(i00, i10, i11);
          mesh.addTriangle(i00, i11, i01);
        } else {
          mesh.addTriangle(i01, i00, i10);
          mesh.addTriangle(i01, i10, i11);
        }
      } else if (h00 != INVALID_HEIGHT && h11 != INVALID_HEIGHT) {
        if (h10 != INVALID_HEIGHT) {
          mesh.addTriangle(i00, i10, i11);
        } 
        if (h01 != INVALID_HEIGHT) {
          mesh.addTriangle(i00, i11, i01);
        }
      } else if (h10 != INVALID_HEIGHT && h01 != INVALID_HEIGHT) {
        if (h00 != INVALID_HEIGHT) {
          mesh.addTriangle(i01, i00, i10);
        } 
        if (h11 != INVALID_HEIGHT) {
          mesh.addTriangle(i01, i10, i11);
        }
      }
      
    }
  }

}

/*
template <class real>
Box3_t<real> HeightMap_t<real>::computeBBox(const Split& s) const {
  const vec2& o = this->_origin;
  const real& cs = this->_cellSize;
  return Box3( vec3( o + vec2(s.s0.x(), s.s0.y()) * cs, s.hmin ),
	       vec3( o + vec2(s.s1.x(), s.s1.y()) * cs, s.hmax ) );
}
*/

template <class real>
void HeightMap_t<real>::subdivide(SplitArray& splits) const {
  _subdivide(splits, size_t(-1), 0, vec2u(0,0), vec2u(this->nx(), this->ny()));
}

template <class real>
void HeightMap_t<real>::_subdivide(SplitArray& splits,
				   size_t parent_index,
				   size_t which_child,
				   const vec2u& s0,
				   const vec2u& s1) const {

  assert( s1.x() > s0.x() );
  assert( s1.y() > s0.y() );
  assert( s0.x() < this->nx() );
  assert( s0.y() < this->ny() );
  assert( s1.x() <= this->nx() );
  assert( s1.y() <= this->ny() );
  
  vec2u size = s1-s0;
  
  size_t cur_index = splits.size();
  Split cur;
  
  cur.s0 = s0;
  cur.s1 = s1;

  real hmin, hmax;
  computeExtents(s0, s1, hmin, hmax);

  if (hmin == INVALID_HEIGHT) { 
    return; 
  }

  cur.box.p0 = vec3(this->cellCenter(s0),            hmin);
  cur.box.p1 = vec3(this->cellCenter(s1-vec2u(1,1)), hmax);

  cur.parent_index = parent_index;
  cur.child_index[0] = size_t(-1);
  cur.child_index[1] = size_t(-1);

  splits.push_back(cur);

  if (parent_index != size_t(-1)) {
    splits[parent_index].child_index[which_child] = cur_index;
  }

  if (hmin == hmax) {
    return;
  }

  if (size.x() >= size.y()) { 
    // split along x
    size_t xmid = s0[0] + size[0]/2;
    _subdivide(splits, cur_index, 0, s0, vec2u(xmid, s1.y()));
    _subdivide(splits, cur_index, 1, vec2u(xmid, s0.y()), s1);
  } else {
    // split along y
    size_t ymid = s0[1] + size[1]/2;
    _subdivide(splits, cur_index, 0, s0, vec2u(s1.x(), ymid));
    _subdivide(splits, cur_index, 1, vec2u(s0.x(), ymid), s1);
  }

  const Split& cur_post = splits[cur_index];
  assert( cur_post.parent_index == parent_index );
  assert( cur_post.child_index[0] == size_t(-1) ||
	  splits[cur_post.child_index[0]].parent_index == cur_index );
  assert( cur_post.child_index[1] == size_t(-1) ||
	  splits[cur_post.child_index[1]].parent_index == cur_index );

}



template class HeightMap_t<float>;
template class HeightMap_t<double>;
