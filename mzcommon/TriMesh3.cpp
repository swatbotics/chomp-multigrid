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

#include "TriMesh3.h"
#include "glstuff.h"
#include "strutils.h"
#include <stdio.h>
#include <stdexcept>
#include <string>
//#include <util/gzstream.h>
#include <fstream>
#include <stdlib.h>
#include <map>
#include <assert.h>

// TODO: pretty normals!

const size_t npos = size_t(-1);

template <class real>
void TriMesh3_t<real>::swap(TriMesh3_t<real>& g) {
  verts.swap(g.verts);
  normals.swap(g.normals);
  faces.swap(g.faces);
}

template <class real>
TriMesh3_t<real>::TriMesh3_t() { }

template <class real>
void TriMesh3_t<real>::clear() {
  verts.clear();
  normals.clear();
  faces.clear();
}

template <class real>
vec3_t<real> TriMesh3_t<real>::centroid(size_t faceidx) const {
  const Tri& f = faces[faceidx];
  return (verts[f.vidx[0]] + verts[f.vidx[1]] + verts[f.vidx[2]]) / 3;
}


template <class real>
bool TriMesh3_t<real>::empty() const {
  return (verts.empty() || faces.empty());
}

template <class real>
size_t TriMesh3_t<real>::addVertex(real x, real y, real z) {
  verts.push_back(vec3(x,y,z));
  return verts.size() - 1;
}

template <class real>
size_t TriMesh3_t<real>::addNormal(real x, real y, real z) {
  normals.push_back(vec3(x,y,z));
  return normals.size() - 1;
}


template <class real>
size_t TriMesh3_t<real>::addVertex(const vec3& v) {
  verts.push_back(v);
  return verts.size() - 1;
}

template <class real>
size_t TriMesh3_t<real>::addNormal(const vec3& v) {
  normals.push_back(v);
  return normals.size() - 1;
}

template <class real>
void TriMesh3_t<real>::addMesh(const TriMesh3_t<real>& g) {

  size_t voffs = verts.size();
  size_t noffs = normals.size();

  verts.resize(voffs + g.verts.size());
  normals.resize(noffs + g.normals.size());
               
  for (size_t i=0; i<g.verts.size(); ++i) {
    verts[voffs+i] = g.verts[i];
  }
  for (size_t i=0; i<g.normals.size(); ++i) {
    normals[noffs+i] = g.normals[i];
  }

  for (size_t i=0; i<g.faces.size(); ++i) {
    faces.push_back(g.faces[i]);
    Tri& f = faces.back();
    for (int j=0; j<3; ++j) {
      if (f.vidx[j] != npos) { f.vidx[j] += voffs; }
      if (f.nidx[j] != npos) { f.nidx[j] += noffs; }
    }
  }

}

template <class real>
void TriMesh3_t<real>::addMesh(const TriMesh3_t<real>& g,
                                   const Transform3& xform) {

  size_t voffs = verts.size();
  size_t noffs = normals.size();

  verts.resize(voffs + g.verts.size());
  normals.resize(noffs + g.normals.size());
               
  for (size_t i=0; i<g.verts.size(); ++i) {
    verts[voffs+i] = xform.transformFwd(g.verts[i]);
  }
  for (size_t i=0; i<g.normals.size(); ++i) {
    normals[noffs+i] = xform.rotFwd() * g.normals[i];
  }

  for (size_t i=0; i<g.faces.size(); ++i) {
    faces.push_back(g.faces[i]);
    Tri& f = faces.back();
    for (int j=0; j<3; ++j) {
      if (f.vidx[j] != npos) { f.vidx[j] += voffs; }
      if (f.nidx[j] != npos) { f.nidx[j] += noffs; }
    }
  }

}

template <class real>
void TriMesh3_t<real>::addMesh(const TriMesh3_t<real>& g,
                                   const mat4& m) {

  size_t voffs = verts.size();
  size_t noffs = normals.size();

  verts.resize(voffs + g.verts.size());
  normals.resize(noffs + g.normals.size());
               
  mat4 m2(m);
  m2(0,3) = m2(1,3) = m2(2,3) = 0.0f;

  for (size_t i=0; i<g.verts.size(); ++i) {
    verts[voffs+i] = m.mult3D(g.verts[i]);
  }
  for (size_t i=0; i<g.normals.size(); ++i) {
    normals[noffs+i] = m2.mult3D(g.normals[i]);
  }

  for (size_t i=0; i<g.faces.size(); ++i) {
    faces.push_back(g.faces[i]);
    Tri& f = faces.back();
    for (int j=0; j<3; ++j) {
      if (f.vidx[j] != npos) { f.vidx[j] += voffs; }
      if (f.nidx[j] != npos) { f.nidx[j] += noffs; }
    }
  }

}



template <class real>
size_t TriMesh3_t<real>::addTriangle(size_t i,
                                     size_t j,
                                     size_t k,
                                     bool flip,
                                     size_t fdata) {
  
  return addTriangle(i, j, k, npos, npos, npos, flip, fdata);

}

template <class real>
size_t TriMesh3_t<real>::addTriangle(size_t i, size_t j, size_t k, 
                                     size_t a, size_t b, size_t c,
                                     bool flip, size_t fdata) {
  
  if (flip) {
    faces.push_back(Tri(k,j,i,c,b,a,fdata));
  } else {
    faces.push_back(Tri(i,j,k,a,b,c,fdata));
  }

  return faces.size() - 1;

}

TriMesh3Base::Tri::Tri(const size_t v[3], const size_t n[3], size_t f) {

  vidx[0] = v[0];
  vidx[1] = v[1];
  vidx[2] = v[2];

  nidx[0] = n[0];
  nidx[1] = n[1];
  nidx[2] = n[2];

  fdata = f;

}

TriMesh3Base::Tri::Tri(size_t i, size_t j, size_t k,
                       size_t a, size_t b, size_t c, size_t f) {

  vidx[0] = i;
  vidx[1] = j;
  vidx[2] = k;

  nidx[0] = a;
  nidx[1] = b;
  nidx[2] = c;

  fdata = f;

}

/*
TriMesh3_t<real>::Tri::Tri() { }
*/

template <class real>
TriMesh3_t<real> TriMesh3_t<real>::box(real length,
			       real width,
			       real height) {

  TriMesh3_t<real> g;
  
  real l2 = length * 0.5;
  real w2 = width * 0.5;
  real h2 = height * 0.5;
  
  g.addVertex( l2,  w2,  h2); // 0
  g.addVertex(-l2,  w2,  h2); // 1
  g.addVertex( l2,  w2, -h2); // 2
  g.addVertex(-l2,  w2, -h2); // 3
  g.addVertex( l2, -w2,  h2); // 4
  g.addVertex(-l2, -w2,  h2); // 5
  g.addVertex( l2, -w2, -h2); // 6
  g.addVertex(-l2, -w2, -h2); // 7

  g.addTriangle(1, 0, 2);
  g.addTriangle(1, 2, 3);
  g.addTriangle(0, 4, 6);
  g.addTriangle(0, 6, 2);
  g.addTriangle(7, 4, 5);
  g.addTriangle(7, 6, 4);
  g.addTriangle(3, 5, 1);
  g.addTriangle(3, 7, 5);
  g.addTriangle(5, 0, 1);
  g.addTriangle(5, 4, 0);
  g.addTriangle(6, 7, 3);
  g.addTriangle(6, 3, 2);

  return g;

}
		       

template <class real>
TriMesh3_t<real> TriMesh3_t<real>::cappedBox(real length,
		     real width,
		     real height,
		     int slices) {

  TriMesh3_t<real> g;

  real l2 = length * 0.5;
  real w2 = width * 0.5;
  real h2 = height * 0.5;

  if (slices < 2) { slices = 2; }

  for (int sy=1; sy>=-1; sy-=2) {

    // center edge
    size_t mf = g.addVertex(l2, sy*w2, 0);
    size_t mr = g.addVertex(-l2, sy*w2, 0);

    // top edge
    size_t uf = g.addVertex(l2, sy*w2, h2);
    size_t ur = g.addVertex(-l2, sy*w2, h2);
    
    // rear fan
    for (int i=1; i<slices; ++i) {
      real angle = i * M_PI / slices;
      g.addVertex(-l2-h2*sin(angle), sy*w2, h2*cos(angle));
    }

    // bottom edge
    size_t lr = g.addVertex(-l2, sy*w2, -h2);
    size_t lf = g.addVertex(l2, sy*w2, -h2);

    // front fan
    for (int i=1; i<slices; ++i) {
      real angle = i * M_PI / slices;
      g.addVertex(l2+h2*sin(angle), sy*w2, -h2*cos(angle));
    }

    bool flip = (sy == -1);

    g.addTriangle(ur, uf, mf, flip);
    g.addTriangle(ur, mf, mr, flip);
    g.addTriangle(mr, mf, lf, flip);
    g.addTriangle(mr, lf, lr, flip);
    
    for (int i=0; i<slices; ++i) {
      g.addTriangle(mr, ur+i+1, ur+i, flip);
      if (i == slices-1) {
	g.addTriangle(mf, uf, lf+i, flip);
      } else {
	g.addTriangle(mf, lf+i+1, lf+i, flip);
      }
    }
    
  }

  // go around top
  size_t nv = (g.verts.size() / 2);
  
  for (size_t i=0; i<nv-2; ++i) {
    size_t next = (i+1)%(nv-2);
    g.addTriangle(2+i, 2+next, 2+i+nv);
    g.addTriangle(2+next, 2+next+nv, 2+i+nv);
  }

  return g;

}

template <class real>
TriMesh3_t<real> TriMesh3_t<real>::turnBox(real length, 
		   real extra,
		   real width,
		   real height,
		   int slices1,
		   int slices2) {
  
  TriMesh3_t<real> g;

  real l2 = length * 0.5;
  real w2 = width * 0.5;
  real h2 = height * 0.5;

  if (slices1 < 2) { slices1 = 2; }
  if (slices2 < 2) { slices2 = 2; }

  for (int sy=1; sy>=-1; sy-=2) {

    // center front
    size_t mf = g.addVertex(l2, sy*w2, 0);

    // top edge
    size_t uf = g.addVertex(l2, sy*w2, h2);
    size_t ur = g.addVertex(-l2, sy*w2, h2);
    
    // rear fan
    for (int i=1; i<=slices2; ++i) {
      real angle = i * M_PI * 0.5 / slices2;
      g.addVertex(-l2-extra*sin(angle), sy*w2, -h2 + height*cos(angle));
    }

    // bottom edge
    size_t lr = g.addVertex(-l2, sy*w2, -h2);
    size_t lf = g.addVertex(l2, sy*w2, -h2);

    // front fan
    for (int i=1; i<slices1; ++i) {
      real angle = i * M_PI / slices1;
      g.addVertex(l2+h2*sin(angle), sy*w2, -h2*cos(angle));
    }

    bool flip = (sy == -1);

    g.addTriangle(mf, ur, uf, flip);
    g.addTriangle(mf, lr, ur, flip);
    g.addTriangle(mf, lf, lr, flip);

    for (int i=0; i<slices2; ++i) {
      g.addTriangle(lr, ur+i+1, ur+i, flip);
    }
    
    // front fan
    for (int i=0; i<slices1; ++i) {
      if (i == slices1-1) {
	g.addTriangle(mf, uf, lf+i, flip);
      } else {
	g.addTriangle(mf, lf+i+1, lf+i, flip);
      }
    }
    
  }

  // go around top
  size_t nv = (g.verts.size() / 2);
  
  for (size_t i=0; i<nv-1; ++i) {
    size_t next = (i+1)%(nv-1);
    g.addTriangle(1+i, 1+next, 1+i+nv);
    g.addTriangle(1+next, 1+next+nv, 1+i+nv);
  }

  return g;

}

template <class real>
TriMesh3_t<real> TriMesh3_t<real>::cylinder(real r,
		    real h,
		    int slices) {

  TriMesh3_t<real> g;

  if (slices < 4) { slices = 4; }

  for (int sz=0; sz<2; ++sz) {
    size_t ctr = g.addVertex(0, 0, sz*h);
    for (int i=0; i<slices; ++i) {
      real angle = i * 2 * M_PI / slices;
      size_t vidx = g.addVertex(r*cos(angle), r*sin(angle), sz*h);
      if (i == slices - 1) {
	g.addTriangle(ctr, vidx, ctr + 1, !sz);
      } else {
	g.addTriangle(ctr, vidx, vidx + 1, !sz);
      }
    }
  }

  // go around sides
  size_t nv = (g.verts.size() / 2);
  
  for (size_t i=0; i<nv-1; ++i) {
    size_t next = (i+1)%(nv-1);
    g.addTriangle(1+i, 1+next, 1+i+nv);
    g.addTriangle(1+next, 1+next+nv, 1+i+nv);
  }

  return g;

}

template <class real>
TriMesh3_t<real> TriMesh3_t<real>::sphere(real r,
		  int slices,
		  int stacks) {

  TriMesh3_t<real> g;
  if (slices < 3) { slices = 3; }
  if (stacks < 2) { stacks = 2; }

  for (int i=0; i<=stacks; ++i) {
    if (i == 0) {

      g.addVertex(0, 0, r);
      //std::cerr << "added vert " << 0 << " = " << g.verts[0] << "\n";

    } else if (i == stacks) {

      size_t end = g.addVertex(0, 0, -r);
      //std::cerr << "added vert " << end << " = " << g.verts[end] << "\n";

      for (int j=0; j<slices; ++j) {
        size_t v0 = end - slices + j;
        int offs = (j < slices-1) ? 1 : -(slices-1);
        g.addTriangle(end, v0+offs, v0);
      }

    } else {

      real phi = i * M_PI / stacks;
      real cp = cos(phi);
      real sp = sin(phi);

      for (int j=0; j<slices; ++j) {

        real theta = j * 2.0 * M_PI / slices;
        //std::cerr << "theta = " << theta << "\n";

        real ct = cos(theta);
        real st = sin(theta);
        size_t vi = g.addVertex(r*sp*ct, r*sp*st, r*cp);
        //std::cerr << "added vert " << vi << " = " << g.verts[vi] << "\n";

        int offs = (j < slices-1) ? 1 : -(slices-1);

        if (i == 1) {
          // fan with start
          g.addTriangle(0, vi, vi+offs);
        } else {
          // strip with prev
          g.addTriangle(vi, vi-slices+offs, vi-slices);
          g.addTriangle(vi, vi+offs, vi-slices+offs);
        }

      }

    }

  }


  return g;

}


template <class real>
void TriMesh3_t<real>::transform(const Transform3& xform, TriMesh3_t& dst) const {
  dst.clear();
  dst.addMesh(*this, xform);
}

template <class real>
void TriMesh3_t<real>::transform(const mat4& xform, TriMesh3_t& dst) const {
  dst.clear();
  dst.addMesh(*this, xform);
}

template <class real>
void TriMesh3_t<real>::applyTransform(const Transform3& xform) {

  for (size_t i=0; i<verts.size(); ++i) {
    verts[i] = xform.transformFwd(verts[i]);
  }

  Transform3 t2(xform.rotation(), vec3(0.0));

  for (size_t j=0; j<normals.size(); ++j) {
    normals[j] = t2 * normals[j];
  }

}

template <class real>
void TriMesh3_t<real>::applyTransform(const mat4& m) {

  for (size_t i=0; i<verts.size(); ++i) {
    verts[i] = m.mult3D(verts[i]);
  }

  mat4 m2(m);
  m2(0,3) = m2(1,3) = m2(2,3) = 0.0f;
  
  for (size_t j=0; j<normals.size(); ++j) {
    normals[j] = m2.mult3D(normals[j]);
    normals[j] = vec3::normalize(normals[j]);
  }

}

template <class real> static void glVertex_t(const real*);
template <class real> static void glNormal_t(const real*);

template <> inline void glVertex_t(const double* v) {
  glVertex3dv(v);
}

template <> inline void glVertex_t(const float* v) {
  glVertex3fv(v);
}

template <> inline void glNormal_t(const double* v) {
  glNormal3dv(v);
}

template <> inline void glNormal_t(const float* v) {
  glNormal3fv(v);
}

  
template <class real>
void TriMesh3_t<real>::renderGL(const SplitInfo* split, size_t which) const {

  glBegin(GL_TRIANGLES);

  for (size_t i=0; i<faces.size(); ++i) {

    if (split && split->fmask[i] != which) { continue; }

    const Tri& f = faces[i];

    if (f.nidx[0] == npos) {

      const vec3& v0 = verts[f.vidx[0]];
      const vec3& v1 = verts[f.vidx[1]];
      const vec3& v2 = verts[f.vidx[2]];

      vec3 normal = vec3::cross((v1-v0), (v2 - v0));
      normal = vec3::normalize(normal);

      glNormal_t(normal.v);

      /*
      std::cerr << "normal is " << normal << "\n";
      std::cerr << "vidx[0] is " << f.vidx[0] << "\n";
      std::cerr << "vidx[1] is " << f.vidx[1] << "\n";
      std::cerr << "vidx[2] is " << f.vidx[2] << "\n";
      std::cerr << "v0 is " << v0 << "\n";
      std::cerr << "v1 is " << v1 << "\n";
      std::cerr << "v2 is " << v2 << "\n";
      */

      glVertex_t(v0.v);
      glVertex_t(v1.v);
      glVertex_t(v2.v);

    } else {

      glNormal_t(normals[f.nidx[0]].v);
      glVertex_t(verts[f.vidx[0]].v); 

      glNormal_t(normals[f.nidx[1]].v);
      glVertex_t(verts[f.vidx[1]].v);

      glNormal_t(normals[f.nidx[2]].v); 
      glVertex_t(verts[f.vidx[2]].v);

   }

  }

  glEnd();

}

template <class real>
void TriMesh3_t<real>::renderGL(const Transform3& xform, const SplitInfo* split, size_t which) const {
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glstuff::mult_transform(xform);
  renderGL(split, which);
  glPopMatrix();
}


template <class real>
Box3_t<real> TriMesh3_t<real>::computeBBox() const {
  Box3 bbox;
  for (size_t i=0; i<verts.size(); ++i) {
    bbox.addPoint(verts[i]);
  }
  return bbox;
}

template <class real>
void TriMesh3_t<real>::splitXY(size_t max_tris,
		   size_t max_splits,
		   SplitInfo& info,
		   Box3 bbox) const {

  info.fmask.clear();
  info.ngroups = 0;

  if (empty()) { return; }

  if (bbox.empty()) {
    bbox = computeBBox();
  }

  info.fmask.resize(faces.size(), 0);
  
  std::vector<size_t> counts;
  counts.push_back(faces.size());

  std::vector<vec3> centroids;

  for (size_t i=0; i<faces.size(); ++i) {
    const Tri& f = faces[i];
    const vec3& v0 = verts[f.vidx[0]];
    const vec3& v1 = verts[f.vidx[1]];
    const vec3& v2 = verts[f.vidx[2]];
    centroids.push_back( (v0 + v1 + v2)/real(3) );
  }

  _split(bbox, centroids, 0, 0, 0, max_tris, max_splits, info.fmask, counts);

  //for (size_t which=0; which<counts.size(); ++which) {
    //std::cerr << "counts[" << which << "] = " << counts[which] << "\n";
  //}

  info.ngroups = counts.size();

}

template <class real>
void TriMesh3_t<real>::_split(Box3 bbox,
		  const std::vector<vec3>& centroids,
		  int axis,
		  size_t depth,
		  size_t which,
		  size_t max_tris,
		  size_t max_splits,
		  std::vector<size_t>& fmask,
		  std::vector<size_t>& counts) const {
  
  // check base cases
  if (depth >= max_splits) {
    return;
  } else if (counts[which] <= max_tris) {
    return;
  }

  // now do split
  Box3 left = bbox;
  Box3 right = bbox;

  right.p0[axis] = 0.5 * (bbox.p0[axis] + bbox.p1[axis]);
  left.p1[axis] = 0.5 * (bbox.p0[axis] + bbox.p1[axis]);

  size_t nw = counts.size();
  counts.push_back(0);

  // place every triangle whose centroid falls in the new bbox in the new group!
  for (size_t i=0; i<faces.size(); ++i) {
    if (fmask[i] != which) { continue; }
    if (right.contains(centroids[i])) {
      --counts[which];
      ++counts[nw];
      fmask[i] = nw;
    }
  }

  // call split twice again
  _split(left, centroids, 1-axis, depth+1, which, max_tris, max_splits, fmask, counts);
  _split(right, centroids, 1-axis, depth+1, nw, max_tris, max_splits, fmask, counts);

}

template <class real>
GLuint TriMesh3_t<real>::compileDisplayList() const {

  GLuint list = glGenLists(1);
  glNewList(list, GL_COMPILE);
  renderGL();
  glEndList();
  return list;
  
}

template <class real>
void TriMesh3_t<real>::compileDisplayLists(const SplitInfo& split,
				       std::vector<GLuint>& lists) const {
  

  /*
  const GLubyte rgb[24] = { 
    128, 0,   255, // purple
    0,   0,   255, // blue
    0,   128, 255, // bc
    0,   255, 255, // cyan
    0,   255, 0,   // green
    255, 255, 0,   // yellow
    255, 128, 0,   // orange
    255, 0,   0,   // red
  };
  */

  for (size_t i=0; i<split.ngroups; ++i) {
    lists.push_back(glGenLists(1));
    glNewList(lists.back(), GL_COMPILE);
    //glColor3ubv(rgb + 3*(i%8));
    renderGL(&split, i);
    glEndList();
  }

}


template <class real>
void TriMesh3_t<real>::renderPoints() const {
  glBegin(GL_POINTS);
  for (size_t i=0; i<verts.size(); ++i) {
    glVertex_t(verts[i].v);
  }
  glEnd();
}



// points are defined CCW on xy plane

template <class real>
TriMesh3_t<real> TriMesh3_t<real>::extrude(const Vec2Array& points, 
				   real dist) {

  size_t n_ext = points.size();

  vec2 centroid(0,0);
  for (size_t i=0; i<n_ext; ++i) {
    centroid = centroid + points[i];
  }
  centroid = centroid / points.size();

  return extrude(points, centroid, dist);

}

template <class real>
TriMesh3_t<real> TriMesh3_t<real>::extrude(Vec2Array points, 
                                   const vec2& centroid,
				   real dist) {

  size_t n_ext = points.size();

  points.push_back(centroid);

  IndexArray indices;
  IndexArray exterior;
  for (size_t i=0; i<n_ext; ++i) {
    indices.push_back(i);
    indices.push_back((i+1) % n_ext);
    indices.push_back(n_ext);
    exterior.push_back(i);
  }

  return extrude(points, indices, exterior, dist);
  
}

template <class real>
TriMesh3_t<real> TriMesh3_t<real>::extrude(const Vec2Array& points, 
				   const IndexArray& indices,
				   const IndexArray& exterior,
				   real dist) {

  dist = fabs(dist);

  TriMesh3_t<real> rval;

  for (int i=0; i<2; ++i) {

    for (size_t j=0; j<points.size(); ++j) {
      const vec2& p = points[j];
      rval.addVertex(p[0], p[1], i==0 ? -0.5*dist : 0.5*dist);
    }

    bool flip = (i==0);

    for (size_t k=0; k<indices.size(); k+=3) {
      rval.addTriangle(indices[k]+i*points.size(), 
		       indices[k+1]+i*points.size(), 
		       indices[k+2]+i*points.size(), 
		       flip);
    }

  }

  size_t n_ext = exterior.size();

  for (size_t e=0; e<n_ext; ++e) {
    size_t i00 = e; 
    size_t i01 = (e+1) % n_ext;
    size_t i10 = i00 + points.size();
    size_t i11 = i01 + points.size();
    rval.addTriangle(i00, i01, i10);
    rval.addTriangle(i01, i11, i10);
  }
  
  return rval;

}


// try to parse the specified thing based on the extension
template <class real>
void TriMesh3_t<real>::parse(const std::string& filename) {

  size_t pos = filename.rfind('.');
  if (pos == std::string::npos) {
    std::cerr << "no extension in TriMesh3_t::parse()\n";
    exit(1);
  }

  size_t len = filename.length() - pos - 1;

  std::ifstream istr(filename.c_str());
  
  return parse(istr, filename.substr(pos+1, len));

}

// try to parse the specified thing based on the extension
template <class real>
void TriMesh3_t<real>::parse(std::istream& istr, const std::string& e) {

  std::string extension = lower(e);

  if (extension == "obj") {
    parseObj(istr);
  } else if (extension == "wrl") {
    parseWrl(istr);
  } else if (extension == "stl") {
    parseStl(istr);
  } else {
    std::cerr << "bad extension '" << extension << "' in TriMesh3_t::parse()\n";
    exit(1);
  }


}


template <class real>
void TriMesh3_t<real>::parseObj(const std::string& filename) {

    /*

  // disabled gz for now
  if (0 && filename.find(".gz") == filename.length() - 3) {
    gzistream istr(filename.c_str());
    parseObj(istr);

    } else */ {
    std::ifstream istr(filename.c_str());
    parseObj(istr);
  }

}

template <class real>
void TriMesh3_t<real>::parseObj(const char* filename) {
  parseObj(std::string(filename));
}


template <class real>
void TriMesh3_t<real>::parseObj(std::istream& istr) {

  if (!istr.good()) {
    throw std::runtime_error("bad stream in parseObj");
  }

  clear();

  float v[3];
  char tag;
  char buf[1024];

  while (istr) {
      
    tag = istr.get();
    if (tag == EOF) { break; }
      
    switch (tag) {
    case 'v': {

      buf[0] = tag;
      tag = istr.get();
      bool normal;

      if (tag == ' ') {
        normal = false;
      } else if (tag == 'n') {
        tag = istr.get();
        if (tag != ' ') {
          throw std::runtime_error("error parsing normal!");
        }
        normal = true;
      } else {
        buf[1] = tag;
        istr.getline(buf+2, 1022);
        //if (strlen(buf) > 2) { std::cerr << "OBJ parser: skipping " << buf << "\n"; }
        continue;
        //throw std::runtime_error("error parsing vertex data!");
      }

      istr.getline(buf, 1024);
      if (sscanf(buf, "%f %f %f", &(v[0]), &(v[1]), &(v[2])) != 3) {
        throw std::runtime_error("error parsing vertex/normal data!");
      }

      if (normal) {
        addNormal(v[0], v[1], v[2]);
      } else {
        addVertex(v[0], v[1], v[2]);
      }

      break;
    }
    case 'f': {

      tag = istr.get();
      if (tag != ' ') {
        throw std::runtime_error("error parsing face!");
      }

      istr.getline(buf, 1024);
	
      int i, j, k, a, b, c;
      int s, t, u;
      bool haveNormals;
	
      if (sscanf(buf, "%d//%d %d//%d %d//%d", 
                 &i, &a, &j, &b, &k, &c) == 6) {
        haveNormals = true;
      } else if (sscanf(buf, "%d/%d %d/%d %d/%d", 
                        &i, &s, &j, &t, &k, &u) == 6) {
        haveNormals = false;
      } else if (sscanf(buf, "%d/%d/%d %d/%d/%d %d/%d/%d", 
                        &i, &s, &a, &j, &t, &b, &k, &u, &c) == 9) {
        haveNormals = true;
      } else if (sscanf(buf, "%d %d %d", &i, &j, &k) == 3) {
        haveNormals = false;
      } else {
        throw std::runtime_error("error parsing shape!");
      }
	
      if (haveNormals) {
        addTriangle(i-1, j-1, k-1, a-1, b-1, c-1);
      } else {
        addTriangle(i-1, j-1, k-1);
      }
	
      break;
    }
    default: 
      buf[0] = tag;
      if (tag != '\n') {
        istr.getline(buf+1, 1023);
        //if (strlen(buf) > 2) { std::cerr << "OBJ parser: skipping " << buf << "\n"; }
      }
      break;
    }

  }


}



template <class real>
void TriMesh3_t<real>::saveObj(const std::string& filename) const {
  saveObj(filename.c_str());
}

template <class real>
void TriMesh3_t<real>::saveObj(const char* filename) const {

  FILE* outfile = fopen(filename, "w");
  if (!outfile) { 
    fprintf(stderr, "martin says this is an error!\n");
    return; 
  }

  for (size_t i=0; i<verts.size(); ++i) {
    const vec3& v = verts[i];
    fprintf(outfile, "v %f %f %f\n", v.x(), v.y(), v.z());
  }
  fprintf(outfile, "\n");

  for (size_t i=0; i<normals.size(); ++i) {
    const vec3& n = normals[i];
    fprintf(outfile, "vn %f %f %f\n", n.x(), n.y(), n.z());
  }
  fprintf(outfile, "\n");

  for (size_t i=0; i<faces.size(); ++i) {
    const Tri& t = faces[i];
    if (t.nidx[0] != npos &&
        t.nidx[0] != npos &&
        t.nidx[0] != npos) {
      fprintf(outfile, "f %u//%u %u//%u %u//%u\n", 
              (unsigned int)(t.vidx[0]+1), 
              (unsigned int)(t.nidx[0]+1), 
              (unsigned int)(t.vidx[1]+1), 
              (unsigned int)(t.nidx[1]+1), 
              (unsigned int)(t.vidx[2]+1),
              (unsigned int)(t.nidx[2]+1));
    } else {
      fprintf(outfile, "f %u %u %u\n", 
              (unsigned int)(t.vidx[0]+1), 
              (unsigned int)(t.vidx[1]+1), 
              (unsigned int)(t.vidx[2]+1));
    }
  }

  fclose(outfile);

}

template <class real>
void TriMesh3_t<real>::compact() {

  IndexArray newvidx(verts.size(), npos);
  IndexArray newnidx(normals.size(), npos);

  Vec3Array newverts;
  Vec3Array newnormals;

  for (size_t i=0; i<faces.size(); ++i) {
    Tri& t = faces[i];
    for (int j=0; j<3; ++j) {
      size_t oldvidx = t.vidx[j];
      if (newvidx[oldvidx] == npos) {
        newvidx[oldvidx] = newverts.size();
        newverts.push_back(verts[oldvidx]);
      }
      t.vidx[j] = newvidx[oldvidx];
      size_t oldnidx = t.nidx[j];
      if (oldnidx != npos) {
        if (newnidx[oldnidx] == npos) {
          newnidx[oldnidx] = newnormals.size();
          newnormals.push_back(normals[oldnidx]);
        }
        t.nidx[j] = newnidx[oldnidx];
      }
    }
  }

  verts.swap(newverts);
  normals.swap(newnormals);
  
}

template <class real>
TriMesh3_t<real> TriMesh3_t<real>::lathe(const Vec2Array& points, int slices) {

  TriMesh3_t<real> g;

  if (points.size() < 2) { return g; }

  const bool closed = (points.front()[0] == 0 && points.back()[0] == 0);


  size_t vs = 0;
  size_t js = 0;
  size_t np = points.size();

  if (closed) {
    g.addVertex(vec3(points.front()[0], points.front()[1], 0));
    g.addVertex(vec3(points.back()[0], points.back()[1], 0));
    vs = 2;
    js = 1;
    np -= 2;
  }

  //std::cerr << "\nclosed = " << closed << "\n";
  //std::cerr << "np = " << np << "\n";

  for (int i=0; i<slices; ++i) {

    int ii = (i + 1) % slices;

    real theta = 2*M_PI * real(i) / real(slices);
    Transform3 t = Transform3::ry(theta, vec3(0.0));

    if (closed) {
      // top triangle
      g.addTriangle(0, vs+np*ii, vs+np*i);
    }

    for (size_t j=0; j<np; ++j) {

      size_t jj = (j + 1) % np;

      const vec2& p = points[j+js]; 
      g.addVertex(t * vec3(p[0], p[1], 0));
     
      size_t v00 = vs +  i * np +  j;
      size_t v01 = vs +  i * np + jj;
      size_t v10 = vs + ii * np +  j;
      size_t v11 = vs + ii * np + jj;

      if (closed && j == np-1) { continue; }

      g.addTriangle(v00, v10, v11);
      g.addTriangle(v00, v11, v01);
      
    }

    if (closed) {
      // bottom triangle
      g.addTriangle(1, vs+np*i+np-1, vs+np*ii+np-1);
    }

  }

  /*
  // debug:
  for (size_t i=0; i<g.verts.size(); ++i) {
    std::cerr << "g.verts[" << i << "] = " << g.verts[i] << "\n";
  }
  for (size_t i=0; i<g.faces.size(); ++i) {
    const TriMesh3_t<real>::Tri& t = g.faces[i];
    std::cerr << "g.faces[" << i << "] = " 
              << t.vidx[0] << ","
              << t.vidx[1] << ","
              << t.vidx[2] << "\n";
  }
  */

  return g;
  
}

#define debug if(0) std::cerr

template <class real>
class VRMLParser {
public:

  TriMesh3_t<real>& mesh;
  std::istream& istr;

  VRMLParser(TriMesh3_t<real>& m, std::istream& i): 
    mesh(m), istr(i) 

  {

    parseHeader();
    
    while (parseNode());

  }
  
  void parseHeader() {
    std::string l;
    std::getline(istr, l);
    if (l != "#VRML V2.0 utf8") { 
      throw std::runtime_error("Bad vrml header: " + l);
    }
    debug << "parsed header\n";
  }


  void consumeWhitespace() {
    while (isspace(istr.peek())) { istr.get(); }
  }

  std::string parseIdentifier() {

    char buf[1024];
    buf[0] = 0;

    char* cur = buf;
    const char* end = buf+1023;

    consumeWhitespace();

    int c = istr.peek();
    if (c == EOF) { return ""; }
    if (!isalpha(c)) { 
      return "";
    }

    while (cur < end) {
      int c = istr.peek();
      if (isalnum(c) || c == '_') { 
        *cur++ = istr.get(); *cur = 0; 
      } else { 
        break; 
      }
    }

    if (cur == end) { 
      throw std::runtime_error("identifier too long");
    }

    debug << "parsed identifier " << buf << "\n";
    return buf;
        
  }

  void matchIdentifier(const std::string& id) {
    if (parseIdentifier() != id) {
      throw std::runtime_error("Expected " + id);
    }
  }

  void parseToken(char c) {
    consumeWhitespace();
    char p = istr.peek();
    if (p != c) { 
      throw std::runtime_error(std::string("Expected ") + c + " but got " + p);
    }
    istr.get();
  }

  void parseTransform() {
    parseToken('{');
    matchIdentifier("children");
    parseToken('[');
    while (parseNode());
    parseToken(']');
    parseToken('}');
  }
  
  void ignore(char open, char close) {

    parseToken(open);
    int cnt = 1;

    while (cnt) {
      char c = istr.get();
      if (!istr) { 
        throw std::runtime_error(std::string("EOF waiting for closing ") + close);
      }
      if (c == open) { ++cnt; } 
      else if (c == close) { --cnt; }
    }

  }
  
  void parseShape() {
    parseToken('{');
    while (true) {
      std::string s = parseIdentifier();
      if (s.empty()) { break; }
      if (s == "appearance" || s == "geometry") { 
        parseNode();
      } 
    }
    parseToken('}');
  }
  
  void parseIndexedFaceSet() {

    parseToken('{');

    size_t offset = mesh.verts.size();

    while (true) {

      std::string s = parseIdentifier();
      if (s.empty()) { break; }
      
      if (s == "coord") {

        matchIdentifier("Coordinate");
        parseToken('{');
        matchIdentifier("point");
        parseToken('[');
        

        vec3_t<real> v;

        while (true) {
          consumeWhitespace();
          char c = istr.peek();
          if (c != '-' && !isdigit(c) && c != '.') { break; }
          if (!(istr >> v[0] >> v[1] >> v[2])) { 
            throw std::runtime_error("Error parsing coordinates");
          }
          parseToken(',');
          //debug << "parsed vertex " << v << "\n";
          mesh.verts.push_back(v);
        }

        

        parseToken(']');
        parseToken('}');

        debug << "done parsing coords\n";
        
      } else if (s == "coordIndex") {

        parseToken('[');

        long int indices[4];

        while (true) {

          consumeWhitespace();
          char c = istr.peek();
          if (!isdigit(c)) { break; }

          int i;
          for (i=0; i<4; ++i) {
            if (!(istr >> indices[i])) { break; }
            parseToken(',');
            if (i != 3) { indices[i] += offset; }
          }

          if (i != 4) { 
            throw std::runtime_error("Error parsing indices");
          }

          if (indices[3] != -1) {
            throw std::runtime_error("Non-triangular geometry");
          }

          mesh.addTriangle(indices[0], indices[1], indices[2]);
          
        }


        parseToken(']');

      }

    }


    parseToken('}');
    debug << "done parsing IndexedFaceSet\n";

  }


  bool parseNode() {
    
    std::string s = parseIdentifier();
    if (s.empty()) { return false; }

    if (s == "DEF") { 
      s = parseIdentifier();
      s = parseIdentifier();
    }

    if (s == "Transform") { 
      parseTransform();
    } else if (s == "Shape") {
      parseShape();
    } else if (s == "IndexedFaceSet") {
      parseIndexedFaceSet();
    } else if (s == "Coordinate") {
      parseIndexedFaceSet();
    } else if (s == "Appearance") {
      ignore('{', '}');
    } else {
      throw std::runtime_error("Invalid node type: " + s);
    }

    return true;

  }


};

template <class real>
void TriMesh3_t<real>::parseWrl(const char* filename) {
  parseWrl(std::string(filename));
}

template <class real>
void TriMesh3_t<real>::parseWrl(const std::string& filename) {
  std::ifstream istr(filename.c_str());
  parseWrl(istr);
}


template <class real>
void TriMesh3_t<real>::parseWrl(std::istream& istr) {
  
  try { 
    VRMLParser<real> parser(*this, istr);
  } catch (std::exception& e) {
    std::cerr << e.what() << "\n";
    throw;
  }

}

template <class real>
void TriMesh3_t<real>::parseStl(const char* filename) {
  parseStl(std::string(filename));
}

template <class real>
void TriMesh3_t<real>::parseStl(const std::string& filename) {
  std::ifstream istr(filename.c_str());
  parseStl(istr);
}

template <class real>
static inline bool operator<(const vec3_t<real>& a,
                             const vec3_t<real>& b) {

  return ( (a.x() < b.x()) ||
           (a.x() == b.x() && a.y() < b.y()) ||
           (a.x() == b.x() && a.y() == b.y() && a.z() < b.z()) );

}

template <class real>
void TriMesh3_t<real>::parseStl(std::istream& istr) {

  typedef std::map< vec3_t<real>, size_t > IndexLookup;

  IndexLookup lookup;
  size_t duplicates = 0;
  

  char header[80];
  const char* magic = "solid ";
  if (!istr.read(header, 6)) {
    std::cerr << "error reading first 6 bytes of STL file!\n";
    exit(1);
  }

  if (memcmp(magic, header, 6) == 0) {

    // it's ASCII
    std::cerr << "error: can't parse ASCII STL files yet!\n";
    exit(1);

  }

  if (!istr.read(header+6, 80-6)) { 
    std::cerr << "error reading binary header!\n";
    exit(1);
  }

  uint32_t count;
  if (!istr.read((char*)&count, sizeof(count))) {
    std::cerr << "error reading count!\n";
    exit(1);
  }

  
  for (uint32_t i=0; i<count; ++i) {
    
    vec3 v[4];
    for (int j=0; j<4; ++j) {
      for (int a=0; a<3; ++a) {
        float c;
        if (!istr.read((char*)&c, sizeof(c))) {
          std::cerr << "error reading vertex!\n";
          exit(1);
        }
        v[j][a] = c;
      }
    }

    size_t idx[3];
    for (int p=0; p<3; ++p) {
      const vec3_t<real>& vp = v[p+1];
      typename IndexLookup::const_iterator it = lookup.find(vp);
      if (it == lookup.end()) {
        lookup[vp] = idx[p] = addVertex(vp);
      } else {
        assert( it->first == vp );
        ++duplicates;
        idx[p] = it->second;
      }
    }

    addTriangle(idx[0], idx[1], idx[2]);

    uint16_t ignore;
    if (!istr.read((char*)&ignore, sizeof(ignore))) {
      std::cerr << "error reading ignore!\n";
      exit(1);
    }

  }

}

template class TriMesh3_t<float>;
template class TriMesh3_t<double>;
