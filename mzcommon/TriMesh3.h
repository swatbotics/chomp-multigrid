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

#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include "glstuff.h"
#include "Box3.h"
#include "vec2.h"

class TriMesh3Base {
public:

  typedef std::vector<size_t> IndexArray;

  struct Tri {

    size_t vidx[3];
    size_t nidx[3];
    size_t fdata;
    
    Tri(const size_t v[3], const size_t n[3], size_t fdata=0);

    Tri(size_t i, size_t j, size_t k,
	size_t a, size_t b, size_t c,
        size_t fdata=0);

  };

  struct SplitInfo {
    IndexArray fmask;
    size_t ngroups;
  };

  typedef std::vector<Tri>  TriArray;

};

template <class real> class TriMesh3_t: public TriMesh3Base {
public:

  typedef vec2_t<real> vec2;
  typedef vec3_t<real> vec3;
  typedef mat4_t<real> mat4;
  typedef Transform3_t<real> Transform3;
  typedef Box3_t<real> Box3;

  typedef std::vector<vec3> Vec3Array;
  typedef std::vector<vec2> Vec2Array;

  Vec3Array verts;
  Vec3Array normals;
  TriArray  faces;

  vec3 centroid(size_t faceidx) const;

  TriMesh3_t();

  void swap(TriMesh3_t<real>& g);
  void clear();
  bool empty() const;

  size_t addVertex(real x, real y, real z);
  size_t addNormal(real x, real y, real z);

  size_t addVertex(const vec3& v);
  size_t addNormal(const vec3& v);

  size_t addTriangle(size_t i, 
                     size_t j, 
                     size_t k, 
                     bool flip=false,
                     size_t fdata=0);
  
  size_t addTriangle(size_t i, 
                     size_t j, 
                     size_t k,
                     size_t a,
                     size_t b,
                     size_t c,
                     bool flip=false,
                     size_t fdata=0);

  void addMesh(const TriMesh3_t<real>& mesh);

  void addMesh(const TriMesh3_t<real>& mesh,
                   const Transform3& xform);

  void addMesh(const TriMesh3_t<real>& mesh,
                   const mat4& xform);

  // Load an obj file into this
  void parseObj(const char* filename);

  void parseObj(const std::string& filename);
  
  void parseObj(std::istream& istr);

  // Load an wrl file into this
  void parseWrl(const char* filename);

  void parseWrl(const std::string& filename);
  
  void parseWrl(std::istream& istr);

  // Load an stl file into this
  void parseStl(const char* filename);

  void parseStl(const std::string& filename);
  
  void parseStl(std::istream& istr);

  // try to parse the specified thing based on the extension
  void parse(const std::string& filename);

  // try to parse the specified thing based on the extension
  void parse(std::istream& istr, const std::string& extension);

  // get rid of unused vertices and normals
  void compact();

  // build a box
  static TriMesh3_t box(real length,
                        real width,
                        real height);
  
  // length along x,
  // width along y,
  // height along z
  static TriMesh3_t cappedBox(real length,
                              real width,
                              real height,
                              int slices);

  static TriMesh3_t turnBox(real l1, real l2,
                            real width,
                            real height,
                            int slices1,
                            int slices2);

  // cylinder along z
  static TriMesh3_t cylinder(real r,
                             real h,
                             int slices);

  static TriMesh3_t sphere(real r,
                           int slices,
                           int stacks);

  // lathe along y, points should be defined CCW on xy plane
  static TriMesh3_t lathe(const Vec2Array& points,
                          int slices);

  // extrude along z, points should be defined CCW on xy plane
  // this will automatically compute the centroid of the points
  // and make a triangle fan around it
  static TriMesh3_t extrude(const Vec2Array& points, 
                            real dist);

  // extrude along z, points should be defined CCW on xy plane
  // uses your centroid
  static TriMesh3_t extrude(Vec2Array points, 
                            const vec2& centroid,
                            real dist);

  // extrude along z, triangles are given by indices
  // note that indices must have length divisible by 3
  static TriMesh3_t extrude(const Vec2Array& points, 
                            const IndexArray& indices,
                            const IndexArray& exterior,
                            real dist);

  void applyTransform(const Transform3& xform);

  void applyTransform(const mat4& xform);

  void transform(const Transform3& xform,
                 TriMesh3_t& dst) const;

  void transform(const mat4& xform,
                 TriMesh3_t& dst) const;

  Box3 computeBBox() const;

  void renderPoints() const;
  
  void renderGL(const SplitInfo* split=0, 
		size_t which=0) const;

  void renderGL(const Transform3& transform,
		const SplitInfo* split=0,
		size_t which=0) const;

  void splitXY(size_t max_tris,
	       size_t max_splits,
	       SplitInfo& info,
	       Box3 bbox = Box3()) const;

  GLuint compileDisplayList() const;
  
  void compileDisplayLists(const SplitInfo& split,
			   std::vector<GLuint>& lvec) const;

  void saveObj(const char* filename) const;

  void saveObj(const std::string& filename) const;
  
  void saveObj(std::ostream& ostr) const;

private:

  void _split(Box3 bbox,
	      const std::vector<vec3>& centroids,
	      int axis,
	      size_t depth,
	      size_t which,
	      size_t max_tris,
	      size_t max_splits,
	      IndexArray& fmask,
	      IndexArray& counts) const;

};

typedef TriMesh3_t<double> TriMesh3d;
typedef TriMesh3_t<float>  TriMesh3f;

#endif
