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

#ifndef _MAT4_H_
#define _MAT4_H_

#include "vec4.h"

/** Encapsulate a 4-by-4 matrix used for transforming homogenous 3D
 *  coordinates. */

template <class real> class mat4_t {
public:

  typedef real value_type;

  /** Default constructor initializes to identity matrix. */
  mat4_t() { *this = identity(); }

  /** Construct with optional initialization.
   *
   *  @param init If this is true, then initialize to the identity
   *              matrix, otherwise, leave data uninitialized.
   */
  mat4_t(bool init) { if (init) { *this = identity(); } }

  /** Copy from row-major 4-by-4 matrix. */
  mat4_t(const real src[16]) { memcpy(data, src, 16*sizeof(real)); }

  /** Copy constructor. */
  mat4_t(const mat4_t& src) { memcpy(data, src.data, 16*sizeof(real)); }

  /** Assignment operator. */
  mat4_t& operator=(const mat4_t& src) { 
    memcpy(data, src.data, 16*sizeof(real));
    return *this;
  }

  /** Return the transpose of this matrix. */
  mat4_t transpose() const {
    const unsigned int tidx[16] = {
      0, 4, 8,  12,
      1, 5, 9,  13,
      2, 6, 10, 14,
      3, 7, 11, 15,
    };
    mat4_t rval(false);
    for (unsigned int i=0; i<16; ++i) {
      rval.data[i] = data[tidx[i]];
    }
    return rval;
  };

  /** Return the identity matrix. */
  static mat4_t identity() {
    const real id[16] = {
      1.0, 0.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 0.0,
      0.0, 0.0, 1.0, 0.0,
      0.0, 0.0, 0.0, 1.0,
    };
    return mat4_t(id);
  };


  vec4_t<real> row(size_t i) const {
    return vec4_t<real>(data[4*i+0], 
                        data[4*i+1], 
                        data[4*i+2],
                        data[4*i+3]);
  }

  vec4_t<real> col(size_t i) const {
    return vec4_t<real>(data[i+ 0],
                        data[i+ 4],
                        data[i+ 8],
                        data[i+12]);
  }

  void setRow(size_t i, const vec4_t<real>& v) {
    data[4*i+0] = v[0];
    data[4*i+1] = v[1];
    data[4*i+2] = v[2];
    data[4*i+3] = v[3];
  }

  void setCol(size_t i, const vec4_t<real>& v) {
    data[i+ 0] = v[0];
    data[i+ 4] = v[1];
    data[i+ 8] = v[2];
    data[i+12] = v[3];
  }

  /** Elementwise access. Data is stored in a row-major format. */
  const real& operator[](unsigned int i) const {
    return data[i];
  }

  /** Subscript access. */
  const real& operator()(unsigned int row, unsigned int col) const {
    return data[row*4 + col];
  }

  /** Elementwise access. Data is stored in a row-major format. */
  real& operator[](unsigned int i)  {
    return data[i];
  }
  
  /** Subscript access. */
  real& operator()(unsigned int row, unsigned int col)  {
    return data[row*4 + col];
  }

  /** Return the product of a matrix and a 3D vector, ignoring the
   *  last row of the matrix.
   *
   *  The w (homogenous coordinate) of the input vector is assumed to
   *  be 1, and the resulting vector is normalized to have a w
   *  coordinate of 1.
   */
  vec3_t<real> mult3D(const vec3_t<real>& v) const {

    const mat4_t& m = *this;
    
    vec3_t<real> rval(v[0] * m(0,0) + v[1] * m(0,1) + v[2] * m(0,2) + m(0,3),
                      v[0] * m(1,0) + v[1] * m(1,1) + v[2] * m(1,2) + m(1,3),
                      v[0] * m(2,0) + v[1] * m(2,1) + v[2] * m(2,2) + m(2,3));
    
    return rval;
    
  }

  /** TODO: document. */
  vec3_t<real> inv3D(const vec3_t<real>& v) const {

    const mat4_t& m = *this;
    
    real a = v[0] - m(0,3);
    real b = v[1] - m(1,3);
    real c = v[2] - m(2,3);
    
    vec3_t<real> rval(a * m(0,0) + b * m(1,0) + c * m(2,0),
                      a * m(0,1) + b * m(1,1) + c * m(2,1),
                      a * m(0,2) + b * m(1,2) + c * m(2,2));
    
    return rval;

  }

  /** Return the product of a matrix and a 3D vector.
   *
   *  The w (homogenous coordinate) of the input vector is assumed to
   *  be 1, and the resulting vector is normalized to have a w
   *  coordinate of 1.
   */
  vec3_t<real> proj3D(const vec3_t<real>& v) const {
    
    const mat4_t& m = *this;

    vec3_t<real> rval(v[0] * m(0,0) + v[1] * m(0,1) + v[2] * m(0,2) + m(0,3),
			 v[0] * m(1,0) + v[1] * m(1,1) + v[2] * m(1,2) + m(1,3),
			 v[0] * m(2,0) + v[1] * m(2,1) + v[2] * m(2,2) + m(2,3));
    
    real norm = 
      v[0] * m(3,0) + v[1] * m(3,1) + v[2] * m(3,2) + m(3,3);
    
    return rval / norm;
    
  }


  /** Outer product matrix. */
  static mat4_t outer(const vec4_t<real>& a, 
                      const vec4_t<real>& b) {
    mat4_t rval(false);
    for (int row=0; row<4; ++row) {
      for (int col=0; col<4; ++col) {
        rval(row,col) = a[row]*b[col];
      }
    }
    return rval;
  }



  /** Underlying data. */
  real data[16];

};


/** Return the product a matrix and a scalar. */
template <class real, class real2> inline mat4_t<real> operator*(const mat4_t<real>& A, real2 s) {
  mat4_t<real> rval(false);
  for (int i=0; i<16; ++i) { rval.data[i] = A.data[i]*s; }
  return rval;
}

/** Return the product a matrix and a scalar. */
template <class real, class real2> inline mat4_t<real> operator*(real2 s, const mat4_t<real>& A) {
  return A*s;
}

/** Return the product of two matrices. */
template <class real> inline mat4_t<real> operator*(const mat4_t<real>& A, const mat4_t<real>& B) {
  mat4_t<real> rval(false);
  for (unsigned int i=0; i<4; ++i) {
    for (unsigned int j=0; j<4; ++j) {
      real accum = 0.0;
      // out(i,j) = ith row of A dotted with jth col of B
      for (unsigned int k=0; k<4; ++k) {
	accum += A(i, k) * B(k, j);
      }
      rval(i, j) = accum;
    }
  }
  return rval;
}

/** Return the product of a matrix and a 4D vector. */
template <class real> inline vec4_t<real> operator*(const mat4_t<real>& m, const vec4_t<real>& v) {
  
  vec4_t<real> rval(v[0] * m(0,0) + v[1] * m(0,1) + v[2] * m(0,2) + v[3] * m(0,3),
		       v[0] * m(1,0) + v[1] * m(1,1) + v[2] * m(1,2) + v[3] * m(1,3),
		       v[0] * m(2,0) + v[1] * m(2,1) + v[2] * m(2,2) + v[3] * m(2,3),
		       v[0] * m(3,0) + v[1] * m(3,1) + v[2] * m(3,2) + v[3] * m(3,3));
  
  return rval;

}

/** Matrix stream output. */
template <class real> inline std::ostream& operator<<(std::ostream& ostr, const mat4_t<real>& m) {
  for (unsigned int i=0; i<4; ++i) {
    for (unsigned int j=0; j<4; ++j) {
      ostr << std::setw(10) << m(i,j) << " ";
    }
    ostr << "\n";
  }
  return ostr;
}

typedef mat4_t<double> mat4d;
typedef mat4_t<float> mat4f;

#endif
