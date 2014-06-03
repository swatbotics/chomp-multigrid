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

#ifndef _VEC4_H_
#define _VEC4_H_

#include "vec3.h"

/** Encapsulate a 3D vector. */
template <class real> class vec4_t {
public:

  typedef real value_type;
  enum { size = 4 };

  real v[4]; 

  const real& x() const { return v[0]; }
  real& x() { return v[0]; }

  const real& y() const { return v[1]; }
  real& y() { return v[1]; }

  const real& z() const { return v[2]; }
  real& z() { return v[2]; }

  const real& w() const { return v[3]; }
  real& w() { return v[3]; }

  /** Construct the null vector. */
  vec4_t() {
    v[0] = v[1] = v[2] = v[3] = 0.0;
  }

  /** Construct such that x, y, and z, all take the given value. */
  explicit vec4_t(real d) {
    v[0] = v[1] = v[2] = v[3] = d;
  }

  /** Construct from given coordinates. */
  vec4_t(real x, real y, real z, real w=1.0) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
    v[3] = w;
  }

  /** Construct from vec3 and scalar */
  vec4_t(const vec3_t<real>& v3, real w) {
    v[0] = v3[0];
    v[1] = v3[1];
    v[2] = v3[2];
    v[3] = w;
  }

  template <class other>
  explicit vec4_t(const vec4_t<other>& vv) {
    v[0] = vv[0];
    v[1] = vv[1];
    v[2] = vv[2];
    v[3] = vv[3];
  }

  /** Element-wise access. */
  const real& operator[](unsigned int i) const { return v[i]; }

  /** Element-wise access. */
  real& operator[](unsigned int i) { return v[i]; }

  /** Return the magnitude of this vector. */
  real norm() const {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
  }

  /** Return the squared magnitude of this vector. */
  real norm2() const {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3] * v[3];
  }

  /** Return the element with the largest magnitude */
  real infnorm() const {
    return std::max(std::max(fabs(v[0]), fabs(v[1])), 
                    std::max(fabs(v[2]), fabs(v[2])));
  }

  /** Return a 3D vector equivelent to scaling the first three coordinates by the fourth. */
  vec3_t<real> proj() const {
    return vec3_t<real>(v[0]/v[3], v[1]/v[3], v[2]/v[3]);
  }

  /** Return a 3D vector containing only the first three components of this vector. */
  vec3_t<real> trunc() const {
    return vec3_t<real>(v[0], v[1], v[2]);
  }

  /** Return a copy of the input vector with the same direction, but unit magnitude. */
  static vec4_t normalize(const vec4_t& v) {
    real l = v.norm();
    return vec4_t(v[0]/l, v[1]/l, v[2]/l, v[3]/l);
  };

  /** Vector dot product. */
  static real dot(const vec4_t& v1, const vec4_t& v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3];
  }

  /** Linear interpolation */
  static vec4_t lerp(const vec4_t& v0, const vec4_t& v1, real u) {
    if (u < 0) { return v0; }
    if (u > 1) { return v1; }
    return (1-u)*v0 + u*v1;
  }

  /** Cubic interpolation */
  static vec4_t smoothstep(const vec4_t& v0, const vec4_t& v1, real u) {
    real uu = 2*u*u*u + 3*u*u;
    return lerp(v0, v1, uu);
  }


  vec4_t& operator+=(const vec4_t& v2) {
    v[0] += v2[0];
    v[1] += v2[1];
    v[2] += v2[2];
    v[3] += v2[3];
    return *this;
  }

  vec4_t& operator-=(const vec4_t& v2) {
    v[0] -= v2[0];
    v[1] -= v2[1];
    v[2] -= v2[2];
    v[3] -= v2[3];
    return *this;
  }

  template <class real2>
  vec4_t& operator*=(real2 s) {
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
    v[3] *= s;
    return *this;
  }

  template <class real2>
  vec4_t& operator/=(real2 s) {
    v[0] /= s;
    v[1] /= s;
    v[2] /= s;
    v[3] /= s;
    return *this;
  }

};

/** Vector addition. */
template <class real> inline vec4_t<real> operator+(const vec4_t<real>& u, const vec4_t<real>& v) {
  return vec4_t<real>(u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3]);
}

/** Vector negation. */
template <class real> inline vec4_t<real> operator-(const vec4_t<real>& v) {
  return vec4_t<real>(-v[0], -v[1], -v[2], -v[3]);
}

/** Vector subtraction. */
template <class real> inline vec4_t<real> operator-(const vec4_t<real>& u, const vec4_t<real>& v) {
  return vec4_t<real>(u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]);
}

/** Vector multiplication by scalar. */
template <class real, class real2> inline vec4_t<real> operator*(const vec4_t<real>& v, real2 s) {
  return vec4_t<real>(v[0]*s, v[1]*s, v[2]*s, v[3]*s);
}

/** Vector multiplication by scalar. */
template <class real, class real2> inline vec4_t<real> operator*(real2 s, const vec4_t<real>& v) {
  return vec4_t<real>(v[0]*s, v[1]*s, v[2]*s, v[3]*s);
}

/** Vector division. */
template <class real, class real2> inline vec4_t<real> operator/(const vec4_t<real>& v, real2 s) {
  return vec4_t<real>(v[0]/s, v[1]/s, v[2]/s, v[3]/s);
}

template <class real> inline bool operator==(const vec4_t<real>& v1, const vec4_t<real>& v2) {
  return v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2] && v1[3] == v2[3];
}

template <class real> inline bool operator!=(const vec4_t<real>& v1, const vec4_t<real>& v2) {
  return v1[0] != v2[0] || v1[1] != v2[1] || v1[2] != v2[2] || v1[3] != v2[3];
}

/** Vector stream output. */
template <class real> inline std::ostream& operator<<(std::ostream& ostr, const vec4_t<real>& v) {
  ostr << "(" << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ")";
  return ostr;
}

typedef vec4_t<double> vec4d;
typedef vec4_t<float> vec4f;

#endif
