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

#ifndef _VEC3_H_
#define _VEC3_H_

#include "vec2.h"

/** Encapsulate a 3D vector. */
template <class real> class vec3_t {
public:

  typedef real value_type;
  enum { size = 3 };

  real v[3];

  const real& x() const { return v[0]; }
  real& x() { return v[0]; }

  const real& y() const { return v[1]; }
  real& y() { return v[1]; }

  const real& z() const { return v[2]; }
  real& z() { return v[2]; }

  /** Construct the null vector. */
  vec3_t() {
    v[0] = v[1] = v[2] = 0.0;
  }

  /** Construct such that x, y, and z, all take the given value. */
  explicit vec3_t(real d) {
    v[0] = v[1] = v[2] = d; 
  }

  /** Construct from vec2 and scalar. */
  vec3_t(const vec2_t<real>& xy, real z) {
    v[0] = xy[0];
    v[1] = xy[1];
    v[2] = z;
  }

  /** Construct from given coordinates. */
  vec3_t(real x, real y, real z) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }

  /** Construct from other Tval. */
  template <class Tval2>
  explicit vec3_t(const vec3_t<Tval2>& vv) {
    v[0] = vv.x();
    v[1] = vv.y();
    v[2] = vv.z();
  }
  
  /** Element-wise access. */
  const real& operator[](unsigned int i) const { return v[i]; }

  /** Element-wise access. */
  real& operator[](unsigned int i) { return v[i]; }

  /** Return the magnitude of this vector. */
  real norm() const {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  }

  /** Return the squared magnitude of this vector. */
  real norm2() const {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  }

  /** Return the element with the largest magnitude */
  real infnorm() const {
    return std::max(std::max(fabs(v[0]), fabs(v[1])), fabs(v[2]));
  }

  /** Return a 3D vector equivelent to scaling the first two coordinates by the third. */
  vec2_t<real> proj() const {
    return vec2_t<real>(v[0]/v[2], v[1]/v[2]);
  }

  /** Return a 3D vector containing only the first two components of this vector. */
  vec2_t<real> trunc() const {
    return vec2_t<real>(v[0], v[1]);
  }

  /** Return a copy of the input vector with the same direction, but unit magnitude. */
  static vec3_t normalize(const vec3_t& v) {
    real l = v.norm();
    return vec3_t(v[0]/l, v[1]/l, v[2]/l);
  };

  /** Vector dot product. */
  static real dot(const vec3_t& v1, const vec3_t& v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  }

  /** Vector cross product. */
  static vec3_t cross(const vec3_t& u, const vec3_t& v) {
    return vec3_t(u[1]*v[2] - u[2]*v[1],
                  u[2]*v[0] - u[0]*v[2],
                  u[0]*v[1] - u[1]*v[0]);
  }
  
  /** Linear interpolation */
  template <class real2>
  static vec3_t lerp(const vec3_t& v0, const vec3_t& v1, real2 u);


  /** Cubic interpolation */
  template <class real2>
  static vec3_t smoothstep(const vec3_t& v0, const vec3_t& v1, real2 u) {
    real2 uu = 2*u*u*u + 3*u*u;
    return lerp(v0, v1, uu);
  }

  vec3_t& operator+=(const vec3_t& v2) {
    v[0] += v2[0];
    v[1] += v2[1];
    v[2] += v2[2];
    return *this;
  }

  vec3_t& operator-=(const vec3_t& v2) {
    v[0] -= v2[0];
    v[1] -= v2[1];
    v[2] -= v2[2];
    return *this;
  }

  template <class real2>
  vec3_t& operator*=(real2 s) {
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
    return *this;
  }

  template <class real2>
  vec3_t& operator/=(real2 s) {
    v[0] /= s;
    v[1] /= s;
    v[2] /= s;
    return *this;
  }

};

/** Vector addition. */
template <class real> inline vec3_t<real> operator+(const vec3_t<real>& u, const vec3_t<real>& v) {
  return vec3_t<real>(u[0] + v[0], u[1] + v[1], u[2] + v[2]);
}

/** Vector negation. */
template <class real> inline vec3_t<real> operator-(const vec3_t<real>& v) {
  return vec3_t<real>(-v[0], -v[1], -v[2]);
}

/** Vector subtraction. */
template <class real> inline vec3_t<real> operator-(const vec3_t<real>& u, const vec3_t<real>& v) {
  return vec3_t<real>(u[0] - v[0], u[1] - v[1], u[2] - v[2]);
}

/** Vector multiplication by scalar. */
template <class real, class real2> inline vec3_t<real> operator*(const vec3_t<real>& v, real2 s) {
  return vec3_t<real>(v[0]*s, v[1]*s, v[2]*s);
}

/** Vector multiplication by scalar. */
template <class real, class real2> inline vec3_t<real> operator*(real2 s, const vec3_t<real>& v) {
  return vec3_t<real>(v[0]*s, v[1]*s, v[2]*s);
}

/** Vector division. */
template <class real, class real2> inline vec3_t<real> operator/(const vec3_t<real>& v, real2 s) {
  return vec3_t<real>(v[0]/s, v[1]/s, v[2]/s);
}

/** Vector stream output. */
template <class real> inline std::ostream& operator<<(std::ostream& ostr, const vec3_t<real>& v) {
  ostr << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
  return ostr;
}


template <class real> template<class real2> inline vec3_t<real> vec3_t<real>::lerp(const vec3_t<real>& v0, const vec3_t<real>& v1, real2 u) {
  if (u < 0) { return v0; }
  if (u > 1) { return v1; }
  return (1-u)*v0 + u*v1;
}

template <class real> inline bool operator==(const vec3_t<real>& v1, const vec3_t<real>& v2) {
  return v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2];
}

template <class real> inline bool operator!=(const vec3_t<real>& v1, const vec3_t<real>& v2) {
  return v1[0] != v2[0] || v1[1] != v2[1] || v1[2] != v2[2];
}

typedef vec3_t<double> vec3d;
typedef vec3_t<float> vec3f;

#endif

