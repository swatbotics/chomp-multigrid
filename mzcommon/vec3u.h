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

#ifndef _VEC3U_H_
#define _VEC3U_H_

#include <iostream>

struct vec3u {

  size_t v[3];

  vec3u() {}

  explicit vec3u(size_t k) {
    v[0] = v[1] = v[2] = k;
  }

  vec3u(size_t x, size_t y, size_t z) { 
    v[0] = x; 
    v[1] = y;
    v[2] = z;
  }

  const size_t& x() const { return v[0]; }
  const size_t& y() const { return v[1]; }
  const size_t& z() const { return v[2]; }

  size_t& x() { return v[0]; }
  size_t& y() { return v[1]; }
  size_t& z() { return v[2]; }

  size_t prod() const {
    return v[0]*v[1]*v[2];
  }

  size_t max() const {
    return std::max( v[0], std::max(v[1], v[2]) );
  }

  size_t min() const {
    return std::min( v[0], std::min(v[1], v[2]) );
  }

  static vec3u min(const vec3u& v1, const vec3u& v2) {
    return vec3u( std::min(v1[0], v2[0]),
                  std::min(v1[1], v2[1]),
                  std::min(v1[2], v2[2]) );
  }

  static vec3u max(const vec3u& v1, const vec3u& v2) {
    return vec3u( std::max(v1[0], v2[0]),
                  std::max(v1[1], v2[1]),
                  std::max(v1[2], v2[2]) );
  }


  const size_t& operator[](size_t d) const { return v[d]; }

  size_t& operator[](size_t d) { return v[d]; }

  vec3u& operator+=(const vec3u& v2) {
    v[0] += v2[0];
    v[1] += v2[1];
    v[2] += v2[2];
    return *this;
  }

  vec3u& operator-=(const vec3u& v2) {
    v[0] -= v2[0];
    v[1] -= v2[1];
    v[2] -= v2[2];
    return *this;
  }

  vec3u& operator*=(size_t s) {
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
    return *this;
  }

  static size_t sub2ind(const vec3u& d, const vec3u& s) {
    return ((s[2]*d[1]) + s[1])*d[0] + s[0];
  }

  static size_t sub2ind(const vec3u& d, size_t x, size_t y, size_t z) {
    return ((z*d[1]) + y)*d[0] + x;
  }

  static vec3u ind2sub(const vec3u& d, size_t idx) { 
    size_t xy = d[0]*d[1];
    size_t rem = idx % (xy);
    return vec3u(rem%d[0], rem/d[0], idx/xy);
  }

};

inline vec3u operator+(const vec3u& u, const vec3u& v) {
  return vec3u(u[0] + v[0], u[1] + v[1], u[2] + v[2]);
}

inline vec3u operator-(const vec3u& u, const vec3u& v) {
  return vec3u(u[0] - v[0], u[1] - v[1], u[2] - v[2]);
}

inline vec3u operator*(const vec3u& v, size_t s) {
  return vec3u(v[0]*s, v[1]*s, v[2]*s);
}

inline vec3u operator*(size_t s, const vec3u& v) {
  return vec3u(v[0]*s, v[1]*s, v[2]*s);
}

// lexical sorting for vec3u
inline bool operator==(const vec3u& a, const vec3u& b) {
  return (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]);
}

inline bool operator!=(const vec3u& a, const vec3u& b) {
  return (a[0] != b[0] || a[1] != b[1] || a[2] != b[2]);
}

inline bool operator<(const vec3u& a, const vec3u& b) {
  for (int i=0; i<3; ++i) {
    if (a[i] < b[i]) { return true; }
    else if (a[i] > b[i]) { return false; }
  }
  return false;
}

inline std::ostream& operator<<(std::ostream& ostr, const vec3u& v) {
  return ostr << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
}

#endif
