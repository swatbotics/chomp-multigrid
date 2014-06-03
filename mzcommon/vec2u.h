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

#ifndef _VEC2U_H_
#define _VEC2U_H_

#include <iostream>

struct vec2u {

  size_t v[2];

  vec2u() {}

  explicit vec2u(size_t k) {
    v[0] = v[1] = k;
  }

  vec2u(size_t x, size_t y) {
    v[0] = x; 
    v[1] = y;
  }

  const size_t& x() const { return v[0]; }
  const size_t& y() const { return v[1]; }

  size_t& x() { return v[0]; }
  size_t& y() { return v[1]; }

  size_t prod() const {
    return v[0]*v[1];
  }

  size_t min() const {
    return std::min( v[0], v[1] );
  }

  size_t max() const {
    return std::max( v[0], v[1] );
  }

  static vec2u min(const vec2u& v1, const vec2u& v2) {
    return vec2u( std::min(v1[0], v2[0]),
                  std::min(v1[1], v2[1]) );
  }

  static vec2u max(const vec2u& v1, const vec2u& v2) {
    return vec2u( std::max(v1[0], v2[0]),
                  std::max(v1[1], v2[1]) );
  }

  const size_t& operator[](size_t d) const { return v[d]; }

  size_t& operator[](size_t d) { return v[d]; }

  vec2u& operator+=(const vec2u& v2) {
    v[0] += v2[0];
    v[1] += v2[1];
    return *this;
  }

  vec2u& operator-=(const vec2u& v2) {
    v[0] -= v2[0];
    v[1] -= v2[1];
    return *this;
  }

  vec2u& operator*=(size_t s) {
    v[0] *= s;
    v[1] *= s;
    return *this;
  }

  static size_t sub2ind(const vec2u& d, const vec2u& s) {
    return s[1]*d[0] + s[0];
  }

  static size_t sub2ind(const vec2u& d, size_t x, size_t y) {
    return y*d[0] + x;
  }

  static vec2u ind2sub(const vec2u& d, size_t idx) { 
    return vec2u(idx%d[0], idx/d[0]);
  }

};

inline vec2u operator+(const vec2u& u, const vec2u& v) {
  return vec2u(u[0] + v[0], u[1] + v[1]);
}

inline vec2u operator-(const vec2u& u, const vec2u& v) {
  return vec2u(u[0] - v[0], u[1] - v[1]);
}

inline vec2u operator*(const vec2u& v, size_t s) {
  return vec2u(v[0]*s, v[1]*s);
}

inline vec2u operator*(size_t s, const vec2u& v) {
  return vec2u(v[0]*s, v[1]*s);
}

// lexical sorting for vec2u
inline bool operator==(const vec2u& a, const vec2u& b) {
  return (a[0] == b[0] && a[1] == b[1]);
}

inline bool operator!=(const vec2u& a, const vec2u& b) {
  return (a[0] != b[0] || a[1] != b[1]);
}

inline bool operator<(const vec2u& a, const vec2u& b) {
  for (int i=0; i<2; ++i) {
    if (a[i] < b[i]) { return true; }
    else if (a[i] > b[i]) { return false; }
  }
  return false;
}

inline std::ostream& operator<<(std::ostream& ostr, const vec2u& v) {
  return ostr << "[" << v[0] << ", " << v[1] << "]";
}

#endif
