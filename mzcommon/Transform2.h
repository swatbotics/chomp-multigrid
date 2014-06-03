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

#ifndef _MZTRANSFORM2_H_
#define _MZTRANSFORM2_H_

#include "mat2.h"
#include "AngleUtil.h"

template <class real> class Transform2_t {
public:

  typedef vec2_t<real> vec2;
  typedef mat2_t<real> mat2;

  vec2 translation;

  Transform2_t(): 
    translation(0), _rotation(0), _rotDirty(false), _c(1), _s(0) {}

  Transform2_t(const vec2& tx):
    translation(tx), _rotation(0), _rotDirty(false), _c(1), _s(0) {}

  Transform2_t(const vec2& tx, real rot): 
    translation(tx), _rotation(rot), _rotDirty(true) {}

  Transform2_t(real rot): 
    translation(0), _rotation(rot), _rotDirty(true) {}


  template <class other>
  explicit Transform2_t(const Transform2_t<other>& t):
    translation(t.translation), _rotation(t.rotation()), _rotDirty(true) {}

  static Transform2_t lerp(const Transform2_t& a,
                           const Transform2_t& b,
                           real u) {
    return Transform2_t( a.translation + u*(b.translation-a.translation),
                         interp_angle(a.rotation(), b.rotation(), u) );
  }

  vec2 transformFwd(const vec2& p) const {
    return rotFwd() * p + translation;
  }

  vec2 transformInv(const vec2& p) const {
    return rotInv() * (p - translation);
  }
  
  Transform2_t inverse() const {
    return Transform2_t( transformInv(vec2(0,0)), -rotation() );
  }
  
  mat2 rotFwd() const {
    if (_rotDirty) { 
      _c = cos(_rotation);
      _s = sin(_rotation);
      _rotDirty = false; 
    }
    return mat2(_c, -_s, _s, _c);
  }

  mat2 rotInv() const {
    if (_rotDirty) { 
      _c = cos(_rotation);
      _s = sin(_rotation);
      _rotDirty = false; 
    }
    return mat2(_c, _s, -_s, _c);
  }

  real rotation() const { return _rotation; }

  void setRotation(real r) { _rotation = r; _rotDirty = true; }

private:

  real _rotation;
  mutable bool _rotDirty;
  mutable real _c, _s;

};

template <class real> inline vec2_t<real>
operator*(const Transform2_t<real>& tx, const vec2_t<real>& v) {
  return tx.transformFwd(v);
}

template <class real> inline Transform2_t<real>
operator*(const Transform2_t<real>& t1, const Transform2_t<real>& t2) {
  real r2 = clamp_angle(t1.rotation() + t2.rotation());
  return Transform2_t<real>(t1.transformFwd(t2.translation), r2);
}

typedef Transform2_t<double> Transform2d;
typedef Transform2_t<float> Transform2f;

template <class real> 
inline std::ostream& operator<<(std::ostream& ostr, const Transform2_t<real>& t) {
  return ostr << "< " << t.translation.x() 
              << ", " << t.translation.y() 
              << ", " << t.rotation() * 180 / M_PI 
              << ">";
}

#endif

