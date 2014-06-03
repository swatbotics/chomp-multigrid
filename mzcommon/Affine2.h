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

#ifndef _AFFINE2_H_
#define _AFFINE2_H_

#include "Transform2.h"

/** Represent a 2D affine transformation. */
template <class real>
class Affine2_t {
public:

  typedef mat2_t<real> mat2;
  typedef vec2_t<real> vec2;
  typedef Transform2_t<real> Transform2;

  /** The matrix part of the affine transform. */
  mat2 m;

  /** The vector part of the affine transform. */
  vec2 v;

  /** Constructs the identity transform. */
  Affine2_t(): m(mat2::identity()), v(0) {}

  /** Construct from matrix and vector. */
  Affine2_t(const mat2& matrix, const vec2& vector): m(matrix), v(vector) {}

  /** Construct from Transform2. */
  Affine2_t(const Transform2& t): m(t.rotFwd()), v(t.translation) {}

  /** Construct from matrix only. */
  explicit Affine2_t(const mat2& matrix): m(matrix), v(0) {}

  /** Construct from vector only. */
  explicit Affine2_t(const vec2& vector): m(mat2::identity()), v(vector) {}

  /** Returns the inverse of this affine transform.  The calculation
   *  involves dividing by the determinant returned by det(). If the
   *  determinant is zero then the Affine2_t is not invertible and the
   *  inverse will contain NaN entries.
   */
  Affine2_t inverse() const {
    Affine2_t rval;
    real det = m.determinant();
    m.inverse(rval.m, det);
    rval.v = -(rval.m * v);
    return rval;
  }

  /** Return the identity transform. */
  static Affine2_t identity() { return Affine2_t(); }

  /** Return an Affine2_t which translates points by (@a tx, @a ty) */
  static Affine2_t translation(real tx, real ty) {
    return Affine2_t(vec2(tx, ty));
  }

  /** Return an Affine2_t which translates points by @a t. */
  static Affine2_t translation(const vec2& t) {
    return Affine2_t(t);
  }

  /** Return an Affine2_t which rotates points around the origin by 
   *  an angle of @a theta radians. 
   */
  static Affine2_t rotation(real theta) {
    return Affine2_t(mat2::rotMat(theta));
  }

  /** Return an Affine2_t which scales points by (@a sx, @a sy). */
  static Affine2_t scale(real sx, real sy) {
    Affine2_t rval;
    rval.m(0,0) = sx;  rval.m(0,1) = 0;
    rval.m(1,0) = 0;   rval.m(1,1) = sy;
    return rval;
  }

  /** Return the point which is the product of this Affine2_t and @a p. */
  vec2 transform(const vec2& p) const {
    return m*p + v;
  }

};

/** @relates Affine2_t Compose two affine transformations. */
template <class real>
inline Affine2_t<real> operator*(const Affine2_t<real>& a, 
                                 const Affine2_t<real>& b) {
  return Affine2_t<real>(a.m * b.m, a.transform(b.v));
}

/** @relates Affine2_t Syntactic sugar. */
template <class real>
inline vec2_t<real> operator*(const Affine2_t<real>& a, 
                              const vec2_t<real> p) {
  return a.transform(p);
}

typedef Affine2_t<double> Affine2d;
typedef Affine2_t<float>  Affine2f;

#endif
