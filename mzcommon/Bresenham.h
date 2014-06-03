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

#ifndef _BRESENHAM_H_
#define _BRESENHAM_H_

#include "vec2u.h"
#include "vec3u.h"
#include <cmath>
#include <algorithm>

template <typename Tpixfunc2>
void bresenham2D(const vec2u& s1, 
                 const vec2u& s2,
                 Tpixfunc2 PUT_PIXEL) {

  vec2u pixel = s1;

  int dx = s2[0]-s1[0];
  int dy = s2[1]-s1[1];
  
  int x_inc = (dx < 0) ? -1 : 1;
  int l = abs(dx);

  int y_inc = (dy < 0) ? -1 : 1;
  int m = abs(dy);

  int dx2 = l << 1;
  int dy2 = m << 1;

  if (l >= m) {

    // we are X-major
    int err = dy2 - l;

    for (int i = 0; i < l; i++) {

      PUT_PIXEL(pixel);

      if (err > 0) {
        pixel[1] += y_inc;
        err -= dx2;
      }

      err += dy2;
      pixel[0] += x_inc;
      
    }

  } else {

    // we are y-major

    int err = dx2 - m;

    for (int i = 0; i < m; i++) {

      PUT_PIXEL(pixel);

      if (err > 0) {
        pixel[0] += x_inc;
        err -= dy2;
      }

      err += dx2;

      pixel[1] += y_inc;

    }

  }

  PUT_PIXEL(pixel);
  
}

template <typename Tpixfunc3>
void bresenham3D(const vec3u& s1, const vec3u& s2,
                 Tpixfunc3 PUT_PIXEL) 
{

  vec3u pixel = s1;

  int dx = s2[0]-s1[0];
  int dy = s2[1]-s1[1];
  int dz = s2[2]-s1[2];
  
  int x_inc = (dx < 0) ? -1 : 1;
  int l = abs(dx);

  int y_inc = (dy < 0) ? -1 : 1;
  int m = abs(dy);

  int z_inc = (dz < 0) ? -1 : 1;
  int n = abs(dz);

  int dx2 = l << 1;
  int dy2 = m << 1;
  int dz2 = n << 1;
  
  if ((l >= m) && (l >= n)) {

    // we are X-major

    int err_1 = dy2 - l;
    int err_2 = dz2 - l;

    for (int i = 0; i < l; i++) {

      PUT_PIXEL(pixel);

      if (err_1 > 0) {
        pixel[1] += y_inc;
        err_1 -= dx2;
      }

      if (err_2 > 0) {
        pixel[2] += z_inc;
        err_2 -= dx2;
      }

      err_1 += dy2;
      err_2 += dz2;
      pixel[0] += x_inc;
      
    }

  } else if ((m >= l) && (m >= n)) {

    // we are y-major

    int err_1 = dx2 - m;
    int err_2 = dz2 - m;

    for (int i = 0; i < m; i++) {

      PUT_PIXEL(pixel);

      if (err_1 > 0) {
        pixel[0] += x_inc;
        err_1 -= dy2;
      }

      if (err_2 > 0) {
        pixel[2] += z_inc;
        err_2 -= dy2;
      }

      err_1 += dx2;
      err_2 += dz2;

      pixel[1] += y_inc;

    }

  } else {

    // we are z-major

    int err_1 = dy2 - n;
    int err_2 = dx2 - n;

    for (int i = 0; i < n; i++) {
      
      PUT_PIXEL(pixel);

      if (err_1 > 0) {
        pixel[1] += y_inc;
        err_1 -= dz2;
      }

      if (err_2 > 0) {
        pixel[0] += x_inc;
        err_2 -= dz2;
      }

      err_1 += dy2;
      err_2 += dx2;
      pixel[2] += z_inc;

    }

  }

  PUT_PIXEL(pixel);

}

#endif
