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

#ifndef _ANGLEUTIL_H_
#define _ANGLEUTIL_H_

#include <math.h>

template <class TAng> 
inline TAng clamp_angle(TAng angle, bool twopi=false) {

  const TAng lower = twopi ? 0 : -M_PI;
  const TAng upper = twopi ? 2*M_PI : M_PI;

  while (angle < lower) { angle += 2*M_PI; }
  while (angle >= upper) { angle -= 2*M_PI; }

  return angle;

}

template <class TAng>
inline TAng delta_angle(TAng angle2, TAng angle1) {
  TAng diff = angle2-angle1;
  while (diff > M_PI) { diff -= 2*M_PI; }
  while (diff < -M_PI) { diff += 2*M_PI; }
  return diff;
}

template <class TAng, class TInterp>
inline TAng interp_angle(TAng angle1, TAng angle2, TInterp u, bool twopi=false) {
  TAng diff = delta_angle(angle2, angle1);
  return clamp_angle(angle1+u*diff, twopi);
}

#endif
