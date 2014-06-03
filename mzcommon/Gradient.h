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

#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include "vec3.h"
#include <map>

class Gradient {
public:

  typedef std::map<float, vec3f> StopMap;
  StopMap stops;
  
  vec3f lookup(float u) const {
    
    if (stops.empty()) {
      return vec3f(0,0,0);
    }

    StopMap::const_iterator i = stops.begin();
    StopMap::const_iterator j = i;
    ++j;

    if (u < i->first || j == stops.end()) { return i->second; }

    for ( ; j!=stops.end(); ++i, ++j) {
      float min = i->first;
      float max = j->first;
      if (u >= min && u <= max) {
        if (max == min) { return i->second; }
        u = (u-min)/(max-min);
        return vec3f::lerp(i->second, j->second, u);
      }
    }

    return i->second;

  }

  static Gradient rainbow() {
    Gradient g;
    const float a = 1.0f/7.0f;
    g.stops[0*a] = vec3f(0.5, 0.0, 1.0); // pur
    g.stops[1*a] = vec3f(0.0, 0.0, 1.0); // bl
    g.stops[2*a] = vec3f(0.0, 0.5, 1.0); // bc
    g.stops[3*a] = vec3f(0.0, 1.0, 1.0); // cy
    g.stops[4*a] = vec3f(0.0, 1.0, 0.0); // g
    g.stops[5*a] = vec3f(1.0, 1.0, 0.0); // y
    g.stops[6*a] = vec3f(1.0, 0.5, 0.0); // or
    g.stops[7*a] = vec3f(1.0, 0.0, 0.0); // red
    return g;
  }

  static Gradient jet() {

    Gradient g;

    const float a = 1.0f/15.0f;

    static const float stuff[16*3] = {
      0.0000,   0.0000,   0.7500,
      0.0000,   0.0000,   1.0000,
      0.0000,   0.2500,   1.0000,
      0.0000,   0.5000,   1.0000,
      0.0000,   0.7500,   1.0000,
      0.0000,   1.0000,   1.0000,
      0.2500,   1.0000,   0.7500,
      0.5000,   1.0000,   0.5000,
      0.7500,   1.0000,   0.2500,
      1.0000,   1.0000,   0.0000,
      1.0000,   0.7500,   0.0000,
      1.0000,   0.5000,   0.0000,
      1.0000,   0.2500,   0.0000,
      1.0000,   0.0000,   0.0000,
      0.7500,   0.0000,   0.0000,
      0.5000,   0.0000,   0.0000,
    };

    int offs = 0;

    for (int i=0; i<16; ++i) {
      g.stops[a*i] = vec3f(stuff[offs+0], stuff[offs+1], stuff[offs+2]);
      offs += 3;
    }

    return g;

  }


  static Gradient bone() {
    Gradient g;
    g.stops[0.0f] = vec3f(0.0f, 0.0f, 0.0f);
    g.stops[1.0f/3.0f] = (7*vec3f(1.0f/3.0f) + vec3f(0.0f, 0.0f, 1.0f)) / 8;
    g.stops[2.0f/3.0f] = (7*vec3f(2.0f/3.0f) + vec3f(0.0f, 1.0f, 1.0f)) / 8;
    g.stops[1.0f] = vec3f(1.0f, 1.0f, 1.0f);
    return g;
  }

  static Gradient summer() {
    Gradient g;
    g.stops[0.0f] = vec3f(0.0f, 0.5f, 0.4f);
    g.stops[1.0f] = vec3f(1.0f, 1.0f, 0.4f);
    return g;
  }

  static Gradient hot() {
    Gradient g;
    g.stops[0.0f] = vec3f(0.0f, 0.0f, 0.0f);
    g.stops[1.0f/3.0f] = vec3f(1.0f, 0.0f, 0.0f);
    g.stops[2.0f/3.0f] = vec3f(1.0f, 1.0f, 0.0f);
    g.stops[1.0f] = vec3f(1.0f, 1.0f, 1.0f);
    return g;
  }

  static Gradient bw() {
    Gradient g;
    g.stops[0.0f] = vec3f(0,0,0);
    g.stops[1.0f] = vec3f(1,1,1);
    return g;
  }

};

#endif
