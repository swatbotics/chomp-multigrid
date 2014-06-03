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

#ifndef _SIMPLETOKENIZER_H_
#define _SIMPLETOKENIZER_H_

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>

class SimpleTokenizer {
public:

  std::istream& istr;

  SimpleTokenizer(std::istream& i): istr(i) {}

  void consumeWhitespaceAndComments() {
    
    bool in_comment = false;

    while (1) {
      int c = istr.peek();
      if (c == EOF) {
	break;
      } else if (in_comment) {
	istr.get();
	if (c == '\n') { in_comment = false; }
      } else if (c == '#') { 
	istr.get();
	in_comment = true;
      } else if (isspace(c)) {
	istr.get();
      } else {
	break;
      }
    }

  }

  void parseLiteral(const std::string& literal) {
    consumeWhitespaceAndComments();
    std::string str;
    if (!(istr >> str) || str != literal) {
      std::cerr << "error: expected " << literal << "\n";
      exit(1);
    }
  }

  std::string parseString(bool allowEOF) {
    consumeWhitespaceAndComments();
    std::string str;
    if (!(istr >> str)) {
      if (!allowEOF) {
	std::cerr << "error: unexpected EOF\n";
	exit(1);
      }
    }
    return str;
  }

  template <class Tnumber>
  Tnumber parseNumber() {
    consumeWhitespaceAndComments();
    Tnumber d;
    if (!(istr >> d)) {
      std::cerr << "error parsing number!\n";
      exit(1);
    }
    return d;
  }

  template <class Tnumber>
  void parseNumber(Tnumber& value) {
    value = parseNumber<Tnumber>();
  }

  int parseFromList(const char* const* strs) {
    consumeWhitespaceAndComments();
    std::string s = parseString(false);
    for (int i=0; strs[i]; ++i) {
      if (strs[i] == s) { 
	return i;
      }
    }
    std::cerr << "unrecognized item: " << s << "\n";
    exit(1);
  }



  double parseLengthToMeters() {

    double quantity = parseNumber<double>();

    static const char* const units[] =  { "m", "cm", "mm", "ft", "in", 0 };
    static const double scales[] = { 1, 0.01, 0.001, 0.3048, 0.0254, 0 };

    int idx = parseFromList(units);

    return quantity * scales[idx];

  }

  double parseAngleToRadians() {

    double quantity = parseNumber<double>();

    static const char* const units[] =  { "rad", "deg" };
    static const double scales[] = { 1, M_PI/180.0 };

    int idx = parseFromList(units);

    return quantity * scales[idx];

  }

};

#endif
