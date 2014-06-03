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

#ifndef _STRUTILS_H_
#define _STRUTILS_H_

#include <string>
#include <ctype.h>

inline std::string directoryOf(const std::string& filename) {

  size_t pos = filename.rfind('/');

  if (pos == size_t(-1)) { 
    return ".";
  } else {
    return filename.substr(0, pos);
  }

}

inline std::string extensionOf(const std::string& filename) {
  size_t pos = filename.rfind('.');
  if (pos == size_t(-1)) {
    return "";
  } else {
    return filename.substr(pos+1, filename.length()-pos-1);
  }
}

inline std::string filenameOf(const std::string& filename) {
  size_t pos = filename.rfind('/');
  if (pos == size_t(-1)) {
    return filename;
  } else {
    return filename.substr(pos+1, filename.length()-pos-1);
  }
}

inline std::string combineDir(const std::string& dir, const std::string& rel) {
  if (rel.length() && rel[0] == '/') {
    return rel;
  } else {
    return dir + "/" + rel;
  }
}

inline std::string lower(const std::string& s) {
  std::string rval = s;
  for (size_t i=0; i<s.length(); ++i) { rval[i] = tolower(rval[i]); }
  return rval;
}

inline std::string upper(const std::string& s) {
  std::string rval = s;
  for (size_t i=0; i<s.length(); ++i) { rval[i] = toupper(rval[i]); }
  return rval;
}


inline std::string trimws(const std::string& s) {

  if (s.empty()) { return ""; }
 
  size_t p0 = 0; 
  while (p0 < s.length() && isspace(s[p0])) { ++p0; } 

  size_t p1 = s.length()-1; 
  while (p1 > p0 && isspace(s[p1])) { --p1; } 

  return s.substr(p0, p1-p0+1);

}

inline bool split(const std::string& s, char c, 
                  std::string& s1, std::string& s2) {
  
  size_t p = s.find(c);
  if (p == std::string::npos) { return false; }

  s1 = trimws(s.substr(0, p));
  s2 = trimws(s.substr(p+1, s.length()-p-1));

  return true;

}

#endif
