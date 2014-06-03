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

#include "GlTexture2D.h"  
#include <iostream>

/*
static const char* typestr(GLenum type) {
  return (type == GL_TEXTURE_RECTANGLE_ARB) ? 
    "GL_TEXTURE_RECTANGLE_ARB" : "GL_TEXTURE_2D";
}
*/

class Tex2DBinder {
public:

  GlTexture2D& t;
  GLint prev;
  bool restore;

  Tex2DBinder(GlTexture2D& tex): 
    t(tex), restore(false) 
  {
    if (t.handle()) {
      GLenum binding = (t.type() == GL_TEXTURE_RECTANGLE_ARB) ?
        GL_TEXTURE_BINDING_RECTANGLE_ARB : GL_TEXTURE_BINDING_2D;
      glGetIntegerv(binding, &prev);
      if (GLuint(prev) != t.handle()) {
        restore = true;
        //std::cerr << "calling bind " << typestr(t.type()) << ", " << t.handle() << " in ctor\n";
        glBindTexture(t.type(), t.handle());
      }
    } 
  }

  ~Tex2DBinder() { 
    if (restore) {
      //std::cerr << "calling bind " << typestr(t.type()) << ", " << prev << " in dtor\n";
      glBindTexture(t.type(), prev); 
    }
  }

};

GlTexture2D::GlTexture2D(): _texname(0), _type(Texture2D) {}

GlTexture2D::GlTexture2D(GLenum type): _texname(0), _type(Texture2D) {
  createTexture(type);
}

GlTexture2D::~GlTexture2D() {
  deleteTexture();
}


void GlTexture2D::texImage(GLint level,
                           GLint internalFormat,
                           GLsizei width, GLsizei height,
                           GLint border, 
                           GLenum format,
                           GLenum type,
                           const GLvoid* pixels) {
  Tex2DBinder b(*this);
  glTexImage2D(_type, level, internalFormat, width, height, border, format, type, pixels);
}

void  GlTexture2D::setParameter(GLenum pname, GLint pvalue) {
  Tex2DBinder b(*this);
  glTexParameteri(_type, pname, pvalue);
}

GLint GlTexture2D::getParameter(GLenum pname) {
  Tex2DBinder b(*this);
  GLint pvalue;
  glGetTexParameteriv(_type, pname, &pvalue);
  return pvalue;
}

/*
void GlTexture2D::generateMipmaps() {
  Tex2DBinder b(*this);
  glGenerateMipmapEXT(_type);
}
*/

void GlTexture2D::createTexture(GLenum type) {
  deleteTexture();
  glGenTextures(1, &_texname);
  _type = type;
}

void GlTexture2D::deleteTexture() {
  if (_texname && glIsTexture(_texname)) {
    glDeleteTextures(1, &_texname);
  }
  _texname = 0;
}

void GlTexture2D::release() {
  _texname = 0;
}

GLuint GlTexture2D::handle() {
  return _texname;
}

GLenum GlTexture2D::type() const {
  return _type;
}

void GlTexture2D::enable() {
  //std::cerr << "enable " << typestr(_type) << "\n";
  glEnable(_type);
}

void GlTexture2D::disable() {
  //std::cerr << "disable " << typestr(_type) << "\n";
  glDisable(_type);
}

void GlTexture2D::bind() {
  //std::cerr << "bind " << typestr(_type) << ", " << _texname << "\n";
  glBindTexture(_type, _texname);
}

void GlTexture2D::unbind() {
  //std::cerr << "bind " << typestr(_type) << ", " << 0 << "\n";
  glBindTexture(_type, 0);
}

GLint GlTexture2D::getMaxTextureUnits() {
  GLint rval;
  glGetIntegerv(GL_MAX_TEXTURE_UNITS, &rval);
  return rval;
}

