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

#ifndef _GLTEXTURE2D_H_
#define _GLTEXTURE2D_H_


#include "glstuff.h"

class GlTexture2D {
public:

  enum TargetType {
    Texture2D=GL_TEXTURE_2D,
    TextureRectangle=GL_TEXTURE_RECTANGLE_ARB
  };

  enum ParamName {
    WrapS=GL_TEXTURE_WRAP_S,
    WrapT=GL_TEXTURE_WRAP_T,
    MagFilter=GL_TEXTURE_MAG_FILTER,
    MinFilter=GL_TEXTURE_MIN_FILTER
  };

  enum FilterType {
    Nearest=GL_NEAREST,
    Linear=GL_LINEAR,
    NearestMipmapNearest=GL_NEAREST_MIPMAP_NEAREST,
    LinearMipmapNearest=GL_LINEAR_MIPMAP_NEAREST,
    LinearMipmapLinear=GL_LINEAR_MIPMAP_LINEAR
  };

  enum WrapType {
    Clamp=GL_CLAMP,
    Repeat=GL_REPEAT,
    ClampToEdge=GL_CLAMP_TO_EDGE,
  };

  GlTexture2D();
  GlTexture2D(GLenum type);
  ~GlTexture2D();

  void loadPNG(const char* filename);
  void loadPNG(const std::string& filename);

  void texImage(GLint level,
                GLint internalFormat,
                GLsizei width, GLsizei height,
                GLint border, 
                GLenum format,
                GLenum type,
                const GLvoid* pixels);

  void  setParameter(GLenum pname, GLint pvalue);
  GLint getParameter(GLenum pname);

  /*
  void generateMipmaps();
  */

  void createTexture(GLenum type);
  void deleteTexture();
  void release();
  GLuint handle();
  GLenum type() const;

  void enable();
  void disable();

  void bind();
  void unbind();

  static GLint getMaxTextureUnits();

private:

  GLuint _texname;
  GLenum _type;

};

#endif
