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

#ifndef _GLFRAMEBUFFEROBJECT_H_
#define _GLFRAMEBUFFEROBJECT_H_

#include <GL/glew.h>
#include "GlTexture2D.h"

//////////////////////////////////////////////////////////////////////

class GlRenderbufferObject {
public:

  GlRenderbufferObject();
  GlRenderbufferObject(GLint internalFormat,
                       GLsizei width, GLsizei height);
  ~GlRenderbufferObject();

  void createObject(GLint internalFormat,
                    GLsizei width, GLsizei height);
  void deleteObject();
  void release();
  GLuint handle();

  void bind();
  static void unbind();
  
private:

  GLuint _handle;

};

//////////////////////////////////////////////////////////////////////

class GlFramebufferObject {
public:

  GlFramebufferObject();
  GlFramebufferObject(bool init);

  ~GlFramebufferObject();

  void createObject();
  void deleteObject();
  void release();
  GLuint handle();

  void bind();
  
  void attach(GlTexture2D& texture,
              GLenum attachment=GL_COLOR_ATTACHMENT0_EXT,
              GLint miplevel=0);

  void attach(GlRenderbufferObject& rbuf,
              GLenum attachment=GL_DEPTH_ATTACHMENT_EXT);

  void attachTexture2D(GLenum target, GLuint texname,
                       GLenum attachment = GL_COLOR_ATTACHMENT0_EXT,
                       GLint mipLevel=0);

  void attachRenderBuffer(GLuint id, GLenum attachment=GL_DEPTH_ATTACHMENT_EXT);

  void detach(GLenum attachment);
  GLenum status() const;
  bool valid() const;

  static const char* statusString(GLenum status);

  GLenum getAttachmentType(GLenum attachment);
  GLuint getAttachmentID(GLenum attachment);

  static GLint getMaxColorAttachments();
  static void disable();
  
private:

  GLuint _handle;

};

#endif
