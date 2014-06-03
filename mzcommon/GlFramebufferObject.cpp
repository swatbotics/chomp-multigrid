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

#include "GlFramebufferObject.h"


GlRenderbufferObject::GlRenderbufferObject(): _handle(0) {}

GlRenderbufferObject::GlRenderbufferObject(GLint internalFormat,
                                           GLsizei width, GLsizei height): _handle(0) {
  createObject(internalFormat, width, height);
}

GlRenderbufferObject::~GlRenderbufferObject() {
  deleteObject();
}

void GlRenderbufferObject::createObject(GLint internalFormat,
                                        GLsizei width, GLsizei height) {
  deleteObject();
  glGenRenderbuffersEXT(1, &_handle);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _handle);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, internalFormat, width, height);
}

void GlRenderbufferObject::deleteObject() {
  if (_handle && glIsRenderbufferEXT(_handle)) {
    glDeleteRenderbuffersEXT(1, &_handle);
  }
  _handle = 0;
}

void GlRenderbufferObject::release() {
  _handle = 0; 
}

GLuint GlRenderbufferObject::handle() {
  return _handle;
}

void GlRenderbufferObject::bind() {
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _handle);
}

void GlRenderbufferObject::unbind() {
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
}

//////////////////////////////////////////////////////////////////////

class FboBinder {
public:

  GLint prev;
  bool restore;

  FboBinder(GLuint handle): prev(0), restore(false) {
    glGetIntegerv( GL_FRAMEBUFFER_BINDING_EXT, &prev );
    if ( GLuint(prev) != handle ) {
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, handle );
      restore = true;
    }
  }

  ~FboBinder() {
    if (restore) {
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, prev);
    }
  }
    

};

GlFramebufferObject::GlFramebufferObject(): _handle(0) { }

GlFramebufferObject::GlFramebufferObject(bool init): _handle(0) {
  if (init) { createObject(); }
}

GlFramebufferObject::~GlFramebufferObject() {
  deleteObject();
}

void GlFramebufferObject::createObject() {
  deleteObject();
  glGenFramebuffersEXT(1, &_handle);
}

void GlFramebufferObject::deleteObject() {
  if (_handle && glIsFramebufferEXT(_handle)) {
    glDeleteFramebuffersEXT(1, &_handle);
  }
  _handle = 0;
}

void GlFramebufferObject::release() {
  _handle = 0;
}

GLuint GlFramebufferObject::handle() {
  return _handle;
}

void GlFramebufferObject::bind() {
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _handle);
}
  
void GlFramebufferObject::attach(GlTexture2D& texture,
                                 GLenum attachment,
                                 GLint miplevel) {
  attachTexture2D(texture.type(), texture.handle(), attachment, miplevel);
}

void GlFramebufferObject::attachTexture2D(GLenum target, GLuint texname,
                                          GLenum attachment,
                                          GLint mipLevel) {
  FboBinder b(_handle);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, attachment,
                            target, texname, mipLevel);
}

void GlFramebufferObject::attach(GlRenderbufferObject& rbuf,
                                 GLenum attachment) {
  attachRenderBuffer(rbuf.handle(), attachment);
}

void GlFramebufferObject::attachRenderBuffer(GLuint id, GLenum attachment) {
  FboBinder b(_handle);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, attachment,
                               GL_RENDERBUFFER_EXT, id);
}

void GlFramebufferObject::detach(GLenum attachment) {

  FboBinder b(_handle);
  GLenum type = getAttachmentType(attachment);

  switch (type) {
  case GL_RENDERBUFFER_EXT:
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, attachment, GL_RENDERBUFFER_EXT, 0);
    break;
  case GL_TEXTURE:
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, attachment, GL_TEXTURE_2D, 0, 0);
    break;
  default:
    break;
  }

}

GLenum GlFramebufferObject::status() const {
  FboBinder b(_handle);
  GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  return status;
}

const char* GlFramebufferObject::statusString(GLenum status) {

  switch (status) {
  case GL_FRAMEBUFFER_COMPLETE_EXT: // Everything's OK
    return "complete";
  case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
    return "incomplete attachment";
  case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
    return "no attachments present";
  case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
    return "attachment dimension mismatch";
  case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
    return "color attachment format mismatch";
  case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
    return "draw buffer incomplete";
  case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
    return "read buffer incomplete";
  case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
    return "unsupported combination of internal formats";
  default:
    break;
  };

  return "UNKNOWN";

}

bool GlFramebufferObject::valid() const{
  return status() == GL_FRAMEBUFFER_COMPLETE_EXT;
}

GLenum GlFramebufferObject::getAttachmentType(GLenum attachment) {

  FboBinder b(_handle);

  GLint type = 0;
  glGetFramebufferAttachmentParameterivEXT(GL_FRAMEBUFFER_EXT, attachment,
                                           GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE_EXT, 
                                           &type);
  return GLenum(type);

}

GLuint GlFramebufferObject::getAttachmentID(GLenum attachment) {

  FboBinder b(_handle);
  
  GLint id = 0;
  glGetFramebufferAttachmentParameterivEXT(GL_FRAMEBUFFER_EXT, attachment,
                                           GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT,
                                           &id);
  return GLuint(id);
  
}

GLint GlFramebufferObject::getMaxColorAttachments() {
  GLint maxAttach = 0;
  glGetIntegerv( GL_MAX_COLOR_ATTACHMENTS_EXT, &maxAttach );
  return maxAttach;
}

void GlFramebufferObject::disable() {
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}
