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

#ifndef _MZGLUTAPP_H_
#define _MZGLUTAPP_H_

#include "GlCamera.h"

#ifdef __APPLE__
#include <Glut/glut.h>
#else
#include <GL/glut.h>
#endif

#include <cstdlib>

class MzGlutApp {
public:

  enum {
    DEFAULT_PARAMS = GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB,
    NICE_PARAMS = DEFAULT_PARAMS | GLUT_MULTISAMPLE
  };

  int width;
  int height;
  GlCamera camera;
  GlCamera::MouseMode mmode;

  virtual ~MzGlutApp();

  void initWindowPosition(int x, int y);
  void initWindowSize(int w, int h);
  void createWindow(const char* name);
  void setTimer(unsigned int msec, int value);
  void setIdle();
  void clearIdle();

  void setupBasicLight(const vec4f& position=vec4f(-1,1,1,0));

  virtual void run();

  virtual void display();
  virtual void reshape(int width, int height);
  virtual void keyboard(unsigned char key, int x, int y);
  virtual void special(int key, int x, int y);
  virtual void mouse(int button, int state, int x, int y);
  virtual void motion(int x, int y);
  virtual void timer(int value);
  virtual void idle();
  virtual void viewChanged(bool interactive);

  void drawString(int rasterx, int rastery, 
                  const std::string& string,
                  void* font = GLUT_BITMAP_8_BY_13,
                  int linespacing = 13);

  virtual GlCamera::MouseMode btn2cam(int button) const;

  static MzGlutApp& getInstance();

protected:
  MzGlutApp(int argc, char** argv, int displaymode=DEFAULT_PARAMS);

private:

  static MzGlutApp* _instance;

};

#endif
