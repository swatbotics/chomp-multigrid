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

#include "MzGlutApp.h"
#include <stdexcept>

MzGlutApp* MzGlutApp::_instance = 0;

MzGlutApp& MzGlutApp::getInstance() { return *_instance; }

static void _reshape(int w, int h) { 
  MzGlutApp::getInstance().reshape(w, h); 
}

static void _display() { 
  MzGlutApp::getInstance().display(); 
}

static void _keyboard(unsigned char key, int x, int y) { 
  MzGlutApp::getInstance().keyboard(key, x, y); 
}

static void _special(int key, int x, int y) { 
  MzGlutApp::getInstance().special(key, x, y); 
}

static void _mouse(int button, int state, int x, int y) {
  MzGlutApp::getInstance().mouse(button, state, x, y); 
}

static void _motion(int x, int y) { 
  MzGlutApp::getInstance().motion(x,y); 
}

static void _timer(int value) {
  MzGlutApp::getInstance().timer(value);
}

static void _idle() {
  MzGlutApp::getInstance().idle();
}

GlCamera::MouseMode MzGlutApp::btn2cam(int button) const {
  switch (button) {
  case GLUT_LEFT_BUTTON: {
    int mod = glutGetModifiers();
    switch (mod) {
    case 0:
      return GlCamera::MOUSE_ROTATE;
    case GLUT_ACTIVE_CTRL:
      return GlCamera::MOUSE_ZOOM;
    case GLUT_ACTIVE_SHIFT:
      return GlCamera::MOUSE_PAN_XY;
    default:
      break;
    }
    return GlCamera::MOUSE_NONE;
  }
  case GLUT_MIDDLE_BUTTON:
    return GlCamera::MOUSE_ZOOM;
  case GLUT_RIGHT_BUTTON:
    return GlCamera::MOUSE_PAN_XY;
  default:
    break;
  }
  return GlCamera::MOUSE_NONE;
}

MzGlutApp::MzGlutApp(int argc, char** argv, int displaymode) {
  if (_instance) { throw std::runtime_error("GLUT already initialized!"); }
  _instance = this;
  glutInit(&argc, argv);
  if (displaymode) { glutInitDisplayMode(displaymode); }
  width = 640;
  height = 480;
  mmode = GlCamera::MOUSE_NONE;
}

MzGlutApp::~MzGlutApp() { 
  _instance = 0;
}

void MzGlutApp::initWindowPosition(int x, int y) {
  glutInitWindowPosition(x, y);
}

void MzGlutApp::initWindowSize(int w, int h) {
  width = w;
  height = h;
  glutInitWindowSize(w, h);
}

void MzGlutApp::createWindow(const char* name) {

  glutCreateWindow(name);

  glutReshapeFunc(_reshape);
  glutDisplayFunc(_display);
  glutKeyboardFunc(_keyboard);
  glutSpecialFunc(_special);
  glutMouseFunc(_mouse);
  glutMotionFunc(_motion);

  

  glClearColor(1,1,1,1);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glFrontFace(GL_CCW);

}

void MzGlutApp::setTimer(unsigned int msec, int value) {
  glutTimerFunc(msec, _timer, value);
}

void MzGlutApp::setIdle() { glutIdleFunc(_idle); }

void MzGlutApp::clearIdle() { glutIdleFunc(0); }

void MzGlutApp::setupBasicLight(const vec4f& position) {

  glMatrixMode(GL_PROJECTION);
  camera.loadMatrix(GlCamera::MATRIX_PROJECTION);

  glMatrixMode(GL_MODELVIEW);
  camera.loadMatrix(GlCamera::MATRIX_MODELVIEW);

  const GLfloat lightdim[4] = { 0.3, 0.3, 0.3, 1.0 };
  const GLfloat lightbrt[4] = { 0.7, 0.7, 0.7, 1.0 };
  const GLfloat white[4] = { 1, 1, 1, 1 };

  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_LIGHT0);
  
  glLightfv(GL_LIGHT0, GL_POSITION, position.v);
  glLightfv(GL_LIGHT0, GL_AMBIENT,  lightdim);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  lightbrt);
  glLightfv(GL_LIGHT0, GL_SPECULAR, white);
  
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, white);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 120.0f);
  
  glEnable(GL_LIGHTING);


}

void MzGlutApp::display() {

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


}

void MzGlutApp::reshape(int w, int h) {

  glClearColor(1,1,1,1);
  
  width = w;
  height = h;

  glViewport(0, 0, width, height);
  camera.setViewport(0, 0, width, height);

  glMatrixMode(GL_PROJECTION);
  camera.loadMatrix(GlCamera::MATRIX_PROJECTION);

  viewChanged(true);

}

void MzGlutApp::keyboard(unsigned char key, int x, int y) { 

  switch (key) { 
  case 27: // ESC
    exit(0);
    break;
  case ' ': // SPACE
    camera.recallHomePosition();
    viewChanged(false);
    break;
  case 'w':
  case 'W':
    camera.pan(0, 0, -25, true);
    viewChanged(false);
    break;
  case 's':
  case 'S':
    camera.pan(0, 0, 25, true);
    viewChanged(false);
    break;
  case 'a':
  case 'A':
    camera.pan(-25, 0, 0, true);
    viewChanged(false);
    break;
  case 'd':
  case 'D':
    camera.pan(25, 0, 0, true);
    viewChanged(false);
    break;
  }

}

void MzGlutApp::special(int key, int x, int y) { }

void MzGlutApp::mouse(int button, int state, int x, int y) {

  if (state == GLUT_DOWN) {
    mmode = btn2cam(button);
    camera.mousePress(x,y,mmode);
  } else {
    camera.mouseRelease(x,y,mmode);
    mmode = GlCamera::MOUSE_NONE;
  }

  viewChanged(state == GLUT_DOWN);

  
}

void MzGlutApp::motion(int x, int y) {

  if (mmode != GlCamera::MOUSE_NONE) {
    camera.mouseMove(x,y,mmode);
    viewChanged(true);
  }
  
}

void MzGlutApp::timer(int value) {
}

void MzGlutApp::idle() {
  
}

void MzGlutApp::run() {

  // run main loop
  glutMainLoop();
  
}

void MzGlutApp::viewChanged(bool interactive) {

  glMatrixMode(GL_MODELVIEW);
  camera.loadMatrix(GlCamera::MATRIX_MODELVIEW);

  glutPostRedisplay();
  
}

void MzGlutApp::drawString(int rasterx, int rastery, 
                           const std::string& str,
                           void* font,
                           int linespacing) {

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, width, height, 0);
  
  glRasterPos2f(rasterx, rastery);
  for (std::string::const_iterator i=str.begin(); i!=str.end(); ++i) {
    if (*i == '\n') {
      rastery += linespacing;
      glRasterPos2f(rasterx, rastery);
    } else {
      glutBitmapCharacter(font, *i);
    }
  }

  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

}
