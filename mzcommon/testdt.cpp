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

#include "DtGrid.h"
#include "MzGlutApp.h"
#include <assert.h>

class DtDemo: public MzGlutApp {
public: 

  const DtGridf& grid;
  float miniso, isoval;
  int pointSize;

  DtDemo(int argc, char** argv, const DtGridf& g): 
    MzGlutApp(argc, argv),
    grid(g),
    pointSize(3)

  {

    if (grid.minDist() < 0) {
      miniso = 0;
    } else {
      miniso = 0.5*grid.cellSize();
    }
    isoval = miniso;

    initWindowSize(640, 480);
    createWindow("DtGrid test");
    setupBasicLight();
    glClearColor(0,0,0,1);

    camera.aim(vec3f(3, 0, 0),
               vec3f(0, 0, 0),
               vec3f(0, 0, 1));

    camera.setPerspective();

    camera.setRotateType(GlCamera::ROTATE_2_AXIS);

    camera.setHomePosition();

  }

  void renderDisc(const vec3f& pos,
                  const vec3f& normal) const {
    glNormal3fv(normal.v);
    glVertex3fv(pos.v);
  }

  void renderGrid() const {

    float tol = 0.87*grid.cellSize();

    glPointSize(pointSize);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glBegin(GL_POINTS);
    
    for (size_t z=0; z<grid.nz(); ++z) {
      for (size_t y=0; y<grid.ny(); ++y) {
        for (size_t x=0; x<grid.nx(); ++x) {
          
          size_t i = grid.sub2ind(x,y,z);
          
          float val = grid[i];
          float dist = isoval - val;
          
          if (fabs(dist) < tol) {
            vec3f n = grid.normal(i);
            if (n.norm2() < 0.9) { continue; }
            vec3f pos = grid.cellCenter(i);
            pos += dist * n;
            renderDisc(pos, n);
          }
          
        }
      }
    }
    

    glEnd();
    glPopMatrix();

  }

  virtual void display() {

    MzGlutApp::display();

    glColor3ub(127,191,255);
    
    renderGrid();

    glutSwapBuffers();
    
  }

  virtual void keyboard(unsigned char key, int x, int y) {
    switch (key) {
    case 27:
      exit(0);
      break;
    case 'a':
      isoval += 0.5 * grid.cellSize();
      if (isoval > grid.maxDist()) { 
        isoval = grid.maxDist();
      }
      std::cerr << "isoval = " << isoval << "\n";
      glutPostRedisplay();
      break;
    case 'z':
      isoval -= 0.5 * grid.cellSize();
      if (isoval < grid.minDist()) { 
        isoval = grid.minDist();
      }
      std::cerr << "isoval = " << isoval << "\n";
      glutPostRedisplay();
      break;
    case 'x':
      isoval = miniso;
      std::cerr << "isoval = " << isoval << "\n";
      glutPostRedisplay();
      break;
    case '+':
    case '=':
      pointSize = std::min(pointSize+1, 10);
      glutPostRedisplay();
      break;
    case '-':
      pointSize = std::max(pointSize-1, 1);
      glutPostRedisplay();
      break;
    }
  }

  virtual void special(int key, int x, int y) {
    switch (key) {
    case GLUT_KEY_HOME:
      break;
    case GLUT_KEY_END:
      break;
    }
  }

};

int main(int argc, char** argv) {

  DtGridf g;

  float cellSize = 0.01;
  float border = 0.25;
  float amp = 0.125;
  float freq = 4*M_PI;
  float dx = 1;
  float dy = 1;
  float dz = 0.5;

  vec3f p0(-dx, -dy, -dz);
  vec3f p1( dx,  dy,  dz);

  g.resize(p0, p1, DtGridf::AXIS_Z, cellSize);

  HeightMapf& h = g.heightMap();
  assert(h.nx() == g.nx());
  assert(h.ny() == g.ny()); 
        

  for (size_t y=0; y<g.ny(); ++y) {
    for (size_t x=0; x<g.nx(); ++x) {
      vec2f p = h.cellCenter(x,y);
      if (p.x() > p0.x() + border && 
          p.x() < p1.x() - border && 
          p.y() > p0.y() + border &&
          p.y() < p1.y() - border) {
        h(x,y) = amp*cos(freq*p.x()) + amp*cos(freq*p.y());
      } else {
        h(x,y) = HeightMapf::INVALID_HEIGHT;
      }
    }
  }

  g.computeDistsFromHeightMap();

  DtDemo demo(argc, argv, g);
  demo.run();

  return 0;

}
