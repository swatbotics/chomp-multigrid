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

#include "GlCamera.h"    
#include "glstuff.h"
#include <iostream>

GlCamera::PerspectiveInfo::PerspectiveInfo():
  fovy(45), near(0.001f), far(100.0f) {}

GlCamera::PerspectiveInfo::PerspectiveInfo(float a, float b, float c):
  fovy(a), near(b), far(c) {}

GlCamera::AimInfo::AimInfo():
  position(0,0,0),
  target(0,0,-1),
  up(0,1,0) {}


GlCamera::AimInfo::AimInfo(const vec3f& p, const vec3f& t, const vec3f& u):
  position(p), target(t), up(u) {}

vec3f GlCamera::unproject(int x, int y, float z) const {
  return unproject(vec3f(x, _viewport[3]-y-1, z));
}

vec3f GlCamera::unproject(const vec3f& wincoords) const {

  GLint vp[4];

  for (int i=0; i<4; ++i) { vp[i] = _viewport[i]; }

  GLdouble mmat[16], pmat[16];
  
  int offs = 0;
  for (int col=0; col<4; ++col) {
    for (int row=0; row<4; ++row) {
      mmat[offs] = _modelview(row,col);
      pmat[offs] = _projection(row,col);
      ++offs;
    }
  }

  GLdouble x, y, z;

  gluUnProject(wincoords[0], wincoords[1], wincoords[2],
               mmat, pmat, vp, &x, &y, &z);

  return vec3f(x,y,z);
    
}

void GlCamera::AimInfo::translate(const vec3f& v) {

  Transform3f t = xform();
  vec3f vt = t.rotInv() * v;
  position += vt;
  target += vt;

}

void GlCamera::AimInfo::rotate(const quatf& q) {

  float d = (position-target).norm();

  Transform3f told = xform();

  quatf qnew = q * told.rotation();
  mat3f Rnew = qnew.toMat3();
    
  vec3f fwd(Rnew(2,0), Rnew(2,1), Rnew(2,2));
  position = target + fwd * d;

  up = vec3f(Rnew(1,0), Rnew(1,1), Rnew(1,2));

}

Transform3f GlCamera::AimInfo::xform() const {

  const vec3f& p = position;

  quatf q = quatf::fromTwoVectors( p-target, 2, up, 1 );
  
  return Transform3f(q, p).inverse();

}

GlCamera::PoseInfo::PoseInfo(): 
  xrot(0), yrot(0), orthoViewDims(1, 1, 100) {}

void GlCamera::setZoomRange(float z0, float z1) {
  _minsz = z0;
  _maxsz = z1;
}

void GlCamera::setXRotRange(float x0, float x1) {
  _minxrot = x0;
  _maxxrot = x1;
}

void GlCamera::setYRotRange(float y0, float y1) {
  _minyrot = y0;
  _maxyrot = y1;
}

void GlCamera::zoom(float factor) {

  float oldsz = _winsz;
  float newsz = oldsz * factor;

  if (_minsz < _maxsz) {
    if (newsz > _maxsz) {
      newsz = _maxsz;
    } else if (newsz < _minsz) {
      newsz = _minsz;
    }
  }
  
  _current.orthoViewDims *= newsz / _current.orthoViewDims.y();
  _updateProjection();
  
  vec3f diff = _current.aim.position - _current.aim.target;
  diff /= diff.norm();
  float f = tan(0.5 * _pinfo.fovy * M_PI / 180.0);
  float d = 0.5 * newsz / f;
  _current.aim.position = _current.aim.target + diff * d;
  _updateModelview();

}

void GlCamera::setRotateType(RotateType t) {
  // TODO: handle switching from ROTATE_2_AXIS
  _rotateType = t;
  _updateModelview();
}

void GlCamera::setCameraType(CameraType t) {
  _cameraType = t;
  _updateProjection();
  _updateModelview();
}

GlCamera::CameraType GlCamera::cameraType() const {
  return _cameraType;
}

GlCamera::RotateType GlCamera::rotateType() const {
  return _rotateType;
}

GlCamera::GlCamera():
  _cameraType(CAMERA_PERSPECTIVE),
  _rotateType(ROTATE_TRACKBALL),
  _minsz(1),
  _maxsz(0),
  _minxrot(1),
  _maxxrot(0),
  _minyrot(1),
  _maxyrot(0),
  _mouseFlags(MOUSE_ALLOW_ALL),
  _mouseDown(false),
  _mouseMode(MOUSE_NONE)
{
  setViewport(0,0,100,100);
  _updateModelview();
}

void GlCamera::setViewport(int x, int y, int w, int h) {
  _viewport[0] = x;
  _viewport[1] = y;
  _viewport[2] = w;
  _viewport[3] = h;
  _updateProjection();
}

void GlCamera::setOrtho(const vec3f& viewDims) {
  _cameraType = CAMERA_ORTHO;
  _current.orthoViewDims = viewDims;
  // TODO: set camera distance to match current
  _updateProjection();
}

void GlCamera::setPerspective(float fovy, float near, float far) {
  setPerspective(PerspectiveInfo(fovy, near, far));
}

void GlCamera::setPerspective(const PerspectiveInfo& pi) {
  _cameraType = CAMERA_PERSPECTIVE;
  _pinfo = pi;
  // TODO: set ortho view size to match current
  _updateProjection();
}

void GlCamera::aim(const vec3f& position,
                   const vec3f& target,
                   const vec3f& up) {

  aim(AimInfo(position, target, up));

}

void GlCamera::aim(const AimInfo& aim) {

  _current.aim = aim;
  _updateModelview();

}

const vec3f& GlCamera::orthoViewDims() const {

  return _current.orthoViewDims;

}

const GlCamera::PerspectiveInfo& GlCamera::perspectiveInfo() const {
  return _pinfo;
}

const GlCamera::AimInfo& GlCamera::aimInfo() const {
  return _current.aim;
}

void GlCamera::get2AxisAngles(float& xr, float& yr) const {
  xr = _current.xrot;
  yr = _current.yrot;
}

void GlCamera::set2AxisAngles(float xr, float yr) {
  _current.xrot = xr;
  _current.yrot = yr;
  if (_minyrot < _maxyrot) {
    _current.yrot = std::min(std::max(_current.yrot, _minyrot), _maxyrot);
  }
  if (_minxrot < _maxxrot) {
    _current.xrot = std::min(std::max(_current.xrot, _minxrot), _maxxrot);
  }
  _updateModelview();
}

void GlCamera::loadMatrix(MatrixType m) const {
  glstuff::load_mat4(getMatrix(m));
}

void GlCamera::multMatrix(MatrixType m) const {
  glstuff::mult_mat4(getMatrix(m));
}

const mat4f& GlCamera::getMatrix(MatrixType m) const {
  switch (m) {
  case MATRIX_PROJECTION_INVERSE:
    return _projection_inv;
  case MATRIX_MODELVIEW_INVERSE:
    return _modelview_inv;
  case MATRIX_PROJECTION:
    return _projection;
  case MATRIX_MODELVIEW:
  default:
    return _modelview;
  }
}

const mat4f& GlCamera::projection() const { return _projection; }

const mat4f& GlCamera::modelview() const { return _modelview; }

const mat4f& GlCamera::projectionInverse() const { return _projection_inv; }

const mat4f& GlCamera::modelviewInverse() const { return _modelview_inv; }

const int* GlCamera::viewport() const { return _viewport; }

void GlCamera::mousePress(int x, int y, MouseMode m) {
  bool allowed = true;
  switch (m) {
  case MOUSE_PAN_XY: 
  case MOUSE_PAN_Z:
    allowed = (_mouseFlags & MOUSE_ALLOW_PAN);
    break;
  case MOUSE_ROTATE:
    allowed = (_mouseFlags & MOUSE_ALLOW_ROTATE);
    break;
  case MOUSE_ZOOM:
    allowed = (_mouseFlags & MOUSE_ALLOW_ZOOM);
    break;
  default:
    allowed = false;
    break;
  }
  if (_mouseDown || !allowed) {
    mouseReset();
    return;
  }
  _mouseDown = true;
  _mouseMode = m;
  _mouseOrigPose = _current;
  _mouseX = x;
  _mouseY = y;
  _mouseWinSz = _winsz;
  _mouseHemi = _hemiCoords(x, y);
}

vec3f GlCamera::_hemiCoords(int x, int y) const {

  float xc = x - _viewport[0] - 0.5f * _viewport[2];
  float yc = y - _viewport[1] - 0.5f * _viewport[3];

  float scl = 0.5f * std::min(_viewport[2], _viewport[3]);

  float ux = xc / scl;
  float uy = -yc / scl;

  float d2 = ux * ux + uy * uy;
  float uz = 0;

  if (d2 < 1) {
    uz = sqrt(1 - d2);
  } else {
    float d = sqrt(d2);
    ux /= d;
    uy /= d;
  }

  vec3f rval(ux, uy, uz);

  return rval;
  
}

void GlCamera::mouseMove(int x, int y, MouseMode m) {

  if (!_mouseDown || m != _mouseMode) {
    mouseReset();
    return;
  } 

  int dx = x - _mouseX;
  int dy = y - _mouseY;

  switch (m) {
  case MOUSE_ZOOM: {
    float factor = pow(1.05, dy);
    _winsz = _mouseWinSz;
    _current = _mouseOrigPose;
    zoom(factor);
    break;
  }
  case MOUSE_ROTATE:
    if (_rotateType == ROTATE_2_AXIS) {
      _current = _mouseOrigPose;
      set2AxisAngles(_current.xrot+dy, _current.yrot+dx);
    } else { // trackball
      vec3f h1 = _hemiCoords(x, y);
      vec3f axis = vec3f::cross(_mouseHemi, h1);
      float st = axis.norm();
      float ct = vec3f::dot(_mouseHemi, h1);
      if (st  > 1e-5) {
        axis /= st;
        float angle = acos(std::min(1.0f, std::max(ct, -1.0f)));
        _current = _mouseOrigPose;
        _current.aim.rotate(quatf::fromAxisAngle(axis, angle));
        _updateModelview();
      }
    }
    break;
  case MOUSE_PAN_XY:
  case MOUSE_PAN_Z:
  default:  {
    _current = _mouseOrigPose;
    int psz = std::min(_viewport[2], _viewport[3]);
    float scl = _winsz / psz;
    if (m == MOUSE_PAN_Z) {
      pan(0, 0, 4*dy*scl);
    } else {
      pan(-dx*scl, dy*scl, 0);
    }
    break;
  }
  }

}

void GlCamera::pan(float px, float py, float pz, bool autoscale) {
  pan(vec3f(px,py,pz),autoscale);
}

void GlCamera::pan(const vec3f& p, bool autoscale) {
  float s = 1;
  if (autoscale) {
    int psz = std::min(_viewport[2], _viewport[3]);
    s = _winsz / psz;
  }
  vec3f offs = _cxform.rotInv() * (s * p);
  _current.aim.position += offs;
  _current.aim.target += offs;
  _updateModelview();
}



void GlCamera::mouseRelease(int x, int y, MouseMode m) {
  if (!_mouseDown || m != _mouseMode) {
    mouseReset();
    return;
  } 
  mouseMove(x, y, m);
  _mouseDown = false;
}

void GlCamera::mouseReset() {
  if (_mouseDown) {
    _mouseDown = false;
    _current = _mouseOrigPose;
    _updateModelview();
  }
  _mouseMode = MOUSE_NONE;
}

void GlCamera::setHomePosition() {
  _home = _current;
}

void GlCamera::recallHomePosition() {
  _current = _home;
  _updateModelview();
  _updateProjection();
}

void GlCamera::_updateProjection() {

  int w = _viewport[2];
  int h = _viewport[3];
  if (!h) { ++h; }

  _aspect = float(w)/float(h);

  if (_cameraType == CAMERA_PERSPECTIVE) {

    float f = 1.0f / tan(0.5 * _pinfo.fovy * M_PI / 180.0);
    
    float zn = _pinfo.near;
    float zf = _pinfo.far;
    
    float a = (zf + zn) / (zn - zf);
    float b = (2 * zf * zn) / (zn - zf);
    
    _projection = mat4f::identity();
    _projection(0,0) = f / _aspect;
    _projection(1,1) = f;
    _projection(2,2) = a;
    _projection(2,3) = b;
    _projection(3,2) = -1;
    _projection(3,3) = 0;

    _projection_inv = _projection;

    _projection_inv(0,0) = _aspect / f;
    _projection_inv(1,1) = 1 / f;
    _projection_inv(2,2) = 0;
    _projection_inv(2,3) = -1;
    _projection_inv(3,2) = 1/b;
    _projection_inv(3,3) = a/b;

  } else {

    const vec3f& vdims = _current.orthoViewDims;
    
    float rw = vdims[0];
    float rh = vdims[1];
    float rd = vdims[2];

    float ra = rw / rh;

    float va = _aspect;

    float vw, vh;

    if (ra > va) {
      vw = rw;
      vh = rh*ra/va;
    } else {
      vw = rw*va/ra;
      vh = rh;
    }
    
    _projection = mat4f::identity();

    _projection(0,0) = 2 / vw;
    _projection(1,1) = 2 / vh;
    _projection(2,2) = -2 / rd;
    _projection(2,3) = -0.5;

    _projection_inv = _projection;

    _projection(0,0) = 0.5 * vw;
    _projection(1,1) = 0.5 * vh;
    _projection(2,2) = -0.5 * rd;
    _projection(2,3) = 0.25 * rd;

  }
  
  
}

void GlCamera::_updateModelview() {

  if (_cameraType == CAMERA_PERSPECTIVE) {
    float d = (_current.aim.position - _current.aim.target).norm();
    float f = tan(0.5 * _pinfo.fovy * M_PI / 180.0);
    _winsz = 2 * d * f;
  } else {
    _winsz = _current.orthoViewDims.y();
  }
  
  if (_rotateType == ROTATE_2_AXIS) {
    AimInfo tmp = _current.aim;
    quatf rx = quatf::fromAxisAngle(vec3f(1,0,0), _current.xrot * M_PI / 180);
    quatf ry = quatf::fromAxisAngle(vec3f(0,1,0), _current.yrot * M_PI / 180);
    tmp.rotate(rx * ry);
    _cxform = tmp.xform();
  } else {
    _cxform = _current.aim.xform();
  };

  _modelview = _cxform.matrix();
  _modelview_inv = _cxform.inverse().matrix();

}

void GlCamera::setMouseFlags(int flags) {
  _mouseFlags = flags;
}

int GlCamera::mouseFlags() const {
  return _mouseFlags;
}

const Transform3f& GlCamera::xform() const {
  return _cxform;
}
