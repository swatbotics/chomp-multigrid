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

#ifndef _GLCAMERA_H_
#define _GLCAMERA_H_

#include "Transform3.h"

class GlCamera {
public:

  enum CameraType {
    CAMERA_ORTHO,
    CAMERA_PERSPECTIVE
  };

  enum RotateType {
    ROTATE_2_AXIS,
    ROTATE_TRACKBALL
  };

  enum MatrixType {
    MATRIX_PROJECTION,
    MATRIX_MODELVIEW,
    MATRIX_PROJECTION_INVERSE,
    MATRIX_MODELVIEW_INVERSE,
  };

  enum MouseMode {
    MOUSE_NONE,
    MOUSE_PAN_XY,
    MOUSE_PAN_Z,
    MOUSE_ROTATE,
    MOUSE_ZOOM
  };

  enum MouseFlags {
    MOUSE_ALLOW_PAN    = 0x01,
    MOUSE_ALLOW_ROTATE = 0x02,
    MOUSE_ALLOW_ZOOM   = 0x04,
    MOUSE_ALLOW_ALL    = 0x07,
    MOUSE_ALLOW_NONE   = 0x00
  };

  //////////////////////////////////////////////////

  class PerspectiveInfo {
  public:
    float fovy;
    float near;
    float far;

    PerspectiveInfo();
    PerspectiveInfo(float fovy, float near, float far);

  };

  //////////////////////////////////////////////////

  class AimInfo {
  public:

    vec3f position;
    vec3f target;
    vec3f up;

    AimInfo();
    AimInfo(const vec3f& p, const vec3f& t, const vec3f& u);

    // translates in the current frame
    void translate(const vec3f& t);

    // rotates around the target position
    void rotate(const quatf& q);

    Transform3f xform() const;

  };

  //////////////////////////////////////////////////

  GlCamera();

  //////////////////////////////////////////////////

  void setViewport(int x, int y, int w, int h);
  
  void setOrtho(const vec3f& viewDims);

  void setPerspective(const PerspectiveInfo& pi);

  void setPerspective(float fovy=45.0f, 
                      float near=0.1,
                      float far=100.0);


  void setCameraType(CameraType t);
  void setRotateType(RotateType t);
  
  void aim(const AimInfo& aim);

  void aim(const vec3f& position,
           const vec3f& target,
           const vec3f& up);

  void setZoomRange(float minSize, float maxSize);
  void setXRotRange(float xmin, float xmax);
  void setYRotRange(float ymin, float ymax);
  void set2AxisAngles(float xr, float yr);
  void setMouseFlags(int flags);

  //////////////////////////////////////////////////

  CameraType             cameraType() const;
  RotateType             rotateType() const;
  const vec3f&           orthoViewDims() const;
  const PerspectiveInfo& perspectiveInfo() const;
  const AimInfo&         aimInfo() const;
  void                   get2AxisAngles(float& xr, float& yr) const;
  int                    mouseFlags() const;

  //////////////////////////////////////////////////

  void   loadMatrix(MatrixType m) const;
  void   multMatrix(MatrixType m) const;
  const mat4f& getMatrix(MatrixType m) const;
  
  const Transform3f& xform() const;

  const mat4f& projection() const;
  const mat4f& modelview() const;

  const mat4f& projectionInverse() const;
  const mat4f& modelviewInverse() const;

  const int* viewport() const;

  //////////////////////////////////////////////////

  void setHomePosition();
  void recallHomePosition();

  //////////////////////////////////////////////////

  void pan(float px, float py, float pz=0, bool autoscale=false);
  void pan(const vec3f& p, bool autoscale=false);
  void zoom(float factor);

  //////////////////////////////////////////////////

  void mousePress(int x, int y, MouseMode m);
  void mouseMove(int x, int y, MouseMode m);
  void mouseRelease(int x, int y, MouseMode m);
  void mouseReset();

  //////////////////////////////////////////////////

  vec3f unproject(int x, int y, float z=0.5) const;
  vec3f unproject(const vec3f& wincoords) const;

private:

  vec3f _hemiCoords(int x, int y) const;
  void _updateProjection();
  void _updateModelview();

  CameraType _cameraType;
  RotateType  _rotateType;

  int   _viewport[4];
  float _aspect; // derived from _viewport

  // rotating moves position but not target for trackball, just changes xrot/yrot for tilt
  // panning moves position AND target
  // zooming moves position only AND changes view dims

  class PoseInfo {
  public:
    AimInfo aim;
    float   xrot, yrot;
    vec3f   orthoViewDims;
    PoseInfo();
  };

  PerspectiveInfo _pinfo;

  PoseInfo _home;
  PoseInfo _current;
  
  mat4f       _projection; // derived from _current, _pinfo, _viewport
  float       _winsz; // derived from orthoViewDims or pinfo + aim
  float       _minsz;
  float       _maxsz;
  float       _minxrot;
  float       _maxxrot;
  float       _minyrot;
  float       _maxyrot;

  Transform3f _cxform; // derived from _current
  mat4f       _modelview; // derived from _current

  mat4f       _projection_inv;
  mat4f       _modelview_inv;
  

  int _mouseFlags;
  bool _mouseDown;
  int  _mouseMode;
  PoseInfo _mouseOrigPose;
  int _mouseX, _mouseY;
  float _mouseWinSz;
  vec3f _mouseHemi;
  
    
};

#endif
