# CMakeLists.txt for hubomz/mzcommon

set(mzcommon_srcs
  mersenne.c
  gauss.c
  GlCamera.cpp
  GlFrameGrabber.cpp
  GlTexture2D.cpp
  MzGlutApp.cpp
  TinyDom.cpp
  TriMesh3.cpp
  HeightMap.cpp
  DtGrid.cpp
  Geom2.cpp
  glstuff.cpp
  CCDWrapper.cpp)

if (GLEW_FOUND)
  set(mzcommon_srcs ${mzcommon_srcs} GlFramebufferObject.cpp)
endif()

add_library( mzcommon ${mzcommon_srcs} )
target_link_libraries(mzcommon ${OPENGL_LIBRARY} ${GLUT_LIBRARY} ${EXPAT_LIBRARY} ${PNG_LIBRARY} ${CCD_LDFLAGS} )

if (GLEW_FOUND)
  target_link_libraries(mzcommon ${GLEW_LIBRARIES})
endif()

add_executable(testwrl testwrl.cpp)
target_link_libraries(testwrl mzcommon ${OPENGL_LIBRARY})

add_gui_app(testdt testdt.cpp)
target_link_libraries(testdt mzcommon ${OPENGL_LIBRARY} ${GLUT_LIBRARY})

add_executable(testgeom2 testgeom2.cpp)

if (${CAIRO_FOUND})
target_link_libraries(testgeom2 mzcommon ${CAIRO_LIBRARIES})
else()
target_link_libraries(testgeom2 mzcommon)
endif()

add_gui_app( ccddemo ccddemo.cpp )
target_link_libraries( ccddemo mzcommon 
  ${CCD_LDFLAGS} ${GLUT_LIBRARY} ${OPENGL_LIBRARY} expat )
