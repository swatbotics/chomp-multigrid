set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

if(APPLE)
  include_directories(/opt/local/include) # MacPorts
  link_directories(/opt/local/lib)
  find_library(OPENGL_LIBRARY OpenGL)
else()
  find_library(OPENGL_LIBRARY GL)
  find_library(GLU_LIBRARY GLU)
  set(OPENGL_LIBRARY ${OPENGL_LIBRARY} ${GLU_LIBRARY})
endif()

include(FindPkgConfig)

find_library(PNG_LIBRARY png)
find_library(EXPAT_LIBRARY expat)
find_library(GLUT_LIBRARY glut)

pkg_search_module(EIGEN3 REQUIRED eigen3>=3)
pkg_search_module(CCD libccd ccd)
pkg_search_module(CAIRO cairo)

if (CAIRO_FOUND)
add_definitions(-DMZ_HAVE_CAIRO)
endif ()

include_directories(${EIGEN3_INCLUDE_DIRS})


set(CMAKE_C_FLAGS "-Wall -g")
set(CMAKE_CXX_FLAGS "-Wall -g")

if(APPLE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wsign-compare -Wno-deprecated")
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wsign-compare -Wno-deprecated")
endif()

set(CMAKE_C_FLAGS_DEBUG "-O")
set(CMAKE_CXX_FLAGS_DEBUG "-O")

set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

macro(add_gui_app name)
  if(APPLE)
    add_executable(${name} MACOSX_BUNDLE ${ARGN})
  else()
    add_executable(${name} ${ARGN})
  endif()
endmacro(add_gui_app)

include_directories(${PROJECT_SOURCE_DIR})
