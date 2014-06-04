chomp-multigrid
===============

Demonstration of various aspects of CHOMP, with & without constraints.

See the paper "Multigrid CHOMP with Local Smoothing" by He, Martin, &
Zucker, 2013 for details. (In the repository under the `paper`
directory).

Please contact Matt Zucker <mzucker1@swarthmore.edu> with questions,
comments, or rants.

Building
========

To build the software, you will need these libraries:
 
  - GLUT
  - eigen3
  - expat
  - cairo 
  - ccd (optional for collision detection, use ros-groovy-libccd)

You must install cmake to build the software as well.  To build:

    cd /path/to/chomp-multigrid
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make

Programs
========
 
From the build directory, run

    ../chomp/map2d_tests.sh
    
or

    ../chomp/circle_tests.sh

