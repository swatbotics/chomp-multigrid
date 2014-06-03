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

#ifndef _CHOMPUTIL_H_
#define _CHOMPUTIL_H_

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>

namespace chomp {

  typedef Eigen::MatrixXd MatX;

  template <class Derived>
  static inline double mydot(const Eigen::MatrixBase<Derived>& a,
                             const Eigen::MatrixBase<Derived>& b) {

    return a.cwiseProduct(b).sum();

  }

  template <class Derived1, class Derived2, class Derived3>
  inline void diagMul(const Eigen::MatrixBase<Derived1>& coeffs, // e.g. [1, -4, 6]
                      const Eigen::MatrixBase<Derived2>& x,
                      const Eigen::MatrixBase<Derived3>& Ax_const) {

    assert( Ax_const.rows() == x.rows() && Ax_const.cols() == x.cols() );

    Eigen::MatrixBase<Derived3>& Ax = 
      const_cast<Eigen::MatrixBase<Derived3>&>(Ax_const);

    int n = x.rows();
    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;

    for (int i=0; i<x.rows(); ++i) {

      int j0 = std::max(i-o, int(0));
      int j1 = std::min(i+nc, n);

      Ax.row(i) = x.row(j0) * coeffs(j0-i+o);

      for (int j=j0+1; j<=i; ++j) {
        Ax.row(i) += x.row(j) * coeffs(j-i+o);
      }
    
      for (int j=i+1; j<j1; ++j) {
        Ax.row(i) += x.row(j) * coeffs(i-j+o);
      }

    }

  }

  //////////////////////////////////////////////////////////////////////

  template <class Derived1, class Derived2>
  void skylineChol(int n,
                   const Eigen::MatrixBase<Derived1>& coeffs, // e.g. [-1,2] or [1,-4,6]
                   Eigen::PlainObjectBase<Derived2>& L) {

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;

    L.resize(n, nc);
  
    for (int j=0; j<n; ++j) {
      int i1 = std::min(j+nc, n);
      for (int i=j; i<i1; ++i) {
        double sum = 0;
        int k0 = std::max(0,i-o);
        for (int k=k0; k<j; ++k) {
          sum += L(i,k-i+o) * L(j,k-j+o); // k < j < i
        }
        if (i == j) {
          L(j,o) = sqrt(coeffs(o) - sum);
        } else {
          L(i,j-i+o) = (coeffs(j-i+o) - sum) / L(j,o);
        }
      }
    }

  }

  //////////////////////////////////////////////////////////////////////

  template <class Derived1, class Derived2>
  inline void skylineCholSolve(const Eigen::MatrixBase<Derived1>& L,
                               const Eigen::MatrixBase<Derived2>& x_const) {

    int n = L.rows();
    assert(x_const.rows() == n);

    Eigen::MatrixBase<Derived2>& x = 
      const_cast<Eigen::MatrixBase<Derived2>&>(x_const);

    int nc = L.cols();
    int o = nc-1;

    for (int i=0; i<n; ++i) {
      int j0 = std::max(0, i-o);
      for (int j=j0; j<i; ++j) {
        x.row(i) -= L(i,j-i+o)*x.row(j); // here j < i so col < row
      }
      x.row(i) /= L(i,o);
    }

    for (int i=n-1; i>=0; --i) {
      int j1 = std::min(i+nc, n);
      for (int j=i+1; j<j1; ++j) {
        x.row(i) -= L(j,i-j+o) * x.row(j); // here j > i so col < row
      }
      x.row(i) /= L(i,o);
    }

  }

  //////////////////////////////////////////////////////////////////////

  template <class Derived1, class Derived2>
  inline void skylineCholSolveMulti(const Eigen::MatrixBase<Derived1>& L, 
                                    const Eigen::MatrixBase<Derived2>& xx_const) {


    int n = L.rows();
    int m = xx_const.rows() / n;
    assert(xx_const.rows() == m*n);

    Eigen::MatrixBase<Derived2>& xx = 
      const_cast<Eigen::MatrixBase<Derived2>&>(xx_const);

    for (int i=0; i<m; ++i) {
      skylineCholSolve(L, xx.block(n*i, 0, n, xx.cols()));
    }
    
  }


  template <class Derived>
  inline MatX getPos(const Eigen::MatrixBase<Derived>& x,
                     double h) {

    double fac = 1;
    double hn = 1;

    MatX rval = x.row(0);
    for (int i=1; i<x.rows(); ++i) {
      fac *= i;
      hn *= h;
      rval += (hn / fac) * x.row(i);
    }

    return rval;
    
  }

  template <class Derived1, class Derived2, class Derived3>
  inline double createBMatrix(int n, 
                              const Eigen::MatrixBase<Derived1>& coeffs,
                              const Eigen::MatrixBase<Derived2>& x0,
                              const Eigen::MatrixBase<Derived2>& x1,
                              const Eigen::MatrixBase<Derived3>& b_const,
                              double dt) {

    int nc = coeffs.rows() * coeffs.cols();
    int o = nc-1;
    
    assert(b_const.rows() == n);
    assert(b_const.cols() == x0.cols());
    assert(b_const.cols() == x1.cols());

    Eigen::MatrixBase<Derived3>& b = 
      const_cast<Eigen::MatrixBase<Derived3>&>(b_const);

    b.setZero();

    double c = 0;

    for (int i0=0; i0<o; ++i0) {
      int nn = nc-i0-1;
      int i1 = n-i0-1;
      for (int j=0; j<nn; ++j) {
        // so when j=o, we want t0=0, and when j<o we want t0<0
        //    when j=o, we want t1=0, and when j<o we want t1>0
        int t0 = j-o;
        int t1 = -t0;
        b.row(i0) += coeffs(j)*getPos(x0, t0*dt);
        b.row(i1) += coeffs(j)*getPos(x1, t1*dt);
      }
      c += mydot(b.row(i0), b.row(i0));
      if (i0 != i1) { c += mydot(b.row(i1), b.row(i1)); }
    }

    return 0.5*c;

  }

}


  //////////////////////////////////////////////////////////////////////


  /*

//////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////

inline fr::Mat safe_reshape(const fr::Mat& mat, int newrows, int newcols) {
const cv::Mat& orig = mat;
fr::Mat rval = orig.reshape(1, newrows);
assert(rval.rows == newrows && rval.cols == newcols);
return rval;
}

  */

#endif
