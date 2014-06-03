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

#include "ConstraintFactory.h"
#include "Constraint.h"

namespace chomp {


  void ConstraintFactory::getAll(size_t total, std::vector<Constraint*>& constraints) {
    constraints.clear(); 
    for (size_t i=0; i<total; i++){
      constraints.push_back(getConstraint(i, total));
    }
    assert(constraints.size() == total);
  }

  void ConstraintFactory::evaluate(const std::vector<Constraint*>& constraints, 
                                   const MatX& xi, 
                                   MatX& h_tot, 
                                   MatX& H_tot, 
                                   int step)
  {

    size_t DoF = xi.cols();

    assert(size_t(xi.rows()) == constraints.size());

    size_t numCons = 0;

    size_t timesteps = 0;

    for (size_t t=0; t<constraints.size(); t+=step) {
      Constraint* c = constraints[t];
      if (c) {
        numCons += c->numOutputs();
      }
      ++timesteps;
    }

    if (size_t(h_tot.rows()) != numCons || size_t(h_tot.cols()) != 1) {
      h_tot.resize(numCons,1);
    }

    if (size_t(H_tot.rows()) != numCons || size_t(H_tot.cols()) != DoF*timesteps) {
      H_tot.resize(numCons,DoF*timesteps);
    }

    H_tot.setZero(); 

    MatX h, H; // individual constraints at time t
    size_t row = 0; // starting row for the current bundle o constraints

    for (size_t t=0, i=0; t<constraints.size(); t+=step, ++i) {

      Constraint* c = constraints[t];

      //get individual h and H from qt
      if (!c || c->numOutputs()==0){ continue; }
      c->evaluateConstraints(xi.row(t), h, H);

      if (h.rows() == 0) { 
        assert(H.rows() == 0);
        continue;
      }
    
      //stick all the h and H into h_tot and H_tot
      assert((size_t)H.cols()==DoF);
      assert(h.cols() == 1);
      assert((size_t)H.rows()<=constraints[t]->numOutputs());
      assert(H.rows() == h.rows());

      for (size_t r = 0; r<(size_t)h.rows(); r++){
        h_tot(row) = h(r);
        for (size_t j = 0; j<DoF; j++){
          assert( row < size_t(H_tot.rows()) );
          assert( j*timesteps+i < size_t(H_tot.cols()) );
          assert( j < size_t(H.cols()) );
          H_tot(row, j*timesteps+i) = H(r,j);
        }
        row++;
      }

    }  

    assert( row <= (size_t)numCons );
    if ( row < numCons ) {
      h_tot = MatX(h_tot.block(0, 0, row, 1));
      H_tot = MatX(H_tot.block(0, 0, row, DoF*timesteps));
    }

    assert(h_tot.rows() == H_tot.rows());
    assert(h_tot.cols() == 1);
    assert(size_t(H_tot.cols()) == DoF*timesteps);

  }

}
