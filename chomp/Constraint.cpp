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

#include "Constraint.h"

namespace chomp {

  Constraint::~Constraint() {}

  NullConstraint::~NullConstraint() {}

  size_t NullConstraint::numOutputs() {
    return 0;
  }

  void NullConstraint::evaluateConstraints(const MatX& qt, 
                                           MatX& h, 
                                           MatX& H){
  }

  ConstantConstraint::~ConstantConstraint() {}



  ConstantConstraint::ConstantConstraint(const std::vector<size_t>& consIndex, 
                                         const std::vector<double>& consValue){
    // Check to see if the input is valid
    assert(consIndex.size()==consValue.size()); 
    // Make sure at most 1 constraint per DoF
    index = consIndex;
    std::sort(index.begin(), index.end());
    std::vector<size_t>::iterator it;
    it = std::adjacent_find(index.begin(), index.end());
    assert(it==index.end());
    // Assign values
    index = consIndex;
    value = consValue;
    numCons = consIndex.size();
  }

  size_t ConstantConstraint::numOutputs() {
    return numCons;
  }

  void ConstantConstraint::evaluateConstraints(const MatX& qt, 
                                               MatX& h, 
                                               MatX& H)

  {

    assert(qt.rows() == 1 || qt.cols() == 1);

    size_t DoF = std::max(qt.rows(), qt.cols());

    if (size_t(h.rows()) != numCons || h.cols() != 1) {
      h.resize(numCons, 1);
    }
    
    if (size_t(H.rows()) != numCons || size_t(H.cols()) != DoF) {
      H.resize(numCons, DoF);
    }

    H.setZero();

    for(size_t i=0; i<numCons; i++){
      assert(index[i] < DoF);
      h(i) = qt(index[i])-value[i];
      H(i,index[i]) = 1;
    }

  }

}
