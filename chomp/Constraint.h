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

#ifndef _CONSTRAINT_H_
#define _CONSTRAINT_H_

#include "chomputil.h"
#include <vector>

namespace chomp {

  class Constraint {
  public:
  
    virtual ~Constraint();

    virtual size_t numOutputs() =0;

    virtual void evaluateConstraints(const MatX& qt, 
                                     MatX& h, 
                                     MatX& H) =0;

  };

  class NullConstraint: public Constraint {
  public:

    virtual ~NullConstraint();
  
    virtual size_t numOutputs();

    virtual void evaluateConstraints(const MatX& qt, 
                                     MatX& h, 
                                     MatX& H);

  };

  class ConstantConstraint: public Constraint {
  public:

    size_t numCons;
    std::vector<size_t> index;
    std::vector<double> value;

    ConstantConstraint(const std::vector<size_t>& consIndex, 
                       const std::vector<double>& consValue);
  
    virtual ~ConstantConstraint();
  
    virtual size_t numOutputs();

    virtual void evaluateConstraints(const MatX& qt, 
                                     MatX& h, 
                                     MatX& H);

  };

}

#endif
