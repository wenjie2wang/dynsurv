//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef DYNAMIC_COX_MODEL_V2_H_
#define DYNAMIC_COX_MODEL_V2_H_

#include "DynamicCoxModel.h"

/*------------------------------------------------------------------*/
/* Declaration of dynamic Cox model (treat log(lambda) as intercept)*/
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior>
class DynamicCoxModel_v2:
  public DynamicCoxModel<Prior>
{
public :
  DynamicCoxModel_v2<Prior>(const boost::shared_ptr<IntRegData>& pd)
    : IntRegModel<Prior, DynamicCoxPar>(pd),
      DynamicCoxModel<Prior>(pd) {}

  /* Overload virtual functions */
  DynamicCoxPar initPar() const;
};
}

/*------------------------------------------------------------------*/
/* Implementation of dynamic Cox model (treat log(lambda) as intercept)*/
/*------------------------------------------------------------------*/
namespace ir {
/* Initial dynamic Cox model parameters */
template <typename Prior>
DynamicCoxPar DynamicCoxModel_v2<Prior>::
initPar() const
{
  /* Initalize jump indicator matrix */
  ublas::matrix<int> jumpMat(this->K_, this->nBeta_, 1);
  // start at 1
  for (Size j = 1; j < this->nBeta_; ++j) {

    // place 1 jump for every 4 intervals
    for (Size k = 0; k < this->K_; ++k)
      jumpMat(k, j) = (k % 4 == 3) ? 1 : 0;

    jumpMat(this->K_ - 1, j) = 1;
  }

  // Workaround: the first column of betaMat is log(lambda)
  ublas::matrix<double> betaMat(this->K_, this->nBeta_, 0);
  ublas::column(betaMat, 0) = ublas::log(this->initLambda() +
                                         ublas::vector<double>(this->K_, 0.001));

  return DynamicCoxPar(ublas::vector<double>(this->K_, 1.0),
                       betaMat, ublas::vector<double>(this->nBeta_, 1.0), jumpMat);
}
}

#endif // DYNAMIC_COX_MODEL_V2_H_
