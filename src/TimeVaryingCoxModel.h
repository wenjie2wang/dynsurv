//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef TIME_VARYING_COX_MODEL_H_
#define TIME_VARYING_COX_MODEL_H_

#include "IntRegPrior.h"
#include "IntRegPar.h"
#include "TimeVaryingModel.h"
#include "CoxModel.h"

/*------------------------------------------------------------------*/
/* Declaration of time varying coef Cox model class                 */
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior>
class TimeVaryingCoxModel:
  public TimeVaryingModel<Prior, TimeVaryingCoxPar>,
  public CoxModel<Prior, TimeVaryingCoxPar>
{
public :
  TimeVaryingCoxModel(const boost::shared_ptr<IntRegData>& pd)
    : IntRegModel<Prior, TimeVaryingCoxPar>(pd),
      TimeVaryingModel<Prior, TimeVaryingCoxPar>(pd),
      CoxModel<Prior, TimeVaryingCoxPar>(pd) {}

  /* Overload virtual functions */
  TimeVaryingCoxPar initPar() const;

  void gibbsKernel(const Prior& prior,
                   TimeVaryingCoxPar& par);

  ublas::vector<double> likeVec(const TimeVaryingCoxPar& par) const;
};
}

/*------------------------------------------------------------------*/
/* Implementation of time varying coef Cox model                    */
/*------------------------------------------------------------------*/
namespace ir {
/* Initial time varying Bayes Cox model parameters */
template <typename Prior>
TimeVaryingCoxPar TimeVaryingCoxModel<Prior>::
initPar() const
{
  return TimeVaryingCoxPar(this->initLambda(),
                           ublas::matrix<double>(this->K_, this->nBeta_, 0),
                           ublas::vector<double>(this->nBeta_, 1.0));
}

/* Gibbs kernel implementation */
template <typename Prior>
void TimeVaryingCoxModel<Prior>::
gibbsKernel(const Prior& prior,
            TimeVaryingCoxPar& par)
{
  ublas::matrix<double> expXb(ublas::exp(ublas::prod(this->pd_->X(),
                                         ublas::trans(par.beta))));

  /* Sample dNMat and YMat */
  ublas::matrix<int> dNMat(this->N_, this->K_, 0);
  ublas::matrix<double> YMat(this->N_, this->K_, 1.0);
  this->sampleMat(par.lambda, expXb, dNMat, YMat);

  /* Sample lambda */
  this->sampleLambda(expXb, dNMat, YMat, prior.base_prior, par.lambda);

  /* Sample beta */
  ublas::vector<double> omega(this->N_, 1.0);
  this->sampleBeta(par.lambda, dNMat, YMat, omega, prior.coef_prior, par.beta, par.nu);
}

/* Likelihood function */
template <typename Prior>
ublas::vector<double> TimeVaryingCoxModel<Prior>::
likeVec(const TimeVaryingCoxPar& par) const
{
  return this->coxLikeVec(par.lambda, par.beta);
}
}
#endif // TIME_VARYING_COX_MODEL_H_
