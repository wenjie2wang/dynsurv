//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef TIME_VARYING_GORH_MODEL_H_
#define TIME_VARYING_GORH_MODEL_H_

#include "IntRegPrior.h"
#include "IntRegPar.h"
#include "TimeVaryingModel.h"
#include "GORHModel.h"

/*------------------------------------------------------------------*/
/* Declaration of time varying coef GORH model class                */
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior>
class TimeVaryingGORHModel:
  public TimeVaryingModel<Prior, TimeVaryingGORHPar>,
  public GORHModel<Prior, TimeVaryingGORHPar>
{
public :
  TimeVaryingGORHModel<Prior>(const boost::shared_ptr<IntRegData>& pd)
    : IntRegModel<Prior, TimeVaryingGORHPar>(pd),
      TimeVaryingModel<Prior, TimeVaryingGORHPar>(pd),
      GORHModel<Prior, TimeVaryingGORHPar>(pd) {}

  /* Overload virtual functions */
  TimeVaryingGORHPar initPar() const;

  void gibbsKernel(const Prior& prior,
                   TimeVaryingGORHPar& par);

  ublas::vector<double> likeVec(const TimeVaryingGORHPar& par) const;
};
}

/*------------------------------------------------------------------*/
/* Implementation of time varying coef GORH model                   */
/*------------------------------------------------------------------*/
namespace ir {
/* Initial time varying Bayes GORH model parameters */
template <typename Prior>
TimeVaryingGORHPar TimeVaryingGORHModel<Prior>::
initPar() const
{
  return TimeVaryingGORHPar(initLambda(), ublas::matrix<double>(this->K_, this->nBeta_, 0),
                            ublas::vector<double>(this->nBeta_, 1.0), 0.5);
}

/* Gibbs kernel implementation */
template <typename Prior>
void TimeVaryingGORHModel<Prior>::
gibbsKernel(const Prior& prior,
            TimeVaryingGORHPar& par)
{
  ublas::matrix<double> expXb(ublas::exp(ublas::prod(this->pd_->X(),
                                         ublas::trans(par.beta))));

  /* Sample dNMat and YMat */
  ublas::matrix<int> dNMat(this->N_, this->K_, 0);
  ublas::matrix<double> YMat(this->N_, this->K_, 1.0);
  ublas::vector<double> tLambda(par.theta * par.lambda);

  sampleMat(tLambda, expXb, par.theta, dNMat, YMat);

  /* Sample theta */
  sampleTheta(expXb, dNMat, YMat, tLambda, prior.base_prior, prior.theta_prior,
              par.theta);
  tLambda = par.theta * par.lambda;

  /* Sample omega */
  ublas::vector<double> omega(this->N_, 1.0);
  sampleOmega(expXb, dNMat, YMat, tLambda, par.theta, omega);

  /* Sample lambda */
  sampleTLambda(expXb, dNMat, YMat, omega, par.theta, prior.base_prior, tLambda);
  par.lambda = tLambda / par.theta;

  /* Sample beta */
  sampleBeta(tLambda, dNMat, YMat, omega, prior.coef_prior, par.beta, par.nu);
}

/* Likelihood function */
template <typename Prior>
ublas::vector<double> TimeVaryingGORHModel<Prior>::
likeVec(const TimeVaryingGORHPar& par) const
{
  return gorhLikeVec(par.theta * par.lambda, par.beta, par.theta);
}
}
#endif // TIME_VARYING_GORH_MODEL_H_
