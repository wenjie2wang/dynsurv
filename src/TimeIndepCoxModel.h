//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef TIME_INDEP_COX_MODEL_H_
#define TIME_INDEP_COX_MODEL_H_

#include "IntRegPrior.h"
#include "IntRegPar.h"
#include "TimeIndepModel.h"
#include "CoxModel.h"

/*------------------------------------------------------------------*/
/* Declaration of time independent coef Cox model class             */
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior>
class TimeIndepCoxModel:
  public TimeIndepModel<Prior, TimeIndepCoxPar>,
  public CoxModel<Prior, TimeIndepCoxPar>
{
public :
  TimeIndepCoxModel<Prior>(const boost::shared_ptr<IntRegData>& pd)
    : IntRegModel<Prior, TimeIndepCoxPar>(pd),
      TimeIndepModel<Prior, TimeIndepCoxPar>(pd),
      CoxModel<Prior, TimeIndepCoxPar>(pd) {}

  /* Overload virtual functions */
  TimeIndepCoxPar initPar() const;

  void gibbsKernel(const Prior& prior,
                   TimeIndepCoxPar& par);

  ublas::vector<double> likeVec(const TimeIndepCoxPar& par) const;
};
}

/*------------------------------------------------------------------*/
/* Implementation of time indepedent coef Cox model                 */
/*------------------------------------------------------------------*/
namespace ir {
/* Initial model parameters */
template <typename Prior>
TimeIndepCoxPar TimeIndepCoxModel<Prior>::
initPar() const
{
  return TimeIndepCoxPar(this->initLambda(),
                         ublas::vector<double>(this->nBeta_, 0));
}

/* Gibbs kernel implementation */
template <typename Prior>
void TimeIndepCoxModel<Prior>::
gibbsKernel(const Prior& prior,
            TimeIndepCoxPar& par)
{
  ublas::matrix<double> tBetaMat(ublas::outer_prod(par.beta,
                                 ublas::vector<double>(this->K_, 1.0)));
  ublas::matrix<double> expXb(ublas::exp(ublas::prod(this->pd_->X(), tBetaMat)));

  /* Sample dNMat and YMat */
  ublas::matrix<int> dNMat(this->N_, this->K_, 0);
  ublas::matrix<double> YMat(this->N_, this->K_, 1.0);
  this->sampleMat(par.lambda, expXb, dNMat, YMat);

  /* Sample lambda */
  this->sampleLambda(expXb, dNMat, YMat, prior.base_prior, par.lambda);

  /* Sample beta */
  ublas::vector<double> omega(this->N_, 1.0);
  this->sampleBeta(par.lambda, dNMat, YMat, omega, prior.coef_prior, par.beta);
}

/* Likelihood function */
template <typename Prior>
ublas::vector<double> TimeIndepCoxModel<Prior>::
likeVec(const TimeIndepCoxPar& par) const
{
  ublas::matrix<double> betaMat(ublas::outer_prod(
                                  ublas::vector<double>(this->K_, 1.0), par.beta));
  return this->coxLikeVec(par.lambda, betaMat);
}
}

#endif // TIME_INDEP_COX_MODEL_H_
