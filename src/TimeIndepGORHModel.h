//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef TIME_INDEP_GORH_MODEL_H_
#define TIME_INDEP_GORH_MODEL_H_

#include "IntRegPrior.h"
#include "IntRegPar.h"
#include "TimeIndepModel.h"
#include "GORHModel.h"

/*------------------------------------------------------------------*/
/* Declaration of time independent coef GORH model class            */
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior>
class TimeIndepGORHModel:
  public TimeIndepModel<Prior, TimeIndepGORHPar>,
  public GORHModel<Prior, TimeIndepGORHPar>
{
public :
  TimeIndepGORHModel<Prior>(const boost::shared_ptr<IntRegData>& pd)
    : IntRegModel<Prior, TimeIndepGORHPar>(pd),
      TimeIndepModel<Prior, TimeIndepGORHPar>(pd),
      GORHModel<Prior, TimeIndepGORHPar>(pd) {}

  /* Overload load virtual functions */
  TimeIndepGORHPar initPar() const;

  void gibbsKernel(const Prior& prior,
                   TimeIndepGORHPar& par);

  ublas::vector<double> likeVec(const TimeIndepGORHPar& par) const;
};
}

/*------------------------------------------------------------------*/
/* Implementation of time indepedent coef GORH model                */
/*------------------------------------------------------------------*/
namespace ir {
/* Initial model parameters */
template <typename Prior>
TimeIndepGORHPar TimeIndepGORHModel<Prior>::
initPar() const
{
  return TimeIndepGORHPar(initLambda(), ublas::vector<double>(this->nBeta_, 0), 0.5);
}

/* Gibbs kernel implementation */
template <typename Prior>
void TimeIndepGORHModel<Prior>::
gibbsKernel(const Prior& prior,
            TimeIndepGORHPar& par)
{
  ublas::matrix<double> tBetaMat(ublas::outer_prod(par.beta,
                                 ublas::vector<double>(this->K_, 1.0)));
  ublas::matrix<double> expXb(ublas::exp(ublas::prod(this->pd_->X(), tBetaMat)));

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
  sampleBeta(tLambda, dNMat, YMat, omega, prior.coef_prior, par.beta);
}

/* Likelihood function */
template <typename Prior>
ublas::vector<double> TimeIndepGORHModel<Prior>::
likeVec(const TimeIndepGORHPar& par) const
{
  ublas::matrix<double> betaMat(ublas::outer_prod(
                                  ublas::vector<double>(this->K_, 1.0), par.beta));
  return gorhLikeVec(par.theta * par.lambda, betaMat, par.theta);
}
}
#endif // TIME_INDEP_GORH_MODEL_H_
