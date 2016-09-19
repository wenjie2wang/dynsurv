//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef DYNAMIC_GORH_MODEL_H_
#define DYNAMIC_GORH_MODEL_H_

#include "IntRegPrior.h"
#include "IntRegPar.h"
#include "DynamicModel.h"
#include "GORHModel.h"

/*------------------------------------------------------------------*/
/* Declaration of dynamic coef GORH model                           */
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior>
class DynamicGORHModel:
  public DynamicModel<Prior, DynamicGORHPar>,
  public GORHModel<Prior, DynamicGORHPar>
{
public :
  /* Class constructor */
  DynamicGORHModel<Prior>(const boost::shared_ptr<IntRegData>& pd)
    : IntRegModel<Prior, DynamicGORHPar>(pd),
      DynamicModel<Prior, DynamicGORHPar>(pd),
      GORHModel<Prior, DynamicGORHPar>(pd) {}

  /* Overload virtual functions */
  DynamicGORHPar initPar() const;

  void gibbsKernel(const Prior& prior,
                   DynamicGORHPar& par);

  ublas::vector<double> likeVec(const DynamicGORHPar& par) const;
};
}

/*------------------------------------------------------------------*/
/* Implementation of dynamic GORH model                             */
/*------------------------------------------------------------------*/
namespace ir {
/* Initial dynamic GORH model parameters */
template <typename Prior>
DynamicGORHPar DynamicGORHModel<Prior>::
initPar() const
{
  /* Initalize jump indicator matrix */
  ublas::matrix<int> jumpMat(this->K_, this->nBeta_);
  for (Size j = 0; j < this->nBeta_; ++j) {

    // place 1 jump for every 4 intervals
    for (Size k = 0; k < this->K_; ++k)
      jumpMat(k, j) = (k % 4 == 3) ? 1 : 0;

    jumpMat(this->K_ - 1, j) = 1;
  }

  return DynamicGORHPar(initLambda(), ublas::matrix<double>(this->K_, this->nBeta_, 0),
                        ublas::vector<double>(this->nBeta_, 1.0), jumpMat, 0.5);
}

/* Gibbs kernel implementation */
template <typename Prior>
void DynamicGORHModel<Prior>::
gibbsKernel(const Prior& prior,
            DynamicGORHPar& par)
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

  /* Loop through covariates */
  for (Size j = 0; j < this->nBeta_; ++j) {
    double u = rmath::unif_rand();
    double nJump = ublas::sum(ublas::column(par.jump, j));

    bool prop_flag = false;
    double prop_ratio = 1.0;
    DynamicGORHPar prop_par(par);

    /* Birth move */
    if (u < this->prob_(0) && nJump < this->K_) {
      prop_ratio = propBirth(j, par.beta, par.jump, prop_par.beta,
                             prop_par.jump);

      prop_flag = true;
    }

    /* Death move */
    if (u > this->prob_(0) && u < this->prob_(0) + this->prob_(1) && nJump > 1) {
      prop_ratio = propDeath(j, par.beta, par.jump, prop_par.beta,
                             prop_par.jump);

      prop_flag = true;
    }

    /* Decide whether to accept */
    if (prop_flag == true) {
      double prior_ratio = std::exp(logCoefPrior(ublas::column(prop_par.jump, j),
                                                 ublas::column(prop_par.beta, j),
                                                 prior.coef_prior) -
                                    logCoefPrior(ublas::column(par.jump, j),
                                                 ublas::column(par.beta, j),
                                                 prior.coef_prior));

      double like_ratio = std::exp(ublas::sum(ublas::log(gorhLikeVec(tLambda,
                                                                     prop_par.beta,
                                                                     prop_par.theta))) -
                                   ublas::sum(ublas::log(gorhLikeVec(tLambda,
                                                                     par.beta,
                                                                     par.theta))));

      double ratio = prop_ratio * prior_ratio * like_ratio;

      if (rmath::unif_rand() < std::min<double>(1.0, ratio))
        par = prop_par;
    }

    if (prop_flag == false)
      sampleBeta(j, dNMat, YMat, tLambda, par.jump, omega, prior.coef_prior,
                 par.beta, par.nu);
  }
}

/* Likelihood function */
template <typename Prior>
ublas::vector<double> DynamicGORHModel<Prior>::
likeVec(const DynamicGORHPar& par) const
{
  return gorhLikeVec(par.theta * par.lambda, par.beta, par.theta);
}
}

#endif // DYNAMIC_GORH_MODEL_H_
