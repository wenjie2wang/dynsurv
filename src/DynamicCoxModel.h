//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef DYNAMIC_COX_MODEL_H_
#define DYNAMIC_COX_MODEL_H_

#include "IntRegPrior.h"
#include "IntRegPar.h"
#include "DynamicModel.h"
#include "CoxModel.h"

/*------------------------------------------------------------------*/
/* Declaration of dynamic Cox model                                 */
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior>
class DynamicCoxModel:
  public DynamicModel<Prior, DynamicCoxPar>,
  public CoxModel<Prior, DynamicCoxPar>
{
public :
  DynamicCoxModel<Prior>(const boost::shared_ptr<IntRegData>& pd)
    : IntRegModel<Prior, DynamicCoxPar>(pd),
      DynamicModel<Prior, DynamicCoxPar>(pd),
      CoxModel<Prior, DynamicCoxPar>(pd) {}

  /* Overload virtual functions */
  DynamicCoxPar initPar() const;

  void gibbsKernel(const Prior& prior,
                   DynamicCoxPar& par);

  ublas::vector<double> likeVec(const DynamicCoxPar& par) const;
};
}

/*------------------------------------------------------------------*/
/* Implementation of dynamic Cox model                              */
/*------------------------------------------------------------------*/
namespace ir {
/* Initial dynamic Cox model parameters */
template <typename Prior>
DynamicCoxPar DynamicCoxModel<Prior>::
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

  return DynamicCoxPar(this->initLambda(),
                       ublas::matrix<double>(this->K_, this->nBeta_, 0),
                       ublas::vector<double>(this->nBeta_, 1.0), jumpMat);
}

/* Gibbs kernel implementation */
template <typename Prior>
void DynamicCoxModel<Prior>::
gibbsKernel(const Prior& prior,
            DynamicCoxPar& par)
{
  ublas::matrix<double> expXb(ublas::exp(ublas::prod(this->pd_->X(),
                                         ublas::trans(par.beta))));

  /* Sample dNMat and YMat */
  ublas::matrix<int> dNMat(this->N_, this->K_, 0);
  ublas::matrix<double> YMat(this->N_, this->K_, 1.0);
  this->sampleMat(par.lambda, expXb, dNMat, YMat);

  /* Sample lambda */
  this->sampleLambda(expXb, dNMat, YMat, prior.base_prior, par.lambda);

  /* Reversible jump MCMC for beta */
  ublas::vector<double> omega(this->N_, 1.0);

  /* Loop through covariates */
  for (Size j = 0; j < this->nBeta_; ++j) {
    double u = rmath::unif_rand();
    double nJump = ublas::sum(ublas::column(par.jump, j));

    bool prop_flag = false;
    double prop_ratio = 1.0;
    DynamicCoxPar prop_par(par);

    /* Birth move */
    if (u < this->prob_(0) && nJump < this->K_) {
      prop_ratio = this->propBirth(j, par.beta, par.jump, prop_par.beta,
                                   prop_par.jump);

      prop_flag = true;
    }

    /* Death move */
    if (u > this->prob_(0) && u < this->prob_(0) + this->prob_(1) && nJump > 1) {
      prop_ratio = this->propDeath(j, par.beta, par.jump, prop_par.beta,
                                   prop_par.jump);

      prop_flag = true;
    }

    /* Decide whether to accept */
    if (prop_flag == true) {
      double prior_ratio = std::exp(this->logCoefPrior(ublas::column(prop_par.jump, j),
                                                       ublas::column(prop_par.beta, j),
                                                       prior.coef_prior) -
                                    this->logCoefPrior(ublas::column(par.jump, j),
                                                       ublas::column(par.beta, j),
                                                       prior.coef_prior));

      double like_ratio = std::exp(ublas::sum(ublas::log(this->coxLikeVec(prop_par.lambda,
                                                                          prop_par.beta))) -
                                   ublas::sum(ublas::log(this->coxLikeVec(par.lambda,
                                                                          par.beta))));

      double ratio = prop_ratio * prior_ratio * like_ratio;

      if (rmath::unif_rand() < std::min<double>(1.0, ratio))
        par = prop_par;
    }

    if (prop_flag == false)
      this->sampleBeta(j, dNMat, YMat, par.lambda, par.jump, omega,
                 prior.coef_prior, par.beta, par.nu);
  }
}

/* Likelihood function */
template <typename Prior>
ublas::vector<double> DynamicCoxModel<Prior>::
likeVec(const DynamicCoxPar& par) const
{
  return this->coxLikeVec(par.lambda, par.beta);
}
}

#endif // DYNAMIC_COX_MODEL_H_
