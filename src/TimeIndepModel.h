//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef TIME_INDEP_MODEL_H_
#define TIME_INDEP_MODEL_H_

#include "IntRegModel.h"
#include "rng.h"

/*------------------------------------------------------------------*/
/* Declaration of time independent coef model class                 */
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior, typename Par>
class TimeIndepModel: virtual public IntRegModel<Prior, Par>
{
public :
  TimeIndepModel<Prior, Par>(const boost::shared_ptr<IntRegData>& pd)
    : IntRegModel<Prior, Par>(pd) {}

protected :
  /* Overload this function to sample with alternative prior */
  void sampleBeta(const ublas::vector<double>& lambda,
                  const ublas::matrix<int>& dNMat,
                  const ublas::matrix<double>& YMat,
                  const ublas::vector<double>& omega,
                  const NormalPrior& coef_prior,
                  ublas::vector<double>& beta);
};
}

/*------------------------------------------------------------------*/
/* Implementation of time indepedent coef model                     */
/*------------------------------------------------------------------*/
namespace ir {
/* Sample beta */
template <typename Prior, typename Par>
void TimeIndepModel<Prior, Par>::
sampleBeta(const ublas::vector<double>& lambda,
           const ublas::matrix<int>& dNMat,
           const ublas::matrix<double>& YMat,
           const ublas::vector<double>& omega,
           const NormalPrior& coef_prior,
           ublas::vector<double>& beta)
{
  /* Declare array paramters of log density function for beta */
  double *ldp_X = new double[this->N_];
  double *ldp_dleY = new double[this->N_];

  /* Set up parameters of arms_simple */
  int ninit = 4;
  double xl = -15.0, xr = 15.0;
  int dometrop = 0;
  double xprev = 0.0;

  /* error code from arms */
  int err;

  ublas::vector<double> DL(ublas::element_prod(this->delta_, lambda));
  ublas::vector<double> dNMat_rs(ublas::row_sum(dNMat));

  double ldp_sg2 = coef_prior.sd * coef_prior.sd;

  /* Loop through covariates */
  for (Size j = 0; j < this->nBeta_; ++j) {

    /* Prepare data for LogDenPar_type */
    double ldp_mu = ldp_sg2 * ublas::inner_prod(ublas::column(this->pd_->X(), j),
                    dNMat_rs);

    ublas::vector<double> temp_beta(beta);
    temp_beta(j) = 0;

    for (Size i = 0; i < this->N_; ++i) {
      ldp_X[i] = this->pd_->X()(i, j);

      ldp_dleY[i] = omega(i) * ublas::inner_prod(ublas::row(YMat, i), DL) *
        std::exp(ublas::inner_prod(ublas::row(this->pd_->X(), i), temp_beta));
    }

    /* Construct object of LogDenPar_type */
    struct IntRegModel<Prior, Par>::LogDenPar data = {
	ldp_mu, ldp_sg2,
	static_cast<int> (this->N_),
	ldp_X, ldp_dleY
    };
    double xsamp = 0.0;

    err = arms_simple(ninit, &xl, &xr, &IntRegModel<Prior, Par>::logDen,
                      &data, dometrop, &xprev, &xsamp);

    beta(j) = xsamp;
  }

  delete [] ldp_X;
  delete [] ldp_dleY;
}
}
#endif // TIME_INDEP_MODEL_H_
