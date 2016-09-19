//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once
#ifndef TIME_VARYING_MODEL_H_
#define TIME_VARYING_MODEL_H_

#include "IntRegPrior.h"
#include "rng.h"

#include <limits>

/*------------------------------------------------------------------*/
/* Declaration of time varying coef model                           */
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior, typename Par>
class TimeVaryingModel: virtual public IntRegModel<Prior, Par>
{
public :
  TimeVaryingModel<Prior, Par>(const boost::shared_ptr<IntRegData>& pd)
    : IntRegModel<Prior, Par>(pd),
      a0_(100.0) {}

  /* Set functions */
  void set_a0(const double a) {
    a0_ = a;
  }

protected :
  /* Amplifying factor for the variance of beta_1 */
  double a0_;

  /* Overload this function to sample with alternative prior */
  void sampleBeta(const ublas::vector<double>& tLambda,
                  const ublas::matrix<int>& dNMat,
                  const ublas::matrix<double>& YMat,
                  const ublas::vector<double>& omega,
                  const NormalProcessPrior& coef_prior,
                  ublas::matrix<double>& betaMat,
                  ublas::vector<double>& nu);

  void sampleBeta(const ublas::vector<double>& tLambda,
                  const ublas::matrix<int>& dNMat,
                  const ublas::matrix<double>& YMat,
                  const ublas::vector<double>& omega,
                  const NormalInvGammaProcessPrior& coef_prior,
                  ublas::matrix<double>& betaMat,
                  ublas::vector<double>& nu);
};
}

/*------------------------------------------------------------------*/
/* Implementation of time varying coef model                        */
/*------------------------------------------------------------------*/
namespace ir {
/* Sample beta with NormalProcess prior */
template <typename Prior, typename Par>
void TimeVaryingModel<Prior, Par>::
sampleBeta(const ublas::vector<double>& tLambda,
           const ublas::matrix<int>& dNMat,
           const ublas::matrix<double>& YMat,
           const ublas::vector<double>& omega,
           const NormalProcessPrior& coef_prior,
           ublas::matrix<double>& betaMat,
           ublas::vector<double>& nu)
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

  const double INF = std::numeric_limits<double>::max();

  nu = ublas::vector<double>(this->nBeta_, coef_prior.sd * coef_prior.sd);
  ublas::matrix<double> sg2Mat(this->K_, this->nBeta_, coef_prior.sd * coef_prior.sd);
  ublas::row(sg2Mat, 0) *= a0_;

  /* Loop through time and covariate */
  for (Size k = 0; k < this->K_; ++k) {
    for (Size j = 0; j < this->nBeta_; ++j) {

      /* Prepare data for LogDenPar_type */
      ublas::vector<double> temp_beta(ublas::row(betaMat, k));
      temp_beta(j) = 0;

      for (Size i = 0; i < this->N_; ++i) {
        ldp_X[i] = this->pd_->X()(i, j);

        ldp_dleY[i] = omega(i) * this->delta_(k) * tLambda(k) *
          std::exp(ublas::inner_prod(ublas::row(this->pd_->X(), i), temp_beta)) *
          YMat(i, k);
      }

      double cur_sg2 = sg2Mat(k, j);
      double next_sg2 = (k < this->K_ - 1) ? sg2Mat(k + 1, j) : INF;

      double ldp_sg2 = 1.0 / (1.0 / cur_sg2 + 1.0 / next_sg2);

      double prev_b = (k == 0) ? 0.0 : betaMat(k - 1, j);
      double next_b = (k < this->K_ - 1) ? betaMat(k + 1, j) : 1.0;

      double ldp_mu = ldp_sg2 * (ublas::inner_prod(ublas::column(this->pd_->X(), j),
                                 ublas::column(dNMat, k)) + prev_b / cur_sg2 +
                                 next_b / next_sg2);

      /* Construct object of LogDenPar_type */
      struct IntRegModel<Prior, Par>::LogDenPar data = {
	  ldp_mu, ldp_sg2, 
	  static_cast<int> (this->N_),
	  ldp_X, ldp_dleY
      };
      double xsamp = 0.0;

      err = arms_simple(ninit, &xl, &xr, &IntRegModel<Prior, Par>::logDen,
                        &data, dometrop, &xprev, &xsamp);

      betaMat(k, j) = xsamp;
    }
  }

  delete [] ldp_X;
  delete [] ldp_dleY;
}

/* Sample beta with NormalInvGammaProcess prior */
template <typename Prior, typename Par>
void TimeVaryingModel<Prior, Par>::
sampleBeta(const ublas::vector<double>& tLambda,
           const ublas::matrix<int>& dNMat,
           const ublas::matrix<double>& YMat,
           const ublas::vector<double>& omega,
           const NormalInvGammaProcessPrior& coef_prior,
           ublas::matrix<double>& betaMat,
           ublas::vector<double>& nu)
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

  const double INF = std::numeric_limits<double>::max();

  /* Sample nv (scale of the variance) */
  ublas::matrix<double> sg2Mat(this->K_, this->nBeta_, 1.0);
  ublas::row(sg2Mat, 0) *= a0_;

  double shape = coef_prior.shape + this->K_ / 2;

  for (Size j = 0; j < this->nBeta_; ++j) {
    double scale = coef_prior.scale;

    double prev_b = 0.0;
    for (Size k = 0; k < this->K_; ++ k) {
      scale += (betaMat(k, j) - prev_b) * (betaMat(k, j) - prev_b) /
               (2 * sg2Mat(k, j));
      prev_b = betaMat(k, j);
    }

    nu(j) = 1.0 / rmath::rgamma(shape, 1.0 / scale);
    ublas::column(sg2Mat, j) *= nu(j);
  }

  /* Loop through time and covariate */
  for (Size k = 0; k < this->K_; ++k) {
    for (Size j = 0; j < this->nBeta_; ++j) {

      /* Prepare data for LogDenPar_type */
      ublas::vector<double> temp_beta(ublas::row(betaMat, k));
      temp_beta(j) = 0;

      for (Size i = 0; i < this->N_; ++i) {
        ldp_X[i] = this->pd_->X()(i, j);

        ldp_dleY[i] = omega(i) * this->delta_(k) * tLambda(k) *
          std::exp(ublas::inner_prod(ublas::row(this->pd_->X(), i), temp_beta)) *
          YMat(i, k);
      }

      double cur_sg2 = sg2Mat(k, j);
      double next_sg2 = (k < this->K_ - 1) ? sg2Mat(k + 1, j) : INF;

      double ldp_sg2 = 1.0 / (1.0 / cur_sg2 + 1.0 / next_sg2);

      double prev_b = (k == 0) ? 0.0 : betaMat(k - 1, j);
      double next_b = (k < this->K_ - 1) ? betaMat(k + 1, j) : 1.0;

      double ldp_mu = ldp_sg2 * (ublas::inner_prod(ublas::column(this->pd_->X(), j),
                                 ublas::column(dNMat, k)) + prev_b / cur_sg2 +
                                 next_b / next_sg2);

      /* Construct object of LogDenPar_type */
      struct IntRegModel<Prior, Par>::LogDenPar data = {
	  ldp_mu, ldp_sg2, 
	  static_cast<int> (this->N_),
	  ldp_X, ldp_dleY
      };
      double xsamp = 0.0;

      err = arms_simple(ninit, &xl, &xr, &IntRegModel<Prior, Par>::logDen,
                        &data, dometrop, &xprev, &xsamp);

      betaMat(k, j) = xsamp;
    }
  }

  delete [] ldp_X;
  delete [] ldp_dleY;
}
}

#endif // TIME_VARYING_MODEL_H_
