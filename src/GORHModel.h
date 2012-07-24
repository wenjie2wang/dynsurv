//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once
#ifndef GORH_MODEL_H_
#define GORH_MODEL_H_

#include "IntRegModel.h"
#include "rng.h"

/*------------------------------------------------------------------*/
/* Declaration of GORH model                                        */
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior, typename Par>
class GORHModel: virtual public IntRegModel<Prior, Par>
{
public :
  GORHModel<Prior, Par>(const boost::shared_ptr<IntRegData>& pd)
    : IntRegModel<Prior, Par>(pd) {}

  typedef typename Prior::BasePrior_type BasePrior_type;
  typedef typename Prior::CoefPrior_type CoefPrior_type;
  typedef typename Prior::ThetaPrior_type ThetaPrior_type;

protected :
  /* Helper functions for derived class */
  // tLambda = theta * lambda
  void sampleMat(const ublas::vector<double>& tLambda,
                 const ublas::matrix<double>& expXb,
                 const double theta,
                 ublas::matrix<int>& dNMat,
                 ublas::matrix<double>& YMat);

  /* Overload this function to sample with alternative prior */
  void sampleTheta(const ublas::matrix<double>& expXb,
                   const ublas::matrix<int>& dNMat,
                   const ublas::matrix<double>& YMat,
                   const ublas::vector<double>& tLambda,
                   const GammaPrior& base_prior,
                   const InvGammaPrior& theta_prior,
                   double& theta);

  /* Overload this function to sample with alternative prior */
  void sampleTheta(const ublas::matrix<double>& expXb,
                   const ublas::matrix<int>& dNMat,
                   const ublas::matrix<double>& YMat,
                   const ublas::vector<double>& tLambda,
                   const BasePrior_type& base_prior,
                   const ConstValuePrior& theta_prior,
                   double& theta);

  void sampleOmega(const ublas::matrix<double>& expXb,
                   const ublas::matrix<int>& dNMat,
                   const ublas::matrix<double>& YMat,
                   const ublas::vector<double>& tLambda,
                   const double theta,
                   ublas::vector<double>& omega);

  /* Overload this function to sample with alternative prior */
  void sampleTLambda(const ublas::matrix<double>& expXb,
                     const ublas::matrix<int>& dNMat,
                     const ublas::matrix<double>& YMat,
                     const ublas::vector<double>& omega,
                     const double theta,
                     const GammaPrior& base_prior,
                     ublas::vector<double>& tLambda);

  ublas::vector<double> gorhLikeVec(const ublas::vector<double>& tLambda,
                                    const ublas::matrix<double>& betaMat,
                                    const double theta) const;
};
}

/*------------------------------------------------------------------*/
/* Implementation of GORH model                                     */
/*------------------------------------------------------------------*/
namespace ir {
/* Sample dNMat and YMat in GORH model */
template <typename Prior, typename Par>
void GORHModel<Prior, Par>::
sampleMat(const ublas::vector<double>& tLambda,
          const ublas::matrix<double>& expXb,
          const double theta,
          ublas::matrix<int>& dNMat,
          ublas::matrix<double>& YMat)
{
  ublas::matrix<double> expXbDL(expXb);
  for (Size k = 0; k < this->K_; ++k)
    ublas::column(expXbDL, k) *= this->delta_(k) * tLambda(k);

  ublas::matrix<double> survMat(ublas::exp(- 1 / theta * ublas::log(
                                  ublas::matrix<double>(this->N_, this->K_, 1.0) +
                                  ublas::prod(expXbDL, this->csMat_))));

  ublas::vector<double> eMat_i(this->K_);

  for (Size i = 0; i < this->N_; ++i) {

    /* Calculate eMat_i */
    for (Size k = 0; k < this->K_; ++k) {
      eMat_i(k) = (k == 0) ? (1 - survMat(i, k)) :
                  (survMat(i, k - 1) - survMat(i, k));
      eMat_i(k) *= this->iMat_(i, k) * this->isIC_(i);
    }

    /* Sample dNMat_i for finite interval censored subject */
    Size k1;
    if (this->isIC_(i) == 1) {
      eMat_i /= ublas::sum(eMat_i);

      double cumSum = 0;
      double u = rmath::unif_rand();
      for (Size k = 0; k < this->K_; ++k) {
        cumSum += eMat_i(k);
        if (cumSum > u) {
          dNMat(i, k) = 1;
          k1 = k;
          break;
        }
      }
    }

    /* Sample YMat_i for finite interval censored subject */
    ublas::row(YMat, i) = ublas::prod(ublas::row(dNMat, i),
                                      ublas::trans(this->csMat_)) * this->isIC_(i) + ublas::row(this->YMatRC_, i);
    if (this->isIC_(i) == 1) {
      double a = expXbDL(i, k1) / (1 + ublas::sum(ublas::project(
                                     ublas::row(expXbDL, i), ublas::range(0, k1))));
      double u = rmath::unif_rand();

      YMat(i, k1) = (pow(1 - u + u * pow(1 + a, - 1 / theta), - theta) - 1) / a;
    }
  }
}

/* Sample theta in GORH model */
template <typename Prior, typename Par>
void GORHModel<Prior, Par>::
sampleTheta(const ublas::matrix<double>& expXb,
            const ublas::matrix<int>& dNMat,
            const ublas::matrix<double>& YMat,
            const ublas::vector<double>& tLambda,
            const GammaPrior& base_prior,
            const InvGammaPrior& theta_prior,
            double& theta)
{
  ublas::vector<double> DL(ublas::element_prod(this->delta_, tLambda));
  ublas::matrix<double> eY(ublas::element_prod(expXb, YMat));

  /* Shape and scale parameters for inverse gamma*/
  double shape = theta_prior.shape + this->K_ * base_prior.shape +
                 ublas::sum(ublas::row_sum(dNMat));
  double scale = theta_prior.scale + base_prior.rate * ublas::sum(tLambda) +
                 ublas::sum(ublas::log(ublas::vector<double>(this->N_, 1.0) + ublas::prod(eY, DL)));

  theta = 1.0 / rmath::rgamma(shape, 1.0 / scale);
}

template <typename Prior, typename Par>
void GORHModel<Prior, Par>::
sampleTheta(const ublas::matrix<double>& expXb,
            const ublas::matrix<int>& dNMat,
            const ublas::matrix<double>& YMat,
            const ublas::vector<double>& tLambda,
            const BasePrior_type& base_prior,
            const ConstValuePrior& theta_prior,
            double& theta)
{
  theta = theta_prior.value0;
}

/* Sample omega in GORH model */
template <typename Prior, typename Par>
void GORHModel<Prior, Par>::
sampleOmega(const ublas::matrix<double>& expXb,
            const ublas::matrix<int>& dNMat,
            const ublas::matrix<double>& YMat,
            const ublas::vector<double>& tLambda,
            const double theta,
            ublas::vector<double>& omega)
{
  ublas::vector<double> DL(ublas::element_prod(this->delta_, tLambda));
  ublas::matrix<double> eY(ublas::element_prod(expXb, YMat));

  for (Size i = 0; i < this->N_; ++i) {
    double shape = 1.0 / theta + ublas::sum(ublas::row(dNMat, i));
    double rate = 1.0 + ublas::inner_prod(DL, ublas::row(eY, i));

    omega(i) = rmath::rgamma(shape, 1.0 / rate);
  }
}

/* Sample lambda in GORH model */
template <typename Prior, typename Par>
void GORHModel<Prior, Par>::
sampleTLambda(const ublas::matrix<double>& expXb,
              const ublas::matrix<int>& dNMat,
              const ublas::matrix<double>& YMat,
              const ublas::vector<double>& omega,
              const double theta,
              const GammaPrior& base_prior,
              ublas::vector<double>& tLambda)
{
  for (Size k = 0; k < this->K_; ++k) {
    double shape = base_prior.shape + ublas::sum(ublas::column(dNMat, k));

    double doeY = 0.0;
    for (Size i = 0; i < this->N_; ++i)
      doeY += this->delta_(k) * omega(i) * expXb(i, k) * YMat(i, k);

    double rate = base_prior.rate / theta + doeY;

    tLambda(k) = rmath::rgamma(shape, 1.0 / rate);
  }
}

/* GORH model likelihood */
template <typename Prior, typename Par>
ublas::vector<double> GORHModel<Prior, Par>::
gorhLikeVec(const ublas::vector<double>& tLambda,
            const ublas::matrix<double>& betaMat,
            const double theta) const
{
  ublas::vector<double> DL(ublas::element_prod(this->delta_, tLambda));
  ublas::vector<double> like(this->N_, 1.0), expXb_i(this->K_, 1.0),
        expXbDL_i(this->K_, 1.0);

  for (Size i = 0; i < this->N_; ++i) {
    expXb_i = ublas::exp(ublas::prod(betaMat, ublas::row(this->pd_->X(), i)));
    expXbDL_i = ublas::element_prod(expXb_i, DL);
    like(i) = pow(1 + ublas::inner_prod(expXbDL_i, ublas::column(this->lcsMat_, i)),
                  - 1 / theta) - pow(1 + ublas::inner_prod(expXbDL_i,
                                     ublas::column(this->rcsMat_, i)), - 1 / theta) * this->isIC_(i);
  }

  return like;
}
}
#endif // GORH_MODEL_H_
