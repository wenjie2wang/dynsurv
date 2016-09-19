//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef COX_MODEL_H_
#define COX_MODEL_H_

#include "IntRegModel.h"
#include "rng.h"

/*------------------------------------------------------------------*/
/* Declaration of Cox model                                         */
/*------------------------------------------------------------------*/
namespace ir {
  template <typename Prior, typename Par>
    class CoxModel: virtual public IntRegModel<Prior, Par>
    {
    public :
      CoxModel<Prior, Par>(const boost::shared_ptr<IntRegData>& pd)
        : IntRegModel<Prior, Par>(pd) {}

      typedef typename Prior::BasePrior_type BasePrior_type;
      typedef typename Prior::CoefPrior_type CoefPrior_type;

    protected :
      /* Helper functions for derived class */
      void sampleMat(const ublas::vector<double>& lambda,
                     const ublas::matrix<double>& expXb,
                     ublas::matrix<int>& dNMat,
                     ublas::matrix<double>& YMat);

      /* Overload this function to sample with alternative prior */
      void sampleLambda(const ublas::matrix<double>& expXb,
                        const ublas::matrix<int>& dNMat,
                        const ublas::matrix<double>& YMat,
                        const GammaPrior& base_prior,
                        ublas::vector<double>& lambda);

      void sampleLambda(const ublas::matrix<double>& expXb,
                        const ublas::matrix<int>& dNMat,
                        const ublas::matrix<double>& YMat,
                        const GammaProcessPrior& base_prior,
                        ublas::vector<double>& lambda);

      void sampleLambda(const ublas::matrix<double>& expXb,
                        const ublas::matrix<int>& dNMat,
                        const ublas::matrix<double>& YMat,
                        const ConstValuePrior& base_prior,
                        ublas::vector<double>& lambda);

      ublas::vector<double> coxLikeVec(const ublas::vector<double>& lambda,
                                       const ublas::matrix<double>& betaMat)
        const;
    };
}

/*------------------------------------------------------------------*/
/* Implementation of Cox model                                      */
/*------------------------------------------------------------------*/
namespace ir {
  /* Sample dNMat and YMat in Cox model */
  template <typename Prior, typename Par>
    void CoxModel<Prior, Par>::
    sampleMat(const ublas::vector<double>& lambda,
              const ublas::matrix<double>& expXb,
              ublas::matrix<int>& dNMat,
              ublas::matrix<double>& YMat)
    {
      ublas::matrix<double> expXbDL(expXb);
      for (Size k = 0; k < this->K_; ++k)
        ublas::column(expXbDL, k) *= this->delta_(k) * lambda(k);

      ublas::matrix<double> survMat(ublas::exp(- ublas::prod(expXbDL,
                                                             this->csMat_)));
      ublas::vector<double> eMat_i(this->K_);

      for (Size i = 0; i < this->N_; ++i) {

        /* Calculate eMat_i */
        for (Size k = 0; k < this->K_; ++k) {
          eMat_i(k) = (k == 0) ? (1 - survMat(i, k)) :
            (survMat(i, k - 1) - survMat(i, k));
          eMat_i(k) *= this->iMat_(i, k) * this->isIC_(i);
        }

        /* Sample dNMat_i for finite interval censored subject */
        Size k1 = 0;
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
                                          ublas::trans(this->csMat_)) *
          this->isIC_(i) + ublas::row(this->YMatRC_, i);
        if (this->isIC_(i) == 1 && this->isObs_(i) != 1) {
          double u = rmath::unif_rand();
          YMat(i, k1) = - log(1 - u + u * std::exp(- expXbDL(i, k1))) /
            expXbDL(i, k1);
        }
      }
    }

  /* Sample lambda in Cox model */
  template <typename Prior, typename Par>
    void CoxModel<Prior, Par>::
    sampleLambda(const ublas::matrix<double>& expXb,
                 const ublas::matrix<int>& dNMat,
                 const ublas::matrix<double>& YMat,
                 const GammaPrior& base_prior,
                 ublas::vector<double>& lambda)
    {
      for (Size k = 0; k < this->K_; ++k) {
        double shape = base_prior.shape + ublas::sum(ublas::column(dNMat, k));
        double rate = base_prior.rate +
          ublas::inner_prod(ublas::column(expXb, k), ublas::column(YMat, k)) *
          this->delta_(k);

        lambda(k) = rmath::rgamma(shape, 1.0 / rate);
      }
    }

  template <typename Prior, typename Par>
    void CoxModel<Prior, Par>::
    sampleLambda(const ublas::matrix<double>& expXb,
                 const ublas::matrix<int>& dNMat,
                 const ublas::matrix<double>& YMat,
                 const GammaProcessPrior& base_prior,
                 ublas::vector<double>& lambda)
    {
      for (Size k = 0; k < this->K_; ++k) {
        double shape = base_prior.mean * base_prior.control * this->delta_(k) +
          ublas::sum(ublas::column(dNMat, k));
        double rate = base_prior.control * this->delta_(k) +
          ublas::inner_prod(ublas::column(expXb, k),
                            ublas::column(YMat, k)) * this->delta_(k);

        lambda(k) = rmath::rgamma(shape, 1.0 / rate);
      }
    }

  template <typename Prior, typename Par>
    void CoxModel<Prior, Par>::
    sampleLambda(const ublas::matrix<double>& expXb,
                 const ublas::matrix<int>& dNMat,
                 const ublas::matrix<double>& YMat,
                 const ConstValuePrior& base_prior,
                 ublas::vector<double>& lambda)
    {
      lambda = ublas::vector<double>(lambda.size(), base_prior.value0);
    }

  /* Cox model likelihood */
  template <typename Prior, typename Par>
    ublas::vector<double> CoxModel<Prior, Par>::
    coxLikeVec(const ublas::vector<double>& lambda,
               const ublas::matrix<double>& betaMat) const
    {
      ublas::vector<double> DL(ublas::element_prod(this->delta_, lambda));
      ublas::vector<double> like(this->N_, 1.0), expXb_i(this->K_, 1.0),
        expXbDL_i(this->K_, 1.0);

      for (Size i = 0; i < this->N_; ++i) {
        expXb_i = ublas::exp(ublas::prod(betaMat,
                                         ublas::row(this->pd_->X(), i)));
        expXbDL_i = ublas::element_prod(expXb_i, DL);

        Size sumLeft = ublas::sum(ublas::column(this->lcsMat_, i));
        Size sumRight = ublas::sum(ublas::column(this->rcsMat_, i));

        // Event time is observed
        if (sumLeft == sumRight && sumLeft > 0) { //if (this->isObs_(i)) {
          Size k = sumLeft - 1;
          like(i) = lambda(k) * expXb_i(k) *
            std::exp(- ublas::inner_prod(expXbDL_i,
                                         ublas::column(this->lcsMat_, i)));
        }
        // Interval is observed
        else {
          like(i) =
            std::exp(- ublas::inner_prod(expXbDL_i,
                                         ublas::column(this->lcsMat_, i))) -
            std::exp(- ublas::inner_prod(expXbDL_i,
                                         ublas::column(this->rcsMat_, i))) *
            this->isIC_(i);
        }

        //if (like(i) <= 0)
        //std::cerr << "Non-positive likelihood evaluation!" << std::endl;
        //REprintf("Non-positive likelihood evaluation!")
      }

      return like;
    }
}

#endif // COX_MODEL_H_
