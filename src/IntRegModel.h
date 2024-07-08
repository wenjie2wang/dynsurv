//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef INTREG_MODEL_H_
#define INTREG_MODEL_H_

#include "IntRegData.h"
#include "ublas_ext.h"

#include <boost/smart_ptr/shared_ptr.hpp>

/*------------------------------------------------------------------*/
/* Declaration of abstract base model class                         */
/*------------------------------------------------------------------*/
namespace ir {
template <typename Prior, typename Par>
class IntRegModel
{
public:
  IntRegModel(const boost::shared_ptr<IntRegData>& pd);

  //~IntRegModel() {};

  typedef Prior Prior_type;
  typedef Par Par_type;

  /* Virtual functions to be defined in derived class */
  virtual Par_type initPar() const = 0;

  virtual void gibbsKernel(const Prior_type& prior,
                           Par_type& par) = 0;

  virtual ublas::vector<double> likeVec(const Par_type& par) const = 0;

  /* Accessor */
  Size N() const {
    return pd_->N();
  }

protected:
  /* Shared pointer to data */
  const boost::shared_ptr<IntRegData> pd_;

  ublas::matrix<int> iMat_;
  ublas::matrix<int> YMatRC_;
  ublas::matrix<int> csMat_;
  ublas::matrix<int> lcsMat_;
  ublas::matrix<int> rcsMat_;
  ublas::vector<int> isRC_;
  ublas::vector<int> isIC_;
  ublas::vector<int> isObs_;
  ublas::vector<double> delta_;
  const Size N_;
  const Size K_;
  const Size nBeta_;

  /* Helper functions for derived class */
  ublas::vector<double> initLambda() const;

  /* Parameters in the log density function of beta */
  struct LogDenPar
  {
    double ldp_mu;
    double ldp_sg2;
    int    ldp_N;
    double *ldp_X;
    double *ldp_dleY;
  };

  typedef LogDenPar LogDenPar_type;

  /* Log density function pointer */
  static double logDen(double x, void *data)
  {
    struct LogDenPar *d = static_cast<struct LogDenPar *>(data);

    double ld = - (x - d->ldp_mu) * (x - d->ldp_mu) / (2.0 * d->ldp_sg2);
    for (int i = 0; i < d->ldp_N; ++i)
      ld += - std::exp(d->ldp_X[i] * x) * d->ldp_dleY[i];

    return ld;
  }
};
}

/*------------------------------------------------------------------*/
/* Member function definitions of abstract base model class         */
/*------------------------------------------------------------------*/
namespace ir {
/* Constructor, prepare protect members to be used in derived class */
template <typename Prior, typename Par>
IntRegModel<Prior, Par>::
IntRegModel(const boost::shared_ptr<IntRegData>& pd)
  : pd_  (pd),
    iMat_  (pd->N(), pd->K(), 0),
    YMatRC_(pd->N(), pd->K(), 0),
    csMat_ (pd->K(), pd->K(), 0),
    lcsMat_(pd->K(), pd->N(), 0),
    rcsMat_(pd->K(), pd->N(), 0),
    isRC_  (pd->N(), 0),
    isIC_  (pd->N(), 0),
    isObs_ (pd->N(), 0),
    delta_ (pd->K(), 0.0),
    N_     (pd->N()),
    K_     (pd->K()),
    nBeta_ (pd->nBeta())
{
  double tau = pd->grid()(K_ - 1);

  /* Initialize indicator matrix and vector */
  for (Size i = 0; i < N_; ++i) {
    isRC_(i) = (pd->right()(i) > tau);
    isIC_(i) = (pd->right()(i) <= tau);
    if (isIC_(i) == 1 && pd->right()(i) - pd->left()(i) < 1e-8)
      isObs_(i) = 1;

    for (Size j = 0; j < K_; ++j) {
      iMat_(i, j) = ((pd->left()(i) < pd->grid()(j)) &&
                     (pd->right()(i) >= pd->grid()(j)));
      if (isObs_(i) == 1 && pd->right()(i) - pd->grid()(j) < 1e-8)
        iMat_(i, j) = 1;

      YMatRC_(i, j) = (pd->left()(i) >= pd->grid()(j)) * isRC_(i);
    }
  }

  for (Size i = 0; i < K_; ++i)
    for (Size j = 0; j < K_; ++j)
      csMat_(i, j) = (i <= j) ? 1 : 0;

  for (Size i = 0; i < K_; ++i) {
    for (Size j = 0; j < N_; ++j) {
      lcsMat_(i, j) = (pd->grid()(i) <= pd->left()(j)) ? 1 : 0;
      rcsMat_(i, j) = (pd->grid()(i) <= pd->right()(j)) ? 1 : 0;
    }
  }

  /* Initialize delta */
  this->delta_(0) = pd->grid()(0);
  for (Size j = 1; j < K_; ++j)
    this->delta_(j) = pd->grid()(j) - pd->grid()(j - 1);
}


/* Helper function to initialize lambda */
template <typename Prior, typename Par>
ublas::vector<double> IntRegModel<Prior, Par>::
initLambda() const
{
  ublas::matrix<double> eMat(iMat_);
  ublas::matrix<double> rMat(YMatRC_);

  for (Size i = 0; i < eMat.size1(); ++i) {
    double rowSum = ublas::sum(ublas::row(eMat, i));
    if (rowSum > 0) {
      ublas::row(eMat, i) *= (1 / rowSum) * isIC_(i);
    }

    ublas::row(rMat, i) += ublas::prod(csMat_, ublas::row(eMat, i)) *
                           isIC_(i);
  }

  ublas::vector<double> lambda(ublas::element_div(ublas::col_sum(eMat),
                               ublas::col_sum(rMat)));
  lambda = ublas::element_div(lambda, delta_);

  return lambda;
}
}
#endif // INTREG_MODEL_H_
