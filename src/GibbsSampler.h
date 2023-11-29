//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once
#ifndef GIBBS_SAMPLER_H_
#define GIBBS_SAMPLER_H_

//extern "C" {
#include <R_ext/Print.h>
//}
#include "ublas.h"
#include "ublas_ext.h"

#include <iostream>
#include <vector>
#include <boost/smart_ptr/shared_ptr.hpp>

// 3 public member function of M required: initPar, gibbsKernel, likeVec
// 1 template function required: mean

/*------------------------------------------------------------------*/
/* Declaration of Gibbs sampler class                               */
/*------------------------------------------------------------------*/
namespace ir {
// M: interval regression model class
template<class M>
class GibbsSampler
{
public :
  typedef typename M::Prior_type Prior_type;
  typedef typename M::Par_type Par_type;
  typedef std::vector<Par_type> VectorPar_type;

  GibbsSampler(const boost::shared_ptr<M>& pm,
               const Size iter)
    : pm_(pm),
      iter_(iter),
      N_   (pm->N()) {}

  void runGibbs(const Prior_type& prior,
                bool trace = true,
                Size nReport = 1);

  std::ostream& summaryFit(std::ostream& out,
                           Size burn = 25,
                           Size thin = 1) const;

  void summaryFitR(Size burn,
                   Size thin,
                   Par_type& par,
                   double& LPML,
                   double& DHat,
                   double& DBar,
                   double& pD,
                   double& DIC);

  std::ofstream& outputSample(std::ofstream& out) const;

  ~GibbsSampler() {}

private :
  const boost::shared_ptr<M> pm_;
  const Size iter_;
  const Size N_;
  VectorPar_type samples_;
};
}

/*------------------------------------------------------------------*/
/* Member function definitions of Gibbs sampler class               */
/*------------------------------------------------------------------*/
namespace ir {
/* Run Gibbs sampler */
template<class M>
void GibbsSampler<M>::runGibbs(const Prior_type& prior,
                               bool trace,
                               Size nReport)
{
  Par_type par(pm_->initPar());

  for (Size i = 0; i < iter_; ++i) {

    if (trace && (i % nReport == 0)) {
      //std::cout << "Iteration(" << i << ")" << std::endl;
      Rprintf("Iteration(%u)\n", static_cast<unsigned int>(i));
    }

    pm_->gibbsKernel(prior, par);

    samples_.push_back(par);
  }
}

/* Summary fit statistics */
template<class M>
std::ostream& GibbsSampler<M>::summaryFit(std::ostream& out,
    Size burn,
    Size thin) const
{
  /* Select a subset of samples_*/
  if (burn >= iter_)
    //std::cerr << "burn must be smaller than iter!" << std::endl;
    REprintf("burn must be smaller than iter!\n");

  const Size n = (iter_ - burn) / thin;
  std::vector<int> sq(n);

  for (Size i = 0; i < n; ++i)
    sq[i] = burn + i * thin;

  /* Calculate observed likelihood matrix */
  ublas::matrix<double> obsLikeMat(n, N_);

  for (Size i = 0; i < n; ++i)
    ublas::row(obsLikeMat, i) = pm_->likeVec(samples_[sq[i]]);

  /* CPO */
  ublas::matrix<double> invObsLikeMat(ublas::element_div(
                                        ublas::matrix<double>(n, N_, 1.0), obsLikeMat));
  ublas::vector<double> CPO(ublas::element_div(
                              ublas::vector<double>(N_, 1.0), ublas::col_mean(invObsLikeMat)));

  /* log pseudo likelihood*/
  double LPML = ublas::sum(ublas::log(CPO));

  /* DHat */
  Par_type par(mean<Par_type>(samples_));
  double DHat = - 2 * ublas::sum(ublas::log(pm_->likeVec(par)));

  /* DBar */
  ublas::matrix<double> logObsLikeMat(ublas::log(obsLikeMat));
  double DBar = - 2 * ublas::sum(ublas::col_mean(logObsLikeMat));

  double pD = DBar - DHat;
  double DIC = DBar + pD;

  out	<< "\nGibbs samples"
      << "\n start      : " << burn + 1
      << "\n end        : " << iter_
      << "\n thin       : " << thin
      << "\n summarized : " << n << '\n';

  out << "\nBayesian estimates\n";
  par.print(out);

  out	<< "\nStatistics"
      << "\n LPML : " << LPML
      << "\n DHat : " << DHat
      << "\n DBar : " << DBar
      << "\n pD   : " << pD
      << "\n DIC  : " << DIC << '\n';

  return out;
}

/* Summary fit statistics, R version */
template<class M>
void GibbsSampler<M>::summaryFitR(Size burn,
                                  Size thin,
                                  Par_type& par,
                                  double& LPML,
                                  double& DHat,
                                  double& DBar,
                                  double& pD,
                                  double& DIC)
{
  /* Select a subset of samples_*/
  if (burn >= iter_)
    //std::cerr << "burn must be smaller than iter!" << std::endl;
    REprintf("burn must be smaller than iter!\n");

  const Size n = (iter_ - burn) / thin;
  std::vector<int> sq(n);

  for (Size i = 0; i < n; ++i)
    sq[i] = burn + i * thin;

  /* Calculate observed likelihood matrix */
  ublas::matrix<double> obsLikeMat(n, N_);

  for (Size i = 0; i < n; ++i) {
    ublas::row(obsLikeMat, i) = pm_->likeVec(samples_[sq[i]]);
  }

  /* CPO */
  ublas::matrix<double> invObsLikeMat(ublas::element_div(
                                        ublas::matrix<double>(n, N_, 1.0), obsLikeMat));
  ublas::vector<double> CPO(ublas::element_div(
                              ublas::vector<double>(N_, 1.0), ublas::col_mean(invObsLikeMat)));

  /* log pseudo likelihood*/
  LPML = ublas::sum(ublas::log(CPO));

  /* DHat */
  par = mean<Par_type>(samples_);
  DHat = - 2 * ublas::sum(ublas::log(pm_->likeVec(par)));

  /* DBar */
  ublas::matrix<double> logObsLikeMat(ublas::log(obsLikeMat));
  DBar = - 2 * ublas::sum(ublas::col_mean(logObsLikeMat));

  pD = DBar - DHat;
  DIC = DBar + pD;
}

/* Output samples to file */
template<class M>
std::ofstream& GibbsSampler<M>::outputSample(std::ofstream& out) const
{
  for (Size i = 0; i < samples_.size(); ++i) {
    samples_[i].output(out);
    out << '\n';
  }

  return out;
}
}
#endif // GIBBS_SAMPLER_H_
