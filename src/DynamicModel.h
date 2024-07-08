//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once
#ifndef DYNAMIC_MODEL_H_
#define DYNAMIC_MODEL_H_

#include "IntRegModel.h"
#include "rng.h"

#include <limits>

/*------------------------------------------------------------------*/
/* Declaration of dynamic coef model                                */
/*------------------------------------------------------------------*/
namespace ir {
  template <typename Prior, typename Par>
    class DynamicModel: virtual public IntRegModel<Prior, Par>
    {
    public :
      /* Class constructor */
      DynamicModel(const boost::shared_ptr<IntRegData>& pd)
        : IntRegModel<Prior, Par>(pd),
        a0_  (100.0),
        prob_(2, 0.35),
        eps0_(0.5) {}

      /* Set functions */
      void set_a0(const double a) {
        a0_ = a;
      }
      void set_prob(const ublas::vector<double>& p) {
        prob_ = p;
      }
      void set_eps0(const double e) {
        eps0_ = e;
      }

    protected :
      /* Amplifying factor for the variance of beta_1 */
      double a0_;

      /* probability of birth and death */
      ublas::vector<double> prob_;

      /* u ~ uniform(-eps0, eps0) */
      double eps0_;

      /* Log of coef prior, to be overloaded with alternative coef prior */
      double logCoefPrior(const ublas::vector<int>& jump,
                          const ublas::vector<double>& beta,
                          const NormalProcessPrior& coef_prior) const;

      double logCoefPrior(const ublas::vector<int>& jump,
                          const ublas::vector<double>& beta,
                          const NormalInvGammaProcessPrior& coef_prior) const;

      /* Propose beta and jump */
      double propBirth(const Size j,
                       const ublas::matrix<double>& betaMat,
                       const ublas::matrix<int>& jumpMat,
                       ublas::matrix<double>& prop_betaMat,
                       ublas::matrix<int>& prop_jumpMat);

      double propDeath(const Size j,
                       const ublas::matrix<double>& betaMat,
                       const ublas::matrix<int>& jumpMat,
                       ublas::matrix<double>& prop_betaMat,
                       ublas::matrix<int>& prop_jumpMat);

      /* Overload this function to sample with alternative prior */
      void sampleBeta(const Size j,
                      const ublas::matrix<int>& dNMat,
                      const ublas::matrix<double>& YMat,
                      const ublas::vector<double>& lambda,
                      const ublas::matrix<int>& jumpMat,
                      const ublas::vector<double>& omega,
                      const NormalProcessPrior& coef_prior,
                      ublas::matrix<double>& betaMat,
                      ublas::vector<double>& nu);

      void sampleBeta(const Size j,
                      const ublas::matrix<int>& dNMat,
                      const ublas::matrix<double>& YMat,
                      const ublas::vector<double>& lambda,
                      const ublas::matrix<int>& jumpMat,
                      const ublas::vector<double>& omega,
                      const NormalInvGammaProcessPrior& coef_prior,
                      ublas::matrix<double>& betaMat,
                      ublas::vector<double>& nu);
    };
}

/*------------------------------------------------------------------*/
/* Implementation of dynamic model                                  */
/*------------------------------------------------------------------*/
namespace ir {
  /* Log value of NormalProcess prior */
  template <typename Prior, typename Par>
    double DynamicModel<Prior, Par>::
    logCoefPrior(const ublas::vector<int>& jump,
                 const ublas::vector<double>& beta,
                 const NormalProcessPrior& coef_prior) const
    {
      ublas::vector<double> sg2Vec(this->K_, coef_prior.sd * coef_prior.sd);
      for (Size k = 0; k < this->K_; ++k) {
        sg2Vec(k) *= a0_;
        if (jump(k) == 1)
          break;
      }

      double res = 0.0;
      double prev_b = 0.0;
      for (Size k = 0; k < this->K_; ++k) {
        if (jump(k) == 1) {
          res += - (beta(k) - prev_b) * (beta(k) - prev_b) / (2 * sg2Vec(k)) -
            0.5 * log(2 * M_PI * sg2Vec(k));

          prev_b = beta(k);
        }
      }

      return res;
    }

  /* Log value of NormalInvGammaProcess prior */
  template <typename Prior, typename Par>
    double DynamicModel<Prior, Par>::
    logCoefPrior(const ublas::vector<int>& jump,
                 const ublas::vector<double>& beta,
                 const NormalInvGammaProcessPrior& coef_prior) const
    {
      ublas::vector<double> sg2Vec(this->K_, coef_prior.scale);
      for (Size k = 0; k < this->K_; ++k) {
        sg2Vec(k) *= a0_;
        if (jump(k) == 1)
          break;
      }

      double res = 0.0;
      double prev_b = 0.0;
      for (Size k = 0; k < this->K_; ++k) {
        if (jump(k) == 1) {
          res += - (coef_prior.shape + 0.5) *
            log(1 +  (beta(k) - prev_b) *
                (beta(k) - prev_b) / (2 * sg2Vec(k))) -
            //log(rmath::beta(coef_prior.shape, 0.5)) -
            log(rmath::gammafn(coef_prior.shape) * rmath::gammafn(0.5) /
                rmath::gammafn(coef_prior.shape + 0.5)) -
            0.5 * log(2 * coef_prior.scale);

          prev_b = beta(k);
        }
      }

      return res;
    }

  /* Propose beta and jump in birth move */
  template <typename Prior, typename Par>
    double DynamicModel<Prior, Par>::
    propBirth(const Size j,
              const ublas::matrix<double>& betaMat,
              const ublas::matrix<int>& jumpMat,
              ublas::matrix<double>& prop_betaMat,
              ublas::matrix<int>& prop_jumpMat)
    {
      Size nJump = ublas::sum(ublas::column(jumpMat, j));

      /* Randomly select one non jump point */
      Size j0 = static_cast<Size>((this->K_ - nJump) * rmath::unif_rand()) + 1;
      Size cur = 0, counter = 0;
      for (Size k = 0; k < this->K_; ++k) {
        if (jumpMat(k, j) == 0)
          counter++;

        if (counter == j0) {
          cur = k;
          break;
        }
      }

      prop_jumpMat(cur, j) = 1;

      /* Subinterval range */
      Size end = this->K_ - 1;
      for (Size k = cur + 1; k < this->K_; ++k) {
        if (jumpMat(k, j) == 1) {
          end = k;
          break;
        }
      }

      Size sta = 0;
      for (Size k = cur; k > 0; --k) {
        if (jumpMat(k - 1, j) == 1) {
          sta = k;
          break;
        }
      }

      ublas::range rj(j, j + 1), r1(sta, cur + 1), r2(cur + 1, end + 1),
        r12(sta, end + 1);

      /* Weight proportional to width */
      double wt = ublas::sum(ublas::project(this->delta_, r1)) /
        ublas::sum(ublas::project(this->delta_, r12));


      /* Get beta from previous and next interval */
      double prev_b = (sta == 0) ? betaMat(sta, j) : betaMat(sta - 1, j);
      double next_b = (end == this->K_ - 1) ?
        betaMat(end, j) : betaMat(end + 1, j);

      /* Propose beta */
      double u = rmath::runif(- eps0_, eps0_);
      ublas::project(prop_betaMat, r1, rj) =
        ublas::matrix<double>(r1.size(), 1, prev_b * wt +
                              (betaMat(cur, j) + u) * (1 - wt));

      ublas::project(prop_betaMat, r2, rj) =
        ublas::matrix<double>(r2.size(), 1, (betaMat(cur, j) - u) *
                              wt + next_b * (1 - wt));

      /* Jacobian */
      double jacob = 2 * wt * (1 - wt);
      if (sta == 0)
        jacob += wt * wt;

      if (end == this->K_ - 1)
        jacob += (1 - wt) * (1 - wt);

      return jacob * 2 * eps0_;
    }

  /* Propose beta and jump in death move */
  template <typename Prior, typename Par>
    double DynamicModel<Prior, Par>::
    propDeath(const Size j,
              const ublas::matrix<double>& betaMat,
              const ublas::matrix<int>& jumpMat,
              ublas::matrix<double>& prop_betaMat,
              ublas::matrix<int>& prop_jumpMat)
    {
      Size nJump = ublas::sum(ublas::column(jumpMat, j));

      /* Randomly select one jump point */
      // No death move for the last jump point
      Size j1 = static_cast<Size>((nJump - 1) * rmath::unif_rand()) + 1;
      Size cur = 0, counter = 0;
      for (Size k = 0; k < this->K_; ++k) {
        if (jumpMat(k, j) == 1)
          counter++;

        if (counter == j1) {
          cur = k;
          break;
        }
      }

      prop_jumpMat(cur, j) = 0;

      /* Subinterval range */
      Size end = this->K_ - 1;
      for (Size k = cur + 1; k < this->K_; ++k) {
        if (jumpMat(k, j) == 1) {
          end = k;
          break;
        }
      }

      Size sta = 0;
      for (Size k = cur; k > 0; --k) {
        if (jumpMat(k - 1, j) == 1) {
          sta = k;
          break;
        }
      }

      ublas::range rj(j, j + 1), r1(sta, cur + 1), r2(cur + 1, end + 1),
        r12(sta, end + 1);

      /* Weight proportional to width */
      double wt = ublas::sum(ublas::project(this->delta_, r1)) /
        ublas::sum(ublas::project(this->delta_, r12));

      /* Get beta from previous and next interval */
      double prev_b = (sta == 0) ? betaMat(sta, j) : betaMat(sta - 1, j);
      double next_b = (end == this->K_ - 1) ?
        betaMat(end, j) : betaMat(end + 1, j);

      /* Propose beta */
      ublas::project(prop_betaMat, r12, rj) =
        ublas::matrix<double>(r12.size(), 1, 0.5 *
                              (prev_b * (- wt / (1 - wt)) + betaMat(cur, j) *
                               (1.0 / (1- wt)) + betaMat(end, j) * (1.0 / wt) +
                               next_b * (- (1 - wt) / wt)));

      /* Jacobian */
      double jacob = 1.0 / (2 * wt * (1 - wt));
      if (sta == 0)
        jacob *= 1 - wt;

      if (end == this->K_ - 1)
        jacob *= wt;

      return jacob / (2 * eps0_);
    }

  /* Sample beta with NormalProcess prior */
  template <typename Prior, typename Par>
    void DynamicModel<Prior, Par>::
    sampleBeta(const Size j,
               const ublas::matrix<int>& dNMat,
               const ublas::matrix<double>& YMat,
               const ublas::vector<double>& lambda,
               const ublas::matrix<int>& jumpMat,
               const ublas::vector<double>& omega,
               const NormalProcessPrior& coef_prior,
               ublas::matrix<double>& betaMat,
               ublas::vector<double>& nu)
    {
      Size nJump = ublas::sum(ublas::column(jumpMat, j));

      std::vector<Size> sta, end;
      std::vector<ublas::range> rg;
      sta.push_back(0);
      for (Size k = 0; k < this->K_; ++k) {
        if (jumpMat(k, j) == 1) {
          end.push_back(k);
          rg.push_back(ublas::range(sta.back(), k + 1));

          if (k < this->K_ - 1)
            sta.push_back(k + 1);
        }
      }

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

      nu(j) = coef_prior.sd * coef_prior.sd;
      ublas::vector<double> shr_sg2Vec(nJump, coef_prior.sd * coef_prior.sd);
      shr_sg2Vec(0) *= a0_;

      /* Loop through time intervals */
      for (Size s = 0; s < nJump; ++s) {

        /* Prepare data for LogDenPar_type */
        for (Size i = 0; i < this->N_; ++i) {
          ldp_X[i] = this->pd_->X()(i, j);
          ldp_dleY[i] = 0;

          for (Size k = sta[s]; k < end[s] + 1; ++k) {
            ublas::vector<double> temp_beta(ublas::row(betaMat, k));
            temp_beta(j) = 0;

            ldp_dleY[i] += omega(i) * this->delta_(k) * lambda(k) *
              std::exp(ublas::inner_prod(ublas::row(this->pd_->X(), i),
                                         temp_beta)) * YMat(i, k);
          }
        }

        double cur_sg2 = shr_sg2Vec(s);
        double next_sg2 = (s < nJump - 1) ? shr_sg2Vec(s + 1) : INF;

        double prev_b = (s == 0) ? 0.0 : betaMat(sta[s] - 1, j);
        double next_b = (s < nJump - 1) ? betaMat(end[s] + 1, j) : 1.0;

        double ldp_sg2 = 1.0 / (1.0 / cur_sg2 + 1.0 / next_sg2);
        double ldp_mu = ldp_sg2 *
          (ublas::sum(ublas::prod(ublas::column(this->pd_->X(), j),
                                  project(dNMat, ublas::range(0, this->N_),
                                          rg[s]))) +
           prev_b / cur_sg2 + next_b / next_sg2);

        /* Construct object of LogDenPar_type */
        struct IntRegModel<Prior, Par>::LogDenPar data = {
          ldp_mu, ldp_sg2,
          static_cast<int> (this->N_),
          ldp_X, ldp_dleY
        };
        double xsamp = 0.0;

        err = arms_simple(ninit, &xl, &xr, &IntRegModel<Prior, Par>::logDen,
                          &data, dometrop, &xprev, &xsamp);

        ublas::project(betaMat, rg[s], ublas::range(j, j + 1)) =
          ublas::matrix<double>(rg[s].size(), 1, xsamp);
      }

      delete [] ldp_X;
      delete [] ldp_dleY;
    }

  /* Sample beta with NormalInvGammaProcess prior */
  template <typename Prior, typename Par>
    void DynamicModel<Prior, Par>::
    sampleBeta(const Size j,
               const ublas::matrix<int>& dNMat,
               const ublas::matrix<double>& YMat,
               const ublas::vector<double>& lambda,
               const ublas::matrix<int>& jumpMat,
               const ublas::vector<double>& omega,
               const NormalInvGammaProcessPrior& coef_prior,
               ublas::matrix<double>& betaMat,
               ublas::vector<double>& nu)
    {
      Size nJump = ublas::sum(ublas::column(jumpMat, j));

      std::vector<Size> sta, end;
      std::vector<ublas::range> rg;
      sta.push_back(0);
      for (Size k = 0; k < this->K_; ++k) {
        if (jumpMat(k, j) == 1) {
          end.push_back(k);
          rg.push_back(ublas::range(sta.back(), k + 1));

          if (k < this->K_ - 1)
            sta.push_back(k + 1);
        }
      }

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

      /* Sample tau (scale of the variance) */
      ublas::vector<double> shr_sg2Vec(nJump, 1.0);
      shr_sg2Vec(0) *= a0_;

      double shape = coef_prior.shape + nJump / 2;
      double scale = coef_prior.scale;

      double prev_b = 0.0;
      for (Size s = 0; s < nJump; ++s) {
        scale += (betaMat(sta[s], j) - prev_b) * (betaMat(sta[s], j) - prev_b) /
          (2 * shr_sg2Vec(s));
        prev_b = betaMat(sta[s], j);
      }

      nu(j) = 1.0 / rmath::rgamma(shape, 1.0 / scale);
      shr_sg2Vec *= nu(j);

      /* Loop through time intervals */
      for (Size s = 0; s < nJump; ++s) {

        /* Prepare data for LogDenPar_type */
        for (Size i = 0; i < this->N_; ++i) {
          ldp_X[i] = this->pd_->X()(i, j);
          ldp_dleY[i] = 0;

          for (Size k = sta[s]; k < end[s] + 1; ++k) {
            ublas::vector<double> temp_beta(ublas::row(betaMat, k));
            temp_beta(j) = 0;

            ldp_dleY[i] += omega(i) * this->delta_(k) * lambda(k) *
              std::exp(ublas::inner_prod(ublas::row(this->pd_->X(), i),
                                         temp_beta)) * YMat(i, k);
          }
        }

        double cur_sg2 = shr_sg2Vec(s);
        double next_sg2 = (s < nJump - 1) ? shr_sg2Vec(s + 1) : INF;

        double prev_b = (s == 0) ? 0.0 : betaMat(sta[s] - 1, j);
        double next_b = (s < nJump - 1) ? betaMat(end[s] + 1, j) : 1.0;

        double ldp_sg2 = 1.0 / (1.0 / cur_sg2 + 1.0 / next_sg2);
        double ldp_mu = ldp_sg2 *
          (ublas::sum(ublas::prod(ublas::column(this->pd_->X(), j),
                                  project(dNMat, ublas::range(0, this->N_),
                                          rg[s]))) +
           prev_b / cur_sg2 + next_b / next_sg2);

        /* Construct object of LogDenPar_type */
        struct IntRegModel<Prior, Par>::LogDenPar data = {
          ldp_mu, ldp_sg2,
          static_cast<int> (this->N_),
          ldp_X, ldp_dleY
        };
        double xsamp = 0.0;

        err = arms_simple(ninit, &xl, &xr, &IntRegModel<Prior, Par>::logDen,
                          &data, dometrop, &xprev, &xsamp);

        ublas::project(betaMat, rg[s], ublas::range(j, j + 1)) =
          ublas::matrix<double>(rg[s].size(), 1, xsamp);
      }

      delete [] ldp_X;
      delete [] ldp_dleY;
    }
}

#endif // DYNAMIC_MODEL_H_
