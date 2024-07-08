//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once
#ifndef INTREG_PRIOR_H_
#define INTREG_PRIOR_H_

#include <iostream>

/*------------------------------------------------------------------*/
/* Definition of basic prior class                                  */
/*------------------------------------------------------------------*/
namespace ir {
struct IntRegPrior
{
  virtual std::ostream& print(std::ostream& out) const = 0;
};

struct ConstValuePrior: public IntRegPrior
{
  ConstValuePrior(double v = 1.0)
    : value0(v) {}

  double value0;

  std::ostream& print(std::ostream& out) const
  {
    out	<< "ConstValue " << "value0=" << value0;
    return out;
  }
};

struct GammaPrior: public IntRegPrior
{
  GammaPrior(double s = 0.2,
             double r = 0.4)
    : shape(s),
      rate(r) {}

  double shape;	// eta0
  double rate;	// gamma0

  std::ostream& print(std::ostream& out) const
  {
    out	<< "Gamma shape=" << shape << " rate=" << rate;
    return out;
  }
};

struct GammaProcessPrior: public IntRegPrior
{
  GammaProcessPrior(double m = 0.1,
                    double c = 1)
    : mean(m),
      control(c) {}

  double mean;	//
  double control;	//

  std::ostream& print(std::ostream& out) const
  {
    out	<< "GammaProcess mean=" << mean << " control=" << control;
    return out;
  }
};

struct InvGammaPrior: public IntRegPrior
{
  InvGammaPrior(double sh = 1.0,
                double sc = 1.0)
    : shape(sh),
      scale(sc) {}

  double shape;	// shape
  double scale;	// scale

  std::ostream& print(std::ostream& out) const
  {
    out	<< "InvGamma shape=" << shape << " scale=" << scale;
    return out;
  }
};

struct NormalPrior: public IntRegPrior
{
  NormalPrior(double m = 0.0,
              double s = 1.0)
    : mean(m),
      sd(s) {}

  double mean;
  double sd;

  std::ostream& print(std::ostream& out) const
  {
    out	<< "Normal mean=" << mean << " sd=" << sd;
    return out;
  }
};

struct NormalProcessPrior: public IntRegPrior
{
  NormalProcessPrior(double s = 1.0)
    : sd(s) {}

  double sd;

  std::ostream& print(std::ostream& out) const
  {
    out	<< "NormalProcess sd=" << sd;
    return out;
  }
};

struct NormalInvGammaProcessPrior: public IntRegPrior
{
  NormalInvGammaProcessPrior(double sh = 2.0,
                             double sc = 1.0)
    : shape(sh),
      scale(sc) {}

  double shape;
  double scale;

  std::ostream& print(std::ostream& out) const
  {
    out	<< "NormalInvGammaProcess shape=" << shape << " scale=" << scale;
    return out;
  }
};
}

/*------------------------------------------------------------------*/
/* Definition of prior class-template                               */
/*------------------------------------------------------------------*/
namespace ir {
template <typename BasePrior, typename CoefPrior>
struct CoxPrior
{
  CoxPrior(const BasePrior& bp,
           const CoefPrior& cp)
    : base_prior(bp),
      coef_prior(cp) {}

  BasePrior base_prior;
  CoefPrior coef_prior;

  typedef BasePrior BasePrior_type;
  typedef CoefPrior CoefPrior_type;
};

template <typename BasePrior, typename CoefPrior, typename ThetaPrior>
struct GORHPrior
{
  GORHPrior(const BasePrior& bp,
            const CoefPrior& cp,
            const ThetaPrior& tp)
    : base_prior(bp),
      coef_prior(cp),
      theta_prior(tp) {}

  BasePrior base_prior;
  CoefPrior coef_prior;
  ThetaPrior theta_prior;

  typedef BasePrior BasePrior_type;
  typedef CoefPrior CoefPrior_type;
  typedef ThetaPrior ThetaPrior_type;
};
}
#endif // INTREG_PRIOR_H_
