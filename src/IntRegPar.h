//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef INTREG_PAR_H_
#define INTREG_PAR_H_

#include "ublas.h"

#include <iostream>
#include <fstream>
#include <vector>

/*------------------------------------------------------------------*/
/* Definition of model parameter struct                              */
/*------------------------------------------------------------------*/
namespace ir {
/* Baseline hazard */
struct BaseHazPar
{
  BaseHazPar(ublas::vector<double> l)
    : lambda(l) {}

  ublas::vector<double> lambda;

  virtual std::ostream& print(std::ostream& out) const;
  virtual std::ofstream& output(std::ofstream& out) const;
};

/* Time Indep Cox coef */
struct TimeIndepCoxPar: public BaseHazPar
{
  TimeIndepCoxPar(ublas::vector<double> l,
                  ublas::vector<double> b)
    : BaseHazPar(l),
      beta(b) {}

  ublas::vector<double> beta;

  virtual std::ostream& print(std::ostream& out) const;
  virtual std::ofstream& output(std::ofstream& out) const;
};

/* Time Indep CORH coef */
struct TimeIndepGORHPar: public TimeIndepCoxPar
{
  TimeIndepGORHPar(ublas::vector<double> l,
                   ublas::vector<double> b,
                   double t)
    : TimeIndepCoxPar(l, b),
      theta(t) {}

  double theta;

  virtual std::ostream& print(std::ostream& out) const;
  virtual std::ofstream& output(std::ofstream& out) const;
};

/* Time varying Cox coef */
struct TimeVaryingCoxPar: public BaseHazPar
{
  TimeVaryingCoxPar(ublas::vector<double> l,
                    ublas::matrix<double> b,
                    ublas::vector<double> n)
    : BaseHazPar(l),
      beta(b),
      nu(n) {}

  ublas::matrix<double> beta;
  ublas::vector<double> nu;

  virtual std::ostream& print(std::ostream& out) const;
  virtual std::ofstream& output(std::ofstream& out) const;
};

/* Time varying CORH coef */
struct TimeVaryingGORHPar: public TimeVaryingCoxPar
{
  TimeVaryingGORHPar(ublas::vector<double> l,
                     ublas::matrix<double> b,
                     ublas::vector<double> n,
                     double t)
    : TimeVaryingCoxPar(l, b, n),
      theta(t) {}

  double theta;

  virtual std::ostream& print(std::ostream& out) const;
  virtual std::ofstream& output(std::ofstream& out) const;
};

/* Time varying Cox coef with jumps */
struct DynamicCoxPar: public TimeVaryingCoxPar
{
  DynamicCoxPar(ublas::vector<double> l,
                ublas::matrix<double> b,
                ublas::vector<double> n,
                ublas::matrix<int> jp)
    : TimeVaryingCoxPar(l, b, n),
      jump(jp) {}

  ublas::matrix<int> jump;

  virtual std::ostream& print(std::ostream& out) const;
  virtual std::ofstream& output(std::ofstream& out) const;
};

/* Time varying CORH coef with jumps*/
struct DynamicGORHPar: public DynamicCoxPar
{
  DynamicGORHPar(ublas::vector<double> l,
                 ublas::matrix<double> b,
                 ublas::vector<double> n,
                 ublas::matrix<int> jp,
                 double t)
    : DynamicCoxPar(l, b, n, jp),
      theta(t) {}

  double theta;

  virtual std::ostream& print(std::ostream& out) const;
  virtual std::ofstream& output(std::ofstream& out) const;
};
}

/*------------------------------------------------------------------*/
/* Template function, mean of vector of model parameters            */
/*------------------------------------------------------------------*/
namespace ir {
template <typename P>
inline P mean(const std::vector<P>& vp);

/* P = BaseHazPar */
template <>
inline BaseHazPar mean(const std::vector<BaseHazPar>& vp)
{
  Size n = vp.size();
  BaseHazPar par(vp[0]);
  for (Size i = 1; i < n; ++i)
    par.lambda += vp[i].lambda;

  par.lambda /= n;
  return par;
}

/* P = TimeIndepCoxPar */
template <>
inline TimeIndepCoxPar mean(const std::vector<TimeIndepCoxPar>& vp)
{
  Size n = vp.size();
  TimeIndepCoxPar par(vp[0]);
  for (Size i = 1; i < n; ++i) {
    par.lambda += vp[i].lambda;
    par.beta += vp[i].beta;
  }

  par.lambda /= n;
  par.beta /= n;
  return par;
}

/* P = TimeIndepGORHPar */
template <>
inline TimeIndepGORHPar mean(const std::vector<TimeIndepGORHPar>& vp)
{
  Size n = vp.size();
  TimeIndepGORHPar par(vp[0]);
  for (Size i = 1; i < n; ++i) {
    par.lambda += vp[i].lambda;
    par.beta += vp[i].beta;
    par.theta += vp[i].theta;
  }

  par.lambda /= n;
  par.beta /= n;
  par.theta /= n;
  return par;
}

/* P = TimeVaryingCoxPar */
template <>
inline TimeVaryingCoxPar mean(const std::vector<TimeVaryingCoxPar>& vp)
{
  Size n = vp.size();
  TimeVaryingCoxPar par(vp[0]);
  for (Size i = 1; i < n; ++i) {
    par.lambda += vp[i].lambda;
    par.beta += vp[i].beta;
    par.nu += vp[i].nu;
  }

  par.lambda /= n;
  par.beta /= n;
  par.nu /= n;
  return par;
}

/* P = TimeVaryingGORHPar */
template <>
inline TimeVaryingGORHPar mean(const std::vector<TimeVaryingGORHPar>& vp)
{
  Size n = vp.size();
  TimeVaryingGORHPar par(vp[0]);
  for (Size i = 1; i < n; ++i) {
    par.lambda += vp[i].lambda;
    par.beta += vp[i].beta;
    par.nu += vp[i].nu;
    par.theta += vp[i].theta;
  }

  par.lambda /= n;
  par.beta /= n;
  par.nu /= n;
  par.theta /= n;
  return par;
}

/* P = DynamicCoxPar */
template <>
inline DynamicCoxPar mean(const std::vector<DynamicCoxPar>& vp)
{
  Size n = vp.size();
  DynamicCoxPar par(vp[0]);
  for (Size i = 1; i < n; ++i) {
    par.lambda += vp[i].lambda;
    par.beta += vp[i].beta;
    par.nu += vp[i].nu;
    par.jump += vp[i].jump;
  }

  par.lambda /= n;
  par.beta /= n;
  par.nu /= n;
  return par;
}

/* P = DynamicGORHPar */
template <>
inline DynamicGORHPar mean(const std::vector<DynamicGORHPar>& vp)
{
  Size n = vp.size();
  DynamicGORHPar par(vp[0]);
  for (Size i = 1; i < n; ++i) {
    par.lambda += vp[i].lambda;
    par.beta += vp[i].beta;
    par.nu += vp[i].nu;
    par.jump += vp[i].jump;
    par.theta += vp[i].theta;
  }

  par.lambda /= n;
  par.beta /= n;
  par.nu /= n;
  par.theta /= n;
  return par;
}
}
#endif // INTREG_PAR_H_
