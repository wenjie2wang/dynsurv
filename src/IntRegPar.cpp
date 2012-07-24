//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#include "IntRegPar.h"

/*------------------------------------------------------------------*/
/* Definition of model parameter struct                              */
/*------------------------------------------------------------------*/
namespace ir {
/* Baseline hazard */
std::ostream& BaseHazPar::print(std::ostream& out) const
{
  out << "lambda = " << lambda << std::endl;
  return out;
}

std::ofstream& BaseHazPar::output(std::ofstream& out) const
{
  for (Size k = 0; k < lambda.size(); ++k)
    out << lambda(k) << ' ';

  return out;
}

/* Time Indep Cox coef */
std::ostream& TimeIndepCoxPar::print(std::ostream& out) const
{
  BaseHazPar::print(out);
  out	<< "beta = " << beta << std::endl;
  return out;
}

std::ofstream& TimeIndepCoxPar::output(std::ofstream& out) const
{
  BaseHazPar::output(out);
  for (Size j = 0; j < beta.size(); ++j)
    out << beta(j) << ' ';

  return out;
}

/* Time Indep GORH coef */
std::ostream& TimeIndepGORHPar::print(std::ostream& out) const
{
  TimeIndepCoxPar::print(out);
  out	<< "theta = " << theta << std::endl;
  return out;
}

std::ofstream& TimeIndepGORHPar::output(std::ofstream& out) const
{
  TimeIndepCoxPar::output(out);
  out << theta << ' ';
  return out;
}


/* Time varying Cox coef */
std::ostream& TimeVaryingCoxPar::print(std::ostream& out) const
{
  BaseHazPar::print(out);
  out	<< "beta = " << beta << '\n'
      << "nu = " << nu << std::endl;
  return out;
}

std::ofstream& TimeVaryingCoxPar::output(std::ofstream& out) const
{
  BaseHazPar::output(out);
  for (Size j = 0; j < beta.size2(); ++j)
    for (Size k = 0; k < beta.size1(); ++k)
      out << beta(k, j) << ' ';

  for (Size j = 0; j < nu.size(); ++j)
    out << nu(j) << ' ';

  return out;
}

/* Time varying GORH coef */
std::ostream& TimeVaryingGORHPar::print(std::ostream& out) const
{
  TimeVaryingCoxPar::print(out);
  out	<< "theta = " << theta << std::endl;
  return out;
}

std::ofstream& TimeVaryingGORHPar::output(std::ofstream& out) const
{
  TimeVaryingCoxPar::output(out);
  out << theta << ' ';
  return out;
}

/* Time varying Cox coef with jumps */
std::ostream& DynamicCoxPar::print(std::ostream& out) const
{
  TimeVaryingCoxPar::print(out);
  out	<< "jump = " << jump << std::endl;
  return out;
}

std::ofstream& DynamicCoxPar::output(std::ofstream& out) const
{
  TimeVaryingCoxPar::output(out);
  for (Size j = 0; j < jump.size2(); ++j)
    for (Size k = 0; k < jump.size1(); ++k)
      out << jump(k, j) << ' ';

  return out;
}

/* Time varying GORH coef with jumps */
std::ostream& DynamicGORHPar::print(std::ostream& out) const
{
  DynamicCoxPar::print(out);
  out	<< "theta = " << theta << std::endl;
  return out;
}

std::ofstream& DynamicGORHPar::output(std::ofstream& out) const
{
  DynamicCoxPar::output(out);
  out << theta << ' ';
  return out;
}
}
