//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once
#ifndef INTREG_DATA_H_
#define INTREG_DATA_H_

#include "ublas.h"

/*------------------------------------------------------------------*/
/* Definition of interval censored data class                       */
/*------------------------------------------------------------------*/
namespace ir {
class IntRegData
{
public:
  IntRegData(const ublas::matrix<double>& LRX,
             const ublas::vector<double>& grid)
    : X_  (ublas::subrange(LRX, 0, LRX.size1(), 2, LRX.size2())),
      oneX_ (ublas::subrange(LRX, 0, LRX.size1(), 1, LRX.size2())),
      left_ (ublas::column(LRX, 0)),
      right_(ublas::column(LRX, 1)),
      grid_ (grid),
      N_    (LRX.size1()),
      K_    (grid.size()),
      nBeta_(LRX.size2() - 2)
  {
    ublas::column(oneX_, 0) = ublas::vector<double>(N_, 1.0);
  }

  /* Public accessor */
  const ublas::matrix<double>& X() const {
    return X_;
  }
  const ublas::matrix<double>& oneX() const {
    return oneX_;
  }
  const ublas::vector<double>& left() const {
    return left_;
  }
  const ublas::vector<double>& right() const {
    return right_;
  }
  const ublas::vector<double>& grid() const {
    return grid_;
  }
  Size N() {
    return N_;
  }
  Size K() {
    return K_;
  }
  Size nBeta() {
    return nBeta_;
  }

private:
  ublas::matrix<double> X_;
  ublas::matrix<double> oneX_;
  ublas::vector<double> left_;
  ublas::vector<double> right_;
  ublas::vector<double> grid_;
  const Size N_;
  const Size K_;
  const Size nBeta_;
};
}
#endif // INTREG_DATA_H_
