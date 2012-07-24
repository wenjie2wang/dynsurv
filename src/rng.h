//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once
#ifndef RNG_H_
#define RNG_H_

/* Rmath */
#ifdef R_PACKAGE
extern "C" {
#include <R.h>
}
extern "C" {
  namespace rmath {
#include <Rmath.h>
  }
}
# ifdef beta
#  undef beta
# endif
#else
# define MATHLIB_STANDALONE
extern "C" {
  namespace rmath {
# include <Rmath.h>
  }
}
#endif

/* ARMS */
#include <math.h>
#include <stdlib.h>
#ifdef R_PACKAGE
extern "C" {
# include "arms.h"
}
#else
extern "C" {
# include <arms.h>
}
#endif

#endif // RNG_H_
