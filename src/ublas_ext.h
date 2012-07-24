//
//  Copyright (c) 2010
//  Xiaojing Wang
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//#pragma once

#ifndef UBLAS_EXT_H_
#define UBLAS_EXT_H_

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

// Functional
namespace boost {
namespace numeric {
namespace ublas {
template<class T>
struct scalar_exp:
  public scalar_unary_functor<T> {
  typedef typename scalar_unary_functor<T>::argument_type argument_type;
  typedef typename scalar_unary_functor<T>::result_type result_type;

  static BOOST_UBLAS_INLINE
  result_type apply (argument_type t) {
    return std::exp(t);
  }
};

template<class T>
struct scalar_log:
  public scalar_unary_functor<T> {
  typedef typename scalar_unary_functor<T>::argument_type argument_type;
  typedef typename scalar_unary_functor<T>::result_type result_type;

  static BOOST_UBLAS_INLINE
  result_type apply (argument_type t) {
    return std::log(t);
  }
};
}
}
}

// Expression
namespace boost {
namespace numeric {
namespace ublas {

/* Matrix element wise exponential */
template<class E>
BOOST_UBLAS_INLINE
typename vector_unary_traits<E, scalar_exp<typename E::value_type> >::result_type
exp (const vector_expression<E> &e) {
  typedef typename vector_unary_traits<E, scalar_exp<typename E::value_type> >::expression_type expression_type;
  return expression_type (e ());
}

/* Matrix element wise log */
template<class E>
BOOST_UBLAS_INLINE
typename vector_unary_traits<E, scalar_log<typename E::value_type> >::result_type
log (const vector_expression<E> &e) {
  typedef typename vector_unary_traits<E, scalar_log<typename E::value_type> >::expression_type expression_type;
  return expression_type (e ());
}

/* Vector element wise exponential */
template<class E>
BOOST_UBLAS_INLINE
typename matrix_unary1_traits<E, scalar_exp<typename E::value_type> >::result_type
exp (const matrix_expression<E> &e) {
  typedef typename matrix_unary1_traits<E, scalar_exp<typename E::value_type> >::expression_type expression_type;
  return expression_type (e ());
}

/* Vector element wise log */
template<class E>
BOOST_UBLAS_INLINE
typename matrix_unary1_traits<E, scalar_log<typename E::value_type> >::result_type
log (const matrix_expression<E> &e) {
  typedef typename matrix_unary1_traits<E, scalar_log<typename E::value_type> >::expression_type expression_type;
  return expression_type (e ());
}

}
}
}

namespace boost {
namespace numeric {
namespace ublas {
/* Column sum */
template<class T>
BOOST_UBLAS_INLINE
vector<T> col_sum (const matrix<T> &m) {
  return prod(vector<int>(m.size1(), 1), m);
}

/* Column mean */
template<class T>
BOOST_UBLAS_INLINE
vector<double> col_mean (const matrix<T> &m) {
  return prod(vector<double>(m.size1(), 1.0), m) / m.size1();
}

/* Row sum */
template<class T>
BOOST_UBLAS_INLINE
vector<T> row_sum (const matrix<T> &m) {
  return prod(m, vector<int>(m.size2(), 1));
}

/* Row mean */
template<class T>
BOOST_UBLAS_INLINE
vector<double> row_mean (const matrix<T> &m) {
  return prod(m, vector<double>(m.size2(), 1.0)) / m.size2();
}
}
}
}

#endif // UBLAS_EXT_H_
