/*
 * complexfun.h
 *
 *  Created on: Jan 28, 2009
 *      Author: eh7
 */

#ifndef COMPLEXFUN_H_
#define COMPLEXFUN_H_

#include <mra/mra.h>

namespace madness {

//***************************************************************************
double abs(double x) {return x;}
//***************************************************************************

//***************************************************************************
double real(double x) {return x;}
//***************************************************************************

//***************************************************************************
double imag(double x) {return 0.0;}
//***************************************************************************

//***************************************************************************
template <typename Q>
void tensor_real2complex(const Tensor<Q>& r, Tensor< std::complex<Q> >& c)
{
  BINARY_OPTIMIZED_ITERATOR(Q, r, std::complex<Q>, c, *_p1 = std::complex<Q>(0.0,*_p0););
}
//***************************************************************************

//***************************************************************************
template <typename Q>
Tensor<Q> tensor_real(const Tensor< std::complex<Q> >& c)
{
  Tensor<Q> r(c.dim[0], c.dim[1]);
  BINARY_OPTIMIZED_ITERATOR(Q, r, std::complex<Q>, c, *_p0 = real(*_p1););
  return r;
}
//***************************************************************************

//***************************************************************************
template <typename Q>
Tensor<Q> tensor_imag(const Tensor< std::complex<Q> >& c)
{
  Tensor<Q> r(c.dim[0], c.dim[1]);
  BINARY_OPTIMIZED_ITERATOR(Q, r, std::complex<Q>, c, *_p0 = imag(*_p1););
  return r;
}
//***************************************************************************

//***************************************************************************
template <typename Q>
Tensor<Q> tensor_abs(const Tensor< std::complex<Q> >& c)
{
  Tensor<Q> r(c.dim[0], c.dim[1]);
  BINARY_OPTIMIZED_ITERATOR(Q, r, std::complex<Q>, c, *_p0 = abs(*_p1););
  return r;
}
//***************************************************************************

//***************************************************************************
template <typename Q, int NDIM>
struct abs_square_op
{
  typedef typename TensorTypeData<Q>::scalar_type resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim, t.dim);
    BINARY_OPTIMIZED_ITERATOR(Q, t, resultT, result, resultT d = abs(*_p0); *_p1 = d*d);
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};
//***************************************************************************

//***************************************************************************
template<typename Q, int NDIM>
Function<typename TensorTypeData<Q>::scalar_type,NDIM> abs_square(const Function<Q,NDIM>& func)
{
  return unary_op(func, abs_square_op<Q,NDIM>());
}
//***************************************************************************

}
#endif /* COMPLEXFUN_H_ */
