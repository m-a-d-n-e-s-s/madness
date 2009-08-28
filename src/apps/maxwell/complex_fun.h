#ifndef __complex_fun__
#define __complex_fun__

template <typename Q, int NDIM>
struct real_op
{
  typedef typename TensorTypeData<Q>::scalar_type resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim, t.dim);
    BINARY_OPTIMIZED_ITERATOR(Q, t, resultT, result, *_p1 = real(*_p0););
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};

template<typename Q, int NDIM>
Function<typename TensorTypeData<Q>::scalar_type,NDIM> real(const Function<Q,NDIM>& func)
{
  return unary_op_coeffs(func, real_op<Q,NDIM>());
}

inline void ln(const Key<3> &key, Tensor<std::complex<double> > &t) {
	UNARY_OPTIMIZED_ITERATOR(std::complex<double>, t,
		*_p0 = log(*_p0));
}

#endif
