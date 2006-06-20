#ifndef TENSOR_LAPACK_H
#define TENSOR_LAPACK_H

#include <tensor/tensor.h>

/// \file tensor_lapack.h
/// \brief Prototypes for a partial interface from Tensor to LAPACK

namespace madness {

    template <typename T>
    void svd(const Tensor<T>& a, Tensor<T>* U,
             Tensor< typename Tensor<T>::scalar_type >* s, Tensor<T>* VT);

    template <typename T>
    void gesv(const Tensor<T>& a, const Tensor<T>& b, Tensor<T>* x);

    template <typename T>
    void gelss(const Tensor<T>& a, const Tensor<T>& b, double rcond,
               Tensor<T>* x, Tensor< typename Tensor<T>::scalar_type >* s, long *rank);
    template <typename T>
    void syev(const Tensor<T>& A,
              Tensor<T>* V, Tensor< typename Tensor<T>::scalar_type >* e);

    template <typename T>
    void sygv(const Tensor<T>& A, const Tensor<T>& B, int itype,
              Tensor<T>* V, Tensor< typename Tensor<T>::scalar_type >* e);

    bool test_tensor_lapack();

}

#endif
