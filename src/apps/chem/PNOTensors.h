/*
 * PNOmadness::Tensors.h
 *
 *  Created on: Oct 22, 2018
 *      Author: kottmanj
 */

#ifndef PAPER_CODE_PNOTENSORS_H_
#define PAPER_CODE_PNOTENSORS_H_

#include <stdio.h>
#include <valarray>
#include "madness.h"


namespace PNOTensors {
  /// index in packed lower triangle matrix, row >= col
  template <typename Int>
    size_t tridx(Int row, Int col) {
    return (row < col) ? (row + col * (col + 1) / 2)
      : (col + row * (row + 1) / 2);
  }
  
  /// # of elements in lower triangle matrix
  template <typename Int>
    size_t ntri(Int n) {
    return n * (n + 1) / 2;
  }
  
  /// stores pairs of pairs that share first index (i.e. {ij} and {ik}, where i >= j) as a 3-index tensor:
  /// case 1 -- i >= j, i >= k: only need j >= k, or ij >= ik; store as {{ij},k}
  /// case 2 -- i >= j, k > i: can also be found at Tensor_IJ_JK(k,i,j)
  template <typename T>
    class Tensor_IJ_IK {
  public:
  Tensor_IJ_IK(size_t n) : n_(n), data_(ntri(n) * n) {
      reset();
    }
    ~Tensor_IJ_IK() = default;
    
    std::tuple<size_t,bool> ijk(size_t i, size_t j, size_t k) const {
      assert(i < size_t(n_));
      assert(j < size_t(n_));
      assert(k < size_t(n_));
      assert(i >= j);
      if (i >= k) {
	return j >= k ? std::make_tuple(tridx(i,j) * n_ + k,false) : std::make_tuple(tridx(i,k) * n_ + j, true);
      }
      else  // i < k
	return std::make_tuple(tridx(i,j) * n_ + k,false);
    }
    
    bool is_unique(size_t i, size_t j, size_t k) const {
      assert(i < size_t(n_));
      assert(j < size_t(n_));
      assert(k < size_t(n_));
      assert(i >= j);
      if (i < k)
	return true;
      else
	return j >= k;
    }
    
    bool is_initialized(size_t i, size_t j, size_t k) const {
      assert(i < size_t(n_));
      assert(j < size_t(n_));
      assert(k < size_t(n_));
      assert(i >= j);
      size_t ijk;
      std::tie(ijk,std::ignore) = this->ijk(i,j,k);
      return data_[ijk].size() != 0;
    }
    madness::Tensor<T> get(size_t i, size_t j, size_t k) const {
      assert(i < size_t(n_));
      assert(j < size_t(n_));
      assert(k < size_t(n_));
      assert(i >= j);
      size_t ijk; bool swap;
      std::tie(ijk, swap) = this->ijk(i,j,k);
      return swap ? data_[ijk].swapdim(0,1) : data_[ijk];
    }
    void set(size_t i, size_t j, size_t k, const madness::Tensor<T>& t) {
      assert(i < size_t(n_));
      assert(j < size_t(n_));
      assert(k < size_t(n_));
      assert(i >= j);
      size_t ijk; bool swap;
      std::tie(ijk, swap) = this->ijk(i,j,k);
      data_[ijk] = swap ? t.swapdim(0,1) : t;
    }
    
    void reset() {
      std::fill(std::begin(data_), std::end(data_), madness::Tensor<T>());
    }
    
  private:
    int n_;
    std::valarray<madness::Tensor<T>> data_;
  };
  
  /// stores pairs of pairs that share second index (i.e. {ij} and {kj}, where i >= j) as a 3-index tensor:
  /// case 1 -- i >= j, j >= k: store as {{ij},k}
  /// case 2 -- i >= j, k > j: only need i >= k, or ij >= kj; store as {{ij},k}
  template <typename T>
    class Tensor_IJ_KJ {
  public:
  Tensor_IJ_KJ(size_t n) : n_(n), data_(ntri(n) * n) {
      reset();
    }
    ~Tensor_IJ_KJ() = default;
    
    std::tuple<size_t,bool> ijk(size_t i, size_t j, size_t k) const {
      assert(i < size_t(n_));
      assert(j < size_t(n_));
      assert(k < size_t(n_));
      assert(i >= j);
      if (k > j)
	return i >= k ? std::make_tuple(tridx(i,j) * n_ + k,false) : std::make_tuple(tridx(k,j) * n_ + i, true);
      else
	return std::make_tuple(tridx(i,j) * n_ + k,false);
    }
    bool is_unique(size_t i, size_t j, size_t k) const {
      assert(i < size_t(n_));
      assert(j < size_t(n_));
      assert(k < size_t(n_));
      assert(i >= j);
      if (k <= j)
	return true;
      else
	return (i >= k);
    }
    bool is_initialized(size_t i, size_t j, size_t k) const {
      assert(i < size_t(n_));
      assert(j < size_t(n_));
      assert(k < size_t(n_));
      assert(i >= j);
      
      size_t ijk;
      std::tie(ijk, std::ignore) = this->ijk(i, j, k);
      return data_[ijk].size() != 0;
    }
    madness::Tensor<T> get(size_t i, size_t j, size_t k) const {
      assert(i < size_t(n_));
      assert(j < size_t(n_));
      assert(k < size_t(n_));
      assert(i >= j);
      size_t ijk;
      bool swap;
      std::tie(ijk, swap) = this->ijk(i,j,k);
      return swap ? data_[ijk].swapdim(0,1) : data_[ijk];
    }
    void set(size_t i, size_t j, size_t k, const madness::Tensor<T>& t) {
      assert(i < size_t(n_));
      assert(j < size_t(n_));
      assert(k < size_t(n_));
      assert(i >= j);
      size_t ijk;
      bool swap;
      std::tie(ijk, swap) = this->ijk(i,j,k);
      data_[ijk] = swap ? t.swapdim(0,1) : t;
    }
    
    void reset() {
      std::fill(std::begin(data_), std::end(data_), madness::Tensor<T>());
    }
    
  private:
    int n_;
    std::valarray<madness::Tensor<T>> data_;
  };
  
}  // was anonymous namespace --- now PNOTensors

#endif /* PAPER_CODE_PNOTENSORS_H_ */
