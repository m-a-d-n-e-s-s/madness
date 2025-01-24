/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#ifndef MADNESS_WORLD_VECTOR_H__INCLUDED
#define MADNESS_WORLD_VECTOR_H__INCLUDED

/**
 \file vector.h
 \brief Implement the \c madness:Vector class, an extension of \c std::array
    that supports some mathematical operations.
 \ingroup containers
*/

#include <madness/madness_config.h>
#include <madness/world/madness_exception.h>
#include <madness/world/worldhash.h>
#include <madness/world/array_addons.h>
#include <madness/world/archive.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>

namespace madness {

    /// A simple, fixed dimension vector.

    /// This class eliminates memory allocation cost, is just POD (it can be
    /// copied easily and allocated on the stack), and the known
    /// dimension permits aggressive compiler optimizations.
    ///
    /// Provides additional mathematical and I/O operations.
    /// \tparam T The type of data stored in the vector.
    /// \tparam N The size of the vector.
    template <typename T, std::size_t N>
    class Vector {
    public:
        using arrayT = std::array<T,N>; ///< The underlying array type.

        template <typename Q,std::size_t M>
        friend class Vector;

    private:
        arrayT data_; ///< The underlying array.

    public:
        // type defs... these are just wrappers to the underlying array types
        using value_type = typename arrayT::value_type; ///< The data value type.
        using iterator = typename arrayT::iterator; ///< Iterator type.
        using const_iterator = typename arrayT::const_iterator; ///< Const iterator type.
        using reverse_iterator = typename arrayT::reverse_iterator; ///< Reverse iterator type.
        using const_reverse_iterator = typename arrayT::const_reverse_iterator; ///< Const reverse iterator type.
        using reference = typename arrayT::reference; ///< Reference type.
        using const_reference = typename arrayT::const_reference; ///< Const reference type.
        using size_type = typename arrayT::size_type; ///< Size type.
        using difference_type = typename arrayT::difference_type; ///< Difference type.

        /// The size of the \c Vector.
        static inline constexpr size_type static_size = N;

        /// Default constructor; does not initialize vector contents.
        constexpr Vector() = default;

        /// Initialize all elements to value \c t.

        /// \tparam Q The type of \c t.
        /// \param[in] t The value used to initialized the \c Vector.
        template <typename Q>
        constexpr explicit Vector(Q t) : data_{} {
            fill(t);
        }

        /// Construct from a C-style array of the same dimension.

        /// \tparam Q The type of data in \c t.
        /// \param[in] t The C-style array.
        template <typename Q>
        constexpr explicit Vector(const Q (&t)[N]) {
#if __cplusplus >= 202002L
          std::copy(t, t+N, data_.begin());
#else
          for(std::size_t i=0; i!=N; ++i) data_[i] = t[i];
#endif
        }

        /// Construct from an STL vector of equal or greater length.

        /// \tparam Q Type of data stored in the \c std::vector.
        /// \tparam A Allocator type for the \c std::vector.
        /// \param[in] t The \c std::vector.
        template <typename Q, typename A>
        constexpr explicit Vector(const std::vector<Q, A>& t) {
            operator=(t);
        }

        /// Construct from a \c std::array of equal length.

        /// \tparam Q Type of data stored in the original \c std::array.
        /// \param[in] t The \c std::array.
        template <typename Q>
        constexpr explicit Vector(const std::array<Q, N>& t) {
            data_ = t;
        }

        /// Copy constructor is deep (because \c Vector is POD).

        /// \param[in] other The \c Vector to copy.
        constexpr Vector(const Vector<T,N>& other) {
MADNESS_PRAGMA_GCC(diagnostic push)
MADNESS_PRAGMA_GCC(diagnostic ignored "-Wuninitialized")
MADNESS_PRAGMA_GCC(diagnostic ignored "-Wmaybe-uninitialized")
            data_ = other.data_;
MADNESS_PRAGMA_GCC(diagnostic pop)
        }

        /// Copy constructor is deep (because \c Vector is POD).

        /// \tparam Q Type of the \c Vector to copy.
        /// \param[in] other The \c Vector to copy.
        template <typename Q>
        constexpr Vector(const Vector<Q,N>& other) {
            std::copy(other.data_.begin(), other.data_.end(), data_.begin());
        }

        /// List initialization constructor (deep copy because \c Vector is POD).

        /// This constructor allows initialization using, e.g.,
        /// \code
        /// Vector<double, 3> v{ 1.5, 2.4, -1.9 };
        /// \endcode
        /// \throw MadnessException if the list does not contain exactly \c N
        ///    elements.
        /// \param[in] list The initializer list; elements are copied to the
        ///    \c Vector.
        constexpr Vector(std::initializer_list<T> list) :
            data_()
        {
            MADNESS_ASSERT(list.size() == N);
#if __cplusplus >= 202002L
            std::copy(list.begin(), list.end(), data_.begin());
#else
            for(std::size_t i=0; i!=N; ++i) data_[i] = *(list.begin()+i);
#endif
        }

        /// Assignment is deep (because a \c Vector is POD).

        /// \param[in] other The \c Vector to copy.
        /// \return This \c Vector.
        constexpr Vector<T,N>& operator=(const Vector<T,N>& other) {
            data_ = other.data_;
            return *this;
        }

        /// Assignment is deep (because \c Vector is POD).

        /// \tparam Q The type of the \c Vector to copy.
        /// \param[in] other The \c Vector to copy.
        /// \return This \c Vector.
        template <typename Q>
        constexpr Vector<T,N>& operator=(const Vector<Q,N>& other) {
            std::copy(other.data_.begin(), other.data_.end(), data_.begin());
            return *this;
        }

        /// Assignment is deep (because \c Vector is POD).

        /// Make sure the size of \c other is at least \c N.
        /// \tparam Q The type of data in the \c std::vector.
        /// \tparam A The allocator type for the \c std::vector.
        /// \param[in] other The \c std::vector to copy.
        /// \return This \c Vector.
        template <typename Q, typename A>
        constexpr Vector<T,N>& operator=(const std::vector<Q, A>& other) {
            MADNESS_ASSERT(other.size() >= N);
#if __cplusplus >= 202002L
            std::copy(other.begin(), other.begin() + N, data_.begin());
#else
            for(std::size_t i=0; i!=N; ++i) data_[i] = other[i];
#endif
            return *this;
        }

        /// List initialization assignment (deep copy because \c Vector is POD).

        /// This assignment operator allows initialization using, e.g.,
        /// \code
        /// v = { 1.5, 2.4, -1.9 };
        /// \endcode
        /// \throw MadnessException if the list does not contain exactly \c N
        ///    elements.
        /// \param[in] list The initializer list; elements are copied to the
        ///    \c Vector.
        constexpr Vector<T,N>& operator=(std::initializer_list<T> list) {
            MADNESS_ASSERT(list.size() == N);
#if __cplusplus >= 202002L
            std::copy(list.begin(), list.end(), data_.begin());
#else
            for(std::size_t i=0; i!=N; ++i) data_[i] = *(list.begin() + i);
#endif
            return *this;
        }

        /// Fill from a scalar value.

        /// \param[in] t The scalar to use for filling.
        /// \return This \c Vector.
        constexpr Vector<T,N>& operator=(const T& t) {
            fill(t);
            return *this;
        }

        /// Type conversion to a \c std::array.

        /// \return The underlying \c std::array.
        constexpr operator std::array<T,N> () { return data_; }

        // iterator support
        /// Iterator starting at the first element.

        /// \return Iterator to the starting element.
        constexpr iterator begin() { return data_.begin(); }

        /// Const iterator starting at the first element.

        /// \return Const iterator to the starting element.
        constexpr const_iterator begin() const { return data_.begin(); }

        /// Iterator to the end (past the last element).

        /// \return Iterator to the end.
        constexpr iterator end() { return data_.end(); }

        /// Const iterator to the end (past the last element).

        /// \return Const iterator to the end.
        constexpr const_iterator end() const { return data_.end(); }

        // reverse iterator support
        /// Reverse iterator starting at the last element.

        /// \return Reverse iterator to the last element.
        constexpr reverse_iterator rbegin() { return data_.rbegin(); }

        /// Const reverse iterator starting at the last element.

        /// \return Const reverse iterator to the last element.
        constexpr const_reverse_iterator rbegin() const { return data_.rbegin(); }

        /// Reverse iterator to the beginning (before the first element).

        /// \return Reverse iterator to the beginning.
        constexpr reverse_iterator rend() { return data_.rend(); }

        /// Const reverse iterator to the beginning (before the first element).

        /// \return Const reverse iterator to the beginning.
        constexpr const_reverse_iterator rend() const { return data_.rend(); }

        // capacity
        /// Accessor for the number of elements in the \c Vector.

        /// \return The number of elements.
        constexpr size_type size() const { return data_.size(); }

        /// Check if the \c Vector is empty.

        /// \return True if the \c Vector is empty; false otherwise. This
        ///    should be false unless `N == 0`.
        constexpr  bool empty() const { return data_.empty(); }

        /// Get the maximum size of the \c Vector.

        /// \return The maximum size, \c N.
        constexpr size_type max_size() const { return data_.max_size(); }

        // element access
        /// Access element \c i of the \c Vector.

        /// Bounds checking is not performed.
        /// \param[in] i The index.
        /// \return A reference to element \c i.
        constexpr reference operator[](size_type i) { return data_[i]; }

        /// Access element \c i of the \c Vector.

        /// Bounds checking is not performed.
        /// \param[in] i The index.
        /// \return A const reference to element \c i.
        constexpr const_reference operator[](size_type i) const { return data_[i]; }

        /// Access element \c i of the \c Vector with bounds checking.

        /// \param[in] i The index.
        /// \return A reference to element \c i.
        constexpr reference at(size_type i) { return data_.at(i); }

        /// Access element \c i of the \c Vector with bounds checking.

        /// \param[in] i The index.
        /// \return A const reference to element \c i.
        constexpr const_reference at(size_type i) const { return data_.at(i); }

        /// Access the first element.

        /// \return A reference to the first element.
        constexpr reference front() { return data_.front(); }

        /// Access the first element.

        /// \return A const reference to the first element.
        constexpr const_reference front() const { return data_.front(); }

        /// Access the last element.

        /// \return A reference to the last element.
        constexpr reference back() { return data_.back(); }

        /// Access the last element.

        /// \return A const reference to the last element.
        constexpr const_reference back() const { return data_.back(); }

        /// Direct access to the underlying array.

        /// \return Pointer to the underlying array.
        constexpr T* data() { return data_.data(); }

        /// Direct access to the underlying array.

        /// \return Const pointer to the underlying array.
        constexpr const T* data() const { return data_.data(); }

        // modifiers
        /// Swap the contents with another \c Vector.

        /// \param[in] other The other vector.
        constexpr void swap(Vector<T, N>& other) { data_.swap(other.data_); }

        /// Fill the \c Vector with the specified value.

        /// \param[in] t The value used to fill the \c Vector.
        constexpr void fill(const T& t) {
#if __cplusplus >= 202002L
          data_.fill(t);
#else
          for (std::size_t i = 0; i < N; i++)
            data_[i] = t;
#endif
        }

        /// In-place, element-wise multiplcation by a scalar.

        /// \tparam Q Type of the scalar.
        /// \param[in] q The scalar.
        /// \return A reference to this for chaining operations.
        /// \todo Do we want a similar division operation?
        template <typename Q>
        constexpr Vector<T,N>& operator*=(Q q) {
            for(size_type i = 0; i < N; ++i)
                data_[i] *= q;
            return *this;
        }

        /// In-place, element-wise addition of another \c Vector.

        /// \tparam Q Type stored in the other \c Vector.
        /// \param[in] q The other \c Vector.
        /// \return A reference to this for chaining operations.
        template <typename Q>
        constexpr Vector<T,N>& operator+=(const Vector<Q,N>& q) {
            for(size_type i = 0; i < N; ++i)
                data_[i] += q[i];
            return *this;
        }

        /// In-place, element-wise subtraction of another \c Vector.

        /// \tparam Q Type stored in the other \c Vector.
        /// \param[in] q The other \c Vector.
        /// \returns A reference to this for chaining operations.
        template <typename Q>
        constexpr Vector<T,N>& operator-=(const Vector<Q,N>& q) {
            for(size_type i = 0; i < N; ++i)
                data_[i] -= q[i];
            return *this;
        }

        /// Calculate the 2-norm of the vector elements.

        /// \return The 2-norm.
        /// \todo Is there a reason this is "normf" and not "norm2"?
        constexpr T normf() const {
        	T d = 0.;
        	for(std::size_t i=0; i<N; ++i)
                d += (data_[i])*(data_[i]);
        	return sqrt(d);
        }

        /// Support for MADNESS serialization.

        /// \tparam Archive The archive type.
        /// \param[in,out] ar The archive.
        template <typename Archive>
        constexpr void serialize(Archive& ar) {
            ar & data_;
        }

        /// Support for MADNESS hashing.

        /// \return The hash.
        constexpr hashT hash() const {
            return hash_value(data_);
        }

        // comparisons
        /// Check if each element is equal to its partner in the other \c Vector.

        /// \param[in] l One \c Vector.
        /// \param[in] r The other \c Vector.
        /// \return True if each element is equal to its partner; false otherwise.
        friend constexpr bool operator==(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ == r.data_;
        }

        /// Check if any element is not equal to its partner in the other \c Vector.

        /// \param[in] l One \c Vector.
        /// \param[in] r The other \c Vector.
        /// \return True if any element is not equal to its partner; false otherwise.
        friend constexpr bool operator!=(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ != r.data_;
        }

        /// Compare \c l and \c r lexicographically.

        /// \param[in] l One \c Vector.
        /// \param[in] r The other \c Vector.
        /// \return True if the contents of \c l are lexicographically less than the contents of \c r; false otherwise.
        friend constexpr bool operator<(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ < r.data_;
        }

        /// Compare \c l and \c r lexicographically.

        /// \param[in] l One \c Vector.
        /// \param[in] r The other \c Vector.
        /// \return True if the contents of \c l are lexicographically greater than the contents of \c r; false otherwise.
        friend constexpr bool operator>(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ > r.data_;
        }

        /// Compare \c l and \c r lexicographically.

        /// \param[in] l One \c Vector.
        /// \param[in] r The other \c Vector.
        /// \return True if the contents of \c l are lexicographically less than or equal to the contents of \c r; false otherwise.
        friend constexpr bool operator<=(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ <= r.data_;
        }

        /// Compare \c l and \c r lexicographically.

        /// \param[in] l One \c Vector.
        /// \param[in] r The other \c Vector.
        /// \return True if the contents of \c l are lexicographically greater than or equal to the contents of \c r; false otherwise.
        friend constexpr bool operator>=(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ >= r.data_;
        }

        /// Output a \c Vector to stream.

        /// \param[in] s The output stream.
        /// \param[in] v The \c Vector to output.
        /// \return The output stream.
        friend std::ostream& operator<<(std::ostream& s, const Vector<T,N>& v) {
            using madness::operators::operator<<;
            s << v.data_;
            return s;
        }
    }; // class Vector

    /// Swap the contents of two `Vector`s.

    /// \tparam T The type of data stored in the `Vector`s.
    /// \tparam N The size of the `Vector`s.
    /// \param[in,out] l One \c Vector.
    /// \param[in,out] r The other \c Vector.
    template <typename T, std::size_t N>
    constexpr void swap(Vector<T,N>& l, Vector<T,N>& r) {
        l.swap(r);
    }


    // Arithmetic operators

    /// Scale a \c Vector.

    /// Multiply each \c Vector element by the scalar value \c r.
    /// \tparam T The left-hand \c Vector element type.
    /// \tparam N The \c Vector size.
    /// \tparam U The right-hand scalar type.
    /// \param[in] l The left-hand \c Vector.
    /// \param[in] r The right-hand scalar value to be multiplied by
    ///    the \c Vector elements.
    /// \return A new \c Vector, \c c, where `c[i]==(l[i]*r)`.
    template <typename T, std::size_t N, typename U>
    constexpr Vector<T,N> operator*(Vector<T,N> l, U r) {
        // coordinate passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] *= r;
        return l;
    }


    /// Scale a \c Vector.

    /// Multiply the scalar value \c l by each \c Vector element.
    /// \tparam T The left-hand scalar type.
    /// \tparam U The right-hand \c Vector element type.
    /// \tparam N The \c Vector size.
    /// \param[in] l The left-hand scalar value to be multiplied by
    ///    the \c Vector elements.
    /// \param[in] r The right-hand \c Vector.
    /// \return A new \c Vector, \c c, where `c[i]==(l*r[i])`.
    template <typename T, typename U, std::size_t N>
    constexpr Vector<T,N> operator*(T l, Vector<U,N> r) {
        // coordinate passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            r[i] *= l;
        return r;
    }

    /// Multiply (element-wise) two `Vector`s.

    /// Do an element-wise multiplication of \c l and \c r and return the
    /// result in a new \c Vector.
    /// \tparam T The left-hand \c Vector element type.
    /// \tparam N The \c Vector size.
    /// \tparam U The right-hand \c Vector element type.
    /// \param[in] l The left-hand \c Vector.
    /// \param[in] r The right-hand \c Vector.
    /// \return A new \c Vector, \c c, where `c[i]==(l[i]*r[i])`.
    template <typename T, std::size_t N, typename U>
    constexpr Vector<T,N> operator*(Vector<T,N> l, const Vector<U,N>& r) {
        // coordinate r passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] *= r[i];
        return l;
    }

    /// Divide (element-wise) two `Vector`s.

    /// Do an element-wise division of \c l by \c r and return the
    /// result in a new \c Vector.
    /// \tparam T The left-hand \c Vector element type.
    /// \tparam N The \c Vector size.
    /// \tparam U The right-hand \c Vector element type.
    /// \param[in] l The left-hand \c Vector.
    /// \param[in] r The right-hand \c Vector.
    /// \return A new \c Vector, \c c, where `c[i]==(l[i]/r[i])`.
    template <typename T, std::size_t N, typename U>
    constexpr Vector<T,N> operator/(Vector<T,N> l, const Vector<U,N>& r) {
      for (std::size_t i = 0; i < N; ++i)
        l[i] /= r[i];
      return l;
    }

    /// Add a scalar to a \c Vector.

    /// Add the scalar value \c r to each \c Vector element.
    /// \tparam T The left-hand \c Vector element type.
    /// \tparam N The \c Vector size.
    /// \tparam U The right-hand scalar type.
    /// \param[in] l The left-hand \c Vector.
    /// \param[in] r The right-hand scalar value to be added to the \c Vector.
    /// \return A new \c Vector, \c c, where `c[i]==(l[i]+r)`.
    template <typename T, std::size_t N, typename U>
    constexpr Vector<T,N> operator+(Vector<T,N> l, U r) {
        // coordinate passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] += r;
        return l;
    }

    /// Add (element-wise) two `Vector`s.

    /// Do an element-wise addition of \c l and \c r and return the
    /// result in a new \c Vector.
    /// \tparam T The left-hand \c Vector element type.
    /// \tparam N The \c Vector size.
    /// \tparam U The right-hand \c Vector element type.
    /// \param[in] l The left-hand \c Vector.
    /// \param[in] r The right-hand \c Vector.
    /// \return A new \c Vector, \c c, where `c[i]==(l[i]+r[i])`.
    template <typename T, std::size_t N, typename U>
    constexpr Vector<T,N> operator+(Vector<T,N> l, const Vector<U,N>& r) {
        // l passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] += r[i];
        return l;
    }

    /// Subtract a scalar from a \c Vector.

    /// Subtract the scalar value \c r from the \c Vector elements `l[i]`.
    /// \tparam T The left-hand \c Vector element type.
    /// \tparam N The \c Vector size.
    /// \tparam U The right-hand scalar type.
    /// \param[in] l The left-hand \c Vector.
    /// \param[in] r The right-hand scalar value to be added to the \c Vector.
    /// \return A new \c Vector, \c c, where `c[i]==(l[i]-r)`.
    template <typename T, std::size_t N, typename U>
    constexpr Vector<T,N> operator-(Vector<T,N> l, U r) {
        // l passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] -= r;
        return l;
    }

    /// Subtract (element-wise) two `Vector`s.

    /// Do an element-wise subtraction of \c l and \c r and return the
    /// result in a new \c Vector.
    /// \tparam T The left-hand \c Vector element type.
    /// \tparam N The \c Vector size.
    /// \tparam U The right-hand \c Vector element type.
    /// \param[in] l The left-hand \c Vector.
    /// \param[in] r The right-hand \c Vector.
    /// \return A new \c Vector, \c c, where `c[i]==(l[i]-r[i])`.
    template <typename T, std::size_t N, typename U>
    constexpr Vector<T,N> operator-(Vector<T,N> l, const Vector<U,N>& r) {
        // l passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] -= r[i];
        return l;
    }


    /// compute the inner product two `Vector`s.

    /// Do an inner product of \c l and \c r and return the
    /// result in a new \c Vector.
    /// \tparam T The \c Vector element type.
    /// \tparam N The \c Vector size.
    /// \param[in] l The left-hand \c Vector.
    /// \param[in] r The right-hand \c Vector.
    /// \return the inner product, where `result==\sum_i(l[i]*r[i])`.
    template <typename T, std::size_t N>
    constexpr T inner(const Vector<T,N>& l, const Vector<T,N>& r) {
        static_assert(std::is_arithmetic_v<T>);  // for complex need to conjugate l
    	T result=0.0;
        for (std::size_t i = 0; i < N; ++i)
            result+=l[i]*r[i];
        return result;
    }

    /// compute the cross product two `Vector`s of dimension 3

    /// Do a cross product of \c l and \c r and return the result in a new \c Vector.
    /// \tparam T The \c Vector element type.
    /// \tparam N The \c Vector size.
    /// \param[in] l The left-hand \c Vector.
    /// \param[in] r The right-hand \c Vector.
    /// \return the cross product
    template <typename T, std::size_t N>
    constexpr typename std::enable_if<N==3, Vector<T,N> >::type
	cross(const Vector<T,N>& l, const Vector<T,N>& r) {
    	Vector<T,N> result;
    	result[0]=l[1]*r[2] - r[1]*l[2];
    	result[1]=l[2]*r[0] - r[2]*l[0];
    	result[2]=l[0]*r[1] - r[0]*l[1];
        return result;
    }


    /// Factory function for creating a \c madness::Vector.

    /// Variadic templates are used to create factories that mimic
    /// \code
    /// inline madness::Vector<T, N> vec(T t1, ..., T tN) {
    ///     std::Vector<T, N> ret;
    ///     ret[0] = t1;
    ///     ...
    ///     ret[N-1] = tN;
    ///     return ret;
    /// }
    /// \endcode
    ///
    /// This function counts the number of arguments passed in through the
    /// argument pack, creates a \c std::array of the appropriate size,
    /// forwards the arguments to the `std::array`'s constructor, and passes
    /// the \c std::array to \c madness::Vector.
    ///
    /// \deprecated This function has been replaced by the list-initialization
    ///    constructor and assignment operator. Rather than
    ///    \code
    ///    Vector<double, 3> v = vec(1.4, 2.5, 4.0);
    ///    \endcode
    ///    use
    ///    \code
    ///    Vector<double, 3> v{ 1.4, 2.5, 4.0 };
    ///    \endcode
    ///    or
    ///    \code
    ///    Vector<double, 3> v = { 1.4, 2.5, 4.0 };
    ///    \endcode
    ///
    /// \note The first argument is separated from the pack to prevent 0-size
    ///    arrays and also so that the caller doesn't have to explicitly
    ///    specify \c T. It is assumed that all arguments are of type \c T or
    ///    are convertible to type \c T.
    ///
    /// \tparam T The data type for the array.
    /// \tparam Ts The argument pack; that is, the list of arguments. The
    ///    size of the resulting \c Vector is directly determined from the
    ///    size of the argument pack.
    /// \param[in] t The first argument.
    /// \param[in] ts The rest of the arguments.
    /// \return The \c Vector with the arguments put into it.
    template<typename T, typename... Ts>
    constexpr inline Vector<T, sizeof...(Ts) + 1> vec(T t, Ts... ts) {
        return Vector<T, sizeof...(Ts) + 1> {
            std::array<T, sizeof...(Ts) + 1>
                {{ t, static_cast<T>(ts)... }}
        };
    }


	/// Construct a unit-`Vector` that has the same direction as \c r.

    /// \tparam T The type of data stored in the \c Vector.
    /// \tparam N The size of the \c Vector.
    /// \param[in] r The \c Vector.
    /// \param[in] eps A precision. If the norm of \c r is less than \c eps,
    ///    the zero vector is returned.
    /// \return The desired unit-`Vector` (unless \c r is numerically the zero
    ///    \c Vector).
    template<typename T, std::size_t N>
    constexpr Vector<T,N> unitvec(const Vector<T,N>& r, const double eps=1.e-6) {
		const double norm=r.normf();
		if(norm < eps)
            return Vector<T,N>(0.0);
		return r * (1.0/norm);
	}

} // namespace madness

#endif // MADNESS_WORLD_VECTOR_H__INCLUDED
