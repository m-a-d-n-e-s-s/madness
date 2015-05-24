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
#include <array>
#include <madness/world/array_addons.h>
#include <vector>
#include <algorithm>
#include <iostream>

namespace madness {

    /// A simple, fixed dimension vector.

    /// This class eliminates memory allocation cost, is just POD (it can be
    /// copied easily and allocated on the stack), and the known
    /// dimension permits aggressive compiler optimizations.
    ///
    /// Provides additional mathematical and I/O operations.
    /// \tparam T The type of data stored in the vector.
    /// \tparam N The size of the vector.
    /// \todo Could this class be derived from std::array (to prevent
    ///    duplication of many functions), rather than simply encapsulating
    ///    an array? Not documenting duplicated functions until this question
    ///    is resolved.
    template <typename T, std::size_t N>
    class Vector {
    public:
        typedef std::array<T,N> arrayT;

    private:
        arrayT data_;

    public:
        // type defs
        typedef typename arrayT::value_type             value_type;
        typedef typename arrayT::iterator               iterator;
        typedef typename arrayT::const_iterator         const_iterator;
        typedef typename arrayT::reverse_iterator       reverse_iterator;
        typedef typename arrayT::const_reverse_iterator const_reverse_iterator;
        typedef typename arrayT::reference              reference;
        typedef typename arrayT::const_reference        const_reference;
        typedef typename arrayT::size_type              size_type;
        typedef typename arrayT::difference_type        difference_type;

        /// The size of the vector.
        static const size_type static_size = N;

        /// Default constructor; does not initialize vector contents.
        Vector() = default;

        /// Initialize all elements to value \c t.

        /// \tparam Q The type of \c t.
        /// \param[in] t The value used to initialized the \c Vector.
        template <typename Q>
        explicit Vector(Q t) {
            fill(t);
        }

        /// Construct from a C-style array of the same dimension.

        /// \tparam Q The type of data in \c t.
        /// \param[in] t The C-style array.
        template <typename Q>
        explicit Vector(const Q (&t)[N]) {
            std::copy(t, t + N, data_.begin());
        }

        /// Construct from an STL vector of equal or greater length.

        /// \tparam Q Type of data stored in the \c std::vector.
        /// \tparam A Allocator type for the \c std::vector.
        /// \param[in] t The \c std::vector.
        template <typename Q, typename A>
        explicit Vector(const std::vector<Q, A>& t) {
            operator=(t);
        }

        /// Construct from a \c std::array of equal length.

        /// \tparam Q Type of data stored in the original \c std::array.
        /// \param[in] t The \c std::array.
        template <typename Q>
        explicit Vector(const std::array<Q, N>& t) {
            data_ = t;
        }

        /// Copy constructor is deep (because \c Vector is POD).

        /// \param[in] other The \c Vector to copy.
        Vector(const Vector<T,N>& other) {
            data_ = other.data_;
        }

        /// Copy constructor is deep (because \c Vector is POD).

        /// \tparam Q Type of the \c Vector to copy.
        /// \param[in] other The \c Vector to copy.
        template <typename Q>
        Vector(const Vector<Q,N>& other) {
            data_ = other.data_;
        }

        /// Assignment is deep (because a \c Vector is POD).

        /// \param[in] other The \c Vector to copy.
        /// \return This \c Vector.
        Vector<T,N>& operator=(const Vector<T,N>& other) {
            data_ = other.data_;
            return *this;
        }

        /// Assignment is deep (because \c Vector is POD).

        /// \tparam Q The type of the \c Vector to copy.
        /// \param[in] other The \c Vector to copy.
        /// \return This \c Vector.
        template <typename Q>
        Vector<T,N>& operator=(const Vector<Q,N>& other) {
            data_ = other.data_;
            return *this;
        }

        /// Assignment is deep (because \c Vector is POD).

        /// Make sure the size of \c other is at least \c N.
        /// \tparam Q The type of data in the \c std::vector.
        /// \tparam A The allocator type for the \c std::vector.
        /// \param[in] other The \c std::vector to copy.
        /// \return This \c Vector.
        template <typename Q, typename A>
        Vector<T,N>& operator=(const std::vector<Q, A>& other) {
            MADNESS_ASSERT(other.size() >= N);
            std::copy(other.begin(), other.begin() + N, data_.begin());
            return *this;
        }

        /// Fill from a scalar value.

        /// \param[in] t The scalar to use for filling.
        /// \return This \c Vector.
        Vector<T,N>& operator=(const T& t) {
            fill(t);
            return *this;
        }

        // type conversion
        operator std::array<T,N> () { return data_; }

         // iterator support
         iterator begin() { return data_.begin(); }
         const_iterator begin() const { return data_.begin(); }
         iterator end() { return data_.end(); }
         const_iterator end() const { return data_.end(); }

         // reverse iterator support
         reverse_iterator rbegin() { return data_.rbegin(); }
         const_reverse_iterator rbegin() const { return data_.rbegin(); }
         reverse_iterator rend() { return data_.rend(); }
         const_reverse_iterator rend() const { return data_.rend(); }

         // capacity
         size_type size() const { return data_.size(); }
         bool empty() const { return data_.empty(); }
         size_type max_size() const { return data_.max_size(); }

         // element access
         reference operator[](size_type i) { return data_[i]; }
         const_reference operator[](size_type i) const { return data_[i]; }
         reference at(size_type i) { return data_.at(i); }
         const_reference at(size_type i) const { return data_.at(i); }
         reference front() { return data_.front(); }
         const_reference front() const { return data_.front(); }
         reference back() { return data_.back(); }
         const_reference back() const { return data_.back(); }
         const T* data() const { return data_.data(); }
         T* c_array() { return data_.data(); }

         // modifiers
         void swap(Vector<T, N>& other) { data_.swap(other.data_); }
         void fill(const T& t) {
             data_.fill(t);
         }

        /// In-place, element-wise multiplcation by a scalar.

        /// \tparam Q Type of the scalar.
        /// \param[in] q The scalar.
        /// \return A reference to this for chaining operations.
        /// \todo Do we want a similar division operation?
        template <typename Q>
        Vector<T,N>& operator*=(Q q) {
            for(size_type i = 0; i < N; ++i)
                data_[i] *= q;
            return *this;
        }

        /// In-place, element-wise addition of another \c Vector.

        /// \tparam Q Type stored in the other \c Vector.
        /// \param[in] q The other \c Vector.
        /// \return A reference to this for chaining operations.
        template <typename Q>
        Vector<T,N>& operator+=(const Vector<Q,N>& q) {
            for(size_type i = 0; i < N; ++i)
                data_[i] += q[i];
            return *this;
        }

        /// In-place, element-wise subtraction of another \c Vector.

        /// \tparam Q Type stored in the other \c Vector.
        /// \param[in] q The other \c Vector.
        /// \returns A reference to this for chaining operations.
        template <typename Q>
        Vector<T,N>& operator-=(const Vector<Q,N>& q) {
            for(size_type i = 0; i < N; ++i)
                data_[i] -= q[i];
            return *this;
        }

        /// Calculate the 2-norm of the vector elements.

        /// \return The 2-norm.
        T normf() const {
        	T d=0;
        	for (std::size_t i=0; i<N; ++i) d+=(data_[i])*(data_[i]);
        	return sqrt(d);
        }

        /// Support for MADNESS serialization.

        /// \tparam Archive The archive type.
        /// \param[in,out] ar The archive.
        template <typename Archive>
        void serialize(Archive& ar) {
            ar & data_;
        }

        /// Support for MADNESS hashing.

        /// \return The hash.
        hashT hash() const {
            return hash_value(data_);
        }

        // comparisons
        friend bool operator==(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ == r.data_;
        }

        friend bool operator!=(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ != r.data_;
        }

        friend bool operator<(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ < r.data_;
        }

        friend bool operator>(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ > r.data_;
        }

        friend bool operator<=(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ <= r.data_;
        }

        friend bool operator>=(const Vector<T, N>& l, const Vector<T, N>& r) {
            return l.data_ >= r.data_;
        }

        /// Output a \c Vector to stream.

        /// \param[in] s The output stream.
        /// \param[in] v The \c Vector to output.
        /// \return The output stream.
        friend std::ostream& operator<<(std::ostream& s, const Vector<T,N>& v) {
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
    void swap(Vector<T,N>& l, Vector<T,N>& r) {
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
    Vector<T,N> operator*(Vector<T,N> l, U r) {
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
    Vector<T,N> operator*(T l, Vector<U,N> r) {
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
    Vector<T,N> operator*(Vector<T,N> l, const Vector<U,N>& r) {
        // coordinate r passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] *= r[i];
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
    Vector<T,N> operator+(Vector<T,N> l, U r) {
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
    Vector<T,N> operator+(Vector<T,N> l, const Vector<U,N>& r) {
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
    Vector<T,N> operator-(Vector<T,N> l, U r) {
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
    Vector<T,N> operator-(Vector<T,N> l, const Vector<U,N>& r) {
        // l passed by value to allow compiler optimization
        for (std::size_t i = 0; i < N; ++i)
            l[i] -= r[i];
        return l;
    }

    /// Compute norm of a \c Vector

    /// \tparam T The \c Vector element type
    /// \tparam N The \c Vector size
    /// \param v The \c Vector
    /// \return The vector norm, \f$ ||v||_2 = \sqrt{\sum_{k=1}^N v_i^2} \f$
    /// \todo Duplicate function with the normf member function. Delete soon.
    template <typename T, std::size_t N>
    T norm(Vector<T,N> v) {
      T norm2 = 0.0;
      for (std::size_t i = 0; i < N; ++i)
        norm2 += v[i] * v[i];
      return sqrt(norm2);
    }

    /// Your friendly neighborhood factory function

    /// \todo Replace this set of factory functions by a variadic template.
    template <typename T>
    Vector<T,1> vec(T x) {
        Vector<T,1> r; r[0] = x;
        return r;
    }

    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,2> vec(T x, T y) {
        Vector<T,2> r; r[0] = x; r[1] = y;
        return r;
    }

    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,3> vec(T x, T y, T z) {
        Vector<T,3> r; r[0] = x; r[1] = y; r[2] = z;
        return r;
    }

    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,4> vec(T x, T y, T z, T xx) {
        Vector<T,4> r; r[0] = x; r[1] = y; r[2] = z; r[3] = xx;
        return r;
    }

    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,5> vec(T x, T y, T z, T xx, T yy) {
        Vector<T,5> r; r[0] = x; r[1] = y; r[2] = z; r[3] = xx; r[4] = yy;
        return r;
    }

    /// Your friendly neighborhood factory function
    template <typename T>
    Vector<T,6> vec(T x, T y, T z, T xx, T yy, T zz) {
        Vector<T,6> r; r[0] = x; r[1] = y; r[2] = z; r[3] = xx; r[4] = yy; r[5] = zz;
        return r;
    }


	/// Construct a unit-`Vector` that has the same direction as \c r.

    /// \tparam T The type of data stored in the \c Vector.
    /// \tparam NDIM The size of the \c Vector.
    /// \param[in] r The \c Vector.
    /// \param[in] eps A precision. If the norm of \c r is less than \c eps,
    ///    the zero vector is returned.
    /// \return The desired unit-`Vector` (unless \c r is numerically the zero
    ///    \c Vector).
    /// \todo Make NDIM into N for consistency. Also, give this function a
    ///    much more useful name.
    template<typename T, std::size_t NDIM>
	Vector<T,NDIM> n12(const Vector<T,NDIM>& r, const double eps=1.e-6) {
		const double norm=r.normf();
		if (norm<1.e-6) return Vector<T,NDIM>(0.0);
		return r*(1.0/norm);
	}

} // namespace madness

#endif // MADNESS_WORLD_VECTOR_H__INCLUDED
