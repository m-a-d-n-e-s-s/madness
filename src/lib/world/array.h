/*
  This file is part of MADNESS.

  Copyright (C) <2007> <Oak Ridge National Laboratory>

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


  $Id$
*/

#ifndef MADNESS_WORLD_ARRAY_H__INCLUDED
#define MADNESS_WORLD_ARRAY_H__INCLUDED

#include <vector>
#include <world/worldexc.h>
#include <world/worldhash.h>

namespace madness {

    /// A simple, fixed dimension Vector

    /// Eliminates memory allocation cost, is just POD so can be
    /// copied easily and allocated on the stack, and the known
    /// dimension permits agressive compiler optimizations.
    template <typename T, int N>
    class Vector {
    private:
        T v[N];
    public:
        typedef T* iterator;
        typedef const T* const_iterator;

        /// Default constructor does not initialize vector contents
        Vector() {}

        /// Initialize all elements to value t
        Vector(T t) {
            for (int i=0; i<N; i++) v[i] = t;
        }

        /// Construct from a C++ array of the same dimension
        Vector(const T(&t)[N])  {
            for (int i=0; i<N; i++) v[i] = t[i];
        }

        /// Construct from an STL vector of equal or greater length
        template <typename Q>
        Vector(const std::vector<Q> t) {
            MADNESS_ASSERT(t.size() >= N);
            for (int i=0; i<N; i++) v[i] = t[i];
        }

        /// Copy constructor is deep (since a vector is POD)
        Vector(const Vector<T,N>& other) {
            for (int i=0; i<N; i++) v[i] = other.v[i];
        }

        /// Assignment is deep (since a vector is POD)
        Vector& operator=(const Vector<T,N>& other) {
            for (int i=0; i<N; i++) v[i] = other.v[i];
            return *this;
        }

        /// Assignment is deep (since a vector is POD)
        Vector& operator=(const std::vector<T>& other) {
            for (int i=0; i<N; i++) v[i] = other[i];
            return *this;
        }

        /// Fill from scalar value
        Vector& operator=(T t) {
            for (int i=0; i<N; i++) v[i] = t;
            return *this;
        }

        /// Test for element-wise equality
        bool operator==(const Vector<T,N>& other) const {
            for (int i=0; i<N; i++)
                if (v[i] != other.v[i]) return false;
            return true;
        }

        /// Same as !(*this==other)
        bool operator!=(const Vector<T,N>& other) const {
            return !(*this==other);
        }

        /// Tests a<b in sense of lexically ordered index

        /// Can be used to impose a standard ordering on vectors.
        bool operator<(const Vector<T,N>& other) const {
            for (int i=0; i<N; i++) {
                if (v[i] < other.v[i]) return true;
                else if (v[i] > other.v[i]) return false;
            }
            if (v[N-1] == other.v[N-1]) return false; // equality
            return true;
        }

        /// Indexing
        T& operator[](int i) {
            return v[i];
        }

        /// Indexing
        const T& operator[](int i) const {
            return v[i];
        }

        /// Element-wise multiplcation by a scalar

        /// Returns a new vector
        template <typename Q>
        Vector<T,N> operator*(Q q) const {
            Vector<T,N> r;
            for (int i=0; i<N; i++) r[i] = v[i] * q;
            return r;
        }

        /// In-place element-wise multiplcation by a scalar

        /// Returns a reference to this for chaining operations
        template <typename Q>
        Vector<T,N>& operator*=(Q q) {
            for (int i=0; i<N; i++) v[i] *= q;
            return *this;
        }

        /// Element-wise multiplcation by another vector

        /// Returns a new vector
        template <typename Q>
        Vector<T,N> operator*(const Vector<Q,N>& q) const {
            Vector<T,N> r;
            for (int i=0; i<N; i++) r[i] = v[i]*q[i];
            return r;
        }

        /// Element-wise addition of a scalar

        /// Returns a new vector
        template <typename Q>
        Vector<T,N> operator+(Q q) const {
            Vector<T,N> r;
            for (int i=0; i<N; i++) r[i] = v[i] + q;
            return r;
        }

        /// In-place element-wise addition of a scalar

        /// Returns a reference to this for chaining operations
        template <typename Q>
        Vector<T,N>& operator+=(Q q) {
            for (int i=0; i<N; i++) v[i] += q;
            return *this;
        }

        /// Element-wise addition of another vector

        /// Returns a new vector
        template <typename Q>
        Vector<T,N> operator+(const Vector<Q,N>& q) const {
            Vector<T,N> r;
            for (int i=0; i<N; i++) r[i] = v[i] + q[i];
            return r;
        }

        /// Element-wise subtraction of a scalar

        /// Returns a new vector
        template <typename Q>
        Vector<T,N> operator-(Q q) const {
            Vector<T,N> r;
            for (int i=0; i<N; i++) r[i] = v[i] - q;
            return r;
        }

        /// Element-wise subtraction of another vector

        /// Returns a new vector
        template <typename Q>
        Vector<T,N> operator-(const Vector<Q,N>& q) const {
            Vector<T,N> r;
            for (int i=0; i<N; i++) r[i] = v[i] - q[i];
            return r;
        }

        /// STL iterator support
        iterator begin() {
            return v;
        }

        /// STL iterator support
        const_iterator begin() const {
            return v;
        }

        /// STL iterator support
        iterator end() {
            return v+N;
        }

        /// STL iterator support
        const_iterator end() const {
            return v+N;
        }

        /// Length of the vector
        int size() const {
            return N;
        }

        /// Support for MADNESS serialization
        template <typename Archive>
        void serialize(Archive& ar) {
            ar & v;
        }

        /// Support for MADNESS hashing
        hashT hash() const {
            return madness::hash(v);
        }
    };

    /// Output vector to stream for human consumption
    template <typename T,int N>
    std::ostream& operator<<(std::ostream& s, const Vector<T,N>& a) {
        s << "[";
        for (int i=0; i<N; i++) {
            s << a[i];
            if (i != (N-1)) s << ",";
        }
        s << "]";
        return s;
    }


    /// A simple, fixed-size, stack
    template <typename T,int N>
    class Stack {
    private:
        T t[N];
//         Vector<T,N> t;
        int n;
        
    public:
        Stack() : n(0) {}

        void push(const T& value) {
            MADNESS_ASSERT(n < N);
            t[n++] = value;
        }

        T& pop() {
            MADNESS_ASSERT(n > 0);
            return t[--n];
        }

        T& front() {
            MADNESS_ASSERT(n > 0);
            return t[n-1];
        }

        int size() const {
            return n;
        }

        void clear() {
            n = 0;
        }
    };


    // Sigh ...vector_factory already in use by tensor

    /// Returns a Vector<T,1> initialized from the arguments
    template <typename T>
    inline Vector<T,1> VectorFactory(const T& v0) {
        Vector<T,1> v;
        v[0] = v0;
        return v;
    }

    /// Returns a Vector<T,2> initialized from the arguments
    template <typename T>
    inline Vector<T,2> VectorFactory(const T& v0, const T& v1) {
        Vector<T,2> v;
        v[0] = v0;
        v[1] = v1;
        return v;
    }

    /// Returns a Vector<T,3> initialized from the arguments
    template <typename T>
    inline Vector<T,3> VectorFactory(const T& v0, const T& v1,
                                     const T& v2) {
        Vector<T,3> v;
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        return v;
    }

    /// Returns a Vector<T,4> initialized from the arguments
    template <typename T>
    inline Vector<T,4> VectorFactory(const T& v0, const T& v1,
                                     const T& v2, const T& v3) {
        Vector<T,4> v;
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        return v;
    }

    /// Returns a Vector<T,5> initialized from the arguments
    template <typename T>
    inline Vector<T,5> VectorFactory(const T& v0, const T& v1,
                                     const T& v2, const T& v3,
                                     const T& v4) {
        Vector<T,5> v;
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
        return v;
    }

    /// Returns a Vector<T,6> initialized from the arguments
    template <typename T>
    inline Vector<T,6> VectorFactory(const T& v0, const T& v1,
                                     const T& v2, const T& v3,
                                     const T& v4, const T& v5) {
        Vector<T,6> v;
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
        v[5] = v5;
        return v;
    }
}
#endif // MADNESS_WORLD_ARRAY_H__INCLUDED
