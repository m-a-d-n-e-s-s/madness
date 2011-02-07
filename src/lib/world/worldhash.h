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


  $Id$
*/


#ifndef MADNESS_WORLD_WORLDHASH_H__INCLUDED
#define MADNESS_WORLD_WORLDHASH_H__INCLUDED

/// \file worldhash.h
/// \brief Defines hash functions for use in distributed containers

#include <madness_config.h>
#include <world/typestuff.h>
#include <cstddef>

namespace madness {

    typedef unsigned int hashT;

    extern "C" hashT hashlittle(const void *key, size_t length, int initval);
    extern "C" hashT hashword(const hashT *k, size_t length, hashT initval);

    static inline hashT hashulong(const unsigned long* k, size_t length, unsigned long initval) {
        return hashword((const hashT*) k, length*sizeof(long)/sizeof(hashT),
                        (hashT)(initval));
    }

    template <typename T> struct Hash;

    /// Hash a single instance
    template <class T>
#ifndef HAVE_IBMBGP
    static
#endif
    inline
    typename madness::enable_if<std::is_fundamental<T>, hashT>::type
    hash(const T& t, hashT initval=0) {
        // Use heavily optimized hashword when sizeof(T) is multiple
        // of sizeof(hashT) and presumably correctly aligned.
        if (((sizeof(T)/sizeof(hashT))*sizeof(hashT)) == sizeof(T)) {
            return hashword((const hashT *) &t, sizeof(T)/sizeof(hashT), initval);
        }
        else {
            return hashlittle((const void *) &t, sizeof(T), initval);
        }
    }

    template <class T>
#ifndef HAVE_IBMBGP
    static
#endif
    inline
    typename madness::disable_if<std::is_fundamental<T>, hashT>::type
    hash(const T& t, hashT initval=0) {
        hashT h = Hash<T>::hash(t);
        if (initval) h = hashword(&h, 1, initval);
        return h;
    }

    /// Hash a variable sized array
    template <class T>
#ifndef HAVE_IBMBGP
    static
#endif
    inline
    typename madness::enable_if<std::is_fundamental<T>, hashT>::type
    hash(const T* t, std::size_t n, hashT initval=0) {
        // Use heavily optimized hashword when sizeof(T) is multiple
        // of sizeof(hashT)
        if (((sizeof(T)/sizeof(hashT))*sizeof(hashT)) == sizeof(T)) {
            //std::cout << "hashing words ";
            //for (int i=0; i<n; i++) std::cout << t[i] << " ";
            hashT result = hashword((const hashT *) t, n*sizeof(T)/sizeof(hashT), initval);
            //std::cout << " ---> " << result << std::endl;
            return result;
        }
        else {
            return hashlittle((void *) t, n*sizeof(T), initval);
        }
    }

    template <class T>
#ifndef HAVE_IBMBGP
    static
#endif
    inline
    typename madness::disable_if<std::is_fundamental<T>, hashT>::type
    hash(const T* t, std::size_t n, hashT initval=0) {
        hashT sum=0;
        for (std::size_t i=0; i<n; i++) sum = hash(t[i],sum);
        return sum;
    }

    /// Default \c Hash<T>::hash(t) invokes t.hash()
    template <typename T>
    struct Hash {
        static hashT hash(const T& t) {
            return t.hash();
        };
    };

    /// Specialization for fixed dim arrays invokes hash(t,n)
    template <class T, std::size_t n>
    struct Hash<T[n]> {
        static hashT hash(const T(&t)[n], hashT initval=0) {
            return madness::hash(t, n, initval);
        }
    };

    /// Default \c Hash<T>::hash(t) invokes t.hash()
    template <typename T>
    struct Hash<T*> {
        static hashT hash(void* p) {
            unsigned long n = reinterpret_cast<unsigned long>(p);
            return madness::hash(n);
        };
    };

}


#endif // MADNESS_WORLD_WORLDHASH_H__INCLUDED
