#ifndef WORLD_HASH_H
#define WORLD_HASH_H

/// \file worldhash.h
/// \brief Defines hash functions for use in distributed containers


#define WORLDDC_USES_GNU_HASH_MAP
#define HASH_MAP_NAMESPACE __gnu_cxx

#ifdef WORLDDC_USES_GNU_HASH_MAP
#  ifdef __GNUG__
#    include <ext/hash_map>
#  else
#    include <hash_map>
#  endif
#endif


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
    static
    inline
    typename madness::enable_if<madness::is_fundamental<T>, hashT>::type
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
    static
    inline
    typename madness::disable_if<madness::is_fundamental<T>, hashT>::type
    hash(const T& t, hashT initval=0) {
        hashT h = Hash<T>::hash(t);
        if (initval) h = hashword(&h, 1, initval);
        return h;
    }

    /// Hash a variable sized array
    template <class T>
    static
    inline
    typename madness::enable_if<madness::is_fundamental<T>, hashT>::type
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
    static
    inline
    typename madness::disable_if<madness::is_fundamental<T>, hashT>::type
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
        static hashT hash(const T (&t)[n], hashT initval=0) {
            return madness::hash(t, n, initval);
        };
    };
}


#endif
