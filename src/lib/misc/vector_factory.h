#ifndef VECTOR_FACTORY_H
#define VECTOR_FACTORY_H

#include <vector>

/// \file vector_factory.h
/// \brief Declares and implements factories for short vectors

namespace madness {
    
    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
        inline std::vector<T> vector_factory(const T& v0) {
        std::vector<T> v(1);
        v[0] = v0;
        return v;
    }
    
    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
        inline std::vector<T> vector_factory(const T& v0, const T& v1) {
        std::vector<T> v(2);
        v[0] = v0;
        v[1] = v1;
        return v;
    }
    
    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
        inline std::vector<T> vector_factory(const T& v0, const T& v1, 
                                             const T& v2) {
        std::vector<T> v(3);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        return v;
    }
    
    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
        inline std::vector<T> vector_factory(const T& v0, const T& v1, 
                                             const T& v2, const T& v3) {
        std::vector<T> v(4);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        return v;
    }
    
    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
        inline std::vector<T> vector_factory(const T& v0, const T& v1, 
                                             const T& v2, const T& v3, 
                                             const T& v4) {
        std::vector<T> v(5);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
        return v;
    }
    
    /// Returns a std::vector<T> initialized from the arguments
    template <typename T>
        inline std::vector<T> vector_factory(const T& v0, const T& v1, 
                                             const T& v2, const T& v3, 
                                             const T& v4, const T& v5) {
        std::vector<T> v(6);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
        v[5] = v5;
        return v;
    }
}

#endif
