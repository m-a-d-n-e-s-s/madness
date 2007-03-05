#ifndef MAD_VECAR_H
#define MAD_VECAR_H

/// \file vecar.h
/// \brief Implements archive wrapping an STL vector

#include <vector>
#include <world/archive.h>


namespace madness {
namespace archive {

// With a bit of thought this could be generalized to several STL containers


        /// Wraps an archive around an STL vector for output
class VectorOutputArchive : public BaseOutputArchive {
    mutable std::vector<unsigned char>& v;
public:
    VectorOutputArchive(std::vector<unsigned char>& v, std::size_t hint=262144) : v(v) {
        open(hint);
    };

    template <class T>
    inline
    typename madness::enable_if< madness::is_fundamental<T>, void >::type
    store(const T* t, long n) const {
        const unsigned char* ptr = (unsigned char*) t;
        v.insert(v.end(),ptr,ptr+n*sizeof(T));
    }

    void open(std::size_t hint=262144) {
        v.clear();
        v.reserve(hint);
    };

    void close() {};

    void flush() {};
};


        /// Wraps an archive around an STL vector for input
class VectorInputArchive : public BaseInputArchive {
    mutable std::vector<unsigned char>& v;
    mutable std::size_t i;
public:
    VectorInputArchive(std::vector<unsigned char>& v) : v(v) , i(0) {}

    template <class T>
    inline
    typename madness::enable_if< madness::is_fundamental<T>, void >::type
    load(T* t, long n) const {
        std::size_t m = n*sizeof(T);
        if (m+i >  v.size()) throw "VectorInputArchive: reading past end";
        memcpy((unsigned char*) t, &v[i], m);
        i += m;
    }

    void open() {};

    void rewind() const {i=0;};

    std::size_t nbyte_avail() const {return v.size()-i;};

    void close() {}
};
}
}
#endif
