#ifndef MAD_VECAR_H
#define MAD_VECAR_H

#include <vector>
#include <serialize/archive.h>


namespace madness {
    namespace archive {

        // With a bit of thought this could be generalized to several STL containers
        class VectorOutputArchive : public BaseOutputArchive {
            mutable std::vector<unsigned char>& v;
        public:
            VectorOutputArchive(std::vector<unsigned char>& v) : v(v) {
                open();
            };

            template <class T>
            inline
            typename madness::enable_if< madness::is_fundamental<T>, void >::type
            store(const T* t, long n) const {
                const unsigned char* ptr = (unsigned char*) t;
                v.insert(v.end(),ptr,ptr+n*sizeof(T));
            }

            void open() {
                v.clear();
                v.reserve(262144);
            };

            void close() {};

            void flush() {};
        };


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
