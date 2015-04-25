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

#ifndef MADNESS_WORLD_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_ARCHIVE_H__INCLUDED

/**
 \file world/archive.h
 \brief Interface templates for the archives (serialization).
 \ingroup serialization
*/

#include <complex>
#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
//#include <madness/world/worldprofile.h>
#include <madness/world/enable_if.h>
#include <madness/world/type_traits.h>
#include <madness/world/madness_exception.h>

#define ARCHIVE_COOKIE "archive"
#define ARCHIVE_MAJOR_VERSION 0
#define ARCHIVE_MINOR_VERSION 1

//#define MAD_ARCHIVE_DEBUG_ENABLE

#ifdef MAD_ARCHIVE_DEBUG_ENABLE
#define MAD_ARCHIVE_DEBUG(s) s
//using std::endl;
#else
#define MAD_ARCHIVE_DEBUG(s)
#endif

namespace madness {


    // Forward declarations
    template <typename T> class Tensor;

    namespace archive {

        // Forward declarations
        template <class>
        class archive_array;
        template <class T>
        inline archive_array<T> wrap(const T*, unsigned int);
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T*, unsigned int);
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T&);

        // There are 64 empty slots for user types.  Free space for
        // registering user types begins at cookie=128.

#ifdef MAD_ARCHIVE_TYPE_NAMES_CC
        const char *archive_type_names[256];
#else
        extern const char *archive_type_names[256];
#endif
        void archive_initialize_type_names();

        // Used to enable type checking inside archives
        template <typename T>
        struct archive_typeinfo {
            static const unsigned char cookie = 255; ///< 255 indicates unknown type
        };

        // Returns the name of the type, or unknown if not registered.
        template <typename T>
        const char* get_type_name() {
            return archive_type_names[archive_typeinfo<T>::cookie];
        }


#if defined(ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE) && defined(ARCHIVE_REGISTER_TYPE_IBMBUG)
#define ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T) \
        ; const unsigned char archive_typeinfo< T >::cookie
#else
#define ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T)
#endif

        /// \def ARCHIVE_REGISTER_TYPE(T, cooky)
        /// \brief Used to associate type with cookie value inside archive
        ///
        /// \def  ARCHIVE_REGISTER_TYPE_AND_PTR(T, cooky)
        /// \brief Used to associate type and ptr to type with cookie value inside archive
        ///
        /// \def ARCHIVE_REGISTER_TYPE_NAME(T)
        /// \brief Used to associate names with types
        ///
        /// \def ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(T)
        /// \brief Used to associate names with types and pointers

#define ARCHIVE_REGISTER_TYPE(T, cooky) \
	template <> struct archive_typeinfo< T > { \
		static const unsigned char cookie = cooky; \
	} \
        ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T)

#define ARCHIVE_REGISTER_TYPE_AND_PTR(T, cooky) \
        ARCHIVE_REGISTER_TYPE(T, cooky); \
        ARCHIVE_REGISTER_TYPE(T*, cooky+64)

#define ATN ::madness::archive::archive_type_names
#define ATI ::madness::archive::archive_typeinfo
#define ARCHIVE_REGISTER_TYPE_NAME(T) \
     if (strcmp(ATN[ATI< T >::cookie],"invalid")) {\
        std::cout << "archive_register_type_name: slot/cookie already in use! "<< #T << " " << ATN[ATI< T >::cookie]<< std::endl; \
        MADNESS_EXCEPTION("archive_register_type_name: slot/cookie already in use!", 0); \
     } \
     ATN[ATI< T >::cookie] = #T

#define ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(T) \
     ARCHIVE_REGISTER_TYPE_NAME(T); \
     ARCHIVE_REGISTER_TYPE_NAME(T*)

        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned char,0);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned short,1);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned int,2);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned long,3);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned long long,4);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed char,5);
        ARCHIVE_REGISTER_TYPE_AND_PTR(char,5);	// Needed, but why?
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed short,6);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed int,7);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed long,8);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed long long,9);
        ARCHIVE_REGISTER_TYPE_AND_PTR(bool,10);
        ARCHIVE_REGISTER_TYPE_AND_PTR(float,11);
        ARCHIVE_REGISTER_TYPE_AND_PTR(double,12);
        ARCHIVE_REGISTER_TYPE_AND_PTR(long double,13);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::complex<float>,14);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::complex<double>,15);

        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<char>,20);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned char>,21);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<short>,22);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned short>,23);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<int>,24);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned int>,25);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<long>,26);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned long>,27);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<bool>,28);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<float>,29);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<double>,30);

        ARCHIVE_REGISTER_TYPE_AND_PTR(std::string,31);

        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<int>,32);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<long>,33);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<float>,34);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<double>,35);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor< std::complex<float> >,36);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor< std::complex<double> >,37);

        /// Base class for all archives classes
        class BaseArchive {
        public:
            static const bool is_archive = true;
            static const bool is_input_archive = false;
            static const bool is_output_archive = false;
            static const bool is_parallel_archive = false;
            BaseArchive() {
                archive_initialize_type_names();
            }
        }; // class BaseArchive

        /// Base class for input archives classes
        class BaseInputArchive : public BaseArchive {
        public:
            static const bool is_input_archive = true;
        }; // class BaseInputArchive


        /// Base class for output archives classes
        class BaseOutputArchive : public BaseArchive {
        public:
            static const bool is_output_archive = true;
        }; // class BaseOutputArchive

        /// Checks that \c T is an archive type

        /// If \c T is an archive type then \c is_archive will be inherited from
        /// \c std::true_type , otherwise it is inherited from
        /// \c std::false_type .
        template <typename T>
        struct is_archive : public std::is_base_of<BaseArchive, T> {};

        /// Checks that \c T is an input archive type

        /// If \c T is an input archive type then \c is_archive will be
        /// inherited from \c std::true_type , otherwise it is inherited from
        /// \c std::false_type .
        template <typename T>
        struct is_input_archive : public std::is_base_of<BaseInputArchive, T> {};

        /// Checks that \c T is an output archive type

        /// If \c T is an output archive type then \c is_archive will be
        /// inherited from \c std::true_type , otherwise it is inherited from
        /// \c std::false_type .
        template <typename T>
        struct is_output_archive : public std::is_base_of<BaseOutputArchive, T> {};

        // Serialize an array of fundamental stuff
        template <class Archive, class T>
        typename enable_if_c< is_serializable<T>::value && is_output_archive<Archive>::value >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize fund array" << std::endl);
            ar.store(t,n);
        }


        // Deserialize an array of fundamental stuff
        template <class Archive, class T>
        typename enable_if_c< is_serializable<T>::value && is_input_archive<Archive>::value >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "deserialize fund array" << std::endl);
            ar.load((T*) t,n);
        }


        // (de)Serialize an array of non-fundamental stuff
        template <class Archive, class T>
        typename enable_if_c< ! is_serializable<T>::value && is_archive<Archive>::value >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize non-fund array" << std::endl);
            for (unsigned int i=0; i<n; ++i) ar & t[i];
        }


        /// Default implementation of pre/postamble
        template <class Archive, class T>
        struct ArchivePrePostImpl {
            static inline void preamble_load(const Archive& ar) {
                unsigned char ck = archive_typeinfo<T>::cookie;
                unsigned char cookie;
                ar.load(&cookie, 1); // cannot use >>
                if (cookie != ck) {
                    char msg[255];
                    std::sprintf(msg,"InputArchive type mismatch: expected cookie "
                                 "%u (%s) but got %u (%s) instead",
                                 ck, archive_type_names[ck],
                                 cookie,archive_type_names[cookie]);
                    std::cerr << msg << std::endl;
                    MADNESS_EXCEPTION(msg, static_cast<int>(cookie));
                }
                else {
                    MAD_ARCHIVE_DEBUG(std::cout << "read cookie " << archive_type_names[cookie] << std::endl);
                }
            }


            /// Serialize a cookie for type checking
            static inline void preamble_store(const Archive& ar) {
                unsigned char ck = archive_typeinfo<T>::cookie;
                ar.store(&ck, 1); // cannot use <<
                MAD_ARCHIVE_DEBUG(std::cout << "wrote cookie " << archive_type_names[ck] << std::endl);
            }


            /// By default there is no postamble
            static inline void postamble_load(const Archive& /*ar*/) {}

            /// By default there is no postamble
            static inline void postamble_store(const Archive& /*ar*/) {}
        };


        /// Default symmetric serialization of a non-fundamental thingy
        template <class Archive, class T>
        struct ArchiveSerializeImpl {
            static inline void serialize(const Archive& ar, T& t) {
                t.serialize(ar);
            }
        };


        // Redirect \c serialize(ar,t) to \c serialize(ar,&t,1) for fundamental types
        template <class Archive, class T>
        inline
        typename enable_if_c< is_serializable<T>::value && is_archive<Archive>::value >::type
        serialize(const Archive& ar, const T& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize(ar,t) -> serialize(ar,&t,1)" << std::endl);
            serialize(ar,&t,1);
        }


        // Redirect \c serialize(ar,t) to \c ArchiveSerializeImpl for non-fundamental types
        template <class Archive, class T>
        inline
        typename enable_if_c< !is_serializable<T>::value && is_archive<Archive>::value >::type
        serialize(const Archive& ar, const T& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize(ar,t) -> ArchiveSerializeImpl" << std::endl);
            ArchiveSerializeImpl<Archive,T>::serialize(ar,(T&) t);
        }


        /// Default store of a thingy via serialize(ar,t)
        template <class Archive, class T>
        struct ArchiveStoreImpl {
            static inline void store(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "store(ar,t) default" << std::endl);
                serialize(ar,t);
            }
        };


        /// Default load of a thingy via serialize(ar,t)
        template <class Archive, class T>
        struct ArchiveLoadImpl {
            static inline void load(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "load(ar,t) default" << std::endl);
                serialize(ar,t);
            }
        };


        /// Default implementation of wrap_store and wrap_load
        template <class Archive, class T>
        struct ArchiveImpl {
            static inline const Archive& wrap_store(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for default" << std::endl);
                ArchivePrePostImpl<Archive,T>::preamble_store(ar);
                ArchiveStoreImpl<Archive,T>::store(ar,t);
                ArchivePrePostImpl<Archive,T>::postamble_store(ar);
                return ar;
            }

            static inline const Archive& wrap_load(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for default" << std::endl);
                ArchivePrePostImpl<Archive,T>::preamble_load(ar);
                ArchiveLoadImpl<Archive,T>::load(ar,(T&) t);  // Loses constness here!
                ArchivePrePostImpl<Archive,T>::postamble_load(ar);
                return ar;
            }
        };


        // Redirect \c << to ArchiveImpl::wrap_store for output archives
        template <class Archive, class T>
        inline
        typename enable_if<is_output_archive<Archive>, const Archive&>::type
        operator<<(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_store(ar,t);
        }

        // Redirect \c >> to ArchiveImpl::wrap_load for input archives
        template <class Archive, class T>
        inline
        typename enable_if<is_input_archive<Archive>, const Archive&>::type
        operator>>(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_load(ar,t);
        }

        // Redirect \c & to ArchiveImpl::wrap_store for output archives
        template <class Archive, class T>
        inline
        typename enable_if<is_output_archive<Archive>, const Archive&>::type
        operator&(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_store(ar,t);
        }

        // Redirect \c & to ArchiveImpl::wrap_load for input archives
        template <class Archive, class T>
        inline
        typename enable_if<is_input_archive<Archive>, const Archive&>::type
        operator&(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_load(ar,t);
        }


        ///////////////////////////////////////////////////////////////

        /// Wrapper for opaque pointer ... bitwise copy of the pointer ... no remapping performed
        template <class T>
        class archive_ptr {
        public:
            T* ptr;

            archive_ptr(T* t = 0) : ptr(t) {}

            T& operator*(){return *ptr;}

            template <class Archive>
            void serialize(const Archive& ar) {ar & wrap_opaque(&ptr, 1);}
        };

	template <class T>
	inline
	archive_ptr<T> wrap_ptr(T* p) {return archive_ptr<T>(p);}

        /// Wrapper for dynamic arrays and pointers
        template <class T>
        class archive_array {
        public:
            const T* ptr;
            unsigned int n;

            archive_array(const T *ptr, unsigned int n) : ptr(ptr), n(n) {}

            archive_array() : ptr(0), n(0) {}
        };


        /// Factory function to wrap dynamically allocated pointer as typed archive_array
        template <class T>
        inline
        archive_array<T> wrap(const T* ptr, unsigned int n) {
            return archive_array<T>(ptr,n);
        }

        /// Factory function to wrap pointer to contiguous data as opaque (uchar) archive_array
        template <class T>
        inline
        archive_array<unsigned char> wrap_opaque(const T* ptr, unsigned int n) {
            return archive_array<unsigned char>((unsigned char*) ptr, n*sizeof(T));
        }

        /// Factory function to wrap contiguous scalar as opaque (uchar) archive_array
        template <class T>
        inline
        archive_array<unsigned char> wrap_opaque(const T& t) {
            return archive_array<unsigned char>((unsigned char*) &t,sizeof(t));
        }

        /// Serialize function pointer
        template <class Archive, typename resT, typename... paramT>
        struct ArchiveSerializeImpl<Archive, resT(*)(paramT...)> {
            static inline void serialize(const Archive& ar, resT(*fn)(paramT...)) {
                ar & wrap_opaque(fn);
            }
        };

        /// Serialize member function pointer
        template <class Archive, typename resT, typename objT, typename... paramT>
        struct ArchiveSerializeImpl<Archive, resT(objT::*)(paramT...)> {
            static inline void serialize(const Archive& ar, resT(objT::*memfn)(paramT...)) {
                ar & wrap_opaque(memfn);
            }
        };

        /// Serialize const member function pointer
        template <class Archive, typename resT, typename objT, typename... paramT>
        struct ArchiveSerializeImpl<Archive, resT(objT::*)(paramT...) const> {
            static inline void serialize(const Archive& ar, resT(objT::*memfn)(paramT...) const) {
                ar & wrap_opaque(memfn);
            }
        };

        /// Partial specialization for archive_array

        /// This makes use of stuff that user specializations need not
        template <class Archive, class T>
        struct ArchiveImpl< Archive, archive_array<T> > {
            static inline const Archive& wrap_store(const Archive& ar, const archive_array<T>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for archive_array" << std::endl);
                ArchivePrePostImpl<Archive,T*>::preamble_store(ar);
                //ar << t.n;
                //ArchivePrePostImpl<Archive,T>::preamble_store(ar);
                serialize(ar,(T *) t.ptr,t.n);
                //ArchivePrePostImpl<Archive,T>::postamble_store(ar);
                ArchivePrePostImpl<Archive,T*>::postamble_store(ar);
                return ar;
            }

            static inline const Archive& wrap_load(const Archive& ar, const archive_array<T>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for archive_array" << std::endl);
                ArchivePrePostImpl<Archive,T*>::preamble_load(ar);
                //unsigned int n;
                //ar >> n;
                //if (n != t.n)
                //    MADNESS_EXCEPTION("deserializing archive_array: dimension mismatch", n);
                //ArchivePrePostImpl<Archive,T>::preamble_load(ar);
                serialize(ar,(T *) t.ptr,t.n);
                //ArchivePrePostImpl<Archive,T>::postamble_load(ar);
                ArchivePrePostImpl<Archive,T*>::postamble_load(ar);
                return ar;
            }
        };


        /// Partial specialization for fixed dimension array redirects to archive_array
        template <class Archive, class T, std::size_t n>
        struct ArchiveImpl<Archive, T[n]> {
            static inline const Archive& wrap_store(const Archive& ar, const T(&t)[n]) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for array" << std::endl);
                ar << wrap(&t[0],n);
                return ar;
            }

            static inline const Archive& wrap_load(const Archive& ar, const T(&t)[n]) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for array" << std::endl);
                ar >> wrap(&t[0],n);
                return ar;
            }
        };


        /// Serialize a complex number
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, std::complex<T> > {
            static inline void store(const Archive& ar, const std::complex<T>& c) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize complex number" << std::endl);
                ar & c.real() & c.imag();
            }
        };


        /// Deserialize a complex number
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, std::complex<T> > {
            static inline void load(const Archive& ar, std::complex<T>& c) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize complex number" << std::endl);
                T r = 0, i = 0;
                ar & r & i;
                c = std::complex<T>(r,i);
            }
        };


        /// Serialize STL vector.
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, std::vector<T> > {
            static inline void store(const Archive& ar, const std::vector<T>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL vector" << std::endl);
                ar & v.size();
                ar & wrap(&v[0],v.size());
            }
        };


        /// Deserialize STL vector. Clears & resizes as necessary.
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, std::vector<T> > {
            static void load(const Archive& ar, std::vector<T>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL vector" << std::endl);
                std::size_t n = 0ul;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                ar & wrap((T *) &v[0],n);
            }
        };

        /// Serialize STL vector<bool> (as plain array of bool)
        template <class Archive>
        struct ArchiveStoreImpl< Archive, std::vector<bool> > {
            static inline void store(const Archive& ar, const std::vector<bool>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL vector<bool>" << std::endl);
                std::size_t n = v.size();
                bool* b = new bool[n];
                for (std::size_t i=0; i<n; ++i) b[i] = v[i];
                ar & n & wrap(b,v.size());
                delete [] b;
            }
        };


        /// Deserialize STL vector<bool>. Clears & resizes as necessary.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, std::vector<bool> > {
            static void load(const Archive& ar, std::vector<bool>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL vector" << std::endl);
                std::size_t n = 0ul;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                bool* b = new bool[n];
                ar & wrap(b,v.size());
                for (std::size_t i=0; i<n; ++i) v[i] = b[i];
                delete [] b;
            }
        };

        /// Serialize STL string
        template <class Archive>
        struct ArchiveStoreImpl< Archive, std::string > {
            static void store(const Archive& ar, const std::string& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL string" << std::endl);
                ar & v.size();
                ar & wrap((const char*) &v[0],v.size());
            }
        };


        /// Deserialize STL string. Clears & resizes as necessary.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, std::string > {
            static void load(const Archive& ar, std::string& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL string" << std::endl);
                std::size_t n = 0ul;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                ar & wrap((char*) &v[0],n);
            }
        };


        /// (de)Serialize an STL pair.
        template <class Archive, typename T, typename Q>
        struct ArchiveSerializeImpl< Archive, std::pair<T,Q> > {
            static inline void serialize(const Archive& ar, std::pair<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize STL pair" << std::endl);
                ar & t.first & t.second;
            }
        };


        /// Serialize an STL map (crudely).
        template <class Archive, typename T, typename Q>
        struct ArchiveStoreImpl< Archive, std::map<T,Q> > {
            static void store(const Archive& ar, const std::map<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL map" << std::endl);
                ar << t.size();
                for (typename std::map<T,Q>::const_iterator p = t.begin();
                        p != t.end(); ++p) {
                    // Fun and games here since IBM's iterator (const or
                    // otherwise) gives us a const qualified key
                    // (p->first) which buggers up the type matching
                    // unless the user defines pair(T,Q) and pair(const
                    // T,Q) to have cookie (which is tedious).
                    std::pair<T,Q> pp = *p;
                    ar & pp;
                }
            }
        };


        /// Deserialize an STL map.  Map is NOT cleared; duplicate elements are replaced.
        template <class Archive, typename T, typename Q>
        struct ArchiveLoadImpl< Archive, std::map<T,Q> > {
            static void load(const Archive& ar, std::map<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL map" << std::endl);
                std::size_t n = 0;
                ar & n;
                while (n--) {
                    std::pair<T,Q> p;
                    ar & p;
                    t[p.first] = p.second;
                }
            }
        };



    }
}

#endif // MADNESS_WORLD_ARCHIVE_H__INCLUDED
