#ifndef MAD_ARCHIVE_H
#define MAD_ARCHIVE_H

/// \file archive.h
/// \brief Interface templates for the archives (serialization)

/// The user should not need to include this file directly.  Instead,
/// include the header file for the actual archive (binary file, text/xml
/// file, vector in memory, ...) that you want to use.
///
/// The interface and implementation are deliberately modelled, albeit
/// loosely, upon the Boost serialization class (thanks boost!).  At
/// some point we might switch, but for now we cannot.  The major
/// differences are that this archive class does \em not break cycles
/// and does \em not automatically store unique copies of data
/// referenced by multiple objects.  Also, classes are responsbible
/// for managing their own version information.  At the lowest level,
/// the interface to an archive also differs to facilitate
/// vectorization and high-bandwidth data transfer.  The
/// implementation employs templates that are almost entirely inlined.
/// This should enable low-overhead use of archives in applications
/// such as interprocess communication.
///
/// How to use an archive?  An archive is a uni-directional stream of
/// typed data usually to disk, memory, or another process.  Whether
/// the stream is for input or for output, you can use the \c &
/// operator to transfer data to the stream.  If you really want, you
/// can also use the \c << and \c >> for output or input,
/// respectively, but there is no reason to do so (I don't know why
/// Boost did this).  The \c & operator chains just like \c << for \c
/// cout or \c >> for \c cin.  Unless type checking has not been
/// implemented by an archive for reasons of efficiency (e.g., message
/// passing) a C-string exception will be thrown on a type-mismatch
/// when deserializing.  End-of-file, out-of-memory and other others
/// also generate string exceptions.
///
/// You may discover in \c archive.h other interfaces but you should
/// \em not use them --- use the \& operator!  The lower level
/// interfaces will probably not, or only inconsistently, incorpoate
/// type information and may even appear to work when they are not.
///
/// Fundamental types (see below), STL complex, vector, strings, pairs
/// and maps, and madness tensors (int, long, float, double,
/// float_complex, double_complex) all just work without you doing
/// anything, as do fixed dimension arrays of the same.  E.g.,
/// \code
/// map<int,double> fred;
/// int info[3];
/// bool finished;
/// map[0]=55.0; map[1]=99.0;
/// info[0]=1; info[1]=33; info[2]=-1;
/// BinaryFstreamOutputArchive ar('restart.dat');
/// ar & map & info & finished;
/// \endcode
/// Deserializing is identical, except that you need to use an input archive.
/// \code
/// map<int,double> fred;
/// int info[3];
/// bool finished;
/// BinaryFstreamInputArchive ar('restart.dat');
/// ar & map & info & finished;
/// \endcode
///
/// Variable dimension and dynamically allocated arrays do not have
/// their dimension encoded in their type.  The best way to 
/// (de)serialize them is to wrap them in an \c archive_array as follows.
/// \code
/// int a[n]; // n is not known at compile time
/// double *p = new double[n];
/// ar & wrap(a,n) & wrap(p,n);
/// \endcode 
/// The \c wrap function template is a factory function to simplify
/// instantiation of a correctly typed \c archive_array template.
/// Note that when deserializing you must have first allocated the
/// array --- the above code can be used for both serializing and
/// deserializing.  If you want the memory to be automatically allocated
/// consider using either an STL vector or a madness tensor.
///
/// User-defined types require a little more effort.  Three
/// cases are distinguished.  
/// - symmetric load and store
///    - intrusive
///    - non-intrusive
/// - non-symmetric load and store
/// We will examine each in turn, but we first need to discuss a little
/// about the implementation.
///
/// When transfering an object \c obj to/from an archive \c ar with \c ar&obj,
/// you are invoking the templated function 
/// \code
/// template <class Archive, class T>
/// inline const Archive& operator&(const Archive& ar, T& obj);
/// \endcode
/// which then invokes other templated functions to redirect to input
/// or output streams as appropriate, manage type checking, etc..
/// What we'd now like to do is essentially overload the behaviour of
/// some of these functions in order accomodate your fancy object.
/// However, function templates are tricky buggers when it comes to
/// overloading.  Rather than do this (or more precisely we did and
/// got burnt), we are adopting the technique recommend in
/// http://www.gotw.ca/publications/mill17.htm (look for moral#2).
/// Each of the templated functions directly calls a member of a
/// templated class.  Classes, unlike functions, can be partially
/// specialized so it is easy to control and predict what is
/// happening.  Thus, in order to change the behaviour of all archives
/// for an object you just have to provide a partial specialization of
/// the appropriate class(es).  Do \em not overload any of the
/// function templates.
///
/// Many classes can use the same code for serializing and
/// deserializing.  If such a class can be modified, the cleanest way
/// of enabling serialization is to add a templated method as follows.
/// \code
/// class A {
///     float a;
/// public:
///     A(float a = 0.0) : a(a) {};
///  
///     template <class Archive>
///     inline void serialize(const Archive& ar) {ar & a;}
/// };
/// \endcode
/// 
/// If a class with symmetric serialization cannot be modified, then
/// you can define an external class template with this signature in
/// the \c madness::archive namespace (where \c Obj is the name of
/// your type).
/// \code
/// template <class Archive>
/// struct ArchiveSerializeImpl<Archive,Obj> {
///       static inline void serialize(const Archive& ar, Obj& obj);
/// };
/// \endcode
///
/// For example,
/// \code
/// class B {
/// public:
///     bool b;
///     B(bool b = false) : b(b) {};
/// };
///
/// namespace madness {
///     namespace archive {
///         template <class Archive>
///         struct ArchiveSerializeImpl<Archive,B> {
///             static inline void serialize(const Archive& ar, B& b) {ar & b.b;};
///         };
///     }
/// }
/// \endcode
///
/// For classes that do not have symmetric (de)serialization you must
/// define separate partial templates for load and store,
/// respectively, with these signatures and again in the \c
/// madness::archive namespace.
/// \code
/// template <class Archive> 
/// struct ArchiveLoadImpl<Archive,Obj> {
///     static inline void load(const Archive& ar, Obj& obj);
/// };
/// 
/// template <class Archive> 
/// struct ArchiveStoreImpl<Archive,Obj> {
///     static inline void store(const Archive& ar, Obj& obj);
/// };
/// \endcode
///
/// First a simple, but artificial example.
/// \code
/// class C {
/// public:
///   long c;
///   C(long c = 0) : c(c) {};
/// };
/// 
/// template <class Archive>
/// inline void store(const Archive& ar, const C& c) {ar & c.c;}
/// 
/// template <class Archive>
/// inline void load(const Archive& ar, C& c) {ar & c.c;}
/// \endcode
///
/// Now a more complicated example that genuinely requires asymmetric
/// load and store.  First, a class definition for a simple linked list.
/// \code
/// class linked_list {
///   int value;
///   linked_list *next;
/// public:
///   linked_list(int value = 0) : value(value), next(0) {};
/// 
///   void append(int value) {
///     if (next) next->append(value);
///     else next = new linked_list(value);
///   };
/// 
///   void set_value(int val) {value = val;};
/// 
///   int get_value() const {return value;};
/// 
///   linked_list* get_next() const {return next;};
/// };
/// \endcode
/// And this is how you (de)serialize it.
/// \code
/// namespace madness {
///     namespace archive {
///         template <class Archive>
///         struct ArchiveStoreImpl<Archive,linked_list> {
///             static void store(const Archive& ar, const linked_list& c) {
///                 ar & c.get_value() & bool(c.get_next());
///                 if (c.get_next()) ar & *c.get_next();
///             };
///         };
///         
///         template <class Archive>
///         struct ArchiveLoadImpl<Archive,linked_list> {
///             static void load(const Archive& ar, linked_list& c) {
///                 int value;  bool flag;
///                 ar & value & flag;
///                 c.set_value(value);
///                 if (flag) {
///                     c.append(0);
///                     ar & *c.get_next();
///                 }
///             };
///         };
///     }
/// }
/// \endcode
///
/// Given the above implementation of a linked list, you can 
/// (de)serialize an entire list using a single statement.
/// \code
/// linked_list list(0);
/// for (int i=1; i<=10; i++) list.append(i);
/// BinaryFstreamOutputArchive ar('list.dat'); 
/// ar & list;
/// \endcode;
///
/// Presently provided are archives mapped to binary or text (XML)
/// files by wrapping an instance of the \c fstream class, and an STL
/// \c vector<unsigned char> for storing an archive in memory.  The \c
/// vector archive is bitwise identical to the binary file archive.
/// Soon to come will be an efficient mapping to MPI for passing data
/// between processes, but this can presently be acomplished by
/// serializing to a vector, sending the vector, and then deserializing.
/// 
/// Minimally, an archive must derive from either BaseInputArchive or
/// BaseOutputArchive and define for arrays of fundamental types
/// either a \c load or \c store method, as appropriate.  Additional
/// methods can be provided to manipulate the target stream.  Here is
/// a simple, but functional, implementation of a binary file archive.
/// \code
/// #include <fstream>
/// using namespace std;
/// class OutputArchive : public BaseOutputArchive {
///   mutable ofstream os;
/// public:
///   OutputArchive(const char* filename) 
///     : os(filename, ios_base::binary | ios_base::out | ios_base::trunc) {};
///   
///   template <class T>
///   void store(const T* t, long n) const {
///     os.write((const char *) t, n*sizeof(T));
///   }
/// };
/// 
/// class InputArchive : public BaseInputArchive {
///   mutable ifstream is;
/// public:
///   InputArchive(const char* filename) 
///     : is(filename, ios_base::binary | ios_base::in) {};
///   
///   template <class T>
///   void load(T* t, long n) const {
///     is.read((char *) t, n*sizeof(T));
///   }
/// };
/// \endcode

#include <complex>
#include <iostream>
#include <vector>
#include <map>
#include <typestuff.h>
#include <tensor/tensor.h>

//#ifndef OCTTREE_H
#include <octtree/octtree.h>
//#endif
//#include <mra/mra.h>

#define ARCHIVE_COOKIE "archive"
#define ARCHIVE_MAJOR_VERSION 0
#define ARCHIVE_MINOR_VERSION 1

//#define MAD_ARCHIVE_DEBUG_ENABLE

#ifdef MAD_ARCHIVE_DEBUG_ENABLE
#define MAD_ARCHIVE_DEBUG(s) s
using std::endl;
#else
#define MAD_ARCHIVE_DEBUG(s) 
#endif

namespace madness {
    
    namespace archive {
        
      // There are 64 empty slots for user types.  Free space for
      // registering user types begins at cookie=128.
        
#ifdef MAD_ARCHIVE_TYPE_NAMES_CC
        char *archive_type_names[256];
#else
        extern char *archive_type_names[256];
#endif
        void archive_initialize_type_names();
        
        template <typename T>
        struct archive_typeinfo {
            static const unsigned char cookie = 255;
        };
        
#if defined(ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE) && defined(ARCHIVE_REGISTER_TYPE_IBMBUG)
#define ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T) \
        ; const unsigned char archive_typeinfo< T >::cookie
#else
#define ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T)
#endif

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
        throw "archive_register_type_name: slot/cookie already in use!"; \
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
        
        ARCHIVE_REGISTER_TYPE_AND_PTR(OctTree<double>,38);
        
        /// Base class for all archives
        class BaseArchive {
        public:
            static const bool is_archive = true;
            static const bool is_input_archive = false;
            static const bool is_output_archive = false;
            BaseArchive() {
                archive_initialize_type_names();
            };
        };
        
        
        /// Base class for input archives 
        class BaseInputArchive : public BaseArchive {
        public:
            static const bool is_input_archive = true;
        };
        
        
        /// Base class for output archives 
        class BaseOutputArchive : public BaseArchive {
        public:
            static const bool is_output_archive = true;
        };
        

        // Serialize an array of fundamental stuff
        template <class Archive, class T>
        typename madness::enable_if< madness::type_and_c< madness::is_fundamental<T>::value, 
                                                          Archive::is_output_archive >, void >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize fund array" << std::endl);
            ar.store(t,n);
        }
        
        
        // Deserialize an array of fundamental stuff
        template <class Archive, class T>
        typename madness::enable_if< madness::type_and_c< madness::is_fundamental<T>::value, 
                                                          Archive::is_input_archive >, void >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "deserialize fund array" << std::endl);
            ar.load((T*) t,n);
        }
        
        
        // (de)Serialize an array of non-fundamental stuff
        template <class Archive, class T>
        typename madness::enable_if< madness::type_and_c< !madness::is_fundamental<T>::value, 
                                                          Archive::is_archive >, void >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize non-fund array" << std::endl);
            for (unsigned int i=0; i<n; i++) ar & t[i];
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
                    throw msg;
                }
                else {
                    MAD_ARCHIVE_DEBUG(std::cout << "read cookie " << archive_type_names[cookie] << std::endl);
                }
            };
        
        
            /// Serialize a cookie for type checking
            static inline void preamble_store(const Archive& ar) {
                unsigned char ck = archive_typeinfo<T>::cookie;
                ar.store(&ck, 1); // cannot use <<
                MAD_ARCHIVE_DEBUG(std::cout << "wrote cookie " << archive_type_names[ck] << std::endl);
            };
        
        
            /// By default there is no postamble
            static inline void postamble_load(const Archive& ar) {};

            /// By default there is no postamble
            static inline void postamble_store(const Archive& ar) {};
        };


        /// Default symmetric serialization of a non-fundamental thingy
        template <class Archive, class T>
        struct ArchiveSerializeImpl {
            static inline void serialize(const Archive& ar, T& t) {t.serialize(ar);};
        };


        // Redirect \c serialize(ar,t) to \c serialize(ar,&t,1) for fundamental types
        template <class Archive, class T>
        inline
        typename madness::enable_if< madness::type_and_c< madness::is_fundamental<T>::value, 
                                                          Archive::is_archive >, 
                                     void >::type
        serialize(const Archive& ar, const T& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize(ar,t) -> serialize(ar,&t,1)" << std::endl);
            serialize(ar,&t,1);
        }
        

        // Redirect \c serialize(ar,t) to \c ArchiveSerializeImpl for non-fundamental types
        template <class Archive, class T>
        inline
        typename madness::enable_if< madness::type_and_c< !madness::is_fundamental<T>::value, 
                                                          Archive::is_archive >, 
                                     void >::type
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
            };
        };

/*
        /// Default load of a thingy via serialize(ar,t)
        template <class Archive, class T> 
        struct ArchiveLoadImpl {
	    template <class U>
            static inline void load(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "load(ar,t) default" << std::endl);
                serialize(ar,t);
            };
        };
*/

        /// Default load of a thingy via serialize(ar,t)
        template <class Archive, class T> 
        struct ArchiveLoadImpl {
            static inline void load(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "load(ar,t) default" << std::endl);
                serialize(ar,t);
            };
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
            };
            
            static inline const Archive& wrap_load(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for default" << std::endl);
                ArchivePrePostImpl<Archive,T>::preamble_load(ar);
                ArchiveLoadImpl<Archive,T>::load(ar,(T&) t);  // Loses constness here!
                ArchivePrePostImpl<Archive,T>::postamble_load(ar);
                return ar;
            };
        };        


        // Redirect \c << to ArchiveImpl::wrap_store for output archives
        template <class Archive, class T>
        inline
        typename madness::enable_if_c<Archive::is_output_archive, const Archive&>::type
        operator<<(const Archive& ar, const T& t) {
            return ArchiveImpl<Archive,T>::wrap_store(ar,t);
        }
        
        // Redirect \c >> to ArchiveImpl::wrap_load for input archives
        template <class Archive, class T>
        inline
        typename madness::enable_if_c<Archive::is_input_archive, const Archive&>::type
        operator>>(const Archive& ar, const T& t) {
            return ArchiveImpl<Archive,T>::wrap_load(ar,t);
        }

        // Redirect \c & to ArchiveImpl::wrap_store for output archives
        template <class Archive, class T>
        inline
        typename madness::enable_if_c<Archive::is_output_archive, const Archive&>::type
        operator&(const Archive& ar, const T& t) {
            return ArchiveImpl<Archive,T>::wrap_store(ar,t);
        }
        
        // Redirect \c & to ArchiveImpl::wrap_load for input archives
        template <class Archive, class T>
        inline
        typename madness::enable_if_c<Archive::is_input_archive, const Archive&>::type
        operator&(const Archive& ar, const T& t) {
            return ArchiveImpl<Archive,T>::wrap_load(ar,t);
        }

        
        ///////////////////////////////////////////////////////////////
        
        
        
        /// Wrapper for dynamic arrays and pointers
        template <class T>
        class archive_array {
        public:
            const T* ptr;
            unsigned int n;
            inline
            archive_array(const T *ptr, unsigned int n) : ptr(ptr), n(n) {};
        };
        
        
        /// Factory function to wrap pointer as archive_array
        template <class T>
        inline
        archive_array<T> wrap(const T* ptr, unsigned int n) {
            return archive_array<T>(ptr,n);
        }
        
        
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
            };
            
            static inline const Archive& wrap_load(const Archive& ar, const archive_array<T>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for archive_array" << std::endl);
                ArchivePrePostImpl<Archive,T*>::preamble_load(ar);
                //unsigned int n;
                //ar >> n;
                //if (n != t.n) 
                //    throw "deserializing archive_array: dimension mismatch";
                //ArchivePrePostImpl<Archive,T>::preamble_load(ar);
                serialize(ar,(T *) t.ptr,t.n);
                //ArchivePrePostImpl<Archive,T>::postamble_load(ar);
                ArchivePrePostImpl<Archive,T*>::postamble_load(ar);
                return ar;
            };
        };
        
        
        /// Partial specialization for fixed dimension array redirects to archive_array
        template <class Archive, class T, std::size_t n> 
        struct ArchiveImpl<Archive, T[n]> {
            static inline const Archive& wrap_store(const Archive& ar, const T (&t)[n]) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for array" << std::endl);
                ar << wrap(&t[0],n);
                return ar;
            };
            
            static inline const Archive& wrap_load(const Archive& ar, const T (&t)[n]) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for array" << std::endl);
                ar >> wrap(&t[0],n);
                return ar;
            };
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
		T r, i;
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
            };
        };
        
        
        /// Deserialize STL vector. Clears & resizes as necessary.
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, std::vector<T> > {
            static void load(const Archive& ar, std::vector<T>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL vector" << std::endl);
                std::size_t n;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                ar & wrap((T *) &v[0],n);
            };
        };
        
        /// Serialize STL vector<bool> (as plain array of bool)
        template <class Archive>
        struct ArchiveStoreImpl< Archive, std::vector<bool> > {
            static inline void store(const Archive& ar, const std::vector<bool>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL vector<bool>" << std::endl);
                std::size_t n = v.size();
                bool* b = new bool[n];
                for (int i=0; i<n; i++) b[i] = v[i];
                ar & n & wrap(b,v.size());
                delete [] b;
            };
        };
        
        
        /// Deserialize STL vector<bool>. Clears & resizes as necessary.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, std::vector<bool> > {
            static void load(const Archive& ar, std::vector<bool>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL vector" << std::endl);
                std::size_t n;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                bool* b = new bool[n];
                ar & wrap(b,v.size());
                for (int i=0; i<n; i++) v[i] = b[i];
                delete [] b;
            };
        };
        
        /// Serialize STL string
        template <class Archive>
        struct ArchiveStoreImpl< Archive, std::string > {
            static void store(const Archive& ar, const std::string& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL string" << std::endl);
                ar & v.size();
                ar & wrap((const char*) &v[0],v.size());
            };
        };
        
        
        /// Deserialize STL string. Clears & resizes as necessary.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, std::string > {
            static void load(const Archive& ar, std::string& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL string" << std::endl);
                std::size_t n;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                ar & wrap((char*) &v[0],n);
            };
        };
        
        
        /// (de)Serialize an STL pair.
        template <class Archive, typename T, typename Q>
        struct ArchiveSerializeImpl< Archive, std::pair<T,Q> > {
            static inline void serialize(const Archive& ar, std::pair<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize STL pair" << std::endl);
                ar & t.first & t.second;
            };
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
            };
        };
        
        
        /// Deserialize an STL map.  Map is NOT cleared; duplicate elements are replaced.
        template <class Archive, typename T, typename Q>
        struct ArchiveLoadImpl< Archive, std::map<T,Q> > {
            static void load(const Archive& ar, std::map<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL map" << std::endl);
                std::size_t n;  ar & n;
                while (n--) {
                    std::pair<T,Q> p; ar & p;
                    t[p.first] = p.second;
                }
            };
        };
        
        
        /// Serialize a tensor
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, Tensor<T> > {
            static inline void store(const Archive& s, const Tensor<T>& t) {
                if (t.iscontiguous()) s & t.id & t.ndim & t.dim & wrap(t.ptr(),t.size); 
                else s & copy(t);
            };
        };
        
        
        /// Deserialize a tensor ... existing tensor is replaced
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, Tensor<T> > {
            static inline void load(const Archive& s, Tensor<T>& t) { 
                long id;
                s & id;
		if (id != t.id) throw "type mismatch deserializing a tensor";
                long ndim, dim[TENSOR_MAXDIM];  // Uncool reference to this macro
                s & ndim & dim;
                t = Tensor<T>(ndim, dim, false);
		t.fillrandom();
                s & wrap(t.ptr(), t.size);
            };
        };

// This disaster brought to you by hqi
	/// Serialize an OctTree<T>
	template <class Archive, class T>
	struct ArchiveStoreImpl< Archive, OctTree<T> > {
	    static inline void store(const Archive& ar, const OctTree<T>& t) {
//		std::cout << "serializing OctTree" << std::endl;
		ar & t._x & t._y & t._z & t._n & t._remote & t._rank & t._cost;
//		std::cout << "sent x y z = (" << t._x << "," << t._y << "," << t._z
//			<< ") ... cost" << std::endl;
		if (t.islocal())
		{
		    ar & t._data;
//		    std::cout << "sent data" << std::endl;
		    if (t._c[0][0][0])
		    {
		    	ar & 1;
//			std::cout << "t is a parent" << std::endl;

		        FORIJK(
		    	    if (t._c[i][j][k])
		    	    {
			    	store(ar, *(t._c[i][j][k]));
		    	    }
		        );

		    }
		    else
		    {
		    	ar & 0;
//			std::cout << "t is not a parent" << std::endl;
		    }
		}
	    };
	};


	/// Deserialize an OctTree<T>
	template <class Archive, class T>
	struct ArchiveLoadImpl< Archive, OctTree<T> > {
	    static inline void load(const Archive& ar, OctTree<T>& t) {
//		std::cout << "deserializing OctTree" << std::endl;
		ar & t._x & t._y & t._z & t._n & t._remote & t._rank & t._cost;
//		std::cout << "received x y z = (" << t._x << "," << t._y << "," << t._z
//			<< ") ... cost" << std::endl;
		if (!(t._remote))
		{
		    ar & t._data;
//		    std::cout << "received data" << std::endl;
		    int hasChildren;
//		    std::cout << "about to receive hasChildren" << std::endl;
		    ar & hasChildren;
//		    std::cout << "has children: " << hasChildren << std::endl;

		    if (hasChildren)
		    {
		    	FORIJK(
//			    std::cout << "load c[" << i << "," << j << "," << k << "]" << std::endl;
			    OctTree<T> *child = new OctTree<T>();
			    t._c[i][j][k] = shared_ptr<OctTree<T> >(child);
			    load(ar, *child);
//			    std::cout << "loaded c[" << i << "," << j << "," << k << "], (" <<
//				(t._c[i][j][k])->x() << ", " << (t._c[i][j][k])->y() << ", " <<
//				(t._c[i][j][k])->z() << ")"<< std::endl;
		    	);
		    }

		}
//		std::cout << "end of load (prolly won't see this)" << std::endl;
	    };
	};

	/// Serialize a Tensor thru a BaseTensor pointer
        template <class Archive>
        struct ArchiveStoreImpl< Archive, BaseTensor* > {
            static inline void store(const Archive& s, const BaseTensor* t) {
//		std::cout << "serizialing thru bt\n";
		s & t->id;
		if (t->id == TensorTypeData<double>::id) {
//		    std::cout << "serizialing thru bt ... it's a double!\n";
		    s & *((const Tensor<double>*) t);
		}
		else {
		    throw "not yet";
		}
            };
        };

	/// Deserialize a Tensor thru a BaseTensor pointer

	/// It allocates a NEW tensor ... the original point is stomped on
	/// and so should not be preallocated.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, BaseTensor* > {
            static inline void load(const Archive& s, BaseTensor*& t) { 
                long id;
                s & id;
//		std::cout << "deserizialing thru bt\n";
		if (id == TensorTypeData<double>::id) {
//		    std::cout << "deserizialing thru bt ... it's a double!\n";
		    Tensor<double>* d = new Tensor<double>();
		    s & *d;
		    t = d;
		}
		else {
		    throw "not yet";
		}
            };
	};
        
    }
}

#endif
