#include <madness/tensor/clapack.h>
#include <madness/fortran_ctypes.h>

/*!
  \file linalg_wrappers.h
  \brief Template wrappers for LAPACK routines 
  \ingroup linalg
@{
*/

namespace madness {
namespace detail {

    template <typename T>
    struct real_type {
      using type = T;
    };
  
    template <typename T>
    struct real_type< std::complex<T> > {
      using type = T;
    };

}

  
    /// Compute the SVD via LAPACK
    template <typename T>
    void svd( char jobu, char jobvt, integer m, integer n, T* A, integer lda,
              typename detail::real_type<T>::type* S, T* U, integer ldu,
              T* VT, integer ldvt ); 

    /// Solve the EVP via LAPACK
    template <typename T>
    void hereig( char jobz, char uplo, integer n, T* A, integer lda, 
                 typename detail::real_type<T>::type* W );

    /// Solve the Generalized EVP via LAPACK
    template <typename T>
    void hereig_gen( integer itype, char jobz, char uplo, integer n, T* A, 
                     integer lda, T* B, integer ldb, 
                     typename detail::real_type<T>::type* W );

    /// Compute the Cholesky Factorization via LAPACK
    template <typename T>
    void cholesky( char uplo, integer n, T* A, integer lda );





    /// Linear algebra Exception
    class LinAlgException : public std::exception {
        const char* msg;
        const char* assertion;
        int value;
        int line;
        const char *function;
        const char *filename;
  
public:
        LinAlgException(const char* s, const char *a, int err, 
                        int lin, const char *func, const char *file)
                : msg(s)
                , assertion(a)
                , value(err)
                , line(lin)
                , function(func)
                , filename(file) { }

        virtual const char* what() const throw() {
            return msg;
        }

        virtual ~LinAlgException() throw() {}

    friend std::ostream& operator <<(std::ostream& out, const LinAlgException& e) {
        out << "LinAlgException: msg='";
        if (e.msg) out << e.msg;
        out << "'\n";
        if (e.assertion) out << "                 failed assertion='" <<
            e.assertion << "'\n";
        out << "                 value=" << e.value << "\n";
        if (e.line) out << "                 line=" << e.line << "\n";
        if (e.function) out << "                 function='" <<
            e.function << "'\n";
        if (e.filename) out << "                 filename='" <<
            e.filename << "'\n";

        return out;
    }

    };



#define LINALG_STRINGIZE(X) #X
#define LINALG_EXCEPTION_AT(F, L) LINALG_STRINGIZE(F) "(" LINALG_STRINGIZE(L) ")"

#define LINALG_EXCEPTION(msg,value) \
    throw ::madness::LinAlgException("LINALG EXCEPTION: " LINALG_EXCEPTION_AT( __FILE__, __LINE__ ) ": " msg , \
    0,value,__LINE__,__FUNCTION__,__FILE__)

#define LINALG_ASSERT(condition,msg,value) \
do {if (!(condition)) \
        throw ::madness::LinAlgException("LINALG ASSERTION FAILED: " LINALG_EXCEPTION_AT( __FILE__, __LINE__ ) ": " msg , \
        #condition,value,__LINE__,__FUNCTION__,__FILE__); \
   } while (0)
}

/* @} */
