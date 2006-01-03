#ifndef SEPOP_H
#define SEPOP_H

#include <map>

#include <madness_config.h>
#include <tensor/tensor.h>

namespace madness {
    /// Simplified interface to STL map for caching matrix elements
    class Cache {
        // Default copy, assignment and destructor should be satisfactory
    private:
        typedef std::map< unsigned long, Tensor<double> > mapT;
        mapT cache;
        const mapT::iterator end;
        
        /// Private:  turns (n,lx) into key
        inline unsigned long key1(long n, long lx) const {
            return unsigned( (n<<8) + (lx+127));
        };
        
        /// Private:  turns (n,lx,ly,lz) into key
        inline unsigned long key3(long n, long lx, long ly,  long lz) const {
            return unsigned( (n<<24) + ((lx+127)<<16) + ((ly+127)<<8) + (lz+127));
        };
        
        /// Private: If key is present return pointer to cached value, otherwise NULL
        Tensor<double>* getptr(unsigned long key) {
            mapT::iterator test = cache.find(key);
            if (test == end) return false;
            return &((*test).second);
        };
        
        /// Private: Set value associated with key
        inline void set(unsigned long key, const Tensor<double>& val) {
            cache.insert(std::pair< unsigned long, Tensor<double> >(key,val));
        };
        
    public:
        Cache() : cache(), end(cache.end()) {};
        
        Cache(const Cache& c) : cache(c.cache), end(cache.end()) {};
        Cache& operator=(const Cache& c) {
            if (this != &c) {
                cache.clear();
                cache = c.cache;
            }
            return *this;
        };
        
        /// If (n,lx) is present return pointer to cached value, otherwise return NULL
        
        /// Assumes that the translations are in range |lx| < 128
        inline Tensor<double>* getptr(long n,  long lx) {
            return getptr(key1(n,lx));
        };
        
        /// If (n,lx,ly,lz) is present return pointer to cached value, otherwise return NULL
        
        /// Assumes that the translations are in range |lx| < 128
        inline Tensor<double>* getptr(long n,  long lx, long ly,  long lz) {
            return getptr(key3(n,lx,ly,lz));
        };
        
        /// Set value associated with (n,lx)
        
        /// Assumes that the translations are in range |lx| < 128
        inline void set(long n, long lx, const Tensor<double>& val) {
            set(key1(n,lx),val);
        };
        
        /// Set value associated with (n,lx,ly,lz)
        
        /// Assumes that the translations are in range |lx| < 128
        inline void set(long n, long lx, long ly, long lz, const Tensor<double>& val) {
            set(key3(n,lx,ly,lz),val);
        };
    };
    
    class GaussianConvolution {
    public:
        double coeff;
        double expnt;
        int k;
        int npt;
        Tensor<double> c;
        Tensor<double> hgT;
        Tensor<double> quad_x;
        Tensor<double> quad_w;
        Cache rnlij_cache;
        Cache ns_cache;
        Cache ns_T_cache;
        
        GaussianConvolution() {};
        GaussianConvolution(int k, double coeff, double expnt);
        Tensor<double> rnlp(long n, long l);
        const Tensor<double>& rnlij(long n, long lx);
        const Tensor<double>& nonstandard(long n, long lx);
        Tensor<double>& nonstandard_T(long n, long lx);
        bool issmall(long n, long lx);
    };
    
    class SeparatedConvolution {
    public:
        const long k;
        const Tensor<double> coeff;
        const Tensor<double> expnt;
        long rank;
        std::vector<GaussianConvolution> ops;
        std::vector<double> signs;
        
        SeparatedConvolution(long k, 
                             const Tensor<double>& coeffs, 
                             const Tensor<double>& expnts);
        
        double munorm(long mu, long n, long x, long y, long z);
        double norm(long n, long x, long y, long z);
        void opxv(long n, long x, long y, long z,
                  const Tensor<double>& f, Tensor<double>& result,
                  double tol);
        void opxvt(long n, long x, long y, long z,
                   const Tensor<double>& f, Tensor<double>& result,
                   double tol);
        void muopxv(const long mu, 
                    const long n, const long x, const long y, const long z,
                    const Tensor<double>& f, Tensor<double>& result,
                    const double tol);
        void muopxvt(long mu, long n, long x, long y, long z,
                     const Tensor<double>& f, Tensor<double>& result,
                     double tol);
    };
    
}

#endif
