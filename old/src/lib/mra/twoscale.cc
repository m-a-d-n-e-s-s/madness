#include <iostream>
using std::cout;
using std::endl;
using std::ios;

#include <cstdio>
#include <cstdlib>

#include <cmath>
using std::abs;

#include <mra/twoscale.h>
#include <tensor/tensor.h>
#include <misc/misc.h>

/// \file twoscale.cc
/// \brief Routines to provide twoscale & correlation coeffs for Legendre basis

namespace madness {

    static const int kmax = 34;
    static const char *asciifile = "coeffs";
    static const double cksum[] = {
                                      1.57015282218922890,1.48602374574943450
                                  };

    static class twoscale_cache_class {
        /// This caches the two-scale coefficients
    public:
        Tensor<double> h0;
        Tensor<double> h1;
        Tensor<double> g0;
        Tensor<double> g1;
    }
    cache[kmax+1];

    static bool loaded = 0;

    static void checksum(int kmax, double *s0, double *s1) {
        double sum0=0, sum1=0;
        for (int k=0; k<=kmax; k++) {
            Tensor<double> h0 = cache[k].h0;
            Tensor<double> h1 = cache[k].h1;
            Tensor<double> g0 = cache[k].g0;
            Tensor<double> g1 = cache[k].g1;
            for (int i=0; i<k; i++) {
                for (int j=0; j<k; j++) {
                    double ij = double(j+k)/(i+k);
                    sum0 += ij*(std::abs(h0(i,j)) + std::abs(h1(i,j)) +
                                std::abs(g0(i,j)) + std::abs(g1(i,j)));
                    sum1 += ij*(h0(i,j) + h1(i,j) +
                                g0(i,j) + g1(i,j));
                    while (sum0 > 2.0) sum0 *= 0.25;
                    while (std::abs(sum1) > 2.0) sum1 *= 0.25;
                }
            }
        }
        *s0 = sum0;
        *s1 = sum1;
    }

    static Tensor<double> readmat(int k, FILE* file) {
        Tensor<double> a(k,k);
        for (int i=0; i<k; i++) {
            for (int j=0; j<k; j++) {
                double c;
                if (fscanf(file,"%lf",&c) != 1) {
                    cout << "readmat: twoscale missing coeff?\n";
                    throw "readmat";
                }
                a(i,j) = c;
            }
        }
        return a;
    }

    static inline double phase(long i) {
        return (i&1) ? -1.0 : 1.0;
    }

    static bool readascii(int kmax,const char *filename) {
        FILE* file = fopen(filename,"r");
        if (!file) {
            cout << "twoscale: failed opening file with twoscale coefficients\n";
            return false;
        }
        for (int k=0; k<kmax+1; k++) {
            Tensor <double> h0, g0;
            try {
                h0 = readmat(k,file);
                g0 = readmat(k,file);
            } catch (char *e) {
                fclose(file);
                return false;
            }

            Tensor<double> h1(k,k);
            Tensor<double> g1(k,k);

            for (int i=0; i<k; i++) {
                for (int j=0; j<k; j++) {
                    h1(i,j) = h0(i,j)*phase(i+j);
                    g1(i,j) = g0(i,j)*phase(i+j+k);
                }
            }

            cache[k].h0 = h0;
            cache[k].h1 = h1;
            cache[k].g0 = g0;
            cache[k].g1 = g1;
        }
        fclose(file);

        double sum0, sum1;
        checksum(kmax,&sum0,&sum1);
        if (std::abs(cksum[0]-sum0) > 1e-14 || std::abs(cksum[1]-sum1) > 3e-14) {
            cout.setf(ios::scientific);
            cout.precision(17);
            cout.width(30);
            cout << cksum[0] << " " << cksum[1] << endl;
            cout << sum0 << " " << sum1 << endl;
            return false;
        }
        //cout << "\ntwoscale: read valid twoscale coeffients from disk\n";

        loaded = true;
        return true;
    }


    /// Return the two scale coefficients in the Legendre basis

    /// Returns true on success, false on failure (and prints error message)
    bool two_scale_coefficients(int k,
                                Tensor<double>* h0, Tensor<double>* h1,
                                Tensor<double>* g0, Tensor<double>* g1) {
        if (!loaded) {
            if (!readascii(kmax,asciifile)) return false;
        }

        if (k < 1 || k > kmax) return false;

        *h0 = copy(cache[k].h0);
        *h1 = copy(cache[k].h1);
        *g0 = copy(cache[k].g0);
        *g1 = copy(cache[k].g1);

        return true;
    }

    bool two_scale_hg(int k, Tensor<double>* hg) {
        Tensor<double> h0(k,k), h1(k,k), g0(k,k), g1(k,k);

        if (!two_scale_coefficients(k, &h0, &h1, &g0, &g1)) return false;

        *hg = Tensor<double>(2*k,2*k);

        Slice sk(0,k-1), sk2(k,-1);
        (*hg) = Tensor<double>(2*k,2*k);
        (*hg)(sk,sk)   = h0;
        (*hg)(sk,sk2)  = h1;
        (*hg)(sk2,sk)  = g0;
        (*hg)(sk2,sk2) = g1;

        return true;
    }

    bool test_two_scale_coefficients() {
        /// Test the two scale coefficients for orthogonality
        for (int k=1; k<kmax; k++) {
            Tensor<double> hg;
            if (!two_scale_hg(k,&hg)) return false;
            Tensor<double> ident(2*k,2*k);
            for (int i=0; i<2*k; i++) ident(i,i) = 1.0;

            double err0 = (inner(hg,hg,0,0)-ident).absmax();
            if (err0 > 7e-16) {
                std::cout << "twoscale failed 0: " << k << " " << err0 << std::endl;
                std::cout << (inner(hg,hg,0,0)-ident);
                return false;
            }

            double err1 = (inner(hg,hg,1,1)-ident).absmax();
            if (err1 > 7e-16) {
                std::cout << "twoscale failed 1: " << k << " " << err1 << std::endl;
                std::cout << (inner(hg,hg,1,1)-ident);
                return false;
            }
        }
        return true;
    }


    // BELOW HERE THE AUTOCORRELATION ROUTINES

    static const char *filename = "autocorr";
    static const int kmax_autoc = 27;
    static int kread = -1;
    static Tensor<double> _c;

    static bool read_data(int k);

    /// Return the autocorrelation coefficients for scaling functions of given order

    /// Returned is view of a cached copy of the data ... do not modify.
    ///
    /// The autocorrelation functions are defined as
    /// \code
    ///    Phi_ij(z) = int(0,z+1) phi_i(x) phi_j(x-z)  z<=0
    ///    Phi_ij(z) = int(z,1) phi_i(x) phi_j(x-z)    z>=0
    /// \endcode
    /// and are expanded in the double order Legendre scaling functions on either
    /// side of the origin
    /// \code
    ///    Phi_ij(z) = sum(p) [phi_p(z+1)*cminus_ijp + phi_p(z)*cplus_ijp]
    ///
    ///    cplus_ijp  = int(-1,0) Phi_ij(z) phi_p(z+1)
    ///    cminus_ijp = int(0,1)  Phi_ij(z) phi_p(z)
    /// \endcode
    ///
    /// The expansion coefficients \c cminus and \c cplus have been precomputed
    /// with Maple to a precision of about 1e-30.  Only \c cplus is stored and we use
    /// \code
    ///    cplus(i,j,p) = (-1)**(i+j+p) * cminus(i,j,p)
    ///    cplus(i,j,p) = (-1)**(i+j)   * cplus(j,i,p)
    /// \endcode
    ///
    /// The returned tensor concatenates cminus and cplus.
    ///
    /// Return true on success, false on failure (which also prints error
    /// message to stdout).
    bool autoc(int k, Tensor<double>* c) {
        if (k < 1 || k > kmax_autoc) {
            cout << "autoc: invalid k " << k << endl;
            return false;
        }

        // The data is cached in _c ... reread for each new value of k
        if (kread != k) {
            if (!read_data(k))  return false;
        }

        *c = _c; //copy(_c(Slice(0,k-1),Slice(0,k-1),Slice(0,4*k-1)));
        return true;
    }

    /// Return +1 if i is even, -1 if i is odd ... perverse no if-test form
    static inline int parity(int i) {
        return 1 - ((i&1)<<1);
    }

    bool test_autoc() {
        unsigned long correct = 0x638a9b;
        unsigned long computed = checksum_file(filename);
        if (correct != computed)
            cout << "test_autoc: file checksum invalid: correct="
            << correct << " computed=" << computed << endl;

        return (correct == computed);
    }

    static bool read_data(int k) {
        if (!test_autoc()) return false;
        kread = -1;
        FILE *file = fopen(filename,"r");
        if (!file) {
            cout << "autoc: failed opening file with autocorrelation coefficients" << endl;
            return false;
        }

        _c = Tensor<double>(k,k,4*k);

        long twok = 2*k;
        while (1) {
            long i, j, p;
            double val;
            if (fscanf(file,"%ld %ld %ld %lf",&i,&j,&p,&val) != 4) {
                cout<<"autoc: failed reading file " << endl;
                fclose(file);
                return false;
            }
            if (i >= k) break;
            double ij = parity(i+j);
            double ijp= parity(i+j+p);

            _c(i,j,p)      = val*ijp;	// c-
            _c(j,i,p)      = val*ij*ijp;
            _c(i,j,p+twok) = val;	// c+
            _c(j,i,p+twok) = val*ij;
        }
        fclose(file);
        kread = k;
        return true;
    }

    /// Collective routine to load and cache twoscale & autorrelation coefficients

    /// Only process rank 0 will access the files.
    void load_coeffs(Communicator& comm) {
        if (!loaded) {
            int autoc_k = 15;   // Plausible maximum value
            if (comm.rank() == 0) {
                if (!readascii(kmax,asciifile))
                    throw "load_coeffs: failed reading twoscale coeffs";
                if (!read_data(15))
                    throw "load_coeffs: failed reading coeffs";
            } else {
                for (int k=0; k<=kmax; k++) {
                    cache[k].h0 = Tensor<double>(k,k);
                    cache[k].h1 = Tensor<double>(k,k);
                    cache[k].g0 = Tensor<double>(k,k);
                    cache[k].g1 = Tensor<double>(k,k);
                }
                _c = Tensor<double>(autoc_k,autoc_k,4*autoc_k);
            }
            for (int k=0; k<=kmax; k++) {
                comm.Bcast(cache[k].h0.ptr(), k*k, 0);
                comm.Bcast(cache[k].h1.ptr(), k*k, 0);
                comm.Bcast(cache[k].g0.ptr(), k*k, 0);
                comm.Bcast(cache[k].g1.ptr(), k*k, 0);
            }

            comm.Bcast(_c.ptr(), autoc_k*autoc_k*4*autoc_k, 0);

            loaded = true;
        }
    }

}



