#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/adquad.h>
#include <misc/cfft.h>
#include <misc/interpolation_1d.h>
#include <cmath>
#include <complex>
#include <constants.h>
#include <iostream>
#include <vector>

/// \file qmprop.cc
/// \brief Implements BandlimitedPropagator and qm_free_particle_propagator


namespace madness {
    /// Class to evaluate the filtered Schrodinger free-particle propagator in real space

    /// Follows the corresponding Maple worksheet and the implementation notes.
    class BandlimitedPropagator {
    private:
        double width;
        double xmax;
        CubicInterpolationTable<double_complex> fit;

        double_complex g0_filtered(double k, double c, double t) {
            double r = k/c;
            if (fabs(r) > 4.0) return 0.0;
            double_complex arg(0.0, -k*k*t*0.5);
            return exp(arg)/(1.0 + pow(r,30.0));
        }

    public:
        typedef double_complex returnT;

        BandlimitedPropagator(double c, double t, double width) : width(width) {
            // 4.0 is range of window in Fourier space, 32x for
            // oversampling in real space for accurate cubic
            // interpolation
            const double kmax = 32.0*4.0*c;
            const double fac = kmax/constants::pi;
            const int N = 131072;
            const double hk = 2.0*kmax/N;
            const double hx = constants::pi/kmax;
            std::vector<double_complex> s(N);

            for (int i=0; i<N/2; i++) {
                double k = i*hk;
                s[i] = g0_filtered(k, c, t)*fac;
                if (i) s[N-i] = g0_filtered(-k, c, t)*fac;
            }

            CFFT::Inverse(&s[0], N);

            int n;
            for (n=N/3; n>=0 && std::abs(s[n])<1e-14; n--);
            n++;

            s.resize(n);

            xmax = (n-1)*hx;
            fit = CubicInterpolationTable<double_complex>(0.0, xmax, n, s);

            print("QM", c, t, width);
            std::cout.precision(12);
            print(0.001,(*this)(0.001));
            print(0.002,(*this)(0.002));
//             for (int i=0; i<10001; i++) {
//                 double x = i*xmax/10000.0/width;
//                 double_complex value = (*this)(x);
//                 print(x,value.real(),value.imag());
//             }
        }

        std::complex<double> operator()(double x) const {
            x = fabs(x)*width;
            if (x >= xmax) return 0.0;
            else return fit(x)*width;
        }

        static void test() {
            std::complex<double> maple(1.138514411208581,-0.986104972271240);
            BandlimitedPropagator bp(31.4, 0.07, 1.0);
            if (std::abs(bp(0.1)-maple) > 1e-11) {
                std::cout.precision(14);
                std::cout << bp(0.1) << " " << maple << " " << bp(0.1)-maple << std::endl;
                throw "BandlimitedPropagator: failed test";
            }
            return;
        }
    };

    Convolution1D<double_complex>*
    qm_1d_free_particle_propagator(int k, double bandlimit, double timestep, double width) {
        return new GenericConvolution1D<double_complex,BandlimitedPropagator>(k,BandlimitedPropagator(bandlimit,timestep,width));
    }

    template <int NDIM>
    SeparatedConvolution<double_complex,NDIM>
    qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep) {
        double width = FunctionDefaults<NDIM>::get_cell_min_width(); // Assuming cubic so all dim equal
        std::vector< SharedPtr< Convolution1D<double_complex> > > q(1);
        q[0] = SharedPtr< Convolution1D<double_complex> >(qm_1d_free_particle_propagator(k, bandlimit, timestep, width));
        return SeparatedConvolution<double_complex,NDIM>(world, k, q, true);
    }		
		
    template <int NDIM>
    SeparatedConvolution<double_complex,NDIM>*
    qm_free_particle_propagatorPtr(World& world, int k, double bandlimit, double timestep) {
        double width = FunctionDefaults<NDIM>::get_cell_min_width(); // Assuming cubic so all dim equal
        std::vector< SharedPtr< Convolution1D<double_complex> > > q(1);
        q[0] = SharedPtr< Convolution1D<double_complex> >(qm_1d_free_particle_propagator(k, bandlimit, timestep, width));
        return new SeparatedConvolution<double_complex,NDIM>(world, k, q, true);
    }

#ifdef FUNCTION_INSTANTIATE_1
    template SeparatedConvolution<double_complex,1> qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep);
	template SeparatedConvolution<double_complex,1>* qm_free_particle_propagatorPtr(World& world, int k, double bandlimit, double timestep);
#endif

#ifdef FUNCTION_INSTANTIATE_2
    template SeparatedConvolution<double_complex,2> qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep);
#endif

#ifdef FUNCTION_INSTANTIATE_3
    template SeparatedConvolution<double_complex,3> qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep);
#endif

#ifdef FUNCTION_INSTANTIATE_4
    template SeparatedConvolution<double_complex,4> qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep);
#endif
}
