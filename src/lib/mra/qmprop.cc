#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/adquad.h>
#include <cmath>
#include <complex>
#include <constants.h>
#include <iostream>

/// \file qmprop.cc
/// \brief Implements BandlimitedPropagator and qm_free_particle_propagator


namespace madness {
    /// Class to evaluate the filtered Schrodinger free-particle propagator in real space

    /// Follows the corresponding Maple worksheet and the implementation notes.
    class BandlimitedPropagator {
    private:
        const double c;
        const double t;
        const double width;
        const double ctop;
        const double L;
        const double h;
        const int n;
        const double dc;

        std::complex<double> ff(double k) const {
            if (k>2.54*c) return std::complex<double>(0.0,0.0);
            const std::complex<double> arg(0,-k*k*t*0.5);
            return std::exp(arg)/(1.0+std::pow(k/c, 30.0));
        }

    public:
        typedef double_complex returnT;

        BandlimitedPropagator(double c, double t, double width)
                : c(c)
                , t(t)
                , width(width)
                , ctop(3.0*c)
                , L(sqrt(6.1*c*c*t*t + 82000.0/c/c))
                , h(3.14/ctop)
                , n(2.0*L/h+1)
                , dc(2*ctop/(n-1)) {
            //         std::cout << " c " << c << std::endl;
            //         std::cout << " t " << t << std::endl;
            //         std::cout << " ctop " << ctop << std::endl;
            //         std::cout << " L " << L << std::endl;
            //         std::cout << " h " << h << std::endl;
            //         std::cout << " n " << n << std::endl;
            //         std::cout << " dc " << dc << std::endl;
        }

        std::complex<double> operator()(double x) const {
            x = width*x;  // <<<<<<<< from simulation to user coords
            if (fabs(x) > L) return std::complex<double>(0.0,0.0);
            std::complex<double> base = exp(std::complex<double>(0.0,-x*ctop));
            std::complex<double>  fac = exp(std::complex<double>(0.0,x*dc));
            std::complex<double> sum(0.0,0.0);
            double W = -ctop;
            for (int i=0; i<n; i++,W+=dc) {
                //std::complex<double> arg(0.0,x*W);
                //std::complex<double> base = std::exp(arg);
                sum += ff(W)*base;
                base *= fac;
            }
            return width*sum*dc*0.5/madness::constants::pi;
        }

        static void test() {
            std::complex<double> maple(1.13851441120840,-.986104972277800);
            BandlimitedPropagator bp(31.4, 0.07, 1.0);
            //std::cout << bp(0.1) << " " << maple << std::endl;
            if (std::abs(bp(0.1)-maple) > 1e-11) throw "BandlimitedPropagator: failed test";
            return;
        }
    };

    template <int NDIM>
    SeparatedConvolution<double_complex,NDIM>
    qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep, double width) {
        BandlimitedPropagator::test();
        std::vector< SharedPtr< Convolution1D<double_complex> > > q(1);
        q[0] = SharedPtr< Convolution1D<double_complex> >(new GenericConvolution1D<double_complex,BandlimitedPropagator>(k,BandlimitedPropagator(bandlimit,timestep,width)));
        return SeparatedConvolution<double_complex,NDIM>(world, k, q, true);
	}		
		
	template <int NDIM>
    SeparatedConvolution<double_complex,NDIM>*
    qm_free_particle_propagatorPtr(World& world, int k, double bandlimit, double timestep, double width) {
        BandlimitedPropagator::test();
        std::vector< SharedPtr< Convolution1D<double_complex> > > q(1);
        q[0] = SharedPtr< Convolution1D<double_complex> >(new GenericConvolution1D<double_complex,BandlimitedPropagator>(k,BandlimitedPropagator(bandlimit,timestep,width)));
        return new SeparatedConvolution<double_complex,NDIM>(world, k, q, true);
    }

#ifdef FUNCTION_INSTANTIATE_1
    template SeparatedConvolution<double_complex,1> qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep, double width);    
	template SeparatedConvolution<double_complex,1>* qm_free_particle_propagatorPtr(World& world, int k, double bandlimit, double timestep, double width);
#endif

#ifdef FUNCTION_INSTANTIATE_2
    template SeparatedConvolution<double_complex,2> qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep, double width);
#endif

#ifdef FUNCTION_INSTANTIATE_3
    template SeparatedConvolution<double_complex,3> qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep, double width);
#endif

#ifdef FUNCTION_INSTANTIATE_4
    template SeparatedConvolution<double_complex,4> qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep, double width);
#endif
}
