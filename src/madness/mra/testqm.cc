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


  $Id: test.cc 257 2007-06-25 19:09:38Z HartmanBaker $
*/

/// \file testsuite.cc
/// \brief The QA/test suite for Function

#include <madness/mra/mra.h>
#include <unistd.h>
#include <cstdio>
#include <madness/constants.h>
#include <madness/mra/qmprop.h>

using namespace madness;

typedef Convolution1D<double_complex> complex_operatorT;

double_complex psi0(const coord_1d& r) {
    return double_complex(exp(-r[0]*r[0]*0.5),0.0);
}


double mask(const coord_1d& r) {
    return 1.0/(1.0 + std::pow(r[0]*0.05,30.0));
}

double_complex V(const coord_1d& r) {
    double v0 = 0.5*r[0]*r[0];
    double m = mask(r);
    return v0*m + (1.0-m)*255.0;
}

double_complex dVsq(const coord_1d& r) {
    return double_complex(r[0]*r[0],0.0);
}

double energy(World& world, const complex_function_1d& v, const complex_function_1d& psi) {
    complex_derivative_1d D(world, 0);
    complex_function_1d du = D(psi);
    double_complex ke = 0.5*du.inner(du);
    double_complex pe = psi.inner(v*psi);

    print("  ke", real(ke), "pe", real(pe), "total", real(ke+pe), "norm", psi.norm2());

    return real(ke+pe);
}

complex_function_1d chin_chen(const complex_function_1d& expV,
                    const complex_function_1d& expVtilde,
                    const complex_operatorT* G,
                    const complex_function_1d& psi0) {

    // psi(t) = exp(-i*V*t/6) exp(-i*T*t/2) exp(-i*2*Vtilde*t/3) exp(-i*T*t/2) exp(-i*V*t/6)

    complex_function_1d psi1;

    psi1 = expV*psi0;
    psi1.reconstruct();
    psi1.broaden();
    psi1.broaden();
    psi1 = apply_1d_realspace_push(*G, psi1, 0);
    psi1.sum_down();
    psi1 = expVtilde*psi1;
    psi1 = apply_1d_realspace_push(*G, psi1, 0);
    psi1.sum_down();
    psi1 = expV*psi1;
    psi1.truncate();

    return psi1;
}

complex_function_1d trotter(const complex_function_1d& expV,
                  const complex_operatorT* G,
                  const complex_function_1d& psi0) {
    //    psi(t) = exp(-i*T*t/2) exp(-i*V*t) exp(-i*T*t/2) psi(0)

    psi0.reconstruct();
    psi0.broaden();
    psi0.broaden();
    complex_function_1d psi1 = apply_1d_realspace_push(*G, psi0, 0);
    psi1.sum_down();
    psi1 = expV*psi1;
    psi1 = apply_1d_realspace_push(*G, psi1, 0);
    psi1.sum_down();
    psi1.truncate();

    return psi1;
}

template<typename T, std::size_t NDIM>
struct unaryexp {
    void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
        UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};


/// Returns exp(-I*t*V) with truncation
complex_function_1d make_exp(double t, const complex_function_1d& v) {
    v.reconstruct();
    complex_function_1d expV = double_complex(0.0,-t)*v;
    expV.unaryop(unaryexp<double_complex,1>());
    expV.truncate();
    return expV;
}


void test_trotter(World& world) {

    const double L = 30.0;
    const int k = 16;
    const double thresh = 1e-12;
    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);

    std::cout.precision(8);

    complex_function_1d psi = complex_factory_1d(world).f(psi0);
    psi.scale(1.0/psi.norm2());
    complex_function_1d psi0 = copy(psi);

    complex_function_1d v = complex_factory_1d(world).f(V);

    energy(world,v,psi);

    double c = 10.0*sqrt(0.5) * 1.86;
    //double c = 10.0*sqrt(0.5) * 1.86 * 2.0;
    double tcrit = 2*constants::pi/(c*c);

    double tstep = tcrit * 0.125;

    print("The time step is", tstep, "\n");

    complex_operatorT* G0 = qm_1d_free_particle_propagator(k, c, tstep, 2*L);

    complex_function_1d expV = make_exp(tstep,v);

    for (int step=0; step<40; ++step) {
        double time = step * tstep;
        double_complex phase = psi0.inner(psi);
        double radius = abs(phase);
        double theta = arg(phase);
        double theta_exact = -time*0.5;
        while (theta_exact > constants::pi) theta_exact -= 2.0*constants::pi;

        print("step", step, "time", time, "radius", radius, "arg", theta, "exact", theta_exact, "phase err", theta_exact-theta);
        energy(world, v, psi);

        psi = trotter(expV, G0, psi);
    }

}

void test_chin_chen(World& world) {

    const double L = 20.0;
    const int k = 16;
    const double thresh = 1e-12;
    double c = 10.0*sqrt(0.5) * 1.86;
    double tcrit = 2*constants::pi/(c*c);
    double tstep = tcrit;
    print("The time step is", tstep, "\n");

    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);

    std::cout.precision(8);

    complex_function_1d psi = complex_factory_1d(world).f(psi0);
    psi.scale(1.0/psi.norm2());
    complex_function_1d psi0 = copy(psi);

    complex_function_1d v = complex_factory_1d(world).f(V);
    complex_function_1d dvsq = complex_factory_1d(world).f(dVsq);

    complex_function_1d vtilde = v - dvsq*(tstep*tstep/48.0);

    complex_operatorT* G0 = qm_1d_free_particle_propagator(k, c, tstep, 1400.0);
    //G.doleaves = true;
    complex_function_1d expV = make_exp(tstep/6.0, v);
    complex_function_1d expVtilde = make_exp(2.0*tstep/3.0, vtilde);

    for (int step=0; step<100; ++step) {
        double time = step * tstep;
        double_complex phase = psi0.inner(psi);
        double radius = abs(phase);
        double theta = arg(phase);
        double theta_exact = -time*0.5;
        while (theta_exact > constants::pi) theta_exact -= 2.0*constants::pi;
        while (theta_exact < -constants::pi) theta_exact += 2.0*constants::pi;

        print("step", step, "time", time, "radius", radius, "arg", theta, "exact", theta_exact, "phase err", theta_exact-theta);
        energy(world, v, psi);

        psi = chin_chen(expV, expVtilde, G0, psi);
    }

}

int main(int argc, char**argv) {
    const int required = MADNESS_MPI_THREAD_LEVEL;
    SafeMPI::Init_thread(argc, argv, required);
    World world(SafeMPI::COMM_WORLD);
    try {
        startup(world,argc,argv);

        bandlimited_propagator_plot();
        //test_trotter(world);
        //test_chin_chen(world);

    }
    catch (const SafeMPI::Exception& e) {
        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    }
//    catch (char* s) {
//        print(s);
//        error("caught a c-string exception");
//    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    SafeMPI::Finalize();

    return 0;
}

