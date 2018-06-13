#ifndef MADNESS_APPS_MOLDFT_FCWF_H_INCLUDED
#define MADNESS_APPS_MOLDFT_FCWF_H_INCLUDED

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <vector>
#include <math.h>
#include <complex>

using namespace madness;

class Fcwf{
     std::vector<complex_function_3d> m_psi;
     bool m_initialized;

public:
     
     Fcwf();

     Fcwf(const complex_function_3d& wf1,
          const complex_function_3d& wf2,
          const complex_function_3d& wf3,
          const complex_function_3d& wf4);

     Fcwf(World& world);

     complex_function_3d& operator[](const int i);

     const complex_function_3d& operator[](const int i) const ;
     
     Fcwf(std::vector<complex_function_3d>& phi);

     bool getinitialize();

     bool getinitialize() const ;

     unsigned int size();

     unsigned int size() const ;

     Fcwf(const Fcwf& phi);

     Fcwf operator=(const Fcwf& phi);

     Fcwf operator-(const Fcwf& phi) const ;

     Fcwf operator+(const Fcwf& phi);

     Fcwf operator*(std::complex<double> a);
     
     void scale(std::complex<double> a);
     
     Fcwf operator+=(const Fcwf& phi);

     Fcwf operator-=(const Fcwf& phi);

     double norm2();

     void normalize();

     Fcwf operator*(madness::complex_function_3d phi);

     Fcwf operator*(madness::real_function_3d phi);

     void truncate();

};

std::complex<double> inner(const Fcwf& psi, const Fcwf& phi);

Fcwf apply(real_convolution_3d& op, const Fcwf& psi);

real_function_3d squaremod(Fcwf psi);
real_function_3d squaremod_small(Fcwf psi);
real_function_3d squaremod_large(Fcwf psi);

complex_function_3d inner_func(Fcwf psi, Fcwf phi);

Fcwf copy(Fcwf psi);

std::complex<double> inner(std::vector<Fcwf> a, std::vector<Fcwf> b);

std::vector<Fcwf> operator*(std::vector<Fcwf> psis, std::complex<double> a);

std::vector<Fcwf> operator*(std::complex<double> a, std::vector<Fcwf> psis);

void operator+=(std::vector<Fcwf>& phi, std::vector<Fcwf> psi);

std::vector<Fcwf> operator-(std::vector<Fcwf> phi, std::vector<Fcwf> psi);

class Fcwf_vector_allocator {
     World& world;
     unsigned int m_size;
     public:
          //Constructor
          Fcwf_vector_allocator(World& world, unsigned int m_size);

          //Overloading () operator
          std::vector<Fcwf> operator()();

          //Copy Constructor
          //according to Bryan, according to Jakob, this is necessary for KAIN?
          Fcwf_vector_allocator operator=(const Fcwf_vector_allocator& other);

          void set_size(int size);
};

#endif

//kthxbye
