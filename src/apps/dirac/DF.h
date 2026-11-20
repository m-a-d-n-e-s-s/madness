#ifndef MADNESS_APPS_DF_H_INCLUDED
#define MADNESS_APPS_DF_H_INCLUDED

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/nonlinsol.h>  // The kain solver
#include <tuple>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <complex>
#include <cmath>
#include <string>
#include <algorithm> 
#include <madness/chem/molecule.h>
#include "DFConvergence.h"
#include "DFKramers.h"
#include "DFParameters.h"
#include "InitParameters.h"


using namespace madness;

/// Real-Hilbert-space, gauge-consistent KAIN nonlinear subspace solver for 4-component spinors.
class DFKainSolver {
private:
     Fcwf_vector_allocator alloc_;
     unsigned int maxsub_ = 10;
     std::vector<std::vector<Fcwf>> ulist_;
     std::vector<std::vector<Fcwf>> rlist_;
     Tensor<double> Q_;

     static std::vector<Fcwf> deep_copy_vec(const std::vector<Fcwf>& src) {
          std::vector<Fcwf> dst;
          dst.reserve(src.size());
          for (const auto& f : src) {
               dst.push_back(copy(f));
          }
          return dst;
     }

     static double real_inner_vec(std::vector<Fcwf>& a, std::vector<Fcwf>& b) {
          return std::real(inner(a, b));
     }

public:
     explicit DFKainSolver(const Fcwf_vector_allocator& alloc, unsigned int maxsub = 10)
          : alloc_(alloc), maxsub_(maxsub) {}

     void set_maxsub(int maxsub) {
          maxsub_ = (maxsub > 0) ? static_cast<unsigned int>(maxsub) : 1u;
     }

     bool empty() const noexcept { return ulist_.empty(); }
     std::size_t size() const noexcept { return ulist_.size(); }

     std::vector<std::vector<Fcwf>>& get_ulist() noexcept { return ulist_; }
     std::vector<std::vector<Fcwf>>& get_rlist() noexcept { return rlist_; }

     void clear_subspace() {
          ulist_.clear();
          rlist_.clear();
          Q_ = Tensor<double>();
     }

     void rebuild_Q() {
          const long m = static_cast<long>(ulist_.size());
          if (m == 0) {
               Q_ = Tensor<double>();
               return;
          }
          Tensor<double> Qnew(m, m);
          for (long i = 0; i < m; ++i) {
               for (long j = 0; j < m; ++j) {
                    Qnew(i, j) = real_inner_vec(ulist_[i], rlist_[j]);
               }
          }
          Q_ = Qnew;
     }

     std::vector<Fcwf> update(const std::vector<Fcwf>& u,
                              const std::vector<Fcwf>& r,
                              double cabsmax = 3.0,
                              double rcond_max = 0.01) {
          if (maxsub_ <= 1) {
               return u - r;
          }

          ulist_.push_back(deep_copy_vec(u));
          rlist_.push_back(deep_copy_vec(r));
          const long m = static_cast<long>(ulist_.size());

          Tensor<double> Qnew(m, m);
          if (m > 1) {
               Qnew(Slice(0, -2), Slice(0, -2)) = Q_;
          }
          for (long i = 0; i < m; ++i) {
               Qnew(i, m - 1) = real_inner_vec(ulist_[i], rlist_[m - 1]);
               if (i < m - 1) {
                    Qnew(m - 1, i) = real_inner_vec(ulist_[m - 1], rlist_[i]);
               }
          }
          Q_ = Qnew;

          auto [c, reset] = solve_real_kain_weights(Q_, cabsmax, rcond_max);
          if (reset) {
               std::vector<Fcwf> u_last = std::move(ulist_.back());
               std::vector<Fcwf> r_last = std::move(rlist_.back());
               const double q_last = Q_(m - 1, m - 1);
               ulist_.clear();
               rlist_.clear();
               ulist_.push_back(std::move(u_last));
               rlist_.push_back(std::move(r_last));
               Q_ = Tensor<double>(1L, 1L);
               Q_(0L, 0L) = q_last;
          }

          std::vector<Fcwf> unew = alloc_();
          const std::size_t n_cur = ulist_.size();
          for (std::size_t s = 0; s < n_cur; ++s) {
               const std::complex<double> cs(c(static_cast<long>(s)), 0.0);
               unew += (ulist_[s] - rlist_[s]) * cs;
          }

          if (ulist_.size() >= maxsub_) {
               ulist_.erase(ulist_.begin());
               rlist_.erase(rlist_.begin());
               Q_ = copy(Q_(Slice(1, -1), Slice(1, -1)));
          }
          return unew;
     }
};


/// Given a molecule and nonrelativistic ground state orbitals, solve the Dirac-Hartree-Fock equations
class DF {
     private:
          // Member variables

          // DFParameter object to hold all user input variables
          DFParameters DFparams;

          // InitParametesr object to hold all variables needed from
          // nonrelativistic ground state calculation. Read from an archive
          InitParameters Init_params;

          // Timer variables
          std::vector<double> sss, ttt;

          // Tensor for holding energies
          Tensor<double> energies;        

          //Vector of DF Fcwf occupied orbitals
          std::vector<Fcwf> occupieds;

          //Total energy of the system
          double total_energy;

          //Whether or not the calculation to be done is closed shell
          bool closed_shell;

     public:

          // Start a timer
          void start_timer(World & world);

          // Needed to do timers correctly
          double pop(std::vector<double> & v);

          // Stop a timer
          Tensor<double> end_timer(World & world);

          //Find current time (relative to job start)
          Tensor<double> get_times(World& world);

          // Collective constructor uses contents of file \c filename and broadcasts to all nodes
          DF(World & world,            // MADNESS world object
              const char* input_file);  // Input file 

          // Collective constructor uses contents of stream \c input and broadcasts to all nodes
          DF(World & world,                        // MADNESS world object
              std::shared_ptr<std::istream> input); // Pointer to input stream

          static void help() {
              print_header2("help page for DIRAC ");
              print("The DIRAC code computes Dirac-Fock energies");
          }

          static void print_parameters() {
              print("no parameter help is available for the dirac code");
              print("\nYou need a dft-block and a DiracFock-block in an input file named 'input'");
          }



        //Calculates the kinetic+rest energy expectation value of psi
          double rele(World& world, Fcwf& psi);

          //Applies the exchange operator to all of occupieds
          void exchange(World& world, real_convolution_3d& op, std::vector<Fcwf>& Kpsis);

          //diagonalizes occupieds in the Fock space. Transforms occupieds, Kpsis, and optional KAIN subspace.
          void diagonalize(World& world, real_function_3d& myV, real_convolution_3d& op,
                           std::vector<Fcwf>& Kpsis, DFKainSolver* kainsolver = nullptr);

          // Small function to print geometry of a molecule nicely
          // Straight up stolen from Bryan
          void print_molecule(World &world);

          //Saves the state of a DF job so that it can be restarted
          void saveDF(World& world);

          //Creates the gaussian nuclear potential from the molecule object
          void make_gaussian_potential(World& world, real_function_3d& potential);

          //Creates the gaussian nuclear potential from the molecule object. Also calculates the nuclear repulsion energy
          void make_gaussian_potential(World& world, real_function_3d& potential, double& nuclear_repulsion_energy);

          //Creates the fermi nuclear potential from the molecule object
          void make_fermi_potential(World& world, real_convolution_3d& op, real_function_3d& potential);

          //Creates the fermi nuclear potential from the molecule object. Also calculates the nuclear repulsion energy
          void make_fermi_potential(World& world, real_convolution_3d& op, real_function_3d& potential, double& nuclear_repulsion_energy);

          //Creates the point nuclear potential from the molecule object
          void make_point_potential(World& world, real_function_3d& potential);

          //Creates the point nuclear potential from the molecule object. Also calculates the nuclear repulsion energy
          void make_point_potential(World& world, real_function_3d& potential, double& nuclear_repulsion_energy);

          //Load balancing function
          void DF_load_balance(World& world, real_function_3d& Vnuc);

          //Does one full SCF iteration
          std::tuple<bool, double, real_function_3d>
          iterate(World &world, real_function_3d &V, real_convolution_3d &op,
                  real_function_3d &JandV, std::vector<Fcwf> &Kpsis,
                  DFKainSolver &kainsolver,
                  double &tolerance, int &iteration_number,
                  double &nuclear_repulsion_energy, real_function_3d &prev_rho);

          //Runs the job specified in the input parameters
          void solve(World& world);

          //solves the Dirac Fock equation for the occupied orbitals   
          void solve_occupied(World & world);

          //Lineplot the densities. Currently only along x axis from 0 to L
          void make_density_lineplots(World& world, const char* filename, int npt, double endpnt);

          //Lineplot the densities of the large and small component separately. only along x axis from 0 to L
          void make_component_lineplots(World& world, const char* filename1, const char* filename2, int npt, double endpnt);

          //orthogonormalize occupieds, overwriting the ones in memory
          void orthogonalize_inplace(World& world);

          //Lineplot the densities of the large and small component separately. only along x axis on log scale from 10^-startpnt to 10^endpnt with pts evenly spaced in log space
          void make_component_logplots(World& world, const char* filename1, const char* filename2, int npt, int startpnt, int endpnt);

          //print the number of coefficients used to represent all occupied orbitals. primarily used when debugging.
          void print_sizes(World& world, bool individual);

};


#endif

//kthxbye
