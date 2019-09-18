#ifndef MADNESS_APPS_DF_H_INCLUDED
#define MADNESS_APPS_DF_H_INCLUDED

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/nonlinsol.h>  // The kain solver
#include <vector>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <complex>
#include <cmath>
#include <string>
#include <algorithm> 
#include "../chem/molecule.h"
#include "DFParameters.h"
#include "InitParameters.h"


using namespace madness;


/// Given a molecule and nonrelativistic ground state orbitals, solve the Dirac-Hartree-Fock equations
class DF {
     private:
          // Member variables

          // DFParameter object to hold all user input variables
          DFParameters DFparams;

          // NRParameter object to hold all variables needed from
          // nonrelativistic ground state calculation. Read from an archive
          InitParameters Init_params;

          // Timer variables
          std::vector<double> sss, ttt;

          // Tensor for holding energies
          Tensor<double> energies;        

          //Tensor for holding energies of virtual orbitals
          Tensor<double> v_energies;

          //Vector of DF Fcwf occupied orbitals
          std::vector<Fcwf> occupieds;

          //Vector of DF Fcwf virtual orbitals
          std::vector<Fcwf> virtuals;

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

          // Collective constructor for response uses contents of file \c filename and broadcasts to all nodes
          DF(World & world,            // MADNESS world object
              const char* input_file);  // Input file 

          // Collective constructor for DF uses contents of stream \c input and broadcasts to all nodes
          DF(World & world,                        // MADNESS world object
              std::shared_ptr<std::istream> input); // Pointer to input stream

          //Calculates the kinetic+rest energy expectation value of psi
          double rele(World& world, Fcwf& psi);

          //Applies the exchange operator to all of psis
          void exchange(World& world, real_convolution_3d& op, std::vector<Fcwf>& Kpsis);

          //Applies the exchange operator defined by psis to a single generic input fcwf
          Fcwf apply_K(World& world, real_convolution_3d& op, Fcwf& phi);
              
          //diagonalizes psis in the Fock space. Transforms psis and Kpsis.
          void diagonalize(World& world, real_function_3d& myV,real_convolution_3d& op, std::vector<Fcwf>& Kpsis);

          //diagonalizes virtual orbitals. Transforms virtuals and Kpsis
          Tensor<double>  diagonalize_virtuals(World& world, real_function_3d& JandV,real_convolution_3d& op, std::vector<Fcwf>& Kpsis);
              
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

          //Load balancing function
          void DF_load_balance(World& world, real_function_3d& Vnuc);

          //Does one full SCF iteration
          bool iterate(World& world, real_function_3d& V, real_convolution_3d& op, real_function_3d& JandV, std::vector<Fcwf>& Kpsis, XNonlinearSolver<std::vector<Fcwf>, std::complex<double>, Fcwf_vector_allocator>& kainsolver, double& tolerance, int& iteration_number, double& nuclear_repulsion_energy);
                

          //Runs the job specified in the input parameters
          void solve(World& world);

          //solves the Dirac Fock equation for the occupied orbitals   
          void solve_occupied(World & world);

          //solves for n-1 occupieds, then the last occupied, then some number of virtuals
          void solve_virtuals1(World& world);

          //Lineplot the densities. Currently only along x axis from 0 to L
          void make_density_lineplots(World& world, const char* filename, int npt, double endpnt);

          //Lineplot the densities of the large and small component separately. only along x axis from 0 to L
          void make_component_lineplots(World& world, const char* filename1, const char* filename2, int npt, double endpnt);

          void orthogonalize_inplace(World& world);

          //Lineplot the densities of the large and small component separately. only along x axis on log scale from 10^-startpnt to 10^endpnt with pts evenly spaced in log space
          void make_component_logplots(World& world, const char* filename1, const char* filename2, int npt, int startpnt, int endpnt);

          void print_sizes(World& world);

};


#endif

//kthxbye

