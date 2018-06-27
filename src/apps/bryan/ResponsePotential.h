/*
 *
 *   Small class to hold all the functionality needed in creating/manipulating
 *   the  term in the response equations. (Gamma is the perturbed two 
 *   electron piece.)
 *
 */

#ifndef MADNESS_APPS_TDHF_RESPONSEPOTENTIAL_INCLUDE
#define MADNESS_APPS_TDHF_RESPONSEPOTENTIAL_INCLUDE

#include "ResponseFunction2.h"
#include <vector>
#include <memory>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include "../chem/projector.h"

namespace madness {

class ResponsePotential {
   // Member Variables
   private:
      World& world;                                                    // MPI communicator
      std::vector<real_function_3d> ground_orbitals;                   // Pointer to ground state orbitals
      QProjector<double,3> projector;                                  // Will project out ground state
      const double small;                                              // Smallest lengthscale for coulomb op.
      const double thresh;                                             // Truncation threshold for coulomb op.
      real_convolution_3d op;                                          // Coulomb operator
      ResponseFunction potential;                                      // For storing the ground state
      const bool store_potential;                                      // Flag for storing the potential (default = true)
      const bool TDA;                                                  // Flag indicating whether the TDA is being used (default = true)

   // Member Functions
   public:
   
      // ResponsePotential constructor
      ResponsePotential(World& world, std::vector<real_function_3d>& ground, double small, double thresh, 
            int r_states, int g_states, bool store = true, bool TDA = true)
            : world(world)
            , ground_orbitals(ground)
            , projector(QProjector<double,3>(world, ground_orbitals))
            , small(small)
            , thresh(thresh)
            , op(CoulombOperator(world, small, thresh))
            , potential(ResponseFunction(world, r_states, g_states)) 
            , store_potential(store)
            , TDA(TDA) {}

      /// Calculates coulomb like terms in the potential
      /// @param[in]    f        current perturbed orbitals
      /// @param[inout] gammaK   perturbed potential applied to ground state orbitals
      /// @param[inout] groundJ  ground state potential applied to perturbed orbitals
      void coulomb_terms(ResponseFunction& f,
                         ResponseFunction& gammaK,
                         ResponseFunction& groundJ) {
          // Calculate intermediaries
          // If store_potential is true, only calculate it once
          if(!store_potential or potential.r_states == 0) {
             // Clear any old data
             potential.clear();

             // Calculate potential
             for(unsigned int i = 0; i < f.r_states; i++) {
                std::vector<real_function_3d> psif = mul_sparse(world, ground_orbitals[i], ground_orbitals, thresh);
                truncate(world, psif);
                psif = apply(world, op, psif);
                truncate(world, psif);
                potential.push_back(psif);
             }            
          }

          // Calc. gammaK
          gammaK.clear();
          for(unsigned int k = 0; k < f.r_states; k++) {
             for(unsigned int p = 0; p < f.g_states; p++) {
                std::vector<real_function_3d> temp = std::vector<real_function_3d>(f.g_states);
                for(unsigned int i = 0; i < f.g_states; i++) {
                   temp[i] += potential[i][p] * f[k][i];
                }
                gammaK.push_back(temp);
             }
          } 

          // Calc. groundJ
          real_function_3d coulomb = potential[0][0];
          for(unsigned int i = 1; i < f.g_states; i++) {
             coulomb += potential[i][i];
          }
          groundJ = f * coulomb;

          // Project out ground states from gammaK
          for(unsigned int i = 0; i < gammaK.r_states; i++) {
             gammaK[i] = projector(gammaK[i]);
          }
          // Project out ground states from groundJ 
          for(unsigned int i = 0; i < groundJ.r_states; i++) {
             groundJ[i] = projector(groundJ[i]);
          }
      } 


      /// Calculates exchange like terms in the potential
      /// @param[in]    f        current perturbed orbitals
      /// @param[inout] gammaJ   perturbed potential applied to ground state orbitals
      /// @param[inout] groundK  ground state potential applied to perturbed orbitals
      void exchange_terms(ResponseFunction& f,
                          ResponseFunction& gammaJ,
                          ResponseFunction& groundK) {
          // Clear inputs
          gammaJ.zero();
          groundK.zero();

          // Going to calculate each transition density
          // and then use it to calculare correct pieces of
          // gammaJ and groundK before moving on          
          for(unsigned int k = 0; k < f.r_states; k++) {
             for(unsigned int p = 0; p < f.g_states; p++) {
                for(unsigned int i = 0; i < f.g_states; i++) {
                   // Get the transition density
                   real_function_3d rho = f[k][p] * ground_orbitals[i];

                   // Apply coulomb op
                   rho = apply(op, rho);

                   // Add pieces to gammaJ
                   gammaJ[k][p] += rho * ground_orbitals[p];

                   // Add pieces to groundK
                   groundK[k][p] += rho * ground_orbitals[i];
                }
             }
          } 

          // Project out ground states from gammaK
          for(unsigned int i = 0; i < gammaJ.r_states; i++) {
             gammaJ[i] = projector(gammaJ[i]);
          }
          // Project out ground states from groundJ 
          for(unsigned int i = 0; i < groundK.r_states; i++) {
             groundK[i] = projector(groundK[i]);
          }
      
     }
};

} // End namespace madness

#endif

// Dueces
