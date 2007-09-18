#include "hartreefock.h"

namespace madness
{

  //***************************************************************************
  HartreeFock::HartreeFock(World& world, funcT V, std::vector<funcT> phis)
    : _world(world), _V(V), _phis(phis)
  {
    ones = functorT(new OnesFunctor());
    zeros = functorT(new ZerosFunctor());
  }
  //***************************************************************************
  
  //***************************************************************************
  HartreeFock::~HartreeFock()
  {
  }
  //***************************************************************************
  
  //***************************************************************************
  void HartreeFock::hartree_fock(int maxits)
  {
    for (int it = 0; it < maxits; it++)
    {
      for (std::vector<funcT>::iterator pi = _phis.begin();
            pi != _phis.end(); ++pi)
      {
        // Get psi from collection
        funcT psi = (*pi);
        // Calculate nuclear contribution
        funcT pnuclear = _V*psi;
        // Calculate the Coulomb contribution to the Fock operator (J)
        funcT pcoulomb = calculate_coulomb(psi);
        // Calculate the Exchange contribution to the Fock operator (K)
        funcT pexchange = calculate_exchange(psi);
        // Get new wavefunction
        //funcT pfunc = pnuclear + 2.0 * pcoulomb - pexchange;
        funcT pfunc = FunctionFactory<double,3>(_world).functor(zeros);
        // Create the free-particle Green's function operator
        
        // Apply the Green's function operator (stubbed)
        funcT tmp = FunctionFactory<double,3>(_world).functor(zeros);
        // (Not sure whether we have to do this mask thing or not!)
        for (std::vector<funcT>::iterator pj = _phis.begin();
              pj != pi; ++pj)
        {
          // Project out the lower states
          if (pi != pj)
          {
            // Get other wavefunction
            funcT psij = (*pj);
            double overlap = 0.0;
            //tmp -= overlap*psij;
          }
        }
      }
    }
  }
  //***************************************************************************

  //***************************************************************************
  funcT HartreeFock::calculate_coulomb(funcT psi)
  {
    // Return value
    funcT rfunc = FunctionFactory<double,3>(_world).functor(zeros);
    if (include_coulomb())
    {
      // Create Coulomb operator
      
      for (std::vector<funcT>::iterator pj = _phis.begin(); pj != _phis.end(); ++pj)
      {
        // Get phi(j) from iterator
        funcT& phij = (*pj);
        // Compute the density
        funcT prod = phij*phij;
        // Transform Coulomb operator into a function (stubbed)
        funcT Vc = FunctionFactory<double,3>(_world).functor(zeros);
        // Note that we are not using psi
        // The density is built from all of the wavefunctions. The contribution
        // psi will be subtracted out later during the exchange.
        //rfunc += Vc*psi;
      }
    }
    return rfunc;
  }
  //***************************************************************************

  //***************************************************************************
  funcT HartreeFock::calculate_exchange(funcT psi)
  {
    // Return value
    funcT rfunc = FunctionFactory<double,3>(_world).functor(zeros);
    if (include_exchange())
    {
      // Create Coulomb operator
      
      // Use the psi and pj wavefunctions to build a product so that the K 
      // operator can be applied to the wavefunction indexed by pj, NOT PSI.
      for (std::vector<funcT>::iterator pi = _phis.begin(); pi != _phis.end(); ++pi)
      {
        for (std::vector<funcT>::iterator pj = _phis.begin(); pj != _phis.end(); ++pj)
        {
          // Get phi(j) from iterator
          funcT& phij = (*pj);
          // NOTE that psi is involved in this calculation
          funcT prod = phij*psi;
          // Transform Coulomb operator into a function (stubbed)
          funcT Vex = FunctionFactory<double,3>(_world).functor(zeros);
          // NOTE that the index is j.
          //rfunc += Vex*phij;
        }
      }
    }
    return rfunc;
  }
}
