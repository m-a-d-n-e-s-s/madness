#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include "dft.h"
#include "util.h"

namespace madness
{
  //*************************************************************************
  DFTNuclearPotentialOp::DFTNuclearPotentialOp(World& world, funcT V, 
      double coeff, double thresh) : EigSolverOp(world, coeff, thresh)
  {
    _V = V;
    ones = functorT(new OnesFunctor());
    zeros = functorT(new ZerosFunctor());
  }
  //*************************************************************************
  
  //*************************************************************************
  funcT DFTNuclearPotentialOp::op_r(const funcT& rho, const funcT& psi)
  {
    funcT rfunc = _V*psi;
    return rfunc;
  }
  //*************************************************************************

  //***************************************************************************
  DFT::DFT(World& world, funcT V, funcT phi, double eig, double thresh)
  : _world(world), _V(V), _thresh(thresh)
  {
    // Create temporary list for eigsolver
    std::vector<funcT> phis;
    std::vector<double> eigs;
    phis.push_back(phi);
    eigs.push_back(eig);
    // Create ops list 
    std::vector<EigSolverOp*> ops;
    // Add nuclear potential to ops list
    ops.push_back(new DFTNuclearPotentialOp(world, V, 1.0, thresh));

    // Create solver
    _solver = new EigSolver(world, phis, eigs, ops, thresh);
    _solver->addObserver(this);

    // Misc.
    ones = functorT(new OnesFunctor());
    zeros = functorT(new ZerosFunctor());
  }
  //***************************************************************************
    
  //*****************************************************************************
  DFT::~DFT()
  {
    delete _solver;
  }
  //*****************************************************************************
}
