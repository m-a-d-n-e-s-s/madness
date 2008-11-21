#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include "dft.h"
#include "util.h"
#include <moldft/xc/f2c.h>
#include <vector>
#include "poperator.h"
#include "lda.h"

typedef madness::Vector<double,3> coordT;

namespace madness
{
  //***************************************************************************
  template <typename T, int NDIM>
  DFTNuclearChargeDensityOp<T,NDIM>::DFTNuclearChargeDensityOp(World& world, funcT rhon,
      double coeff, double thresh, bool periodic) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("NuclearChargeDensityOp");
    _rhon = rhon;
    SeparatedConvolution<T,NDIM>* cop;
    if (periodic)
    {
      Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
      cop = PeriodicCoulombOpPtr<T,NDIM>(world, FunctionDefaults<NDIM>::get_k(),
          1e-4, thresh, L);
    }
    else
    {
      cop =
        CoulombOperatorPtr<T,NDIM>(world,
            FunctionDefaults<NDIM>::get_k(), 1e-8, thresh);
    }
    // Apply operator to get potential
    _Vnuc = apply(*cop, rhon);
    delete cop;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  DFTNuclearPotentialOp<T,NDIM>::DFTNuclearPotentialOp(World& world, funcT V,
      double coeff, double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("NuclearPotentialOp");
    _V = V;
  }
  //***************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  DFTCoulombOp<T,NDIM>::DFTCoulombOp(World& world, double coeff,
      double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("CoulombOp");
    // For now, no spin polarized
    _spinpol = false;
    // Create Coulomb operator
    _cop = CoulombOperatorPtr<T,NDIM>(world, FunctionDefaults<NDIM>::get_k(),
        1e-4, thresh);
    // Initialize potential
    _Vc = FunctionFactory<T,NDIM>(world);
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  DFTCoulombPeriodicOp<T,NDIM>::DFTCoulombPeriodicOp(World& world, funcT rhon, double coeff,
      double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Nuclear charge density
    _rhon = rhon;
    // Message for the matrix element output
    this->messageME("PeriodicCoulombOp");
    // For now, no spin polarized
    _spinpol = false;
    // Create Coulomb operator
    Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
    _cop = PeriodicCoulombOpPtr<T,NDIM>(world, FunctionDefaults<NDIM>::get_k(),
        1e-4, thresh, L);
  }
  //*************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFTCoulombOp<T,NDIM>::prepare_op(Function<double,NDIM> rho)
  {
    _Vc = apply(*_cop, rho);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFTCoulombPeriodicOp<T,NDIM>::prepare_op(Function<double,NDIM> rho)
  {
    _Vc = apply(*_cop, rho + _rhon);
////    _Vc.reconstruct();
////    rho.reconstruct();
////    _rhon.reconstruct();
//    Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
//    printf("\n");
//    double bstep = L[0] / 100;
//    for (int i=0; i<101; i++)
//    {
//      coordT p(-L[0]/2 + i*bstep);
//      printf("%.2f\t\t%.8f\t%.8f\t%.8f\n", p[0], _rhon(p), rho(p), _Vc(p));
//    }
//    printf("\n");
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTNuclearChargeDensityOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    funcT rfunc = _Vnuc * psi;
    return rfunc;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTNuclearPotentialOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    funcT rfunc = _V * psi;
    return rfunc;
  }
  //***************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTCoulombOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
//    if (_world.rank() == 0) printf("Applying Coulomb operator ...\n\n");
    //printf("Applying Coulomb operator ...\n\n");
    funcT rfunc = _Vc * psi;
    return  rfunc;
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTCoulombPeriodicOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
//    if (_world.rank() == 0) printf("Applying Coulomb operator ...\n\n");
    //printf("Applying Coulomb operator ...\n\n");
    funcT rfunc = _Vc * psi;
    return  rfunc;
  }
  //*************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  XCFunctionalLDA<T,NDIM>::XCFunctionalLDA(World& world, double coeff, double thresh)
    : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("XCFunctionalLDA");
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> XCFunctionalLDA<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    funcT V_rho = 0.5 * copy(rho);
    V_rho.reconstruct();
    V_rho.unaryop(&::dft_xc_lda_V<NDIM>);
//    funcT V_rho = copy(rho);
//    V_rho.scale(0.5);
//    V_rho.unaryop(&::ldaop);
    funcT rfunc = V_rho * psi;

      if (this->_world.rank() == 0)  printf("\n");
      double L = 30.0;
      double bstep = L / 100.0;
      rho.reconstruct();
      V_rho.reconstruct();
      for (int i = 0; i < 101; i++)
      {
        coordT p(-L / 2 + i * bstep);
        if (this->_world.rank() == 0)
          printf("%.2f\t\t%.8f\t%.8f\n", p[0], rho(p), V_rho(p));
      }
      if (this->_world.rank() == 0) printf("\n");

    return rfunc;
  }
  //***************************************************************************
}

namespace madness
{
  //***************************************************************************
  template <typename T, int NDIM>
  DFT<T,NDIM>::DFT(World& world, funcT rhon, std::vector<funcT> phis,
      std::vector<double> eigs, ElectronicStructureParams params)
  : _world(world), _rhon(rhon), _params(params)
  {

    if (world.rank() == 0 && !params.periodic) printf("DFT constructor (non-peridic) ...\n\n");
    if (world.rank() == 0 && params.periodic) printf("DFT constructor (periodic) ...\n\n");

    // Create ops list
    std::vector<EigSolverOp<T,NDIM>*> ops;
    // Add nuclear potential to ops list
//    ops.push_back(new DFTNuclearPotentialOp<T,NDIM>(world, V, 1.0, thresh));
    if (params.periodic)
    {
      ops.push_back(new DFTCoulombPeriodicOp<T,NDIM>(world, rhon, 1.0, params.thresh));
    }
    else
    {
      ops.push_back(new DFTNuclearChargeDensityOp<T,NDIM>(world, rhon, 1.0, params.thresh, false));
      ops.push_back(new DFTCoulombOp<T,NDIM>(world, 1.0, params.thresh));
    }
    _xcfunc = new XCFunctionalLDA<T,NDIM>(world, 1.0, params.thresh);
    ops.push_back(_xcfunc);

    // Create solver
    if (params.periodic)
    {
      std::vector<kvecT> kpoints;
      kvecT gammap(0.0);
      kpoints.push_back(gammap);
      _solver = new EigSolver<T,NDIM>(world, rhon, phis, eigs, ops, kpoints, params);
    }
    else
    {
      _solver = new EigSolver<T,NDIM>(world, rhon, phis, eigs, ops, params);
    }
    _solver->addObserver(this);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  DFT<T,NDIM>::~DFT()
  {
    delete _solver;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFT<T,NDIM>::solve(int maxits)
  {
    _solver->solve(maxits);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_ke_sp(funcT psi, bool periodic)
  {
    // Do calculation
    double kenergy = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
      funcT dpsi = diff(psi, axis);
      kenergy += 0.5 * inner(dpsi, dpsi);
    }
    return kenergy;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_ke_sp(const std::vector<funcT>& phis, bool spinpol,
      bool periodic)
  {
    double tot_ke = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate kinetic energy contribution from psi
      tot_ke += calculate_ke_sp(psi);
    }
    if (!spinpol) tot_ke *= 2.0;
    return tot_ke;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_pe_sp(const World& world, const Function<double,NDIM>& rho,
      const Function<double,NDIM>& rhon, bool spinpol, const double thresh, bool periodic)
  {
    // Create Coulomb operator
    SeparatedConvolution<T,NDIM>* op = 0;
    if (periodic)
    {
      Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
      op = CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
          FunctionDefaults<NDIM>::get_k(), 1e-4, thresh);
    }
    else
    {
      op = CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
          FunctionDefaults<NDIM>::get_k(), 1e-4, thresh);
    }
    // Apply Coulomb operator and trace with the density
    funcT tmp = rhon + rho;
    funcT Vnuc = apply(*op, tmp);

//    // DEBUG ******************************************************************
//    if (world.rank() == 0) cout << "Printing out the electronic charge density and potential ..." << endl;
//    if (world.rank() == 0) printf("\n");
//    //Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
//    double LLL = L[0];
//    double bstep = LLL / 100.0;
//    Vnuc.reconstruct();
//    Vnuc2.reconstruct();
//    for (int i=0; i<101; i++)
//    {
//      coordT p(-LLL/2 + i*bstep);
//      if (world.rank() == 0) printf("%.2f\t\t%.8f\t%.8f\n", p[0], Vnuc(p), Vnuc2(p));
//    }
//    if (world.rank() == 0) printf("\n");
//    // DEBUG ******************************************************************

    double tot_pe = inner(Vnuc, rho);
    delete op;
    return tot_pe;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_coulomb_energy(const World& world, const Function<double, NDIM>& rho,
      bool spinpol, const double thresh, bool periodic)
  {
    // Create Coulomb operator
    SeparatedConvolution<T,NDIM>* op;
    if (periodic)
    {
      Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
      op = PeriodicCoulombOpPtr<T,NDIM>(const_cast<World&>(world),
        FunctionDefaults<NDIM>::get_k(), 1e-4, thresh, L);
    }
    else
    {
      op = CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
        FunctionDefaults<NDIM>::get_k(), 1e-4, thresh);
    }
    // Apply Coulomb operator and trace with the density
    funcT Vc = apply(*op, rho);

    double tot_ce = 2*Vc.inner(rho);
    delete op;
    return tot_ce;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_xc_energy(const Function<double, NDIM>& rho)
  {
    funcT enefunc = 0.5 * copy(rho);
    enefunc.reconstruct();
    enefunc.unaryop(&dft_xc_lda_ene<NDIM>);
    return enefunc.trace();
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFT<T,NDIM>::iterateOutput(const std::vector<funcT>& phis,
      const std::vector<double>& eigs, const Function<double, NDIM>& rho,
      const int& iter, bool periodic)
  {
    if (iter%3 == 0)
    {
      if (world().rank() == 0) printf("Calculating energies ...\n");
      if (world().rank() == 0) printf("Calculating KE ...\n");
      double ke = DFT::calculate_tot_ke_sp(phis, false, periodic);
      if (world().rank() == 0) printf("Calculating PE and CE...\n");
      double pece = DFT::calculate_tot_pe_sp(_world, rho, _rhon, false, _params.thresh, periodic);
//      if (world().rank() == 0) printf("Calculating CE ...\n");
//      double ce = DFT::calculate_tot_coulomb_energy(_world, rho, false, _thresh, periodic);
      if (world().rank() == 0) printf("Calculating EE ...\n");
      double xce = DFT::calculate_tot_xc_energy(rho);
      if (world().rank() == 0) printf("Calculating NE ...\n");
      double ne = 0.0;
      if (world().rank() == 0) printf("Kinetic energy:\t\t\t\t %.8f\n", ke);
      if (world().rank() == 0) printf("Potential and Coulomb energy:\t\t %.8f\n", pece);
//      if (world().rank() == 0) printf("Coulomb energy:\t\t\t %.8f\n", ce);
      if (world().rank() == 0) printf("XC energy:\t\t\t\t %.8f\n", xce);
      if (world().rank() == 0) printf("Total energy:\t\t\t\t %.8f\n", ke + pece + xce + ne);
      if (world().rank() == 0) printf("gs ene\t\t\t\t\t%.4f\n", eigs[0]);
      if (world().rank() == 0) printf("1st es ene\t\t\t\t%.4f\n", eigs[1]);
      T mtxe = matrix_element(phis[0], phis[0]);
      if (world().rank() == 0) printf("\nKS matrix element:\t\t\t%.8f\n\n", mtxe);
      print_matrix_elements(phis[0], phis[0]);
    }
  }
  //***************************************************************************

  //***************************************************************************
//  template class DFT<double, 1>;
//  template class DFT<double, 2>;
  template class DFT<double, 3>;


//  template class DFT< std::complex<double>, 3>;
  //***************************************************************************
}
