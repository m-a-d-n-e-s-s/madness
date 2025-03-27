
#include "ResponsePreliminariesCalculator.hpp"
#include "functypedefs.h"
#include "molresponse_v2/GroundStateData.hpp"
#include "molresponse_v2/ResponsePreliminaries.hpp"
#include "xcfunctional.h"
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/tensor/tensor.h>
typedef std::shared_ptr<operatorT> poperatorT;
using namespace madness;

ResponsePreliminariesCalculator::ResponsePreliminariesCalculator(
    World &world, const GroundStateData &ground_state, double coulomb_lo_thresh,
    double potential_thresh)
    : world_(world), ground_state_(ground_state),
      coulomb_lo_thresh_(coulomb_lo_thresh),
      potential_thresh_(potential_thresh) {

  xcf_.initialize(ground_state_.xc, ground_state_.spinrestricted, world, true);
  potential_manager_ =
      std::make_shared<PotentialManager>(ground_state_.molecule, "a");
  coulomb_operator_ = poperatorT(
      CoulombOperatorPtr(world, coulomb_lo_thresh_, potential_thresh_));

  density_ = computeDensity();
}

Tensor<double> ResponsePreliminariesCalculator::computeKineticEnergy() {
  auto phi = copy(world_, ground_state_.orbitals, true);
  reconstruct(world_, phi);

  real_derivative_3d Dx(world_, 0), Dy(world_, 1), Dz(world_, 2);
  auto fx = apply(world_, Dx, phi);
  auto fy = apply(world_, Dy, phi);
  auto fz = apply(world_, Dz, phi);

  compress(world_, fx, true);
  compress(world_, fy, true);
  compress(world_, fz, true);
  compress(world_, phi, true);
  world_.gop.fence();

  Tensor<double> T =
      0.5 * (matrix_inner(world_, fx, fx) + matrix_inner(world_, fy, fy) +
             matrix_inner(world_, fz, fz));

  phi.clear();
  fx.clear();
  fy.clear();
  fz.clear();

  return T;
}

real_function_3d ResponsePreliminariesCalculator::computeDensity() {
  auto phi_squared = square(world_, ground_state_.orbitals);
  compress(world_, phi_squared);
  real_function_3d rho = sum(world_, phi_squared);
  rho.truncate();
  phi_squared.clear();
  return rho;
}

real_function_3d ResponsePreliminariesCalculator::computeNuclearPotential() {
  potential_manager_->make_nuclear_potential(world_);
  auto v_nuc = potential_manager_->vnuclear();
  v_nuc.truncate(potential_thresh_);
  return v_nuc;
}

real_function_3d ResponsePreliminariesCalculator::computeCoulombPotential() {
  auto v_coul = 2.0 * apply(*coulomb_operator_, density_);
  v_coul.truncate(potential_thresh_);
  return v_coul;
}

real_function_3d ResponsePreliminariesCalculator::computeXCPotential() {
  XCOperator<double, 3> xc_operator{world_, ground_state_.xc,
                                    ground_state_.spinrestricted, density_,
                                    density_};
  auto v_xc = xc_operator.make_xc_potential();
  return v_xc;
}

vector_real_function_3d
ResponsePreliminariesCalculator::computeHFExchangeEnergy() {
  const double lo = 1.e-10;
  Exchange<double, 3> k{world_, lo};

  auto phi = ground_state_.orbitals;
  k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
  k.set_bra_and_ket(phi, phi);
  // k.set_symmetric(true).set_printlevel(r_params.print_level());
  auto kphi = k(phi);
  auto excv = inner(world_, kphi, phi);
  double exchf = 0.0;
  for (size_t i = 0; i < phi.size(); ++i) {
    exchf -= 0.5 * excv[static_cast<long>(i)];
  }
  auto Vpsi = -xcf_.hf_exchange_coefficient() * kphi;

  return Vpsi;
}

ResponsePreliminaries ResponsePreliminariesCalculator::computeAll() {
  ResponsePreliminaries prelim;

  auto phi = ground_state_.orbitals;

  // Compute kinetic energy once
  auto T = computeKineticEnergy();

  // Compute and combine nuclear and Coulomb potentials once
  auto V_nuc = computeNuclearPotential();
  auto V_coul = computeCoulombPotential();

  prelim.V_local = V_nuc + V_coul;

  // Conditionally compute XC potential for DFT
  if (xcf_.is_dft() && (xcf_.hf_exchange_coefficient() != 1.0)) {
    auto V_xc = computeXCPotential();
    prelim.V_local += V_xc;
  }

  // Conditionally compute HF exchange potential
  vector_real_function_3d V_total =
      zero_functions<double, 3>(world_, phi.size());
  if (xcf_.hf_exchange_coefficient() > 0.0) {
    V_total = computeHFExchangeEnergy();
  }

  auto vtol = FunctionDefaults<3>::get_thresh() * 0.1;
  gaxpy(world_, 1.0, V_total, 1.0,
        mul_sparse(world_, prelim.V_local, phi, vtol));
  truncate(world_, V_total);
  auto phi_V_phi = matrix_inner(world_, phi, V_total);

  prelim.Hamiltonian = T + phi_V_phi;
  auto new_hamiltonian_no_diag = copy(prelim.Hamiltonian);
  for (size_t i = 0; i < phi.size(); i++)
    new_hamiltonian_no_diag(long(i), long(i)) = 0.0;
  prelim.Hamiltonian_no_diag = new_hamiltonian_no_diag;

  return prelim;
}
