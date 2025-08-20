#include "GroundStateData.hpp"

using namespace madness;

GroundStateData::GroundStateData(World &world, std::string archiveFile, Molecule mol)
    : molecule(std::move(mol)), spinrestricted(false), num_orbitals(0), L(0.0), k(0), xc(), localize_method(),
      converged_for_thresh(0.0), original_k(0), archive_File(std::move(archiveFile)) {

  if (world.rank() == 0) {
    print("Constructing GroundStateData for molecule: ");
    molecule.print();
  }
  load(world);
}

GroundStateData::GroundStateData(World &world, std::shared_ptr<SCF> scf_calc)
    : molecule(scf_calc->molecule), spinrestricted(scf_calc->is_spin_restricted()),
      num_orbitals(scf_calc->param.nalpha()), L(scf_calc->param.L()), k(FunctionDefaults<3>::get_k()), xc(scf_calc->param.xc()),
      localize_method(scf_calc->param.localize_method()), converged_for_thresh(FunctionDefaults<3>::get_thresh()),
      original_k(FunctionDefaults<3>::get_k()) {

  if (world.rank() == 0) {
    print("Constructing GroundStateData from SCF calculation for molecule: ");
    molecule.print();
  }
  orbitals = scf_calc->get_amo();
}

void GroundStateData::load(World &world) {
  archive::ParallelInputArchive input(world, archive_File.c_str());

  unsigned int version;
  double dummy_energy;
  std::vector<int> dummy_sets;

  input & version;
  input & dummy_energy;
  input & spinrestricted;
  input & L;
  input & k;
  input & molecule;
  input & xc;
  input & localize_method;

  if (version > 3) {
    input & converged_for_thresh;
  }

  input & num_orbitals;
  input & energies;
  input & occ;
  input & dummy_sets;

  if (k < 1 || k > 30) {
    if (world.rank() == 0)
      madness::print("Invalid wavelet order read from archive, setting to default k=8.");
    k = 8;
  }

  // save the current k
  int current_k = FunctionDefaults<3>::get_k();

  // change Default k to the one read from the archive so that read orbitals
  // is done correctly
  FunctionDefaults<3>::set_k(k);

  world.gop.fence();
  orbitals.clear();
  world.gop.fence();

  for (unsigned int i = 0; i < num_orbitals; ++i) {
    real_function_3d orbital;
    input & orbital;
    orbitals.push_back(orbital);
  }

  world.gop.fence();
  truncate(world, orbitals);

  original_k = k;
  FunctionDefaults<3>::set_truncate_mode(1);
  // Set the default k back to the original value
  FunctionDefaults<3>::set_k(current_k);
}

void GroundStateData::print_info() const {
  madness::print("\nGround State Orbital Information:");
  madness::print("------------------------");
  madness::print("XC Functional:", xc);
  madness::print("Localization Method:", localize_method);
  madness::print("Spin Restricted:", spinrestricted);
  madness::print("Number of Orbitals:", num_orbitals);
  madness::print("Box Size L:", L);
  madness::print("Wavelet Order k:", k);
  madness::print("Converged for Threshold:", converged_for_thresh);
  madness::print("Orbital Energies:", energies);
}

void GroundStateData::prepareOrbitals(World &world, int current_k, double thresh) {

  // Check if the current orbitals' polynomial order is higher (projectable)
  if (original_k < current_k) {
    if (world.rank() == 0) {
      print("Cannot project orbitals: current k (", k, ") is lower than requested k (", current_k, ").");
    }
    MADNESS_EXCEPTION("Orbital polynomial order too low for projection.", k);
  }
  bool needs_reload = false;

  // Check if orbitals need reloading due to polynomial order mismatch or
  // convergence threshold
  if (k != current_k || converged_for_thresh > thresh) {
    if (world.rank() == 0) {
      print("Reloading orbitals: current k (", k, ") is different from requested k (", current_k, ").");
    }
    needs_reload = true;
  }

  if (needs_reload) {
    load(world);
    reconstruct(world, orbitals);

    for (auto &orbital : orbitals) {
      orbital = project(orbital, current_k, thresh, true);
    }
    truncate(world, orbitals, thresh);

    // Update orbital metadata
    k = current_k;
    converged_for_thresh = thresh;
  } else {
    truncate(world, orbitals, thresh);
  }
  if (world.rank() == 0) {
    print_info();
  }
  Qhat = QProjector<double, 3>(orbitals);
}
// TODO: Need to add functionality for spin restricted
Tensor<double> GroundStateData::tryLoadHamiltonianFromJson(World &world, const json &fock_json, double thresh, int k) {
  std::pair<Tensor<double>, Tensor<double>> fock_data;
  auto protocol = std::string("thresh: ") + std::to_string(thresh) + std::string(" k: ") + std::to_string(k);
  // Check if the protocol exists in the JSON
  if (fock_json.find(protocol) == fock_json.end()) {
    std::cerr << "Protocol not found in JSON: " << protocol << std::endl;
    return {};
  } else {
    fock_data.first = tensor_from_json<double>(fock_json[protocol]["focka"]);
    if (!spinrestricted) {
      fock_data.second = tensor_from_json<double>(fock_json[protocol]["fockb"]);
    }
    return fock_data.first;
  }
}

// TODO: Need to add functionality for spin restricted
// TODO: use fock_a and fock_b
void GroundStateData::computePreliminaries(World &world, const operatorT &coulop, double vtol,
                                           const std::string &fock_json_file) {

  auto thresh = FunctionDefaults<3>::get_thresh();
  auto k = FunctionDefaults<3>::get_k();
  density_ = computeDensity(world);

  potential_manager_ = std::make_shared<PotentialManager>(molecule, "a");
  potential_manager_->make_nuclear_potential(world);
  world.gop.fence();
  world.gop.fence();
  xcf_.initialize(xc, spinrestricted, world, true);

  auto V_nuc = computeNuclearPotential(vtol);
  auto V_coul = computeCoulombPotential(coulop, vtol);
  V_local = V_nuc + V_coul;
  if (xcf_.is_dft() && (xcf_.hf_exchange_coefficient() != 1.0)) {
    auto V_xc = computeXCPotential(world);
    V_local += V_xc;
  }

  V_nuc.clear();
  V_coul.clear();

  bool json_exists = std::filesystem::exists(fock_json_file);
  bool loaded_from_file = false;

  if (json_exists) {
    auto fock_json = broadcast_json_file(world, fock_json_file);
    world.gop.fence();
    Hamiltonian = tryLoadHamiltonianFromJson(world, fock_json, thresh, k);
  }
  loaded_from_file = Hamiltonian.has_data();

  if (!loaded_from_file) {
    auto T = computeKineticEnergy(world);
    vector_real_function_3d V_hf_phi = zero_functions<double, 3>(world, orbitals.size());
    if (xcf_.hf_exchange_coefficient() > 0.0) {
      V_hf_phi = computeHFExchangeEnergy(world);
    }

    auto V_local_phi = mul_sparse(world, V_local, orbitals, vtol);
    auto V_phi = gaxpy_oop(1.0, V_local_phi, 1.0, V_hf_phi);
    truncate(world, V_phi);
    auto phi_V_phi = matrix_inner(world, orbitals, V_phi);

    Hamiltonian = T + phi_V_phi;
  }
  auto new_hamiltonian_no_diag = copy(Hamiltonian);
  for (size_t i = 0; i < orbitals.size(); i++)
    new_hamiltonian_no_diag(long(i), long(i)) = 0.0;
  Hamiltonian_no_diag = new_hamiltonian_no_diag;
}

real_function_3d GroundStateData::computeDensity(World &world) const {
  auto phi_squared = square(world, orbitals);
  compress(world, phi_squared);
  real_function_3d rho = sum(world, phi_squared);
  rho.truncate();
  phi_squared.clear();
  return rho;
}

real_function_3d GroundStateData::computeNuclearPotential(double vtol) const {
  auto v_nuc = potential_manager_->vnuclear();
  v_nuc.truncate(vtol);
  return v_nuc;
}

[[nodiscard]] real_function_3d GroundStateData::computeCoulombPotential(const operatorT &coulop, double vtol) const {
  auto v_coul = 2.0 * apply(coulop, density_, true);
  v_coul.truncate(vtol);
  return v_coul;
}

real_function_3d GroundStateData::computeXCPotential(World &world) const {
  XCOperator<double, 3> xc_operator{world, xc, spinrestricted, density_, density_};
  auto v_xc = xc_operator.make_xc_potential();
  return v_xc;
}

vector_real_function_3d GroundStateData::computeHFExchangeEnergy(World &world) const {
  const double lo = 1.e-10;
  Exchange<double, 3> k{world, lo};

  auto phi = copy(world, orbitals, true);
  k.set_algorithm(Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
  k.set_bra_and_ket(phi, phi);
  // k.set_symmetric(true).set_printlevel(r_params.print_level());
  auto kphi = k(phi);
  auto Vpsi = -xcf_.hf_exchange_coefficient() * kphi;

  return Vpsi;
}
Tensor<double> GroundStateData::computeKineticEnergy(World &world) const {
  auto phi = copy(world, orbitals, true);
  reconstruct(world, phi);

  real_derivative_3d Dx(world, 0);
  real_derivative_3d Dy(world, 1);
  real_derivative_3d Dz(world, 2);

  auto fx = apply(world, Dx, phi);
  auto fy = apply(world, Dy, phi);
  auto fz = apply(world, Dz, phi);

  compress(world, fx, true);
  compress(world, fy, true);
  compress(world, fz, true);
  compress(world, phi, true);
  world.gop.fence();

  Tensor<double> T = 0.5 * (matrix_inner(world, fx, fx) + matrix_inner(world, fy, fy) + matrix_inner(world, fz, fz));

  phi.clear();
  fx.clear();
  fy.clear();
  fz.clear();

  return T;
}
