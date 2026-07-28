# src/apps/molresponse/CMakeLists.txt

if(NOT TARGET MADchem)
  message(STATUS "molresponse: MADchem not available, skipping")
  return()
endif()

# Optional VTK for native HyperTreeGrid (.htg) export from dump_mra_trees (--htg).
# OFF by default so a stock build needs no VTK; turn on with -DMADNESS_ENABLE_VTK=ON
# (point CMake at the cluster VTK via -DVTK_DIR=...). Only dump_mra_trees links it.
option(MADNESS_ENABLE_VTK "Build dump_mra_trees with native VTK HyperTreeGrid (.htg) export" OFF)
if(MADNESS_ENABLE_VTK)
  find_package(VTK 9.0 REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)
  if(VTK_FOUND)
    message(STATUS "molresponse: VTK ${VTK_VERSION} found -> dump_mra_trees --htg enabled")
  endif()
endif()

add_mad_executable(molresponse "main.cpp;GroundState.cpp" "MADchem")
# Legacy direct-FD test binaries (test_v3_solver / test_es_solver) retired in R0b
# along with the top-level FDSolver.hpp / ESSolver.hpp they exercised; the solver
# core is covered by test_calc_manager_run / test_run_response / test_v3_*_skeleton.
add_mad_executable(test_v3_es_skeleton
                   "tests/test_v3_es_skeleton.cpp;GroundState.cpp"
                   "MADchem")
add_mad_executable(test_v3_fd_skeleton
                   "tests/test_v3_fd_skeleton.cpp;GroundState.cpp"
                   "MADchem")
add_mad_executable(test_protocol_key
                   "tests/test_protocol_key.cpp"
                   "MADchem")
add_mad_executable(test_response_metadata
                   "tests/test_response_metadata.cpp"
                   "MADchem")
add_mad_executable(test_response_property_planner
                   "tests/test_response_property_planner.cpp"
                   "MADchem")
add_mad_executable(test_state_metrics
                   "tests/test_state_metrics.cpp"
                   "MADchem")
add_mad_executable(test_calc_manager
                   "tests/test_calc_manager.cpp"
                   "MADchem")
add_mad_executable(test_calc_manager_run
                   "tests/test_calc_manager_run.cpp;GroundState.cpp"
                   "MADchem")
add_mad_executable(test_run_response
                   "tests/test_run_response.cpp;GroundState.cpp"
                   "MADchem")
add_mad_executable(test_kernel_equivalence
                   "tests/test_kernel_equivalence.cpp;GroundState.cpp"
                   "MADchem")
# Pins the three algebraic identities (M1/M2/M3) behind the state-grouped
# reformulation of the 2PA two-electron residue. Research-branch only.
add_mad_executable(test_tpa_regroup_identities
                   "tests/test_tpa_regroup_identities.cpp;GroundState.cpp"
                   "MADchem")
add_mad_executable(test_step_restriction
                   "tests/test_step_restriction.cpp"
                   "MADchem")
# k-CONSISTENCY regression (review fix, confirmed HIGH): try_load_fd_state must
# re-project a coarser-rung source to the active (k, thresh) inside the loader,
# so beta assembly can never contract mixed-k functions ("functions have
# different k" abort). Synthetic state, no ground-state fixture needed.
add_mad_executable(test_fd_kfallback_load
                   "tests/test_fd_kfallback_load.cpp"
                   "MADchem")
# Slot-identity contract of rs::diagonalize (dominance order, phase fix,
# rank-deficient dim padding) — guards against re-introducing the
# ascending-omega re-sort that caused the coarse-rung slot churn.
add_mad_executable(test_diagonalize_slots
                   "tests/test_diagonalize_slots.cpp"
                   "MADchem")
# Inc 1 of the ES node-φ state-parallel design (docs/21): node-aligned subworld
# creation (one subworld per physical node via Split_type SHARED) + map test.
add_mad_executable(test_node_subworlds
                   "tests/test_node_subworlds.cpp"
                   "MADchem")
# Inc 2 (docs/21): prove ground state shipped into node-aligned subworlds is
# ONE distributed copy per node (not per-rank replicated); de-risk cross-world copy.
add_mad_executable(test_node_phi
                   "tests/test_node_phi.cpp"
                   "MADchem")
# Keystone of the persistent-subworld recompute design (docs/22 §8): assemble the
# subspace matrix from per-node partial columns + universe allreduce (no function
# movement), == direct universe matrix.
add_mad_executable(test_subspace_allreduce
                   "tests/test_subspace_allreduce.cpp"
                   "MADchem")
# Inc 1 of the FD exchange tensor layer (docs 26/27/28): the tensor-layer θ
# (build_ctx + assemble_theta; shared Tx + cached g0, dot-contracted) matches the
# reference compute_V0x − compute_E0x + compute_gamma (Static ClosedShell).
add_mad_executable(test_exchange_ctx
                   "tests/test_exchange_ctx.cpp"
                   "MADchem")
# S1 of the persistent-subworld ES design (docs/23): prove the GS-dependent build
# kernels (V0x/γ/E0x/T0x → Λ) reproduce bit-for-bit in a node-aligned subworld,
# the one piece the φ-count + A/S-scalar proofs did not cover.
add_mad_executable(test_es_build_subworld
                   "tests/test_es_build_subworld.cpp"
                   "MADchem")
# F1 of the FD state-parallel design (docs/31,32): independent dipole FD states
# solved one-per-node-subworld give an alpha tensor bit-identical to the single-
# World path (pure fan-out -> solve -> filesystem gather). Needs GroundState.cpp
# (uses a real moldft archive, unlike the synthetic-GS subworld proofs above).
add_mad_executable(test_fd_subworld_fanout
                   "tests/test_fd_subworld_fanout.cpp;GroundState.cpp"
                   "MADchem")
# Offline MRA tree-skeleton dumper (R2/L4 structural export for viz/partitioning).
# Loads converged orbitals and dumps the adaptive octree per orbital; never solves.
add_mad_executable(dump_mra_trees
                   "tools/dump_mra_trees.cpp;GroundState.cpp"
                   "MADchem")
# tpa_from_dalton: contraction-only 2PA harness — reads DALTON's own converged
# vectors (RSPVEC: EXCITLAB X_f + DIPLEN N^A(w_f/2); molden.inp: occupied MOs),
# projects them to MRA, and passes them through tpa_sources/tpa_moment with ZERO
# solves. Equation-vs-vector discriminator against DALTON's printed S tensors
# (TPA_SCOPING.md §5k). No HDF5 needed (text molden + Fortran-binary RSPVEC).
add_mad_executable(tpa_from_dalton
                   "tools/tpa_from_dalton.cpp;GroundState.cpp"
                   "MADchem")
# Native HTG export is optional: only when -DMADNESS_ENABLE_VTK=ON found VTK.
# Guarded by MADNESS_HAS_VTK in the source; a no-VTK build is unchanged.
if(MADNESS_ENABLE_VTK AND VTK_FOUND)
  target_compile_definitions(dump_mra_trees PRIVATE MADNESS_HAS_VTK)
  target_link_libraries(dump_mra_trees PRIVATE ${VTK_LIBRARIES})
  vtk_module_autoinit(TARGETS dump_mra_trees MODULES ${VTK_LIBRARIES})
endif()
# --- io-hdf5 thread: HDF5-backed I/O (opt-in; docs/30) -----------------------
# Native MADNESS archive stays the default. No-op unless configured with
# -DMADNESS_ENABLE_HDF5=ON, so a normal build is byte-for-byte unchanged.
option(MADNESS_ENABLE_HDF5 "Enable HDF5-backed function/checkpoint I/O in molresponse (opt-in; native archive remains the default)" OFF)
if(MADNESS_ENABLE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C)
  message(STATUS "molresponse: HDF5 enabled — ${HDF5_LIBRARIES}")
  # P1 build-wiring smoke: pure HDF5 C-API round-trip, no MADNESS function code.
  add_mad_executable(test_hdf5_smoke "tests/test_hdf5_smoke.cpp" "MADworld")
  target_include_directories(test_hdf5_smoke PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_hdf5_smoke PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(test_hdf5_smoke PRIVATE MADNESS_HAS_HDF5)
  # P2 Layer A: Function raw-coeff round-trip vs the native archive.
  add_mad_executable(test_function_hdf5 "tests/test_function_hdf5.cpp;GroundState.cpp" "MADchem")
  target_include_directories(test_function_hdf5 PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_function_hdf5 PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(test_function_hdf5 PRIVATE MADNESS_HAS_HDF5)
  # P2a: ResponseStateX round-trip through the opt-in restart path.
  add_mad_executable(test_state_archive_hdf5 "tests/test_state_archive_hdf5.cpp" "MADchem")
  target_include_directories(test_state_archive_hdf5 PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_state_archive_hdf5 PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(test_state_archive_hdf5 PRIVATE MADNESS_HAS_HDF5)
  # Make the real solver binaries HDF5-aware so cm_run / madqc can use the opt-in
  # restart path in response_state.hpp. MADNESS_HAS_HDF5 only enables the path;
  # default (env MADRESPONSE_IO_HDF5 unset) still writes the legacy archive.
  # dump_mra_trees joins the list so its --coeffs native-coefficient export can call
  # save_function_hdf5 (the structured /meta+/keys+/coeffs writer in
  # solvers/function_hdf5_io.hpp); guarded by MADNESS_HAS_HDF5 in the source.
  foreach(_hdf5_tgt molresponse test_calc_manager_run dump_mra_trees)
    if(TARGET ${_hdf5_tgt})
      target_include_directories(${_hdf5_tgt} PRIVATE ${HDF5_INCLUDE_DIRS})
      target_link_libraries(${_hdf5_tgt} PRIVATE ${HDF5_LIBRARIES})
      target_compile_definitions(${_hdf5_tgt} PRIVATE MADNESS_HAS_HDF5)
    endif()
  endforeach()
  # dalton-cis-seed thread: converter that turns pymra per-occ-orbital .mad.h5
  # response functions into a TDA excited-state seed bundle (save_es_roots,
  # converged=false) the ES node Resumes from. HDF5-only (uses load_function_hdf5).
  add_mad_executable(seed_es_from_hdf5 "tools/seed_es_from_hdf5.cpp;GroundState.cpp" "MADchem")
  target_include_directories(seed_es_from_hdf5 PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(seed_es_from_hdf5 PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(seed_es_from_hdf5 PRIVATE MADNESS_HAS_HDF5)
  # dalton-cis-seed (b): overlap-based pinning — identify which converged ES root
  # matches the operator-selected seed by |<seed|root>| (character, not energy).
  add_mad_executable(es_overlap "tools/es_overlap.cpp;GroundState.cpp" "MADchem")
  target_include_directories(es_overlap PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(es_overlap PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(es_overlap PRIVATE MADNESS_HAS_HDF5)
endif()

# dalton-cis-seed (c): native DALTON reader — reads RSPVEC + molden.inp directly,
# projects x_i onto MRA, and emits the ES seed bundle without a Python/HDF5
# intermediate. Writes the NATIVE save_es_roots bundle (converged=false) that the
# ES node Resumes from — no HDF5 at runtime, and every HDF5 path in
# response_state.hpp is #ifdef MADNESS_HAS_HDF5 guarded, so this compiles in a
# plain (non-HDF5) build. Kept OUT of the MADNESS_ENABLE_HDF5 block so it is
# always available (enabling HDF5 here drags parallel-HDF5's Intel-MPI include
# path in and breaks the xeonmax toolchain — see es_bench MPI notes).
add_mad_executable(seed_from_dalton "tools/seed_from_dalton.cpp;GroundState.cpp" "MADchem")

# dalton-cis-seed (GS): project the occupied DALTON/molden MOs and write a
# moldft version-4 .restartdata (native path only, no HDF5). Ground-state analogue
# of seed_from_dalton; moldft --restart resumes from it (SCF::load_mos).
add_mad_executable(seed_moldft_from_dalton "tools/seed_moldft_from_dalton.cpp;GroundState.cpp" "MADchem")

add_dependencies(applications-madness molresponse)

install(TARGETS molresponse
        DESTINATION "${MADNESS_INSTALL_BINDIR}")

# --- CTest registration (M3): tiered so a bare `ctest` is fast and green on ---
# --- any machine, and heavier tiers are opt-in by label / cache variable. ----
if (BUILD_TESTING)
  # Tier 1 — "unit": pure logic, no fixtures, no MPI launcher required.
  foreach(_t test_calc_manager test_protocol_key test_response_metadata
             test_response_property_planner test_state_metrics
             test_diagonalize_slots)
    add_test(NAME madness/test/molresponse/${_t}/run COMMAND ${_t})
    set_property(TEST madness/test/molresponse/${_t}/run
                 PROPERTY LABELS "unit;response;short" APPEND)
  endforeach()

  # Tier 2 — "response-small": self-contained solver tests (SYNTHETIC ground
  # state; no fixture archives). Registered SERIAL on purpose: launcher choice
  # is cluster-specific (mpirun hangs in batch on some systems) and the
  # multi-rank coverage of these paths lives in the es_bench state-parallel
  # gates. Serial runs still validate the numerics end-to-end.
  foreach(_t test_exchange_ctx test_node_subworlds test_node_phi
             test_subspace_allreduce test_es_build_subworld test_step_restriction
             test_fd_kfallback_load)
    add_test(NAME madness/test/molresponse/${_t}/run COMMAND ${_t})
    set_property(TEST madness/test/molresponse/${_t}/run
                 PROPERTY LABELS "response-small;response;medium" APPEND)
  endforeach()

  # Tier 3 — "response-fixture": need a moldft ground-state archive; registered
  # only when the fixture tree is available so CI without fixtures stays green.
  # (test_fd_subworld_fanout is intentionally NOT here — it requires multi-rank
  # subworlds and is covered by the state-parallel gate.)
  set(MADNESS_RESPONSE_FIXTURE_DIR "" CACHE PATH
      "fixture tree for response regression tests (<dir>/h2o/mad.restartdata.*)")
  if (MADNESS_RESPONSE_FIXTURE_DIR)
    add_test(NAME madness/test/molresponse/test_kernel_equivalence/run
             COMMAND test_kernel_equivalence
                     --archive=${MADNESS_RESPONSE_FIXTURE_DIR}/h2o/mad.restartdata)
    add_test(NAME madness/test/molresponse/test_calc_manager_run/run
             COMMAND test_calc_manager_run
                     --archive=${MADNESS_RESPONSE_FIXTURE_DIR}/h2o/mad.restartdata
                     --omega=0.0 --axes=z --protocol=1e-4 --maxiter=25
                     --calc-dir=ctest_calc_manager_run)
    add_test(NAME madness/test/molresponse/test_run_response/run
             COMMAND test_run_response
                     --archive=${MADNESS_RESPONSE_FIXTURE_DIR}/h2o/mad.restartdata)
    foreach(_t test_kernel_equivalence test_calc_manager_run test_run_response)
      set_property(TEST madness/test/molresponse/${_t}/run
                   PROPERTY LABELS "response-fixture;response;long" APPEND)
    endforeach()
  endif()
endif (BUILD_TESTING)
