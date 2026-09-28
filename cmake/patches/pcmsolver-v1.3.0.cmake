# Patch script for the PCMSolver source tree fetched by FindOrFetchPCMSolver.cmake.
#
# Run as `cmake -DPCMSOLVER_SOURCE_DIR=<dir> -P pcmsolver-v1.3.0.cmake`, from
# FetchContent's PATCH_COMMAND. It is idempotent: a hunk already applied is
# skipped, and a hunk whose context is missing is a hard error -- so bumping
# MADNESS_TRACKED_PCMSOLVER_TAG past the point where upstream moves any of this
# fails loudly instead of silently building something untested.
#
# PCMSolver v1.3.0 is from 2020 and receives no maintenance. Hunks 1-4 are
# toolchain-compatibility fixes; the later hunks make the library usable at
# 10^4 tesserae, the cavity of a protein, and leave its surface charges unchanged.

cmake_minimum_required(VERSION 3.12.0)

if (NOT PCMSOLVER_SOURCE_DIR)
  message(FATAL_ERROR "PCMSOLVER_SOURCE_DIR must be set")
endif ()

# Replace `old` with `new` in `relpath`, exactly once.
function(pcm_patch relpath old new what)
  set(_file "${PCMSOLVER_SOURCE_DIR}/${relpath}")
  if (NOT EXISTS "${_file}")
    message(FATAL_ERROR "pcmsolver patch: ${relpath} does not exist -- "
                        "the fetched PCMSolver is not the pinned v1.3.0")
  endif ()
  file(READ "${_file}" _contents)
  string(FIND "${_contents}" "${new}" _already)
  if (NOT _already EQUAL -1)
    return()  # already patched
  endif ()
  string(FIND "${_contents}" "${old}" _found)
  if (_found EQUAL -1)
    message(FATAL_ERROR "pcmsolver patch: could not apply '${what}' to ${relpath} -- "
                        "upstream text has changed; revisit this patch script")
  endif ()
  string(REPLACE "${old}" "${new}" _contents "${_contents}")
  file(WRITE "${_file}" "${_contents}")
endfunction()

# 1. CMake >= 4.0 refuses a project whose cmake_minimum_required is below 3.5.
pcm_patch(CMakeLists.txt
    "cmake_minimum_required(VERSION 3.3 FATAL_ERROR)"
    "cmake_minimum_required(VERSION 3.5 FATAL_ERROR)"
    "raise cmake_minimum_required to 3.5")

# 2. Boost >= 1.87 requires C++14 (boost::numeric::odeint uses generic lambdas),
#    and SphericalDiffuse.cpp pulls odeint in.
pcm_patch(cmake/custom/compilers/CXXFlags.cmake
    "set(CMAKE_CXX_STANDARD 11)"
    "set(CMAKE_CXX_STANDARD 14)"
    "build PCMSolver as C++14")

# 3. The bundled Eigen 3.3.2 does not compile with clang >= 11: `trt` is a
#    Transpose<TranspositionsBase<...>>, which has no derived(). Fixed upstream
#    in Eigen 3.3.5; applied here because PCMSolver vendors the header.
pcm_patch(external/eigen3/include/eigen3/Eigen/src/Core/Transpositions.h
    "Product<OtherDerived, Transpose, AliasFreeProduct>(matrix.derived(), trt.derived())"
    "Product<OtherDerived, Transpose, AliasFreeProduct>(matrix.derived(), trt)"
    "Eigen 3.3.2 Transpositions.h clang fix")

# 4. `update_version` is far too generic a target name to add to a build tree
#    that MADNESS may itself be a subproject of (see the CMake harness notes in
#    CLAUDE.md). Nothing else in PCMSolver refers to it.
pcm_patch(cmake/custom/pcmsolver.cmake
    "add_custom_target(update_version"
    "add_custom_target(pcmsolver-update-version"
    "rename the update_version target")

# 5. Both solvers refactorize their N x N block matrix at every surface-charge
#    evaluation (`.lu().solve(...)` in computeCharge_impl), an O(N^3) LU per SCF
#    iteration: 4 s at 3374 tesserae, minutes at 10^4. Factorize once when the
#    matrices are built and keep the factorizations.
pcm_patch(src/solver/IEFSolver.hpp
    "#include <Eigen/Core>"
    "#include <Eigen/Core>
#include <Eigen/LU>"
    "IEFSolver.hpp: include Eigen/LU")
pcm_patch(src/solver/IEFSolver.hpp
    "  /*! R_infinity matrix, symmetry blocked form */
  std::vector<Eigen::MatrixXd> blockRinfinity_;"
    "  /*! R_infinity matrix, symmetry blocked form */
  std::vector<Eigen::MatrixXd> blockRinfinity_;
  /*! LU factorizations of the T(epsilon) blocks, and of their adjoints when
   *  hermitivitize_ is set: built once with the matrices, reused by every
   *  computeCharge_impl call instead of refactorizing there */
  std::vector<Eigen::PartialPivLU<Eigen::MatrixXd> > blockTepsilonLU_;
  std::vector<Eigen::PartialPivLU<Eigen::MatrixXd> > blockTepsilonAdjLU_;"
    "IEFSolver.hpp: cached LU factorizations")
pcm_patch(src/solver/IEFSolver.cpp
    "  utils::symmetryPacking(blockRinfinity_, Rinfinity_, dimBlock, nrBlocks);

  built_ = true;"
    "  utils::symmetryPacking(blockRinfinity_, Rinfinity_, dimBlock, nrBlocks);

  // factorize once; computeCharge_impl then only does the triangular solves
  blockTepsilonLU_.clear();
  blockTepsilonAdjLU_.clear();
  for (size_t i = 0; i < blockTepsilon_.size(); ++i) {
    blockTepsilonLU_.push_back(Eigen::PartialPivLU<Eigen::MatrixXd>(blockTepsilon_[i]));
    if (hermitivitize_)
      blockTepsilonAdjLU_.push_back(Eigen::PartialPivLU<Eigen::MatrixXd>(blockTepsilon_[i].adjoint()));
  }

  built_ = true;"
    "IEFSolver.cpp: factorize the T(epsilon) blocks at build time")
pcm_patch(src/solver/IEFSolver.cpp
    "-blockTepsilon_[irrep].lu().solve("
    "-blockTepsilonLU_[irrep].solve("
    "IEFSolver.cpp: reuse the LU in computeCharge_impl")
pcm_patch(src/solver/IEFSolver.cpp
    "blockTepsilon_[irrep].adjoint().lu().solve("
    "blockTepsilonAdjLU_[irrep].solve("
    "IEFSolver.cpp: reuse the adjoint LU in computeCharge_impl")
pcm_patch(src/solver/CPCMSolver.hpp
    "#include <Eigen/Core>"
    "#include <Eigen/Core>
#include <Eigen/LU>"
    "CPCMSolver.hpp: include Eigen/LU")
pcm_patch(src/solver/CPCMSolver.hpp
    "  /*! S matrix, symmetry blocked form */
  std::vector<Eigen::MatrixXd> blockS_;"
    "  /*! S matrix, symmetry blocked form */
  std::vector<Eigen::MatrixXd> blockS_;
  /*! LU factorizations of the S blocks, built once, reused by computeCharge_impl */
  std::vector<Eigen::PartialPivLU<Eigen::MatrixXd> > blockSLU_;"
    "CPCMSolver.hpp: cached LU factorizations")
pcm_patch(src/solver/CPCMSolver.cpp
    "  utils::symmetryPacking(blockS_, S_, dimBlock, nrBlocks);

  built_ = true;"
    "  utils::symmetryPacking(blockS_, S_, dimBlock, nrBlocks);

  // factorize once; computeCharge_impl then only does the triangular solves
  blockSLU_.clear();
  for (size_t i = 0; i < blockS_.size(); ++i)
    blockSLU_.push_back(Eigen::PartialPivLU<Eigen::MatrixXd>(blockS_[i]));

  built_ = true;"
    "CPCMSolver.cpp: factorize the S blocks at build time")
pcm_patch(src/solver/CPCMSolver.cpp
    "-blockS_[irrep].lu().solve("
    "-blockSLU_[irrep].solve("
    "CPCMSolver.cpp: reuse the LU in computeCharge_impl")

# 6. The cavity restart file cavity.npz is written at every construction, two
#    zip entries per tessera in append mode: O(N^2), 2.1 GB and the 2 GB zip
#    limit at 24k tesserae. Only `cavity_type restart` reads it, so write it
#    only when PCMSOLVER_SAVE_CAVITY is set.
pcm_patch(src/interface/Meddle.cpp
    "#include <string>
#include <vector>"
    "#include <cstdlib>
#include <string>
#include <vector>"
    "Meddle.cpp: include cstdlib")
pcm_patch(src/interface/Meddle.cpp
    "  cavity_->saveCavity();"
    "  // opt-in: the .npz is only read back by cavity_type = restart, and writing it
  // (append mode, two entries per tessera) is O(N^2) and fails past 2 GB
  if (std::getenv(\"PCMSOLVER_SAVE_CAVITY\")) cavity_->saveCavity();"
    "Meddle.cpp: cavity.npz on request only")

# 7. GePol's tesserae and vertex work arrays are sized at compile time (50,000
#    tesserae, 100,000 vertices) and the C++ side passes the same numbers; a
#    787-atom cavity at 0.3 A^2 has 29,467 tesserae and 107,375 vertices. Four
#    times the room. The sphere and centre limits (1000 atoms) stay, since the
#    arrays behind them are quadratic (DERCEN(MXSP,MXCENT,3,3)).
pcm_patch(src/pedra/pcm_pcmdef.inc
    "PARAMETER (MXTS=50000, MXSP=1000, MXTSPT = 2*MXTS)"
    "PARAMETER (MXTS=200000, MXSP=1000, MXTSPT = 2*MXTS)"
    "GePol: 200,000 tesserae")
pcm_patch(src/pedra/pcm_pcmdef.inc
    "PARAMETER (MXVER = 100000)"
    "PARAMETER (MXVER = 400000)"
    "GePol: 400,000 vertices")
pcm_patch(src/cavity/GePolCavity.cpp
    "build(suffix, 50000, 1000, 100000);"
    "build(suffix, 200000, 1000, 400000);"
    "GePolCavity.cpp: array sizes in step with pcm_pcmdef.inc")

# 8. Collocation::computeS_impl and computeD_impl copy an Element (two dynamic
#    Eigen matrices) for every pair (i, j): ~10^9 heap allocations for a 13.9k
#    tessera cavity, minutes per matrix. Take references. The loops are then
#    independent over i and run under OpenMP where the build has it.
pcm_patch(src/bi_operators/Collocation.cpp
    "    Element source = elems[i];"
    "    const Element & source = elems[i];"
    "Collocation.cpp: no Element copy per row")
pcm_patch(src/bi_operators/Collocation.cpp
    "      Element probe = elems[j];"
    "      const Element & probe = elems[j];"
    "Collocation.cpp: no Element copy per pair")
pcm_patch(src/bi_operators/Collocation.cpp
    "  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {"
    "  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
#pragma omp parallel for schedule(dynamic, 32)
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {"
    "Collocation.cpp: S fill in parallel")
pcm_patch(src/bi_operators/Collocation.cpp
    "  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {"
    "  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(cavitySize, cavitySize);
#pragma omp parallel for schedule(dynamic, 32)
  for (PCMSolverIndex i = 0; i < cavitySize; ++i) {"
    "Collocation.cpp: D fill in parallel")

# 9. The positive-definiteness check of S is an Eigen LDLT, which Eigen 3.3
#    implements unblocked and serially, O(N^3). The LLT fails exactly when S is
#    not positive-definite and is blocked on the parallel products.
pcm_patch(src/bi_operators/IBoundaryIntegralOperator.cpp
    "  Eigen::LDLT<Eigen::MatrixXd> Sldlt(biop);
  if (!Sldlt.isPositive()) {"
    "  Eigen::LLT<Eigen::MatrixXd> Sldlt(biop);
  if (Sldlt.info() != Eigen::Success) {"
    "IBoundaryIntegralOperator.cpp: LLT instead of LDLT for the S check")

# 10. The tessera-area diagonal is materialized as a dense N x N matrix and
#     multiplied as one: a dense O(N^3) product and N^2 doubles in every solver
#     build. Keep it diagonal.
pcm_patch(src/solver/SolverImpl.cpp
    "  Eigen::MatrixXd a = cav.elementArea().asDiagonal();"
    "  const Eigen::DiagonalMatrix<double, Eigen::Dynamic> a = cav.elementArea().asDiagonal();"
    "SolverImpl.cpp: keep the area matrix diagonal")
