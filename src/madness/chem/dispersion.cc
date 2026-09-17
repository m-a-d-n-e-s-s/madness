/*
 This file is part of MADNESS.

 Copyright (C) 2007,2010 Oak Ridge National Laboratory

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

 For more information please contact:

 Robert J. Harrison
 Oak Ridge National Laboratory
 One Bethel Valley Road
 P.O. Box 2008, MS-6367

 email: harrisonrj@ornl.gov
 tel:   865-241-3937
 fax:   865-572-0680
 */

#include <madness/madness_config.h>
#include <madness/chem/dispersion.h>

#include <madness/world/worldgop.h>

#include <algorithm>
#include <cctype>
#include <sstream>
#include <vector>

#ifdef MADNESS_HAS_DFTD3
// s-dftd3.h, not dftd3.h: the latter is a one-line wrapper and shares its name
// with Grimme's older standalone dftd3 program.
//
// Included through an absolute path baked in by CMake, not through an -I: the
// library usually lives in a conda prefix whose include/ also holds headers
// MADNESS resolves elsewhere (nlohmann/json.hpp), and an -I here would shadow
// them in this translation unit alone. Same trick as MADNESS_MPI_HEADER.
#ifdef MADNESS_DFTD3_HEADER
#include MADNESS_DFTD3_HEADER
#else
#include <s-dftd3.h>
#endif
#endif

namespace madness {

namespace {

std::string lowercase(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

/// the first whitespace-separated token of a line, lowercased ("" if none)
std::string first_token(const std::string& line) {
    std::istringstream is(line);
    std::string token;
    is >> token;
    return lowercase(token);
}

#ifdef MADNESS_HAS_DFTD3

/// turn a set error handle into a MADNESS exception

/// The message is printed rather than handed to MADNESS_EXCEPTION because that
/// macro stores the `const char*` it is given without copying it.
void throw_on_error(dftd3_error error, const char* what) {
    if (not dftd3_check_error(error)) return;
    char buffer[512];
    const int buflen = sizeof(buffer);
    dftd3_get_error(error, buffer, &buflen);
    print("\nsimple-dftd3 failed while", what);
    print(" ", buffer);
    MADNESS_EXCEPTION("simple-dftd3 reported an error -- see the message above", 1);
}

/// load the damping parameters for `method`, or explain why it could not be done
dftd3_param load_damping(dftd3_error error, const bool rational,
                         const std::string& method, const bool atm) {
    // the C API takes a mutable char*
    std::vector<char> name(method.begin(), method.end());
    name.push_back('\0');
    dftd3_param param = rational
        ? dftd3_load_rational_damping(error, name.data(), atm)
        : dftd3_load_zero_damping(error, name.data(), atm);
    if (dftd3_check_error(error)) {
        char buffer[512];
        const int buflen = sizeof(buffer);
        dftd3_get_error(error, buffer, &buflen);
        print("\nsimple-dftd3 has no D3 damping parameters for the functional '"
              + method + "'");
        print(" ", buffer);
        print(" set `dispersion_functional` in the dft group to a name it knows,");
        print(" e.g. pbe, pbe0, b3lyp, blyp, bp86, tpss, revpbe, hf");
        MADNESS_EXCEPTION("no D3 damping parameters for the requested functional", 1);
    }
    return param;
}

#endif // MADNESS_HAS_DFTD3

} // anonymous namespace


bool DispersionCorrection::available() {
#ifdef MADNESS_HAS_DFTD3
    return true;
#else
    return false;
#endif
}


std::string DispersionCorrection::library_version() {
#ifdef MADNESS_HAS_DFTD3
    const int v = dftd3_get_version();
    return std::to_string(v / 10000) + "." + std::to_string((v / 100) % 100)
           + "." + std::to_string(v % 100);
#else
    return "";
#endif
}


DispersionCorrection::DispersionCorrection(const std::string& spec,
                                           const std::string& xc_line,
                                           const std::string& functional,
                                           const bool atm1)
    : atm(atm1) {

    const std::string s = first_token(spec);
    if (s.empty() or s == "none") {
        damping = none;
        atm = false;
        return;
    } else if (s == "d3bj") {
        damping = rational;
    } else if (s == "d3zero") {
        damping = zero;
    } else {
        print("\nunknown dispersion correction '" + spec + "'");
        print(" known values are: none, d3bj, d3zero");
        MADNESS_EXCEPTION("unknown value for the `dispersion` keyword", 1);
    }

    // The `xc` line is free-form (an alias like `pbe0`, or a raw libxc spec like
    // `GGA_X_PBE 1.0 GGA_C_PBE 1.0`), so its first token is a guess. It is only
    // a fallback; `dispersion_functional` overrides it, and an unusable guess
    // fails loudly below rather than silently picking the wrong parameters.
    // "none" is the parameter default, spelled the way pcm_data/ac_data/nwfile
    // spell "unset" -- an empty default value cannot be read back out of
    // QCCalculationParametersBase.
    const std::string requested = first_token(functional);
    method = (requested.empty() or requested == "none") ? first_token(xc_line)
                                                        : requested;
    if (method.empty()) {
        print("\nno functional name to look up D3 damping parameters with");
        print(" set `dispersion_functional` in the dft group");
        MADNESS_EXCEPTION("no functional name for the D3 damping parameters", 1);
    }

#ifdef MADNESS_HAS_DFTD3
    // fail fast: resolve the parameters now, at construction, rather than at the
    // first energy evaluation. This runs on every rank and touches no MPI, so a
    // bad name throws collectively.
    dftd3_error error = dftd3_new_error();
    dftd3_param param = load_damping(error, damping == rational, method, atm);
    dftd3_delete_param(&param);
    dftd3_delete_error(&error);
#else
    print("\n`dispersion " + spec + "` was requested, but this MADNESS was built");
    print(" without simple-dftd3. Install it (e.g. `conda install -c conda-forge");
    print(" simple-dftd3`) and reconfigure with -DENABLE_DFTD3=ON, pointing");
    print(" -DDFTD3_ROOT_DIR at the install prefix.");
    MADNESS_EXCEPTION("dispersion correction requested, but MADNESS was built without simple-dftd3", 1);
#endif
}


void DispersionCorrection::reject(const std::string& spec, const char* engine) {
    const std::string s = first_token(spec);
    if (s.empty() or s == "none") return;
    print("\n`dispersion " + spec + "` was requested, but", engine,
          "has no dispersion term in its energy expression.");
    print(" Remove `dispersion` from the dft group, or use a workflow that");
    print(" supports it (scf, nemo).");
    MADNESS_EXCEPTION("dispersion correction requested on an engine that does not support it", 1);
}


std::string DispersionCorrection::description() const {
    if (damping == none) return "none";
    std::string d = (damping == rational) ? "D3(BJ)" : "D3(0)";
    if (atm) d += "-ATM";
    return d + "/" + method;
}


void DispersionCorrection::print_citation(World& world) const {
    if (damping == none or citation_printed) return;
    citation_printed = true;
    if (world.rank() != 0) return;

    print("");
    print(" Empirical dispersion correction:", description());
    print("   S. Grimme, J. Antony, S. Ehrlich, H. Krieg,");
    print("     J. Chem. Phys. 132, 154104 (2010)      doi:10.1063/1.3382344");
    if (damping == rational) {
        print("   S. Grimme, S. Ehrlich, L. Goerigk,          [Becke-Johnson damping]");
        print("     J. Comput. Chem. 32, 1456 (2011)       doi:10.1002/jcc.21759");
    }
    print("   implementation: simple-dftd3 " + library_version()
          + ", https://github.com/dftd3/simple-dftd3");
    print("");
}


void DispersionCorrection::evaluate(const std::vector<int>& numbers,
                                    const Tensor<double>& coords,
                                    double& e, Tensor<double>& g) const {
#ifdef MADNESS_HAS_DFTD3
    dftd3_error error = dftd3_new_error();

    // no lattice, no periodicity -- molecular
    dftd3_structure structure = dftd3_new_structure(
            error, int(numbers.size()), numbers.data(), coords.ptr(), nullptr, nullptr);
    throw_on_error(error, "constructing the molecular structure");

    dftd3_model model = dftd3_new_d3_model(error, structure);
    throw_on_error(error, "constructing the D3 model");

    dftd3_param param = load_damping(error, damping == rational, method, atm);

    // energy in Hartree, gradient in Hartree/bohr laid out [3*atom + axis],
    // which is exactly what SCF::derivatives uses -- no reordering, no unit
    // conversion. The virial is meaningless for a molecular system but is
    // still written to: passing a null pointer segfaults inside the library.
    double sigma[9];
    dftd3_get_dispersion(error, structure, model, param, &e, g.ptr(), sigma);
    throw_on_error(error, "evaluating the dispersion energy and gradient");

    dftd3_delete_param(&param);
    dftd3_delete_model(&model);
    dftd3_delete_structure(&structure);
    dftd3_delete_error(&error);
#else
    (void)numbers; (void)coords; (void)e; (void)g;
    MADNESS_EXCEPTION("dispersion correction requested, but MADNESS was built without simple-dftd3", 1);
#endif
}


/// the atomic numbers, in the order simple-dftd3 expects the coordinates
std::vector<int> DispersionCorrection::atomic_numbers(const Molecule& mol) {
    std::vector<int> numbers(mol.natom());
    for (size_t i = 0; i < mol.natom(); ++i)
        numbers[i] = int(mol.get_atomic_number(static_cast<unsigned int>(i)));
    return numbers;
}


void DispersionCorrection::compute(World& world, const Molecule& mol) const {

    const Tensor<double> coords = mol.get_all_coords();  // natom x 3, in bohr
    if (cache_valid and cached_coords.size() == coords.size()
        and (cached_coords - coords).absmax() == 0.0) return;

    const long natom = static_cast<long>(mol.natom());
    double e = 0.0;
    Tensor<double> g(3 * natom);

    // Evaluated on one rank and broadcast rather than redundantly everywhere:
    // simple-dftd3 may be OpenMP-parallel, and its reduction order -- hence the
    // last bits of the result -- is not guaranteed to agree across ranks running
    // with different thread counts.
    if (world.rank() == 0) evaluate(atomic_numbers(mol), coords, e, g);

    world.gop.broadcast(e, 0);
    world.gop.broadcast(g.ptr(), g.size(), 0);

    cached_coords = copy(coords);
    cached_gradient = g;
    cached_energy = e;
    cache_valid = true;
}


Tensor<double> DispersionCorrection::hessian(World& world, const Molecule& mol) const {

    const long natom = static_cast<long>(mol.natom());
    const long ndim = 3 * natom;
    Tensor<double> h(ndim, ndim);
    if (damping == none) return h;

    // The whole finite-difference loop runs on rank 0 and the matrix is
    // broadcast once, rather than going through compute() per displacement:
    // 6*natom collectives inside a loop for a quantity this cheap is the wrong
    // trade, and it would evict the energy/gradient cache on every step.
    if (world.rank() == 0) {
        const std::vector<int> numbers = atomic_numbers(mol);
        Tensor<double> coords = copy(mol.get_all_coords());

        // 1e-4 bohr: the gradient being differenced is analytic, so the only
        // errors are the O(h^2) truncation and the O(eps/h) cancellation, which
        // cross near this value and leave ~1e-9 Ha/bohr^2.
        const double step = 1.e-4;
        double e = 0.0;
        Tensor<double> gplus(ndim), gminus(ndim);

        for (long i = 0; i < ndim; ++i) {
            const long iatom = i / 3;
            const long iaxis = i % 3;
            const double x0 = coords(iatom, iaxis);

            coords(iatom, iaxis) = x0 + step;
            evaluate(numbers, coords, e, gplus);
            coords(iatom, iaxis) = x0 - step;
            evaluate(numbers, coords, e, gminus);
            coords(iatom, iaxis) = x0;  // exact restore, not x0-step+step

            for (long j = 0; j < ndim; ++j)
                h(i, j) = (gplus[j] - gminus[j]) / (2.0 * step);
        }

        // The exact Hessian is symmetric; the finite-difference one is only
        // nearly so, and averaging the two off-diagonal estimates is free
        // accuracy. Done explicitly rather than via transpose() to keep the
        // in-place update obvious.
        for (long i = 0; i < ndim; ++i) {
            for (long j = 0; j < i; ++j) {
                const double avg = 0.5 * (h(i, j) + h(j, i));
                h(i, j) = avg;
                h(j, i) = avg;
            }
        }
    }

    world.gop.broadcast(h.ptr(), h.size(), 0);
    return h;
}


double DispersionCorrection::energy(World& world, const Molecule& mol) const {
    if (damping == none) return 0.0;
    compute(world, mol);
    return cached_energy;
}


Tensor<double> DispersionCorrection::gradient(World& world, const Molecule& mol) const {
    if (damping == none) return Tensor<double>(static_cast<long>(3 * mol.natom()));
    compute(world, mol);
    return copy(cached_gradient);
}

} // namespace madness
