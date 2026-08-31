#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Standalone GTO evaluation and Molden parser for DALTON output.
// No MADNESS dependencies. Supports spherical shells l=0..4 (s,p,d,f,g).
// Normalization convention: each spherical harmonic component is normalized
// to unit self-overlap via combo_overlap (Schlegel & Frisch 1995).

// --- Harmonic tables ---------------------------------------------------------

struct GtoTerm { double coeff; int lx, ly, lz; };
using GtoHarmonic = std::vector<GtoTerm>;

// d: d0, d+1, d-1, d+2, d-2  (Molden spherical order, [5D] convention)
inline const std::vector<GtoHarmonic>& d_harmonics() {
    static const std::vector<GtoHarmonic> H = {
        // d0:  2z^2 - x^2 - y^2
        {{-0.5, 2, 0, 0}, {-0.5, 0, 2, 0}, {1.0, 0, 0, 2}},
        // d+1: xz
        {{1.0, 1, 0, 1}},
        // d-1: yz
        {{1.0, 0, 1, 1}},
        // d+2: x^2 - y^2
        {{0.5 * 1.7320508075688772, 2, 0, 0}, {-0.5 * 1.7320508075688772, 0, 2, 0}},
        // d-2: xy
        {{1.0, 1, 1, 0}},
    };
    return H;
}

// f: f0, f+1, f-1, f+2, f-2, f+3, f-3  (Schlegel & Frisch 1995, Table I)
inline const std::vector<GtoHarmonic>& f_harmonics() {
    static const std::vector<GtoHarmonic> H = {
        // f0
        {{2.0, 0, 0, 3}, {-3.0, 2, 0, 1}, {-3.0, 0, 2, 1}},
        // f+1
        {{4.0, 1, 0, 2}, {-1.0, 3, 0, 0}, {-1.0, 1, 2, 0}},
        // f-1
        {{4.0, 0, 1, 2}, {-1.0, 2, 1, 0}, {-1.0, 0, 3, 0}},
        // f+2
        {{1.0, 2, 0, 1}, {-1.0, 0, 2, 1}},
        // f-2
        {{1.0, 1, 1, 1}},
        // f+3
        {{1.0, 3, 0, 0}, {-3.0, 1, 2, 0}},
        // f-3
        {{3.0, 2, 1, 0}, {-1.0, 0, 3, 0}},
    };
    return H;
}

// g: g0, g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
inline const std::vector<GtoHarmonic>& g_harmonics() {
    static const std::vector<GtoHarmonic> H = {
        // g0
        {{3.0, 4, 0, 0}, {3.0, 0, 4, 0}, {8.0, 0, 0, 4},
         {6.0, 2, 2, 0}, {-24.0, 2, 0, 2}, {-24.0, 0, 2, 2}},
        // g+1
        {{4.0, 1, 0, 3}, {-3.0, 3, 0, 1}, {-3.0, 1, 2, 1}},
        // g-1
        {{4.0, 0, 1, 3}, {-3.0, 2, 1, 1}, {-3.0, 0, 3, 1}},
        // g+2
        {{-1.0, 4, 0, 0}, {1.0, 0, 4, 0}, {6.0, 2, 0, 2}, {-6.0, 0, 2, 2}},
        // g-2
        {{6.0, 1, 1, 2}, {-1.0, 3, 1, 0}, {-1.0, 1, 3, 0}},
        // g+3
        {{1.0, 3, 0, 1}, {-3.0, 1, 2, 1}},
        // g-3
        {{3.0, 2, 1, 1}, {-1.0, 0, 3, 1}},
        // g+4
        {{1.0, 4, 0, 0}, {-6.0, 2, 2, 0}, {1.0, 0, 4, 0}},
        // g-4
        {{1.0, 3, 1, 0}, {-1.0, 1, 3, 0}},
    };
    return H;
}

inline const std::vector<GtoHarmonic>& harmonics_for_l(int l) {
    if (l == 2) return d_harmonics();
    if (l == 3) return f_harmonics();
    if (l == 4) return g_harmonics();
    throw std::runtime_error("harmonics_for_l: l=" + std::to_string(l) +
                             " has no harmonic table (use s/p path)");
}

// --- Normalization helpers ---------------------------------------------------

static inline double double_factorial(int n) {
    if (n <= 0) return 1.0;
    double r = 1.0;
    while (n > 1) { r *= n; n -= 2; }
    return r;
}

// <terms_a|terms_b> for contracted Gaussian with shared exponents/raw_coeffs.
// Analytic closed-form integral; odd (lx,ly,lz) sums vanish by symmetry.
static inline double combo_overlap(const std::vector<double>& exps,
                                   const std::vector<double>& raw_coeffs,
                                   const GtoHarmonic& terms_a,
                                   const GtoHarmonic& terms_b) {
    double total = 0.0;
    const int np = static_cast<int>(exps.size());
    for (const auto& ta : terms_a) {
        for (const auto& tb : terms_b) {
            int lx = ta.lx + tb.lx;
            int ly = ta.ly + tb.ly;
            int lz = ta.lz + tb.lz;
            if (lx % 2 || ly % 2 || lz % 2) continue;
            double df = double_factorial(lx - 1) *
                        double_factorial(ly - 1) *
                        double_factorial(lz - 1);
            int ltot2 = (lx + ly + lz) / 2;
            for (int p = 0; p < np; p++) {
                for (int q = 0; q < np; q++) {
                    double ab = exps[p] + exps[q];
                    double mom = std::pow(M_PI / ab, 1.5) * df /
                                 std::pow(2.0 * ab, static_cast<double>(ltot2));
                    total += ta.coeff * tb.coeff *
                             raw_coeffs[p] * raw_coeffs[q] * mom;
                }
            }
        }
    }
    return total;
}

// --- DaltonSphericalShell ---------------------------------------------------

struct DaltonSphericalShell {
    double cx, cy, cz;
    int l;
    std::vector<double> exponents;
    std::vector<double> raw_coeffs;
    std::vector<double> norms;   // one per harmonic component (1/sqrt(combo_overlap))
    int n_ao;                    // 2l+1 for spherical (1 for s, 3 for p, ...)

    DaltonSphericalShell(double cx_, double cy_, double cz_, int l_,
                         std::vector<double> exps, std::vector<double> coeffs)
        : cx(cx_), cy(cy_), cz(cz_), l(l_),
          exponents(std::move(exps)), raw_coeffs(std::move(coeffs)) {
        static const int n_ao_table[5] = {1, 3, 5, 7, 9};
        assert(l >= 0 && l <= 4);
        n_ao = n_ao_table[l];
        norms.resize(static_cast<size_t>(n_ao));

        // Molden contraction coefficients multiply UNIT-NORMALIZED primitives,
        // whose norm depends on the exponent: N(a) ∝ (2a/pi)^{3/4} (4a)^{l/2}.
        // Fold that a-dependence into the coefficients BEFORE contracting —
        // otherwise every contracted shell (nprim>1) has wrong RELATIVE
        // primitive weights, which the per-component renormalization below
        // cannot repair (it only fixes the overall scale). This bug broke MO
        // orthonormality at the 0.3 level on aug-cc-pVXZ water (caught by
        // tpa_from_dalton's C^T·S·C import-fidelity check, 2026-07-21); the
        // constant (2l-1)!!-type factor is per-shell and absorbed by norms[].
        for (size_t p = 0; p < exponents.size(); ++p)
            raw_coeffs[p] *= std::pow(2.0 * exponents[p] / M_PI, 0.75) *
                             std::pow(4.0 * exponents[p], 0.5 * l);

        if (l == 0) {
            GtoHarmonic h = {{1.0, 0, 0, 0}};
            double ov = combo_overlap(exponents, raw_coeffs, h, h);
            norms[0] = 1.0 / std::sqrt(ov);
        } else if (l == 1) {
            // px, py, pz
            const int comps[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
            for (int k = 0; k < 3; k++) {
                GtoHarmonic h = {{1.0, comps[k][0], comps[k][1], comps[k][2]}};
                double ov = combo_overlap(exponents, raw_coeffs, h, h);
                norms[static_cast<size_t>(k)] = 1.0 / std::sqrt(ov);
            }
        } else {
            const auto& harmon = harmonics_for_l(l);
            for (int k = 0; k < n_ao; k++) {
                double ov = combo_overlap(exponents, raw_coeffs,
                                          harmon[static_cast<size_t>(k)],
                                          harmon[static_cast<size_t>(k)]);
                norms[static_cast<size_t>(k)] = 1.0 / std::sqrt(ov);
            }
        }
    }

    // Evaluate all n_ao harmonic components at Cartesian point (x, y, z).
    // Writes n_ao values into bf[]; caller must provide at least n_ao doubles.
    void evaluate(double x, double y, double z, double* bf) const {
        double dx = x - cx, dy = y - cy, dz = z - cz;
        double r2 = dx*dx + dy*dy + dz*dz;

        // Shared radial: sum_p raw_coeffs[p] * exp(-exps[p]*r2)
        double radial = 0.0;
        for (size_t p = 0; p < exponents.size(); p++)
            radial += raw_coeffs[p] * std::exp(-exponents[p] * r2);

        if (l == 0) {
            bf[0] = norms[0] * radial;
        } else if (l == 1) {
            bf[0] = norms[0] * dx * radial;
            bf[1] = norms[1] * dy * radial;
            bf[2] = norms[2] * dz * radial;
        } else {
            const auto& harmon = harmonics_for_l(l);
            for (int k = 0; k < n_ao; k++) {
                double poly = 0.0;
                for (const auto& t : harmon[static_cast<size_t>(k)]) {
                    double mono = t.coeff;
                    for (int ix = 0; ix < t.lx; ix++) mono *= dx;
                    for (int iy = 0; iy < t.ly; iy++) mono *= dy;
                    for (int iz = 0; iz < t.lz; iz++) mono *= dz;
                    poly += mono;
                }
                bf[static_cast<size_t>(k)] = norms[static_cast<size_t>(k)] * poly * radial;
            }
        }
    }
};

// --- Molden data structures -------------------------------------------------

struct DaltonMoldenBasis {
    std::vector<DaltonSphericalShell> shells;
    std::vector<int> ao_offsets;   // first AO index for each shell
    int n_ao;
};

struct DaltonMoldenResult {
    std::vector<std::string> atom_symbols;
    std::vector<int> atomic_numbers;             // [Atoms] column 3
    std::vector<std::array<double, 3>> coords;   // in Bohr
    DaltonMoldenBasis basis;
    std::vector<double> mo_coeffs;     // column-major: C[mu + n_ao*mo_idx]
    std::vector<double> mo_energies;
    std::vector<double> mo_occ;
    int n_ao;
    int n_mo;
};

// --- Molden parser ----------------------------------------------------------

namespace dalton_gto_detail {

inline std::string to_upper(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

inline std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Find (start_line, end_line) of a named block (e.g. "GTO", "Atoms", "MO").
// Block starts at a line matching ^\s*\[name(\]|\s) (case-insensitive).
// Ends at the next line starting with '[' or end of file.
// Returns {-1,-1} if not found.
inline std::pair<int,int> block_bounds(const std::vector<std::string>& lines,
                                        const std::string& name) {
    std::regex pat("^\\s*\\[" + name + "(\\]|\\s)", std::regex::icase);
    int start = -1;
    for (int i = 0; i < static_cast<int>(lines.size()); i++) {
        if (std::regex_search(lines[static_cast<size_t>(i)], pat)) {
            start = i;
            break;
        }
    }
    if (start < 0) return {-1, -1};
    int end = static_cast<int>(lines.size());
    for (int i = start + 1; i < static_cast<int>(lines.size()); i++) {
        const std::string t = trim(lines[static_cast<size_t>(i)]);
        if (!t.empty() && t[0] == '[') { end = i; break; }
    }
    return {start, end};
}

} // namespace dalton_gto_detail

// Read a DALTON Molden file (standard [Molden Format] / [Atoms] / [GTO] / [MO]).
// Handles [5D7F] and [9G] markers (pure spherical; they are recorded and skipped).
// MO coefficients in the [MO] block are sparse (1-based index, only non-zero listed).
inline DaltonMoldenResult read_molden(const std::string& path) {
    using namespace dalton_gto_detail;
    std::ifstream fh(path);
    if (!fh.is_open())
        throw std::runtime_error("read_molden: cannot open file: " + path);

    std::vector<std::string> lines;
    {
        std::string line;
        while (std::getline(fh, line)) lines.push_back(line);
    }

    // ---- [Atoms] block ----
    auto [ab_start, ab_end] = block_bounds(lines, "Atoms");
    if (ab_start < 0)
        throw std::runtime_error("read_molden: no [Atoms] block in " + path);
    {
        std::string hdr = to_upper(trim(lines[static_cast<size_t>(ab_start)]));
        // unit_angstrom: header contains "ANGS" (not just "AU")
    }
    bool unit_angstrom = false;
    {
        std::string hdr = to_upper(trim(lines[static_cast<size_t>(ab_start)]));
        unit_angstrom = (hdr.find("ANGS") != std::string::npos);
        // If only "AU" present and no "ANGS", it's Bohr (DALTON default).
    }

    std::vector<std::string> atom_symbols;
    std::vector<int> atomic_numbers;
    std::vector<std::array<double, 3>> coords;
    for (int i = ab_start + 1; i < ab_end; i++) {
        std::istringstream ss(lines[static_cast<size_t>(i)]);
        std::string sym;
        int idx, anum;
        double x, y, z;
        // Format: symbol  idx  atomic_num  x  y  z
        if (!(ss >> sym >> idx >> anum >> x >> y >> z)) continue;
        atom_symbols.push_back(sym);
        atomic_numbers.push_back(anum);
        if (unit_angstrom) {
            const double bohr_per_ang = 1.8897259886;
            coords.push_back({x * bohr_per_ang, y * bohr_per_ang, z * bohr_per_ang});
        } else {
            coords.push_back({x, y, z});
        }
    }

    // ---- [GTO] block ----
    auto [gb_start, gb_end] = block_bounds(lines, "GTO");
    if (gb_start < 0)
        throw std::runtime_error("read_molden: no [GTO] block in " + path);

    std::vector<DaltonSphericalShell> shells;
    std::vector<int> ao_offsets;
    int n_ao = 0;
    int atom_idx = 0;
    const std::string ang_labels = "spdfg";

    for (int i = gb_start + 1; i < gb_end; ) {
        const std::string line = trim(lines[static_cast<size_t>(i)]);
        if (line.empty()) { i++; continue; }
        // Skip [5D7F], [9G], etc. block markers that may appear inside [GTO]
        if (line[0] == '[') { i++; continue; }

        std::istringstream ss(line);
        std::string tok1, tok2;
        ss >> tok1;

        // Atom index line: first token is a digit string
        bool is_digit = !tok1.empty();
        for (char c : tok1) if (!std::isdigit(static_cast<unsigned char>(c))) { is_digit = false; break; }
        if (is_digit) {
            atom_idx = std::stoi(tok1) - 1;  // 1-based -> 0-based
            i++;
            continue;
        }

        // Shell line: "s 3 1.00" (angular label, nprim, scale)
        std::string ang_str = tok1;
        // Lowercase
        for (auto& c : ang_str) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (ang_str.size() != 1 || ang_labels.find(ang_str[0]) == std::string::npos) {
            i++;
            continue;
        }
        int l = static_cast<int>(ang_labels.find(ang_str[0]));
        int nprim = 0;
        ss >> nprim;

        std::vector<double> exps, coeffs;
        for (int j = 0; j < nprim; j++) {
            std::string pline = lines[static_cast<size_t>(i + 1 + j)];
            // Replace Fortran D/d exponents with E
            for (auto& c : pline) if (c == 'D' || c == 'd') c = 'E';
            std::istringstream pss(pline);
            double ex, co;
            pss >> ex >> co;
            exps.push_back(ex);
            coeffs.push_back(co);
        }

        const auto& crd = coords[static_cast<size_t>(atom_idx)];
        shells.emplace_back(crd[0], crd[1], crd[2], l, exps, coeffs);
        ao_offsets.push_back(n_ao);
        n_ao += shells.back().n_ao;
        i += 1 + nprim;
    }

    // ---- [MO] block ----
    auto [mb_start, mb_end] = block_bounds(lines, "MO");
    if (mb_start < 0)
        throw std::runtime_error("read_molden: no [MO] block in " + path);

    std::vector<double> mo_coeffs_all;  // column-major
    std::vector<double> mo_energies;
    std::vector<double> mo_occ;

    for (int i = mb_start + 1; i < mb_end; ) {
        const std::string line = trim(lines[static_cast<size_t>(i)]);
        if (line.empty()) { i++; continue; }

        std::string lup = to_upper(line);
        if (lup.substr(0, 4) == "SYM=") {
            // Start of a new MO
            double ene = 0.0, occ = 0.0;
            i++;
            // Read scalar header fields (lines containing '=')
            while (i < mb_end) {
                const std::string fl = trim(lines[static_cast<size_t>(i)]);
                if (fl.find('=') == std::string::npos) break;
                std::string flup = to_upper(fl);
                if (flup.substr(0, 4) == "ENE=") {
                    ene = std::stod(fl.substr(fl.find('=') + 1));
                } else if (flup.substr(0, 6) == "OCCUP=") {
                    occ = std::stod(fl.substr(fl.find('=') + 1));
                }
                // Spin= lines are skipped (contain '=', not parsed further)
                i++;
            }
            // Read sparse coefficient lines (idx coeff, no '=')
            std::vector<double> col(static_cast<size_t>(n_ao), 0.0);
            while (i < mb_end) {
                const std::string fl = trim(lines[static_cast<size_t>(i)]);
                if (fl.empty() || fl.find('=') != std::string::npos ||
                    to_upper(fl).substr(0, 4) == "SYM=") break;
                std::istringstream css(fl);
                int idx;
                double coeff;
                if (!(css >> idx >> coeff)) { i++; continue; }
                if (idx >= 1 && idx <= n_ao)
                    col[static_cast<size_t>(idx - 1)] = coeff;
                i++;
            }
            for (double v : col) mo_coeffs_all.push_back(v);
            mo_energies.push_back(ene);
            mo_occ.push_back(occ);
        } else {
            i++;
        }
    }

    int n_mo = static_cast<int>(mo_energies.size());

    DaltonMoldenBasis basis;
    basis.shells = std::move(shells);
    basis.ao_offsets = std::move(ao_offsets);
    basis.n_ao = n_ao;

    DaltonMoldenResult res;
    res.atom_symbols = std::move(atom_symbols);
    res.atomic_numbers = std::move(atomic_numbers);
    res.coords = std::move(coords);
    res.basis = std::move(basis);
    res.mo_coeffs = std::move(mo_coeffs_all);
    res.mo_energies = std::move(mo_energies);
    res.mo_occ = std::move(mo_occ);
    res.n_ao = n_ao;
    res.n_mo = n_mo;
    return res;
}
