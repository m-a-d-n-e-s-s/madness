// =========================================================================
// Unit test for the dalton.dir import contract (solvers/dalton_import.hpp,
// showcase W3). Pure filesystem + parsing — no MPI, no MADNESS World, no
// fixtures: writes a synthetic molden.inp / RSPVEC / .out (and a tarball)
// into a scratch dir. Covers:
//   1. locate by convention: loose molden.inp + RSPVEC
//   2. locate: .out discovery (unique), multiple *.out -> hard error,
//      explicit override wins
//   3. locate: multiple *.molden.inp with no molden.inp -> hard error
//   4. locate: missing RSPVEC -> hard error naming what was looked for
//   5. tarball fallback: RSPVEC/molden.inp extracted from a unique *.tar.gz
//   6. geometry fingerprint: match / coordinate mismatch (max_dev) /
//      element mismatch / atom-count mismatch
//   7. geometry hash: stable across re-locate; differs for a moved atom
//   8. basis + method parsed from the .out (HF and KS-DFT variants)
//   9. exact 1-to-1 frequency matching: all-present -> "", one absent ->
//      message lists BOTH the missing request and the available ladder
// =========================================================================

#include "../solvers/dalton_import.hpp"

#include <unistd.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

namespace {

int failed = 0;

#define EXPECT(cond, label)                                                \
  do {                                                                     \
    if (cond) { std::printf("  [PASS]  %s\n", label); }                    \
    else      { std::printf("  [FAIL]  %s\n", label); ++failed; }          \
  } while (0)

namespace fs = std::filesystem;

// A minimal-but-valid DALTON molden file: H2-like, two s-only atoms (bohr).
const char *kMolden =
    "[Molden Format]\n"
    "[Atoms] AU\n"
    "O_a        1     8         0.0000000000         0.0000000000         0.2216648690\n"
    "H_a        2     1         0.0000000000         1.4309006196        -0.8866594990\n"
    "[GTO]\n"
    "  1 0\n"
    " s    1 1.00\n"
    "     1.0000000000e+00  1.0000000000e+00\n"
    "\n"
    "  2 0\n"
    " s    1 1.00\n"
    "     8.0000000000e-01  1.0000000000e+00\n"
    "\n"
    "[MO]\n"
    " Sym= 1a\n"
    " Ene= -0.5\n"
    " Spin= Alpha\n"
    " Occup= 2.0\n"
    "   1   0.70\n"
    "   2   0.30\n"
    " Sym= 2a\n"
    " Ene= 0.2\n"
    " Spin= Alpha\n"
    " Occup= 0.0\n"
    "   1   0.30\n"
    "   2  -0.70\n";

void write_file(const fs::path &p, const std::string &s) {
  std::ofstream(p) << s;
}

// Fortran unformatted sequential record (4-byte little-endian markers).
void frec(std::ofstream &f, const void *data, std::int32_t n) {
  f.write(reinterpret_cast<const char *>(&n), 4);
  f.write(reinterpret_cast<const char *>(data), n);
  f.write(reinterpret_cast<const char *>(&n), 4);
}

// Synthetic RSPVEC: n_occ=1, n_vir=1 (norb=2), entries per (label, freq).
void write_rspvec(const fs::path &p,
                  const std::vector<std::pair<std::string, double>> &recs) {
  std::ofstream f(p, std::ios::binary);
  std::int32_t rec1[33] = {0};
  rec1[0]  = 1;   // NISH
  rec1[16] = 2;   // NORB
  rec1[24] = 2;   // NBAS
  rec1[32] = 1;   // NSYM
  frec(f, rec1, sizeof rec1);
  for (const auto &[lab, freq] : recs) {
    char hdr[76];
    std::memset(hdr, ' ', sizeof hdr);
    std::memcpy(hdr, lab.c_str(), lab.size());          // lab1 (8, blank-pad)
    double d;
    d = freq; std::memcpy(hdr + 16, &d, 8);             // freq1
    d = 0.0;  std::memcpy(hdr + 24, &d, 8);             // freq2
    std::int32_t i32;
    i32 = 1;  std::memcpy(hdr + 32, &i32, 4);           // isym1
    i32 = 1;  std::memcpy(hdr + 36, &i32, 4);           // isym2
    d = -1.0; std::memcpy(hdr + 40, &d, 8);             // antsym
    d = 0.0;  std::memcpy(hdr + 48, &d, 8);             // rsd
    i32 = 2;  std::memcpy(hdr + 56, &i32, 4);           // length (Z+Y, 1x1)
    d = 0.0;  std::memcpy(hdr + 60, &d, 8);             // emcscf
    i32 = 2;  std::memcpy(hdr + 68, &i32, 4);           // nbast
    i32 = 2;  std::memcpy(hdr + 72, &i32, 4);           // norbt
    frec(f, hdr, sizeof hdr);
    double vec[2] = {0.5, -0.5};
    frec(f, vec, sizeof vec);
  }
  const char eof[8] = {'E', 'O', 'F', 'L', 'A', 'B', 'E', 'L'};
  frec(f, eof, 8);
}

const char *kOutHF =
    "     Basis set used is \"aug-cc-pVTZ\" from the basis set library.\n"
    "@    Restricted, closed shell Hartree-Fock calculation.\n"
    "@    Wave function type        --- HF ---\n";

const char *kOutDFT =
    "     Basis set used is \"aug-cc-pVTZ\" from the basis set library.\n"
    "@    Wave function type        --- KS-DFT ---\n";

std::string err_of(const std::function<void()> &f) {
  try { f(); } catch (const std::exception &e) { return e.what(); }
  return {};
}

} // namespace

int main() {
  using namespace molresponse_v3;

  const fs::path tmp =
      fs::temp_directory_path() /
      ("test_dalton_import." + std::to_string(::getpid()));
  fs::remove_all(tmp);
  const fs::path ex = tmp / "extract";

  // ---- 1. locate by convention (loose files) ------------------------------
  std::printf("\n-- locate (loose) --\n");
  const fs::path d1 = tmp / "d1";
  fs::create_directories(d1);
  write_file(d1 / "molden.inp", kMolden);
  write_rspvec(d1 / "RSPVEC", {{"XDIPLEN", 0.0}, {"ZDIPLEN", 0.0596000625}});
  write_file(d1 / "run.out", kOutHF);
  {
    auto m = locate_dalton_dir(d1.string(), ex.string());
    EXPECT(m.molden_path == (d1 / "molden.inp").string(), "molden by convention");
    EXPECT(m.rspvec_path == (d1 / "RSPVEC").string(),     "RSPVEC by convention");
    EXPECT(m.out_path    == (d1 / "run.out").string(),    "unique .out discovered");
    EXPECT(m.basis  == "aug-cc-pVTZ", "basis parsed from .out");
    EXPECT(m.method == "HF",          "method parsed from .out (HF)");
    EXPECT(m.atomic_numbers.size() == 2 && m.atomic_numbers[0] == 8,
           "geometry parsed (2 atoms, Z=8 first)");
    EXPECT(m.geometry_hash.size() == 16, "geometry hash present (fnv1a64 hex)");

    auto m2 = locate_dalton_dir(d1.string(), ex.string());
    EXPECT(m2.geometry_hash == m.geometry_hash, "geometry hash stable");
  }

  // ---- 2. multiple .out -> hard error; override wins -----------------------
  std::printf("\n-- .out ambiguity --\n");
  write_file(d1 / "other.out", kOutHF);
  {
    const std::string e =
        err_of([&] { (void)locate_dalton_dir(d1.string(), ex.string()); });
    EXPECT(e.find("multiple *.out") != std::string::npos,
           "two .out files -> hard error naming the ambiguity");
    auto m = locate_dalton_dir(d1.string(), ex.string(), "", "",
                               (d1 / "run.out").string());
    EXPECT(m.out_path == (d1 / "run.out").string(), "explicit --dalton-out wins");
  }
  fs::remove(d1 / "other.out");

  // ---- 3. multiple *.molden.inp, no molden.inp -> hard error ---------------
  std::printf("\n-- molden ambiguity --\n");
  const fs::path d2 = tmp / "d2";
  fs::create_directories(d2);
  write_file(d2 / "a.molden.inp", kMolden);
  write_file(d2 / "b.molden.inp", kMolden);
  write_rspvec(d2 / "RSPVEC", {{"XDIPLEN", 0.0}});
  {
    const std::string e =
        err_of([&] { (void)locate_dalton_dir(d2.string(), ex.string()); });
    EXPECT(e.find("multiple molden") != std::string::npos,
           "two *.molden.inp -> hard error");
    auto m = locate_dalton_dir(d2.string(), ex.string(),
                               (d2 / "a.molden.inp").string());
    EXPECT(m.molden_path == (d2 / "a.molden.inp").string(),
           "explicit molden override wins");
  }

  // ---- 4. missing RSPVEC -> hard error -------------------------------------
  std::printf("\n-- missing RSPVEC --\n");
  const fs::path d3 = tmp / "d3";
  fs::create_directories(d3);
  write_file(d3 / "molden.inp", kMolden);
  {
    const std::string e =
        err_of([&] { (void)locate_dalton_dir(d3.string(), ex.string()); });
    EXPECT(e.find("no RSPVEC") != std::string::npos,
           "missing RSPVEC -> hard error naming what was looked for");
  }

  // ---- 5. tarball fallback --------------------------------------------------
  std::printf("\n-- tarball extraction --\n");
  const fs::path d4 = tmp / "d4";
  fs::create_directories(d4);
  {
    // Stage members elsewhere, tar them into d4 (dir itself has no loose files).
    const fs::path stage = tmp / "stage";
    fs::create_directories(stage);
    write_file(stage / "molden.inp", kMolden);
    write_rspvec(stage / "RSPVEC", {{"ZDIPLEN", 0.0596000625}});
    const std::string cmd = "tar -czf '" + (d4 / "quad_test.tar.gz").string() +
                            "' -C '" + stage.string() + "' RSPVEC molden.inp";
    EXPECT(std::system(cmd.c_str()) == 0, "staging tarball created");
    const fs::path ex4 = tmp / "extract4";
    auto m = locate_dalton_dir(d4.string(), ex4.string());
    EXPECT(m.molden_path == (ex4 / "molden.inp").string(),
           "molden extracted from unique *.tar.gz");
    EXPECT(m.rspvec_path == (ex4 / "RSPVEC").string(),
           "RSPVEC extracted from unique *.tar.gz");
  }

  // ---- 6. geometry fingerprint ----------------------------------------------
  std::printf("\n-- geometry fingerprint --\n");
  {
    auto m = locate_dalton_dir(d1.string(), ex.string(), "", "",
                               (d1 / "run.out").string());
    madness::Molecule good;
    good.add_atom(0.0, 0.0, 0.2216648690, 8.0, 8);
    good.add_atom(0.0, 1.4309006196, -0.8866594990, 1.0, 1);
    auto c = fingerprint_dalton_geometry(m, good);
    EXPECT(c.ok && c.max_dev < 1e-10, "identical geometry -> MATCH");
    EXPECT(c.report.find("MATCH") != std::string::npos, "report says MATCH");

    madness::Molecule shifted;
    shifted.add_atom(0.0, 0.0, 0.2226648690, 8.0, 8);   // +1e-3 bohr in z
    shifted.add_atom(0.0, 1.4309006196, -0.8866594990, 1.0, 1);
    auto cs = fingerprint_dalton_geometry(m, shifted);
    EXPECT(!cs.ok && std::abs(cs.max_dev - 1e-3) < 1e-9,
           "1e-3 bohr shift -> MISMATCH with max_dev = 1e-3");
    EXPECT(cs.report.find("MISMATCH") != std::string::npos,
           "report says MISMATCH");

    madness::Molecule wrong_el;
    wrong_el.add_atom(0.0, 0.0, 0.2216648690, 6.0, 6);   // C instead of O
    wrong_el.add_atom(0.0, 1.4309006196, -0.8866594990, 1.0, 1);
    auto ce = fingerprint_dalton_geometry(m, wrong_el);
    EXPECT(!ce.ok, "element mismatch -> MISMATCH");
    EXPECT(ce.report.find("ELEMENT MISMATCH") != std::string::npos,
           "report flags the element");

    madness::Molecule fewer;
    fewer.add_atom(0.0, 0.0, 0.2216648690, 8.0, 8);
    auto cn = fingerprint_dalton_geometry(m, fewer);
    EXPECT(!cn.ok, "atom-count mismatch -> MISMATCH");
    EXPECT(cn.report.find("ATOM-COUNT MISMATCH") != std::string::npos,
           "report flags the atom count");
  }

  // ---- 7. geometry hash sensitivity -----------------------------------------
  std::printf("\n-- geometry hash --\n");
  {
    const fs::path d5 = tmp / "d5";
    fs::create_directories(d5);
    std::string moved = kMolden;
    const auto pos = moved.find("0.2216648690");
    moved.replace(pos, std::strlen("0.2216648690"), "0.3216648690");
    write_file(d5 / "molden.inp", moved);
    write_rspvec(d5 / "RSPVEC", {{"XDIPLEN", 0.0}});
    auto ma = locate_dalton_dir(d1.string(), ex.string(), "", "",
                                (d1 / "run.out").string());
    auto mb = locate_dalton_dir(d5.string(), ex.string());
    EXPECT(ma.geometry_hash != mb.geometry_hash,
           "moved atom -> different geometry hash");
  }

  // ---- 8. DFT method parse ---------------------------------------------------
  std::printf("\n-- method parse --\n");
  {
    const fs::path d6 = tmp / "d6";
    fs::create_directories(d6);
    write_file(d6 / "molden.inp", kMolden);
    write_rspvec(d6 / "RSPVEC", {{"XDIPLEN", 0.0}});
    write_file(d6 / "dft.out", kOutDFT);
    auto m = locate_dalton_dir(d6.string(), ex.string());
    EXPECT(m.method == "KS-DFT",
           "KS-DFT wave function type parsed (rejected later by the gate)");
  }

  // ---- 9. exact 1-to-1 frequency matching ------------------------------------
  std::printf("\n-- frequency matching --\n");
  {
    auto rsp = read_rspvec((d1 / "RSPVEC").string());
    const auto &entries = rsp.second;
    EXPECT(match_dalton_frequencies(
               entries, {{0, 0.0}, {2, 0.0596000625}}, "RSPVEC").empty(),
           "all requested present -> no error");
    // representation noise within 1e-9 still matches
    EXPECT(match_dalton_frequencies(
               entries, {{2, 0.0596000625 + 5e-10}}, "RSPVEC").empty(),
           "sub-1e-9 representation noise tolerated");
    const std::string e =
        match_dalton_frequencies(entries, {{2, 0.057}}, "RSPVEC");
    EXPECT(!e.empty(), "absent frequency -> error");
    EXPECT(e.find("ZDIPLEN @ 0.057000000") != std::string::npos,
           "error names the missing request");
    EXPECT(e.find("0.059600063") != std::string::npos,
           "error lists the available ladder");
    EXPECT(e.find("nearest") != std::string::npos,
           "error states the no-nearest-fallback contract");
    // a frequency present for Z but requested for Y is still a mismatch
    const std::string ey =
        match_dalton_frequencies(entries, {{1, 0.0596000625}}, "RSPVEC");
    EXPECT(!ey.empty() && ey.find("YDIPLEN") != std::string::npos,
           "right frequency, wrong operator label -> error");
  }

  fs::remove_all(tmp);
  std::printf("\n%s: %d failure(s)\n",
              failed == 0 ? "ALL PASS" : "FAILED", failed);
  return failed == 0 ? 0 : 1;
}
