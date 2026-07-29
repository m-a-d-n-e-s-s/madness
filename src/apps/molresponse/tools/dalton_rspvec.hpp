#pragma once

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// RSPVEC binary reader: Fortran unformatted sequential (gfortran, little-endian,
// 4-byte record markers). Layout pinned from DALTON/rsp/rspcr5.F (WRTRSP/REARSP)
// and DALTON/sirius/sirset.F (JWOP ordering for the C1 case).

struct RspVecInfo {
    int32_t nish[8];
    int32_t nash[8];
    int32_t norb[8];
    int32_t nbas[8];
    int32_t nsym;
};

struct RspVecEntry {
    char lab1[9];   // null-terminated, 8 chars from file
    char lab2[9];
    double freq1;
    double freq2;
    int32_t isym1;
    int32_t isym2;
    double antsym;
    double rsd;
    int32_t length;
    double emcscf;
    int32_t nbast;
    int32_t norbt;
    std::vector<double> vec;
};

// Reads one Fortran unformatted record (4-byte leading + trailing markers).
// Returns payload bytes. Returns empty vector at clean EOF (zero bytes for
// leading marker). Throws on truncation or marker mismatch.
inline std::vector<char> read_fortran_record(std::ifstream& f) {
    int32_t n_lead = 0;
    f.read(reinterpret_cast<char*>(&n_lead), 4);
    if (f.gcount() == 0) return {};
    if (f.gcount() != 4)
        throw std::runtime_error("RSPVEC: truncated leading record marker");
    std::vector<char> payload(static_cast<size_t>(n_lead));
    if (n_lead > 0) {
        f.read(payload.data(), n_lead);
        if (static_cast<int32_t>(f.gcount()) != n_lead)
            throw std::runtime_error("RSPVEC: truncated record payload");
    }
    int32_t n_trail = 0;
    f.read(reinterpret_cast<char*>(&n_trail), 4);
    if (f.gcount() != 4)
        throw std::runtime_error("RSPVEC: missing trailing record marker");
    if (n_trail != n_lead)
        throw std::runtime_error("RSPVEC: record marker mismatch");
    return payload;
}

// Read all entries from a DALTON RSPVEC file.
// Returns (RspVecInfo, entries). Stops at EOFLABEL or physical EOF.
inline std::pair<RspVecInfo, std::vector<RspVecEntry>>
read_rspvec(const std::string& path) {
    std::ifstream fh(path, std::ios::binary);
    if (!fh.is_open())
        throw std::runtime_error("RSPVEC: cannot open file: " + path);

    // Record 1: 33 little-endian int32_t
    auto rec1 = read_fortran_record(fh);
    if (rec1.empty())
        throw std::runtime_error("RSPVEC: empty file (no record 1)");
    if (static_cast<int>(rec1.size()) != 33 * 4)
        throw std::runtime_error("RSPVEC: record 1 size mismatch");

    RspVecInfo info;
    const int32_t* ints = reinterpret_cast<const int32_t*>(rec1.data());
    for (int i = 0; i < 8; i++) info.nish[i] = ints[i];
    for (int i = 0; i < 8; i++) info.nash[i] = ints[8 + i];
    for (int i = 0; i < 8; i++) info.norb[i] = ints[16 + i];
    for (int i = 0; i < 8; i++) info.nbas[i] = ints[24 + i];
    info.nsym = ints[32];

    std::vector<RspVecEntry> entries;

    // Header record layout (packed, 76 bytes total, offsets shown):
    //   0:  lab1     8 chars
    //   8:  lab2     8 chars
    //  16:  freq1    double
    //  24:  freq2    double
    //  32:  isym1    int32_t
    //  36:  isym2    int32_t
    //  40:  antsym   double
    //  48:  rsd      double
    //  56:  length   int32_t
    //  60:  emcscf   double
    //  68:  nbast    int32_t
    //  72:  norbt    int32_t
    //  76: (end)
    constexpr int HEADER_SIZE = 76;

    while (true) {
        auto header = read_fortran_record(fh);
        if (header.empty()) break;

        // EOF sentinel: exactly 8 bytes containing "EOFLABEL"
        if (header.size() == 8) {
            if (std::string(header.data(), 8) == "EOFLABEL") break;
            throw std::runtime_error("RSPVEC: unexpected 8-byte record (not EOFLABEL)");
        }
        if (static_cast<int>(header.size()) != HEADER_SIZE)
            throw std::runtime_error("RSPVEC: header record has wrong size");

        RspVecEntry e;
        const char* p = header.data();

        std::memcpy(e.lab1, p + 0, 8);  e.lab1[8] = '\0';
        std::memcpy(e.lab2, p + 8, 8);  e.lab2[8] = '\0';
        // Trim trailing spaces
        for (int i = 7; i >= 0 && e.lab1[i] == ' '; i--) e.lab1[i] = '\0';
        for (int i = 7; i >= 0 && e.lab2[i] == ' '; i--) e.lab2[i] = '\0';

        std::memcpy(&e.freq1,  p + 16, 8);
        std::memcpy(&e.freq2,  p + 24, 8);
        std::memcpy(&e.isym1,  p + 32, 4);
        std::memcpy(&e.isym2,  p + 36, 4);
        std::memcpy(&e.antsym, p + 40, 8);
        std::memcpy(&e.rsd,    p + 48, 8);
        std::memcpy(&e.length, p + 56, 4);
        std::memcpy(&e.emcscf, p + 60, 8);
        std::memcpy(&e.nbast,  p + 68, 4);
        std::memcpy(&e.norbt,  p + 72, 4);

        auto data = read_fortran_record(fh);
        if (e.length > 0) {
            if (static_cast<int>(data.size()) != e.length * 8)
                throw std::runtime_error("RSPVEC: data record size mismatch");
            e.vec.resize(static_cast<size_t>(e.length));
            std::memcpy(e.vec.data(), data.data(), static_cast<size_t>(e.length) * 8);
        }
        entries.push_back(std::move(e));
    }
    return {info, entries};
}

// Reshape a flat RSPVEC vector into (X, Y) in row-major (n_occ, n_vir) order.
// TDA/CIS vectors have length n_occ*n_vir -> Y is empty.
// Full RPA vectors have length 2*n_occ*n_vir -> first half is X, second is Y.
// Raises on unexpected length (e.g. NSYM>1 multi-block file).
inline std::pair<std::vector<double>, std::vector<double>>
split_ov(const std::vector<double>& vec, int n_occ, int n_vir) {
    const int n_ov = n_occ * n_vir;
    if (static_cast<int>(vec.size()) == n_ov) {
        return {std::vector<double>(vec.begin(), vec.end()), {}};
    }
    if (static_cast<int>(vec.size()) == 2 * n_ov) {
        return {
            std::vector<double>(vec.begin(), vec.begin() + n_ov),
            std::vector<double>(vec.begin() + n_ov, vec.end())
        };
    }
    throw std::runtime_error(
        "split_ov: vector length " + std::to_string(vec.size()) +
        " matches neither n_occ*n_vir (" + std::to_string(n_ov) +
        ") nor 2*n_occ*n_vir (" + std::to_string(2 * n_ov) +
        ") -- multi-symmetry (NSYM>1) file not supported");
}
