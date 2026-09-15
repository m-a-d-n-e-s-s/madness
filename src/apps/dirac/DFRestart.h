/// \file DFRestart.h
/// \brief Format version for DFdriver restart archives, and its validation.

#ifndef MADNESS_APPS_DIRAC_DFRESTART_H_INCLUDED
#define MADNESS_APPS_DIRAC_DFRESTART_H_INCLUDED

#include <madness/world/madness_exception.h>

namespace madness {

/// Format version for DFdriver restart archives.
///
/// Earlier archives carry no version field, and this code rejects them. The read
/// goes through the type cookie that the BinaryFstream archives write.
///
/// Current Version:
/// 0) DF restart format version (unsigned int)
/// 1) Total energy (double)
/// 2) Krestricted (boolean)
/// 3) closed_shell (boolean)
/// 3) number of occupied orbitals (int)
/// 4) orbital energies (vector of doubles)
/// 5) box size (double)
/// 6) wavelet order (int)
/// 7) molecule (molecule)
/// 8) occupied orbitals as complex functions
///
/// /note v1 introduced format version
inline constexpr unsigned int DF_RESTART_VERSION = 1;

/// Writes the current DF restart format version as the first datum in the archive.
template <typename Archive>
void write_df_restart_version(const Archive& ar) {
    const unsigned int version = DF_RESTART_VERSION;
    ar & version;
}

/// Returns the format version, or 0 for an unversioned archive.
///
/// Reads from a *local* archive, so the caller decides which rank reads.
template <typename Archive>
unsigned int read_df_restart_version(const Archive& ar) {
    unsigned int version = 0;
    try {
        ar & version;
    } catch (const MadnessException&) {
        // Type cookie mismatch: the first datum is a double, so the archive
        // has no version. The archive prints its own message to stderr
        // before it throws. The what() of that exception points at a stack
        // buffer, and this code must not read it.
        return 0;
    }
    return version;
}

/// Throws unless \c version is a DF restart format this build can read.
inline void require_supported_df_restart(const unsigned int version) {
    MADNESS_CHECK_THROW(version != 0,
        "unversioned DF restart archive! start from a moldft archive or "
        "rerun with the older DFdriver");
    MADNESS_CHECK_THROW(version == DF_RESTART_VERSION,
        "unsupported DF restart archive version! DF restart archive was written"
        " by a newer DFdriver; use that DFdriver or start from a moldft "
        "archive");
}

}  // namespace madness

#endif
