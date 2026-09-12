/// \file DFRestart.h
/// \brief Format version for DFdriver restart archives, and its validation.

#ifndef MADNESS_APPS_DIRAC_DFRESTART_H_INCLUDED
#define MADNESS_APPS_DIRAC_DFRESTART_H_INCLUDED

#include <madness/world/madness_exception.h>

namespace madness {

/// Format version for DFdriver restart archives.
///
/// Earlier archives carry no version field at all and are rejected. The field
/// is read back through the type cookie the BinaryFstream archives write, so a
/// legacy archive — whose first datum is a double total energy — fails the read
/// by construction.
inline constexpr unsigned int DF_RESTART_VERSION = 1;

/// Writes the current DF restart format version as the archive's first datum.
template <typename Archive>
void write_df_restart_version(const Archive& ar) {
    const unsigned int version = DF_RESTART_VERSION;
    ar & version;
}

/// Returns the archive's format version, or 0 for an unversioned archive.
///
/// Reads from a *local* archive, so the caller decides which rank reads.
template <typename Archive>
unsigned int read_df_restart_version(const Archive& ar) {
    unsigned int version = 0;
    try {
        ar & version;
    } catch (const MadnessException&) {
        // Type cookie mismatch: the first datum is a double, so this is a
        // pre-versioning archive. The archive prints its own message to
        // stderr before throwing; its what() points at a stack buffer and
        // must not be read here.
        return 0;
    }
    return version;
}

/// Throws unless \c version is a DF restart format this build can read.
inline void require_supported_df_restart(const unsigned int version) {
    MADNESS_CHECK_THROW(version != 0,
        "unversioned DF restart archive; start from a moldft archive or "
        "rerun with the older DFdriver");
    MADNESS_CHECK_THROW(version == DF_RESTART_VERSION,
        "unsupported DF restart archive version; start from a moldft archive "
        "or regenerate it with this DFdriver");
}

}  // namespace madness

#endif
