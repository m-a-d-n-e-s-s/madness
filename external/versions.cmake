set(MADNESS_TRACKED_PARSEC_TAG 064918a398327371789999591c224aff3df8077b)
set(MADNESS_TRACKED_NLOHMANN_JSON_VERSION 3.12.0)
set(MADNESS_TRACKED_NLOHMANN_JSON_TAG v3.12.0)
set(MADNESS_TRACKED_LIBXSMM_VERSION 2.1.0)
set(MADNESS_TRACKED_LIBXSMM_TAG 2.1.0)
set(MADNESS_TRACKED_PCMSOLVER_VERSION 1.3.0)
set(MADNESS_TRACKED_PCMSOLVER_TAG v1.3.0)

# Boost headers for the PCMSolver source build, fetched only when the host has
# none. The b2-nodocs release tarball is the smallest upstream artifact that
# still carries the merged boost/ header tree (~50 MB, vs ~190 MB for the
# classic source tarball); PCMSolver needs boost/any.hpp and
# boost/numeric/odeint.hpp, and no compiled Boost libraries at all.
set(MADNESS_TRACKED_BOOST_VERSION 1.92.0)
set(MADNESS_TRACKED_BOOST_URL_HASH
    SHA256=ea7b982002cc9dfbe59b0b217b206f470dc75f3de0bb2973d844118934d82411)
