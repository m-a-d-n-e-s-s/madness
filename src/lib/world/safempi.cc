#include <world/safempi.h>

#ifdef SERIALIZE_MPI
madness::SCALABLE_MUTEX_TYPE SafeMPI::charon;
#endif

