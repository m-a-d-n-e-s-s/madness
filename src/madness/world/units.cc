//
// Created by Eduard Valeyev on 2/5/24.
//

#include <madness/world/units.h>

#include <madness/world/madness_exception.h>


#include <sstream>

namespace madness {

std::uint64_t cstr_to_memory_size(const char *str) {
  // Convert the string into bytes
  std::stringstream ss(str);
  std::uint64_t memory = 0;
  if (ss >> memory) {
    if (memory > 0) {
      std::string unit;
      if (ss >> unit) { // Failure == assume bytes
        if (unit == "KB" || unit == "kB" || unit == "KiB") {
          MADNESS_ASSERT(memory <= (1ul << 54));
          memory *= 1024ul;
        } else if (unit == "MB" || unit == "MiB") {
          MADNESS_ASSERT(memory <= (1ul << 44));
          memory *= 1024ul * 1024;
        } else if (unit == "GB" || unit == "GiB") {
          MADNESS_ASSERT(memory <= (1ul << 34));
          memory *= 1024ul * 1024 * 1024;
        } else if (unit == "TB" || unit == "TiB") {
          MADNESS_ASSERT(memory <= (1ul << 24));
          memory *= 1024ul * 1024 * 1024 * 1024;
        } else if (unit == "PB" || unit == "PiB") {
          MADNESS_ASSERT(memory <= (1ul << 14));
          memory *= 1024ul * 1024 * 1024 * 1024 * 1024;
        } else if (unit == "EB" || unit == "EiB") {
          MADNESS_ASSERT(memory <= (1ul << 4));
          memory *= 1024ul * 1024 * 1024 * 1024 * 1024 * 1024;
        }
      }
    }
  }

  return memory;
}

}  // namespace madness