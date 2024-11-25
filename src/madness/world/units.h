//
// Created by Eduard Valeyev on 2/5/24.
//

#ifndef MADNESS_UNITS_H
#define MADNESS_UNITS_H

#include <cstddef>
#include <cstdint>

namespace madness {

namespace units::literals {

constexpr std::size_t operator""_KiB(unsigned long long int x) {
  return 1024ULL * x;
}

constexpr std::size_t operator""_MiB(unsigned long long int x) {
  return 1024_KiB * x;
}

constexpr std::size_t operator""_GiB(unsigned long long int x) {
  return 1024_MiB * x;
}

constexpr std::size_t operator""_TiB(unsigned long long int x) {
  return 1024_GiB * x;
}

constexpr std::size_t operator""_PiB(unsigned long long int x) {
  return 1024_TiB * x;
}

constexpr std::size_t operator""_EiB(unsigned long long int x) {
  return 1024_PiB * x;
}

constexpr std::size_t operator""_kB(unsigned long long int x) {
  return 1024ULL * x;
}

constexpr std::size_t operator""_KB(unsigned long long int x) {
  return 1024ULL * x;
}

constexpr std::size_t operator""_MB(unsigned long long int x) {
  return 1024_kB * x;
}

constexpr std::size_t operator""_GB(unsigned long long int x) {
  return 1024_MB * x;
}

constexpr std::size_t operator""_TB(unsigned long long int x) {
  return 1024_GB * x;
}

constexpr std::size_t operator""_PB(unsigned long long int x) {
  return 1024_TB * x;
}

constexpr std::size_t operator""_EB(unsigned long long int x) {
  return 1024_PB * x;
}

} // namespace units::literals

/// Unit-aware conversion of a C string to a size_t

/// See https://en.wikipedia.org/wiki/Kilobyte for units. This assumes
/// memory units throughout, i.e. `MiB` and `MB` both mean 1024*1024 bytes, etc.
/// To reduce confusion, `KiB`, `kB`, and `KB` all mean 1024 bytes.
/// @param str The C string to convert.
/// @return The memory size, in bytes
std::uint64_t cstr_to_memory_size(const char* str);

} // namespace madness

#endif // MADNESS_UNITS_H
