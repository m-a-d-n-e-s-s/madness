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

#ifndef MADNESS_WORLD_CEREAL_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_CEREAL_ARCHIVE_H__INCLUDED

#if __has_include(<cereal/cereal.hpp>)

#ifndef MADNESS_HAS_CEREAL
#define MADNESS_HAS_CEREAL 1
#endif

#include <memory>
#include <cereal/cereal.hpp>
#include <cereal/details/traits.hpp>
#include <madness/world/archive.h>

namespace madness {
namespace archive {
/// Wraps an output archive around a Cereal archive
template <typename Muesli>
class CerealOutputArchive : public ::madness::archive::BaseOutputArchive {
  mutable std::shared_ptr<Muesli>
      muesli; ///< The cereal archive being wrapped, deleter determines whether this is an owning ptr

public:
  CerealOutputArchive(Muesli &muesli) : muesli(&muesli, [](Muesli *) {}) {}
  template <typename Arg, typename... RestOfArgs,
            typename = std::enable_if_t<
                !std::is_same<Muesli, std::decay_t<Arg>>::value>>
  CerealOutputArchive(Arg &&arg, RestOfArgs &&... rest_of_args)
      : muesli(new Muesli(std::forward<Arg>(arg),
                          std::forward<RestOfArgs>(rest_of_args)...),
               std::default_delete<Muesli>{}) {}

  template <class T, class Cereal = Muesli>
  inline std::enable_if_t<
      madness::is_trivially_serializable<T>::value &&
          !cereal::traits::is_text_archive<Cereal>::value,
      void>
  store(const T *t, long n) const {
    const unsigned char *ptr = (unsigned char *)t;
    (*muesli)(cereal::binary_data(ptr, sizeof(T) * n));
  }

  template <class T, class Cereal = Muesli>
  inline std::enable_if_t<
      !madness::is_trivially_serializable<T>::value ||
          cereal::traits::is_text_archive<Cereal>::value,void>
  store(const T *t, long n) const {
    for (long i = 0; i != n; ++i)
      *muesli & t[i];
  }

  void open(std::size_t hint) {}
  void close(){};
  void flush(){};
};

/// Wraps an input archive around a Cereal archive
template <typename Muesli> class CerealInputArchive : public BaseInputArchive {
  std::shared_ptr<Muesli> muesli; ///< The cereal archive being wrapped, deleter determines whether this is an owning ptr

public:
  CerealInputArchive(Muesli &muesli) : muesli(&muesli, [](Muesli *) {}) {}
  template <typename Arg, typename... RestOfArgs,
            typename = std::enable_if_t<
                !std::is_same<Muesli, std::decay_t<Arg>>::value>>
  CerealInputArchive(Arg &&arg, RestOfArgs &&... rest_of_args)
      : muesli(new Muesli(std::forward<Arg>(arg),
                          std::forward<RestOfArgs>(rest_of_args)...),
               std::default_delete<Muesli>{}) {}

  template <class T, class Cereal = Muesli>
  inline std::enable_if_t<
      madness::is_trivially_serializable<T>::value &&
          !cereal::traits::is_text_archive<Cereal>::value,
      void>
  load(T *t, long n) const {
    (*muesli)(cereal::binary_data(t, sizeof(T) * n));
  }

  template <class T, class Cereal = Muesli>
  inline std::enable_if_t<
      !madness::is_trivially_serializable<T>::value ||
          cereal::traits::is_text_archive<Cereal>::value,
      void>
  load(T *t, long n) const {
    for (long i = 0; i != n; ++i)
      *muesli & t[i];
  }

  void open(std::size_t hint) {}
  void rewind() const {}
  void close(){};
};
} // namespace archive

template <typename Muesli>
struct is_text_archive<
    archive::CerealInputArchive<Muesli>,
    std::enable_if_t<cereal::traits::is_text_archive<Muesli>::value>>
    : std::true_type {};
template <typename Muesli>
struct is_text_archive<
    archive::CerealOutputArchive<Muesli>,
    std::enable_if_t<cereal::traits::is_text_archive<Muesli>::value>>
    : std::true_type {};

template <typename Muesli, typename T>
struct is_serializable<
    archive::CerealOutputArchive<Muesli>, T,
    std::enable_if_t<(is_trivially_serializable<T>::value &&
        !cereal::traits::is_text_archive<Muesli>::value) ||
        (cereal::traits::detail::count_output_serializers<T, Muesli>::value != 0 &&
         cereal::traits::is_text_archive<Muesli>::value)>>
    : std::true_type {};
template <typename Muesli, typename T>
struct is_serializable<
    archive::CerealInputArchive<Muesli>, T,
    std::enable_if_t<
        (is_trivially_serializable<T>::value &&
            !cereal::traits::is_text_archive<Muesli>::value) ||
            (cereal::traits::detail::count_input_serializers<T, Muesli>::value != 0 &&
             cereal::traits::is_text_archive<Muesli>::value)>>
    : std::true_type {};

}  // namespace madness

#endif  // have cereal/cereal.hpp

#endif  // MADNESS_WORLD_CEREAL_ARCHIVE_H__INCLUDED