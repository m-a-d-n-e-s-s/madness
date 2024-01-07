/*
  This file is part of MADNESS.

  Copyright (C) 2013  Virginia Tech

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

#include <madness/world/dist_cache.h>
#include <madness/world/worldgop.h>
#include <iostream>

namespace madness {
namespace detail {

std::unique_ptr<std::array<std::pair<std::atomic<std::size_t>, std::atomic<std::size_t>>, 1000000>> dist_caches_stats;

struct DistCacheStatsInitializer {
  DistCacheStatsInitializer() {
    dist_caches_stats = std::make_unique<std::array<std::pair<std::atomic<std::size_t>, std::atomic<std::size_t>>, 1000000>>();
    for (auto&& [nreads, nwrites] : *dist_caches_stats) {
      nreads.store(0);
      nwrites.store(0);
    }
  }
};
DistCacheStatsInitializer dist_cache_stats_initializer;

void dump_dist_caches_and_stats(std::ostream &os) {
  os << "DistCache stats:\n";
  std::size_t obj_id = 0;
  for (auto&& [nreads, nwrites] : *dist_caches_stats) {
    if (nreads || nwrites)
      os << "  " << obj_id << ": " << nreads << " reads != " << nwrites
         << " writes\n";
    ++obj_id;
  }
  DistCache<ProcessKey<DistributedID, WorldGopInterface::PointToPointTag>>::dump(os);
}

void stdout_dump_dist_caches_and_stats() {
  dump_dist_caches_and_stats(std::cout);
}

} // namespace detail
} // namespace madness
