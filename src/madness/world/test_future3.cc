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

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/MADworld.h>
#include <madness/world/world_object.h>

using namespace madness;

class Foo : public WorldObject<Foo> {
    std::vector<int> bar;
public:
    Foo(World& world, int bar) : WorldObject<Foo>(world), bar(1, bar) {
        process_pending();
    }

    virtual ~Foo() { }

    decltype(auto) get() const {
        return bar;
    }
};


int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    madness::World world(SafeMPI::COMM_WORLD);

    Foo a(world,world.rank()), b(world,world.rank()*10);

    for (ProcessID p=0; p<world.size(); ++p) {
        auto futa = a.send(p,&Foo::get);
        const auto futb = b.send(p,&Foo::get);
        // Could work here until the results are available
        MADNESS_CHECK(futa.get() == std::vector<int>{p});
        MADNESS_CHECK(futb.get() == std::vector<int>{p*10});

      // test various flavors of Future::get:
        {
          // 1) Future &::get()
          static_assert(
              std::is_same_v<decltype(futa.get()), std::vector<int> &>);
          std::vector<int> &vala_ref = futa.get();
          MADNESS_CHECK(vala_ref == std::vector<int>{p});
          // 2) Future const &::get()
          static_assert(
              std::is_same_v<decltype(futb.get()), std::vector<int> const &>);
          std::vector<int> const &valb_ref = futb.get();
          MADNESS_CHECK(valb_ref == std::vector<int>{p * 10});
          // 2) Future &&::get()
          static_assert(std::is_same_v<decltype(std::move(futa).get()),
                                       std::vector<int>>);
          std::vector<int> vala_copy = a.send(p, &Foo::get).get();
          MADNESS_CHECK(vala_copy == std::vector<int>{p});
        }

      // test various flavors of Future<T>::operator T:
        {
          // 1) Future &::operator T&()
          std::vector<int> &vala_ref = futa;
          MADNESS_CHECK(vala_ref == std::vector<int>{p});
          // 2) Future const &::operator T const&()
          std::vector<int> const &valb_ref = futb;
          MADNESS_CHECK(valb_ref == std::vector<int>{p * 10});
          // 2) Future &&::operator T()
          std::vector<int> vala_copy = a.send(p, &Foo::get);
          MADNESS_CHECK(vala_copy == std::vector<int>{p});
        }
    }
    world.gop.fence();
    if (world.rank() == 0) print("OK!");

    madness::finalize();
    return 0;
}
