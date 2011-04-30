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

  $Id$
*/

#include <madness_config.h>

#ifdef MADNESS_HAS_GOOGLE_TEST

#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <world/worldreduce.h>
#include <world/functional.h>
#include <world/functional.h>
#include <world/deferred_deleter.h>
#include <gtest/gtest.h>

madness::World* pworld;

namespace {

    class Reducer : public madness::WorldReduce<Reducer, std::size_t> {
        typedef madness::WorldReduce<Reducer, std::size_t> WorldReducer_;
    public:

        Reducer(madness::World& w) :
            WorldReducer_(w)
        { process_pending(); }

        virtual ~Reducer() { }
    };

    class WorldReduceTest : public ::testing::Test {
    public:
        WorldReduceTest() : reducer(new Reducer(*pworld), madness::DeferredDeleter<Reducer>(*pworld)) {
            for(ProcessID r = 0; r < pworld->size(); ++r) {
                if((r % 2) == 0)
                    even.push_back(r);

                all.push_back(r);
            }
        }

        virtual ~WorldReduceTest() { }

        std::vector<ProcessID> all;
        std::vector<ProcessID> even;
        std::shared_ptr<Reducer> reducer;
    };

    TEST_F(WorldReduceTest, Construct) {
        EXPECT_NO_THROW(Reducer(*pworld));
    }

    TEST_F(WorldReduceTest, ReduceAll) {
        madness::Future<ProcessID> result = reducer->reduce(0, pworld->rank(),
                std::plus<ProcessID>(), all.begin(), all.end(), 0);

        if(pworld->rank() == 0) {
            ProcessID sum = 0;
            for(std::vector<ProcessID>::const_iterator it = all.begin(); it != all.end(); ++it)
                sum += *it;
            EXPECT_EQ(sum, result.get());
        }
    }

    TEST_F(WorldReduceTest, ReduceEven) {
        if((pworld->rank() % 2) == 0) {
            madness::Future<ProcessID> result = reducer->reduce(1, pworld->rank(),
                    std::plus<ProcessID>(), even.begin(), even.end(), 0);

            if(pworld->rank() == 0) {
                ProcessID sum = 0;
                for(std::vector<ProcessID>::const_iterator it = even.begin(); it != even.end(); ++it)
                    sum += *it;
                EXPECT_EQ(sum, result.get());
            }
        }
    }
}

int main(int argc, char **argv) {
    madness::initialize(argc,argv);
    int status = 0;
    {
        madness::World world(MPI::COMM_WORLD);
        pworld = &world;

        if (world.rank()) madness::redirectio(world);
        world.args(argc,argv);
        world.gop.fence();

        ::testing::InitGoogleTest(&argc, argv);
        status = RUN_ALL_TESTS();

        world.gop.fence();
    }
    madness::finalize();
    return status;
}


#else

#include <iostream>
int main() {
    std::cout << "U need to build with Google test to enable the world test code\n";
    return 0;
}

#endif
