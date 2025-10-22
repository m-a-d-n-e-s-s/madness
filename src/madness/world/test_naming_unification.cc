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
*/

#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>

using namespace madness;

/// Test DistributionType from_string converter
void test_distribution_type_converter() {
    // Test Distributed
    MADNESS_CHECK(distribution_type_from_string("Distributed") == Distributed);
    MADNESS_CHECK(distribution_type_from_string("distributed") == Distributed);
    MADNESS_CHECK(distribution_type_from_string("DISTRIBUTED") == Distributed);
    
    // Test RankReplicated
    MADNESS_CHECK(distribution_type_from_string("RankReplicated") == RankReplicated);
    MADNESS_CHECK(distribution_type_from_string("rank_replicated") == RankReplicated);
    MADNESS_CHECK(distribution_type_from_string("rank-replicated") == RankReplicated);
    MADNESS_CHECK(distribution_type_from_string("rankreplicated") == RankReplicated);
    MADNESS_CHECK(distribution_type_from_string("rank") == RankReplicated);
    
    // Test NodeReplicated
    MADNESS_CHECK(distribution_type_from_string("NodeReplicated") == NodeReplicated);
    MADNESS_CHECK(distribution_type_from_string("node_replicated") == NodeReplicated);
    MADNESS_CHECK(distribution_type_from_string("node-replicated") == NodeReplicated);
    MADNESS_CHECK(distribution_type_from_string("nodereplicated") == NodeReplicated);
    MADNESS_CHECK(distribution_type_from_string("node") == NodeReplicated);
    MADNESS_CHECK(distribution_type_from_string("host_replicated") == NodeReplicated);
    MADNESS_CHECK(distribution_type_from_string("host-replicated") == NodeReplicated);
    MADNESS_CHECK(distribution_type_from_string("hostreplicated") == NodeReplicated);
    MADNESS_CHECK(distribution_type_from_string("host") == NodeReplicated);
    
    // Test DistributionTypeFromString wrapper
    DistributionTypeFromString wrapper1("distributed");
    MADNESS_CHECK(wrapper1.value == Distributed);
    MADNESS_CHECK(static_cast<DistributionType>(wrapper1) == Distributed);
    
    DistributionTypeFromString wrapper2("rank");
    MADNESS_CHECK(wrapper2.value == RankReplicated);
    
    DistributionTypeFromString wrapper3("node");
    MADNESS_CHECK(wrapper3.value == NodeReplicated);
    
    print("test_distribution_type_converter passed");
}

/// Test WorldContainer naming consistency
void test_worldcontainer_naming(World& world) {
    WorldContainer<int, int> c(world);
    
    // Test initial state
    MADNESS_CHECK(c.is_distributed());
    MADNESS_CHECK(!c.rank_replication());
    MADNESS_CHECK(!c.node_replication());
    
    // Test deprecated aliases still work
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    MADNESS_CHECK(!c.is_replicated());
    MADNESS_CHECK(!c.is_host_replicated());
    #pragma GCC diagnostic pop
    
    // Add some data
    if (world.rank() == 0) {
        for (int i = 0; i < 10; ++i) {
            c.replace(i, i * 10);
        }
    }
    world.gop.fence();
    
    // Test rank replication with new canonical name
    c.replicate_on_ranks();
    MADNESS_CHECK(c.rank_replication());
    MADNESS_CHECK(!c.node_replication());
    MADNESS_CHECK(!c.is_distributed());
    
    // Verify data is replicated to all ranks
    for (int i = 0; i < 10; ++i) {
        auto it = c.find(i);
        MADNESS_CHECK(it.get() != c.end());
        MADNESS_CHECK(it.get()->second == i * 10);
    }
    
    // Test deprecated replicate method still works
    WorldContainer<int, int> c2(world);
    if (world.rank() == 0) {
        for (int i = 0; i < 5; ++i) {
            c2.replace(i, i * 5);
        }
    }
    world.gop.fence();
    
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    c2.replicate();
    MADNESS_CHECK(c2.is_replicated());
    #pragma GCC diagnostic pop
    
    print("test_worldcontainer_naming passed on rank", world.rank());
}

int main(int argc, char** argv) {
    try {
        World& world = initialize(argc, argv);
        
        if (world.rank() == 0) {
            test_distribution_type_converter();
        }
        
        test_worldcontainer_naming(world);
        
        world.gop.fence();
        if (world.rank() == 0) {
            print("\nAll naming unification tests passed!");
        }
        
        finalize();
        return 0;
    }
    catch (const SafeMPI::Exception& e) {
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a string exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }
    
    return 1;
}
