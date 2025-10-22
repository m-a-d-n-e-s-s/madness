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

#include <madness.h>
#include <madness/chem/SCFOperators.h>

using namespace madness;

/// Test Exchange::Algorithm from_string converter
void test_algorithm_converter() {
    using AlgType = Exchange<double, 3>::Algorithm;
    
    // Test small_memory
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("small_memory") == AlgType::small_memory);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("smallmemory") == AlgType::small_memory);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("small") == AlgType::small_memory);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("Small_Memory") == AlgType::small_memory);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("SMALL-MEMORY") == AlgType::small_memory);
    
    // Test large_memory
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("large_memory") == AlgType::large_memory);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("largememory") == AlgType::large_memory);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("large") == AlgType::large_memory);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("Large Memory") == AlgType::large_memory);
    
    // Test multiworld_efficient
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("multiworld_efficient") == AlgType::multiworld_efficient);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("multiworldefficient") == AlgType::multiworld_efficient);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("multiworld") == AlgType::multiworld_efficient);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("MultiWorld-Efficient") == AlgType::multiworld_efficient);
    
    // Test multiworld_efficient_row
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("multiworld_efficient_row") == AlgType::multiworld_efficient_row);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("multiworld_efficientrow") == AlgType::multiworld_efficient_row);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("MultiWorld Efficient Row") == AlgType::multiworld_efficient_row);
    
    // Test fetch_compute
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("fetch_compute") == AlgType::fetch_compute);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("fetchcompute") == AlgType::fetch_compute);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("fetch") == AlgType::fetch_compute);
    MADNESS_CHECK(Exchange<double, 3>::from_string_algorithm("Fetch-Compute") == AlgType::fetch_compute);
    
    // Test deprecated from_string wrapper
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    MADNESS_CHECK(Exchange<double, 3>::from_string("small") == AlgType::small_memory);
    MADNESS_CHECK(Exchange<double, 3>::from_string("large") == AlgType::large_memory);
    #pragma GCC diagnostic pop
    
    // Test free function
    MADNESS_CHECK(algorithm_from_string<double, 3>("small") == AlgType::small_memory);
    MADNESS_CHECK(algorithm_from_string<double, 3>("large") == AlgType::large_memory);
    
    // Test AlgorithmFromString wrapper
    AlgorithmFromString<double, 3> wrapper1("small");
    MADNESS_CHECK(wrapper1.value == AlgType::small_memory);
    MADNESS_CHECK(static_cast<AlgType>(wrapper1) == AlgType::small_memory);
    
    AlgorithmFromString<double, 3> wrapper2("large");
    MADNESS_CHECK(wrapper2.value == AlgType::large_memory);
    
    // Test user-defined literal
    MADNESS_CHECK("small"_alg == AlgType::small_memory);
    MADNESS_CHECK("large"_alg == AlgType::large_memory);
    MADNESS_CHECK("multiworld"_alg == AlgType::multiworld_efficient);
    
    // Test invalid string throws exception
    bool caught = false;
    try {
        Exchange<double, 3>::from_string_algorithm("invalid_algorithm");
    } catch (const std::invalid_argument& e) {
        caught = true;
    }
    MADNESS_CHECK(caught);
    
    print("test_algorithm_converter passed");
}

int main(int argc, char** argv) {
    try {
        World& world = initialize(argc, argv);
        
        if (world.rank() == 0) {
            test_algorithm_converter();
            print("\nAll Exchange::Algorithm converter tests passed!");
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
