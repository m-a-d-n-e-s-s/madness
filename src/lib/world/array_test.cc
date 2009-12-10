#include <world/array.h>
#include <madness_config.h>

#ifdef HAVE_GOOGLE_TEST

#include <gtest/gtest.h>

namespace {
    TEST(VectorTest, SizeWorks) {
        madness::Vector<double,3> v;
        EXPECT_EQ(3, v.size());
    }

    TEST(VectorTest, InitializerWorks) {
        madness::Vector<double,33> v(1.0);
        for (int i=0; i<v.size(); i++) {
            EXPECT_EQ(1.0, v[i]);
        }
    } 
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#else

#include <iostream>
int main() {
    std::cout << "U need to build with Google test to enable the array test code\n";
    return 0;
}

#endif
