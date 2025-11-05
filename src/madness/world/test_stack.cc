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

#include <madness/madness_config.h>
#ifdef MADNESS_HAS_GOOGLE_TEST

#define MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE 0
#include <madness/world/stack.h>
#include <gtest/gtest.h>


namespace {

    using namespace madness;


    class C {
    public:
        int value;

        C(int v) : value(v) { }
    };


    template <typename T>
    struct StackTestBase {
      static T make(int i) { return i; }
      static int value(const T& i) { return i; }
    };

    template <>
    struct StackTestBase<int*> {

        union type_converter {
            int* p;
            int v;
        };

        static int* make(int i) {
            type_converter tc;
            tc.p = nullptr;
            tc.v = i;
            return tc.p;
        }

        static int value(int* i) {
            type_converter tc;
            tc.p = i;
            return tc.v;
        }
    };

    template <>
    struct StackTestBase<std::shared_ptr<int> > {
        static std::shared_ptr<int> make(int i)
        { return std::shared_ptr<int>(new int(i)); }

        static int value(const std::shared_ptr<int>& i) { return *i; }
    };

    template <>
    struct StackTestBase<C> {
      static C make(int i) { return C(i); }
      static int value(const C& i) { return i.value; }
    };

    template <typename T>
    class StackTest : public ::testing::Test, public StackTestBase<T> {
    public:
        typedef Stack<T,4> Stack_;

        StackTest() { }

        virtual ~StackTest() { }

        virtual void SetUp() { }

        virtual void TearDown() { }

    };

    // Set types use in tests
    typedef ::testing::Types<int, int*, std::shared_ptr<int>, C> MyTypes;
    TYPED_TEST_CASE(StackTest, MyTypes);

    TYPED_TEST(StackTest, DefaultConstructor) {
        typedef typename StackTest<gtest_TypeParam_>::Stack_ Stack_;
        Stack_ s;

        EXPECT_EQ(0u, s.size());
        EXPECT_EQ(4u, s.capacity());
        EXPECT_TRUE(s.empty());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(s.top(), madness::MadnessException);
        EXPECT_THROW(s.pop(), madness::MadnessException);
#endif
    }

    TYPED_TEST(StackTest, PushPop) {
        // Add names from the test fixture to this scope
        typedef typename StackTest<gtest_TypeParam_>::Stack_ Stack_;
        const auto& make = StackTest<gtest_TypeParam_>::make;
        const auto& value = StackTest<gtest_TypeParam_>::value;
        Stack_ s;

        unsigned int i = 0u;
        unsigned int cap = s.capacity();
        const typename Stack_::value_type* const small_buffer = s.data();

        // Test pushing data up to the small buffer size
        for(; i < cap; ++i) {
            EXPECT_EQ(i, s.size());
            EXPECT_NO_THROW(s.push(make(i)));
            EXPECT_EQ(i, value(s.top()));
            EXPECT_EQ(i + 1u, s.size());
            EXPECT_EQ(cap, s.capacity());
            EXPECT_FALSE(s.empty());
        }

        // Test popping data to empty
        for(; i > 0; --i) {
            EXPECT_FALSE(s.empty());
            EXPECT_EQ(i, s.size());
            EXPECT_EQ(i - 1, value(s.top()));
            EXPECT_NO_THROW(s.pop());
            EXPECT_EQ(i - 1, s.size());
            EXPECT_EQ(cap, s.capacity());
        }

        // Check that all elements have been popped
        EXPECT_EQ(small_buffer, s.data());
        EXPECT_EQ(0u, s.size());
        EXPECT_TRUE(s.empty());

        // Test pushing data with large buffer size
        for(unsigned int n = 0u; n < 3; ++n) {
            for(; i < cap; ++i) {
                EXPECT_EQ(i, s.size());
                EXPECT_NO_THROW(s.push(make(i)));
                EXPECT_EQ(i, value(s.top()));
                EXPECT_EQ(i + 1u, s.size());
                EXPECT_EQ(cap, s.capacity());
                EXPECT_FALSE(s.empty());
            }

            // Test that the buffer is reallocated and increases in size when
            // the stack is at capacity
            typename Stack_::value_type* const buffer = s.data();
            EXPECT_NO_THROW(s.push(make(i)));
            EXPECT_LT(cap, s.capacity());
            EXPECT_NE(buffer, s.data());
            cap = s.capacity();
            ++i;
        }

        // Test popping data to empty
        for(; i > 0; --i) {
            EXPECT_FALSE(s.empty());
            EXPECT_EQ(i, s.size());
            EXPECT_EQ(i - 1, value(s.top()));
            EXPECT_NO_THROW(s.pop());
            EXPECT_EQ(i - 1, s.size());
            EXPECT_EQ(cap, s.capacity());
        }

        // Check that all elements have been popped
        EXPECT_NE(small_buffer, s.data());
        EXPECT_EQ(0u, s.size());
        EXPECT_TRUE(s.empty());
    }


    TYPED_TEST(StackTest, CopyCtor) {
        // Add names from the test fixture to this scope
        typedef typename StackTest<gtest_TypeParam_>::Stack_ Stack_;
        const auto& make = StackTest<gtest_TypeParam_>::make;
        const auto& value = StackTest<gtest_TypeParam_>::value;
        Stack_ s;

        unsigned int i = 0u;
        unsigned int cap = s.capacity();

        for(; i < cap * 0.5; ++i) {
            s.push(make(i));
        }

        // Test moving a small stack

        Stack_ css(s);

        // Check that the size and capacity has been correctly copied
        EXPECT_EQ(s.size(), css.size());
        EXPECT_EQ(s.capacity(), css.capacity());
        EXPECT_NE(s.data(), css.data());
        EXPECT_FALSE(s.empty());
        EXPECT_FALSE(css.empty());

        // Test that the move target stack contains the data that was held by s.
        for(; i > 0u; --i) {
            EXPECT_EQ(i, css.size());
            EXPECT_EQ(value(s.top()), value(css.top()));
            EXPECT_NO_THROW(css.pop());
            EXPECT_NO_THROW(s.pop());
            EXPECT_EQ(i - 1, css.size());
            EXPECT_EQ(cap, css.capacity());
        }

        for(; i < cap * 1.5; ++i) {
            EXPECT_NO_THROW(s.push(make(i)));
        }
        cap = s.size();

        // Test moving a large stack

        Stack_ cls(s);

        // Check that the size and capacity has been correctly copied
        EXPECT_EQ(s.size(), cls.size());
        EXPECT_EQ(s.size(), cls.capacity());
        EXPECT_NE(s.data(), cls.data());
        EXPECT_FALSE(s.empty());
        EXPECT_FALSE(cls.empty());

        // Test that the move target stack contains the data that was held by s.
        for(; i > 0; --i) {
            EXPECT_EQ(i, cls.size());
            EXPECT_EQ(value(s.top()), value(cls.top()));
            EXPECT_NO_THROW(cls.pop());
            EXPECT_NO_THROW(s.pop());
            EXPECT_EQ(i - 1, cls.size());
            EXPECT_EQ(cap, cls.capacity());
        }

    }

    TYPED_TEST(StackTest, MoveCtor) {
        // Add names from the test fixture to this scope
        typedef typename StackTest<gtest_TypeParam_>::Stack_ Stack_;
        const auto& make = StackTest<gtest_TypeParam_>::make;
        const auto& value = StackTest<gtest_TypeParam_>::value;
        Stack_ s;

        unsigned int i = 0u;
        unsigned int cap = s.capacity();
        const typename Stack_::value_type* const small_buffer = s.data();

        for(; i < cap * 0.5; ++i) {
            s.push(make(i));
        }

        // Test moving a small stack

        Stack_ mss(std::move(s));

        // Check that the size and capacity has been correctly moved
        EXPECT_EQ(small_buffer, s.data());
        EXPECT_EQ(0u, s.size());
        EXPECT_EQ(cap, s.capacity());
        EXPECT_NE(small_buffer, mss.data());
        EXPECT_EQ(i, mss.size());
        EXPECT_EQ(cap, mss.capacity());
        EXPECT_TRUE(s.empty());
        EXPECT_FALSE(mss.empty());

        // Test that the move target stack contains the data that was held by s.
        for(; i > 0u; --i) {
            EXPECT_EQ(i, mss.size());
            EXPECT_EQ(int(i - 1), value(mss.top()));
            EXPECT_NO_THROW(mss.pop());
            EXPECT_EQ(i - 1, mss.size());
            EXPECT_EQ(cap, mss.capacity());
        }

        for(; i < cap * 1.5; ++i) {
            EXPECT_NO_THROW(s.push(make(i)));
        }
        cap = s.capacity();
        typename Stack_::value_type* buffer = s.data();

        // Test moving a large stack

        Stack_ mls(std::move(s));

        // Check that the size and capacity has been correctly moved
        EXPECT_EQ(small_buffer, s.data());
        EXPECT_EQ(0u, s.size());
        EXPECT_EQ(size_t(4), s.capacity());
        EXPECT_EQ(buffer, mls.data());
        EXPECT_EQ(i, mls.size());
        EXPECT_EQ(cap, mls.capacity());
        EXPECT_TRUE(s.empty());
        EXPECT_FALSE(mls.empty());

        // Test that the move target stack contains the data that was held by s.
        for(; i > 0; --i) {
            EXPECT_EQ(i, mls.size());
            EXPECT_EQ(i - 1, value(mls.top()));
            EXPECT_NO_THROW(mls.pop());
            EXPECT_EQ(i - 1, mls.size());
            EXPECT_EQ(cap, mls.capacity());
        }

    }

    TYPED_TEST(StackTest, CopyAssign) {
        // Add names from the test fixture to this scope
        typedef typename StackTest<gtest_TypeParam_>::Stack_ Stack_;
        const auto& make = StackTest<gtest_TypeParam_>::make;
        const auto& value = StackTest<gtest_TypeParam_>::value;
        Stack_ s;
        const unsigned int small_cap = s.capacity();

        for(double fill = 0.0; fill <= 3.0; fill += 0.75) {

            unsigned int i = 0u;
            for(; i < small_cap * 0.5; ++i) {
                s.push(make(i));
            }

            // Test moving a small stack

            Stack_ css;
            for(unsigned int x = 0u; x < small_cap * fill; ++x) {
                css.push(make(x + 20));
            }
            unsigned int cap = css.capacity();
            if(cap < s.size())
                cap = s.size();
            css = s;

            // Check that the size and capacity has been correctly copied
            EXPECT_EQ(s.size(), css.size());
            EXPECT_EQ(cap, css.capacity());
            EXPECT_NE(s.data(), css.data());
            EXPECT_FALSE(s.empty());
            EXPECT_FALSE(css.empty());

            // Test that the move target stack contains the data that was held by s.
            for(; i > 0u; --i) {
                EXPECT_EQ(i, css.size());
                EXPECT_EQ(value(s.top()), value(css.top()));
                EXPECT_NO_THROW(css.pop());
                EXPECT_NO_THROW(s.pop());
                EXPECT_EQ(i - 1, css.size());
                EXPECT_EQ(cap, css.capacity());
            }

            for(; i < small_cap * 1.5; ++i) {
                EXPECT_NO_THROW(s.push(make(i)));
            }

            // Test moving a large stack

            Stack_ cls;
            for(unsigned int x = 0u; x < small_cap * fill; ++x) {
                cls.push(make(x + 20));
            }
            cap = cls.capacity();
            if(cap < s.size())
                cap = s.size();
            cls = s;

            // Check that the size and capacity has been correctly copied
            EXPECT_EQ(s.size(), cls.size());
            EXPECT_EQ(cap, cls.capacity());
            EXPECT_NE(s.data(), cls.data());
            EXPECT_FALSE(s.empty());
            EXPECT_FALSE(cls.empty());

            // Test that the move target stack contains the data that was held by s.
            for(; i > 0; --i) {
                EXPECT_EQ(i, cls.size());
                EXPECT_EQ(value(s.top()), value(cls.top()));
                EXPECT_NO_THROW(cls.pop());
                EXPECT_NO_THROW(s.pop());
                EXPECT_EQ(i - 1, cls.size());
                EXPECT_EQ(cap, cls.capacity());
            }
        }

    }

    TYPED_TEST(StackTest, MoveAssign) {
        // Add names from the test fixture to this scope
        typedef typename StackTest<gtest_TypeParam_>::Stack_ Stack_;
        const auto& make = StackTest<gtest_TypeParam_>::make;
        const auto& value = StackTest<gtest_TypeParam_>::value;
        Stack_ s;
        const unsigned int small_cap = s.capacity();

        for(double fill = 0.0; fill <= 3.0; fill += 0.75) {

            unsigned int i = 0u;
            unsigned int cap = s.capacity();
            const typename Stack_::value_type* const small_buffer = s.data();

            for(; i < small_cap * 0.5; ++i) {
                s.push(make(i));
            }

            // Test moving a small stack

            Stack_ mss;
            for(unsigned int x = 0u; x < small_cap * fill; ++x) {
                mss.push(make(x + 20));
            }
            mss = std::move(s);

            // Check that the size and capacity has been correctly moved
            EXPECT_EQ(small_buffer, s.data());
            EXPECT_EQ(0u, s.size());
            EXPECT_EQ(cap, s.capacity());
            EXPECT_NE(small_buffer, mss.data());
            EXPECT_EQ(i, mss.size());
            EXPECT_EQ(cap, mss.capacity());
            EXPECT_TRUE(s.empty());
            EXPECT_FALSE(mss.empty());

            // Test that the move target stack contains the data that was held by s.
            for(; i > 0u; --i) {
                EXPECT_EQ(i, mss.size());
                EXPECT_EQ(i - 1, value(mss.top()));
                EXPECT_NO_THROW(mss.pop());
                EXPECT_EQ(i - 1, mss.size());
                EXPECT_EQ(cap, mss.capacity());
            }

            for(; i < small_cap * 1.5; ++i) {
                EXPECT_NO_THROW(s.push(make(i)));
            }
            cap = s.capacity();
            typename Stack_::value_type* buffer = s.data();

            // Test moving a large stack

            Stack_ mls;
            for(unsigned int x = 0u; x < small_cap * fill; ++x) {
                mls.push(make(x + 20));
            }
            mls = std::move(s);

            // Check that the size and capacity has been correctly moved
            EXPECT_EQ(small_buffer, s.data());
            EXPECT_EQ(0u, s.size());
            EXPECT_EQ(small_cap, s.capacity());
            EXPECT_EQ(buffer, mls.data());
            EXPECT_EQ(i, mls.size());
            EXPECT_EQ(cap, mls.capacity());
            EXPECT_TRUE(s.empty());
            EXPECT_FALSE(mls.empty());

            // Test that the move target stack contains the data that was held by s.
            for(; i > 0; --i) {
                EXPECT_EQ(i, mls.size());
                EXPECT_EQ(i - 1, value(mls.top()));
                EXPECT_NO_THROW(mls.pop());
                EXPECT_EQ(i - 1, mls.size());
                EXPECT_EQ(cap, mls.capacity());
            }
        }
    }

} // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int status = RUN_ALL_TESTS();

    return status;
}


#else

#include <iostream>
int main() {
    std::cout << "!!! Error: You need to build with Google test to enable WorldRef test code\n";
    return 1;
}

#endif
