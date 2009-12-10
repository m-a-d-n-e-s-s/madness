#include <tensor/tensor.h>
#include <world/print.h>

#ifdef HAVE_GOOGLE_TEST

using madness::_;
using madness::___;
using madness::_reverse;

#include <iostream>
#include <gtest/gtest.h>

namespace {
    template <typename T>
    class TensorTest : public ::testing::Test {
    public:
        TensorTest() {}

        virtual ~TensorTest() {}

        virtual void SetUp() {}

        virtual void TearDown() {}

        static T* move_test_ptr;        
        
        static madness::Tensor<T> f() { // For testing move semantics
            madness::Tensor<T> result(3);
            move_test_ptr = result.ptr();
            return result;
        }
    };

    template <typename T> T* TensorTest<T>::move_test_ptr;

    typedef ::testing::Types<int, long, float, double, float_complex, double_complex> TensorTestTypes;
    TYPED_TEST_CASE(TensorTest, TensorTestTypes);

    TYPED_TEST(TensorTest, Basic) {
        for (int ndim=0; ndim<=TENSOR_MAXDIM; ndim++) {
            try {
                std::vector<long> dim(ndim);
                for (int pass=0; pass<10; pass++) {
                    long nelem = 1;
                    for (int i=0; i<ndim; i++) {
                        dim[i] = (madness::RandomValue<long>()&0x5) + 1;
                        nelem *= dim[i];
                    }

                    // Use ASSERT not EXPECT otherwise we get a huge volume of output
                    // from failing tests.  Also, subsequent tests usually rely on 
                    // succcess of previous ones.
                    
                    madness::Tensor<TypeParam> empty; // Verify default constructor
                    ASSERT_EQ(empty.size(),0);
                    ASSERT_EQ(empty.ptr(),(TypeParam*)0);
                    ASSERT_EQ(empty.ndim(),-1);
                    ASSERT_EQ(empty.id(),madness::TensorTypeData<TypeParam>::id);
                    
                    madness::Tensor<TypeParam> d(dim);
                    ASSERT_EQ(d.id(),madness::TensorTypeData<TypeParam>::id);
                    ASSERT_EQ(d.ndim(),ndim);
                    ASSERT_EQ(d.size(),nelem);
                    ASSERT_NE(d.ptr(),(TypeParam*)0);
                    for (int i=0; i<ndim; i++) {
                        ASSERT_EQ(d.dim(i),dim[i]);
                    }
                    // Should be initialized to zero
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(0)));
                    
                    d.fillindex();   // Check fillindex and indexing
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(_index)));

                    d.fill(TypeParam(33)); // Check fill
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(33)));

                    d = TypeParam(21);     // Check fill
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(21)));
                 
                    madness::Tensor<TypeParam> q(d); // Copy constructor
                    ASSERT_EQ(q.id(),madness::TensorTypeData<TypeParam>::id);
                    ASSERT_EQ(q.ptr(),d.ptr()); // Was it shallow?
                    ASSERT_EQ(q.id(),madness::TensorTypeData<TypeParam>::id);
                    ASSERT_EQ(q.ndim(),ndim);
                    ASSERT_EQ(q.size(),nelem);
                    ITERATOR(d, ASSERT_EQ(q(IND),d(IND)));

                    q.clear();  // Check that clear is working
                    ASSERT_EQ(q.size(),0);
                    ASSERT_EQ(q.ptr(),(TypeParam*)0);
                    ASSERT_EQ(q.ndim(),-1);
                    ASSERT_EQ(q.id(),madness::TensorTypeData<TypeParam>::id);

                    q = d;    // Assignment
                    ASSERT_EQ(q.ptr(),d.ptr()); // Was it shallow?
                    ITERATOR(d, ASSERT_EQ(q(IND),d(IND)));

                    q = copy(d); // Deep copy
                    ASSERT_NE(q.ptr(),d.ptr()); // Was it deep?
                    ITERATOR(d, ASSERT_EQ(q(IND),d(IND)));

                    madness::Tensor<TypeParam> r = TensorTest<TypeParam>::f(); // Should invoke move semantics
                    ASSERT_EQ(r.ptr(),TensorTest<TypeParam>::move_test_ptr);

                    r.clear();  // Does clear work?
                    ASSERT_EQ(r.size(),0);
                    ASSERT_EQ(r.ptr(),(TypeParam*)0);
                    ASSERT_EQ(r.ndim(),-1);
                    ASSERT_EQ(r.id(),madness::TensorTypeData<TypeParam>::id);

                    r = TensorTest<TypeParam>::f(); // Should invoke move semantics
                    ASSERT_EQ(r.ptr(),TensorTest<TypeParam>::move_test_ptr);

                    q.fillindex();

                    d.fillindex();
                    
                    q *= TypeParam(2); // Check inplace scaling
                    d.scale(TypeParam(3));
                    ITERATOR(q, ASSERT_EQ(q(IND),TypeParam(_index*2)));
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(_index*3)));
                    
                    d += TypeParam(2); // Check inplace scalar addition
                    ITERATOR(q, ASSERT_EQ(q(IND),TypeParam(_index*2)));
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(_index*3+2)));
                    
                    d -= q;     // Check inplace tensor subtraction
                    ITERATOR(q, ASSERT_EQ(q(IND),TypeParam(_index*2)));
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(_index+2)));

                    d += q;     // Check inplace tensor addition
                    ITERATOR(q, ASSERT_EQ(q(IND),TypeParam(_index*2)));
                    ITERATOR(d, ASSERT_EQ(d(IND),TypeParam(_index*3+2)));
                }
            }
            catch (const madness::TensorException& e) {
                if (ndim != 0) std::cout << e;
                EXPECT_EQ(ndim,0);
            }
            catch(...) {
                std::cout << "Caught unknown exception" << std::endl;
                EXPECT_EQ(1,0);
            }
        }
    }

    TYPED_TEST(TensorTest, Slicing1d) {
        madness::Tensor<TypeParam> t(5L);
        t(_) = 3;
        ITERATOR(t,EXPECT_EQ(t(IND),TypeParam(3)));
    }

    TYPED_TEST(TensorTest, Transpose) {
        for (long n=1; n<=21; n++) {
            madness::Tensor<TypeParam> a(n,n+3);
            ASSERT_TRUE(a.iscontiguous());
            ASSERT_EQ(a.size(),n*(n+3));
            a.fillrandom();
            ITERATOR(a, ASSERT_NE(a(IND),TypeParam(0)));

            madness::Tensor<TypeParam> aT = transpose(a);
            ASSERT_TRUE(aT.iscontiguous());
            ASSERT_EQ(aT.size(),n*(n+3));
            ASSERT_NE(aT.ptr(),a.ptr()); // Was it deep?
            ITERATOR2(aT, ASSERT_EQ(aT(_i,_j),a(_j,_i)));

            a.fillrandom();
            aT.clear();
            aT = transpose(a);
            ASSERT_EQ(aT.size(),n*(n+3));
            ASSERT_TRUE(a.iscontiguous());
            ASSERT_TRUE(aT.iscontiguous());
            ASSERT_NE(aT.ptr(),a.ptr()); // Was it deep?
            ITERATOR2(aT, ASSERT_EQ(aT(_i,_j),a(_j,_i)));

            aT = a.swapdim(0,1);
            ASSERT_EQ(aT.size(),n*(n+3));
            ASSERT_TRUE(a.iscontiguous());
            ASSERT_FALSE(aT.iscontiguous());
            ASSERT_EQ(aT.ptr(),a.ptr()); // Was it shallow?
            ITERATOR2(aT, ASSERT_EQ(aT(_i,_j),a(_j,_i)));
        }
    }            

//     TYPED_TEST(TensorTest, Container) {
//         typedef madness::ConcurrentHashMap< int, Tensor<TypeParam> > containerT;
//         static const int N = 100;
//         containerT c;
//         Tensor<TypeParam> a[N];
//         for (int i=0; i<N; i++) {
            
//         }
//     }

}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#else

#include <iostream>
int main() {
    std::cout << "U need to build with Google test to enable the tensor test code\n";
    return 0;
}

#endif
