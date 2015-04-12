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

/// \file tensor/test.cc
/// \brief New test code for Tensor class using Google unit test

#include <madness/tensor/tensor.h>
#include <madness/tensor/gentensor.h>
#include <madness/world/print.h>


#ifdef MADNESS_HAS_GOOGLE_TEST

// The test code deliberately uses only the dumb ITERATOR macros
// in order to test the optimized iterators used by the implementation.

using namespace madness;
using madness::_;
using madness::___;
using madness::_reverse;

#include <iostream>
#include <gtest/gtest.h>
#include <tr1/tuple>

namespace {

    template <typename T>
    class GenTensorTest : public ::testing::Test {
    public:
        GenTensorTest() {}

        virtual ~GenTensorTest() {}

        virtual void SetUp() {}

        virtual void TearDown() {}
    };

    double eps=1.e-4;

    /// will yield different ranks (null: 0, index: 2, random: full)
    enum TensorFill {null, index, random};

    /// prepare a full rank tensor according to target rank
    template<typename T>
    Tensor<T> prep_tensor(const std::vector<long>& dim, const TensorFill& f) {
    	Tensor<T> t0=Tensor<double>(dim);
		if (f==index) t0.fillindex();
		if (f==random) t0.fillrandom();
		if (f!=null) {
			t0+=1;
			t0.scale(1.0/t0.normf());
		}
		return t0;
    }


    class UnaryGenTest : public ::testing::TestWithParam<std::tr1::tuple<
    		::madness::TensorType, long, double, TensorFill> > {
    public:
    	std::vector<long> dim;
        TensorType tt;
        double eps;

        TensorFill fill0;

        Tensor<double> t0;

		GenTensor<double> g0;

		UnaryGenTest() {
			madness::Random::Random(100);
            int alldim= (madness::RandomValue<long>()&0x5) + 3;

            tt=std::tr1::get<0>(this->GetParam());
            long ndim=std::tr1::get<1>(this->GetParam());
            eps=std::tr1::get<2>(this->GetParam());

            fill0=std::tr1::get<3>(this->GetParam());
            dim=std::vector<long>(ndim,alldim);

    		// fill the full tensors t0 and t1 with numbers that will yield certain ranks
    		t0=prep_tensor<double>(dim,fill0);

    		g0=GenTensor<double>(t0,TensorArgs(eps,tt));
    	}
    };

    class BinaryGenTest : public ::testing::TestWithParam<std::tr1::tuple<
    		::madness::TensorType, long, double, TensorFill, TensorFill> > {
    public:
    	std::vector<long> dim;
        TensorType tt;
        double eps;

        TensorFill fill0;
        TensorFill fill1;

        Tensor<double> t0;
		Tensor<double> t1;

		GenTensor<double> g0;
		GenTensor<double> g1;


        BinaryGenTest() {
            int alldim= (madness::RandomValue<long>()&0x5) + 3;

            tt=std::tr1::get<0>(this->GetParam());
            long ndim=std::tr1::get<1>(this->GetParam());
            eps=std::tr1::get<2>(this->GetParam());

            fill0=std::tr1::get<3>(this->GetParam());
            fill1=std::tr1::get<4>(this->GetParam());
            dim=std::vector<long>(ndim,alldim);

    		// fill the full tensors t0 and t1 with numbers that will yield certain ranks
    		t0=prep_tensor<double>(dim,fill0);
    		t1=prep_tensor<double>(dim,fill1);

    		g0=GenTensor<double>(t0,TensorArgs(eps,tt));
    		g1=GenTensor<double>(t1,TensorArgs(eps,tt));


    	}
    };

    INSTANTIATE_TEST_CASE_P(UnaryGenTest,UnaryGenTest,
    				::testing::Combine(::testing::Values(::madness::TT_FULL, ::madness::TT_2D),
    								   ::testing::Values(2,4,6),
    								   ::testing::Values(1.e-3, 1.e-4, 1.e-5),
     								   ::testing::Values(null,index,random)
    ));

    INSTANTIATE_TEST_CASE_P(BinaryGenTest,BinaryGenTest,
    				::testing::Combine(::testing::Values(::madness::TT_FULL, ::madness::TT_2D),
    								   ::testing::Values(2l,4l),
    								   ::testing::Values(1.e-3, 1.e-4, 1.e-5),
     								   ::testing::Values(null,index,random),
      								   ::testing::Values(null,index,random)
    ));

    TEST_P(UnaryGenTest, Norms_etc) {
    	try {

    		// should not change contents of the GenTensor
    		g0.normalize();
    		ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

    		// Frobenius norm
    		ASSERT_LT(g0.normf()-t0.normf(),eps);
    		ASSERT_LT(g0.svd_normf()-t0.normf(),eps);

    		// svd_normf computes the norm of the weights only
    		g0.scale(1.1);
    		ASSERT_LT(g0.svd_normf()-1.1*(t0.normf()),eps);

    	} catch (const madness::TensorException& e) {
    		if (dim.size() != 0) std::cout << e;
    		EXPECT_EQ(dim.size(),0);
    	} catch(...) {
    		std::cout << "Caught unknown exception" << std::endl;
    		EXPECT_EQ(1,0);
    	}
    }

    TEST_P(BinaryGenTest, Addition) {
    	for (int ipass=0; ipass<10; ++ipass) {
    	try {
    		// check for addition
    		t0+=t1;
    		g0+=g1;
    		ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

    		// check for subtraction
    		t0-=t1;
    		g0-=g1;
    		ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

    		// check for inplace stuff
    		t0.gaxpy(-0.7, t1, 0.1);
    		g0.gaxpy(-0.7, g1, 0.1);
    		ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

    		/*
    		 *  this is somewhat special, but nevertheless important to test
    		 *
    		 *  also: note the fac_reduce which is important for accuracy.
    		 *  Many subsequent operations will deteriorate the error.
    		 */
    		g0=GenTensor<double>(t0,TensorArgs(eps,tt));
    		g1=GenTensor<double>(t1,TensorArgs(eps,tt));

    		g0.config().orthonormalize(eps*GenTensor<double>::fac_reduce());	// this line break it!
    		ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

    		g1.config().orthonormalize(eps*GenTensor<double>::fac_reduce());
    		ASSERT_LT((g1.full_tensor_copy()-t1).normf(),eps);

    		t0+=t1;
    		g0.add_SVD(g1,eps);
    		ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);


    	} catch (const madness::TensorException& e) {
    		if (dim.size() != 0) std::cout << e;
    		EXPECT_EQ(dim.size(),0);
    	} catch(...) {
    		std::cout << "Caught unknown exception" << std::endl;
    		EXPECT_EQ(1,0);
    	}
    	}
    }

    TEST_P(BinaryGenTest, ScalarMultiplication) {
    	try {
    		// check for multiplication
    		Tensor<double> t1=.337*t0;
    		GenTensor<double> g1=0.3370*g0;
    		ASSERT_LT((g1.full_tensor_copy()-t1).normf(),eps);

    		// check for scaling
    		t0.scale(4.84);
    		g0.scale(4.84);
    		ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

    	} catch (const madness::TensorException& e) {
    		if (dim.size() != 0) std::cout << e;
    		EXPECT_EQ(dim.size(),0);
    	} catch(...) {
    		std::cout << "Caught unknown exception" << std::endl;
    		EXPECT_EQ(1,0);
    	}
    }

    // test for inner and outer product
    TEST_P(BinaryGenTest, InnerOuter) {
    	try {
    		// check for overlap
    		double inner_t=t0.trace(t1);
    		double inner_g=g0.trace_conj(g1);
    		ASSERT_LT(inner_t-inner_g,eps);

			// check for outer product
    		if ((t0.ndim()+t1.ndim()<=TENSOR_MAXDIM) and (tt==TT_FULL)) {
    			Tensor<double> t2=outer(t0,t1);
    			GenTensor<double> g2=outer(g0,g1);
    			ASSERT_LT((g2.full_tensor_copy()-t2).normf(),eps);
    		}

    	} catch (const madness::TensorException& e) {
    		if (dim.size() != 0) std::cout << e;
    		EXPECT_EQ(dim.size(),0);
    	} catch(...) {
    		std::cout << "Caught unknown exception" << std::endl;
    		EXPECT_EQ(1,0);
    	}
    }

    // checks for transform, general_transform, transform_dir
    TEST_P(UnaryGenTest, Transform) {
    	try {

    		// this is the transform tensor
    		Tensor<double> c(dim[0],dim[0]);
    		Tensor<double> cc[TENSOR_MAXDIM];
    		for (unsigned int idim=0; idim<dim.size(); idim++) {
    			cc[idim]=Tensor<double>(dim[0],dim[0]);
    			cc[idim].fillrandom();
    		}

    		c.fillindex();
    		c.scale(1.0/c.normf());

    		// check for transform
    		Tensor<double> t1=transform(t0,c);
    		GenTensor<double> g1=transform(g0,c);
			ASSERT_LT((g1.full_tensor_copy()-t1).normf(),eps);

    		// check for general transform
    		t1=general_transform(t0,cc);
    		g1=general_transform(g0,cc);
			ASSERT_LT((g1.full_tensor_copy()-t1).normf(),eps);

    		// check for transform in one direction
			for (int idim=0; idim<dim.size(); ++idim) {
				t1=transform_dir(t0,c,idim);
				g1=transform_dir(g0,c,idim);
				ASSERT_LT((g1.full_tensor_copy()-t1).normf(),eps);
			}


    	} catch (const madness::TensorException& e) {
    		if (dim.size() != 0) std::cout << e;
    		EXPECT_EQ(dim.size(),0);
    	} catch(...) {
    		std::cout << "Caught unknown exception" << std::endl;
    		EXPECT_EQ(1,0);
    	}
    }

    // checks for rank reduction, relies on Addition working properly
    TEST_P(BinaryGenTest, RankReduction) {
    	try {
    		// check for addition
    		t0+=t1;
    		g0+=g1;
    		g0.reduce_rank(eps);
    		ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

    		// check for subtraction
    		t0-=t1;
    		g0-=g1;
    		g0.reduce_rank(eps);
    		ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

    		// check for inplace stuff
    		t0.gaxpy(-0.7, t1, 0.1);
    		g0.gaxpy(-0.7, g1, 0.1);
    		g0.reduce_rank(eps);
    		ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

    	} catch (const madness::TensorException& e) {
    		if (dim.size() != 0) std::cout << e;
    		EXPECT_EQ(dim.size(),0);
    	} catch(...) {
    		std::cout << "Caught unknown exception" << std::endl;
    		EXPECT_EQ(1,0);
    	}
    }

    // checks for addition with slices
    TEST_P(BinaryGenTest, SliceAddition) {
        Tensor<double> t0_save=copy(t0);
		Tensor<double> t1_save=copy(t1);

		GenTensor<double> g0_save=copy(g0);
		GenTensor<double> g1_save=copy(g1);

		// need multiple passes to to random number generation
    	for (int ipass=0; ipass<10; ++ipass) {
			try {

				// reset for each pass
				t0=copy(t0_save);
				t1=copy(t1_save);
				g0=copy(g0_save);
				g1=copy(g1_save);

				long mindim= std::min(dim.front(), dim.back());

				long lo=(madness::RandomValue<long>() % 2);
				long hi=(madness::RandomValue<long>() % mindim-1);
				if (hi<lo) hi=lo;
				std::vector<Slice> s(dim.size(),Slice(lo,hi));

				// check for addition g0(s)+=g1(s)
				t0(s)+=t1(s);
				g0(s)+=g1(s);
				ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

				// special setup: shrink one of the tensors
				t1=copy(t1(s));
				g1=copy(g1(s));

				// check for addition g0(s)+=g1
				t0(s)+=t1;
				g0(s)+=g1;
				ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

				// check for subtraction g0(s)-=g1
				t0(s)-=t1;
				g0(s)-=g1;
				ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

				// check for addition g0+=g1(s)
				t1+=t0(s);
				g1+=g0(s);
				ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

				// check for subtraction g0-=g1(s)
				t1-=t0(s);
				g1-=g0(s);
				ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

				// check for zeroing slices, aka subtraction of that slice
				t0(s)=0.0;
				g0(s)=0.0;
				ASSERT_LT((g0.full_tensor_copy()-t0).normf(),eps);

			} catch (const madness::TensorException& e) {
				if (dim.size() != 0) std::cout << e;
				EXPECT_EQ(dim.size(),0);
			} catch(...) {
				std::cout << "Caught unknown exception" << std::endl;
				EXPECT_EQ(1,0);
			}
    	}
    }

    // let's keep construction as a typed test so that at least compiling will always work
    typedef ::testing::Types<float, double, float_complex, double_complex> GenTensorTestTypes;
    TYPED_TEST_CASE(GenTensorTest, GenTensorTestTypes);

    TYPED_TEST(GenTensorTest, Construction) {
        for (int ndim=0; ndim<=TENSOR_MAXDIM; ndim+=2) {	// even number of dimensions only
            try {
                std::vector<long> dim(ndim);
                for (int pass=0; pass<10; ++pass) {
                    long nelem = 1;
                    int alldim= (madness::RandomValue<long>()&0x5) + 1;
                    for (int i=0; i<ndim; ++i) {
                        dim[i] = alldim;
                        nelem *= dim[i];
                    }

                    // this is the full rank reference tensor
                    madness::Tensor<TypeParam> fullrank(dim);
                    fullrank.fillindex();
                    fullrank+=1.0;
                    fullrank.scale(1.0/fullrank.normf());

                    // Use ASSERT not EXPECT otherwise we get a huge volume of output
                    // from failing tests.  Also, subsequent tests usually rely on
                    // success of previous ones.

                    // verify default constructor
                    madness::GenTensor<TypeParam> empty;
                    ASSERT_EQ(empty.size(),0);
                    ASSERT_EQ(empty.rank(),0);
                    ASSERT_EQ(empty.ndim(),-1);
                    ASSERT_TRUE(empty.has_no_data());

                    // verify various "empty" constructors a given dimension but no tensor;
                    // Will result in the construction of an empty SRConf, which means that
                    // the gentensor "has data", which is a matrix of zeros, and rank is 0
                    {
                    	// TT_FULL always returns rank -1, irrespective of the content
						GenTensor<TypeParam> d2(dim, TT_FULL);
						ASSERT_TRUE(d2.has_data());
						ASSERT_EQ(d2.rank(),-1);

						GenTensor<TypeParam> d3(dim, TensorArgs(eps,TT_FULL));
						ASSERT_TRUE(d3.has_data());
						ASSERT_EQ(d3.rank(),-1);

						GenTensor<TypeParam> d4(TT_FULL, alldim, ndim);
						ASSERT_TRUE(d4.has_data());
						ASSERT_EQ(d4.rank(),-1);
                    } {
						GenTensor<TypeParam> d2(dim, TT_2D);
						ASSERT_TRUE(d2.has_data());
						ASSERT_EQ(d2.rank(),0);

						GenTensor<TypeParam> d3(dim, TensorArgs(eps,TT_2D));
						ASSERT_TRUE(d3.has_data());
						ASSERT_EQ(d3.rank(),0);

						GenTensor<TypeParam> d4(TT_2D, alldim, ndim);
						ASSERT_TRUE(d4.has_data());
						ASSERT_EQ(d4.rank(),0);
                    }

                    // verify various constructors with actual data
                    // test relies on the correct reconstruction to a full rank tensor
                    {
						// full rank -> TT_FULL
						madness::GenTensor<TypeParam> d1(fullrank,TensorArgs(eps,TT_FULL));
						ASSERT_EQ(d1.rank(),-1);
						ASSERT_NE(fullrank.ptr(),d1.config().vector_[0].ptr()); // Was it shallow?
						ASSERT_LT((fullrank-d1.full_tensor_copy()).normf(),eps);

						// full rank -> TT_2D
						madness::GenTensor<TypeParam> d_svd(fullrank,TensorArgs(eps,TT_2D));
						ASSERT_LT((fullrank-d_svd.full_tensor_copy()).normf(),eps);

						// TT_2D -> TT_2D
						madness::GenTensor<TypeParam> d_svd2(d_svd);
						ASSERT_EQ(d_svd.config().vector_[0].ptr(),d_svd2.config().vector_[0].ptr()); // shallow?
						ASSERT_LT((d_svd2.full_tensor_copy()-d_svd.full_tensor_copy()).normf(),eps);

                    }

                    // verify assigment operator
                    {
                    	madness::GenTensor<TypeParam> d0(fullrank,TensorArgs(eps,TT_2D));

                    	madness::GenTensor<TypeParam> d1(fullrank,TensorArgs(eps,TT_2D));
                    	d1=d0;
						ASSERT_LT((d0.full_tensor_copy()-d1.full_tensor_copy()).normf(),eps);
						ASSERT_EQ(d0.config().vector_[0].ptr(),d1.config().vector_[0].ptr()); // shallow?

                    	madness::GenTensor<TypeParam> t0(fullrank,TensorArgs(eps,TT_FULL));
                    	madness::GenTensor<TypeParam> t1(fullrank,TensorArgs(eps,TT_FULL));
						t1=t0;
						ASSERT_LT((t0.full_tensor_copy()-t1.full_tensor_copy()).normf(),eps);
						ASSERT_EQ(t0.config().vector_[0].ptr(),t1.config().vector_[0].ptr()); // shallow?

                    }

                    // verify copy construction and deep copy
                    {
                    	madness::GenTensor<TypeParam> d0(fullrank,TensorArgs(eps,TT_2D));

                    	madness::GenTensor<TypeParam> d1=d0;
						ASSERT_LT((d0.full_tensor_copy()-d1.full_tensor_copy()).normf(),eps);
						ASSERT_EQ(d0.config().vector_[0].ptr(),d1.config().vector_[0].ptr()); // shallow?

                    	d1=copy(d0);
						ASSERT_LT((d0.full_tensor_copy()-d1.full_tensor_copy()).normf(),eps);
						ASSERT_NE(d0.config().vector_[0].ptr(),d1.config().vector_[0].ptr()); // deep?

                    }
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



    /// test slices of GenTensors
    TYPED_TEST(GenTensorTest, SlicingConstruction) {
        for (int ndim=0; ndim<=TENSOR_MAXDIM; ndim+=2) {	// even number of dimensions only
            try {
                std::vector<long> dim(ndim,5);

				// this is the full rank reference tensor
				madness::Tensor<TypeParam> fullrank(dim);
				fullrank.fillindex();
				fullrank+=1.0;
				fullrank.scale(1.0/fullrank.normf());

				std::vector<Slice> s(ndim,Slice(0,3));

				// Use ASSERT not EXPECT otherwise we get a huge volume of output
				// from failing tests.  Also, subsequent tests usually rely on
				// success of previous ones.

				// g0=g1(s)
				Tensor<TypeParam> t0=copy(fullrank);
				Tensor<TypeParam> t1=t0(s);
				GenTensor<TypeParam> g0_svd(t0,TensorArgs(eps,TT_2D));
				GenTensor<TypeParam> g1_svd=g0_svd(s);
				GenTensor<TypeParam> g0_full(t0,TensorArgs(eps,TT_FULL));
				GenTensor<TypeParam> g1_full=g0_full(s);

				// check for shallowness -- should be deep(!)
				ASSERT_NE(g0_svd.config().vector_[0].ptr(),g1_svd.config().vector_[0].ptr());
				ASSERT_NE(g0_full.config().vector_[0].ptr(),g1_full.config().vector_[0].ptr());

				// check for numerical correctness
				ASSERT_LT((g1_svd.full_tensor_copy()-t1).normf(),eps);
				ASSERT_LT((g1_full.full_tensor_copy()-t1).normf(),eps);

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

}

int main(int argc, char** argv) {
	madness::default_random_generator.setstate(3149);

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
