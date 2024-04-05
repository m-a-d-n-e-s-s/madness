#define NO_GENTENSOR
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/misc/ran.h>
#include <madness/world/test_utilities.h>

const double PI = 3.1415926535897932384;

using namespace madness;

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt)

template <typename T>
T complexify(T c) {
    return c;
}

template <> double_complex complexify<double_complex>(double_complex c) {
    return double_complex(c.real(),c.real()*c.real());
}

template <> float_complex complexify<float_complex>(float_complex c) {
    return c*float_complex(c.real(),c.real()*c.real());
}

/// struct to test the multi_to_multi_op_values
template<std::size_t NDIM>
struct many_to_many_op {

    many_to_many_op() : result_size(2) {}

    std::size_t result_size;
    std::size_t get_result_size() const {return result_size;}

    std::vector<madness::Tensor<double> > operator()(const madness::Key<NDIM> & key,
            const std::vector< madness::Tensor<double> >& t) const {
        std::vector<madness::Tensor<double> > result(result_size);
        result[0]=2.0*t[0];
        result[1]=t[1]+t[2];
        return result;
    }
};


template <typename T, std::size_t NDIM>
class Gaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;

    Gaussian(const coordT& center, double exponent, T coefficient)
            : center(center), exponent(exponent), coefficient(complexify(coefficient)) {};

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    };
};

/// Makes a square-normalized Gaussian with random origin and exponent
template <typename T, std::size_t NDIM>
Gaussian<T,NDIM>*
RandomGaussian(const Tensor<double> cell, double expntmax=1e5) {
    typedef Vector<double,NDIM> coordT;
    coordT origin;
    for (std::size_t i=0; i<NDIM; ++i) {
        origin[i] = RandomValue<double>()*(cell(i,1)-cell(i,0)) + cell(i,0);
    }
    double lo = log(0.01);
    double hi = log(expntmax);
    double expnt = exp(RandomValue<double>()*(hi-lo) + lo);
    T coeff = pow(2.0*expnt/PI,0.25*NDIM);
    //print("RandomGaussian: origin", origin, "expnt", expnt, "coeff", coeff);
    return new Gaussian<T,NDIM>(origin,expnt,coeff);
}


template <typename T, std::size_t NDIM>
void test_add(World& world) {

    const double thresh=1.e-7;
    Tensor<double> cell(NDIM,2);
    for (std::size_t i=0; i<NDIM; ++i) {
        cell(i,0) = -11.0-2*i;  // Deliberately asymmetric bounding box
        cell(i,1) =  10.0+i;
    }
    FunctionDefaults<NDIM>::set_cell(cell);
    FunctionDefaults<NDIM>::set_k(8);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_truncate_mode(1);

    std::size_t nvec=5;
    std::vector<Function<T,NDIM> > add1(nvec), add2(nvec), sum(nvec), diff(nvec);

    for (int i=0; i<nvec; ++i) {
        add1[i]=FunctionFactory<T,NDIM>(world).functor([] (const Vector<double,3>& r) {return r.normf();});
        add2[i]=FunctionFactory<T,NDIM>(world).functor([] (const Vector<double,3>& r) {return 2.0*r.normf();});
        sum[i]=FunctionFactory<T,NDIM>(world).functor([] (const Vector<double,3>& r) {return 3.0*r.normf();});
        diff[i]=FunctionFactory<T,NDIM>(world).functor([] (const Vector<double,3>& r) {return -1.0*r.normf();});
    }

    std::vector<Function<T,NDIM> > r1=add1+add2;
    std::vector<Function<T,NDIM> > r3=add1+add2[0];
    std::vector<Function<T,NDIM> > r5=add1[0]+add2;

    std::vector<Function<T,NDIM> > r2=add1-add2;
    std::vector<Function<T,NDIM> > r4=add1-add2[0];
    std::vector<Function<T,NDIM> > r6=add1[0]-add2;

    double error1=0.0,error2=0.0,error3=0.0,error4=0.0,error5=0.0,error6=0.0;
    for (int i=0; i<nvec; ++i) {
    	error1+=(r1[i]-sum[i]).norm2();
    	error3+=(r3[i]-sum[i]).norm2();
    	error5+=(r5[i]-sum[i]).norm2();

    	error2+=(r2[i]-diff[i]).norm2();
    	error4+=(r4[i]-diff[i]).norm2();
    	error6+=(r6[i]-diff[i]).norm2();
    }
    print("errors in add", error1,error3,error5,error2,error4,error6);

}



template <typename T, typename R, int NDIM, bool sym>
void test_inner(World& world) {
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > ffunctorT;
    typedef std::shared_ptr< FunctionFunctorInterface<R,NDIM> > gfunctorT;

    const double thresh=1.e-7;
    Tensor<double> cell(NDIM,2);
    for (std::size_t i=0; i<NDIM; ++i) {
        cell(i,0) = -11.0-2*i;  // Deliberately asymmetric bounding box
        cell(i,1) =  10.0+i;
    }
    FunctionDefaults<NDIM>::set_cell(cell);
    FunctionDefaults<NDIM>::set_k(8);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_truncate_mode(1);

    const int nleft=95, nright=sym ? nleft : 94;

    if (world.rank() == 0) 
        print("testing matrix_inner<",archive::get_type_name<T>(),",",archive::get_type_name<R>(),">","sym =",sym);

    START_TIMER;
    std::vector< Function<T,NDIM> > left(nleft);
    for (int i=0; i<nleft; ++i) {
        ffunctorT f(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),0.5));
        left[i] = FunctionFactory<T,NDIM>(world).functor(f);
    }
    std::vector< Function<R,NDIM> > right(nright);
    std::vector< Function<R,NDIM> >* pright = &right;
    if (sym) {
        pright = (std::vector< Function<R,NDIM> >*)(&left);
    }
    else {
        for (int i=0; i<nright; ++i) {
            gfunctorT f(RandomGaussian<R,NDIM>(FunctionDefaults<NDIM>::get_cell(),0.5));
            right[i] = FunctionFactory<R,NDIM>(world).functor(f);
        }
    }
    END_TIMER("project");

    START_TIMER;
    compress(world,left);
    compress(world,right);
    END_TIMER("compress");
    
    START_TIMER;
    Tensor<TENSOR_RESULT_TYPE(T,R)> rnew = matrix_inner(world,left,*pright,sym);
    END_TIMER("new");
    START_TIMER;
    Tensor<TENSOR_RESULT_TYPE(T,R)> rold = matrix_inner_old(world,left,*pright,sym);
    END_TIMER("old");

    if (world.rank() == 0) 
        print("error norm",(rold-rnew).normf(),"\n");
}

template<typename T, std::size_t NDIM>
int test_transform(World& world) {
    test_output to("testing transform");
    to.set_cout_to_terminal();
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > ffunctorT;


    const double thresh=1.e-7;
    Tensor<double> cell(NDIM,2);
    for (std::size_t i=0; i<NDIM; ++i) {
        cell(i,0) = -11.0-2*i;  // Deliberately asymmetric bounding box
        cell(i,1) =  10.0+i;
    }
    FunctionDefaults<NDIM>::set_cell(cell);
    FunctionDefaults<NDIM>::set_k(8);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_truncate_mode(1);


    const int nleft=RandomValue<int>()%10;
    const int nright=RandomValue<int>()%10;
    print("nleft, nright",nleft,nright);

    START_TIMER;
    std::vector< Function<T,NDIM> > left(nleft);
    for (int i=0; i<nleft; ++i) {
        ffunctorT f(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),0.5));
        left[i] = FunctionFactory<T,NDIM>(world).functor(f);
    }
    END_TIMER("project");
    to.checkpoint(true,"initial projection");

    Tensor<T> c(nleft,nright);
    auto result1=transform(world,left,c);
    to.checkpoint(true,"reference transform");

    change_tree_state(left,reconstructed);
    auto result2=transform_reconstructed(world,left,c,false);
    world.gop.fence();
    for (auto& r : result2) MADNESS_CHECK(r.get_impl()->get_tree_state()==redundant_after_merge);
    change_tree_state(result2,reconstructed);
    double err1=norm2(world, result1-result2);
    to.checkpoint(err1,thresh,"reference transform");

    change_tree_state(left,compressed);
    auto result3=transform_reconstructed(world,left,c);
    double err2=norm2(world, result1-result3);
    to.checkpoint(err2,thresh,"reference transform");


    return to.end();
}





template <typename T, typename R, int NDIM>
void test_cross(World& world) {
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > ffunctorT;
    typedef std::shared_ptr< FunctionFunctorInterface<R,NDIM> > gfunctorT;

    const double thresh=1.e-7;
    Tensor<double> cell(NDIM,2);
    for (std::size_t i=0; i<NDIM; ++i) {
        cell(i,0) = -11.0-2*i;  // Deliberately asymmetric bounding box
        cell(i,1) =  10.0+i;
    }
    FunctionDefaults<NDIM>::set_cell(cell);
    FunctionDefaults<NDIM>::set_k(8);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_truncate_mode(1);

    const int nleft=3, nright=3;

    if (world.rank() == 0)
        print("testing  cross product<",archive::get_type_name<T>(),",",archive::get_type_name<R>(),">");

    START_TIMER;
    std::vector< Function<T,NDIM> > left(nleft);
    for (int i=0; i<nleft; ++i) {
        ffunctorT f(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),0.5));
        left[i] = FunctionFactory<T,NDIM>(world).functor(f);
    }
    std::vector< Function<R,NDIM> > right(nright);
    for (int i=0; i<nright; ++i) {
    	gfunctorT f(RandomGaussian<R,NDIM>(FunctionDefaults<NDIM>::get_cell(),0.5));
    	right[i] = FunctionFactory<R,NDIM>(world).functor(f);
    }
    END_TIMER("project");

    START_TIMER;
    compress(world,left);
    compress(world,right);
    END_TIMER("compress");

    // cross product is anti-commuting
    std::vector<Function<TENSOR_RESULT_TYPE(T,R),NDIM> > result1=cross(left,right)+cross(right,left);
    double err1=norm2(world,result1);
    if (world.rank() == 0) print("error norm1",err1,"\n");

    // cross product with self vanishes
    std::vector<Function<TENSOR_RESULT_TYPE(T,R),NDIM> > result2=cross(left,left);
    double err2=norm2(world,result2);
    if (world.rank() == 0) print("error norm2",err2,"\n");

    // cross product with self vanishes
    std::vector<Function<TENSOR_RESULT_TYPE(T,R),NDIM> > result3=cross(left,right);
    Function<TENSOR_RESULT_TYPE(T,R), NDIM> reference0=left[1]*right[2] - left[2]*right[1];

    double err3=(reference0-result3[0]).norm2();
    if (world.rank() == 0) print("error norm3",err3,"\n");

    double err4=(result3[0]).norm2();
    if (world.rank() == 0) print("error norm4",err4," (should not be zero)\n");

}


template <typename T, int NDIM>
void test_rot(World& world) {
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > ffunctorT;

    const double thresh=1.e-7;
    Tensor<double> cell(NDIM,2);
    for (std::size_t i=0; i<NDIM; ++i) {
        cell(i,0) = -11.0-2*i;  // Deliberately asymmetric bounding box
        cell(i,1) =  10.0+i;
    }
    FunctionDefaults<NDIM>::set_cell(cell);
    FunctionDefaults<NDIM>::set_k(8);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_truncate_mode(1);

    const int nleft=3;

    if (world.rank() == 0)
        print("testing rot operator<",archive::get_type_name<T>(),",",archive::get_type_name<T>(),">");

    START_TIMER;
    std::vector< Function<T,NDIM> > left(nleft);
    for (int i=0; i<nleft; ++i) {
        ffunctorT f(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::get_cell(),0.5));
        left[i] = FunctionFactory<T,NDIM>(world).functor(f);
    }
    END_TIMER("project");

    START_TIMER;
    compress(world,left);
    END_TIMER("compress");

    // rot grad v vanishes
    std::vector<Function<T,NDIM> > result1=rot(grad(left[0]));
    double err1=norm2(world,result1);
    if (world.rank() == 0) print("error norm1",err1,"\n");

    // div rot v vanishes
    Function<T,NDIM> result2=div(rot(left));
    double err2=result2.norm2();
    if (world.rank() == 0) print("error norm2",err2,"\n");


}

template<typename T, int NDIM>
void test_matrix_mul_sparse(World &world) {
    typedef std::shared_ptr<FunctionFunctorInterface<T, NDIM> > ffunctorT;

    if (world.rank()==0) print("entering test_mul_sparse");
    const double thresh = 1.e-7;
    Tensor<double> cell(NDIM, 2);
    for (std::size_t i = 0; i < NDIM; ++i) {
        cell(i, 0) = -11.0 - 2 * i;  // Deliberately asymmetric bounding box
        cell(i, 1) = 10.0 + i;
    }
    FunctionDefaults<NDIM>::set_cell(cell);
    FunctionDefaults<NDIM>::set_k(8);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_truncate_mode(1);

    START_TIMER;
    const long nrow = 3;
    const long ncolumn = 4;
    std::vector<Function<T, NDIM> > row(nrow);
    for (long i = 0; i < nrow; ++i) {
        ffunctorT f(RandomGaussian<T, NDIM>(FunctionDefaults<NDIM>::get_cell(), 0.5));
        row[i] = FunctionFactory<T, NDIM>(world).functor(f);
    }
    std::vector<Function<T, NDIM> > column(ncolumn);
    for (long i = 0; i < ncolumn; ++i) {
        ffunctorT f(RandomGaussian<T, NDIM>(FunctionDefaults<NDIM>::get_cell(), 0.5));
        column[i] = FunctionFactory<T, NDIM>(world).functor(f);
    }
    END_TIMER("project");

    START_TIMER;
    auto result = matrix_mul_sparse<T, T, NDIM>(world, row, column, 0.0);
    END_TIMER("matrix_mul_sparse");



    MADNESS_CHECK(result.size()==nrow);
    for (long i = 0; i < nrow; ++i) {
        MADNESS_CHECK(result[i].size()==ncolumn);
        for (long j = 0; j < ncolumn; ++j) {
            Function<T, NDIM> tmp = row[i] * column[j];
            double err=(tmp-result[i][j]).norm2();
            MADNESS_CHECK(err<1.e-10);
        }
    }

    if (world.rank()==0) print("leaving test_mul_sparse");
}




template <std::size_t NDIM>
void test_multi_to_multi_op(World& world) {

    typedef Function<double,NDIM> functionT;
    typedef std::vector<Function<double,NDIM> > vecfuncT;
    typedef std::shared_ptr< FunctionFunctorInterface<double,NDIM> > ffunctorT;

    const double thresh=1.e-7;
    Tensor<double> cell(NDIM,2);
    for (std::size_t i=0; i<NDIM; ++i) {
        cell(i,0) = -11.0-2*i;  // Deliberately asymmetric bounding box
        cell(i,1) =  10.0+i;
    }
    FunctionDefaults<NDIM>::set_cell(cell);
    FunctionDefaults<NDIM>::set_k(8);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_truncate_mode(1);


    many_to_many_op<NDIM> op;
    vecfuncT vin(3);

    for (functionT& in : vin) {
        ffunctorT f(RandomGaussian<double,NDIM>(FunctionDefaults<NDIM>::get_cell(),0.5));
        in=FunctionFactory<double,NDIM>(world).functor(f);
    }
    refine_to_common_level(world,vin);

    vecfuncT vout=multi_to_multi_op_values(op, vin);

    std::vector<functionT> result(2);
    result[0]=2.0*vin[0]-vout[0];
    result[1]=vout[1]-vin[1]-vin[2];

    double norm_in=norm2(world,vin);
    double norm_out=norm2(world,vout);
    double error=norm2(world,result);
    if (world.rank()==0) print("error in multi_to_multi ",error,norm_in,norm_out);

}

int main(int argc, char**argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    bool smalltest = false;
    if (getenv("MAD_SMALL_TESTS")) smalltest=true;
    for (int iarg=1; iarg<argc; iarg++) if (strcmp(argv[iarg],"--small")==0) smalltest=true;
    std::cout << "small test : " << smalltest << std::endl;
    if (smalltest) return 0;
    madness::default_random_generator.setstate(int(cpu_time()*1000)%4149);


    try {
        startup(world,argc,argv);

        // test_add<double,3>(world);
        // test_add<std::complex<double>,3 >(world);

        test_inner<double,double,1,false>(world);
        test_inner<double,double,1,true>(world);
        test_multi_to_multi_op<1>(world);
        test_multi_to_multi_op<2>(world);

        test_cross<double,double,2>(world);
        test_cross<std::complex<double>,double,2>(world);
        test_cross<std::complex<double>,std::complex<double>,2>(world);

        test_transform<double,2>(world);
        test_transform<double,2>(world);
        test_transform<double,2>(world);
        test_transform<double,2>(world);

        test_rot<double,3>(world);
        test_rot<std::complex<double>,3>(world);

        test_matrix_mul_sparse<double,2>(world);
        test_matrix_mul_sparse<double,3>(world);

        if (!smalltest) test_multi_to_multi_op<3>(world);
#if !HAVE_GENTENSOR
        test_inner<double,std::complex<double>,1,false>(world);
        if (!smalltest) {
            test_inner<std::complex<double>,double,1,false>(world);
            test_inner<std::complex<double>,std::complex<double>,1,false>(world);
            test_inner<std::complex<double>,std::complex<double>,1,true>(world);
        }
#endif
    }
    catch (const SafeMPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }
    finalize();

    return 0;
}
