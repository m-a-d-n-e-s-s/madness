#define NO_GENTENSOR
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/misc/ran.h>

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

    try {
        World world(SafeMPI::COMM_WORLD);
        startup(world,argc,argv);

        test_inner<double,double,1,false>(world);
        test_inner<double,double,1,true>(world);
        test_multi_to_multi_op<1>(world);
        test_multi_to_multi_op<2>(world);
        test_multi_to_multi_op<3>(world);
#if !HAVE_GENTENSOR
        test_inner<double,std::complex<double>,1,false>(world);
        test_inner<std::complex<double>,double,1,false>(world);
        test_inner<std::complex<double>,std::complex<double>,1,false>(world);
        test_inner<std::complex<double>,std::complex<double>,1,true>(world);
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
    catch (const char* s) {
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
