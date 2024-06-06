//
// Created by Florian Bischoff on 12/19/23.
//

#include <madness/mra/mra.h>
#include <madness/chem/projector.h>
#include <madness/world/test_utilities.h>


using namespace madness;

template<typename T, std::size_t NDIM>
int test_projector(World& world) {
    test_output t1("testing projector for dimension " + std::to_string(NDIM));
//    t1.set_cout_to_terminal();

    FunctionDefaults<NDIM>::set_tensor_type(TT_FULL);

    auto g1=[](const Vector<double,NDIM>& r){return exp(-inner(r,r));};
    auto g2=[](const Vector<double,NDIM>& r){return exp(-2.0*inner(r,r));};
    auto g3=[](const Vector<double,NDIM>& r){return exp(-3.0*inner(r,r));};

    // set up an orthonormal basis for projection
    std::vector<Function<T,NDIM>> vphi;
    vphi.push_back(FunctionFactory<T,NDIM>(world).functor(g1));
    vphi.push_back(FunctionFactory<T,NDIM>(world).functor(g2));
    vphi.push_back(FunctionFactory<T,NDIM>(world).functor(g3));
    vphi=orthonormalize_symmetric(vphi);
    auto S=matrix_inner(world,vphi,vphi);
    print("overlap");
    print(S);

    std::vector<Function<T,NDIM>> pspace1={vphi[0],vphi[1]};
    std::vector<Function<T,NDIM>> pspace2={vphi[2]};

    // P1 and P2 project onto orthogonal subspaces
    Projector<T,NDIM> P1(pspace1);
    Projector<T,NDIM> P2(pspace2);
    QProjector<T,NDIM> Q1(world,pspace1);

    // set up trial function
    auto trial1=[](const Vector<double,NDIM>& r){return inner(r,r)*exp(-inner(r,r));};
    auto trial2=[](const Vector<double,NDIM>& r){return r.normf()*exp(-inner(r,r));};
    Function<T,NDIM> f1=FunctionFactory<T,NDIM>(world).functor(trial1);
    Function<T,NDIM> f2=FunctionFactory<T,NDIM>(world).functor(trial2);
    std::vector<Function<T,NDIM>> f={f1,f2};


    // check that the projector is indeed a projector: idempotency
    std::vector<Function<T,NDIM>> pf=P1(f);
    std::vector<Function<T,NDIM>> ppf=P1(pf);
    double err=norm2(world,(pf-ppf))/norm2(world,pf);
    print("err",err);
    t1.checkpoint(fabs(err)<FunctionDefaults<NDIM>::get_thresh(),"P1 projector is a projector");

    // check that P1 and P2 are orthogonal
    std::vector<Function<T,NDIM>> p2f=P2(f);
    double err2=matrix_inner(world,p2f,pf).normf();
    print("err2",err2);
    t1.checkpoint(fabs(err2)<FunctionDefaults<NDIM>::get_thresh(),"P1 and P2 are orthogonal");

    // check that P2 projects into a subspace of Q1
    std::vector<Function<T,NDIM>> q1p2f=Q1(p2f);
    double err3=norm2(world,(p2f-q1p2f))/norm2(world,p2f);
    print("err3",err3);
    t1.checkpoint(fabs(err3)<FunctionDefaults<NDIM>::get_thresh(),"P2 projects into a subspace of Q1");






    return t1.end();
}

template<typename T, std::size_t NDIM>
int test_projector_outer(World& world) {
    test_output t1("testing projector_outer for dimension " + std::to_string(NDIM));
    constexpr std::size_t LDIM=NDIM/2;
    static_assert(2*LDIM==NDIM);

    auto g1=[](const Vector<double,LDIM>& r){return exp(-inner(r,r));};
    auto g_hidim=[](const Vector<double,NDIM>& r){return 2.0*exp(-3.0*inner(r,r));};
    Function<double,LDIM> f1=FunctionFactory<double,LDIM>(world).f(g1);
    Function<double,NDIM> f_hidim=FunctionFactory<double,NDIM>(world).f(g_hidim);


    // compare explicit SO projector Q12 and outer product projector Q1Q2
    StrongOrthogonalityProjector<double,LDIM> Q1(world);
    Q1.set_spaces({f1});

    QProjector<double,LDIM> q(world,{f1});
    auto Q2=outer(q,q);

    auto Q1f=Q1(f_hidim);
    auto Q2f=Q2(f_hidim);
    double err=(Q1f-Q2f).norm2();
    print("error",err);
    double norm1=Q1f.norm2();
    double norm2=Q2f.norm2();
    print("norm1/2",norm1,norm2);
    double trace1=Q1f.trace();
    double trace2=Q2f.trace();
    print("trace1/2",trace1,trace2);

    t1.checkpoint(norm1-norm2,FunctionDefaults<NDIM>::get_thresh(),"Q1 direct and Q2 outer are the same");
    t1.checkpoint(trace1-trace2,FunctionDefaults<NDIM>::get_thresh(),"Q1 direct and Q2 outer are the same");
    // loosen threshold due to outer product
    t1.checkpoint(err,FunctionDefaults<NDIM>::get_thresh()*3.0,"Q1 direct and Q2 outer are the same");

    return t1.end();
}

template<typename T, std::size_t NDIM>
int test_Q12_projector(World& world) {
    test_output t1("testing Q12 projector for dimension "+std::to_string(NDIM));
    t1.set_cout_to_terminal();
    constexpr std::size_t LDIM=NDIM/2;
    static_assert(NDIM==LDIM*2);
    double thresh=FunctionDefaults<NDIM>::get_thresh();
    FunctionDefaults<NDIM>::set_tensor_type(TT_2D);
    FunctionDefaults<LDIM>::set_tensor_type(TT_FULL);

    auto g1=[](const Vector<double,LDIM>& r){return exp(-inner(r,r));};
    auto g2=[](const Vector<double,LDIM>& r){return exp(-2.0*inner(r,r));};
    auto g_hidim=[](const Vector<double,NDIM>& r){return exp(-3.0*inner(r,r));};

    // set up an orthonormal basis for projection
    std::vector<Function<T,LDIM>> vphi;
    vphi.push_back(FunctionFactory<T,LDIM>(world).functor(g1));
    vphi.push_back(FunctionFactory<T,LDIM>(world).functor(g2));
    vphi=orthonormalize_symmetric(vphi);
//    auto S=matrix_inner(world,vphi,vphi);
//    print("overlap");
//    print(S);

    StrongOrthogonalityProjector<T,LDIM> SO(world);
    SO.set_spaces(vphi);

    Function<T,NDIM> f=FunctionFactory<T,NDIM>(world).functor(g_hidim);
    double fnorm=f.norm2();
    print("fnorm",fnorm);

    // check that the projector is indeed a projector
    Function<T,NDIM> f1=SO(f);
    double sonorm=f1.norm2();
    print("Q12(f) norm",sonorm);
    double refnorm;
    if (NDIM==2) {
        refnorm=0.0028346312885398958; // according to mathematica
        print("err0",fabs(sonorm-refnorm));
        t1.checkpoint(fabs(sonorm-refnorm)<thresh,"SO projector is correct");
    } else if (NDIM==4) {
        refnorm=5.258329e-03;      // according to madness ..
        print("err0",fabs(sonorm-refnorm));
        t1.checkpoint(fabs(sonorm-refnorm)<thresh,"SO projector is correct");
    }

    Function<T,NDIM> f2=SO(f1);
    double err=(f1-f2).norm2()/f.norm2();
    print("err",err);
    t1.checkpoint(fabs(err)<thresh,"SO projector is a projector");


    // check projector being consistent
    // SO(f) = f - O1(f) - O2(f) + O1O2(f)
    Projector<T,LDIM> O1(vphi);
    Projector<T,LDIM> O2(vphi);
    O1.set_particle(0);
    O2.set_particle(1);
    Function<T,NDIM> f3=f-O1(f)-O2(f)+O1(O2(f));
    double err1=(f1-f3).norm2()/f.norm2();
    print("err1",err1);
    t1.checkpoint(fabs(err1)<thresh*20.0,"SO projector is consistent with 1-O1-O2+O1O2");


    return t1.end();
}

int main(int argc, char**argv) {

    World& world=initialize(argc,argv);
    commandlineparser parser(argc,argv);

    // the parameters
    double L=40;
    long k=5;
    double thresh=1.e-4;

    if (parser.key_exists("k")) k=atoi(parser.value("k").c_str());
    if (parser.key_exists("thresh")) k=atof(parser.value("thresh").c_str());
    print("k, thresh", k, thresh);

    srand(time(nullptr));
    startup(world,argc,argv);

    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<2>::set_thresh(thresh);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<4>::set_thresh(thresh);
    FunctionDefaults<5>::set_thresh(thresh);
    FunctionDefaults<6>::set_thresh(thresh);

    FunctionDefaults<1>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<2>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<4>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<5>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<6>::set_cubic_cell(-L/2,L/2);

    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<2>::set_k(k);
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<4>::set_k(k);
    FunctionDefaults<5>::set_k(k);
    FunctionDefaults<6>::set_k(k);

    print("entering testsuite for 6-dimensional functions\n");
    print("k            ",k);
    print("thresh       ",thresh);
    print("boxsize      ",L);
    print("tensor type: ", FunctionDefaults<6>::get_tensor_type());
    print("");


    int error=0;
    {
        error+=test_projector<double,2>(world);
        error+=test_projector<double,3>(world);
        error+=test_projector<double,4>(world);

        error+=test_projector_outer<double,2>(world);

        if (HAVE_GENTENSOR) {
            error+=test_Q12_projector<double,2>(world);
            error+=test_Q12_projector<double,4>(world);
//        error+=test_Q12_projector<double,6>(world);
        }

    }


    world.gop.fence();
    finalize();

    return error;
}

