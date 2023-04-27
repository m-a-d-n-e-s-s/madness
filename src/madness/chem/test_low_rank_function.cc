//
// Created by Florian Bischoff on 4/24/23.
//


#include<madness.h>
#include<madness/chem/ccpairfunction.h>
#include<madness/chem/correlationfactor.h>
#include<madness/chem/electronic_correlation_factor.h>
#include<madness/chem/CCStructures.h>
#include<madness/chem/CCPotentials.h>
#include<madness/chem/projector.h>
#include<madness/chem/localizer.h>

#include<madness/world/test_utilities.h>
#include <random>

template<std::size_t NDIM>
struct gauss {
    double a=2.0;
    gauss() = default;
    double operator()(const Vector<double,NDIM>& r) const {
        return exp(-a*inner(r,r));
    }
};

struct hidim_orbitals {
    std::shared_ptr<real_convolution_3d> f12;
    real_function_3d phi1, phi2;
    hidim_orbitals(std::shared_ptr<real_convolution_3d>& f, const real_function_3d& phi1, const real_function_3d& phi2)
            : f12(f), phi1(phi1), phi2(phi2) {
    }

    double operator()(const Vector<double,6>& coord) const {

        Vector<double,3> left,right;
        for (int i=0; i<3; ++i) {
            left[i]=coord[i];
            right[i]=coord[i+3];
        }

        auto f12_functor=[](const Vector<double,6>& r) {
            Vector<double,3> r1,r2;
            for (int i=0; i<3; ++i) {
                r1[i]=r[i];
                r2[i]=r[i+3];
            }
            return exp(-(r1-r2).normf());
        };

        return phi1(left)*phi2(right)*f12_functor(coord);
    }

};

template<typename T, std::size_t LDIM>
void orthonormalize_svd(World& world, std::vector<Function<T,LDIM>>& g, std::vector<Function<T,LDIM>>& h, Tensor<double>& s) {
    /**
     *  |g >s< h| = |g_ortho><g_ortho | g> s < h | h_ortho ><h_ortho |
     *           = |g_ortho> gg s hh <h_ortho |
     *           = |g_ortho> U s VT <h_ortho |
     */
    timer t(world);
    std::vector<Function<T,LDIM>> g_ortho=orthonormalize_canonical(g,1.e-8);
    std::vector<Function<T,LDIM>> h_ortho=orthonormalize_canonical(h,1.e-8);
    t.tag("orthonormalize_canonical");
    auto gg=matrix_inner(world,g_ortho,g);
    auto hh=matrix_inner(world,h,h_ortho);
    t.tag("matrix_inner");
    for (int i=0; i<gg.dim(0); ++i) gg(_,i)*=s(i);
    auto ovlp=inner(gg,hh);
    Tensor<double> U,VT;
    svd(ovlp,U,s,VT);
    auto V=transpose(VT);

    // truncate
//    for (int i=1; i<s.size(); ++i) {
//        if (s[i]<1.e-2) {
//            s=s(Slice(0,i-1));
//            U=U(_,Slice(0,i-1));
//            V=V(_,Slice(0,i-1));
//            print("truncating svd at i",i);
//            break;
//        }
//    }
    g=transform(world,g_ortho,U);
    h=transform(world,h_ortho,V);
    t.tag("transform");

    // test
//    auto gg1=matrix_inner(world,g,g);
//    auto hh1=matrix_inner(world,h,h);
//    for (int i=0; i<gg1.dim(0); ++i) {
//        gg1(i,i)-=1.0;
//        hh1(i,i)-=1.0;
//    }
//    double gmatnorm=gg1.normf()/gg1.size();
//    double hmatnorm=hh1.normf()/hh1.size();
//    print("g/h identity",gmatnorm,hmatnorm);
//    print("singular values",s);

}


template<typename T, std::size_t LDIM>
std::vector<Function<T,LDIM>> form_inner(const hidim_orbitals& f, const std::vector<Function<T,LDIM>>& vec) {
    double cpu0=cpu_time();
    World& world=vec.front().world();
    auto tmp1=vec*f.phi1;
    std::vector<Function<T,LDIM>> result=apply(world,(*f.f12),tmp1);
    result=result*f.phi2;
    double cpu1=cpu_time();
    std::printf("form inner finished after %4.2fs \n",cpu1-cpu0);
    return result;
}

// note that we assume f is symmetric!!!!!
template<typename T, std::size_t LDIM>
std::vector<Function<T,LDIM>> form_inner(const Function<T,2*LDIM>& f, const std::vector<Function<T,LDIM>>& vec) {
    double cpu0=cpu_time();
    std::vector<Function<T,LDIM>> result(vec.size());
    if constexpr (std::is_same<Function<T,LDIM>,Function<T,1>>::value) {
        for (int i=0; i<vec.size(); ++i) result[i]=inner(f,vec[i],{1},{0});
    } else {
        for (int i=0; i<vec.size(); ++i) result[i]=inner(f,vec[i],{2,3},{0,1});
    }
    double cpu1=cpu_time();
    std::printf("form inner finished after %4.2fs \n",cpu1-cpu0);
    return result;
}

// note that we assume f is symmetric!!!!!
template<typename T, std::size_t LDIM>
typename std::enable_if<LDIM==3,std::vector<Function<T,LDIM>>>::type form_inner(const real_convolution_3d& f, const std::vector<Function<T,LDIM>>& vec) {
    double cpu0=cpu_time();
    World& world=vec.front().world();
    std::vector<Function<T,LDIM>> result=apply(world,f,vec);
    double cpu1=cpu_time();
    std::printf("form inner finished after %4.2fs \n",cpu1-cpu0);
    return result;
}

template<typename T, std::size_t LDIM>
class LowRank {
public:

    Tensor<double> s;
    std::vector<Function<T,LDIM>> g,h;

    LowRank(Tensor<double>& s, std::vector<Function<T,LDIM>> g, std::vector<Function<T,LDIM>> h)
            : s(s), g(g), h(h) {}

//    LowRank(World& world, long n) : world(world) {
//        g= zero_functions_compressed<T,LDIM>(world,n);
//        h= zero_functions_compressed<T,LDIM>(world,n);
//    }

    LowRank() =default;      // Default constructor necessary for storage in vector

    LowRank(const LowRank& a) : s(a.s), g(copy(a.g.front().world(),a.g)), h(copy(a.h.front().world(),a.h)) {} // Copy constructor necessary

    LowRank& operator=(const LowRank& f) { // Assignment required for storage in vector
        LowRank ff(f);
        std::swap(ff.s,s);
        std::swap(ff.g,g);
        std::swap(ff.h,h);
        return *this;
    }

    template<typename hidimT>
    void project(World& world, const hidimT& f, const std::vector<Function<T,LDIM>>& regular_grid_Y, const long ntrial) {

        timer t1(world);
        auto Y=random_vectors(world, regular_grid_Y, ntrial);
        t1.tag("Yforming");
        print("Y.size()",Y.size());


        Localizer localizer;

//    std::vector<Function<double,LDIM>> g=orthonormalize_canonical(Y,1.e-8);
        g=orthonormalize_cd(Y);
        print("g.size()",g.size());
        t1.tag("Y orthonormalizing");

        bool localize=false;
        if (localize) {
            localizer.set_method("boys");
            MolecularOrbitals<double,LDIM> mos;
            mos.set_mos(g);
            std::vector<int> set(g.size(),0);
            mos.update_localize_set(set);
            localizer.localize(mos,true);
            g=mos.get_mos();
            t1.tag("Y localizing");
        }

        h=form_inner(f,g);
        s=Tensor<double>(g.size());
        s=1.0;
        t1.tag("Y backprojection");
    }

    void orthonormalize_svd() {
        ::orthonormalize_svd(g.front().world(),g,h,s);
    }

    template<typename hidimT>
    void optimize(const hidimT& f, const std::vector<Function<T,LDIM>>& regular_grid_Y) {
        World& world=g.front().world();
        for (int iopt=0; iopt<8; ++iopt) {
            timer t(world);
            std::vector<Function<double,LDIM>> htmp(g.size()), gtmp(g.size());

            for (int i=0; i<h.size(); ++i) h[i]*=1.0/s[i];
            gtmp=form_inner(f,h);       // only correct if f is symmetric!!!
            for (int i=0; i<g.size(); ++i) g[i]*=1.0/s[i];
            htmp=form_inner(f,g);

            g=gtmp;
            h=htmp;

            if (g.size()>1) orthonormalize_svd();
            plot_plane<2*LDIM>(world,*this,"lrf_iter"+std::to_string(iopt));

            if (iopt%2==0) compute_error(f,regular_grid_Y);
            t.end("finished optimization iteration "+std::to_string(iopt));
        }
    }

    long rank() const {return s.size();}

    T operator()(const Vector<double,2*LDIM>& coord) const {
        Vector<double,LDIM> left,right;
        for (int i=0; i<LDIM; ++i) {
            left[i]=coord[i];
            right[i]=coord[i+LDIM];
        }
        T result=0.0;
        for (int i=0; i<rank(); ++i) result+=s(i)*g[i](left)*h[i](right);
        return result;
    }

    // check accuracy by applying hi/lo-dim function on a set of trial functions
    // cf Halko, sec 4.3
    template<typename hidimT>
    double compute_error(const hidimT& f, const std::vector<Function<T,LDIM>>& regular_grid_Y) const {
        const long ntrial=10;
        World& world=g.front().world();
        timer t(world);
        auto trial= random_vectors(world,regular_grid_Y,ntrial);

        // result[j,1] = \sum_2 f(1,2) t(2,j) - \sum_2r g(1,r)s(r)h(r,2) t(2,j)
        // ref[j,1] = \sum_2 f(1,2) t(2,j)
        const std::vector<Function<T,LDIM>> ref=form_inner(f,trial);

        // check[j,1] = \sum_jr g(1,r)s(r)h(r,2) t(2,j)
        Tensor<double> ht=matrix_inner(world,h,trial);     // ht(r,j)
        for (int r=0; r<s.size(); ++r) ht(r,_)*=s(r);
        std::vector<Function<T,LDIM>> check=transform(world,g,ht);

        std::vector<double> norms=norm2s(world,ref-check);
        double max=*std::max_element(norms.begin(),norms.end());
        double err=10.0*sqrt(2.0/constants::pi)*max;
        std::stringstream ss;

        ss << "rank and error in f_approx " <<  g.size() << " " << std::scientific << err;
        t.tag(ss.str());
        return err;
    }


//    LowRank operator-(const LowRank& b) const { // Operator- necessary
//        return LowRank(g-b.g,h-b.h);
//    }
//
//    LowRank& operator+=(const LowRank& b) { // Operator+= necessary
//        g+=b.g;
//        h+=b.h;
//        return *this;
//    }

    LowRank operator*(double a) const { // Scale by a constant necessary
        return LowRank(s*a,g,h);
    }

//    double get() const {return x;}
};

// This interface is necessary to compute inner products
template<typename T, std::size_t LDIM>
double inner(const LowRank<T,LDIM>& a, const LowRank<T,LDIM>& b) {
    World& world=a.world;
    return (matrix_inner(world,a.g,b.g).emul(matrix_inner(world,a.h,b.h))).sum();
}



Tensor<double> gaussian_random_distribution(double mean, double variance, long n) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{mean, variance};
    Tensor<double> result(n);
    for (int i = 0; i < n; ++i) result(i)=d(gen);
    return result;
}

template<typename T, std::size_t LDIM>
typename std::enable_if<(LDIM==3), real_convolution_3d>::type hidim(World& world) {
    auto f_op= SlaterOperator(world,1.0,1.e-6,FunctionDefaults<3>::get_thresh());
    auto f12_functor=[](const Vector<double,2*LDIM>& r) {
        Vector<double,LDIM> r1,r2;
        for (int i=0; i<LDIM; ++i) {
            r1[i]=r[i];
            r2[i]=r[i+LDIM];
        }
        return exp(-(r1-r2).normf());
    };
    plot_plane<2*LDIM>(world,f12_functor,"f12_functor");
    return f_op;
}

template<typename T, std::size_t LDIM>
typename std::enable_if<(LDIM<3), Function<T,2*LDIM>>::type hidim(World& world) {

    auto f12_functor=[](const Vector<double,2*LDIM>& r) {
        Vector<double,LDIM> r1,r2;
        for (int i=0; i<LDIM; ++i) {
            r1[i]=r[i];
            r2[i]=r[i+LDIM];
        }
        return exp(-(r1-r2).normf())* exp(-0.2*inner(r,r));
    };
    auto f=FunctionFactory<double,2*LDIM>(world).functor(f12_functor);
    plot_plane<2*LDIM>(world,f12_functor,"f12_functor");
//    Function<double,2*LDIM> f;
    return f;
}

// nran random vectors on a regular grid
template<typename T, std::size_t LDIM>
std::vector<Function<T,LDIM>> random_vectors(World &world,
                                             const std::vector<Function<T,LDIM>>& regular_grid_Y,
                                             const long nran) {
    const long n = regular_grid_Y.size();
    Tensor<double> factors(n,nran);
    for (long i=0; i<nran; ++i) factors(_,i) = gaussian_random_distribution(0.0, 1.0, n);
    /// Transforms a vector of functions according to new[i] = sum[j] old[j]*c[j,i]
    return truncate(transform(world,regular_grid_Y,factors));
};

// check accuracy by applying hi/lo-dim function on a set of trial functions
// cf Halko, sec 4.3
template<typename T, std::size_t LDIM, typename hidimT>
double compute_error_implicitly(const Tensor<T>& s,
                                const std::vector<Function<T,LDIM>>& g,
                                const std::vector<Function<T,LDIM>>& h,
                                const hidimT& f,
                                const std::vector<Function<T,LDIM>>& regular_grid_Y) {
    const long ntrial=10;
    World& world=g.front().world();
    auto trial= random_vectors(world,regular_grid_Y,ntrial);

    // result[j,1] = \sum_2 f(1,2) t(2,j) - \sum_2r g(1,r)s(r)h(r,2) t(2,j)
    // ref[j,1] = \sum_2 f(1,2) t(2,j)
    const std::vector<Function<T,LDIM>> ref=form_inner(f,trial);

    // check[j,1] = \sum_jr g(1,r)s(r)h(r,2) t(2,j)
    Tensor<double> ht=matrix_inner(world,h,trial);     // ht(r,j)
    for (int r=0; r<s.size(); ++r) ht(r,_)*=s(r);
    std::vector<Function<T,LDIM>> check=transform(world,g,ht);

    std::vector<double> norms=norm2s(world,ref-check);
    double max=*std::max_element(norms.begin(),norms.end());
    double err=10.0*sqrt(2.0/constants::pi)*max;
    return err;

}

template<typename T, std::size_t LDIM>
typename std::enable_if<LDIM<=2,double>::type compute_error (const Tensor<T>& s,
                                                             const std::vector<Function<T,LDIM>>& g,
                                                             const std::vector<Function<T,LDIM>>& h,
                                                             const Function<double,2*LDIM>& f,
                                                             const std::vector<Function<T,LDIM>>& regular_grid_Y,
                                                             std::string name="") {
    auto fapprox=s[0]*hartree_product(g[0],h[0]);
    for (int i=1; i<g.size(); ++i) fapprox+=s[i]*hartree_product(g[i],h[i]);
    double err=(f-fapprox).norm2();
    double err_implicit= compute_error_implicitly(s,g,h,f,regular_grid_Y);
    print("rank and error in f_approx_opt", g.size(), err,err_implicit);
    if (not name.empty())
        plot<2*LDIM>({f,fapprox,f-fapprox},name,std::vector<std::string>({"adsf","asdf","diff"}));
    return err;
};


template<typename T, std::size_t LDIM>
typename std::enable_if<LDIM==3,double>::type compute_error (const Tensor<T>& s,
                                                             const std::vector<Function<T,LDIM>>& g,
                                                             const std::vector<Function<T,LDIM>>& h,
                                                             const SeparatedConvolution<T,LDIM>& f,
                                                             const std::vector<Function<T,LDIM>>& regular_grid_Y,
                                                             std::string name="") {
    double cpu0=cpu_time();
    double err_implicit= compute_error_implicitly(s,g,h,f,regular_grid_Y);
    double err=-1.0;
    double cpu1=cpu_time();
    print("rank and error in f_approx_opt", g.size(), err,err_implicit, "after ",cpu1-cpu0);
    return err;
};

template<std::size_t NDIM>
struct cartesian_grid {
    Vector<double,NDIM> lovec,hivec;
    std::vector<long> stride;
    long index=0;
    long n_per_dim;
    long total_n;
    Vector<double,NDIM> increment;

    cartesian_grid(const long n_per_dim, const double lo, const double hi)
            : n_per_dim(n_per_dim) {
        lovec.fill(lo);
        hivec.fill(hi);
        increment=(hivec-lovec)*(1.0/double(n_per_dim-1));
        stride=std::vector<long>(NDIM,1l);
        total_n=std::pow(n_per_dim,NDIM);
        for (long i=NDIM-2; i>=0; --i) stride[i]=n_per_dim*stride[i+1];
    }

    void operator++() {
        index++;
    }

    bool operator()() const {
        return index < total_n;
    }

    Vector<double,NDIM> get_coordinates() const {
        Vector<double,NDIM> tmp(NDIM);
        for (int idim=0; idim<NDIM; ++idim) {
            tmp[idim]=(index/stride[idim])%n_per_dim;
        }
        return lovec+tmp*increment;
    }

};

template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM>> uniformly_distributed_slater(World& world,
                     const double lo, const double hi, const long n_per_dim) {
    Vector<double,NDIM> R;
    auto sl=[&R](const Vector<double,NDIM>& r) {
        return exp(-sqrt(inner(r-R,r-R)+1.e-3));
    };

    std::vector<Function<T,NDIM>> result;
    cartesian_grid<NDIM> cg(n_per_dim,lo,hi);
    for (auto c=cg; c(); ++c) {
        R=c.get_coordinates();
        result.push_back(FunctionFactory<T,NDIM>(world).functor(sl));
    }
    return result;

}

int test_lowrank_function(World& world) {
    madness::default_random_generator.setstate(int(cpu_time())%4149);

    print("");
    constexpr std::size_t LDIM=3;
    long n_per_dim=10;
    long ntrial=50;



    print("LDIM, n_per_dim, ntrial",LDIM,n_per_dim,ntrial);
    constexpr std::size_t NDIM=2*LDIM;

    Function<double,LDIM> phi=FunctionFactory<double,LDIM>(world)
            .functor([](const Vector<double,LDIM>& r) {return exp(-2.0*inner(r,r));});

//    auto f=hidim<double,LDIM>(world);
    std::shared_ptr<real_convolution_3d> f12(SlaterOperatorPtr(world,1.0,1.e-4,FunctionDefaults<NDIM>::get_thresh()));
    auto f=hidim_orbitals(f12,phi,phi);

    plot_plane<2*LDIM>(world,f,"hidim_orbitals");


    double radius=2.0;
    double lo=-radius;
    double hi=radius;


    auto regular_grid_Y=uniformly_distributed_slater<double,LDIM>(world,lo,hi,n_per_dim);

    std::vector<double> phi1_values(regular_grid_Y.size());
    cartesian_grid<LDIM> cg(n_per_dim,lo,hi);
    for (auto c=cg; c(); ++c) {
        auto R=c.get_coordinates();
        phi1_values[c.index]=phi(R);
    }
    scale(world,regular_grid_Y,phi1_values);
    regular_grid_Y=regular_grid_Y*phi;

    print("regular_grid_Y.size()",regular_grid_Y.size());

    LowRank<double,LDIM> lrf;
    lrf.project(world,f,regular_grid_Y,ntrial);

    lrf.compute_error(f,regular_grid_Y);
    lrf.orthonormalize_svd();
    lrf.compute_error(f,regular_grid_Y);

    plot_plane<2*LDIM>(world,lrf,"lrf0");

    print("\nstarting optimization\n");

    lrf.optimize(f,regular_grid_Y);

    lrf.compute_error(f,regular_grid_Y);
    plot_plane<2*LDIM>(world,lrf,"lrf1");


    return 0;
}

int test_6d_convolution(World& world) {
    timer t1(world);
    real_function_3d one=real_factory_3d(world).functor([](const coord_3d& r){return 1.0;});
    real_function_3d phi=real_factory_3d(world).functor([](const coord_3d& r){return exp(-2.0*inner(r,r));});
//    auto f12= SlaterF12Operator(world,1.0,1.e-4,1.e-3);
    auto g12=CoulombOperator(world,1.e-4,1.e-4);
    g12.particle()=1;
    auto f12=TwoElectronFactory(world).dcut(1.e-7).f12().gamma(1.0);
    {
        real_function_6d f = CompositeFactory<double, 6, 3>(world).g12(f12).particle1(copy(phi)).particle2(copy(phi));
        f.fill_cuspy_tree();
        f.print_size("cuspy tree phi(1)phi(2)");
        t1.tag("fill cuspy tree");
        real_function_6d X = (multiply(f, copy(phi), 1)).truncate();
        real_function_6d Y = g12(X);     // overwrite X to save space
        real_function_6d Z = (multiply(f, copy(phi), 1)).truncate();
        t1.tag("apply exchange");


    }

    {
        real_function_6d f = CompositeFactory<double, 6, 3>(world).g12(f12).particle1(copy(phi)).particle2(copy(one));
        f.fill_cuspy_tree();
        f.print_size("cuspy tree phi(1)one(2)");
        t1.tag("fill cuspy tree");
        real_function_6d X = (multiply(f, copy(phi), 1)).truncate();
        real_function_6d Y = g12(X);     // overwrite X to save space
        real_function_6d Z = (multiply(f, copy(phi), 1)).truncate();
        t1.tag("apply exchange");
    }
    t1.tag("fill cuspy tree");
    return  1.0;

}


int main(int argc, char **argv) {

    madness::World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    commandlineparser parser(argc, argv);
    FunctionDefaults<6>::set_tensor_type(TT_2D);
    FunctionDefaults<3>::set_thresh(1.e-5);
    FunctionDefaults<3>::set_cubic_cell(-1.0,1.0);
    FunctionDefaults<1>::set_thresh(1.e-5);
    FunctionDefaults<1>::set_cubic_cell(-10.,10.);
    FunctionDefaults<2>::set_thresh(1.e-4);
    FunctionDefaults<2>::set_cubic_cell(-10.,10.);
    FunctionDefaults<3>::set_thresh(1.e-4);
    FunctionDefaults<3>::set_cubic_cell(-10.,10.);
    FunctionDefaults<4>::set_thresh(1.e-4);
    FunctionDefaults<4>::set_cubic_cell(-10.,10.);
    FunctionDefaults<6>::set_thresh(1.e-3);
    FunctionDefaults<6>::set_cubic_cell(-10.,10.);
    print("numerical parameters: k, eps(3D), eps(6D)", FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(),
          FunctionDefaults<6>::get_thresh());
    int isuccess=0;
#ifdef USE_GENTENSOR

    try {
        parser.set_keyval("geometry", "he");
        parser.print_map();
        Molecule mol(world, parser);
        mol.print();
        CCParameters ccparam(world, parser);


//        isuccess+= test_6d_convolution(world);

//        std::shared_ptr<NuclearCorrelationFactor> ncf = create_nuclear_correlation_factor(world,
//                mol, nullptr, std::make_pair("slater", 2.0));
        isuccess+=test_lowrank_function(world);
    } catch (std::exception& e) {
        madness::print("an error occured");
        madness::print(e.what());
    }
#else
    print("could not run test_ccpairfunction: U need to compile with ENABLE_GENTENSOR=1");
#endif
    finalize();

    return isuccess;
}
