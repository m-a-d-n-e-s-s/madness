//
// Created by Florian Bischoff on 4/11/24.
//

#include <madness.h>
#include<madness/world/test_utilities.h>
#include<madness/chem/IntegratorXX.h>

using namespace madness;


// \int s*s dV = (pi/2)^(3/2)
double s_func(const Vector<double,3>& r) { return exp(-inner(r,r)); }
// \int px*px dV = pi^(3/2)/(8 Sqrt(2))
double px_func(const Vector<double,3>& r) { return r[0]*exp(-inner(r,r)); }
double py_func(const Vector<double,3>& r) { return r[1]*exp(-inner(r,r)); }
double pz_func(const Vector<double,3>& r) { return r[2]*exp(-inner(r,r)); }
double dxy_func(const Vector<double,3>& r) { return r[0]*r[1]*exp(-inner(r,r)); }
double dxz_func(const Vector<double,3>& r) { return r[0]*r[2]*exp(-inner(r,r)); }
double dyz_func(const Vector<double,3>& r) { return r[1]*r[2]*exp(-inner(r,r)); }
double dx2y2_func(const Vector<double,3>& r) { return (r[0]*r[0]-r[1]*r[1])*exp(-inner(r,r)); }
double dz2_func(const Vector<double,3>& r) { return (r[2]*r[2]-0.5)*exp(-inner(r,r)); }

struct spherical_harmonic {
    int l,m;
    Vector<double,3> origin;
    spherical_harmonic(int l1, int m1, Vector<double,3> origin1=Vector<double,3>()) : l(l1), m(m1), origin(origin1) {
        if (l==0) func=s_func;
        else if (l==1) {
            if (m==-1) func=px_func;
            else if (m==0) func=py_func;
            else if (m==1) func=pz_func;
        } else if (l==2) {
            if (m==-2) func=dxy_func;
            else if (m==-1) func=dxz_func;
            else if (m==0) func=dyz_func;
            else if (m==1) func=dx2y2_func;
            else if (m==2) func=dz2_func;
        } else {
            MADNESS_EXCEPTION("unsupported spherical harmonic",1);
        }
    }

    double operator()(const Vector<double,3>& r) const {
        return func(r-origin);
    }

    std::function<double(const Vector<double,3>&)> func;
};



int test_construction() {
    test_output t("DFT grid construction");
    GridBuilder builder;
    builder.make_grid();
    auto points=builder.get_points();
    return t.end();
}

int test_integration(World& world) {

    test_output t("DFT grid integration");

    // center the grid off the origin
    Vector<double,3> gridcenter={1,0,-0.5};


    auto error=[&world, &gridcenter](const int l, const int m, const auto points, const auto weights) {
        auto f=spherical_harmonic(l,m,gridcenter);
        double integral=0.0;
        for (size_t i=0; i<points.size(); i++) {
            double r2=inner(points[i],points[i]);
            integral+=weights[i]*f(points[i])*f(points[i]);
        }
        // print("integral by grid",integral);

        real_function_3d ff=real_factory_3d(world).functor(f);
        double integral2=ff.norm2();
        // print("integral by MRA ",integral2*integral2);
        double exact=0.0;
        if (l==0) exact=std::pow(M_PI*0.5,1.5);
        else if (l==1) exact=std::pow(M_PI,1.5)/8.0/std::sqrt(2.0);
        else exact=integral2*integral2;
        return std::abs(integral-exact)/std::abs(exact);
    };

    // expected error(nrad,angular order,l)
    std::map<std::tuple<int,int,int>,double> errors;
    errors[std::tuple(10l,3l,0l)]=1.0e-2;
    errors[std::tuple(30,3,0)]=1.0e-4;
    errors[std::tuple(50,3,0)]=1.0e-6;

    errors[std::tuple(10,3,1)]=1.0e-1;
    errors[std::tuple(30,3,1)]=1.0e-4;
    errors[std::tuple(50,3,1)]=1.0e-5;

    errors[std::tuple(10,3,2)]=1.0e1;       // l=2 will fail by construction
    errors[std::tuple(30,3,2)]=1.0e1;       // l=2 will fail by construction
    errors[std::tuple(50,3,2)]=1.0e1;       // l=2 will fail by construction

    errors[std::tuple(10,4,0)]=1.0e-2;
    errors[std::tuple(30,4,0)]=1.0e-4;
    errors[std::tuple(50,4,0)]=1.0e-6;

    errors[std::tuple(10,4,1)]=1.0e-1;
    errors[std::tuple(30,4,1)]=1.0e-4;
    errors[std::tuple(50,4,1)]=1.0e-5;

    errors[std::tuple(10,4,2)]=2.0e-1;
    errors[std::tuple(30,4,2)]=2.0e-3;
    errors[std::tuple(50,4,2)]=1.0e-4;

    for (std::size_t order=3; order<5; ++order) {
        for (std::size_t nrad=10; nrad<61; nrad+=20) {
            GridBuilder builder;
            builder.set_nradial(nrad);
            builder.set_origin(gridcenter);
            builder.set_angular_order(order);
            builder.make_grid();
            auto points=builder.get_points();
            auto weights=builder.get_weights();
            for (int l=0; l<3; ++l) {
                for (int m=-l; m<l+1; ++m) {
                    double err=error(l,m,points,weights);
                    std::stringstream ss;
                    ss << std::setprecision(1) << std::scientific <<  err;
                    print("error, l,m, nrad, angular_order",l,m,builder.get_nradial(),builder.get_angular_order(),err);
                    std::string msg="nrad, ang order, l, rel. error: " +std::to_string(nrad)
                        + " " + std::to_string(builder.get_angular_order())
                        + " " + std::to_string(l)
                        + " " + ss.str();
                    t.checkpoint(err<errors[std::tie(nrad,order,l)],msg);
                }
            }
        }
    }
    return t.end();
}



int main(int argc, char **argv) {

    madness::World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    FunctionDefaults<3>::set_cubic_cell(-10,10);
    FunctionDefaults<3>::set_thresh(1.e-6);
    FunctionDefaults<3>::set_k(7);

    int success=0;
    success+=test_construction();
    success+=test_integration(world);

    world.gop.fence();
    madness::finalize();

    return success;
}
