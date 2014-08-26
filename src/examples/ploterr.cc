//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>

using namespace madness;

double myf(const coord_1d& r) {
    return 1.0/sqrt(fabs(r[0])+1e-6);
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world,argc,argv);
    FunctionDefaults<1>::set_cubic_cell(-3,3);

    for (int t=1; t<6; t++) {
        double thresh = pow(0.1,t);
        FunctionDefaults<1>::set_thresh(thresh);
        for (int k=1; k<=6; k++) {
            FunctionDefaults<1>::set_k(k);

            real_function_1d f = real_factory_1d(world).f(myf);
            f.truncate();
            char fname[32];
            sprintf(fname, "err_t%1.1d_k%1.1d.dat", t, k);
            FILE* file = fopen(fname, "w");
            for (int i=0; i<=10000; i++) {
                double x = 6.0*i/10000 - 3.0;
                coord_1d r(x);
                fprintf(file, "%12.6f %12.6f %12.6f %.3e\n",
                        x, fabs(myf(r)-f(r)), myf(r), f(r));
            }
        }
    }

    finalize();
    return 0;
}
