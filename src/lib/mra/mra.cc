#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/loadbal.h>

/// \file mra.cc
/// \file Declaration and initialization of static data and some implementation

namespace madness {

    // Definition and initialization of FunctionDefaults static members
    // It cannot be an instance of FunctionFactory since we want to
    // set the defaults independent of the data type.  
    
    template <int NDIM> int FunctionDefaults<NDIM>::k = 7;
    template <int NDIM> double FunctionDefaults<NDIM>::thresh = 1e-5;
    template <int NDIM> int FunctionDefaults<NDIM>::initial_level = 2;
    template <int NDIM> int FunctionDefaults<NDIM>::max_refine_level = 30;
    template <int NDIM> int FunctionDefaults<NDIM>::truncate_method = 0;
    template <int NDIM> bool FunctionDefaults<NDIM>::compress = true;
    template <int NDIM> bool FunctionDefaults<NDIM>::refine = true;
    template <int NDIM> bool FunctionDefaults<NDIM>::autorefine = false;
    template <int NDIM> bool FunctionDefaults<NDIM>::debug = false;

    template <> int FunctionDefaults<1>::bc[1][2] = {{0,0}};
    template <> double FunctionDefaults<1>::cell[1][2] = {{0.0,1.0}};
    
    template <> int FunctionDefaults<2>::bc[2][2] = {{0,0}, {0,0}};
    template <> double FunctionDefaults<2>::cell[2][2] = {{0.0,1.0}, {0.0,1.0}};
    
    template <> int FunctionDefaults<3>::bc[3][2] = {{0,0}, {0,0}, {0,0}};
    template <> double FunctionDefaults<3>::cell[3][2] = {{0.0,1.0}, {0.0,1.0}, {0.0,1.0}};
    
    template <> int FunctionDefaults<4>::bc[4][2] = {{0,0}, {0,0}, {0,0},{0,0}};
    template <> double FunctionDefaults<4>::cell[4][2] = {{0.0,1.0}, {0.0,1.0}, {0.0,1.0}, {0.0,1.0}};
    
    template <> int FunctionDefaults<5>::bc[5][2] = {{0,0}, {0,0}, {0,0}, {0,0}, {0,0}};
    template <> double FunctionDefaults<5>::cell[5][2] = {{0.0,1.0}, {0.0,1.0}, {0.0,1.0}, {0.0,1.0}, {0.0,1.0}};
    
    template <> int FunctionDefaults<6>::bc[6][2] = {{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}};
    template <> double FunctionDefaults<6>::cell[6][2] = {{0.0,1.0}, {0.0,1.0}, {0.0,1.0}, {0.0,1.0}, {0.0,1.0}, {0.0,1.0}};
    
    template class FunctionDefaults<1>;
    template class FunctionDefaults<2>;
    template class FunctionDefaults<3>;
    template class FunctionDefaults<4>;
    template class FunctionDefaults<5>;
    template class FunctionDefaults<6>;
    
    template <typename T, int NDIM, typename Pmap> FunctionCommonData<T,NDIM> FunctionImpl<T,NDIM,Pmap>::commondata[FunctionImpl<T,NDIM,Pmap>::MAXK+1];
    template <typename T, int NDIM, typename Pmap> bool FunctionImpl<T,NDIM,Pmap>::initialized = false;

    template <typename T, int NDIM>
    void FunctionCommonData<T,NDIM>::_make_dc_periodic() {
        // See ABGV for details
        r0 = Tensor<double>(k,k);
        rp = Tensor<double>(k,k);
        rm = Tensor<double>(k,k);

        double iphase = 1.0;
        for (int i=0; i<k; i++) {
            double jphase = 1.0;
            for (int j=0; j<k; j++) {
                double gammaij = sqrt(double((2*i+1)*(2*j+1)));
                double Kij;
                if (((i-j)>0) && (((i-j)%2)==1))
                    Kij = 2.0;
                else
                    Kij = 0.0;

                r0(i,j) = 0.5*(1.0 - iphase*jphase - 2.0*Kij)*gammaij;
                rm(i,j) = 0.5*jphase*gammaij;
                rp(i,j) =-0.5*iphase*gammaij;
                jphase = -jphase;
            }
            iphase = -iphase;
        }

        // Make the rank-1 forms of rm and rp
        rm_left = Tensor<double>(k);
        rm_right = Tensor<double>(k);
        rp_left = Tensor<double>(k);
        rp_right = Tensor<double>(k);

        iphase = 1.0;
        for (int i=0; i<k; i++) {
            double gamma = sqrt(0.5*(2*i+1));
            rm_left(i)  = rp_right(i) = gamma;
            rm_right(i) = rp_left(i)  = gamma*iphase;
            iphase *= -1.0;
        }
        rp_left.scale(-1.0);

//         Tensor<double> rm_test = outer(rm_left,rm_right);
//         Tensor<double> rp_test = outer(rp_left,rp_right);
    }

    template <typename T, int NDIM>
    void FunctionCommonData<T,NDIM>::_init_twoscale() {
        if (! two_scale_hg(k, &hg)) throw "failed to get twoscale coefficients";
        hgT = transpose(hg);
        hgsonly = copy(hg(Slice(0,k-1),_));
    }

    template <typename T, int NDIM>
    void FunctionCommonData<T,NDIM>::_init_quadrature() {
        quad_x = Tensor<double>(npt);
        quad_w = Tensor<double>(npt);
        quad_phi = Tensor<double>(npt,k);
        quad_phiw = Tensor<double>(npt,k);

        gauss_legendre(npt,0.0,1.0,quad_x.ptr(),quad_w.ptr());
        for (int i=0; i<npt; i++) {
            double phi[200];
            legendre_scaling_functions(quad_x(i),k,phi);
            for (int j=0; j<k; j++) {
                quad_phi(i,j) = phi[j];
                quad_phiw(i,j) = quad_w(i)*phi[j];
            }
        }
        quad_phit = transpose(quad_phi);
    }


    template class Function<double, 1>;
    template class Function<std::complex<double>, 1>;
    template class FunctionImpl<double, 1>;
    template class FunctionImpl<std::complex<double>, 1>;
    template class Function<double, 1, DClass<1>::MyProcMap>;
    template class Function<std::complex<double>, 1, DClass<1>::MyProcMap>;
    template class FunctionImpl<double, 1, DClass<1>::MyProcMap>;
    template class FunctionImpl<std::complex<double>, 1, DClass<1>::MyProcMap>;
    template class FunctionCommonData<double, 1>;
    template class FunctionCommonData<double_complex, 1>;

    template class Function<double, 2>;
    template class Function<std::complex<double>, 2>;
    template class FunctionImpl<double, 2>;
    template class FunctionImpl<std::complex<double>, 2>;
    template class Function<double, 2, DClass<2>::MyProcMap>;
    template class Function<std::complex<double>, 2, DClass<2>::MyProcMap>;
    template class FunctionImpl<double, 2, DClass<2>::MyProcMap>;
    template class FunctionImpl<std::complex<double>, 2, DClass<2>::MyProcMap>;
    template class FunctionCommonData<double, 2>;
    template class FunctionCommonData<double_complex, 2>;

    template class Function<double, 3>;
    template class Function<std::complex<double>, 3>;
    template class FunctionImpl<double, 3>;
    template class FunctionImpl<std::complex<double>, 3>;
    template class Function<double, 3, DClass<3>::MyProcMap>;
    template class Function<std::complex<double>, 3, DClass<3>::MyProcMap>;
    template class FunctionImpl<double, 3, DClass<3>::MyProcMap>;
    template class FunctionImpl<std::complex<double>, 3, DClass<3>::MyProcMap>;
    template class FunctionCommonData<double, 3>;
    template class FunctionCommonData<double_complex, 3>;

    template class Function<double, 4>;
    template class Function<std::complex<double>, 4>;
    template class FunctionImpl<double, 4>;
    template class FunctionImpl<std::complex<double>, 4>;
    template class Function<double, 4, DClass<4>::MyProcMap>;
    template class Function<std::complex<double>, 4, DClass<4>::MyProcMap>;
    template class FunctionImpl<double, 4, DClass<4>::MyProcMap>;
    template class FunctionImpl<std::complex<double>, 4, DClass<4>::MyProcMap>;
    template class FunctionCommonData<double, 4>;
    template class FunctionCommonData<double_complex, 4>;

    template class Function<double, 5>;
    template class Function<std::complex<double>, 5>;
    template class FunctionImpl<double, 5>;
    template class FunctionImpl<std::complex<double>, 5>;
    template class Function<double, 5, DClass<5>::MyProcMap>;
    template class Function<std::complex<double>, 5, DClass<5>::MyProcMap>;
    template class FunctionImpl<double, 5, DClass<5>::MyProcMap>;
    template class FunctionImpl<std::complex<double>, 5, DClass<5>::MyProcMap>;
    template class FunctionCommonData<double, 5>;
    template class FunctionCommonData<double_complex, 5>;

    template class Function<double, 6>;
    template class Function<std::complex<double>, 6>;
    template class FunctionImpl<double, 6>;
    template class FunctionImpl<std::complex<double>, 6>;
    template class Function<double, 6, DClass<6>::MyProcMap>;
    template class Function<std::complex<double>, 6, DClass<6>::MyProcMap>;
    template class FunctionImpl<double, 6, DClass<6>::MyProcMap>;
    template class FunctionImpl<std::complex<double>, 6, DClass<6>::MyProcMap>;
    template class FunctionCommonData<double, 6>;
    template class FunctionCommonData<double_complex, 6>;

}
