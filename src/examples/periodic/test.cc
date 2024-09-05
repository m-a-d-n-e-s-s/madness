#include "testpc.h"
#include <fftw3.h>

/*****************

Use FFT to compute the electrostatic potential due to some test
charge densities.  Note that FFT produces a result with

  a) mean value that is zero because it projects out the constant (k=0) mode.

  b) an opposing electric field to cancel that arising from the periodic
     sum of the dipole moment potential

Test cases --- select by number
1) Gaussian spheropole
2) Gaussian dipole
3) Gaussian quadrupole
4) Cosine function

Compile with "g++ -O3 -Wall test.cc -I.. -lfftw3 -llapacke"

 *****************/

int main() {
    std::cout.precision(15);

    set_test_case(2); // select test case
    
    const size_t n = 50*L; // number of lattice points in each direction
    const size_t nh = n/2+1; // for real to complex transform last dimension
    const size_t mid = n/2;  // for mapping hermitian symmetry in transform (and testing and plotting)

    const double lo = -0.5*L;
    const double hi = 0.5*L;
    const double h = (hi - lo) / (n - 1);
    std::vector<double> x = linspace(lo,hi,n,false); // exclude right endpoint
    std::vector<double> y = x; // Presently assuming cubic domain with equal spacing in x, y, and z
    std::vector<double> z = x;

    std::vector<double> F = tabulate(f, x, y, z);
    {
        madness::Gnuplot g("set style data l");
        // Extract the middle z column?
        std::vector<double> col(&F[mid*n*n + mid*n], &F[mid*n*n + mid*n + n]);
        g.plot(z, col);
    }

    // sum the values in F to get the total charge
    double norm = std::reduce(F.begin(), F.end(), 0.0)*h*h*h;
    std::cout << "norm = " << norm << std::endl;
    
    // Perform a 3D Fourier transform of F and store the result in G
    std::vector<fftw_complex> G(n*n*nh);
    {
        fftw_plan plan = fftw_plan_dft_r2c_3d(n, n, n, F.data() , G.data(), FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    }

    // Navigate the way FFTW stores the Hermitian symmetric data
    auto fudge = [=](size_t j) {return j<nh ? j : n-j;};

    double kscale = 4.0*pi*L*L/(n*n);

    // Apply Coulomb GF in fourier space
    for (size_t jx=0; jx<n; jx++) {
        double kx = 2.0*pi*fudge(jx)/n;
        for (size_t jy=0; jy<n; jy++) {
            double ky = 2.0*pi*fudge(jy)/n;
            for (size_t jz=0; jz<nh; jz++) {
                double kz = 2.0*pi*jz/n;
                double ksq = kx*kx + ky*ky + kz*kz;

                if (ksq > 1.0e-20) {
                    double kfac = kscale/ksq;
                    size_t j = jx*n*nh + jy*nh + jz;
                    //if (std::abs(G[j][0]) > 1e-6) {
                        //std::cout << jx << " (" << fudge(jx) << ") " << jy << " (" << fudge(jy) << ") " << jz << " " << ksq << " " << G[j][0] << " " << G[j][1] << std::endl;
                    //}
                    G[j][0] *= kfac;
                    G[j][1] *= kfac;
                }
            }
        }
    }
    
    // Do the reverse transform into FF
    std::vector<double> FF(n*n*n);
    {
        fftw_plan plan = fftw_plan_dft_c2r_3d(n, n, n, G.data() , FF.data(), FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        const double scale = 1.0/(n*n*n); // FFTW transforms are not normalized
        for (double& u : FF) u *= scale;
    }

    // Check the mean potential --- expected to be zero
    double FFmean = std::reduce(FF.begin(), FF.end(), 0.0)*h*h*h / (L*L*L);
    std::cout << "FFmean = " << FFmean << std::endl;
    
    size_t j = mid; //mid+3; // test point (avoid origin) -- don't need to avoid origin now we are treating erf(0)/0 correctly
    {
        madness::Gnuplot g("set style data l");
        std::vector<double> col(&FF[j*n*n + j*n], &FF[j*n*n + j*n + n]);
        std::vector<double> col2(n);
        std::vector<double> col3(n);
        std::cout << " i   z        FFT               analytic       FFT-analytic\n";
        
        for (size_t i = 0; i < n; i++) {
            double v = exact(x[j],y[j],z[i]);
            std::cout << i << " " << z[i] << " " << col[i] << " " << v << " " << col[i]-v << std::endl;
            col2[i] = v;
            col3[i] = col[i] - v;
        }

        //g.plot(z, col, col2);
        g.plot(z, col3);
        
    }

    // std::vector<double> sums = periodic_sum_partial(x[j],y[j],z[0],exact_quadrupoleX,50);
    // //std::vector<double> sums = periodic_sum_partial(x[j],y[j],z[0],exact_dipoleX,50);
    // for (size_t i = 0; i < sums.size(); i++) {
    //     std::cout << i << " " << sums[i] << " " << exact(x[j],y[j],z[0]) << std::endl;
    // }
    // std::cout << exact(x[j],y[j],z[0]) << " " << FF[j*n*n + j*n + 0] <<  " !! " << periodic_sum_accelerated(x[j],y[j],z[0],exact_quadrupoleX) << std::endl;

    // size_t p = 7; // no. of terms in the polynomial fit --- 7 is best for quadrupole
    // for (size_t i=p+1; i<sums.size()-p; i++) {
    //     // Extract the last p terms of the series in v starting at position i
    //     std::vector<double> f(&sums[i], &sums[i+p]);
    //     std::vector<double> N(p);
    //     for (size_t j=0; j<p; j++) {
    //         N[j] = i+j;
    //     }
    //     std::vector<double> c = fit(p, p, N, f);
    //     std::cout << i << " " << c[0] << " " << sums[i]-FF[j*n*n + j*n + 0] << " " << c[0]-FF[j*n*n + j*n + 0] << std::endl;
    // }

    return 0;
}
