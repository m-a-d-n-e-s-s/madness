#include "input.h"
#include "densities.h"
#include "loadbalance.h"
#include "states.h"
#include "potential.h"
#include "energies.h"
#include "output.h"

double tt, ss;

#define START_TIMER world.gop.fence(); tt=wall_time(); ss=cpu_time()
#define END_TIMER(msg) tt=wall_time()-tt; ss=cpu_time()-ss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, ss, tt)

// Calculates effective-mass potential 
void Umeff(World& world, comp_vecfuncT& Upsi_qu, 
                         comp_vecfuncT& Upsi_qd, 
                         const real_functionT& B_q,
                         const comp_vecfuncT& dpsi_qu_dx, 
                         const comp_vecfuncT& dpsi_qu_dy, 
                         const comp_vecfuncT& dpsi_qu_dz,
                         const comp_vecfuncT& dpsi_qd_dx, 
                         const comp_vecfuncT& dpsi_qd_dy, 
                         const comp_vecfuncT& dpsi_qd_dz,
                         const double prec)

{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;
    const double_complex  one( 1.0,0.0);
    const double_complex mone(-1.0,0.0);

    B_q.reconstruct();
    world.gop.fence();

    real_functionT B_q1 = 1.0 * B_q;
    real_functionT B_q2 = 1.0 * B_q;
    real_functionT B_q3 = 1.0 * B_q;
    world.gop.fence();

    compress(world,  Upsi_qu, false);
    compress(world,  Upsi_qd, false);
    world.gop.fence();

    complex_derivative_3d Dx(world, 0);
    complex_derivative_3d Dy(world, 1);
    complex_derivative_3d Dz(world, 2);

    comp_vecfuncT mul1 = mul(world, B_q1, dpsi_qu_dx, false);
    comp_vecfuncT mul2 = mul(world, B_q2, dpsi_qu_dy, false);
    comp_vecfuncT mul3 = mul(world, B_q3, dpsi_qu_dz, false);
    world.gop.fence();

    comp_vecfuncT Ueff1psi_qu = apply(world, Dx, mul1, false);
    comp_vecfuncT Ueff2psi_qu = apply(world, Dy, mul2, false);
    comp_vecfuncT Ueff3psi_qu = apply(world, Dz, mul3, false);
    world.gop.fence();

    mul1.clear(), mul2.clear(), mul3.clear();
    world.gop.fence();

    truncate2(world, Ueff1psi_qu, Ueff2psi_qu, prec);
    world.gop.fence();
    truncate(world,  Ueff3psi_qu, ntol);
    world.gop.fence();

    compress(world, Ueff1psi_qu, false);
    compress(world, Ueff2psi_qu, false);
    compress(world, Ueff3psi_qu, false);
    world.gop.fence();

    gaxpy(world, one, Upsi_qu, mone, Ueff1psi_qu, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qu, mone, Ueff2psi_qu, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qu, mone, Ueff3psi_qu, false);

    world.gop.fence();
    Ueff1psi_qu.clear(); Ueff2psi_qu.clear(); Ueff3psi_qu.clear();
    world.gop.fence();

    comp_vecfuncT mul4 = mul(world, B_q1, dpsi_qd_dx, false);
    comp_vecfuncT mul5 = mul(world, B_q2, dpsi_qd_dy, false);
    comp_vecfuncT mul6 = mul(world, B_q3, dpsi_qd_dz, false);
    world.gop.fence();

    comp_vecfuncT Ueff1psi_qd = apply(world, Dx, mul4, false);
    comp_vecfuncT Ueff2psi_qd = apply(world, Dy, mul5, false);
    comp_vecfuncT Ueff3psi_qd = apply(world, Dz, mul6, false);

    world.gop.fence();
    mul4.clear(), mul5.clear(), mul6.clear();

    truncate2(world, Ueff1psi_qd, Ueff2psi_qd, prec);
    world.gop.fence();
    truncate(world,  Ueff3psi_qd, ntol);
    world.gop.fence();

    compress(world, Ueff1psi_qd, false);
    compress(world, Ueff2psi_qd, false);
    compress(world, Ueff3psi_qd, false);
    world.gop.fence();

    gaxpy(world, one, Upsi_qd, mone, Ueff1psi_qd, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qd, mone, Ueff2psi_qd, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qd, mone, Ueff3psi_qd, false);

    world.gop.fence();
    Ueff1psi_qd.clear(), Ueff2psi_qd.clear(), Ueff3psi_qd.clear();
    B_q1.clear(), B_q2.clear(); B_q3.clear();
}



// Calculates spin-orbit potential
void Uso(World& world, comp_vecfuncT& Upsi_qu,
                       comp_vecfuncT& Upsi_qd,
                       const real_functionT& W_q,
                       const comp_vecfuncT& dpsi_qu_dx,
                       const comp_vecfuncT& dpsi_qd_dx,
                       const comp_vecfuncT& dpsi_qu_dy,
                       const comp_vecfuncT& dpsi_qd_dy,
                       const comp_vecfuncT& dpsi_qu_dz,
                       const comp_vecfuncT& dpsi_qd_dz)
{
    const double_complex one(1.0, 0.0);
    const double lambda = -1.0;
    const double_complex lam(lambda, 0.0);
    //const double_complex I = double_complex(0.0, 1.0);

    W_q.reconstruct();
    world.gop.fence();

    real_functionT dW_q_dx = real_derivative_3d(world, 0)(W_q);
    world.gop.fence();
    real_functionT dW_q_dy = real_derivative_3d(world, 1)(W_q);
    world.gop.fence();
    real_functionT dW_q_dz = real_derivative_3d(world, 2)(W_q);

    compress(world, Upsi_qu, false);
    compress(world, Upsi_qd, false);
    world.gop.fence();

    comp_vecfuncT mul1 = mul(world, dW_q_dy, dpsi_qu_dx, false);
    world.gop.fence();
    comp_vecfuncT mul2 = mul(world, dW_q_dy, dpsi_qd_dx, false);
    world.gop.fence();
    comp_vecfuncT mul3 = mul(world, dW_q_dz, dpsi_qu_dx, false);
    world.gop.fence();
    comp_vecfuncT mul4 = mul(world, dW_q_dz, dpsi_qd_dx, false);
    world.gop.fence();

    compress(world, mul1, false);
    compress(world, mul2, false);
    compress(world, mul3, false);
    compress(world, mul4, false);
    world.gop.fence();

    gaxpy(world, one, Upsi_qu, -I*lam, mul1, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qd,  I*lam, mul2, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qd,   -lam, mul3, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qu,    lam, mul4, false);
    world.gop.fence();

    mul1.clear(), mul2.clear(), mul3.clear(), mul4.clear();
    world.gop.fence();

    comp_vecfuncT mul5 = mul(world, dW_q_dx, dpsi_qu_dy, false);
    world.gop.fence();
    comp_vecfuncT mul6 = mul(world, dW_q_dx, dpsi_qd_dy, false);
    world.gop.fence();
    comp_vecfuncT mul7 = mul(world, dW_q_dz, dpsi_qu_dy, false);
    world.gop.fence();
    comp_vecfuncT mul8 = mul(world, dW_q_dz, dpsi_qd_dy, false);
    world.gop.fence();

    compress(world, mul5, false);
    compress(world, mul6, false);
    compress(world, mul7, false);
    compress(world, mul8, false);
    world.gop.fence();

    gaxpy(world, one, Upsi_qu,  I*lam, mul5, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qd, -I*lam, mul6, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qd, -I*lam, mul7, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qu, -I*lam, mul8, false);
    world.gop.fence();

    mul5.clear(), mul6.clear(), mul7.clear(), mul8.clear();
    world.gop.fence();

    comp_vecfuncT mul9  = mul(world, dW_q_dx, dpsi_qd_dz, false);
    world.gop.fence();
    comp_vecfuncT mul10 = mul(world, dW_q_dx, dpsi_qu_dz, false);
    world.gop.fence();
    comp_vecfuncT mul11 = mul(world, dW_q_dy, dpsi_qd_dz, false);
    world.gop.fence();
    comp_vecfuncT mul12 = mul(world, dW_q_dy, dpsi_qu_dz, false);
    world.gop.fence();

    gaxpy(world, one, Upsi_qu,   -lam, mul9,  false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qd,    lam, mul10, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qu,  I*lam, mul11, false);
    world.gop.fence();
    gaxpy(world, one, Upsi_qd,  I*lam, mul12, false);
    world.gop.fence();

    mul9.clear(), mul10.clear(), mul11.clear(), mul12.clear();
}


// Makes Skyrme Potential
void Potential(World& world, const comp_vecfuncT& psi_pu,
                             const comp_vecfuncT& psi_pd,
                             const comp_vecfuncT& psi_nu,
                             const comp_vecfuncT& psi_nd,
                             real_tensorT& energy_p,
                             real_tensorT& energy_n,
                             real_functionT& U_p,
                             real_functionT& U_n,
                             real_functionT& rho_p,
                             real_functionT& rho_n,
                             real_functionT& rho,
                             const int iter,
                             double& BE,
                             const double delta_psi,
                             const double tol,
                             const double brad,
                             real_functionT& U_p_old,
                             real_functionT& U_n_old,
                             real_functionT& lap_p_old,
                             real_functionT& lap_n_old,
                             const double L,
                             const int jellium,
                             const int spinorbit,
                             const int screening,
                             const double screenl,
                             const double avg_pot,
                             const double avg_lap,
                             const int lap_comp,
                             const double prec,
                             const int vtk_output,
                             const int txt_output,
                             const real_tensorT b,
                             const real_tensorT bp,
                             const double alpha,
                             const double k_fn,
                             const double k_fp,
                             const int timing)
{

    double ntol = FunctionDefaults<3>::get_thresh() * prec;
    int np = psi_pu.size();
    int nn = psi_nu.size();

    /////////////////////////////////////////////////////////////////////////////
    /////////                         Densities                          ////////
    /////////////////////////////////////////////////////////////////////////////

    // Make number density and kinetic densities
    if (world.rank() == 0) {print("   Make densities");}
    if (BE == 0.0) {
        if (timing == 1) {START_TIMER;}
        rho_p = ndensity(world, psi_pu, psi_pd);
        rho_n = ndensity(world, psi_nu, psi_nd);
        if (timing == 1) {END_TIMER("ndensity");}
    }

    double p_number = rho_p.trace();
    double n_number = rho_n.trace();
    if (world.rank() == 0) {
        print("   Protons  =", p_number);
        print("   Neutrons =", n_number);
    }
    world.gop.fence();

    real_functionT dJ_p  = 0.0 * rho_p;
    real_functionT dJ_n  = 0.0 * rho_n;
    comp_tensorT K_p(np,np);
    comp_tensorT K_n(nn,nn);

    complex_derivative_3d Dx(world, 0);
    complex_derivative_3d Dy(world, 1);
    complex_derivative_3d Dz(world, 2);
    world.gop.fence();

    comp_vecfuncT dpsi_pu_dx = apply(world, Dx, psi_pu, false);
    comp_vecfuncT dpsi_pd_dy = apply(world, Dy, psi_pd, false);
    world.gop.fence();

    comp_vecfuncT dpsi_pd_dx = apply(world, Dx, psi_pd, false);
    comp_vecfuncT dpsi_pu_dz = apply(world, Dz, psi_pu, false);
    world.gop.fence();

    comp_vecfuncT dpsi_pu_dy = apply(world, Dy, psi_pu, false);
    comp_vecfuncT dpsi_pd_dz = apply(world, Dz, psi_pd, false);
    world.gop.fence();

    if (world.rank() == 0) {print("   Kdensities ");}
    if (timing == 1) {START_TIMER;}
    real_functionT tau_p = kdensity(world, dpsi_pu_dx, dpsi_pd_dx, dpsi_pu_dy, dpsi_pd_dy, dpsi_pu_dz, dpsi_pd_dz);
    if (timing == 1) {END_TIMER("kdensity");}

    real_functionT lap_p = 0.0 * rho_p;
    world.gop.fence();

    // Gaussian smoothing 
    const coordT origin(0.0);
    Tensor<double> exponent(1), coeff(1);
    exponent = 1.0/(2.0 * brad * brad);
    coeff = pow(1.0/(brad * sqrt(2.0 * M_PI)), 3.0);
    operatorT G(world, coeff, exponent);

    if (timing == 1) {START_TIMER;}
    if (lap_comp == 0) {
        if (world.rank() == 0) {print("   Laplacian method 0");}
        real_functionT rho_pt = G(rho_p);
        lap_p = laplacian(world, rho_p, brad, prec);
        world.gop.fence();
    }
    else if (lap_comp == 1) {
        if (world.rank() == 0) {print("   Laplacian method 1");}
        lap_p = laplacian1(world, dpsi_pu_dx, dpsi_pd_dx, dpsi_pu_dy, dpsi_pd_dy, dpsi_pu_dz, dpsi_pd_dz, psi_pu, psi_pd);
    }
    else if (lap_comp == 2) {
        if (world.rank() == 0) {print("   Laplacian method 2");}
        lap_p = laplacian2(world, dpsi_pu_dx, dpsi_pd_dx, dpsi_pu_dy, dpsi_pd_dy, dpsi_pu_dz, dpsi_pd_dz, psi_pu, psi_pd, prec);
    }
    else {if (world.rank() == 0) {print("   No method to calculate laplacian selected");}}
    if (timing == 1) {END_TIMER("Laplacian for p");}
    world.gop.fence();

    if(spinorbit != 0) {
        if (timing == 1) {START_TIMER;}
        if (world.rank() == 0) {print("   Spin-orbit densities ");}
        dJ_p  = so_density(world, dpsi_pu_dx, dpsi_pd_dx, dpsi_pu_dy, dpsi_pd_dy, dpsi_pu_dz, dpsi_pd_dz, prec);
        if (timing == 1) {END_TIMER("so_density");}
        world.gop.fence();
    }

    // Matrix of proton single particle kinetic energies 
    if (iter%100 == 0 || iter == 0) {
        if (world.rank() == 0) {print("   Kmatrix for protons ");}
        K_p = Kmatrix(world, psi_pu, psi_pd, dpsi_pu_dx, dpsi_pd_dx, dpsi_pu_dy, dpsi_pd_dy, dpsi_pu_dz, dpsi_pd_dz, prec, k_fp);
        world.gop.fence();
    }
    if (world.rank() == 0) {print("   done ");}

    dpsi_pu_dx.clear(); dpsi_pd_dx.clear();
    dpsi_pu_dy.clear(); dpsi_pd_dy.clear();
    dpsi_pu_dz.clear(); dpsi_pd_dz.clear();
    world.gop.fence();

    comp_vecfuncT dpsi_nu_dx = apply(world, Dx, psi_nu, false);
    comp_vecfuncT dpsi_nd_dy = apply(world, Dy, psi_nd, false);
    world.gop.fence();

    comp_vecfuncT dpsi_nd_dx = apply(world, Dx, psi_nd, false);
    comp_vecfuncT dpsi_nu_dz = apply(world, Dz, psi_nu, false);
    world.gop.fence();

    comp_vecfuncT dpsi_nu_dy = apply(world, Dy, psi_nu, false);
    comp_vecfuncT dpsi_nd_dz = apply(world, Dz, psi_nd, false);
    world.gop.fence();

        if (timing == 1) {START_TIMER;}
    real_functionT tau_n = kdensity(world, dpsi_nu_dx, dpsi_nd_dx, dpsi_nu_dy, dpsi_nd_dy, dpsi_nu_dz, dpsi_nd_dz);
    real_functionT lap_n = 0.0 * rho_n;
    world.gop.fence();

    if (lap_comp == 0) {
        if (world.rank() == 0) {print("   Laplacian method 0");}
        real_functionT rho_nt = G(rho_n);
        lap_n = laplacian(world, rho_nt, brad, prec);
    }
    else if (lap_comp == 1) {
        if (world.rank() == 0) {print("   Laplacian method 1");}
        lap_n = laplacian1(world, dpsi_nu_dx, dpsi_nd_dx, dpsi_nu_dy, dpsi_nd_dy, dpsi_nu_dz, dpsi_nd_dz, psi_nu, psi_nd);
    }
    else if (lap_comp == 2) {
        if (world.rank() == 0) {print("   Laplacian method 2");}
        lap_n = laplacian2(world, dpsi_nu_dx, dpsi_nd_dx, dpsi_nu_dy, dpsi_nd_dy, dpsi_nu_dz, dpsi_nd_dz, psi_nu, psi_nd, prec);
    }
    else {
        if (world.rank() == 0) {print("   Error: No method to calculate density laplacian selected");}
        assert(0);
    }
    if (timing == 1) {END_TIMER("Laplacian for n");}
    world.gop.fence();

    if (spinorbit != 0) {
        if (world.rank() == 0) {print("   Spin-orbit densities ");}
        dJ_n  = so_density(world, dpsi_nu_dx, dpsi_nd_dx, dpsi_nu_dy, dpsi_nd_dy, dpsi_nu_dz, dpsi_nd_dz, prec);
        world.gop.fence();
    }

    // Matrix of neutron single particle kinetic energies 
    if (iter%100 == 0 || iter == 0) {
        K_n = Kmatrix(world, psi_nu, psi_nd, dpsi_nu_dx, dpsi_nd_dx, dpsi_nu_dy, dpsi_nd_dy, dpsi_nu_dz, dpsi_nd_dz, prec, k_fn);
        world.gop.fence();
    }

    dpsi_nu_dx.clear(); dpsi_nd_dx.clear();
    dpsi_nu_dy.clear(); dpsi_nd_dy.clear();
    dpsi_nu_dz.clear(); dpsi_nd_dz.clear();
    world.gop.fence();

    // Mixing of density laplacians
    if (avg_lap == 1 && brad < 0.0 && iter > 0) {
        if (world.rank() == 0) {print("   Mix laplacians");}
        double lap_avg = 0.4;
        lap_p = lap_avg * lap_p + (1.0 - lap_avg) * lap_p_old;
        lap_n = lap_avg * lap_n + (1.0 - lap_avg) * lap_n_old;
    }

    // Make total densities
    rho = rho_p + rho_n;
    real_functionT tau = tau_p + tau_n;
    real_functionT lap = lap_p + lap_n;

    lap_p_old = copy(lap_p);
    lap_n_old = copy(lap_n);
    world.gop.fence();

    if (iter >= 0 && brad > 0.0) {
        if (world.rank() == 0) print("   Gaussian smoothing");
        lap_p = G(lap_p); lap_n = G(lap_n); lap = G(lap);
        tau_p = G(tau_p); tau_n = G(tau_n); tau = G(tau);
        world.gop.fence();
    }

    real_functionT power = real_factoryT(world);
    power = 0.0 * power + alpha;
    real_functionT rho_1   = binary_op(rho, power + 1.0, rho_power());

    real_functionT rho_2_n = rho_n * binary_op(rho, power, rho_power());
    real_functionT rho_2_p = rho_p * binary_op(rho, power, rho_power());

    real_functionT rho_3_n = rho_n * rho_n * binary_op(rho, power - 1.0, rho_power());
    real_functionT rho_3_p = rho_p * rho_p * binary_op(rho, power - 1.0, rho_power());
    world.gop.fence();

    /////////////////////////////////////////////////////////////////////////////
    /////////               Local and Coulomb potentials                 ////////
    /////////////////////////////////////////////////////////////////////////////

    real_convolution_3d cop = CoulombOperator(world, ntol*0.001, ntol*0.01); //gc
    real_functionT U_c = e2 * cop(rho_p);

    // Make Exchange potential, no screening
    real_functionT U_ex = binary_op(rho_p, rho_p, uex());
    U_ex = e2 * rho_p * U_ex;

    real_functionT rho_e = real_factoryT(world);
    rho_e = 0.0 * rho_e + (1.0 * np/(8.0 * L * L * L));
    double mu = 1.0/screenl;
    real_convolution_3d bop = BSHOperator3D(world, mu, tol, ntol);

    // Make Exchange potential, yes screening
        if (screening == 1) {
        real_functionT lambda  = 0.0 * rho_p + mu;
        real_functionT Fg = binary_op(rho_p, lambda, Fgamma());
        U_ex = e2 * rho_p * Fg * binary_op(rho_p, rho_p, uex());
    }

    if (jellium == 1 && screening == 0) {
        if (world.rank() == 0) {print("   Jellium, no screening");}
        U_c = e2 * cop(rho_p - rho_e);
    }

    if (jellium == 1 && screening == 1) {
        if (world.rank() == 0) {print("   Jellium, screening");}
        U_c = e2 * bop(4.0 * M_PI * (rho_p - rho_e));
    }

    if (jellium == 0 && screening == 1) {
        if (world.rank() == 0) {print("   No jellium, screening");}
        U_c = e2 * bop(4.0 * M_PI * rho_p);
    }
    world.gop.fence();

    if (world.rank() == 0) {print("   Local potential");}
    U_n  = b[0] * rho - bp[0] * rho_n + b[3] * ((alpha + 2.0)/3.0) * rho_1
         - bp[3] * (2.0/3.0) * rho_2_n - bp[3] * (alpha/3.0) * (rho_3_p + rho_3_n)
         + b[1] * tau_p + (b[1] - bp[1]) * tau_n - b[2] * lap + bp[2] * lap_n
         - b[4] * (dJ_p + dJ_n) - bp[4] * dJ_n;

    U_p  = b[0] * rho - bp[0] * rho_p + b[3] * ((alpha + 2.0)/3.0) * rho_1
         - bp[3] * (2.0/3.0) * rho_2_p - bp[3] * (alpha/3.0) * (rho_3_p + rho_3_n)
         + (b[1] - bp[1]) * tau_p + b[1] * tau_n - b[2] * lap + bp[2] * lap_p
         - b[4] * (dJ_p + dJ_n) - bp[4] * dJ_p
         + U_c + U_ex;
    world.gop.fence();

    if (iter > 0 && avg_pot == 1) {
        if (world.rank() == 0) {print("   Mix potentials");}
        double pot_avg = 0.4;
        U_p = pot_avg * U_p + (1.0 - pot_avg) * U_p_old;
        U_n = pot_avg * U_n + (1.0 - pot_avg) * U_n_old;
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////             Energy output and plotting data                ////////
    /////////////////////////////////////////////////////////////////////////////
    
    // Calculate all energy contributions (just for info - not used in iterations)
    //if (iter%10 == 0 || BE == 0.0)
    {
        if (world.rank() == 0) {print("   Energies");}
        double Etotal = Ebinding(world, rho_p, rho_n, tau_p, tau_n, lap_p, lap_n, dJ_p, dJ_n,
                        U_c, psi_pu, psi_pd, iter, spinorbit, b, bp, alpha, k_fn, k_fp);
        BE = Etotal / (1.0 * (nn + np));
    }

    // Calculate single particle state energies (for info)
    if (iter%100 == 0 || iter == 0) {
        Energies(world, U_n, psi_nu, psi_nd, rho, rho_n, energy_n, iter, K_n);
        Energies(world, U_p, psi_pu, psi_pd, rho, rho_p, energy_p, iter, K_p);
    }

    // Make output
    if( iter%10 == 0 || iter == 0 ) {
        output(world, rho_p, rho_n, tau, lap_p, lap_n, U_p, U_c, U_ex, iter, L,
               vtk_output, txt_output);
    }

    world.gop.fence();
    tau_p.clear(); tau_n.clear(); tau.clear(); lap_p.clear(); lap_n.clear(); lap.clear();
    rho_1.clear(); rho_2_p.clear(); rho_3_p.clear(); rho_2_n.clear(); rho_3_n.clear();
    dJ_p.clear(); dJ_n.clear(); U_c.clear(); U_ex.clear();
    world.gop.fence();
}

