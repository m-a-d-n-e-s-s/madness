#include "densities.h"
#include "loadbalance.h"
#include "states.h"


// Calculates the laplacian of rho from
//  2nd derivatives (fastes laplacian but also noisiest)
real_functionT laplacian(World& world, real_functionT& rho_q,
                                       const double brad,
                                       const double prec)
{   
    double ntol = FunctionDefaults<3>::get_thresh() * prec;
    
    real_functionT drho_q_dx    = real_derivative_3d(world,0)(rho_q);
    real_functionT drho_q_dy    = real_derivative_3d(world,1)(rho_q);
    real_functionT drho_q_dz    = real_derivative_3d(world,2)(rho_q);
    world.gop.fence();
    
    real_functionT ddrho_q_dxdx = real_derivative_3d(world,0)(drho_q_dx);
    real_functionT ddrho_q_dydy = real_derivative_3d(world,1)(drho_q_dy);
    real_functionT ddrho_q_dzdz = real_derivative_3d(world,2)(drho_q_dz);
    world.gop.fence();    

    real_functionT lap_q = ddrho_q_dxdx + ddrho_q_dydy + ddrho_q_dzdz;
    lap_q.truncate(ntol);
    
    world.gop.fence(); 
    drho_q_dx.clear(); drho_q_dy.clear(); drho_q_dz.clear();
    ddrho_q_dxdx.clear(); ddrho_q_dydy.clear(); ddrho_q_dzdz.clear();
    world.gop.fence();
    
    return lap_q;
}


// Calculates laplacian of rho from single particle states                                       
// (this laplacian is a bit faster than laplacian2 but results in more 
// noise)
real_functionT laplacian1(World& world, const comp_vecfuncT& dpsi_qu_dx,
                                        const comp_vecfuncT& dpsi_qd_dx,
                                        const comp_vecfuncT& dpsi_qu_dy,
                                        const comp_vecfuncT& dpsi_qd_dy,
                                        const comp_vecfuncT& dpsi_qu_dz,
                                        const comp_vecfuncT& dpsi_qd_dz,
                                        const comp_vecfuncT& psi_qu,
                                        const comp_vecfuncT& psi_qd)
{
    comp_functionT laplx = comp_factoryT(world);
    comp_functionT laply = comp_factoryT(world);
    comp_functionT laplz = comp_factoryT(world);

    complex_derivative_3d Dx(world, 0);

    comp_vecfuncT conj0a = conj(world, dpsi_qu_dx, false);
    comp_vecfuncT conj1  = conj(world, dpsi_qd_dx, false);
    comp_vecfuncT conj2  = conj(world, psi_qu, false);
    comp_vecfuncT conj3  = conj(world, psi_qd, false);
    world.gop.fence();

    comp_vecfuncT termx = mul(world, conj0a, psi_qu, false);
    comp_vecfuncT mul1  = mul(world, conj1, psi_qd, false);
    comp_vecfuncT mul2  = mul(world, conj2, dpsi_qu_dx, false);
    comp_vecfuncT mul3  = mul(world, conj3, dpsi_qd_dx, false);
    world.gop.fence();

    conj0a.clear(), conj1.clear(), conj2.clear(), conj3.clear();
    compress(world, mul1,  false);
    compress(world, mul2,  false);
    compress(world, mul3,  false);
    compress(world, termx, false);
    world.gop.fence();

    //gaxpy(world, 1.0, termx, 1.0, mul1, false);
    //world.gop.fence();
    //gaxpy(world, 1.0, termx, 1.0, mul2, false);
    //world.gop.fence();
    //gaxpy(world, 1.0, termx, 1.0, mul3, false);

    gaxpy(world, 1.0, mul2,  1.0, mul3, false);
    gaxpy(world, 1.0, termx, 1.0, mul1, false);
    world.gop.fence();
    gaxpy(world, 1.0, termx, 1.0, mul2, false);

    world.gop.fence();
    mul1.clear(), mul2.clear(), mul3.clear();

    comp_vecfuncT Dtermx = apply(world, Dx, termx, false);
    world.gop.fence();

    termx.clear();
    complex_derivative_3d Dy(world, 1);

    comp_vecfuncT conj0b = conj(world, dpsi_qu_dy, false);
    comp_vecfuncT conj4  = conj(world, dpsi_qd_dy, false);
    comp_vecfuncT conj5  = conj(world, psi_qu, false);
    comp_vecfuncT conj6  = conj(world, psi_qd, false);
    world.gop.fence();

    comp_vecfuncT termy = mul(world, conj0b, psi_qu, false);
    comp_vecfuncT mul4  = mul(world, conj4, psi_qd, false);
    comp_vecfuncT mul5  = mul(world, conj5, dpsi_qu_dy, false);
    comp_vecfuncT mul6  = mul(world, conj6, dpsi_qd_dy, false);
    world.gop.fence();

    conj0b.clear(), conj4.clear(), conj5.clear(), conj6.clear();
    compress(world, mul4,  false);
    compress(world, mul5,  false);
    compress(world, mul6,  false);
    compress(world, termy, false);
    world.gop.fence();

    //gaxpy(world, 1.0, termy, 1.0, mul4, false);
    //world.gop.fence();
    //gaxpy(world, 1.0, termy, 1.0, mul5, false);
    //world.gop.fence();
    //gaxpy(world, 1.0, termy, 1.0, mul6, false);

    gaxpy(world, 1.0, mul5,  1.0, mul6, false);
    gaxpy(world, 1.0, termy, 1.0, mul4, false);
    world.gop.fence();
    gaxpy(world, 1.0, termy, 1.0, mul5, false);

    world.gop.fence();
    mul4.clear(), mul5.clear(), mul6.clear();

    comp_vecfuncT Dtermy = apply(world, Dy, termy, true);
    world.gop.fence();

    termy.clear();
    complex_derivative_3d Dz(world, 2);

    comp_vecfuncT conj0c = conj(world, dpsi_qu_dz, false);
    comp_vecfuncT conj7  = conj(world, dpsi_qd_dz, false);
    comp_vecfuncT conj8  = conj(world, psi_qu, false);
    comp_vecfuncT conj9  = conj(world, psi_qd, false);
    world.gop.fence();

    comp_vecfuncT termz = mul(world, conj0c, psi_qu, false);
    comp_vecfuncT mul7  = mul(world, conj7, psi_qd, false);
    comp_vecfuncT mul8  = mul(world, conj8, dpsi_qu_dz, false);
    comp_vecfuncT mul9  = mul(world, conj9, dpsi_qd_dz, false);
    world.gop.fence();

    conj0c.clear(), conj7.clear(), conj8.clear(), conj9.clear();
    compress(world, termz, false);
    compress(world, mul7, false);
    compress(world, mul8, false);
    compress(world, mul9, false);
    world.gop.fence();

    //gaxpy(world, 1.0, termz, 1.0, mul7, false);
    //world.gop.fence();
    //gaxpy(world, 1.0, termz, 1.0, mul8, false);
    //world.gop.fence();
    //gaxpy(world, 1.0, termz, 1.0, mul9, false);

    gaxpy(world, 1.0, mul8,  1.0, mul9, false);
    gaxpy(world, 1.0, termz, 1.0, mul7, false);
    world.gop.fence();
    gaxpy(world, 1.0, termz, 1.0, mul8, false);

    world.gop.fence();
    mul7.clear(), mul8.clear(), mul9.clear();

    comp_vecfuncT Dtermz = apply(world, Dz, termz, false);
    world.gop.fence();

    termz.clear();

    compress(world, Dtermx, false);
    compress(world, Dtermy, false);
    compress(world, Dtermz, false);
    world.gop.fence();

    laplx = Dtermx[0];
    laply = Dtermy[0];
    laplz = Dtermz[0];

    for (unsigned int i=1; i < Dtermx.size(); i++) {
        laplx.gaxpy(1.0, Dtermx[i], 1.0, true);
        laply.gaxpy(1.0, Dtermy[i], 1.0, true);
        laplz.gaxpy(1.0, Dtermz[i], 1.0, true);
    }

    world.gop.fence();
    return real(laplx) + real(laply) + real(laplz);
}


// Calculates laplacian of rho from single particle states
// (least noisy laplacian, but also most expensive to calculate)
real_functionT laplacian2(World& world, const comp_vecfuncT& dpsi_qu_dx,
                                        const comp_vecfuncT& dpsi_qd_dx,
                                        const comp_vecfuncT& dpsi_qu_dy,
                                        const comp_vecfuncT& dpsi_qd_dy,
                                        const comp_vecfuncT& dpsi_qu_dz,
                                        const comp_vecfuncT& dpsi_qd_dz,
                                        const comp_vecfuncT& psi_qu,
                                        const comp_vecfuncT& psi_qd,
                                        const double prec)
{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;
    comp_functionT lapl = comp_factoryT(world);
    complex_derivative_3d Dx(world, 0);

    comp_vecfuncT c_dpsi_qu_d = conj(world, dpsi_qu_dx);
    comp_vecfuncT c_dpsi_qd_d = conj(world, dpsi_qd_dx);

    comp_vecfuncT   c_ddpsi_qu_dxdx = apply(world, Dx, c_dpsi_qu_d);
    truncate(world, c_ddpsi_qu_dxdx, ntol);
    compress(world, c_ddpsi_qu_dxdx);

    comp_vecfuncT   ddpsi_qu_dxdx = apply(world, Dx, dpsi_qu_dx);
    truncate(world, ddpsi_qu_dxdx, ntol);
    compress(world, ddpsi_qu_dxdx);

    comp_vecfuncT   ddpsi_qd_dxdx = apply(world, Dx, dpsi_qd_dx);
    truncate(world, ddpsi_qd_dxdx, ntol);
    compress(world, ddpsi_qd_dxdx);

    comp_vecfuncT   c_ddpsi_qd_dxdx = apply(world, Dx, c_dpsi_qd_d);
    truncate(world, c_ddpsi_qd_dxdx, ntol);
    compress(world, c_ddpsi_qd_dxdx);


    comp_vecfuncT term = mul(world, c_ddpsi_qu_dxdx, psi_qu);
    truncate(world, term, ntol);
    compress(world, term);
    compress(world, c_dpsi_qu_d);
    compress(world, c_dpsi_qd_d);

    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, psi_qu), ddpsi_qu_dxdx));
    gaxpy(world, 1.0, term, 2.0, mul(world, c_dpsi_qu_d, dpsi_qu_dx));
    gaxpy(world, 1.0, term, 1.0, mul(world, c_ddpsi_qd_dxdx, psi_qd));
    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, psi_qd), ddpsi_qd_dxdx));
    gaxpy(world, 1.0, term, 2.0, mul(world, c_dpsi_qd_d, dpsi_qd_dx));
    truncate(world, term, ntol);

    world.gop.fence();
    c_ddpsi_qu_dxdx.clear();
    c_ddpsi_qd_dxdx.clear();
    ddpsi_qu_dxdx.clear();
    ddpsi_qd_dxdx.clear();
    world.gop.fence();

    truncate(world, term, ntol);
    compress(world, term);

    complex_derivative_3d Dy(world, 1);
    c_dpsi_qu_d = conj(world, dpsi_qu_dy);
    c_dpsi_qd_d = conj(world, dpsi_qd_dy);
    compress(world, c_dpsi_qu_d);
    compress(world, c_dpsi_qd_d);

    gaxpy(world, 1.0, term, 1.0, mul(world, apply(world, Dy, c_dpsi_qu_d), psi_qu));
    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, psi_qu), apply(world, Dy, dpsi_qu_dy)));
    gaxpy(world, 1.0, term, 2.0, mul(world, c_dpsi_qu_d, dpsi_qu_dy));
    gaxpy(world, 1.0, term, 1.0, mul(world, apply(world, Dy, c_dpsi_qd_d), psi_qd));
    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, psi_qd), apply(world, Dy, dpsi_qd_dy)));
    gaxpy(world, 1.0, term, 2.0, mul(world, c_dpsi_qd_d, dpsi_qd_dy));

    truncate(world, term, ntol);
    compress(world, term);

    complex_derivative_3d Dz(world, 2);
    c_dpsi_qu_d = conj(world, dpsi_qu_dz);
    c_dpsi_qd_d = conj(world, dpsi_qd_dz);
    compress(world, c_dpsi_qu_d);
    compress(world, c_dpsi_qd_d);

    gaxpy(world, 1.0, term, 1.0, mul(world, apply(world, Dz, c_dpsi_qu_d), psi_qu));
    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, psi_qu), apply(world, Dz, dpsi_qu_dz)));
    gaxpy(world, 1.0, term, 2.0, mul(world, c_dpsi_qu_d, dpsi_qu_dz));
    gaxpy(world, 1.0, term, 1.0, mul(world, apply(world, Dz, c_dpsi_qd_d), psi_qd));
    gaxpy(world, 1.0, term, 1.0, mul(world, conj(world, psi_qd), apply(world, Dz, dpsi_qd_dz)));
    gaxpy(world, 1.0, term, 2.0, mul(world, c_dpsi_qd_d, dpsi_qd_dz));

    world.gop.fence();
    c_dpsi_qu_d.clear();
    c_dpsi_qd_d.clear();
    truncate(world, term, ntol);
    compress(world, term);
    world.gop.fence();

    lapl = term[0];
    lapl.compress();
    for (unsigned int i=1; i<term.size(); i++) {
        lapl.gaxpy(double_complex(1.0, 0.), term[i], double_complex(1.0,0.), true);
    }

    world.gop.fence();
    return real(lapl);
}





// Calculates number density rho
real_functionT ndensity(World& world, const comp_vecfuncT& psi_qu,
                                      const comp_vecfuncT& psi_qd)
{
    comp_functionT comp_density;

    compress(world, psi_qu, false);
    compress(world, psi_qd, false);
    world.gop.fence();

    comp_vecfuncT conj1 = conj(world, psi_qd, false);
    comp_vecfuncT conj2 = conj(world, psi_qu, false);
    world.gop.fence();

    comp_vecfuncT mul1 = mul(world, conj1, psi_qd, false);
    comp_vecfuncT mul2 = mul(world, conj2, psi_qu, false);
    world.gop.fence();

    // calculate single particle state norms
    comp_vecfuncT sps_norm = add(world, mul2, mul1, false);
    world.gop.fence();

    compress(world, sps_norm, false);
    world.gop.fence();

    comp_density = sps_norm[0];
    comp_density.compress();

    // clean up
    world.gop.fence();
    conj1.clear(), conj2.clear();
    mul1.clear(),  mul2.clear();

    for(unsigned int i = 1; i < psi_qu.size(); i++) {
        comp_density.gaxpy(double_complex(1.0,0.), sps_norm[i], double_complex(1.0,0.), true);
    }

    // clean up
    world.gop.fence();
    sps_norm.clear();

    comp_density.compress();
    world.gop.fence();

    // make sure that the number density is a real function 
    real_functionT density = real(comp_density);

    return binary_op(density, density, abs_func());
}


// Calculates kinetic density tau
real_functionT kdensity(World& world, const comp_vecfuncT& dpsi_qu_dx,
                                        const comp_vecfuncT& dpsi_qd_dx,
                                        const comp_vecfuncT& dpsi_qu_dy,
                                        const comp_vecfuncT& dpsi_qd_dy,
                                        const comp_vecfuncT& dpsi_qu_dz,
                                        const comp_vecfuncT& dpsi_qd_dz)
{
    comp_functionT kdensity = comp_factoryT(world);

    comp_vecfuncT conj1a = conj(world, dpsi_qu_dx, false);
    comp_vecfuncT conj1b = conj(world, dpsi_qd_dx, false);
    world.gop.fence();

    comp_vecfuncT mul1a = mul(world, conj1a, dpsi_qu_dx, false);
    comp_vecfuncT mul1b = mul(world, conj1b, dpsi_qd_dx, false);
    world.gop.fence();

    conj1a.clear(), conj1b.clear();
    world.gop.fence();

    comp_vecfuncT kabs = add(world, mul1a, mul1b, false);
    world.gop.fence();

    mul1a.clear(), mul1b.clear();
    world.gop.fence();

    compress(world, kabs, false);
    compress(world, dpsi_qu_dy, false);
    compress(world, dpsi_qd_dy, false);
    world.gop.fence();

    comp_vecfuncT conj2 = conj(world, dpsi_qu_dy, false);
    comp_vecfuncT conj3 = conj(world, dpsi_qd_dy, false);
    world.gop.fence();

    comp_vecfuncT mul2  = mul(world, conj2, dpsi_qu_dy, false);
    comp_vecfuncT mul3  = mul(world, conj3, dpsi_qd_dy, false);
    world.gop.fence();

    conj2.clear(), conj3.clear();
    world.gop.fence();

    comp_vecfuncT mul4 = add(world, mul2, mul3, true);
    world.gop.fence();

    mul2.clear(), mul3.clear();
    compress(world, mul4, false);
    world.gop.fence();

    gaxpy(world, double_complex(1.0, 0.0), kabs, double_complex(1.0, 0.0), mul4);
    world.gop.fence();

    mul4.clear();
    world.gop.fence();

    compress(world, dpsi_qu_dz, false);
    compress(world, dpsi_qd_dz, false);
    world.gop.fence();

    comp_vecfuncT conj5 = conj(world, dpsi_qu_dz, false);
    comp_vecfuncT conj6 = conj(world, dpsi_qd_dz, false);
    world.gop.fence();

    comp_vecfuncT mul5  = mul(world, conj5, dpsi_qu_dz, false);
    comp_vecfuncT mul6  = mul(world, conj6, dpsi_qd_dz, false);
    world.gop.fence();

    conj5.clear(), conj6.clear();
    world.gop.fence();

    comp_vecfuncT mul7 = add(world, mul5, mul6, false);
    world.gop.fence();

    mul5.clear(), mul6.clear();
    compress(world, mul7, false);
    world.gop.fence();

    gaxpy(world, double_complex(1.0, 0.0), kabs, double_complex(1.0, 0.0), mul7);
    world.gop.fence();

    mul7.clear();
    world.gop.fence();

    compress(world, kabs);
    kdensity.compress();
    world.gop.fence();

    for (unsigned int i = 0; i < kabs.size(); i++) {
        kdensity.gaxpy(double_complex(1.0, 0.0), kabs[i], double_complex(1.0,0.0), true);   // false
    }
    world.gop.fence();

    kdensity.compress();
    real_functionT rkdensity = real(kdensity);
    world.gop.fence();
    return binary_op(rkdensity, rkdensity, abs_func());
}


// Calculates spin-orbit density
real_functionT so_density(World& world, const comp_vecfuncT& dpsi_qu_dx,
                                        const comp_vecfuncT& dpsi_qd_dx,
                                        const comp_vecfuncT& dpsi_qu_dy,
                                        const comp_vecfuncT& dpsi_qd_dy,
                                        const comp_vecfuncT& dpsi_qu_dz,
                                        const comp_vecfuncT& dpsi_qd_dz,
                                        const double prec)
{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;
    int nst = dpsi_qu_dx.size();
    comp_functionT dj_q = comp_factoryT(world);

    const double_complex one(1.0, 0.0);

    comp_vecfuncT conj1 = conj(world, dpsi_qu_dx, false);
    world.gop.fence();
    comp_vecfuncT DJ_q_part = mul(world, conj1, dpsi_qu_dy, false);
    world.gop.fence();
    truncate(world, DJ_q_part, ntol, false);
    world.gop.fence();
    compress(world, DJ_q_part, false);
    world.gop.fence();


    comp_vecfuncT conj2 = conj(world, dpsi_qu_dy, false);
    comp_vecfuncT conj4 = conj(world, dpsi_qu_dz, false);
    world.gop.fence();


    comp_vecfuncT mul1  = mul(world, conj1, dpsi_qd_dz, false);
    comp_vecfuncT mul2  = mul(world, conj2, dpsi_qu_dx, false);
    comp_vecfuncT mul5  = mul(world, conj4, dpsi_qd_dy, false);
    world.gop.fence();

    comp_vecfuncT mul3  = mul(world, conj2, dpsi_qd_dz, false);
    comp_vecfuncT mul4  = mul(world, conj4, dpsi_qd_dx, false);
    world.gop.fence();

    world.gop.fence();
    conj1.clear(), conj2.clear(), conj4.clear();
    world.gop.fence();

    compress(world, mul1, false);
    compress(world, mul2, false);
    compress(world, mul3, false);
    compress(world, mul4, false);
    compress(world, mul5, false);
    world.gop.fence();

    //gaxpy(world, one, DJ_q_part, I, mul1);
    //world.gop.fence();
    //gaxpy(world, one, DJ_q_part, -one, mul2);
    //world.gop.fence();
    //gaxpy(world, one, DJ_q_part, one, mul3);
    //world.gop.fence();
    //gaxpy(world, one, DJ_q_part, -I, mul4);
    //world.gop.fence();
    //gaxpy(world, one, DJ_q_part, -one, mul5);
    
    gaxpy(world, double_complex(1.0,0.0), mul2,  double_complex(1.0,0.0), mul5);
    gaxpy(world, double_complex(1.0,0.0), mul1, -double_complex(1.0,0.0), mul4);
    world.gop.fence();
    gaxpy(world, double_complex(1.0,0.0), mul3,     -double_complex(1.0,0.0), mul2);
    gaxpy(world, double_complex(1.0,0.0), DJ_q_part, double_complex(0.0,1.0), mul1);
    world.gop.fence();
    gaxpy(world, double_complex(1.0,0.0), DJ_q_part, double_complex(1.0,0.0), mul3);

    world.gop.fence();
    mul1.clear(), mul2.clear(), mul3.clear(), mul4.clear(), mul5.clear();
    //loadbalance_v1(world, DJ_q_part);
    world.gop.fence();

    comp_vecfuncT conj6  = conj(world, dpsi_qd_dx, false);
    comp_vecfuncT conj8  = conj(world, dpsi_qd_dy, false);
    comp_vecfuncT conj10 = conj(world, dpsi_qd_dz, false);
    world.gop.fence();


    comp_vecfuncT mul6  = mul(world, conj6, dpsi_qu_dz, false);
    comp_vecfuncT mul11 = mul(world, conj10, dpsi_qu_dy, false);
    world.gop.fence();

    comp_vecfuncT mul7  = mul(world, conj6, dpsi_qd_dy, false);
    comp_vecfuncT mul9  = mul(world, conj8, dpsi_qd_dx, false);
    world.gop.fence();

    comp_vecfuncT mul8  = mul(world, conj8, dpsi_qu_dz, false);
    comp_vecfuncT mul10 = mul(world, conj10, dpsi_qu_dx, false);
    world.gop.fence();

    world.gop.fence();
    conj6.clear(), conj8.clear(), conj10.clear();
    world.gop.fence();

    compress(world, mul6,  false);
    compress(world, mul7,  false);
    compress(world, mul8,  false);
    compress(world, mul9,  false);
    compress(world, mul10, false);
    compress(world, mul11, false);
    world.gop.fence();

    //gaxpy(world, one, DJ_q_part, -I , mul6, false);
    //world.gop.fence();
    //gaxpy(world, one, DJ_q_part, -one, mul7, false);
    //world.gop.fence();
    //gaxpy(world, one, DJ_q_part, one, mul8, false);
    //world.gop.fence();
    //gaxpy(world, one, DJ_q_part, one, mul9, false);
    //world.gop.fence();
    //gaxpy(world, one, DJ_q_part, I, mul10, false);
    //world.gop.fence();
    //gaxpy(world, one, DJ_q_part, -one, mul11, false);

    gaxpy(world, double_complex(1.0,0.0), mul10, -double_complex(1.0,0.0), mul6, false);
    gaxpy(world, double_complex(1.0,0.0), mul8,   double_complex(1.0,0.0), mul9, false);    
    gaxpy(world, double_complex(1.0,0.0), mul7,   double_complex(1.0,0.0), mul11, false);
    world.gop.fence();
    gaxpy(world, double_complex(1.0,0.0), mul8,    -double_complex(1.0,0.0), mul7, false);
    gaxpy(world, double_complex(1.0,0.0), DJ_q_part, double_complex(0.0,1.0), mul10, false);
    world.gop.fence();
    gaxpy(world, double_complex(1.0,0.0), DJ_q_part, double_complex(0.0,1.0), mul8, false);

    world.gop.fence();
    mul6.clear(), mul7.clear(), mul8.clear(), mul9.clear(), mul10.clear(), mul11.clear();
    //loadbalance_v1(world, DJ_q_part);
    world.gop.fence();

    compress(world, DJ_q_part, false);
    world.gop.fence();

    dj_q = DJ_q_part[0];
    world.gop.fence();
    dj_q.compress();
    world.gop.fence();

    for(int i = 1; i < nst; i++) {
        dj_q.gaxpy(double_complex(1.0, 0.0), DJ_q_part[i], double_complex(1.0, 0.0), true); 
    }
    world.gop.fence();

    return real(-1.0 * I * dj_q);
}

