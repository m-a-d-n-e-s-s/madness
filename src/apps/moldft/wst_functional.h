//#define WORLD_INSTANTIATE_STATIC_TEMPLATES  

#include <madness/mra/mra.h>
#include <chem/xcfunctional.h>

using namespace madness;

class WSTFunctional
{
  public:

  real_function_3d make_dft_potential(World& world, const XCfunctional& xc, const vector_real_function_3d& vf, int ispin, int what)
  {
    return multiop_values<double, xc_potential, 3>(xc_potential(xc, ispin, what), vf);
  }
  
  double make_dft_energy(World& world, const XCfunctional& xc, const vector_real_function_3d& vf, int ispin)
  {
    real_function_3d vlda = multiop_values<double, xc_functional, 3>(xc_functional(xc, ispin), vf);
    return vlda.trace();
  }

  std::pair<real_function_3d, double> apply_xc(World& world, const XCfunctional& xc, real_function_3d& rho)
  {
    vector_real_function_3d delrho;
    vector_real_function_3d vf;

    rho.reconstruct();
    vf.push_back(rho);
    if (xc.is_gga()) 
    {
      for(int axis = 0; axis < 3; ++axis)
      {
        Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
        delrho.push_back(D(rho));
      }
      real_function_3d saa = delrho[0]*delrho[0]+delrho[1]*delrho[1]+delrho[2]*delrho[2];
      vf.push_back(saa); // sigma_aa
      if (vf.size()) 
      {
        reconstruct(world, vf);
//        rho.refine_to_common_level(vf); // Ugly but temporary (I hope!)
        refine_to_common_level(world,vf); // Ugly but temporary (I hope!)
      }
    }
    double exc = make_dft_energy(world, xc, vf, 0);
    real_function_3d vxc =  make_dft_potential(world, xc, vf, 0, 0);

    if (xc.is_gga()) 
    {
      real_function_3d vsigaa = make_dft_potential(world, xc, vf, 0, 1).truncate();
      for (int axis=0; axis<1; axis++)
      {
        Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
        real_function_3d  gradn = D(rho);
        real_function_3d  ddel = vsigaa*gradn;
        ddel.scale(2.0);

        real_function_3d vxc2=D(ddel).truncate();
        vxc = vxc - vxc2;
      }
    } 
    return std::pair<real_function_3d, double>(vxc, exc);
  }
};

