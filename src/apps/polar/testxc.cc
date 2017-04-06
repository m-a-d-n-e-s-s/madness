/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>
#include <fstream>
#include "xcfunctional.h"

using namespace madness;

static std::string df_repo_functionals[] = {
"lda_x",
"lda_c_vwn_rpa",
"lda_c_vwn",
"lda_c_pz",
"lda_c_pw",
"hyb_gga_xc_b3lyp",
"gga_xc_hcth_93",
"gga_xc_hcth_407",
"gga_xc_hcth_147",
"gga_xc_hcth_120",
"gga_xc_edf1",
"gga_xc_b97_2",
"gga_xc_b97_1",
"gga_xc_b97",
"gga_x_pw91",
"gga_x_pbe",
"gga_x_ft97_b",
"gga_x_b88",
"gga_c_pw91",
"gga_c_pbe",
"gga_c_p86",
"gga_c_lyp"};


struct xcfunc_data_point
{
  double rhoa, rhob;
  double sigmaaa, sigmaab, sigmabb;
  double zk;
  double vrhoa, vrhob;
  double vsigmaaa, vsigmaab, vsigmabb;
  double v2rhoa2, v2rhoab, v2rhob2;
};

std::vector<xcfunc_data_point> read_test_data(const std::string& dfname,
                                              bool spin_polarized)
{
  std::ifstream fstr(dfname.c_str());
  char buffer[120];
  fstr.getline(buffer, 120);
  fstr.getline(buffer, 120);

  std::string tmpstr;

  std::vector<xcfunc_data_point> dps;
  while(!fstr.eof())
  {
    xcfunc_data_point dp;
    fstr >> tmpstr; fstr >> dp.rhoa;
    fstr >> tmpstr; fstr >> dp.rhob;
    fstr >> tmpstr; fstr >> dp.sigmaaa;
    fstr >> tmpstr; fstr >> dp.sigmaab;
    fstr >> tmpstr; fstr >> dp.sigmabb;

    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.zk;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.vrhoa;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.vrhob;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.vsigmaaa;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.vsigmaab;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.vsigmabb;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.v2rhoa2;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.v2rhoab;
    fstr >> tmpstr; fstr >> tmpstr; fstr >> dp.v2rhob2;

    // skip unwanted lines for now
    for (int iskip = 0; iskip < 15; iskip++)
      fstr.getline(buffer,120);

    if (!spin_polarized)
    {
      if (std::abs(dp.vrhoa-dp.vrhob) <= 1e-10)
        dps.push_back(dp);
    }
    else
    {
      dps.push_back(dp);
    }
  }

  return dps;

}

void test_xcfunctional()
{
  bool spin_polarized = false;

  for (int istr = 6; istr < 22; istr++)
  {
    XCfunctional xcfunc;
    std::string xcfuncstr = df_repo_functionals[istr];
    std::cout << "Testing exchange-correlation functional:  "<< xcfuncstr << std::endl;

    xcfuncstr += " 1.0";
    xcfunc.initialize(xcfuncstr,spin_polarized);

    std::string fpath("df_repo/");
    fpath += df_repo_functionals[istr];
    fpath += ".data";
    std::vector<xcfunc_data_point> dps = read_test_data(fpath.c_str(),spin_polarized);

    Tensor<double> rhoa_t((long)dps.size());
    Tensor<double> rhob_t((long)dps.size());
    Tensor<double> sigmaaa_t((long)dps.size());
    Tensor<double> sigmaab_t((long)dps.size());
    Tensor<double> sigmabb_t((long)dps.size());
    std::vector<Tensor<double> > xc_args;
    for (unsigned int idp = 0; idp < dps.size(); idp++)
    {
      rhoa_t(idp) = dps[idp].rhoa;
      rhob_t(idp) = dps[idp].rhob;
      sigmaaa_t(idp) = dps[idp].sigmaaa;
      sigmaab_t(idp) = dps[idp].sigmaab;
      sigmabb_t(idp) = dps[idp].sigmabb;
    }
    if (spin_polarized)
    {
      xc_args.push_back(rhoa_t);
      xc_args.push_back(rhob_t);
      xc_args.push_back(sigmaaa_t);
      xc_args.push_back(sigmaab_t);
      xc_args.push_back(sigmabb_t);
    }
    else
    {
      xc_args.push_back(rhoa_t);
      xc_args.push_back(sigmaaa_t);
    }
    Tensor<double> vr = xcfunc.vxc(xc_args, 0, 0);

    for (unsigned int idp = 0; idp < dps.size(); idp++)
    {
      printf("%25.12e %25.12e  %25.12e %25.12e   %25.12e\n",
          rhoa_t[idp], sigmaaa_t[idp], dps[idp].vrhoa, vr[idp],
          std::abs(dps[idp].vrhoa - vr[idp]));
    }
    print("\n\n");
  }

  spin_polarized = true;

  std::cout << "Testing spin-polarized case: " << std::endl << std::endl;

  for (int istr = 6; istr < 22; istr++)
  {
    XCfunctional xcfunc;
    std::string xcfuncstr = df_repo_functionals[istr];
    std::cout << "Testing exchange-correlation functional:  "<< xcfuncstr << std::endl;

    xcfuncstr += " 1.0";
    xcfunc.initialize(xcfuncstr,spin_polarized);

    std::string fpath("df_repo/");
    fpath += df_repo_functionals[istr];
    fpath += ".data";
    std::vector<xcfunc_data_point> dps = read_test_data(fpath.c_str(),spin_polarized);

    Tensor<double> rhoa_t((long)dps.size());
    Tensor<double> rhob_t((long)dps.size());
    Tensor<double> sigmaaa_t((long)dps.size());
    Tensor<double> sigmaab_t((long)dps.size());
    Tensor<double> sigmabb_t((long)dps.size());
    std::vector<Tensor<double> > xc_args;
    for (unsigned int idp = 0; idp < dps.size(); idp++)
    {
      rhoa_t(idp) = dps[idp].rhoa;
      rhob_t(idp) = dps[idp].rhob;
      sigmaaa_t(idp) = dps[idp].sigmaaa;
      sigmaab_t(idp) = dps[idp].sigmaab;
      sigmabb_t(idp) = dps[idp].sigmabb;
    }
    if (spin_polarized)
    {
      xc_args.push_back(rhoa_t);
      xc_args.push_back(rhob_t);
      if (xcfunc.is_gga())
      {
        xc_args.push_back(sigmaaa_t);
        xc_args.push_back(sigmaab_t);
        xc_args.push_back(sigmabb_t);
      }
    }
    else
    {
      xc_args.push_back(rhoa_t);
      if (xcfunc.is_gga())
      {
        xc_args.push_back(sigmaaa_t);
      }
    }
    Tensor<double> vr = xcfunc.vxc(xc_args,0, 0);

    if (xcfunc.is_spin_polarized())
    {
      for (unsigned int idp = 0; idp < dps.size(); idp++)
      {
        printf("%25.12e %25.12e %25.12e %25.12e %25.12e %25.12e %25.12e %25.12e\n",
            rhoa_t[idp], rhob_t[idp], sigmaaa_t[idp], sigmaab_t[idp], sigmabb_t[idp],
            dps[idp].vrhoa, vr[idp],
            std::abs(dps[idp].vrhoa - vr[idp]));
      }
    }
    else
    {
      for (unsigned int idp = 0; idp < dps.size(); idp++)
      {
        printf("%25.12e %25.12e  %25.12e %25.12e   %25.12e\n",
            rhoa_t[idp], sigmaaa_t[idp], dps[idp].vrhoa, vr[idp],
            std::abs(dps[idp].vrhoa - vr[idp]));
      }
    }
    print("\n\n");
  }

}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();

#ifdef MADNESS_HAS_LIBXC
    test_xcfunctional();
#else
    std::cout << "WARNING: To run this program you need libXC. All tests will be skipped." << std::endl;
#endif

    madness::finalize();
    return 0;
}
