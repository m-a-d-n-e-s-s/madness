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
  
  $Id$
*/
#include <madness/mra/mra.h>

using std::cout;
using std::endl;
using std::max;

namespace madness {
    using std::abs;

    bool test_rnlp(const bool log_errors) {
        long i, n, l;
        Tensor<double> r;
        double maxerr = 0.0;
        const double pi = 3.14159265358979323846264338328;
        double exact, err;
        const double err_tolerance = 1e-13;
        {
            if (log_errors && err >= err_tolerance) cout << "Testing accuracy of rnlp against results from Maple" << endl;
            double a = 1e-4;
            GaussianConvolution1D<double> g(8,sqrt(a/pi),a,0,false);
            n = 0;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.64057956799766339e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -3.19741057080945913e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.81799864246383464e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.64170777792480805e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.06615886053515213e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.82194856723151343e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.64170777792480805e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.06615886053515213e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.82194856723151343e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.64057956799766339e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 3.19741057080945913e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.81799864246383464e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.63832382494524470e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 5.32546574080415871e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.81010590087371315e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.41047267328181405e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.22018791972360980e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.41047377521394878e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -4.06729836170390096e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.41047377521394878e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.06729836170390096e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.41047267328181405e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.22018791972360980e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.41047046942012737e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.03364123692283645e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52618488461890419e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52618489537996637e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52618489537996637e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52618488461890419e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52618486309677983e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546224281108936e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546224291617701e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546224291617701e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546224281108936e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546224260091135e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556073330313e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556073340579e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556073340579e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556073330313e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556073309815e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183355599e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183355769e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183355769e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183355599e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183355430e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
        }
        {
            double a = 1e-2;
            GaussianConvolution1D<double> g(8,sqrt(a/pi),a,0,false);
            n = 0;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.51198365960967815e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -3.07636514759748877e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.39004237222238520e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 3.05711391835261338e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.62314580091424487e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.06031221284974586e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.77004995102873103e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 1.09765925398799329e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.62314580091424487e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.06031221284974586e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.77004995102873103e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -1.09765925398799329e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.51198365960967815e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 3.07636514759748877e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.39004237222238520e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -3.05711391835261338e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.29620851243244872e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.79304691694024367e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.69661200992196493e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -4.36848047760689648e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.41034540671661907e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.22000442699401760e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.20320972136612724e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.41045559353827048e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -4.06721097934550646e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.20550942666829470e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.41045559353827048e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.06721097934550646e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.20550942666829470e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.41034540671661907e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.22000442699401760e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.20320972136612724e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.41012505889650305e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.03281124249789191e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.19861192757157094e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52618364171646210e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -4.65465182972746039e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52618471782244636e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.55155139906883680e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52618471782244636e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.55155139906883680e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52618364171646210e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.65465182972746039e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52618148950547801e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 7.75774515795399146e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546223067336716e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546224118221681e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546224118221681e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546223067336716e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546220965566786e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556072144996e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556073171247e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556073171247e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556072144996e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556070092493e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183344046e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183354075e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183354075e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183344046e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183323988e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
        }
        {
            double a = 1e0;
            GaussianConvolution1D<double> g(8,sqrt(a/pi),a,0,false);
            n = 0;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 7.63107360346189367e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.39326380508023582e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -3.96841677858563145e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 8.75661865077648466e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 8.82410737907618627e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -1.30949951473492010e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 4.21350396474857447e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -6.17983307932157278e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 4.43199670152667480e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 3.35647464179302069e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -3.57967506112678285e-11;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -1.36118274898628887e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 4.21350396474857447e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 6.17983307932157278e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 4.43199670152667480e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -3.35647464179302069e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -3.57967506112678285e-11;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 1.36118274898628887e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 7.63107360346189367e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.39326380508023582e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -3.96841677858563145e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -8.75661865077648466e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 8.82410737907618627e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 1.30949951473492010e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.32782224202434003e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -3.68552701680425695e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 4.13775061144892418e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 1.60042507022451686e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.96088571198545215e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 1.41367897661257123e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.39768454157893490e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.20177712692920616e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -3.94868219678417840e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 7.31566593570543275e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.40863955444774147e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -4.05848287082555859e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.17520807696470182e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 2.51168860115838423e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.40863955444774147e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.05848287082555859e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.17520807696470182e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -2.51168860115838423e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.39768454157893490e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.20177712692920616e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -3.94868219678417840e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -7.32150795271814975e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.37602944516084069e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.95125817624941723e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -3.51139113424138104e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -1.15110485207894889e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52605935399162088e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -4.65437839301828836e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52616696215172984e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.55153837798979186e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52616696215172984e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.55153837798979186e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52605935399162088e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.65437839301828836e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52584414752294961e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 7.75650823576525769e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546101690129993e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.77560940279273474e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546206778619640e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -5.91869918524889379e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546206778619640e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 5.91869918524889379e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81546101690129993e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.77560940279273474e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81545891513188171e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.95934782871507479e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386555953612550e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556056238053e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386556056238053e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386555953612550e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386555748361586e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390182186538e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183188666e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390183188666e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390182186538e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390180182065e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
        }
        {
            double a = 1e2;
            GaussianConvolution1D<double> g(8,sqrt(a/pi),a,0,false);
            n = 0;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.00000000000000000e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 6.10756752342258014e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.93454443422216120e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -6.74867032786444360e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.24739852592635793e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 1.43294933024935476e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.00000000000000000e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -6.10756752342258014e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.93454443422216120e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 6.74867032786444360e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.24739852592635793e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -1.43294933024935476e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.99318492136080461e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.08696368094917623e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 2.60447168323340892e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -9.04034504985510627e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -2.14998633928712021e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 5.58835531816329497e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.24648176437683600e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -3.27908149234218621e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.67276950422990212e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 1.64624909472826329e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 6.70688497777513430e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -9.05569878193988709e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.24648176437683600e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 3.27908149234218621e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.67276950422990212e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -1.64624909472826329e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 6.70688497777513430e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 9.05569878193988709e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.99318492136080461e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.08696368094917623e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 2.60447168323340892e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 9.04034504985510627e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -2.14998633928712021e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -5.58835531816329497e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.38179858827323476e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.04763936483394230e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -8.32772220224538394e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -1.89437614141419528e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 2.88769508210870525e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 2.21803339122174113e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.51365573776570017e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -4.62710590578847031e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -6.11623077961541036e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52439220749447890e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.55023685987324124e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -6.24931708968097842e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.52439220749447890e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.55023685987324124e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -6.24931708968097842e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.51365573776570017e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.62710590578847031e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -6.11623077961541036e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.49228081923241285e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 7.63366360227468385e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -5.85369846273186721e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81533964065561448e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.77556865689395017e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81544472821516045e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -5.91867978227500295e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81544472821516045e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 5.91867978227500295e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81533964065561448e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.77556865689395017e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81512946929470659e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.95916350513851956e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386544100369154e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -6.77341329040006257e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386554362917535e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.25780460536204789e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386554362917535e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.25780460536204789e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386544100369154e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 6.77341329040006257e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386523575273816e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.12890203983799996e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390066432148e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390166652414e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390166652414e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966390066432148e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966389865991790e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
        }
        {
            double a = 1e4;
            GaussianConvolution1D<double> g(8,sqrt(a/pi),a,0,false);
            n = 0;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.00000000000000000e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.23528272214066237e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 1.41178217654244809e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 1.26724337383508057e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 9.24808666820860048e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 5.16058662037806704e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.00000000000000000e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.23528272214066237e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 1.41178217654244809e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -1.26724337383508057e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 9.24808666820860048e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -5.16058662037806704e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.93406055705484376e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.39497184635340733e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 4.08935561473207502e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 2.66300430413872028e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.33508873895376489e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 5.08640731141814058e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.00000000000000000e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.34688631911635626e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -3.48488099568184218e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 2.62343309407099060e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.42392358397174803e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -5.74390741499019122e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.00000000000000000e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.34688631911635626e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -3.48488099568184218e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -2.62343309407099060e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.42392358397174803e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 5.74390741499019122e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.93406055705484376e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -4.39497184635340733e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 4.08935561473207502e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -2.66300430413872028e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.33508873895567407e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -5.08640731150656195e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.49146408254436480e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.50647102869495280e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 2.64278820712894470e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 3.15088994971621284e-11;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -3.70274884997262088e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 2.02622223000723244e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.35475425674241334e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.42579697710634181e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.61160183056987646e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 4.87550411609523399e-11;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 2.47120513789010578e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -1.86594916521298455e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.35475425674241334e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.42579697710634181e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.61160183056987646e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -4.87550411609523399e-11;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 2.47120513789010578e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 1.86594916521298455e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.49146408254436480e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.50647102869495280e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 2.64278820712894470e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -3.15088994971621284e-11;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -3.70274884997262088e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -2.02622223000727306e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.37408108427214515e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.09038475936154631e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 3.68884880267772498e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 3.66150600926587787e-11;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 3.77722151485774017e-16;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -1.03592207658699442e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.80321162122826961e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.77149821469458344e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -9.25065232468101948e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81371108109231360e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -5.91673982825609260e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -9.32840028258958169e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.81371108109231360e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 5.91673982825609260e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -9.32840028258958169e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.80321162122826961e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.77149821469458344e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -9.25065232468101948e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.78225020949705004e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.94078061354701635e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -9.09598921287351880e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20385358779692841e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -6.77335257389880792e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386385030982762e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.25780171409048271e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20386385030982762e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.25780171409048271e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20385358779692841e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 6.77335257389880792e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20383306291449504e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.12887457303027408e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966378490998915e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.58385225835037511e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966388513019073e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -8.61284112227903377e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966388513019073e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 8.61284112227903377e-20;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966378490998915e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.58385225835037511e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966358446959223e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.30642016947285981e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
        }
        {
            double a = 1e6;
            GaussianConvolution1D<double> g(8,sqrt(a/pi),a,0,false);
            n = 0;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.00000000000000000e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.31393925175791226e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 1.76043401220533680e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 2.07091850787566534e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 2.28733404941548857e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 2.42623389685600399e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=0;
            r = g.rnlp(n,l);
            i = 9;
            exact = -2.07091850787566534e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 2.28733404941548857e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -2.42623389685600399e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=0;
            r = g.rnlp(n,l);
            i = 9;
            exact = -3.50842272140536116e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.68215019539905630e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -1.77563550367174916e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.64684064677916715e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.77372027143066234e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 2.48546774012851273e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 6.34469829004409155e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 8.20364207845620860e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 5.02838600880059851e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 7.99999973531593511e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.28136988133445895e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.97337646440325010e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 2.42595127151247827e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -1.01342125713415232e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 1.19548791588270655e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 7.99999973531593511e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.28136988133445895e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.97337646440325010e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -2.42595127151247827e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -1.01342125713415232e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -1.19548791588270655e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.64684064677916715e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -4.77372027143066234e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 2.48546774012851273e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -6.34469829004409155e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 8.20364207845620860e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -5.02838600880059851e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 7.68112666770073726e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.40335952063720378e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.16852188536233089e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 2.27739586876326983e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -9.64560459462787522e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.64340268954690671e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -5.72613518237629813e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -8.32996156019435581e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 1.22452367554029254e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 3.28001683556767755e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 8.64340268954690671e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 5.72613518237629813e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -8.32996156019435581e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -1.22452367554029254e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 3.28001683556767755e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 7.68112666770073726e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.40335952063720378e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.16852188536233089e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -2.27739586876326983e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -9.64560459462787522e-18;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 6.06603017048240378e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.52348873733158482e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 4.50962719862293780e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -9.54505770868048384e-13;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -3.34357728315444496e-17;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20266863365880416e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -6.76728333900617373e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.38645113252605220e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20369453020092898e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.25751260692224813e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.39098307532258949e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20369453020092898e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.25751260692224813e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.39098307532258949e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20266863365880416e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 6.76728333900617373e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.38645113252605220e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.20061827312391323e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.12613077452750488e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.37740622815052894e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50965220949072187e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.58384321086449390e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966223149730427e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -8.61283681394682842e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50966223149730427e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 8.61283681394682842e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50965220949072187e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.58384321086449390e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50963216553224777e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.30637924047538213e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
        }
        {
            double a = 1e8;
            GaussianConvolution1D<double> g(8,sqrt(a/pi),a,0,false);
            n = 0;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.00000000000000000e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.32198023075050020e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 1.79850757146013485e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 2.16840442140282663e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 2.47807154022197462e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 2.74638507267385013e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 0;
            l=0;
            r = g.rnlp(n,l);
            i = 9;
            exact = -2.16840442140282663e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 2.47807154022197462e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -2.74638507267385013e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=0;
            r = g.rnlp(n,l);
            i = 12;
            exact = 8.66639007453691868e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -8.91729870909936828e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=0;
            r = g.rnlp(n,l);
            i = 9;
            exact = -6.97748517721321448e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.97852113142792607e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 2.44271302079553987e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.77621741053728029e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.80561296694321666e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 2.14965049420876670e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 2.87273668988151739e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -5.78715292978526717e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 3.75446225456850381e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.19822378257343587e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.59910355189205866e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 1.88940266256636957e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -6.48224605439561008e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 8.95940657127647480e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -2.01590163575141197e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.19822378257343587e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.59910355189205866e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 1.88940266256636957e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 6.48224605439561008e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 8.95940657127647480e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 2.01590163575141197e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.77621741053728029e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.80561296694321666e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 2.14965049420876670e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -2.87273668988151739e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -5.78715292978526717e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -3.75446225456850381e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.60267034903721741e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.56247390565564865e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 1.00204874077533328e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -1.68527185127306691e-11;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.32633156186227036e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -4.97857828210297353e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.08775684768582366e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -6.18391032576240027e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -9.17718829164970323e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 7.27641771748153090e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.45035821991845720e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.18688012530858096e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.22880079087442982e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.33180134100535050e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 2.89311654858050329e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 3.32835978333160121e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.18688012530858096e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.22880079087442982e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.33180134100535050e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -2.89311654858050329e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 3.32835978333160121e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.08775684768582366e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 6.18391032576240027e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -9.17718829164970323e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -7.27641771748153090e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 1.45035821991845720e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.90278525861057943e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 8.78775137729619093e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.52594114788140973e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -8.24285673816764725e-14;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -1.10123909261290486e-19;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50849480740676434e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.58293860288335944e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.07063600559659463e-16;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50949687272003708e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -8.61240599236203343e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.07327540364256062e-16;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50949687272003708e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 8.61240599236203343e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.07327540364256062e-16;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50849480740676434e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.58293860288335944e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.07063600559659463e-16;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.50649122351290821e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 4.30228801910451485e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.06536152937712491e-16;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
        }
        {
            double a = 1e10;
            GaussianConvolution1D<double> g(8,sqrt(a/pi),a,0,false);
            n = 0;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.00000000000000000e-01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.32278609519677226e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 1.80234849054933588e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 2.17834302709225502e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = 2.49780041123299235e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 2.78011463893526400e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 4;
            l=0;
            r = g.rnlp(n,l);
            i = 12;
            exact = 9.85994473913733493e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -1.08963047031019240e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 8;
            l=0;
            r = g.rnlp(n,l);
            i = 12;
            exact = 3.17378652019804868e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -3.10371707224047384e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 12;
            l=0;
            r = g.rnlp(n,l);
            i = 9;
            exact = -3.55193478516153949e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -8.76624415465518325e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 4.45211896310375366e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.95753699701407591e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 1.01817522582965636e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.26716664222736262e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 9.05830429656977076e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -4.09146442457520497e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 2.46678858892003107e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.24040427814888830e+02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -6.51252735045721654e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 7.18690903438790468e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -4.13584998317234614e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -1.74068139575386427e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 5.03038032706968060e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 1.24040427814888830e+02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 6.51252735045721654e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 7.18690903438790468e-02;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 4.13584998317234614e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -1.74068139575386427e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -5.03038032706968060e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 3.95753699701407591e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.01817522582965636e+00;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -4.26716664222736262e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -9.05830429656977076e-06;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -4.09146442457520497e-08;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = -2.46678858892003107e-10;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 16;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 2.03517586076878889e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -1.49274934514448811e-03;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = 9.66085723513747846e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -9.25718189118173804e-07;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 12;
            exact = -1.51883153358735018e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 15;
            exact = 2.59891477668870922e-12;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.39414070815156137e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -2.49387003254647314e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.78376303140302747e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 1.87233252628364406e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=-1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.49300601109782676e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = -8.56943996085143253e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.03858147916271143e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = 6.67700051743989307e-16;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=0;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.49300601109782676e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 8.56943996085143253e-05;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -2.03858147916271143e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -6.67700051743989307e-16;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=1;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.39414070815156137e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 2.49387003254647314e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.78376303140302747e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -1.87233252628364406e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            n = 20;
            l=2;
            r = g.rnlp(n,l);
            i = 0;
            exact = 5.20171632178656012e+01;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 3;
            exact = 3.90952885279176158e-04;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 6;
            exact = -1.31483311173371317e-09;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
            i = 9;
            exact = -2.71474799023069306e-15;
            err = abs(r(i)-exact);
            if (log_errors && err >= err_tolerance) cout << a << " " << n << " " << l << " " << i << " " << exact << " " << r(i) << " " << err << endl;
            maxerr = max(maxerr,err);
        }
        if (log_errors && maxerr >= err_tolerance) cout << "MAXERR " << maxerr << endl;
        return (maxerr < err_tolerance);
    }

    bool test_rnlp_rangelimited(const bool log_errors) {
      double maxerr = 0.0;
      const double pi = 3.14159265358979323846264338328;
      double err;
      const auto err_tolerance = 1e-13;
      {
        const int alog10s[] = {-4, -2, 0, 2, 4, 6, 8, 10};
        const int ps[] = {0, 3};
        const int ns[] = {0, 1, 2};
        const int ls[] = {-2, -1, 0, 1, 2};
        /// these do not include 2^{n/2} factor baked into the GaussianConvolution1D (via the multiscale adjustment of the normalization constant?)
        const double exact_tr_values
            [sizeof(alog10s) / sizeof(int)][sizeof(ns) / sizeof(int)]
            [sizeof(ls) / sizeof(int)][sizeof(ps) / sizeof(int)] = {
                {{{0, 0},
                  {0.002820924410015775, -0.000932925282557223},
                  {0.002820924410015775, 0.000932925282557223},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.002820924410015775, -3.331885264499784e-14},
                  {0.002820924410015775, 3.331885264499784e-14},
                  {0, 0},
                  {0, 0}},
                 {{0.001410453389628289, -3.123609897563486e-15},
                  {0.001410471020387486, -1.041224991186846e-15},
                  {0.001410471020387486, 1.041224991186846e-15},
                  {0.001410453389628289, 3.123609897563486e-15},
                  {0, 0}}},
                {{{0, 0},
                  {0.02818598889850831, -0.00931386887557233},
                  {0.02818598889850831, 0.00931386887557233},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.02818598889850831, -3.327307292171618e-9},
                  {0.02818598889850831, 3.327307292171618e-9},
                  {0, 0},
                  {0, 0}},
                 {{0.0140841872463443, -3.116101632173124e-10},
                  {0.01410180165216401, -1.040867135941334e-10},
                  {0.01410180165216401, 1.040867135941334e-10},
                  {0.0140841872463443, 3.116101632173124e-10},
                  {0, 0}}},
                {{{0, 0},
                  {0.2602499389065233, -0.07865853072930514},
                  {0.2602499389065233, 0.07865853072930514},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.2602499389065233, -0.0002901773850574695},
                  {0.2602499389065233, 0.0002901773850574695},
                  {0, 0},
                  {0, 0}},
                 {{0.1220867438224048, -0.00002440297110419346},
                  {0.1381631950841185, -0.00001005737075613944},
                  {0.1381631950841185, 0.00001005737075613944},
                  {0.1220867438224048, 0.00002440297110419346},
                  {0, 0}}},
                {{{0, 0},
                  {0.4999999999992313, 0.6107567523421995},
                  {0.4999999999992313, -0.6107567523421995},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.4999999999992313, 0.205938225992832},
                  {0.4999999999992313, -0.205938225992832},
                  {0, 0},
                  {0, 0}},
                 {{0.0002034760079537496, 0.0002159002695817145},
                  {0.4997965239912775, -0.03888049185695126},
                  {0.4997965239912775, 0.03888049185695126},
                  {0.0002034760079537496, -0.0002159002695817145},
                  {0, 0}}},
                {{{0, 0},
                  {0.5, 1.235282722140662},
                  {0.5, -1.235282722140662},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.5, 1.151568853395805},
                  {0.5, -1.151568853395805},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.5, 0.995420060567129},
                  {0.5, -0.995420060567129},
                  {0, 0},
                  {0, 0}}},
                {{{0, 0},
                  {0.5, 1.313939251757912},
                  {0.5, -1.313939251757912},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.5, 1.305042444690875},
                  {0.5, -1.305042444690875},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.5, 1.287367262429561},
                  {0.5, -1.287367262429561},
                  {0, 0},
                  {0, 0}}},
                {{{0, 0},
                  {0.5, 1.3219802307505},
                  {0.5, -1.3219802307505},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.5, 1.32108520274184},
                  {0.5, -1.32108520274184},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.5, 1.319296336685672},
                  {0.5, -1.319296336685672},
                  {0, 0},
                  {0, 0}}},
                {{{0, 0},
                  {0.5, 1.322786095196772},
                  {0.5, -1.322786095196772},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.5, 1.322696538829787},
                  {0.5, -1.322696538829787},
                  {0, 0},
                  {0, 0}},
                 {{0, 0},
                  {0.5, 1.322517438001069},
                  {0.5, -1.322517438001069},
                  {0, 0},
                  {0, 0}}}};
        auto ia = 0;
        for (auto alog10 : alog10s) {
          const auto a = pow(10., alog10);
          GaussianConvolution1D<double> tg(8, sqrt(a / pi), a, 0, false, 0., 1);
          GaussianConvolution1D<double> g(8, sqrt(a / pi), a, 0, false, 0.);

          auto in = 0;
          for (auto n : ns) {

            auto il = 0;
            for (auto l : ls) {

              const auto tr = tg.rnlp(n, l);

              auto ip = 0;
              for (auto p : ps) {
                const auto exact_tr = std::pow(M_SQRT2,n) * exact_tr_values[ia][in][il][ip];
                err = abs(tr(p) - exact_tr);
                if (log_errors && err >= err_tolerance) {
                  cout << a << " " << n << " " << l << " " << p << " "
                       << exact_tr << " " << tr(p) << " " << err << endl;
                }
                maxerr = max(maxerr, err);
                ip++;
              }
              il++;
            }
            in++;
          }
          ia++;
        }
      }

      return (maxerr < err_tolerance);
    }

    bool test_rnlp_rangelimited_erf(const bool log_errors) {
      double maxerr = 0.0;
      const double pi = 3.14159265358979323846264338328;
      double err;
      const auto err_tolerance = 3.3e-13; /* I suspect the reference values are not any more precise than this */
      {
        const int alog10s[] = {-4, 0, 4};
        const int ps[] = {0, 3, 7};
        const int ns[] = {0, 1, 2};
        const int ls[] = {-2, -1, 0, 1, 2};
        const int σlog10s[] = {-2, -1, 0};
        /// these do not include 2^{n/2} factor baked into the GaussianConvolution1D (via the multiscale adjustment of the normalization constant?)
        const double exact_tr_values
            [sizeof(alog10s) / sizeof(int)][sizeof(σlog10s) / sizeof(int)]
            [sizeof(ns) / sizeof(int)][sizeof(ls) / sizeof(int)]
            [sizeof(ps) / sizeof(int)] = {{{{{0,0,0},{0.002820924395911391,-0.0009318063413062203,
                                                         -0.0004243830414805509},{0.002820924395911391,0.0009318063413062203,
                                              0.00042438304148055957},{0,0,0},{0,0,0}},
                                            {{0.00001591508932156427,0.00003779555818054598,0.000036617002143814076},
                                             {0.002805009306589827,0.0000377955890572059,0.00003661702159918728},
                                             {0.002805009306589827,-0.00003779558905720578,-0.00003661702159918547},
                                             {0.00001591508932156427,-0.00003779555818054597,-0.00003661700214380724},
                                             {0,0,0}},{{0.0013945382862023397,0.000033807197824626645,
                                              0.000020477638782848537},{0.0014104710203874877,-1.0411890544416563e-15,
                                              -1.6957271727974724e-17},{0.0014104710203874877,1.041264904692606e-15,
                                              -1.6967441166862152e-17},{0.0013945382862023397,-0.00003380719782462656,
                                              -0.000020477638782846738},{0.000015915089321564267,-0.00003380717268513069,
                                              -0.000020477638054641078}}},{{{4.178170993418692e-17,
                                              9.855895855395598e-17,9.459945406185167e-17},
                                             {0.0028209229995781386,-0.000826572761246539,-0.00024345581006846996},
                                             {0.0028209229995781386,0.0008265727612465391,0.0002434558100684785},
                                             {4.178170993418692e-17,-9.855895855395597e-17,-9.45994540618344e-17},
                                             {0,0,0}},{{0.00015915020600155503,0.00011931312084461618,
                                              -6.289260984209179e-6},{0.0026617727935765835,0.00011931293482868884,
                                              -6.289332344783605e-6},{0.0026617727935765835,-0.0001193129348286887,
                                              6.289332344784805e-6},{0.00015915020600155503,-0.00011931312084461618,
                                              6.289260984216599e-6},{4.1781731504516154e-17,-8.767838933949522e-17,
                                              -5.4754345758947715e-17}},{{0.0012513220167850525,
                                              0.000020183026420882527,1.5444010562782382e-7},
                                             {0.0014104507767915308,2.2571114123512386e-8,1.332154115046817e-9},
                                             {0.0014104507767915308,-2.2571114123453666e-8,-1.332154110762501e-9},
                                             {0.0012513220167850525,-0.000020183026420882463,-1.5444010562660777e-7},
                                             {0.0001591299634856514,-0.000020183640686752766,-1.5443132799552644e-7}}},
                                           {{{0.0005387571200970742,6.406585944812713e-6,3.00861499978693e-9},
                                             {0.0028208791553762202,-0.000018476538773579174,-5.150075146259928e-9},
                                             {0.0028208791553762202,0.000018476538773579303,5.1500751548823895e-9},
                                             {0.0005387571200970742,-6.4065859448126955e-6,-3.0086149956949922e-9},
                                             {0.000024110182048810063,-2.7521271452166193e-6,7.186304650539819e-10}},
                                            {{0.0010283165458855514,-1.0119950509696495e-6,-1.3119534020791191e-11},
                                             {0.0017925626094906692,-1.0122164982684448e-6,-1.3126711167049264e-11},
                                             {0.0017925626094906692,1.012216498268539e-6,1.3126674287158598e-11},
                                             {0.0010283165458855514,1.0119950509696984e-6,1.3119511421621816e-11},
                                             {0.00042135005068331396,-7.7356481427139e-8,-1.4861211009784566e-11}},
                                            {{0.0008036740477723163,-7.435135138641168e-8,-7.727319706994729e-14},
                                             {0.0009888885617183529,-4.877000030532881e-8,-1.76940419171849e-14},
                                             {0.0009888885617183529,4.8770000305373456e-8,1.7673093494330113e-14},
                                             {0.0008036740477723163,7.435135138643921e-8,7.725576168166772e-14},
                                             {0.000606764986163338,7.434375817569097e-8,7.724179426320363e-14}}}},
                                          {{{{0,0,0},{0.26023895481078885,-0.07857138858370427,
                                                         -0.03322010711931271},{0.26023895481078885,0.07857138858370429,
                                              0.03322010711931377},{0,0,0},{0,0,0}},
                                            {{0.0012339872371034185,0.0029315174860420064,0.0028442154352782863},
                                             {0.2590049675736854,0.0026654120913752544,0.0028593564280339307},
                                             {0.2590049675736854,-0.0026654120913752444,-0.002859356428033791},
                                             {0.0012339872371034183,-0.0029315174860420055,-0.002844215435277754},
                                             {0,0,0}},{{0.1208417724895668,0.0026183278995025785,
                                              0.0015951365688988814},{0.13816319508411862,-0.000010057370756136057,
                                              -2.553462005667635e-11},{0.13816319508411862,0.000010057370756143887,
                                              2.5531308515578387e-11},{0.1208417724895668,-0.002618327899502572,
                                              -0.0015951365688987138},{0.001233987237103418,-0.002623150195491068,
                                              -0.0015945687961760704}}},{{{1.5083891790033054e-15,
                                              3.565692611787655e-15,3.449480227249289e-15},
                                             {0.2591582795849656,-0.07040459977553397,-0.019165813275565966},
                                             {0.2591582795849656,0.07040459977553398,0.019165813275567017},
                                             {1.5083891790033054e-15,-3.565692611787655e-15,-3.449480227248656e-15},
                                             {0,0,0}},{{0.011828532546311768,0.009373979519692642,
                                              -0.0004654553738560623},{0.24732974703865382,0.008944386407754032,
                                              -0.0005229143887150141},{0.24732974703865382,-0.00894438640775402,
                                              0.0005229143887151022},{0.011828532546311768,-0.009373979519692642,
                                              0.00046545537385663395},{1.5083905197015853e-15,-3.1790544183974636e-15,
                                              -2.017435862297216e-15}},{{0.1091684687049073,0.001316879321155076,
                                              0.000015609438479641545},{0.13816127833374653,-7.934123402749891e-6,
                                              1.2136172114235441e-7},{0.13816127833374653,7.934123402755763e-6,
                                              -1.2136172073900235e-7},{0.1091684687049073,-0.0013168793211550702,
                                              -0.000015609438479517983},{0.01182740800819827,-0.001818558870972423,
                                              -8.588759174373388e-6}}},{{{0.01095315211141335,
                                              0.0017592303634754408,-1.773472055278654e-6},
                                             {0.22808502180327883,-0.006607516600575764,-2.3092403293847193e-6},
                                             {0.22808502180327883,0.006607516600575772,2.3092403302455585e-6},
                                             {0.01095315211141335,-0.0017592303634754406,1.7734720554923817e-6},
                                             {0.00002199254764344912,-0.000011981723363715384,-7.815515267869049e-8}},
                                            {{0.06130217827291403,-0.00009508077383050712,2.495651154633274e-8},
                                             {0.1667828435303648,-0.00059386372114806,-4.0952904983310634e-8},
                                             {0.1667828435303648,0.0005938637211480676,4.095290262082e-8},
                                             {0.06130217827291403,0.00009508077383050979,-2.4956511278195987e-8},
                                             {0.010220196785229523,-0.00017809647929935645,2.966496945184784e-9}},
                                            {{0.06982777471182787,-0.000036495266669874715,-1.2148723556440645e-10},
                                             {0.09695506881853691,-0.000032925481535526365,-1.4600710290998759e-10},
                                             {0.09695506881853691,0.000032925481535530153,1.4600499331341376e-10},
                                             {0.06982777471182787,0.00003649526666987838,1.2148561889888637e-10},
                                             {0.04135189809576318,0.000017218932338011736,-5.061041038431077e-11}}}},
                                          {{{{0,0,0},{0.5000000000000006,1.2352827221406644,
                                                         1.3934393915350205},{0.5000000000000006,-1.235282722140664,
                                              -1.3934393915347183},{0,0,0},{0,0,0}},
                                            {{0,0,0},{0.5000000000000006,1.1515688533958064,0.971495788445294},
                                             {0.5000000000000006,-1.1515688533958062,-0.971495788445133},{0,0,0},
                                             {0,0,0}},{{0,0,0},
                                             {0.5000000000000006,0.9954200605671301,0.406093458003402},
                                             {0.5000000000000006,-0.9954200605671301,-0.40609345800337426},{0,0,0},
                                             {0,0,0}}},{{{0,0,0},
                                             {0.49999999999924527,1.2352827221388498,1.393439391533205},
                                             {0.49999999999924527,-1.2352827221388496,-1.3934393915329026},{0,0,0},
                                             {0,0,0}},{{0,0,0},
                                             {0.49999999999924527,1.1515688533941655,0.9714957884442748},
                                             {0.49999999999924527,-1.1515688533941653,-0.9714957884441134},{0,0,0},
                                             {0,0,0}},{{0,0,0},
                                             {0.49999999999924527,0.995420060565806,0.40609345800330227},
                                             {0.49999999999924527,-0.9954200605658062,-0.4060934580032746},{0,0,0},
                                             {0,0,0}}},{{{0,0,0},
                                             {0.3788799981204239,0.9361706288484614,1.0565968861847341},
                                             {0.3788799981204239,-0.9361706288484609,-1.0565968861845048},{0,0,0},
                                             {0,0,0}},{{0,0,0},
                                             {0.3788799981204239,0.8728487850758117,0.7372872173052638},
                                             {0.3788799981204239,-0.8728487850758114,-0.7372872173051412},{0,0,0},
                                             {0,0,0}},{{0,0,0},
                                             {0.3788799981204239,0.7547241841524247,0.30911695286978946},
                                             {0.3788799981204239,-0.7547241841524247,-0.30911695286976826},{0,0,0},
                                             {0,0,0}}}}};
        auto ia = 0;
        for (auto alog10 : alog10s) {
          const auto a = pow(10., alog10);
          auto iσ = 0;
          for (auto σlog10 : σlog10s) {
            const auto σ = pow(10., σlog10);
            GaussianConvolution1D<double> tg(8, sqrt(a / pi), a, 0, false, 0.,
                                             KernelRange(1, σ));
            GaussianConvolution1D<double> g(8, sqrt(a / pi), a, 0, false, 0.);

            auto in = 0;
            for (auto n : ns) {

              auto il = 0;
              for (auto l : ls) {

                const auto tr = tg.rnlp(n, l);

                auto ip = 0;
                for (auto p : ps) {
                  const auto exact_tr = std::pow(M_SQRT2, n) *
                                        exact_tr_values[ia][iσ][in][il][ip];
                  err = abs(tr(p) - exact_tr);
                  if (log_errors && err >= err_tolerance) {
                    cout << a << " " << σ << " " << n << " " << l << " " << p
                         << " " << exact_tr << " " << tr(p) << " " << err
                         << endl;
                  }
                  maxerr = max(maxerr, err);
                  ip++;
                }
                il++;
              }
              in++;
            }
            iσ++;
          }
          ia++;
        }
      }

      return (maxerr < err_tolerance);
    }

    bool test_rnlij_rangelimited(const bool log_errors) {
      double maxerr = 0.0;
      const double pi = 3.14159265358979323846264338328;
      double err;
      const auto err_tolerance = 1e-13;
      {
        const int alog10s[] = {0, 2};
        const int ns[] = {0, 1, 2};
        const int ls[] = {-4, -3, -2, -1, 0, 1, 2, 3, 4};
        const int kmax = 3;
        const int D = 1;

        const double exact_tr_values
            [sizeof(alog10s) / sizeof(int)][sizeof(ns) / sizeof(int)]
            [sizeof(ls) / sizeof(int)][kmax + 1][kmax + 1] = {
                {{{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0.06239914704001694, -0.07295744403619789,
                    0.0369357677412754, -0.001904406234260277},
                   {0.07295744403619789, -0.08048503031516791,
                    0.03190461853341899, 0.00910210942647432},
                   {0.0369357677412754, -0.03190461853341899,
                    -0.003852872191759357, 0.02302140423433085},
                   {0.001904406234260277, 0.00910210942647432,
                    -0.02302140423433085, 0.01862442139435793}},
                  {{0.3957015837330126, 0, -0.07387153548255079, 0},
                   {0, 0.1760059389991514, 0, -0.03139833227690751},
                   {-0.07387153548255079, 0, 0.03137919753706641, 0},
                   {0, -0.03139833227690751, 0, -0.02788296971843104}},
                  {{0.06239914704001694, 0.07295744403619789,
                    0.0369357677412754, 0.001904406234260277},
                   {-0.07295744403619789, -0.08048503031516791,
                    -0.03190461853341899, 0.00910210942647432},
                   {0.0369357677412754, 0.03190461853341899,
                    -0.003852872191759357, -0.02302140423433085},
                   {-0.001904406234260277, 0.00910210942647432,
                    0.02302140423433085, 0.01862442139435793}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}},
                 {{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0.1247982940800339, -0.07567279010025063,
                    0.002406686110882992, 0.0003586857987265181},
                   {0.07567279010025063, -0.007347125872983151,
                    -0.02977841036946986, 0.001329880617009164},
                   {0.002406686110882992, 0.02977841036946986,
                    -0.002836223356072624, -0.01871527459360532},
                   {-0.0003586857987265181, 0.001329880617009164,
                    0.01871527459360532, -0.001259489520761049}},
                  {{0.2709032896529788, 0, -0.004813372221765983, 0},
                   {0, 0.01091766074249585, 0, -0.000174713124627314},
                   {-0.004813372221765983, 0, 0.0002690098419760488, 0},
                   {0, -0.000174713124627314, 0, 4.764945643053651e-6}},
                  {{0.1247982940800339, 0.07567279010025063,
                    0.002406686110882992, -0.0003586857987265181},
                   {-0.07567279010025063, -0.007347125872983151,
                    0.02977841036946986, 0.001329880617009164},
                   {0.002406686110882992, -0.02977841036946986,
                    -0.002836223356072624, 0.01871527459360532},
                   {0.0003586857987265181, 0.001329880617009164,
                    -0.01871527459360532, -0.001259489520761049}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}},
                 {{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0.05914480661826465, -0.0352961498075315,
                    0.0008511266384317965, 0.00003456132305199318},
                   {0.0352961498075315, -0.002370553564141238,
                    -0.01432084200059843, 0.000309382263802944},
                   {0.0008511266384317965, 0.01432084200059843,
                    -0.0006699147777973966, -0.00930029157727551},
                   {-0.00003456132305199318, 0.000309382263802944,
                    0.00930029157727551, -0.000307372220005017}},
                  {{0.1313069749235385, -0.004668508830198368,
                    -0.0005298213407802086, 0.00001813531878048097},
                   {0.004668508830198368, 0.001188374184419877,
                    -0.00006241841553979117, -4.178039760209366e-6},
                   {-0.0005298213407802086, 0.00006241841553979117,
                    6.392819570634216e-6, -4.072217041620477e-7},
                   {-0.00001813531878048097, -4.178039760209366e-6,
                    4.072217041620477e-7, 2.417479793162679e-8}},
                  {{0.1395963147294403, 0, -0.0006426105953031759, 0},
                   {0, 0.001442060388883122, 0, -5.867533905635802e-6},
                   {-0.0006426105953031759, 0, 8.980617332941763e-6, 0},
                   {0, -5.867533905635802e-6, 0, 4.00121333360345e-8}},
                  {{0.1313069749235385, 0.004668508830198368,
                    -0.0005298213407802086, -0.00001813531878048097},
                   {-0.004668508830198368, 0.001188374184419877,
                    0.00006241841553979117, -4.178039760209366e-6},
                   {-0.0005298213407802086, -0.00006241841553979117,
                    6.392819570634216e-6, 4.072217041620477e-7},
                   {0.00001813531878048097, -4.178039760209366e-6,
                    -4.072217041620477e-7, 2.417479793162679e-8}},
                  {{0.05914480661826465, 0.0352961498075315,
                    0.0008511266384317965, -0.00003456132305199318},
                   {-0.0352961498075315, -0.002370553564141238,
                    0.01432084200059843, 0.000309382263802944},
                   {0.0008511266384317965, -0.01432084200059843,
                    -0.0006699147777973966, 0.00930029157727551},
                   {0.00003456132305199318, 0.000309382263802944,
                    -0.00930029157727551, -0.000307372220005017}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}}},
                {{{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0.02820947917699604, -0.04453012417103717,
                    0.04756936948027392, -0.04191644513315848},
                   {0.04453012417103717, -0.07019262711553054,
                    0.0747499664442194, -0.06550833223839031},
                   {0.04756936948027392, -0.0747499664442194,
                    0.07906098685098572, -0.06845248583481641},
                   {0.04191644513315848, -0.06550833223839031,
                    0.06845248583481641, -0.05798132187123867}},
                  {{0.94358104164447, 0, -0.0951387389605478, 0},
                   {0, 0.8318715041031744, 0, -0.1593889860964341},
                   {-0.0951387389605478, 0, 0.723479401312235, 0},
                   {0, -0.1593889860964341, 0, 0.6203974508809246}},
                  {{0.02820947917699604, 0.04453012417103717,
                    0.04756936948027392, 0.04191644513315848},
                   {-0.04453012417103717, -0.07019262711553054,
                    -0.0747499664442194, -0.06550833223839031},
                   {0.04756936948027392, 0.0747499664442194,
                    0.07906098685098572, 0.06845248583481641},
                   {-0.04191644513315848, -0.06550833223839031,
                    -0.06845248583481641, -0.05798132187123867}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}},
                 {{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0.05641895835399209, -0.08039999430492184,
                    0.06916711686405792, -0.04229641365344437},
                   {0.08039999430492184, -0.113770391733522, 0.0962281120412795,
                    -0.05670779950056223},
                   {0.06916711686405792, -0.0962281120412795,
                    0.07801595248299861, -0.04157611664901071},
                   {0.04229641365344437, -0.05670779950056223,
                    0.04157611664901071, -0.01629194504123798}},
                  {{0.8871620832904784, 0, -0.1383342337281158, 0},
                   {0, 0.6705132832080146, 0, -0.179622550701738},
                   {-0.1383342337281158, 0, 0.4787790951354259, 0},
                   {0, -0.179622550701738, 0, 0.3222142289420809}},
                  {{0.05641895835399209, 0.08039999430492184,
                    0.06916711686405792, 0.04229641365344437},
                   {-0.08039999430492184, -0.113770391733522,
                    -0.0962281120412795, -0.05670779950056223},
                   {0.06916711686405792, 0.0962281120412795,
                    0.07801595248299861, 0.04157611664901071},
                   {-0.04229641365344437, -0.05670779950056223,
                    -0.04157611664901071, -0.01629194504123798}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}},
                 {{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0.00001435241351443475, -0.00002152375883897738,
                    0.00002083312096219267, -0.00001600260181591134},
                   {0.00002152375883897738, -0.00003224468529552152,
                    0.00003114338849154711, -0.00002384226488442202},
                   {0.00002083312096219267, -0.00003114338849154711,
                    0.00002994801558399589, -0.00002276961064552485},
                   {0.00001600260181591134, -0.00002384226488442202,
                    0.00002276961064552485, -0.00001712431182868865}},
                  {{0.1128092118809553, -0.1261656431623922,
                    0.06465783383904638, -0.0142939514887932},
                   {0.1261656431623922, -0.1347065150820696,
                    0.05982467094752078, -0.00541888437694819},
                   {0.06465783383904638, -0.05982467094752078,
                    0.01342028920843553, 0.01068448196651819},
                   {0.0142939514887932, -0.00541888437694819,
                    -0.01068448196651819, 0.01161773908828892}},
                  {{0.7743528714095231, 0, -0.1293573339200171, 0},
                   {0, 0.3950780610746883, 0, -0.08836768614958088},
                   {-0.1293573339200171, 0, 0.1636764991511842, 0},
                   {0, -0.08836768614958088, 0, 0.05606202228640211}},
                  {{0.1128092118809553, 0.1261656431623922, 0.06465783383904638,
                    0.0142939514887932},
                   {-0.1261656431623922, -0.1347065150820696,
                    -0.05982467094752078, -0.00541888437694819},
                   {0.06465783383904638, 0.05982467094752078,
                    0.01342028920843553, -0.01068448196651819},
                   {-0.0142939514887932, -0.00541888437694819,
                    0.01068448196651819, 0.01161773908828892}},
                  {{0.00001435241351443475, 0.00002152375883897738,
                    0.00002083312096219267, 0.00001600260181591134},
                   {-0.00002152375883897738, -0.00003224468529552152,
                    -0.00003114338849154711, -0.00002384226488442202},
                   {0.00002083312096219267, 0.00003114338849154711,
                    0.00002994801558399589, 0.00002276961064552485},
                   {-0.00001600260181591134, -0.00002384226488442202,
                    -0.00002276961064552485, -0.00001712431182868865}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},
                  {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}}}};

        auto ia = 0;
        for (auto alog10 : alog10s) {
          const auto a = pow(10., alog10);
          GaussianConvolution1D<double> tg(8, sqrt(a / pi), a, 0, false, 0., D);

          auto in = 0;
          for (auto n : ns) {

            auto il = 0;
            for (auto l : ls) {

              const auto tr = tg.rnlij(n, l);

              for (int ii = 0; ii != kmax; ++ii) {
                for (int ij = 0; ij != kmax; ++ij) {
                  const auto exact_tr = exact_tr_values[ia][in][il][ii][ij];
                  err = abs(tr(ii, ij) - exact_tr);
                  if (log_errors && err >= err_tolerance) {
                    cout << a << " " << n << " " << l << " " << ii << " " << ij
                         << " " << exact_tr << " " << tr(ii, ij) << " " << err
                         << endl;
                  }
                  maxerr = max(maxerr, err);
                }
              }
              il++;
            }
            in++;
          }
          ia++;
        }
      }

      return (maxerr < err_tolerance);
    }

}  // namespace madness
