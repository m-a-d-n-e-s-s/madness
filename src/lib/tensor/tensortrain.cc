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


  $Id: test.cc 2816 2012-03-23 14:59:52Z 3ru6ruWu $
*/

/// \file tensor/test.cc
/// \brief New test code for Tensor class using Google unit test

#include <tensor/tensor.h>
#include <tensor/gentensor.h>
#include <world/print.h>
#include <tensor/tensortrain.h>

using namespace madness;

template<typename T>
int tt_svd(const Tensor<T>& t, const double eps) {
  abort();
  return -1;
}

int main(int argc, char** argv) {

	const long d=5;
	const long k=9;
	const double eps=1.e-5;
	print("entering tensortrain, eps=",eps);

	const std::vector<long> dim(d,k);
	Tensor<double> t(dim);
	print("constructed tensor with d=",d," and k=",k,"; size=",t.size());
	t.fillrandom();


	double cpu0=cpu_time();
	TensorTrain<double> tt(t,eps);
	double cpu1=cpu_time();
	print("time 1",cpu1-cpu0);

	cpu0=cpu1;
	Tensor<double> t2=tt.reconstruct();
	cpu1=cpu_time();
	print("time reconstruct",cpu1-cpu0);

	cpu0=cpu1;
	print("norms", t.normf()," " ,t2.normf(), (t-t2).normf());
	cpu1=cpu_time();
	print("time normf",cpu1-cpu0);


	Tensor<double> U, VT, s;
	tt.two_mode_representation(U,VT,s);
	Tensor<double> t3=inner(U,VT);
	print("norms 2-mode", t.normf()," " ,t3.normf(), (t-t3).normf());



	cpu0=cpu1;
	GenTensor<double> g(t,eps,TT_2D);
	cpu1=cpu_time();
	print("time 2",cpu1-cpu0);

	print("norms gentensor", t.normf()," " ,g.normf(), (t-g.full_tensor_copy()).normf());

	print("leaving tensortrain");
	return 0;
}

