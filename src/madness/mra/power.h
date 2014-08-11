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


#ifndef MADNESS_MRA_POWER_H__INCLUDED
#define MADNESS_MRA_POWER_H__INCLUDED


namespace madness {

    template <int D>
    inline int power(int base = 2) {
        return (int) std::pow((double) base, (int) D);
    }

    template <>
    inline int power<0>(int base) {
        return 1;
    }

    template <>
    inline int power<1>(int base) {
        return base;
    }

    template <>
    inline int power<2>(int base) {
        return (int)(base*base);
    }

    template <>
    inline int power<3>(int base) {
        return (int)(base*base*base);
    }

    template <>
    inline int power<4>(int base) {
        return power<2>(power<2>(base));
    }

    template <>
    inline int power<5>(int base) {
        return (power<2>(base)*power<3>(base));
    }

    template <>
    inline int power<6>(int base) {
        return power<2>(power<3>(base));
    }

    template <>
    inline int power<7>(int base) {
        return (power<4>(base)*power<3>(base));
    }

    template <>
    inline int power<8>(int base) {
        return (power<2>(power<4>(base)));
    }

    template <>
    inline int power<9>(int base) {
        return (power<3>(power<3>(base)));
    }

    template <>
    inline int power<10>(int base) {
        return (power<2>(power<5>(base)));
    }

    template <>
    inline int power<12>(int base) {
        return (power<3>(power<4>(base)));
    }


    // Florian: following code from the internet (stackoverflow)
    namespace detail {
		template<class T, int N>
		struct helper {
			static T pow(const T x){
				return helper<T, N-1>::pow(x) * x;
			}
		};

		template<class T>
		struct helper<T, 1> {
			static T pow(const T x){
				return x;
			}
		};

		template<class T>
		struct helper<T, 0> {
			static T pow(const T x){
				return 1;
			}
		};
    }

    template<int N, class T>
    T power(T const x) {
    	return detail::helper<T, N>::pow(x);
    }

}

#endif
