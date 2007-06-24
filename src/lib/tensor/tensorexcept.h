/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

  $LastChangedDate$
  $Rev$
*/

  
#ifndef TENSOREXCEPT_H
#define TENSOREXCEPT_H

/// \file tensorexcept.h
/// \brief Declares and implements TensorException

namespace madness {
/// Tensor is intended to throws only TensorExceptions
    class TensorException {
    public:
        const char* msg;
        const char* assertion;
        const int value;
        const BaseTensor* t;
        const int line;
        const char *function;
        const char *filename;

        // Capturing the line/function/filename info is best done with the macros below
        TensorException(const char* s, const char *a, int err, const BaseTensor* tp,
                        int lin, const char *func, const char *file)
                : msg(s)
                , assertion(a)
                , value(err)
                , t(tp)
                , line(lin)
                , function(func)
        , filename(file) {};
    };

// implemented in tensor.cc
    std::ostream& operator <<(std::ostream& out, const TensorException& e);

#define TENSOR_EXCEPTION(msg,value,t) \
throw TensorException(msg,0,value,t,__LINE__,__FUNCTION__,__FILE__)

#define TENSOR_ASSERT(condition,msg,value,t) \
do {if (!(condition)) \
    throw TensorException(msg,#condition,value,t,__LINE__,__FUNCTION__,__FILE__); \
   } while (0)

}

#endif
