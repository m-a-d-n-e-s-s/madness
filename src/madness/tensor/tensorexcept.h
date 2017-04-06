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

  
#ifndef MADNESS_TENSOR_TENSOREXCPT_H__INCLUDED
#define MADNESS_TENSOR_TENSOREXCPT_H__INCLUDED

/// \file tensorexcept.h
/// \brief Declares and implements TensorException

/// \ingroup tensor

#include <iosfwd>
#include <exception>
#include <madness/tensor/basetensor.h>

namespace madness {
/// Tensor is intended to throw only TensorExceptions
    class TensorException : public std::exception {
        const char* msg;
        const char* assertion;
        int value;
        BaseTensor t;
        const BaseTensor *tp;
        int line;
        const char *function;
        const char *filename;

public:
        // Capturing the line/function/filename info is best done with the macros below
        TensorException(const char* s, const char *a, int err, const BaseTensor* tp,
                        int lin, const char *func, const char *file)
                : msg(s)
                , assertion(a)
                , value(err)
                , tp(tp)
                , line(lin)
                , function(func)
                , filename(file) {
            // We must copy the pointed-to tensor because during unwinding it
            // might go out of scope and dereferencing its pointer is invalid
	    if (tp != nullptr)
              new(&t) BaseTensor(*tp);
        }

        virtual const char* what() const throw() {
            return msg;
        }

        virtual ~TensorException() throw() {}

    /// Print a TensorException to the stream (for human consumption)
    friend std::ostream& operator <<(std::ostream& out, const TensorException& e) {
        out << "TensorException: msg='";
        if (e.msg) out << e.msg;
        out << "'\n";
        if (e.assertion) out << "                 failed assertion='" <<
            e.assertion << "'\n";
        out << "                 value=" << e.value << "\n";
        if (e.line) out << "                 line=" << e.line << "\n";
        if (e.function) out << "                 function='" <<
            e.function << "'\n";
        if (e.filename) out << "                 filename='" <<
            e.filename << "'\n";
        if (e.tp != nullptr) {
            out << "                 tensor=Tensor<";
            if (e.t.id()>=0 && e.t.id()<=TENSOR_MAX_TYPE_ID) {
                out << tensor_type_names[e.t.id()] << ">(";
            }
            else {
                out << "invalid_type_id>(";
            }
            if (e.t.ndim()>=0 && e.t.ndim()<TENSOR_MAXDIM) {
                for (int i=0; i<e.t.ndim(); ++i) {
                    out << e.t.dim(i);
                    if (i != (e.t.ndim()-1)) out << ",";
                }
                out << ")";
            }
            else {
                out << "invalid_dimensions)";
            }
            out << " at 0x" << (void *)(e.tp) << "\n";
        }

        return out;
    }

 };


#define TENSOR_STRINGIZE(X) #X
#define TENSOR_EXCEPTION_AT(F, L) TENSOR_STRINGIZE(F) "(" TENSOR_STRINGIZE(L) ")"

#define TENSOR_EXCEPTION(msg,value,t) \
    throw ::madness::TensorException("TENSOR EXCEPTION: " TENSOR_EXCEPTION_AT( __FILE__, __LINE__ ) ": " msg , \
    0,value,t,__LINE__,__FUNCTION__,__FILE__)

#define TENSOR_ASSERT(condition,msg,value,t) \
do {if (!(condition)) \
        throw ::madness::TensorException("TENSOR ASSERTION FAILED: " TENSOR_EXCEPTION_AT( __FILE__, __LINE__ ) ": " msg , \
        #condition,value,t,__LINE__,__FUNCTION__,__FILE__); \
   } while (0)

}

#endif // MADNESS_TENSOR_TENSOREXCPT_H__INCLUDED
