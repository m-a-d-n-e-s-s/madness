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

  
  $Id$
*/

  
#ifndef MAD_BUFAR_H
#define MAD_BUFAR_H

/// \file bufar.h
/// \brief Implements an archive wrapping a memory buffer


#include <world/archive.h>
#include <world/worldexc.h>
#include <world/print.h>


namespace madness {
    namespace archive {
        
        /// Wraps an archive around a memory buffer for output
        
        /// Type checking is disabled for efficiency. 
        ///
        /// Throws MadnessException in case of buffer overflow
        class BufferOutputArchive : public BaseOutputArchive {
        private:
            unsigned char * const ptr;
            const std::size_t nbyte;
            mutable std::size_t i;
        public:
            BufferOutputArchive(void* ptr, std::size_t nbyte) 
                : ptr((unsigned char *) ptr), nbyte(nbyte), i(0) {};

            template <class T>
            inline
            typename madness::enable_if< madness::is_serializable<T>, void >::type
            store(const T* t, long n) const {
                std::size_t m = n*sizeof(T);
                if (i+m > nbyte) {
                    madness::print("BufferOutputArchive:ptr,nbyte,i,n,m,i+m:",(void *)ptr,nbyte,i,n,m,i+m);
                    MADNESS_ASSERT(i+m<=nbyte);
                }
                memcpy(ptr+i, t, m);
                i += m;
            }
            
            void open(std::size_t hint) {};
            
            void close() {};
            
            void flush() {};
            
            inline std::size_t size() const {return i;};
        };
        
        
        /// Wraps an archive around a memory buffer for input
        
        /// Type checking is disabled for efficiency. 
        ///
        /// Throws MadnessException in case of buffer overrun
        class BufferInputArchive : public BaseInputArchive {
        private:
            const unsigned char* const ptr;
            const std::size_t nbyte;
            mutable std::size_t i;
            
        public:
            BufferInputArchive(const void* ptr, std::size_t nbyte) 
                : ptr((const unsigned char *) ptr), nbyte(nbyte), i(0) {};
            
            template <class T>
            inline
            typename madness::enable_if< madness::is_serializable<T>, void >::type
            load(T* t, long n) const {
                std::size_t m = n*sizeof(T);
                MADNESS_ASSERT(m+i <=  nbyte);
                memcpy((unsigned char*) t, ptr+i, m);
                i += m;
            }
            
            void open() {};
            
            void rewind() const {i=0;};
            
            std::size_t nbyte_avail() const {return nbyte-i;};
            
            void close() {}
        };


        // No type checking over Buffer stream for efficiency
        template <class T>
        struct ArchivePrePostImpl<BufferOutputArchive,T> {
            static inline void preamble_store(const BufferOutputArchive& ar) {};
            static inline void postamble_store(const BufferOutputArchive& ar) {};
        };
        
        // No type checking over Buffer stream for efficiency
        template <class T>
        struct ArchivePrePostImpl<BufferInputArchive,T> {
            static inline void preamble_load(const BufferInputArchive& ar) {};
            static inline void postamble_load(const BufferInputArchive& ar) {};
        };
    }
}
#endif
