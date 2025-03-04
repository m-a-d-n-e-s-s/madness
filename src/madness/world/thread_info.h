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

#ifndef MADNESS_WORLD_THREAD_INFO_H__INCLUDED
#define MADNESS_WORLD_THREAD_INFO_H__INCLUDED

/**
 \file thread_info.h
 \brief Implements thread introspection for Pthreads backend
 \ingroup threads
*/

namespace madness {

    /// bitfield describing thread type
    enum ThreadTag {
        ThreadTag_NonMADNESS = 0b0000,
        ThreadTag_MADNESS = 0b0001,
        ThreadTag_Main = 0b0010
    };

    namespace detail {
        inline ThreadTag& thread_tag_accessor() {
            thread_local ThreadTag value = ThreadTag_NonMADNESS;
            return value;
        }
    }   // namespace madness::detail

    /// sets the thread tag for the thread invoking this function
    /// @param tag thread tag for the calling thread
    inline void set_thread_tag(ThreadTag tag) {
        detail::thread_tag_accessor() = tag;
    }

    /// sets the thread tag for the thread invoking this function
    /// @param tag thread tag for the calling thread
    inline void set_thread_tag(int tag) {
        detail::thread_tag_accessor() = static_cast<ThreadTag>(tag);
    }

    /// @return true if this thread is used to execute MADNESS tasks; this is true for the main thread but only true for a non-main thread if the Pthreads backend is used and the thread is part of the ThreadPool
    inline bool is_madness_thread() {
        return detail::thread_tag_accessor() & ThreadTag_MADNESS;
    }

}  // namespace madness

#endif // MADNESS_WORLD_THREAD_INFO_H__INCLUDED
