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


#ifndef WORLD_DEP_H
#define WORLD_DEP_H

/// \file worlddep.h
/// \brief Defines DependencyInterface and CallbackInterface

#include <world/array.h>
#include <world/print.h>
#include <world/worldthread.h>
#include <world/atomicint.h>

#include <typeinfo>

namespace madness {

    /// This class used for callbacks (e.g., for dependency tracking)
    class CallbackInterface {
    public:
        virtual void notify() = 0;

        virtual ~CallbackInterface() {}
    };


    /// Provides interface for tracking dependencies
    class DependencyInterface : public CallbackInterface, private Spinlock {
    private:
        AtomicInt ndepend;   ///< Counts dependencies
        
        static const int MAXCALLBACKS = 8;
        typedef Stack<CallbackInterface*,MAXCALLBACKS> callbackT;
        mutable volatile callbackT callbacks; ///< Called ONCE by dec() when ndepend==0

        void do_callbacks() const {
            // ASSUME WE HAVE THE LOCK IN HERE
            callbackT& cb = const_cast<callbackT&>(callbacks);

            // Note that we now execute the callback BEFORE popping
            // the stack.  This is so the destructor only has to get
            // the lock if callbacks.size() is non-zero.  We are
            // accessing callbacks thru a non-volatile reference but
            // that is OK since we only rely up on the size in memory
            // being updated after the callbacks are executed.
            while (cb.size()) {
                cb.front()->notify();
                cb.pop();
            }
        };

    public:
        DependencyInterface(int ndep = 0) {
            ndepend = ndep;
        }


        /// Returns the number of unsatisfied dependencies
        int ndep() const {
            return ndepend;
        }

        /// Returns true if ndepend == 0
        bool probe() const {
            return ndep() == 0;
        }


        /// Invoked by callbacks to notifiy of dependencies being satisfied
        void notify() {
            dec();
        }


        /// Registers a callback for when ndepend=0

        /// If ndepend == 0, the callback is immediately invoked.
        void register_callback(CallbackInterface* callback) {
            ScopedMutex<Spinlock> hold(this);
            const_cast<callbackT&>(this->callbacks).push(callback);
            if (ndep() == 0) do_callbacks();
        }


        /// Increment the number of dependencies
        void inc() {
            ndepend++;
        }


        /// Decrement the number of dependencies and invoke callback if ndepend=0
        void dec() {
            if (ndepend.dec_and_test()) {
                ScopedMutex<Spinlock> hold(this);
                do_callbacks();
            }
        }


        virtual ~DependencyInterface() {
            // How to avoid this lock?  It is here because execution of a
            // callback might destroy the object before all callbacks have
            // been performed.  Hence, we only need the lock if the number of
            // unexecuted callbacks is non-zero.

            // But if we whack the object while some poor bugger is holding
            // the lock (probably attempting to free it) bad things will happen.
            // Thus, there is no choice but to get the damn lock.

            ScopedMutex<Spinlock> hold(this);

            // Paranoia is good
            if (ndepend) {
                print("DependencyInterface: destructor with ndepend =",ndepend,"?", typeid(*this).name());
            }
        }
    };
}
#endif
