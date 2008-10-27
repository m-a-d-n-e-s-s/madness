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

namespace madness {

    /// This class used for callbacks (e.g., for dependency tracking)
    class CallbackInterface {
    public:
        virtual void notify() = 0;

	virtual ~CallbackInterface(){};
    };
        

    /// Provides interface for tracking dependencies
    class DependencyInterface : public CallbackInterface, private Mutex {
    private:
        volatile int ndepend;   ///< Counts dependencies
        
        static const int MAXCALLBACKS = 8;
        typedef Stack<CallbackInterface*,MAXCALLBACKS> callbackT;
        mutable volatile callbackT callbacks; ///< Called ONCE by dec() when ndepend==0

        void do_callbacks() const {
            // ASSUME THAT WE HAVE THE LOCK WHILE IN HERE
            callbackT& cb = const_cast<callbackT&>(callbacks);
            if (cb.size()) {
                while (cb.size()) {
                    CallbackInterface* p = cb.pop();
                    p->notify();
                }
            }
        };

    public:
        DependencyInterface(int ndep = 0) : ndepend(ndep) {}; 


        /// Returns the number of unsatisfied dependencies
        int ndep() const {return ndepend;}

        /// Returns true if ndepend == 0
        bool probe() const {return ndep() == 0;};


        /// Invoked by callbacks to notifiy of dependencies being satisfied
        void notify() {dec();};


        /// Registers a callback for when ndepend=0

        /// If ndepend == 0, the callback is immediately invoked.
        void register_callback(CallbackInterface* callback) {
            ScopedMutex<Mutex> hold(this);
            const_cast<callbackT&>(this->callbacks).push(callback);
            if (ndepend == 0) do_callbacks();
        };
        

        /// Increment the number of dependencies
        void inc() {
            ScopedMutex<Mutex> hold(this);
            MADNESS_ASSERT(ndepend>=0);
            ndepend++;
        };
        

        /// Decrement the number of dependencies and invoke callback if ndepend=0
        void dec() {
            ScopedMutex<Mutex> hold(this);
            if (--ndepend == 0) do_callbacks();
        };


        virtual ~DependencyInterface() {
            ScopedMutex<Mutex> hold(this); 
            if (ndepend) {
                print("DependencyInterface: destructor with ndepend =",ndepend,"?");
            }
        };
    };
}
#endif
