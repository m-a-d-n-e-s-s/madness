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

  
#ifndef WORLD_DEP_H
#define WORLD_DEP_H

/// \file worlddep.h
/// \brief Defines DependencyInterface and CallbackInterface

namespace madness {

    /// This class used for callbacks (e.g., for dependency tracking)
    class CallbackInterface {
    public:
        virtual void notify() = 0;

	virtual ~CallbackInterface(){};
    };
        

    /// Provides interface for tracking dependencies
    class DependencyInterface : public CallbackInterface {
    private:
        int ndepend;        //< Counts dependencies
        
        //mutable std::vector<CallbackInterface*> callbacks; //< Called ONCE by dec() when ndepend==0
        static const int MAXCALLBACKS = 8;
        mutable Stack<CallbackInterface*,MAXCALLBACKS> callbacks; //< Called ONCE by dec() when ndepend==0

        inline void do_callbacks() const {
//             for (int i=0; i<(int)callbacks.size(); i++)
//                 callbacks[i]->notify();
//             callbacks.clear();
            while (callbacks.size()) callbacks.pop()->notify();
        };

    public:
        DependencyInterface(int ndepend = 0) : 
            ndepend(ndepend), callbacks() {};

        /// Returns the number of unsatisfied dependencies
        inline int ndep() const {
            return ndepend;
        };

        /// Returns true if ndepend == 0
        inline bool probe() const {return ndep() == 0;};


        /// Invoked by callbacks to notifiy of dependencies being satisfied
        void notify() {dec();};


        /// Registers a callback for when ndepend=0

        /// If ndepend == 0, the callback is immediately invoked.
        inline void register_callback(CallbackInterface* callback) {
            //this->callbacks.push_back(callback);
            this->callbacks.push(callback);
            if (ndepend == 0) do_callbacks();
        };
        
        /// Increment the number of dependencies
        inline void inc() {
            MADNESS_ASSERT(ndepend>=0);
            ndepend++;
        };
        
        /// Decrement the number of dependencies and invoke callback if ndepend=0
        inline void dec() {
            ndepend--;
            MADNESS_ASSERT(ndepend>=0);
            if (ndepend == 0) do_callbacks();
        };

        virtual ~DependencyInterface() {
            if (ndepend) {
                print("DependencyInterface: destructor with ndepend =",ndepend,"?");
            }
        };
    };
}
#endif
