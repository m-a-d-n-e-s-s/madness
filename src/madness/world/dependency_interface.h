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

/**
 \file dependency_interface.h
 \brief Defines \c DependencyInterface and \c CallbackInterface.
 \ingroup world
*/

#ifndef MADNESS_WORLD_DEPENDENCY_INTERFACE_H__INCLUDED
#define MADNESS_WORLD_DEPENDENCY_INTERFACE_H__INCLUDED

#include <atomic>
#include <cassert>
#include <typeinfo>

#include <madness/world/stack.h>
#include <madness/world/worldmutex.h>
#include <madness/world/atomicint.h>
#include <madness/world/world.h>
#include <madness/world/print.h>

#include <cassert>
#include <iterator>
#include <numeric>
#include <mutex>
#include <typeinfo>
#include <unordered_map>

namespace madness {

    /// The class used for callbacks (e.g., dependency tracking).
    class CallbackInterface {
    public:
        /// Invoked by the callback to notify when a dependency is satisfied.
        virtual void notify() = 0;

        /// Same as notify(), but tracks how many times called from each \c caller
        virtual void notify_debug(const char* caller) {
          notify_debug_impl(caller);
          this->notify();
        }

//        virtual ~CallbackInterface() = default;
        virtual ~CallbackInterface() {}

    protected:
      virtual void notify_debug_impl(const char* caller) {
#if !defined(NDEBUG)
        {
          mx_.lock();
          const auto caller_str = caller;
          auto it = callers_.find(caller_str);
          if (it != callers_.end())
            it->second += 1;
          else
            callers_[caller] = 1;
          mx_.unlock();
        }
#endif
      }

    private:
#if !defined(NDEBUG)
        std::unordered_map<std::string,int> callers_;
        std::mutex mx_;
#endif
    };


    /// Provides an interface for tracking dependencies.
    class DependencyInterface : public CallbackInterface, private Spinlock {
    private:
        // Dependency counter is still atomic to allow fast(er) probing from outside critical sections
        // note that all *updates* to this counter will occur from critical sections
        // thus this is for lock-free probing only
        std::atomic<int> ndepend; ///< Counts dependencies.

        static const int MAXCALLBACKS = 8; ///< Maximum number of callbacks.
        using callbackT = Stack<CallbackInterface*,MAXCALLBACKS>; ///< \todo Brief description needed.
        mutable volatile callbackT callbacks; ///< Called ONCE by \c dec() when `ndepend==0`.
#if !defined(NDEBUG)
        mutable volatile bool finalized = false;  ///< Set to true when register_final_callback is called; not safe to execute any callbacks after the final callback
        volatile int max_ndepend; ///< max value of \c ndepend
        constexpr static const bool print_debug = false;  // change to true to log state changes in {inc,dec,notify}_debug, (debug) ctor, and dtor
#endif

        /// \todo Brief description needed.

        /// Main design point is that, because a callback might destroy this
        /// object, when callbacks are invoked we cannot be holding the lock
        /// and all necessary data must be on the stack (i.e., not from the
        /// object state).
        /// \todo Parameter description needed.
        /// \param[in,out] cb Description needed.
        void do_callbacks(callbackT& cb) const {
            while (!cb.empty()) {
                cb.top()->notify();
                cb.pop();
            }
        }

    public:
        /// \todo Constructor that...

        /// \param[in] ndep The number of unsatisfied dependencies.
        DependencyInterface(int ndep = 0) : ndepend(ndep) {}

        /// Same as ctor above, but keeps track of \c caller . If a given object was constructed
        /// via this ctor, or DependencyInterface::inc_debug() had been called once, debugging variants
        /// of mutating calls ( \c {inc/dec/notify}_debug ) must be used for the rest of this object's lifetime.
        /// \param[in] ndep The number of unsatisfied dependencies.
        DependencyInterface(int ndep, const char* caller) : ndepend(ndep)
#if !defined(NDEBUG)
            , max_ndepend(ndep)
#endif
        {
#if !defined(NDEBUG)
          callers_[caller] = ndep;
          if (print_debug)
            print("DependencyInterface debug ctor: this=", this, " caller=", caller, " ndepend=", ndepend);
#endif
        }

        /// Returns the number of unsatisfied dependencies.

        /// \return The number of unsatisfied dependencies.
        int ndep() const { return ndepend; }

#if !defined(NDEBUG)
        /// Returns the number of unsatisfied dependencies tracked by the debugging instrumentation.

        /// \return The number of unsatisfied dependencies tracked via _debug calls (or debug ctor).
        int ndep_debug() const {
          return std::accumulate(begin(callers_), end(callers_), 0,
                                 [](int partial_sum, auto& x) { return partial_sum + x.second; });
        }
#endif

        /// Returns true if `ndepend == 0` (no unsatisfied dependencies).

        /// \return True if there are no unsatisfied dependencies.
        bool probe() const {
          return ndep() == 0;
        }

        /// Invoked by callbacks to notify of dependencies being satisfied.
        void notify() { dec(); }

        /// Overload of CallbackInterface::notify_debug(), updates dec()
        void notify_debug(const char* caller) {
#if !defined(NDEBUG)
          CallbackInterface::notify_debug_impl(caller);
#endif
          this->dec_debug(caller);
        }

        /// \brief Registers a callback that will be executed when `ndepend==0`; immediately invoked
        ///    if `ndepend==0`.

        /// \param[in] callback The callback to use.
        void register_callback(CallbackInterface* callback) {
            callbackT cb;
            {
                ScopedMutex<Spinlock> obolus(this);
#if !defined(NDEBUG)
                if (print_debug && !callers_.empty())
                  print("DependencyInterface::register_callback: this=", this, " ndepend=", ndepend);
                if (finalized)
                  error("DependencyInterface::register_callback() cannot be called after register_final_callback");
#endif
                const_cast<callbackT&>(callbacks).push(callback);
                if (probe()) {
                    cb = std::move(const_cast<callbackT&>(callbacks));
                }
            }
            do_callbacks(cb);
        }


        /// \brief Registers the final callback to be executed when `ndepend==0`; immediately invoked
        ///    if `ndepend==0`.

        /// No additional callbacks can be registered after this call since execution
        ///  of the final callback can cause destruction of this object.
        /// \param[in] callback The callback to use.
        void register_final_callback(CallbackInterface* callback) {
          callbackT cb;
          {
            ScopedMutex<Spinlock> obolus(this);
#if !defined(NDEBUG)
            MADNESS_ASSERT(finalized == false);
            if (print_debug && !callers_.empty())
              print("DependencyInterface::register_final_callback: this=", this, " ndepend=", ndepend);
            if (finalized)
              error("DependencyInterface::register_final_callback() called more than once");
#endif
            const_cast<callbackT&>(callbacks).push(callback);
            if (probe()) {
              cb = std::move(const_cast<callbackT&>(callbacks));
#if !defined(NDEBUG)
              finalized = true;
#endif
            }
          }
          do_callbacks(cb);
        }

        /// Increment the number of dependencies.
        void inc() {
            ScopedMutex<Spinlock> obolus(this);
#if !defined(NDEBUG)
            if (finalized && callbacks.empty())
              error("DependencyInterface::inc() called after the final callback had been executed");
            if (!callers_.empty())
              error("DependencyInterface::inc() called for an object that is being debugged", "");
#endif
            ++ndepend;
        }

        /// Decrement the number of dependencies and invoke the callback if `ndepend==0`.
        void dec() {
            callbackT cb;
            {
                ScopedMutex<Spinlock> obolus(this);
                MADNESS_ASSERT(ndepend > 0);
#if !defined(NDEBUG)
                if (finalized && callbacks.empty())
                  error("DependencyInterface::dec() called after the final callback had been executed");
                if (!callers_.empty())
                  error("DependencyInterface::dec() called for an object that is being debugged", "");
#endif
                if (ndepend == 1) {
                    cb = std::move(const_cast<callbackT&>(callbacks));
                }
                // NB safe to update ndepend now, was not safe to do that before since that makes it observable and
                //    e.g. result in its destruction before all changes to state are done
                --ndepend;
            }
            do_callbacks(cb);
        }

        /// Same as inc(), but keeps track of \c caller; calling dec_debug() will signal error if no matching inc_debug() had been invoked        
        void inc_debug(const char* caller) {
          ScopedMutex<Spinlock> obolus(this);
#if !defined(NDEBUG)
          if (finalized && callbacks.empty())
              error("DependencyInterface::inc_debug() called after the final callback had been executed: caller =", caller);
#endif
          ++ndepend;
#if !defined(NDEBUG)
          const auto caller_str = caller;
          auto it = callers_.find(caller_str);
          if (it != callers_.end())
            it->second += 1;
          else
            callers_[caller] = 1;
          max_ndepend = std::max(static_cast<int>(max_ndepend), static_cast<int>(ndepend));
          if (ndep() != ndep_debug())
            error("DependencyInterface::inc_debug(): ndepend != ndepend_debug, caller = ", caller);
          if (print_debug)
            print("DependencyInterface::inc_debug: this=", this, " caller=", caller, " ndep=", callers_[caller], " ndepend=", ndepend);
#endif
        }

        void dec_debug(const char* caller) {
            callbackT cb;
            {
                ScopedMutex<Spinlock> obolus(this);
                MADNESS_ASSERT(ndepend > 0);
#if !defined(NDEBUG)
                const auto caller_str = caller;
                auto it = callers_.find(caller_str);
                if (it != callers_.end()) {
                  MADNESS_ASSERT(it->second > 0);
                }
                else {
                  assert(false && "DependencyInterface::dec_debug() called without matching inc_debug()");
                }
#endif
                if (ndepend == 1) {
                    cb = std::move(const_cast<callbackT&>(callbacks));
#if !defined(NDEBUG)
                    if (ndep() != ndep_debug())
                      error("DependencyInterface::dec_debug(): ndepend != ndepend_debug, caller = ", caller);
                    if (print_debug)
                      print("DependencyInterface::dec_debug: callback spawned, this=", this, " caller=", caller, " ndep=", it->second-1, " ndepend=", ndepend-1);
#endif
                }
#if !defined(NDEBUG)
                else {
                  if (print_debug)
                    print("DependencyInterface::dec_debug: this=", this, " caller=", caller, " ndep=", it->second-1, " ndepend=", ndepend-1);
                }
#endif
                // NB safe to update ndepend now, was not safe to do that before since that makes it observable and
                //    e.g. result in its destruction before all changes to state are done
#if !defined(NDEBUG)
                it->second -= 1;
#endif
                --ndepend;
            }
            do_callbacks(cb);
        }

        /// Destructor.
        virtual ~DependencyInterface() {
#ifdef MADNESS_ASSERTIONS_THROW
            if(ndepend != 0)
                error("DependencyInterface::~DependencyInterface(): ndepend =", ndepend);
#else
            MADNESS_ASSERT(ndepend == 0);
#endif
#if !defined(NDEBUG)
            if (!callers_.empty()) {
              MADNESS_ASSERT(max_ndepend > 0);
              if (print_debug)
                print("DependencyInterface dtor: this=", this, " max_ndepend=", max_ndepend);
            }
#endif
#if !defined(NDEBUG)
          for(const auto& c: callers_) {
            if (c.second != 0) {
              const auto error_msg = std::string("ndepend=") + std::to_string(c.second) + " for caller " + c.first;
              error("DependencyInterface::~DependencyInterface(): ", error_msg);
            }
          }
#endif
        }

private:
#if !defined(NDEBUG)
        std::unordered_map<std::string, int> callers_;
#endif
    };
}
#endif // MADNESS_WORLD_DEPENDENCY_INTERFACE_H__INCLUDED
