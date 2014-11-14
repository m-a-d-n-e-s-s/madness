/*
  This file is a part of MADNESS.
  Copyright (C) 2014  Virginia Tech

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

#ifndef MADNESS_WORLD_MEM_FUNC_WRAPPER_H__INCLUDED
#define MADNESS_WORLD_MEM_FUNC_WRAPPER_H__INCLUDED

#include <madness/world/worldexc.h>
#include <madness/world/shared_ptr.h>
#include <madness/world/typestuff.h>

namespace madness {
    namespace detail {


        template <typename T>
        struct DefaultInitPtr {
            static T init() { return NULL; }
        }; // struct DefaultInitPtr

        template <typename T>
        struct DefaultInitPtr<std::shared_ptr<T> > {
            static std::shared_ptr<T> init() { return std::shared_ptr<T>(); }
        }; // struct DefaultInitPtr<std::shared_ptr<T> >

        /// Functor wrapper for object and member function pointers

        /// \tparam ptrT Pointer type
        /// \tparam memfnT Member function pointer type
        /// \tparam resT result type of the member function
        template <typename ptrT, typename memfnT, typename resT>
        class MemFuncWrapper {
        private:
            ptrT ptr_;
            memfnT memfn_;

            friend memfnT get_mem_func_ptr<ptrT, memfnT, resT>(const MemFuncWrapper<ptrT, memfnT, resT>&);

        public:
            typedef MemFuncWrapper<ptrT, memfnT, void> MemFuncWrapper_;
            typedef resT result_type;
            typedef memfnT memfn_type;

            MemFuncWrapper() : ptr_(DefaultInitPtr<ptrT>::init()), memfn_() { }

            MemFuncWrapper(const MemFuncWrapper_& other) :
                ptr_(other.ptr_), memfn_(other.memfn_)
            { }

            MemFuncWrapper(ptrT ptr, memfnT memfn) :
                ptr_(ptr), memfn_(memfn)
            { }

            MemFuncWrapper_& operator=(const MemFuncWrapper_& other) {
                ptr_ = other.ptr_;
                memfn_ = other.memfn_;
                return *this;
            }

            resT operator()() const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)();
            }

            template <typename a1T>
            resT operator()(a1T& a1) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1);
            }

            template <typename a1T, typename a2T>
            resT operator()(a1T& a1, a2T& a2) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2);
            }

            template <typename a1T, typename a2T, typename a3T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4, a5);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8, a9T& a9) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8, a9);
            }

            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & ptr_ & memfn_;
            }

            friend memfnT get_mem_func_ptr(const MemFuncWrapper_& wrapper) {
                return wrapper.memfn_;
            }
        }; // class MemFuncWrapper


        /// Functor wrapper for object and member function pointers

        /// \tparam ptrT Pointer type
        /// \tparam memfnT Member function pointer type
        template <typename ptrT, typename memfnT>
        class MemFuncWrapper<ptrT, memfnT, void> {
        private:
            ptrT ptr_;
            memfnT memfn_;

            friend memfnT get_mem_func_ptr<ptrT, memfnT, void>(const MemFuncWrapper<ptrT, memfnT, void>&);

        public:
            typedef MemFuncWrapper<ptrT, memfnT, void> MemFuncWrapper_;
            typedef void result_type;
            typedef memfnT memfn_type;

            MemFuncWrapper() : ptr_(DefaultInitPtr<ptrT>::init()), memfn_() { }

            MemFuncWrapper(const MemFuncWrapper_& other) :
                ptr_(other.ptr_), memfn_(other.memfn_)
            { }

            MemFuncWrapper(ptrT ptr, memfnT memfn) :
                ptr_(ptr), memfn_(memfn)
            { }

            MemFuncWrapper_& operator=(const MemFuncWrapper_& other) {
                ptr_ = other.ptr_;
                memfn_ = other.memfn_;
                return *this;
            }

            void operator()() const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)();
            }

            template <typename a1T>
            void operator()(a1T& a1) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1);
            }

            template <typename a1T, typename a2T>
            void operator()(a1T& a1, a2T& a2) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2);
            }

            template <typename a1T, typename a2T, typename a3T>
            void operator()(a1T& a1, a2T& a2, a3T& a3) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4, a5);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8, a9T& a9) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8, a9);
            }

            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & ptr_ & memfn_;
            }

        }; // class MemFuncWrapper<ptrT, memfnT, void>

        template <typename objT, typename memfnT>
        MemFuncWrapper<objT*, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(objT& obj, memfnT memfn) {
            return MemFuncWrapper<objT*, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(& obj, memfn);
        }

        template <typename objT, typename memfnT>
        MemFuncWrapper<const objT*, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(const objT& obj, memfnT memfn) {
            return MemFuncWrapper<const objT*, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(& obj, memfn);
        }

        template <typename objT, typename memfnT>
        MemFuncWrapper<objT*, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(objT* obj, memfnT memfn) {
            return MemFuncWrapper<objT*, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(obj, memfn);
        }

        template <typename objT, typename memfnT>
        MemFuncWrapper<const objT*, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(const objT* obj, memfnT memfn) {
            return MemFuncWrapper<const objT*, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(obj, memfn);
        }

        template <typename objT, typename memfnT>
        MemFuncWrapper<std::shared_ptr<objT>, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(const std::shared_ptr<objT>& obj, memfnT memfn) {
            return MemFuncWrapper<std::shared_ptr<objT>, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(obj, memfn);
        }

        template <typename objT, typename memfnT>
        MemFuncWrapper<std::shared_ptr<objT>, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(std::shared_ptr<objT>& obj, memfnT memfn) {
            return MemFuncWrapper<std::shared_ptr<objT>, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(obj, memfn);
        }

        template <typename ptrT, typename memfnT, typename resT>
        memfnT get_mem_func_ptr(const MemFuncWrapper<ptrT, memfnT, resT>& wrapper) {
            return wrapper.memfn_;
        }

    }  // namespace detail
} // namespace madness

#endif // MADNESS_WORLD_MEM_FUNC_WRAPPER_H__INCLUDED
