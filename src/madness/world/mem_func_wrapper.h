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

/**
 \file mem_func_wrapper.h
 \brief Defines tools for encapsulating a pointer to a member function.
 \ingroup libraries

 The member pointer is stored, along with a pointer to the object with which
 it should be dereferenced.

 \todo If we're moving to C++11, could we replace this with builtin features
      (\c bind)?
*/

#ifndef MADNESS_WORLD_MEM_FUNC_WRAPPER_H__INCLUDED
#define MADNESS_WORLD_MEM_FUNC_WRAPPER_H__INCLUDED

#include <madness/world/madexception.h>

namespace madness {
    namespace detail {

        /// Default pointer to a. object of type \c T.

        /// Returns a \c NULL pointer, by default.
        /// \tparam T The pointer type.
        template <typename T>
        struct DefaultInitPtr {

            /// Get a default pointer.

            /// \return \c NULL.
            static T init() {
                return NULL;
            }
        }; // struct DefaultInitPtr


        /// Default shared pointer to an object of type \c T.

        /// Returns a \"NULL\" shared pointer, by default. Specialization for
        /// shared pointers.
        /// \tparam T The type pointed to by the shared pointer.
        template <typename T>
        struct DefaultInitPtr<std::shared_ptr<T> > {

            /// Get a default shared pointer.

            /// \return Default shared pointer.
            static std::shared_ptr<T> init() {
                return std::shared_ptr<T>();
            }
        }; // struct DefaultInitPtr<std::shared_ptr<T> >

        /// Functor wrapper for object and member function pointers.

        /// This class encapsulates a pointer to an object and a member
        /// pointer to a function of that object's type.
        ///
        /// \tparam ptrT Pointer type.
        /// \tparam memfnT Member function pointer type.
        /// \tparam resT Result type of the member function.
        template <typename ptrT, typename memfnT, typename resT>
        class MemFuncWrapper {
        private:
            ptrT ptr_; ///< Pointer to the object.
            memfnT memfn_; ///< Member function of the object's type.

            friend memfnT get_mem_func_ptr<ptrT, memfnT, resT>(const MemFuncWrapper<ptrT, memfnT, resT>&);

        public:
            /// Alias for a wrapper that returns \c void.
            typedef MemFuncWrapper<ptrT, memfnT, void> MemFuncWrapper_;

            /// Alias for the function's result type.
            typedef resT result_type;

            /// Alias for the type of the member function pointer.
            typedef memfnT memfn_type;

            /// \brief Constructor that sets the pointer to the default value
            ///     from `DefaultInitPtr<ptrT>`.
            MemFuncWrapper()
                : ptr_(DefaultInitPtr<ptrT>::init()), memfn_()
            { }

            /// \brief Copy constructor that copies the object pointer and
            ///     member function pointer.
            MemFuncWrapper(const MemFuncWrapper_& other) :
                ptr_(other.ptr_), memfn_(other.memfn_)
            { }

            /// \brief Constructs a wrapper from an object pointer and a member
            ///     function pointer.
            ///
            /// \param[in] ptr The object pointer.
            /// \param[in] memfn The member function pointer.
            MemFuncWrapper(ptrT ptr, memfnT memfn) :
                ptr_(ptr), memfn_(memfn)
            { }

            /// Copy assignment operator.
            ///
            /// \param[in] other The wrapper to be copied.
            /// \return This wrapper, which is now a copy of \c other.
            MemFuncWrapper_& operator=(const MemFuncWrapper_& other) {
                ptr_ = other.ptr_;
                memfn_ = other.memfn_;
                return *this;
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// no arguments.
            /// \return The member function's return value.
            resT operator()() const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)();
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 1 argument.
            /// \tparam a1T The type of argument 1.
            /// \param[in,out] a1 Argument 1.
            /// \return The member function's return value.
            template <typename a1T>
            resT operator()(a1T& a1) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 2 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \return The member function's return value.
            template <typename a1T, typename a2T>
            resT operator()(a1T& a1, a2T& a2) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 3 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \return The member function's return value.
            template <typename a1T, typename a2T, typename a3T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 4 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \return The member function's return value.
            template <typename a1T, typename a2T, typename a3T, typename a4T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 5 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \tparam a5T The type of argument 5.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \param[in,out] a5 Argument 5.
            /// \return The member function's return value.
            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4, a5);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 6 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \tparam a5T The type of argument 5.
            /// \tparam a6T The type of argument 6.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \param[in,out] a5 Argument 5.
            /// \param[in,out] a6 Argument 6.
            /// \return The member function's return value.
            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 7 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \tparam a5T The type of argument 5.
            /// \tparam a6T The type of argument 6.
            /// \tparam a7T The type of argument 7.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \param[in,out] a5 Argument 5.
            /// \param[in,out] a6 Argument 6.
            /// \param[in,out] a7 Argument 7.
            /// \return The member function's return value.
            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 8 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \tparam a5T The type of argument 5.
            /// \tparam a6T The type of argument 6.
            /// \tparam a7T The type of argument 7.
            /// \tparam a8T The type of argument 8.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \param[in,out] a5 Argument 5.
            /// \param[in,out] a6 Argument 6.
            /// \param[in,out] a7 Argument 7.
            /// \param[in,out] a8 Argument 8.
            /// \return The member function's return value.
            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 9 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \tparam a5T The type of argument 5.
            /// \tparam a6T The type of argument 6.
            /// \tparam a7T The type of argument 7.
            /// \tparam a8T The type of argument 8.
            /// \tparam a9T The type of argument 9.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \param[in,out] a5 Argument 5.
            /// \param[in,out] a6 Argument 6.
            /// \param[in,out] a7 Argument 7.
            /// \param[in,out] a8 Argument 8.
            /// \param[in,out] a9 Argument 9.
            /// \return The member function's return value.
            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8, a9T& a9) const {
                MADNESS_ASSERT(ptr_);
                return ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8, a9);
            }

            /// Serializes a \c MemFuncWrapper.

            /// \tparam Archive The archive type.
            /// \param[in,out] ar The archive.
            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & ptr_ & memfn_;
            }

            friend memfnT get_mem_func_ptr(const MemFuncWrapper_& wrapper) {
                return wrapper.memfn_;
            }
        }; // class MemFuncWrapper


        /// Functor wrapper for object and member function pointers.

        /// This is a specialization of \c MemFuncWrapper for the case where
        /// the member function has \c void return type.
        /// \tparam ptrT Pointer type.
        /// \tparam memfnT Member function pointer type.
        template <typename ptrT, typename memfnT>
        class MemFuncWrapper<ptrT, memfnT, void> {
        private:
            ptrT ptr_; ///< Pointer to the object.
            memfnT memfn_; ///< Pointer to the desired member function.

            friend memfnT get_mem_func_ptr<ptrT, memfnT, void>(const MemFuncWrapper<ptrT, memfnT, void>&);

        public:
            /// Alias for a member function with \c void return type.
            typedef MemFuncWrapper<ptrT, memfnT, void> MemFuncWrapper_;

            /// Alias for the function's return type.
            typedef void result_type;

            /// Alias for the member function pointer type.
            typedef memfnT memfn_type;

            /// Construct a wrapper with the default pointer for this type.
            MemFuncWrapper()
                : ptr_(DefaultInitPtr<ptrT>::init()), memfn_()
            { }

            /// Copy constructor.

            /// \param[in] other The wrapper to copy.
            MemFuncWrapper(const MemFuncWrapper_& other) :
                ptr_(other.ptr_), memfn_(other.memfn_)
            { }

            /// \brief Construct a wrapper from an object pointer and a member
            ///     function pointer.

            /// \param[in] ptr The object pointer.
            /// \param[in] memfn The member function pointer.
            MemFuncWrapper(ptrT ptr, memfnT memfn) :
                ptr_(ptr), memfn_(memfn)
            { }

            /// Copy assignment operator.
            ///
            /// \param[in] other The wrapper to be copied.
            /// \return This wrapper, which is now a copy of \c other.
            MemFuncWrapper_& operator=(const MemFuncWrapper_& other) {
                ptr_ = other.ptr_;
                memfn_ = other.memfn_;
                return *this;
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// no arguments.
            void operator()() const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)();
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 1 argument.
            /// \tparam a1T The type of argument 1.
            /// \param[in,out] a1 Argument 1.
            template <typename a1T>
            void operator()(a1T& a1) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 2 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            template <typename a1T, typename a2T>
            void operator()(a1T& a1, a2T& a2) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 3 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            template <typename a1T, typename a2T, typename a3T>
            void operator()(a1T& a1, a2T& a2, a3T& a3) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 4 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            template <typename a1T, typename a2T, typename a3T, typename a4T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 5 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \tparam a5T The type of argument 5.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \param[in,out] a5 Argument 5.
            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4, a5);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 6 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \tparam a5T The type of argument 5.
            /// \tparam a6T The type of argument 6.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \param[in,out] a5 Argument 5.
            /// \param[in,out] a6 Argument 6.
            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 7 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \tparam a5T The type of argument 5.
            /// \tparam a6T The type of argument 6.
            /// \tparam a7T The type of argument 7.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \param[in,out] a5 Argument 5.
            /// \param[in,out] a6 Argument 6.
            /// \param[in,out] a7 Argument 7.
            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 8 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \tparam a5T The type of argument 5.
            /// \tparam a6T The type of argument 6.
            /// \tparam a7T The type of argument 7.
            /// \tparam a8T The type of argument 8.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \param[in,out] a5 Argument 5.
            /// \param[in,out] a6 Argument 6.
            /// \param[in,out] a7 Argument 7.
            /// \param[in,out] a8 Argument 8.
            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8);
            }

            /// \brief Evaluates the member function, when dereferenced from the
            ///     object pointer.

            /// This variant should only be called for a member function with
            /// 9 arguments.
            /// \tparam a1T The type of argument 1.
            /// \tparam a2T The type of argument 2.
            /// \tparam a3T The type of argument 3.
            /// \tparam a4T The type of argument 4.
            /// \tparam a5T The type of argument 5.
            /// \tparam a6T The type of argument 6.
            /// \tparam a7T The type of argument 7.
            /// \tparam a8T The type of argument 8.
            /// \tparam a9T The type of argument 9.
            /// \param[in,out] a1 Argument 1.
            /// \param[in,out] a2 Argument 2.
            /// \param[in,out] a3 Argument 3.
            /// \param[in,out] a4 Argument 4.
            /// \param[in,out] a5 Argument 5.
            /// \param[in,out] a6 Argument 6.
            /// \param[in,out] a7 Argument 7.
            /// \param[in,out] a8 Argument 8.
            /// \param[in,out] a9 Argument 9.
            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8, a9T& a9) const {
                MADNESS_ASSERT(ptr_);
                ((*ptr_).*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8, a9);
            }

            /// Serializes a \c MemFuncWrapper.

            /// \tparam Archive The archive type.
            /// \param[in,out] ar The archive.
            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & ptr_ & memfn_;
            }

        }; // class MemFuncWrapper<ptrT, memfnT, void>

        /// \brief Create a member function wrapper (\c MemFuncWrapper) from an
        ///     object and a member function pointer.

        /// \tparam objT The object type.
        /// \tparam memfnT The member function pointer type.
        /// \param[in] obj The object.
        /// \param[in] memfn The member function pointer.
        /// \return A wrapped member function pointer.
        template <typename objT, typename memfnT>
        MemFuncWrapper<objT*, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(objT& obj, memfnT memfn) {
            return MemFuncWrapper<objT*, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(& obj, memfn);
        }

        /// \brief Create a member function wrapper (\c MemFuncWrapper) from a
        ///     const object and a member function pointer.

        /// \tparam objT The object type.
        /// \tparam memfnT The member function pointer type.
        /// \param[in] obj The object.
        /// \param[in] memfn The member function pointer.
        /// \return A wrapped member function pointer.
        template <typename objT, typename memfnT>
        MemFuncWrapper<const objT*, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(const objT& obj, memfnT memfn) {
            return MemFuncWrapper<const objT*, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(& obj, memfn);
        }

        /// \brief Create a member function wrapper (\c MemFuncWrapper) from a
        ///     pointer and a member function pointer.

        /// \tparam objT The object type.
        /// \tparam memfnT The member function pointer type.
        /// \param[in] obj Pointer to the object.
        /// \param[in] memfn The member function pointer.
        /// \return A wrapped member function pointer.
        template <typename objT, typename memfnT>
        MemFuncWrapper<objT*, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(objT* obj, memfnT memfn) {
            return MemFuncWrapper<objT*, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(obj, memfn);
        }

        /// \brief Create a member function wrapper (\c MemFuncWrapper) from a
        ///     const pointer and a member function pointer.

        /// \tparam objT The object type.
        /// \tparam memfnT The member function pointer type.
        /// \param[in] obj Pointer to the object.
        /// \param[in] memfn The member function pointer.
        /// \return A wrapped member function pointer.
        template <typename objT, typename memfnT>
        MemFuncWrapper<const objT*, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(const objT* obj, memfnT memfn) {
            return MemFuncWrapper<const objT*, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(obj, memfn);
        }

        /// \brief Create a member function wrapper (\c MemFuncWrapper) from a
        ///     shared pointer and a member function pointer.

        /// \tparam objT The object type.
        /// \tparam memfnT The member function pointer type.
        /// \param[in] obj Shared pointer to the object.
        /// \param[in] memfn The member function pointer.
        /// \return A wrapped member function pointer.
        template <typename objT, typename memfnT>
        MemFuncWrapper<std::shared_ptr<objT>, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(std::shared_ptr<objT>& obj, memfnT memfn) {
            return MemFuncWrapper<std::shared_ptr<objT>, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(obj, memfn);
        }

        /// \brief Create a member function wrapper (\c MemFuncWrapper) from a
        ///     const shared pointer and a member function pointer.

        /// \tparam objT The object type.
        /// \tparam memfnT The member function pointer type.
        /// \param[in] obj Shared pointer to the object.
        /// \param[in] memfn The member function pointer.
        /// \return A wrapped member function pointer.
        template <typename objT, typename memfnT>
        MemFuncWrapper<std::shared_ptr<objT>, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(const std::shared_ptr<objT>& obj, memfnT memfn) {
            return MemFuncWrapper<std::shared_ptr<objT>, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(obj, memfn);
        }

        /// Returns the member function pointer from a wrapper.

        /// \tparam ptrT Pointer type.
        /// \tparam memfnT Member function pointer type.
        /// \tparam resT Member function return type.
        /// \param wrapper Wrapper to the member function.
        /// \return The member function pointer.
        template <typename ptrT, typename memfnT, typename resT>
        memfnT get_mem_func_ptr(const MemFuncWrapper<ptrT, memfnT, resT>& wrapper) {
            return wrapper.memfn_;
        }

    }  // namespace detail
} // namespace madness

#endif // MADNESS_WORLD_MEM_FUNC_WRAPPER_H__INCLUDED
