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
 \file worldptr.h
 \brief The \c madness::detail::WorldPtr class for global pointers.
 \ingroup world
*/

#ifndef MADNESS_WORLD_WORLDPTR_H__INCLUDED
#define MADNESS_WORLD_WORLDPTR_H__INCLUDED

#include <madness/world/madness_exception.h>
#include <madness/world/worldtypes.h>   // for ProcessID
#include <madness/world/archive.h>      // for wrap_opaque
#include <madness/world/world.h>
#include <algorithm>            // for std::swap
#include <iostream>             // for std::iostream

/// \addtogroup world
/// @{

namespace madness {

    namespace detail {

        /// \todo Brief description needed.

        /// \todo Description needed.
        /// \tparam U Description needed.
        template<typename U>
        struct ptr_traits {
            /// \todo Brief description needed.
            typedef U & reference;
        };

        /// Specialization of \c ptr_traits for type \c void.
        template<>
        struct ptr_traits<void> {
            /// \todo Brief description needed.
            typedef void reference;
        };

        /// A global pointer address, valid anywhere in the world.

        /// Stores a globally addressable pointer. It can be sent to any
        /// process in the world.
        /// \tparam T The pointer type.
        template <typename T>
        class WorldPtr {
        public:
            typedef unsigned long worldidT; ///< World ID type.

        private:
            World* world_; ///< A pointer to the world.
            worldidT worldid_; ///< The world ID.
            ProcessID rank_; ///< The rank of the node that the pointer belongs to.
            T* pointer_; ///< The pointer being referenced.

            template<typename>
            friend class WorldPtr;

            /// Current local rank.

            /// \return The rank of the current node. If the pointer is not
            ///     set, then -2.
            /// \note -2 is returned so it is not equal to the null value of -1.
            ProcessID local_rank() const {
                return (world_ != nullptr ? world_->rank() : -2);
            }

        public:

            /// Alias for the pointer type.
            typedef T* pointer;

            /// \todo Brief description needed.
            typedef typename ptr_traits<T>::reference reference;


            /// Default constructor

            /// Creates a \c NULL pointer. There is no owner; i.e. the owner is
            /// set to -1.
            /// \todo Would it be worth adding a static constant \c ProcessID equal to -1 to signify an unowned pointer?
            WorldPtr() :
                world_(nullptr),
                worldid_(0),
                rank_(-1),
                pointer_(nullptr)
            { }


            /// %World pointer constructor.

            /// Construct a world pointer form a local pointer.
            /// \param[in] w A reference to the local world.
            /// \param[in] p The local pointer.
            WorldPtr(World& w, T* p) :
                world_(&w),
                worldid_(w.id() + 1),
                rank_(w.rank()),
                pointer_(p)
            { }


            /// Copy constructor

            /// \param[in] other The world pointer to be copied.
            WorldPtr(const WorldPtr<T>& other) :
                world_(other.world_),
                worldid_(other.worldid_),
                rank_(other.rank_),
                pointer_(other.pointer_)
            { }


            /// Copy conversion constructor.

            /// Copy and convert a pointer from \c U* to \c T* type.
            /// \tparam U The pointer type of the \c other pointer.
            /// \param[in] other The world pointer to be copied.
            /// \note \c U* must be implicitly convertible to \c T* type.
            template <typename U>
            WorldPtr(const WorldPtr<U>& other) :
                world_(other.world_),
                worldid_(other.worldid_),
                rank_(other.rank_),
                pointer_(other.pointer_)
            { }


            /// Copy assignment operator.

            /// \param[in] other The world pointer to be copied.
            /// \return A reference to this object.
            WorldPtr<T>& operator=(const WorldPtr<T>& other) {
                world_ = other.world_;
                worldid_ = other.worldid_;
                rank_ = other.rank_;
                pointer_ = other.pointer_;

                return *this;
            }


            /// Copy conversion assignment operator.

            /// Copy and convert a pointer from \c U* to \c T* type.
            /// \tparam U The pointer type of the \c other pointer.
            /// \param[in] other The world pointer to be copied.
            /// \return A reference to this object.
            /// \note \c U* must be implicitly convertible to \c T* type.
            template <typename U>
            WorldPtr<T>& operator=(const WorldPtr<U>& other) {
                world_ = other.world_;
                worldid_ = other.worldid_;
                rank_ = other.rank_;
                pointer_ = other.pointer_;

                return *this;
            }


            /// Check that the world pointer references a local pointer.

            /// \return True if the pointer points to a local address; false
            ///     if it points to a remote address or is NULL.
            bool is_local() const {
                return local_rank() == rank_;
            }


            /// Check that the world pointer has an owner.

            /// \return True if the pointer has a valid owner; false otherwise.
            bool has_owner() const {
                return (rank_ != -1) && (world_ != nullptr);
            }


            /// Pointer accessor.

            /// Get the pointer from the world pointer.
            /// \note A default initialized pointer is not considered to be
            /// local because it is not associated with a world.
            /// \return The local pointer.
            /// \throw MadnessException When the pointer references a remote
            ///     address.
            pointer get() const {
                // It is not safe to access this pointer remotely unless null.
                MADNESS_ASSERT(is_local());
                return pointer_;
            }


            /// Dereference operator.

            /// Dereference the local pointer.
            /// \return A reference to the local pointer.
            /// \throw MadnessException If the pointer references a remote
            ///     address, or if the pointer is \c NULL.
            reference operator*() const {
                // It is not safe to access a NULL pointer with this operator.
                MADNESS_ASSERT(pointer_ != nullptr);
                // It is not safe to access this pointer remotely.
                MADNESS_ASSERT(is_local());
                return *pointer_;
            }


            /// Pointer arrow operator.

            /// Access members of the pointer.
            /// \return The local pointer.
            /// \throw MadnessException If the pointer references a remote
            ///     address, or if the pointer is \c NULL.
            pointer operator->() const {
                // It is not safe to access a NULL pointer with this operator.
                MADNESS_ASSERT(pointer_ != nullptr);
                // It is not safe to access this pointer remotely.
                MADNESS_ASSERT(is_local());
                return pointer_;
            }


            /// Boolean conversion operator.

            /// \return True if the pointer is non-null; false otherwise.
            operator bool () const {
                return pointer_;
            }


            /// Boolean conversion (not) operator.

            /// \return True if the pointer is null; false otherwise.
            bool operator ! () const {
                return !pointer_;
            }


            /// Equality comparison operator.

            /// \tparam U Another pointer type.
            /// \param[in] other The pointer to compare with.
            /// \return True if the pointers refer to the same address from
            ///     the same node in the same world; false otherwise.
            template <typename U>
            bool operator==(const WorldPtr<U>& other) const {
                return (pointer_ == other.pointer_) && (rank_ == other.rank_)
                    && (worldid_ == other.worldid_);
            }


            /// Inequality comparison operator.

            /// \tparam U Another pointer type.
            /// \param[in] other The other pointer to compare with.
            /// \return True if the pointers refer to different addresses or
            ///     different nodes or different worlds; false otherwise.
            template <typename U>
            bool operator!=(const WorldPtr<U>& other) const {
                return (pointer_ != other.pointer_) || (rank_ != other.rank_)
                    || (worldid_ != other.worldid_);
            }


            /// Less-than comparison operator.

            /// This operator does a lexicographical comparison of world ID,
            /// rank, and pointer (in that order).
            /// \tparam U Another pointer type.
            /// \param[in] other The other pointer to compare with.
            /// \return True if the lexicographical comparison of world ID,
            ///     rank, and pointer is true; false otherwise.
            template <typename U>
            bool operator<(const WorldPtr<U>& other) const {
                return (worldid_ < other.worldid_) ||
                        ((worldid_ == other.worldid_) && ((rank_ < other.rank_) ||
                        ((rank_ == other.rank_) && (pointer_ < other.pointer_))));
            }


            /// World accessor.

            /// \return A reference to the world (may be \c NULL).
            /// \throw MadnessException When the pointer world has not been set
            ///     (i.e. when \c has_owner()==false).
            World& get_world() const {
                MADNESS_ASSERT(world_ != nullptr);
                return *world_;
            }


            /// %World ID accessor.

            /// \return The world ID of the world that the pointer belongs to.
            worldidT get_worldid() const {
                return worldid_ - 1;
            }


            /// Rank accessor.

            /// \todo Finish this sentence: If the pointer is not associated with
            /// \return The rank of the process that owns the pointer.
            ProcessID owner() const {
                return rank_;
            }


            /// Swap the content of \c this with \c other.

            /// \tparam U The other world pointer type.
            /// \param[in,out] other The other world pointer.
            /// \note \c U* must be implicitly convertible to \c T* type.
            template <typename U>
            void swap(WorldPtr<U>& other) {
                std::swap(world_, other.world_);
                std::swap(worldid_, other.worldid_);
                std::swap(rank_, other.rank_);
                std::swap(pointer_, other.pointer_);
            }


            /// Deserialize the world pointer.

            /// Deserialize the world pointer for remote communication or write
            /// to disk.
            /// \tparam Archive The archive object type.
            /// \param[in] ar The archive.
            template <class Archive>
            inline void load_internal_(const Archive& ar) {
                ar & worldid_ & rank_ & archive::wrap_opaque(pointer_);
                world_ = (worldid_ != 0 ? World::world_from_id(get_worldid()) : nullptr);
            }

            /// Serialize the world pointer.

            /// Serialize the world pointer for remote communication or write
            /// to disk.
            /// \tparam Archive The archive object type.
            /// \param[in] ar The archive.
            template <class Archive>
            inline void store_internal_(const Archive& ar) const {
                ar & worldid_ & rank_ & archive::wrap_opaque(pointer_);
            }

            /// Output stream insertion operator for world pointers.

            /// \param[in,out] out The output stream.
            /// \param[in] p The world pointer.
            /// \return The output stream.
            /// \todo Does this \c friend function need to be implemented in the class or can we move it to a \c *.cc file?
            friend std::ostream& operator<<(std::ostream& out, const WorldPtr<T>& p) {
                out << "WorldPointer(ptr=" << p.pointer_ << ", rank=";
                if(p.rank_ >= 0)
                    out << p.rank_;
                else
                    out << "none";
                out << ", worldid=";
                if(p.worldid_ != 0)
                    out << (p.worldid_ - 1);
                else
                    out << "none";
                out << ")";
                return out;
            }
        }; // class WorldPtr


        /// Swap the content of \c l with \c r.

        /// \tparam T The world pointer type.
        /// \param[in,out] l One world pointer.
        /// \param[in,out] r The other world pointer.
        template <typename T>
        void swap(WorldPtr<T>& l, WorldPtr<T>& r) {
            l.swap(r);
        }


        /// Greater-than comparison operator.

        /// This operator does a lexicographical comparison of world ID,
        /// rank, and pointer (in that order).
        /// \tparam T One pointer type.
        /// \tparam U Another pointer type.
        /// \param[in] left One pointer.
        /// \param[in] right The other pointer.
        /// \return True if the lexicographical comparison of world ID,
        ///     rank, and pointer is true; false otherwise.
        template <typename T, typename U>
        bool operator>(const WorldPtr<T>& left, const WorldPtr<U>& right) {
            return right < left;
        }


        /// Less-than-equal-to comparison operator.

        /// This operator does a lexicographical comparison of world ID,
        /// rank, and pointer (in that order).
        /// \tparam T One pointer type.
        /// \tparam U Another pointer type.
        /// \param[in] left One pointer.
        /// \param[in] right The other pointer.
        /// \return True if the lexicographical comparison of world ID,
        ///     rank, and pointer is true; false otherwise.
        template <typename T, typename U>
        bool operator<=(const WorldPtr<T>& left, const WorldPtr<U>& right) {
            return !(right < left);
        }


        /// Greater-than-equal-to comparison operator.

        /// This operator does a lexicographical comparison of world ID,
        /// rank, and pointer (in that order).
        /// \tparam T One pointer type.
        /// \tparam U Another pointer type.
        /// \param[in] left One pointer.
        /// \param[in] right The other pointer.
        /// \return True if the lexicographical comparison of world ID,
        ///     rank, and pointer is true; false otherwise.
        template <typename T, typename U>
        bool operator>=(const WorldPtr<T>& left, const WorldPtr<U>& right) {
            return !(left < right);
        }

    } // namespace detail

    namespace archive {

        /// Specialization of \c ArchiveLoadImpl for world pointers.

        /// \tparam Archive The archive type.
        /// \tparam T The world pointer type.
        template <typename Archive, typename T>
        struct ArchiveLoadImpl<Archive, madness::detail::WorldPtr<T> > {

            /// Load a world pointer from the archive.

            /// \param[in] ar The archive.
            /// \param[out] p The world pointer.
            static inline void load(const Archive& ar, madness::detail::WorldPtr<T>& p) {
                p.load_internal_(ar);
            }
        };

        /// Specialization of \c ArchiveStoreImpl for world pointers.

        /// \tparam Archive The archive type.
        /// \tparam T The world pointer type.
        template <typename Archive, typename T>
        struct ArchiveStoreImpl<Archive, madness::detail::WorldPtr<T> > {

            /// Store a world pointer from the archive.

            /// \param[in] ar The archive.
            /// \param[in] p The world pointer.
            static inline void store(const Archive& ar, const madness::detail::WorldPtr<T>& p) {
                p.store_internal_(ar);
            }
        };

    } // namespace archive

} // namespace madness

/// @}

#endif // MADNESS_WORLD_WORLDPTR_H__INCLUDED
