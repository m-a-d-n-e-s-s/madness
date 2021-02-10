/*
  This file is part of MADNESS.

  Copyright (C) 2013  Virginia Tech

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

#ifndef MADNESS_WORLD_DISTRIBUTED_ID_H__INCLUDED
#define MADNESS_WORLD_DISTRIBUTED_ID_H__INCLUDED

/**
 \file distributed_id.h
 \brief \todo Description needed.
 \ingroup parallel_runtime
*/

#include <madness/world/uniqueid.h>
#include <madness/world/worldtypes.h>
#include <madness/world/print.h>

namespace madness {

    /// Distributed ID which is used to identify objects.
    typedef std::pair<uniqueidT, std::size_t> DistributedID;

    /// Distributed ID equality comparison operator.

    /// \param[in] left The first key to compare.
    /// \param[in] right The second key to compare.
    /// \return True when \c first and \c second of \c left and \c right are
    ///    equal; otherwise false.
    inline bool operator==(const DistributedID& left, const DistributedID& right) {
        return (left.first == right.first) && (left.second == right.second);
    }

    /// Distributed ID inequality comparison operator.

    /// \param[in] left The first key to compare.
    /// \param[in] right The second key to compare.
    /// \return True when \c first or \c second of \c left and \c right are
    ///    not equal; otherwise false.
    inline bool operator!=(const DistributedID& left, const DistributedID& right) {
        return (left.first != right.first) || (left.second != right.second);
    }

    /// Overload redirect to std::ostream to be discoverable via ADL
    inline std::ostream& operator<<(std::ostream& os, const DistributedID& did) {
      return madness::operators::operator<<(os, did);
    }

    /// Key object that includes the process information.

    /// \tparam Key The base key type.
    /// \tparam Tag A type to differentiate key types.
    template <typename Key, typename Tag = void>
    class ProcessKey {
    private:
        Key key_; ///< The base key type.
        ProcessID proc_; ///< The process that generated the key.

    public:

        /// Default constructor.
        ProcessKey() : key_(), proc_(-1) { }

        /// Constructor.

        /// \param[in] key The base key.
        /// \param[in] proc The process that generated the key.
        ProcessKey(const Key& key, const ProcessID proc) :
            key_(key), proc_(proc)
        { }

        /// Copy constructor.

        /// \param[in] other The key to be copied.
        ProcessKey(const ProcessKey<Key, Tag>& other) :
            key_(other.key_), proc_(other.proc_)
        { }

        /// Copy assignment operator.

        /// \param[in] other The key to be copied.
        /// \return A reference to this object.
        ProcessKey<Key, Tag>& operator=(const ProcessKey<Key, Tag>& other) {
            key_ = other.key_;
            proc_ = other.proc_;
            return *this;
        }

        /// Base key accessor.

        /// \return The base key.
        const Key& key() const { return key_; }

        /// Process ID accessor.

        /// \return The process ID.
        ProcessID proc() const { return proc_; }

        /// Equality comparison.

        /// \param[in] other The key to be compared to this.
        /// \return True when other key and other process are equal to that of
        ///    this key; otherwise false.
        bool operator==(const ProcessKey<Key, Tag>& other) const {
            return ((key_ == other.key_) && (proc_ == other.proc_));
        }

        /// Inequality comparison.

        /// \param[in] other The key to be compared to this.
        /// \return True when other key or other process are not equal to that
        ///    of this key; otherwise false.
        bool operator!=(const ProcessKey<Key, Tag>& other) const {
            return ((key_ != other.key_) || (proc_ != other.proc_));
        }

        /// Serialize this key.

        /// \tparam Archive The archive type.
        /// \param[in,out] ar The archive object that will serialize this object.
        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & key_ & proc_;
        }

        /// Hashing function.

        /// \param[in] key The key to be hashed.
        /// \return The hashed key value.
        friend hashT hash_value(const ProcessKey<Key, Tag>& key) {
            Hash<Key> hasher;
            hashT seed = hasher(key.key_);
            madness::detail::combine_hash(seed, key.proc_);
            return seed;
        }

    }; // class ProcessKey

    /// Key object that uses a tag to differentiate keys.

    /// \tparam Key The base key type.
    /// \tparam Tag A type to differentiate key types.
    template <typename Key, typename Tag>
    class TaggedKey {
    private:
        Key key_; ///< The base key type.

    public:

        /// Default constructor.
        TaggedKey() : key_() { }

        /// Constructor.

        /// \param[in] key The base key.
        TaggedKey(const Key& key) : key_(key) { }

        /// Copy constructor.

        /// \param[in] other The key to be copied.
        TaggedKey(const TaggedKey<Key, Tag>& other) : key_(other.key_) { }

        /// Copy assignment operator.

        /// \param[in] other The key to be copied.
        /// \return A reference to this object.
        TaggedKey<Key, Tag>& operator=(const TaggedKey<Key, Tag>& other) {
            key_ = other.key_;
            return *this;
        }

        /// Base key accessor.

        /// \return The base key.
        const Key& key() const { return key_; }

        /// Equality comparison.

        /// \param[in] other The key to be compared to this.
        /// \return True when other key and other process are equal to that of
        ///    this key; otherwise false.
        bool operator==(const TaggedKey<Key, Tag>& other) const {
            return (key_ == other.key_);
        }

        /// Inequality comparison.

        /// \param[in] other The key to be compared to this.
        /// \return True when other key or other process are not equal to that
        ///    of this key; otherwise false.
        bool operator!=(const TaggedKey<Key, Tag>& other) const {
            return (key_ != other.key_);
        }

        /// Serialize this key.

        /// \tparam Archive The archive type.
        /// \param[in,out] ar The archive object that will serialize this object.
        template <typename Archive>
        void serialize(const Archive& ar) { ar & key_; }

        /// Hashing function.

        /// \param[in] key The key to be hashed.
        /// \return The hashed key value.
        friend hashT hash_value(const TaggedKey<Key, Tag>& key) {
            Hash<Key> hasher;
            return hasher(key.key_);
        }

    }; // class TagKey

} // namespace madness

namespace std {

    /// Hash a \c DistributedID.

    /// \param[in] id The distributed ID to be hashed.
    /// \return The hash value of \c id.
    inline madness::hashT hash_value(const madness::DistributedID& id) {
        madness::hashT seed = madness::hash_value(id.first);
        madness::detail::combine_hash(seed, madness::hash_value(id.second));

        return seed;
    }

} // namespace std

#endif // MADNESS_WORLD_DISTRIBUTED_ID_H__INCLUDED
