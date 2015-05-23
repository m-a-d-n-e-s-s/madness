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

#ifndef MADNESS_WORLD_DIST_CACHE_H__INCLUDED
#define MADNESS_WORLD_DIST_CACHE_H__INCLUDED

#include <madness/world/madness_exception.h>
#include <madness/world/worldhashmap.h>
#include <madness/world/future.h>

namespace madness {
    namespace detail {

        /// Distributed caching utility

        /// This object implements a local, key-value caching mechanism that can
        /// be used by remote tasks or active messages to move data between
        /// processes without synchronization. Because cache values are
        /// retrieved via a \c Future, you can get or set the cache in any
        /// order. The first call \c get_cache_value or \c set_cache_value or
        /// will insert the cache element, and the second call to these
        /// functions will remove it. Therefore, \c set_cache_value and
        /// get_cache_value can only be called once each per cache value.
        /// \tparam keyT The key type of the cache
        template <typename keyT>
        class DistCache {
        private:

            // Forward declarations
            class Cache;
            template <typename> class CacheData;

            typedef madness::ConcurrentHashMap<keyT, Cache*> cache_container;
            ///< The container that holds cache values
            typedef typename cache_container::datumT datum_type;
            ///< Cache container datum type

            static cache_container caches_; ///< Cache container

            /// Cache interface class

            /// This base class is used to access derived class values
            class Cache {
            public:

                /// Virtual destructor
                virtual ~Cache() { }

                /// Cache data accessor

                /// \tparam valueT The cached data type
                /// \return A const reference to the cached future
                template <typename valueT>
                const madness::Future<valueT>& get() const {
                    MADNESS_ASSERT(this->get_type_info() == typeid(CacheData<valueT>));
                    return static_cast<const CacheData<valueT>*>(this)->value();
                }

            private:

                /// Typeid accessor of the derived class

                /// \return The std::type_info of the derived class
                virtual const std::type_info& get_type_info() const = 0;

            }; // class Cache

            /// Cache value container

            /// \tparam valueT The data type stored in the cache
            template <typename valueT>
            class CacheData : public Cache {
            private:
                madness::Future<valueT> value_; ///< Local cached data

            public:

                /// Default constructor
                CacheData() : value_() { }

                /// Constructor with future initialization
                CacheData(const madness::Future<valueT>& value) : value_(value) { }

                /// Constructor with data initialization
                CacheData(const valueT& value) : value_(value) { }

                /// Virtual destructor
                virtual ~CacheData() { }

                /// Data accessor

                /// \return A const reference to the data
                const madness::Future<valueT>& value() const { return value_; }

            private:

                /// Typeid accessor of the derived class

                /// \return The std::type_info of the derived class
                virtual const std::type_info& get_type_info() const {
                    return typeid(CacheData<valueT>);
                }

            }; // class CacheData

            public:

            /// Set the cache value accosted with \c key

            /// This will set the value associated with \c key to \c value. If
            /// the cache element does not exist, it is inserted into the cache.
            /// Otherwise, it is removed from the cache.
            /// \tparam valueT The object type that will be used to set the
            /// cache (may be a \c madness::Future type).
            /// \param key The key associated with \c value
            /// \param value The data that will be cached
            template <typename valueT>
            static void set_cache_value(const keyT& key, const valueT& value) {
                typedef typename madness::remove_future<valueT>::type value_type;

                // Retrieve the cached future
                typename cache_container::accessor acc;
                if(caches_.insert(acc, datum_type(key, static_cast<Cache*>(nullptr)))) {

                    // A new element was inserted, so create a new cache object.
                    acc->second = new CacheData<value_type>(value);
                    acc.release();

                } else {

                    // The element already existed, so retrieve the data
                    Cache* cache = acc->second;
                    caches_.erase(acc);

                    // Set the cache value
                    madness::Future<value_type> f =
                            cache->template get<value_type>();
                    MADNESS_ASSERT(! f.probe());
                    f.set(value);

                    // Cleanup cache
                    delete cache;
                }
            }

            /// Get the cache value accosted with \c key

            /// This will get the value associated with \c key to \c value. The
            /// value is given in the form of a \c Future, which is set by a
            /// call to \c set_cache_value. If the cache element does not exist,
            /// it is inserted into the cache. Otherwise, it is removed from the
            /// cache.
            /// \tparam valueT The object type that will be used to set the cache
            /// \param[in] key The key associated with \c value
            /// \param[out] value The data that will be cached
            template <typename valueT>
            static void get_cache_value(const keyT& key, madness::Future<valueT>& value) {
                // Retrieve the cached future
                typename cache_container::accessor acc;
                if(caches_.insert(acc, datum_type(key, static_cast<Cache*>(nullptr)))) {
                    // A new element was inserted, so create a new cache object.
                    acc->second = new CacheData<valueT>(value);
                    acc.release();
                } else {
                    // The element already existed, so retrieve the data and
                    // remove the cache element.
                    Cache* cache = acc->second;
                    caches_.erase(acc);

                    // Get the result
                    value.set(cache->template get<valueT>());
                    delete cache;
                }
            }

            /// Get the cache value accosted with \c key

            /// This will get the value associated with \c key to \c value. If
            /// the cache element does not exist, it is inserted into the cache.
            /// Otherwise, it is removed from the cache.
            /// \tparam valueT The object type that will be used to set the cache
            /// \param[in] key The key associated with \c value
            /// \return A \c Future that holds/will hold the cache value, which
            /// will be set by a call to \c set_cache_value.
            template <typename valueT>
            static madness::Future<valueT> get_cache_value(const keyT& key) {
                // Retrieve the cached future
                typename cache_container::accessor acc;
                if(caches_.insert(acc, datum_type(key, static_cast<Cache*>(nullptr)))) {
                    // A new element was inserted, so create a new cache object.
                    acc->second = new CacheData<valueT>();
                    madness::Future<valueT> value(acc->second->template get<valueT>());
                    acc.release();

                    return value;
                } else {
                    // The element already existed, so retrieve the data and
                    // remove the cache element.
                    Cache* cache = acc->second;
                    caches_.erase(acc);

                    // Get the result
                    madness::Future<valueT> value(cache->template get<valueT>());
                    delete cache;

                    return value;
                }
            }

        }; // class DistCache

        template <typename keyT>
        typename DistCache<keyT>::cache_container DistCache<keyT>::caches_;

    }  // namespace detail
} // namespace madness

#endif // MADNESS_WORLD_DIST_CACHE_H__INCLUDED
