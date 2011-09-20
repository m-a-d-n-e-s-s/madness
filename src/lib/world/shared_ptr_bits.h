#ifndef MADNESS_SHARED_PTR_BITS_H
#define MADNESS_SHARED_PTR_BITS_H

#include <world/atomicint.h>
#include <world/worldexc.h>
#include <world/type_traits.h>

//  shared_ptr.hpp
//
//  (C) Copyright Greg Colvin and Beman Dawes 1998, 1999.
//  Copyright (c) 2001-2008 Peter Dimov
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
//  See http://www.boost.org/libs/smart_ptr/shared_ptr.htm for documentation.
//
//  Ported by JC into madness 9/20/2011
//  X86 assembly replaced by madness::AtomicInt by RJH 9/20/2011

#include <algorithm>            // for std::swap
#include <functional>           // for std::less
#include <typeinfo>             // for std::bad_cast
#include <cstddef>              // for std::size_t
#include <iosfwd>               // for std::basic_ostream

namespace madness {

    // hash_value

    template <class T> struct Hash;

    namespace tr1 {

        namespace shptr {

        template <class T> class shared_ptr;
        template <class T> class weak_ptr;
        template <class T> class enable_shared_from_this;
        template <class T> class enable_shared_from_this2;

        struct bad_weak_ptr : public std::exception { };

        namespace detail {

            template <class Y, class T> struct sp_convertible {
                typedef char (&yes)[1];
                typedef char (&no)[2];

                static yes f(T*);
                static no f(... );

                enum _vt {
                    value = sizeof((f)(static_cast<Y*>(0))) == sizeof(yes)
                };
            };

            struct sp_empty { };

            template <bool> struct enable_if_convertible_impl;

            template <> struct enable_if_convertible_impl<true> { typedef sp_empty type; };

            template <> struct enable_if_convertible_impl<false> { };

            template <class Y, class T> struct enable_if_convertible :
                public enable_if_convertible_impl< sp_convertible<Y, T>::value>
            { };

            struct static_cast_tag { };
            struct const_cast_tag { };
            struct dynamic_cast_tag { };

            struct sp_nothrow_tag { };

            // Shared pointer reference trait types for handling void pointers
            template <class T> struct shared_ptr_traits { typedef T & reference; };
            template <> struct shared_ptr_traits<void> { typedef void reference; };
            template <> struct shared_ptr_traits<void const> { typedef void reference; };
            template <> struct shared_ptr_traits<void volatile> { typedef void reference; };
            template <> struct shared_ptr_traits<void const volatile> { typedef void reference; };

            // enable_shared_from_this support
            template <class X, class Y, class T>
            inline void sp_enable_shared_from_this(shared_ptr<X> const * ppx, Y const * py, enable_shared_from_this<T> const * pe) {
                if(pe != 0)
                    pe->_internal_accept_owner(ppx, const_cast<Y*>(py));
            }

            template <class X, class Y, class T>
            inline void sp_enable_shared_from_this( shared_ptr<X> * ppx, Y const * py, enable_shared_from_this2<T> const * pe) {
                if(pe != 0)
                    pe->_internal_accept_owner(ppx, const_cast<Y*>(py));
            }

            inline void sp_enable_shared_from_this(... ) { }


            class sp_counted_base
            {
            private:

                sp_counted_base( sp_counted_base const & );
                sp_counted_base & operator= ( sp_counted_base const & );

                madness::AtomicInt use_count_;        // #shared
                madness::AtomicInt weak_count_;       // #weak + (#shared != 0)

            public:

                sp_counted_base() 
                {
                    use_count_ = 1;
                    weak_count_ = 1;
                }

                virtual ~sp_counted_base() // nothrow
                {
                }

                // dispose() is called when use_count_ drops to zero, to release
                // the resources managed by *this.

                virtual void dispose() = 0; // nothrow

                // destroy() is called when weak_count_ drops to zero.

                virtual void destroy() // nothrow
                {
                    delete this;
                }

                virtual void * get_deleter( std::type_info const & ti ) = 0;

                void add_ref_copy()
                {
                    use_count_++;
                }

                bool add_ref_lock() // true on success
                {
                    //return atomic_conditional_increment( &use_count_ ) != 0;
                    throw "shared_ptr::add_ref_lock - not supported";
                }

                void release() // nothrow
                {
                    // if (atomic_exchange_and_add( &use_count_, -1 ) == 1 )
                    if (use_count_.dec_and_test())
                    {
                        dispose();
                        weak_release();
                    }
                }

                void weak_add_ref() { weak_count_++; }

                void weak_release() {
                    //if( atomic_exchange_and_add( &weak_count_, -1 ) == 1 )
                    if( weak_count_.dec_and_test())
                        destroy();
                }

                long use_count() const { return use_count_; }
            };

            template<class X>
            class sp_counted_impl_p : public sp_counted_base {
            private:

                X * px_;

                sp_counted_impl_p( sp_counted_impl_p const & );
                sp_counted_impl_p & operator= ( sp_counted_impl_p const & );

                typedef sp_counted_impl_p<X> this_type;

            public:

                explicit sp_counted_impl_p( X * px ): px_( px ) { }

                virtual void dispose() { madness::detail::checked_delete( px_ ); }

                virtual void * get_deleter( std::type_info const & ) { return 0; }
            };

            template<class Y, class D>
            class sp_counted_impl_pd : public sp_counted_base {
            private:

                Y* ptr; // copy constructor must not throw
                D del; // copy constructor must not throw

                sp_counted_impl_pd( sp_counted_impl_pd const & );
                sp_counted_impl_pd & operator= ( sp_counted_impl_pd const & );

                typedef sp_counted_impl_pd<Y, D> this_type;

            public:

                // pre: d(p) must not throw

                sp_counted_impl_pd( Y* p, D d ): ptr( p ), del( d ) { }

                sp_counted_impl_pd( Y* p ): ptr( p ), del() { }

                virtual void dispose()  { del( ptr ); }

                virtual void * get_deleter( std::type_info const & ti ) {
                    return ti == typeid(D)? &reinterpret_cast<char&>( del ): 0;
                }
            };

            class weak_count;

            class shared_count
            {
            private:

                sp_counted_base * pi_;

                friend class weak_count;

            public:

                shared_count(): pi_(0) { }

                template<class Y> explicit shared_count( Y * p ): pi_( 0 )
                {
                    try {
                        pi_ = new sp_counted_impl_p<Y>( p );
                    } catch(...) {
                        madness::detail::checked_delete( p );
                        throw;
                    }
                }

                template<class Y, class D> shared_count( Y * p, D d ): pi_(0)
                {
                    try {
                        pi_ = new sp_counted_impl_pd<Y, D>(p, d);
                    } catch(...) {
                        d(p); // delete p
                        throw;
                    }
                }

                ~shared_count()  {
                    if( pi_ != 0 ) pi_->release();
                }

                shared_count(shared_count const & r): pi_(r.pi_) {
                    if( pi_ != 0 ) pi_->add_ref_copy();
                }

                explicit shared_count(weak_count const & r); // throws bad_weak_ptr when r.use_count() == 0

                shared_count & operator= (shared_count const & r) {
                    sp_counted_base * tmp = r.pi_;

                    if( tmp != pi_ )
                    {
                        if( tmp != 0 ) tmp->add_ref_copy();
                        if( pi_ != 0 ) pi_->release();
                        pi_ = tmp;
                    }

                    return *this;
                }

                void swap(shared_count & r) {
                    sp_counted_base * tmp = r.pi_;
                    r.pi_ = pi_;
                    pi_ = tmp;
                }

                long use_count() const { return pi_ != 0? pi_->use_count(): 0; }

                bool unique() const  { return use_count() == 1; }

                bool empty() const { return pi_ == 0; }

                friend inline bool operator==(shared_count const & a, shared_count const & b) {
                    return a.pi_ == b.pi_;
                }

                friend inline bool operator<(shared_count const & a, shared_count const & b) {
                    return std::less<sp_counted_base *>()( a.pi_, b.pi_ );
                }

                void * get_deleter( std::type_info const & ti ) const {
                    return pi_? pi_->get_deleter( ti ): 0;
                }
            };


            class weak_count {
            private:

                sp_counted_base * pi_;

                friend class shared_count;

            public:

                weak_count() : pi_(0) { }

                weak_count(shared_count const & r) : pi_(r.pi_) {
                    if(pi_ != 0) pi_->weak_add_ref();
                }

                weak_count(weak_count const & r) : pi_(r.pi_) {
                    if(pi_ != 0) pi_->weak_add_ref();
                }

                ~weak_count() { if(pi_ != 0) pi_->weak_release(); }

                weak_count & operator= (shared_count const & r) {
                    sp_counted_base * tmp = r.pi_;

                    if( tmp != pi_ ) {
                        if(tmp != 0) tmp->weak_add_ref();
                        if(pi_ != 0) pi_->weak_release();
                        pi_ = tmp;
                    }

                    return *this;
                }

                weak_count & operator= (weak_count const & r) {
                    sp_counted_base * tmp = r.pi_;

                    if( tmp != pi_ ) {
                        if(tmp != 0) tmp->weak_add_ref();
                        if(pi_ != 0) pi_->weak_release();
                        pi_ = tmp;
                    }

                    return *this;
                }

                void swap(weak_count & r) {
                    sp_counted_base * tmp = r.pi_;
                    r.pi_ = pi_;
                    pi_ = tmp;
                }

                long use_count() const { return pi_ != 0 ? pi_->use_count() : 0; }

                bool empty() const { return pi_ == 0; }

                friend inline bool operator==(weak_count const & a, weak_count const & b) {
                    return a.pi_ == b.pi_;
                }

                friend inline bool operator<(weak_count const & a, weak_count const & b) {
                    return std::less<sp_counted_base *>()(a.pi_, b.pi_);
                }
            };

            inline shared_count::shared_count( weak_count const & r ) : pi_( r.pi_ ) {
                if( pi_ == 0 || !pi_->add_ref_lock() )
                    throw bad_weak_ptr();
            }


        } // namespace detail

        //
        //  shared_ptr
        //
        //  An enhanced relative of scoped_ptr with reference counted copy semantics.
        //  The object pointed to is deleted when the last shared_ptr pointing to it
        //  is destroyed or reset.
        //
        template <class T> class shared_ptr {
        private:
            typedef shared_ptr<T> this_type;

        public:

            typedef T element_type;
            typedef T value_type;
            typedef T * pointer;
            typedef typename detail::shared_ptr_traits<T>::reference reference;

            shared_ptr() :
                    px(0), pn() // never throws in 1.30+
            { }

            template <class Y>
            explicit shared_ptr(Y * p) :
                    px(p), pn(p) // Y must be complete
            {
                detail::sp_enable_shared_from_this(this, p, p);
            }

            //
            // Requirements: D's copy constructor must not throw
            //
            // shared_ptr will release p by calling d(p)
            //

            template <class Y, class D> shared_ptr(Y * p, D d) :
                    px(p), pn(p, d)
            {
                detail::sp_enable_shared_from_this(this, p, p);
            }

            //  generated copy constructor, destructor are fine

            template <class Y>
            explicit shared_ptr(weak_ptr<Y> const & r) :
                    pn(r.pn) // may throw
            {
                // it is now safe to copy r.px, as pn(r.pn) did not throw
                px = r.px;
            }

            template <class Y>
            shared_ptr(shared_ptr<Y> const & r,
                    typename detail::enable_if_convertible<Y, T>::type =
                             detail::sp_empty()) :
                    px(r.px), pn(r.pn) // never throws
            { }

            // aliasing
            template <class Y>
            shared_ptr(shared_ptr<Y> const & r, T * p) :
                    px(p), pn(r.pn) // never throws
            { }

            template <class Y>
            shared_ptr(shared_ptr<Y> const & r, detail::static_cast_tag) :
                    px(static_cast<element_type *>(r.px)), pn(r.pn)
            { }

            template <class Y>
            shared_ptr(shared_ptr<Y> const & r, detail::const_cast_tag) :
                    px(const_cast<element_type *>(r.px)), pn(r.pn)
            { }

            template <class Y>
            shared_ptr(shared_ptr<Y> const & r, detail::dynamic_cast_tag) :
                    px(dynamic_cast<element_type *>(r.px)), pn(r.pn)
            {
                if(px == 0) // need to allocate new counter -- the cast failed
                    pn = detail::shared_count();
            }

            // assignment

            shared_ptr & operator=(shared_ptr const & r) {
                this_type(r).swap(*this);
                return *this;
            }

            void reset() { this_type().swap(*this); }

            template <class Y>
            void reset(Y * p) { // Y must be complete
                MADNESS_ASSERT(p == 0 || p != px); // catch self-reset errors
                this_type(p).swap(*this);
            }

            template <class Y, class D>
            void reset(Y * p, D d) { this_type(p, d).swap(*this); }

            template <class Y>
            void reset(shared_ptr<Y> const & r, T * p) {
                this_type(r, p).swap(*this);
            }

            reference operator*() const {
                MADNESS_ASSERT(px != 0);
                return *px;
            }

            T * operator->() const {
                MADNESS_ASSERT(px != 0);
                return px;
            }

            T * get() const { return px; }

            operator bool() const { return px != NULL; }


            // operator! is redundant, but some compilers need it
            bool operator! () const { return px == 0; }

            bool unique() const  { return pn.unique(); }

            long use_count() const { return pn.use_count(); }

            void swap(shared_ptr<T> & other) {
                std::swap(px, other.px);
                pn.swap(other.pn);
            }

            template <class Y> bool _internal_less(shared_ptr<Y> const & rhs) const {
                return pn < rhs.pn;
            }

            void * _internal_get_deleter(std::type_info const & ti) const {
                return pn.get_deleter(ti);
            }

            bool _internal_equiv(shared_ptr const & r) const {
                return px == r.px && pn == r.pn;
            }

        private:

            template <class Y> friend class shared_ptr;
            template <class Y> friend class weak_ptr;

            T * px; // contained pointer
            detail::shared_count pn; // reference counter

        };
        // shared_ptr

        template <class T, class U> inline bool operator==(shared_ptr<T> const & a,
                shared_ptr<U> const & b) {
            return a.get() == b.get();
        }

        template <class T, class U> inline bool operator!=(shared_ptr<T> const & a,
                shared_ptr<U> const & b) {
            return a.get() != b.get();
        }

        template <class T, class U> inline bool operator<(shared_ptr<T> const & a,
                shared_ptr<U> const & b) {
            return a._internal_less(b);
        }

        template <class T> inline void swap(shared_ptr<T> & a, shared_ptr<T> & b) {
            a.swap(b);
        }

        template <class T, class U> shared_ptr<T> static_pointer_cast(shared_ptr<U> const & r) {
            return shared_ptr<T>(r, detail::static_cast_tag());
        }

        template <class T, class U> shared_ptr<T> const_pointer_cast(shared_ptr<U> const & r) {
            return shared_ptr<T>(r, detail::const_cast_tag());
        }

        template <class T, class U> shared_ptr<T> dynamic_pointer_cast(shared_ptr<U> const & r) {
            return shared_ptr<T>(r, detail::dynamic_cast_tag());
        }

        // get_pointer() enables boost::mem_fn to recognize shared_ptr

        template <class T>
        inline T * get_pointer(shared_ptr<T> const & p) { return p.get(); }

        // operator<<

        template <class Y>
        std::ostream & operator<<(std::ostream & os, shared_ptr<Y> const & p) {
            os << p.get();
            return os;
        }

        // get_deleter

        template <class D, class T>
        D * get_deleter(shared_ptr<T> const & p) {
            return static_cast<D *>(p._internal_get_deleter(typeid(D)));
        }

        template <class T> std::size_t hash_value(shared_ptr<T> const & p) {
            return madness::Hash<T*>()(p.get());
        }

        template<class T>
        class weak_ptr {
        private:

            typedef weak_ptr<T> this_type;

        public:

            typedef T element_type;

            weak_ptr(): px(0), pn() { }

            //  generated copy constructor, assignment, destructor are fine


            //
            //  The "obvious" converting constructor implementation:
            //
            //  template<class Y>
            //  weak_ptr(weak_ptr<Y> const & r): px(r.px), pn(r.pn) // never throws
            //  {
            //  }
            //
            //  has a serious problem.
            //
            //  r.px may already have been invalidated. The px(r.px)
            //  conversion may require access to *r.px (virtual inheritance).
            //
            //  It is not possible to avoid spurious access violations since
            //  in multithreaded programs r.px may be invalidated at any point.
            //

            template<class Y>
            weak_ptr( weak_ptr<Y> const & r, typename detail::enable_if_convertible<Y,T>::type = detail::sp_empty() )
                : px(r.lock().get()), pn(r.pn) // never throws
            { }

            template<class Y>
            weak_ptr( shared_ptr<Y> const & r, typename detail::enable_if_convertible<Y,T>::type = detail::sp_empty() )
            : px( r.px ), pn( r.pn ) // never throws
            { }

            shared_ptr<T> lock() const { return shared_ptr<element_type>( *this, detail::sp_nothrow_tag() ); }

            long use_count() const { return pn.use_count(); }

            bool expired() const { return pn.use_count() == 0; }

            void reset() { this_type().swap(*this); }

            void swap(this_type & other) {
                std::swap(px, other.px);
                pn.swap(other.pn);
            }

            void _internal_assign(T * px2, detail::shared_count const & pn2) {
                px = px2;
                pn = pn2;
            }

            template<class Y> bool _internal_less(weak_ptr<Y> const & rhs) const {
                return pn < rhs.pn;
            }

        private:

            template<class Y> friend class weak_ptr;
            template<class Y> friend class shared_ptr;


            T * px;                       // contained pointer
            detail::weak_count pn; // reference counter

        };  // weak_ptr

        template<class T, class U>
        inline bool operator<(weak_ptr<T> const & a, weak_ptr<U> const & b) {
            return a._internal_less(b);
        }

        template<class T>
        void swap(weak_ptr<T> & a, weak_ptr<T> & b) {
            a.swap(b);
        }

        template <class T>
        class enable_shared_from_this {
        protected:

            enable_shared_from_this() { }

            enable_shared_from_this(enable_shared_from_this const &) { }

            enable_shared_from_this &
            operator=(enable_shared_from_this const &) { return *this; }

            ~enable_shared_from_this() { }

        public:

            shared_ptr<T> shared_from_this() {
                shared_ptr<T> p( weak_this_ );
                MADNESS_ASSERT( p.get() == this );
                return p;
            }

            shared_ptr<T const> shared_from_this() const {
                shared_ptr<T const> p( weak_this_ );
                MADNESS_ASSERT( p.get() == this );
                return p;
            }

        public: // actually private, but avoids compiler template friendship issues

            // Note: invoked automatically by shared_ptr; do not call
            template <class X, class Y>
            void _internal_accept_owner( shared_ptr<X> const * ppx, Y * py ) const {
                if( weak_this_.expired() )
                    weak_this_ = shared_ptr<T>( *ppx, py );
            }

        private:

            mutable weak_ptr<T> weak_this_;
        };

        namespace detail {

            class esft2_deleter_wrapper
            {
            private:

                shared_ptr<void> deleter_;

            public:

                esft2_deleter_wrapper() { }

                template< class T >
                void set_deleter( shared_ptr<T> const & deleter ) { deleter_ = deleter; }

                template< class T>
                void operator()( T* ) {
                    MADNESS_ASSERT( deleter_.use_count() <= 1 );
                    deleter_.reset();
                }
            };

        } // namespace detail

        template< class T >
        class enable_shared_from_this2 {
        protected:

            enable_shared_from_this2() { }

            enable_shared_from_this2( enable_shared_from_this2 const & ) { }

            enable_shared_from_this2 & operator=( enable_shared_from_this2 const & ) { return *this; }

            ~enable_shared_from_this2() {
                MADNESS_ASSERT( shared_this_.use_count() <= 1 ); // make sure no dangling shared_ptr objects exist
            }

        private:

            mutable weak_ptr<T> weak_this_;
            mutable shared_ptr<T> shared_this_;

        public:

            shared_ptr<T> shared_from_this() {
                init_weak_once();
                return shared_ptr<T>( weak_this_ );
            }

            shared_ptr<T const> shared_from_this() const {
                init_weak_once();
                return shared_ptr<T>( weak_this_ );
            }

        private:

            void init_weak_once() const {
                if( weak_this_._empty() ) {
                    shared_this_.reset( static_cast< T* >( 0 ), detail::esft2_deleter_wrapper() );
                    weak_this_ = shared_this_;
                }
            }

        public: // actually private, but avoids compiler template friendship issues

            // Note: invoked automatically by shared_ptr; do not call
            template<class X, class Y>
            void _internal_accept_owner( shared_ptr<X> * ppx, Y * py ) const {
                MADNESS_ASSERT( ppx != 0 );

                if( weak_this_.use_count() == 0 ) {
                    weak_this_ = shared_ptr<T>( *ppx, py );
                } else if( shared_this_.use_count() != 0 ) {
                    MADNESS_ASSERT( ppx->unique() ); // no weak_ptrs should exist either, but there's no way to check that

                    detail::esft2_deleter_wrapper * pd = get_deleter<detail::esft2_deleter_wrapper>( shared_this_ );
                    MADNESS_ASSERT( pd != 0 );

                    pd->set_deleter( *ppx );

                    ppx->reset( shared_this_, ppx->get() );
                    shared_this_.reset();
                }
            }
        }; // class enable_shared_from_this2

        } // namespace shptr
    } // namespace tr1
} // namespace madness

#endif
