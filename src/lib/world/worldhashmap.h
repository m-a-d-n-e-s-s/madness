#ifndef WORLD_HASH_HEADER_INCLUDED
#define WORLD_HASH_HEADER_INCLUDED

/// \file worldhashmap.h
/// \brief Defines and implements a concurrent hashmap


// Why does this exist?  It's a bridge from where we are to where we
// want to be, which is a mutlthreaded environment probably
// based upon the Intel TBB.  Don't have the resources right now to
// bite off the entire TBB but we probably must in the future.
// This is a basic, functional-enough, fast-enough hash map with
// vague compatibility with the TBB API.

#include <world/worldthread.h>
#include <world/worldexc.h>
#include <world/worldhash.h>
#include <new>
#include <iostream>

namespace madness {

    template <class keyT, class valueT, class hashT> class ConcurrentHashMap;

    template <class keyT, class valueT, class hashfunT>
    class ConcurrentHashMap;

    namespace Hash_private {

        // A hashtable is an array of nbin bins.
        // Each bin is a linked list of entries protected by a spinlock.
        // Each entry holds a key+value pair, a read-write mutex, and a link to the next entry.

        template <typename keyT, typename valueT>
        class entry : public madness::MutexReaderWriter {
        public:
            typedef std::pair<const keyT, valueT> datumT;
            datumT datum;

            class entry<keyT,valueT> * volatile next;

            entry(const datumT& datum, entry<keyT,valueT>* next)
                    : datum(datum), next(next) {}
        };

        template <class keyT, class valueT>
        class bin : private madness::Spinlock {
        private:
            typedef entry<keyT,valueT> entryT;
            typedef std::pair<const keyT, valueT> datumT;
            // Could pad here to avoid false sharing of cache line but
            // perhaps better to just use more bins.
        public:

            entryT* volatile p;

            bin() : p(0) {}

            ~bin() {
                clear();
            }

            void clear() {
                lock();             // BEGIN CRITICAL SECTION
                while (p) {
                    entryT* n=p->next;
                    delete p;
                    p=n;
                }
                unlock();           // END CRITICAL SECTION
            }

            entryT* find(const keyT& key, const int lockmode) const {
                bool gotlock;
                entryT* result;
                madness::MutexWaiter waiter;
                do {
                    lock();             // BEGIN CRITICAL SECTION
                    result = match(key);
                    if (result) {
                        gotlock = result->try_lock(lockmode);
                    }
                    else {
                        gotlock = true;
                    }
                    unlock();           // END CRITICAL SECTION
                    if (!gotlock) waiter.wait(); //cpu_relax(); 
                }
                while (!gotlock);

                return result;
            }

            std::pair<entryT*,bool> insert(const datumT& datum, int lockmode) {
                bool gotlock;
                entryT* result;
                bool notfound;
                madness::MutexWaiter waiter;
                do {
                    lock();             // BEGIN CRITICAL SECTION
                    result = match(datum.first);
                    notfound = !result;
                    if (notfound) result = p = new entryT(datum,p);
                    gotlock = result->try_lock(lockmode);
                    unlock();           // END CRITICAL SECTION
                    if (!gotlock) waiter.wait(); //cpu_relax(); 
                }
                while (!gotlock);

                return std::pair<entryT*,bool>(result,notfound);
            }

            bool del(const keyT& key, int lockmode) {
                bool status = false;
                lock();             // BEGIN CRITICAL SECTION
                for (entryT *t=p,*prev=0; t; prev=t,t=t->next) {
                    if (t->datum.first == key) {
                        if (prev) {
                            prev->next = t->next;
                        }
                        else {
                            p = t->next;
                        }
                        t->unlock(lockmode);
                        delete t;
                        status = true;
                        break;
                    }
                }
                unlock();           // END CRITICAL SECTION
                return status;
            }

            std::size_t size() const {
                lock();             // BEGIN CRITICAL SECTION
                std::size_t sum = 0;
                for (entryT *t=p; t; t=t->next) sum++;
                unlock();           // END CRITICAL SECTION
                return sum;
            };

        private:
            entryT* match(const keyT& key) const {
                entryT* t;
                for (t=p; t; t=t->next)
                    if (t->datum.first == key) break;
                return t;
            }

        };

        template <typename T>
        class defhashT {
        public:
            hashT operator()(const T& t) const {
                return hash(t);
            }
        };

        /// iterator for hash
        template <class hashT, class entryT, class datumT> class HashIterator {
        private:
            hashT* h;               // Associated hash table
            int bin;                // Current bin
            entryT* entry;          // Current entry in bin ... zero means at end

            /// If the entry is null (end of current bin) finds next non-empty bin
            void next_non_null_entry() {
                while (!entry) {
                    bin++;
                    if ((unsigned) bin == h->nbins) {
                        entry = 0;
                        return;
                    }
                    entry = h->bins[bin].p;
                }
                return;
            }

        public:
            typedef std::forward_iterator_tag iterator_category;
            typedef datumT value_type;
            typedef std::ptrdiff_t difference_type;
            typedef datumT* pointer;
            typedef datumT& reference;

            /// Makes invalid iterator
            HashIterator() : h(0), bin(-1), entry(0) {}

            /// Makes begin/end iterator
            HashIterator(hashT* h, bool begin)
                    : h(h), bin(-1), entry(0) {
                if (begin) next_non_null_entry();
            }

            /// Makes iterator to specific entry
            HashIterator(hashT* h, int bin, entryT* entry)
                    : h(h), bin(bin), entry(entry) {}

            HashIterator& operator++() {
                if (!entry) return *this;
                entry = entry->next;
                next_non_null_entry();
                return *this;
            }

            HashIterator operator++(int) {
                HashIterator old(*this);
                ++(*this);
                return old;
            }


            bool operator==(const HashIterator& a) const {
                return entry==a.entry;
            }

            bool operator!=(const HashIterator& a) const {
                return entry!=a.entry;
            }

            reference operator*() const {
                MADNESS_ASSERT(entry);
                //if (!entry) throw "Hash iterator: operator*: at end";
                return entry->datum;
            }

            pointer operator->() const {
                MADNESS_ASSERT(entry);
                //if (!entry) throw "Hash iterator: operator->: at end";
                return &entry->datum;
            }
        };

        template <class hashT, typename entryT, typename datumT, int lockmode>
        class HashAccessor : NO_DEFAULTS {
            template <class a,class b,class c> friend class ConcurrentHashMap;
        private:
            entryT* entry;
            bool gotlock;

            /// Used by Hash to set entry (assumed that it has the lock already)
            void set(entryT* entry) {
                release();
                this->entry = entry;
                gotlock = true;
            }

            /// Used by Hash after having already released lock and deleted entry
            void unset() {
                gotlock = false;
                entry = 0;
            }

            void convert_read_lock_to_write_lock() {
                if (entry) entry->convert_read_lock_to_write_lock();
            }


        public:
            HashAccessor() : entry(0), gotlock(false) {}

            HashAccessor(entryT* entry) : entry(entry), gotlock(true) {}

            datumT& operator*() const {
                if (!entry) throw "Hash accessor: operator*: no value";
                return entry->datum;
            }

            datumT* operator->() const {
                if (!entry) throw "Hash accessor: operator->: no value";
                return &entry->datum;
            }

            void release() {
                if (gotlock) {
                    entry->unlock(lockmode);
                    entry=0;
                    gotlock = false;
                }
            }

            ~HashAccessor() {
                release();
            }
        };

    } // End of namespace Hash_private

    template < class keyT, class valueT, class hashfunT = Hash_private::defhashT<keyT> >
    class ConcurrentHashMap {
    public:
        typedef ConcurrentHashMap<keyT,valueT,hashfunT> hashT;
        typedef std::pair<const keyT,valueT> datumT;
        typedef Hash_private::entry<keyT,valueT> entryT;
        typedef Hash_private::bin<keyT,valueT> binT;
        typedef Hash_private::HashIterator<hashT,entryT,datumT> iterator;
        typedef Hash_private::HashIterator<const hashT,const entryT,const datumT> const_iterator;
        typedef Hash_private::HashAccessor<hashT,entryT,datumT,entryT::WRITELOCK> accessor;
        typedef Hash_private::HashAccessor<const hashT,const entryT,const datumT,entryT::READLOCK> const_accessor;

        friend class Hash_private::HashIterator<hashT,entryT,datumT>;
        friend class Hash_private::HashIterator<const hashT,const entryT,const datumT>;

    protected:
        const size_t nbins;         // Number of bins
        binT* bins;                 // Array of bins

    private:
        const iterator _end;
        const const_iterator _const_end;
        hashfunT hashfun;

        //unsigned int hash(const keyT& key) const {return hashfunT::hash(key)%nbins;}

        static int nbins_prime(int n) {
            static const int primes[] = {11, 23, 31, 41, 53, 61, 71, 83, 101, 131, 181, 239, 293, 359, 421, 557, 673, 821, 953, 1021, 1231, 1531, 1747, 2069, 2543, 3011, 4003, 5011, 6073, 7013, 8053, 9029, 9907, 17401, 27479, 37847, 48623, 59377, 70667, 81839, 93199, 104759, 224759, 350411, 479951, 611969, 746791, 882391, 1299743, 2750171, 4256257, 5800159, 7368811, 8960477, 10570871, 12195269, 13834133};
            static const int nprimes = sizeof(primes)/sizeof(int);
            // n is a user provided estimate of the no. of elements to be put
            // in the table.  Want to make the number of bins a prime number
            // larger than this.
            for (int i=0; i<nprimes; i++) if (n<=primes[i]) return primes[i];
            return primes[nprimes-1];
        }

        unsigned int hash_to_bin(const keyT& key) const {
            return hashfun(key)%nbins;
        }

    public:
        ConcurrentHashMap(int n=1021)
                : nbins(hashT::nbins_prime(n))
                , bins(new binT[nbins])
                , _end(this,false)
                , _const_end(this,false) {}

        ConcurrentHashMap(const  hashT& h)
                : nbins(h.nbins)
                , bins(new binT[nbins])
                , _end(this,false)
                , _const_end(this,false) {
            *this = h;
        }

        virtual ~ConcurrentHashMap() {
            delete [] bins;
        }

        hashT& operator=(const  hashT& h) {
            if (this != &h) {
                this->clear();
                for (const_iterator p=h.begin(); p!=h.end(); ++p) {
                    insert(*p);
                }
            }
            return *this;
        }

        std::pair<iterator,bool> insert(const datumT& datum) {
            int bin = hash_to_bin(datum.first);
            std::pair<entryT*,bool> result = bins[bin].insert(datum,entryT::NOLOCK);
            return std::pair<iterator,bool>(iterator(this,bin,result.first),result.second);
        }

        /// Returns true if new pair was inserted; false if key is already in the map
        bool insert(accessor& result, const keyT& key) {
            result.release();
            int bin = hash_to_bin(key);
            std::pair<entryT*,bool> r = bins[bin].insert(datumT(key,valueT()),entryT::WRITELOCK);
            result.set(r.first);
            return r.second;
        }

        /// Returns true if new pair was inserted; false if key is already in the map
        bool insert(const_accessor& result, const keyT& key) {
            result.release();
            int bin = hash_to_bin(key);
            std::pair<entryT*,bool> r = bins[bin].insert(datumT(key,valueT()),entryT::READLOCK);
            result.set(r.first);
            return r.second;
        }

        std::size_t erase(const keyT& key) {
            if (bins[hash_to_bin(key)].del(key,entryT::NOLOCK)) return 1;
            else return 0;
        }

        void erase(const iterator& it) {
            if (it == end()) throw "ConcurrentHashMap: erase(iterator): at end";
            erase(it->first);
        }

        void erase(accessor& item) {
            bins[hash_to_bin(item->first)].del(item->first,entryT::WRITELOCK);
            item.unset();
        }

        void erase(const_accessor& item) {
            item.convert_read_lock_to_write_lock();
            bins[hash_to_bin(item->first)].del(item->first,entryT::WRITELOCK);
            item.unset();
        }

        iterator find(const keyT& key) {
            int bin = hash_to_bin(key);
            entryT* entry = bins[bin].find(key,entryT::NOLOCK);
            if (!entry) return end();
            else return iterator(this,bin,entry);
        }

        const_iterator find(const keyT& key) const {
            int bin = hash_to_bin(key);
            const entryT* entry = bins[bin].find(key,entryT::NOLOCK);
            if (!entry) return end();
            else return const_iterator(this,bin,entry);
        }

        bool find(accessor& result, const keyT& key) {
            result.release();
            int bin = hash_to_bin(key);
            entryT* entry = bins[bin].find(key,entryT::WRITELOCK);
            bool foundit = entry;
            if (foundit) result.set(entry);
            return foundit;
        }

        bool find(const_accessor& result, const keyT& key) const {
            result.release();
            int bin = hash_to_bin(key);
            entryT* entry = bins[bin].find(key,entryT::READLOCK);
            bool foundit = entry;
            if (foundit) result.set(entry);
            return foundit;
        }

        void clear() {
            for (unsigned int i=0; i<nbins; i++) bins[i].clear();
        }

        size_t size() const {
            size_t sum = 0;
            for (size_t i=0; i<nbins; i++) sum += bins[i].size();
            return sum;
        }

        valueT& operator[](const keyT& key) {
            std::pair<iterator,bool> it = insert(datumT(key,valueT()));
            return it.first->second;
        }

        iterator begin() {
            return iterator(this,true);
        }

        const_iterator begin() const {
            return const_iterator(this,true);
        }

        iterator end() {
            return _end;
        }

        const_iterator end() const {
            return _const_end;
        }

    };
}

#endif
