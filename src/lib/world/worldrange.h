#ifndef MADNESS_WORLD_WORLDRANGE_H__INCLUDED
#define MADNESS_WORLD_WORLDRANGE_H__INCLUDED

/// \file worldrange.h
/// \brief Implement Range class for parallel iteration

namespace madness {

    /// Dummy class a la Intel TBB used to distinguish splitting constructor
    class Split {};

    /// Range vaguely a la Intel TBB encapsulates STL-like start and end iterators with chunksize

    /// Any more sophisticated will wait for the jump to TBB
    template <typename iteratorT>
    class Range {
        long n;
        iteratorT start;
        iteratorT finish;
        int chunksize;
    public:
        typedef iteratorT iterator;

        /// Makes the range [start,finish) ... cost is O(n) due to dumb, linear counting of items

        /// The motivated reader should look at the Intel TBB range,
        /// partitioner, split, concepts, etc..
        ///
        /// Default chunksize is to make 10 tasks per thread to
        /// facilitate dynamic load balancing.
        Range(const iterator& start, const iterator& finish, int chunk=-1)
                : n(0), start(start), finish(finish), chunksize(chunk) {
            for (iterator it=start; it!=finish; ++it) n++;
            if (chunksize == -1) chunksize = n / (10*(ThreadPool::size()+1));
            if (chunksize < 1) chunksize = 1;
        }

        /// Copy constructor ... cost is O(1)
        Range(const Range& r)
                : n(r.n), start(r.start), finish(r.finish), chunksize(r.chunksize) {}

        /// Splits range between new and old (r) objects ... cost is O(n/2)

        /// Presently only bisection down to given chunksize and
        /// executes iterator circa Nlog(N) times so it had better be cheap
        /// compared to the operation being performed.
        Range(Range& left, const Split& split)
                : n(0), start(left.finish), finish(left.finish), chunksize(left.chunksize) {
            //print("SPLITTING: input", left.n, left.chunksize);
            if (left.n > chunksize) {
                start = left.start;
                long nhalf = left.n/2;
                left.n -= nhalf;
                n = nhalf;
                while (nhalf--) {
                    ++start;
                }
                left.finish = start;
            }
            //print("SPLITTING: output: left", left.n, "right", n);
        }

        /// Returns number of items in the range (cost is O(1))
        size_t size() const {
            return n;
        }

        /// Returns true if size=0
        bool empty() const {
            return n==0;
        }

        const iterator& begin() const {
            return start;
        }

        const iterator& end() const {
            return finish;
        }

        unsigned int get_chunksize() const {
            return chunksize;
        }
    };


//         std::vector<iterator> iterators() {
//             size_t n = size();
//             std::vector<iterator> r(n);
//             unsigned int i=0;
//             for (iterator it=begin(); it!=end(); ++it) r[i++] = it;
//             if (i != n) throw "ConcurrentHashMap: count wrong in iterators";
//             return r;
//         }


//     template <typename T>
//     class Range<std::vector<T>::iterator> {



//     }
}

#endif // MADNESS_WORLD_WORLDRANGE_H__INCLUDED
