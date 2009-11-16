#ifndef MADNESS_INDEXIT_H
#define MADNESS_INDEXIT_H

/// \file indexit.h
/// \brief Provides IndexIterator

#include <vector>
#include <world/print.h>

namespace madness {

    /// Facilitates iteration through a multidimension index space.
    /// Since there are multiple ways for navigating index space (column- vs.
    /// row-major, etc.), this class should be abstract, with an abstract ++
    /// operator.  The original IndexIterator assumed the highest dimension
    /// (NDIM-1) would be iterated quickest (this index changes each time ++
    /// is called), however, sometimes a different order is desired.
    ///
    /// For legacy purposes, operator++ is thus NOT abstract, but has the
    /// implementation of the HighDimIndexIterator defined below.  Eventually,
    /// the IndexIterator::operator++ should be deprecated, and such
    /// instances of IndexIterator replaced with HighDimIndexIterator.
    class IndexIterator {
    private:
        /// Bury the default constructor
        IndexIterator() {}

    protected:
        std::vector<long> n; ///< User specified upper limits for each dimension
        std::vector<long> i; ///< Current index
        bool finished;

    public:
        /// Iterates dimension d from 0 to limts[d]-1 inclusive
        template<typename V>
        IndexIterator(const V& limits) :
                n(limits.size()), i(limits.size(), 0), finished(false) {
            for (unsigned int d = 0; d < n.size(); d++)
                n[d] = limits[d];
        }

        /// Iterates dimension d from 0 to limts[d]-1 inclusive
        IndexIterator(int ndim, const long limits[]) :
                n(ndim), i(ndim, 0), finished(false) {
            for (unsigned int d = 0; d < n.size(); d++)
                n[d] = limits[d];
        }

        /// Iterates all dimensions from 0 to top-1 inclusive
        IndexIterator(int ndim, long top) :
                n(ndim,top), i(ndim, 0), finished(false) {
        }

        virtual ~IndexIterator() {}

        IndexIterator&
        reset() {
            for (unsigned int d = 0; d < n.size(); d++)
                i[d] = 0;
            finished = false;
            return *this;
        }

        long
        operator[](int d) const {
            MADNESS_ASSERT(!finished);
            return i[d];
        }

        const std::vector<long>&
        operator*() const {
            MADNESS_ASSERT(!finished);
            return i;
        }

        operator bool() const {
            return !finished;
        }

        /// this function should be abstracted and deprecated
        virtual IndexIterator&
        operator++() {
            for (int d = n.size() - 1; d >= 0; d--) {
                i[d]++;
                if (i[d] < n[d])
                    return *this;
                else
                    i[d] = 0;
            }
            finished = true;
            return *this;
        }

        /// this function should also be deprecated
        static void
        test() {
            Vector<int, 4> n(3);
            for (IndexIterator it(n); it; ++it) {
                print(*it);
            }
        }
    };

    /// The inherited IndexIterator for iterating over the high dimensions
    /// quickly and the low dimensions slowly (all elements of dimension 0
    /// at index i will be visited before any element of index i+1 in dim. 0).
    ///
    /// This is equivalent to the original implementation of IndexIterator.
    class HighDimIndexIterator : public IndexIterator {
    private:
        /// Bury the default constructor
        HighDimIndexIterator() : IndexIterator(0, 0l) {}

    public:
        /// Iterates dimension d from 0 to limts[d]-1 inclusive
        template<typename V>
        HighDimIndexIterator(const V& limits) : IndexIterator(limits) {}

        /// Iterates dimension d from 0 to limts[d]-1 inclusive
        HighDimIndexIterator(int ndim, const long limits[]) :
                IndexIterator(ndim, limits) {}

        /// Iterates all dimensions from 0 to top-1 inclusive
        HighDimIndexIterator(int ndim, long top) : IndexIterator(ndim, top) {}

        virtual HighDimIndexIterator() {}

        /// increment the highest dimension first and check for overflows
        /// up through dimension 0
        virtual IndexIterator&
        operator++() {
            for (int d = n.size() - 1; d >= 0; d--) {
                i[d]++;
                if (i[d] < n[d])
                    return *this;
                else
                    i[d] = 0;
            }
            finished = true;
            return *this;
        }
    };

    /// The inherited IndexIterator for iterating over the low dimensions
    /// quickly and the high dimensions slowly (all elements of dimension
    /// ndim-1 at index i will be visited before any element of index i+1 in
    /// dim. ndim-1).
    class LowDimIndexIterator : public IndexIterator {
    private:
        /// Bury the default constructor
        LowDimIndexIterator() : IndexIterator(0, 0l) {}

    public:
        /// Iterates dimension d from 0 to limts[d]-1 inclusive
        template<typename V>
        LowDimIndexIterator(const V& limits) : IndexIterator(limits) {}

        /// Iterates dimension d from 0 to limts[d]-1 inclusive
        LowDimIndexIterator(int ndim, const long limits[]) :
                IndexIterator(ndim, limits) {}

        /// Iterates all dimensions from 0 to top-1 inclusive
        LowDimIndexIterator(int ndim, long top) : IndexIterator(ndim, top) {}

        virtual ~LowDimIndexIterator() {}

        /// increment the lowest dimension first and check for overflows
        /// up through dimension 0
        virtual IndexIterator&
        operator++() {
            int ndim = n.size();
            for (int d = 0; d < ndim; d++) {
                i[d]++;
                if (i[d] < n[d])
                    return *this;
                else
                    i[d] = 0;
            }
            finished = true;
            return *this;
        }
    };

    /// The inherited IndexIterator for iterating over the dimensions in a
    /// specified order.
    ///
    /// NOTE: if iterating quickly over the high dimensions and slowly over
    /// the low dimensions (in dimensional order), use HighDimIndexIterator.
    ///
    /// NOTE: if iterating quickly over the low dimensions and slowly over
    /// the high dimensions (in dimensional order), use LowDimIndexIterator.
    class NonstandardIndexIterator : public IndexIterator {
    private:
        /// Bury the default constructor
        NonstandardIndexIterator() : IndexIterator(0, 0l) {}

    protected:
        /// the array storing the dimensional order for iteration
        /// dim[0] is the quickest dimension over which to iterate
        /// dim[ndim-1] is the slowest dimension
        std::vector<int> dim;

    public:
        /// Iterates dimension d from 0 to limts[d]-1 inclusive
        ///
        /// order[0] is the dimension to be iterated over quickest
        /// ...
        /// order[d-1] is the dimension to be iterated over slowest
        template<typename V, typename D>
        NonstandardIndexIterator(const V& limits, const D& order) :
                IndexIterator(limits), dim(order.size()) {
            int i, j, ndim = order.size();

            MADNESS_ASSERT(limits.size() == ndim);
            for(i = 0; i < ndim; ++i) {
                MADNESS_ASSERT(order[i] >= 0 && order[i] < ndim);

                // make sure we haven't seen this dimension before
                for(j = 0; j < i; ++j)
                    MADNESS_ASSERT(order[i] != order[j]);

                dim[i] = order[i];
            }
        }

        /// Iterates dimension d from 0 to limts[d]-1 inclusive
        ///
        /// order[0] is the dimension to be iterated over quickest
        /// ...
        /// order[d-1] is the dimension to be iterated over slowest
        NonstandardIndexIterator(int ndim, const long limits[],
                const int order[]) : IndexIterator(ndim, limits) {
            int i, j;

            for(i = 0; i < ndim; ++i) {
                MADNESS_ASSERT(order[i] >= 0 && order[i] < ndim);

                // make sure we haven't seen this dimension before
                for(j = 0; j < i; ++j)
                    MADNESS_ASSERT(order[i] != order[j]);

                dim[i] = order[i];
            }
        }

        /// Iterates all dimensions from 0 to top-1 inclusive
        ///
        /// order[0] is the dimension to be iterated over quickest
        /// ...
        /// order[d-1] is the dimension to be iterated over slowest
        template<typename D>
        NonstandardIndexIterator(int ndim, long top, const D &order) :
                IndexIterator(ndim, top), dim(order.size()) {
            int i, j;

            MADNESS_ASSERT(order.size() == ndim);
            for(i = 0; i < ndim; ++i) {
                MADNESS_ASSERT(order[i] >= 0 && order[i] < ndim);

                // make sure we haven't seen this dimension before
                for(j = 0; j < i; ++j)
                    MADNESS_ASSERT(order[i] != order[j]);

                dim[i] = order[i];
            }
        }

        /// Iterates all dimensions from 0 to top-1 inclusive
        ///
        /// order[0] is the dimension to be iterated over quickest
        /// ...
        /// order[d-1] is the dimension to be iterated over slowest
        NonstandardIndexIterator(int ndim, long top, const int order[]) :
                IndexIterator(ndim, top), dim(ndim) {
            int i, j;

            for(i = 0; i < ndim; ++i) {
                MADNESS_ASSERT(order[i] >= 0 && order[i] < ndim);

                // make sure we haven't seen this dimension before
                for(j = 0; j < i; ++j)
                    MADNESS_ASSERT(order[i] != order[j]);

                dim[i] = order[i];
            }
        }

        virtual ~NonstandardIndexIterator() {}

        /// increment the dimensions in the order detailed in dim
        virtual IndexIterator&
        operator++() {
            int ndim = n.size();
            for (int d = 0; d < ndim; d++) {
                i[dim[d]]++;
                if (i[dim[d]] < n[dim[d]])
                    return *this;
                else
                    i[dim[d]] = 0;
            }
            finished = true;
            return *this;
        }
    };
}

#endif
