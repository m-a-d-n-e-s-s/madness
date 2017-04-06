#ifndef MADNESS_DISTRIBUTED_MATRIX_H
#define MADNESS_DISTRIBUTED_MATRIX_H

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

// THE STUFF IN THIS FILE IS IN TRANSITION!  THE API AND
// IMPLEMENTATION WILL BOTH SHIFT RAPIDLY AS WE TRANSITION FROM
// REPLICATED TO DISTRIBUTED MATRIX ALGORITHMS, AND SUBSEQUENTLY
// REFINE THE DESIGN AND INTERFACE TO 3RD PARTY PACKAGES.

#include <madness/world/MADworld.h>
#include <utility>
#include <madness/tensor/tensor.h>

namespace madness {

    // If in a fit of misplaced enthusiasm you desire to change
    // int64_t to either long or std::size_t you should be aware that
    // some uses below may contain quantities greater than may be
    // represented in a 32-bit integer and may also be negative.
    // I.e., a simple global replace will fail, though the existing
    // test suite may not detect that.  Also, large skinny matrices
    // could easily need more than 32 bit integers to address.


    // Forward declarations for friends
    class DistributedMatrixDistribution;
    template <typename T> class DistributedMatrix;
    
    static inline DistributedMatrixDistribution column_distributed_matrix_distribution(World& world, int64_t n, int64_t m, int64_t coltile=0);
    static inline DistributedMatrixDistribution row_distributed_matrix_distribution(World& world, int64_t n, int64_t m, int64_t rowtile=0);

    template <typename T>
    DistributedMatrix<T> concatenate_rows(const DistributedMatrix<T>& a, const DistributedMatrix<T>& b);

    template <typename T>
    DistributedMatrix<T> interleave_rows(const DistributedMatrix<T>& a, const DistributedMatrix<T>& b);

    class DistributedMatrixDistribution {
        friend DistributedMatrixDistribution column_distributed_matrix_distribution(World& world, int64_t n, int64_t m, int64_t coltile);
        friend DistributedMatrixDistribution row_distributed_matrix_distribution(World& world, int64_t n, int64_t m, int64_t rowtile);
        template <typename T> friend class DistributedMatrix;

    protected:
        World* pworld;
        int64_t P;                ///< No. of processors
        ProcessID rank;           ///< My processor rank
        int64_t n;                ///< Column dimension of A(n,m)
        int64_t m;                ///< Row dimension of A(n,m)
        int64_t tilen;            ///< Tile size for column
        int64_t tilem;            ///< Tile size for row
        int64_t Pcoldim;          ///< Column dimension of processor grid
        int64_t Prowdim;          ///< Row dimension of processor grid
        int64_t Pcol;             ///< Column of processor grid for this processor
        int64_t Prow;             ///< Row of processor grid for this processor
        int64_t ilo,ihi;          ///< Range of column indices on this processor
        int64_t jlo,jhi;          ///< Range of row indices on this processor
        int64_t idim,jdim;        ///< Dimension of data on this processor


        /// Constructs distribution and size info for a matrix (for use by factory functions only)

        /// This routine is dumb and just copies the given arguments,
        /// hence it can easily make an invalid matrix.  The smarts
        /// are in the factory functions, hence this constructor is
        /// not for general use.
        ///
        /// The matrix is tiled over a grid of processes as specified by the tile sizes.
        /// @param[in] world The world
        /// @param[in] n The matrix column dimension 
        /// @param[in] m The matrix row dimension
        /// @param[in] coltile Tile size for the columns
        /// @param[in] rowtile Tile size for the rows
        DistributedMatrixDistribution(World& world, int64_t n, int64_t m, int64_t coltile, int64_t rowtile)
            : pworld(&world)
            , P(world.size())
            , rank(world.rank())
            , n(n)
            , m(m)
            , tilen(coltile)
            , tilem(rowtile)
            , Pcoldim((n-1)/tilen+1)
            , Prowdim((m-1)/tilem+1)
            , Pcol(rank/Prowdim)
            , Prow(rank - Pcol*Prowdim)
            , ilo(Pcol*tilen)
            , ihi(std::min(ilo+tilen-1,n-1))
            , jlo(Prow*tilem)
            , jhi(std::min(jlo+tilem-1,m-1))
            , idim(std::max(ihi-ilo+1,int64_t(0)))
            , jdim(std::max(jhi-jlo+1,int64_t(0)))
        {
            if (ilo > ihi || jlo > jhi) {
                ilo = jlo = 0;
                ihi = jhi = -1;
            }
        }


    public:

        /// Default constructor makes an invalid distribution
        DistributedMatrixDistribution() 
            : pworld(0)
            , P(0)
            , rank(0)
            , n(0)
            , m(0)
            , tilen(0)
            , tilem(0)
            , Pcoldim(0)
            , Prowdim(0)
            , Pcol(0)
            , Prow(0)
            , ilo(0)
            , ihi(-1)
            , jlo(0)
            , jhi(-1)
            , idim(0)
            , jdim(0)
        {}


        /// Resets state to same as default constructor
        void clear() {
            pworld = (World*)(0);
            P = rank = n = m = tilen = tilem = Pcoldim = Prowdim = Pcol = Prow = ilo = ihi = jlo = jhi = idim = jdim = 0;
        }

        bool operator==(const DistributedMatrixDistribution& d) const {
            return 
                pworld  == d.pworld &&
                P       == d.P      &&
                rank    == d.rank   &&
                n       == d.n      &&
                m       == d.m      &&
                tilen   == d.tilen  &&
                tilem   == d.tilem  &&
                Pcoldim == d.Pcoldim &&
                Prowdim == d.Prowdim &&
                Pcol    == d.Pcol   &&
                Prow    == d.Prow   &&
                ilo     == d.ilo    &&
                ihi     == d.ihi    &&
                jlo     == d.jlo    &&
                jhi     == d.jhi    &&
                idim    == d.idim   &&
                jdim    == d.jdim;
        }


        /// Returns the column dimension of the matrix ... i.e., n for A(n,m)

        /// @return Column dimension of the matrix ... i.e., n for A(n,m)
        int64_t coldim() const {
            return n;
        }

        /// Returns the row dimension of the matrix ... i.e., m for A(n,m)


        /// @return Row dimension of the matrix ... i.e., m for A(n,m)
        int64_t rowdim() const {
            return m;
        }


        /// Returns the column tile size

        /// @return Column tile size
        int64_t coltile() const {
            return tilen;
        }


        /// Returns the row tile size

        /// @return Row tile size
        int64_t rowtile() const {
            return tilem;
        }


        /// Returns the no. of processors in the column dimension

        /// @return No. of processors in the column dimension
        int64_t process_coldim() const {return Pcoldim;}


        /// Returns the no. of processors in the row dimension

        /// @return No. of processors in the rown dimension
        int64_t process_rowdim() const {return Prowdim;}


        /// Returns the total no. of elements stored on this processor

        /// @return Total no. of elements stored on this processor
        int64_t local_size() const {return idim*jdim;}


        /// Returns the no. of column elements stored on this processor

        /// @return No. of column elements stored on this processor (may be zero)
        int64_t local_coldim() const {return idim;}


        /// Returns the no. of row elements stored on this processor

        /// @return No. of row elements stored on this processor
        int64_t local_rowdim() const {return jdim;}

        /// Returns the inclusive range of column indices on this processor

        /// If there is no data on this processor it returns ilow=0 and ihigh=-1
        /// @param[out] ilow First column index on this processor (0 if no data)
        /// @param[out] ihigh Last column index on this processor (-1 if no data)
        void local_colrange(int64_t& ilow, int64_t& ihigh) const {
            ilow = ilo;
            ihigh = ihi;
        }


        /// Returns the inclusive range of row indices on this processor

        /// If there is no data on this processor it returns jlow=0 and jhigh=-1
        /// @param[out] jlow First row index on this processor (0 if no data)
        /// @param[out] jhigh Last row index on this processor (-1 if no data)
        void local_rowrange(int64_t& jlow, int64_t& jhigh) const {
            jlow = jlo;
            jhigh = jhi;
        }


        /// Returns the first column index on this processor (0 if no data present)
        int64_t local_ilow() const {
            return ilo;
        }

        /// Returns the last column index on this processor (-1 if no data present)
        int64_t local_ihigh() const {
            return ihi;
        }


        /// Returns the first row index on this processor (0 if no data present)
        int64_t local_jlow() const  {
            return jlo;
        }


        /// Returns the last row index on this processor (0 if no data present)
        int64_t local_jhigh() const  {
            return jhi;
        }


        /// Returns the inclusive ranges of column and row indicies on processor p

        /// If is no data on processor p it returns ilow=jlow=0 and ihigh=jhigh=-1
        /// @param[in] p The processor p of interest
        /// @param[out] ilow The first column index on the processor
        /// @param[out] ihigh The last column index on the processor (-1 if none)
        /// @param[out] jlow The first row index on the processor
        /// @param[out] jhigh The last row index on the processor (-1 if none)
        void get_range(int p, int64_t& ilow, int64_t& ihigh, int64_t& jlow, int64_t& jhigh) const {
            int pi = p/Prowdim;
            int pj = p - pi*Prowdim;
            if (pi >= process_coldim() || pj >= process_rowdim()) {
                ilow = jlow = 0;
                ihigh = jhigh = -1;
            }
            else {
                ilow = pi*tilen;
                jlow = pj*tilem;
                ihigh= std::min(ilow+tilen-1,n-1);
                jhigh= std::min(jlow+tilem-1,m-1);
            }

            return;
        }


        /// Returns the inclusive range of column indices on processor p

        /// If is no data on processor p it returns ilow=0 and ihigh=-1
        /// @param[in] p The processor p of interest
        /// @param[out] ilow The first column index on the processor
        /// @param[out] ihigh The last column index on the processor (-1 if none)
        void get_colrange(int p, int64_t& ilow, int64_t& ihigh) const {
            int64_t jlow, jhigh;
            get_range(p, ilow, ihigh, jlow, jhigh);

            return;
        }


        /// Returns the inclusive range of row indices on processor p

        /// If is no data on processor p it returns jlow=0 and jhigh=-1
        /// @param[in] p The processor p of interest
        /// @param[out] jlow The first row index on the processor
        /// @param[out] jhigh The last row index on the processor (-1 if none)
        void get_rowrange(int p, int64_t& jlow, int64_t& jhigh) const {
            int64_t ilow, ihigh;
            get_range(p, ilow, ihigh, jlow, jhigh);

            return;
        }


        /// Returns the associated world

        /// @return The world
        World& get_world() const {return *pworld;}


        /// Returns true if the matrix is column distributed (i.e., row dimension not distributed)

        /// @return True if the matrix is column distributed (i.e., row dimension not distributed)
        bool is_column_distributed() const {return process_rowdim()==1;}


        /// Returns true if the matrix is row distributed (i.e., column dimension not distributed)

        /// @return True if the matrix is row distributed (i.e., column dimension not distributed)
        bool is_row_distributed() const {return process_coldim()==1;}


        /// Returns the distribution (aka *this)
        const DistributedMatrixDistribution& distribution() const {return *this;}


        /// Returns the number of the process that owns element (i,j)
        ProcessID owner(int64_t i, int64_t j) const {
            int pcol = i/coltile();
            int prow = j/rowtile();
            
            return pcol*process_rowdim() + prow;
        }

        virtual ~DistributedMatrixDistribution() {}
    };


    /// Manages data associated with a row/column/block distributed array

    /// The class itself provides limited functionality for accessing the
    /// data and is primarily intended to provide base functionality for
    /// use by matrix algorithms and other matrix classes.
    ///
    /// The constructor is deliberately simple.  Factory functions
    /// are expected to be the main construction tool.
    ///
    /// Assignment and copy are shallow just like for tensor and for the same reasons.  
    ///
    /// To get a deep copy use the copy function (again just like for tensors).
    template <typename T>
    class DistributedMatrix : public DistributedMatrixDistribution {
        friend DistributedMatrix<T> interleave_rows<T>(const DistributedMatrix<T>& a, const DistributedMatrix<T>& b);
        friend DistributedMatrix<T> concatenate_rows<T>(const DistributedMatrix<T>& a, const DistributedMatrix<T>& b);

        Tensor<T> t;            ///< The data

        static T idij(const int64_t i, const int64_t j) {return (i==j) ?  T(1) : T(0);}

    protected:

        /// Constructs a distributed matrix dimension (n,m) with specified tile sizes and initialized to zero

        /// [deprecated ... use factory functions instead]
        /// 
        /// The matrix is tiled over a grid of processes as specified by the tile sizes.
        /// @param[in] world The world
        /// @param[in] n The matrix column dimension 
        /// @param[in] m The matrix row dimension
        /// @param[in] coltile Tile size for the columns
        /// @param[in] rowtile Tile size for the rows
        DistributedMatrix(World& world, int64_t n, int64_t m, int64_t coltile, int64_t rowtile)
            : DistributedMatrixDistribution(world, n, m, coltile, rowtile)
        {
            if (idim>0 && jdim>0) t = Tensor<T>(idim,jdim);
        }

    public:

        /// Default constructor makes an empty matrix that cannot be used except as a target for assignemnt
        DistributedMatrix()
            : DistributedMatrixDistribution()
            , t()
        {}


        /// Constructs a distributed matrix with given distribution info
        DistributedMatrix(const DistributedMatrixDistribution& d)
            : DistributedMatrixDistribution(d)
        {
            if (idim>0 && jdim>0) t = Tensor<T>(idim,jdim);
        }


        /// Copy constructor copies dimensions, distribution, and shallow copy of content (unless deepcopy=true)
        DistributedMatrix(const DistributedMatrix<T>& A, bool deepcopy=false)
            : DistributedMatrixDistribution(A)
            , t(deepcopy ? copy(A.t) : A.t) 
        {}


        /// Assigment copies dimensions, distribution, and shallow copy of content
        DistributedMatrix<T>& operator=(const DistributedMatrix<T>& A) {
            if (this != &A) {
                DistributedMatrixDistribution::operator=(A);
                t = A.t;
            }
            return *this;
        }

        virtual ~DistributedMatrix() {}


        /// Frees memory and resets state to same as default constructor
        void clear() {
            DistributedMatrixDistribution::clear();
            t.clear();
        }


        /// Fills the matrix with the provided function of the indices

        /// @param[in] f The matrix is filled using \c a[i,j]=f(i,j)
        template <typename funcT>
        void fill(const funcT& f) {
            for (int64_t i=ilo; i<=ihi; i++) {
                for (int64_t j=jlo; j<=jhi; j++) {
                    t(i-ilo,j-jlo) = f(i,j);
                }
            }
        }


        /// Fills the matrix with a scalar

        /// @param[in] value The matrix is filled using \c a[i,j]=value
        void fill(T value) {
            t.fill(value);
        }


        void fill_identity() {
            fill(DistributedMatrix<T>::idij);
        }


        /// Returns reference to the local data 

        /// The local data is a tensor dimension \c (local_coldim,local_rowdim) and if either of the dimensions
        /// are zero there is no data.  A natural way to loop thru the data that gives you the actual row and column indices is
        /// \code
        /// const Tensor<T>& t = A.data();
        /// for (int64_t i=A.get_ilow(); i<=A.get_ihigh(); i++) {
        ///    for (int64_t j=A.get_jlow(); j<=A.get_jhigh(); j++) {
        ///        the ij'th element is t(i-ilo, j-jlo)
        ///    }
        /// }
        /// \endcode
        /// @return Reference to non-constant tensor holding the local data
        Tensor<T>& data() {return t;}

        /// Returns const reference to data

        /// The local data is a tensor dimension \c (local_coldim,local_rowdim) and if either of the dimensions
        /// are zero there is no data.  A natural way to loop thru the data that gives you the actual row and column indices is
        /// \code
        /// const Tensor<T>& t = A.data();
        /// for (int64_t i=A.get_ilow(); i<=A.get_ihigh(); i++) {
        ///    for (int64_t j=A.get_jlow(); j<=A.get_jhigh(); j++) {
        ///        the ij'th element is t(i-ilo, j-jlo)
        ///    }
        /// }
        /// \endcode
        /// @return Reference to constant tensor holding the local data
        const Tensor<T>& data() const {return t;}


        /// Copy from the replicated \c (m,n) matrix into the distributed matrix

        /// @param[in] s The input tensor
        void copy_from_replicated(const Tensor<T>& s) {
            if (local_size() > 0) t(___) = s(Slice(ilo,ihi),Slice(jlo,jhi));
        }

        /// Copy from the distributed \c (m,n) matrix into the replicated matrix (collective call)

        /// The entire output tensor is zeroed, the local data copied
        /// into it, and then a global sum performed to replicate the
        /// data.
        /// @param[out] s The output tensor (assumed already allocated
        /// to be at least large enough in both dimensions)
        void copy_to_replicated(Tensor<T>& s) const {
            MADNESS_ASSERT(s.iscontiguous());
            s = 0.0;
            if (local_size() > 0) s(Slice(ilo,ihi),Slice(jlo,jhi)) = t(___);
            get_world().gop.sum(s.ptr(), s.size());
        }

        /// Copy from replicated patch (inclusive index range) into the distributed matrix

        /// @param[in] ilow First \c i index in patch
        /// @param[in] ihigh Last \c i index in patch
        /// @param[in] jlow First \c j index in patch
        /// @param[in] jhigh Last \c j index in patch
        void copy_from_replicated_patch(int64_t ilow, int64_t ihigh, int64_t jlow, int64_t jhigh, const Tensor<T>& s) {
            int64_t i0 = std::max(ilo,ilow);
            int64_t j0 = std::max(jlo,jlow);
            int64_t i1 = std::min(ihi,ihigh);
            int64_t j1 = std::min(jhi,jhigh);
            if (i0<=i1 && j0<=j1) {
                t(Slice(i0-ilo,i1-ilo),Slice(j0-jlo,j1-jlo)) = s(Slice(i0-ilow,i1-ilow),Slice(j0-jlow,j1-jlow));
            }
        }

        /// Copy from distributed matrix into replicated patch (inclusive index range; collective call)

        /// The entire output tensor is zeroed, relevant local data
        /// copied into it, and then a global sum performed to
        /// replicate the data.
        /// @param[in] ilow First \c i index in patch
        /// @param[in] ihigh Last \c i index in patch
        /// @param[in] jlow First \c j index in patch
        /// @param[in] jhigh Last \c j index in patch
        void copy_to_replicated_patch(int64_t ilow, int64_t ihigh, int64_t jlow, int64_t jhigh, Tensor<T>& s) const {
            MADNESS_ASSERT(s.iscontiguous());
            s = 0;
            int64_t i0 = std::max(ilo,ilow);
            int64_t j0 = std::max(jlo,jlow);
            int64_t i1 = std::min(ihi,ihigh);
            int64_t j1 = std::min(jhi,jhigh);
            if (i0<=i1 && j0<=j1) {
                t(Slice(i0-ilo,i1-ilo),Slice(j0-jlo,j1-jlo)) = s(Slice(i0-ilow,i1-ilow),Slice(j0-jlow,j1-jlow));
            }
            get_world().gop.sum(s.ptr(), s.size());
        }

        void extract_columns(int64_t jlow, int64_t jhigh, DistributedMatrix<T>& U) const {
            int newrowdim = jhigh - jlow + 1;
            MADNESS_ASSERT(jlow >= 0);
            MADNESS_ASSERT(jhigh < rowdim());
            MADNESS_ASSERT(newrowdim == U.rowdim());
            MADNESS_ASSERT(coldim() == U.coldim());
            MADNESS_ASSERT(is_column_distributed());
            MADNESS_ASSERT(U.is_column_distributed());
            MADNESS_ASSERT(coltile() == U.coltile());

            int64_t i0 = ilo;
            int64_t j0 = std::max(jlo,jlow);
            int64_t i1 = ihi;
            int64_t j1 = std::min(jhi,jhigh);
            if (i0<=i1 && j0<=j1) {
                U.data()(___) = t(Slice(i0-ilo,i1-ilo),Slice(j0-jlo,j1-jlo));
            }
        }

        template <typename R> 
        bool has_same_dimension_and_distribution(const DistributedMatrix<R>& A) {
            return DistributedMatrixDistribution::operator==(A);
        }

        /// Inplace addition --- dimensions and distribution must be identical

        /// @param[in] A The matrix to add to the current matrix
        /// @return A reference to the current matrix
        DistributedMatrix<T>& operator+=(const DistributedMatrix<T>& A) {
            MADNESS_ASSERT(has_same_dimension_and_distribution(A));
            t += A.t;
            return *this;
        }


        /// Out of place addition  --- dimensions and distribution must be identical

        /// @param[in] A The matrix to add to the current matrix
        /// @return A new matrix with the same dimensions and distribution as the inputs
        DistributedMatrix<T> operator+(const DistributedMatrix<T>& A) const {
            MADNESS_ASSERT(has_same_dimension_and_distribution(A));
            return copy(*this)+=A;
        }


        /// Inplace scale by a constant

        /// @param[in] s The scaling factor
        /// @return A reference to the current matrix
        DistributedMatrix<T>& operator*=(const T s) {
            t.scale(s);
            return *this;
        }

        /// Sets element (i,j) to v if (i,j) is local, otherwise throws MadnessException
        void set(int64_t i, int64_t j, const T x) {
            MADNESS_ASSERT(i>=ilo && i<=ihi && j>=jlo && j<=jhi);
            t(i-ilo,j-jlo) = x;
        }

        /// Gets element (i,j) if (i,j) is local, otherwise throws MadnessException
        T get(int64_t i, int64_t j) const {
            MADNESS_ASSERT(i>=ilo && i<=ihi && j>=jlo && j<=jhi);
            return t(i-ilo,j-jlo);
        }
    };


    /// Deep copy of content

    /// @param[in] A The matrix to be copied
    /// @return A new matrix with identical dimensions, distribution and content (deep copy)
    template <typename T>
    DistributedMatrix<T> copy(const DistributedMatrix<T>& A) {
        return DistributedMatrix<T>(A,true);
    }
    

    /// Generates distribution for an (n,m) matrix distributed by columns (row dimension is not distributed)

    /// Quietly forces an even column tile size for ease of use in the systolic matrix algorithms
    /// @param[in] world The world
    /// @param[in] n The column (first) dimension
    /// @param[in] m The row (second) dimension
    /// @param[in] coltile Tile size for columns forced to be even (default is to use all processes)
    /// @return An object encoding the dimension and distribution information
    static inline DistributedMatrixDistribution
    column_distributed_matrix_distribution(World& world, int64_t n, int64_t m, int64_t coltile) { // default coltile=0 above
        if (world.size()*coltile < n) coltile = (n-1)/world.size() + 1;
        coltile = std::min(coltile,n);
        if ((coltile&0x1)) ++coltile; // ??? Was before the previous statement

        return DistributedMatrixDistribution(world, n, m, coltile, m);
    }


    /// Generates an (n,m) matrix distributed by columns (row dimension is not distributed)

    /// Quietly forces an even column tile size for ease of use in the systolic matrix algorithms
    /// @param[in] world The world
    /// @param[in] n The column (first) dimension
    /// @param[in] m The row (second) dimension
    /// @param[in] coltile Tile size for columns forced to be even (default is to use all processes)
    /// @return A new zero matrix with the requested dimensions and distribution
    template <typename T>
    DistributedMatrix<T> column_distributed_matrix(World& world, int64_t n, int64_t m, int64_t coltile=0) {
        return DistributedMatrix<T>(column_distributed_matrix_distribution(world, n, m, coltile));
    }


    /// Generates an (n,m) matrix distribution distributed by rows (column dimension is not distributed)

    /// @param[in] world The world
    /// @param[in] n The column (first) dimension
    /// @param[in] m The row (second) dimension
    /// @param[in] rowtile Tile size for row (default is to use all processes)
    /// @return An object encoding the dimension and distribution information
    static inline DistributedMatrixDistribution
    row_distributed_matrix_distribution(World& world, int64_t n, int64_t m, int64_t rowtile) { // default rowtile=0 above
        if (world.size()*rowtile < m) rowtile = (m-1)/world.size() + 1;
        rowtile = std::min(rowtile,m);

        return DistributedMatrixDistribution(world, n, m, n, rowtile);
    }

    /// Generates an (n,m) matrix distributed by rows (column dimension is not distributed)

    /// @param[in] world The world
    /// @param[in] n The column (first) dimension
    /// @param[in] m The row (second) dimension
    /// @param[in] rowtile Tile size for row (default is to use all processes)
    /// @return A new zero matrix with the requested dimensions and distribution
    template <typename T>
    DistributedMatrix<T> row_distributed_matrix(World& world, int64_t n, int64_t m, int64_t rowtile=0) {
        return DistributedMatrix<T>(row_distributed_matrix_distribution(world, n, m, rowtile));
    }


    /// Generates a distributed matrix with rows of \c a and \c b interleaved

    /// I.e., the even rows of the result will be rows of \c a , and the
    /// odd rows those of \c b .
    ///
    /// The matrices a and b must have the same dimensions and be
    /// identically distributed.  The result will have a doubled column
    /// dimension and column tile size.  The row dimension is unchanged.
    /// @param[in] a The matrix providing even rows of the result
    /// @param[in] b The matrix providing odd rows of the result
    /// @return The result matrix
    template <typename T>
    DistributedMatrix<T> interleave_rows(const DistributedMatrix<T>& a, const DistributedMatrix<T>& b) {
        MADNESS_ASSERT(a.rowdim()==b.rowdim() && a.coldim()==b.coldim() && a.coltile()==b.coltile() && a.rowtile()==b.rowtile());

        DistributedMatrix<T> c(a.get_world(), a.coldim()*2, a.rowdim(), a.coltile()*2, a.rowtile());
        c.data()(Slice(0,-1,2),_) = a.data()(___);
        c.data()(Slice(1,-1,2),_) = b.data()(___);
    }


    /// Generates a column-distributed matrix with rows of \c a and \c b contatenated

    /// I.e., c[i,j] = a[i,j] if j<na or b[i,j-na] if j>=na
    ///
    /// The matrices a and b must have the same column size (i.e., the
    /// same number of rows) and be column distributed with the same
    /// column tilesze.  The result is also column distributed with
    /// the same column tilesize as the input matrices.
    /// @param[in] a The first matrix
    /// @param[in] b The second matrix
    /// @return The result matrix
    template <typename T>
    DistributedMatrix<T> concatenate_rows(const DistributedMatrix<T>& a, const DistributedMatrix<T>& b) {
        MADNESS_ASSERT(a.coldim()==b.coldim() && a.coltile()==b.coltile() && a.is_column_distributed() && b.is_column_distributed());

        int64_t ma = a.rowdim();
        int64_t mb = b.rowdim();

        DistributedMatrix<T> c(a.get_world(), a.coldim(), ma+mb, a.coltile(), ma+mb);
        
        int64_t ilow, ihigh;
        a.local_colrange(ilow, ihigh);
        if (ilow <= ihigh) {
            c.data()(_,Slice(0,ma-1)) = a.data()(___);
            c.data()(_,Slice(ma,-1))  = b.data()(___);
        }

        return c;
    }

    /// Generates a column-distributed matrix with rows of \c a, \c b, \c c, and \c d contatenated in order

    /// I.e., c[i,j] = a[i,j] if j in [0,na), b[i,j-na] if j in [na,na+nb), c[i,j-na-nb] if j in [na+nb,na+nb+nc), etc.
    /// The matrices must have the same column size (i.e., the same
    /// number of rows) and be column distributed with the same column
    /// tilesze.  The result is also column distributed with the same
    /// column tilesize as the input matrices.
    /// @param[in] a The first matrix
    /// @param[in] b The second matrix
    /// @param[in] c The third matrix
    /// @param[in] d The fourth matrix
    /// @return The result matrix
    template <typename T>
    DistributedMatrix<T> concatenate_rows( const DistributedMatrix<T>& a, const DistributedMatrix<T>& b, const DistributedMatrix<T>& c, const DistributedMatrix<T>& d) {
        MADNESS_ASSERT(a.coldim()==b.coldim() && b.coldim()==c.coldim() && c.coldim()==d.coldim());
        MADNESS_ASSERT(a.coltile()==b.coltile() && b.coltile()==c.coltile() && c.coltile()==d.coltile());
        MADNESS_ASSERT(a.is_column_distributed() && b.is_column_distributed() && c.is_column_distributed() && d.is_column_distributed());
        
        int64_t ma = a.rowdim();
        int64_t mb = b.rowdim();
        int64_t mc = c.rowdim();
        int64_t md = d.rowdim();
        
        DistributedMatrix<T> result(a.get_world(), a.coldim(), ma+mb+mc+md, a.coltile(), ma+mb+mc+md);

        if(a.local_size() > 0) result.data()( _ , Slice(0,ma-1) ) = a.data()(___);
        if(b.local_size() > 0) result.data()( _ , Slice(ma, ma+mb-1) ) = b.data()(___);
        if(c.local_size() > 0) result.data()( _ , Slice(ma+mb, ma+mb+mc-1) ) = c.data()(___);
        if(d.local_size() > 0) result.data()( _ , Slice(ma+mb+mc, -1) ) = d.data()(___);

        return result;
    }


    /// Generates a row-distributed matrix with rows of \c a and \c b contatenated

    /// I.e., c[i,j] = a[i,j] if i<ma or b[i-ma,j] if i>=ma
    ///
    /// The matrices a and b must have the same row size (i.e., the
    /// same number of columns) and be row distributed with the same
    /// row tilesze.  The result is also row distributed with
    /// the same row tilesize as the input matrices.
    /// @param[in] a The first matrix
    /// @param[in] b The second matrix
    /// @return The result matrix
    template <typename T>
    DistributedMatrix<T> concatenate_columns(const DistributedMatrix<T>& a, const DistributedMatrix<T>& b) {
        MADNESS_ASSERT(a.rowdim()==b.rowdim() && a.rowtile()==b.rowtile() && a.is_row_distributed() && b.is_row_distributed());

        int64_t ma = a.coldim();
        int64_t mt = ma + b.coldim();

        DistributedMatrix<T> c(a.get_world(), mt, a.rowdim(), b.rowtile(), mt);

        if(a.local_size() > 0) c.data()( Slice(0,ma-1), _ ) = a.data()(___);
        if(a.local_size() > 0) c.data()( Slice(ma,-1), _ ) = b.data()(___);

        return c;
    }
}

#endif
