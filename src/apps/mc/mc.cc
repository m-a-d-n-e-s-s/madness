#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>

#include <madness/misc/ran.h>
#include <madness/world/print.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>

#include "mc.h"

using namespace madness;

namespace std {
    template <> struct hash<std::pair<size_t,size_t>> {
        size_t operator()(const std::pair<size_t,size_t>& ij) const {
            return (ij.first<<32) | ij.second;
        }
    };
}

namespace madness {
    
    // Simple sparse vector class 
    template <typename T=double>
    class sparse_vector {
        size_t M;
        std::vector<size_t> indices;
        std::vector<T> data;
    public:

        // Make an empty sparse vector of length M
        sparse_vector(size_t M) : M(M) {}
        
        // Make a sparse vector from a dense vector with screening
        template <typename vecT>
        sparse_vector(const vecT& x, T tol=0) : M(x.size()) {
            for (size_t i=0; i<M; i++) {
                if (std::abs(x[i]) > tol) {
                    indices.push_back(i);
                    data.push_back(x[i]);
                }
            }
        }

        // Make a dense tensor 
        Tensor<T> dense() const {
            Tensor<T> v(M);
            for (auto [i,x] : *this) v[i] = x;
            return v;
        }
        
        // Return the indices of the non-zero elements
        const std::vector<size_t>& get_indices() const {return indices;}
        std::vector<size_t>& get_indices() {return indices;}

        // Return the values of the non-zero elements
        const std::vector<T>& get_data() const {return data;}
        std::vector<T>& get_data() {return data;}

        // Returns the length (if dense) of the vector
        size_t size() const {return M;}

        // Returns the number of nonzero elements of the vector
        size_t nnz() const {return indices.size();}
        
        // Push (i,value) onto the end of the indices and data
        void push_back(size_t i, T value) {
            MADNESS_ASSERT(i>=0 && i<M);
            indices.push_back(i);
            data.push_back(value);
        }
        
        // Sort the index vector into increasing order
        void sort_indices() {
            using pairT = std::pair<size_t,T>;
            std::vector<pairT> x; x.reserve(nnz());
            for (size_t i=0; i<data.size(); i++) x.push_back({indices[i],data[i]});
            std::sort(x.begin(), x.end(), [](const pairT& a, const pairT& b){return a.first<b.first;});
            indices = std::vector<size_t>(indices.size());
            data = std::vector<T>(data.size());
            for (size_t i=0; i<data.size(); i++) {indices[i]=x[i].first,data[i]=x[i].second;}
        }

        // Compress data keeping entries with absolute value greater than tol
        void compress(T tol=0.0) {
            std::vector<size_t> ind; ind.reserve(nnz());
            std::vector<T> dat; dat.reserve(nnz());
            for (auto [i,x] : *this) {
                if (std::abs(x)>tol) {
                    ind.push_back(i);
                    dat.push_back(x);
                }
            }
            indices = ind; // Assign/copy should shrink capacity to match size but perhaps move is done instead?
            data = dat;
        }

        // return (true,value) if index i exists, (false,0.0) otherwise --- this is presently implemented with dumb linear search and so is slooooow
        std::pair<bool,T> find(size_t i) const {
            for (auto [j,x] : *this) if (j == i) return std::make_pair(true,x);
            return {false,0.0};
        }

        // Non-const iterator deferences into pair(index&,value&)
        struct iterator {
            std::vector<size_t>& indices;
            std::vector<T>& data;
            size_t i;
            iterator(const iterator& other) : indices(other.indices), data(other.data), i(other.i) {}
            iterator(std::vector<size_t>& indices, std::vector<T>& data, size_t start=0) : indices(indices), data(data), i(start) {}
            iterator& operator++() {++i; return *this;}
            bool operator!=(const iterator& other) const {return i!=other.i;}
            std::pair<size_t&,T&> operator*() {return std::pair<size_t&,T&>(indices[i],data[i]);}
        };

        // Const iterator deferences into pair(index,value)
        struct const_iterator {
            const std::vector<size_t>& indices;
            const std::vector<T>& data;
            size_t i;
            const_iterator(const iterator& other) : indices(other.indices), data(other.data), i(other.i) {}
            const_iterator(const std::vector<size_t>& indices, const std::vector<T>& data, size_t start=0) : indices(indices), data(data), i(start) {}
            const_iterator& operator++() {++i; return *this;}
            bool operator!=(const const_iterator& other) const {return i!=other.i;}
            std::pair<size_t,T> operator*() {return std::pair<size_t,T>(indices[i],data[i]);}
        };

        const_iterator begin() const {return const_iterator(indices,data,0);}
        const_iterator end() const {return const_iterator(indices,data,indices.size());}
        iterator begin() {return iterator(indices,data,0);}
        iterator end() {return iterator(indices,data,indices.size());}
    };
    
    // Print sparse vector
    template <typename T>
    std::ostream& operator<<(std::ostream& s, const sparse_vector<T>& v) {
        using madness::operators::operator<<;        
        s<<"(";
        s<<v.get_indices();
        s<<", ";
        s<<v.get_data();
        s<<")";
        return s;
    }
    
    // Simple class storing sparse matrix as map {i,j}->value --- inherits most of interface from std::unordered_map
    template <typename T=double>
    class sparse_matrix_coo : public  std::unordered_map<std::pair<size_t,size_t>,T> {
        size_t N,M;
    public:
        sparse_matrix_coo(size_t N, size_t M) : N(N), M(M) {}

        // Return the number of rows or columns in the matrix
        size_t nrow() const {return N;}
        size_t ncol() const {return M;}    
    };
    
    // Simple class storing sparse matrix as compressed sparse row
    template <typename T=double>
    class sparse_matrix_csr {
        size_t N,M;
        std::vector<sparse_vector<T>> data;
    public:
        typedef sparse_vector<T>* iterator;
        typedef const sparse_vector<T>* const_iterator;
        
        sparse_matrix_csr(size_t N, size_t M) : N(N), M(M), data(N,M) {}
        
        sparse_matrix_csr(const sparse_matrix_coo<T>& A) : N(A.nrow()), M(A.ncol()), data(N,M) { 
            for (const auto& [ij,value] : A) {
                const auto [i,j] = ij;
                if (i<0 || i>=N) print("bad", i, j, value);
                get_row(i).push_back(j,value);
            }
            for (auto& row : data) row.sort_indices();
        }

        // Returns the sparse vector representing row i
        const sparse_vector<T>& get_row(size_t i) const {return data.at(i);}
        sparse_vector<T>& get_row(size_t i) {return data.at(i);}

        const_iterator begin() const {return &data[0];}
        const_iterator end() const {return (&data[0])+N;}
        iterator begin() {return &data[0];}
        iterator end() {return (&data[0])+N;}
        
        // Return the number of rows or columns in the matrix
        size_t nrow() const {return N;}
        size_t ncol() const {return M;}

        // Compress out zero/small elements
        void compress() {
            for (auto& row : *this) row.compress();
        }

        // Return a new matrix with elements in the specified slices.
        // Note a matrix is always returned and zero size slice will
        // throw.  Also, MADNESS slices are inclusive ranges whereas
        // Python slice have an exclusive endpoint.  Look at slice.h
        // and tensor.h for more info.
        sparse_matrix_csr<T> operator()(const Slice& rows, const Slice& cols) const {
            long rowstart = rows.start;  if (rowstart < 0) rowstart += N;
            long rowend = rows.end;  if (rowend < 0) rowend += N;
            long rowstep = rows.step;
            long colstart = cols.start;  if (colstart < 0) colstart += M;
            long colend = cols.end;  if (colend < 0) colend += M;
            long colstep = cols.step;

            if (rowstart == rowend && rowstep==0) throw "empty row slice not permitted";
            if (colstart == colend && colstep==0) throw "empty col slice not permitted";
            if (colstep!=1 || rowstep!=1) throw "non-unit step not supported yet";
            MADNESS_ASSERT(rowstart >= rowend);
            MADNESS_ASSERT(colstart >= colend);

            sparse_matrix_csr<T> A(rowend-rowstart+1, colend-colstart+1);
            for (long i=rowstart; i<=rowend; i++) {
                const auto& oldrow = get_row(i);
                auto& newrow = A.get_row(i-rowstart);
                for (auto [j,x] : oldrow) {
                    if (long(j)>=colstart && long(j)<=colend)
                        newrow.push_back(j-colstart,x);
                }
            }
            return A;
        }
    };

    // Read A test matrix from the configuration interaction (CI) code
    sparse_matrix_csr<double> load_ci_matrix(const std::string& filename) {
        std::fstream file(filename,std::ios::in);
        size_t N;
        file >> N;
        //sparse_matrix_coo<double> A(N,N);
        sparse_matrix_csr<double> A(N,N);
        size_t i, j;
        double value;
        size_t count=0;
        while (file >> i >> j >> value) {
            i--; j--; // fortran indexes from 1, c++ from 1
            MADNESS_ASSERT(i>=0 && i<N && j>=0 && j<N);
            //A[{i,j}] = value; A[{j,i}] = value; // read full square
            A.get_row(i).push_back(j,value);
            if (i != j) A.get_row(j).push_back(i,value);
            count++;
        }
        print("load_ci_matrix:",N,count);
        
        //sparse_matrix_csr<double> B(A);
        //print(B.nrow());
        for (auto& row : A) row.sort_indices();
        return A;
    }

    // Converts a sparse_matrix_csr into a MADNESS tensor
    template <typename T>
    Tensor<T> sparse_matrix_to_tensor(const sparse_matrix_csr<T>& A) {
        Tensor<T> B(A.nrow(),A.ncol());
        size_t i=0;
        for (const auto& row: A) {
            for (const auto& [j,value] : row) B(i,j) = value;
            i++;
        }
        return B;
    }

    // Converts a MADNESS tensor into a sparse_matrix_csr
    template <typename T>
    sparse_matrix_csr<T> tensor_to_sparse_matrix(const Tensor<T>& A, T tol=0.0) {
        size_t N=A.dim(0), M=A.dim(1);
        sparse_matrix_csr<T> B(N,M);

        for (size_t i=0; i<N; i++) {
            sparse_vector<T>& row = B.get_row(i);
            for (size_t j=0; j<M; j++) {
                T value = A(i,j);
                if (std::abs(value) > tol) row.push_back(j,value);
            }
            row.get_indices().shrink_to_fit();
            row.get_data().shrink_to_fit();
        }
        return B;
    }
    
    // Multiplies a sparse matrix (CSR) onto a vector (that can be either std::vector or madness::Tensor)
    template <typename T, typename vecT>
    vecT sparse_matrix_times_vector(const sparse_matrix_csr<T>& A, const vecT& v) {
        MADNESS_ASSERT(A.ncol() == size_t(v.size()));
        vecT Av(v.size());
        size_t i = 0;
        for (const auto& row: A) {
            T sum = 0.0;
            for (const auto& [j,aij] : row) sum += aij * v[j];
            Av[i++] = sum;
        }
        return Av;
    }
    
    // Power iteration with shift intended to crudely estimate extreme eigenvalues
    template <typename T>
    T simple_power_iteration(const sparse_matrix_csr<T>& A, T shift, double tol=0.001, bool ranguess=true) {
        Tensor<double> x(A.nrow());
        if (ranguess) x.fillrandom();
        else x[0] = 1.0;
        
        T eprev = 0.0;
        for (int iter=0; iter<1000; iter++) {
            x.scale(1.0/x.normf());
            Tensor<T> Ax = sparse_matrix_times_vector(A, x);
            T eval = x.trace(Ax);
            //print(iter,eval);
            if ((iter%10)==1) {
                if (std::abs((eval-eprev)/eval) < tol) return eval;
                eprev = eval;
            }
            x = Ax - shift*x;
        }
        return eprev;
    }
    
    // More efficient power iteration for lowest eigenvalue using Chebyshev nodes
    template <typename T>
    T power_iteration(const sparse_matrix_csr<T>& A, T e0, T e1, T ehi, double tol=0.0001) {
        Tensor<double> x(A.nrow());
        x[0] = 1.0;
        
        T xlo = (e1-e0)/(ehi-e1);
        T xhi = 1.0;
        //T pts[] = {0.0};
        //T pts[] = {std::sqrt(2.0)/2.0,-std::sqrt(2.0)/2.0};
        //T pts[] = {std::sqrt(3.0)/2.0,0.0,-std::sqrt(3.0)/2.0};
        T pts[] = {-std::sqrt(2.0+std::sqrt(2.0))/2.0,-std::sqrt(2.0-std::sqrt(2.0))/2.0,std::sqrt(2.0-std::sqrt(2.0))/2.0,std::sqrt(2.0+std::sqrt(2.0))/2.0};
        T eprev = 0.0;
        for (int iter=0; iter<1000; iter++) {
            T eval;
            for (T pt : pts) {
                T omega = 1.0/((0.5*(pt+1.0)*(xhi-xlo)+ xlo));
                //omega = 1.0;
                x.scale(1.0/x.normf());
                Tensor<T> Ax = sparse_matrix_times_vector(A, x);
                eval = x.trace(Ax);
                x = x + omega*((-1.0/Ax.normf())*Ax - x);
            }
            print(iter,eval);
            if ((iter%10)==1) {
                if (std::abs((eval-eprev)/eval) < tol) return eval;
                eprev = eval;
            }
        }
        return eprev;
    }
}

// Computes the Shur energy expression using dense matrix algebra
template <typename T>
T ShurDense(const Tensor<T>& H, T E) {
    size_t N = H.dim(0);

    // Reference state energy
    T E0 = H(0,0);

    // Extract Shur complement matrix and vector
    Tensor<T> A = copy(H(Slice(1,-1),Slice(1,-1)));
    Tensor<T> x = copy(H(0,Slice(1,-1)));

    // Shift and extract diagonal, and scale x
    Tensor<double> D(N-1);
    for (size_t i=0; i<N-1; i++) {
        T d = A(i,i) - E;
        A(i,i) = d;
        d = 1.0/std::sqrt(d);
        D(i) = d;
        x(i) *= d;
    }

    // Scale A
    for (size_t i=0; i<N-1; i++) {
        for (size_t j=0; j<N-1; j++) {
            A(i,j) *= D(i)*D(j);
        }
    }

    // Tensor<double> evals,evecs;
    // syev(A, evecs, evals);
    // print("dense Aevals",evals[0],evals[A.dim(0)-1]);

    // Solve the linear equation
    Tensor<T> y;
    gesv(A, x, y);

    // // Project x and y onto Aevecs out of curiosity
    // Tensor<T> xA = inner(x,evecs);
    // Tensor<T> yA = inner(y,evecs);
    // print("");
    // print("xA",xA);
    // print("yA",yA);
    // std::exit(0);

    // Compute new energy
    return E0 - x.trace(y);
}

// Computes the Shur energy expression using (mostly) sparse matrix algebra
// (we need to extend the sparse matrix class a little to avoid the dense operations)
template <typename T>
T ShurSparse(const Tensor<T>& H, T E) {
    size_t N = H.dim(0);

    // Reference state energy
    T E0 = H(0,0);

    // Extract Shur complement matrix and vector
    Tensor<T> A = copy(H(Slice(1,-1),Slice(1,-1)));
    Tensor<T> x = copy(H(0,Slice(1,-1)));

    // Shift and extract diagonal, and scale x
    Tensor<double> D(N-1);
    Tensor<double> s(N-1);
    for (size_t i=0; i<N-1; i++) {
        T d = A(i,i) - E;
        A(i,i) = d;
        d = 1.0/std::sqrt(d);
        D(i) = d;
        x(i) *= d;
        if (x(i)<0) {
            s(i) = -1.0;
            x(i) = -x(i);
        }
        else {
            s(i) = 1.0;
        }
    }

    // Scale A
    for (size_t i=0; i<N-1; i++) {
        for (size_t j=0; j<N-1; j++) {
            A(i,j) *= D(i)*D(j)*s(i)*s(j);
        }
    }

    // Convert to sparse
    sparse_matrix_csr<T> Asp = tensor_to_sparse_matrix(A,1e-6);

    // Iteratively solve the linear equation --- for now just bute force iteration --- really need to estimate the range of the spectrum using power method
    // then switch to faster iteration to solve.

    Tensor<T> y(N-1);
    double omega = 0.773;
    for (int iter=0; iter<1000; iter++) {
        Tensor<T> r = x - sparse_matrix_times_vector(Asp,y);
        double err = r.normf();
        //print(iter,err);
        y += omega*r;
        if (err < 1e-6) break;
    }

    // Compute new energy
    return E0 - x.trace(y);
}

// Computes the Shur energy expression using fully sparse matrix algebra
template <typename T>
T ShurSparse(const sparse_matrix_csr<T>& H, T E) {
    size_t N = H.nrow();

    // Reference state energy
    T E0 = H.get_row(0).get_data()[0];

    // Extract Shur complement matrix and vector
    sparse_matrix_csr<T> A = H(Slice(1,-1),Slice(1,-1));
    Tensor<T> x = copy(H.get_row(0).dense()(Slice(1,-1)));

    // Shift and extract diagonal, and scale x
    Tensor<double> D(N-1);
    for (size_t i=0; i<N-1; i++) {
        auto [found, d] = A.get_row(i).find(i);
        if (!found) throw "diagonal element is missing?";
        d = 1.0/std::sqrt(d-E);
        D(i) = d;
        x(i) *= d;
    }

    // Scale A and shift its diagonal (which becomes 1)
    for (size_t i=0; i<N-1; i++) {
        auto& rowi = A.get_row(i);
        auto& ind = rowi.get_indices();
        auto& dat = rowi.get_data();
        for (size_t jj=0; jj<ind.size(); jj++) {
            size_t j = ind[jj];
            if (i == j) {
                dat[jj] = 1.0;
            }
            else {
                dat[jj] *= D(i)*D(j);
            }
        }
    }

    // Iteratively solve the linear equation --- for now just bute force iteration --- really need to estimate the range of the spectrum using power method
    // then switch to faster iteration to solve.

    Tensor<T> y(N-1);
    double omega = 0.773;
    for (int iter=0; iter<1000; iter++) {
        Tensor<T> r = x - sparse_matrix_times_vector(A,y);
        double err = r.normf();
        //print(iter,err);
        y += omega*r;
        if (err < 1e-6) break;
    }

    // Compute new energy
    return E0 - x.trace(y);
}


// Computes the Shur energy expression using (mostly) sparse matrix algebra and Monte Carlo
// (we need to extend the sparse matrix class a little to avoid the dense operations)
template <typename T>
T ShurSparseMC(const Tensor<T>& H, T E) {
    size_t N = H.dim(0);

    // Reference state energy
    T E0 = H(0,0);

    // Extract Shur complement matrix and vector
    Tensor<T> A = copy(H(Slice(1,-1),Slice(1,-1)));
    Tensor<T> x = copy(H(0,Slice(1,-1)));

    // Shift and extract diagonal, and scale x
    Tensor<T> D(N-1);
    Tensor<T> s(N-1);
    for (size_t i=0; i<N-1; i++) {
        T d = A(i,i) - E;
        A(i,i) = d;
        d = 1.0/std::sqrt(d);
        D(i) = d;
        x(i) *= d;
        // if (x(i)<0) {  // experiment with using sign information of x
        //     s(i) = -1.0;
        //     x(i) = -x(i);
        // }
        // else {
           s(i) = 1.0;
        // }
    }

    // Find highest index non-zero element of x so we can optimize truncation etc.
    int ximax = N-2;
    while (x[ximax] == 0.0) ximax--;
    print("ximax", ximax);

    // Scale A
    for (size_t i=0; i<N-1; i++) {
        for (size_t j=0; j<N-1; j++) {
            A(i,j) *= D(i)*D(j);
        }
    }

    // Tensor<double> evals,evecs;
    // syev(A, evecs, evals);
    // print("Aevals",evals[0],evals[1],evals[2],evals[3],evals[A.dim(0)-1]);


    // For MC algorithm zero the diagonal entries of A which we know are all 1
    for (size_t i=0; i<N-1; i++) {
        if (std::abs(A(i,i)-1.0) > 1e-12) print("bad diag", i, A(i,i));
        A(i,i) = 0.0;
    }

    // // Examine spectra of +ve and -ve parts of A
    // Tensor<T> Apos=copy(A), Aneg=copy(A);
    // for (size_t i=0; i<N-1; i++) {
    //     for (size_t j=0; j<N-1; j++) {
    //         if (A(i,j)<0) {
    //             Apos(i,j) = 0.0;
    //         }
    //         else {
    //             Aneg(i,j) = 0.0;
    //         }
    //     }
    // }
    // Tensor<double> evals,evecs;
    // {syev(Aneg, evecs, evals); print("Aneg",evals[0],evals[A.dim(0)-1]);}
    // {syev(Apos, evecs, evals); print("Apos",evals[0],evals[A.dim(0)-1]);}
    // {syev(Apos+Aneg, evecs, evals); print("Apos+Aneg",evals[0],evals[A.dim(0)-1]);}
    // {syev(Apos-Aneg, evecs, evals); print("Apos-Aneg",evals[0],evals[A.dim(0)-1]);}
    
    // Convert to sparse
    sparse_matrix_csr<T> Asp = tensor_to_sparse_matrix(A,1e-6);

    double omega = 0.1; //0.773;
    Tensor<T> y = omega*copy(x); // can generate a better initial guess

    double sum = 0.0;
    double sumsq = 0.0;
    double sum2 = 0.0;
    double sumsq2 = 0.0;
    size_t count = 0;

    const double cut = 0.01;
    const double pkill = 0.75;

    {
        T eval0 = simple_power_iteration(Asp, -1.0);
        T eval1 = simple_power_iteration(Asp, eval0);
        print("Aspevals", eval0, eval1);
    }

    for (int iter=0; iter<3500; iter++) {
        Tensor<T> ynew(N-1);
        sparse_matrix_coo<double> Asample(N,N);
        for (size_t i=0; i<N-1; i++) {
            if (y[i]) {
                int nrep = 3;
                for (int rep=0; rep<nrep; rep++) {
                    // Randomly sample non-zero element in row A[i,*]
                    auto row = Asp.get_row(i);
                    size_t n = row.get_indices().size();
                    size_t j = n*RandomValue<double>();
                    size_t jj = row.get_indices()[j];
                    double aji = s(i)*s(jj)*row.get_data()[j];
                    Asample[{i,jj}] = aji;
                    Asample[{jj,i}] = aji;
                    double gji = double(nrep)/double(n);
                    ynew[jj] -= 0.5*omega*aji*y[i] / gji; 
                    ynew[i]  -= 0.5*omega*aji*y[jj] / gji; // Enforcing symmetry is important for stability (part of detailed balance?)
                }
                //for (auto [j,aji] : Asp.get_row(i)) ynew[j] -= omega*aji*y[i];
            }
            //ynew[i] += (1.0-omega)*y[i] + omega*x[i];
        }
        sparse_matrix_csr<double> B(Asample);
        T eval0 = simple_power_iteration(B, -1.0);
        T eval1 = simple_power_iteration(B, eval0);
        print("evals", eval0, eval1);

        //ynew -= omega*sparse_matrix_times_vector(Asp,y);
        double xzy = x.trace(x) + x.trace(ynew)/(omega); // using (1-z)^-1 = 1 + z (1-z)^-1
        ynew += (1.0-omega)*y + omega*x;
        double ecor = ynew.trace(x);
        //print(iter,ecor,(ynew-y).normf());
        y = ynew;
        if (iter >500) {
            sum += ecor;
            sumsq += ecor*ecor;
            sum2 += xzy;
            sumsq2 += xzy*xzy;
            count += 1;
        }

        size_t nkill = 0.0;
        for (size_t i=0; i<N-1; i++) {
            if (std::abs(y[i]) < cut) {
                if (RandomValue<double>() > pkill) {
                    y[i] *= 1.0/(1.0-pkill);
                }
                else {
                    nkill ++;
                    y[i] = 0.0;
                }
            }
        }
        //print("nkill", nkill);
    }
    sum /= count;
    sumsq /= count;
    double err = std::sqrt((sumsq - sum*sum)/count);
    sum2 /= count;
    sumsq2 /= count;
    double err2 = std::sqrt((sumsq2 - sum2*sum2)/count);
    print("Correlation energy ",-sum,"stderr",err,"total energy", E0 - sum, "F(E)", E0-sum-E);
    print("Correlation energy2",-sum2,"stderr",err2,"total energy", E0 - sum2, "F(E)", E0-sum2-E);
    //print(y);

    // Compute new energy
    return E0 - sum;
}

// Computes the Shur energy expression using fully sparse matrix algebra and Monte Carlo
template <typename T>
T ShurSparseMC(const sparse_matrix_csr<T>& H, T E) {
    size_t N = H.nrow();

    // Reference state energy
    T E0 = H.get_row(0).get_data()[0];

    // Extract Shur complement matrix and vector
    sparse_matrix_csr<T> A = H(Slice(1,-1),Slice(1,-1));
    Tensor<T> x = copy(H.get_row(0).dense()(Slice(1,-1)));

    // Shift and extract diagonal, and scale x
    Tensor<double> D(N-1);
    for (size_t i=0; i<N-1; i++) {
        auto [found, d] = A.get_row(i).find(i);
        if (!found) throw "diagonal element is missing?";
        d = 1.0/std::sqrt(d-E);
        D(i) = d;
        x(i) *= d;
    }

    // Scale A and shift its diagonal (which becomes 1 and here we set to zero since it is handled exactly)
    for (size_t i=0; i<N-1; i++) {
        auto& rowi = A.get_row(i);
        auto& ind = rowi.get_indices();
        auto& dat = rowi.get_data();
        for (size_t jj=0; jj<ind.size(); jj++) {
            size_t j = ind[jj];
            if (i == j) {
                dat[jj] = 0.0;
            }
            else {
                dat[jj] *= D(i)*D(j);
            }
        }
    }

    // Compress out the zeros (diagonal elements)
    A.compress();

    double omega = 0.1; //0.773;
    Tensor<T> y = omega*copy(x); // can generate a better initial guess

    double sum = 0.0;
    double sumsq = 0.0;
    double sum2 = 0.0;
    double sumsq2 = 0.0;
    size_t count = 0;

    const double cut = 0.001;
    const double pkill = 0.75;

    {
        T eval0 = simple_power_iteration(A, -1.0);
        T eval1 = simple_power_iteration(A, eval0);
        print("Aevals", eval0, eval1);
    }

    for (int iter=0; iter<6800; iter++) {
        Tensor<T> ynew(N-1);
        sparse_matrix_coo<double> Asample(N,N);
        for (size_t i=0; i<N-1; i++) {
            if (y[i]) {
                int nrep = 10;
                for (int rep=0; rep<nrep; rep++) {
                    // Randomly sample non-zero element in row A[i,*]
                    auto row = A.get_row(i);
                    size_t n = row.get_indices().size();
                    size_t j = n*RandomValue<double>();
                    size_t jj = row.get_indices()[j];
                    double aji = row.get_data()[j];
                    Asample[{i,jj}] = aji;
                    Asample[{jj,i}] = aji;
                    double gji = double(nrep)/double(n);
                    ynew[jj] -= 0.5*omega*aji*y[i] / gji; 
                    ynew[i]  -= 0.5*omega*aji*y[jj] / gji; // Enforcing symmetry is important for stability (part of detailed balance?)
                }
                //for (auto [j,aji] : A.get_row(i)) ynew[j] -= omega*aji*y[i];
            }
            //ynew[i] += (1.0-omega)*y[i] + omega*x[i];
        }
        // sparse_matrix_csr<double> B(Asample);
        // T eval0 = simple_power_iteration(B, -1.0);
        // T eval1 = simple_power_iteration(B, eval0);
        // print("evals", eval0, eval1);

        //ynew -= omega*sparse_matrix_times_vector(A,y);
        double xzy = x.trace(x) + x.trace(ynew)/(omega); // using (1-z)^-1 = 1 + z (1-z)^-1
        ynew += (1.0-omega)*y + omega*x;
        double ecor = ynew.trace(x);
        y = ynew;
        if (iter > 800) {
            sum += ecor;
            sumsq += ecor*ecor;
            sum2 += xzy;
            sumsq2 += xzy*xzy;
            count += 1;
        }

        size_t nkill = 0.0;
        for (size_t i=0; i<N-1; i++) {
            if (std::abs(y[i]) < cut) {
                if (RandomValue<double>() > pkill) {
                    y[i] *= 1.0/(1.0-pkill);
                }
                else {
                    nkill ++;
                    y[i] = 0.0;
                }
            }
        }
        if ((iter%10) == 1) print(iter,ecor,xzy,(ynew-y).normf(),nkill);
    }
    sum /= count;
    sumsq /= count;
    double err = std::sqrt((sumsq - sum*sum)/count);
    sum2 /= count;
    sumsq2 /= count;
    double err2 = std::sqrt((sumsq2 - sum2*sum2)/count);
    print("Correlation energy ",-sum,"stderr",err,"total energy", E0 - sum, "F(E)", E0-sum-E);
    print("Correlation energy2",-sum2,"stderr",err2,"total energy", E0 - sum2, "F(E)", E0-sum2-E);
    //print(y);

    // Compute new energy
    return E0 - sum;
}

int main() {
    std::cout.precision(10);    
    //sparse_matrix_csr A = load_ci_matrix("sparse.txt");
    sparse_matrix_csr A = load_ci_matrix("hooh-15252.txt");
    //print(A.get_row(0));
    
    //Tensor<double> B = sparse_matrix_to_tensor(A);
    // Tensor<double> evals,evecs;
    // syev(B, evecs, evals);
    // print("Hevals",evals[0],evals[B.dim(0)-1]);

    // Tensor<double> x(A.nrow());
    // x.fillrandom();
    // Tensor<double> Bx = inner(B,x);
    // Tensor<double> Ax = sparse_matrix_times_vector(A, x);
    // print("error in sparse mxv", (Ax - Bx).normf());

    // double elo = simple_power_iteration(A,-65.0,1e-4,false);
    // double ehi = simple_power_iteration(A,elo,1e-4);
    // print(elo,ehi);
    // elo -= 0.1;
    // ehi += 0.1;
    // print(elo,ehi);

    // power_iteration(A, elo, elo+0.2, ehi, 1e-7);

    //    print("sparse1 ", -76.1, ShurSparse(B,-76.1));
    double E0 = -150.807 - 0.4;
    print("sparse2 ", E0, ShurSparse(A,E0));
    //    //print(" dense", -76.1,ShurDense(B,-76.1));

    ShurSparseMC(A,E0);
    
    // Use secant method to find the root ... better would be to use all points computed and fit a low-order polyn to manage noise
    //double E0 = -76.1; // Initial guess
    double F0 = ShurSparse(A,E0) - E0;
    double E1 = E0 + 0.1;
    while (std::abs(E0-E1)>1e-4) {
        double F1 = ShurSparse(A,E1) - E1;
        print(E0,F0,E1,F1);
        double E2 = E1 - F1*(E1-E0)/(F1-F0);
        // Save the old point closest to the new point
        if (std::abs(E2-E1) < std::abs(E2-E0)) {
            E0 = E1; F0 = F1; E1 = E2;
        }
        else {
            E1 = E2;
        }
    }
    
    return 0;
}
