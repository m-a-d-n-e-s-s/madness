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
        
        // Return the indices of the non-zero elements
        const std::vector<size_t>& get_indices() const {return indices;}
        std::vector<size_t>& get_indices() {return indices;}

        // Return the values of the non-zero elements
        const std::vector<T>& get_data() const {return data;}
        std::vector<T>& get_data() {return data;}

        // Returns the length (if dense) of the vector
        size_t size() const {return M;}
        
        // Push (i,value) onto the end of the indices and data
        void push_back(size_t i, T value) {
            MADNESS_ASSERT(i>=0 && i<M);
            indices.push_back(i);
            data.push_back(value);
        }
        
        // Sort the index vector into increasing order
        void sort_indices() {
            using pairT = std::pair<size_t,T>;
            std::vector<pairT> x;
            for (size_t i=0; i<data.size(); i++) x.push_back({indices[i],data[i]});
            std::sort(x.begin(), x.end(), [](const pairT& a, const pairT& b){return a.first<b.first;});
            indices = std::vector<size_t>(indices.size());
            data = std::vector<T>(data.size());
            for (size_t i=0; i<data.size(); i++) {indices[i]=x[i].first,data[i]=x[i].second;}
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
    };

    // Read A test matrix from the configuration interaction (CI) code
    sparse_matrix_coo<double> load_ci_matrix(const std::string& filename) {
        std::fstream file(filename,std::ios::in);
        size_t N;
        file >> N;
        sparse_matrix_coo<double> A(N,N);
        size_t i, j;
        double value;
        size_t count=0;
        while (file >> i >> j >> value) {
            i--; j--; // fortran indexes from 1, c++ from 1
            MADNESS_ASSERT(i>=0 && i<N && j>=0 && j<N);
            A[{i,j}] = value; A[{j,i}] = value; // read full square
            count++;
        }
        print("load_ci_matrix:",N,count);
        
        sparse_matrix_csr<double> B(A);
        print(B.nrow());
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

    // Convert to sparse
    sparse_matrix_csr<T> Asp = tensor_to_sparse_matrix(A,1e-6);

    // Iteratively solve the linear equation --- for now just bute force iteration --- really need to estimate the range of the spectrum using power method
    // then switch to faster iteration to solve.

    Tensor<T> y(N-1);
    double omega = 0.773;
    for (int iter=0; iter<1000; iter++) {
        Tensor<T> r = x - sparse_matrix_times_vector(Asp,y);
        double err = r.normf();
        print(iter,err);
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
    // print("Aevals",evals[0],evals[A.dim(0)-1]);

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
    size_t count = 0;
    
    for (int iter=0; iter<500; iter++) {
        Tensor<T> ynew(N-1);
        for (size_t i=0; i<N-1; i++) {
            if (y[i]) {
                int nrep = 10;
                for (int rep=0; rep<nrep; rep++) {
                    // Randomly sample non-zero element in row A[i,*]
                    auto row = Asp.get_row(i);
                    size_t n = row.get_indices().size();
                    size_t j = n*RandomValue<double>();
                    double aji = row.get_data()[j];
                    double gji = double(nrep)/double(n);
                    ynew[row.get_indices()[j]] -= omega*aji*y[i] / gji;
                }
                //for (auto [j,aji] : Asp.get_row(i)) ynew[j] -= omega*aji*y[i];
            }
            //ynew[i] += (1.0-omega)*y[i] + omega*x[i];
        }

        //ynew -= omega*sparse_matrix_times_vector(Asp,y);
        ynew += (1.0-omega)*y + omega*x;
        double ecor = ynew.trace(x);
        print(iter,ecor,(ynew-y).normf());
        y = ynew;
        if (iter >250) {
            sum += ecor;
            sumsq += ecor*ecor;
            count += 1;
        }

        const double cut = 0.01;
        const double pkill = 0.75;
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
        print("nkill", nkill);
    }
    sum /= count;
    sumsq /= count;
    double err = std::sqrt((sumsq - sum*sum)/count);
    print("Correlation energy",-sum,"stderr",err,"total energy", E0 - sum, "F(E)", E0-sum-E);
    //print(y);

    // Compute new energy
    return E0 - sum;
}

int main() {
    std::cout.precision(10);    
    sparse_matrix_csr A = load_ci_matrix("sparse.txt");
    //print(A.get_row(0));
    
    Tensor<double> B = sparse_matrix_to_tensor(A);
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

    // print("sparse", -76.1,ShurSparse(B,-76.1));
    // print(" dense", -76.1,ShurDense(B,-76.1));

    ShurSparseMC(B,-76.1);
    
    // exact evals [-75.95318203 -75.78742742  ... -65.77660219]
    // for (int i=0; i<20; i++) {
    //     double E = -75.95318203 - i*0.01;
    //     print(E,ShurDense(B,E));
    // }

    // Use secant method to find the root ... better would be to use all points computed and fit a low-order polyn to manage noise
    double E0 = -76.1; // Initial guess
    double F0 = ShurDense(B,E0) - E0;
    double E1 = E0 + 0.1;
    while (std::abs(E0-E1)>1e-4) {
        double F1 = ShurDense(B,E1) - E1;
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
