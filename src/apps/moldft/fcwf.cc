#include "fcwf.h"

using namespace madness;

     
Fcwf::Fcwf(){
     m_initialized = false;
}


Fcwf::Fcwf(const complex_function_3d& wf1,
     const complex_function_3d& wf2,
     const complex_function_3d& wf3,
     const complex_function_3d& wf4){
     MADNESS_ASSERT(m_psi.size() == 0);
     m_psi.push_back(wf1);
     m_psi.push_back(wf2);
     m_psi.push_back(wf3);
     m_psi.push_back(wf4);
     m_initialized = true;
}

Fcwf::Fcwf(World& world){
     MADNESS_ASSERT(m_psi.size() == 0);
     for(int i = 0 ; i < 4 ; i ++){
          m_psi.push_back(complex_factory_3d(world));
     }
     m_initialized = true;
}

complex_function_3d& Fcwf::operator[](const int i){
     MADNESS_ASSERT(i >= 0 && i <= 3);
     MADNESS_ASSERT(m_initialized);
     return m_psi[i];
}

const complex_function_3d& Fcwf::operator[](const int i) const {
     MADNESS_ASSERT(i >= 0 && i <= 3);
     MADNESS_ASSERT(m_initialized);
     return m_psi[i];
}

Fcwf::Fcwf(std::vector<complex_function_3d>& phi){
     MADNESS_ASSERT(m_psi.size() == 0);
     MADNESS_ASSERT(phi.size() == 4);
     for(int i = 0 ; i < 4 ; i++){
          m_psi.push_back(phi[i]);
     }
     m_initialized=true;
}

bool Fcwf::getinitialize(){
     return m_initialized;
}

bool Fcwf::getinitialize() const {
     return m_initialized;
}

unsigned int Fcwf::size(){
     MADNESS_ASSERT(m_initialized);
     return m_psi.size();
}

unsigned int Fcwf::size() const {
     MADNESS_ASSERT(m_initialized);
     return m_psi.size();
}

//copy contructor defaults to deep copy
//if this ever changes, you will need to change the copy() function, as it calls this
Fcwf::Fcwf(const Fcwf& phi){
     MADNESS_ASSERT(m_psi.size() == 0);
     MADNESS_ASSERT(phi.size() == 4);
     for(int i = 0 ; i < 4 ; i++){
          m_psi.push_back(copy(phi[i]));
     }
     m_initialized = true;
}

//Assignment operator defaults to shallow copy
Fcwf Fcwf::operator=(const Fcwf& phi){
     //MADNESS_ASSERT(phi.getinitialize());
     //if (this != &phi) {
     //     if(m_psi.size() == 4){
     //          for(int i = 0 ; i < 4 ; i++){
     //               m_psi[i] = copy(phi[i]);
     //          }
     //     }
     //     else {
     //          MADNESS_ASSERT(m_psi.size() == 0);
     //          for(int i = 0 ; i < 4 ; i++){
     //               m_psi.push_back(copy(phi[i]));
     //          }
     //     }
     //}
     //m_initialized = true;
     m_psi = phi.m_psi;
     m_initialized = phi.m_initialized;
     return *this;
}

Fcwf Fcwf::operator-(const Fcwf& phi) const {
     MADNESS_ASSERT(phi.getinitialize());
     std::vector<complex_function_3d> temp;
     if(m_initialized){
          //temp = madness::sub(m_psi[0].world(), m_psi, phi.m_psi);
          for(int i = 0 ; i < 4 ; i++){
               temp.push_back(m_psi[i] - phi[i]);
          }
     }
     else {
          for(int i = 0 ; i < 4 ; i++){
               temp.push_back(copy(phi[i]));
               temp[i].scale(-1.0);
          }
     }
     return Fcwf(temp);
}

Fcwf Fcwf::operator+(const Fcwf& phi){
     MADNESS_ASSERT(phi.getinitialize());
     std::vector<complex_function_3d> temp;
     if(m_initialized){
          //temp = madness::add(m_psi[0].world(), m_psi, phi.m_psi);
          for(int i = 0 ; i < 4 ; i++){
               temp.push_back(m_psi[i] + phi[i]);
          }
     }
     else {
          for(int i = 0 ; i < 4 ; i++){
               temp.push_back(copy(phi[i]));
          }
     }
     return Fcwf(temp);
}

Fcwf Fcwf::operator*(std::complex<double> a) const {
     MADNESS_ASSERT(m_initialized);
     std::vector<complex_function_3d> temp(4);
     for(int i = 0 ; i < 4 ; i++){
          temp[i] = a*m_psi[i];    
     }
     return Fcwf(temp);
}

void Fcwf::scale(std::complex<double> a){
     MADNESS_ASSERT(m_initialized);
     for(int i = 0 ; i < 4 ; i++){
          //m_psi[i] = a*m_psi[i];
          m_psi[i].scale(a);
     }
}

Fcwf Fcwf::operator+=(const Fcwf& phi){
     if(m_initialized){
          //m_psi = madness::add(m_psi[0].world(), m_psi, phi.m_psi);
          for(int i = 0 ; i < 4 ; i++){
               m_psi[i] += phi[i];
          }
     }
     else {
          MADNESS_ASSERT(m_psi.size()==0);
          for(int i = 0 ; i < 4 ; i++){
               m_psi.push_back(copy(phi[i]));
          }
          m_initialized = true;
     }
     return *this;
}

Fcwf Fcwf::operator-=(const Fcwf& phi){
     if(m_initialized){
          //m_psi = madness::sub(m_psi[0].world(), m_psi, phi.m_psi);
          for(int i = 0 ; i < 4 ; i++){
               m_psi[i] -= phi[i];
          }
     }
     else {
          MADNESS_ASSERT(m_psi.size()==0);
          for(int i = 0 ; i < 4 ; i++){
               m_psi.push_back(copy(phi[i]));
               m_psi[i].scale(-1.0);
          }
          m_initialized = true;
     }
     return *this;
}


double Fcwf::norm2(){
     MADNESS_ASSERT(m_initialized);
     std::complex<double> temp = madness::inner(m_psi[0].world(), m_psi, m_psi).sum();
     //std::complex<double> temp(0,0);
     //for(int i = 0 ; i < 4 ; i++){
     //     temp += madness::inner(m_psi[i],m_psi[i]);
     //}
     return std::sqrt(std::real(temp));
}

void Fcwf::normalize(){
     MADNESS_ASSERT(m_initialized);
     double norm = norm2();
     MADNESS_ASSERT(norm != 0.0);
     for(int i = 0 ; i < 4 ; i++){
          m_psi[i].scale(1.0/norm);// = (1.0/norm)*m_psi[i];
     }
}

Fcwf Fcwf::operator*(madness::complex_function_3d& phi){
     MADNESS_ASSERT(m_initialized);
     std::vector<complex_function_3d> temp(4);// = mul(m_psi[0].world(), phi, m_psi);
     for(int i = 0 ; i < 4 ; i++){
          temp[i] = phi*m_psi[i];
     }
     return Fcwf(temp);
}

Fcwf Fcwf::operator*(madness::real_function_3d& phi){
     MADNESS_ASSERT(m_initialized);
     std::vector<complex_function_3d> temp(4);// = madness::mul(m_psi[0].world(), phi, m_psi);
     for(int i = 0 ; i < 4 ; i++){
          temp[i] = phi*m_psi[i];
     }
     return Fcwf(temp);
}

void Fcwf::truncate(){
     MADNESS_ASSERT(m_initialized);
     for(int i = 0 ; i < 4 ; i++){
          m_psi[i].truncate();    
     }
}

std::complex<double> Fcwf::inner(World& world, const Fcwf& phi) const{
     MADNESS_ASSERT(m_initialized && phi.getinitialize());
     return madness::inner(world, m_psi, phi.m_psi).sum();
}

void Fcwf::apply(World& world, real_convolution_3d& op){
     m_psi = madness::apply(world, op, m_psi);
}

void Fcwf::apply(World& world, complex_derivative_3d& D){
     m_psi = madness::apply(world, D, m_psi);
}

std::complex<double> inner(const Fcwf& psi, const Fcwf& phi){
     //std::complex<double> result(0,0);
     //for(int i = 0 ; i < 4 ; i++){
     //     result += madness::inner(psi[i],phi[i]);
     //}
     //return result;
     MADNESS_ASSERT(psi.getinitialize() && phi.getinitialize());
     return psi.inner(psi[0].world(), phi);
}

//Fcwf apply(real_convolution_3d& op, const Fcwf& psi){
//     std::vector<complex_function_3d> temp;
//     for(int i = 0 ; i < 4 ; i++){
//          temp.push_back(madness::apply(op, psi[i]));
//     }
//     return Fcwf(temp);
//}
Fcwf apply(World& world, real_convolution_3d& op, const Fcwf& psi){
     Fcwf temp = copy(psi);
     temp.apply(world, op);
     return temp;
}

Fcwf apply(World& world, complex_derivative_3d& D, const Fcwf& psi){
     Fcwf temp = copy(psi);
     temp.apply(world, D);
     return temp;
}

real_function_3d squaremod(Fcwf& psi){
     MADNESS_ASSERT(psi.getinitialize());
     real_function_3d temp = abssq(psi[0]) + abssq(psi[1]) + abssq(psi[2]) + abssq(psi[3]);
     return temp;
}

real_function_3d squaremod_small(Fcwf& psi){
     MADNESS_ASSERT(psi.getinitialize());
     real_function_3d temp = abssq(psi[2]) + abssq(psi[3]);
     return temp;
}

real_function_3d squaremod_large(Fcwf& psi){
     MADNESS_ASSERT(psi.getinitialize());
     real_function_3d temp = abssq(psi[0]) + abssq(psi[1]);
     return temp;
}

complex_function_3d inner_func(World& world, Fcwf& psi, Fcwf& phi){
     MADNESS_ASSERT(psi.getinitialize() && phi.getinitialize());
     std::vector<complex_function_3d> a(4);
     std::vector<complex_function_3d> b(4);
     for(unsigned int i = 0; i < 4; i++){
          a[i] = psi[i];
          b[i] = phi[i];
     }
     complex_function_3d result = sum(world, mul(world, conj(world, a), b)); 
     return result;
}

Fcwf copy(Fcwf psi){
     return Fcwf(psi);

}

std::complex<double> inner(std::vector<Fcwf>& a, std::vector<Fcwf>& b){
     MADNESS_ASSERT(a.size() == b.size());
     std::complex<double> result(0,0);
     for(int i = 0; i < a.size(); i++){
          result += inner(a[i],b[i]);    
     }
     return result;
}

std::vector<Fcwf> operator*(const std::vector<Fcwf>& psis, std::complex<double> a){
     std::vector<Fcwf> result;
     if(psis.size() != 0){
          for(int i = 0; i < psis.size(); i++){
               result.push_back(psis[i]*a);
          }
     }
     return result;
}

std::vector<Fcwf> operator*(std::complex<double> a, const std::vector<Fcwf>& psis){
     std::vector<Fcwf> result;
     if(psis.size() != 0){
          for(int i = 0; i < psis.size(); i++){
               result.push_back(psis[i]*a);
          }
     }
     return result;
}

void operator+=(std::vector<Fcwf>& phi, const std::vector<Fcwf>& psi){
     if(phi.size()==0){
          phi = psi;
     }
     else if(psi.size() != 0){
          MADNESS_ASSERT(phi.size()==psi.size());
          for(int i=0; i < psi.size(); i++){
               phi[i]+=psi[i];
          }
     }
}

std::vector<Fcwf> operator-(const std::vector<Fcwf>& phi, const std::vector<Fcwf>& psi){
     std::vector<Fcwf> result;
     if(psi.size()==0){
          result = phi;
     }
     else if(phi.size()==0){
          result = -1.0*psi;
     }
     else{
          MADNESS_ASSERT(phi.size()==psi.size());
          for(int i=0; i < psi.size(); i++){
               result.push_back(phi[i]-psi[i]);
          }
     }
     return result;
}

//Constructor
Fcwf_vector_allocator::Fcwf_vector_allocator(World& world, unsigned int m_size)
: world(world)
, m_size(m_size)
{}

//Overloading () operator
std::vector<Fcwf> Fcwf_vector_allocator::operator()(){
     std::vector<Fcwf> result;
     for(int i=0; i < m_size; i++){
          result.push_back(Fcwf(world));
     }
     return result;
}

//Copy Constructor
//according to Bryan, according to Jakob, this is necessary for KAIN?
Fcwf_vector_allocator Fcwf_vector_allocator::operator=(const Fcwf_vector_allocator& other){
     Fcwf_vector_allocator tmp(world, other.m_size);
     return tmp;
}

void Fcwf_vector_allocator::set_size(int size){
     m_size = size;
}

//Forms the outer product between two vectors of Fcwfs, where each matrix element is the inner product of the two contributing Fcwfs.
Tensor<std::complex<double>> matrix_inner(World& world, std::vector<Fcwf>& a, std::vector<Fcwf>& b){
     unsigned int n = a.size();
     unsigned int m = b.size();
     MADNESS_ASSERT(n==m);

     std::vector<complex_function_3d> a_1(n);
     std::vector<complex_function_3d> a_2(n);
     std::vector<complex_function_3d> a_3(n);
     std::vector<complex_function_3d> a_4(n);
     std::vector<complex_function_3d> b_1(n);
     std::vector<complex_function_3d> b_2(n);
     std::vector<complex_function_3d> b_3(n);
     std::vector<complex_function_3d> b_4(n);

     for(unsigned int i = 0; i < n; i++){
          a_1[i] = a[i][0];
          a_2[i] = a[i][1];
          a_3[i] = a[i][2];
          a_4[i] = a[i][3];
          b_1[i] = b[i][0];
          b_2[i] = b[i][1];
          b_3[i] = b[i][2];
          b_4[i] = b[i][3];
     }

     Tensor<std::complex<double>> component1 = matrix_inner(world, a_1, b_1);
     Tensor<std::complex<double>> component2 = matrix_inner(world, a_2, b_2);
     Tensor<std::complex<double>> component3 = matrix_inner(world, a_3, b_3);
     Tensor<std::complex<double>> component4 = matrix_inner(world, a_4, b_4);
     component1=component1+component2+component3+component4;
     return component1;
}

void transform(World& world, std::vector<Fcwf>& a, Tensor<std::complex<double>> U){
     unsigned int n = a.size();
     unsigned int m = U.dim(0);
     unsigned int k = U.dim(1);
     MADNESS_ASSERT(n==m);
     MADNESS_ASSERT(m==k); //for now only support square transformation

     std::vector<complex_function_3d> a_1(n);
     std::vector<complex_function_3d> a_2(n);
     std::vector<complex_function_3d> a_3(n);
     std::vector<complex_function_3d> a_4(n);
     for(unsigned int i = 0; i < n; i++){
          a_1[i] = a[i][0];
          a_2[i] = a[i][1];
          a_3[i] = a[i][2];
          a_4[i] = a[i][3];
     }
     a_1 = transform(world, a_1, U);
     a_2 = transform(world, a_2, U);
     a_3 = transform(world, a_3, U);
     a_4 = transform(world, a_4, U);

     for(unsigned int i = 0; i < n; i++){
          a[i][0] = a_1[i];
          a[i][1] = a_2[i];
          a[i][2] = a_3[i];
          a[i][3] = a_4[i];
     }
     
}

//loop through fcwf and reconstruct each function
void Fcwf::reconstruct(){
     MADNESS_ASSERT(m_initialized); //needed?
     for(unsigned int i = 0; i < 4; i++) m_psi[i].reconstruct();
}

//loop through fcwf and compress each function
void Fcwf::compress(){
     MADNESS_ASSERT(m_initialized); //needed?
     for(unsigned int i = 0; i < 4; i++) m_psi[i].compress();
}

//kthxbye
