#include "fcwf.h"

//Fcwf (Four component wavefunction) implementation file
//NOTE: The small component of an Fcwf is scaled by c and then stored
//This means that when computing properties and some other operations,
//care must be taken to remove the extra factors of c in the result

using namespace madness;

//Default constructor can't do much other than state that the fcwf hasn't been initialized 
Fcwf::Fcwf(){
     m_initialized = false;
}

//Constructor that takes four complex functions and uses them to make the Fcwf
Fcwf::Fcwf(const complex_function_3d& wf1,
     const complex_function_3d& wf2,
     const complex_function_3d& wf3,
     const complex_function_3d& wf4,
     const double myc){
     MADNESS_ASSERT(m_psi.size() == 0);
     m_psi.push_back(wf1);
     m_psi.push_back(wf2);
     m_psi.push_back(wf3);
     m_psi.push_back(wf4);
     m_initialized = true;
     speed_of_light = myc;
}

//This constructor creates a zero Fcwf
Fcwf::Fcwf(World& world, const double myc){
     MADNESS_ASSERT(m_psi.size() == 0);
     for(int i = 0 ; i < 4 ; i ++){
          m_psi.push_back(complex_factory_3d(world));
     }
     m_initialized = true;
     speed_of_light = myc;
}

//give access to the individual components
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

//Can also initialize from a vector of complex functions
Fcwf::Fcwf(std::vector<complex_function_3d>& phi, const double myc){
     MADNESS_ASSERT(m_psi.size() == 0);
     MADNESS_ASSERT(phi.size() == 4);
     for(int i = 0 ; i < 4 ; i++){
          m_psi.push_back(phi[i]);
     }
     m_initialized=true;
     speed_of_light = myc;
}

//learn whether an Fcwf is initialized
bool Fcwf::getinitialize(){
     return m_initialized;
}

bool Fcwf::getinitialize() const {
     return m_initialized;
}

// get speed of light
double Fcwf::get_myc(){
     MADNESS_ASSERT(m_initialized);
     return speed_of_light;
}

double Fcwf::get_myc() const {
     MADNESS_ASSERT(m_initialized);
     return speed_of_light;
}

//Probably don't need this, but get the size of m_psi, which should only be 0 or 4
unsigned int Fcwf::size(){
     MADNESS_ASSERT(m_initialized);
     return m_psi.size();
}

unsigned int Fcwf::size() const {
     MADNESS_ASSERT(m_initialized);
     return m_psi.size();
}

//copy constructor defaults to deep copy
//if this ever changes, you will need to change the copy() function, as it calls this
Fcwf::Fcwf(const Fcwf& phi, const double myc){
     MADNESS_ASSERT(m_psi.size() == 0);
     MADNESS_ASSERT(phi.size() == 4);
     for(int i = 0 ; i < 4 ; i++){
          m_psi.push_back(copy(phi[i]));
     }
     m_initialized = true;
     speed_of_light = myc;
}

//Assignment operator defaults to shallow copy
Fcwf Fcwf::operator=(const Fcwf& phi){
     m_psi = phi.m_psi;
     m_initialized = phi.m_initialized;
     speed_of_light = phi.speed_of_light;
     return *this;
}

//subtract two Fcwfs. The Fcwf being subtracted must be initialized
Fcwf Fcwf::operator-(const Fcwf& phi) const {
     MADNESS_ASSERT(phi.getinitialize());
     std::vector<complex_function_3d> temp;
     if(m_initialized){
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
     return Fcwf(temp, phi.speed_of_light);
}

//add two Fcwfs. The Fcwf being added must be initialized
Fcwf Fcwf::operator+(const Fcwf& phi){
     MADNESS_ASSERT(phi.getinitialize());
     std::vector<complex_function_3d> temp;
     if(m_initialized){
          for(int i = 0 ; i < 4 ; i++){
               temp.push_back(m_psi[i] + phi[i]);
          }
     }
     else {
          for(int i = 0 ; i < 4 ; i++){
               temp.push_back(copy(phi[i]));
          }
     }
     return Fcwf(temp, phi.speed_of_light);
}

//multiply an Fcwf by a complex number a
Fcwf Fcwf::operator*(std::complex<double> a) const {
     MADNESS_ASSERT(m_initialized);
     std::vector<complex_function_3d> temp(4);
     for(int i = 0 ; i < 4 ; i++){
          temp[i] = a*m_psi[i];    
     }
     return Fcwf(temp, speed_of_light);
}

//scale an Fcwf in place by a complex number a
void Fcwf::scale(std::complex<double> a){
     MADNESS_ASSERT(m_initialized);
     for(int i = 0 ; i < 4 ; i++){
          m_psi[i].scale(a);
     }
}

//in place addition of Fcwfs. Fcwf on right must be initialized
Fcwf Fcwf::operator+=(const Fcwf& phi){
     MADNESS_ASSERT(phi.getinitialize());
     if(m_initialized){
          for(int i = 0 ; i < 4 ; i++){
               m_psi[i] += phi[i];
          }
     }
     else {
          for(int i = 0 ; i < 4 ; i++){
               m_psi.push_back(copy(phi[i]));
          }
          m_initialized = true;
          speed_of_light = phi.speed_of_light;
     }
     return *this;
}

//in place subtraction of Fcwfs. Fcwf on right must be initialized
Fcwf Fcwf::operator-=(const Fcwf& phi){
     MADNESS_ASSERT(phi.getinitialize());
     if(m_initialized){
          for(int i = 0 ; i < 4 ; i++){
               m_psi[i] -= phi[i];
          }
     }
     else {
          for(int i = 0 ; i < 4 ; i++){
               m_psi.push_back(copy(phi[i]));
               m_psi[i].scale(-1.0);
          }
          m_initialized = true;
          speed_of_light = phi.speed_of_light;
     }
     return *this;
}

//Returns the 2-norm of an initialized Fcwf
double Fcwf::norm2(){
     MADNESS_ASSERT(m_initialized);
     double c2 = speed_of_light * speed_of_light; //speed of light in atomic units from CODATA 2022
     std::complex<double> temp(0,0);

     temp += madness::inner(m_psi[0],m_psi[0]);
     temp += madness::inner(m_psi[1],m_psi[1]);

     //Small component is stored with an additional factor of c, so remove it here
     temp += madness::inner(m_psi[2],m_psi[2])/c2;
     temp += madness::inner(m_psi[3],m_psi[3])/c2;

     return std::sqrt(std::real(temp));
}

//Normalize the input Fcwf
void Fcwf::normalize(){
     MADNESS_ASSERT(m_initialized);
     double norm = norm2();
     MADNESS_ASSERT(norm != 0.0);
     for(int i = 0 ; i < 4 ; i++){
          m_psi[i].scale(1.0/norm);
     }
}

//multiply an Fcwf by a complex function
Fcwf Fcwf::operator*(madness::complex_function_3d& phi){
     MADNESS_ASSERT(m_initialized);
     std::vector<complex_function_3d> temp(4);
     for(int i = 0 ; i < 4 ; i++){
          temp[i] = phi*m_psi[i];
     }
     return Fcwf(temp, speed_of_light);
}

Fcwf Fcwf::operator*(madness::real_function_3d& phi){
     MADNESS_ASSERT(m_initialized);
     std::vector<complex_function_3d> temp(4);
     for(int i = 0 ; i < 4 ; i++){
          temp[i] = phi*m_psi[i];
     }
     return Fcwf(temp, speed_of_light);
}

//truncate
void Fcwf::truncate(){
     MADNESS_ASSERT(m_initialized);
     for(int i = 0 ; i < 4 ; i++){
          m_psi[i].truncate();    
     }
}

//Returns the inner product of two Fcwfs
std::complex<double> Fcwf::inner(World& world, const Fcwf& phi) const{
     MADNESS_ASSERT(m_initialized && phi.getinitialize());
     double c2 = speed_of_light * speed_of_light; //speed of light in atomic units from CODATA 2022
     std::complex<double> temp(0,0);

     temp += madness::inner(m_psi[0],phi.m_psi[0]);
     temp += madness::inner(m_psi[1],phi.m_psi[1]);

     //Small component is stored with an extra factor of c, so remove that here
     temp += madness::inner(m_psi[2],phi.m_psi[2])/c2;
     temp += madness::inner(m_psi[3],phi.m_psi[3])/c2;

     return temp;
}

//Apply an integral operator to an Fcwf
void Fcwf::apply(World& world, real_convolution_3d& op){
     m_psi = madness::apply(world, op, m_psi);
}

//Apply a derivative operator to an Fcwf
void Fcwf::apply(World& world, complex_derivative_3d& D){
     m_psi = madness::apply(world, D, m_psi);
}

//Return result of time-reversal of the input Fcwf
Fcwf Fcwf::KramersPair(){
     MADNESS_ASSERT(m_initialized);
     complex_function_3d phi0 = -1.0*conj(m_psi[1]);
     complex_function_3d phi1 = conj(m_psi[0]);
     complex_function_3d phi2 = -1.0*conj(m_psi[3]);
     complex_function_3d phi3 = conj(m_psi[2]);
     return Fcwf(phi0,phi1,phi2,phi3,speed_of_light);
}

//function for computing the inner product of two Fcwfs
std::complex<double> inner(const Fcwf& psi, const Fcwf& phi){
     MADNESS_ASSERT(psi.getinitialize() && phi.getinitialize());
     return psi.inner(psi[0].world(), phi);
}

//function for applying an integral operator fo an Fcwf
Fcwf apply(World& world, real_convolution_3d& op, const Fcwf& psi){
     Fcwf temp = copy(psi);
     temp.apply(world, op);
     return temp;
}

//function for applying an integral operator fo an Fcwf
Fcwf apply(World& world, complex_derivative_3d& D, const Fcwf& psi){
     Fcwf temp = copy(psi);
     temp.apply(world, D);
     return temp;
}

//Returns the square modulus of an Fcwf, which is a real function
real_function_3d squaremod(Fcwf& psi){
     MADNESS_ASSERT(psi.getinitialize());
     double c2 = psi.get_myc() * psi.get_myc(); //speed of light in atomic units from CODATA 2022
     real_function_3d temp = abssq(psi[0]) + abssq(psi[1]) + abssq(psi[2]).scale(1.0/c2) + abssq(psi[3]).scale(1.0/c2);
     return temp;
}

//Returns the square modulus of the small component of an Fcwf, which is a real function
real_function_3d squaremod_small(Fcwf& psi){
     MADNESS_ASSERT(psi.getinitialize());
     double c2 = psi.get_myc() * psi.get_myc(); //speed of light in atomic units from CODATA 2022
     real_function_3d temp = (abssq(psi[2]) + abssq(psi[3])).scale(1.0/c2); //don't forget the factor of c^2
     return temp;
}

//Returns the square modulus of the large component of an Fcwf, which is a real function
real_function_3d squaremod_large(Fcwf& psi){
     MADNESS_ASSERT(psi.getinitialize());
     real_function_3d temp = abssq(psi[0]) + abssq(psi[1]);
     return temp;
}

//compute the function inner product between two Fcwfs. Result is a complex function.
complex_function_3d inner_func(World& world, Fcwf& psi, Fcwf& phi){
     MADNESS_ASSERT(psi.getinitialize() && phi.getinitialize());
     double c = psi.get_myc(); //speed of light in atomic units from CODATA 2022
     std::vector<complex_function_3d> a(4);
     std::vector<complex_function_3d> b(4);
     for(unsigned int i = 0; i < 2; i++){
          a[i] = psi[i];
          b[i] = phi[i];
     }
     //Small components are stored with an extra factor of c, so remove those here
     for(unsigned int i = 2; i < 4; i++){
          a[i] = copy(psi[i]).scale(1.0/c);
          b[i] = copy(phi[i]).scale(1.0/c);
     }
     //vmra function call takes care of the rest
     complex_function_3d result = sum(world, mul(world, conj(world, a), b)); 
     return result;
}

//deep copy of an Fcwf
Fcwf copy(Fcwf psi){
     MADNESS_ASSERT(psi.getinitialize());
     return Fcwf(psi, psi.get_myc());

}

//takes the inner product between two vectors of Fcwfs. Result is a complex number
std::complex<double> inner(std::vector<Fcwf>& a, std::vector<Fcwf>& b){
     MADNESS_ASSERT(a.size() == b.size());
     std::complex<double> result(0,0);
     for(int i = 0; i < a.size(); i++){
          result += inner(a[i],b[i]);    
     }
     return result;
}

//Multiply a vector of Fcwfs by a scalar
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

//In-place addition for a vector of Fcwfs
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

//Subtraction for two vectors of Fcwfs
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

//Constructor for allocator for vector of Fcwfs
Fcwf_vector_allocator::Fcwf_vector_allocator(World& world, unsigned int m_size, const double& myc)
: world(world)
, m_size(m_size)
, speed_of_light(myc)
{}

//Overloading () operator
std::vector<Fcwf> Fcwf_vector_allocator::operator()(){
     std::vector<Fcwf> result;
     for(int i=0; i < m_size; i++){
          result.push_back(Fcwf(world, speed_of_light));
     }
     return result;
}

//Copy Constructor for allocator. Necessary for KAIN
Fcwf_vector_allocator Fcwf_vector_allocator::operator=(const Fcwf_vector_allocator& other){
     Fcwf_vector_allocator tmp(world, other.m_size, speed_of_light);
     return tmp;
}

void Fcwf_vector_allocator::set_size(int size){
     m_size = size;
}

//Forms the outer product between two vectors of Fcwfs, where each matrix element is the inner product of the two contributing Fcwfs.
Tensor<std::complex<double>> matrix_inner(World& world, std::vector<Fcwf>& a, std::vector<Fcwf>& b){
     unsigned int n = a.size();
     unsigned int m = b.size();
     //MADNESS_ASSERT(n==m);

     double c2 = a[0].get_myc() * a[0].get_myc(); //speed of light in atomic units from CODATA 2022

     //Reassign the vectors of Fcwfs to vectors of complex functions to facilitate use of vmra functions
     std::vector<complex_function_3d> a_1(n); //all first components of Fcwfs in input a
     std::vector<complex_function_3d> a_2(n); //all second components of Fcwfs in input a
     std::vector<complex_function_3d> a_3(n);
     std::vector<complex_function_3d> a_4(n);
     std::vector<complex_function_3d> b_1(m); //all first components of Fcwfs in input b
     std::vector<complex_function_3d> b_2(m);
     std::vector<complex_function_3d> b_3(m);
     std::vector<complex_function_3d> b_4(m);

     for(unsigned int i = 0; i < n; i++){
          a_1[i] = a[i][0];
          a_2[i] = a[i][1];
          a_3[i] = a[i][2];
          a_4[i] = a[i][3];
     }
     for(unsigned int i = 0; i < m; i++){
          b_1[i] = b[i][0];
          b_2[i] = b[i][1];
          b_3[i] = b[i][2];
          b_4[i] = b[i][3];
     }

     //create 4 matrices which are the inner products of the components and add them
     Tensor<std::complex<double>> component1 = matrix_inner(world, a_1, b_1);
     Tensor<std::complex<double>> component2 = matrix_inner(world, a_2, b_2);

     //don't forget that small components are scaled by c
     Tensor<std::complex<double>> component3 = (1.0/c2)*matrix_inner(world, a_3, b_3);
     Tensor<std::complex<double>> component4 = (1.0/c2)*matrix_inner(world, a_4, b_4);

     //add
     component1=component1+component2+component3+component4;

     //return
     return component1;
}

//Transform a (row) vector of Fcwfs by right-multiplying by a transformation matrix U
std::vector<Fcwf> transform(World& world, std::vector<Fcwf>& a, Tensor<std::complex<double>> U){
     unsigned int n = a.size();
     unsigned int m = U.dim(0);
     unsigned int k = U.dim(1);
     MADNESS_ASSERT(n==m); //make sure dimensions align

     //Make vectors of the individual components to facilitate use of vmra
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

     //transform the component vectors
     a_1 = transform(world, a_1, U);
     a_2 = transform(world, a_2, U);
     a_3 = transform(world, a_3, U);
     a_4 = transform(world, a_4, U);

     //Now put the components back into Fcwf form and push back into a vector
     std::vector<Fcwf> result;
     Fcwf reader(world, a[0].get_myc()); 
     for(unsigned int i = 0; i < k; i++){
          reader[0] = a_1[i];
          reader[1] = a_2[i];
          reader[2] = a_3[i];
          reader[3] = a_4[i];
          result.push_back(reader); 
     }

     return result;
     
}

//loop through fcwf and reconstruct each function
void Fcwf::reconstruct(){
     MADNESS_ASSERT(m_initialized); 
     for(unsigned int i = 0; i < 4; i++) m_psi[i].reconstruct();
}

//loop through fcwf and compress each function
void Fcwf::compress(){
     MADNESS_ASSERT(m_initialized); 
     for(unsigned int i = 0; i < 4; i++) m_psi[i].compress();
}

//kthxbye
