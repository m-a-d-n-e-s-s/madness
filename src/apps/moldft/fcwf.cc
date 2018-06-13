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
          m_psi.push_back(copy(wf1));
          m_psi.push_back(copy(wf2));
          m_psi.push_back(copy(wf3));
          m_psi.push_back(copy(wf4));
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
               m_psi.push_back(copy(phi[i]));
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
     Fcwf::Fcwf(const Fcwf& phi){
          MADNESS_ASSERT(m_psi.size() == 0);
          MADNESS_ASSERT(phi.size() == 4);
          for(int i = 0 ; i < 4 ; i++){
               m_psi.push_back(copy(phi[i]));
          }
          m_initialized = true;
     }

     //Assignment operator defaults to deep copy
     Fcwf Fcwf::operator=(const Fcwf& phi){
          MADNESS_ASSERT(phi.getinitialize());
          if (this != &phi) {
               if(m_psi.size() == 4){
                    for(int i = 0 ; i < 4 ; i++){
                         m_psi[i] = copy(phi[i]);
                    }
               }
               else {
                    MADNESS_ASSERT(m_psi.size() == 0);
                    for(int i = 0 ; i < 4 ; i++){
                         m_psi.push_back(copy(phi[i]));
                    }
               }
          }
          m_initialized = true;
          return *this;
     }

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
          return Fcwf(temp);
     }

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
          return Fcwf(temp);
     }

     Fcwf Fcwf::operator*(std::complex<double> a){
          //print("multiply1");
          MADNESS_ASSERT(m_initialized);
          std::vector<complex_function_3d> temp(4);
          //print("multiply2");
          for(int i = 0 ; i < 4 ; i++){
               //print("multiply3, ",i);
               temp[i] = a*m_psi[i];    
          }
          //print("multiply4");
          return Fcwf(temp);
     }
     
     void Fcwf::scale(std::complex<double> a){
          MADNESS_ASSERT(m_initialized);
          for(int i = 0 ; i < 4 ; i++){
               m_psi[i] = a*m_psi[i];
          }
     }
     
     Fcwf Fcwf::operator+=(const Fcwf& phi){
          if(m_initialized){
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
               //print("inpl substract 1");
               for(int i = 0 ; i < 4 ; i++){
                    //print("inpl subtract 2, ", i);
                    m_psi[i] -= phi[i];
               }
          }
          else {
               //print("inpl subtract 3");
               MADNESS_ASSERT(m_psi.size()==0);
               for(int i = 0 ; i < 4 ; i++){
                    //print("inpl subtract 4, ",i);
                    m_psi.push_back(copy(phi[i]));
                    //print("inpl subtract 5, ",i);
                    m_psi[i].scale(-1.0);
               }
               //print("inpl subtract 6");
               m_initialized = true;
          }
          //print("inpl subtract 7");
          return *this;
     }


     double Fcwf::norm2(){
          MADNESS_ASSERT(m_initialized);
          std::complex<double> temp(0,0);
          for(int i = 0 ; i < 4 ; i++){
               temp += madness::inner(m_psi[i],m_psi[i]);
          }
          return std::sqrt(std::real(temp));
     }

     void Fcwf::normalize(){
          MADNESS_ASSERT(m_initialized);
          //print("normalize1");
          double norm = norm2();
          //print("normalize2", norm);
          MADNESS_ASSERT(norm != 0.0);
          for(int i = 0 ; i < 4 ; i++){
               //print("normalize3, ", i);
               m_psi[i] = (1.0/norm)*m_psi[i];
          }
     }

     Fcwf Fcwf::operator*(madness::complex_function_3d phi){
          MADNESS_ASSERT(m_initialized);
          std::vector<complex_function_3d> temp(4);
          for(int i = 0 ; i < 4 ; i++){
               temp[i] = phi*m_psi[i];
          }
          return Fcwf(temp);
     }

     Fcwf Fcwf::operator*(madness::real_function_3d phi){
          MADNESS_ASSERT(m_initialized);
          std::vector<complex_function_3d> temp(4);
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


std::complex<double> inner(const Fcwf& psi, const Fcwf& phi){
     //print("inner 1");
     std::complex<double> result(0,0);
     for(int i = 0 ; i < 4 ; i++){
          //print("inner 2, ",i);
          result += madness::inner(psi[i],phi[i]);
     }
     //print("inner 3");
     return result;
}

Fcwf apply(real_convolution_3d& op, const Fcwf& psi){
     std::vector<complex_function_3d> temp;
     for(int i = 0 ; i < 4 ; i++){
          temp.push_back(madness::apply(op, psi[i]));
     }
     return Fcwf(temp);
}

real_function_3d squaremod(Fcwf psi){
     MADNESS_ASSERT(psi.getinitialize());
     real_function_3d temp = abssq(psi[0]) + abssq(psi[1]) + abssq(psi[2]) + abssq(psi[3]);
     return temp;
}

real_function_3d squaremod_small(Fcwf psi){
     MADNESS_ASSERT(psi.getinitialize());
     real_function_3d temp = abssq(psi[2]) + abssq(psi[3]);
     return temp;
}

real_function_3d squaremod_large(Fcwf psi){
     MADNESS_ASSERT(psi.getinitialize());
     real_function_3d temp = abssq(psi[0]) + abssq(psi[1]);
     return temp;
}

complex_function_3d inner_func(Fcwf psi, Fcwf phi){
     MADNESS_ASSERT(psi.getinitialize() && phi.getinitialize());
     complex_function_3d result = conj(psi[0])*phi[0];
     result += conj(psi[1])*phi[1];
     result += conj(psi[2])*phi[2];
     result += conj(psi[3])*phi[3];
     return result;
}

Fcwf copy(Fcwf psi){
     return Fcwf(psi);

}

std::complex<double> inner(std::vector<Fcwf> a, std::vector<Fcwf> b){
     //print("inner 1");
     MADNESS_ASSERT(a.size() == b.size());
     std::complex<double> result(0,0);
     for(int i = 0; i < a.size(); i++){
          //print("inner 2, i =", i);
          result += inner(a[i],b[i]);    
     }
     //print("inner 3");
     return result;
}

std::vector<Fcwf> operator*(std::vector<Fcwf> psis, std::complex<double> a){
     std::vector<Fcwf> result;
     //print("multiply 1");
     if(psis.size() != 0){
          for(int i = 0; i < psis.size(); i++){
               //print("multiply 2, ", i);
               result.push_back(psis[i]*a);
          }
     }
     //print("multiply 3");
     return result;
}

std::vector<Fcwf> operator*(std::complex<double> a, std::vector<Fcwf> psis){
     std::vector<Fcwf> result;
     //print("multiply 1");
     if(psis.size() != 0){
          for(int i = 0; i < psis.size(); i++){
               //print("multiply 2, ", i);
               result.push_back(psis[i]*a);
          }
     }
     //print("multiply 3");
     return result;
}

void operator+=(std::vector<Fcwf>& phi, std::vector<Fcwf> psi){
     std::vector<Fcwf> result;
     //print("inplace add 1");
     if(phi.size()==0){
          //print("inplace add 2");
          phi = psi;
     }
     else if(psi.size() != 0){
          //print("inplace add 3");
          MADNESS_ASSERT(phi.size()==psi.size());
          for(int i=0; i < psi.size(); i++){
               //print("inplace add 4, i =", i);
               phi[i]+=psi[i];
          }
     }
     //print("inplace add 5");
     return;
}

std::vector<Fcwf> operator-(std::vector<Fcwf> phi, std::vector<Fcwf> psi){
     std::vector<Fcwf> result;
     //print("subtract 1");
     if(phi.size()==0){
          //print("subtract 2");
          std::vector<Fcwf> temp = -1.0*psi;
          return temp;
     }
     else if(psi.size()==0){
          //print("subtract 3");
          return phi;
     }
     else{
          //print("subtract 4");
          MADNESS_ASSERT(phi.size()==psi.size());
          for(int i=0; i < psi.size(); i++){
               //print("subtract 5, i =", i );
               result.push_back(phi[i]-psi[i]);
          }
     }
     //print("subtract 6");
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


//kthxbye
