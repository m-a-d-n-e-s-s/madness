/* Program to solve SCF Equations and Compute SCF and MP2 energies
 *  Program outline
 *  1. Read and print out all data from integral.dat
 *  2. Develop some notes in latex that details the equations
 *  3. Use integrals and given MO coefficients to compute SCF energy/ compare with NWCHEM
 *
 *
 *
 *
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/print.h>
#include <cmath>

using namespace madness;

void readTensor(std::ifstream& fs, Tensor<double>& T,int dim, int nbf,bool sym);

int main()
{
    /********* Variables **********************************************************
     * nbf: Number of Basis Functions
     * nocc: ?
     * enrep: effective nuclear repulsion
     * overlap: OverlapIntegrals | Tensor<double> S(nbf,nbf)
     * ke: Kinetic Energy Operator | Tensor<double> KE(nbf,nbf)
     * pe: Potential Energy Operator Integrals | Tensor<double> PE(nbf,nbf)
     * mux: DipoleX | Tensor<double> mux(nbf,nbf)
     * muy: Dipoley | Tensor<double> muy(nbf,nbf)
     * muz: Dipolez | Tensor<double> muz(nbf,nbf)
     * mos: ???? | Tensor<double> mos(nbf,nbf)
     * 2-electron | Tensor<double> C(nbf,nbf,nbf,nbf)
     *****************************************************************************/
    // make ifstream object and read in the first key word and read in the nbf
    std::ifstream fs;
    fs.open("integrals.dat");
    int nbf;
    int nocc;
    double enrep;
    std::string keyword;
    fs>>keyword;//read the first key word
   if (keyword == "nbf") {
        fs>>nbf;
    }
    else {
        throw "input file should start with nbf";
    }

    Tensor<double> S(nbf,nbf);
    Tensor<double> KE(nbf,nbf);
    Tensor<double> PE(nbf,nbf);
    Tensor<double> MUX(nbf,nbf);
    Tensor<double> MUY(nbf,nbf);
    Tensor<double> MUZ(nbf,nbf);
    Tensor<double> MOS(nbf,nbf);
    Tensor<double> Electron(nbf,nbf,nbf,nbf);


    while (fs >> keyword) {

        if (keyword =="nocc") {
            fs >> nocc;
        }
        else if ( keyword == "enrep") {
            fs >> enrep;
        }
        else if ( keyword == "overlap") {
            readTensor(fs,S,2,nbf,true);
        }
        else if ( keyword == "ke") {
            readTensor(fs,KE,2,nbf,true);
        }
        else if ( keyword == "pe") {
            readTensor(fs,PE,2,nbf,true);
        }
        else if ( keyword == "mux") {
            readTensor(fs,MUX,2,nbf,true);
        }
        else if ( keyword == "muy") {
            readTensor(fs,MUY,2,nbf,true);
        }
        else if ( keyword == "muz") {
            readTensor(fs,MUZ,2,nbf,true);
        }
        else if ( keyword == "mos") {
            readTensor(fs,MOS,2,nbf,false);
        }
        else if ( keyword == "2-electron") {
            print("entering");
            readTensor(fs,Electron,4,nbf,true);
        }
        else {
            print("Unknown key word",keyword);
            throw "bad";
        }
    }

/* The first stage of the calculation after having all the integrals define is to diagonalize the 
** the overlap matrix S **/

    Tensor<double> s(nbf,1);//holds the eigenvalues of S
    Tensor<double> U(nbf,nbf);//holds the unitary matrix that transforms S
    syev(S,U,s);// solve SU=sU
    Tensor<double> s12(nbf,nbf);

// Create s^(-1/2) from eigenvalues in s
    for(int mu = 0;mu<nbf;mu++){
        s12(mu,mu)=1/pow(s(mu),.5);
    }
// Create X canoncial orthogonalization of S
/*  X_{mu mu} = sum(mu') Umu mu' * s12 mu' nu' */
    Tensor<double> X = inner(U,s12);
// Compute the Core Hamiltonian with KE and PE 
    Tensor<double> Hcore=KE+PE;




    
// ******************************************************** 
// Everything below is the has the answers
bool answer =false
if(answer == true){
    Tensor<double> Cocc=MOS(_,Slice(0,nocc-1));
    Tensor<double> D=inner(Cocc,Cocc,1,1);

    double pe=PE.trace(D);

    double ke=KE.trace(D);
    print(ke+pe);
    double E1=Hcore.trace(D);
    print(E1);
    double twoEE(0);

    for( int mu =0; mu <nbf; mu++) {
        for (int nu =0; nu <nbf; nu++) {
            for( int lambda =0; lambda <nbf; lambda ++) {
                for ( int sigma=0; sigma <nbf; sigma ++) {
                    twoEE+=Electron(mu,nu,lambda,sigma)*(2*D(mu,nu)*D(lambda,sigma)-D(mu,lambda)*D(nu,sigma));
                }
            }
        }
    }

}


    fs.close();
    return 0;
}

void readTensor(std::ifstream& fs, Tensor<double>& T, int dim, int nbf,
                bool sym) {
  int count = 1;
  int a, b, c, d;
  int* indices;
  double val;              // value
  indices = new int[dim];  // dynamic array of indicies
  std::string line;
  do {
    std::getline(fs, line);
    std::istringstream iss(line);
    // gather indices

    for (int i = 0; i < dim; i++) {
      iss >> indices[i];
    }
    iss >> val;
    if (dim == 4 && indices[0] == 1) {
      print("before if(indices[0]!=-1");
      for (int i = 0; i < dim; i++) {
        std::cout << " " << indices[i];
      }
      std::cout << " " << val << std::endl;
    }

    if (indices[0] != -1 &&
        indices[0] != 0) {  // make sure it passes a blank line
      if (dim == 2) {
        a = indices[0] - 1;
        b = indices[1] - 1;
        T(a, b) = val;

        if (sym) {
          T(b, a) = val;
        }
      } else if (dim == 4) {
        a = indices[0] - 1;
        b = indices[1] - 1;
        c = indices[2] - 1;
        d = indices[3] - 1;

        T(a, b, c, d) = val;
        T(a, b, d, c) = val;
        T(b, a, c, d) = val;
        T(b, a, d, c) = val;
        T(c, d, a, b) = val;
        T(c, d, b, a) = val;
        T(d, c, a, b) = val;
        T(d, c, b, a) = val;

        if (count == 1) print("huh", a, b, c, d, indices[0]);
        count++;

      } else {
        throw "bad";
      }
    }
  } while (indices[0] != -1);

  delete[] indices;  // releases the memory allocated for arrays of elements
                     // using new and a size in brackets([])
}

// 1. Read in the input file "integrals.dat" Save the data to corresponding variables
// I need to figure out what every variable represents.