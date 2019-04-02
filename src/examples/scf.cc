/** Program to solve SCF Equations and Compute SCF and MP2 energies
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

using namespace madness;

void readTensor(std::ifstream& fs, Tensor<double>& T,int dim, int nbf);

int main()
{
    std::ifstream fs;
    fs.open("integrals.dat");
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
    int nbf;
    int nocc;
    double enrep;
    /*************** Read and printout data from integrals.dat *********
     *
     * Current task: Write a function that reads array data
     * Note:First Thing I read in has to be nbf that way I can I can create
     * Tensors with nbf sizes
     ***********************************************************************/
    std::string keyword;
    // I can define Tensors now i think
    Tensor<double> S(nbf,nbf);
    Tensor<double> KE(nbf,nbf);
    Tensor<double> PE(nbf,nbf);
    Tensor<double> MUX(nbf,nbf);
    Tensor<double> MUY(nbf,nbf);
    Tensor<double> MUZ(nbf,nbf);
    Tensor<double> Electron(nbf,nbf,nbf,nbf);


    while (fs >> keyword) {

        if (keyword =="nocc") {
            fs >> nocc;
        }
        else if ( keyword == "enrep") {
            fs >> enrep;
        }
        else if ( keyword == "overlap") {
            readTensor(fs,S,2,nbf);
        }
        else {
            print("Unknown key word",keyword);
            //throw "bad";
        }
    }

    std::cout << nbf <<"this is nbf"<< std::endl;
    fs.close();
    /**********************************************************************************
     *
     * ********************************************************************************/
    Tensor<double> f(2, 2);
    print(f);
    return 0;
}

void readTensor(std::ifstream& fs, Tensor<double>& T,int dim, int nbf) {

    int i;
    int j;
    
    std::string token;
    while (fs >>token ){
      
    }
     

}

// 1. Read in the input file "integrals.dat" Save the data to corresponding variables
// I need to figure out what every variable represents.
