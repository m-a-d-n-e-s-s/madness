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

    // I can define Tensors now i think
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
            readTensor(fs,S,2,nbf);
        }
        else if ( keyword == "ke") {
            readTensor(fs,KE,2,nbf);
        }
        else if ( keyword == "pe") {
            readTensor(fs,PE,2,nbf);
        }
        else if ( keyword == "mux") {
            readTensor(fs,MUX,2,nbf);
        }
        else if ( keyword == "muy") {
            readTensor(fs,MUY,2,nbf);
        }
        else if ( keyword == "muz") {
            readTensor(fs,MUZ,2,nbf);
        }
        else if ( keyword == "mos") {
            readTensor(fs,MOS,2,nbf);
        }
        else if ( keyword == "2-electron") {
            readTensor(fs,Electron,4,nbf);
        }
        else {
            print("Unknown key word",keyword);
            throw "bad";
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

    int * indices;
    double val;//value
    indices= new int [dim];// dynamic array of indicies
    std::string token;
    std::string line;
    do {
        std::getline(fs,line);
        std::istringstream iss(line);
        // gather indices
        for (int i=0; i<dim; i++) {
            iss>>indices[i];
        }
        //
        iss>>val;
        for (int i=0; i<dim; i++) {
            std::cout<<" "<<indices[i];
        }
        std::cout<<" "<<val<<std::endl;
    } while ( indices[0]!=-1);

    delete[] indices;//releases the memory allocated for arrays of elements using new and a size in brackets([])
}

// 1. Read in the input file "integrals.dat" Save the data to corresponding variables
// I need to figure out what every variable represents.
