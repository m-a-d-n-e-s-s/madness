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

void readTensor(std::ifstream &fs, Tensor<double> &T, int dim, int nbf, bool sym);
void computeDensity(Tensor<double> &P, Tensor<double> &Cocc, Tensor<double> &C, int Nocc);
Tensor<double> computeG(Tensor<double> &twoE, Tensor<double> P, int nbf);
double computeTwoEE(Tensor<double> twoE, Tensor<double> P, int nbf);
double computeEtot(Tensor<double> P, Tensor<double> H, Tensor<double> F, int nbf);
double computeMP2(Tensor<double> twoEE, Tensor<double> epsilons, int nocc, int K);

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
    fs >> keyword; //read the first key word
    if (keyword == "nbf")
    {
        fs >> nbf;
    }
    else
    {
        throw "input file should start with nbf";
    }

    Tensor<double> S(nbf, nbf);
    Tensor<double> KE(nbf, nbf);
    Tensor<double> PE(nbf, nbf);
    Tensor<double> MUX(nbf, nbf);
    Tensor<double> MUY(nbf, nbf);
    Tensor<double> MUZ(nbf, nbf);
    Tensor<double> MOS(nbf, nbf);
    Tensor<double> Electron(nbf, nbf, nbf, nbf);

    while (fs >> keyword)
    {

        if (keyword == "nocc")
        {
            fs >> nocc;
        }
        else if (keyword == "enrep")
        {
            fs >> enrep;
        }
        else if (keyword == "overlap")
        {
            readTensor(fs, S, 2, nbf, true);
        }
        else if (keyword == "ke")
        {
            readTensor(fs, KE, 2, nbf, true);
        }
        else if (keyword == "pe")
        {
            readTensor(fs, PE, 2, nbf, true);
        }
        else if (keyword == "mux")
        {
            readTensor(fs, MUX, 2, nbf, true);
        }
        else if (keyword == "muy")
        {
            readTensor(fs, MUY, 2, nbf, true);
        }
        else if (keyword == "muz")
        {
            readTensor(fs, MUZ, 2, nbf, true);
        }
        else if (keyword == "mos")
        {
            readTensor(fs, MOS, 2, nbf, false);
        }
        else if (keyword == "2-electron")
        {
            print("entering");
            readTensor(fs, Electron, 4, nbf, true);
        }
        else
        {
            print("Unknown key word", keyword);
            throw "bad";
        }
    }

    /* The first stage of the calculation after having all the integrals define is to diagonalize the 
** the overlap matrix S **/

    Tensor<double> s(nbf, 1);   //holds the eigenvalues of S
    Tensor<double> U(nbf, nbf); //holds the unitary matrix that transforms S
    syev(S, U, s);              // solve SU=sU
    Tensor<double> s12(nbf, nbf);

    // Create s^(-1/2) from eigenvalues in s
    for (int mu = 0; mu < nbf; mu++)
    {
        s12(mu, mu) = 1 / pow(s(mu), .5);
    }
    // Create X canoncial orthogonalization of S
    /*  X_{mu mu} = sum(mu') Umu mu' * s12 mu' nu' */
    Tensor<double> X = inner(U, s12);
    // Compute the Core Hamiltonian with KE and PE
    Tensor<double> Hcore = KE + PE;
    // Compute the Intial Guess to the Core to Density Matrix
    Tensor<double> C(nbf, nbf);
    Tensor<double> Cocc = C(_, Slice(0, nocc - 1));
    Tensor<double> P = 2 * inner(Cocc, Cocc, 1, 1);

    // Convergence criteria Standard deviation of successive density matrix elements
    int iter = 1;
    double Etot = 100001;
    double deltaE = 10;
    double E0(0);
    // Fock Matrix
    Tensor<double> F(nbf, nbf);
    F = Hcore;
    Tensor<double> G(nbf, nbf);
    Tensor<double> Fprime(nbf, nbf);
    Tensor<double> Cprime(nbf, nbf);
    Tensor<double> epsilons(nbf, 1);
    Tensor<double> CprimeOcc(nbf, nocc);
    double twoEE(0);
    double E1(0);

    double del = 1e-10;
    int maxIter = 100;
    while (std::abs(deltaE) > del && iter < maxIter)
    {

        E0 = Etot;

        // Fprime = transform(F, X);
        // syev(Fprime, Cprime, epsilons); // solve SU=sU
        sygv(F, S, 1, C, epsilons);
        // CprimeOcc = Cprime(_, Slice(0, nocc - 1));
        Cocc = C(_, Slice(0, nocc - 1));
        // Cocc=inner(X,CprimeOcc);
        P = 2 * inner(Cocc, Cocc, 1, 1);
        // compute two electron integral portion of Fock matrix
        G = computeG(Electron, P, nbf);
        F = Hcore + G;

        E1 = Hcore.trace(P);
        twoEE = computeTwoEE(Electron, P, nbf);
        Etot = E1 + twoEE + enrep;

        std::cout << "One-electron energy =    " << E1 << std::endl;
        std::cout << "Two-electron energy =   " << twoEE << std::endl;
        std::cout << "Total SCF energy  =     " << Etot << std::endl;

        std::cout << " Iteration " << iter << std::endl;
        deltaE = std::abs((Etot - E0));
        std::cout << "Detla E =" << deltaE << std::endl;
        std::cout << std::endl;

        iter++;
    }

    // calculation of dipole moments
    Tensor<double> dipole(3, 1);

    Tensor<double> dipole2(3, 1);

    for (int mu = 0; mu < nbf; mu++)
    {
        for (int nu = 0; nu < nbf; nu++)
        {
            dipole(0, 0) += P(mu, nu) * MUX(nu, mu);
            dipole(1, 0) += P(mu, nu) * MUY(nu, mu);
            dipole(2, 0) += P(mu, nu) * MUZ(nu, mu);
        }
    }

    dipole2(0, 0) = P.trace(MUX);
    dipole2(1, 0) = P.trace(MUY);
    dipole2(2, 0) = P.trace(MUZ);

    print(dipole);
    print("Electric Dipole");
    print(dipole2);


    //MP2 Correlation Energy
    double EMP2 = computeMP2(Electron,epsilons,nocc,nbf);
    print("MP2 Correlation Energy");
    print(EMP2);
    print("Total MP2 Energy " ,Etot+EMP2);
    print(epsilons);
    // ********************************************************
    // Everything below is the has the answers
    //This is for Hideo's education
    bool answer = true;
    if (answer == true)
    {
        Tensor<double> Cocc = MOS(_, Slice(0, nocc - 1));
        Tensor<double> D = 2 * inner(Cocc, Cocc, 1, 1);
        double E1 = Hcore.trace(D);
        double twoEE(0);
        twoEE = computeTwoEE(Electron, D, nbf);
        Etot = E1 + twoEE + enrep;
        //
        std::cout << "Answers **************************" << std::endl;
        ;
        std::cout << "One-electron energy =    " << E1 << std::endl;
        std::cout << "Two-electron energy =   " << twoEE << std::endl;
        std::cout << "Total SCF energy  =     " << Etot << std::endl;
    }

    fs.close();
    return 0;
}

void readTensor(std::ifstream &fs, Tensor<double> &T, int dim, int nbf, bool sym)
{
    int count = 1;
    int a, b, c, d;
    int *indices;
    double val;             // value
    indices = new int[dim]; // dynamic array of indicies
    std::string line;
    do
    {
        std::getline(fs, line);
        std::istringstream iss(line);
        // gather indices

        for (int i = 0; i < dim; i++)
        {
            iss >> indices[i];
        }
        iss >> val;
        if (dim == 4 && indices[0] == 1)
        {
        }

        if (indices[0] != -1 &&
            indices[0] != 0)
        { // make sure it passes a blank line
            if (dim == 2)
            {
                a = indices[0] - 1;
                b = indices[1] - 1;
                T(a, b) = val;

                if (sym)
                {
                    T(b, a) = val;
                }
            }
            else if (dim == 4)
            {
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

                if (count == 1)
                    //                   print("huh", a, b, c, d, indices[0]);
                    count++;
            }
            else
            {
                throw "bad";
            }
        }
    } while (indices[0] != -1);

    delete[] indices; // releases the memory allocated for arrays of elements
                      // using new and a size in brackets([])
}

void computeDensity(Tensor<double> &P, Tensor<double> &Cocc, Tensor<double> &C, int Nocc)
{
    Cocc = C(_, Slice(0, Nocc - 1));
    P = 2 * inner(Cocc, Cocc, 1, 1);
}
Tensor<double> computeG(Tensor<double> &twoE, Tensor<double> P, int nbf)
{
    Tensor<double> G(nbf, nbf);
    double Gmn = 0;
    for (int mu = 0; mu < nbf; mu++)
    {
        for (int nu = 0; nu < nbf; nu++)
        {
            for (int lambda = 0; lambda < nbf; lambda++)
            {
                for (int sigma = 0; sigma < nbf; sigma++)
                {
                    Gmn = P(lambda, sigma) * (twoE(mu, nu, lambda, sigma) - 0.5 * twoE(mu, lambda, sigma, nu));
                    G(mu, nu) = G(mu, nu) + Gmn;
                }
            }
        }
    }
    return G;
}
double computeTwoEE(Tensor<double> twoE, Tensor<double> P, int nbf)
{
    double twoEE(0);
    for (int mu = 0; mu < nbf; mu++)
    {
        for (int nu = 0; nu < nbf; nu++)
        {
            for (int lambda = 0; lambda < nbf; lambda++)
            {
                for (int sigma = 0; sigma < nbf; sigma++)
                {
                    twoEE += twoE(mu, nu, lambda, sigma) * (0.5 * P(mu, nu) * P(lambda, sigma) - .25 * P(mu, lambda) * P(nu, sigma));
                }
            }
        }
    }
    return twoEE;
}
double computeEtot(Tensor<double> P, Tensor<double> H, Tensor<double> F, int nbf)
{
    double E0(0);
    for (int mu = 0; mu < nbf; mu++)
    {
        for (int nu = 0; nu < nbf; nu++)
        {
            E0 += P(nu, mu) * (H(mu, nu) + F(mu, nu));
        }
    }
    return E0 / 2;
}

double computeMP2(Tensor<double> twoEE, Tensor<double> epsilons, int nocc, int K)
{
    print(K);
    print(nocc);
    K=nocc*2;
    double E = 0;
    for (int a = 0; a < nocc; a++)
    {

        for (int b = 0; b < nocc; b++)
        {

            for (int r = nocc; r < K; r++)
            {

                for (int s = nocc; s < K; s++)
                {
                    E+=2*twoEE(a,r,b,s)*twoEE(r,a,s,b)/(epsilons(a)+epsilons(b)-epsilons(r)-epsilons(s));
                    E-=twoEE(a,r,b,s)*twoEE(r,b,s,a)/(epsilons(a)+epsilons(b)-epsilons(r)-epsilons(s));
                }
            }
        }
    }
    return E;
}

    // 1. Read in the input file "integrals.dat" Save the data to corresponding variables
    // I need to figure out what every variable represents.
