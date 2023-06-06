# Chemistry API

There are several interfaces and methods for writing quantum chemistry codes and we will write a simple HF code as an example.

## SCFOperators
Operators occuring in SCF methods are implemented in SCFOperators.h. Currently there 
are: 
 - ```Kinetic```
 - ```Laplacian```
 - ```Derivative```
 - ```Coulomb```
 - ```Exchange```
 - ```XCOperator``` (using LibXC)
 - ```Nuclear```
 - ```DNuclear``` (derivative of the nuclear potential)
 - ```Lz``` (angular momentum)
 - ```Local``` (any local potential)
 
All of the above can be put into a ```Fock``` operator, which is simply a linear combinaton 
of ```SCFOperator```, to simplify notation. Construction follows the scheme
```c++
auto T=Kinetic<double,3>(world);
auto Vnuc=Nuclear<double,3>(world,molecule);
```
SCFOperators act on MRA functions or vector of functions, yielding either the a result (vector of) functions, 
or a matrix representation of the operator
```c++
std::vector<Function<double,3>> Vphi=Vnuc(orbitals);
Tensor<double> tmat=T(orbitals,orbitals);
```

## BSHApply
Solving the Schroedinger equation in MRA always looks like
```math
\displaylines{
(T + V)\psi  = E \psi \\
\psi = (T - E)^{-1}V\psi
}
```
If there are more than one wave functions (e.g. several orbitals) there might appear
a coupling term 
```math
\displaylines{F\psi_i  = f_{ij} \psi_j \\
\psi_i = (T - f_{ii})^{-1}V\psi_i + \sum_{j\neq i}\psi_j}
```
Since this operator occurs very often, and might include further techniques, 
such as a level shift if the energy is positive, a class named ```BSHApply``` 
has been introduced that applies the BSHOperator, includes couping and performs 
a level shift, if necessary
```c++
auto [residual, eps_update] = bsh_apply(orbitals, fock, Vorbitals);
```
New orbitals are the taken directly from the residual, or by using a convergence accelerator
```c++
orbitals-=residual;
```

## KAIN Solver
Similar to DIIS, the KAIN solver stabilizes and accelerates the solution of the nonlinear equations.

```c++
auto functionsolver=NonlinearSolver;
auto vectorfunctionsolver=nonlinear_vector_solver<double,3>(world, orbitals.size());
```
The new solution is a linear combination of the old solutions and the new update and can be obtained through
```c++
orbitals=solver.update(orbitals,residual);
```

## Hartree-Fock example codes

* [Simple HF code for helium in C++](https://github.com/m-a-d-n-e-s-s/madness/blob/master/src/examples/hehf.cc) and associated [documentation](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__examplehehf.html)
* [Simple molecular HF code using the chemistry API](https://github.com/m-a-d-n-e-s-s/madness/blob/master/doc/tutorial/simple_hf.cpp)


