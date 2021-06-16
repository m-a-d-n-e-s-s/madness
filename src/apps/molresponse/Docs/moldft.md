# Notes on moldft

# Solve function in SCF.cc

## Objective: For given protocol, solve the DFT/HF/response equations

### Pre iteration

Before iteration we create the necessary variables

  - `functionT arho_old, brho_old;`
  - `const double dconv = std::max(FunctionDefaults<3>::get_thresh(), param.dconv());`
  - `const double trantol = vtol`
  - `const double tolloc = 1e-6;`

### Per iteration

For each iteration

#### localize orbitals and compute densities

- localize orbitals based on some localization method
- compute the density
	- arho
	- brho
- `if (iter < 2|| iter %10) loadbalance`
- compute da and db. density differences
- compute total density rho=arho+brho

#### compute local potentials vlocal

- vcoul=J[rho]
	- ecoul=.5 *inner(rho,vcoul)
- vnuc
- vlocal = vcoul + vnuc
- if pcm add vpcm

#### compute Vpsia and Vpsib

we pass in orbitals and vlocal

- if dft
	- create XCOperator
	- compute exc
	- add :vlocal += vxc// `vloc += xcoperator.make_xc_potential();`

- Create Vpsi
	- vlocal * orbitals

- If HF exchange
	- Create Exchangeoperator
	- compute exchf
	- if not spin_polarized
		- exchf *=2

- Add Kamo to Vpsi using gaxpy
- Vpsi***=Vpsi-Kamo

#### make fock matrix focka and fockb

#### Diagonalization of Fock matrix

#### Compute energies

- enrep
- ekinetic
- enonlocal
- exc
- etot

print

#### Before iteration end

- if iter>0
- if converge or last iter
	- diagonlize for eigen values
	- localize one last time
	- break
```cpp
			if (da < dconv * std::max(size_t(5), molecule.natom()) &&
          db < dconv * std::max(size_t(5), molecule.natom()) &&
          (param.get<bool>("conv_only_dens") || bsh_residual < 5.0 * dconv))
        converged = true;

```
- if does not converge
#### update_subspace

- zero out off diagon lagrange multipliers
- Compute residual

##### Compute residual

- set eps(i) as the min of -.05 or fock(i,i)
- fock(i,i) -= eps(i)

- transform psi using fock matrix
- fpsi=fock*** * psi

- undo the damage from above

- Vpsi=Vpsi-fpsi
- scale Vpsi
- make bsh operator with eps tensor

- newpsi= apply(world,ops,Vpsi)
- r=sub(psi,new_psi)

## Update subspace

`typedef std::pair<vecfuncT, vecfuncT> pairvecfuncT;`

- world
- Vpsia
- Vpsib
- focka
- fockb
- subspace typedef std::vector<pairvecfuncT> subspaceT;
- Q
- bsh_residual
- update_residual

1. First zero off diagonal lagrange mulitpliers of focka
2. Compute_residual
	- set diagonal to min of 0.05 or fock(i,i)
	- fpsi= psi*fock
	- Vpsi=Vpsi-fpsi
	- scale: 2*Vpsi
	- create bsh_operators
	- set thresh for Vpsi
	- newpsi=apply op to Vpsi
	- clear ops and Vpsi and fence
	- truncate newpsi
	- r=psi-newpsi
	- rnorm=nrom2s(world,r)
	- BSH residual: rms and max
	- r is the residual in wavefunctions
3. If spin beta
	- compute residual
	- insert beta to vm=amo
	- insert residuals
4. BSH residual is max (aerr,berr)
5. broadcast(bsh_residual)
6. compress vm and rm
7. Subspace = pair of vm and residuals





        /// Transforms a vector of functions left[i] = sum[j] right[j]*c[j,i] using sparsity
        /// @param[in] vright vector of functions (impl's) on which to be transformed
        /// @param[in] c the tensor (matrix) transformer
        /// @param[in] vleft vector of of the *newly* transformed functions (impl's)










