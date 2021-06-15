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
- if (iter < 2|| iter %10) loadbalance
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












