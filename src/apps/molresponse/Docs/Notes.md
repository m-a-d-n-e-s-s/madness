# Notes on response solver

Here we will get an overview of how the solve function from `void TDHF::solve(World & world)`

## The main iteration

We iterate for a number set by by the size of the vector `r_params.protocol_data.size();`
So it looks like we iterate with a set of threshold values for truncation.
This way we can increase the accuracy progressively per iteration.

`std::vector<double> protocol_data; ///< Different thresholds for truncation`

For each iteration, we check that each function has the same k value... The k value sets the number of polynomials per box

Q? Why would this not be the case?
Q? Any how would this change between iterations?

### Check_k

Returns a bool to redo the groundstate hamiltonian calculation if groundstate orbitals change

If the default wavelet order does not equal to the wavelet order of the GS orbitals then we re-read the orbitals from the archive.  We then reconstruct the orbitals with `reconstruct(world, ground_orbitals);`

`reconstruct(world, ground_orbitals);`

``` reconstruct a vector of functions
set fence to false
for f in vF:
    If f is compressed:


## Calculation Parameters

call initialize with `key value, comment, and vector of allowed values`
