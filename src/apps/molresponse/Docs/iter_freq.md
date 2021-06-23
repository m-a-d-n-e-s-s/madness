# Notes on Response state solver


We are solving for the frequency or static response states.
These are a set of functions held in a `X_space` object
called `Chi` .
`Chi` holds both `x` and `y` response states.

We begin with zero_functions for each of the orbitals in x and y response states

## Iteration

- We compute the response density given the current x and y
	- multiply and add orbitals and response states
- Compute the density residual vectors
- Copy newChi and new_rho into old
- Next we move into update_x_space_response

### Update X Space

First we need to compute the X_space vector $\Theta \Chi$



