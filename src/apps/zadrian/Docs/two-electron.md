# Electron Interaction in Madness Response

## Summary TL;DR

In our response calculations we need to compute and the conjugate,

$$
\hat{\Gamma}^{0,b}=\frac{\partial\hat{g}}{{\partial\rho}}\left[\rho^{0}\right]\star\rho^{(b)}(r,r')

$$

in order it iterate. Considering the **J,K,W** terms we simplify to the final expressions which are implemnted in the code.

### Working Equations

$$
\hat{\Gamma}^{0,b}=H_{x,p}^{(b)}(r)+G_{y,p}^{(b)}(r)

$$
$$
\hat{\Gamma}^{\dagger0,b}=H_{y,p}^{(b)}(r)+G_{x,p}^{(b)}(r)

$$
$$
H_{f,p}^{(k)}(r)=\int dr'\frac{f^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)-\gamma\sum_{i}f_{i}^{(b)}(r)\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)+W_{f,p}^{(k)}(r)

$$
$$
G_{f,p}^{(k)}(r)=\int dr'\frac{f^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)-\gamma\sum_{i}\varphi_{i}(r)\left(\int dr'\frac{y_{i}^{\dagger(b)}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)+\hat{W_{y}}^{0,b}\varphi_{k}(r)

$$

### Implementation

Faster implementation

Every iteration we need to compute these electron interaction terms

$$
\hat{\Gamma}^{0,b}\varphi_{k}(r)=\int dr'\frac{\rho^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)-\gamma\sum_{i}x_{i}^{(b)}(r)\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)-\gamma\sum_{i}\varphi_{i}(r)\left(\int dr'\frac{y_{i}^{\dagger(b)}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)+\hat{W}^{0,b}\varphi_{k}(r)

$$
$$
\hat{\Gamma}^{\dagger0,b}\varphi_{k}(r)=\int dr'\frac{\rho^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)-\gamma\sum_{i}\varphi_{i}(r)\left(\int dr'\frac{x_{i}^{\dagger(b)}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)-\gamma\sum_{i}y_{i}^{(b)}(r)\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)+\hat{W}^{\dagger0,b}\varphi_{k}(r)

$$

1.  Before entering the iteration we compute

$$
\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)

$$

2.  At the start of each iteration we also compute the vector of density resopnses

$$
\rho_{H}^{b_{N}}(x,x')=\sum_{i}\left(x_{i}^{(b_{N})}(x)\varphi_{i}^{\dagger}(x')+\varphi_{i}(x)y_{i}^{\dagger(b_{N})}(x')\right)

$$

3.  During each iteration for each density response we compute the n coloumb terms

$$
\int dr'\frac{\rho^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)

$$

4.  Then for each response density we compute the n exchange terms.

* * *

## Introduction

In our response calculations we need to compute,

$$
\hat{\Gamma}^{0,b}=\frac{\partial\hat{g}}{{\partial\rho}}\left[\rho^{0}\right]\star\rho^{(b)}(r,r')

$$

Here $\hat{g}$ is the electron interaction operator,

$$
\hat{g}=\hat{J}-\gamma\hat{K}+\hat{W}

$$

where $\gamma$ is the scaling factor for the exchange energy.

- $\gamma=0$ is pure DFT
- $\gamma$ is non zero is a hybrid theory
- $\gamma=1$ and remomving exchange-correlation is HF

### Coulomb Derivation

In the coulomb case we start with,

$$
\hat{j}(r')f(r')=\int dx_{1}\frac{\rho(r,r)}{\left|r-r'\right|}f(r')

$$

which can be rewritten as,

$$
\hat{j}[\rho](r)f(r')=\int drds\frac{\rho(r,s)\delta(r-s)}{\left|r-r'\right|}f(r')

$$

taking the functional derivative,

$$
\frac{\partial\hat{J}}{{\partial\rho}}\left[\rho^{0}\right]f(r')=\frac{\delta(r-s)}{\left|r-r'\right|}f(r')

$$

and convolving with the density response derivative,

$$
\frac{\partial\hat{J}}{{\partial\rho}}\left[\rho^{0}\right]\star\rho^{(b)}(r,r')f(r')=\int drds\frac{\delta(r-s)\rho^{(b)}(r,s)}{\left|r-r'\right|}f(r')

$$
$$
\frac{\partial\hat{J}}{\partial\rho}\left[\rho^{0}\right]\star\rho^{(b)}(r,r')f(r')=\int dr\frac{\rho(r,r)}{\left|r-r'\right|}f(r')

$$

### Exchange Derivation

The echange operator is given as,

$$
\hat{k}(r')f(r')=-\int dr\frac{\rho(r,r')}{\left|r-r'\right|}f(r)

$$

which can be rewritten as,

$$
\hat{k}(r')f(r')=-\int\int dsdr\frac{\rho(r,s)\delta(r'-s)}{\left|r-r'\right|}f(r)

$$

Taking the functional derivative,

$$
\frac{\partial\hat{K}}{\partial\rho}\left[\rho^{0}\right]f(r')=-\frac{\delta(r'-s)}{\left|r-r'\right|}f(r')

$$

We convolve this term with the density response,

$$
\frac{\partial\hat{k}}{\partial\rho}\left[\rho^{0}\right]\star\rho^{(b)}(r,r')f(r')=-\int\int dsdr\frac{\delta(r'-s)\rho^{(b)}(r,s)}{\left|r-r'\right|}f(r)

$$

Which leaves us with the final expression:

$$
\frac{\partial\hat{k}}{\partial\rho}\left[\rho^{0}\right]\star\rho^{(b)}(r,r')f(r')=-\int dr\frac{\rho^{(b)}(r,r')}{\left|r-r'\right|}f(r)

$$

* * *

## Deriving final equations

Given the electron density is defined as

$$
\rho_{H}^{b_{N}}(r,r')=\sum_{i}\left(x_{i}^{(b_{N})}(r)\varphi_{i}^{\dagger}(r')+\varphi_{i}(r)y_{i}^{\dagger(b_{N})}(r')\right)

$$

we can write out the expressions for each electron interaction

### Coulomb Terms

Given,

$$
\hat{\Gamma}^{0,b}\varphi_{k}(r)=\int dr'\frac{\rho^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)

$$

### Exchange Terms

$$
\hat{\Gamma}^{0,b}\varphi_{k}(r)=-\int dr\frac{\rho^{(b)}(r,r')}{\left|r-r'\right|}\varphi_{k}(r')

$$

Expanding out,

$$
\hat{\Gamma}^{0,b}\varphi_{k}(r)=-\sum_{i}x_{i}^{(b)}(r)\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)-\sum_{i}\varphi_{i}(r)\left(\int dr'\frac{y_{i}^{\dagger(b)}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)

$$

and the conjugate terms,

$$
\hat{\Gamma}^{\dagger0,b}\varphi_{k}(r)=-\sum_{i}\varphi_{i}(r)\left(\int dr'\frac{x_{i}^{\dagger(b)}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)-\sum_{i}y_{i}^{(b)}(r)\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)

$$

### All together

Putting eveything together,

$$
\hat{\Gamma}^{0,b}\varphi_{k}(r)=\int dr'\frac{\rho^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)-\gamma\sum_{i}x_{i}^{(b)}(r)\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)-\gamma\sum_{i}\varphi_{i}(r)\left(\int dr'\frac{y_{i}^{\dagger(b)}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)+\hat{W}^{0,b}\varphi_{k}(r)

$$
$$
\hat{\Gamma}^{\dagger0,b}\varphi_{k}(r)=\int dr'\frac{\rho^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)-\gamma\sum_{i}\varphi_{i}(r)\left(\int dr'\frac{x_{i}^{\dagger(b)}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)-\gamma\sum_{i}y_{i}^{(b)}(r)\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)+\hat{W}^{\dagger0,b}\varphi_{k}(r)

$$

##### (working equations) Hf and Gf

It is useful to seperato **x** and **y** terms,

Here we simplify notation,

We re-write,

$$
\hat{\Gamma}^{0,b}\varphi_{k}(r)=\int dr'\frac{x^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)+\int dr'\frac{y^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)-\gamma\sum_{i}x_{i}^{(b)}(r)\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)-\gamma\sum_{i}\varphi_{i}(r)\left(\int dr'\frac{y_{i}^{\dagger(b)}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)+\hat{W_{x}}^{0,b}\varphi_{k}(r)+\hat{W_{y}}^{0,b}\varphi_{k}(r)

$$
$$
\hat{\Gamma}^{\dagger0,b}\varphi_{k}(r)=\int dr'\frac{x^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)+\int dr'\frac{x^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)-\gamma\sum_{i}\varphi_{i}(r)\left(\int dr'\frac{x_{i}^{\dagger(b)}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)-\gamma\sum_{i}y_{i}^{(b)}(r)\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)+\hat{W_{x}}^{\dagger0,b}\varphi_{k}(r)+\hat{W_{y}}^{\dagger0,b}\varphi_{k}(r)

$$

We can use functions,

$$
\hat{\Gamma}^{0,b}=H_{x,p}^{(b)}(r)+G_{y,p}^{(b)}(r)

$$
$$
\hat{\Gamma}^{\dagger0,b}=H_{y,p}^{(b)}(r)+G_{x,p}^{(b)}(r)

$$

in order to seperate the electron interaction terms.

$$
G_{f,p}^{(k)}(r)=\int dr'\frac{f^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)+-\gamma\sum_{i}\varphi_{i}(r)\left(\int dr'\frac{y_{i}^{\dagger(b)}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)+\hat{W_{y}}^{0,b}\varphi_{k}(r)

$$
$$
H_{f,p}^{(k)}(r)=\int dr'\frac{f^{(b)}(r',r')}{\left|r-r'\right|}\varphi_{k}(r)-\gamma\sum_{i}f_{i}^{(b)}(r)\left(\int dr'\frac{\varphi_{i}^{\dagger}(r')\varphi_{k}(r')}{\left|r-r'\right|}\right)+W_{f,p}^{(k)}(r)

$$

and G is the conjugate of H