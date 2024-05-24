=======================
Response Properties
=======================

A response function is a measure of how a property of a system changes in 
the presence of one or more perturbations.  We use the same notation used by the developers of Dalton :cite:t:`aidasDaltonQuantumChemistry2014a` for response properties.

.. list-table:: Response Function Notation
   :widths: 30 20 50
   :header-rows: 1

   * - Response Function
     - Notation
     - Description
   * - Linear
     - :math:`\langle\langle A;B \rangle\rangle_{\omega_b}`
     - Correction to expectation value of A due to perturbations in B.
   * - Quadratic
     - :math:`\langle\langle A;B,C \rangle\rangle_{\omega_b,\omega_c}`
     - Correction to expectation value of A due to perturbations in B and C.
   * - Cubic
     - :math:`\langle\langle A;B,C,D \rangle\rangle_{\omega_b,\omega_c,\omega_d}`
     - Correction to expectation value of A due to perturbations in B, C and D.

In other words, response functions are defined as the change in expectation values 
due to perturbations B, C, D, etc. where each perturbation is associated with a
frequency :math:`\omega_b`, :math:`\omega_c`, :math:`\omega_d`.
The perturbations can be considered as external monochromatic fields, or static
perturbations, where the frequency is zero.
In general, the perturbations represent Fourier components of an arbitrary time-dependent perturbation.

Linear Response
===============


The polarizability is a linear response function, and is defined as the change in the expectation value of the dipole moment due to a perturbation in the electric field. :math:`\alpha_{ij} = -\langle \langle x_i; x_j \rangle \rangle_{\omega}`
which for TDDFT is given by the trace of the dipole with the response density :math:`\alpha_{ij} = tr(x_i \rho^{(1)}(\omega))`.

As an example, to compute the Hartree-Fock polarizability the input file to `maddft` will be:

.. code-block:: none

    dft
        xc hf
        econv 0.01
        protocol [0.0001]
    end

    response
        dipole true
        first_order true
        freq_range [0.0,0.056,0.1]
        kain true
        maxiter 10
        maxsub 10
        protocol [0.0001]
    end

or in JSON format:

.. code-block:: json

    {
        "dft": {
            "econv": 0.01,
            "protocol": [
                0.0001
            ]
        },
        "response": {
            "first_order": true,
            "dipole": true,
            "freq_range": [
                0.0,
                0.056,
                0.1
            ],
            "kain": true,
            "maxiter": 10,
            "maxsub": 10,
            "protocol": [
                0.0001
            ],
        }
    }

The corresponding calculation computes the linear response in the density with respect to dipole operator at frequency 0.0, 0.056 and 0.1 a.u.

From the linear response function, one can also compute information about the spectrum of TDDFT Hamiltonian, such as the excitation energies and oscillator strengths.
To calculate the three lowest excitation energies and dipole transition moments for the three lowest excited states one can modify the input file as follows:
 

.. code-block:: none

    dft
        xc hf
        econv 0.01
        protocol [0.0001]
    end

    response
        first_order true
        excited_states true
        states 3
    end

or in JSON format:

Quadratic Response
==================

An example of a quadratic response function is the first hyperpolarizability.
To compute the first hyperpolarizability the input file to `maddft` will be:

.. code-block:: none

    dft
        xc hf
        econv 0.01
        protocol [0.0001]
    end

    response
        dipole true
        first_order true
        freq_range [0.0,0.056,0.1]
        qudratic true
    end

where the :code:`quadratic` keyword is set to true to compute the quadratic response function.
In this case, the first hyperpolarizability will be computed from all mix of all possible frequencies in the :code:`freq_range`.


Two-absorption amplitude can be computed from second-order response :math:`\gamma^{\alpha\beta}`.
From pole-analysis the two-photon absorption amplitude :math:`\sigma_{n\alpha\beta}` is

.. math:: 

   \sigma_{n\alpha\beta} = \langle X_n,Y_n \vert P^{(\alpha\beta)},Q^{(\alpha\beta)} \rangle

where :math:`\omega_{\alpha} +\omega_{\beta} = \omega_n` and :math:`P^{(\alpha\beta)}` and :math:`Q^{(\alpha\beta)}` are the
second order perturbation operators defined from corresponding first-order densities.

To define such a calculation in `maddft` the input file will be:

.. code-block:: none

    dft
        xc hf
        econv 0.01
        protocol [0.0001]
    end

    response
        dipole true
        qudratic true
        excited_states true
    end

In this case, the two-photon absorption is computed by first computing the excited states and the first-order response with 
respect to perturbation in the dipole operator at frequencies :math:`\omega_{n} = \omega_{\beta} + \omega_{\alpha}`.



On beta.json file
=================

beta.json prints the quadratic response at all frequency non-redudant combinations of the frequencies in the freq_range. f

For each frequency, we print 10 components of the first hyperpolarizability tensor :math:`\beta_{ijk}`, the minimal
number of components needed to fully determine each component by symmetry.  

.. code-block:: none

   A B C

1   x y z -> 6
2   x x x -> 1
3   x y y -> 3
4   x z z -> 3
5   y x x -> 3
6   y y y -> 1
7   y z z -> 3
8   z x x -> 3
9   z y y -> 3
10  z z z -> 1



To compute a single component \beta_{ABC}(\omega_a;\omega_b,\omega_c) we use the following formula:

.. math::
   \braket{\braket{A;B,C}} = \braket{X^{(A)}| V^{(BC)}} + \braket{\zeta^{(BC)}_x | v^{(A)} | \zeta^{(BC)}_y} + \braket{\zeta^{(CB)}_x | v^{(A)} | \zeta^{(CB)}_y}
    

In order to compute this we need to define :math:`\zeta^{(BC)}_x` and :math:`\zeta^{(BC)}_y` as well as second order perturbation operators :math:`V^{(BC)}`, 
for pairs XY, XX, YY, ZZ.  From there, we can compute the 10 above components of the first hyperpolarizability tensor.


.. code-block:: none

    X ;XX
    X ;YY
    X ;ZZ
    Y ;XX
    Y ;YY
    Y ;ZZ
    Z ;XX
    Z ;YY
    Z ;ZZ
    X ;YZ


Second order perturbation operators are defined as:

.. math::

    \begin{align*}
    V_p^{(BC)}(r) = & - \hat{Q} \hat{g_1} \bqty{ \hat{\zeta}^{(BC)}} \ket{\phi_p} - \hat{Q} \hat{g_1} \bqty{ \hat{\zeta}^{(BC)}} \ket{\phi_p} \\
    &-\hat{Q} \hat{g_2} \bqty{ \gamma^{(B)}\gamma^{(C)}+\gamma^{(C)}\gamma^{(B)}} \ket{\phi_p} \\
    &-\hat{Q} \hat{F}^{(B)} \ket{x_p^(C)} - \hat{Q} \hat{F}^{(C)} \ket{x_p^(B)}
    \end{align*}







