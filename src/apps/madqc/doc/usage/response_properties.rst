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

Hyperpolarizability
---------------------

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

In this calculation the components of the first hyperpolarizability will be computed by first computing the 
linear response vectors in all direction at all combinations of the frequencies in the freq_range.
In order to compute the quadratic response function from linear response functions the second-order perturbation
operators need to be computed. The second-order perturbation operators are defined as: 

.. math:: 
    \begin{equation}
	    \begin{aligned}
		V_p^{BC}(r) =        & -g^{'}_p[\hat{\zeta}^{BC}](r) -  g^{''}_p[\hat{\chi}^{B} \hat{\chi}^{C}](r)                     \\
		                     & -\hat{Q}^{0} \hat{F}^{B} x_p^{C}(r) + \sum_{k} x_k^{C}(r)F_{kp}^{B} \\
		                     & -\hat{Q}^{0} \hat{F}^{C} x_p^{B}(r) + \sum_{k} x_k^{B}(r)F_{kp}^{C}                                          \\
		V_p^{BC\dagger}(r) = & -g^{'}_p[\hat{\zeta}^{BC\dagger}](r)-  g^{''}_p[\hat{\chi}^{B\dagger} \hat{\chi}^{C\dagger}](r) \\
		                     & -\hat{Q}^0 \hat{F}^{B\dagger} y_p^{C}(r) + \sum_{k} y_k^{C}(r)F_{kp}^{*B}  \\
		                     & -\hat{Q}^0 \hat{F}^{C\dagger} y_p^{B}(r) + \sum_{k} y_k^{B}(r)F_{kp}^{*C}
	    \end{aligned}
    \end{equation}


In general, there are 27 unique elements to the hyperpolarizability tensor.  These are computed by constructing 9 unique second-order perturbation operators 
$XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ$ and tracing with the linear response vectors $X, Y, Z$. 
The hyperpolarizability tensor is outputed in beta.json file, which contains the 27 unique elements of the hyperpolarizability tensor at each frequency in the freq_range,
in dictionary format which can be easily read into a pandas dataframe.

Two-absorption Amplitude
-------------------------

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











