# MADNESS-QCEngine Documentation

Purpose: To provide a single executable for running chemistry applications. This will be 
the main executable used within QCEngine. 

## Current Functionality

- hf/dft 
    - geometry optimization
- linear response
    - excited states
    - dipole
    - nuclear
- response properties
    - polarizability
    - hyperpolarizability
    - gradients of excited states and response states - optimized for excited states


The executable is designed to be a single executable that can be used to run all of the above from either 
an input file in the MADNESS format or JSON input.

Examples:


Mixed response

```
dft
    econv 0.01
    protocol [0.0001]
end

excited_states
    states 5
    protocol [0.0001]
end

response
    perturbation dipole nuclear excited
    freq_range [0.0,0.056,0.1]
    kain true
    maxiter 10
    maxsub 10
    protocol [0.0001]
end

quadratic 
    perturbation dipole nuclear # maybe?
    freq_range [0.0,0.056,0.1]
    protocol [0.0001]
end

```

