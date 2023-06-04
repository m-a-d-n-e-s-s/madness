# Chemistry in Madness

## Running a calculation -- quick and dirty

Running a quantum chemical calculations requires a molecule -- and not much more:

```shell
moldft --geometry=h2o.xyz 
```

This will run a HF calculation on the molecule specified in the `h2o.xyz`-file.
A geometry optimization is performed by

```shell
moldft --geometry=h2o.xyz  --optimize
```

## Getting help

All quantum chemical codes (e.g. `nemo`, `moldft`, etc) have two options

```shell
nemo --help
nemo --print_parameters
```

The `--help` option gives a short description of what the program does.
The `--print_parameters` option lists all parameters available to that program, together with
the current value, a marker if this parameter has been set by the user (`defined`), derived by the
program (`derived`) or left to its default value (`default`), and a brief description of that
parameter.
In some cases there is also a list of the available options in square brackets.

You can copy/paste the printout directly into your input file.

Parameters are case-insensitive.

## The input file -- General structure

Input parameters are passed to the program via the command line or through the input file.
If no input file name is given, the default name is `input`.
In the input file there are groups of parameters, enclosed by the data group key and the word `end`.

Each program uses its own data group, sometimes several, sometimes shared with other programs.
For instance, for a `cis` or an `mp2` calculation a HF reference is needed, so the `dft` data group
will also be used by `cis` and `mp2` and `mp2`.
The data groups will be listed by calling the program with `--print_parameters`.

Example for an HF calculation using `moldft` or `nemo`:

```
dft
  L 20
  k 8
  xc hf
end
```

Example for an excited state calculations using `cis`:

```
response
  freeze 1
  irrep b2
end
```

The input file may also be passed through the command line:

```
cis --dft="L=20; k=8; xc=hf" --response="freeze=1; irrep=b2"
```

will have the same effect as an input file the with two data groups from above.

There are 2 convenience command line options outside the general data group structure:

```
 --geometry=xxx.xyz
 --optimize
```

will use the geometry from the `xyz`-file and optimize the molecule.

### moldft and nemo

These two codes compute HF and DFT ground state wave functions.
Many of the parameters defined here will propagate to post-HF calculations, e.g.
the box size `L`, the polynomial order `k`, the orbital set (localized or canonical),
or the nuclear correlation factor (ncf):

Like `moldft`, the `nemo` code solves SCF equations, but uses a nuclear correlation factor (ncf)
for regularizing  the nuclear potential.

### CIS and CC2

`cis` and `cc2` use a converged `moldft` or `nemo` reference for computing excited states
and MP2/CC2 correlation energies, respectively.

### molresponse

`molresponse` computes frquency-dependent molecular response properties and excited-states.  We will focus on the frquency-dependent response here.


1. **Prepare Your Environment:** Before you can run a `molresponse` calculation,
   ensure that you've successfully installed the `molresponse` software and its
   dependencies. Refer to the official documentation for the installation guide.

2. **Generate Ground State Orbitals:** Before you compute molecular response
   properties, you need to generate the ground state orbitals. This is typically
   done with the `moldft` code and results in a `.restartdata` file.

3. **Configure Your Input File:** Create an input file (usually named
   `response.in`) specifying the parameters of your calculation. Here are the
   bare minimum parameters you'll need to define:

   - **Perturbation Operator:** Set by the `dipole` parameter. When set to
     `True`, the dipole operator is used.
   - **Perturbation Frequency:** Defined by the `omega` parameter.
   - **Ground-State Restart File:** By default, `molresponse` looks for
     `../moldft.restartdata`. You can specify a different path using the
     `archive` parameter.

   For a static response calculation, your `omega` value would be zero. For a
   frequency-dependent response, you would set `omega` to your desired
   frequency.

4. **Run Your Calculation:** Execute the `molresponse` program with your input
   file. The command is usually something like 
   
   ```cpp
   molresponse response.in
   ```

5. **Interpret Your Results:** The output is saved in a `response_base.json`
   file, which you can analyze to interpret your results.  

This is a very brief overview. I highly recommend reading the full tutorial [here](../../src/apps/molresponse/molresponse_tutorial.md) to fully understand how to use `molresponse`
effectively and accurately.

### MADNESS + MPI
To run a MADNESS application in parallel with MPI, you can simply run 

- `mpirun -n #procs qccode`

which will execute the given application with the specified amount of MPI processes. In addition, you will need to set the number of MADNESS threads by exporting the following variable

- `export MAD_NUM_THREADS=#threads`

By default, this is set to the number of CPUs available. When using MPI, each process will spawn the specified number of threads plus one additional communication thread **per process**. Don't use too many threads or the performance will be poor, it is a good idea to set the number to be lower than the total number of CPUs available and leave some capacity to the OS.

#### example

Consider a compute node with two CPUs with 14 cores each. Then

- `export MAD_NUM_THREADS=12`
-  `mpirun -n 2 moldft`

will result in 26 threads in total, leaving two cores to the OS. In general, if $n$ is the number of MPI processes and $m$ is the number of threads, the total number of threads will be $(m+1) \cdot n$. 

#### multiple nodes
If the calculation is distributed over multiple compute nodes, the number of MPI processes *per node* need to be specified via the `-ppn #procs` option, as well as the *total* number of processes. Using the same example as earlier but with ten compute nodes instead of one, we would get

- `mpirun -n 20 -ppn 2 moldft`

to run moldft on ten compute nodes with two MPI processes each and should result in a faster calculation. Performance can be further enhanced by pinning the processes to the socket, in order to avoid data being used by the process being stored on another node.

- `export I_MPI_PIN_DOMAIN=socket`
