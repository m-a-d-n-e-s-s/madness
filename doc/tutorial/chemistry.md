# Chemistry in Madness
## Running a calculation -- quick and dirty
Running a quantum chemical calculations requires a molecule -- and not much more:
```
moldft --geometry=h2o.xyz 
```
This will run a HF calculation on the molecule specified in the `h2o.xyz`-file.
A geometry optimization is performed by 
```
moldft --geometry=h2o.xyz  --optimize
```


## Getting help
All quantum chemical codes (e.g. `nemo`, `moldft`, etc) have two options
```
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
