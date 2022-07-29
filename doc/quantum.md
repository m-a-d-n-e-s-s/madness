# MADNESS quantum chemistry


## What you can do
Quantum chemical capabilities are:
 * SCF: Hartree-Fock and DFT, also in strong magnetic fields: `moldft, nemo, znemo`
 * analytical gradients for SCF: `moldft, nemo, znemo`
 * second nuclear derivatives for SCF (field-free case only), `nemo`
 * CIS, `cis, zcis`
 * response properties, `molresponse`
 * MP2, `cc2`
 * ADC(2), CIS(D), `cc2`
 * CC2 ground and excited states, `cc2`
 * PNO-MP2, `pno`
 * OEP/RKS, `oep`

## Quickstart
A quick help and an overview over all available codes can be
obtained by 

`madqc --help`

All programs can read commandline options or an input file (by default this is named "input").
A full list of all available calculation parameters can be obtained by writing

`qccode --help`
 
where `qccode` stands for any of the qc codes (e.g. moldft, cc2, nemo, ..)


## Calculation parameters
All calculations require parameters, specifying the molecule, the quantum chemical model, charge etc.
If no parameters are given, default parameters will be used. Only the molecule itself must be specified by the user.
Parameters can be specified through an input file and/or command line options.

Some parameters depend on each other and are set automatically, e.g. certain numerical parameters, or the use of 
point group symmetry. Parameters set by the user are always respected.

In the output files of the calculations the complete set of input parameters are printed out, 
together with a short description and further information.

You can see the full list of parameters by typing
> `qccode --print_parameters`

where, again, `qccode` stands for any of the qc codes

## Input file
The input file consists of data groups, starting with the relevant keyword, e.g. "dft" and ending with "end".
All parameters in a data group are given as key/value pairs, where the value can be an integer, a double, a string
even a vector or pair. 

A sample input file looks like
>dft\
>  charge 1          # comment\
>  ncf (slater,2.0) # value is a pair of string and double\
>end
> 
> geometry\
>  O 0 0 0\
>  H 0 1 1\
>  H 0 1 -1\
> end

Blank lines are ignored, as is everything after a hashtag. 
The input file is case-insensitive.
The key/value pairs can be separated by blanks or by the equal sign.
Pairs and vectors must be encased in parantheses or brackets, their entries must be separated by commas.

All programs will output the complete list of input options. You can always run 
> `qccode --help` 

which will output the input parameters and the copy/paste the options verbatim.
Note the once an option appears in the input file it will be considered user-defined and will override all default or derived values.

## Command line options
The data groups in the input file can also be set or augmented through the command line, e.g. the following
line will pass the same calculation parameters as the input file above.

`nemo --dft="charge=1; ncf=(slater,2.0)"`

Different key/value pairs are separated by a semicolon to indicate a newline.
If a given parameter is specified both in the input file and the command line, the command line parameters have 
priority over input file parameters.

The name of the input file can be changed by
> `nemo --input=customfile`


## Numerical parameters
Numerical parameters are specified in the dft block, the most important ones that will
also be used is subsequent calculations (e.g. mp2) are
> dft\
>   k 5\
>   L 20\
> end
 
The polynomial order k can be chosed between 2 and 30, for SCF calculations a value of 7 to 9 is advised, for correlated
calculations it should be chosen 5 or 6.

The calculation box size L should be large enough to hold all relevant phyics. All wavefunctions should have decayed
to practically zero at the box boundary. Note that excited states or anions can be quite large. 
The box size is given in atomic units (1 a.u. = 52 pm).

Generally it is advisable to use as few numerical parameters as possible, as they are usually dependent on each other

## Geometry input
The geometry of the molecule is given in the geometry data group. 
By default atomic units are used, but angstrom can be switched on by adding the line
> geometry\
> units angs\
> ..\
> end

The following example will read an external xyz file, using angstrom by default
>`geometry`\
> `  source_type xyz          # optional `\
> `  source_name h2o.xyz`\
> `end`

or you can use the command line options using the convenience short option
> `nemo --geometry=h2o.xyz`

A small number of geometries are stored in a library, accessible through
> `nemo --geometry="source_type=library; source_name=h2o"`
 
If no source type is given it will be deduced from the file name, if the source is ambiguous,
e.g. a structure in the library has the same same as an input file, the code will stop.

## Geometry optimization
For the following codes/methods there are gradients implemented:
> `nemo`, `moldft`, `znemo`

The optimization is started with the `gopt` flag in the `dft` block, geometry optimization
parameters are set in the `geoopt` block.
> `nemo --dft="k=8; econv=1.e-5; gopt=1"  --geoopt="maxiter=10" --geometry="source=library, h2o"`

## Convenience short options
`--optimize` optimize the geometry\
`--geometry=file.xyz` find the geometry in the xyz file (note Angstrom units!)
