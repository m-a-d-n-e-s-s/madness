# MADNESS quantum chemistry


## What you can do
Quantum chemical capabilities are:
 * SCF: Hartree-Fock and DFT, also in strong magnetic fields
 * analytical gradients for SCF 
 * second nuclear derivatives for SCF (field-free case only)
 * CIS
 * response properties
 * MP2 
 * ADC(2), CIS(D)
 * CC2 ground and excited states
 * PNO-MP2

## Quickstart
A quick help and an overview over all available codes can be
obtained by 
> madqc --help

All programs can read commandline options or an input file (by default this is named "input").
A full list of all available calculation parameters can be obtained by writing
> qccode --help
 
where qccode stands for any of the qc codes (e.g. moldft, cc2, nemo, ..)


## Calculation parameters
All calculations require parameters, specifying the molecule, the quantum chemical model, charge etc.
If no parameters are given, default parameters will be used. Only the molecule itself must be specified by the user.
Paramters can be specified through an input file and/or command line options

Some parameters depend on each other and are set automatically, e.g. certain numerical parameters, or the use of 
point group symmetry. 

In the output files of the calculations the complete set of input parameters are printed out, 
together with a short description and further information.

## Input file
The input file consists of data groups, starting with the relevant keyword, e.g. "dft" and ending with "end".
All parameters in a data group are given as key/value pairs, where the value can be an integer, a double, a string
even a vector or pair
>dft\
>  charge 1          # comment\
>  ncf (slater,2.0) # value is a pair of string and double\
>end

Blank lines are ignored, as is everything after a hashtag. 
The input file is case-insensitive.
The key/value pairs can be separated by blanks or by the equal sign.
Pairs and vectors must be encased in parantheses or brackets, their entries must be separated by commas.

## Command line options
The data groups in the input file can also be set or augmented through the command line, e.g. the following
line will pass the same calculation parameters as the input file above.
> nemo --dft="charge=1; ncf=(slater,2.0)"

Different key/value pairs are separated by a semicolon to indicate a newline.
If a given parameter is specified both in the input file and the command line, the command line parameters have 
priority over input file parameters.


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
