MADNESS Skyrme-Hartree-Fock code
--------------------------------


For general information, see:

Sagert, I.; Fann, G. I.; Fattoyev, F. J.; Postnikov, S.; Horowitz, C. J.</br>
Physical Review C, Volume 93, Issue 5, id.055801, 2016


To compile the code after building MADNESS:
```
make clean
make mshf
```

When running the code, make sure that the input files ```mshf_skyrme``` and ```mshf_input``` are in the same directory as the executable. Otherwise, mshf will run a default calculation of Oxygen-16. For more information on the code, see the manual in the ```mshf_manual``` folder. 
