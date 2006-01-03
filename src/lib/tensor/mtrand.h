#ifndef MTRAND_H
#define MTRAND_H

namespace madness {
    
    /// Initialize Mersenne twister random number generator
    extern void init_genrand(unsigned long s);
    
    /// Initialize Mersenne twister random number generator
    
    /// Initialize by an array with array-length 
    /// init_key is the array for initializing keys
    /// key_length is its length 
    extern void init_by_array(unsigned long init_key[], int key_length);
    
    /// Generates a random number on [0,0xffffffff]-interval
    extern unsigned long genrand_int32(void);
    
    /// Generates a random number on [0,0x7fffffff]-interval
    extern long genrand_int31(void);
    
    /// Generates a random number on [0,1]-real-interval
    extern double genrand_real1(void);
    
    /// Generates a random number on [0,1)-real-interval
    extern double genrand_real2(void);
    
    /// Generates a random number on (0,1)-real-interval 
    extern double genrand_real3(void);
    
    /// Generates a random number on [0,1) with 53-bit resolution
    double genrand_res53(void) ;

}

#endif
