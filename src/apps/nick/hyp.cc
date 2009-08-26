//\file hyp.cc
//\brief Computes 1F1(a,b,z) internally using extended precision

//By: Robert Harrison
#include <nick/hyp.h>

/// Computes 1F1(a,b,z) internally using extended precision

/// If result is larger than 1.0, result should be accurate to
/// full double precision, otherwise result is accurate to
/// somewhat better than 1e-17.  However, the termination
/// test is not very smart so there may be failures.
complexd conhyp(const complexd& a_arg,
                const complexd& b_arg,
                const complexd& z_arg) {
    const double tol = 1e-15;

    // Save input precision so we can reset it on exit
    int nbits_save = extended_real::get_default_prec();
    if (nbits_save < 128) extended_real::set_default_prec(128);

    for (int attempts=0; attempts<3; attempts++) {
        // Convert input arguments to extended precision representation
        extended_complex a(a_arg);
        extended_complex b(b_arg);
        extended_complex z(z_arg);
        
        extended_complex sum(0.0,0.0);   // Accumulates the result
        extended_complex term(1.0,0.0);  // Current term in sum
        
        extended_real absprevterm = 999.0;      // Norm of previous term for convergence test
        extended_real maxsum = 0.0;             // Size of largest intermediate for precision check
        
        bool converged = false;
        for (int n=0; n<20000; n++) {
            sum += term;
            maxsum = max(maxsum,abs(sum));
            
            extended_complex cn(1.0*n,0.0);
            extended_complex cn1(1.0*(n+1),0.0);
            extended_complex ratio = ((a+cn)*z) / ((b+cn)*cn1);
            term*=ratio;
            
            double absterm = abs(term);

            if (absterm<tol && absprevterm<tol) {
                //std::cout << "Converged in " << n << " " << extended_real::get_default_prec() << " " << a_arg << " " << b_arg << " " << z_arg << std::endl;
                converged = true;
                break;
            }
            
            absprevterm = absterm;
        }
        
        if (!converged) 
            throw "no convergence in 20000 terms!";

        // Did we have enough precision? 
        int nbits = extended_real::get_default_prec();
        extended_real twon = extended_real(1) << nbits;
        extended_real test = maxsum*(100.0/tol);
        if (test < twon) {
            // Convert back to standard precision and return
            complexd result(sum.real(),sum.imag());

            // Restore input precision
            extended_real::set_default_prec(nbits_save);

            return result;
        }
        
        // Could get more intelligent here
        std::cout << "EXTENDING PRECISION " << nbits+128 << std::endl;
        extended_real::set_default_prec(nbits+128);
    }
    throw "insufficient precision";
}

// int main() {
//     cout.precision(15);

//     double k = 1.0;
//     double r = 20.0;

//     complexd a(0.0, -1.0/k);
//     complexd b(1.0, 0.0);
//     complexd z(0.0, -r*k);

//     complexd result = conhyp(a, b, z);
//     cout << result.real() << " " << result.imag() << endl;

//     return 0;
// }

    


