#ifndef MADNESS_CONSTANTS_H
#define MADNESS_CONSTANTS_H

/// \file constants.h
/// \brief Defines common mathematical and physical constants


namespace madness {
    namespace constants {

        // Mathematical constants
        const double pi = 3.14159265358979323846264338328; //< Mathematical constant pi


        // Misc physical constants and units
        const double calorie_joule_relationship = 4.184000000 ; // J  ... i.e.,  1 kcal/mol = 4.184 kJ/mol

        const double Debye = 3.335640035 e-30 ; // Cm  (coulomb metre)


        // The following physical constants and units were obtained from NIST http://physics.nist.gov/constants
        // on 2/2/2008.  The comments contain the uncertainty and units.

        const double atomic_mass_constant = 1.660 538 782 e-27 ; //    0.000 000 083 e-27    kg

        const double atomic_unit_of_1st_hyperpolarizablity = 3.206 361 533 e-53 ; //   0.000 000 081 e-53    C^3 m^3 J^-2

        const double atomic_unit_of_2nd_hyperpolarizablity = 6.235 380 95 e-65  ; //    0.000 000 31 e-65     C^4 m^4 J^-3

        const double atomic_unit_of_action = 1.054 571 628 e-34 ; //   0.000 000 053 e-34    J s

        const double atomic_unit_of_charge = 1.602 176 487 e-19 ; //   0.000 000 040 e-19    C

        const double atomic_unit_of_charge_density = 1.081 202 300 e12 ; //     0.000 000 027 e12     C m^-3

        const double atomic_unit_of_current = 6.623 617 63 e-3 ; //     0.000 000 17 e-3      A

        const double atomic_unit_of_electric_dipole_moment = 8.478 352 81 e-30 ; //    0.000 000 21 e-30     C m

        const double atomic_unit_of_electric_quadrupole_moment = 4.486 551 07 e-40 ; //    0.000 000 11 e-40     C m^2

        const double atomic_unit_of_electric_field = 5.142 206 32 e11 ; //      0.000 000 13 e11      V m^-1

        const double atomic_unit_of_electric_field_gradient = 9.717 361 66 e21 ; //  0.000 000 24 e21      V m^-2

        const double atomic_unit_of_electric_polarizablity = 1.648 777 2536 e-41 ; //   0.000 000 0034 e-41   C^2 m^2 J^-1

        const double atomic_unit_of_electric_potential = 27.211 383 86 ; //        0.000 000 68          V

        const double atomic_unit_of_energy = 4.359 743 94 e-18 ; //     0.000 000 22 e-18     J

        const double atomic_unit_of_force = 8.238 722 06 e-8 ; //     0.000 000 41 e-8      N

        const double atomic_unit_of_length = 0.529 177 208 59 e-10 ; // 0.000 000 000 36 e-10 m

        const double atomic_unit_of_magnetic_dipole_moment = 1.854 801 830 e-23 ; //   0.000 000 046 e-23    J T^-1

        const double atomic_unit_of_magnetic_flux_density = 2.350 517 382 e5 ; //      0.000 000 059 e5      T

        const double atomic_unit_of_magnetizability = 7.891 036 433 e-29 ; //    0.000 000 027 e-29    J T^-2

        const double atomic_unit_of_mass = 9.109 382 15 e-31 ; //    0.000 000 45 e-31     kg

        const double atomic_unit_of_momentum = 1.992 851 565 e-24 ; //   0.000 000 099 e-24    kg m s^-1

        const double atomic_unit_of_permittivity = 1.112 650 056 e-10 ; // (exact)               F m^-1

        const double atomic_unit_of_time = 2.418 884 326 505 e-17 ; // 0.000 000 000 016 e-17 s

        const double atomic_unit_of_velocity = 2.187 691 2541 e6 ; //  0.000 000 0015 e6     m s^-1

        const double Avogadro_constant =  6.022 141 79 e23 ; //     0.000 000 30 e23      mol^-1

        const double Bohr_magneton = 927.400 915 e-26  ; //    0.000 023 e-26        J T^-1

        const double Bohr_radius = 0.529 177 208 59 e-10 ; // 0.000 000 000 36 e-10 m

        const double Boltzmann_constant = 1.380 6504 e-23 ; //       0.000 0024 e-23       J K^-1

        const double Compton_wavelength = 2.426 310 2175 e-12 ; //  0.000 000 0033 e-12   m

        const double conductance_quantum = 7.748 091 7004 e-5 ; //  0.000 000 0053 e-5    S

        const double electron_g_factor -2.002 319 304 3622 ; //  0.000 000 000 0015    

        const double electron_gyromagnetic_ratio = 1.760 859 770 e11 ; //   0.000 000 044 e11     s^-1 T^-1

        const double electron_magnetic_moment = -928.476 377 e-26 ; //    0.000 023 e-26        J T^-1

        const double electron_magnetic_moment_anomaly =  1.159 652 181 11 e-3 ; // 0.000 000 000 74 e-3  

        const double electron_magnetic_moment_to_Bohr_magneton_ratio =  -1.001 159 652 181 11 ; // 0.000 000 000 000 74  

        const double electron_magnetic_moment_to_nuclear_magneton_ratio = -1838.281 970 92 ; //     0.000 000 80          

        const double electron mass = 9.109 382 15 e-31 ; //    0.000 000 45 e-31     kg

        const double electron_proton_mass_ratio = 5.446 170 2177 e-4 ; //   0.000 000 0024 e-4    

        const double electron_volt = 1.602 176 487 e-19 ; //   0.000 000 040 e-19    J

        const double electron_volt_hartree_relationship = 3.674 932 540 e-2 ; //     0.000 000 092 e-2     E_h

        const double electron_volt_hertz_relationship = 2.417 989 454 e14 ; //     0.000 000 060 e14     Hz

        const double electron_volt_joule_relationship = 1.602 176 487 e-19 ; //   0.000 000 040 e-19    J

        const double elementary_charge = 1.602 176 487 e-19 ; //    0.000 000 040 e-19    C

        const double Faraday_constant = 96 485.3399 ; //           0.0024                C mol^-1

        const double Fermi_coupling_constant = 1.166 37 e-5 ; // 0.000 01 e-5          GeV^-2

        const double fine_structure_constant = 7.297 352 5376 e-3 ;  //   0.000 000 0050 e-3    

        const double hartree_electron_volt_relationship = 27.211 383 86 ; //        0.000 000 68          eV

        const double Hartree_energy = 4.359 743 94 e-18 ; //     0.000 000 22 e-18     J

        const double hartree_hertz_relationship = 6.579 683 920 722 e15 ; // 0.000 000 000 044 e15 Hz

        const double hartree_inverse_meter_relationship = 2.194 746 313 705 e7 ; //  0.000 000 000 015 e7  m^-1

        const double hartree_joule_relationship = 4.359 743 94 e-18 ; //    0.000 000 22 e-18     J

        const double hartree_kelvin_relationship = 3.157 7465 e5 ; //        0.000 0055 e5         K

        const double hertz_electron_volt_relationship = 4.135 667 33 e-15 ; //     0.000 000 10 e-15     eV

        const double hertz_hartree_relationship = 1.519 829 846 006 e-16 ; // 0.000 000 000 010 e-16 E_h

        const double hertz_joule_relationship = 6.626 068 96 e-34 ; //  0.000 000 33 e-34     J

        const double hertz_kelvin_relationship = 4.799 2374 e-11 ; //      0.000 0084 e-11       K

        const double nuclear_magneton = 5.050 783 24 e-27 ; //     0.000 000 13 e-27     J T^-1

        const double nuclear_magneton_in_eV_per_T = 3.152 451 2326 e-8 ; //   0.000 000 0045 e-8    eV T^-1

        const double Planck_constant = 6.626 068 96 e-34 ; //    0.000 000 33 e-34     J s

        const double Planck_constant_over_2_pi = 1.054 571 628 e-34 ; //   0.000 000 053 e-34    J s

        const double proton_electron_mass_ratio = 1836.152 672 47 ; //      0.000 000 80          

        const double proton_g_factor = 5.585 694 713 ; //        0.000 000 046         

        const double proton_gyromagnetic_ratio = 2.675 222 099 e8 ; //     0.000 000 070 e8      s^-1 T^-1

        const double proton_magnetic_moment = 1.410 606 662 e-26 ; //   0.000 000 037 e-26    J T^-1

        const double proton_mass = 1.672 621 637 e-27 ; //   0.000 000 083 e-27    kg

        const double Rydberg_constant = 10 973 731.568 527 ; //    0.000 073             m^-1

        const double speed_of_light_in_vacuum = 299 792 458 ; //           (exact)               m s^-1

        const double Stefan_Boltzmann_constant = 5.670 400 e-8 ; //      0.000 040 e-8         W m^-2 K^-4

        const double unified_atomic_mass_unit = 1.660 538 782 e-27 ; //   0.000 000 083 e-27    kg

    }

}







#endif
