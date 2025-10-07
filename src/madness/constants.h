/*
  This file is part of MADNESS.
  
  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:
  
  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367
  
  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

/**
 \file constants.h
 \brief Defines common mathematical and physical constants.
 \ingroup libraries

 \todo We should consider a more uniform naming scheme for all of the conversion ratios.
*/

#ifndef MADNESS_CONSTANTS_H
#define MADNESS_CONSTANTS_H

namespace madness {

    namespace constants {

        // Mathematical constants
        constexpr double pi = 3.14159265358979323846264338328; ///< Mathematical constant \f$\pi\f$.
        constexpr double sqrt_pi = 1.77245385090551602729816748334; ///< Mathematical constant \f$\sqrt{\pi}\f$.
        constexpr double inv_sqrt_pi = 0.564189583547756286948079451561; ///< Mathematical constant \f$\pi^{-1/2}\f$.

        // Misc physical constants and units

        constexpr double calorie_joule_relationship = 4.184000000 ; ///< 1 kcal/mol = 4.184 kJ/mol.

        constexpr double Debye = 3.3356409519815204e-30 ; ///< Cm (coulomb metre).

        // The following physical constants and units were obtained from NIST http://physics.nist.gov/constants
        // on 29/09/2025.  The comments contain the uncertainty and units.

        /// Mass constant in atomic units.
        constexpr double atomic_mass_constant = 1.66053906892e-27 ; //    0.00000000052e-27    kg

        /// First hyperpolarizability in atomic units.
        constexpr double atomic_unit_of_1st_hyperpolarizablity = 3.2063612996e-53 ; //   0.0000000015e-53    C^3 m^3 J^-2

        /// Second hyperpolarizability in atomic units.
        constexpr double atomic_unit_of_2nd_hyperpolarizablity = 6.2353799735e-65  ; //    0.0000000039e-65     C^4 m^4 J^-3

        /// Action in atomic units.
        constexpr double atomic_unit_of_action = 1.054571817e-34 ; //   (exact)    J s

        /// Charge in atomic units.
        constexpr double atomic_unit_of_charge = 1.602176634e-19 ; //   (exact)    C

        /// Charge density in atomic units.
        constexpr double atomic_unit_of_charge_density = 1.08120238677e12 ; //     0.00000000051 e12     C m^-3

        /// Current in atomic units.
        constexpr double atomic_unit_of_current = 6.6236182375082e-3 ; //     0.0000000000072e-3      A

        /// Electric dipole moment in atomic units.
        constexpr double atomic_unit_of_electric_dipole_moment = 8.4783536198e-30 ; //    0.0000000013e-30     C m

        /// Electric quadrupole moment in atomic units.
        constexpr double atomic_unit_of_electric_quadrupole_moment = 4.4865515185e-40 ; //    0.0000000014e-40     C m^2

        /// Electric field in atomic units.
        constexpr double atomic_unit_of_electric_field = 5.14220675112e11 ; //      0.00000000080e11      V m^-1

        /// Electric field gradient in atomic units.
        constexpr double atomic_unit_of_electric_field_gradient = 9.7173624424e21 ; //  0.0000000030e21      V m^-2

        /// Electric polarizability in atomic units.
        constexpr double atomic_unit_of_electric_polarizablity = 1.64877727212e-41 ; //   0.00000000051e-41   C^2 m^2 J^-1

        /// Electric potential in atomic units.
        constexpr double atomic_unit_of_electric_potential = 27.211386245981 ; //        0.000000000030          V

        /// Energy in atomic units.
        constexpr double atomic_unit_of_energy = 4.3597447222060e-18 ; //     0.0000000000048e-18     J

        /// Force in atomic units.
        constexpr double atomic_unit_of_force = 8.2387235038e-8 ; //     0.0000000013e-8      N

        /// Length in atomic units.
        constexpr double atomic_unit_of_length = 0.529177210544e-10 ; // 0.000000000082e-10 m

        /// Magnetic dipole moment in atomic units.
        constexpr double atomic_unit_of_magnetic_dipole_moment = 1.85480201315e-23 ; //   0.00000000058e-23    J T^-1

        /// Magnetic flux density in atomic units.
        constexpr double atomic_unit_of_magnetic_flux_density = 2.35051757077e5 ; //      0.00000000073 e5      T

        /// Magnetizability in atomic units.
        constexpr double atomic_unit_of_magnetizability = 7.8910365794e-29 ; //    0.0000000049e-29    J T^-2

        /// Mass in atomic units.
        constexpr double atomic_unit_of_mass = 9.1093837139e-31 ; //    0.0000000028e-31     kg

        /// Momentum in atomic units.
        constexpr double atomic_unit_of_momentum = 1.99285191545e-24 ; //   0.00000000031e-24    kg m s^-1

        /// Permittivity in atomic units.
        constexpr double atomic_unit_of_permittivity = 1.11265005620e-10 ; //  0.00000000017e-10            F m^-1

        /// Time in atomic units.
        constexpr double atomic_unit_of_time = 2.4188843265864e-17 ; // 0.0000000000026e-17 s

        /// Velocity in atomic units.
        constexpr double atomic_unit_of_velocity = 2.18769126216e6 ; //  0.00000000034 e6     m s^-1

        /// Avogadro's number.
        constexpr double Avogadro_constant =  6.02214076e23 ; //     (exact)      mol^-1

        /// Bohr magneton.
        constexpr double Bohr_magneton = 927.40100657e-26  ; //    0.00000029e-26        J T^-1

        /// Bohr radius.
        constexpr double Bohr_radius = 0.529177210544e-10 ; // 0.000000000082e-10 m

        /// Boltzmann constant.
        constexpr double Boltzmann_constant = 1.380649e-23 ; //       (exact)       J K^-1

        /// Compton wavelength.
        constexpr double Compton_wavelength = 2.42631023538e-12 ; //  0.00000000076e-12   m

        /// Quantum of conductance, \f$ 2e^2/h \f$.
        constexpr double conductance_quantum = 7.748091729e-5 ; //  (exact)    S

        /// Electron \f$g\f$ factor.
        constexpr double electron_g_factor = -2.00231930436092 ; //  0.00000000000036

        /// Electron gyromagnetic ratio.
        constexpr double electron_gyromagnetic_ratio = 1.76085962784e11 ; //   0.00000000055 e11     s^-1 T^-1

        /// Electron magnetic moment.
        constexpr double electron_magnetic_moment = -928.47646917e-26 ; //    0.00000029e-26        J T^-1

        /// Electron magnetic moment anomaly.
        constexpr double electron_magnetic_moment_anomaly =  1.15965218046e-3 ; // 0.00000000018e-3

        /// Ratio between the electron magnetic moment and Bohr magneton.
        constexpr double electron_magnetic_moment_to_Bohr_magneton_ratio =  -1.00115965218046 ; // 0.00000000000018

        /// Ratio between the electron magnetic moment and nuclear magneton.
        constexpr double electron_magnetic_moment_to_nuclear_magneton_ratio = -1838.281971877 ; //     0.000000032

        /// Electron mass.
        constexpr double electron_mass = 9.1093837139e-31 ; //    0.0000000028e-31     kg

        /// Ratio of the electron to proton mass.
        constexpr double electron_proton_mass_ratio = 5.446170214889e-4 ; //   0.000000000094e-4

        /// Electron volt.
        constexpr double electron_volt = 1.602176634e-19 ; //   (exact)    J

        /// Electron volt to Hartree conversion.
        constexpr double electron_volt_hartree_relationship = 3.6749322175665e-2 ; //     0.0000000000040e-2     E_h

        /// Electron volt to Hertz relationship.
        constexpr double electron_volt_hertz_relationship = 2.417989242e14 ; //     (exact)     Hz

        /// Electron volt to Joule relationship.
        constexpr double electron_volt_joule_relationship = 1.602176634e-19 ; //   (exact)    J

        /// Elementary charge.
        constexpr double elementary_charge = 1.602176634e-19 ; //    (exact)    C

        /// Faraday constant.
        constexpr double Faraday_constant = 96485.33212 ; //      (exact)     C mol^-1

        /// Fermi coupling constant.
        constexpr double Fermi_coupling_constant = 1.1663787e-5 ; // 0.0000006e-5          GeV^-2

        /// Fine structure constant.
        constexpr double fine_structure_constant = 7.2973525643e-3 ;  //   0.0000000011e-3

        /// Hartree to electron volt relationship.
        constexpr double hartree_electron_volt_relationship = 27.211386245981 ; //        0.000000000030          eV

        /// Hartree energy in Joules.
        constexpr double Hartree_energy = 4.3597447222060e-18 ; //     0.0000000000048e-18     J

        /// Hartree energy in Hertz.
        constexpr double hartree_hertz_relationship = 6.5796839204999e15 ; // 0.0000000000072e15 Hz

        /// Hartree energy in inverse meters.
        constexpr double hartree_inverse_meter_relationship = 2.1947463136314e7 ; //  0.0000000000024 e7  m^-1

        /// Hartree energy in Joules.
        constexpr double hartree_joule_relationship = 4.3597447222060e-18 ; //    0.0000000000048e-18     J

        /// Hartree energy in Kelvin.
        constexpr double hartree_kelvin_relationship = 3.1577502480398e5 ; //        0.0000000000034 e5         K

        /// Hertz energy in electron volts.
        constexpr double hertz_electron_volt_relationship = 4.135667696e-15 ; //     (exact)     eV

        /// Hertz energy in Hartrees.
        constexpr double hertz_hartree_relationship = 1.5198298460574e-16 ; // 0.0000000000017e-16 E_h

        /// Hertz energy in Joules.
        constexpr double hertz_joule_relationship = 6.62607015e-34 ; //  (exact)     J

        /// Hertz energy in Kelvin.
        constexpr double hertz_kelvin_relationship = 4.799243073e-11 ; //      (exact)       K

        /// Nuclear magneton.
        constexpr double nuclear_magneton = 5.0507837393e-27 ; //     0.0000000016e-27     J T^-1

        /// Nuclear magneton in electron volts per Tesla.
        constexpr double nuclear_magneton_in_eV_per_T = 3.15245125417e-8 ; //   0.00000000098e-8    eV T^-1

        /// Planck's constant.
        constexpr double Planck_constant = 6.62607015e-34 ; //    (exact)     J s

        /// Reduced Planck's constant.
        constexpr double Planck_constant_over_2_pi = 1.054571817e-34 ; //   (exact)    J s

        /// Ratio of proton to electron mass.
        constexpr double proton_electron_mass_ratio = 1836.152673426 ; //      0.000000032

        /// Proton \f$ g \f$ factor.
        constexpr double proton_g_factor = 5.5856946893 ; //        0.0000000016

        /// Proton gyromagnetic ratio.
        constexpr double proton_gyromagnetic_ratio = 2.6752218708e8 ; //     0.0000000011 e8      s^-1 T^-1

        /// Proton magnetic moment.
        constexpr double proton_magnetic_moment = 1.41060679545e-26 ; //   0.00000000060e-26    J T^-1

        /// Proton mass.
        constexpr double proton_mass = 1.67262192595e-27 ; //   0.00000000052e-27    kg

        /// Rydberg constant.
        constexpr double Rydberg_constant = 10973731.568157 ; //    0.000012             m^-1

        /// Speed of light in a vacuum.
        constexpr double speed_of_light_in_vacuum = 299792458 ; //           (exact)               m s^-1

        /// Stefan-Boltzmann constant.
        constexpr double Stefan_Boltzmann_constant = 5.670374419e-8 ; //      (exact)         W m^-2 K^-4

        /// Unified atomic mass unit.
        constexpr double unified_atomic_mass_unit = 1.66053906892e-27 ; //   0.00000000052e-27    kg

        /// Atomic mass in atomic units
        constexpr double atomic_mass_in_au = atomic_mass_constant/atomic_unit_of_mass;

        /// conversion from atomic units in reciprocal centimeter
        constexpr double au2invcm = 219474.63136314;  // 0.00000024  cm^-1

        /// the dielectric constant \f$\epsilon_0\f$, or the permittivity of vacuum
        constexpr double dielectric_constant = 8.8541878188e-12;     //   0.0000000014   F m^{-1}

        /// speed of light in vacuum in au
        constexpr double speed_of_light_in_vacuum_in_au = 1.0/fine_structure_constant;   // \approx 137

    }

}

#endif
