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

#include "madness.h"
#include <chem/pcm.h>
#include <chem/atomutil.h>
#include <chem/molecule.h>


#ifdef MADNESS_HAS_PCM

#include "PCMSolver/pcmsolver.h"
#include "PCMSolver/PCMInput.h"

namespace madness {

namespace detail {

/// to convert the apparent surface charge potential to a real_function_3d
struct asc_potential : public FunctionFunctorInterface<double,3> {

    /// the coordinates of the apparent surface charges
    std::vector<coord_3d> asc_coord;

    /// the charges of the apparent surface charges
    std::vector<double> asc_value;

    /// regularization parameter for the (singular) surface charges
    double rcut;

    /// constructor

    /// @param[in]  size    the grid size (number of apparent surface charges, ASC)
    /// @param[in]  coord   the coordinates of the ASC
    /// @param[in]  val     the charges of the ASC
    asc_potential(const int size, const Tensor<double> coord, const Tensor<double> val)
        : rcut(3.0) {

        asc_coord.resize(size);
        asc_value.resize(size);
        for (int i=0; i<size; ++i) {
            coord_3d c;
            c[0]=coord(3*i);
            c[1]=coord(3*i+1);
            c[2]=coord(3*i+2);
            asc_coord[i]=c;
            asc_value[i]=val(i);
        }

        // the smoothing parameter is based on the assumption that
        // electron density is high near the singularity, but the ASC
        // are at the van-der-Waals distance. We can savely use a smaller rcut.
//        double Z=*std::max_element(asc_value.begin(),asc_value.end());
//        double eprec=1.e-4;
//        rcut=smoothing_parameter(Z,eprec);

    }

    /// return the value of the ASC potential
    double operator()(const coord_3d& r) const {

        double sum = 0.0;
        for (unsigned int i=0; i<asc_value.size(); ++i) {
            double rr=(r-asc_coord[i]).normf();
            sum += asc_value[i] * smoothed_potential(rr*rcut)*rcut;
        }
        return sum;

    }
};

void host_writer(const char * message) {
    std::string msg(message);
    print(msg);
}


}

PCM::PCM(World& world, const Molecule& mol, const std::string pcm_data,
        const bool verbose) : pcm_context(0), mep_lbl("TotMEP"), asc_lbl("totASC") {
    if (!pcmsolver_is_compatible_library()) {
        fprintf(stderr, "%s\n", "PCMSolver library not compatible");
        exit(EXIT_FAILURE);
    }

    // convert molecule to pcm format
    const int natom=mol.natom();
    charges=Tensor<double>(mol.natom());
    coordinates=Tensor<double>(3*mol.natom());
    double* ch=charges.ptr();
    double* c=coordinates.ptr();
    for (int iatom=0; iatom<mol.natom(); ++iatom) {
        ch[iatom]=mol.get_atom_charge(iatom);
        c[3*iatom   ]=mol.get_atom(iatom).x;
        c[3*iatom +1]=mol.get_atom(iatom).y;
        c[3*iatom +2]=mol.get_atom(iatom).z;
    }

    // default input
    struct PCMInput host_input = pcmsolver_input();

    // use the explicit input reader provided by PCMSolver
    std::stringstream spcm_data(pcm_data);
    std::string word;
    spcm_data >> word;

    pcmsolver_reader_t pcmsolver_reader=PCMSOLVER_READER_HOST;
    if (word=="PCMSOLVER_READER_OWN") {
        pcmsolver_reader=PCMSOLVER_READER_OWN;
    } else {
        // convert madness input to PCM input
        if (not word.empty()) strcpy(host_input.solvent, pcm_data.c_str());
    }

    // there is only one option allowed
    if (not spcm_data.eof()) {
        print("the pcm flag had more than one option");
        MADNESS_EXCEPTION("only one pcm option is allowed",0);
    }

    // This means the molecular point group has three generators:
    // the Oxy, Oxz and Oyz planes
    // we don't use symmetry
    symmetry_info=Tensor<int>(4);

    pcm_context =std::shared_ptr<pcmsolver_context_t> (
            pcmsolver_new(pcmsolver_reader, natom, charges.ptr(), coordinates.ptr(),
                    symmetry_info.ptr(), &host_input, detail::host_writer),
                    pcmsolver_delete);

    if (verbose and (world.rank()==0)) pcmsolver_print(pcm_context.get());

}

PCMInput PCM::pcmsolver_input() const {
    PCMInput host_input;

    // These parameters would be set by the host input reading
    // Length and area parameters are all assumed to be in Angstrom,
    // the module will convert to Bohr internally
    strcpy(host_input.cavity_type, "gepol");
    host_input.patch_level = 2;
    host_input.coarsity = 0.5;
    host_input.area = 0.2;
    host_input.min_distance = 0.1;
    host_input.der_order = 4;
    host_input.scaling = true;
    strcpy(host_input.radii_set, "bondi");
    strcpy(host_input.restart_name, "cavity.npz");
    host_input.min_radius = 100.0;

    strcpy(host_input.solver_type, "iefpcm");
    strcpy(host_input.solvent, "water");
    strcpy(host_input.equation_type, "secondkind");
    host_input.correction = 0.0;
    host_input.probe_radius = 1.0;

    strcpy(host_input.inside_type, "vacuum");
    host_input.outside_epsilon = 1.0;
    strcpy(host_input.outside_type, "uniformdielectric");

    return host_input;
}

Tensor<double> PCM::nuclear_mep(int nr_nuclei, const Tensor<double>& charges,
        const Tensor<double>& coordinates, const int grid_size,
        const Tensor<double>& grid) const {
    Tensor<double> mep(grid_size);
    for (int i = 0; i < nr_nuclei; i++) {
        for (int j = 0; j < grid_size; j++) {
            // Column-major ordering. Offsets: col_idx * nr_rows + row_idx
            double dist = pow((coordinates(i * 3) - grid(j * 3)), 2) +
                    pow((coordinates(i * 3 + 1) - grid(j * 3 + 1)), 2) +
                    pow((coordinates(i * 3 + 2) - grid(j * 3 + 2)), 2);
            mep(j) += charges(i) / sqrt(dist);
        }
    }
    return mep;
}



real_function_3d PCM::compute_pcm_potential(const real_function_3d& coulomb_potential,
        const bool dynamic) const {

    MADNESS_ASSERT(coulomb_potential.is_initialized());
    const int grid_size = pcmsolver_get_cavity_size(pcm_context.get());

    Tensor<double> grid(3*grid_size);
    pcmsolver_get_centers(pcm_context.get(), grid.ptr());

    // compute the molecular electrostatic potential (mep) from the nuclei
    Tensor<double> mep = nuclear_mep(charges.size(), charges, coordinates, grid_size, grid);

    // nuclear potential is density independent
    if(dynamic) mep = 0.0;

    // add the electronic contribution to the mep
    for (int i=0; i<grid_size; ++i) {
        coord_3d evalpoint={grid(3*i),grid(3*i+1),grid(3*i+2)};
        mep[i]-=coulomb_potential(evalpoint);
    }

    // This is the Ag irreducible representation (totally symmetric)
    int irrep = 0;

    Tensor<double> asc(grid_size);

    // compute the contribution to the response kernel
    if (dynamic) {

        const std::string mep_neq_lbl="mep_neq_lbl";
        const std::string asc_neq_lbl="asc_neq_lbl";

        pcmsolver_set_surface_function(pcm_context.get(), grid_size, mep.ptr(), mep_neq_lbl.c_str());
        pcmsolver_compute_response_asc(pcm_context.get(), mep_neq_lbl.c_str(),asc_neq_lbl.c_str(), irrep);
        pcmsolver_get_surface_function(pcm_context.get(), grid_size, asc.ptr(), asc_neq_lbl.c_str());

    } else {
        pcmsolver_set_surface_function(pcm_context.get(), grid_size, mep.ptr(), mep_lbl.c_str());
        pcmsolver_compute_asc(pcm_context.get(), mep_lbl.c_str(), asc_lbl.c_str(), irrep);
        pcmsolver_get_surface_function(pcm_context.get(), grid_size, asc.ptr(), asc_lbl.c_str());
    }

    detail::asc_potential ascpot(grid_size,grid,asc);
    World& world = coulomb_potential.world();
    real_function_3d v=real_factory_3d(world).functor(ascpot);

    return -1.0*v;

}

double PCM::compute_pcm_energy() const {
    double pcm_energy =
            pcmsolver_compute_polarization_energy(pcm_context.get(), mep_lbl.c_str(), asc_lbl.c_str());
    return pcm_energy;
}


} /* namespace madness */
#else // MADNESS_HAS_PCM

namespace madness {

PCM::PCM(World& world, const Molecule& mol, const std::string pcm_data,
            const bool verbose) {
    MADNESS_EXCEPTION("no PCMSolver configured and available in MADNESS",1);
}

real_function_3d PCM::compute_pcm_potential(const real_function_3d& coulomb_potential,
        const bool dynamic) const {
    MADNESS_EXCEPTION("no PCMSolver configured and available in MADNESS",1);
    real_function_3d dummy;
    return dummy;
}

double PCM::compute_pcm_energy() const {
    MADNESS_EXCEPTION("no PCMSolver configured and available in MADNESS",1);
    return 0.0;
}

Tensor<double> PCM::nuclear_mep(int nr_nuclei, const Tensor<double>& charges,
                         const Tensor<double>& coordinates, const int grid_size,
                         const Tensor<double>& grid) const {
    MADNESS_EXCEPTION("no PCMSolver configured and available in MADNESS",1);
    return Tensor<double>();
}

} // namespace madness
#endif // MADNESS_HAS_PCM
