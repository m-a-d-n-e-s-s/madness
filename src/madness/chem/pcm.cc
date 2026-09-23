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

#include <cstring>

#include "madness.h"
#include<madness/chem/pcm.h>
#include<madness/chem/atomutil.h>
#include<madness/chem/molecule.h>


// The PCMParameters group is compiled whether or not PCMSolver is available:
// madqc prints it either way. Only to_pcmsolver_input() needs the library.
namespace madness {

const std::map<std::string, std::pair<double, double> >&
PCMParameters::solvent_data() {
    // PCMSolver's twenty canonical solvents as {epsilon_static, probe radius/Angstrom},
    // transcribed from its src/utils/Solvent.cpp (the optical permittivity is not
    // exposed by PCMInput and therefore not kept here). Sorted by permittivity, as
    // upstream is. Six of the names are 16 characters or longer and so cannot be
    // written into PCMInput::solvent[16] with their terminating NUL -- those are
    // rejected in set_derived_values(), which uses this table to say what to use
    // instead.
    static const std::map<std::string, std::pair<double, double> > table = {
        {"n-heptane",            { 1.920, 3.125}},
        {"cyclohexane",          { 2.023, 2.815}},
        {"carbon tetrachloride", { 2.228, 2.685}},
        {"benzene",              { 2.247, 2.630}},
        {"1,4-dioxane",          { 2.250, 2.630}},
        {"toluene",              { 2.379, 2.820}},
        {"chloroform",           { 4.900, 2.480}},
        {"chlorobenzene",        { 5.621, 2.805}},
        {"aniline",              { 6.890, 2.800}},
        {"tetrahydrofurane",     { 7.580, 2.900}},
        {"methylenechloride",    { 8.930, 2.270}},
        {"1,2-dichloroethane",   {10.360, 2.505}},
        {"acetone",              {20.700, 2.380}},
        {"ethanol",              {24.550, 2.180}},
        {"methanol",             {32.630, 1.855}},
        {"acetonitrile",         {36.640, 2.155}},
        {"nitromethane",         {38.200, 2.155}},
        {"dimethylsulfoxide",    {46.700, 2.455}},
        {"propylenecarbonate",   {64.960, 1.385}},
        {"water",                {78.390, 1.385}}
    };
    return table;
}

PCMParameters::PCMParameters() {

    // All twenty canonical names are allowed here so that a misspelling is caught by
    // the framework with an accurate message; the six that do not fit PCMInput are
    // rejected separately in set_derived_values(), which can say what to use instead.
    // "explicit" means "no tabulated solvent": epsilon/probe_radius are used directly.
    initialize<std::string>("solvent","none","solvent surrounding the molecule",
            {"none","explicit",
             "n-heptane","cyclohexane","carbon tetrachloride","benzene","1,4-dioxane",
             "toluene","chloroform","chlorobenzene","aniline","tetrahydrofurane",
             "methylenechloride","1,2-dichloroethane","acetone","ethanol","methanol",
             "acetonitrile","nitromethane","dimethylsulfoxide","propylenecarbonate",
             "water"});
    initialize<std::string>("solver_type","iefpcm","the integral equation solver",
            {"iefpcm","cpcm"});
    initialize<double>("correction",0.0,"CPCM apparent-surface-charge scaling k; cpcm only");
    initialize<std::string>("cavity_type","gepol","cavity model; `restart` reads restart_name",
            {"gepol","restart"});
    initialize<double>("area",0.2,"average area of a tessera [Angstrom^2]");
    // PCMSolver also knows Allinger's MM3 radii, but "allinger" is eight characters and
    // PCMInput::radii_set is char[8], so the name cannot be passed through the C API.
    initialize<std::string>("radii_set","bondi","set of atomic radii building the cavity",
            {"bondi","uff"});
    initialize<bool>("scaling",true,"scale the atomic radii by 1.2");
    initialize<double>("min_radius",100.0,"minimal radius of an added sphere [Angstrom]");
    initialize<std::string>("restart_name","cavity.npz","the .npz cavity restart file");
    initialize<double>("probe_radius",1.0,"radius of the solvent probe [Angstrom]; "
            "derived from `solvent`, to be set only for solvent `explicit`");
    initialize<double>("epsilon",1.0,"static permittivity outside the cavity; "
            "derived from `solvent`, to be set only for solvent `explicit`");
}

void PCMParameters::set_derived_values(const CalculationParameters& cparam) {

    const std::string pcm_data = cparam.pcm_data();

    // The dft group's `pcm_data` is the trigger and carries the solvent name. It is a
    // *derived* value here, so an explicit `solvent` in the `pcm` group overrides it.
    reader_own_ = (pcm_data == "pcmsolver_reader_own");
    if (pcm_data != "none" and not reader_own_) set_derived_value<std::string>("solvent",pcm_data);

    // In reader-own mode PCMSolver parses its own input file and none of the keys of
    // this group are consulted, so there is nothing left to validate.
    if (reader_own_) return;

    // A tabulated solvent fixes the permittivity and the probe radius, so derive both
    // rather than leaving the placeholder defaults in the printed group. PCMSolver
    // looks them up in its own database and ignores what we pass unless the solvent is
    // `explicit` -- which is exactly the case where nothing is derived and the values
    // the user gave stand.
    const auto it = solvent_data().find(solvent());
    if (it != solvent_data().end()) {
        set_derived_value<double>("epsilon",it->second.first);
        set_derived_value<double>("probe_radius",it->second.second);
    }

    // Checks that allowed_values cannot express: fixed buffer capacities, and
    // consistency between keys.
    if (solvent().size() > max_solvent_name_length) {
        MADNESS_CHECK(it != solvent_data().end());    // allowed_values guarantees this
        // Qualified: QCCalculationParametersBase::print(header,footer) is in scope here
        // and would swallow these -- the one-argument calls silently, as a header.
        //
        // No World here to gate the output on rank 0, but this aborts the run on every
        // rank anyway, and saying what to use instead is worth the duplication.
        ::madness::print("solvent '" + it->first + "' cannot be passed to PCMSolver:");
        ::madness::print("   PCMInput::solvent holds", max_solvent_name_length,
                "characters, this name has", it->first.size());
        ::madness::print("   use instead, in the pcm group:");
        ::madness::print("      solvent       explicit");
        ::madness::print("      epsilon      ", it->second.first);
        ::madness::print("      probe_radius ", it->second.second);
        MADNESS_EXCEPTION("solvent name too long for the PCMSolver interface",1);
    }
    if (restart_name().size() >= 20)
        MADNESS_EXCEPTION("pcm restart_name must be shorter than 20 characters "
                "(PCMInput::restart_name is char[20])",1);
    if (correction() != 0.0 and solver_type() != "cpcm")
        MADNESS_EXCEPTION("the pcm `correction` keyword applies to solver_type cpcm only",1);
}

} /* namespace madness */


#ifdef MADNESS_HAS_PCM

#include "PCMSolver/pcmsolver.h"
#include "PCMSolver/PCMInput.h"

namespace madness {

namespace detail {

/// to convert the apparent surface charge potential to a real_function_3d
struct asc_potential : public FunctionFunctorInterface<double,3> {

    using FunctionFunctorInterface<double,3>::operator();

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

PCM::PCM(World& world, const Molecule& mol, const PCMParameters& param,
        const bool verbose) : mep_lbl("TotMEP"), asc_lbl("totASC") {
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

    // Either PCMSolver reads its own input file (`pcm_data pcmsolver_reader_own` in the
    // dft group), or -- the normal case -- we hand it the `pcm` group as its host input
    // struct. Validation of that struct happened in PCMParameters::set_derived_values.
    struct PCMInput host_input = param.to_pcmsolver_input();
    const pcmsolver_reader_t pcmsolver_reader =
            param.reader_own() ? PCMSOLVER_READER_OWN : PCMSOLVER_READER_HOST;

    if (verbose and (world.rank()==0)) param.print("pcm","end");

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

namespace detail {

/// copy a parameter into one of PCMInput's fixed char arrays, refusing to overrun it

/// PCMInput's string fields are small fixed arrays (solvent is char[16], radii_set and
/// cavity_type char[8], ...), and the previous code strcpy'd into them unchecked. The
/// allowed-value lists already keep oversized names out, so reaching this is a bug --
/// but a named exception beats a smashed struct if a future PCMSolver shrinks a buffer.
static void copy_field(char* dst, const std::size_t capacity, const std::string& src,
        const char* what) {
    if (src.size() + 1 > capacity) {
        ::madness::print("pcm: value '"+src+"' does not fit PCMSolver's", what,
                "field of", capacity, "bytes");
        MADNESS_EXCEPTION("pcm parameter too long for the PCMSolver interface",1);
    }
    std::strncpy(dst, src.c_str(), capacity);
    dst[capacity-1] = '\0';
}

}

PCMInput PCMParameters::to_pcmsolver_input() const {

    // Start from PCMSolver's own defaults rather than restating them: that keeps the
    // fields this group deliberately does not expose (patch_level, coarsity,
    // min_distance, der_order, equation_type, and the two Green's function types
    // inside_type/outside_type -- `vacuum` is the only string that fits
    // PCMInput::inside_type[7], and on the host-side path `uniformdielectric` is the
    // only outside type whose interface parameters PCMSolver initializes) at whatever
    // upstream considers right, instead of pinning them to values copied here once in
    // 2016.
    PCMInput host_input = pcmsolver_default_input();

    // Lengths and areas go in as Angstrom; PCMSolver converts to Bohr internally.
    // (The coordinates handed to pcmsolver_new are the exception -- those are au.)
    detail::copy_field(host_input.cavity_type, sizeof(host_input.cavity_type),
            cavity_type(), "cavity_type");
    host_input.area = area();
    host_input.scaling = scaling();
    detail::copy_field(host_input.radii_set, sizeof(host_input.radii_set),
            radii_set(), "radii_set");
    detail::copy_field(host_input.restart_name, sizeof(host_input.restart_name),
            restart_name(), "restart_name");
    host_input.min_radius = min_radius();

    detail::copy_field(host_input.solver_type, sizeof(host_input.solver_type),
            solver_type(), "solver_type");
    host_input.correction = correction();
    host_input.probe_radius = probe_radius();

    // "explicit" is PCMSolver's sentinel for "no tabulated solvent": it then takes the
    // permittivity and probe radius from the fields above instead of a database entry.
    // "none" cannot reach here -- PCM is only constructed when pcm_data is set.
    const std::string solv = (solvent()=="none") ? std::string("explicit") : solvent();
    detail::copy_field(host_input.solvent, sizeof(host_input.solvent), solv, "solvent");

    host_input.outside_epsilon = epsilon();

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

PCM::PCM(World& world, const Molecule& mol, const PCMParameters& param,
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
