/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "elements/Source.H"

#include <AMReX_REAL.H>

#ifdef ImpactX_USE_OPENPMD
#   include "elements/diagnostics/openPMD.H"
#   include <openPMD/openPMD.hpp>
namespace io = openPMD;
#endif

#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <string>


namespace impactx::elements
{
    void
    SourceGPU::check_step_selection () const
    {
        if (m_load_step.has_value() && m_load_step_index.has_value()) {
            throw std::invalid_argument(
                "Source: set either load_step or load_step_index, not both."
            );
        }
        if (m_load_step.has_value() && *m_load_step < 0) {
            throw std::invalid_argument(
                "Source: load_step " + std::to_string(*m_load_step) + " is negative, but it is "
                "the ImpactX step stored in the file; use load_step_index to count back from "
                "the last step in the file."
            );
        }
    }

    void
    SourceGPU::operator() (
        ImpactXParticleContainer & pc,
        int,
        int period
    )
    {
        // Check if only active for lattice period zero (0), e.g., in rings
        if (m_active_once && period > 0) {
            return;
        }

        check_step_selection();

#ifdef ImpactX_USE_OPENPMD
        auto series = io::Series(m_series_name, io::Access::READ_ONLY
#   if openPMD_HAVE_MPI==1
            , amrex::ParallelDescriptor::Communicator()
#   endif
        );

        // select the step (openPMD iteration) to read
        auto const & iterations = series.iterations;
        if (iterations.empty()) {
            throw std::runtime_error("Source: no step found in " + m_series_name);
        }

        // the steps stored in the file, for error messages
        // a long list, e.g. one step per turn in a ring, is elided in the middle
        auto const available_steps = [&iterations] () {
            constexpr std::size_t max_listed = 100; //! list at most this many steps
            std::size_t const num_steps = iterations.size();
            bool const elide = num_steps > max_listed;
            std::size_t const num_head = elide ? max_listed / 2 : num_steps;
            std::size_t const num_tail = elide ? max_listed - max_listed / 2 : 0;

            std::string list;
            auto append_step = [&list] (auto const & iteration_entry) {
                if (!list.empty()) { list += ", "; }
                list += std::to_string(iteration_entry->first);
            };

            auto it = iterations.begin();
            for (std::size_t i = 0; i < num_head; ++i, ++it) { append_step(it); }
            if (elide) {
                list += ", ...";
                std::advance(it, static_cast<long>(num_steps - num_head - num_tail));
                for (std::size_t i = 0; i < num_tail; ++i, ++it) { append_step(it); }
            }

            return std::to_string(num_steps) +
                   (num_steps == 1 ? " available step: " : " available steps: ") + list;
        };

        io::Iteration::IterationIndex_t read_iteration = 0;
        if (m_load_step.has_value())
        {
            // an ImpactX step, which the beam monitor wrote as the openPMD iteration
            auto const requested = static_cast<io::Iteration::IterationIndex_t>(*m_load_step);
            if (!iterations.contains(requested)) {
                throw std::runtime_error(
                    "Source: load_step " + std::to_string(*m_load_step) + " not found in " +
                    m_series_name + " (" + available_steps() + ")"
                );
            }
            read_iteration = requested;
        }
        else
        {
            // a position in the file: 0 is the first step, -1 the last, as in Python
            constexpr int last_step_index = -1; //! read the last step in the file by default
            int const step_index = m_load_step_index.value_or(last_step_index);

            auto const num_steps = static_cast<long>(iterations.size());
            long const index = step_index < 0
                ? num_steps + static_cast<long>(step_index)
                : static_cast<long>(step_index);
            if (index < 0 || index >= num_steps) {
                throw std::runtime_error(
                    "Source: load_step_index " + std::to_string(step_index) +
                    " is out of range in " + m_series_name + " (" + available_steps() + ")"
                );
            }
            read_iteration = std::next(iterations.begin(), index)->first;
        }

        io::Iteration iteration = series.iterations[read_iteration];

        // TODO: later we can make the particle species name an option
        std::string const species_name = "beam";
        io::ParticleSpecies beam = iteration.particles[species_name];
        // TODO: later we can make the bunch charge an option (i.e., allow rescaling a distribution)
        // amrex::ParticleReal bunch_charge = beam.getAttribute("charge_C").get<amrex::ParticleReal>();

        auto const scalar = openPMD::RecordComponent::SCALAR;
        auto const getComponentRecord = [&beam](std::string comp_name) {
            return diagnostics::detail::get_component_record(beam, std::move(comp_name));
        };

        int const npart = beam["id"][scalar].getExtent()[0];  // how many particles to read total

        // restore the reference particle from the species metadata
        if (m_load_ref_particle)
        {
            RefPart & ref = pc.GetRefParticle();
            auto const read_attr = [&beam, this] (std::string const & attr_name) {
                if (!beam.containsAttribute(attr_name)) {
                    throw std::runtime_error(
                        "Source: attribute '" + attr_name + "' not found in " + m_series_name +
                        " (load_ref_particle requires a file written by an ImpactX beam_monitor;"
                        " set load_ref_particle=false to keep the current reference particle)"
                    );
                }
                return beam.getAttribute(attr_name).get<amrex::ParticleReal>();
            };

            ref.s = read_attr("s_ref");
            ref.x = read_attr("x_ref");
            ref.y = read_attr("y_ref");
            ref.z = read_attr("z_ref");
            ref.t = read_attr("t_ref");
            ref.px = read_attr("px_ref");
            ref.py = read_attr("py_ref");
            ref.pz = read_attr("pz_ref");
            ref.pt = read_attr("pt_ref");
            ref.gyromagnetic_anomaly = read_attr("gyromagnetic_anomaly_ref");
            ref.mass = read_attr("mass_ref");
            ref.charge = read_attr("charge_ref");
            // ref.sedge: could be set to ref.s, but should be set anyway in tracking loops on element entry.
        }

        // read the particles

        // Logic: We initialize 1/Nth of particles, independent of their
        // position, per MPI rank. We then measure the distribution's spatial
        // extent, create a grid, resize it to fit the beam, and then
        // redistribute particles so that they reside on the correct MPI rank.
        int const myproc = amrex::ParallelDescriptor::MyProc();
        int const nprocs = amrex::ParallelDescriptor::NProcs();
        int const navg = npart / nprocs;  // note: integer division
        int const nleft = npart - navg * nprocs;
        std::uint64_t const npart_this_proc = (myproc < nleft) ? navg+1 : navg;  // add 1 to each proc until distributed
        std::uint64_t npart_before_this_proc = (myproc < nleft) ? (navg+1) * myproc : navg * myproc + nleft;

        // alloc data for particle attributes
        std::map<std::string, amrex::Gpu::PinnedVector<amrex::ParticleReal>> pinned_SoA;
        std::vector<std::string> real_soa_names = pc.GetRealSoANames();
        for (auto real_idx = 0; real_idx < pc.NumRealComps(); real_idx++) {
            pinned_SoA[real_soa_names.at(real_idx)].resize(npart_this_proc);
        }

        // read from file
        // idcpu: TODO
        //amrex::Gpu::PinnedVector<std::uint64_t> pinned_idcpu(npart_this_proc);
        //beam["id"][scalar].loadChunkRaw(pinned_idcpu.data(), {npart_before_this_proc}, {npart_this_proc});
        // SoA: Real
        {
            for (auto real_idx = 0; real_idx < pc.NumRealComps(); real_idx++) {
                auto const component_name = real_soa_names.at(real_idx);
                getComponentRecord(component_name).loadChunkRaw(
                    pinned_SoA[component_name].data(),
                    {npart_before_this_proc},
                    {npart_this_proc}
                );
            }
        }
        // SoA: Int
        std::vector<std::string> int_soa_names = pc.GetIntSoANames();
        static_assert(IntSoA::nattribs == 0); // not yet used
        if (!int_soa_names.empty())
            throw std::runtime_error("BeamMonitor: int_soa_names output not yet implemented!");

        series.flush();
        series.close();

        // copy to device
        std::map<std::string, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_SoA;
        for (auto real_idx = 0; real_idx < pc.NumRealComps(); real_idx++) {
            auto const component_name = real_soa_names.at(real_idx);
            d_SoA[component_name].resize(npart_this_proc);
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, pinned_SoA[component_name].begin(), pinned_SoA[component_name].end(), d_SoA[component_name].begin());
        }

        // finalize distributions and deallocate temporary device global memory
        amrex::Gpu::streamSynchronize();

        // TODO: at this point, we ignore the "id" and "qm" in the file. Could be improved

        pc.AddNParticles(d_SoA["position_x"], d_SoA["position_y"], d_SoA["position_t"],
                         d_SoA["momentum_x"], d_SoA["momentum_y"], d_SoA["momentum_t"],
                         pc.GetRefParticle().qm_ratio_SI(),
                         std::nullopt,
                         d_SoA["weighting"],
                         d_SoA["spin_x"], d_SoA["spin_y"], d_SoA["spin_z"]);

#else  // ImpactX_USE_OPENPMD
        amrex::ignore_unused(pc);
        throw std::runtime_error("BeamMonitor: openPMD not compiled");
#endif  // ImpactX_USE_OPENPMD
    }

} // namespace impactx
