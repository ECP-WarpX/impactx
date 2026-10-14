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

#include "particles/transformation/ToFixedS.H"

#include <ablastr/constant.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_REAL.H>

#ifdef ImpactX_USE_OPENPMD
#   include "elements/diagnostics/openPMD.H"
#   include <openPMD/openPMD.hpp>
namespace io = openPMD;
#endif

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>


namespace impactx::elements
{
#ifdef ImpactX_USE_OPENPMD
namespace
{
    //! particle attributes on the host, by name
    using PinnedSoA = std::map<std::string, amrex::Gpu::PinnedVector<amrex::ParticleReal>>;

    /** A chunk of a floating point openPMD record component
     *
     * Read the values after the series was flushed.
     */
    struct RealChunk
    {
        /** Request a chunk of a record component.
         *
         * @param[in] component openPMD record component to read
         * @param[in] offset index of the first value to read
         * @param[in] extent number of values to read
         * @param[in] where name of the component and its file, for error messages
         * @throws std::runtime_error if the component is not stored as float or double
         */
        RealChunk (
            io::RecordComponent component,
            std::uint64_t offset,
            std::uint64_t extent,
            std::string const & where
        )
        : scale(component.unitSI())
        {
            // an MPI rank without particles
            if (extent == 0) { return; }

            switch (component.getDatatype())
            {
                case io::Datatype::FLOAT:
                    m_float = component.loadChunk<float>({offset}, {extent});
                    break;
                case io::Datatype::DOUBLE:
                    m_double = component.loadChunk<double>({offset}, {extent});
                    break;
                default:
                    throw std::runtime_error(
                        "Source: " + where + " is not stored as float or double.");
            }
        }

        /** The value at an index of the chunk, times scale */
        double operator[] (std::uint64_t i) const
        {
            return scale * (m_double ? m_double.get()[i] : static_cast<double>(m_float.get()[i]));
        }

        double scale; //! factor on the stored values: starts as unitSI, which converts them to SI

    private:
        std::shared_ptr<float> m_float; //! the values, if the component is stored as float
        std::shared_ptr<double> m_double; //! the values, if the component is stored as double
    };

    /** Read particles at fixed t and transform them to fixed s.
     *
     * Positions and momenta are read as in the openPMD standard: after multiplication with
     * unitSI, the positions plus their positionOffset are in meters and the momenta are in
     * kg*m/s. Both are taken in the frame of the reference particle, whose z axis points
     * along the reference momentum. Each particle then drifts to the position s of the
     * reference particle.
     *
     * @param[in,out] series openPMD series, which is flushed to read the particles
     * @param[in] beam openPMD particle species in the series
     * @param[in] series_name path of the series, for error messages
     * @param[in] ref reference particle
     * @param[in] offset index of the first particle to read on this MPI rank
     * @param[in] extent number of particles to read on this MPI rank
     * @param[out] pinned_SoA particles at fixed s: positions, momenta, weighting and spin
     * @throws std::runtime_error if the reference particle has no momentum, if positions or
     *         momenta lack their openPMD unit or if particles do not move forward along z
     */
    void
    read_fixed_t (
        io::Series & series,
        io::ParticleSpecies & beam,
        std::string const & series_name,
        RefPart const & ref,
        std::uint64_t offset,
        std::uint64_t extent,
        PinnedSoA & pinned_SoA
    )
    {
        auto const scalar = io::RecordComponent::SCALAR;
        std::array<std::string, 3> const axes = {"x", "y", "z"};

        // momentum of the reference particle, which normalizes the momenta
        double const p_ref_SI =
            static_cast<double>(ref.beta_gamma()) *
            static_cast<double>(ref.mass) *
            ablastr::constant::SI::c_v<double>;
        if (!(p_ref_SI > 0.0)) {
            throw std::runtime_error(
                "Source: the reference particle has no momentum, but the particles at fixed t "
                "in " + series_name + " are transformed to fixed s relative to it. "
                "Set its species and energy or use load_ref_particle.");
        }

        // positions and momenta carry the units of the openPMD standard
        //   powers of length, mass, time, current, temperature, amount and luminous intensity
        using UnitDimension = std::array<double, 7>;
        auto const check_unit = [&beam, &series_name] (
            std::string const & record_name,
            UnitDimension const & unit,
            std::string const & unit_name
        ) {
            if (beam[record_name].unitDimension() != unit) {
                throw std::runtime_error(
                    "Source: " + record_name + " in " + series_name + " is not " + unit_name +
                    " (openPMD unitDimension). Particles at fixed t (coord_system=t) are read "
                    "as in the openPMD standard, in SI after multiplication with unitSI.");
            }
        };
        check_unit("position", {1., 0., 0., 0., 0., 0., 0.}, "a length");
        check_unit("momentum", {1., 1., -1., 0., 0., 0., 0.}, "a momentum");
        bool const has_position_offset = beam.contains("positionOffset");
        if (has_position_offset) {
            check_unit("positionOffset", {1., 0., 0., 0., 0., 0., 0.}, "a length");
        }

        // momenta of macroparticles are divided by a power of the weighting
        io::Record momentum_record = beam["momentum"];
        bool const macro_weighted =
            momentum_record.containsAttribute("macroWeighted") &&
            momentum_record.getAttribute("macroWeighted").get<std::uint32_t>() != 0u;
        double const weighting_power = momentum_record.containsAttribute("weightingPower")
            ? momentum_record.getAttribute("weightingPower").get<double>()
            : 0.0;
        bool const per_macroparticle = macro_weighted && weighting_power != 0.0;

        // request the particle attributes
        auto const request = [&beam, &series_name, offset, extent] (
            std::string const & record_name,
            std::string const & component_name
        ) {
            std::string const where = component_name == io::RecordComponent::SCALAR
                ? record_name + " in " + series_name
                : record_name + "/" + component_name + " in " + series_name;
            return RealChunk(beam[record_name][component_name], offset, extent, where);
        };

        std::vector<RealChunk> position;
        std::vector<RealChunk> position_offset;
        std::vector<RealChunk> momentum;
        for (auto const & axis : axes)
        {
            position.push_back(request("position", axis));
            if (has_position_offset) {
                position_offset.push_back(request("positionOffset", axis));
            }
            // in units of the reference momentum
            momentum.push_back(request("momentum", axis));
            momentum.back().scale /= p_ref_SI;
        }

        if (!beam.contains("weighting")) {
            throw std::runtime_error("Source: no weighting found in " + series_name + ".");
        }
        RealChunk const weighting = request("weighting", scalar);

        // spin is optional: other codes usually do not write it
        bool const has_spin = beam.contains("spin");
        std::vector<RealChunk> spin;
        if (has_spin) {
            for (auto const & axis : axes) { spin.push_back(request("spin", axis)); }
        }

        series.flush();

        // particles at fixed s
        amrex::ParticleReal * const AMREX_RESTRICT part_x = pinned_SoA["position_x"].data();
        amrex::ParticleReal * const AMREX_RESTRICT part_y = pinned_SoA["position_y"].data();
        amrex::ParticleReal * const AMREX_RESTRICT part_t = pinned_SoA["position_t"].data();
        amrex::ParticleReal * const AMREX_RESTRICT part_px = pinned_SoA["momentum_x"].data();
        amrex::ParticleReal * const AMREX_RESTRICT part_py = pinned_SoA["momentum_y"].data();
        amrex::ParticleReal * const AMREX_RESTRICT part_pt = pinned_SoA["momentum_t"].data();
        amrex::ParticleReal * const AMREX_RESTRICT part_w = pinned_SoA["weighting"].data();
        amrex::ParticleReal * const AMREX_RESTRICT part_sx = pinned_SoA["spin_x"].data();
        amrex::ParticleReal * const AMREX_RESTRICT part_sy = pinned_SoA["spin_y"].data();
        amrex::ParticleReal * const AMREX_RESTRICT part_sz = pinned_SoA["spin_z"].data();

        // position of the reference particle
        std::array<double, 3> const ref_position = {
            static_cast<double>(ref.x),
            static_cast<double>(ref.y),
            static_cast<double>(ref.z)
        };

        impactx::transformation::ToFixedS const to_s(ref.beta_gamma());

        // particles that do not move forward along z never pass the reference position s
        amrex::Long num_backward = 0;

        for (std::uint64_t i = 0; i < extent; ++i)
        {
            double const w = weighting[i];
            double const per_real_particle =
                per_macroparticle ? std::pow(w, -weighting_power) : 1.0;

            // relative to the reference particle: first subtract the two large numbers
            std::array<double, 3> pos{};
            std::array<double, 3> mom{};
            for (std::size_t d = 0; d < 3; ++d)
            {
                double const coarse = has_position_offset ? position_offset[d][i] : 0.0;
                pos[d] = position[d][i] + (coarse - ref_position[d]);
                mom[d] = momentum[d][i] * per_real_particle;
            }
            if (!(mom[2] > 0.0)) { ++num_backward; }
            // the longitudinal momentum is a deviation from the reference momentum
            mom[2] -= 1.0;

            auto x = static_cast<amrex::ParticleReal>(pos[0]);
            auto y = static_cast<amrex::ParticleReal>(pos[1]);
            auto z = static_cast<amrex::ParticleReal>(pos[2]);
            auto px = static_cast<amrex::ParticleReal>(mom[0]);
            auto py = static_cast<amrex::ParticleReal>(mom[1]);
            auto pz = static_cast<amrex::ParticleReal>(mom[2]);

            // z and pz become t and pt
            to_s(x, y, z, px, py, pz);

            part_x[i] = x;
            part_y[i] = y;
            part_t[i] = z;
            part_px[i] = px;
            part_py[i] = py;
            part_pt[i] = pz;
            part_w[i] = static_cast<amrex::ParticleReal>(w);
            part_sx[i] = static_cast<amrex::ParticleReal>(has_spin ? spin[0][i] : 0.0);
            part_sy[i] = static_cast<amrex::ParticleReal>(has_spin ? spin[1][i] : 0.0);
            part_sz[i] = static_cast<amrex::ParticleReal>(has_spin ? spin[2][i] : 0.0);
        }

        // all MPI ranks need to agree on the error
        amrex::ParallelDescriptor::ReduceLongSum(num_backward);
        if (num_backward > 0) {
            throw std::runtime_error(
                "Source: " + std::to_string(num_backward) + " particles in " + series_name +
                " do not move forward along z. They cannot be transformed from fixed t to "
                "fixed s.");
        }
    }
} // namespace
#endif // ImpactX_USE_OPENPMD

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

        // the particles are added at fixed s
        if (pc.GetCoordSystem() != CoordSystem::s) {
            throw std::runtime_error(
                "Source: the particles in the beam must be at fixed s to add particles.");
        }

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

        // the longitudinal components tell the coordinate system of the file:
        // t at fixed s, z at fixed t
        auto const has_longitudinal = [&beam] (std::string const & component_name) {
            return beam.contains("position") && beam["position"].contains(component_name) &&
                   beam.contains("momentum") && beam["momentum"].contains(component_name);
        };
        bool const file_at_fixed_s = has_longitudinal("t");
        bool const file_at_fixed_t = has_longitudinal("z");
        if (m_coord_system == CoordSystem::s && !file_at_fixed_s) {
            throw std::runtime_error(
                "Source: no particles at fixed s (position t, momentum t) found in " +
                m_series_name + "." +
                (file_at_fixed_t ? " It holds particles at fixed t: set coord_system to t." : ""));
        }
        if (m_coord_system == CoordSystem::t && !file_at_fixed_t) {
            throw std::runtime_error(
                "Source: no particles at fixed t (position z, momentum z) found in " +
                m_series_name + "." +
                (file_at_fixed_s ? " It holds particles at fixed s: set coord_system to s." : ""));
        }

        // how many particles to read total
        int const npart = m_coord_system == CoordSystem::s
            ? beam["id"][scalar].getExtent()[0]
            : beam["position"]["x"].getExtent()[0];

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
        PinnedSoA pinned_SoA;
        std::vector<std::string> real_soa_names = pc.GetRealSoANames();
        for (auto real_idx = 0; real_idx < pc.NumRealComps(); real_idx++) {
            pinned_SoA[real_soa_names.at(real_idx)].resize(npart_this_proc);
        }

        // read from file
        // idcpu: TODO
        //amrex::Gpu::PinnedVector<std::uint64_t> pinned_idcpu(npart_this_proc);
        //beam["id"][scalar].loadChunkRaw(pinned_idcpu.data(), {npart_before_this_proc}, {npart_this_proc});
        // SoA: Real
        if (m_coord_system == CoordSystem::t)
        {
            read_fixed_t(series, beam, m_series_name, pc.GetRefParticle(),
                         npart_before_this_proc, npart_this_proc, pinned_SoA);
        }
        else
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
