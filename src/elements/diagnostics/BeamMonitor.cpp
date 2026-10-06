/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */

#include "BeamMonitor.H"

#include "ImpactXVersion.H"
#include "diagnostics/FilePrefix.H"
#include "diagnostics/ReducedBeamCharacteristics.H"
#include "particles/ImpactXParticleContainer.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

#ifdef ImpactX_USE_OPENPMD
#   include "elements/diagnostics/openPMD.H"
#   include <openPMD/openPMD.hpp>
namespace io = openPMD;
#endif

#include <filesystem>
#include <fstream>  // for std::ofstream
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace impactx::elements::diagnostics
{
namespace detail {
    ImpactXParticleCounter::ImpactXParticleCounter (ParticleContainer & pc)
    {
        m_MPISize = amrex::ParallelDescriptor::NProcs();
        m_MPIRank = amrex::ParallelDescriptor::MyProc();

        m_ParticleCounterByLevel.resize(pc.finestLevel()+1);
        m_ParticleOffsetAtRank.resize(pc.finestLevel()+1);
        m_ParticleSizeAtRank.resize(pc.finestLevel()+1);

        for (auto currentLevel = 0; currentLevel <= pc.finestLevel(); currentLevel++)
        {
            amrex::Long numParticles = 0;  // numParticles on this MPI rank

            for (ParticleIter pti(pc, currentLevel); pti.isValid(); ++pti) {
                amrex::Long const numParticleOnTile = pti.numParticles();
                numParticles += numParticleOnTile;
            }

            unsigned long long offset=0; // offset of this level
            unsigned long long sum=0; // numParticles in this level (sum from all processors)

            GetParticleOffsetOfProcessor(numParticles, offset,  sum);

            m_ParticleCounterByLevel[currentLevel] = sum;
            m_ParticleOffsetAtRank[currentLevel] = offset;
            m_ParticleSizeAtRank[currentLevel] = numParticles;

            // adjust offset, it should be numbered after particles from previous levels
            for (auto lv=0; lv<currentLevel; lv++)
                m_ParticleOffsetAtRank[currentLevel] += m_ParticleCounterByLevel[lv];

            m_Total += sum;
        }
    }

    // get the offset in the overall particle id collection
    //
    // note: this is a MPI-collective operation
    //
    // input: num of particles  of from each   processor
    //
    // output:
    //     offset within <all> the particles in the comm
    //     sum of all particles in the comm
    //
    void
    ImpactXParticleCounter::GetParticleOffsetOfProcessor (
            amrex::Long const & numParticles,
            unsigned long long& offset,
            unsigned long long& sum
    ) const
    {
        offset = 0;
#if defined(AMREX_USE_MPI)
        std::vector<amrex::Long> result(m_MPISize, 0);
    amrex::ParallelGather::Gather (numParticles, result.data(), -1, amrex::ParallelDescriptor::Communicator());

    sum = 0;
    int const num_results = result.size();
    for (int i=0; i<num_results; i++) {
        sum += result[i];
        if (i<m_MPIRank)
            offset += result[i];
    }
#else
        sum = numParticles;
#endif
    }

} // namespace detail

    void BeamMonitor::finalize ()
    {
#ifdef ImpactX_USE_OPENPMD
        // close shared series alias
        if (m_series.has_value())
        {
            auto series = std::any_cast<io::Series>(m_series);
            series.close();
            m_series.reset();
        }

        // remove from unique series map
        if (!m_series_key.empty() && m_unique_series.count(m_series_key) != 0u) {
            m_unique_series.erase(m_series_key);
            m_series_key.clear();
        }
#endif // ImpactX_USE_OPENPMD
    }

    BeamMonitor::BeamMonitor (
        std::string series_name,
        std::string backend,
        std::string encoding,
        int period_sample_intervals,
        bool particles,
        bool beam_moments
    ) :
        m_series_name(std::move(series_name)),
        m_OpenPMDFileType(std::move(backend)),
        m_encoding(std::move(encoding)),
        m_period_sample_intervals(period_sample_intervals),
        m_particles(particles),
        m_beam_moments(beam_moments)
    {
        if (period_sample_intervals < 1) {
            throw std::invalid_argument("BeamMonitor: period_sample_intervals must be >= 1");
        }
    }

    void BeamMonitor::open ()
    {
#ifdef ImpactX_USE_OPENPMD
        // pick first available backend if default is chosen
        if (m_OpenPMDFileType == "default")
#   if openPMD_HAVE_ADIOS2==1
        m_OpenPMDFileType = "bp4";
#   elif openPMD_HAVE_ADIOS1==1
        m_OpenPMDFileType = "bp";  // bp3
#   elif openPMD_HAVE_HDF5==1
        m_OpenPMDFileType = "h5";
#   else
        m_OpenPMDFileType = "json";
#   endif

        // encoding of iterations in the series
        openPMD::IterationEncoding series_encoding = openPMD::IterationEncoding::groupBased;
        if ("v" == m_encoding)
            series_encoding = openPMD::IterationEncoding::variableBased;
        else if ("g" == m_encoding)
            series_encoding = openPMD::IterationEncoding::groupBased;
        else if ("f" == m_encoding)
            series_encoding = openPMD::IterationEncoding::fileBased;

        // BP5 does not support groupBased (metadata explosion)
        if ((m_OpenPMDFileType == "bp5" || m_OpenPMDFileType == "bp") &&
            (series_encoding == openPMD::IterationEncoding::groupBased))
        {
            throw std::runtime_error("BeamMonitor: groupBased encoding not supported for BP5.");
        }

        amrex::ParmParse pp_diag("diag");
        // legacy options from other diagnostics
        pp_diag.queryAddWithParser("file_min_digits", m_file_min_digits);

        std::filesystem::path filepath = impactx::diagnostics::FilePrefixPath("openPMD");
        std::string filename = m_series_name;

        if (series_encoding == openPMD::IterationEncoding::fileBased)
        {
            std::string const fileSuffix =
                std::string("_%0") + std::to_string(m_file_min_digits) + std::string("T");
            filename.append(fileSuffix);
        }
        filename.append(".").append(m_OpenPMDFileType);

        std::string series_path = (filepath / filename).string();

        // transform paths for Windows
#   ifdef _WIN32
        series_path = openPMD::auxiliary::replace_all(series_path, "/", "\\");
#   endif
        m_series_key = series_path;

        // Ensure m_series is the same for the same full output path.
        if (m_unique_series.count(m_series_key) == 0u) {
            auto series = io::Series(m_series_key, io::Access::CREATE
#   if openPMD_HAVE_MPI==1
                , amrex::ParallelDescriptor::Communicator()
#   endif
            );
            series.setSoftware("ImpactX", IMPACTX_VERSION);
            series.setIterationEncoding( series_encoding );
            m_series = series;
            m_unique_series[m_series_key] = series;

            // create a little helper file for ParaView 5.9+
            if (amrex::ParallelDescriptor::IOProcessor())
            {
                std::filesystem::create_directories(filepath);
                std::ofstream pv_helper_file((filepath / "paraview.pmd").string());
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pv_helper_file.is_open(), "Could not open paraview.pmd file.");
                pv_helper_file << filename << "\n";
                pv_helper_file.close();
            }
        }
        else {
            m_series = m_unique_series[m_series_key];
        }
#else
        amrex::AllPrint() << "Warning: openPMD output requested but not compiled for series=" << m_series_name << "\n";
#endif
    }

    void BeamMonitor::prepare (
        PinnedContainer & pc,
        std::vector<std::string> const & real_soa_names,
        std::vector<std::string> const & int_soa_names,
        RefPart const & ref_part,
        int step
    ) {
#ifdef ImpactX_USE_OPENPMD
        m_step = step;

        // series & iteration
        //   note: the ADIOS2 BP4 engine drops attribute-only iterations with
        //         the steps-based writeIterations() API, thus use random
        //         access when no particle data sets are written
        auto series = std::any_cast<io::Series>(m_series);
        io::Iteration iteration = use_io_steps() ?
                                  io::Iteration(series.writeIterations()[m_step]) :
                                  series.iterations[m_step];
        io::ParticleSpecies beam = iteration.particles["beam"];

        // openPMD 1.* needs "seconds" here, but we fake it as "s"
        iteration.setTime(ref_part.s);

        // reference particle information
        beam.setAttribute( "beta_ref", ref_part.beta() );
        beam.setAttribute( "gamma_ref", ref_part.gamma() );
        beam.setAttribute( "beta_gamma_ref", ref_part.beta_gamma() );
        beam.setAttribute( "s_ref", ref_part.s );
        beam.setAttribute( "x_ref", ref_part.x );
        beam.setAttribute( "y_ref", ref_part.y );
        beam.setAttribute( "z_ref", ref_part.z );
        beam.setAttribute( "t_ref", ref_part.t );
        beam.setAttribute( "px_ref", ref_part.px );
        beam.setAttribute( "py_ref", ref_part.py );
        beam.setAttribute( "pz_ref", ref_part.pz );
        beam.setAttribute( "pt_ref", ref_part.pt );
        beam.setAttribute( "gyromagnetic_anomaly_ref", ref_part.gyromagnetic_anomaly );
        beam.setAttribute( "mass_ref", ref_part.mass );
        beam.setAttribute( "charge_ref", ref_part.charge );

        // total particle bunch information
        //   @see impactx::diagnostics::reduced_beam_characteristics
        for (const auto &kv : m_rbc) {
            beam.setAttribute(kv.first, kv.second);
        }

        // individual particle data
        if (m_particles)
        {
            // calculate & update particle offset in MPI-global particle array, per level
            auto const num_levels = pc.finestLevel() + 1;
            m_offset = std::vector<uint64_t>(num_levels);
            auto counter = detail::ImpactXParticleCounter(pc);
            auto const np = counter.GetTotalNumParticles();
            for (auto currentLevel = 0; currentLevel < num_levels; currentLevel++) {
                m_offset.at(currentLevel) = static_cast<uint64_t>( counter.m_ParticleOffsetAtRank[currentLevel] );
            }

            // helpers to parse strings to openPMD
            auto const scalar = openPMD::RecordComponent::SCALAR;
            auto const getComponentRecord = [&beam](std::string comp_name) {
                return detail::get_component_record(beam, std::move(comp_name));
            };

            // define data set and metadata
            io::Datatype const dtype_fl = io::determineDatatype<amrex::ParticleReal>();
            io::Datatype const dtype_ui = io::determineDatatype<uint64_t>();
            auto d_fl = io::Dataset(dtype_fl, {np});
            auto d_ui = io::Dataset(dtype_ui, {np});

            // openPMD coarse position: for global coordinates
            {

                beam["positionOffset"]["x"].resetDataset(d_fl);
                beam["positionOffset"]["x"].makeConstant(ref_part.x);
                beam["positionOffset"]["y"].resetDataset(d_fl);
                beam["positionOffset"]["y"].makeConstant(ref_part.y);
                beam["positionOffset"]["t"].resetDataset(d_fl);
                beam["positionOffset"]["t"].makeConstant(ref_part.t);
            }

            // unique, global particle index
            beam["id"][scalar].resetDataset(d_ui);

            // SoA: Real
            {
                for (auto real_idx = 0; real_idx < pc.NumRealComps(); real_idx++) {
                    auto const component_name = real_soa_names.at(real_idx);
                    getComponentRecord(component_name).resetDataset(d_fl);
                }
            }
            // SoA: Int
            static_assert(IntSoA::nattribs == 0); // not yet used
            if (!int_soa_names.empty())
                throw std::runtime_error("BeamMonitor: int_soa_names output not yet implemented!");
        }
        else
        {
            // zero-extent constant record: openPMD-api 0.17+ readers reject
            // particle species without any component extent by default
            auto const scalar = openPMD::RecordComponent::SCALAR;
            io::Datatype const dtype_fl = io::determineDatatype<amrex::ParticleReal>();
            auto d_none = io::Dataset(dtype_fl, {0});

            beam["empty"][scalar].resetDataset(d_none);
            beam["empty"][scalar].makeConstant(amrex::ParticleReal(0.0));
        }
#else
        amrex::ignore_unused(pc, real_soa_names, int_soa_names, ref_part, step);
#endif // ImpactX_USE_OPENPMD
    }

    void
    BeamMonitor::operator() (
        ImpactXParticleContainer & pc,
        int step,
        int period
    )
    {
        // filter out this turn?
        if (period % m_period_sample_intervals != 0)
            return;

        this->open();

#ifdef ImpactX_USE_OPENPMD
        std::string profile_name = "impactx::push::" + std::string(BeamMonitor::type);
        BL_PROFILE(profile_name);

        // the output schema describes particles at fixed s
        if (pc.GetCoordSystem() != CoordSystem::s)
        {
            throw std::runtime_error(
                "BeamMonitor (" + m_series_name + "): the particles must be at fixed s.");
        }

        // preparing to access reference particle data: RefPart
        RefPart & ref_part = pc.GetRefParticle();

        // optional: add and calculate additional particle properties
        if (m_particles)
        {
            add_optional_properties(m_series_name, pc);
        }
        else
        {
            // H & I are per-particle columns and cannot be written without particle output
            bool add_nll_invariants = false;
            amrex::ParmParse pp_element(m_series_name);
            pp_element.queryAdd("nonlinear_lens_invariants", add_nll_invariants);
            if (add_nll_invariants)
            {
                throw std::runtime_error(
                    "BeamMonitor (" + m_series_name +
                    "): nonlinear_lens_invariants=true requires particles=true");
            }
        }

        // optional: calculate total particle bunch information
        m_rbc.clear();
        if (m_beam_moments)
        {
            m_rbc = impactx::diagnostics::reduced_beam_characteristics(pc);
            m_rbc["period"] = period;  // current period (turn or cycle)
        }

        // component names
        std::vector<std::string> real_soa_names = pc.GetRealSoANames();
        std::vector<std::string> int_soa_names = pc.GetIntSoANames();

        // pinned memory copy
        PinnedContainer pinned_pc = pc.make_alike<amrex::PolymorphicArenaAllocator>();
        pinned_pc.SetArena(amrex::The_Pinned_Arena());
        if (m_particles)
        {
            pinned_pc.copyParticles(pc, true);  // no filtering
        }

        // TODO: filtering
        /*
        using SrcData = WarpXParticleContainer::ParticleTileType::ConstParticleTileDataType;
        tmp.copyParticles(*pc,
                          [=] AMREX_GPU_HOST_DEVICE (const SrcData& src, int ip, const amrex::RandomEngine& engine)
                          {
                              const SuperParticleType& p = src.getSuperParticle(ip);
                              return random_filter(p, engine) * uniform_filter(p, engine)
                                     * parser_filter(p, engine) * geometry_filter(p, engine);
                          }, true);
        */

        // prepare element access & write reference particle
        this->prepare(pinned_pc, real_soa_names, int_soa_names, ref_part, step);

        if (m_particles)
        {
            // loop over refinement levels
            int const nLevel = pinned_pc.finestLevel();
            for (int lev = 0; lev <= nLevel; ++lev)
            {
                // loop over all particle boxes
                //using ParIt = ImpactXParticleContainer::iterator;
                using ParIt = PinnedContainer::ParIterType;
                // note: openPMD-api is not thread-safe, so do not run OMP parallel here
                for (ParIt pti(pinned_pc, lev); pti.isValid(); ++pti) {
                    // write beam particles relative to reference particle
                    this->operator()(pti, real_soa_names, int_soa_names, ref_part);
                } // end loop over all particle boxes
            } // end mesh-refinement level loop
        }

        //   note: as in prepare(), use random access when no particle data
        //         sets are written
        auto series = std::any_cast<io::Series>(m_series);
        io::Iteration iteration = use_io_steps() ?
                                  io::Iteration(series.writeIterations()[m_step]) :
                                  series.iterations[m_step];

        // close iteration: writes all particle tiles (MPI-collective) while the
        // pinned particle copy is still alive
        iteration.close();
#else
        amrex::ignore_unused(pc, step, period);
#endif // ImpactX_USE_OPENPMD
    }

    void
    BeamMonitor::operator() (
        PinnedContainer::ParIterType & pti,
        std::vector<std::string> const & real_soa_names,
        std::vector<std::string> const & int_soa_names,
        RefPart const & ref_part
    )
    {
#ifdef ImpactX_USE_OPENPMD
        int const currentLevel = pti.GetLevel();

        auto & offset = m_offset.at(currentLevel); // ...

        // series & iteration
        auto series = std::any_cast<io::Series>(m_series);
        io::WriteIterations iterations = series.writeIterations();
        io::Iteration iteration = iterations[m_step];

        // writing
        io::ParticleSpecies beam = iteration.particles["beam"];

        auto const numParticleOnTile = pti.numParticles();
        uint64_t const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

        // Do not call storeChunk() with zero-sized particle tiles:
        //   https://github.com/openPMD/openPMD-api/issues/1147
        //if (numParticleOnTile == 0) continue;

        auto const scalar = openPMD::RecordComponent::SCALAR;
        auto const getComponentRecord = [&beam](std::string comp_name) {
            return detail::get_component_record(beam, std::move(comp_name));
        };

        // SoA
        auto const& soa = pti.GetStructOfArrays();
        //   particle id arrays
        {
            beam["id"][scalar].storeChunkRaw(soa.GetIdCPUData().data(), {offset}, {numParticleOnTile64});
        }
        //   SoA floating point (ParticleReal) properties
        {
            for (auto real_idx=0; real_idx < soa.NumRealComps(); real_idx++) {
                auto const component_name = real_soa_names.at(real_idx);
                getComponentRecord(component_name).storeChunkRaw(
                soa.GetRealData(real_idx).data(), {offset}, {numParticleOnTile64});
            }
        }
        //   SoA integer (int) properties (not yet used)
        {
            static_assert(IntSoA::nattribs == 0); // not yet used
            if (!int_soa_names.empty())
                throw std::runtime_error("BeamMonitor: int_soa_names output not yet implemented!");
            /*
            // comment this in once IntSoA::nattribs is > 0

            std::copy(IntSoA::names_s.begin(), IntSoA::names_s.end(), int_soa_names.begin());

            for (auto int_idx=0; int_idx < RealSoA::nattribs; int_idx++) {
                auto const component_name = int_soa_names.at(int_idx);
                getComponentRecord(component_name).storeChunkRaw(
                    soa.GetIntData(int_idx).data(), {offset}, {numParticleOnTile64});
            }
            */
        }

        // TODO
        amrex::ignore_unused(ref_part);

        // needs to be higher for next pti; must be reset for next step via prepare
        offset += numParticleOnTile64;
#else
        amrex::ignore_unused(pti, real_soa_names, int_soa_names, ref_part);
#endif   // ImpactX_USE_OPENPMD
    }

} // namespace impactx::diagnostics
