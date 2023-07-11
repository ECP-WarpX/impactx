/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */

#include "openPMD.H"
#include "ImpactXVersion.H"
#include "particles/ImpactXParticleContainer.H"

#include <ablastr/particles/IndexHandling.H>

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

#ifdef ImpactX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
namespace io = openPMD;
#endif

namespace impactx::diagnostics
{
namespace detail
{
    ImpactXParticleCounter::ImpactXParticleCounter (ParticleContainer & pc)
    {
        m_MPISize = amrex::ParallelDescriptor::NProcs();
        m_MPIRank = amrex::ParallelDescriptor::MyProc();

        m_ParticleCounterByLevel.resize(pc.finestLevel()+1);
        m_ParticleOffsetAtRank.resize(pc.finestLevel()+1);
        m_ParticleSizeAtRank.resize(pc.finestLevel()+1);

        for (auto currentLevel = 0; currentLevel <= pc.finestLevel(); currentLevel++)
        {
            long numParticles = 0; // numParticles in this processor

            for (ParticleIter pti(pc, currentLevel); pti.isValid(); ++pti) {
                auto numParticleOnTile = pti.numParticles();
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
            const long& numParticles,
            unsigned long long& offset,
            unsigned long long& sum
    ) const
    {
        offset = 0;
#if defined(AMREX_USE_MPI)
        std::vector<long> result(m_MPISize, 0);
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

#ifdef ImpactX_USE_OPENPMD
    /** Unclutter a real_names to openPMD record
     *
     * TODO: move to ABLASTR
     *
     * @param fullName name as in real_names variable
     * @return pair of openPMD record and component name
     */
    inline std::pair< std::string, std::string >
    name2openPMD ( std::string const& fullName )
    {
        std::string record_name = fullName;
        std::string component_name = io::RecordComponent::SCALAR;

        // we use "_" as separator in names to group vector records
        std::size_t startComp = fullName.find_last_of("_");
        if( startComp != std::string::npos ) {  // non-scalar
            record_name = fullName.substr(0, startComp);
            component_name = fullName.substr(startComp + 1u);
        }
        return make_pair(record_name, component_name);
    }

    // TODO: move to ablastr
    io::RecordComponent get_component_record (
        io::ParticleSpecies & species,
        std::string const comp_name
    ) {
        // handle scalar and non-scalar records by name
        const auto [record_name, component_name] = name2openPMD(comp_name);
        return species[record_name][component_name];
    }
#endif
} // namespace detail

    void BeamMonitor::finalize ()
    {
        // close shared series alias
        if (m_series.has_value())
        {
            auto series = std::any_cast<io::Series>(m_series);
            series.close();
            m_series.reset();
        }

        // remove from unique series map
        if (m_unique_series.count(m_series_name) != 0u)
            m_unique_series.erase(m_series_name);
    }

    BeamMonitor::BeamMonitor (std::string series_name, std::string backend, std::string encoding) :
        m_series_name(series_name), m_OpenPMDFileType(backend)
    {
#ifdef ImpactX_USE_OPENPMD
        // pick first available backend if default is chosen
        if( m_OpenPMDFileType == "default" )
#   if openPMD_HAVE_ADIOS2==1
        m_OpenPMDFileType = "bp";
#   elif openPMD_HAVE_ADIOS1==1
        m_OpenPMDFileType = "bp";
#   elif openPMD_HAVE_HDF5==1
        m_OpenPMDFileType = "h5";
#   else
        m_OpenPMDFileType = "json";
#   endif

        // encoding of iterations in the series
        openPMD::IterationEncoding series_encoding = openPMD::IterationEncoding::groupBased;
        if ( 0 == encoding.compare("v") )
            series_encoding = openPMD::IterationEncoding::variableBased;
        else if ( 0 == encoding.compare("g") )
            series_encoding = openPMD::IterationEncoding::groupBased;
        else if ( 0 == encoding.compare("f") )
            series_encoding = openPMD::IterationEncoding::fileBased;

        // legacy options from other diagnostics
        amrex::ParmParse pp_diag("diag");
        pp_diag.queryAdd("file_min_digits", m_file_min_digits);

        // Ensure m_series is the same for the same names.
        if (m_unique_series.count(m_series_name) == 0u) {
            std::string filepath = "diags/openPMD/";
            filepath.append(m_series_name);

            if (series_encoding == openPMD::IterationEncoding::fileBased)
            {
                std::string fileSuffix = std::string("_%0") + std::to_string(m_file_min_digits) + std::string("T");
                filepath.append(fileSuffix);
            }
            filepath.append(".").append(m_OpenPMDFileType);

            // transform paths for Windows
#   ifdef _WIN32
            filepath = openPMD::auxiliary::replace_all(filepath, "/", "\\");
#   endif

            auto series = io::Series(filepath, io::Access::CREATE
#   if openPMD_HAVE_MPI==1
                , amrex::ParallelDescriptor::Communicator()
#   endif
                , "adios2.engine.usesteps = true"
            );
            series.setSoftware("ImpactX", IMPACTX_VERSION);
            series.setIterationEncoding( series_encoding );
            m_series = series;
            m_unique_series[m_series_name] = series;
        }
        else {
            m_series = m_unique_series[m_series_name];
        }
#else
        amrex::AllPrint() << "Warning: openPMD output requested but not compiled for series=" << m_series_name << "\n";
#endif
    }

    void BeamMonitor::prepare (
        PinnedContainer & pc,
        RefPart const & ref_part,
        int step
    ) {
#ifdef ImpactX_USE_OPENPMD
        m_step = step;

        // series & iteration
        auto series = std::any_cast<io::Series>(m_series);
        io::WriteIterations iterations = series.writeIterations();
        io::Iteration iteration = iterations[m_step];
        io::ParticleSpecies beam = iteration.particles["beam"];

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
        auto const getComponentRecord = [&beam](std::string const comp_name) {
            return detail::get_component_record(beam, comp_name);
        };

        // define data set and metadata
        io::Datatype dtype_fl = io::determineDatatype<amrex::ParticleReal>();
        io::Datatype dtype_ui = io::determineDatatype<uint64_t>();
        auto d_fl = io::Dataset(dtype_fl, {np});
        auto d_ui = io::Dataset(dtype_ui, {np});

        // AoS: Real
        {
            std::vector<std::string> real_aos_names(RealAoS::names_s.size());
            std::copy(RealAoS::names_s.begin(), RealAoS::names_s.end(), real_aos_names.begin());
            for (auto real_idx=0; real_idx < RealAoS::nattribs; real_idx++) {
                auto const component_name = real_aos_names.at(real_idx);
                getComponentRecord(component_name).resetDataset(d_fl);
            }
        }

        // beam mass
        beam.setAttribute( "mass", ref_part.mass );
        beam.setAttribute( "beta_ref", ref_part.beta() );
        beam.setAttribute( "gamma_ref", ref_part.gamma() );

        // openPMD coarse position
        {
            beam["positionOffset"]["x"].resetDataset(d_fl);
            beam["positionOffset"]["x"].makeConstant(ref_part.x);
            beam["positionOffset"]["y"].resetDataset(d_fl);
            beam["positionOffset"]["y"].makeConstant(ref_part.y);
            beam["positionOffset"]["t"].resetDataset(d_fl);
            beam["positionOffset"]["t"].makeConstant(ref_part.t);
        }

        // AoS: Int
        beam["id"][scalar].resetDataset(d_ui);

        // SoA: Real
        {
            std::vector<std::string> real_soa_names(RealSoA::names_s.size());
            std::copy(RealSoA::names_s.begin(), RealSoA::names_s.end(), real_soa_names.begin());
            for (auto real_idx = 0; real_idx < RealSoA::nattribs; real_idx++) {
                auto const component_name = real_soa_names.at(real_idx);
                getComponentRecord(component_name).resetDataset(d_fl);
            }
        }
        // SoA: Int
        static_assert(IntSoA::nattribs == 0); // not yet used
#else
        amrex::ignore_unused(pc, step);
#endif
    }

    void
    BeamMonitor::operator() (
        ImpactXParticleContainer & pc,
        int step
    )
    {
        // preparing to access reference particle data: RefPart
        RefPart & ref_part = pc.GetRefParticle();

        // pinned memory copy
        PinnedContainer pinned_pc = pc.make_alike<amrex::PinnedArenaAllocator>();
        pinned_pc.copyParticles(pc, true);  // no filtering

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

        // prepare element access
        this->prepare(pinned_pc, ref_part, step);

        // loop over refinement levels
        int const nLevel = pinned_pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            // loop over all particle boxes
            //using ParIt = ImpactXParticleContainer::iterator;
            using ParIt = PinnedContainer::ParIterType;
            // note: openPMD-api is not thread-safe, so do not run OMP parallel here
            for (ParIt pti(pinned_pc, lev); pti.isValid(); ++pti) {
                // push reference particle in global coordinates
                this->operator()(ref_part);

                // push beam particles relative to reference particle
                this->operator()(pti, ref_part);
            } // end loop over all particle boxes
        } // end mesh-refinement level loop

        auto series = std::any_cast<io::Series>(m_series);
        io::WriteIterations iterations = series.writeIterations();
        io::Iteration iteration = iterations[m_step];

        // close iteration
        iteration.close();
    }

    void
    BeamMonitor::operator() (
        PinnedContainer::ParIterType & pti,
        RefPart const & ref_part
    )
    {
#ifdef ImpactX_USE_OPENPMD
        int const currentLevel = pti.GetLevel();

        auto & offset = m_offset.at(currentLevel); // ...

        // preparing access to particle data: AoS
        auto& aos = pti.GetArrayOfStructs();

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
        auto const getComponentRecord = [&beam](std::string const comp_name) {
            return detail::get_component_record(beam, comp_name);
        };

        // AoS: position and particle ID
        {
            using vs = std::vector<std::string>;
            vs const positionComponents{"x", "y", "t"}; // TODO: generalize
            for (auto currDim = 0; currDim < AMREX_SPACEDIM; currDim++) {
                std::shared_ptr<amrex::ParticleReal> curr(
                    new amrex::ParticleReal[numParticleOnTile],
                    [](amrex::ParticleReal const *p) { delete[] p; }
                );
                for (auto i = 0; i < numParticleOnTile; i++) {
                    curr.get()[i] = aos[i].pos(currDim);
                }
                std::string const positionComponent = positionComponents[currDim];
                beam["position"][positionComponent].storeChunk(curr, {offset},
                                                               {numParticleOnTile64});
            }

            // save particle ID after converting it to a globally unique ID
            std::shared_ptr<uint64_t> ids(
                new uint64_t[numParticleOnTile],
                [](uint64_t const *p) { delete[] p; }
            );
            for (auto i = 0; i < numParticleOnTile; i++) {
                ids.get()[i] = ablastr::particles::localIDtoGlobal(aos[i].id(), aos[i].cpu());
            }
            beam["id"][scalar].storeChunk(ids, {offset}, {numParticleOnTile64});
        }

        // SoA: everything else
        auto const& soa = pti.GetStructOfArrays();
        //   SoA floating point (ParticleReal) properties
        {
            std::vector<std::string> real_soa_names(RealSoA::names_s.size());
            std::copy(RealSoA::names_s.begin(), RealSoA::names_s.end(), real_soa_names.begin());

            for (auto real_idx=0; real_idx < RealSoA::nattribs; real_idx++) {
                auto const component_name = real_soa_names.at(real_idx);
                getComponentRecord(component_name).storeChunkRaw(
                soa.GetRealData(real_idx).data(), {offset}, {numParticleOnTile64});
            }
        }
        //   SoA integer (int) properties (not yet used)
        {
            static_assert(IntSoA::nattribs == 0); // not yet used
            /*
            // comment this in once IntSoA::nattribs is > 0

            std::vector<std::string> int_soa_names(IntSoA::names_s.size);
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

        // TODO could be done once after all pti are processed
        // TODO at that point, we could also close the iteration/step
        series.flush();
#else
        amrex::ignore_unused(pti, ref_part);
#endif
    }

} // namespace impactx::diagnostics
