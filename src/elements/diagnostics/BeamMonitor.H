/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_DIAGS_BEAMMONITOR_H
#define IMPACTX_ELEMENTS_DIAGS_BEAMMONITOR_H

#include "elements/mixin/thin.H"
#include "particles/CovarianceMatrix.H"
#include "particles/ImpactXParticleContainer.H"
#include "elements/mixin/lineartransport.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <any>
#include <string>
#include <unordered_map>
#include <vector>


namespace impactx::elements::diagnostics
{
namespace detail
{
    class ImpactXParticleCounter {
    public:
        using ParticleContainer = typename ImpactXParticleContainer::ContainerLike<amrex::PolymorphicArenaAllocator>;
        using ParticleIter = typename ParticleContainer::ParIterType;

        ImpactXParticleCounter (ParticleContainer & pc);

        unsigned long long GetTotalNumParticles () { return m_Total; }

        std::vector<unsigned long long> m_ParticleOffsetAtRank;
        std::vector<unsigned long long> m_ParticleSizeAtRank;
    private:
        /** get the offset in the overall particle id collection
        *
        * @param[out] numParticles particles on this processor  / amrex fab
        * @param[out] offset particle offset over all, mpi-global amrex fabs
        * @param[out] sum number of all particles from all amrex fabs
        */
        void GetParticleOffsetOfProcessor (amrex::Long const &numParticles,
                                           unsigned long long &offset,
                                           unsigned long long &sum) const;

        int m_MPIRank = 0;
        int m_MPISize = 1;

        unsigned long long m_Total = 0;

        std::vector<unsigned long long> m_ParticleCounterByLevel;
    };
} // namespace detail

    /** This element writes the particle beam out to openPMD data.
     *
     * This class behaves like a singleton if constructed with the
     * same series name as an existing instance.
     */
    struct BeamMonitor
    : public mixin::LinearTransport<BeamMonitor>,
      public mixin::Thin
    {
        static constexpr auto type = "BeamMonitor";
        using PType = typename ImpactXParticleContainer::ParticleType;
        using PinnedContainer = typename ImpactXParticleContainer::ContainerLike<amrex::PolymorphicArenaAllocator>;

        static constexpr char const * DEFAULT_backend = "default";
        static constexpr char const * DEFAULT_encoding = "g";
        static constexpr int DEFAULT_period_sample_intervals = 1;
        static constexpr bool DEFAULT_particles = true;
        static constexpr bool DEFAULT_beam_moments = true;

        /** This element writes the particle beam out to openPMD data.
         *
         * Elements with the same series name are identical.
         *
         * @param series_name name of the data series, usually the element name
         * @param backend file format backend for openPMD, e.g., "bp5", "bp4", or "h5"
         * @param encoding openPMD iteration encoding: "v"ariable based, "f"ile based, "g"roup based (default)
         * @param period_sample_intervals for periodic lattice, only output every Nth period (turn)
         * @param particles output the individual particles of the beam (default: true)
         * @param beam_moments output the beam moments (reduced beam characteristics) as attributes (default: true)
         */
        BeamMonitor (
            std::string series_name,
            std::string backend = DEFAULT_backend,
            std::string encoding = DEFAULT_encoding,
            int period_sample_intervals = DEFAULT_period_sample_intervals,
            bool particles = DEFAULT_particles,
            bool beam_moments = DEFAULT_beam_moments
        );

        BeamMonitor (BeamMonitor const & other) = default;
        BeamMonitor (BeamMonitor && other) = default;
        BeamMonitor& operator= (BeamMonitor const & other) = default;
        BeamMonitor& operator= (BeamMonitor && other) = default;

        /** A distinct monitor with the same output configuration
         *
         * The open output series is deliberately not carried over: two monitors writing
         * through one series handle would interleave into the same iterations. The copy
         * opens its own on first use.
         */
        [[nodiscard]] BeamMonitor copy () const
        {
            BeamMonitor out(*this);
            out.m_series.reset();
            out.m_series_key.clear();
            out.m_step = 0;
            out.m_offset.clear();
            out.m_rbc.clear();
            return out;
        }

        /** Reverse the element (no-op for diagnostics) */
        void reverse () { }

        /** Prepare entering the element before starting push logic.
         *
         * And write reference particle.
         *
         * @param[in] pc particle container
         * @param[in] real_soa_names ParticleReal component names
         * @param[in] int_soa_names integer component names
         * @param[in] ref_part reference particle
         * @param[in] step global step for diagnostics
         */
        void prepare (
            PinnedContainer & pc,
            std::vector<std::string> const & real_soa_names,
            std::vector<std::string> const & int_soa_names,
            RefPart const & ref_part,
            int step
        );

        /** Dump all particles.
         *
         * Particles are relative to the reference particle.
         *
         * @param[in,out] pc particle container to push
         * @param[in] step global step for diagnostics
         * @param[in] period for periodic lattices, this is the current period (turn or cycle)
         */
        void operator() (
            ImpactXParticleContainer & pc,
            int step,
            int period
        );

        /** Write a tile of particles
         *
         * @param pti particle tile iterator
         * @param[in] real_soa_names ParticleReal component names
         * @param[in] int_soa_names integer component names
         * @param ref_part reference particle
         */
        void operator() (
            PinnedContainer::ParIterType & pti,
            std::vector<std::string> const & real_soa_names,
            std::vector<std::string> const & int_soa_names,
            RefPart const & ref_part
        );

        /** This does nothing to the reference particle. */
        using Thin::operator();

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * @param[in] refpart reference particle
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            // nothing to do
            return Map6x6::Identity();
        }

        /** Get the name of the series (and element)
         *
         * Elements with the same series name are identical.
         */
        AMREX_FORCE_INLINE
        std::string name () const { return m_series_name; }

        /** Provide impactx::elements::mixin::Named::has_name */
        AMREX_FORCE_INLINE
        bool has_name () const { return true; }

        /** openPMD backend string passed at construction (e.g. "default", "bp4", "h5") */
        AMREX_FORCE_INLINE
        std::string backend () const { return m_OpenPMDFileType; }

        /** openPMD iteration encoding: "v", "f", or "g" */
        AMREX_FORCE_INLINE
        std::string encoding () const { return m_encoding; }

        /** For periodic lattices, only output every Nth period (turn or cycle) */
        AMREX_FORCE_INLINE
        int period_sample_intervals () const { return m_period_sample_intervals; }

        /** Output the individual particles of the beam */
        AMREX_FORCE_INLINE
        bool particles () const { return m_particles; }

        /** Toggle output of the individual particles of the beam */
        void set_particles (bool particles) { m_particles = particles; }

        /** Output the beam moments (reduced beam characteristics) as openPMD attributes */
        AMREX_FORCE_INLINE
        bool beam_moments () const { return m_beam_moments; }

        /** Toggle output of the beam moments */
        void set_beam_moments (bool beam_moments) { m_beam_moments = beam_moments; }

        /** track all openPMD series instances
         *
         * Ensure m_series is the same for the same output path.
         */
        static inline std::map<std::string, std::any> m_unique_series = {};

        /** Close and deallocate all data series and backends.
         */
        void
        finalize ();

    private:
        /** Open the openPMD series (on first access) */
        void open ();

        /** Whether to write through the steps-based openPMD writeIterations() API
         *
         * I/O steps are the recommended way to write and are fully supported
         * by ADIOS2 BP5, HDF5 and JSON, with one exception: the ADIOS2 BP4
         * engine drops attribute-only iterations written through I/O steps
         * (verified with openPMD-api 0.16.1 and 0.17.1). Thus, when no
         * particle data sets are written, use random access for "bp4". The
         * generic "bp" ending resolves to ADIOS2's default engine, which is
         * BP5 for ADIOS2 >= 2.9. Variable-based encoding strictly requires
         * I/O steps (random access overwrites attributes of prior
         * iterations).
         * Note: openPMD-api cannot report the selected ADIOS2
         *   engine (Series::backend() only returns "ADIOS2"), so this checks
         *   the resolved file ending.
         *   https://github.com/openPMD/openPMD-api/issues/1913
         */
        bool use_io_steps () const
        {
            return m_particles ||
                   m_encoding == "v" ||
                   m_OpenPMDFileType != "bp4";
        }

        std::string m_series_name; //! openPMD filename
        std::string m_OpenPMDFileType; //! openPMD backend: usually HDF5 (h5) or ADIOS2 (bp/bp4/bp5) or ADIOS2 SST (sst)
        std::string m_encoding; //! openPMD iteration encoding: "v"ariable based, "f"ile based, "g"roup based (default)
        std::any m_series; //! openPMD::Series that holds potentially multiple outputs
        std::string m_series_key; //! full path used as key in m_unique_series
        int m_step = 0; //! global step for output

        int m_file_min_digits = 6; //! minimum number of digits to iteration number in file name

        int m_period_sample_intervals = DEFAULT_period_sample_intervals; //! only output every Nth period (turn or cycle) of periodic lattices

        bool m_particles = DEFAULT_particles; //! output the individual particles of the beam
        bool m_beam_moments = DEFAULT_beam_moments; //! output the beam moments as openPMD attributes

        /** This rank's offset in the MPI-global particle array, by level
         *
         * This MUST be updated by prepare() before the next step's output.
         */
        std::vector<uint64_t> m_offset;

        /** Reduced Beam Characteristics
         *
         * In situ calculated particle bunch moments.
         */
        std::unordered_map<std::string, amrex::ParticleReal> m_rbc;

    };

    /** Calculate additional particle properties.
     *
     * @param[in] element_name name of the element (for input queries)
     * @param[inout] pc particle container
     */
    void
    add_optional_properties (
        std::string const & element_name,
        ImpactXParticleContainer & pc
    );

} // namespace impactx::diagnostics

#endif // IMPACTX_ELEMENTS_DIAGS_BEAMMONITOR_H
