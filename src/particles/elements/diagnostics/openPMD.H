/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_DIAGS_OPENPMD_H
#define IMPACTX_ELEMENTS_DIAGS_OPENPMD_H

#include "particles/elements/mixin/thin.H"
#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <any>
#include <string>
#include <vector>


namespace impactx::diagnostics
{
namespace detail
{
    class ImpactXParticleCounter {
    public:
        using ParticleContainer = typename ImpactXParticleContainer::ContainerLike<amrex::PinnedArenaAllocator>;
        using ParticleIter = typename ParticleContainer::ParIterType;

        ImpactXParticleCounter (ParticleContainer & pc);

        unsigned long GetTotalNumParticles () { return m_Total; }

        std::vector<unsigned long long> m_ParticleOffsetAtRank;
        std::vector<unsigned long long> m_ParticleSizeAtRank;
    private:
        /** get the offset in the overall particle id collection
        *
        * @param[out] numParticles particles on this processor  / amrex fab
        * @param[out] offset particle offset over all, mpi-global amrex fabs
        * @param[out] sum number of all particles from all amrex fabs
        */
        void GetParticleOffsetOfProcessor (const long &numParticles,
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
    : public elements::Thin
    {
        static constexpr auto name = "BeamMonitor";
        using PType = typename ImpactXParticleContainer::ParticleType;
        using PinnedContainer = typename ImpactXParticleContainer::ContainerLike<amrex::PinnedArenaAllocator>;

        /** This element writes the particle beam out to openPMD data.
         *
         * Elements with the same series name are identical.
         *
         * @param series_name name of the data series, usually the element name
         * @param backend file format backend for openPMD, e.g., "bp" or "h5"
         * @param encoding openPMD iteration encoding: "v"ariable based, "f"ile based, "g"roup based (default)
         */
        BeamMonitor (std::string series_name, std::string backend="default", std::string encoding="g");

        BeamMonitor (BeamMonitor const & other) = default;
        BeamMonitor (BeamMonitor && other) = default;
        BeamMonitor& operator= (BeamMonitor const & other) = default;
        BeamMonitor& operator= (BeamMonitor && other) = default;

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
         */
        void operator() (
            ImpactXParticleContainer & pc,
            int step
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

        /** Get the name of the series
         *
         * Elements with the same series name are identical.
         */
        std::string series_name () const { return m_series_name; }

        /** track all m_series_name instances
         *
         * Ensure m_series is the same for the same name.
         */
        static inline std::map<std::string, std::any> m_unique_series = {};

        /** Close and deallocate all data series and backends.
         */
        void
        finalize ();

    private:
        std::string m_series_name; //! ...
        std::string m_OpenPMDFileType; //! ...
        std::any m_series; //! openPMD::Series; ...
        int m_step = 0; //! ...

        int m_file_min_digits = 6; //! minimum number of digits to iteration number in file name

        /** This rank's offset in the MPI-global particle array, by level
         *
         * This MUST be updated by prepare() before the next step's output.
         */
        std::vector<uint64_t> m_offset;

    };

} // namespace impactx::diagnostics

#endif // IMPACTX_ELEMENTS_DIAGS_OPENPMD_H
