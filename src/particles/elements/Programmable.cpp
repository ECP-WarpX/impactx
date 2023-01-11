/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */

#include "Programmable.H"


namespace impactx
{
    void
    Programmable::operator() (
        ImpactXParticleContainer::iterator & pti,
        RefPart & ref_part
    ) const
    {
        if (m_beam_particles == nullptr)
            amrex::AllPrint() << "Programmable element - all particles: NO HOOK\n";
        else
            m_beam_particles(&pti, ref_part);
    }

    void
    Programmable::operator() (RefPart & ref_part) const
    {
        if (m_ref_particle == nullptr)
            amrex::AllPrint() << "Programmable element - ref particles: NO HOOK\n";
        else
            m_ref_particle(ref_part);
    }

} // namespace impactx
