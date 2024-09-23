/* Copyright 2022-2024 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_NAMED_H
#define IMPACTX_ELEMENTS_MIXIN_NAMED_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Extension.H>

#include <optional>
#include <string>


namespace impactx::elements
{
    /** This is a helper class for lattice elements with a user-defined name
     */
    struct Named
    {
        /** A user-named element
         *
         * @param name a user defined and not necessarily unique name of the element
         */
        Named (
            std::optional<std::string> name
        )
        : m_name(name)
        {
        }

        /** Return the user provided name of the element
         *
         * @return user defined and not necessarily unique name
         */
        AMREX_FORCE_INLINE
        std::string name () const
        {
            return m_name.value();
        }

        /** Return true if the user provided a name for this element
         *
         * @return true if the user provided a name
         */
        AMREX_FORCE_INLINE
        bool has_name () const
        {
            return m_name.has_value();
        }

        std::optional<std::string> m_name; //! a user defined and not necessarily unique name of the element
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_NAMED_H
