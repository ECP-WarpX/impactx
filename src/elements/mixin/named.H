/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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

#include "hostonly.H"

#include <AMReX_Extension.H>

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>


namespace impactx::elements::mixin
{
    /** This is a helper class for lattice elements with a user-defined name
     *
     * A name is host-only metadata: it describes the element to a human, and takes no
     * part in the physics. It is carried by the public element envelope
     * (@see impactx::elements::Element), which keeps it out of the execution snapshot
     * that is copied into a particle-push kernel.
     *
     * That separation is what lets this be an ordinary rule-of-zero value type.
     */
    struct Named : HostOnly
    {
        static constexpr std::nullopt_t DEFAULT_name = std::nullopt;

        /** An unnamed element */
        Named () = default;

        /** A user-named element
         *
         * An empty name is treated as no name.
         *
         * @param name a user defined and not necessarily unique name of the element
         */
        Named (
            std::optional<std::string> name
        )
          : m_name(std::move(name))
        {
            if (m_name.has_value() && m_name->empty())
            {
                m_name.reset();
            }
        }

        /** Overwrite the name of the element
         *
         * An empty name unsets the name.
         *
         * @param new_name set a new name on the element
         */
        void set_name (
            std::string const & new_name
        )
        {
            if (new_name.empty())
            {
                m_name.reset();
            }
            else
            {
                m_name = new_name;
            }
        }

        /** Return the user provided name of the element
         *
         * @return user defined and not necessarily unique name
         */
        [[nodiscard]] AMREX_FORCE_INLINE
        std::string name () const
        {
            if (!has_name())
            {
                throw std::runtime_error("Name not set on element!");
            }
            return *m_name;
        }

        /** Return true if the user provided a name for this element
         *
         * @return true if the user provided a name
         */
        [[nodiscard]] AMREX_FORCE_INLINE
        bool has_name () const noexcept
        {
            return m_name.has_value();
        }

    private:
        std::optional<std::string> m_name; //! a user defined and not necessarily unique name of the element
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_NAMED_H
