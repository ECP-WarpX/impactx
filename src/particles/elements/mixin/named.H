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
#include <AMReX_GpuQualifiers.H>

#include <cstring>
#include <optional>
#include <stdexcept>
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
        AMREX_GPU_HOST
        Named (
            std::optional<std::string> name
        )
        {
            if (name.has_value()) {
                m_name = new char[name->size() + 1];
                std::strcpy(m_name, name->c_str());
            }
        }

        AMREX_GPU_HOST_DEVICE
        ~Named ()
        {
            AMREX_IF_ON_HOST((
                if (m_name != nullptr) {
                    delete[] m_name;
                    m_name = nullptr;
                }
            ))
        }

        AMREX_GPU_HOST
        Named (Named const & other)
        {
            if (other.has_name()) {
                m_name = new char[std::strlen(other.m_name) + 1];
                std::strcpy(m_name, other.m_name);
            }
        }

        AMREX_GPU_HOST
        Named& operator=(Named const & other)
        {
            if (&other != this) {
                if (other.has_name()) {
                    m_name = new char[std::strlen(other.m_name) + 1];
                    std::strcpy(m_name, other.m_name);
                }
            }
            return *this;
        }

        AMREX_GPU_HOST_DEVICE
        Named ([[maybe_unused]] Named && other) noexcept
        {
            AMREX_IF_ON_HOST((
                std::swap(this->m_name, other.m_name);
                other.m_name = nullptr;
            ))
        }

        AMREX_GPU_HOST_DEVICE
        Named& operator=([[maybe_unused]] Named && other) noexcept
        {
            AMREX_IF_ON_HOST((
                std::swap(this->m_name, other.m_name);
                other.m_name = nullptr;
            ))
            return *this;
        }

        /** Return the user provided name of the element
         *
         * @return user defined and not necessarily unique name
         */
        AMREX_FORCE_INLINE
        std::string name () const
        {
            if (!has_name()) {
                throw std::runtime_error("Name not set on element!");
            }
            return std::string(m_name);
        }

        /** Return true if the user provided a name for this element
         *
         * @return true if the user provided a name
         */
        AMREX_FORCE_INLINE
        bool has_name () const
        {
            return m_name != nullptr;
        }

        // Implementation note:
        // This is used as a mixin class in elements that are copied to GPU. GPU compilers copy
        // a whole element by its sizeof(T).
        // To save on this copy operation at kernel start, this is implemented
        // as a simple C pointer (8 bytes), contrary to a std::string (32 bytes) or
        // a std::optional<std::string> (40 bytes). m_name points to host-side memory and
        // must not be dereferenced (used) on GPU.
        char * m_name = nullptr; //! a user defined and not necessarily unique name of the element
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_NAMED_H
