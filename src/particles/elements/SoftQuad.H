/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SOFTQUAD_H
#define IMPACTX_SOFTQUAD_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/integrators/Integrators.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"

#include <ablastr/constant.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <array>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <vector>


namespace impactx
{
   /** Default Fourier coefficients
    *
    * Fourier coefficients for the on-axis quadrupole field gradient
    * k = 1/(Brho)*dBy/dx as a function of longitudinal position z.
    * The default on-axis gradient is described by the function:
    * k(z) = gscale/2*(tanh((z-zi)/g) - tanh((z-zf)/g)),
    * where zf = -zi = L/2, L is the magnetic length of the
    * quadrupole, and g is the full gap parameter.  Here the
    * parameters are chosen so that g/L = 1/10.  This is
    * based on a similar model used for soft-edge elements
    * in MaryLie and ML/I, eg.
    * http://www.physics.umd.edu/dsat/docs/MaryLieMan.pdf
    *
    * This is equivalent to a single-parameter Enge model for
    * the fields in the fringe region.
    *
    * Usually superseded by user-provided input values.
    */
    struct Quad_field_data
    {
       amrex::Vector<amrex::ParticleReal> default_cos_coef = {
             0.834166514794446,
             0.598104328994702,
             0.141852844428785,
            -0.118211272182381,
            -9.056149864743113E-002,
             1.803476331179615E-002,
             4.464887700797893E-002,
             7.364410636252136E-003,
            -1.697541023436736E-002,
            -9.012679515542771E-003,
             4.367667630047725E-003,
             5.444030542119803E-003,
            -5.889959910931886E-005,
            -2.409098101409192E-003,
            -7.962712154165590E-004,
             7.855814707106538E-004,
             6.174930463182168E-004,
            -1.340154094301854E-004,
            -3.167213724698439E-004,
            -4.925292460592617E-005,
             1.221580597451921E-004,
             6.331025910961789E-005,
            -3.202416719002774E-005,
            -3.872103053895529E-005,
             8.212882937116278E-007
            };

       amrex::Vector<amrex::ParticleReal> default_sin_coef = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0
            };
    };

/** Dynamic data for the SoftQuadrupole elements
 *
 * Since we copy the element to the device, we cannot store this data on the element itself.
 * But we can store pointers to this data with the element and keep a lookup table here,
 * which we clean up in the end.
 */
namespace SoftQuadrupoleData
{
    //! last used id for a created soft quad
    static inline int next_id = 0;

    //! host: cosine coefficients in Fourier expansion of on-axis magnetic field Bz
    static inline std::map<int, std::vector<amrex::ParticleReal>> h_cos_coef = {};
    //! host: sine coefficients in Fourier expansion of on-axis magnetic field Bz
    static inline std::map<int, std::vector<amrex::ParticleReal>> h_sin_coef = {};

    //! device: cosine coefficients in Fourier expansion of on-axis magnetic field Bz
    static inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_cos_coef = {};
    //! device: sine coefficients in Fourier expansion of on-axis magnetic field Bz
    static inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_sin_coef = {};

} // namespace SoftQuadrupoleData

    struct SoftQuadrupole
    : public elements::BeamOptic<SoftQuadrupole>,
      public elements::Thick,
      public elements::Alignment
    {
        static constexpr auto name = "SoftQuadrupole";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A soft-edge quadrupole
         *
         * @param ds Segment length in m
         * @param gscale Scaling factor for on-axis field gradient Bz in 1/m^2
         * @param cos_coef cosine coefficients in Fourier expansion of on-axis magnetic field Bz
         * @param sin_coef sine coefficients in Fourier expansion of on-axis magnetic field Bz
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error around the s axis (in the x-y plane) [degrees]
         * @param mapsteps number of integration steps per slice used for
         *        map and reference particle push in applied fields
         * @param nslice number of slices used for the application of space charge
         */
        SoftQuadrupole (
            amrex::ParticleReal ds,
            amrex::ParticleReal gscale,
            std::vector<amrex::ParticleReal> cos_coef,
            std::vector<amrex::ParticleReal> sin_coef,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            int mapsteps = 1,
            int nslice = 1
        )
          : Thick(ds, nslice),
            Alignment(dx, dy, rotation_degree),
            m_gscale(gscale), m_mapsteps(mapsteps), m_id(SoftQuadrupoleData::next_id)
        {
            // next created soft quad has another id for its data
            SoftQuadrupoleData::next_id++;

            // validate sin and cos coefficients are the same length
            m_ncoef = int(cos_coef.size());
            if (m_ncoef != int(sin_coef.size()))
                throw std::runtime_error("SoftQuadrupole: cos and sin coefficients must have same length!");

            // host data
            SoftQuadrupoleData::h_cos_coef[m_id] = cos_coef;
            SoftQuadrupoleData::h_sin_coef[m_id] = sin_coef;
            m_cos_h_data = SoftQuadrupoleData::h_cos_coef[m_id].data();
            m_sin_h_data = SoftQuadrupoleData::h_sin_coef[m_id].data();

            // device data
            SoftQuadrupoleData::d_cos_coef.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_ncoef));
            SoftQuadrupoleData::d_sin_coef.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_ncoef));
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  cos_coef.begin(), cos_coef.end(),
                                  SoftQuadrupoleData::d_cos_coef[m_id].begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  sin_coef.begin(), sin_coef.end(),
                                  SoftQuadrupoleData::d_sin_coef[m_id].begin());
            amrex::Gpu::streamSynchronize();

            // low-level objects we can use on device
            m_cos_d_data = SoftQuadrupoleData::d_cos_coef[m_id].data();
            m_sin_d_data = SoftQuadrupoleData::d_sin_coef[m_id].data();
       }

        /** Push all particles */
        using BeamOptic::operator();

        /** This is a soft-edge quadrupole functor, so that a variable of this type can be used
         *  like a soft-edge quadrupole function.
         *
         * @param p Particle AoS data for positions and cpu/id
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            PType& AMREX_RESTRICT p,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            [[maybe_unused]] RefPart const & refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);

            // access AoS data such as positions and cpu/id
            amrex::ParticleReal const x = p.pos(RealAoS::x);
            amrex::ParticleReal const y = p.pos(RealAoS::y);
            amrex::ParticleReal const t = p.pos(RealAoS::t);

            // initialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // get the linear map
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;

            // symplectic linear map for a quadrupole is computed using the
            // Hamiltonian formalism as described in:
            // https://uspas.fnal.gov/materials/09UNM/ComputationalMethods.pdf .
            // R denotes the transfer matrix in the basis (x,px,y,py,t,pt),
            // so that, e.g., R(3,4) = dyf/dpyi.

            // push particles using the linear map
            p.pos(RealAoS::x) = R(1,1)*x + R(1,2)*px + R(1,3)*y
                     + R(1,4)*py + R(1,5)*t + R(1,6)*pt;
            pxout = R(2,1)*x + R(2,2)*px + R(2,3)*y
                  + R(2,4)*py + R(2,5)*t + R(2,6)*pt;
            p.pos(RealAoS::y) = R(3,1)*x + R(3,2)*px + R(3,3)*y
                     + R(3,4)*py + R(3,5)*t + R(3,6)*pt;
            pyout = R(4,1)*x + R(4,2)*px + R(4,3)*y
                  + R(4,4)*py + R(4,5)*t + R(4,6)*pt;
            p.pos(RealAoS::t) = R(5,1)*x + R(5,2)*px + R(5,3)*y
                     + R(5,4)*py + R(5,5)*t + R(5,6)*pt;
            ptout = R(6,1)*x + R(6,2)*px + R(6,3)*y
                  + R(6,4)*py + R(6,5)*t + R(6,6)*pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(p.pos(RealAoS::x), p.pos(RealAoS::y), px, py);
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // assign input reference particle values
            amrex::ParticleReal const x = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const s = refpart.s;
            amrex::ParticleReal const sedge = refpart.sedge;

            // initialize linear map (deviation) values
            for (int i=1; i<7; i++) {
               for (int j=1; j<7; j++) {
                  auto const default_value = (i == j) ? 1.0_prt : 0.0_prt;
                  refpart.map(i, j) = default_value;
               }
            }

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // compute intial value of beta*gamma
            amrex::ParticleReal const bgi = sqrt(pow(pt, 2) - 1.0_prt);

            // call integrator to advance (t,pt)
            amrex::ParticleReal const zin = s - sedge;
            amrex::ParticleReal const zout = zin + slice_ds;
            int const nsteps = m_mapsteps;

            integrators::symp2_integrate(refpart,zin,zout,nsteps,*this);
            amrex::ParticleReal const ptf = refpart.pt;

            /*
            // print computed linear map:
               for(int i=1; i<7; ++i){
                 for(int j=1; j<7; ++j){
                    amrex::PrintToFile("QuadMap.txt") << i << " " <<
                    j << " " << refpart.map(i,j) << "\n";
                 }
               }
            //
            */

            // advance position (x,y,z)
            refpart.x = x + slice_ds*px/bgi;
            refpart.y = y + slice_ds*py/bgi;
            refpart.z = z + slice_ds*pz/bgi;

            // compute final value of beta*gamma
            amrex::ParticleReal const bgf = sqrt(pow(ptf, 2) - 1.0_prt);

            // advance momentum (px,py,pz)
            refpart.px = px*bgf/bgi;
            refpart.py = py*bgf/bgi;
            refpart.pz = pz*bgf/bgi;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** This evaluates the on-axis magnetic field Bz at a fixed location
         *  z, together with certain required integrals and derivatives.
         *  The field returned is normalized to a peak value of 1.
         *
         * @param zeval Longitudinal on-axis location in m
         */
        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Quad_Bfield (amrex::ParticleReal const zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // pick the right data depending if we are on the host side
            // (reference particle push) or device side (particles):
#if AMREX_DEVICE_COMPILE
            amrex::ParticleReal* cos_data = m_cos_d_data;
            amrex::ParticleReal* sin_data = m_sin_d_data;
#else
            amrex::ParticleReal* cos_data = m_cos_h_data;
            amrex::ParticleReal* sin_data = m_sin_h_data;
#endif

            // specify constants
            using ablastr::constant::math::pi;
            amrex::ParticleReal const zlen = m_ds;
            amrex::ParticleReal const zmid = zlen / 2.0_prt;

            // compute on-axis magnetic field (z is relative to quadrupole midpoint)
            amrex::ParticleReal bfield = 0.0;
            amrex::ParticleReal bfieldp = 0.0;
            amrex::ParticleReal bfieldint = 0.0;
            amrex::ParticleReal const z = zeval - zmid;

            if (std::abs(z) <= zmid)
            {
               bfield = 0.5_prt*cos_data[0];
               bfieldint = z*bfield;
               for (int j=1; j < m_ncoef; ++j)
               {
                 bfield = bfield + cos_data[j] * cos(j * 2 * pi * z / zlen) +
                         sin_data[j] * sin(j * 2 * pi * z / zlen);
                 bfieldp = bfieldp - j * 2 * pi * cos_data[j] * sin(j * 2 * pi * z / zlen) / zlen +
                           j * 2 * pi * sin_data[j] * cos(j * 2 * pi * z / zlen) / zlen;
                 bfieldint = bfieldint + zlen * cos_data[j] * sin(j * 2 * pi * z / zlen) / (j * 2 * pi) -
                             zlen * sin_data[j] * cos(j * 2 * pi * z / zlen) / (j * 2 * pi);
               }
            }
            return std::make_tuple(bfield, bfieldp, bfieldint);
        }

        /** This pushes the reference particle and the linear map matrix
         *  elements for a quadrupole through the symplectic map associated
         *  with H_1 in the Hamiltonian splitting H = H_1 + H_2.
         *
         * @param tau Map step size in m
         * @param[in,out] refpart reference particle
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map1 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   [[maybe_unused]] amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // push the reference particle
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const z = zeval;

            if (pt < -1.0_prt) {
                refpart.t = t + tau/sqrt(1.0_prt - pow(pt, -2));
                refpart.pt = pt;
            }
            else {
                refpart.t = t;
                refpart.pt = pt;
            }

            zeval = z + tau;

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
            amrex::ParticleReal const betgam = refpart.beta_gamma();

            refpart.map(1,1) = R(1,1) + tau*R(2,1);
            refpart.map(1,2) = R(1,2) + tau*R(2,2);
            refpart.map(1,3) = R(1,3) + tau*R(2,3);
            refpart.map(1,4) = R(1,4) + tau*R(2,4);

            refpart.map(3,1) = R(3,1) + tau*R(4,1);
            refpart.map(3,2) = R(3,2) + tau*R(4,2);
            refpart.map(3,3) = R(3,3) + tau*R(4,3);
            refpart.map(3,4) = R(3,4) + tau*R(4,4);

            refpart.map(5,5) = R(5,5) + tau*R(6,5)/pow(betgam,2);
            refpart.map(5,6) = R(5,6) + tau*R(6,6)/pow(betgam,2);

        }

        /** This pushes the reference particle and the linear map matrix
         *  elements for a quadrupole through the symplectic map associated
         *  with H_2 in the Hamiltonian splitting H = H_1 + H_2.
         *
         * @param tau Map step size in m
         * @param[in,out] refpart reference particle
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map2 (amrex::ParticleReal const tau,
                   RefPart & refpart,
                   amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;

            // Define parameters and intermediate constants
            amrex::ParticleReal const G0 = m_gscale;

            // push the reference particle
            auto [bz, bzp, bzint] = Quad_Bfield(zeval);
            amrex::ignore_unused(bzp, bzint);

            refpart.t = t;
            refpart.pt = pt;

            // push the linear map equations
            amrex::Array2D<amrex::ParticleReal, 1, 6, 1, 6> const R = refpart.map;
            amrex::ParticleReal const alpha = G0*bz;

            refpart.map(2,1) = R(2,1) - tau*alpha*R(1,1);
            refpart.map(2,2) = R(2,2) - tau*alpha*R(1,2);
            refpart.map(2,3) = R(2,3) - tau*alpha*R(1,3);
            refpart.map(2,4) = R(2,4) - tau*alpha*R(1,4);

            refpart.map(4,1) = R(4,1) + tau*alpha*R(3,1);
            refpart.map(4,2) = R(4,2) + tau*alpha*R(3,2);
            refpart.map(4,3) = R(4,3) + tau*alpha*R(3,3);
            refpart.map(4,4) = R(4,4) + tau*alpha*R(3,4);

        }

        /** Close and deallocate all data and handles.
         */
        void
        finalize ()
        {
            // remove from unique data map
            if (SoftQuadrupoleData::h_cos_coef.count(m_id) != 0u)
                SoftQuadrupoleData::h_cos_coef.erase(m_id);
            if (SoftQuadrupoleData::h_sin_coef.count(m_id) != 0u)
                SoftQuadrupoleData::h_sin_coef.erase(m_id);

            if (SoftQuadrupoleData::d_cos_coef.count(m_id) != 0u)
                SoftQuadrupoleData::d_cos_coef.erase(m_id);
            if (SoftQuadrupoleData::d_sin_coef.count(m_id) != 0u)
                SoftQuadrupoleData::d_sin_coef.erase(m_id);
        }

    private:
        amrex::ParticleReal m_gscale; //! scaling factor for quad field gradient
        int m_mapsteps; //! number of map integration steps per slice
        int m_id; //! unique soft quad id used for data lookup map

        int m_ncoef = 0; //! number of Fourier coefficients
        amrex::ParticleReal* m_cos_h_data = nullptr; //! non-owning pointer to host cosine coefficients
        amrex::ParticleReal* m_sin_h_data = nullptr; //! non-owning pointer to host sine coefficients
        amrex::ParticleReal* m_cos_d_data = nullptr; //! non-owning pointer to device cosine coefficients
        amrex::ParticleReal* m_sin_d_data = nullptr; //! non-owning pointer to device sine coefficients
    };

} // namespace impactx

#endif // IMPACTX_SOFTQUAD_H
