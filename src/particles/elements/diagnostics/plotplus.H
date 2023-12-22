/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Junmin Gu
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_DIAGS_OPENPMD_WITH_AMREX_H
#define IMPACTX_ELEMENTS_DIAGS_OPENPMD_WITH_AMREX_H

#include "particles/elements/mixin/thin.H"
#include "particles/ImpactXParticleContainer.H"

#include <AMReX_PlotFileUtil.H> // for AMReXWithOpenPMD
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <any>

#ifdef AMREX_USE_OPENPMD_API
namespace impactx::diagnostics
{
  class AMReX_impactxWriter: public amrex::openpmd_api::AMReX_openPMDWriter
  {
  public:
    AMReX_impactxWriter(openPMD::IterationEncoding ie =  openPMD::IterationEncoding::variableBased,
            std::string filetype = "default",
            std::string openpmdOptions = "{}")
      :AMReX_openPMDWriter("{impactX}", ie, filetype, openpmdOptions)
    {}


    // subclasses
    inline std::vector< std::string > getParticlePositionComponentLabels() const
    {
      using vs = std::vector< std::string >;
      vs const positionComponents{"x", "y", "t"};
      return positionComponents;
    }
    // internal  to impact x
    /** Prepare entering the element before starting push logic.
     *
     * @param[in] pc particle container
     * @param[in] ref_part reference particle
     * @param[in] step global step for diagnostics
     */
    void SetConstantRefPart (openPMD::ParticleSpecies& beam,
                 //const unsigned long long& np,
                 RefPart const & ref_part
                 )
    {
      // beam mass
        beam.setAttribute( "beta_ref", ref_part.beta() );
        beam.setAttribute( "gamma_ref", ref_part.gamma() );
        beam.setAttribute( "s_ref", ref_part.s );
        beam.setAttribute( "x_ref", ref_part.x );
        beam.setAttribute( "y_ref", ref_part.y );
        beam.setAttribute( "z_ref", ref_part.z );
        beam.setAttribute( "t_ref", ref_part.t );
        beam.setAttribute( "px_ref", ref_part.px );
        beam.setAttribute( "py_ref", ref_part.py );
        beam.setAttribute( "pz_ref", ref_part.pz );
        beam.setAttribute( "pt_ref", ref_part.pt );
        beam.setAttribute( "mass", ref_part.mass );
        beam.setAttribute( "charge", ref_part.charge );

    beam["positionOffset"]["x"].makeConstant( ref_part.x );
    beam["positionOffset"]["y"].makeConstant( ref_part.y );
    beam["positionOffset"]["t"].makeConstant( ref_part.t );

    }

    template<class PIt>
    void Save_impactx_PosID(PIt& pti,
                openPMD::ParticleSpecies& currSpecies,
                unsigned long long offset) const
    {
      //  use default is enough I guess?
      SavePosId(pti, currSpecies, offset);
    }

    void
    GetNames(amrex::Vector<std::string>& real_names,
         amrex::Vector<std::string>& int_names,
         amrex::Vector<int>& int_flags,
         amrex::Vector<int>& real_flags) const

    {
      real_names.resize(RealSoA::names_s.size());
      std::copy(RealSoA::names_s.begin(), RealSoA::names_s.end(), real_names.begin());

      for (auto real_idx=0; real_idx < RealSoA::nattribs; real_idx++) {
           auto const component_name = real_names.at(real_idx);
      }
      real_flags.resize(real_names.size(), 1);


      // currently no use of int names
      int_names.resize(0);
      int_flags.resize(0);

      /*
    // enable when IntSoA is in
      amrex::Print() <<" => [check] int names ?  "<<IntSoA::names_s.size()<<std::endl;
      int_names.resize(IntSoA::names_s.size());
      std::copy(IntSoA::names_s.begin(), IntSoA::names_s.end(), int_names.begin());

      for (auto int_idx=0; int_idx < IntSoA::nattribs; int_idx++) {
    auto const component_name = int_names.at(int_idx);
      }
      int_flags.resize(int_names.size(), 1);
      */
    }

  };  // class impactxWriter




class AMReXWithOpenPMD
{
public:
  AMReXWithOpenPMD();

  AMReXWithOpenPMD (AMReXWithOpenPMD && other) = default;
  AMReXWithOpenPMD& operator= (AMReXWithOpenPMD && other) = default;

  ~AMReXWithOpenPMD();

  void SetWriter(amrex::openpmd_api::AMReX_openPMDWriter* w);

  bool InitLocalHandler(const std::string& prefix);

  std::unique_ptr<amrex::openpmd_api::AMReX_openPMDHandler> m_UserHandler;

private:
  std::string m_Prefix = "";
};




    /** This element writes the particle beam out to openPMD data.
     *
     * This class behaves like a singleton if constructed with the
     * same series name as an existing instance.
     */
    struct BeamPlotplus
    : public elements::Thin
    {
        static constexpr auto name = "BeamPlotplus";
        using PinnedContainer = typename ImpactXParticleContainer::ContainerLike<amrex::PinnedArenaAllocator>;

        /** This element writes the particle beam out to openPMD data.
         *
         * Elements with the same series name are identical.
         *
         * @param series_name name of the data series, usually the element name
         * @param backend file format backend for openPMD, e.g., "bp" or "h5"
         * @param encoding openPMD iteration encoding: "v"ariable based, "f"ile based, "g"roup based (default)
         */
        BeamPlotplus (std::string series_name, std::string backend="default", std::string encoding="g");

        BeamPlotplus (BeamPlotplus const & other) = default;
        BeamPlotplus (BeamPlotplus && other) = default;
        BeamPlotplus& operator= (BeamPlotplus const & other) = default;
        BeamPlotplus& operator= (BeamPlotplus && other) = default;


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

        /** This does nothing to the reference particle. */
        using Thin::operator();

        /** Close and deallocate all data series and backends.
         */
        void
        finalize ();

    private:
      std::string m_seriesName;
      void*  m_plotWriter = NULL;
      static inline std::map<std::string, void*> m_uniqueWriter = {};

      //int m_file_min_digits = 6; //! minimum number of digits to iteration number in file name
    };

} // namespace impactx::diagnostics

#endif // AMREX_USE_OPENPMD_API
#endif // IMPACTX_ELEMENTS_DIAGS_OPENPMD_WITH_AMREX_H
