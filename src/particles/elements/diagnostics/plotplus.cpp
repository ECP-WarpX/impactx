/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */

/// #include "openPMD.H"
#include "ImpactXVersion.H"
#include "particles/ImpactXParticleContainer.H"

#include <ablastr/particles/IndexHandling.H>

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

#include <utility>

#ifdef AMREX_USE_OPENPMD_API
#include "plotplus.H"

namespace impactx::diagnostics
{

  class StepMgr
{
public:
  StepMgr(int step, AMReXWithOpenPMD* owner)
    :m_Step(step),
     m_Owner(owner)
  {
    m_Owner = owner;
    m_Owner->m_UserHandler->m_Writer->SetStep(m_Step);
  }
  ~StepMgr()
  {
    m_Owner->m_UserHandler->m_Writer->CloseStep(m_Step);
  }

private:
  int m_Step;
  AMReXWithOpenPMD* m_Owner;
};


  AMReXWithOpenPMD::AMReXWithOpenPMD()
  {
    // warpx has multiple diags, each should maintain its own handler
    m_UserHandler = amrex::openpmd_api::InitUserHandler(m_Prefix);
  }

  void AMReXWithOpenPMD::SetWriter(amrex::openpmd_api::AMReX_openPMDWriter* w)
  {
    BL_ASSERT ( m_UserHandler != nullptr );
    BL_ASSERT ( w != nullptr );

    // so the openpmd filepath assigned from input file is still in use
    w->m_openPMDPrefix = m_UserHandler->m_Writer->m_openPMDPrefix;
    w->m_openPMDEncoding = m_UserHandler->m_Writer->m_openPMDEncoding;
    w->m_openPMDFileType = m_UserHandler->m_Writer->m_openPMDFileType;
    w->m_openPMDSeriesOptions = m_UserHandler->m_Writer->m_openPMDSeriesOptions;

    m_UserHandler->m_Writer.reset(w);
  }


  AMReXWithOpenPMD::~AMReXWithOpenPMD()
{
  amrex::openpmd_api::CloseUserHandler(m_UserHandler);
}

bool AMReXWithOpenPMD::InitLocalHandler(const std::string& prefix)
{
  if (m_Prefix.compare(prefix) == 0)
    return false;

  amrex::Print()<<" openpmd_api::Init handler "<<m_Prefix<<" \n";
  m_Prefix = prefix;
  m_UserHandler = amrex::openpmd_api::InitUserHandler(prefix);
  return true;
}



  //
  // BeamMonitor, with I/O handled by amrex::openpmd-api
  //
    void BeamPlotplus::finalize ()
    {
      BL_PROFILE("BeamPlotplus::finalize()")
      if (m_uniqueWriter.count(m_seriesName) == 0u)
    return;

      if (m_plotWriter != NULL) {
         auto m_Writer = (AMReXWithOpenPMD*)(m_plotWriter);
         delete m_Writer;
      }
      amrex::Print()<<" => [check]: is things in BeamPlotplus::finalize() addressed? \n";
      m_uniqueWriter.erase(m_seriesName);

    }

    BeamPlotplus::BeamPlotplus (std::string series_name, std::string backend, std::string encoding)
      :m_seriesName(std::move(series_name))
    {
      BL_PROFILE("BeamPlotplus::BeamPlotplus()")
      if (m_uniqueWriter.count(m_seriesName) > 0u)
    {
      m_plotWriter = m_uniqueWriter[m_seriesName];
    }
      else
    {
      auto m_Writer =  new AMReXWithOpenPMD();
#ifdef ImpactX_USE_OPENPMD
        // encoding of iterations in the series
        openPMD::IterationEncoding series_encoding = openPMD::IterationEncoding::groupBased;
        if ( "v" == encoding )
            series_encoding = openPMD::IterationEncoding::variableBased;
        else if ( "g" == encoding )
            series_encoding = openPMD::IterationEncoding::groupBased;
        else if ( "f" == encoding )
            series_encoding = openPMD::IterationEncoding::fileBased;

    if ( m_Writer->InitLocalHandler(m_seriesName) )
      {
        AMReX_impactxWriter* testWriter = new AMReX_impactxWriter(series_encoding);
        m_Writer->SetWriter(testWriter);
      }
#else
        amrex::AllPrint() << "Warning: openPMD output requested but not compiled for series=" << m_series_name << "\n";
#endif
    m_plotWriter = m_Writer;
    m_uniqueWriter[m_seriesName] = m_plotWriter;
    }// else
    }

  #ifdef NEVER
    void BeamPlotplus::prepare (
        PinnedContainer & pc,
        RefPart const & ref_part,
        int step
    ) {
#ifdef ImpactX_USE_OPENPMD

    /* should be covered by amrex-openpmd-io
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
    */
#else
        amrex::ignore_unused(pc, step);
#endif
    }
#endif
    void
    BeamPlotplus::operator() (
        ImpactXParticleContainer & pc,
        int step
    )
    {
      BL_PROFILE("BeamPlotplus::(pc, step)")
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

    auto m_Writer = (AMReXWithOpenPMD*)(m_plotWriter);

    StepMgr sm(step, m_Writer);
    pinned_pc.CountParticles();


      AMReX_impactxWriter* impactxWriter =  (AMReX_impactxWriter*) (m_Writer->m_UserHandler->m_Writer.get());
      amrex::Vector<std::string> real_names;
      amrex::Vector<std::string> int_names;
      amrex::Vector<int> int_flags;
      amrex::Vector<int> real_flags;
      impactxWriter->GetNames(real_names, int_names, int_flags, real_flags);
      m_Writer->m_UserHandler->m_Writer->DumpParticles(pinned_pc,
                               "beam",
                               real_flags,
                               int_flags,
                               real_names,
                               int_names,

                               [=] ([[maybe_unused]] auto& ppc, openPMD::ParticleSpecies& currSpecies, [[maybe_unused]]  unsigned long long localTotal)
                               {
                                  impactxWriter->SetConstantRefPart(currSpecies,  ref_part);
                               },
                               [=] (auto& pti, openPMD::ParticleSpecies& currSpecies, unsigned long long offset)
                               {
                                   // use the default
                                   impactxWriter->Save_impactx_PosID(pti, currSpecies, offset);
                               });


    }
} // namespace impactx::diagnostics


#endif //#ifdef AMREX_USE_OPENPMD_API
