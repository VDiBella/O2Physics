// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file treeReaderDplusToPiKPi.cxx
/// \brief Writer of D+ → π+ K- π+ candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug, local optimization of analysis on small samples or ML training.
///        In this file are defined and filled the output tables
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

#include "PWGHF/CutPairDplus/DplusPairCut.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTreeReaderDplusToPiKPi {
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 1, "Selection Flag for Dplus"};

  using SelectedCandidates = soa::Filtered<o2::aod::HfCandDpFulls>;
  Filter filterSelectCandidates = aod::full::candidateSelFlag >= selectionFlagDplus;

  HistogramConfigSpec hTH1NCand{HistType::kTH1F, {{7, -0.5, 6.5}}};
  HistogramRegistry registry{
    "registry",
    {{"hNCand", "Number of D candidates per event;xxxx", hTH1NCand}}};

  void init(InitContext const&)
  {
    registry.add("hMassDplus", "D+ candidates;inv. mass (#pi#pi K) (GeV/#it{c}^{2}))", {HistType::kTH1F, {{120, 1.5848, 2.1848}}});
    registry.add("hMassDplusMatched", "D+ matched candidates;inv. mass (#pi#pi K) (GeV/#it{c}^{2}))", {HistType::kTH1F, {{120, 1.5848, 2.1848}}});
    registry.add("hMassDplusDminus", "D+D- pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});inv. mass (#pi K) (GeV/#it{c}^{2})", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {120, 1.5848, 2.1848}}});
  }

  void processLocalData(o2::aod::HfCandDpFullEvs::iterator const& localCollision,
                        SelectedCandidates const& localCandidates)
  {
    // LOG(info) << "BC: " << localCollision.bcId() << "  #cand: " << localCandidates.size();
    registry.fill(HIST("hNCand"), localCandidates.size());

    for (const auto& cand1 : localCandidates) {
      auto mass1 = cand1.m();
      for (auto cand2 = cand1 + 1; cand2 != localCandidates.end(); ++cand2) {
        auto mass2 = cand2.m();
        registry.fill(HIST("hMassDplusDminus"), mass2, mass1);
      }
    }
  }
  PROCESS_SWITCH(HfTreeReaderDplusToPiKPi, processLocalData, "Process local data", true);

  void processLocalDataMc(o2::aod::HfCandDpFullEvs::iterator const& localCollision,
                          SelectedCandidates const& localCandidates)
  {
    registry.fill(HIST("hNCand"), localCandidates.size());

    for (const auto& cand1 : localCandidates) {
      auto mass1 = cand1.m();
      registry.fill(HIST("hMassDplus"), mass1);
      if (std::abs(cand1.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)
        registry.fill(HIST("hMassDplusMatched"), mass1);
    }
  }
  PROCESS_SWITCH(HfTreeReaderDplusToPiKPi, processLocalDataMc, "Process local MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTreeReaderDplusToPiKPi>(cfgc)};
}
