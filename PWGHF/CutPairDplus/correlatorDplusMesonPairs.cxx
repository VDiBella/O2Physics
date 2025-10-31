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

/// \file correlatorDplusMesonPairs.cxx
/// \brief D+ correlator task
///
/// \author Valerio Di Bella (valerio.di.bella@cern.ch) IPHC Strasbourg

#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/DMesonPairsTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfCorrelatorDplusMesonPairs {

  // Produces<aod::DplusPair> entryDplusPair;
  // Produces<aod::DplusPairMcInfo> entryDplusPairMcInfo;
  Produces<aod::D0PairMcGenInfo> entryDplusPairMcGenInfo;

  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for Dplus"}; // 7 corresponds to topo+PID cuts

  Configurable<float> yCandMax{"yCandMax", 0.8, "maxmum |y| of Dplus candidates"};
  Configurable<float> ptCandMin{"ptCandMin", -1., "minimum pT of Dplus candidates"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots"};
  Configurable<bool> selectSignalRegionOnly{"selectSignalRegionOnly", false, "only use events close to PDG peak"};
  Configurable<float> massCut{"massCut", 0.05, "Maximum deviation from PDG peak allowed for signal region"};
  Configurable<bool> daughterTracksCutFlag{"daughterTracksCutFlag", false, "Flag to add cut on daughter tracks"};

  HfHelper hfHelper;

  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandDplusMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  using McParticles3Prong = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Filter filterDplusFlag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi))) != static_cast<uint8_t>(0);
  Filter filterDplusSelFlag = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  HistogramConfigSpec hTH1Pt{HistType::kTH1F, {{180, 0., 36.}}};
  HistogramConfigSpec hTH1Y{HistType::kTH1F, {{100, -5., 5.}}};
  HistogramConfigSpec hTH1NContrib{HistType::kTH1F, {{120, -0.5, 119.5}}};
  HistogramConfigSpec hTH1Phi{HistType::kTH1F, {{32, 0., o2::constants::math::TwoPI}}};
  HistogramConfigSpec hTH2Pid{HistType::kTH2F, {{500, 0., 10.}, {400, -20., 20.}}};
  HistogramConfigSpec hTH3PtVsYVsNContrib{HistType::kTH3F, {{360, 0., 36.}, {20, -1., 1.}, {120, -0.5, 119.5}}};

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "D meson candidates;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtProng0", "D meson candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtProng1", "D meson candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hEta", "D meson candidates;candidate #it{#eta};entries", hTH1Y},
     {"hPhi", "D meson candidates;candidate #it{#varphi};entries", hTH1Phi},
     {"hY", "D meson candidates;candidate #it{y};entries", hTH1Y},
     {"hPVContrib", "D meson candidates;candidate Number of PV contributors;entries", hTH1NContrib},
     // MC Gen plots
     {"hPtCandMcGen", "D meson candidates MC Gen;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hEtaMcGen", "D meson candidates MC Gen;candidate #it{#eta};entries", hTH1Y},
     {"hPhiMcGen", "D meson candidates MC Gen;candidate #it{#varphi};entries", hTH1Phi},
     {"hPtVsYVsNContribMcGen", "D meson candidates MC Gen;candidate #it{p}_{T} (GeV/#it{c});#it{y};Number of contributors", hTH3PtVsYVsNContrib},
     {"hPtVsYVsNContribMcGenPrompt", "D meson candidates MC Gen Prompt;candidate #it{p}_{T} (GeV/#it{c});#it{y};Number of contributors", hTH3PtVsYVsNContrib},
     {"hPtVsYVsNContribMcGenNonPrompt", "D meson candidates MC Gen Prompt;candidate #it{p}_{T} (GeV/#it{c});#it{y};Number of contributors", hTH3PtVsYVsNContrib},
     {"hNContribMcGen", "D meson candidates MC Gen;Number of PV contributors", hTH1NContrib}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;

    constexpr int kNBinsMatching = 8;
    std::string labelsMatching[kNBinsMatching];
    // Cand1 analysis
    labelsMatching[0] = "total # of Cand 1";
    labelsMatching[1] = "# of matched D Cand 1";
    labelsMatching[2] = "# of matched Dbar Cand 1";
    labelsMatching[3] = "# of unmatched Cand 1";
    // Cand2 analysis
    labelsMatching[4] = "total # of Cand 2";
    labelsMatching[5] = "# of matched D Cand 2";
    labelsMatching[6] = "# of matched Dbar Cand 2";
    labelsMatching[7] = "# of unmatched Cand 2";

    AxisSpec axisMatching = {kNBinsMatching, 0.5, kNBinsMatching + 0.5, ""};
    registry.add("hMatchingMcRec", "D Meson candidates; MC matching status;entries", HistType::kTH1F, {axisMatching});
    registry.add("hMatchingMcGen", "D Meson candidates; MC matching status;entries", HistType::kTH1F, {axisMatching});

    for (int iBin = 0; iBin < kNBinsMatching; iBin++) {
      registry.get<TH1>(HIST("hMatchingMcRec"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMatching[iBin].data());
      registry.get<TH1>(HIST("hMatchingMcGen"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMatching[iBin].data());
    }

    registry.add("hMassDplus", "D+ candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDminus", "D- candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDplusDplus", "D+D+ pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});inv. mass (#pi K) (GeV/#it{c}^{2})", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {120, 1.5848, 2.1848}}});
    registry.add("hMassDminusDminus", "D-D- pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});inv. mass (#pi K) (GeV/#it{c}^{2})", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {120, 1.5848, 2.1848}}});
    registry.add("hMassDplusDminus", "D+D- pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});inv. mass (#pi K) (GeV/#it{c}^{2})", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {120, 1.5848, 2.1848}}});

    // Associated with MC
    registry.add("hMassDplusMc", "D+ associated candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDminusMc", "D- associated candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void processData(aod::Collision const& collision, CandDplusData const& candidates, aod::Tracks const&)
  {
    for (const auto& candidate1 : candidates) {
      if (std::abs(hfHelper.yDplus(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }

      auto prong0Cand1 = candidate1.template prong0_as<aod::Tracks>();
      auto prong1Cand1 = candidate1.template prong1_as<aod::Tracks>();
      auto prong2Cand1 = candidate1.template prong2_as<aod::Tracks>();

      auto mass1 = hfHelper.invMassDplusToPiKPi(candidate1);
      bool isSignalDplusCand1 = std::abs(mass1 - MassDPlus) < massCut;
      if (selectSignalRegionOnly && !isSignalDplusCand1) {
        continue;
      }

      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), candidate1.y(MassDPlus));
      registry.fill(HIST("hPVContrib"), collision.numContrib());

      bool Dplus1 = prong1Cand1.sign() < 0;
      if (Dplus1) { // D+
        registry.fill(HIST("hMassDplus"), mass1, candidate1.pt());
      } else { // D-
        registry.fill(HIST("hMassDminus"), mass1, candidate1.pt());
      }

      for (auto candidate2 = candidate1 + 1; candidate2 != candidates.end(); ++candidate2) {
        if (std::abs(hfHelper.yDplus(candidate2)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate2.pt() < ptCandMin) {
          continue;
        }

        auto prong0Cand2 = candidate2.template prong0_as<aod::Tracks>();
        auto prong1Cand2 = candidate2.template prong1_as<aod::Tracks>();
        auto prong2Cand2 = candidate2.template prong2_as<aod::Tracks>();
        if (daughterTracksCutFlag) {
          if (prong0Cand1 == prong0Cand2)
            continue;
          if (prong0Cand1 == prong1Cand2)
            continue;
          if (prong0Cand1 == prong2Cand2)
            continue;
          if (prong1Cand1 == prong0Cand2)
            continue;
          if (prong1Cand1 == prong1Cand2)
            continue;
          if (prong1Cand1 == prong2Cand2)
            continue;
          if (prong2Cand1 == prong0Cand2)
            continue;
          if (prong2Cand1 == prong1Cand2)
            continue;
          if (prong2Cand1 == prong2Cand2)
            continue;
        }

        auto mass2 = hfHelper.invMassDplusToPiKPi(candidate2);
        bool isSignalDplusCand2 = std::abs(mass2 - MassDPlus) < massCut;
        if (selectSignalRegionOnly && !isSignalDplusCand2) {
          continue;
        }

        bool Dplus2 = prong1Cand2.sign() < 0;
        if (Dplus2) { // D+
          registry.fill(HIST("hMassDplus"), mass2, candidate2.pt());
          if (Dplus1) {
            registry.fill(HIST("hMassDplusDplus"), mass2, mass1);
          } else {
            registry.fill(HIST("hMassDplusDminus"), mass2, mass1);
          }
        } else { // D-
          registry.fill(HIST("hMassDminus"), mass2, candidate2.pt());
          if (!Dplus1) {
            registry.fill(HIST("hMassDminusDminus"), mass2, mass1);
          } else {
            registry.fill(HIST("hMassDplusDminus"), mass2, mass1);
          }
        }

      } // end inner loop (Cand2)
    } // end outer loop (Cand1)
  }
  PROCESS_SWITCH(HfCorrelatorDplusMesonPairs, processData, "Process data mode", true);

  void processMcRec(aod::Collision const&, CandDplusMcReco const& candidates, aod::TracksWMc const&,
                    aod::McParticles const& mcParticles)
  {
    for (const auto& candidate1 : candidates) {
      if (std::abs(hfHelper.yDplus(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }

      auto prong0Cand1 = candidate1.template prong0_as<aod::TracksWMc>();
      auto prong1Cand1 = candidate1.template prong1_as<aod::TracksWMc>();
      auto prong2Cand1 = candidate1.template prong2_as<aod::TracksWMc>();

      auto mass1 = hfHelper.invMassDplusToPiKPi(candidate1);
      bool isSignalDplusCand1 = std::abs(mass1 - MassDPlus) < massCut;
      if (selectSignalRegionOnly && !isSignalDplusCand1) {
        continue;
      }

      bool Dplus1 = prong1Cand1.sign() < 0;
      if (Dplus1) { // D+
        registry.fill(HIST("hMassDplus"), mass1, candidate1.pt());
      } else { // D-
        registry.fill(HIST("hMassDminus"), mass1, candidate1.pt());
      }

      bool isTrueDCand1 = std::abs(candidate1.flagMcMatchRec()) == 1 << o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi;
      if (isTrueDCand1) {
        auto indexMother = RecoDecay::getMother(mcParticles, prong1Cand1.mcParticle(), o2::constants::physics::Pdg::kDPlus, true);
        auto motherD = mcParticles.rawIteratorAt(indexMother);
        // LOG(info) << motherD.pdgCode() << ' ' << motherD.pz() << ' ' << candidate1.pz() << " !!!!!";

        if (Dplus1) { // D+
          registry.fill(HIST("hMassDplusMc"), mass1, motherD.pt());
        } else { // D-
          registry.fill(HIST("hMassDminusMc"), mass1, motherD.pt());
        }
      }

      for (auto candidate2 = candidate1 + 1; candidate2 != candidates.end(); ++candidate2) {
        if (std::abs(hfHelper.yDplus(candidate2)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate2.pt() < ptCandMin) {
          continue;
        }

        auto prong0Cand2 = candidate2.template prong0_as<aod::TracksWMc>();
        auto prong1Cand2 = candidate2.template prong1_as<aod::TracksWMc>();
        auto prong2Cand2 = candidate2.template prong2_as<aod::TracksWMc>();
        if (daughterTracksCutFlag) {
          if (prong0Cand1 == prong0Cand2)
            continue;
          if (prong0Cand1 == prong1Cand2)
            continue;
          if (prong0Cand1 == prong2Cand2)
            continue;
          if (prong1Cand1 == prong0Cand2)
            continue;
          if (prong1Cand1 == prong1Cand2)
            continue;
          if (prong1Cand1 == prong2Cand2)
            continue;
          if (prong2Cand1 == prong0Cand2)
            continue;
          if (prong2Cand1 == prong1Cand2)
            continue;
          if (prong2Cand1 == prong2Cand2)
            continue;
        }

        auto mass2 = hfHelper.invMassDplusToPiKPi(candidate2);
        bool isSignalDplusCand2 = std::abs(mass2 - MassDPlus) < massCut;
        if (selectSignalRegionOnly && !isSignalDplusCand2) {
          continue;
        }

        bool Dplus2 = prong1Cand2.sign() < 0;
        if (Dplus2) { // D+
          registry.fill(HIST("hMassDplus"), mass2, candidate2.pt());
          if (Dplus1) {
            registry.fill(HIST("hMassDplusDplus"), mass2, mass1);
          } else {
            registry.fill(HIST("hMassDplusDminus"), mass2, mass1);
          }
        } else { // D-
          registry.fill(HIST("hMassDminus"), mass2, candidate2.pt());
          if (!Dplus1) {
            registry.fill(HIST("hMassDminusDminus"), mass2, mass1);
          } else {
            registry.fill(HIST("hMassDplusDminus"), mass2, mass1);
          }
        }
        /*
              bool isTrueDCand2 = std::abs(candidate2.flagMcMatchRec()) == 1 << o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi;
              if (isTrueDCand2) {
                if (Dplus2) { // D+
                  registry.fill(HIST("hMassDplusMc"), mass2, candidate2.pt());
                } else { // D-
                  registry.fill(HIST("hMassDminusMc"), mass2, candidate2.pt());
                }
              }
        */
      } // end inner loop (Cand2)
    } // end outer loop (Cand1)
  }
  PROCESS_SWITCH(HfCorrelatorDplusMesonPairs, processMcRec, "Process Mc reco mode", false);

  void processMcGen(aod::McCollision const&, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, McParticles3Prong const& mcParticles)
  {
    int numPvContributorsGen{0};
    for (const auto& collision : collisions) { // loop over reco collisions associated to this gen collision
      int numPvContributors = collision.numContrib();

      if (numPvContributors > numPvContributorsGen) { // we take the associated reconstructed collision with higher number of PV contributors
        numPvContributorsGen = numPvContributors;
      }
    }

    for (const auto& particle1 : mcParticles) {
      // check if the particle is Dplus or Dminus
      if (std::abs(particle1.pdgCode()) != Pdg::kDPlus) {
        continue;
      }

      if (std::abs(particle1.y()) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle1.pt() < ptCandMin) {
        continue;
      }

      registry.fill(HIST("hEtaMcGen"), particle1.eta());
      registry.fill(HIST("hPhiMcGen"), particle1.phi());
      registry.fill(HIST("hPtCandMcGen"), particle1.pt());

      // check if it's MC matched
      int8_t matchedGen1 = particle1.flagMcMatchGen();
      // check origin
      int8_t originGen1 = particle1.originMcGen();

      registry.fill(HIST("hPtVsYVsNContribMcGen"), particle1.pt(), particle1.y(), numPvContributorsGen);
      if (originGen1 == 1) {
        registry.fill(HIST("hPtVsYVsNContribMcGenPrompt"), particle1.pt(), particle1.y(), numPvContributorsGen);
      }
      if (originGen1 == 2) {
        registry.fill(HIST("hPtVsYVsNContribMcGenNonPrompt"), particle1.pt(), particle1.y(), numPvContributorsGen);
      }
      registry.fill(HIST("hNContribMcGen"), numPvContributorsGen);

      for (auto particle2 = particle1 + 1; particle2 != mcParticles.end(); ++particle2) {
        // check if the particle is Dplus or Dminus
        if (std::abs(particle2.pdgCode()) != Pdg::kDPlus) {
          continue;
        }
        if (std::abs(particle2.y()) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle2.pt() < ptCandMin) {
          continue;
        }

        // check if it's MC matched
        int8_t matchedGen2 = particle2.flagMcMatchGen();
        // check origin
        int8_t originGen2 = particle2.originMcGen();

        // Fill hMatchingMcGen - Cand 1
        registry.fill(HIST("hMatchingMcGen"), 1);
        if (matchedGen1 == 1) {
          registry.fill(HIST("hMatchingMcGen"), 2);
        } else if (matchedGen1 == -1) {
          registry.fill(HIST("hMatchingMcGen"), 3);
        } else if (matchedGen1 == 0) {
          registry.fill(HIST("hMatchingMcGen"), 4);
        }
        // Fill hMatchingMcRec - Cand 2
        registry.fill(HIST("hMatchingMcGen"), 5);
        if (matchedGen2 == 1) {
          registry.fill(HIST("hMatchingMcGen"), 6);
        } else if (matchedGen2 == -1) {
          registry.fill(HIST("hMatchingMcGen"), 7);
        } else if (matchedGen2 == 0) {
          registry.fill(HIST("hMatchingMcGen"), 8);
        }

        // Fill tables
        entryDplusPairMcGenInfo(originGen1, originGen2, matchedGen1, matchedGen2);

      } // end inner loop
    } // end outer loop
  }
  PROCESS_SWITCH(HfCorrelatorDplusMesonPairs, processMcGen, "Process Mc Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDplusMesonPairs>(cfgc)};
}
