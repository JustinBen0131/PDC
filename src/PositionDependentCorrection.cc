#include "PositionDependentCorrection.h"
#include "PDCGeo.h"
#include <numeric>
#include <cmath>
#include <fstream>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <TKey.h>
#include <array>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <regex>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerDefs.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <mbd/BbcGeom.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtHit.h>
#include <TFile.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TH2.h>
#include <TH3.h>
#include <TNtuple.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TTree.h>
#include <Event/Event.h>
#include <Event/packet.h>
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include"TRandom3.h"
#include <algorithm>
#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>
#include <map>

#include <g4main/PHG4TruthInfoContainer.h>
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/BEmcRecCEMC.h"

constexpr const char* ANSI_BOLD  = "\033[1m";
constexpr const char* ANSI_RESET = "\033[0m";
constexpr const char* ANSI_GREEN = "\033[32m";
constexpr const char* ANSI_YELLOW= "\033[33m";
constexpr const char* ANSI_CYAN    = "\033[36m";
constexpr const char* ANSI_RED   = "\033[31m";
constexpr const char* ANSI_MAGENTA = "\033[35m";

R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)

using namespace std;
namespace CLHEP { class Hep3Vector; }
std::atomic<uint64_t> PositionDependentCorrection::s_verbosityLevel{0};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PositionDependentCorrection::PositionDependentCorrection(const std::string &name,
                                                         const std::string &filename)
  : SubsysReco(name)
  , detector("HCALIN")
  , outfilename(filename)
  , m_isSimulation(false)
  , g4hitntuple(nullptr)
  , g4cellntuple(nullptr)
  , towerntuple(nullptr)
  , clusterntuple(nullptr)
{
  _eventcounter = 0;
  /* vertex‑Z limits                                                     */
  m_vzTightCut  = 10.f;           // |z| ≤ 10 cm  → “physics” histograms
  m_vzSliceMax  = vzEdge.back();  // |z| ≤ 30 cm  → Δη(E,vz) spectra only
  m_nWinRAW = m_nWinCP = m_nWinBCorr = 0;
  s_verbosityLevel.store( Verbosity() );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PositionDependentCorrection::~PositionDependentCorrection()
{
  delete hm;

  if (g4hitntuple)   delete g4hitntuple;
  if (g4cellntuple)  delete g4cellntuple;
  if (towerntuple)   delete towerntuple;
  if (clusterntuple) delete clusterntuple;
}


// ===================================================================
// EMCal energy resolution  (≈10%/√E ⊕ 2%)
// Returned value has the same units as E.
// ===================================================================
static inline double sigmaE(double E)
{
  const double stoch    = 0.10 / std::sqrt(E);   // 10%/√E
  const double constant = 0.02;                  // 2 %
  return E * std::hypot(stoch, constant);        // quadrature
}

// ------------------------------------------------------------------
// (helper) translate a photon/cluster energy into a slice index
// ------------------------------------------------------------------
int PositionDependentCorrection::getEnergySlice(float E) const
{
  const int nSlices = N_Ebins;          // 8 slices

  /* -------------------------------------------------------------
   * 1)  Regular RANGE mode  (unchanged)
   * ----------------------------------------------------------- */
  if (m_binningMode == EBinningMode::kRange)
  {
    for (int i = 0; i < nSlices; ++i)
      if (E >= eEdge[i] && E < eEdge[i + 1]) return i;
    return -1;                               // outside table
  }

  /* -------------------------------------------------------------
   * 2)  DISCRETE mode  (resolution‑scaled window)
   * ----------------------------------------------------------- */
  int    bestIdx  = -1;
  double bestDiff = 1e9;

  for (int i = 0; i < nSlices; ++i)
  {
    const double diff = std::fabs(E - eEdge[i]);   // centre distance
    if (diff < bestDiff) { bestDiff = diff; bestIdx = i; }
  }

  /* dynamic tolerance:  max( Nσ·σ_E , floor ) */
  const double dynTol =
      std::max(kTolFactor * sigmaE(eEdge[bestIdx]), kTolMinGeV);

  /* optional runtime printout to check efficiency */
  if (Verbosity() > 2 && m_binningMode == EBinningMode::kDiscrete)
    std::cout << "[getEnergySlice]  E=" << E
              << "  centre=" << eEdge[bestIdx]
              << "  diff="   << bestDiff
              << "  tol="    << dynTol
              << "  -> "     << (bestDiff < dynTol ? "ACCEPT" : "REJECT")
              << '\n';

  return (bestDiff < dynTol) ? bestIdx : -1;
}

/* ------------------------------------------------------------------------- */
/*  PositionDependentCorrection  –  histogram booking                        */
/* ------------------------------------------------------------------------- */

/** Book all histograms that are needed for _both_ data and simulation */
void
PositionDependentCorrection::bookCommonHistograms
(const std::function<std::string(int)>& makeLabel)
{
  /* ------------------------------------------------------------------ */
  /* 1)  “Always‑on” histograms (identical to the original listing)      */
  /* ------------------------------------------------------------------ */
  h_mass_eta_lt      = new TH2F("h_mass_eta_lt",      "", 50, 0, 0.5, 96, 0, 96);
  h_mass_eta_lt_rw   = new TH2F("h_mass_eta_lt_rw",   "", 50, 0, 0.5, 96, 0, 96);
  h_pt_eta           = new TH2F("h_pt_eta",           "", 100, 0, 10, 96, 0, 96);
  h_pt_eta_rw        = new TH2F("h_pt_eta_rw",        "", 100, 0, 10, 96, 0, 96);
  h_cemc_etaphi      = new TH2F("h_cemc_etaphi",      "", 96, 0, 96, 256, 0, 256);
  h_InvMass          = new TH1F("h_InvMass",          "Invariant Mass", 500, 0, 1.0);
  h_InvMass_w        = new TH1F("h_InvMass_w",        "Invariant Mass", 500, 0, 1.0);
  h_InvMassMix       = new TH1F("h_InvMassMix",       "Invariant Mass", 120, 0, 1.2);
    
  // ------------------------------------------------------------------
  // π0 invariant‑mass vs photon‑energy slice      (2 × TH2F)
  // Y‑axis uses EXACTLY the edges from eEdge[9].
  // ------------------------------------------------------------------
  h_mE_raw  = new TH2F("h_mE_raw",
                         "π^{0} mass vs E_{#gamma}^{max} – RAW;"
                         "M_{γγ} [GeV];E_{#gamma}^{max} [GeV]",
                         120, 0.00, 0.30,         // X: mass
                         N_Ebins, eEdge);         // Y: slice edges 2…30 GeV

  h_mE_corr = new TH2F("h_mE_corr",
                         "π^{0} mass vs E_{#gamma}^{max} – φ‑corr.;"
                         "M_{γγ} [GeV];E_{#gamma}^{max} [GeV]",
                         120, 0.00, 0.30,
                         N_Ebins, eEdge);

  hm->registerHisto(h_mE_raw );
  hm->registerHisto(h_mE_corr);

  h_m_blk_raw  = new TH3F("h_m_blk_raw",
                            "π^{0} mass – RAW;η_{loc};φ_{loc};M_{γγ} [GeV]",
                            14, -0.5, 1.5,     // X: ηloc
                            14, -0.5, 1.5,     // Y: φloc
                            60, 0.00, 0.30);   // Z: mass

  h_m_blk_corr = new TH3F("h_m_blk_corr",
                            "π^{0} mass – φ‑corr.;η_{loc};φ_{loc};M_{γγ} [GeV]",
                            14, -0.5, 1.5,
                            14, -0.5, 1.5,
                            60, 0.00, 0.30);

  hm->registerHisto(h_m_blk_raw );
  hm->registerHisto(h_m_blk_corr);
    
    
  h_tower_e          = new TH1F("h_tower_e",          "", 1000,-1,5);
  h_etaphi_clus      = new TH2F("h_etaphi_clus",      "", 140,-1.2,1.2,64,-TMath::Pi(),TMath::Pi());
  h_clusE            = new TH1F("h_clusE",            "", 100, 0, 10);
  h_clusE_nTow       = new TH2F("h_clusE_nTow",       "", 20,0,20, 50,0,50);
  h_emcal_e_eta      = new TH1F("h_emcal_e_eta",      "", 96, 0, 96);
  h_pt1              = new TH1F("h_pt1",              "", 100, 0, 5);
  h_pt2              = new TH1F("h_pt2",              "", 100, 0, 5);
  h_nclusters        = new TH1F("h_nclusters",        "", 100, 0, 100);
  h_matched_res      = new TH2F("h_matched_res",      "", 100,0,1.5,20,-1,1);
  h_res_e            = new TH2F("h_res_e",            "", 100,0,1.5,20,0,20);
  h_res_e_phi        = new TH3F("h_res_e_phi",        "", 100,0,1.5,10,0,20,256,0,256);
  h_res_e_eta        = new TH3F("h_res_e_eta",        "", 300,0,1.5,40,0,20,96,0,96);
  h_m_pt_eta         = new TH3F("h_m_pt_eta",         "", 70,0,0.7,10,0,10,96,0,96);
  h_m_ptTr_eta       = new TH3F("h_m_ptTr_eta",       "", 70,0,0.7,10,0,10,96,0,96);
  h_m_ptTr_eta_trKin = new TH3F("h_m_ptTr_eta_trKin", "", 70,0,0.7,10,0,10,96,0,96);
  h_res              = new TH1F("h_res",              "", 50, 0, 1.5);
  h_delEta_e_eta     = new TH3F("h_delEta_e_eta",     "", 100,-0.1,0.1,10,0,20,96,0,96);
  h_delR_e_eta       = new TH3F("h_delR_e_eta",       "", 100,-0.1,0.1,10,0,20,96,0,96);
  h_delPhi_e_eta     = new TH3F("h_delPhi_e_eta",     "", 100,-0.3,0.3,20,0,20,96,0,96);
  h_delPhi_e_phi     = new TH3F("h_delPhi_e_phi",     "", 100,-0.1,0.1,20,0,20,256,0,256);
  pr_eta_shower      = new TProfile("pr_eta_shower",  "", 96,-48.5,47.5,-1,1.5);
  pr_phi_shower      = new TProfile("pr_phi_shower",  "", 256,-128.5,127.5,-1,1.5);
  h_vert_xy          = new TH2F("h_vert_xy",          "", 500,-120,120,500,-120,120);

  /* χ²‑maps (kept verbatim) */
  h2_chi2_tot_etaPhi = new TH2F(
        "h2_chi2_tot_etaPhi",
        "Clusters BEFORE #chi^{2} cut;Tower #eta index;Tower #varphi index;Counts",
        96, 0, 96,
        256, 0, 256);

  h2_chi2_rej_etaPhi = new TH2F(
        "h2_chi2_rej_etaPhi",
        "Clusters REJECTED by #chi^{2} cut;Tower #eta index;Tower #varphi index;Counts",
        96, 0, 96,
        256, 0, 256);

  p_chi2_pass_etaPhi = new TProfile2D(
        "p_chi2_pass_etaPhi",
        "Pass fraction after #chi^{2} cut;Tower #eta index;Tower #varphi index;⟨pass⟩",
        96, 0, 96,
        256, 0, 256,
        -0.1, 1.1);

  hm->registerHisto(h2_chi2_tot_etaPhi);
  hm->registerHisto(h2_chi2_rej_etaPhi);
  hm->registerHisto(p_chi2_pass_etaPhi);

  /* per‑block histograms (unchanged) */
  for (int ie = 0; ie < NBinsBlock; ++ie)
    for (int ip = 0; ip < NBinsBlock; ++ip)
    {
      h_mass_block_pt[ie][ip] = new TH2F(Form("h_mass_block_%d_%d_pt",ie,ip),
                                         "", 100,0,1, 5,0,10);
      h_res_block_E [ie][ip]  = new TH2F(Form("h_res_block_%d_%d_E",ie,ip),
                                         "", 120,0,1.2,5,0,10);
    }

  /* block‑coordinate × energy TH3 (identical code) */
  Double_t xEdges[15], yEdges[15];
  { const double s = 2.0/14.0; for (int i=0;i<=14;++i){ xEdges[i]=-0.5+i*s; yEdges[i]=-0.5+i*s; } }

  static constexpr Double_t eEdges[9] = {2,4,6,8,10,12,15,20,30};
  const char* modeTag = (m_binningMode == EBinningMode::kRange ? "range" : "disc");

  h3_cluster_block_cord_E = new TH3F(
        Form("h3_blockCoord_E_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Uncorrected local block coords vs. E_{slice}" :
           "Uncorrected local block coords vs. E_{centre}"),
        14,xEdges, 14,yEdges, 8,eEdges);

  h3_cluster_block_cord_E_corrected = new TH3F(
        Form("h3_blockCoord_Ecorr_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Corrected local block coords vs. E_{slice}" :
           "Corrected local block coords vs. E_{centre}"),
        14,xEdges, 14,yEdges, 8,eEdges);

  h_block_bin            = new TH1F("h_block_bin","",14,-0.5,1.5);
  pr_phi_vs_blockcoord   = new TProfile("pr_phi_vs_blockcoord","",14,-0.5,1.5,-0.2,0.2);

  /* per‑slice corrected local η/φ (unchanged) */
  for (int i=0;i<N_Ebins;++i)
  {
    const float eLo = eEdge[i], eHi = eEdge[i+1];
    const std::string lab = makeLabel(i);

    h_localPhi_corrected[i] = new TH1F(
          Form("h_localPhi_corr_%s",lab.c_str()),
          (m_binningMode==EBinningMode::kRange ?
             Form("Corrected local #phi;%.0f < E < %.0f GeV",eLo,eHi) :
             Form("Corrected local #phi;E = %.0f GeV",eLo)),
          50,-0.5,0.5);

    h_localEta_corrected[i] = new TH1F(
          Form("h_localEta_corr_%s",lab.c_str()),
          (m_binningMode==EBinningMode::kRange ?
             Form("Corrected local #eta;%.0f < E < %.0f GeV",eLo,eHi) :
             Form("Corrected local #eta;E = %.0f GeV",eLo)),
          50,-0.5,0.5);
  }

  /* vertex‑Z (always booked) */
  h_reco_vz = new TH1F("h_reco_vz","Reco Vertex Z;z_{reco} (cm);Counts",200,-100,100);
  hm->registerHisto(h_reco_vz);
}

/** Extra histograms that are **only** needed when m_isSimulation == true */
void
PositionDependentCorrection::bookSimulationHistograms
(const std::function<std::string(int)>& makeLabel)
{
  /* ---------------------- truth & matching ------------------------- */
  h_truth_eta      = new TH1F("h_truth_eta", "", 100,-1.2,1.2);
  h_truth_e        = new TH1F("h_truth_e",   "", 100, 0, 10);
  h_delR_recTrth   = new TH1F("h_delR_recTrth","",500,0,5);

  h_truth_vz       = new TH1F("h_truth_vz",
                               "Truth Vertex Z;z_{truth} (cm);Counts",200,-100,100);
  hm->registerHisto(h_truth_vz);

  h2_truthReco_vz  = new TH2F("h2_truthReco_vz",
                               "Truth vs Reco Vertex Z;z_{truth} (cm);z_{reco} (cm)",
                               200,-100,100,200,-100,100);
  hm->registerHisto(h2_truthReco_vz);

  h_truthE         = new TH1F("h_truthE","",10000,0,30);

    /* ──────────────────────────────────────────────────────────────────────
       (I)  Per‑E‑slice histograms that existed before
       (II) NEW  Δη(E, vz) histograms
       ──────────────────────────────────────────────────────────────────── */
    for (int i = 0; i < N_Ebins; ++i)
    {
      const float        eLo = eEdge[i], eHi = eEdge[i+1];
      const std::string  lab = makeLabel(i);

      /* -------- RAW / CORR Δφ, Δη (unchanged) --------------------------- */
      h_phi_diff_raw_E[i] = new TH1F(
            Form("h_phi_diff_raw_%s", lab.c_str()),
            (m_binningMode==EBinningMode::kRange ?
               Form("#Delta#phi raw;%.0f < E < %.0f GeV", eLo, eHi) :
               Form("#Delta#phi raw;E = %.0f GeV",        eLo)),
            200, -0.1, 0.1);

      h_eta_diff_raw_E[i] = new TH1F(
            Form("h_eta_diff_raw_%s", lab.c_str()),
            (m_binningMode==EBinningMode::kRange ?
               Form("#Delta#eta raw;%.0f < E < %.0f GeV", eLo, eHi) :
               Form("#Delta#eta raw;E = %.0f GeV",        eLo)),
            200, -0.1, 0.1);

      h_phi_diff_corrected_E[i] = new TH1F(
            Form("h_phi_diff_corr_%s", lab.c_str()),
            (m_binningMode==EBinningMode::kRange ?
               Form("#Delta#phi corrected;%.0f < E < %.0f GeV", eLo, eHi) :
               Form("#Delta#phi corrected;E = %.0f GeV",        eLo)),
            200, -0.1, 0.1);

      h_eta_diff_corrected_E[i] = new TH1F(
            Form("h_eta_diff_corr_%s", lab.c_str()),
            (m_binningMode==EBinningMode::kRange ?
               Form("#Delta#eta corrected;%.0f < E < %.0f GeV", eLo, eHi) :
               Form("#Delta#eta corrected;E = %.0f GeV",        eLo)),
            200, -0.1, 0.1);

      hm->registerHisto(h_phi_diff_raw_E[i]);
      hm->registerHisto(h_eta_diff_raw_E[i]);
      hm->registerHisto(h_phi_diff_corrected_E[i]);
      hm->registerHisto(h_eta_diff_corrected_E[i]);

      /* -------- CP helper histograms (needed by code below) ------------- */
      h_phi_diff_cpRaw_E [i] = new TH1F(
            Form("h_phi_diff_cpRaw_%s",  lab.c_str()),
            "#Delta#phi RAW‑CP;#Delta#phi (rad);Counts",   200, -0.1, 0.1);
      h_phi_diff_cpCorr_E[i] = new TH1F(
            Form("h_phi_diff_cpCorr_%s", lab.c_str()),
            "#Delta#phi CP‑corr;#Delta#phi (rad);Counts",  200, -0.1, 0.1);
      h_phi_diff_cpBcorr_E[i] = new TH1F(
            Form("h_phi_diff_cpBcorr_%s", lab.c_str()),
            "#Delta#phi b‑corr;#Delta#phi (rad);Counts",   200, -0.1, 0.1);

      h_eta_diff_cpRaw_E [i] = new TH1F(
            Form("h_eta_diff_cpRaw_%s",  lab.c_str()),
            "#Delta#eta RAW‑CP;#Delta#eta;Counts",         200, -0.1, 0.1);
      h_eta_diff_cpCorr_E[i] = new TH1F(
            Form("h_eta_diff_cpCorr_%s", lab.c_str()),
            "#Delta#eta CP‑corr;#Delta#eta;Counts",        200, -0.1, 0.1);
      h_eta_diff_cpBcorr_E[i] = new TH1F(
            Form("h_eta_diff_cpBcorr_%s", lab.c_str()),
            "#Delta#eta b‑corr;#Delta#eta;Counts",         200, -0.1, 0.1);

      hm->registerHisto(h_phi_diff_cpRaw_E [i]);
      hm->registerHisto(h_phi_diff_cpCorr_E[i]);
      hm->registerHisto(h_phi_diff_cpBcorr_E[i]);
      hm->registerHisto(h_eta_diff_cpRaw_E [i]);
      hm->registerHisto(h_eta_diff_cpCorr_E[i]);
      hm->registerHisto(h_eta_diff_cpBcorr_E[i]);

      /* ------------------------------------------------------------------ *
       * (II)  NEW  Δη(E, vz) histograms – one extra loop over vz slices
       * ------------------------------------------------------------------ */
      for (int iVz = 0; iVz < N_VzBins; ++iVz)
      {
        const float vzLo = vzEdge[iVz];
        const float vzHi = vzEdge[iVz+1];
        const std::string vzTag = Form("vz%.0f_%.0f", vzLo, vzHi);

        h_eta_diff_raw_E_vz[i][iVz] = new TH1F(
              Form("h_eta_diff_raw_%s_%s",  lab.c_str(), vzTag.c_str()),
              Form("#Delta#eta raw;%.0f<E<%.0f GeV, %.0f<|z_{vtx}|<%.0f cm",
                   eLo, eHi, vzLo, vzHi), 200, -0.1, 0.1);

        h_eta_diff_corrected_E_vz[i][iVz] = new TH1F(
              Form("h_eta_diff_corr_%s_%s", lab.c_str(), vzTag.c_str()),
              Form("#Delta#eta corr;%.0f<E<%.0f GeV, %.0f<|z_{vtx}|<%.0f cm",
                   eLo, eHi, vzLo, vzHi), 200, -0.1, 0.1);

        /* CP variants */
        h_eta_diff_cpRaw_E_vz [i][iVz] = new TH1F(
              Form("h_eta_diff_cpRaw_%s_%s",  lab.c_str(), vzTag.c_str()),
              "#Delta#eta RAW‑CP;#Delta#eta;Counts", 200, -0.1, 0.1);
        h_eta_diff_cpCorr_E_vz[i][iVz] = new TH1F(
              Form("h_eta_diff_cpCorr_%s_%s", lab.c_str(), vzTag.c_str()),
              "#Delta#eta CP‑corr;#Delta#eta;Counts", 200, -0.1, 0.1);
        h_eta_diff_cpBcorr_E_vz[i][iVz] = new TH1F(
              Form("h_eta_diff_cpBcorr_%s_%s", lab.c_str(), vzTag.c_str()),
              "#Delta#eta BCORR;#Delta#eta;Counts", 200, -0.1, 0.1);

        hm->registerHisto(h_eta_diff_raw_E_vz     [i][iVz]);
        hm->registerHisto(h_eta_diff_corrected_E_vz[i][iVz]);
        hm->registerHisto(h_eta_diff_cpRaw_E_vz   [i][iVz]);
        hm->registerHisto(h_eta_diff_cpCorr_E_vz  [i][iVz]);
        hm->registerHisto(h_eta_diff_cpBcorr_E_vz [i][iVz]);
      }
    }


  /* Δy scan & log‑weight scan histograms – exactly as before */
  for (int iE=0;iE<N_Ebins;++iE)
  {
    const std::string lab = makeLabel(iE);
    for (double b : m_bScan)
      hm->registerHisto(new TH1F(Form("h_dy_ash_b%.4f_%s",b,lab.c_str()),
                                 ";y_{reco}-y_{true}  [cm];Counts / 0.1 cm",
                                 240,-12,12));
    for (double w0 : m_w0Scan)
      hm->registerHisto(new TH1F(Form("h_dy_log_w0%.2f_%s",w0,lab.c_str()),
                                 ";y_{reco}-y_{true}  [cm];Counts / 0.1 cm",
                                 240,-12,12));
  }

  /* Δx scan (identical) */
  constexpr double DX_MAX = 12.0, BIN_W = 0.10;
  constexpr int    NBINS  = int(2*DX_MAX/BIN_W + 0.5);

  if (!alreadyDeclaredHistograms)
  {
    for (double b = 0.01; b <= 0.50 + 1e-9; b += 0.01)
      m_bScan.push_back(std::round(b*10000.)/10000.);
    for (double w0 = 1.5;  w0 <= 7.0  + 1e-9; w0+=0.10)
      m_w0Scan.push_back(std::round(w0*100.)/100.);

    for (int iE=0;iE<N_Ebins;++iE)
    {
      const std::string lab = makeLabel(iE);

      for (double bVal : m_bScan)
      {
        const TString hName = Form("h_dx_ash_b%.4f_%s",bVal,lab.c_str());
        hm->registerHisto(new TH1F(hName,
                                   ";x_{reco}-x_{true}  [cm];Counts / 0.1 cm",
                                   NBINS,-DX_MAX,DX_MAX));
      }
      for (double w0 : m_w0Scan)
      {
        const TString hName = Form("h_dx_log_w0%.2f_%s",w0,lab.c_str());
        hm->registerHisto(new TH1F(hName,
                                   ";x_{reco}-x_{true}  [cm];Counts / 0.1 cm",
                                   NBINS,-DX_MAX,DX_MAX));
      }
    }
    alreadyDeclaredHistograms = true;
  }
}


// ------------------------------------------------------------------
//  Reads "bValues.txt" and populates m_bVals{Phi,Eta}.
//
//  Return ‑ true  : tables loaded, flags (isFitDoneForPhi/Eta) updated
//          ‑ false : file missing → both corrections disabled
// ------------------------------------------------------------------
bool PositionDependentCorrection::loadBValues(const std::string& bFilePath)
{
  if (!isFitDoneForPhi && !isFitDoneForEta)
  {
    std::cout << "[INFO]   isFitDoneForPhi && isFitDoneForEta are both false ⇒ "
                 "skip reading bValues.txt\n";
    return true;                                   // nothing to do, not an error
  }

  std::ifstream bfile(bFilePath);
  if (!bfile.is_open())
  {
    std::cout << "[WARN]  bValues.txt NOT found at " << bFilePath
              << "  ➜  disabling φ‑ and η‑corrections.\n";
    isFitDoneForPhi = isFitDoneForEta = false;
    return false;
  }

  std::cout << "[INFO]  Reading b‑values from " << bFilePath << '\n';

  std::vector<bool> gotBinPhi(N_Ebins,false), gotBinEta(N_Ebins,false);
  const std::regex lineRe(
      R"(^\s*(PHI|ETA)\s*\[\s*([0-9]*\.?[0-9]+)\s*,\s*([0-9]*\.?[0-9]+)\s*\)\s*([0-9]*\.?[0-9]+)\s*$)"
  );

  std::string line;
  while (std::getline(bfile,line))
  {
    if (line.empty() || line[0]=='#') continue;
    std::smatch m;
    if (!std::regex_match(line,m,lineRe)) continue;

    const std::string dim = m[1];
    const float eLo   = std::stof(m[2]);
    const float eHi   = std::stof(m[3]);
    const float bVal  = std::stof(m[4]);

    for (int i=0;i<N_Ebins;++i)
      if (std::fabs(eLo-expectedLo[i])<1e-3 && std::fabs(eHi-expectedHi[i])<1e-3)
      {
        if (dim=="PHI") { m_bValsPhi[i]=bVal; gotBinPhi[i]=true; }
        else            { m_bValsEta[i]=bVal; gotBinEta[i]=true; }
        break;
      }
  }
  bfile.close();

  const bool allPhi = std::all_of(gotBinPhi.begin(),gotBinPhi.end(),[](bool v){return v;});
  const bool allEta = std::all_of(gotBinEta.begin(),gotBinEta.end(),[](bool v){return v;});
  isFitDoneForPhi &= allPhi;
  isFitDoneForEta &= allEta;

  auto printTable=[&](const char* tag,const std::vector<bool>& got,
                      const float* arr,bool enabled)
  {
    std::cout << "[INFO]  " << tag << (enabled
                ?"  (enabled)\n":"  (DISABLED — missing slices)\n");
    for(int i=0;i<N_Ebins;++i)
    {
      std::cout << "        ["<<expectedLo[i]<<','<<expectedHi[i]<<")  : ";
      if (got[i]) std::cout << arr[i] << '\n';
      else        std::cout << "-- missing --\n";
    }
  };
  printTable("PHI b‑values",gotBinPhi,m_bValsPhi,isFitDoneForPhi);
  printTable("ETA b‑values",gotBinEta,m_bValsEta,isFitDoneForEta);

  return true;
}


// ------------------------------------------------------------------
//  Dummy parser for the pass‑1 mass‑fit table.
//  Expected columns:  TAG  μ  σ  E_LOW  E_HIGH
// ------------------------------------------------------------------
bool PositionDependentCorrection::loadMassWindowTable(const std::string& path)
{
  std::ifstream fin(path);
  if (!fin.is_open())
  {
    std::cerr << "[MassWin]  table \"" << path << "\" not found – skip\n";
    return false;
  }

  std::cout << "[MassWin]  loading pass‑1 fit results from " << path << '\n';

  std::string tag;
  float mu, sigma, elo, ehi;
  while (fin >> tag >> mu >> sigma >> elo >> ehi)
  {
    int idx = -1;
    for (int i = 0; i < N_Ebins; ++i)
      if (std::fabs(elo-eEdge[i])<1e-3f && std::fabs(ehi-eEdge[i+1])<1e-3f)
      { idx = i; break; }

    if (idx < 0) continue;                           // unknown slice → ignore

    if      (strcasecmp(tag.c_str(),"CORR")   == 0) m_winCorr[idx] = {mu,sigma};
    else if (strcasecmp(tag.c_str(),"UNCORR") == 0) m_winRaw [idx] = {mu,sigma};
  }
  fin.close();
  std::cout << "[MassWin]  table parsed – dummy values stored.\n";
  return true;
}


// -------------------------------------------------------------------------
//  New, compact Init() – functional behaviour identical to the original
// -------------------------------------------------------------------------
int PositionDependentCorrection::Init(PHCompositeNode* /*topNode*/)
{
  /* 0) Basic sanity / I/O managers ----------------------------------- */
  if (Verbosity() > 0)
    std::cout << "[DEBUG] PositionDependentCorrection::Init() called.\n";

  hm = new Fun4AllHistoManager(Name());
  if (!hm)
  {
    std::cerr << "[ERROR] Fun4AllHistoManager allocation failed.\n";
    return -1;
  }

  outfile = new TFile(outfilename.c_str(),"RECREATE");
  std::cout << "[INFO-OUTFILENAME] Writing histograms to: "
            << outfilename << std::endl;
  if (!outfile || outfile->IsZombie())
  {
    std::cerr << "[ERROR] Could not open output file '" << outfilename << "'\n";
    return -1;
  }

  trigAna = new TriggerAnalyzer();

  /* 1) Survey geometry ----------------------------------------------- */
  if (m_useSurveyGeometry)
  {
      if (Verbosity() > 0)
        std::cout << ANSI_BOLD << ANSI_CYAN
                  << "[PDC]  Loading CDB CALO_TOWER_GEOMETRY  "
                  << "(tag = " << m_cdbTag << ",  ts = " << m_timeStamp << ")"
                  << ANSI_RESET << '\n';

      /* make CDB happy */
      recoConsts* rc = recoConsts::instance();
      rc->set_StringFlag("CDB_GLOBALTAG", m_cdbTag);
      rc->set_uint64Flag("TIMESTAMP",     m_timeStamp);


      const std::string url =
          CDBInterface::instance()->getUrl("CALO_TOWER_GEOMETRY");

      std::unique_ptr<CDBTTree> geomTree(new CDBTTree(url));
      geomTree->LoadCalibrations();                    // <— fills the tree

      if (Verbosity() > 0)
        std::cout << ANSI_GREEN
                  << "      … geometry loaded from\n      " << url
                  << ANSI_RESET << "\n";

      // ────────────────────────────────────────────────────────────────
      // Extract the rigid φ-offset calculation from Geometry
      // ────────────────────────────────────────────────────────────────
      if (Verbosity() > 0)
      {
        constexpr int kNPhi = 256;                     // full fine bins
        std::cout << ANSI_BOLD << "[φ-survey summary]" << ANSI_RESET << '\n';

        const double phi_first  = geomTree->GetDoubleValue(0, "cemc_phi_first");
        const double phiOffset  = M_PI/2.0 - phi_first;    // your convention
        std::cout << "      global Δφ₀ = "
                  << ANSI_YELLOW << std::fixed << std::setprecision(6)
                  << phiOffset
                  << "  rad"
                  << ANSI_RESET << "  =  "
                  << ANSI_YELLOW << std::setprecision(4)
                  << phiOffset * 180.0 / M_PI
                  << "°" << ANSI_RESET << '\n';

        std::cout << ANSI_CYAN
                  << "      bin   φ_first [rad]   φ_second [rad]\n"
                  << "      ─────────────────────────────────────\n"
                  << ANSI_RESET;

      /* print either all bins (verbosity>1) or the first few */
      const int nPrint = (Verbosity() > 1) ? kNPhi : 6;

      for (int i = 0; i < nPrint; ++i)
      {
          const double ph1 = geomTree->GetDoubleValue(i, "cemc_phi_first");
          const double ph2 = geomTree->GetDoubleValue(i, "cemc_phi_second");
          std::cout << std::setw(10) << i
                    << std::setw(15) << std::setprecision(6) << ph1
                    << std::setw(15) << ph2 << '\n';
      }
      if (Verbosity() == 1 && kNPhi > nPrint)
          std::cout << "      … (" << kNPhi - nPrint
                    << " more bins – increase Verbosity() for full table)\n";
      }
  }
  else
  {
      if (Verbosity() > 0)
        std::cout << ANSI_BOLD << ANSI_YELLOW
                  << "[PDC]  Survey geometry NOT loaded "
                     "(m_useSurveyGeometry = false)."
                  << ANSI_RESET << '\n';
  }

  /* 2) Histogram booking --------------------------------------------- */
  auto makeLabel = [&](int i)->std::string
  { return (m_binningMode==EBinningMode::kRange)
           ? Form("%.0f_%.0f",eEdge[i],eEdge[i+1])
           : Form("E%.0f",eEdge[i]); };

  bookCommonHistograms(makeLabel);                 // always‑on

  if (m_isSimulation)                              // truth extras
    bookSimulationHistograms(makeLabel);
  else if (Verbosity()>0)
    std::cout << "[PDC]  Data taking: truth histograms are NOT booked.\n";

  /* 3) Optional tables ------------------------------------------------ */

  /* 3a)  b‑parameters */
  const std::string bFilePath =
      "/sphenix/u/patsfan753/scratch/PDCrun24pp/bParameters/bValues.txt";
  loadBValues(bFilePath);

  /* 3b)  π0‑mass window (only in data mode & after pass‑1 fits) */
  if (!m_isSimulation && m_massFitsDone)
  {
    const std::string mwPath =
        "/path/to/massWindowTable.txt";            // adjust if needed
    if (!loadMassWindowTable(mwPath))
    {
      std::cout << "[WARN] Mass‑window file missing – block‑maps will be empty.\n";
      m_massFitsDone = false;                      // graceful fallback
    }
  }

  /* 4) Misc ----------------------------------------------------------- */
  rnd = new TRandom3();

  if (Verbosity()>0)
    std::cout << "[DEBUG] PositionDependentCorrection::Init() completed successfully.\n";

  return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int PositionDependentCorrection::process_event(PHCompositeNode* topNode)
{
    std::cout << ANSI_BOLD
    << "\n---------  processing event " << _eventcounter << "  ---------\n"
    << ANSI_RESET << std::endl;
    
    if (Verbosity() > 0)
    {
        std::cout << "[DEBUG] PositionDependentCorrection::process_event() called. "
        << "Current event counter = " << _eventcounter << std::endl;
    }
    
    ++_eventcounter;
    process_towers(topNode);
    return Fun4AllReturnCodes::EVENT_OK;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float
PositionDependentCorrection::retrieveVertexZ(PHCompositeNode* topNode)
{

  float vtx_z = 0.f;
  /*-----------------------------------------------------------------*/
  /* 2) locate the GlobalVertexMap node                              */
  /*-----------------------------------------------------------------*/
  if (Verbosity() > 0)
    std::cout << "  Searching for node \"GlobalVertexMap\" … ";

  auto* vtxMap =
      findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

  if (!vtxMap)
  {
    if (Verbosity() > 0)
    {
      std::cout << ANSI_RED << "NOT FOUND" << ANSI_RESET << '\n'
                << ANSI_YELLOW
                << "  → cannot fill h_reco_vz.\n"
                << ANSI_RESET;
    }
    return vtx_z;
  }

  if (Verbosity() > 0)
    std::cout << ANSI_GREEN << "OK" << ANSI_RESET << '\n';

  /*-----------------------------------------------------------------*/
  /* 3) map must not be empty                                        */
  /*-----------------------------------------------------------------*/
  if (vtxMap->empty())
  {
    if (Verbosity() > 0)
    {
      std::cout << ANSI_RED
                << "  GlobalVertexMap is EMPTY.\n"
                << ANSI_RESET
                << ANSI_YELLOW
                << "  → h_reco_vz NOT filled.\n"
                << ANSI_RESET;
    }
    return vtx_z;
  }

  /*-----------------------------------------------------------------*/
  /* 4) retrieve the first vertex                                    */
  /*-----------------------------------------------------------------*/
  GlobalVertex* vtx = vtxMap->begin()->second;
  if (!vtx)
  {
    if (Verbosity() > 0)
      std::cout << ANSI_RED
                << "  First vertex pointer is NULL.\n"
                << ANSI_RESET;
    return vtx_z;
  }

  vtx_z = vtx->get_z();
  if (Verbosity() > 0)
    std::cout << "  Raw vtx_z read .......... " << vtx_z << " cm\n";

  /*-----------------------------------------------------------------*/
  /* 5) sanity: finite number?                                       */
  /*-----------------------------------------------------------------*/
  if (!std::isfinite(vtx_z))
  {
    if (Verbosity() > 0)
      std::cout << ANSI_RED
                << "  vtx_z is NOT finite!  Using 0.\n"
                << ANSI_RESET;
    vtx_z = 0.f;
    return vtx_z;
  }

  /*-----------------------------------------------------------------*/
  /* 6) fill histogram if it exists                                  */
  /*-----------------------------------------------------------------*/
  if (h_reco_vz)
  {
    h_reco_vz->Fill(vtx_z);
    if (Verbosity() > 0)
      std::cout << ANSI_GREEN
                << "  h_reco_vz filled.\n"
                << ANSI_RESET;
  }
  else if (Verbosity() > 0)
  {
    std::cout << ANSI_RED
              << "  h_reco_vz pointer is NULL → NOT filled!\n"
              << ANSI_RESET;
  }

  /*-----------------------------------------------------------------*/
  /* 7) normal exit                                                  */
  /*-----------------------------------------------------------------*/
  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD << ANSI_CYAN
              << "[retrieveVertexZ] EXIT  →  vtx_z = "
              << vtx_z << " cm"
              << ANSI_RESET << "\n";
  }
  return vtx_z;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PositionDependentCorrection::fillTowerInfo(
    PHCompositeNode* topNode,
    float emcal_hit_threshold,
    float& tower_tot_e,
    std::vector<float>& ht_eta,
    std::vector<float>& ht_phi)
{
  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  if (towers)
  {
    int size = towers->size();

    for (int channel = 0; channel < size; channel++)
    {
      TowerInfo* tower = towers->get_tower_at_channel(channel);

      float offlineenergy = tower->get_energy();

      unsigned int towerkey = towers->encode_key(channel);
      int ieta = towers->getTowerEtaBin(towerkey);
      int iphi = towers->getTowerPhiBin(towerkey);

      bool isGood = !(tower->get_isBadChi2());

      tower_tot_e += offlineenergy;

      h_tower_e->Fill(offlineenergy);

      if (!isGood && offlineenergy > 0.2)
      {
        ht_eta.push_back(ieta);
        ht_phi.push_back(iphi);
      }
      if (isGood)
      {
        h_emcal_e_eta->Fill(ieta, offlineenergy);
      }

      if (offlineenergy > emcal_hit_threshold)
      {
        h_cemc_etaphi->Fill(ieta, iphi);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PositionDependentCorrection::fillTruthInfo(
    PHCompositeNode* topNode,
    float& vtx_z,
    std::vector<TLorentzVector>& truth_photons,
    std::vector<TLorentzVector>& truth_meson_photons)
{
  float wieght = 1;

  if (Verbosity() > 0)
  {
    std::cout << std::endl
              << ANSI_BOLD << ANSI_CYAN
              << "[fillTruthInfo] START" << ANSI_RESET << "\n"
              << "  Attempting to retrieve PHG4TruthInfoContainer named 'G4TruthInfo'..."
              << std::endl;
  }

  // -------------------------------------------------------------------------
  // 0) Retrieve PHG4TruthInfoContainer
  // -------------------------------------------------------------------------
  PHG4TruthInfoContainer* truthinfo =
      findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!truthinfo)
  {
    if (Verbosity() > 0)
    {
      std::cout << ANSI_BOLD << ANSI_YELLOW
                << "  [WARNING] 'G4TruthInfo' not found in topNode. "
                << "No truth info to fill => returning early..." << ANSI_RESET
                << std::endl
                << ANSI_BOLD << "[fillTruthInfo] END\n" << ANSI_RESET
                << std::endl;
    }
    return;
  }

  if (Verbosity() > 0)
  {
    std::cout << "  --> Successfully retrieved 'G4TruthInfo' container.\n"
              << "  --> Now analyzing primary particles..." << std::endl;
  }
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second; ++iter)
  {
    const PHG4Particle* truth = iter->second;
    if (!truthinfo->is_primary(truth)) continue;

    TLorentzVector myVector;
    myVector.SetXYZT(truth->get_px(), truth->get_py(),
                     truth->get_pz(), truth->get_e());

    float energy = myVector.E();
    h_truth_eta->Fill(myVector.Eta());
    h_truth_e->Fill(energy, wieght);

    if (debug)
    {
      std::cout << ANSI_BOLD << "[Debug]" << ANSI_RESET
                << " Primary pid=" << truth->get_pid()
                << "   E=" << energy
                << "  pt=" << myVector.Pt()
                << "  eta=" << myVector.Eta()
                << "  phi=" << myVector.Phi() << std::endl;
    }
    truth_photons.push_back(myVector);
  }
  if (Verbosity() > 0)
  {
    std::cout << "  --> Finished analyzing primary particles."
              << " Now looking at secondary particles..." << std::endl;
  }

  // -------------------------------------------------------------------------
  // 2) Secondary particles
  // -------------------------------------------------------------------------
  PHG4TruthInfoContainer::Range second_range = truthinfo->GetSecondaryParticleRange();
  float m_g4 = 0;
  for (PHG4TruthInfoContainer::ConstIterator siter = second_range.first;
       siter != second_range.second; ++siter)
  {
    if (m_g4 >= 19999) break; // original guard

    const PHG4Particle* truth = siter->second;
    if (truth->get_pid() == 22) // photon
    {
      PHG4Particle* parent = truthinfo->GetParticle(truth->get_parent_id());
      if (!parent)
      {
        if (Verbosity() > 0)
        {
          std::cout << ANSI_BOLD << ANSI_YELLOW
                    << "  [WARNING] Secondary photon has no valid parent => skipping pi0/eta check..."
                    << ANSI_RESET << std::endl;
        }
      }
      else
      {
        // If the parent is pi0 (111) or eta (221), treat photon as "meson decay photon"
        if (parent->get_pid() == 111 || parent->get_pid() == 221)
        {
          float phot_pt = std::sqrt(truth->get_px() * truth->get_px() +
                                    truth->get_py() * truth->get_py());
          if (phot_pt < 0.1) continue;

          float phot_e   = truth->get_e();
          float phot_phi = std::atan2(truth->get_py(), truth->get_px());
          float phot_eta = std::atanh(
              truth->get_pz() /
              std::sqrt(truth->get_px()*truth->get_px() +
                        truth->get_py()*truth->get_py() +
                        truth->get_pz()*truth->get_pz()));

          // build TLorentzVector
          TLorentzVector myVector;
          myVector.SetXYZT(truth->get_px(), truth->get_py(),
                           truth->get_pz(), truth->get_e());

          truth_meson_photons.push_back(myVector);

          if (debug)
          {
            std::cout << ANSI_BOLD << "[Debug]" << ANSI_RESET
                      << " 2nd photon  pt=" << phot_pt
                      << "  e=" << phot_e
                      << "  phi=" << phot_phi
                      << "  eta=" << phot_eta << std::endl;
          }
        }
      }
    }
  }
  if (Verbosity() > 0)
  {
    std::cout << "  --> Finished analyzing secondary particles."
              << " Now looking at vertices..." << std::endl;
  }

  // -------------------------------------------------------------------------
  // 3) Vertex loop
  // -------------------------------------------------------------------------
  PHG4TruthInfoContainer::VtxRange vtxrange = truthinfo->GetVtxRange();
  int n_vertex = 0;
  for (PHG4TruthInfoContainer::ConstVtxIterator iter = vtxrange.first;
       iter != vtxrange.second; ++iter)
  {
    PHG4VtxPoint* vtx = iter->second;
    h_vert_xy->Fill(vtx->get_x(), vtx->get_y());

    // If vtx id=1 => record primary vertex z
    if (vtx->get_id() == 1)
    {
      vtx_z = vtx->get_z();
      h_truth_vz->Fill(vtx_z);
    }

    // Hidden debug block if 'false' changed to 'true'
    if (vtx->get_id() == 1 && false)
    {
      std::cout << "vx=" << vtx->get_x()
                << "  vy=" << vtx->get_y()
                << "   vz=" << vtx->get_z()
                << "  id=" << vtx->get_id() << std::endl;
    }

    n_vertex++;
    if (n_vertex >= 100000) break; // original guard
  }
  if (Verbosity() > 0)
  {
    std::cout << "  --> Finished analyzing vertices.\n"
              << "  --> Summary:\n"
              << "       # of primary photons (truth_photons):      "
              << truth_photons.size() << "\n"
              << "       # of meson decay photons (truth_meson_photons): "
              << truth_meson_photons.size() << "\n"
              << "       Primary vertex z-position:                 "
              << vtx_z << "\n"
              << ANSI_BOLD << "[fillTruthInfo] END" << ANSI_RESET
              << std::endl << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
RawTowerGeomContainer* PositionDependentCorrection::checkTowerGeometry(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "PositionDependentCorrection::checkTowerGeometry - START" << std::endl;
    std::cout << "  --> Attempting to retrieve 'TOWERGEOM_CEMC' from the node tree." << std::endl;
  }

  std::string towergeomnodename = "TOWERGEOM_CEMC";
  RawTowerGeomContainer* geo = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!geo)
  {
    // Existing error printout, preserved:
    if (Verbosity() > 0)
    {
      std::cout << "  [ERROR] " << Name() << "::CreateNodeTree"
                << ": Could not find node " << towergeomnodename << std::endl;
      std::cout << "  --> Throwing std::runtime_error due to missing TOWERGEOM node." << std::endl;
    }
    throw std::runtime_error("failed to find TOWERGEOM node");
  }

  if (Verbosity() > 0)
  {
    std::cout << "  --> Successfully retrieved 'TOWERGEOM_CEMC' node." << std::endl;
    std::cout << "PositionDependentCorrection::checkTowerGeometry - END" << std::endl << std::endl;
  }
  return geo;
}

// ----------------------------------------------------------------------------

RawClusterContainer* PositionDependentCorrection::retrieveClusterContainer(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "PositionDependentCorrection::retrieveClusterContainer - START" << std::endl;
    std::cout << "  --> Attempting to retrieve 'CLUSTERINFO_CEMC' from the node tree." << std::endl;
  }

  // Attempt to get the cluster container
  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
  if (!clusterContainer)
  {
    // Existing printout, preserved:
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE
                << "funkyCaloStuff::process_event - Fatal Error - CLUSTER_CEMC node is missing."
                << std::endl;
      std::cout << "PositionDependentCorrection::retrieveClusterContainer - END (nullptr returned)"
                << std::endl << std::endl;
    }
    return nullptr;
  }

  if (Verbosity() > 0)
  {
    std::cout << "  --> Successfully retrieved 'CLUSTERINFO_CEMC' container." << std::endl;
    std::cout << "PositionDependentCorrection::retrieveClusterContainer - END" << std::endl << std::endl;
  }
  return clusterContainer;
}

// ----------------------------------------------------------------------------

int PositionDependentCorrection::countClusters(
    RawClusterContainer* clusterContainer,
    float vtx_z,
    float nClus_ptCut,
    float clus_chisq_cut,
    int& nClusContainer)
{
  if (Verbosity() > 0)
  {
    std::cout << "PositionDependentCorrection::countClusters - START" << std::endl;
    std::cout << "  --> vtx_z = " << vtx_z
              << ", ptCut = " << nClus_ptCut
              << ", chisq_cut = " << clus_chisq_cut << std::endl;
    std::cout << "  --> About to loop over clusters in container and apply basic cuts." << std::endl;
  }

  // We do a quick pass to count how many clusters pass basic cuts
  RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
  RawClusterContainer::ConstIterator clusterIter;

  int nClusCount = 0;
  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; ++clusterIter)
  {
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(0, 0, vtx_z);
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetEVec(*recoCluster, vertex);

    float clus_pt    = E_vec_cluster.perp();
    float clus_chisq = recoCluster->get_chi2();

    // Existing debug print
    if (debug && clus_pt > 0.1)
    {
      std::cout << "  [Debug] clus #" << nClusCount
                << "  E=" << E_vec_cluster.mag()
                << "  eta=" << E_vec_cluster.pseudoRapidity()
                << "  chi2=" << clus_chisq << std::endl;
    }

    nClusContainer++;

    if (clus_pt < nClus_ptCut)       continue;
    if (clus_chisq > clus_chisq_cut) continue;

    // If it passes the cut, increment
    nClusCount++;
  }

  if (Verbosity() > 0)
  {
    std::cout << "  --> Total clusters passing cuts: " << nClusCount << std::endl;
    std::cout << "  --> nClusContainer incremented to: " << nClusContainer << std::endl;
    std::cout << "PositionDependentCorrection::countClusters - END" << std::endl << std::endl;
  }
  return nClusCount;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// doLogWeightCoord(...)
// --------------------------------------------------------------------
// • Computes the local block-φ  ∈ [0,1)  using the “log-weight” centroid
//   w_i = max( 0 , w0 + ln(E_i / ΣE) ).
// • Wrap handling and unit convention are now **identical** to the Ash
//   branch, so the two algorithms can be compared one-to-one.
//
//   INPUT:
//     towerPhi   – fine φ-index of every tower hit            (0 … 255)
//     towerE     – corresponding tower energies               (GeV)
//     w0         – logarithmic weight parameter               (unit-less)
//
//   RETURNS:
//     local φ inside its 2×2 block (0 ≤ φloc < 1)             (tower units)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float
PositionDependentCorrection::doLogWeightCoord( const std::vector<int>&   towerPhi,
                                               const std::vector<float>& towerE,
                                               float                    w0 )
{
  if (towerPhi.empty() || towerE.empty())
  {
    if (Verbosity() > 0)
      std::cout << "[doLogWeightCoord]  empty input → return 0\n";
    return 0.f;
  }
  if (towerPhi.size() != towerE.size())
  {
    if (Verbosity() > 0)
      std::cout << ANSI_RED << "[doLogWeightCoord]  size mismatch: "
                << towerPhi.size() << " φ bins vs "
                << towerE.size()   << " energies  → return 0"
                << ANSI_RESET << '\n';
    return 0.f;
  }

  /*──────────────────────────────
   * 1) Total cluster energy
   *──────────────────────────────*/
  const double sumE = std::accumulate( towerE.begin(), towerE.end(), 0.0 );
  if (sumE < 1e-9)
  {
    if (Verbosity() > 0)
      std::cout << "[doLogWeightCoord]  ΣE≈0 → return 0\n";
    return 0.f;
  }

  /*──────────────────────────────
   * 2) Reference tower  (highest-E)
   *──────────────────────────────*/
  const std::size_t iRef =
      std::distance( towerE.begin(),
                     std::max_element( towerE.begin(), towerE.end() ) );
  const int refPhi = towerPhi[iRef];                   // 0 … 255
  constexpr int nPhiBin = 256;

  if (Verbosity() > 0)
  {
    std::cout << "[doLogWeightCoord]  nTower=" << towerE.size()
              << "  ΣE=" << sumE
              << "  refPhi=" << refPhi << '\n';
  }

  /*──────────────────────────────
   * 3) Log-weight centroid
   *──────────────────────────────*/
  double sumW = 0.0, sumWPhi = 0.0;
  for (std::size_t i = 0; i < towerE.size(); ++i)
  {
    const double frac = towerE[i] / sumE;
    double w = w0 + std::log(frac);
    if (w < 0.0) w = 0.0;                               // clamp → 0

    int phi = towerPhi[i];
    /* wrap fix relative to refPhi */
    if (phi - refPhi < -nPhiBin/2) phi += nPhiBin;
    if (phi - refPhi >  nPhiBin/2) phi -= nPhiBin;

    sumW    += w;
    sumWPhi += w * phi;

    if (Verbosity() > 1)
      std::cout << "    tower#" << i << "  φ=" << towerPhi[i]
                << "  φAdj=" << phi << "  E=" << towerE[i]
                << "  w=" << w << '\n';
  }

  if (sumW < 1e-9)
  {
    if (Verbosity() > 0)
      std::cout << "[doLogWeightCoord]  Σw≈0 → return 0\n";
    return 0.f;
  }

  double phiAvg = sumWPhi / sumW;                       // could be negative
  phiAvg = std::fmod( phiAvg + nPhiBin, nPhiBin );      // wrap → 0 … 256

  /*──────────────────────────────
   * 4) Convert to local block-φ
   *──────────────────────────────*/
  const double localPhiRaw = std::fmod( phiAvg + 0.5, 2.0 ); // 0 … 2
  const float  localPhi    =
      static_cast<float>( localPhiRaw < 1.0 ? localPhiRaw
                                            : localPhiRaw - 1.0 ); // 0 … 1

  if (Verbosity() > 0)
  {
    std::cout << "[doLogWeightCoord]  φAvg="   << phiAvg
              << "  → φloc=" << localPhi << '\n';
  }
  return localPhi;
}


/************************************************************************
 *  fillAshLogDx(...) — centimetre‑accurate block‑φ diagnostics
 *
 *  Purpose
 *  ──────────────────────────────────────────────────────────────────
 *  For one cluster and one truth photon
 *   → propagate several candidate local‑φ values (Ash‑b scan + log‑w₀ scan)
 *   → to shower depth (x‑coordinate in cm)
 *   → fill Δx = xReco – xTrue histograms.
 *
 *  Inputs
 *  ──────
 *    • cluster            — pointer to RawCluster under study
 *    • recoPhoton         — 4‑vector obtained from cluster (energy slice tag)
 *    • truthPhoton        — matched truth photon 4‑vector
 *    • blockCord          — (ηloc, φloc) of the cluster’s 2×2 super‑cell
 *    • blockPhiBin        — coarse φ‑block index  (0 … 127)
 *    • tower_phis         — fine φ indices of all towers in the cluster
 *    • tower_energies     — corresponding tower energies
 *
 *  Histogram naming scheme
 *  ───────────────────────
 *    h_dx_ash_b%.4f_<sliceTag>   – Ash‑shift with parameter *b*
 *    h_dx_log_w0%.2f_<sliceTag>  – log‑weight centroid with parameter *w₀*
 *
 *  Notes on implementation
 *  ───────────────────────
 *  • One *outer* constant   xTrueRef   is computed once (front‑face x of truth).
 *  • The inner λ‑function **tryFill** converts a candidate local‑φ
 *    to shower‑depth xReco and fills the requested histogram.
 *  • ‑Wshadow was raised because the original code **redeclared** xTrue
 *    inside the λ.  This fix removes that local redeclaration and reuses
 *    xTrueRef captured by reference.
 ************************************************************************/
void PositionDependentCorrection::fillAshLogDx(
        RawCluster*                     cluster,
        const TLorentzVector&           recoPhoton,
        const TLorentzVector&           truthPhoton,
        const std::pair<float,float>&   blockCord,     // (ηloc , φloc)
        int                             blockPhiBin,   // 0 … 127
        const std::vector<int>&         tower_phis,
        const std::vector<float>&       tower_energies )
{
    // ────────────────────────────────────────────────────────────────
    // 0)  Energy‑slice bookkeeping
    // ────────────────────────────────────────────────────────────────
    const float eReco  = recoPhoton.E();
    const int   iSlice = getEnergySlice(eReco);
    if (iSlice < 0) return;                  // outside configured range
    const std::string lab = sliceTag(iSlice);

    // ────────────────────────────────────────────────────────────────
    // 1)  Sanity guards & lead‑tower geometry
    // ────────────────────────────────────────────────────────────────
    if (!cluster || !m_geometry || !m_bemcRec) return;

    float eMax = -1.F;  RawCluster::TowerConstIterator best{};
    for (auto it = cluster->get_towers().first;
              it != cluster->get_towers().second; ++it)
        if (it->second > eMax) { eMax = it->second; best = it; }
    if (eMax < 1e-6F) return;               // pathological cluster

    auto* leadGeo = m_geometry->get_tower_geometry(best->first);
    if (!leadGeo) return;

    const double rFront = leadGeo->get_center_radius();  // front‑face radius
    const double zFront = leadGeo->get_center_z();       // front‑face z
    const int    ixLead = leadGeo->get_binphi();         // fine φ index (0…255)
    const int    iyLead = leadGeo->get_bineta();         // fine η index (0…95)

    // ────────────────────────────────────────────────────────────────
    // 2)  Truth reference at shower depth (here we keep front‑face x)
    // ────────────────────────────────────────────────────────────────
    const float  phiSDtruth = TVector2::Phi_mpi_pi(truthPhoton.Phi());
    const double xTrueRef   = rFront * std::cos(phiSDtruth);

    if (Verbosity() > 1)
        std::cout << ANSI_CYAN
                  << "[fillAshLogDx] slice=" << iSlice
                  << "  E="        << eReco
                  << " GeV  φSDtruth=" << phiSDtruth
                  << "  → xTrue="  << xTrueRef << " cm"
                  << ANSI_RESET   << '\n';

    // ────────────────────────────────────────────────────────────────
    // 3)  Helper λ  –  propagate φloc → xReco and fill histogram
    // ────────────────────────────────────────────────────────────────
    auto tryFill = [&](float        localPhi,
                       const TString& hName,
                       int&         okCtr,
                       int&         missCtr)
    {
        if (localPhi < 0.F || localPhi >= 1.F) return;   // outside canonical cell

        /* front‑face φ (tower grid) */
        const float phiFront = convertBlockToGlobalPhi(blockPhiBin, localPhi);

        /* front‑face → shower‑depth */
        float xSD, ySD, zSD;
        m_bemcRec->CorrectShowerDepth( ixLead, iyLead,
                                       eReco,
                                       rFront * std::cos(phiFront),
                                       rFront * std::sin(phiFront),
                                       zFront,
                                       xSD, ySD, zSD );

        const double xReco = static_cast<double>(xSD);

        if (auto* h = dynamic_cast<TH1F*>(hm->getHisto(hName.Data())))
        {
            h->Fill(xReco - xTrueRef);      // Δx in centimetres
            ++okCtr;
        }
        else
            ++missCtr;
    };

    // ────────────────────────────────────────────────────────────────
    // 4)  Ash‑b scan (front‑loaded block‑φ variant)
    // ────────────────────────────────────────────────────────────────
    int nAshOK = 0, nAshMiss = 0;
    if (blockCord.second >= 0.F && blockCord.second < 1.F)
        for (double b : m_bScan)
        {
            const double bVal   = std::round(b * 1e4) / 1e4;   // nicer label
            const float  phiLoc = doAshShift(blockCord.second, bVal);
            tryFill(phiLoc,
                    Form("h_dx_ash_b%.4f_%s", bVal, lab.c_str()),
                    nAshOK, nAshMiss);
        }

    // ────────────────────────────────────────────────────────────────
    // 5)  Log‑weight w₀ scan
    // ────────────────────────────────────────────────────────────────
    int nLogOK = 0, nLogMiss = 0;
    if (!tower_phis.empty() && !tower_energies.empty())
        for (double w0 : m_w0Scan)
        {
            const double w0Val  = std::round(w0 * 100.) / 100.;
            const float  phiLoc = doLogWeightCoord(tower_phis, tower_energies,
                                                   w0Val);
            tryFill(phiLoc,
                    Form("h_dx_log_w0%.2f_%s", w0Val, lab.c_str()),
                    nLogOK, nLogMiss);
        }

    // ────────────────────────────────────────────────────────────────
    // 6)  Verbose summary (optional)
    // ────────────────────────────────────────────────────────────────
    if (Verbosity() > 0)
        std::cout << "  ➜ Ash filled/miss = "
                  << nAshOK << '/' << nAshMiss
                  << "   |  Log filled/miss = "
                  << nLogOK << '/' << nLogMiss << '\n';
}



// ============================================================================
//  local helpers (only this .cc needs them – keep them ‘static inline’)
// ============================================================================
namespace
{

/* -------------------------------------------------------------------------- */
/* Verbosity access for helpers (0 = silent)                                  */
/* -------------------------------------------------------------------------- */
static inline uint64_t vbLevel()
{
  return PositionDependentCorrection::getVerbosityLevel();   // ← fixed access
}
/* -------------------------------------------------------------------------- */
/* 1.  CG‑ and block‑level coordinate conversions                             */
/* -------------------------------------------------------------------------- */

// Tower‑grid (xT,yT)  ➜  global φ  (already at shower depth)
static inline float cg2GlobalPhi(BEmcRecCEMC* rec,
                                 float eReco,
                                 float xT,
                                 float yT)
{
  float gx, gy, gz;
  rec->Tower2Global(eReco, xT, yT, gx, gy, gz);
  const float phiSD = std::atan2(gy, gx);
  if (vbLevel() > 2)
    std::cout << "[cg2GlobalPhi] (xT,yT)=(" << xT << ',' << yT
              << ")  →  (gx,gy)=(" << gx << ',' << gy
              << ")  φ_SD=" << phiSD << '\n';
  return phiSD;
}

// Tower‑grid (xT,yT)  ➜  η at shower depth
static inline float cg2ShowerEta(BEmcRecCEMC* rec,
                                 float eReco,
                                 float xT,
                                 float yT,
                                 float vtxZ)
{
  float gx, gy, gz;
  rec->Tower2Global(eReco, xT, yT, gx, gy, gz);
  const float etaSD = std::asinh((gz - vtxZ) / std::hypot(gx, gy));
  if (vbLevel() > 2)
    std::cout << "[cg2ShowerEta] (gx,gy,gz)=(" << gx << ',' << gy << ',' << gz
              << ")  vtxZ=" << vtxZ << "  →  η_SD=" << etaSD << '\n';
  return etaSD;
}

// Front‑face φ  ➜  φ at shower depth
static inline float front2ShowerPhi(BEmcRecCEMC* rec,
                                    float  eReco,
                                    double rFront,
                                    double zFront,
                                    float  phiFront,
                                    int    ixFine,
                                    int    iyFine)
{
  float xSD, ySD, zSD;
  rec->CorrectShowerDepth(ixFine, iyFine, eReco,
                          rFront * std::cos(phiFront),
                          rFront * std::sin(phiFront),
                          zFront,
                          xSD, ySD, zSD);
  const float phiSD = std::atan2(ySD, xSD);
  if (vbLevel() > 2)
    std::cout << "[front2ShowerPhi] φ_front=" << phiFront
              << "  →  (xSD,ySD)=(" << xSD << ',' << ySD
              << ")  φ_SD=" << phiSD << '\n';
  return phiSD;
}

// Front‑face (η,φ)  ➜  η at shower depth
static inline float front2ShowerEta(BEmcRecCEMC* rec,
                                    float  eReco,
                                    double rFront,
                                    int    ixFine,
                                    int    iyFine,
                                    float  etaFront,
                                    float  phiFront,
                                    float  vtxZ)
{
  const double thetaF = 2.0 * std::atan(std::exp(-etaFront));
  const double zF     = rFront / std::tan(thetaF);

  float xSD, ySD, zSD;
  rec->CorrectShowerDepth(ixFine, iyFine, eReco,
                          rFront * std::cos(phiFront),
                          rFront * std::sin(phiFront),
                          zF,
                          xSD, ySD, zSD);
  const float etaSD = std::asinh((zSD - vtxZ) / std::hypot(xSD, ySD));
  if (vbLevel() > 2)
    std::cout << "[front2ShowerEta] (η_front,φ_front)=(" << etaFront << ','
              << phiFront << ")  →  η_SD=" << etaSD << '\n';
  return etaSD;
}

/* -------------------------------------------------------------------------- */
/* 2.  Cluster centre‑of‑gravity (tower‑moment analysis)                      */
/* -------------------------------------------------------------------------- */
static inline bool
clusterCentreOfGravity(const RawCluster*      clus,
                       RawTowerGeomContainer* geom,
                       const BEmcRecCEMC*     bemc,
                       float&                 xCG,
                       float&                 yCG)
{
  if (!clus || !geom || !bemc) return false;

  std::vector<EmcModule> hits;
  hits.reserve(clus->getNTowers());
  const int nPhiRec = bemc->GetNx();

  for (auto it = clus->get_towers().first; it != clus->get_towers().second; ++it)
    if (const auto* tg = geom->get_tower_geometry(it->first))
      hits.emplace_back(tg->get_binphi() + tg->get_bineta() * nPhiRec,
                        it->second,
                        0.f);

  if (hits.empty()) return false;

  float xx, yy, xy, eTot;
  bemc->Momenta(&hits, eTot, xCG, yCG, xx, yy, xy);

  if (vbLevel() > 2)
    std::cout << "[clusterCentreOfGravity] NTowers=" << hits.size()
              << "  ΣE=" << eTot
              << "  →  (xCG,yCG)=(" << xCG << ',' << yCG << ")\n";
  return true;
}

/* -------------------------------------------------------------------------- */
/* 3.  Light‑weight flavour records                                           */
/* -------------------------------------------------------------------------- */
struct PhiRec { const char* tag; float loc; float phi; float d {0.f}; };
struct EtaRec { const char* tag; float loc; float eta; float d {0.f}; };

} // anonymous namespace
// ============================================================================


// ============================================================================
//  Δφ – five flavours
// ============================================================================
void PositionDependentCorrection::fillDPhiAllVariants(
        RawCluster*                   clus,
        const TLorentzVector&         /*recoPhoton*/,
        const TLorentzVector&         truthPhoton,
        const std::pair<float,float>& blkCoord,     // (ηloc , φloc)
        int                           blkPhiCoarse,
        float                         /*vtx_z*/,
        TH1F*                         hRAW [N_Ebins],
        TH1F*                         hCP  [N_Ebins] )
{
  /* ── A) fast guards ───────────────────────────────────────────────── */
  if (!clus || !m_geometry || !m_bemcRec) return;
  const int vb = Verbosity();          // snapshot once – cheaper

  /* ── B) slice lookup ─────────────────────────────────────────────── */
  const float eReco = clus->get_energy();
  const int   iE    = getEnergySlice(eReco);
  if (iE < 0 || iE >= N_Ebins) return;

  if (vb > 0)
    std::cout << "\n[fillDPhiAllVariants]  →  slice=" << iE
              << "  E=" << eReco << " GeV\n";

  /* ── C) cluster centre of gravity (tower grid) ───────────────────── */
  float xCG{}, yCG{};
  if (!clusterCentreOfGravity(clus, m_geometry, m_bemcRec, xCG, yCG))
  {
    if (vb > 1) std::cout << "  ‑ CG calculation failed – abort path\n";
    return;
  }
  if (vb > 2)
    std::cout << "  CG(x,y)=(" << xCG << ',' << yCG << ")\n";

  /* ── D) lead‑tower geometry (needed for PDC path) ────────────────── */
  const auto lead =
        std::max_element(clus->get_towers().first, clus->get_towers().second,
                         [](auto& a, auto& b){ return a.second < b.second; });
  if (lead == clus->get_towers().second) return;
  const auto* tgLead = m_geometry->get_tower_geometry(lead->first);
  if (!tgLead) return;

  const double rFront = tgLead->get_center_radius();
  const double zFront = tgLead->get_center_z();
  const int    ixFine = tgLead->get_binphi();
  const int    iyFine = tgLead->get_bineta();

  if (vb > 2)
    std::cout << "  lead‑tower  rFront=" << rFront
              << "  zFront=" << zFront
              << "  ixFine="  << ixFine
              << "  iyFine="  << iyFine << '\n';

  /* ── E) build the five flavour records ───────────────────────────── */
  PhiRec rec[5];
  const float phiTruth = TVector2::Phi_mpi_pi(truthPhoton.Phi());

  //---------------- (1) CLUS‑RAW -------------------------------------
  rec[0] = {"CLUS-RAW", xCG,
            cg2GlobalPhi(m_bemcRec, eReco, xCG, yCG)};

  //---------------- (2) CLUS‑CP  -------------------------------------
  float xCP = xCG, yCP = yCG;
  m_bemcRec->CorrectPosition(eReco, xCG, yCG, xCP, yCP);
  rec[1] = {"CLUS-CP", xCP,
            cg2GlobalPhi(m_bemcRec, eReco, xCP, yCP)};

  //---------------- (3) CLUS‑BCORR -----------------------------------
  {
    const float nX   = static_cast<float>(m_bemcRec->GetNx());
    const float bΦ   = m_bValsPhi[iE];
    float       locB = xCG;

    if (bΦ > 1e-9F)
    {
      const int   ix0 = int(xCG + .5F);
      const float off = xCG - ix0;
      locB = ix0 + doPhiBlockCorr(off, bΦ);
      while (locB < -0.5F)     locB += nX;
      while (locB >= nX-.5F)   locB -= nX;
    }
    rec[2] = {"CLUS-BCORR", locB,
        cg2GlobalPhi(m_bemcRec, eReco,  locB, yCG)};
  }

  //---------------- (4) PDC‑RAW  -------------------------------------
  {
    const float phiFront = convertBlockToGlobalPhi(blkPhiCoarse,
                                                   blkCoord.second);
    const float phiSD    = front2ShowerPhi(m_bemcRec, eReco,
                                           rFront, zFront,
                                           phiFront,
                                           ixFine, iyFine);
    rec[3] = {"PDC-RAW", blkCoord.second, phiSD};
  }

  //---------------- (5) PDC‑CORR -------------------------------------
  {
    const float bΦ = m_bValsPhi[iE];
    if (isFitDoneForPhi && bΦ > 1e-9F)
    {
      float loc = PDC::Geo::phi::undoAsh(blkCoord.second, bΦ);
      if (loc > 0.5F && loc <= 1.5F) ++g_nCorrPhiRight;
      int   blk = blkPhiCoarse;
      PDC::Geo::foldAndStep<true>(loc, blk, PDC::Geo::kCoarsePhiBins);

      const float phiFront = convertBlockToGlobalPhi(blk, loc);
      const float phiSD    = front2ShowerPhi(m_bemcRec, eReco,
                                             rFront, zFront,
                                             phiFront, ixFine, iyFine);
      rec[4] = {"PDC-CORR", loc, phiSD};
    }
    else
      rec[4] = {"PDC-CORR", blkCoord.second, rec[3].phi};
  }

  /* ── F) residuals and histogram filling ─────────────────────────── */
  auto wrap = [](float d){ return TVector2::Phi_mpi_pi(d); };

  for (auto& r : rec)
  {
    r.d = wrap(r.phi - phiTruth);
    if (!std::isfinite(r.phi)) ++g_nanPhi;

    if (vb > 3)
      std::cout << "    " << r.tag << "  loc=" << r.loc
                << "  φ_SD=" << r.phi
                << "  Δφ="   << r.d << '\n';
  }

#define HFill(H,V)  do{ if((H)&&std::isfinite(V)) (H)->Fill(V); }while(0)
  HFill(hRAW [iE], rec[0].d);
  HFill(hCP  [iE], rec[1].d);
  HFill(h_phi_diff_cpBcorr_E [iE], rec[2].d);
  HFill(h_phi_diff_raw_E     [iE], rec[3].d);
  HFill(h_phi_diff_corrected_E[iE], rec[4].d);
#undef  HFill

  /* ── G) “winner” bookkeeping ────────────────────────────────────── */
  int best5 = 0; for (int i=1;i<5;++i)
    if (std::fabs(rec[i].d) < std::fabs(rec[best5].d)) best5 = i;
  switch(best5){
    case 0: ++m_phiWinCLUSraw;   break;
    case 1: ++m_phiWinCLUScp;    break;
    case 2: ++m_phiWinCLUSbcorr; break;
    case 3: ++m_phiWinPDCraw;    break;
    case 4: ++m_phiWinPDCcorr;   break;
  }

  int best3 = 0; for (int i=1;i<3;++i)
    if (std::fabs(rec[i].d) < std::fabs(rec[best3].d)) best3 = i;
  switch(best3){
    case 0: ++m_nWinRAW;   break;
    case 1: ++m_nWinCP;    break;
    case 2: ++m_nWinBCorr; break;
  }

  /* ── H) per‑event verbose table (verbosity ≥ 2) ─────────────────── */
  if (vb >= 2)
  {
    std::cout << "  ── Δφ summary table ────────────────────────────────\n"
              << "  tag         loc      φ_SD(rad)   Δφ(rad)   |Δφ|\n";
    for (const auto& r : rec)
      std::cout << std::setw(10)<<r.tag<<"  "
                << std::setw(7)<<std::fixed<<std::setprecision(3)<<r.loc<<"  "
                << std::setw(11)<<r.phi<<"  "
                << std::setw(9)<<r.d<<"  "
                << std::setw(7)<<std::fabs(r.d)<<"\n";
    std::cout << "  WINNER: " << rec[best5].tag << '\n'
              << "  ────────────────────────────────────────────────\n";
  }
}


// ============================================================================
//  Δη – five flavours
// ============================================================================
void PositionDependentCorrection::fillDEtaAllVariants(
        RawCluster*                   clus,
        const TLorentzVector&         /*recoPhoton*/,
        const TLorentzVector&         truthPhoton,
        const std::pair<float,float>& blkCoord,      // (ηloc , φloc)
        int                           blkEtaCoarse,  // coarse‑η index
        int                           blkPhiCoarse,  // coarse‑φ index
        float                         vtxZ,          // primary‑vertex z  [cm]
        TH1F*                         hRAW [N_Ebins],
        TH1F*                         hCP  [N_Ebins],
        bool                          fillGlobal )
{
  /* ─────────────────────────────────────────── guards ────────────── */
  if (!clus || !m_geometry || !m_bemcRec) return;
  const int vb = Verbosity();

  const float eReco = clus->get_energy();
  const int   iE    = getEnergySlice(eReco);
  if (iE < 0 || iE >= N_Ebins) return;

  if (vb > 0)
    std::cout << "\n[fillDEtaAllVariants] slice=" << iE
              << "  E=" << eReco << " GeV\n";

  /* ────────────────────────────────── CG of the cluster ──────────── */
  float xCG{}, yCG{};
  if (!clusterCentreOfGravity(clus, m_geometry, m_bemcRec, xCG, yCG))
  {
    if (vb > 1) std::cout << "  ‑ CG calculation failed – abort path\n";
    return;
  }
  if (vb > 2)
    std::cout << "  CG(x,y)=(" << xCG << ',' << yCG << ")\n";

  /* ────────────────────────────────── lead‑tower geom ────────────── */
  const auto lead =
        std::max_element(clus->get_towers().first, clus->get_towers().second,
                         [](auto&a,auto&b){return a.second<b.second;});
  if (lead == clus->get_towers().second) return;
  const auto* tgLead = m_geometry->get_tower_geometry(lead->first);
  if (!tgLead) return;

  const double rFront = tgLead->get_center_radius();
  const int    ixFine = tgLead->get_binphi();

  const float phiFrontRaw =
        convertBlockToGlobalPhi(blkPhiCoarse, blkCoord.second);

  /* ─────────────────────────────── build the five variants ──────── */
  EtaRec rec[5];
  const float etaTruth = truthPhoton.Eta();

  // 1) CLUS‑RAW
  rec[0] = {"CLUS-RAW", yCG,
            cg2ShowerEta(m_bemcRec, eReco, xCG, yCG, vtxZ)};

  // 2) CLUS‑CP
  float xCP = xCG, yCP = yCG;
  m_bemcRec->CorrectPosition(eReco, xCG, yCG, xCP, yCP);
  rec[1] = {"CLUS-CP", yCP,
            cg2ShowerEta(m_bemcRec, eReco, xCP, yCP, vtxZ)};

  // 3) CLUS‑BCORR
  {
    const float bEta = m_bValsEta[iE];
    float loc = yCG;
    if (bEta > 1e-9F)
    {
      const int   iy0 = int(yCG + .5F);
      const float off = yCG - iy0;
      loc = iy0 + doEtaBlockCorr(off, bEta);
    }
    rec[2] = {"CLUS-BCORR", loc,
        cg2ShowerEta(m_bemcRec, eReco, xCG,  loc, vtxZ)}; // consistent pair
  }

  // 4) PDC‑RAW
  {
    const float etaFront =
        convertBlockToGlobalEta(blkEtaCoarse, blkCoord.first);
    const int iyFine =
        PDC::Geo::fineEtaIdx(blkEtaCoarse, blkCoord.first);

    const float etaSD =
        front2ShowerEta(m_bemcRec, eReco,
                        rFront, ixFine, iyFine,
                        etaFront, phiFrontRaw,
                        vtxZ);
    rec[3] = {"PDC-RAW", blkCoord.first, etaSD};
  }

  // 5) PDC‑CORR
  {
    const float bEta = m_bValsEta[iE];
    if (isFitDoneForEta && bEta > 1e-9F)
    {
      float loc = PDC::Geo::eta::undoAsh(blkCoord.first, bEta);
      int   blk = blkEtaCoarse;
      PDC::Geo::foldAndStep<false>(loc, blk, PDC::Geo::kCoarseEtaBins);

      const float etaFront =
          convertBlockToGlobalEta(blk, loc);
      const int iyFine =
          PDC::Geo::fineEtaIdx(blk, loc);

      const float etaSD =
          front2ShowerEta(m_bemcRec, eReco,
                          rFront, ixFine, iyFine,
                          etaFront, phiFrontRaw,
                          vtxZ);
      rec[4] = {"PDC-CORR", loc, etaSD};
    }
    else
      rec[4] = {"PDC-CORR", blkCoord.first, rec[3].eta};
  }

  /* ───────────────── residuals + histogramming ───────────────────── */
  for (auto& r : rec)
  {
    r.d = r.eta - etaTruth;
    if (!std::isfinite(r.eta)) ++g_nanEta;

    if (vb > 3)
      std::cout << "    " << r.tag << "  loc=" << r.loc
                << "  η_SD=" << r.eta
                << "  Δη="   << r.d << '\n';
  }

#define HFill(H,V)  do{ if((H)&&std::isfinite(V)) (H)->Fill(V); }while(0)
  if (fillGlobal)
  {
    HFill(hRAW [iE], rec[0].d);
    HFill(hCP  [iE], rec[1].d);
    HFill(h_eta_diff_cpBcorr_E [iE], rec[2].d);
    HFill(h_eta_diff_raw_E     [iE], rec[3].d);
    HFill(h_eta_diff_corrected_E[iE], rec[4].d);
  }
  const int iVz = getVzSlice(std::fabs(vtxZ));
  if (iVz >= 0 && iVz < N_VzBins)
  {
    HFill(h_eta_diff_raw_E_vz     [iE][iVz], rec[0].d);
    HFill(h_eta_diff_corrected_E_vz[iE][iVz], rec[4].d);
    HFill(h_eta_diff_cpRaw_E_vz   [iE][iVz], rec[1].d - rec[0].d);
    HFill(h_eta_diff_cpCorr_E_vz  [iE][iVz], rec[1].d - rec[4].d);
    HFill(h_eta_diff_cpBcorr_E_vz [iE][iVz], rec[1].d - rec[2].d);
  }
#undef  HFill

  /* ───────────────── winner bookkeeping ──────────────────────────── */
  int best5 = 0; for(int i=1;i<5;++i)
    if (std::fabs(rec[i].d) < std::fabs(rec[best5].d)) best5 = i;
  switch(best5){
    case 0: ++m_etaWinCLUSraw;   break;
    case 1: ++m_etaWinCLUScp;    break;
    case 2: ++m_etaWinCLUSbcorr; break;
    case 3: ++m_etaWinPDCraw;    break;
    case 4: ++m_etaWinPDCcorr;   break;
  }

  int best3 = 0; for(int i=1;i<3;++i)
    if (std::fabs(rec[i].d) < std::fabs(rec[best3].d)) best3 = i;
  switch(best3){
    case 0: ++m_nWinRAW_Eta;   break;
    case 1: ++m_nWinCP_Eta;    break;
    case 2: ++m_nWinBCorr_Eta; break;
  }

  /* ───────────────── per‑event verbose table (verbosity ≥ 2) ─────── */
  if (vb >= 2)
  {
    std::cout << "  ── Δη summary table (|z|="<<std::fabs(vtxZ)<<" cm) ─────\n"
              << "  tag         loc      η_SD       Δη        |Δη|\n";
    for (const auto& r : rec)
      std::cout << std::setw(10)<<r.tag<<"  "
                << std::setw(7)<<std::fixed<<std::setprecision(3)<<r.loc<<"  "
                << std::setw(10)<<r.eta<<"  "
                << std::setw(10)<<r.d<<"  "
                << std::setw(7)<<std::fabs(r.d)<<"\n";
    std::cout << "  WINNER: " << rec[best5].tag << '\n'
              << "  ────────────────────────────────────────────────\n";
  }
}



/* =======================================================================
 *  Fills raw & corrected 2×2‑block‑coordinate histograms
 * ===================================================================== */
void PositionDependentCorrection::fillBlockCoordinateHistograms(
        const std::pair<float,float>& blkCoord,  // (ηloc , φloc)
        int   blkEtaCoarse,                      // coarse‑η index
        int   blkPhiCoarse,                      // coarse‑φ index
        float clusE,
        int   iEbin,
        [[maybe_unused]] std::size_t nTowers )
{
    constexpr double kFillW = 1.0;

    /* ------------ A) raw fill ---------------------------------------- */
    const float locEta = blkCoord.first;
    const float locPhi = blkCoord.second;
    h3_cluster_block_cord_E->Fill(locEta, locPhi, clusE, kFillW);

    /* ------------ B) optional Ash-correction ------------------------- */
    float corrEta = locEta;
    float corrPhi = locPhi;

    /* -------- φ (periodic) ------------------------------------------ */
    if (isFitDoneForPhi && iEbin >= 0 && iEbin < N_Ebins)
    {
        const float bPhi = m_bValsPhi[iEbin];
        if (bPhi > 1e-9F)
        {
            auto p        = PDC::Geo::undoAshAndReindex<true>(corrPhi,
                                                              blkPhiCoarse,
                                                              bPhi);
            corrPhi       = p.loc;      // folded, legacy‑exact
            blkPhiCoarse  = p.blk;      // keep index coherent
        }
    }

    /* -------- η (finite) -------------------------------------------- */
    if (isFitDoneForEta && iEbin >= 0 && iEbin < N_Ebins)
    {
        const float bEta = m_bValsEta[iEbin];
        if (bEta > 1e-9F)
        {
            auto e        = PDC::Geo::undoAshAndReindex<false>(corrEta,
                                                               blkEtaCoarse,
                                                               bEta);
            corrEta       = e.loc;
            blkEtaCoarse  = e.blk;
        }
    }


    h3_cluster_block_cord_E_corrected
        ->Fill(corrEta, corrPhi, clusE, kFillW);

    if (corrPhi > 0.5F) ++g_nCorrPhiRight;   // QA counter
    /* ------------ C) occupancy statistics for End() summary ----------- */
    const int iEta = (corrEta > 0.5F) ? 1 : 0;      // 0 ↔ (-0.5…0.5] , 1 ↔ (0.5…1.5]
    const int iPhi = (corrPhi > 0.5F) ? 1 : 0;
    m_blkLocCount[iEta * 2 + iPhi]++;               // atomic – thread‑safe
}



// ============================================================================
//  Helper: handle all photon–truth matching and truth–QA histograms
//           (was the complete body of the previous “if (m_isSimulation)” block)
// ============================================================================
void PositionDependentCorrection::processSimulationTruthMatches(
        RawCluster*                     clus1,
        const TLorentzVector&           photon1,
        [[maybe_unused]] int                             lt_eta,
        int                             lt_phi,
        const std::pair<float,float>&   blkCoord,
        int                             blkEtaCoarse,
        int                             blkPhiCoarse,
        const std::vector<int>&         towerPhis,
        const std::vector<float>&       towerEs,
        float                           vtx_z,
        const std::vector<TLorentzVector>& truth_photons,
        const std::vector<TLorentzVector>& truth_meson_photons,
        /* out‑params ---------------------------------------------------- */
        bool&                           match1,
        TLorentzVector&                 ph1_trEtaPhi )
{
  const float clusE = photon1.E();   // convenience alias

  /* ------------------------------------------------------------------ */
  /* 1)  Primary‑photon matching  ------------------------------------- */
  /* ------------------------------------------------------------------ */
  for (auto& trPhoton : truth_photons)
  {
    const float dR    = photon1.DeltaR(trPhoton);          // ΔR(reco,truth)
    const float ratio = photon1.E() / trPhoton.E();        // Ereco / Etruth

    if (dR   > 0.03f)                    continue;
    if (ratio < 0.30f || ratio > 1.30f)  continue;

    const float dPhi = TVector2::Phi_mpi_pi(
                         photon1.Phi() - trPhoton.Phi() );

    h_delR_recTrth ->Fill(dR);

    if (Verbosity() > 0)
      std::cout << "[DEBUG] => cluster1 E=" << photon1.E()
                << " matched to truthE="    << trPhoton.E()
                << "  dR="   << dR
                << "  dPhi=" << dPhi
                << "  ratio="<< ratio << '\n';

    /* histogramming -------------------------------------------------- */
    h_matched_res      ->Fill(ratio, photon1.Eta());
    h_res_e            ->Fill(ratio, photon1.E());
    h_res              ->Fill(ratio);

    const int iLTeta = lt_eta;
    const int iLTphi = lt_phi;

    h_res_e_eta        ->Fill(ratio, trPhoton.E(), iLTeta);
    h_res_e_phi        ->Fill(ratio, trPhoton.E(), iLTphi);
    h_delEta_e_eta     ->Fill(photon1.Eta()-trPhoton.Eta(), trPhoton.E(), iLTeta);
    h_delR_e_eta       ->Fill(dR, trPhoton.E(), iLTeta);
    h_delPhi_e_eta     ->Fill(dPhi, trPhoton.E(), iLTeta);
    h_delPhi_e_phi     ->Fill(dPhi, trPhoton.E(), iLTphi);
    h_truthE           ->Fill(trPhoton.E());

    /* specialised QA helpers ---------------------------------------- */
    fillAshLogDx(clus1, photon1, trPhoton,
                 blkCoord, blkPhiCoarse,
                 towerPhis, towerEs);

    /* Δφ – five‑way */
    fillDPhiAllVariants( clus1,            // the RawCluster*
                           photon1,          // reco photon 4‑vector
                           trPhoton,         // truth photon 4‑vector
                           blkCoord,         // (ηloc,φloc)
                           blkPhiCoarse,     // coarse φ‑block
                           /*vtx_z=*/0.f,    // not used in φ path
                           h_phi_diff_cpRaw_E,           // RAW   per‑E slice array
                           h_phi_diff_cpCorr_E );        // CP    per‑E slice array


    /* Δη – five‑way (tight‑vz → global spectra on) */
    fillDEtaAllVariants( clus1,                 // RawCluster*
                           photon1,               // reco four‑vector
                           trPhoton,              // truth photon
                           blkCoord,              // (ηloc , φloc)
                           blkEtaCoarse, blkPhiCoarse,
                           vtx_z,                 // *actual* vertex‑Z
                           h_eta_diff_cpRaw_E,    // RAW  per‑E slice array
                           h_eta_diff_cpCorr_E ); // CP   per‑E slice array
      
    if (blkEtaCoarse >= 0 && blkEtaCoarse < NBinsBlock &&
        blkPhiCoarse >= 0 && blkPhiCoarse < NBinsBlock)
    {
      h_res_block_E[blkEtaCoarse][blkPhiCoarse]
        ->Fill(ratio, trPhoton.E());
    }
  }

  /* ------------------------------------------------------------------ */
  /* 2)  Meson‑decay photon matching  --------------------------------- */
  /* ------------------------------------------------------------------ */
  match1        = false;
  ph1_trEtaPhi  = TLorentzVector(0,0,0,0);

  for (auto& trPhot : truth_meson_photons)
  {
    const float dR2    = photon1.DeltaR(trPhot);
    const float ratio2 = photon1.E() / trPhot.E();

    if (dR2 < 0.03f && ratio2 > 0.7f && ratio2 < 1.5f)
    {
      ph1_trEtaPhi.SetPtEtaPhiE( clusE / TMath::CosH(trPhot.Eta()),
                                 trPhot.Eta(),
                                 trPhot.Phi(),
                                 clusE );

      if (Verbosity() > 0)
        std::cout << "[DEBUG] => meson‑decay match found, E="
                  << ph1_trEtaPhi.E() << "  eta="
                  << ph1_trEtaPhi.Eta() << "\n";

      match1 = true;
      break;    // only need the first successful match
    }
  }
}

// ============================================================================
// Helper • loop over ALL “other” clusters, build π0 candidates, apply cuts
//         and fill (reco + truth) histograms.
//
//    – Always fill h_mE_raw / h_mE_corr  (input for pass‑1 μ/σ fits)
//    – If m_massFitsDone==true apply slice‑specific μ±3σ window
//      and populate the block–space maps h_m_blk_raw / h_m_blk_corr.
// ============================================================================
void PositionDependentCorrection::processClusterPairs(
        RawClusterContainer               *clusterContainer,
        RawClusterContainer::ConstIterator cIt1,
        const CLHEP::Hep3Vector           &vertex,
        const TLorentzVector              &photon1,
        float                              clusE,
        int                                lt_eta,
        int                                lt_phi,
        int                                blkEtaCoarse,
        int                                blkPhiCoarse,
        const std::pair<float,float>      &blkCoord,      // (ηloc , φloc) – RAW
        float                              maxAlpha,
        float                              ptMaxCut,
        float                              pt2ClusCut,
        float                              pi0ptcut,
        float                              weight,
        bool                               match1,
        const TLorentzVector              &ph1_trEtaPhi,
        const std::vector<TLorentzVector> &truth_meson_photons,
        bool                               isSimulation)
{
  const bool vb = (Verbosity() > 4);            // detailed debug switch
  (void) lt_phi; 

  RawClusterContainer::ConstRange    cRange2 = clusterContainer->getClusters();
  for (auto cIt2 = cRange2.first; cIt2 != cRange2.second; ++cIt2)
  {
    // ----------------------------------------------------------------------
    // 0)  Guard: do not pair a cluster with itself
    // ----------------------------------------------------------------------
    if (cIt2 == cIt1) continue;
    RawCluster *clus2 = cIt2->second;
    if (!clus2)       continue;

    // ----------------------------------------------------------------------
    // 1)  Basic 4‑vector of the second cluster and quality cuts
    // ----------------------------------------------------------------------
    const CLHEP::Hep3Vector eVec2 = RawClusterUtility::GetEVec(*clus2, vertex);

    const float clus2E    = eVec2.mag();
    const float clus2Pt   = eVec2.perp();
    const float clus2Eta  = eVec2.pseudoRapidity();
    const float clus2Phi  = eVec2.phi();
    const float clus2Chi2 = clus2->get_chi2();

    if (clus2Pt <  pt2ClusCut || clus2Pt >  ptMaxCut)   continue;
    if (clus2Chi2 > 1.0e4)                              continue;

    const float alpha = std::fabs(clusE - clus2E) / (clusE + clus2E);
    if (alpha   >  maxAlpha)                           continue;

    TLorentzVector photon2;  photon2.SetPtEtaPhiE(clus2Pt, clus2Eta,
                                                  clus2Phi, clus2E);

    // ----------------------------------------------------------------------
    // 2)  Optional meson‑decay truth matching (unchanged)
    // ----------------------------------------------------------------------
    bool           match2       = false;
    TLorentzVector ph2_trEtaPhi(0,0,0,0);

    if (isSimulation)
    {
      for (const auto &trPhot : truth_meson_photons)
      {
        const float dR    = photon2.DeltaR(trPhot);
        const float ratio = photon2.E() / trPhot.E();
        if (dR < 0.02f && ratio > 0.7f && ratio < 1.5f)
        {
          ph2_trEtaPhi.SetPtEtaPhiE( clus2E / TMath::CosH(trPhot.Eta()),
                                     trPhot.Eta(), trPhot.Phi(), clus2E );
          if (match1) match2 = true;
          break;
        }
      }
    }

    // ----------------------------------------------------------------------
    // 3)  Build the π0 candidate (reco‑based) and truth variant
    // ----------------------------------------------------------------------
    const TLorentzVector pi0Reco   = photon1 + photon2;
    const TLorentzVector pi0TruthK = ph1_trEtaPhi + ph2_trEtaPhi;

    if (isSimulation && pi0Reco.Pt() < pi0ptcut)         continue;

    // ----------------------------------------------------------------------
    // 4)  Energy‑slice index (use **leading‑γ** energy!)
    // ----------------------------------------------------------------------
    const float  eLead  = std::max(clusE, clus2E);
    const int    iSlice = getEnergySlice(eLead);     // −1 ⇒ outside table

    // ----------------------------------------------------------------------
    // 5)  -- ALWAYS -- fill slice histograms (pass‑1 input)
    // ----------------------------------------------------------------------
    if (h_mE_raw )
      h_mE_raw ->Fill(pi0Reco.M(), eLead);
    if (h_mE_corr)
      h_mE_corr->Fill(pi0Reco.M(), eLead);

    // ----------------------------------------------------------------------
    // 6)  Conditional (pass‑2) block‑map filling
    // ----------------------------------------------------------------------
    if (m_massFitsDone && iSlice >= 0 && iSlice < N_Ebins)
    {
      /* ---------- 6a) RAW window ------------------------------------- */
      const auto& winR = m_winRaw [iSlice];
      if (winR.sigma > 1e-6f)
      {
        const float lo = winR.mu - 3.f*winR.sigma;
        const float hi = winR.mu + 3.f*winR.sigma;
        if (pi0Reco.M() >= lo && pi0Reco.M() <= hi)
          if (h_m_blk_raw)
            h_m_blk_raw->Fill(blkCoord.first, blkCoord.second, pi0Reco.M());
      }

      /* ---------- 6b) CORR window  (needs corrected coords) ---------- */
      const auto& winC = m_winCorr[iSlice];
      if (winC.sigma > 1e-6f)
      {
        float locEtaCorr = blkCoord.first;
        float locPhiCorr = blkCoord.second;

        if (isFitDoneForEta && m_bValsEta[iSlice] > 1e-9f)
          locEtaCorr = doEtaBlockCorr(locEtaCorr, m_bValsEta[iSlice]);
        if (isFitDoneForPhi && m_bValsPhi[iSlice] > 1e-9f)
          locPhiCorr = doPhiBlockCorr(locPhiCorr, m_bValsPhi[iSlice]);

        const float lo = winC.mu - 3.f*winC.sigma;
        const float hi = winC.mu + 3.f*winC.sigma;
        if (pi0Reco.M() >= lo && pi0Reco.M() <= hi)
          if (h_m_blk_corr)
            h_m_blk_corr->Fill(locEtaCorr, locPhiCorr, pi0Reco.M());
      }
    }

    // ----------------------------------------------------------------------
    // 7)  Legacy histogramming (unchanged)
    // ----------------------------------------------------------------------
    h_pt1            ->Fill(photon1.Pt());
    h_pt2            ->Fill(photon2.Pt());
    h_InvMass        ->Fill(pi0Reco.M());
    h_InvMass_w      ->Fill(pi0Reco.M(), weight);
    h_mass_eta_lt    ->Fill(pi0Reco.M(), lt_eta);
    h_mass_eta_lt_rw ->Fill(pi0Reco.M(), lt_eta, weight);
    h_m_pt_eta       ->Fill(pi0Reco.M(), pi0Reco.E(), lt_eta);

    if (blkEtaCoarse >= 0 && blkEtaCoarse < NBinsBlock &&
        blkPhiCoarse >= 0 && blkPhiCoarse < NBinsBlock)
    {
      h_mass_block_pt[blkEtaCoarse][blkPhiCoarse]
            ->Fill(pi0Reco.M(), pi0Reco.E());
    }

    if (isSimulation && match2 && pi0TruthK.M() > 1e-3f)
    {
      h_m_ptTr_eta      ->Fill(pi0Reco.M(), pi0TruthK.E(), lt_eta);
      h_m_ptTr_eta_trKin->Fill(pi0TruthK.M(), pi0TruthK.E(), lt_eta);
    }

    if (vb)
      std::cout << "[PDC::processClusterPairs]  "
                << "π0 mass=" << pi0Reco.M()
                << "  slice=" << iSlice
                << "  filled (RAW/CORR) = "
                << (m_massFitsDone ? "yes" : "no") << '\n';
  }
}


// ============================================================================
void PositionDependentCorrection::finalClusterLoop(
    PHCompositeNode* topNode,
    RawClusterContainer* clusterContainer,
    float vtx_z,
    const std::vector<TLorentzVector>& truth_photons,
    const std::vector<TLorentzVector>& truth_meson_photons,
    float tower_tot_e,
    float max_nClusCount,
    int nClusCount,
    float maxAlpha,
    float ptMaxCut,
    float pt1ClusCut,
    float pt2ClusCut,
    float pi0ptcut,
    float weight)
{
    
  /* ------------------------------------------------------------------
     * Trigger‑based event pre‑selection   (data only)
  * ------------------------------------------------------------------ */
  bool analyseThisEvent = true;              // default for simulation

  if (!m_isSimulation)                       // running on DATA
  {
      analyseThisEvent = false;                // require an accepted trigger

      /* 1) make sure the TriggerAnalyzer exists */
      if (!trigAna)
      {
        std::cerr << "[PDC]  ERROR – TriggerAnalyzer pointer is NULL\n";
        return;                                // drop event gracefully
      }

      /* 2) decode and test the triggers stored in triggerNameMap */
      trigAna->decodeTriggers(topNode);

      for (const auto& kv : triggerNameMap)        // loop over DB names
        if (trigAna->didTriggerFire(kv.first))
        {
          analyseThisEvent = true;  break;         // any match ⇒ keep event
        }

      if (Verbosity() > 1)
        std::cout << "[PDC]  trigger gate ⇒ "
                  << (analyseThisEvent ? "ACCEPT" : "REJECT")
                  << " event\n";
  }
    
  if (Verbosity() > 0)
  {
    std::cout << "\n\n"
              << ANSI_BOLD << "[DEBUG] finalClusterLoop() ENTER" << ANSI_RESET << "\n"
              << "    tower_tot_e  = " << tower_tot_e << "\n"
              << "    nClusCount   = " << nClusCount << "\n"
              << "    max_nClusCount = " << max_nClusCount << "\n"
              << "    pt1ClusCut   = " << pt1ClusCut << ", pt2ClusCut = " << pt2ClusCut << "\n"
              << "    maxAlpha     = " << maxAlpha << "\n"
              << "    pi0ptcut     = " << pi0ptcut << "\n"
              << "    maxPtCut     = " << ptMaxCut << "\n"
              << "    weight       = " << weight << "\n"
              << "    vtx_z        = " << vtx_z << std::endl;
  }

  h_nclusters->Fill(nClusCount);
    

  // If the number of clusters is too large => skip
  if (nClusCount > max_nClusCount)
  {
    if (Verbosity() > 0)
    {
      std::cout << ANSI_BOLD << ANSI_YELLOW
                << "[DEBUG] finalClusterLoop: nClusCount=" << nClusCount
                << " > max_nClusCount=" << max_nClusCount
                << " => SKIPPING event"
                << ANSI_RESET << std::endl;
    }
    return;
  }

  RawClusterContainer::ConstRange clusterRange = clusterContainer->getClusters();
  RawClusterContainer::ConstIterator cIt1, cIt2;

  // For optional pT smearing
  float smear = 0.00;

  // Will hold booleans to see if first/second matched to meson photons
  bool match1 = false;
  bool match2 [[maybe_unused]] = false;

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] finalClusterLoop: Starting OUTER loop over clusters...\n";
  }
  // -----------------------------------------------------------------
  // 1) Outer loop over clusters
  // -----------------------------------------------------------------
  for (cIt1 = clusterRange.first; cIt1 != clusterRange.second; ++cIt1)
  {
    RawCluster* clus1 = cIt1->second;
    if (!clus1) continue;

    CLHEP::Hep3Vector vertex(0,0,vtx_z);

    CLHEP::Hep3Vector E_vec_1 = RawClusterUtility::GetEVec(*clus1, vertex);
      

    if (std::isnan(E_vec_1.mag()))
    {
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG] finalClusterLoop: WARNING - E_vec_1.mag() is NaN!\n"
                    << " -> cluster ID=" << clus1->get_id()
                    << "  vtx_z=" << vtx_z
                    << std::endl;
        }
        continue;
      }
      else if (E_vec_1.mag() < 1e-9)
      {
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG] finalClusterLoop: E_vec_1.mag() < 1e-9!"
                    << " -> Possibly zero energy cluster?\n"
                    << " -> cluster ID=" << clus1->get_id()
                    << "  vtx_z=" << vtx_z
                    << std::endl;
        }
        continue;
    }
    float clusE   = E_vec_1.mag();
    float clusEta = E_vec_1.pseudoRapidity();
    float clusPhi = E_vec_1.phi();
    float clusPt  = E_vec_1.perp();
    float clusChi2= clus1->get_chi2();

    if (Verbosity() > 10)
    {
        std::cout << "[DEBUG] finalClusterLoop: cluster ID=" << clus1->get_id()
                  << " => E="   << clusE
                  << "  Pt="    << clusPt
                  << "  Eta="   << clusEta
                  << "  Phi="   << clusPhi
                  << "  Chi2="  << clusChi2
                  << std::endl;
    }
    clusPt *= rnd->Gaus(1, smear);

    int lt_eta = clus1->get_lead_tower().first;
    int lt_phi = clus1->get_lead_tower().second;
    h2_chi2_tot_etaPhi->Fill(lt_eta, lt_phi);

    if (clusE < 0.1f)
    {
      if (Verbosity() > 11)
      {
        std::cout << "[DEBUG]  => Skipping cluster with E=" << clusE
                  << " (<0.1). \n";
      }
      continue;
    }
      constexpr float kChi2Cut = 10.0f;   // keep the magic number centralised

      if (clusChi2 > kChi2Cut)
      {
        // ---- χ² QA (numerator & profile) ----
        h2_chi2_rej_etaPhi->Fill(lt_eta, lt_phi);       // rejected count
        p_chi2_pass_etaPhi ->Fill(lt_eta, lt_phi, 0.0); // pass-flag = 0
        if (Verbosity() > 5)
          std::cout << "[χ²-CUT]  cluster Chi2=" << clusChi2
                    << "  (lead η,φ)=" << lt_eta << ',' << lt_phi
                    << "  REJECTED\n";
        continue;
      }
      else
      {
        // cluster is accepted
        p_chi2_pass_etaPhi->Fill(lt_eta, lt_phi, 1.0);  // pass-flag = 1
      }
    // 3) If lead tower index is out of range
    if (lt_eta > 95)
    {
      if (Verbosity() > 5)
      {
        std::cout << "[DEBUG]  => Skipping cluster lead_tower eta="
                  << lt_eta << " (>95). \n";
      }
      continue;
    }
    h_clusE->Fill(clusE);

    RawCluster::TowerConstRange towerCR = clus1->get_towers();
    std::vector<int> towerEtas, towerPhis;
    std::vector<float> towerEs;

    int nTow=0;
    for (auto tIter = towerCR.first; tIter != towerCR.second; ++tIter)
    {
      nTow++;
      float twE = tIter->second;

      RawTowerDefs::keytype tk = tIter->first;
      int iEta = m_geometry->get_tower_geometry(tk)->get_bineta();
      int iPhi = m_geometry->get_tower_geometry(tk)->get_binphi();

      towerEtas.push_back(iEta);
      towerPhis.push_back(iPhi);
      towerEs.push_back(twE);
    }
    h_clusE_nTow->Fill(clusE, nTow);
      
    TLorentzVector photon1;
    photon1.SetPtEtaPhiE(clusPt, clusEta, clusPhi, clusE);

    if (clusPt < pt1ClusCut || clusPt > ptMaxCut)
    {
      if (Verbosity() > 11)
      {
        std::cout << "[DEBUG]  => cluster #1 pT=" << clusPt
                  << " fails cut ( <" << pt1ClusCut << " or >" << ptMaxCut
                  << "). \n";
      }
      continue;
    }

    // -----------------------------------------------------------------
    // (C) Determine block coordinates (coarse & local)
    // -----------------------------------------------------------------
    /* ➊ build LocalCoord in one line – no side‑effects */
    PDC::Geo::LocalCoord loc =
              PDC::Geo::computeLocal(towerEtas, towerPhis, towerEs);

    /* ➋ copy what downstream code still expects … */
    int   blkPhiCoarse = loc.blkPhi;      // 0 … 127
    int   blkEtaCoarse = loc.blkEta;      // 0 …  47
    auto  blkCoord     = std::make_pair(loc.locEta, loc.locPhi);

    /* 4) find energy slice ........................................... */
    const int iEbin = getEnergySlice( clusE );


    if (analyseThisEvent)
        fillBlockCoordinateHistograms( blkCoord,
                                       blkEtaCoarse,
                                       blkPhiCoarse,
                                       clusE,
                                       iEbin,
                                       towerEs.size() );



    // QA: fill pT vs leadTower, fill cluster-level eta/phi
    h_pt_eta->Fill(clusPt, lt_eta);
    h_pt_eta_rw->Fill(clusPt, lt_eta, weight);
    h_etaphi_clus->Fill(clusEta, clusPhi);
      
    match1 = false;
    TLorentzVector ph1_trEtaPhi(0,0,0,0);
    if (m_isSimulation)
    {
      processSimulationTruthMatches( clus1,
                                       photon1,
                                       lt_eta,
                                       lt_phi,
                                       blkCoord,
                                       blkEtaCoarse,
                                       blkPhiCoarse,
                                       towerPhis,
                                       towerEs,
                                       vtx_z,
                                       truth_photons,
                                       truth_meson_photons,
                                       match1,
                                       ph1_trEtaPhi );
    }

    if (Verbosity() > 2)
    {
      std::cout << "[DEBUG] => Starting INNER loop for cluster pairs.\n";
    }

    /* (ii) π0 pair building --------------------------------------- */
    if (analyseThisEvent)
        processClusterPairs( clusterContainer,
                             cIt1,
                             vertex,
                             photon1,
                             clusE,
                             lt_eta,
                             lt_phi,
                             blkEtaCoarse,
                             blkPhiCoarse,
                             blkCoord,
                             maxAlpha,
                             ptMaxCut,
                             pt2ClusCut,
                             pi0ptcut,
                             weight,
                             match1,
                             ph1_trEtaPhi,
                             truth_meson_photons,
                             m_isSimulation );

  } // end cluster1 loop

  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD
              << "[DEBUG] finalClusterLoop: COMPLETED cluster loops.\n"
              << " => # of clusters in event = " << nClusCount << "\n"
              << " => weight=" << weight << "\n"
              << "finalClusterLoop() EXIT"
              << ANSI_RESET << "\n\n";
  }
}


// ============================================================================
// Helper: compute the circular‑mean φ–offset between the ideal tower grid
//         and the actual detector geometry.  Executed once, the first time
//         we have a valid `m_geometry`.
//
// Output:
//   • m_phi0Offset   –  rigid shift  (range −π … +π)
//   • m_hasOffset    –  set to true on success
//
// Return value:
//   • true  – offset successfully measured
//   • false – geometry container empty  ➜ caller may skip the event
// ============================================================================
bool PositionDependentCorrection::computeRigidPhiOffset()
{
  /* ------------------------------------------------------------------ */
  /* 1) Geometry summary (printed once – useful QA)                     */
  /* ------------------------------------------------------------------ */
  const int   nEtaBins   = m_geometry->get_etabins();   // ~96
  const int   nPhiBins   = m_geometry->get_phibins();   // 128 in sPHENIX
  const float radPerBin  = 2.F * M_PI / static_cast<float>(nPhiBins);

  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD << "[φ‑anchor] Geometry summary\n" << ANSI_RESET
              << "      nEtaBins  = " << nEtaBins  << '\n'
              << "      nPhiBins  = " << nPhiBins  << '\n'
              << "      Δφ(bin)   = " << radPerBin << " rad\n";
  }

  /* ------------------------------------------------------------------ */
  /* 2) Circular mean of φ‑residuals over *all* towers                  */
  /* ------------------------------------------------------------------ */
  double      sumSin  = 0.0, sumCos = 0.0;
  std::size_t nTowers = 0;

  const auto [beg, end] = m_geometry->get_tower_geometries();
  for (auto it = beg; it != end; ++it)
  {
    const RawTowerGeom* tg = it->second;
    if (!tg) continue;

    const int   iphi     = tg->get_binphi();
    const float phiIdeal = (iphi + 0.5F) * radPerBin;
    const float delta    = TVector2::Phi_mpi_pi( tg->get_phi() - phiIdeal );

    sumSin += std::sin(delta);
    sumCos += std::cos(delta);
    ++nTowers;

    if (Verbosity() > 1 && nTowers <= 3)
      std::cout << "      tower#" << nTowers
                << "  iphi="    << iphi
                << "  φ(real)=" << tg->get_phi()
                << "  φ(ideal)="<< phiIdeal
                << "  Δφ="      << delta << '\n';
  }

  if (nTowers == 0)
  {
    std::cerr << " [φ‑anchor] WARNING: geometry container is empty – skip event\n";
    return false;                      // helper **failed**
  }

  /* ------------------------------------------------------------------ */
  /* 3) Mean offset  (range −π … +π], stored for later corrections)     */
  /* ------------------------------------------------------------------ */
  m_phi0Offset = std::atan2(sumSin, sumCos);   // circular mean
  m_hasOffset  = true;

  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD << ANSI_CYAN
              << "[φ‑anchor] measured rigid barrel tilt  Δφ₀ = "
              << m_phi0Offset << "  rad"
              << ANSI_RESET  << "\n"
              << "      towers used = " << nTowers << '\n'
              << "      ⟨sinΔφ⟩ = "    << sumSin / nTowers
              << " ,  ⟨cosΔφ⟩ = "     << sumCos / nTowers << '\n';
  }
  return true;                         // helper **succeeded**
}

bool PositionDependentCorrection::computeRigidEtaOffset()
{
    constexpr float kEtaMin  = -1.1f;
    constexpr float kDEtaBin = 2.2f / 96.f;   // full span / bins

    double sum = 0.0;   std::size_t n = 0;

    const auto [beg,end] = m_geometry->get_tower_geometries();
    for (auto it=beg; it!=end; ++it)
    {
        auto *tg = it->second;  if (!tg) continue;
        int   ieta  = tg->get_bineta();
        float ideal = kEtaMin + (ieta + 0.5f)*kDEtaBin;
        sum += tg->get_eta() - ideal;
        ++n;
    }
    if (!n) return false;

    m_eta0Offset   = float(sum / n);
    m_hasEtaOffset = true;

    if (Verbosity()>0)
        std::cout << "[η‑anchor]  rigid offset  Δη₀ = "
                  << m_eta0Offset << '\n';
    return true;
}



// ============================================================================
int PositionDependentCorrection::process_towers(PHCompositeNode* topNode)
{
    if (Verbosity() > 0)
    {
        std::cout << ANSI_BOLD << ANSI_CYAN
        << "[process_towers] START" << ANSI_RESET
        << "  => event counter: " << _eventcounter
        << std::endl;
    }
    if ((_eventcounter % 1000) == 0)
    {
        std::cout << ANSI_BOLD << ANSI_YELLOW
        << "[Info] Processing event " << _eventcounter
        << ANSI_RESET << std::endl;
    }
    
    float emcal_hit_threshold = 0.5;  // GeV
    
    if (debug)
    {
        std::cout << ANSI_BOLD << "[DEBUG] " << ANSI_RESET
        << "-----------------------------------" << std::endl;
    }
    float maxAlpha       = 0.6;
    float clus_chisq_cut = 4.0;
    float nClus_ptCut    = 0.5;
    int   max_nClusCount = 3000000;
    
    // ------------------------------------------------------------------------
    // 1) Retrieve vertex z-position
    // ------------------------------------------------------------------------
    float vtx_z = retrieveVertexZ(topNode);
    float truth_vz = 0.0;
    
    if (Verbosity() > 0)
    {
        std::cout << ANSI_BOLD
        << "[process_towers] Retrieved vertex Z position: "
        << ANSI_RESET << vtx_z << std::endl;
    }
    /* ================================================================ *
     * Two‑tier vertex‑Z acceptance
     * ================================================================ */
    const float absVz = std::fabs(vtx_z);
    /* (B) outside the configured slice table → ignore event completely */
    if (absVz > m_vzSliceMax)
    {
        if (Verbosity() > 1)
            std::cout << "[PDC]  |z_vtx|=" << absVz
                      << " cm > " << m_vzSliceMax
                      << " cm – skip event\n";
        return Fun4AllReturnCodes::EVENT_OK;
    }
    /* (C) tight physics window                                         */
    const bool passTightVz = (absVz <= m_vzTightCut);

  // ------------------------------------------------------------------------
  // 2) Tower info
  // ------------------------------------------------------------------------
  std::vector<float> ht_eta;
  std::vector<float> ht_phi;
  float tower_tot_e = 0;
  fillTowerInfo(topNode, emcal_hit_threshold, tower_tot_e, ht_eta, ht_phi);

  if (Verbosity() > 5)
  {
    std::cout << ANSI_BOLD
              << "[process_towers] Finished fillTowerInfo, total tower energy="
              << ANSI_RESET << tower_tot_e
              << ", #tower etas=" << ht_eta.size()
              << ", #tower phis=" << ht_phi.size() << std::endl;
  }
  float weight = 1;

  // ------------------------------------------------------------------------
  // 3) Retrieve truth info
  // ------------------------------------------------------------------------
  std::vector<TLorentzVector> truth_photons;
  std::vector<TLorentzVector> truth_meson_photons;
    
  if (m_isSimulation)
  {
      fillTruthInfo(topNode, vtx_z, truth_photons, truth_meson_photons);
      h2_truthReco_vz->Fill(truth_vz, vtx_z);
      
      if (Verbosity() > 0)
      {
          std::cout << "[process_towers] => fillTruthInfo done.\n"
          << "    #truth_photons=" << truth_photons.size()
          << ", #truth_meson_photons=" << truth_meson_photons.size()
          << std::endl;
      }
  }
    
    m_geometry = checkTowerGeometry(topNode);
    
    if (Verbosity() > 0 && m_geometry)
    {
        std::cout << ANSI_BOLD
                  << "[process_towers] => Tower geometry retrieved successfully."
                  << ANSI_RESET << '\n';
    }

    // ──────────────────────────────────────────────────────────────────────
    // One‑time measurement of the rigid φ–offset (barrel tilt)
    // ──────────────────────────────────────────────────────────────────────
    if (!m_hasOffset && m_geometry)
        if (!computeRigidPhiOffset())          // returns *false* ⇢ geometry empty
         return Fun4AllReturnCodes::EVENT_OK;
      
      
    if (!m_hasEtaOffset && m_geometry)
        computeRigidEtaOffset();     // non‑fatal if it fails
    
    
    /* --------------------------------------------------------------- *
     * Loose but *not* tight → fill only Δη(E,vz) spectra              *
     * --------------------------------------------------------------- */
    if (!passTightVz)
    {
        RawClusterContainer* cc = retrieveClusterContainer(topNode);
        if (cc && m_geometry)
        {
            RawClusterContainer::ConstRange rng = cc->getClusters();
            for (auto it = rng.first; it != rng.second; ++it)
            {
                RawCluster* cl = it->second;     if (!cl) continue;
                CLHEP::Hep3Vector vtx(0,0,vtx_z);
                CLHEP::Hep3Vector evec = RawClusterUtility::GetEVec(*cl, vtx);
                TLorentzVector    reco;  reco.SetPtEtaPhiE(evec.perp(),
                                                           evec.pseudoRapidity(),
                                                           evec.phi(),
                                                           evec.mag());
                std::vector<int>   tEta, tPhi;   std::vector<float> tE;
                for (auto t = cl->get_towers().first;
                     t != cl->get_towers().second; ++t)
                {
                    const auto* tg = m_geometry->get_tower_geometry(t->first);
                    if (!tg) continue;
                    tEta.push_back(tg->get_bineta());
                    tPhi.push_back(tg->get_binphi());
                    tE  .push_back(t->second);
                }
                PDC::Geo::LocalCoord blk =
                        PDC::Geo::computeLocal(tEta, tPhi, tE);
                
                int blkEta = blk.blkEta;
                int blkPhi = blk.blkPhi;
                /* dummy truth photon – only the residual is used            */
                TLorentzVector truth;
                /* Δη – five‑way */
                fillDEtaAllVariants(cl,
                                    reco,
                                    truth,
                                    {blk.locEta, blk.locPhi},   // (ηloc , φloc)
                                    blkEta,
                                    blkPhi,
                                    vtx_z,
                                    h_eta_diff_cpRaw_E,
                                    h_eta_diff_cpCorr_E,
                                    /*fillGlobal=*/false );
                
            }
        }
     return Fun4AllReturnCodes::EVENT_OK;      // done with loose event
  }

  RawClusterContainer* clusterContainer = retrieveClusterContainer(topNode);
  if (!clusterContainer)
  {
    if (Verbosity() > 0)
    {
      std::cout << ANSI_BOLD << ANSI_YELLOW
                << "[process_towers] => Cluster container missing => returning 0."
                << ANSI_RESET << std::endl
                << "[process_towers] END\n" << std::endl;
    }
    return 0;
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cout << ANSI_BOLD
                << "[process_towers] => Cluster container retrieved OK."
                << ANSI_RESET << std::endl;
    }
  }

  int nClusContainer = 0;
  int nClusCount = countClusters(clusterContainer,
                                 vtx_z,
                                 nClus_ptCut,
                                 clus_chisq_cut,
                                 nClusContainer);

  if (Verbosity() > 0)
  {
    std::cout << "[process_towers] => countClusters => "
              << ANSI_BOLD << "nClusCount=" << nClusCount << ANSI_RESET
              << ", nClusContainer=" << nClusContainer
              << std::endl;
  }

  float ptMaxCut   = 100;
  float pt1ClusCut = 0.9;  // 1.3
  float pt2ClusCut = 0.9;  // 0.7
  float pi0ptcut   = 0.f; // was 0 by default

 if (passTightVz)
 {
     finalClusterLoop(topNode,
                      clusterContainer,
                      vtx_z,
                      truth_photons,
                      truth_meson_photons,
                      tower_tot_e,
                      max_nClusCount,
                      nClusCount,
                      maxAlpha,
                      ptMaxCut,
                      pt1ClusCut,
                      pt2ClusCut,
                      pi0ptcut,
                      weight);
 }

  // Clear vectors for tower bin indices
  ht_phi.clear();
  ht_eta.clear();

  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD << "[process_towers] END" << ANSI_RESET
              << "\n" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

// ============================================================================
int PositionDependentCorrection::End(PHCompositeNode* /*topNode*/)
{
  /* ───────────────────────────────────────── File handling ─────────── */
  if (Verbosity() > 0)
    std::cout << "[DEBUG] PositionDependentCorrection::End() – entering\n";

  if (outfile)
  {
    outfile->cd();
    if (Verbosity() > 0)
      std::cout << "[DEBUG] Writing histograms to " << outfilename << " …\n";

    outfile->Write();

    /* optional inventory printout ------------------------------------- */
    if (Verbosity() > 0)
    {
      std::cout << "\n[INFO] Histograms written – content overview\n"
                << std::left << std::setw(30) << "Histogram"
                << std::setw(14)              << "Type"
                << std::setw(12)              << "#Entries\n"
                << "────────────────────────────────────────────────────────\n";

      TIter nextKey(gDirectory->GetListOfKeys());
      while (TKey* key = static_cast<TKey*>(nextKey()))
        if (TH1* h = dynamic_cast<TH1*>(key->ReadObj()))
          std::cout << std::left << std::setw(30) << key->GetName()
                    << std::setw(14)             << key->GetClassName()
                    << std::right<< std::setw(12)<< (Long64_t)h->GetEntries()
                    << '\n';

      std::cout << "────────────────────────────────────────────────────────\n";
    }

    outfile->Close();
    delete outfile;
    outfile = nullptr;
    if (Verbosity() > 0) std::cout << "[DEBUG] Output file closed.\n";
  }
  else if (Verbosity() > 0)
    std::cerr << "[ERROR] outfile pointer is null – nothing written!\n";

  /* ──────────────────────────────────────── helper: % formatting ───── */
  auto pct = [](std::uint64_t n, std::uint64_t d)->std::string
             { return d ? Form("%6.1f %%", 100.*double(n)/double(d))
                         : "   n/a "; };

  /* ──────────────────────────────────────── 1) Δφ five-way table ───── */
  const std::uint64_t nTotPhi =
        m_phiWinCLUSraw + m_phiWinCLUScp + m_phiWinCLUSbcorr +
        m_phiWinPDCraw  + m_phiWinPDCcorr;

  if (Verbosity() > 0)
  {
    std::cout << '\n' << ANSI_BOLD
              << "╭───────────────────────────────╮\n"
              << "│        Δφ   –   FIVE-WAY      │\n"
              << "╰───────────────────────────────╯" << ANSI_RESET << '\n'
              << std::left << std::setw(14) << "variant"
              << std::right<< std::setw(12)<< "wins"
              << std::setw(12)             << "share\n"
              << "────────────────────────────────────────────\n"
              << std::left << std::setw(14) << "CLUS-RAW"
              << std::right<< std::setw(12)<< m_phiWinCLUSraw
              << std::setw(12)             << pct(m_phiWinCLUSraw,nTotPhi) << '\n'
              << std::left << std::setw(14) << "CLUS-CP"
              << std::right<< std::setw(12)<< m_phiWinCLUScp
              << std::setw(12)             << pct(m_phiWinCLUScp ,nTotPhi) << '\n'
              << std::left << std::setw(14) << "CLUS-BCORR"
              << std::right<< std::setw(12)<< m_phiWinCLUSbcorr
              << std::setw(12)             << pct(m_phiWinCLUSbcorr,nTotPhi)<< '\n'
              << std::left << std::setw(14) << "PDC-RAW"
              << std::right<< std::setw(12)<< m_phiWinPDCraw
              << std::setw(12)             << pct(m_phiWinPDCraw ,nTotPhi) << '\n'
              << std::left << std::setw(14) << "PDC-CORR"
              << std::right<< std::setw(12)<< m_phiWinPDCcorr
              << std::setw(12)             << pct(m_phiWinPDCcorr,nTotPhi)<< '\n'
              << "────────────────────────────────────────────\n"
              << std::left << std::setw(14) << "TOTAL"
              << std::right<< std::setw(12)<< nTotPhi << "\n";

    /* quick verdict */
    if (m_phiWinCLUScp)
    {
      const double rel = 100.*double(m_phiWinCLUSbcorr)/double(m_phiWinCLUScp);
      std::cout << "BCORR / CP (φ) : " << std::fixed << std::setprecision(1)
                << rel << " %  →  "
                << (rel >= 100. ? "✅ matches/beats CP" : "⚠️ below CP")
                << '\n';
    }
  }

  /* ──────────────────────────────────────── 2) Δη five-way table ───── */
  const std::uint64_t nTotEta5 =
        m_etaWinCLUSraw + m_etaWinCLUScp + m_etaWinCLUSbcorr +
        m_etaWinPDCraw  + m_etaWinPDCcorr;

  if (Verbosity() > 0)
  {
    std::cout << '\n' << ANSI_BOLD
              << "╭───────────────────────────────╮\n"
              << "│        Δη   –   FIVE-WAY      │\n"
              << "╰───────────────────────────────╯" << ANSI_RESET << '\n'
              << std::left << std::setw(14) << "variant"
              << std::right<< std::setw(12)<< "wins"
              << std::setw(12)             << "share\n"
              << "────────────────────────────────────────────\n"
              << std::left << std::setw(14) << "CLUS-RAW"
              << std::right<< std::setw(12)<< m_etaWinCLUSraw
              << std::setw(12)             << pct(m_etaWinCLUSraw,nTotEta5)<< '\n'
              << std::left << std::setw(14) << "CLUS-CP"
              << std::right<< std::setw(12)<< m_etaWinCLUScp
              << std::setw(12)             << pct(m_etaWinCLUScp ,nTotEta5)<< '\n'
              << std::left << std::setw(14) << "CLUS-BCORR"
              << std::right<< std::setw(12)<< m_etaWinCLUSbcorr
              << std::setw(12)             << pct(m_etaWinCLUSbcorr,nTotEta5)<< '\n'
              << std::left << std::setw(14) << "PDC-RAW"
              << std::right<< std::setw(12)<< m_etaWinPDCraw
              << std::setw(12)             << pct(m_etaWinPDCraw ,nTotEta5)<< '\n'
              << std::left << std::setw(14) << "PDC-CORR"
              << std::right<< std::setw(12)<< m_etaWinPDCcorr
              << std::setw(12)             << pct(m_etaWinPDCcorr,nTotEta5)<< '\n'
              << "────────────────────────────────────────────\n"
              << std::left << std::setw(14) << "TOTAL"
              << std::right<< std::setw(12)<< nTotEta5 << "\n";

    if (m_etaWinCLUScp)
    {
      const double rel = 100.*double(m_etaWinCLUSbcorr)/double(m_etaWinCLUScp);
      std::cout << "BCORR / CP (η) : " << std::fixed << std::setprecision(1)
                << rel << " %  →  "
                << (rel >= 100. ? "✅ matches/beats CP" : "⚠️ below CP")
                << '\n';
    }

  }
    
    /* ─────────────────────── block‑coordinate occupancy ─────────────── */
    if (Verbosity() > 5)
    {
        static const char* kSlot[4] = {"(0,0)","(0,1)","(1,0)","(1,1)"};
        std::uint64_t nTot = 0;
        for (const auto& v : m_blkLocCount) nTot += v.load();

        std::cout << '\n' << ANSI_BOLD
                  << "╭───────────────────────────────╮\n"
                  << "│  Local‑block occupancy (η,φ) │\n"
                  << "╰───────────────────────────────╯" << ANSI_RESET << '\n'
                  << std::left << std::setw(8)  << "slot"
                  << std::right<< std::setw(14) << "counts"
                  << std::setw(12)              << "share\n"
                  << "────────────────────────────────────────────\n";

        for (int i = 0; i < 4; ++i)
        {
            const auto v = m_blkLocCount[i].load();
            std::cout << std::left << std::setw(8)  << kSlot[i]
                      << std::right<< std::setw(14) << v
                      << std::setw(11) << (nTot ? Form("%6.2f %%",
                             100.0 * double(v) / double(nTot)) : "  n/a")
                      << '\n';
        }
        std::cout << "────────────────────────────────────────────\n"
                  << std::left << std::setw(8)  << "TOTAL"
                  << std::right<< std::setw(14) << nTot << "\n";
    }


  return 0;
}
