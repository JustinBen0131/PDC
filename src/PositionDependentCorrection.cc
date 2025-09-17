#include "PositionDependentCorrection.h"
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
#include <TSystem.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawCluster.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawClusterv2.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawClusterContainer.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawClusterUtility.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawTower.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawTowerContainer.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawTowerGeom.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawTowerGeomContainer.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/TowerInfo.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/TowerInfoContainer.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/TowerInfoDefs.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawCluster.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawClusterContainer.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawTowerDefs.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/BEmcRecCEMC.h"

using namespace PDC_detail;

constexpr const char* ANSI_BOLD  = "\033[1m";
constexpr const char* ANSI_RESET = "\033[0m";
constexpr const char* ANSI_GREEN = "\033[32m";
constexpr const char* ANSI_YELLOW= "\033[33m";
constexpr const char* ANSI_CYAN    = "\033[36m";
constexpr const char* ANSI_RED   = "\033[31m";

R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)

using namespace std;
namespace CLHEP { class Hep3Vector; }

PositionDependentCorrection::PositionDependentCorrection(const std::string &name,
                                                         const std::string &filename)
  : SubsysReco(name)
  , m_isSimulation(false)          // << initialise first
  , detector("HCALIN")
  , outfilename(filename)
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


/** Book all histograms that are needed for _both_ data and simulation */
void PositionDependentCorrection::bookCommonHistograms
(const std::function<std::string(int)>& makeLabel)
{
  // ------------------------------------------------------------------
  // Common binning used by the base 3D block-coordinate histograms
  // ------------------------------------------------------------------
  Double_t xEdges[15], yEdges[15];
  {
    const double s = 2.0/14.0;
    for (int i = 0; i <= 14; ++i)
    {
      xEdges[i] = -0.5 + i*s;
      yEdges[i] = -0.5 + i*s;
    }
  }
  static constexpr Double_t eEdges[9] = {2,4,6,8,10,12,15,20,30};
  const char* modeTag = (m_binningMode == EBinningMode::kRange ? "range" : "disc");

  // ==================================================================
  // (A) ALWAYS BOOK the 5 uncorrected 3-D histograms (both modes)
  // ==================================================================
  if (!h3_cluster_block_cord_E)
  {
    h3_cluster_block_cord_E = new TH3F(
        Form("h3_blockCoord_E_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Uncorrected local block coords vs. E_{slice}" :
           "Uncorrected local block coords vs. E_{centre}"),
        14, xEdges, 14, yEdges, 8, eEdges);
  }

  if (!h3_cluster_block_cord_E_full)
  {
    h3_cluster_block_cord_E_full = new TH3F(
        Form("h3_blockCoord_E_full_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Uncorrected local block coords vs. E_{slice} (full |#eta| #leq 1.10)" :
           "Uncorrected local block coords vs. E_{centre} (full |#eta| #leq 1.10)"),
        14, xEdges, 14, yEdges, 8, eEdges);
  }

  if (!h3_cluster_block_cord_E_etaCore)
  {
    h3_cluster_block_cord_E_etaCore = new TH3F(
        Form("h3_blockCoord_E_etaCore_le0p20_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Uncorrected local block coords vs. E_{slice} (|#eta| #leq 0.20)" :
           "Uncorrected local block coords vs. E_{centre} (|#eta| #leq 0.20)"),
        14, xEdges, 14, yEdges, 8, eEdges);
  }

  if (!h3_cluster_block_cord_E_etaMid)
  {
    h3_cluster_block_cord_E_etaMid = new TH3F(
        Form("h3_blockCoord_E_etaMid_0p20to0p70_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Uncorrected local block coords vs. E_{slice} (0.20 < |#eta| #leq 0.70)" :
           "Uncorrected local block coords vs. E_{centre} (0.20 < |#eta| #leq 0.70)"),
        14, xEdges, 14, yEdges, 8, eEdges);
  }

  if (!h3_cluster_block_cord_E_etaEdge)
  {
    h3_cluster_block_cord_E_etaEdge = new TH3F(
        Form("h3_blockCoord_E_etaEdge_0p70to1p10_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Uncorrected local block coords vs. E_{slice} (0.70 < |#eta| #leq 1.10)" :
           "Uncorrected local block coords vs. E_{centre} (0.70 < |#eta| #leq 1.10)"),
        14, xEdges, 14, yEdges, 8, eEdges);
  }

  // QA: Δ between property raw CoG and recomputed raw CoG vs E
  if (!h_dx_prop_vsE)
  {
    h_dx_prop_vsE = new TH2F("h_dx_prop_vsE",
                             "x_{raw}(prop) - x_{raw}(recalc CoG) vs E;E [GeV];#Delta x_{raw} [tower]",
                             80,0,40, 200,-1e-2,1e-2);
  }
  if (!h_dy_prop_vsE)
  {
    h_dy_prop_vsE = new TH2F("h_dy_prop_vsE",
                             "y_{raw}(prop) - y_{raw}(recalc CoG) vs E;E [GeV];#Delta y_{raw} [tower]",
                             80,0,40, 200,-1e-2,1e-2);
  }
  if (hm)
  {
    hm->registerHisto(h_dx_prop_vsE);
    hm->registerHisto(h_dy_prop_vsE);
  }

  // Minimal QA needed by fillTowerInfo/fillTruthInfo/process_towers
  if (!h_tower_e)     { h_tower_e     = new TH1F("h_tower_e","",1000,-1,5); }
  if (!h_emcal_e_eta) { h_emcal_e_eta = new TH1F("h_emcal_e_eta","",96,0,96); }
  if (!h_cemc_etaphi) { h_cemc_etaphi = new TH2F("h_cemc_etaphi","",96,0,96,256,0,256); }

  if (!h_truth_eta) { h_truth_eta = new TH1F("h_truth_eta","",100,-1.2,1.2); }
  if (!h_truth_e)   { h_truth_e   = new TH1F("h_truth_e","",100,0,10); }
  if (!h_truth_vz)  { h_truth_vz  = new TH1F("h_truth_vz","Truth Vertex Z;z_{truth} (cm);Counts",200,-100,100); }
  if (!h2_truthReco_vz)
  {
    h2_truthReco_vz = new TH2F("h2_truthReco_vz",
                               "Truth vs Reco Vertex Z;z_{truth} (cm);z_{reco} (cm)",
                               200,-100,100,200,-100,100);
  }
  if (!h_vert_xy) { h_vert_xy = new TH2F("h_vert_xy","",500,-120,120,500,-120,120); }
  if (!h_reco_vz) { h_reco_vz = new TH1F("h_reco_vz","Reco Vertex Z;z_{reco} (cm);Counts",200,-100,100); }

  // Register only the QA that needs manager
  if (hm)
  {
    hm->registerHisto(h_reco_vz);
    hm->registerHisto(h_truth_vz);
    hm->registerHisto(h2_truthReco_vz);
  }

  // ==================================================================
  // (C) FULL MODE: book everything else (without rebooking base-5)
  // ==================================================================

  // 1) “Always-on” histograms (original list)
  if (!h_mass_eta_lt)    { h_mass_eta_lt    = new TH2F("h_mass_eta_lt","",50,0,0.5,96,0,96); }
  if (!h_mass_eta_lt_rw) { h_mass_eta_lt_rw = new TH2F("h_mass_eta_lt_rw","",50,0,0.5,96,0,96); }
  if (!h_pt_eta)         { h_pt_eta         = new TH2F("h_pt_eta","",100,0,10,96,0,96); }
  if (!h_pt_eta_rw)      { h_pt_eta_rw      = new TH2F("h_pt_eta_rw","",100,0,10,96,0,96); }
  if (!h_cemc_etaphi)    { h_cemc_etaphi    = new TH2F("h_cemc_etaphi","",96,0,96,256,0,256); }
  if (!h_InvMass)        { h_InvMass        = new TH1F("h_InvMass","Invariant Mass",500,0,1.0); }
  if (!h_InvMass_w)      { h_InvMass_w      = new TH1F("h_InvMass_w","Invariant Mass",500,0,1.0); }
  if (!h_InvMassMix)     { h_InvMassMix     = new TH1F("h_InvMassMix","Invariant Mass",120,0,1.2); }

  if (!h_mE_raw)
  {
    h_mE_raw = new TH2F("h_mE_raw",
                        "π^{0} mass vs E_{#gamma}^{max} – RAW;M_{γγ} [GeV];E_{#gamma}^{max} [GeV]",
                        120,0.00,0.30, N_Ebins, eEdge);
    if (hm) { hm->registerHisto(h_mE_raw); }
  }
  if (!h_mE_corr)
  {
    h_mE_corr = new TH2F("h_mE_corr",
                         "π^{0} mass vs E_{#gamma}^{max} – φ-corr.;M_{γγ} [GeV];E_{#gamma}^{max} [GeV]",
                         120,0.00,0.30, N_Ebins, eEdge);
    if (hm) { hm->registerHisto(h_mE_corr); }
  }

  if (!h_m_blk_raw)
  {
    h_m_blk_raw = new TH3F("h_m_blk_raw",
                           "π^{0} mass – RAW;η_{loc};φ_{loc};M_{γγ} [GeV]",
                           14,-0.5,1.5, 14,-0.5,1.5, 60,0.00,0.30);
    if (hm) { hm->registerHisto(h_m_blk_raw); }
  }
  if (!h_m_blk_corr)
  {
    h_m_blk_corr = new TH3F("h_m_blk_corr",
                            "π^{0} mass – φ-corr.;η_{loc};φ_{loc};M_{γγ} [GeV]",
                            14,-0.5,1.5, 14,-0.5,1.5, 60,0.00,0.30);
    if (hm) { hm->registerHisto(h_m_blk_corr); }
  }

  // Full-mode scalar/2D/3D histos (guard against double-booking)
  if (!h_tower_e)     { h_tower_e     = new TH1F("h_tower_e","",1000,-1,5); }
  if (!h_etaphi_clus) { h_etaphi_clus = new TH2F("h_etaphi_clus","",140,-1.2,1.2,64,-TMath::Pi(),TMath::Pi()); }
  if (!h_clusE)       { h_clusE       = new TH1F("h_clusE","",100,0,10); }
  if (!h_clusE_nTow)  { h_clusE_nTow  = new TH2F("h_clusE_nTow","",20,0,20,50,0,50); }
  if (!h_emcal_e_eta) { h_emcal_e_eta = new TH1F("h_emcal_e_eta","",96,0,96); }
  if (!h_pt1)         { h_pt1         = new TH1F("h_pt1","",100,0,5); }
  if (!h_pt2)         { h_pt2         = new TH1F("h_pt2","",100,0,5); }
  if (!h_nclusters)   { h_nclusters   = new TH1F("h_nclusters","",100,0,100); }
  if (!h_matched_res) { h_matched_res = new TH2F("h_matched_res","",100,0,1.5,20,-1,1); }
  if (!h_res_e)       { h_res_e       = new TH2F("h_res_e","",100,0,1.5,20,0,20); }
  if (!h_res_e_phi)   { h_res_e_phi   = new TH3F("h_res_e_phi","",100,0,1.5,10,0,20,256,0,256); }
  if (!h_res_e_eta)   { h_res_e_eta   = new TH3F("h_res_e_eta","",300,0,1.5,40,0,20,96,0,96); }
  if (!h_m_pt_eta)    { h_m_pt_eta    = new TH3F("h_m_pt_eta","",70,0,0.7,10,0,10,96,0,96); }
  if (!h_m_ptTr_eta)  { h_m_ptTr_eta  = new TH3F("h_m_ptTr_eta","",70,0,0.7,10,0,10,96,0,96); }
  if (!h_m_ptTr_eta_trKin) { h_m_ptTr_eta_trKin = new TH3F("h_m_ptTr_eta_trKin","",70,0,0.7,10,0,10,96,0,96); }
  if (!h_res)         { h_res         = new TH1F("h_res","",50,0,1.5); }
  if (!h_delEta_e_eta){ h_delEta_e_eta= new TH3F("h_delEta_e_eta","",100,-0.1,0.1,10,0,20,96,0,96); }
  if (!h_delR_e_eta)  { h_delR_e_eta  = new TH3F("h_delR_e_eta","",100,-0.1,0.1,10,0,20,96,0,96); }
  if (!h_delPhi_e_eta){ h_delPhi_e_eta= new TH3F("h_delPhi_e_eta","",100,-0.3,0.3,20,0,20,96,0,96); }
  if (!h_delPhi_e_phi){ h_delPhi_e_phi= new TH3F("h_delPhi_e_phi","",100,-0.1,0.1,20,0,20,256,0,256); }
  if (!pr_eta_shower) { pr_eta_shower = new TProfile("pr_eta_shower","",96,-48.5,47.5,-1,1.5); }
  if (!pr_phi_shower) { pr_phi_shower = new TProfile("pr_phi_shower","",256,-128.5,127.5,-1,1.5); }
  if (!h_vert_xy)     { h_vert_xy     = new TH2F("h_vert_xy","",500,-120,120,500,-120,120); }

  // χ² QA maps
  if (!h2_chi2_tot_etaPhi)
  {
    h2_chi2_tot_etaPhi = new TH2F("h2_chi2_tot_etaPhi",
                                  "Clusters BEFORE #chi^{2} cut;Tower #eta index;Tower #varphi index;Counts",
                                  96,0,96, 256,0,256);
    if (hm) { hm->registerHisto(h2_chi2_tot_etaPhi); }
  }
  if (!h2_chi2_rej_etaPhi)
  {
    h2_chi2_rej_etaPhi = new TH2F("h2_chi2_rej_etaPhi",
                                  "Clusters REJECTED by #chi^{2} cut;Tower #eta index;Tower #varphi index;Counts",
                                  96,0,96, 256,0,256);
    if (hm) { hm->registerHisto(h2_chi2_rej_etaPhi); }
  }
  if (!p_chi2_pass_etaPhi)
  {
    p_chi2_pass_etaPhi = new TProfile2D("p_chi2_pass_etaPhi",
                                        "Pass fraction after #chi^{2} cut;Tower #eta index;Tower #varphi index;⟨pass⟩",
                                        96,0,96, 256,0,256, -0.1,1.1);
    if (hm) { hm->registerHisto(p_chi2_pass_etaPhi); }
  }

  // Per-block histograms
  for (int ie = 0; ie < NBinsBlock; ++ie)
  {
    for (int ip = 0; ip < NBinsBlock; ++ip)
    {
      if (!h_mass_block_pt[ie][ip])
      {
        h_mass_block_pt[ie][ip] = new TH2F(Form("h_mass_block_%d_%d_pt",ie,ip), "", 100,0,1, 5,0,10);
      }
      if (!h_res_block_E [ie][ip])
      {
        h_res_block_E [ie][ip] = new TH2F(Form("h_res_block_%d_%d_E" ,ie,ip), "", 120,0,1.2, 5,0,10);
      }
    }
  }

  // ==================================================================
  // (D) Corrected TH3 histos (legacy + per-η views) — full mode only
  //     (Do NOT re-create base-5 uncorrected here)
  // ==================================================================
  if (!h3_cluster_block_cord_E_corrected)
  {
    h3_cluster_block_cord_E_corrected = new TH3F(
        Form("h3_blockCoord_Ecorr_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Corrected local block coords vs. E_{slice}" :
           "Corrected local block coords vs. E_{centre}"),
        14,xEdges,14,yEdges,8,eEdges);
  }

  if (!h3_cluster_block_cord_E_full_corr)
  {
    h3_cluster_block_cord_E_full_corr = new TH3F(
        Form("h3_blockCoord_Ecorr_full_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Corrected local block coords vs. E_{slice} (full |#eta| #leq 1.10)" :
           "Corrected local block coords vs. E_{centre} (full |#eta| #leq 1.10)"),
        14,xEdges,14,yEdges,8,eEdges);
  }

  if (!h3_cluster_block_cord_E_etaCore_corr)
  {
    h3_cluster_block_cord_E_etaCore_corr = new TH3F(
        Form("h3_blockCoord_Ecorr_etaCore_le0p20_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Corrected local block coords vs. E_{slice} (|#eta| #leq 0.20)" :
           "Corrected local block coords vs. E_{centre} (|#eta| #leq 0.20)"),
        14,xEdges,14,yEdges,8,eEdges);
  }

  if (!h3_cluster_block_cord_E_etaMid_corr)
  {
    h3_cluster_block_cord_E_etaMid_corr = new TH3F(
        Form("h3_blockCoord_Ecorr_etaMid_0p20to0p70_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Corrected local block coords vs. E_{slice} (0.20 < |#eta| #leq 0.70)" :
           "Corrected local block coords vs. E_{centre} (0.20 < |#eta| #leq 0.70)"),
        14,xEdges,14,yEdges,8,eEdges);
  }

  if (!h3_cluster_block_cord_E_etaEdge_corr)
  {
    h3_cluster_block_cord_E_etaEdge_corr = new TH3F(
        Form("h3_blockCoord_Ecorr_etaEdge_0p70to1p10_%s",modeTag),
        (m_binningMode==EBinningMode::kRange ?
           "Corrected local block coords vs. E_{slice} (0.70 < |#eta| #leq 1.10)" :
           "Corrected local block coords vs. E_{centre} (0.70 < |#eta| #leq 1.10)"),
        14,xEdges,14,yEdges,8,eEdges);
  }

  // ==================================================================
  // (E) Per-slice corrected local η/φ (unchanged logic; guarded)
  // ==================================================================
  for (int i=0;i<N_Ebins;++i)
  {
    const float eLo = eEdge[i], eHi = eEdge[i+1];
    const std::string lab = makeLabel(i);

    if (!h_localPhi_corrected[i])
    {
      h_localPhi_corrected[i] = new TH1F(
            Form("h_localPhi_corr_%s",lab.c_str()),
            (m_binningMode==EBinningMode::kRange ?
               Form("Corrected local #phi;%.0f < E < %.0f GeV",eLo,eHi) :
               Form("Corrected local #phi;E = %.0f GeV",eLo)),
            50,-0.5,0.5);
    }

    if (!h_localEta_corrected[i])
    {
      h_localEta_corrected[i] = new TH1F(
            Form("h_localEta_corr_%s",lab.c_str()),
            (m_binningMode==EBinningMode::kRange ?
               Form("Corrected local #eta;%.0f < E < %.0f GeV",eLo,eHi) :
               Form("Corrected local #eta;E = %.0f GeV",eLo)),
            50,-0.5,0.5);
    }
  }
}



/** Extra histograms that are **only** needed when m_isSimulation == true */
void
PositionDependentCorrection::bookSimulationHistograms
(const std::function<std::string(int)>& makeLabel)
{
    
  // Skip ALL simulation histogram booking in first-pass b-values-only mode
  if (m_firstPassBvaluesOnly) return;
  /* ---------------------- truth & matching ------------------------- */
  h_delR_recTrth   = new TH1F("h_delR_recTrth","",500,0,5);

  h_truthE         = new TH1F("h_truthE","",10000,0,30);
    
    // Book per-variant π0 mass histograms (8 variants × N_Ebins)
    {
      const char* vtag[] = {
        "CLUSRAW", "CLUSCP",
        "EAgeom", "EAetaE", "EAEonly", "EAmix",
        "PDCraw", "PDCcorr"
      };
      for (int v = 0; v < static_cast<int>(VarPi0::NVAR); ++v) {
        for (int i = 0; i < N_Ebins; ++i) {
          h_m_pi0_var[v][i] = new TH1F(
            Form("h_m_pi0_%s_%g_%g", vtag[v], eEdge[i], eEdge[i+1]),
            Form("#pi^{0} mass (%s);M_{#gamma#gamma} [GeV];Counts", vtag[v]),
            120, 0.00, 0.30
          );
          if (hm) hm->registerHisto(h_m_pi0_var[v][i]);
        }
      }
    }

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

      const char* title2D =
        (m_binningMode==EBinningMode::kRange)
          ? Form("#Delta#phi vs #eta;#eta;#Delta#phi  [rad]  (%.0f<E<%.0f GeV)", eLo, eHi)
          : Form("#Delta#phi vs #eta;#eta;#Delta#phi  [rad]  (E=%.0f GeV)",       eLo);

      h2_phi_diff_vsEta_RAW_E[i]     = new TH2F(
          Form("h2_dphi_vs_eta_RAW_%s",    lab.c_str()), title2D,
          96, -1.1, 1.1,   200, -0.1, 0.1);
      h2_phi_diff_vsEta_CP_E[i]      = new TH2F(
          Form("h2_dphi_vs_eta_CP_%s",     lab.c_str()), title2D,
          96, -1.1, 1.1,   200, -0.1, 0.1);
      h2_phi_diff_vsEta_BCORR_E[i]   = new TH2F(
          Form("h2_dphi_vs_eta_BCORR_%s",  lab.c_str()), title2D,
          96, -1.1, 1.1,   200, -0.1, 0.1);
      h2_phi_diff_vsEta_PDCraw_E[i]  = new TH2F(
          Form("h2_dphi_vs_eta_PDCraw_%s", lab.c_str()), title2D,
          96, -1.1, 1.1,   200, -0.1, 0.1);
      h2_phi_diff_vsEta_PDCcorr_E[i] = new TH2F(
          Form("h2_dphi_vs_eta_PDCCorr_%s",lab.c_str()), title2D,
          96, -1.1, 1.1,   200, -0.1, 0.1);

      hm->registerHisto(h2_phi_diff_vsEta_RAW_E[i]);
      hm->registerHisto(h2_phi_diff_vsEta_CP_E[i]);
      hm->registerHisto(h2_phi_diff_vsEta_BCORR_E[i]);
      hm->registerHisto(h2_phi_diff_vsEta_PDCraw_E[i]);
      hm->registerHisto(h2_phi_diff_vsEta_PDCcorr_E[i]);
      
      /* -------- CP helper histograms (needed by code below) ------------- */
      h_phi_diff_cpRaw_E [i] = new TH1F(
            Form("h_phi_diff_cpRaw_%s",  lab.c_str()),
            "#Delta#phi RAW‑CP;#Delta#phi (rad);Counts",   200, -0.1, 0.1);
      h_phi_diff_cpCorr_E[i] = new TH1F(
            Form("h_phi_diff_cpCorr_%s", lab.c_str()),
            "#Delta#phi CP‑corr;#Delta#phi (rad);Counts",  200, -0.1, 0.1);
      // --- EA φ: four variants, uniquely named
      h_phi_diff_cpCorrEA_geom_E[i] = new TH1F(
            Form("h_phi_diff_cpCorrEA_geom_%s", lab.c_str()),
            "#Delta#phi CP-corr(EA geom);#Delta#phi (rad);Counts", 200, -0.1, 0.1);

      h_phi_diff_cpCorrEA_fitEtaDep_E[i] = new TH1F(
            Form("h_phi_diff_cpCorrEA_fitEtaDep_%s", lab.c_str()),
            "#Delta#phi CP-corr(EA |#eta|+E fits);#Delta#phi (rad);Counts", 200, -0.1, 0.1);

      h_phi_diff_cpCorrEA_fitEnergyOnly_E[i] = new TH1F(
            Form("h_phi_diff_cpCorrEA_fitEnergyOnly_%s", lab.c_str()),
            "#Delta#phi CP-corr(EA E-only fits);#Delta#phi (rad);Counts", 200, -0.1, 0.1);

      h_phi_diff_cpCorrEA_fitPhiEnergy_etaEtaDep_E[i] = new TH1F(
            Form("h_phi_diff_cpCorrEA_fitPhiE_etaEtaDep_%s", lab.c_str()),
            "#Delta#phi CP-corr(EA: #varphi E-only, #eta |#eta|+E);#Delta#phi (rad);Counts", 200, -0.1, 0.1);

      h_phi_diff_cpBcorr_E[i] = new TH1F(
            Form("h_phi_diff_cpBcorr_%s", lab.c_str()),
            "#Delta#phi b‑corr;#Delta#phi (rad);Counts",   200, -0.1, 0.1);

      h_eta_diff_cpRaw_E [i] = new TH1F(
            Form("h_eta_diff_cpRaw_%s",  lab.c_str()),
            "#Delta#eta RAW‑CP;#Delta#eta;Counts",         200, -0.1, 0.1);
      h_eta_diff_cpCorr_E[i] = new TH1F(
            Form("h_eta_diff_cpCorr_%s", lab.c_str()),
            "#Delta#eta CP‑corr;#Delta#eta;Counts",        200, -0.1, 0.1);
      // --- EA η: four variants, uniquely named
      h_eta_diff_cpCorrEA_geom_E[i] = new TH1F(
            Form("h_eta_diff_cpCorrEA_geom_%s", lab.c_str()),
            "#Delta#eta CP-corr(EA geom);#Delta#eta;Counts", 200, -0.1, 0.1);

      h_eta_diff_cpCorrEA_fitEtaDep_E[i] = new TH1F(
            Form("h_eta_diff_cpCorrEA_fitEtaDep_%s", lab.c_str()),
            "#Delta#eta CP-corr(EA |#eta|+E fits);#Delta#eta;Counts", 200, -0.1, 0.1);

      h_eta_diff_cpCorrEA_fitEnergyOnly_E[i] = new TH1F(
            Form("h_eta_diff_cpCorrEA_fitEnergyOnly_%s", lab.c_str()),
            "#Delta#eta CP-corr(EA E-only fits);#Delta#eta;Counts", 200, -0.1, 0.1);

      h_eta_diff_cpCorrEA_fitPhiEnergy_etaEtaDep_E[i] = new TH1F(
            Form("h_eta_diff_cpCorrEA_fitPhiE_etaEtaDep_%s", lab.c_str()),
            "#Delta#eta CP-corr(EA: #varphi E-only, #eta |#eta|+E);#Delta#eta;Counts", 200, -0.1, 0.1);

      h_eta_diff_cpBcorr_E[i] = new TH1F(
            Form("h_eta_diff_cpBcorr_%s", lab.c_str()),
            "#Delta#eta b‑corr;#Delta#eta;Counts",         200, -0.1, 0.1);

      hm->registerHisto(h_phi_diff_cpRaw_E [i]);
      hm->registerHisto(h_phi_diff_cpCorr_E[i]);
      hm->registerHisto(h_phi_diff_cpCorrEA_geom_E[i]);
      hm->registerHisto(h_phi_diff_cpCorrEA_fitEtaDep_E[i]);
      hm->registerHisto(h_phi_diff_cpCorrEA_fitEnergyOnly_E[i]);
      hm->registerHisto(h_phi_diff_cpCorrEA_fitPhiEnergy_etaEtaDep_E[i]);
      hm->registerHisto(h_phi_diff_cpBcorr_E[i]);
      hm->registerHisto(h_eta_diff_cpRaw_E [i]);
      hm->registerHisto(h_eta_diff_cpCorr_E[i]);
      hm->registerHisto(h_eta_diff_cpCorrEA_geom_E[i]);
      hm->registerHisto(h_eta_diff_cpCorrEA_fitEtaDep_E[i]);
      hm->registerHisto(h_eta_diff_cpCorrEA_fitEnergyOnly_E[i]);
      hm->registerHisto(h_eta_diff_cpCorrEA_fitPhiEnergy_etaEtaDep_E[i]);
      hm->registerHisto(h_eta_diff_cpBcorr_E[i]);
      
      /* ------------------------------------------------------------------ *
       * (IIa)  NEW  Δφ(E, vz) histograms – one loop over |z| slices
       * ------------------------------------------------------------------ */
      for (int iVz = 0; iVz < N_VzBins; ++iVz)
        {
          const float vzLo = vzEdge[iVz];
          const float vzHi = vzEdge[iVz+1];
          const std::string vzTag = Form("vz%.0f_%.0f", vzLo, vzHi);

          h_phi_diff_raw_E_vz[i][iVz] = new TH1F(
                Form("h_phi_diff_raw_%s_%s",  lab.c_str(), vzTag.c_str()),
                Form("#Delta#phi raw;%.0f<E<%.0f GeV, %.0f<|z_{vtx}|<%.0f cm",
                     eLo, eHi, vzLo, vzHi), 200, -0.1, 0.1);

          h_phi_diff_corrected_E_vz[i][iVz] = new TH1F(
                Form("h_phi_diff_corr_%s_%s", lab.c_str(), vzTag.c_str()),
                Form("#Delta#phi corr;%.0f<E<%.0f GeV, %.0f<|z_{vtx}|<%.0f cm",
                     eLo, eHi, vzLo, vzHi), 200, -0.1, 0.1);

          /* CP variants */
          h_phi_diff_cpRaw_E_vz [i][iVz] = new TH1F(
                Form("h_phi_diff_cpRaw_%s_%s",  lab.c_str(), vzTag.c_str()),
                "#Delta#phi RAW‑CP;#Delta#phi(rad);Counts", 200,-0.1,0.1);
          h_phi_diff_cpCorr_E_vz[i][iVz] = new TH1F(
                Form("h_phi_diff_cpCorr_%s_%s", lab.c_str(), vzTag.c_str()),
                "#Delta#phi CP‑corr;#Delta#phi(rad);Counts",200,-0.1,0.1);
          h_phi_diff_cpBcorr_E_vz[i][iVz] = new TH1F(
                Form("h_phi_diff_cpBcorr_%s_%s", lab.c_str(), vzTag.c_str()),
                "#Delta#phi BCORR;#Delta#phi(rad);Counts", 200,-0.1,0.1);

          hm->registerHisto(h_phi_diff_raw_E_vz       [i][iVz]);
          hm->registerHisto(h_phi_diff_corrected_E_vz [i][iVz]);
          hm->registerHisto(h_phi_diff_cpRaw_E_vz     [i][iVz]);
          hm->registerHisto(h_phi_diff_cpCorr_E_vz    [i][iVz]);
          hm->registerHisto(h_phi_diff_cpBcorr_E_vz   [i][iVz]);
        }

        /* ---------- OPTIONAL signed‑z booking (north / south) ---------- */
        if (m_useSignedVz)
        {
          for (int iVz = 0; iVz < N_VzBinsSigned; ++iVz)
          {
            const bool  isNorth = (iVz < N_VzBins);
            const int   idxAbs  = isNorth ? iVz : iVz - N_VzBins;
            const float vzLo = vzEdge[idxAbs];
            const float vzHi = vzEdge[idxAbs+1];
            const std::string vzTag = Form("%s_vz%s%.0f_%.0f",
                                           lab.c_str(), isNorth?"P":"N",
                                           vzLo, vzHi);

            h_phi_diff_raw_E_vzsgn[i][iVz] = new TH1F(
                  Form("h_phi_diff_raw_%s", vzTag.c_str()),
                  Form("#Delta#phi raw;%.0f<E<%.0f GeV, %c %.0f<z<%.0f cm",
                       eLo, eHi, isNorth?'+':'-', vzLo, vzHi),
                  200,-0.1,0.1);

            h_phi_diff_corrected_E_vzsgn[i][iVz] = new TH1F(
                  Form("h_phi_diff_corr_%s", vzTag.c_str()),
                  Form("#Delta#phi corr;%.0f<E<%.0f GeV, %c %.0f<z<%.0f cm",
                       eLo, eHi, isNorth?'+':'-', vzLo, vzHi),
                  200,-0.1,0.1);

            h_phi_diff_cpRaw_E_vzsgn [i][iVz] = new TH1F(
                  Form("h_phi_diff_cpRaw_%s",  vzTag.c_str()),
                  "#Delta#phi RAW‑CP;#Delta#phi;Counts", 200,-0.1,0.1);
            h_phi_diff_cpCorr_E_vzsgn[i][iVz] = new TH1F(
                  Form("h_phi_diff_cpCorr_%s", vzTag.c_str()),
                  "#Delta#phi CP‑corr;#Delta#phi;Counts",200,-0.1,0.1);
            h_phi_diff_cpBcorr_E_vzsgn[i][iVz] = new TH1F(
                  Form("h_phi_diff_cpBcorr_%s", vzTag.c_str()),
                  "#Delta#phi BCORR;#Delta#phi;Counts", 200,-0.1,0.1);

            hm->registerHisto(h_phi_diff_raw_E_vzsgn       [i][iVz]);
            hm->registerHisto(h_phi_diff_corrected_E_vzsgn [i][iVz]);
            hm->registerHisto(h_phi_diff_cpRaw_E_vzsgn     [i][iVz]);
            hm->registerHisto(h_phi_diff_cpCorr_E_vzsgn    [i][iVz]);
            hm->registerHisto(h_phi_diff_cpBcorr_E_vzsgn   [i][iVz]);
          }
      }
        
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

      for (int ieEA = 0; ieEA < N_Ebins; ++ieEA)
      {
        const std::string labEA = sliceTag(ieEA);

        if (!h2_phi_diffEA_vs_bphi_E[ieEA])
        {
          h2_phi_diffEA_vs_bphi_E[ieEA] = new TH2F(
            Form("h2_dphiEA_vs_bphi_%s", labEA.c_str()),
            (m_binningMode == EBinningMode::kRange ?
               Form("#Delta#varphi(CLUS-CP(EA)) vs b_{#varphi};b_{#varphi};#Delta#varphi (folded) [rad]  (%.0f<E<%.0f GeV)",
                    eEdge[ieEA], eEdge[ieEA+1]) :
               Form("#Delta#varphi(CLUS-CP(EA)) vs b_{#varphi};b_{#varphi};#Delta#varphi (folded) [rad]  (E=%.0f GeV)",
                    eEdge[ieEA])),
            60, 0.12, 0.25,
            120, -0.030, 0.030);
          if (hm) hm->registerHisto(h2_phi_diffEA_vs_bphi_E[ieEA]);
        }

        if (!h2_eta_diffEA_vs_beta_E[ieEA])
        {
          h2_eta_diffEA_vs_beta_E[ieEA] = new TH2F(
            Form("h2_detaEA_vs_beta_%s", labEA.c_str()),
            (m_binningMode == EBinningMode::kRange ?
               Form("#Delta#eta(CLUS-CP(EA)) vs b_{#eta};b_{#eta};#Delta#eta  (%.0f<E<%.0f GeV)",
                    eEdge[ieEA], eEdge[ieEA+1]) :
               Form("#Delta#eta(CLUS-CP(EA)) vs b_{#eta};b_{#eta};#Delta#eta  (E=%.0f GeV)",
                    eEdge[ieEA])),
            80, 0.10, 0.32,
            120, -0.060, 0.060);
          if (hm) hm->registerHisto(h2_eta_diffEA_vs_beta_E[ieEA]);
        }
      }


      
    /* ---------- OPTIONAL signed‑z booking (north / south) ---------- */
    if (m_useSignedVz)
    {
          for (int iVz = 0; iVz < N_VzBinsSigned; ++iVz)
          {
            const bool  isNorth = (iVz < N_VzBins);
            const int   idxAbs  =  isNorth ?  iVz : iVz - N_VzBins;
            const float vzLo    = vzEdge[idxAbs];
            const float vzHi    = vzEdge[idxAbs+1];
            const std::string vzTag = Form("%s_vz%s%.0f_%.0f",
                                           lab.c_str(),
                                           isNorth ? "P" : "N",
                                           vzLo, vzHi);  // “P” = +z, “N” = –z

            /* book once per slice */
            h_eta_diff_raw_E_vzsgn[i][iVz] = new TH1F(
                   Form("h_eta_diff_raw_%s", vzTag.c_str()),
                   Form("#Delta#eta raw;%.0f<E<%.0f GeV, %c %.0f<z<%.0f cm",
                        eLo, eHi, isNorth?'+':'-', vzLo, vzHi),
                   200, -0.1, 0.1);

            h_eta_diff_corrected_E_vzsgn[i][iVz] = new TH1F(
                   Form("h_eta_diff_corr_%s", vzTag.c_str()),
                   Form("#Delta#eta corr;%.0f<E<%.0f GeV, %c %.0f<z<%.0f cm",
                        eLo, eHi, isNorth?'+':'-', vzLo, vzHi),
                   200, -0.1, 0.1);

            h_eta_diff_cpRaw_E_vzsgn  [i][iVz] = new TH1F(
                   Form("h_eta_diff_cpRaw_%s",  vzTag.c_str()),
                   "#Delta#eta RAW‑CP;#Delta#eta;Counts", 200, -0.1, 0.1);
            h_eta_diff_cpCorr_E_vzsgn [i][iVz] = new TH1F(
                   Form("h_eta_diff_cpCorr_%s", vzTag.c_str()),
                   "#Delta#eta CP‑corr;#Delta#eta;Counts", 200, -0.1, 0.1);
            h_eta_diff_cpBcorr_E_vzsgn[i][iVz] = new TH1F(
                   Form("h_eta_diff_cpBcorr_%s", vzTag.c_str()),
                   "#Delta#eta BCORR;#Delta#eta;Counts", 200, -0.1, 0.1);

            hm->registerHisto(h_eta_diff_raw_E_vzsgn       [i][iVz]);
            hm->registerHisto(h_eta_diff_corrected_E_vzsgn [i][iVz]);
            hm->registerHisto(h_eta_diff_cpRaw_E_vzsgn     [i][iVz]);
            hm->registerHisto(h_eta_diff_cpCorr_E_vzsgn    [i][iVz]);
            hm->registerHisto(h_eta_diff_cpBcorr_E_vzsgn   [i][iVz]);
        } // signed‑vz loop
     }   // if(m_useSignedVz)
  }     // E‑slice loop


    // ────────────────────────────────────────────────────────────────
    // Ash/Log scans (per‑E TProfiles).  Also seed default scan grids
    // if the macro didn't pre-fill m_bScan / m_w0Scan.
    // ────────────────────────────────────────────────────────────────
    auto makeEdges = [](const std::vector<double>& grid) {
      std::vector<double> edges;
      if (grid.size() < 2) return edges;
      const int n = static_cast<int>(grid.size());
      edges.resize(n+1);
      for (int i=1; i<n; ++i) edges[i] = 0.5*(grid[i-1] + grid[i]);
      edges[0] = grid[0] - 0.5*(grid[1] - grid[0]);
      edges[n] = grid[n-1] + 0.5*(grid[n-1] - grid[n-2]);
      return edges;
    };

    // 3a) Sane **defaults** (coarse pass) if not set by the macro
    if (m_bScan.empty()) {
      // Ash b: tower‑space parameter used by BEmcRec; cover common 0.10–0.30
      // Use a moderate step to keep Condor jobs light; we refine later around minima.
      for (double b = 0.10; b <= 0.30 + 1e-12; b += 0.02)
        m_bScan.push_back(std::round(b*10000.)/10000.);
    }
    if (m_w0Scan.empty()) {
      // Log w0: per PHENIX/ATLAS conventions, the optimum is typically 3–5.
      for (double w0 = 2.5; w0 <= 6.5 + 1e-12; w0 += 0.25)
        m_w0Scan.push_back(std::round(w0*100.)/100.);
    }

    const std::vector<double> bEdges  = makeEdges(m_bScan);
    const std::vector<double> w0Edges = makeEdges(m_w0Scan);

    for (int i=0; i<N_Ebins; ++i)
    {
      const std::string lab = makeLabel(i);

      // <(Δφ)^2> vs b  (Ash)
      if (!p_phi_rms2_vs_b_E[i] && !bEdges.empty()) {
        p_phi_rms2_vs_b_E[i] = new TProfile(
          Form("p_phi_rms2_vs_b_%s",lab.c_str()),
          "#LT(#Delta#phi)^{2}#GT vs b; b; #LT(#Delta#phi)^{2}#GT",
          static_cast<int>(bEdges.size()-1), bEdges.data());
        if (hm) hm->registerHisto(p_phi_rms2_vs_b_E[i]);
      }
      // <(Δη)^2> vs b  (Ash)
      if (!p_eta_rms2_vs_b_E[i] && !bEdges.empty()) {
        p_eta_rms2_vs_b_E[i] = new TProfile(
          Form("p_eta_rms2_vs_b_%s",lab.c_str()),
          "#LT(#Delta#eta)^{2}#GT vs b; b; #LT(#Delta#eta)^{2}#GT",
          static_cast<int>(bEdges.size()-1), bEdges.data());
        if (hm) hm->registerHisto(p_eta_rms2_vs_b_E[i]);
      }

      // <(Δφ)^2> vs w0 (Log)
      if (!p_phi_rms2_vs_w0_E[i] && !w0Edges.empty()) {
        p_phi_rms2_vs_w0_E[i] = new TProfile(
          Form("p_phi_rms2_vs_w0_%s",lab.c_str()),
          "#LT(#Delta#phi)^{2}#GT vs w_{0}; w_{0}; #LT(#Delta#phi)^{2}#GT",
          static_cast<int>(w0Edges.size()-1), w0Edges.data());
        if (hm) hm->registerHisto(p_phi_rms2_vs_w0_E[i]);
      }
      // <(Δη)^2> vs w0 (Log)
      if (!p_eta_rms2_vs_w0_E[i] && !w0Edges.empty()) {
        p_eta_rms2_vs_w0_E[i] = new TProfile(
          Form("p_eta_rms2_vs_w0_%s",lab.c_str()),
          "#LT(#Delta#eta)^{2}#GT vs w_{0}; w_{0}; #LT(#Delta#eta)^{2}#GT",
          static_cast<int>(w0Edges.size()-1), w0Edges.data());
        if (hm) hm->registerHisto(p_eta_rms2_vs_w0_E[i]);
      }

      // Agreement diagnostic in tower units: | CP(b) – CP(EA) |
      if (!p_abs_dloc_phi_CPea_vs_b_E[i] && !bEdges.empty()) {
        p_abs_dloc_phi_CPea_vs_b_E[i] = new TProfile(
          Form("p_abs_dloc_phi_CPea_vs_b_%s",lab.c_str()),
          "|x_{CP}(b)-x_{EA}| vs b; b; #LT|#Delta x|#GT  (tower units)",
          static_cast<int>(bEdges.size()-1), bEdges.data(), "S");
        if (hm) hm->registerHisto(p_abs_dloc_phi_CPea_vs_b_E[i]);
      }
      if (!p_abs_dloc_eta_CPea_vs_b_E[i] && !bEdges.empty()) {
        p_abs_dloc_eta_CPea_vs_b_E[i] = new TProfile(
          Form("p_abs_dloc_eta_CPea_vs_b_%s",lab.c_str()),
          "|y_{CP}(b)-y_{EA}| vs b; b; #LT|#Delta y|#GT  (tower units)",
          static_cast<int>(bEdges.size()-1), bEdges.data(), "S");
        if (hm) hm->registerHisto(p_abs_dloc_eta_CPea_vs_b_E[i]);
      }
    }
}


bool PositionDependentCorrection::loadBValues(const std::string& bFilePath)
{
  // If both corrections are disabled up-front, there is nothing to read
  if (!isFitDoneForPhi && !isFitDoneForEta)
  {
    std::cout << "[INFO] isFitDoneForPhi and isFitDoneForEta are both false ⇒ skip reading bValues.txt\n";
    return true;
  }

  std::ifstream bfile(bFilePath);
  if (!bfile.is_open())
  {
    std::cout << "[WARN]  bValues.txt NOT found at " << bFilePath
              << "  ➜  disabling φ- and η-corrections.\n";
    isFitDoneForPhi = isFitDoneForEta = false;
    std::fill(m_bPhiReady_view.begin(), m_bPhiReady_view.end(), false);
    std::fill(m_bEtaReady_view.begin(), m_bEtaReady_view.end(), false);
    return false;
  }

  std::cout << "[INFO]  Reading b-values from " << bFilePath << '\n';

  // Reset per-view readiness and values
  for (int v = 0; v < kNEtaViews; ++v)
  {
    m_bPhiReady_view[v] = false;
    m_bEtaReady_view[v] = false;
    for (int i = 0; i < N_Ebins; ++i)
    {
      m_bValsPhi_view[v][i] = 0.f;
      m_bValsEta_view[v][i] = 0.f;
    }
  }

  // Track which slices we successfully filled
  std::vector<bool> gotLegacyPhi(N_Ebins, false), gotLegacyEta(N_Ebins, false);
  std::array<std::array<bool, N_Ebins>, kNEtaViews> gotPhi{}; // default-initialized false
  std::array<std::array<bool, N_Ebins>, kNEtaViews> gotEta{};

  // Map a label → η-view index
  auto viewIndexOf = [](std::string lab) -> int
  {
    if (lab.empty()) return -1;
    // normalize
    std::transform(lab.begin(), lab.end(), lab.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });

    // legacy/default (apply to m_bVals*): originalEta / default / noeta / no-eta
    if (lab == "originaleta" || lab == "default" || lab == "noeta" || lab == "no-eta")
      return -1;

    // view 0: full |eta| <= 1.10
    if (lab == "fulleta" || lab == "full" || lab == "full-eta" || lab == "full_eta")
      return 0;

    // view 1: |eta| <= 0.20
    if (lab == "etacore" || lab == "core" || lab == "eta-core" || lab == "eta_core")
      return 1;

    // view 2: 0.20 < |eta| <= 0.70
    if (lab == "etamid" || lab == "mid" || lab == "eta-mid" || lab == "eta_mid")
      return 2;

    // view 3: 0.70 < |eta| <= 1.10
    if (lab == "etaedge" || lab == "edge" || lab == "eta-edge" || lab == "eta_edge")
      return 3;

    // unknown → treat as legacy/default
    return -1;
  };

  // Accept lines like:
  //   PHI [2.0,4.0)  0.188029  originalEta
  //   ETA [2.0,4.0)  0.194202  fullEta
  // Both 3-column (legacy) and 4-column (with label) are supported.
  const std::regex lineRe(
      R"(^\s*(PHI|ETA)\s*\[\s*([0-9]*\.?[0-9]+)\s*,\s*([0-9]*\.?[0-9]+)\s*\)\s*([0-9]*\.?[0-9]+)\s*(\S+)?\s*$)");

  std::string line;
  while (std::getline(bfile, line))
  {
    if (line.empty() || line[0] == '#') continue;

    std::smatch m;
    if (!std::regex_match(line, m, lineRe)) continue;

    const std::string dim = m[1];   // "PHI" or "ETA"
    const float eLo = std::stof(m[2]);
    const float eHi = std::stof(m[3]);
    const float bVal = std::stof(m[4]);
    const std::string label = m[5].matched ? m[5].str() : std::string{};
    const int vIdx = viewIndexOf(label);

    // find target energy slice
    int ie = -1;
    for (int i = 0; i < N_Ebins; ++i)
    {
      if (std::fabs(eLo - expectedLo[i]) < 1e-3f &&
          std::fabs(eHi - expectedHi[i]) < 1e-3f)
      { ie = i; break; }
    }
    if (ie < 0) continue; // unknown E-range

    if (dim == "PHI")
    {
      if (vIdx >= 0) { m_bValsPhi_view[vIdx][ie] = bVal; gotPhi[vIdx][ie] = true; }
      else           { m_bValsPhi[ie]            = bVal; gotLegacyPhi[ie]  = true; }
    }
    else // "ETA"
    {
      if (vIdx >= 0) { m_bValsEta_view[vIdx][ie] = bVal; gotEta[vIdx][ie] = true; }
      else           { m_bValsEta[ie]            = bVal; gotLegacyEta[ie]  = true; }
    }
  }
  bfile.close();

  // Per-view readiness: only "OK" if we have all slices for that view
  for (int v = 0; v < kNEtaViews; ++v)
  {
    m_bPhiReady_view[v] = std::all_of(gotPhi[v].begin(), gotPhi[v].end(),
                                      [](bool b){ return b; });
    m_bEtaReady_view[v] = std::all_of(gotEta[v].begin(), gotEta[v].end(),
                                      [](bool b){ return b; });
  }

  // Legacy/default readiness (for the always-present corrected histograms)
  const bool legacyPhiReady = std::all_of(gotLegacyPhi.begin(), gotLegacyPhi.end(),
                                          [](bool b){ return b; });
  const bool legacyEtaReady = std::all_of(gotLegacyEta.begin(), gotLegacyEta.end(),
                                          [](bool b){ return b; });

  // Keep legacy semantics: only apply the *default* correction if the default table is complete
  isFitDoneForPhi = isFitDoneForPhi && legacyPhiReady;
  isFitDoneForEta = isFitDoneForEta && legacyEtaReady;

  // ----- optional human-readable summary -----------------------------------
  auto printTable = [&](const char* tag, const std::vector<bool>& got, const float* arr, bool enabled)
  {
    std::cout << "[INFO]  " << tag << (enabled ? "  (enabled)\n" : "  (DISABLED — missing slices)\n");
    for (int i = 0; i < N_Ebins; ++i)
    {
      std::cout << "        [" << expectedLo[i] << ',' << expectedHi[i] << ")  : ";
      if (got[i]) std::cout << arr[i] << '\n';
      else        std::cout << "-- missing --\n";
    }
  };
  printTable("DEFAULT PHI b-values", gotLegacyPhi, m_bValsPhi, isFitDoneForPhi);
  printTable("DEFAULT ETA b-values", gotLegacyEta, m_bValsEta, isFitDoneForEta);

  auto printView = [&](const char* name, int v)
  {
    std::cout << "[INFO]  view=" << name
              << "  PHI(" << (m_bPhiReady_view[v] ? "OK" : "MISS") << ")"
              << "  ETA(" << (m_bEtaReady_view[v] ? "OK" : "MISS") << ")\n";
    for (int i = 0; i < N_Ebins; ++i)
    {
      std::cout << "        [" << expectedLo[i] << ',' << expectedHi[i] << ")  : "
                << "b_phi=" << (gotPhi[v][i] ? m_bValsPhi_view[v][i] : NAN)
                << " | b_eta=" << (gotEta[v][i] ? m_bValsEta_view[v][i] : NAN) << '\n';
    }
  };
  printView("fullEta", 0);
  printView("etaCore", 1);
  printView("etaMid",  2);
  printView("etaEdge", 3);

  // Return true if we can apply *any* corrections (legacy or per-view)
  const bool anyPhi = isFitDoneForPhi || m_bPhiReady_view[0] || m_bPhiReady_view[1]
                                       || m_bPhiReady_view[2] || m_bPhiReady_view[3];
  const bool anyEta = isFitDoneForEta || m_bEtaReady_view[0] || m_bEtaReady_view[1]
                                       || m_bEtaReady_view[2] || m_bEtaReady_view[3];
  return (anyPhi || anyEta);
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

int PositionDependentCorrection::Init(PHCompositeNode* /*topNode*/)
{
  /* 0) Basic sanity / I/O managers ----------------------------------- */
  if (Verbosity() > 0) {
      std::cout << "[DEBUG] PositionDependentCorrection::Init() called.\n";
  }


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


    if (Verbosity()>0) {
      std::cout << "[DEBUG] PositionDependentCorrection::Init() completed successfully.\n";
    }
    
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
    // Reset per-event truth π0 photon cache
    m_truth_pi0_photons.clear();

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
    if (h_truth_eta) h_truth_eta->Fill(myVector.Eta());
    if (h_truth_e)   h_truth_e->Fill(energy, wieght);

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

            TLorentzVector myVector;
            myVector.SetXYZT(truth->get_px(), truth->get_py(),
                             truth->get_pz(), truth->get_e());

            truth_meson_photons.push_back(myVector);

            // NEW: if it’s from a π0, also keep mother id/pid
            if (parent->get_pid() == 111)
            {
              TruthPhoton tp;
              tp.p4 = myVector;
              tp.mother_id  = parent->get_track_id();
              tp.mother_pid = parent->get_pid(); // 111
              m_truth_pi0_photons.push_back(tp);
            }

            if (debug)
            {
              std::cout << ANSI_BOLD << "[Debug]" << ANSI_RESET
                        << " 2nd photon (meson)  e=" << myVector.E()
                        << "  phi=" << myVector.Phi()
                        << "  eta=" << myVector.Eta()
                        << "  parent=" << parent->get_pid() << std::endl;
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

  std::string towergeomnodename = "TOWERGEOM_CEMC_DETAILED";
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
    std::cout << "  --> Attempting to retrieve 'CLUSTERINFO_CEMC_V2' from the node tree." << std::endl;
  }

  RawClusterContainer* clusterContainer =
      findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC_V2");

  if (!clusterContainer)
  {
    std::cerr << "[FATAL] CLUSTERINFO_CEMC_V2 node is missing.\n";
    std::cerr << "        Ensure the builder calls WriteClusterV2(true)." << std::endl;
    gSystem->Exit(1);
    return nullptr;
  }

  if (Verbosity() > 0)
  {
    std::cout << "  --> Successfully retrieved 'CLUSTERINFO_CEMC_V2' container." << std::endl;
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  φ at shower depth
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float PositionDependentCorrection::phiAtShowerDepth( float  energy,
                                                     double rFront,
                                                     double zFront,
                                                     float  phiFront,
                                                     int    ix,        // NEW
                                                     int    iy         // NEW
                                                   ) const
{
  const double xA = rFront * std::cos(phiFront);
  const double yA = rFront * std::sin(phiFront);

  float xSD, ySD, zSD;
  m_bemcRec->CorrectShowerDepth(ix, iy,          // <─ NEW
                                energy,
                                xA, yA, zFront,
                                xSD, ySD, zSD);

  return std::atan2(ySD, xSD);                   // (−π … +π]
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  η at shower depth  (vertex‑aware)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float PositionDependentCorrection::etaAtShowerDepth( float   energy,
                                                     double  rFront,
                                                     double  zFront,
                                                     float   phiFront,
                                                     int     ix,          // fine‑φ
                                                     int     iy,          // fine‑η
                                                     float   vtx_z )const
{
  /* front‑face Cartesian */
  const double xA = rFront * std::cos(phiFront);
  const double yA = rFront * std::sin(phiFront);

  /* propagate to shower depth */
  float xSD,ySD,zSD;
  m_bemcRec->CorrectShowerDepth(ix,iy,
                                energy,
                                xA,yA,zFront,
                                xSD,ySD,zSD);

  /* shift to event vertex and convert to η */
  const double zRel = zSD - static_cast<double>(vtx_z);
  const double R    = std::hypot(xSD,ySD);
  const double theta= std::atan2(R, zRel);                  // 0 … π
  return -std::log( std::tan( 0.5*theta ) );                // pseudorapidity
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  x at shower depth (cm)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double PositionDependentCorrection::xAtShowerDepth( float  energy,
                                                    double rFront,
                                                    double zFront,
                                                    float  phiFront,
                                                    int    ix,         // NEW
                                                    int    iy          // NEW
                                                  ) const
{
  const double xA = rFront * std::cos(phiFront);
  const double yA = rFront * std::sin(phiFront);

  float xSD, ySD, zSD;
  m_bemcRec->CorrectShowerDepth(ix, iy,          // <─ NEW
                                energy,
                                xA, yA, zFront,
                                xSD, ySD, zSD);

  return static_cast<double>(xSD);
}


/// Inverse‑asinh correction in *tower units* about the nearest tower index.
/// Input/Output: xTower = φ- or η‑like cluster CoG coordinate in tower units.
/// Mirrors the φ/η branch inside BEmcRecCEMC::CorrectPosition(...).
float PositionDependentCorrection::doAshShift(float xTower, float b)
{
  if (!std::isfinite(xTower) || b <= 0.f) return xTower;

  // nearest tower index (same rounding as the clusterizer)
  const int   i0 = EmcCluster::lowint(xTower + 0.5f);
  const float du = xTower - i0;

  // numeric guard: only act inside the owner tower
  if (std::fabs(du) > 0.5f) return xTower;

  const float S = std::sinh(0.5f / b);     // sinh(1/2b)
  return i0 + b * std::asinh(2.f * du * S);
}


float
PositionDependentCorrection::doLogWeightCoord( const std::vector<int>&   towerPhi,
                                               const std::vector<float>& towerE,
                                               float                     w0 )
{
  if (towerPhi.size() != towerE.size() || towerPhi.empty())
    return std::numeric_limits<float>::quiet_NaN();

  const double sumE = std::accumulate(towerE.begin(), towerE.end(), 0.0);
  if (sumE <= 0.0) return std::numeric_limits<float>::quiet_NaN();

  // reference: φ index of max‑E tower (handle wrap on the cylinder)
  const std::size_t iRef =
      std::distance(towerE.begin(), std::max_element(towerE.begin(), towerE.end()));
  const int refPhi = towerPhi[iRef];

  const int nPhi = (m_bemcRec ? m_bemcRec->GetNx() : 256);

  double sw = 0.0, sphi = 0.0;
  for (std::size_t i = 0; i < towerE.size(); ++i)
  {
    const double Ei = towerE[i];
    if (Ei <= 0.0) continue;

    // standard logarithmic weighting: wi = max(0, w0 + ln(Ei / Etot))
    double wi = w0 + std::log(Ei / sumE);
    if (wi < 0.0) wi = 0.0;

    int phi = towerPhi[i];
    if (phi - refPhi < -nPhi/2) phi += nPhi;
    if (phi - refPhi >  nPhi/2) phi -= nPhi;

    sw   += wi;
    sphi += wi * static_cast<double>(phi);
  }
  if (sw <= 0.0) return std::numeric_limits<float>::quiet_NaN();

  float x = static_cast<float>(sphi / sw);
  while (x < -0.5f)              x += nPhi;
  while (x >= nPhi - 0.5f)       x -= nPhi;
  return x;
}



namespace {
  // Prefer RawClusterv2 getters; otherwise fall back to Momenta()
  inline bool getTowerCG_fromPropsOrMomenta(const RawCluster* clus,
                                            RawTowerGeomContainer* geo,
                                            BEmcRecCEMC* rec,
                                            float& xCG, float& yCG)
  {
    xCG = yCG = std::numeric_limits<float>::quiet_NaN();
    if (!clus || !geo || !rec) return false;

    // 1) Try RawClusterv2 tower-space getters
    if (const auto* c2 = dynamic_cast<const RawClusterv2*>(clus))
    {
      const float xr = c2->x_tower_raw();
      const float yr = c2->y_tower_raw();
      if (std::isfinite(xr) && std::isfinite(yr)) { xCG = xr; yCG = yr; return true; }
    }

    // 2) Fallback: rebuild the hit list and call Momenta (same as clusterizer)
    std::vector<EmcModule> hitlist;
    hitlist.reserve(std::distance(clus->get_towers().first, clus->get_towers().second));
    const int Nx = rec->GetNx();

    for (auto it = clus->get_towers().first; it != clus->get_towers().second; ++it) {
      const auto tg = geo->get_tower_geometry(it->first);
      if (!tg) continue;
      EmcModule m;
      m.ich = tg->get_bineta()*Nx + tg->get_binphi();
      m.amp = static_cast<float>(it->second);
      m.tof = 0.0f;
      hitlist.push_back(m);
    }
    if (hitlist.empty()) return false;

    float E=0, px=0, py=0, pxx=0, pyy=0, pyx=0;
    const float thr = rec->GetTowerThreshold();
    rec->Momenta(&hitlist, E, px, py, pxx, pyy, pyx, thr);
    if (E <= 0.0f) return false;

    xCG = px; yCG = py;  // Momenta returns tower-space CoG
    return true;
  }

  inline float foldToOneTwPitch(float dphi, int Nx)
  {
    const float tw = 2.f*static_cast<float>(M_PI)/static_cast<float>(Nx);
    float df = TVector2::Phi_mpi_pi(dphi);
    df = static_cast<float>(std::remainder(static_cast<double>(df), static_cast<double>(tw)));
    if (df <= -0.5f*tw) df += tw;
    if (df >   0.5f*tw) df -= tw;
    return df;
  }
} // anon


void PositionDependentCorrection::fillAshLogDx(
        RawCluster*                   clus,
        const TLorentzVector&         recoPhoton,
        const TLorentzVector&         truthPhoton,
        float                         vtxZ_in,
        const std::pair<float,float>& /*blockCord (unused)*/,
        int                           /*blockPhiBin (unused)*/,
        const std::vector<int>&       /*tower_phis (unused here)*/,
        const std::vector<float>&     /*tower_energies (unused here)*/)
{
  if (!clus || !m_geometry || !m_bemcRec) return;

  const float eReco = recoPhoton.E();
  const int   iE    = getEnergySlice(eReco);
  if (iE < 0 || iE >= N_Ebins) return;

  // (1) Raw tower-space CoG (what CorrectPosition would receive)
  float xCG{}, yCG{};
  if (!getTowerCG_fromPropsOrMomenta(clus, m_geometry, m_bemcRec, xCG, yCG)) return;

  // (2) Truth reference at shower depth
  const float phiTruth = TVector2::Phi_mpi_pi(truthPhoton.Phi());
  const float etaTruth = truthPhoton.Eta();
  const int   Nx       = m_bemcRec->GetNx();
  const float vtxZ     = vtxZ_in;

  // Prepare tower lists once (for log-weight scan)
  std::vector<int>   phiIdx; phiIdx.reserve(32);
  std::vector<int>   etaIdx; etaIdx.reserve(32);
  std::vector<float> twE;    twE.reserve(32);

  for (auto it = clus->get_towers().first; it != clus->get_towers().second; ++it) {
    const auto* tg = m_geometry->get_tower_geometry(it->first);
    if (!tg) continue;
    const float Ei = static_cast<float>(it->second);
    if (Ei <= 0.f) continue;
    phiIdx.push_back(tg->get_binphi());
    etaIdx.push_back(tg->get_bineta());
    twE   .push_back(Ei);
  }

  // Helpers to project tower-space (xT,yT) to global at shower depth
  auto phi_at_SD = [&](float xT, float yT){ return PDC_detail::cg2GlobalPhi(m_bemcRec, eReco, xT, yT); };
  auto eta_at_SD = [&](float xT, float yT){ return PDC_detail::cg2ShowerEta(m_bemcRec, eReco, xT, yT, vtxZ); };

  // ============================================================
  // (A) Ash scan (vary b) – in tower units; φ includes ripple term
  // ============================================================
  for (double bd : m_bScan)
  {
    const float b = static_cast<float>(bd);
    if (b <= 0.f) continue;

    // φ side: inverse-asinh around owner cell + legacy ripple
    float xCorr = doAshShift(xCG, b);
    const int   ix8 = int(xCG + 0.5f) / 8;
    const float x8  = xCG + 0.5f - (ix8 * 8) - 4.f;
    float dx = (std::fabs(x8) > 3.3f) ? 0.f : 0.10f * (x8 / 4.f);
    xCorr -= dx;

    const float phiAsh = phi_at_SD(xCorr, yCG);
    const float dphi   = foldToOneTwPitch(phiAsh - phiTruth, Nx);
    if (p_phi_rms2_vs_b_E[iE]) p_phi_rms2_vs_b_E[iE]->Fill(b, dphi*dphi);

    // η side: inverse-asinh on y only
    const float yCorr = doAshShift(yCG, b);
    const float etaAsh = eta_at_SD(xCG, yCorr);
    const float deta   = etaAsh - etaTruth;
    if (p_eta_rms2_vs_b_E[iE]) p_eta_rms2_vs_b_E[iE]->Fill(b, deta*deta);
  }

  // ============================================================
  // (B) Log-weight scan (vary w0) – recompute raw CoG in tower units
  // ============================================================
  auto computeLogCG = [&](float w0, float& x_log, float& y_log)->bool
  {
    if (twE.empty()) return false;
    const double sumE = std::accumulate(twE.begin(), twE.end(), 0.0);
    if (sumE <= 0.0) return false;

    const std::size_t iRef =
        std::distance(twE.begin(), std::max_element(twE.begin(), twE.end()));
    const int refPhi = phiIdx[iRef];

    double sw=0.0, sphi=0.0, seta=0.0;
    const int nx = Nx;

    for (std::size_t i=0; i<twE.size(); ++i)
    {
      const double Ei = twE[i]; if (Ei <= 0.0) continue;
      double wi = w0 + std::log(Ei / sumE); if (wi < 0.0) wi = 0.0;

      int phi = phiIdx[i];
      if (phi - refPhi < -nx/2) phi += nx;
      if (phi - refPhi >  nx/2) phi -= nx;

      sw   += wi;
      sphi += wi * static_cast<double>(phi);
      seta += wi * static_cast<double>(etaIdx[i]);
    }
    if (sw <= 0.0) return false;

    x_log = static_cast<float>(sphi / sw);
    while (x_log < -0.5f)      x_log += nx;
    while (x_log >= nx - 0.5f) x_log -= nx;

    y_log = static_cast<float>(seta / sw);  // η index linear
    return std::isfinite(x_log) && std::isfinite(y_log);
  };

  for (double w0d : m_w0Scan)
  {
    const float w0 = static_cast<float>(w0d);

    float xLog = std::numeric_limits<float>::quiet_NaN();
    float yLog = std::numeric_limits<float>::quiet_NaN();
    if (!computeLogCG(w0, xLog, yLog)) continue;

    const float phiLog = phi_at_SD(xLog, yLog);
    const float etaLog = eta_at_SD(xLog, yLog);

    const float dphi = foldToOneTwPitch(phiLog - phiTruth, Nx);
    const float deta = etaLog - etaTruth;

    if (p_phi_rms2_vs_w0_E[iE]) p_phi_rms2_vs_w0_E[iE]->Fill(w0, dphi*dphi);
    if (p_eta_rms2_vs_w0_E[iE]) p_eta_rms2_vs_w0_E[iE]->Fill(w0, deta*deta);
  }
}

// ============================================================================
//  Δφ – five flavours
// ============================================================================
void PositionDependentCorrection::fillDPhiAllVariants(
        RawCluster*                   clus,
        const TLorentzVector&         /*recoPhoton*/,
        const TLorentzVector&         truthPhoton,
        const std::pair<float,float>& blkCoord,     // (ηloc , φloc)
        int                           blkPhiCoarse,
        float                         vtxZ,
        TH1F*                         hRAW [N_Ebins],
        TH1F*                         hCP  [N_Ebins],
        bool                          fillGlobal )
{
  /* ── A) fast guards ───────────────────────────────────────────────── */
  if (!clus || !m_geometry || !m_bemcRec) return;
  const int vb = Verbosity();          // snapshot once – cheaper

  /* ── B) slice lookup ─────────────────────────────────────────────── */
  const float eReco = clus->get_energy();
  const int   iE    = getEnergySlice(eReco);
  if (iE < 0 || iE >= N_Ebins) return;

    if (vb > 0) {
      std::cout << "\n[fillDPhiAllVariants]  →  slice=" << iE
                << "  E=" << eReco << " GeV\n";
    }

    // ── C) cluster centre of gravity (tower grid) ─────────────────────
    float xCG{}, yCG{};
    if (!PDC_detail::clusterCentreOfGravity(clus, m_geometry, m_bemcRec, xCG, yCG))
    {
      if (vb > 1) std::cout << "  ‑ CG calculation failed – abort path\n";
      return;
    }
    if (vb > 2)
      std::cout << "  CG(x,y)=(" << xCG << ',' << yCG << ")\n";

    /* First-pass RAW CoG cross-check (tower units) using RawClusterv2. */
    {
      static bool s_checked_raw_phi_once = false;
      if (!s_checked_raw_phi_once)
      {
        constexpr float kTolTw = 5e-4f;
        const auto* c2 = dynamic_cast<const RawClusterv2*>(clus);
        if (!c2)
        {
          if (vb > 0)
            std::cout << "[φ-RAW check] cluster ID=" << clus->get_id()
                      << " is not RawClusterv2 — skipping validation.\n";
        }
        else
        {
          const float x_raw = c2->x_tower_raw();
          const float y_raw = c2->y_tower_raw();
          const float dx = std::fabs(x_raw - xCG);
          const float dy = std::fabs(y_raw - yCG);
          if (!std::isfinite(x_raw) || !std::isfinite(y_raw))
          {
            if (vb > 0)
              std::cout << "[φ-RAW check] non-finite v2 CoG — skipping.\n";
          }
          else if (dx > kTolTw || dy > kTolTw)
          {
            std::cerr << std::setprecision(7)
                      << "[WARN] RAW v2 CoG mismatch (φ path)\n"
                      << "  cluster ID=" << clus->get_id()
                      << "  E=" << eReco << " GeV\n"
                      << "  v2 (x,y)=(" << x_raw << "," << y_raw << ")\n"
                      << "  recomputed(x,y)=(" << xCG  << "," << yCG  << ")\n"
                      << "  |Δ|=(" << dx << "," << dy << ")  tol=" << kTolTw << '\n';
          }
          else if (vb > 0)
          {
            std::cout << "[φ-RAW check] v2 vs recomputed agree within "
                      << kTolTw << " tower units.\n";
          }
        }
        s_checked_raw_phi_once = true;
      }
    }

  /* ── D) lead‑tower geometry (needed for PDC path) ────────────────── */
  const auto lead =
        std::max_element(clus->get_towers().first, clus->get_towers().second,
                         [](auto& a, auto& b){ return a.second < b.second; });
  if (lead == clus->get_towers().second) return;
  const auto* tgLead = m_geometry->get_tower_geometry(lead->first);
  if (!tgLead) return;

  const double rFront = tgLead->get_center_radius();
  (void) rFront;          //  add the cast on its own line
  const double zFront = tgLead->get_center_z();
  const int    ixFine = tgLead->get_binphi();
  const int    iyFine = tgLead->get_bineta();

    if (vb > 2) {
      std::cout << "  lead-tower  rFront=" << rFront
                << "  zFront=" << zFront
                << "  ixFine="  << ixFine
                << "  iyFine="  << iyFine << '\n';
    }

    /* ── E) build the five flavour records ───────────────────────────── */
    PhiRec rec[5];
    const float phiTruth = TVector2::Phi_mpi_pi(truthPhoton.Phi());

    // define folding helpers BEFORE any use in this function
    auto __wrapPhi = [](float d){ return TVector2::Phi_mpi_pi(d); };
    const float __kTw_global = 2.f * static_cast<float>(M_PI) / 256.f; // one fine-tower pitch
    auto foldToTowerPitch = [&](float d)->float {
      float df = __wrapPhi(d);  // (−π, +π]
      df = static_cast<float>(
              std::remainder(static_cast<double>(df), static_cast<double>(__kTw_global))
           );                   // (−tw, +tw]
      if (df <= -0.5f * __kTw_global) df += __kTw_global;  // (−tw/2, +tw/2]
      if (df >   0.5f * __kTw_global) df -= __kTw_global;
      return df;
    };

    // NOTE: reuse the existing 'vb' declared earlier; do not redeclare it here.
    if (vb > 9)
    {
      std::cout << "\n[φ-variants] ENTER"
                << "  |  E=" << eReco
                << "  xCG=" << xCG << "  yCG=" << yCG
                << "  |  blkCoord(ηloc,φloc)=(" << blkCoord.first << "," << blkCoord.second << ")"
                << "  blkPhiCoarse=" << blkPhiCoarse
                << "  |  φ_truth=" << phiTruth
                << "  |  iE=" << iE << '\n';
    }
    /* ----------------------------- (1) CLUS-RAW ----------------------------- */
    m_bemcRec->SetPhiTiltVariant(BEmcRecCEMC::ETiltVariant::CLUS_RAW);
    if (vb > 9) std::cout << "[CLUS-RAW] SetPhiTiltVariant → CLUS_RAW\n";

    rec[0] = {"CLUS-RAW", xCG,
              cg2GlobalPhi(m_bemcRec, eReco, xCG, yCG)};
    if (vb > 9)
    {
      std::cout << "[CLUS-RAW] loc=" << rec[0].loc
                << "  φ_SD=" << rec[0].phi
                << "  (from tower2global at shower depth)\n";
    }

    /* ----------------------------- (2) CLUS-CP ------------------------------ */
    m_bemcRec->SetPhiTiltVariant(BEmcRecCEMC::ETiltVariant::CLUS_CP);
    if (vb > 9) std::cout << "[CLUS-CP] SetPhiTiltVariant → CLUS_CP\n";

    float xCP = xCG, yCP = yCG;
    if (vb > 9)
      std::cout << "[CLUS-CP] CorrectPosition  IN  : x=" << xCG << "  y=" << yCG
                << "  (E=" << eReco << ")\n";

    m_bemcRec->CorrectPosition(eReco, xCG, yCG, xCP, yCP);

    if (vb > 9)
      std::cout << "[CLUS-CP] CorrectPosition  OUT : xC=" << xCP << "  yC=" << yCP << '\n';

    /* First-pass CORR CoG cross-check (tower units) using RawClusterv2. */
    {
      static bool s_checked_corr_phi_once = false;
      if (!s_checked_corr_phi_once)
      {
        constexpr float kTolTw = 5e-4f;
        const auto* c2 = dynamic_cast<const RawClusterv2*>(clus);
        if (!c2)
        {
          if (vb > 0)
            std::cout << "[φ-CORR check] cluster ID=" << clus->get_id()
                      << " is not RawClusterv2 — skipping validation.\n";
        }
        else
        {
          const float x_cor = c2->x_tower_corr();
          const float y_cor = c2->y_tower_corr();
          const float dx = std::fabs(x_cor - xCP);
          const float dy = std::fabs(y_cor - yCP);
          if (!std::isfinite(x_cor) || !std::isfinite(y_cor))
          {
            if (vb > 0)
              std::cout << "[φ-CORR check] non-finite v2 CoG — skipping.\n";
          }
          else if (dx > kTolTw || dy > kTolTw)
          {
            std::cerr << std::setprecision(7)
                      << "[WARN] CORR v2 CoG mismatch (φ path)\n"
                      << "  cluster ID=" << clus->get_id()
                      << "  E=" << eReco << " GeV\n"
                      << "  v2 (xC,yC)=(" << x_cor << "," << y_cor << ")\n"
                      << "  recomputed(xC,yC)=(" << xCP   << "," << yCP   << ")\n"
                      << "  |Δ|=(" << dx << "," << dy << ")  tol=" << kTolTw << '\n';
          }
          else if (vb > 0)
          {
            std::cout << "[φ-CORR check] v2 vs recomputed agree within "
                      << kTolTw << " tower units.\n";
          }
        }
        s_checked_corr_phi_once = true;
      }
    }

    rec[1] = {"CLUS-CP", xCP,
              PDC_detail::cg2GlobalPhi(m_bemcRec, eReco, xCP, yCP)};
    if (vb > 9)
    {
      std::cout << "[CLUS-CP] loc=" << rec[1].loc
                << "  φ_SD=" << rec[1].phi << '\n';
    }


    /* ------------------------ (2b) CLUS-CP(EA) — four variants ---------------- */
    // Variant A: EA from Geometry  (also capture bφ, βη for diagnostics)
    m_bemcRec->SetPhiTiltVariant(BEmcRecCEMC::ETiltVariant::CLUS_CP_EA_GEOM);
    float xEA_geom = xCG, yEA_geom = yCG;
    float bphiEA_geom = std::numeric_limits<float>::quiet_NaN();
    float betaEA_geom = std::numeric_limits<float>::quiet_NaN();
    m_bemcRec->CorrectPositionEnergyAwareFromGeometry(eReco, xCG, yCG,
                                                      xEA_geom, yEA_geom,
                                                      &bphiEA_geom, &betaEA_geom);
    const float phiEA_geom  = cg2GlobalPhi(m_bemcRec, eReco, xEA_geom, yEA_geom);
    const float dphiEA_geom = foldToTowerPitch(phiEA_geom - phiTruth);
    if (vb > 3)
      std::cout << "    CLUS-CP(EA geom)  loc=" << xEA_geom
                << "  φ_SD=" << phiEA_geom
                << "  Δφ(folded)=" << dphiEA_geom
                << "  bφ=" << bphiEA_geom << '\n';

    // Variant B: EA with |η|-dependent + energy-dependent fits
    m_bemcRec->SetPhiTiltVariant(BEmcRecCEMC::ETiltVariant::CLUS_CP_EA_FIT_ETADEP);
    float xEA_fitEta = xCG, yEA_fitEta = yCG;
    m_bemcRec->CorrectPositionEnergyAwareEtaAndEnergyDep(eReco, xCG, yCG, xEA_fitEta, yEA_fitEta);
    const float phiEA_fitEta  = cg2GlobalPhi(m_bemcRec, eReco, xEA_fitEta, yEA_fitEta);
    const float dphiEA_fitEta = foldToTowerPitch(phiEA_fitEta - phiTruth);
    if (vb > 3)
      std::cout << "    CLUS-CP(EA |η|+E fit)  loc=" << xEA_fitEta
                << "  φ_SD=" << phiEA_fitEta
                << "  Δφ(folded)=" << dphiEA_fitEta << '\n';

    // Variant C: EA with energy-only fits (no |η| dependence)
    m_bemcRec->SetPhiTiltVariant(BEmcRecCEMC::ETiltVariant::CLUS_CP_EA_FIT_EONLY);
    float xEA_fitE = xCG, yEA_fitE = yCG;
    m_bemcRec->CorrectPositionEnergyAwareEnergyDepOnly(eReco, xCG, yCG, xEA_fitE, yEA_fitE);
    const float phiEA_fitE  = cg2GlobalPhi(m_bemcRec, eReco, xEA_fitE, yEA_fitE);
    const float dphiEA_fitE = foldToTowerPitch(phiEA_fitE - phiTruth);
    if (vb > 3)
      std::cout << "    CLUS-CP(EA E-only fit)  loc=" << xEA_fitE
                << "  φ_SD=" << phiEA_fitE
                << "  Δφ(folded)=" << dphiEA_fitE << '\n';

    // Variant D: mixed — φ(E-only), η(|η|+E)
    m_bemcRec->SetPhiTiltVariant(BEmcRecCEMC::ETiltVariant::CLUS_CP_EA_MIX);
    float xEA_mix = xCG, yEA_mix = yCG;
    m_bemcRec->CorrectPositionEnergyAwareEtaDepOnlyForEtaEnergyForPhi(eReco, xCG, yCG, xEA_mix, yEA_mix);
    const float phiEA_mix  = cg2GlobalPhi(m_bemcRec, eReco, xEA_mix, yEA_mix);
    const float dphiEA_mix = foldToTowerPitch(phiEA_mix - phiTruth);
    if (vb > 3)
      std::cout << "    CLUS-CP(EA φ:E-only, η:|η|+E)  loc=" << xEA_mix
                << "  φ_SD=" << phiEA_mix
                << "  Δφ(folded)=" << dphiEA_mix << '\n';


    /* ----------------------------- (3) CLUS-BCORR --------------------------- */
    {
      const float nX   = static_cast<float>(m_bemcRec->GetNx());
      const float bPhi = m_bValsPhi[iE];
      float       locB = xCG;

      if (vb > 9)
        std::cout << "[CLUS-BCORR] inputs: nX=" << nX
                  << "  bPhi(iE=" << iE << ")=" << bPhi
                  << "  xCG=" << xCG << '\n';

      if (bPhi > 1e-9F)
      {
        const int   ix0  = static_cast<int>(xCG + 0.5F);
        const float off  = xCG - ix0;
        const float corr = doPhiBlockCorr(off, bPhi);
        locB = ix0 + corr;

        int wrapL = 0, wrapR = 0;
        while (locB < -0.5F)    { locB += nX; ++wrapL; }
        while (locB >= nX-.5F)  { locB -= nX; ++wrapR; }

        if (vb > 9)
          std::cout << "[CLUS-BCORR] ix0=" << ix0 << "  off=" << off
                    << "  corr(off,b)=" << corr
                    << "  → locB(pre-wrap)=" << ix0 + corr
                    << "  | wraps(L/R)=" << wrapL << "/" << wrapR
                    << "  → locB(final)=" << locB << '\n';
      }
      else if (vb > 9)
      {
        std::cout << "[CLUS-BCORR] bPhi ≤ 0 → no corr, locB=xCG=" << locB << '\n';
      }

      m_bemcRec->SetPhiTiltVariant(BEmcRecCEMC::ETiltVariant::CLUS_BCORR);
      if (vb > 9) std::cout << "[CLUS-BCORR] SetPhiTiltVariant → CLUS_BCORR\n";

      rec[2] = {"CLUS-BCORR", locB,
                cg2GlobalPhi(m_bemcRec, eReco, locB, yCG)};

      if (vb > 9)
      {
        std::cout << "[CLUS-BCORR] loc=" << rec[2].loc
                  << "  φ_SD=" << rec[2].phi << '\n';
      }
    }

    /* ----------------------------- (4) PDC-RAW ------------------------------ */
    {
      const float phiFront = convertBlockToGlobalPhi(blkPhiCoarse, blkCoord.second);
      if (vb > 9)
        std::cout << "[PDC-RAW] φFront=" << phiFront
                  << "  | blkΦ=" << blkPhiCoarse
                  << "  φloc(raw)=" << blkCoord.second << '\n';

      auto ixFromBlockCenter = [](int blkCoarse, float local)->int {
        float l = local;
        if (l <= -0.5f || l > 1.5f) l = std::fmod(l + 2.0f, 2.0f);
        int ix = blkCoarse * 2 + static_cast<int>(std::floor(l + 0.5f));
        if (ix >= 256) ix -= 256;
        if (ix < 0)    ix += 256;
        return ix;
      };

      const int ixFineBlk = ixFromBlockCenter(blkPhiCoarse, blkCoord.second);
      const int blkEtaCoarse_forBlock = iyFine / 2; // from lead-tower row
      const int iyFineBlk = blkEtaCoarse_forBlock * 2
                          + ((blkCoord.first < 0.5f) ? 0 : 1);

      if (vb > 9)
        std::cout << "[PDC-RAW] fine indices: ixFineBlk=" << ixFineBlk
                  << "  iyFineBlk=" << iyFineBlk
                  << "  (blkEtaCoarse_fromLead=" << blkEtaCoarse_forBlock << ")\n";

      RawTowerDefs::keytype keyBlk =
          RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, iyFineBlk, ixFineBlk);
      const auto* geomBlk = m_geometry->get_tower_geometry(keyBlk);

      if (!geomBlk)
      {
        if (vb > 0)
          std::cerr << ANSI_RED
                    << "[PDC-RAW] geometry lookup FAILED for (iphi=" << ixFineBlk
                    << ", ieta=" << iyFineBlk << ") → φ point skipped"
                    << ANSI_RESET << '\n';

        rec[3] = {"PDC-RAW", blkCoord.second, std::numeric_limits<float>::quiet_NaN()};
      }
      else
      {
        const double rFrontBlk = geomBlk->get_center_radius();
        const double zFrontBlk = geomBlk->get_center_z();

        if (vb > 9)
          std::cout << "[PDC-RAW] geom: rFront=" << rFrontBlk
                    << "  zFront=" << zFrontBlk << '\n';

        m_bemcRec->SetPhiTiltVariant(BEmcRecCEMC::ETiltVariant::PDC_RAW);
        if (vb > 9) std::cout << "[PDC-RAW] SetPhiTiltVariant → PDC_RAW\n";

        const float phiSD = front2ShowerPhi(m_bemcRec, eReco,
                                            rFrontBlk , zFrontBlk ,
                                            phiFront,
                                            ixFineBlk , iyFineBlk);

        rec[3] = {"PDC-RAW", blkCoord.second, phiSD};

        if (vb > 9)
          std::cout << "[PDC-RAW] φloc=" << rec[3].loc
                    << "  φ_SD=" << rec[3].phi << '\n';
      }
    }

    /* ----------------------------- (5) PDC-CORR ----------------------------- */
    /* use LEGACY (originalEta) bϕ only */
    {
      const float bPhi_legacy = m_bValsPhi[iE];
      if (vb > 9)
        std::cout << "[PDC-CORR] bPhi_legacy(iE=" << iE << ")=" << bPhi_legacy
                  << "  |  φloc(raw)=" << blkCoord.second
                  << "  blkΦ(raw)=" << blkPhiCoarse << '\n';

      if (isFitDoneForPhi && bPhi_legacy > 1e-9F)
      {
        float loc = doPhiBlockCorr(blkCoord.second, bPhi_legacy);

        int blk = blkPhiCoarse;
        int wrapL = 0, wrapR = 0;
        while (loc <= -0.5F) { loc += 2.F; --blk; ++wrapL; }
        while (loc >   1.5F) { loc -= 2.F; ++blk; ++wrapR; }
        constexpr int kCoarsePhiBins = 128;
        if (blk < 0)             blk += kCoarsePhiBins;
        if (blk >= kCoarsePhiBins) blk -= kCoarsePhiBins;

        if (vb > 9)
          std::cout << "[PDC-CORR] loc(corr)=" << loc
                    << "  | blk(after fold)=" << blk
                    << "  | wraps(L/R)=" << wrapL << "/" << wrapR << '\n';

        const float phiFront = convertBlockToGlobalPhi(blk, loc);
        if (vb > 9)
          std::cout << "[PDC-CORR] φFront(corrected)=" << phiFront << '\n';

        auto ixFromBlockCenter = [](int blkCoarse, float local)->int {
          float l = local;
          if (l <= -0.5f || l > 1.5f) l = std::fmod(l + 2.0f, 2.0f);
          int ix = blkCoarse * 2 + static_cast<int>(std::floor(l + 0.5f));
          if (ix >= 256) ix -= 256;
          if (ix < 0)    ix += 256;
          return ix;
        };

        const int ixFineBlk = ixFromBlockCenter(blk, loc);
        const int blkEtaCoarse_forBlock = iyFine / 2;
        const int iyFineBlk = blkEtaCoarse_forBlock * 2
                            + ((blkCoord.first < 0.5f) ? 0 : 1);

        if (vb > 9)
          std::cout << "[PDC-CORR] fine indices: ixFineBlk=" << ixFineBlk
                    << "  iyFineBlk=" << iyFineBlk
                    << "  (blkEtaCoarse_fromLead=" << blkEtaCoarse_forBlock << ")\n";

        RawTowerDefs::keytype keyBlk =
            RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, iyFineBlk, ixFineBlk);
        const auto* geomBlk = m_geometry->get_tower_geometry(keyBlk);

        if (!geomBlk)
        {
          if (vb > 0)
            std::cerr << ANSI_RED
                      << "[PDC-CORR] geometry lookup FAILED for (iphi=" << ixFineBlk
                      << ", ieta=" << iyFineBlk << ") → φ point skipped"
                      << ANSI_RESET << '\n';

          rec[4] = {"PDC-CORR", loc, std::numeric_limits<float>::quiet_NaN()};
        }
        else
        {
          const double rFrontBlk = geomBlk->get_center_radius();
          const double zFrontBlk = geomBlk->get_center_z();

          if (vb > 9)
            std::cout << "[PDC-CORR] geom: rFront=" << rFrontBlk
                      << "  zFront=" << zFrontBlk << '\n';

          m_bemcRec->SetPhiTiltVariant(BEmcRecCEMC::ETiltVariant::PDC_CORR);
          if (vb > 9) std::cout << "[PDC-CORR] SetPhiTiltVariant → PDC_CORR\n";

          const float phiSD = front2ShowerPhi(m_bemcRec, eReco,
                                              rFrontBlk , zFrontBlk ,
                                              phiFront,
                                              ixFineBlk , iyFineBlk);

          rec[4] = {"PDC-CORR", loc, phiSD};

          if (vb > 9)
            std::cout << "[PDC-CORR] loc=" << rec[4].loc
                      << "  φ_SD=" << rec[4].phi << '\n';
        }
      }
      else
      {
        if (vb > 9)
          std::cout << "[PDC-CORR] legacy bΦ not available/enabled → copy PDC-RAW φ\n";
        rec[4] = {"PDC-CORR", blkCoord.second, rec[3].phi};
      }
    }

    if (vb > 9) std::cout << "[φ-variants] EXIT\n";

    for (auto& r : rec)
    {
      r.d = foldToTowerPitch(r.phi - phiTruth);
      if (!std::isfinite(r.phi)) ++g_nanPhi;

      if (vb > 3)
        std::cout << "    " << r.tag << "  loc=" << r.loc
                  << "  φ_SD=" << r.phi
                  << "  Δφ(folded)="   << r.d << '\n';
    }


   //rec[0] -->Clusterizer tower2global NO CP
   //rec[1] --> Clusterizer CP with tower2global
   //rec[2] --> nrg dep b corr with tower2global
   //rec[3] --> PDC global conv no corr
   //rec[4] --> PDC global conv w PDC corr

    float etaCLUSraw = std::numeric_limits<float>::quiet_NaN();
    float etaCLUScp  = std::numeric_limits<float>::quiet_NaN();
    float etaBCORR   = std::numeric_limits<float>::quiet_NaN();
    float etaPDCraw  = std::numeric_limits<float>::quiet_NaN();
    float etaPDCcorr = std::numeric_limits<float>::quiet_NaN();

    // CLUS variants → use Tower2Global → shower-depth η
    etaCLUSraw = cg2ShowerEta(m_bemcRec, eReco, xCG,          yCG,          vtxZ);
    etaCLUScp  = cg2ShowerEta(m_bemcRec, eReco, xCP,          yCP,          vtxZ);
    etaBCORR   = cg2ShowerEta(m_bemcRec, eReco, rec[2].loc,   yCG,          vtxZ);

    // PDC variants → map block (η) to global η (front-face); good axis for binning
    const int blkEtaCoarse = iyFine / 2; // use lead-tower row as the parent block
    etaPDCraw  = convertBlockToGlobalEta(blkEtaCoarse, blkCoord.first);
    etaPDCcorr = etaPDCraw;              // φ-correction doesn’t change η in this study

    /* ---- OUT-OF-WINDOW capture for CLUS & PDC variants (|Δφ| > 0.025 rad),
           gated by the tight |z| cut via fillGlobal ---------------------- */
    if (fillGlobal)
    {
      const float kOut = 0.025f;  // 25 mrad in radians

      // CLUS-branch residuals
      const float dClusRaw = rec[0].d;    // CLUS-RAW
      const float dClusCP  = rec[1].d;    // CLUS-CP

      // PDC-branch residuals
      const float dPdcRaw  = rec[3].d;    // PDC-RAW
      const float dPdcCorr = rec[4].d;    // PDC-CORR

      auto push_row = [&](float dphi,
                          std::vector<float>& vecDphi,
                          std::vector<float>& vecEta,
                          std::vector<float>& vecVz,
                          std::vector<float>& vecE,
                          float eta, float vz, float E)
      {
        vecDphi.push_back(dphi);
        vecEta .push_back(eta);
        vecVz  .push_back(vz);
        vecE   .push_back(E);
      };

      // ===== CLUS-RAW =====
      if (std::isfinite(dClusRaw)) {
        if ( dClusRaw >= +kOut) { ++m_clusRawPos1Tw; push_row(dClusRaw, m_clusRawPos1Tw_dphi, m_clusRawPos1Tw_eta, m_clusRawPos1Tw_vz, m_clusRawPos1Tw_E, etaCLUSraw, vtxZ, eReco); }
        if ( dClusRaw <= -kOut) { ++m_clusRawNeg1Tw; push_row(dClusRaw, m_clusRawNeg1Tw_dphi, m_clusRawNeg1Tw_eta, m_clusRawNeg1Tw_vz, m_clusRawNeg1Tw_E, etaCLUSraw, vtxZ, eReco); }
      }

      // ===== CLUS-CP =====
      if (std::isfinite(dClusCP)) {
        if ( dClusCP >= +kOut) { ++m_clusCPPos1Tw; push_row(dClusCP, m_clusCPPos1Tw_dphi, m_clusCPPos1Tw_eta, m_clusCPPos1Tw_vz, m_clusCPPos1Tw_E, etaCLUScp, vtxZ, eReco); }
        if ( dClusCP <= -kOut) { ++m_clusCPNeg1Tw; push_row(dClusCP, m_clusCPNeg1Tw_dphi, m_clusCPNeg1Tw_eta, m_clusCPNeg1Tw_vz, m_clusCPNeg1Tw_E, etaCLUScp, vtxZ, eReco); }
      }

      // ===== PDC-RAW =====
      if (std::isfinite(dPdcRaw)) {
        if ( dPdcRaw >= +kOut) { ++m_pdcRawPos1Tw; push_row(dPdcRaw, m_pdcRawPos1Tw_dphi, m_pdcRawPos1Tw_eta, m_pdcRawPos1Tw_vz, m_pdcRawPos1Tw_E, etaPDCraw, vtxZ, eReco); }
        if ( dPdcRaw <= -kOut) { ++m_pdcRawNeg1Tw; push_row(dPdcRaw, m_pdcRawNeg1Tw_dphi, m_pdcRawNeg1Tw_eta, m_pdcRawNeg1Tw_vz, m_pdcRawNeg1Tw_E, etaPDCraw, vtxZ, eReco); }
      }

      // ===== PDC-CORR =====
      if (std::isfinite(dPdcCorr)) {
        if ( dPdcCorr >= +kOut) { ++m_pdcCorrPos1Tw; push_row(dPdcCorr, m_pdcCorrPos1Tw_dphi, m_pdcCorrPos1Tw_eta, m_pdcCorrPos1Tw_vz, m_pdcCorrPos1Tw_E, etaPDCcorr, vtxZ, eReco); }
        if ( dPdcCorr <= -kOut) { ++m_pdcCorrNeg1Tw; push_row(dPdcCorr, m_pdcCorrNeg1Tw_dphi, m_pdcCorrNeg1Tw_eta, m_pdcCorrNeg1Tw_vz, m_pdcCorrNeg1Tw_E, etaPDCcorr, vtxZ, eReco); }
      }
    }




  #define HFill(H,V)    do{ if((H) && std::isfinite(V))             (H)->Fill(V); }while(0)
  #define HFill2(H,X,Y) do{ if((H) && std::isfinite(X) && std::isfinite(Y)) (H)->Fill((X),(Y)); }while(0)

    if (fillGlobal) {
          // 1D (existing + EA)
          HFill(hRAW [iE], rec[0].d);
          HFill(hCP  [iE], rec[1].d);
          HFill(h_phi_diff_cpCorrEA_geom_E[iE],                   dphiEA_geom);
          HFill(h_phi_diff_cpCorrEA_fitEtaDep_E[iE],              dphiEA_fitEta);
          HFill(h_phi_diff_cpCorrEA_fitEnergyOnly_E[iE],          dphiEA_fitE);
          HFill(h_phi_diff_cpCorrEA_fitPhiEnergy_etaEtaDep_E[iE], dphiEA_mix);
          HFill(h_phi_diff_cpBcorr_E [iE], rec[2].d);
          HFill(h_phi_diff_raw_E     [iE], rec[3].d);
          HFill(h_phi_diff_corrected_E[iE], rec[4].d);
        
          HFill2(h2_phi_diffEA_vs_bphi_E[iE], bphiEA_geom, dphiEA_geom);

          // 2D Δφ vs η (keep existing only)
          HFill2(h2_phi_diff_vsEta_RAW_E    [iE], etaCLUSraw, rec[0].d);
          HFill2(h2_phi_diff_vsEta_CP_E     [iE], etaCLUScp,  rec[1].d);
          HFill2(h2_phi_diff_vsEta_BCORR_E  [iE], etaBCORR,   rec[2].d);
          HFill2(h2_phi_diff_vsEta_PDCraw_E [iE], etaPDCraw,  rec[3].d);
          HFill2(h2_phi_diff_vsEta_PDCcorr_E[iE], etaPDCcorr, rec[4].d);
    }

    /* ---------- vertex-resolved histogramming (unchanged) ---------- */
    if (!skipVertexDep)
    {
        const float absVz   = std::fabs(vtxZ);
        const int   iVzAbs  = getVzSlice(absVz);     // 0 … N_VzBins-1 or −1
        const int   iVzSgn  = getVzSliceSigned(vtxZ);// 0 … 2*N_VzBins-1 or −1

        if (iVzAbs >= 0)                             /* |z| slices */
        {
            HFill(h_phi_diff_cpRaw_E_vz      [iE][iVzAbs], rec[0].d);
            HFill(h_phi_diff_cpCorr_E_vz     [iE][iVzAbs], rec[1].d);
            HFill(h_phi_diff_cpBcorr_E_vz    [iE][iVzAbs], rec[2].d);
            HFill(h_phi_diff_raw_E_vz        [iE][iVzAbs], rec[3].d);
            HFill(h_phi_diff_corrected_E_vz  [iE][iVzAbs], rec[4].d);
        }

        if (m_useSignedVz && iVzSgn >= 0)            /* +z / –z slices */
        {
            HFill(h_phi_diff_cpRaw_E_vzsgn      [iE][iVzSgn], rec[0].d);
            HFill(h_phi_diff_cpCorr_E_vzsgn     [iE][iVzSgn], rec[1].d);
            HFill(h_phi_diff_cpBcorr_E_vzsgn    [iE][iVzSgn], rec[2].d);
            HFill(h_phi_diff_raw_E_vzsgn        [iE][iVzSgn], rec[3].d);
            HFill(h_phi_diff_corrected_E_vzsgn  [iE][iVzSgn], rec[4].d);
        }
    }

   #undef  HFill
   #undef  HFill2
    
    // ---------- Δφ winner bookkeeping (RAW vs CLUS-CP vs 4×EA variants) ----------
    if (fillGlobal)
    {
      const float aRAW       = std::fabs(rec[0].d);
      const float aCP        = std::fabs(rec[1].d);
      const float aEA_geom   = std::fabs(dphiEA_geom);
      const float aEA_fitEta = std::fabs(dphiEA_fitEta);
      const float aEA_fitE   = std::fabs(dphiEA_fitE);
      const float aEA_mix    = std::fabs(dphiEA_mix);

      float best = aRAW; int which = 0; // 0=RAW, 1=CP, 2=geom, 3=fitEta, 4=fitE, 5=mix
      if (aCP        < best) { best = aCP;        which = 1; }
      if (aEA_geom   < best) { best = aEA_geom;   which = 2; }
      if (aEA_fitEta < best) { best = aEA_fitEta; which = 3; }
      if (aEA_fitE   < best) { best = aEA_fitE;   which = 4; }
      if (aEA_mix    < best) { best = aEA_mix;    which = 5; }

      switch (which)
      {
        case 0: ++m_phi6WayWin_RAW;      ++m_phi6WayWinByE_RAW[iE];      break;
        case 1: ++m_phi6WayWin_CP;       ++m_phi6WayWinByE_CP[iE];       break;
        case 2: ++m_phi6WayWin_EA_geom;  ++m_phi6WayWinByE_EA_geom[iE];  break;
        case 3: ++m_phi6WayWin_EA_fitEta;++m_phi6WayWinByE_EA_fitEta[iE];break;
        case 4: ++m_phi6WayWin_EA_fitE;  ++m_phi6WayWinByE_EA_fitE[iE];  break;
        case 5: ++m_phi6WayWin_EA_mix;   ++m_phi6WayWinByE_EA_mix[iE];   break;
      }

      if (vb >= 2)
      {
        const char* tags[6] = {
          "CLUS-RAW", "CLUS-CP", "CLUS-CP(EA geom)",
          "CLUS-CP(EA |eta|+E fits)", "CLUS-CP(EA E-only fits)",
          "CLUS-CP(EA φ:E-only, η:|η|+E)"
        };
        std::cout << "  WINNER(φ per-event, RAW/CP/EA*): "
                  << tags[which] << "  |Δφ|=" << std::fixed << std::setprecision(6) << best << '\n';
      }
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

    if (vb > 0) {
      std::cout << "\n[fillDEtaAllVariants] slice=" << iE
                << "  E=" << eReco << " GeV\n";
    }

    // ────────────────────────────────── CG of the cluster ────────────
    float xCG{}, yCG{};
    if (!PDC_detail::clusterCentreOfGravity(clus, m_geometry, m_bemcRec, xCG, yCG))
    {
      if (vb > 1) std::cout << "  ‑ CG calculation failed – abort path\n";
      return;
    }
    if (vb > 2)
      std::cout << "  CG(x,y)=(" << xCG << ',' << yCG << ")\n";

    /* First-pass RAW CoG cross-check (η path) using RawClusterv2. */
    {
      static bool s_checked_raw_eta_once = false;
      if (!s_checked_raw_eta_once)
      {
        constexpr float kTolTw = 5e-4f;
        const auto* c2 = dynamic_cast<const RawClusterv2*>(clus);
        if (!c2)
        {
          if (vb > 0)
            std::cout << "[η-RAW check] cluster ID=" << clus->get_id()
                      << " is not RawClusterv2 — skipping validation.\n";
        }
        else
        {
          const float x_raw = c2->x_tower_raw();
          const float y_raw = c2->y_tower_raw();
          const float dx = std::fabs(x_raw - xCG);
          const float dy = std::fabs(y_raw - yCG);
          if (!std::isfinite(x_raw) || !std::isfinite(y_raw))
          {
            if (vb > 0)
              std::cout << "[η-RAW check] non-finite v2 CoG — skipping.\n";
          }
          else if (dx > kTolTw || dy > kTolTw)
          {
            std::cerr << std::setprecision(7)
                      << "[WARN] RAW v2 CoG mismatch (η path)\n"
                      << "  cluster ID=" << clus->get_id()
                      << "  E=" << eReco << " GeV\n"
                      << "  v2 (x,y)=(" << x_raw << "," << y_raw << ")\n"
                      << "  recomputed(x,y)=(" << xCG   << "," << yCG   << ")\n"
                      << "  |Δ|=(" << dx << "," << dy << ")  tol=" << kTolTw << '\n';
          }
          else if (vb > 0)
          {
            std::cout << "[η-RAW check] v2 vs recomputed agree within "
                      << kTolTw << " tower units.\n";
          }
        }
        s_checked_raw_eta_once = true;
      }
    }


  /* ────────────────────────────────── lead‑tower geom ────────────── */
  const auto lead =
        std::max_element(clus->get_towers().first, clus->get_towers().second,
                         [](auto&a,auto&b){return a.second<b.second;});
  if (lead == clus->get_towers().second) return;
  const auto* tgLead = m_geometry->get_tower_geometry(lead->first);
  if (!tgLead) return;

  const double rFront = tgLead->get_center_radius();
  (void) rFront;          //  add the cast on its own line

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

    /* First-pass CORR CoG cross-check (η path) using RawClusterv2. */
    {
      static bool s_checked_corr_eta_once = false;
      if (!s_checked_corr_eta_once)
      {
        constexpr float kTolTw = 5e-4f;
        const auto* c2 = dynamic_cast<const RawClusterv2*>(clus);
        if (!c2)
        {
          if (vb > 0)
            std::cout << "[η-CORR check] cluster ID=" << clus->get_id()
                      << " is not RawClusterv2 — skipping validation.\n";
        }
        else
        {
          const float x_cor = c2->x_tower_corr();
          const float y_cor = c2->y_tower_corr();
          const float dx = std::fabs(x_cor - xCP);
          const float dy = std::fabs(y_cor - yCP);
          if (!std::isfinite(x_cor) || !std::isfinite(y_cor))
          {
            if (vb > 0)
              std::cout << "[η-CORR check] non-finite v2 CoG — skipping.\n";
          }
          else if (dx > kTolTw || dy > kTolTw)
          {
            std::cerr << std::setprecision(7)
                      << "[WARN] CORR v2 CoG mismatch (η path)\n"
                      << "  cluster ID=" << clus->get_id()
                      << "  E=" << eReco << " GeV\n"
                      << "  v2 (xC,yC)=(" << x_cor << "," << y_cor << ")\n"
                      << "  recomputed(xC,yC)=(" << xCP   << "," << yCP   << ")\n"
                      << "  |Δ|=(" << dx << "," << dy << ")  tol=" << kTolTw << '\n';
          }
          else if (vb > 0)
          {
            std::cout << "[η-CORR check] v2 vs recomputed agree within "
                      << kTolTw << " tower units.\n";
          }
        }
        s_checked_corr_eta_once = true;
      }
    }


    rec[1] = {"CLUS-CP", yCP,
              PDC_detail::cg2ShowerEta(m_bemcRec, eReco, xCP, yCP, vtxZ)};


    // 2b) CLUS-CP(EA) — four variants
    // Variant A: geometry
    float xEA_geom = xCG, yEA_geom = yCG;
    float bphiEA_geom = std::numeric_limits<float>::quiet_NaN();
    float betaEA_geom = std::numeric_limits<float>::quiet_NaN();
    m_bemcRec->CorrectPositionEnergyAwareFromGeometry(eReco, xCG, yCG,
                                                      xEA_geom, yEA_geom,
                                                      &bphiEA_geom, &betaEA_geom);
    const float etaEA_geom  = cg2ShowerEta(m_bemcRec, eReco, xEA_geom, yEA_geom, vtxZ);
    const float dEtaEA_geom = etaEA_geom - etaTruth;
    if (vb > 3)
      std::cout << "    CLUS-CP(EA geom)  loc=" << yEA_geom
                << "  η_SD=" << etaEA_geom
                << "  Δη=" << dEtaEA_geom
                << "  βη=" << betaEA_geom << '\n';

    // Variant B: |η|+E fits
    float xEA_fitEta = xCG, yEA_fitEta = yCG;
    m_bemcRec->CorrectPositionEnergyAwareEtaAndEnergyDep(eReco, xCG, yCG, xEA_fitEta, yEA_fitEta);
    const float etaEA_fitEta  = cg2ShowerEta(m_bemcRec, eReco, xEA_fitEta, yEA_fitEta, vtxZ);
    const float dEtaEA_fitEta = etaEA_fitEta - etaTruth;
    if (vb > 3)
      std::cout << "    CLUS-CP(EA |η|+E fit)  loc=" << yEA_fitEta
                << "  η_SD=" << etaEA_fitEta
                << "  Δη=" << dEtaEA_fitEta << '\n';

    // Variant C: E-only fits
    float xEA_fitE = xCG, yEA_fitE = yCG;
    m_bemcRec->CorrectPositionEnergyAwareEnergyDepOnly(eReco, xCG, yCG, xEA_fitE, yEA_fitE);
    const float etaEA_fitE  = cg2ShowerEta(m_bemcRec, eReco, xEA_fitE, yEA_fitE, vtxZ);
    const float dEtaEA_fitE = etaEA_fitE - etaTruth;
    if (vb > 3)
      std::cout << "    CLUS-CP(EA E-only fit)  loc=" << yEA_fitE
                << "  η_SD=" << etaEA_fitE
                << "  Δη=" << dEtaEA_fitE << '\n';

    // Variant D: φ(E-only), η(|η|+E)
    float xEA_mix = xCG, yEA_mix = yCG;
    m_bemcRec->CorrectPositionEnergyAwareEtaDepOnlyForEtaEnergyForPhi(eReco, xCG, yCG, xEA_mix, yEA_mix);
    const float etaEA_mix  = cg2ShowerEta(m_bemcRec, eReco, xEA_mix, yEA_mix, vtxZ);
    const float dEtaEA_mix = etaEA_mix - etaTruth;
    if (vb > 3)
      std::cout << "    CLUS-CP(EA φ:E-only, η:|η|+E)  loc=" << yEA_mix
                << "  η_SD=" << etaEA_mix
                << "  Δη=" << dEtaEA_mix << '\n';



  // 3) CLUS‑BCORR
  {
    const float bEta = m_bValsEta[iE];
    float loc = yCG;
    if (bEta > 1e-9F)
    {
        const int   iy0 = int(yCG + 0.5F);
        const float off = yCG - iy0;

        /* physical η‑pitch of this tower row (centre‑to‑centre) */
        const int iyUp = std::min(iy0 + 1, 95);           // stay inside 0…95
        const int iyDn = std::max(iy0 - 1, 0);

        const RawTowerGeom* gU = m_geometry->get_tower_geometry(
                RawTowerDefs::encode_towerid(RawTowerDefs::CEMC,
                                             iyUp, blkPhiCoarse * 2));
        const RawTowerGeom* gD = m_geometry->get_tower_geometry(
                RawTowerDefs::encode_towerid(RawTowerDefs::CEMC,
                                             iyDn, blkPhiCoarse * 2));

        float dEtaPitch = 2.2f / 96.0f;                   // safe default
        if (gU && gD)                                     // use real pitch if possible
        {
            const float etaU = std::asinh(gU->get_center_z() /
                                          gU->get_center_radius());
            const float etaD = std::asinh(gD->get_center_z() /
                                          gD->get_center_radius());
            dEtaPitch = 0.5f * (etaU - etaD);
        }

        loc = iy0 + doEtaBlockCorr(off, bEta, dEtaPitch);

    }
    rec[2] = {"CLUS-BCORR", loc,
        cg2ShowerEta(m_bemcRec, eReco, xCG,  loc, vtxZ)}; // consistent pair
  }

    // 4) PDC-RAW
    {
      const float etaFront = convertBlockToGlobalEta(blkEtaCoarse, blkCoord.first);

      if (vb > 0)
      {
        std::cout << ANSI_CYAN << "[PDC-RAW η] ENTER" << ANSI_RESET << "\n"
                  << "  blkEtaCoarse=" << blkEtaCoarse
                  << "  blkPhiCoarse=" << blkPhiCoarse
                  << "  localEta(raw)=" << blkCoord.first
                  << "  localPhi(raw)=" << blkCoord.second
                  << "  → etaFront=" << etaFront << "\n";
      }

      /* fine-index of the block coordinate with open-right folding + clamp */
      auto ixFromBlockCenter = [](int blkCoarse, float local)->int {
        float l = local;
        if (l <= -0.5f || l >= 1.5f) {
          l = std::fmod(l + 2.0f, 2.0f);
          if (l < 0.f) l += 2.0f;
        }
        int ix = blkCoarse * 2 + static_cast<int>(std::floor(l + 0.5f));
        if (ix >= 256) ix -= 256;
        if (ix < 0)    ix += 256;
        return ix;
      };
      auto iyFromBlockCenter = [](int blkCoarse, float local)->int {
        float l = local;
        if (l <= -0.5f || l >= 1.5f) {
          l = std::fmod(l + 2.0f, 2.0f);
          if (l < 0.f) l += 2.0f;
        }
        int iy = blkCoarse * 2 + static_cast<int>(std::floor(l + 0.5f));
        if (iy < 0)  iy = 0;
        if (iy > 95) iy = 95;
        return iy;
      };

      const int ixFineBlk = ixFromBlockCenter(blkPhiCoarse, blkCoord.second);
      const int iyFine    = iyFromBlockCenter(blkEtaCoarse, blkCoord.first);

      if (vb > 0)
      {
        std::cout << "  ixFineBlk=" << ixFineBlk
                  << "  iyFine="    << iyFine
                  << "  (Nx=256, Ny=96)\n";
      }

      RawTowerDefs::keytype keyBlkEta =
          RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, iyFine, ixFineBlk);
      const auto* geomBlkEta = m_geometry->get_tower_geometry(keyBlkEta);

      if (!geomBlkEta)
      {
        if (vb > 0)
          std::cerr << ANSI_RED
                    << "  ✘ geometry lookup FAILED  key=" << keyBlkEta
                    << "  (iphi=" << ixFineBlk << ", ieta=" << iyFine << ")\n"
                    << "  → skipping PDC-RAW η point"
                    << ANSI_RESET << "\n";
        rec[3] = {"PDC-RAW", blkCoord.first, std::numeric_limits<float>::quiet_NaN()};
      }
      else
      {
        const double rFrontBlk = geomBlkEta->get_center_radius();
        const double zFrontBlk = geomBlkEta->get_center_z();

        if (vb > 0)
        {
          std::cout << "  geom: rFront=" << rFrontBlk
                    << "  zFront="      << zFrontBlk
                    << "  phiFront(raw via convertBlockToGlobalPhi)=" << phiFrontRaw
                    << "\n";
        }

        const float etaSD = front2ShowerEta(m_bemcRec, eReco,
                                            rFrontBlk , ixFineBlk , iyFine ,
                                            etaFront , phiFrontRaw ,
                                            vtxZ);

        if (vb > 0)
        {
          std::cout << "  η_SD=" << etaSD << "  (from front2ShowerEta)\n"
                    << ANSI_CYAN << "[PDC-RAW η] EXIT" << ANSI_RESET << "\n";
        }

        rec[3] = {"PDC-RAW", blkCoord.first, etaSD};
      }
    }


    // 5) PDC-CORR  — use LEGACY (originalEta) bη only
    {
      const float bEta_legacy = m_bValsEta[iE];      // ← originalEta table
      if (isFitDoneForEta && bEta_legacy > 1e-9F)
      {
        int   blk = blkEtaCoarse;

        float dEtaPitchBlk =
            0.5f * ( convertBlockToGlobalEta(blk,  1.5f)
                   - convertBlockToGlobalEta(blk, -0.5f) );
        if (!std::isfinite(dEtaPitchBlk) || dEtaPitchBlk <= 0.0f)
          dEtaPitchBlk = 2.2f / 96.0f;

        float loc = doEtaBlockCorr(blkCoord.first, bEta_legacy, dEtaPitchBlk);

        if (vb > 0)
        {
          std::cout << ANSI_YELLOW << "[PDC-CORR η] ENTER" << ANSI_RESET << "\n"
                    << "  bEta=" << bEta_legacy
                    << "  dEtaPitchBlk=" << dEtaPitchBlk
                    << "  localEta(raw)=" << blkCoord.first
                    << "  localEta(corr, pre-wrap)=" << loc
                    << "  blkEtaCoarse(raw)=" << blk << "\n";
        }

        while (loc <= -0.5F) { loc += 2.F; --blk; }
        while (loc >   1.5F) { loc -= 2.F; ++blk; }
        constexpr int kCoarseEtaBins = 48;
        if (blk < 0)  blk += kCoarseEtaBins;
        if (blk >= kCoarseEtaBins) blk -= kCoarseEtaBins;

        if (vb > 0)
        {
          std::cout << "  localEta(corr, folded)=" << loc
                    << "  blkEtaCoarse(folded)="  << blk << "\n";
        }

        const float etaFront = convertBlockToGlobalEta(blk, loc);

        auto ixFromBlockCenter = [](int blkCoarse, float local)->int {
          float l = local;
          if (l <= -0.5f || l >= 1.5f) {
            l = std::fmod(l + 2.0f, 2.0f);
            if (l < 0.f) l += 2.0f;
          }
          int ix = blkCoarse * 2 + static_cast<int>(std::floor(l + 0.5f));
          if (ix >= 256) ix -= 256;
          if (ix < 0)    ix += 256;
          return ix;
        };
        auto iyFromBlockCenter = [](int blkCoarse, float local)->int {
          float l = local;
          if (l <= -0.5f || l >= 1.5f) {
            l = std::fmod(l + 2.0f, 2.0f);
            if (l < 0.f) l += 2.0f;
          }
          int iy = blkCoarse * 2 + static_cast<int>(std::floor(l + 0.5f));
          if (iy < 0)  iy = 0;
          if (iy > 95) iy = 95;
          return iy;
        };

        const int ixFineBlk = ixFromBlockCenter(blkPhiCoarse, blkCoord.second);
        const int iyFine    = iyFromBlockCenter(blk,           loc);

        if (vb > 0)
        {
          std::cout << "  ixFineBlk=" << ixFineBlk
                    << "  iyFine="    << iyFine
                    << "  → etaFront="<< etaFront << "\n";
        }

        RawTowerDefs::keytype keyBlkEta =
            RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, iyFine, ixFineBlk);
        const auto* geomBlkEta = m_geometry->get_tower_geometry(keyBlkEta);

        if (!geomBlkEta)
        {
          if (vb > 0)
            std::cerr << ANSI_RED
                      << "  ✘ geometry lookup FAILED  key=" << keyBlkEta
                      << "  (iphi=" << ixFineBlk << ", ieta=" << iyFine << ")\n"
                      << "  → skipping PDC-CORR η point"
                      << ANSI_RESET << "\n";
          rec[4] = {"PDC-CORR", loc, std::numeric_limits<float>::quiet_NaN()};
        }
        else
        {
          const double rFrontBlk = geomBlkEta->get_center_radius();

          if (vb > 0)
          {
            std::cout << "  geom: rFront=" << rFrontBlk
                      << "  φFront(raw via convertBlockToGlobalPhi)=" << phiFrontRaw
                      << "\n";
          }

          const float etaSD = front2ShowerEta(m_bemcRec, eReco,
                                              rFrontBlk , ixFineBlk , iyFine ,
                                              etaFront , phiFrontRaw ,
                                              vtxZ);

          if (vb > 0)
          {
            std::cout << "  η_SD=" << etaSD << "  (from front2ShowerEta)\n"
                      << ANSI_YELLOW << "[PDC-CORR η] EXIT" << ANSI_RESET << "\n";
          }

          rec[4] = {"PDC-CORR", loc, etaSD};
        }
      }
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
    // --- NEW 2D fill macro (like in φ function)
    #define HFill2(H,X,Y) do{ if((H) && std::isfinite(X) && std::isfinite(Y)) (H)->Fill((X),(Y)); }while(0)

    // --- NEW: fill Δη vs βη using the geometry EA variant
    HFill2(h2_eta_diffEA_vs_beta_E[iE], betaEA_geom, dEtaEA_geom);
    
    /* ── book‑keeping: |Δη| outside ±0.04 ───────────────────────────── */
    for (int i = 0; i < 5; ++i)
      if (std::fabs(rec[i].d) > 0.04f)
      {
        switch (i)
        {
          case 0: ++m_etaOutCLUSraw;   break;
          case 1: ++m_etaOutCLUScp;    break;
          case 2: ++m_etaOutCLUSbcorr; break;
          case 3: ++m_etaOutPDCraw;    break;
          case 4: ++m_etaOutPDCcorr;   break;
        }
        if (vb > 5)
          std::cout << "    |Δη|>0.04  (" << rec[i].tag << ")  = "
                    << rec[i].d << '\n';
      }

    //rec[0] -->Clusterizer tower2global NO CP
    //rec[1] --> Clusterizer CP with tower2global
    //rec[2] --> nrg dep b corr with tower2global
    //rec[3] --> PDC global conv no corr
    //rec[4] --> PDC global conv w PDC corr
    
#define HFill(H,V)  do{ if((H)&&std::isfinite(V)) (H)->Fill(V); }while(0)
    if (fillGlobal)
    {
      HFill(hRAW [iE], rec[0].d);
      HFill(hCP  [iE], rec[1].d);
        HFill(h_eta_diff_cpCorrEA_geom_E[iE],                   dEtaEA_geom);
        HFill(h_eta_diff_cpCorrEA_fitEtaDep_E[iE],              dEtaEA_fitEta);
        HFill(h_eta_diff_cpCorrEA_fitEnergyOnly_E[iE],          dEtaEA_fitE);
        HFill(h_eta_diff_cpCorrEA_fitPhiEnergy_etaEtaDep_E[iE], dEtaEA_mix);
      HFill(h_eta_diff_cpBcorr_E [iE], rec[2].d);
      HFill(h_eta_diff_raw_E     [iE], rec[3].d);
      HFill(h_eta_diff_corrected_E[iE], rec[4].d);
    }
    /* ---------- vertex-resolved histogramming ---------- */
    if (!skipVertexDep)
    {
        const float absVz   = std::fabs(vtxZ);
        const int   iVzAbs  = getVzSlice(absVz);      // 0 … N_VzBins-1 or −1
        const int   iVzSgn  = getVzSliceSigned(vtxZ); // 0 … 2*N_VzBins-1 or −1

        if (iVzAbs >= 0)                              /* |z| slices */
        {
            HFill(h_eta_diff_cpRaw_E_vz      [iE][iVzAbs], rec[0].d);
            HFill(h_eta_diff_cpCorr_E_vz     [iE][iVzAbs], rec[1].d);
            HFill(h_eta_diff_cpBcorr_E_vz    [iE][iVzAbs], rec[2].d);
            HFill(h_eta_diff_raw_E_vz        [iE][iVzAbs], rec[3].d);
            HFill(h_eta_diff_corrected_E_vz  [iE][iVzAbs], rec[4].d);
        }

        if (m_useSignedVz && iVzSgn >= 0)             /* +z / –z slices */
        {
            HFill(h_eta_diff_cpRaw_E_vzsgn      [iE][iVzSgn], rec[0].d);
            HFill(h_eta_diff_cpCorr_E_vzsgn     [iE][iVzSgn], rec[1].d);
            HFill(h_eta_diff_cpBcorr_E_vzsgn    [iE][iVzSgn], rec[2].d);
            HFill(h_eta_diff_raw_E_vzsgn        [iE][iVzSgn], rec[3].d);
            HFill(h_eta_diff_corrected_E_vzsgn  [iE][iVzSgn], rec[4].d);
        }
    }

#undef  HFill

    // ----- Δη winner bookkeeping (RAW vs CLUS-CP vs 4×EA variants) -----
    if (fillGlobal)
    {
      const float aRAW       = std::fabs(rec[0].d);
      const float aCP        = std::fabs(rec[1].d);
      const float aEA_geom   = std::fabs(dEtaEA_geom);
      const float aEA_fitEta = std::fabs(dEtaEA_fitEta);
      const float aEA_fitE   = std::fabs(dEtaEA_fitE);
      const float aEA_mix    = std::fabs(dEtaEA_mix);

      float best = aRAW; int which = 0; // 0=RAW, 1=CP, 2=geom, 3=fitEta, 4=fitE, 5=mix
      if (aCP        < best) { best = aCP;        which = 1; }
      if (aEA_geom   < best) { best = aEA_geom;   which = 2; }
      if (aEA_fitEta < best) { best = aEA_fitEta; which = 3; }
      if (aEA_fitE   < best) { best = aEA_fitE;   which = 4; }
      if (aEA_mix    < best) { best = aEA_mix;    which = 5; }

      switch (which)
      {
        case 0: ++m_eta6WayWin_RAW;      ++m_eta6WayWinByE_RAW[iE];      break;
        case 1: ++m_eta6WayWin_CP;       ++m_eta6WayWinByE_CP[iE];       break;
        case 2: ++m_eta6WayWin_EA_geom;  ++m_eta6WayWinByE_EA_geom[iE];  break;
        case 3: ++m_eta6WayWin_EA_fitEta;++m_eta6WayWinByE_EA_fitEta[iE];break;
        case 4: ++m_eta6WayWin_EA_fitE;  ++m_eta6WayWinByE_EA_fitE[iE];  break;
        case 5: ++m_eta6WayWin_EA_mix;   ++m_eta6WayWinByE_EA_mix[iE];   break;
      }

      if (vb >= 2)
      {
        const char* tags[6] = {
          "CLUS-RAW", "CLUS-CP", "CLUS-CP(EA geom)",
          "CLUS-CP(EA |eta|+E fits)", "CLUS-CP(EA E-only fits)",
          "CLUS-CP(EA φ:E-only, η:|η|+E)"
        };
        std::cout << "  WINNER(η per-event, RAW/CP/EA*): "
                  << tags[which] << "  |Δη|=" << std::fixed << std::setprecision(6) << best << '\n';
      }
    }

}



// ============================================================================
// Helper • loop over ALL “other” clusters, build π0 candidates, apply cuts
//         and fill (reco + truth) histograms.
//
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

  RawClusterContainer::ConstRange cRange2 = clusterContainer->getClusters();

    // Early-out for single-photon (0 or 1 cluster) events
    {
        auto it = cRange2.first;
        if (it == cRange2.second) return;       // no clusters
        auto it_next = it; ++it_next;
        if (it_next == cRange2.second) return;  // only one cluster → no pairs to build
    }

    // Materialize the "first" cluster for this pairing pass
  RawCluster* clus1 = cIt1->second;
  if (!clus1) return;

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

    if (m_isSimulation && pi0Reco.Pt() < pi0ptcut)       continue;

    // ----------------------------------------------------------------------
    // 4)  Energy‑slice index (use **leading‑γ** energy!)
    // ----------------------------------------------------------------------
    const float  eLead  = std::max(clusE, clus2E);
    const int    iSlice = getEnergySlice(eLead);     // −1 ⇒ outside table

      // ----------------------------------------------------------------------
      // 5)  -- ALWAYS -- fill slice histograms (pass-1 input)
      // ----------------------------------------------------------------------
      if (h_mE_raw )
        h_mE_raw ->Fill(pi0Reco.M(), eLead);
      if (h_mE_corr)
        h_mE_corr->Fill(pi0Reco.M(), eLead);

      // --- NEW: per-variant π0 mass histograms for truth-matched π0 pairs ---
      const int i = iSlice;
      if (isSimulation && i >= 0 && i < N_Ebins)
      {
        // Vertex Z to mimic your fillDPhi/DEta projections
        const float vtxZ = static_cast<float>(vertex.z());

        // (a) get mother ids for both photons by matching to stored π0 daughters
        auto motherIdFor = [&](const TLorentzVector& pReco)->int {
          for (const auto& tp : m_truth_pi0_photons) {
            const float dR    = pReco.DeltaR(tp.p4);
            const float ratio = pReco.E() / tp.p4.E();
            if (dR < 0.02f && ratio > 0.7f && ratio < 1.5f) return tp.mother_id;
          }
          return -1;
        };
        const int mother1 = motherIdFor(photon1);
        const int mother2 = motherIdFor(photon2);

        if (mother1 >= 0 && mother1 == mother2)
        {
          // (b) compute block coords for clus2 (like you do for clus1)
          std::pair<float,float> blkCoord2{0.5f,0.5f};
          int blkPhiCoarse2 = 0, blkEtaCoarse2 = 0;
          {
            std::vector<int> towerEtas2, towerPhis2; std::vector<float> towerEs2;
            RawCluster::TowerConstRange tcr2 = clus2->get_towers();
            for (auto it = tcr2.first; it != tcr2.second; ++it) {
              const auto tg = m_geometry->get_tower_geometry(it->first);
              if (!tg) continue;
              towerEtas2.push_back(tg->get_bineta());
              towerPhis2.push_back(tg->get_binphi());
              towerEs2  .push_back(static_cast<float>(it->second));
            }
            blkCoord2 = getBlockCord(towerEtas2, towerPhis2, towerEs2, blkPhiCoarse2, blkEtaCoarse2);
          }

            // (c) local builder for TLVs per variant (reuses your geometry transforms)
            auto buildPhotonTLV = [&](RawCluster* clus,
                                     VarPi0 variant,
                                     const std::pair<float,float>& blk,
                                     int blkPhiCoarseArg)->TLorentzVector
            {
              TLorentzVector out(0,0,0,0);
              float xCG=0, yCG=0;
              PDC_detail::clusterCentreOfGravity(clus, m_geometry, m_bemcRec, xCG, yCG);
              const float E = clus->get_energy();

              float etaSD = std::numeric_limits<float>::quiet_NaN();
              float phiSD = std::numeric_limits<float>::quiet_NaN();

              switch (variant) {
                case VarPi0::CLUS_RAW: {
                  etaSD = PDC_detail::cg2ShowerEta(m_bemcRec, E, xCG, yCG, vtxZ);
                  phiSD = PDC_detail::cg2GlobalPhi(m_bemcRec, E, xCG, yCG);
                  break;
                }
                case VarPi0::CLUS_CP: {
                  float xC=xCG, yC=yCG;
                  m_bemcRec->CorrectPosition(E, xCG, yCG, xC, yC);
                  etaSD = PDC_detail::cg2ShowerEta(m_bemcRec, E, xC, yC, vtxZ);
                  phiSD = PDC_detail::cg2GlobalPhi (m_bemcRec, E, xC, yC);
                  break;
                }
                case VarPi0::EA_GEOM: {
                  float xEA=xCG, yEA=yCG; float bphiEA=0.f, betaEA=0.f;
                  m_bemcRec->CorrectPositionEnergyAwareFromGeometry(E, xCG, yCG, xEA, yEA, &bphiEA, &betaEA);
                  etaSD = PDC_detail::cg2ShowerEta(m_bemcRec, E, xEA, yEA, vtxZ);
                  phiSD = PDC_detail::cg2GlobalPhi (m_bemcRec, E, xEA, yEA);
                  break;
                }
                case VarPi0::EA_FIT_ETADEP: {
                  float xEA=xCG, yEA=yCG;
                  m_bemcRec->CorrectPositionEnergyAwareEtaAndEnergyDep(E, xCG, yCG, xEA, yEA);
                  etaSD = PDC_detail::cg2ShowerEta(m_bemcRec, E, xEA, yEA, vtxZ);
                  phiSD = PDC_detail::cg2GlobalPhi (m_bemcRec, E, xEA, yEA);
                  break;
                }
                case VarPi0::EA_FIT_EONLY: {
                  float xEA=xCG, yEA=yCG;
                  m_bemcRec->CorrectPositionEnergyAwareEnergyDepOnly(E, xCG, yCG, xEA, yEA);
                  etaSD = PDC_detail::cg2ShowerEta(m_bemcRec, E, xEA, yEA, vtxZ);
                  phiSD = PDC_detail::cg2GlobalPhi (m_bemcRec, E, xEA, yEA);
                  break;
                }
                case VarPi0::EA_MIX: {
                  float xEA=xCG, yEA=yCG;
                  m_bemcRec->CorrectPositionEnergyAwareEtaDepOnlyForEtaEnergyForPhi(E, xCG, yCG, xEA, yEA);
                  etaSD = PDC_detail::cg2ShowerEta(m_bemcRec, E, xEA, yEA, vtxZ);
                  phiSD = PDC_detail::cg2GlobalPhi (m_bemcRec, E, xEA, yEA);
                  break;
                }
                case VarPi0::PDC_RAW:
                case VarPi0::PDC_CORR: {
                  std::pair<float,float> blkLoc = blk;
                  const int sliceIdx = iSlice; // capture outer slice index
                  if (variant == VarPi0::PDC_CORR) {
                    if (isFitDoneForEta && m_bValsEta[sliceIdx] > 1e-9f) blkLoc.first  = doEtaBlockCorr(blkLoc.first,  m_bValsEta[sliceIdx]);
                    if (isFitDoneForPhi && m_bValsPhi[sliceIdx] > 1e-9f) blkLoc.second = doPhiBlockCorr(blkLoc.second, m_bValsPhi[sliceIdx]);
                  }

                  auto ixFromBlockCenter = [](int blkCoarse, float local)->int {
                    float l = local;
                    if (l <= -0.5f || l > 1.5f) l = std::fmod(l + 2.0f, 2.0f);
                    int ix = blkCoarse * 2 + int(std::floor(l + 0.5f));
                    if (ix >= 256) ix -= 256;
                    if (ix < 0)    ix += 256;
                    return ix;
                  };

                  const auto lead = std::max_element(clus->get_towers().first, clus->get_towers().second,
                                                     [](const auto& lhs, const auto& rhs){ return lhs.second < rhs.second; });
                  if (lead == clus->get_towers().second) return out;
                  const auto* tgLead = m_geometry->get_tower_geometry(lead->first);
                  const int iyFineLead = tgLead ? tgLead->get_bineta() : 0;
                  const int blkEtaCoarseFromLead = iyFineLead / 2;

                  const int ixFineBlk = ixFromBlockCenter(blkPhiCoarseArg, blkLoc.second);
                  const int iyFineBlk = blkEtaCoarseFromLead * 2 + ((blkLoc.first < 0.5f) ? 0 : 1);

                  RawTowerDefs::keytype keyBlk = RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, iyFineBlk, ixFineBlk);
                  const auto* geomBlk = m_geometry->get_tower_geometry(keyBlk);
                  if (!geomBlk) return out;

                  const double rFrontBlk = geomBlk->get_center_radius();
                  const double zFrontBlk = geomBlk->get_center_z();

                  const float phiFront = convertBlockToGlobalPhi(blkPhiCoarseArg, blkLoc.second);
                  const float etaFront = convertBlockToGlobalEta(blkEtaCoarseFromLead,  blkLoc.first);

                  phiSD = PDC_detail::front2ShowerPhi(m_bemcRec, E, rFrontBlk, zFrontBlk, phiFront, ixFineBlk, iyFineBlk);
                  etaSD = PDC_detail::front2ShowerEta(m_bemcRec, E, rFrontBlk, ixFineBlk, iyFineBlk, etaFront, phiFront, vtxZ);
                  break;
                }
                default: return out;
              }

              const double pt = E / std::cosh(etaSD);
              out.SetPtEtaPhiE(pt, etaSD, phiSD, E);
              return out;
            };


            // Build eight variants for each cluster
            TLorentzVector p1[static_cast<int>(VarPi0::NVAR)], p2v[static_cast<int>(VarPi0::NVAR)];
            p1[static_cast<int>(VarPi0::CLUS_RAW)]  = buildPhotonTLV(clus1, VarPi0::CLUS_RAW,  blkCoord,  blkPhiCoarse);
            p1[static_cast<int>(VarPi0::CLUS_CP)]   = buildPhotonTLV(clus1, VarPi0::CLUS_CP,   blkCoord,  blkPhiCoarse);
            p1[static_cast<int>(VarPi0::EA_GEOM)]   = buildPhotonTLV(clus1, VarPi0::EA_GEOM,   blkCoord,  blkPhiCoarse);
            p1[static_cast<int>(VarPi0::EA_FIT_ETADEP)] = buildPhotonTLV(clus1, VarPi0::EA_FIT_ETADEP, blkCoord,  blkPhiCoarse);
            p1[static_cast<int>(VarPi0::EA_FIT_EONLY)]  = buildPhotonTLV(clus1, VarPi0::EA_FIT_EONLY,  blkCoord,  blkPhiCoarse);
            p1[static_cast<int>(VarPi0::EA_MIX)]    = buildPhotonTLV(clus1, VarPi0::EA_MIX,    blkCoord,  blkPhiCoarse);
            p1[static_cast<int>(VarPi0::PDC_RAW)]   = buildPhotonTLV(clus1, VarPi0::PDC_RAW,   blkCoord,  blkPhiCoarse);
            p1[static_cast<int>(VarPi0::PDC_CORR)]  = buildPhotonTLV(clus1, VarPi0::PDC_CORR,  blkCoord,  blkPhiCoarse);

            p2v[static_cast<int>(VarPi0::CLUS_RAW)]  = buildPhotonTLV(clus2, VarPi0::CLUS_RAW,  blkCoord2, blkPhiCoarse2);
            p2v[static_cast<int>(VarPi0::CLUS_CP)]   = buildPhotonTLV(clus2, VarPi0::CLUS_CP,   blkCoord2, blkPhiCoarse2);
            p2v[static_cast<int>(VarPi0::EA_GEOM)]   = buildPhotonTLV(clus2, VarPi0::EA_GEOM,   blkCoord2, blkPhiCoarse2);
            p2v[static_cast<int>(VarPi0::EA_FIT_ETADEP)] = buildPhotonTLV(clus2, VarPi0::EA_FIT_ETADEP, blkCoord2, blkPhiCoarse2);
            p2v[static_cast<int>(VarPi0::EA_FIT_EONLY)]  = buildPhotonTLV(clus2, VarPi0::EA_FIT_EONLY,  blkCoord2, blkPhiCoarse2);
            p2v[static_cast<int>(VarPi0::EA_MIX)]    = buildPhotonTLV(clus2, VarPi0::EA_MIX,    blkCoord2, blkPhiCoarse2);
            p2v[static_cast<int>(VarPi0::PDC_RAW)]   = buildPhotonTLV(clus2, VarPi0::PDC_RAW,   blkCoord2, blkPhiCoarse2);
            p2v[static_cast<int>(VarPi0::PDC_CORR)]  = buildPhotonTLV(clus2, VarPi0::PDC_CORR,  blkCoord2, blkPhiCoarse2);

            for (int v = 0; v < static_cast<int>(VarPi0::NVAR); ++v) {
              const double m = (p1[v] + p2v[v]).M();
              if (std::isfinite(m) && h_m_pi0_var[v][i]) h_m_pi0_var[v][i]->Fill(m);
            }

        }
      }


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
    PHCompositeNode* /*topNode*/,
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
    float weight,
    bool  fillGlobal)
{
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

  if (!m_firstPassBvaluesOnly && h_nclusters)
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
    if (!m_firstPassBvaluesOnly && h2_chi2_tot_etaPhi)
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
        if (!m_firstPassBvaluesOnly)
        {
          if (h2_chi2_rej_etaPhi)  h2_chi2_rej_etaPhi->Fill(lt_eta, lt_phi);
          if (p_chi2_pass_etaPhi)  p_chi2_pass_etaPhi->Fill(lt_eta, lt_phi, 0.0);
        }
        if (Verbosity() > 0)
          std::cout << "[χ²-CUT]  cluster Chi2=" << clusChi2
                    << "  (lead η,φ)=" << lt_eta << ',' << lt_phi
                    << "  REJECTED\n";
        continue;
      }
      else
      {
        if (!m_firstPassBvaluesOnly && p_chi2_pass_etaPhi)
          p_chi2_pass_etaPhi->Fill(lt_eta, lt_phi, 1.0);
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
    if (!m_firstPassBvaluesOnly && h_clusE)
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
    if (!m_firstPassBvaluesOnly && h_clusE_nTow)
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
    int   blkPhiCoarse = -1;             // will be filled by getBlockCord
    int blkEtaCoarse = -1;
    auto blkCoord =
            getBlockCord(towerEtas, towerPhis, towerEs,
                         blkPhiCoarse, blkEtaCoarse);

    /* 4) find energy slice ........................................... */
    const int iEbin = getEnergySlice( clusE );
      
      {
          const auto* c2 = dynamic_cast<const RawClusterv2*>(clus1);
          float x_raw_v2 = std::numeric_limits<float>::quiet_NaN();
          float y_raw_v2 = std::numeric_limits<float>::quiet_NaN();
          float x_cor_v2 = std::numeric_limits<float>::quiet_NaN();
          float y_cor_v2 = std::numeric_limits<float>::quiet_NaN();
          if (c2) {
              x_raw_v2 = c2->x_tower_raw();
              y_raw_v2 = c2->y_tower_raw();
              x_cor_v2 = c2->x_tower_corr();
              y_cor_v2 = c2->y_tower_corr();
          }
          
          // recompute CoG (raw) for a sanity check
          float x_cog = std::numeric_limits<float>::quiet_NaN();
          float y_cog = std::numeric_limits<float>::quiet_NaN();
          bool  got_cog = false;
          if (m_bemcRec && m_geometry)
          {
              std::vector<EmcModule> hitlist;
              hitlist.reserve(std::distance(clus1->get_towers().first, clus1->get_towers().second));
              const int Nx = m_bemcRec->GetNx();
              
              auto range = clus1->get_towers();
              for (auto it = range.first; it != range.second; ++it)
              {
                  const auto tg = m_geometry->get_tower_geometry(it->first);
                  if (!tg) continue;
                  EmcModule m;
                  m.ich = tg->get_bineta() * Nx + tg->get_binphi();
                  m.amp = static_cast<float>(it->second);
                  m.tof = 0.0f;
                  hitlist.push_back(m);
              }
              if (!hitlist.empty())
              {
                  float E=0, px=0, py=0, pxx=0, pyy=0, pyx=0;
                  m_bemcRec->Momenta(&hitlist, E, px, py, pxx, pyy, pyx, 0.0f);
                  if (E > 0.0f) { x_cog = px; y_cog = py; got_cog = true; }
              }
          }
          
          // reuse your existing histos to track v2 − recalc
          if (c2 && got_cog)
          {
              if (h_dx_prop_vsE) h_dx_prop_vsE->Fill(clusE, x_raw_v2 - x_cog);
              if (h_dy_prop_vsE) h_dy_prop_vsE->Fill(clusE, y_raw_v2 - y_cog);
          }
          
          if (m_print_first_N_clusters-- > 0)
          {
              std::cout << ANSI_BOLD << ANSI_GREEN
              << "[V2] clusID=" << clus1->get_id()
              << "  E=" << clusE << " GeV"
              << ANSI_RESET << "\n";
              std::cout << "  raw(v2):    (" << x_raw_v2  << "," << y_raw_v2  << ")  "
              << (c2 ? "[OK]" : "[MISSING]") << "\n";
              if (got_cog)
                  std::cout << "  raw(recalc):(" << x_cog << "," << y_cog << ")\n";
              if (c2 && got_cog)
                  std::cout << "  Δ(raw v2 − recalc): (" << (x_raw_v2 - x_cog) << ","
                  << (y_raw_v2 - y_cog) << ")\n";
              std::cout << "  corr(v2):   (" << x_cor_v2  << "," << y_cor_v2  << ")  "
              << (c2 ? "[OK]" : "[MISSING]") << "\n";
          }
      }


    /* optional detailed print ........................................ */
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] finalClusterLoop →"
                  << "  blockCord=(" << blkCoord.first << ',' << blkCoord.second << ")"
                  << "  |  coarse (η,φ)=(" << blkEtaCoarse << ',' << blkPhiCoarse << ")"
                  << "  |  iEbin=" << iEbin
                  << "  |  E=" << clusE << "  pT=" << clusPt
                  << "  |  #towers=" << towerEs.size()
                  << '\n';
    }


      constexpr double kFillW = 1.0;      // weight per fill
      if (fillGlobal)                       // store only when |z_vtx| ≤ 10 cm
      {
        // 1) Legacy default uncorrected (unchanged)
        h3_cluster_block_cord_E->Fill(blkCoord.first,
                                      blkCoord.second,
                                      clusE,
                                      kFillW);

        // 2) Uncorrected η-slices in parallel
        const float absEta = std::fabs(clusEta);

        if (h3_cluster_block_cord_E_full && absEta <= 1.10f)
          h3_cluster_block_cord_E_full->Fill(blkCoord.first, blkCoord.second, clusE, kFillW);

        if (h3_cluster_block_cord_E_etaCore && absEta <= 0.20f)
          h3_cluster_block_cord_E_etaCore->Fill(blkCoord.first, blkCoord.second, clusE, kFillW);

        if (h3_cluster_block_cord_E_etaMid && absEta > 0.20f && absEta <= 0.70f)
          h3_cluster_block_cord_E_etaMid->Fill(blkCoord.first, blkCoord.second, clusE, kFillW);

        if (h3_cluster_block_cord_E_etaEdge && absEta > 0.70f && absEta <= 1.10f)
          h3_cluster_block_cord_E_etaEdge->Fill(blkCoord.first, blkCoord.second, clusE, kFillW);
      }

      if (Verbosity() > 0)
      {
        const int bx = h3_cluster_block_cord_E->GetXaxis()->FindBin(blkCoord.first);
        const int by = h3_cluster_block_cord_E->GetYaxis()->FindBin(blkCoord.second);
        const int bz = h3_cluster_block_cord_E->GetZaxis()->FindBin(clusE);

        std::cout << "[RAW-FILL]  ηloc="  << blkCoord.first
                  << "  φloc="           << blkCoord.second
                  << "  E="              << clusE
                  << "  (binX,Y,Z = "    << bx << "," << by << "," << bz << ")  "
                  << "→ new content = "
                  << h3_cluster_block_cord_E->GetBinContent(bx,by,bz)
                  << '\n';
      }

      // First-pass mode: fill ONLY uncorrected TH3s, skip everything else
      if (m_firstPassBvaluesOnly)
        continue;

      // ---- raw local coordinates
      const float rawEta = blkCoord.first;    // [0,1]
      const float rawPhi = blkCoord.second;   // [0,1]

      // ---- legacy default corrected (kept)
      float corrEta = rawEta;
      float corrPhi = rawPhi;

      if (isFitDoneForPhi && iEbin >= 0 && iEbin < N_Ebins)
      {
        const float bPhi = m_bValsPhi[iEbin];
        if (bPhi > 1e-9f) corrPhi = doPhiBlockCorr(rawPhi, bPhi);
      }
      if (isFitDoneForEta && iEbin >= 0 && iEbin < N_Ebins)
      {
        const float bEta = m_bValsEta[iEbin];
        if (bEta > 1e-9f) corrEta = doEtaBlockCorr(rawEta, bEta);
      }
      if (fillGlobal)
        h3_cluster_block_cord_E_corrected->Fill(corrEta, corrPhi, clusE, kFillW);

      // ---- per-η-view corrected (second pass; fill only if that view has b-values)
      if (fillGlobal)
      {
        const float absEta = std::fabs(clusEta);

        auto applyCorrView = [&](int vIdx, float inEta, float inPhi, float& outEta, float& outPhi)
        {
          outEta = inEta; outPhi = inPhi;
          if (iEbin < 0 || iEbin >= N_Ebins) return;
          if (m_bPhiReady_view[vIdx]) {
            const float bP = m_bValsPhi_view[vIdx][iEbin];
            if (bP > 1e-9f) outPhi = doPhiBlockCorr(outPhi, bP);
          }
          if (m_bEtaReady_view[vIdx]) {
            const float bE = m_bValsEta_view[vIdx][iEbin];
            if (bE > 1e-9f) outEta = doEtaBlockCorr(outEta, bE);
          }
        };

        // fullEta
        if (h3_cluster_block_cord_E_full_corr && absEta <= 1.10f &&
            (m_bPhiReady_view[0] || m_bEtaReady_view[0]))
        {
          float e, p; applyCorrView(0, rawEta, rawPhi, e, p);
          h3_cluster_block_cord_E_full_corr->Fill(e, p, clusE, kFillW);
        }
        // etaCore
        if (h3_cluster_block_cord_E_etaCore_corr && absEta <= 0.20f &&
            (m_bPhiReady_view[1] || m_bEtaReady_view[1]))
        {
          float e, p; applyCorrView(1, rawEta, rawPhi, e, p);
          h3_cluster_block_cord_E_etaCore_corr->Fill(e, p, clusE, kFillW);
        }
        // etaMid
        if (h3_cluster_block_cord_E_etaMid_corr && absEta > 0.20f && absEta <= 0.70f &&
            (m_bPhiReady_view[2] || m_bEtaReady_view[2]))
        {
          float e, p; applyCorrView(2, rawEta, rawPhi, e, p);
          h3_cluster_block_cord_E_etaMid_corr->Fill(e, p, clusE, kFillW);
        }
        // etaEdge
        if (h3_cluster_block_cord_E_etaEdge_corr && absEta > 0.70f && absEta <= 1.10f &&
            (m_bPhiReady_view[3] || m_bEtaReady_view[3]))
        {
          float e, p; applyCorrView(3, rawEta, rawPhi, e, p);
          h3_cluster_block_cord_E_etaEdge_corr->Fill(e, p, clusE, kFillW);
        }
     }



    if (Verbosity() > 3)
    {
      const int bxC = h3_cluster_block_cord_E_corrected
                        ->GetXaxis()->FindBin(corrEta);
      const int byC = h3_cluster_block_cord_E_corrected
                        ->GetYaxis()->FindBin(corrPhi);
      const int bzC = h3_cluster_block_cord_E_corrected
                        ->GetZaxis()->FindBin(clusE);

      const double rawCnt =
          h3_cluster_block_cord_E->GetBinContent(bxC, byC, bzC);
      const double corCnt =
          h3_cluster_block_cord_E_corrected->GetBinContent(bxC, byC, bzC);

      std::cout << "[CORR-FILL] raw(η,φ)=(" << rawEta  << ',' << rawPhi  << ")  "
                << "→ corr(η,φ)=("          << corrEta << ',' << corrPhi << ")  "
                << "E=" << clusE
                << "  (binX,Y,Z = " << bxC << ',' << byC << ',' << bzC << ")\n"
                << "             rawHist now has " << rawCnt
                << "  |  corrHist now has "        << corCnt << '\n';

      /*  global entry-count check (once per event) */
      static Long64_t prevRawEntries = 0, prevCorEntries = 0;
      const Long64_t totRaw = h3_cluster_block_cord_E            ->GetEntries();
      const Long64_t totCor = h3_cluster_block_cord_E_corrected  ->GetEntries();
      if (totRaw != prevRawEntries || totCor != prevCorEntries)
      {
        std::cout << "[SUMMARY]  total entries – raw: " << totRaw
                  << " | corrected: " << totCor
                  << "  (Δ = " << (totRaw - totCor) << ")\n";
        prevRawEntries = totRaw;
        prevCorEntries = totCor;
      }
    }



    // QA: fill pT vs leadTower, fill cluster-level eta/phi
    h_pt_eta->Fill(clusPt, lt_eta);
    h_pt_eta_rw->Fill(clusPt, lt_eta, weight);
    h_etaphi_clus->Fill(clusEta, clusPhi);

    for (auto& trPhoton : truth_photons)
    {
        const float dR    = photon1.DeltaR(trPhoton);          // ΔR(reco,truth)
        const float ratio = photon1.E() / trPhoton.E();        // Ereco / Etruth

        if (dR   > 0.03f)                    continue;
        if (ratio < 0.30f || ratio > 1.30f)  continue;

        const float dPhi = TVector2::Phi_mpi_pi(
                               photon1.Phi() - trPhoton.Phi() );
        h_delR_recTrth->Fill(dR);

        if (Verbosity() > 0)
        {
            std::cout << "[DEBUG] => cluster1 E=" << photon1.E()
                      << " matched to truthE="    << trPhoton.E()
                      << "  dR="   << dR
                      << "  dPhi=" << dPhi
                      << "  ratio="<< ratio << '\n';
        }
        h_matched_res   ->Fill(ratio, photon1.Eta());
        h_res_e         ->Fill(ratio, photon1.E());
        h_res           ->Fill(ratio);

        int iLTeta = lt_eta, iLTphi = lt_phi;

        h_res_e_eta     ->Fill(ratio, trPhoton.E(), iLTeta);
        h_res_e_phi     ->Fill(ratio, trPhoton.E(), iLTphi);
        h_delEta_e_eta  ->Fill(photon1.Eta()-trPhoton.Eta(), trPhoton.E(), iLTeta);
        h_delR_e_eta    ->Fill(dR, trPhoton.E(), iLTeta);
        h_delPhi_e_eta  ->Fill(dPhi, trPhoton.E(), iLTeta);
        h_delPhi_e_phi  ->Fill(dPhi, trPhoton.E(), iLTphi);
        h_truthE        ->Fill(trPhoton.E());

        fillDPhiAllVariants(
            clus1,
            photon1,                       // reco photon
            trPhoton,                      // truth photon
            blkCoord,                      // (ηloc , φloc)
            blkPhiCoarse,                  // coarse φ index
            vtx_z,                         // vertex‑z (only wrapped through)
            h_phi_diff_cpRaw_E,            // Δφ(CLUS‑RAW)  → truth
            h_phi_diff_cpCorr_E,           // Δφ(CLUS‑CP)   → truth
            fillGlobal                     // NEW: apply 10 cm cut inside
        );
        
        fillDEtaAllVariants(
            clus1,
            photon1,
            trPhoton,
            blkCoord,
            blkEtaCoarse,
            blkPhiCoarse,
            vtx_z,
            h_eta_diff_cpRaw_E,            // Δη(CLUS‑RAW)  → truth
            h_eta_diff_cpCorr_E,           // Δη(CLUS‑CP)   → truth
            fillGlobal
        );
        fillAshLogDx(clus1, photon1, trPhoton,
                       vtx_z,
                       blkCoord, blkPhiCoarse,
                       towerPhis, towerEs);

        if (blkEtaCoarse >= 0 && blkEtaCoarse < NBinsBlock &&
              blkPhiCoarse >= 0 && blkPhiCoarse < NBinsBlock)
        {
              h_res_block_E[blkEtaCoarse][blkPhiCoarse]
                  ->Fill(ratio, trPhoton.E());
        }
    }
      match1 = false;
      TLorentzVector ph1_trEtaPhi(0,0,0,0);

      // Prefer exact π0 daughter match (mother pid=111)
      bool matched_pi0 = false;
      for (const auto& tp : m_truth_pi0_photons)
      {
        const float dR    = photon1.DeltaR(tp.p4);
        const float ratio = photon1.E() / tp.p4.E();
        if (dR < 0.03f && ratio > 0.7f && ratio < 1.5f)
        {
          ph1_trEtaPhi.SetPtEtaPhiE( clusE / TMath::CosH(tp.p4.Eta()),
                                     tp.p4.Eta(), tp.p4.Phi(), clusE );
          match1 = true;
          matched_pi0 = true;
          if (Verbosity() > 0)
          {
            std::cout << "[DEBUG] => π0 daughter match (clus1) e=" << ph1_trEtaPhi.E()
                      << "  eta=" << ph1_trEtaPhi.Eta() << "\n";
          }
          break;
        }
      }
      // Fallback: any meson photon (old behavior) if π0 match not found
      if (!matched_pi0)
      {
        for (const auto& trPhot : truth_meson_photons)
        {
          const float dR    = photon1.DeltaR(trPhot);
          const float ratio = photon1.E() / trPhot.E();
          if (dR < 0.03f && ratio > 0.7f && ratio < 1.5f)
          {
            ph1_trEtaPhi.SetPtEtaPhiE( clusE / TMath::CosH(trPhot.Eta()),
                                       trPhot.Eta(), trPhot.Phi(), clusE );
            match1 = true;
            if (Verbosity() > 0)
            {
              std::cout << "[DEBUG] => meson-decay match (fallback), e=" << ph1_trEtaPhi.E()
                        << "  eta=" << ph1_trEtaPhi.Eta() << "\n";
            }
            break;
          }
        }
      }

    if (Verbosity() > 2)
    {
      std::cout << "[DEBUG] => Starting INNER loop for cluster pairs.\n";
    }

    if (fillGlobal)
        processClusterPairs(clusterContainer,
                                cIt1,
                                vertex,
                                photon1,
                                clusE,
                                lt_eta, lt_phi,
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
                                m_isSimulation);



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
  fillTruthInfo(topNode, vtx_z, truth_photons, truth_meson_photons);
  h2_truthReco_vz->Fill(truth_vz, vtx_z);

  if (Verbosity() > 0)
  {
      std::cout << "[process_towers] => fillTruthInfo done.\n"
                << "    #truth_photons=" << truth_photons.size()
                << ", #truth_meson_photons=" << truth_meson_photons.size()
                << std::endl;
  }

  /* ------------------------------------------------------------------
   * Vertex‑Z acceptance
   *   – discard events only when |vz| exceeds the widest histogram bin
   *   – remember whether “physics” (tight‑|z|) histograms are allowed
  * ------------------------------------------------------------------ */
  if (std::fabs(vtx_z) > m_vzSliceMax)
  {
      if (Verbosity() > 1)
        std::cout << "[process_towers] |vz|=" << vtx_z
                  << " cm > " << m_vzSliceMax
                  << " cm – event skipped.\n";
      return Fun4AllReturnCodes::EVENT_OK;
  }
  const bool allowGlobal = (std::fabs(vtx_z) <= m_vzTightCut);
    
  m_geometry = checkTowerGeometry(topNode);

  if (Verbosity() > 0 && m_geometry)
  {
      std::cout << ANSI_BOLD
                << "[process_towers] => Tower geometry retrieved successfully."
                << ANSI_RESET << '\n';
  }

  // ════════════════════════════════════════════════════════════════════
  // Measure the rigid φ–offset between the ideal grid and the detector
  // ════════════════════════════════════════════════════════════════════
  if (!m_hasOffset && m_geometry)
  {
      /* ------------------------------------------------------------------
       * 1. Geometry summary (useful once per job)
       * ----------------------------------------------------------------*/
      const int nEtaBins = m_geometry->get_etabins();   // usually 96
      const int nPhiBins = m_geometry->get_phibins();   // **128 in sPHENIX**
      const float kRadPerBin = 2.F * M_PI / static_cast<float>(nPhiBins);

      if (Verbosity() > 0)
      {
        std::cout << ANSI_BOLD << "[φ-anchor] Geometry summary\n" << ANSI_RESET
                  << "      nEtaBins  = " << nEtaBins  << '\n'
                  << "      nPhiBins  = " << nPhiBins  << '\n'
                  << "      Δφ(bin)   = " << kRadPerBin << " rad\n";
      }

      /* ------------------------------------------------------------------
       * 2. Circular mean of φ-residuals over *all* towers
       * ----------------------------------------------------------------*/
      double sumSin = 0.0, sumCos = 0.0;
      std::size_t nTowers = 0;

      const auto [beg, end] = m_geometry->get_tower_geometries();
      for (auto it = beg; it != end; ++it)
      {
        const RawTowerGeom* tg = it->second;
        if (!tg) continue;

        const int   iphi = tg->get_binphi();
        const float phiIdeal = (iphi + 0.5F)*kRadPerBin;
        const float delta    = TVector2::Phi_mpi_pi( tg->get_phi() - phiIdeal );
        sumSin += std::sin(delta);
        sumCos += std::cos(delta);
        ++nTowers;

        if (Verbosity() > 1 && nTowers <= 3)
          std::cout << "      tower#" << nTowers
                    << "  iphi=" << iphi
                    << "  φ(real)="  << tg->get_phi()
                    << "  φ(ideal)=" << phiIdeal
                    << "  Δφ="       << TVector2::Phi_mpi_pi(delta) << '\n';
      }

      if (nTowers == 0)
      {
        std::cerr << " [φ-anchor]   WARNING: geometry container is empty – skip\n";
        return Fun4AllReturnCodes::EVENT_OK;
      }

      /* ------------------------------------------------------------------
       * 3. Mean offset  (result in (–π, +π])
       * ----------------------------------------------------------------*/
      m_phi0Offset = std::atan2(sumSin, sumCos);   // circular mean
      m_hasOffset  = true;

      if (Verbosity() > 0)
      {
        std::cout << ANSI_BOLD << ANSI_CYAN
                  << "[φ-anchor]  measured rigid barrel tilt  Δφ₀ = "
                  << m_phi0Offset << "  rad"
                  << ANSI_RESET  << "\n"
                  << "      towers used = " << nTowers << '\n'
                  << "      ⟨sinΔφ⟩ = " << sumSin / nTowers
                  << " , ⟨cosΔφ⟩ = " << sumCos / nTowers << '\n';
      }
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
                     weight,
                     allowGlobal);

  // Clear vectors for tower bin indices
  ht_phi.clear();
  ht_eta.clear();

  if (Verbosity() > 0)
  {
    std::cout << "[process_towers] => Completed finalClusterLoop, "
              << "clearing tower vectors now." << std::endl;
  }

  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD << "[process_towers] END" << ANSI_RESET
              << "\n" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}



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
  {
    std::cerr << "[ERROR] outfile pointer is null – nothing written!\n";
  }

  /* Skip all End()-summary printing in first-pass b-values-only mode */
  if (m_firstPassBvaluesOnly)
  {
    return 0;
  }

  /* ──────────────────────────────────────── helper: % formatting ───── */
  auto pct = [](std::uint64_t n, std::uint64_t d)->std::string
             { return d ? Form("%6.1f %%", 100.*double(n)/double(d))
                         : "   n/a "; };

    // ----- Δφ – CORRECTED 6-WAY (RAW / CP / EA variants) -----
    {
      const std::uint64_t nTot6 =
        m_phi6WayWin_RAW + m_phi6WayWin_CP + m_phi6WayWin_EA_geom +
        m_phi6WayWin_EA_fitEta + m_phi6WayWin_EA_fitE + m_phi6WayWin_EA_mix;

      if (Verbosity() > 0 && nTot6)
      {
        std::cout << '\n' << ANSI_BOLD
                  << "╭──────────────────────────────────────────────╮\n"
                  << "│   Δφ  –  CORRECTED 6-WAY (RAW / CP / EA*)     │\n"
                  << "╰──────────────────────────────────────────────╯" << ANSI_RESET << '\n'
                  << std::left << std::setw(26) << "variant"
                  << std::right<< std::setw(12)<< "wins"
                  << std::setw(12)             << "share\n"
                  << "──────────────────────────────────────────────────────\n"
                  << std::left << std::setw(26) << "CLUS-RAW"
                  << std::right<< std::setw(12)<< m_phi6WayWin_RAW
                  << std::setw(12)             << pct(m_phi6WayWin_RAW, nTot6) << '\n'
                  << std::left << std::setw(26) << "CLUS-CP"
                  << std::right<< std::setw(12)<< m_phi6WayWin_CP
                  << std::setw(12)             << pct(m_phi6WayWin_CP, nTot6) << '\n'
                  << std::left << std::setw(26) << "EA (geometry)"
                  << std::right<< std::setw(12)<< m_phi6WayWin_EA_geom
                  << std::setw(12)             << pct(m_phi6WayWin_EA_geom, nTot6) << '\n'
                  << std::left << std::setw(26) << "EA (|eta|+E fits)"
                  << std::right<< std::setw(12)<< m_phi6WayWin_EA_fitEta
                  << std::setw(12)             << pct(m_phi6WayWin_EA_fitEta, nTot6) << '\n'
                  << std::left << std::setw(26) << "EA (E-only fits)"
                  << std::right<< std::setw(12)<< m_phi6WayWin_EA_fitE
                  << std::setw(12)             << pct(m_phi6WayWin_EA_fitE, nTot6) << '\n'
                  << std::left << std::setw(26) << "EA (φ:E-only, η:|η|+E)"
                  << std::right<< std::setw(12)<< m_phi6WayWin_EA_mix
                  << std::setw(12)             << pct(m_phi6WayWin_EA_mix, nTot6) << '\n';
      }
    }


    // ----- Δη – CORRECTED 6-WAY (RAW / CP / EA variants) -----
    {
      const std::uint64_t nTot6 =
        m_eta6WayWin_RAW + m_eta6WayWin_CP + m_eta6WayWin_EA_geom +
        m_eta6WayWin_EA_fitEta + m_eta6WayWin_EA_fitE + m_eta6WayWin_EA_mix;

      if (Verbosity() > 0 && nTot6)
      {
        std::cout << '\n' << ANSI_BOLD
                  << "╭──────────────────────────────────────────────╮\n"
                  << "│   Δη  –  CORRECTED 6-WAY (RAW / CP / EA*)     │\n"
                  << "╰──────────────────────────────────────────────╯" << ANSI_RESET << '\n'
                  << std::left << std::setw(26) << "variant"
                  << std::right<< std::setw(12)<< "wins"
                  << std::setw(12)             << "share\n"
                  << "──────────────────────────────────────────────────────\n"
                  << std::left << std::setw(26) << "CLUS-RAW"
                  << std::right<< std::setw(12)<< m_eta6WayWin_RAW
                  << std::setw(12)             << pct(m_eta6WayWin_RAW, nTot6) << '\n'
                  << std::left << std::setw(26) << "CLUS-CP"
                  << std::right<< std::setw(12)<< m_eta6WayWin_CP
                  << std::setw(12)             << pct(m_eta6WayWin_CP, nTot6) << '\n'
                  << std::left << std::setw(26) << "EA (geometry)"
                  << std::right<< std::setw(12)<< m_eta6WayWin_EA_geom
                  << std::setw(12)             << pct(m_eta6WayWin_EA_geom, nTot6) << '\n'
                  << std::left << std::setw(26) << "EA (|eta|+E fits)"
                  << std::right<< std::setw(12)<< m_eta6WayWin_EA_fitEta
                  << std::setw(12)             << pct(m_eta6WayWin_EA_fitEta, nTot6) << '\n'
                  << std::left << std::setw(26) << "EA (E-only fits)"
                  << std::right<< std::setw(12)<< m_eta6WayWin_EA_fitE
                  << std::setw(12)             << pct(m_eta6WayWin_EA_fitE, nTot6) << '\n'
                  << std::left << std::setw(26) << "EA (φ:E-only, η:|η|+E)"
                  << std::right<< std::setw(12)<< m_eta6WayWin_EA_mix
                  << std::setw(12)             << pct(m_eta6WayWin_EA_mix, nTot6) << '\n';
      }
    }

    /* ───────────────────────── 3) Energy-sliced 6-way (RAW / CP / EA×4) ─────────── */
    if (Verbosity() > 5)
    {
      std::cout << '\n' << ANSI_BOLD
                << "╭──────────────────────────────────────────────────────────╮\n"
                << "│   ENERGY-SLICED  Δφ & Δη  (RAW / CP / EA_geom / EA_|η|+E / EA_E-only / EA_mix) │\n"
                << "╰──────────────────────────────────────────────────────────╯"
                << ANSI_RESET << '\n';

      for (int ie = 0; ie < N_Ebins; ++ie)
      {
        // φ slice totals (6-way)
        const std::uint64_t phiRAW   = m_phi6WayWinByE_RAW[ie];
        const std::uint64_t phiCP    = m_phi6WayWinByE_CP[ie];
        const std::uint64_t phiGeom  = m_phi6WayWinByE_EA_geom[ie];
        const std::uint64_t phiFitEta= m_phi6WayWinByE_EA_fitEta[ie];
        const std::uint64_t phiFitE  = m_phi6WayWinByE_EA_fitE[ie];
        const std::uint64_t phiMix   = m_phi6WayWinByE_EA_mix[ie];
        const std::uint64_t phiTOT   = phiRAW + phiCP + phiGeom + phiFitEta + phiFitE + phiMix;

        // η slice totals (6-way)
        const std::uint64_t etaRAW   = m_eta6WayWinByE_RAW[ie];
        const std::uint64_t etaCP    = m_eta6WayWinByE_CP[ie];
        const std::uint64_t etaGeom  = m_eta6WayWinByE_EA_geom[ie];
        const std::uint64_t etaFitEta= m_eta6WayWinByE_EA_fitEta[ie];
        const std::uint64_t etaFitE  = m_eta6WayWinByE_EA_fitE[ie];
        const std::uint64_t etaMix   = m_eta6WayWinByE_EA_mix[ie];
        const std::uint64_t etaTOT   = etaRAW + etaCP + etaGeom + etaFitEta + etaFitE + etaMix;

        std::cout << ANSI_BOLD << "  • E-slice " << ie << ANSI_RESET
                  << "   (entries: Δφ=" << phiTOT << ", Δη=" << etaTOT << ")\n";

        std::cout << std::left << std::setw(26) << "variant"
                  << std::right<< std::setw(12)<< "Δφ wins"
                  << std::setw(10)             << "%"
                  << std::setw(12)             << "Δη wins"
                  << std::setw(10)             << "%\n"
                  << "  ───────────────────────────────────────────────────────\n";

        auto P = [&](std::uint64_t n, std::uint64_t d){ return pct(n, d); };

        std::cout << std::left << std::setw(26) << "CLUS-RAW"
                  << std::right<< std::setw(12)<< phiRAW
                  << std::setw(10)             << P(phiRAW, phiTOT)
                  << std::setw(12)             << etaRAW
                  << std::setw(10)             << P(etaRAW, etaTOT) << '\n';

        std::cout << std::left << std::setw(26) << "CLUS-CP"
                  << std::right<< std::setw(12)<< phiCP
                  << std::setw(10)             << P(phiCP, phiTOT)
                  << std::setw(12)             << etaCP
                  << std::setw(10)             << P(etaCP, etaTOT) << '\n';

        std::cout << std::left << std::setw(26) << "EA (geometry)"
                  << std::right<< std::setw(12)<< phiGeom
                  << std::setw(10)             << P(phiGeom, phiTOT)
                  << std::setw(12)             << etaGeom
                  << std::setw(10)             << P(etaGeom, etaTOT) << '\n';

        std::cout << std::left << std::setw(26) << "EA (|eta|+E fits)"
                  << std::right<< std::setw(12)<< phiFitEta
                  << std::setw(10)             << P(phiFitEta, phiTOT)
                  << std::setw(12)             << etaFitEta
                  << std::setw(10)             << P(etaFitEta, etaTOT) << '\n';

        std::cout << std::left << std::setw(26) << "EA (E-only fits)"
                  << std::right<< std::setw(12)<< phiFitE
                  << std::setw(10)             << P(phiFitE, phiTOT)
                  << std::setw(12)             << etaFitE
                  << std::setw(10)             << P(etaFitE, etaTOT) << '\n';

        std::cout << std::left << std::setw(26) << "EA (φ:E-only, η:|η|+E)"
                  << std::right<< std::setw(12)<< phiMix
                  << std::setw(10)             << P(phiMix, phiTOT)
                  << std::setw(12)             << etaMix
                  << std::setw(10)             << P(etaMix, etaTOT) << '\n';

        std::cout << "  ───────────────────────────────────────────────────────\n";
      }
    }

  return 0;
}



// ============================================================================
//  getBlockCord – returns the cluster position *inside* one 2×2 coarse block
// ---------------------------------------------------------------------------
//  • Output local coordinates are in (−0.5 … +1.5]  *before* any further
//    calibration (this matches the convention used by downstream code).
//  • A single symmetric fold is applied if |loc| > 1.5 or ≤ −0.5;
//    the associated coarse index is advanced **once** so that the caller
//    sees a self-consistent (blk,loc) pair.
//  • blkPhi is returned via the extra out-parameter `blockPhiBinOut` –
//    *always* the *possibly-updated* value, so that the caller never works
//    with stale information.
// ============================================================================
std::pair<float, float>
PositionDependentCorrection::getBlockCord(const std::vector<int>   &towerEtas,
                                          const std::vector<int>   &towerPhis,
                                          const std::vector<float> &towerEs,
                                          int                      &blockPhiBinOut,
                                          int                      &blockEtaBinOut)
{
  [[maybe_unused]] constexpr int kNEtaFine        =  96;   // mapping only
  [[maybe_unused]] constexpr int kNPhiFine        = 256;
  constexpr int      kFinePerBlock   = 2;                 // 2×2 super-cell
  constexpr int      kNCoarseBlocks  = kNPhiFine / 2;     // 128

  const std::size_t nT = towerEs.size();

  if (Verbosity() > 0)
    std::cout << ANSI_CYAN << "[getBlockCord] ENTER  nTowers = "
              << nT << ANSI_RESET << '\n';

  if (nT == 0 ||
      towerEtas.size() != nT ||
      towerPhis.size() != nT)
  {
    if (Verbosity() > 0)
      std::cout << ANSI_RED
                << "  ✘ empty vectors or size mismatch – return {0,0}\n"
                << ANSI_RESET;
    blockPhiBinOut = 0;
    return {0.0f, 0.0f};
  }

  /* ─── 1) energy-weighted centroids in *fine* bins ────────────────────── */
  const float etaCoG = getAvgEta(towerEtas, towerEs);   //  0 … <96
  const float phiCoG = getAvgPhi(towerPhis, towerEs);   //  0 … <256

  if (Verbosity() > 1)
    std::cout << "    ⟨η⟩_tower=" << etaCoG
              << "  ⟨φ⟩_tower=" << phiCoG << '\n';

  /* ─── 2) parent coarse-block indices ─────────────────────────────────── */
  int blkEta = static_cast<int>(std::floor(etaCoG)) / kFinePerBlock;   // 0 … 47
  int blkPhi = static_cast<int>(std::floor(phiCoG)) / kFinePerBlock;   // 0 … 127

  /* ─── 3) raw local coords (still in fine-bin units) ──────────────────── */
  float locEta = etaCoG - blkEta * kFinePerBlock;
  float locPhi = phiCoG - blkPhi * kFinePerBlock;

  auto foldOnce = [&](float &loc, int &coarse, const char *tag)
    {
      if (loc <= -0.5f || loc > 1.5f)
      {
        const float before = loc;

        // Symmetric one-fold with matching coarse index update
        if (loc <= -0.5f) { loc += 2.0f; --coarse; }
        else              { loc -= 2.0f; ++coarse; }

        // Wrap coarse index both ways
        if (coarse < 0)                    coarse = kNCoarseBlocks - 1;
        else if (coarse >= kNCoarseBlocks) coarse = 0;

        if (Verbosity() > 0)
          std::cout << ANSI_YELLOW
                    << "    • " << tag
                    << " folded: " << before << " → " << loc
                    << "  |  blk → " << coarse
                    << ANSI_RESET << '\n';
      }
  };
    
//auto foldOnce = [&](float &loc, int &coarse, const char *tag)
//{
//  if (loc <= -0.5f || loc > 1.5f)
//  {
//    const float before = loc;
//    loc = std::fmod(loc + 2.0f, 2.0f);           // (0 … 2)
//    if (loc > 1.5f) { loc -= 2.0f; ++coarse; }   // shift block by +1
//
//    if (coarse == kNCoarseBlocks) coarse = 0;    // φ wrap-around
//
//    if (Verbosity() > 0)
//      std::cout << ANSI_YELLOW
//                << "    • " << tag
//                << " folded: " << before << " → " << loc
//                << "  |  blk+1 → " << coarse
//                << ANSI_RESET << '\n';
//  }
//};

  foldOnce(locEta, blkEta, "η");
  foldOnce(locPhi, blkPhi, "φ");

  /* ─── 5) propagate the *updated* block index and return ─────────────── */
  blockPhiBinOut = blkPhi;
  blockEtaBinOut = blkEta;
              
  if (Verbosity() > 0)
    std::cout << ANSI_CYAN << "[getBlockCord] EXIT\n" << ANSI_RESET << '\n';

  return {locEta, locPhi};
}
// ============================================================================
//  getAvgEta – energy–weighted η–centre of a cluster (no wrapping necessary)
//  ----------------------------------------------------------------------------
//  INPUT  towerEtas       int  vector, raw η–row of every tower hit (0 … 95)
//         towerEnergies   float vector, same length, tower energy in GeV
//
//  OUTPUT float           fractional η–index in the 96-row grid
//                         0.0 ≤  avgη  ≤ 95.0
// ============================================================================
float PositionDependentCorrection::getAvgEta(const std::vector<int>   &towerEtas,
                                             const std::vector<float> &towerEnergies)
{
  const int nTowers = towerEtas.size();

  if (Verbosity() > 0)
    std::cout << ANSI_CYAN << "[getAvgEta] called with nTowers = "
              << nTowers << ANSI_RESET << '\n';

  if (nTowers == 0 || nTowers != static_cast<int>(towerEnergies.size()))
  {
    if (Verbosity() > 0)
      std::cout << ANSI_RED << "  ERROR – empty or size-mismatched vectors!"
                << ANSI_RESET << std::endl;
    return 0.f;
  }

  // Single-tower cluster – fast exit
  if (nTowers == 1) return static_cast<float>( towerEtas[0] );

  // -------------------------------------------------------------------------
  // 1) Energy-weighted sum
  // -------------------------------------------------------------------------
  double sumE      = 0.0;   // ΣE
  double sumEeta   = 0.0;   // Σ(E · η)

  for (int i = 0; i < nTowers; ++i)
  {
    const double E   = towerEnergies[i];
    const int    eta = towerEtas[i];          // 0 … 95

    sumE    += E;
    sumEeta += E * eta;

    if (Verbosity() > 1)
      std::cout << "    tower[" << i << "]  η=" << eta
                << "  E=" << E << '\n';
  }

  if (sumE < 1e-9) return 0.f;                // guard against ΣE ≈ 0

  // -------------------------------------------------------------------------
  // 2) Energy-weighted centroid
  // -------------------------------------------------------------------------
  const double avgEta = sumEeta / sumE;       // guaranteed 0.0 … 95.0

  if (Verbosity() > 0)
  {
    std::cout << ANSI_GREEN << "  RESULT ⟨η⟩ = " << avgEta
              << ANSI_RESET << '\n'
              << ANSI_CYAN  << "[getAvgEta] END"
              << ANSI_RESET << "\n";
  }
  return static_cast<float>(avgEta);
}

// ============================================================================
//  getAvgPhi – energy-weighted φ–centroid  (0 ≤ ⟨φ⟩ < 256)
// ----------------------------------------------------------------------------
float
PositionDependentCorrection::getAvgPhi(const std::vector<int>   &towerPhis,
                                       const std::vector<float> &towerEs)
{
  constexpr int   kNPhiFine  = 256;                // total fine bins
  constexpr float kHalfSpan  = kNPhiFine / 2.0f;   // 128.0
  const std::size_t nT       = towerPhis.size();

  if (Verbosity() > 0)
    std::cout << ANSI_CYAN << "[getAvgPhi] ENTER | nTowers = "
              << nT << ANSI_RESET << '\n';

  if (nT == 0 || towerEs.size() != nT)
  {
    if (Verbosity() > 0)
      std::cout << ANSI_RED
                << "  ✘ empty vectors or size mismatch – return 0\n"
                << ANSI_RESET;
    return 0.0f;
  }

  /* ─── 1) choose reference (highest-energy tower) ─────────────────────── */
  const std::size_t iRef =
      static_cast<std::size_t>(
          std::distance(towerEs.begin(),
                        std::max_element(towerEs.begin(), towerEs.end())));
  const int refBin = towerPhis[iRef];              // 0 … 255

  if (Verbosity() > 1)
    std::cout << "    reference φ-bin (max-E) = " << refBin << '\n';

  /* ─── 2) accumulate energy and moment after *one* wrap ───────────────── */
  double sumE = 0.0;
  double sumEphi = 0.0;
  int    nWrapLow = 0, nWrapHigh = 0;

  for (std::size_t i = 0; i < nT; ++i)
  {
    float E   = towerEs[i];
    int   phi = towerPhis[i];                      // raw fine bin

    int d = phi - refBin;                          // distance to reference

    if (d < -kHalfSpan) { phi += kNPhiFine; ++nWrapLow;  }
    if (d >  kHalfSpan) { phi -= kNPhiFine; ++nWrapHigh; }

    sumE     += E;
    sumEphi  += static_cast<double>(E) * phi;

    if (Verbosity() > 2)
      std::cout << "      tower#" << std::setw(2) << i
                << "  bin(raw)=" << std::setw(3) << towerPhis[i]
                << "  bin(adj)=" << std::setw(3) << phi
                << "  E=" << E << '\n';
  }

  if (sumE < 1e-12)
  {
    if (Verbosity() > 0)
      std::cout << ANSI_RED << "  ✘ ΣE ≈ 0 – return 0\n" << ANSI_RESET;
    return 0.0f;
  }

  /* ─── 3) energy-weighted mean and fold to canonical domain ───────────── */
  double phiAvg = sumEphi / sumE;                  // may be negative/large
  phiAvg = std::fmod(phiAvg + kNPhiFine, kNPhiFine); // 0 … <256

  if ((nWrapLow + nWrapHigh) && Verbosity() > 0)
    std::cout << ANSI_YELLOW
              << "    • wrap corrections: " << nWrapLow
              << " low , " << nWrapHigh << " high\n"
              << ANSI_RESET;

  if (Verbosity() > 0)
  {
    std::cout << ANSI_GREEN
              << "    RETURN ⟨φ⟩ = " << phiAvg << "  (fine-bin units)\n"
              << ANSI_RESET
              << ANSI_CYAN << "[getAvgPhi] EXIT\n" << ANSI_RESET << '\n';
  }
  return static_cast<float>(phiAvg);
}


/*==========================================================================*
 *  φ-direction                                                             *
 *==========================================================================*/
float PositionDependentCorrection::doPhiBlockCorr(float localPhi, float b)
{
  /* —————————————————— banner —————————————————— */
  if (Verbosity() > 0)
    std::cout << ANSI_BOLD << "[doPhiBlockCorr] ENTER  "
              << "localPhi = " << localPhi << " ,  b = " << b
              << ANSI_RESET << '\n';

  /* (0) nothing to undo when |b| ≈ 0 */
  if (std::fabs(b) < 1e-9f) {                        // quick exit
    if (Verbosity() > 0) std::cout << "  b≈0 → return unchanged\n";
    return localPhi;
  }

  /* (1) single fold into (-0.5 … +1.5] */
  if (localPhi <= -0.5f || localPhi > 1.5f) {
    if (Verbosity() > 1) std::cout << "  fold: " << localPhi << " → ";
    localPhi = std::fmod(localPhi + 2.f, 2.f);
    if (Verbosity() > 1) std::cout << localPhi << '\n';
  }

  /* (2) map to (-0.5 … +0.5] */
  const float Xmeas = (localPhi < 0.5f) ? localPhi : localPhi - 1.f;
  if (Verbosity() > 1) std::cout << "  Xmeas = " << Xmeas << '\n';

  /* (3) analytic inverse  ———  asinh! */
  const double S = std::sinh(1.0 / (2.0 * b));          // sinh(1/2b)
  if (Verbosity() > 1) std::cout << "  S = sinh(1/2b) = " << S << '\n';

  const float Xtrue = static_cast<float>( b * std::asinh( 2.0 * S * Xmeas ) );
  if (Verbosity() > 1) std::cout << "  Xtrue = " << Xtrue << '\n';

  /* (4) shift back to (-0.5 … +1.5] */
  float corrected = (localPhi < 0.5f) ? Xtrue : Xtrue + 1.f;
  if (Verbosity() > 1)
    std::cout << "  corrected (pre-wrap) = " << corrected << '\n';

  /* (5) safety fold */
  if (corrected <= -0.5f || corrected > 1.5f) {
    corrected = std::fmod(corrected + 2.f, 2.f);
    if (corrected > 1.5f) corrected -= 2.f;
    if (Verbosity() > 1)
      std::cout << "  corrected (post-wrap) = " << corrected << '\n';
  }

  if (Verbosity() > 0)
    std::cout << ANSI_GREEN << "[doPhiBlockCorr] EXIT  → "
              << corrected << " (rad-fraction)"
              << ANSI_RESET << "\n\n";

  return corrected;
}

/*==========================================================================*
 *  η-direction                                *
 *==========================================================================*/
float PositionDependentCorrection::doEtaBlockCorr(float  localEta,
                                                  float  b,
                                                  float  dEtaTower /* = 1.f */)
{
  if (Verbosity() > 0)
    std::cout << ANSI_BOLD << "[doEtaBlockCorr] ENTER  "
              << "ηloc=" << localEta << "  b=" << b
              << ANSI_RESET << '\n';

  (void) dEtaTower;
  if (std::fabs(b) < 1e-9f)                         // nothing to undo
  {
    if (Verbosity() > 0) std::cout << "  b≈0 → return unchanged\n";
    return localEta;
  }
  /* ------------------------------------------------------------------ *
   * 1) fold once into (‑0.5 … +1.5]                                    *
   * ------------------------------------------------------------------ */
  if (localEta <= -0.5f || localEta > 1.5f)
    localEta = std::fmod(localEta + 2.f, 2.f);

  if (!std::isfinite(localEta))
  {
    if (Verbosity() > 0)
      std::cout << ANSI_RED << "  ηloc not finite → NaN\n" << ANSI_RESET;
    return std::numeric_limits<float>::quiet_NaN();
  }
  /* ------------------------------------------------------------------ *
   * 2) map measured coordinate to centre‑of‑gravity frame              *
   * ------------------------------------------------------------------ */
  const float Xmeas = (localEta < 0.5f) ? localEta : localEta - 1.f;

  /* ------------------------------------------------------------------ *
   * 3) analytic inverse of forward Ash‑distortion                      *
   *    (identical to φ)                                                *
   * ------------------------------------------------------------------ */
  const double  S     = std::sinh( 1.0 / (2.0 * b) );
  const float   Xtrue = static_cast<float>( b * std::asinh( 2.0 * S * Xmeas ) );

  /* ------------------------------------------------------------------ *
   * 4) add back integer tower pitch                                    *
   * ------------------------------------------------------------------ */
  float corrected = (localEta < 0.5f) ? Xtrue : Xtrue + 1.f;

  /* ------------------------------------------------------------------ *
   * 5) safety fold & physical edge guard                               *
   * ------------------------------------------------------------------ */
  if (corrected <= -0.5f || corrected > 1.5f)
  {
    corrected = std::fmod(corrected + 2.f, 2.f);
    if (corrected > 1.5f) corrected -= 2.f;
  }
//  /* If a caller supplied a non‑unity tower pitch (rare) apply it now */
//  corrected *= dEtaTower;

  /* keep clusters inside the barrel acceptance (η ≈ ±1.1) */
  if (std::fabs(corrected) > 1.6f)      // 2×2 block edge + margin
    corrected = std::copysign(1.6f, corrected);

  if (Verbosity() > 0)
    std::cout << ANSI_GREEN << "[doEtaBlockCorr] EXIT → "
              << corrected << ANSI_RESET << "\n\n";

  return corrected;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  convertBlockToGlobalPhi – (blk,loc) → absolute EMCAL φ in (−π , +π]
// ----------------------------------------------------------------------------
//  • Accepts any real local coordinate; performs **exactly one** fold into
//    the canonical domain (−0.5 … +1.5] without ever touching the coarse index.
//  • Returns NaN if the coarse index is illegal or if the folded local
//    coordinate is still out of bounds (should never happen).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float
PositionDependentCorrection::convertBlockToGlobalPhi(int   block_phi_bin,
                                                     float localPhi)
{
  constexpr int   kFinePerBlock   = 2;                       // 2 fine bins
  constexpr int   kNCoarseBlocks  = 128;                     // 256 / 2
  constexpr int   kNTowerBins     = 256;                     // total fine
  constexpr float kRadPerBin      = 2.0f * static_cast<float>(M_PI) / kNTowerBins;

  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD << ANSI_CYAN
              << "[convertBlockToGlobalPhi] ENTER\n"
              << "    ▸ block_phi_bin = " << block_phi_bin
              << "  ,  localPhi(raw) = " << localPhi
              << ANSI_RESET << '\n';
  }

  /* ─── sanity on coarse index ─────────────────────────────────────────── */
  if (block_phi_bin < 0 || block_phi_bin >= kNCoarseBlocks)
  {
    if (Verbosity() > 0)
      std::cout << ANSI_RED
                << "  ✘ coarse index outside 0…" << kNCoarseBlocks-1
                << "   → return NaN\n" << ANSI_RESET;
    return std::numeric_limits<float>::quiet_NaN();
  }

  /* ─── single, non-destructive fold of the local coordinate ───────────── */
  if (localPhi <= -0.5f || localPhi > 1.5f)
  {
    localPhi = std::fmod(localPhi + 2.0f, 2.0f);   // (0 … 2)
    if (Verbosity() > 0)
      std::cout << ANSI_YELLOW
                << "    • localPhi out of band – folded once to "
                << localPhi << ANSI_RESET << '\n';
  }

  if (!std::isfinite(localPhi) ||
      localPhi <= -0.5f || localPhi > 1.5f)
  {
    if (Verbosity() > 0)
      std::cout << ANSI_RED
                << "  ✘ localPhi still invalid – return NaN\n"
                << ANSI_RESET;
    return std::numeric_limits<float>::quiet_NaN();
  }

  /* ─── fine-bin centre of the requested tower ─────────────────────────── */
  float fullPhiIndex =
        static_cast<float>(block_phi_bin) * kFinePerBlock +  // coarse origin
        localPhi +                                           // intra-block
        0.5f;                                                // centre of fine bin

  if (Verbosity() > 1)
    std::cout << "    ▸ fullPhiIndex = " << fullPhiIndex
              << "   (coarse*2 + local + 0.5)\n";
    
//  /* --- residual φ‑bias inside an 8‑tower module (same as CorrectPosition) --- */
//  {
//      int  mod8   = static_cast<int>(fullPhiIndex) / 8;           // module #
//      float x8    = fullPhiIndex - mod8 * 8.0f - 4.0f;            // −4 … +4
//      float dBias = 0.10f * x8 / 4.0f;                            // default factors
//      if (std::fabs(x8) > 3.3f) dBias = 0.0f;                     // edge guard
//      fullPhiIndex -= dBias;                                      // apply bias‑undo
//  }


  /* ─── linear → angular & wrap to (−π , +π] ──────────────────────────── */
  float globalPhi = fullPhiIndex * kRadPerBin;  // ideal grid
  if (m_hasOffset) globalPhi += m_phi0Offset;   // shift to real detector
  globalPhi = TVector2::Phi_mpi_pi(globalPhi);

  if (Verbosity() > 0)
  {
    std::cout << ANSI_GREEN
              << "    ▸ OUTPUT  global φ = " << globalPhi << " rad\n"
              << ANSI_RESET
              << ANSI_BOLD << ANSI_CYAN
              << "[convertBlockToGlobalPhi] EXIT\n"
              << ANSI_RESET << '\n';
  }
  return globalPhi;
}


/**********************************************************************
 * convertBlockToGlobalEta – diagnostic edition
 * ---------------------------------------------------------------
 * • Maps a (coarse-η block, local-η) pair to the tower-centre η
 *   in the lab frame, with one-fold logic identical to the
 *   production routine.
 * • Every decision point prints a line when Verbosity() ≥ 4.
 * • On *any* fatal condition the function returns NaN and
 *   prints a highlighted error banner (Verbosity() ≥ 1).
 *********************************************************************/
float
PositionDependentCorrection::convertBlockToGlobalEta(int   block_eta_bin,
                                                     float localEta)  // local ∈ (-∞,∞)
{
  /* ───────────── compile-time constants ───────────── */
  constexpr int   kFinePerBlock   = 2;                 // 2 fine η-towers / coarse block
  constexpr int   kNFineEtaBins   = 96;                // full barrel
  constexpr int   kNCoarseBlocks  = kNFineEtaBins / kFinePerBlock;           // 48
  constexpr float kEtaMin         = -1.1f;             // η centre of iyFine = 0
  constexpr float kDEtaPerFine    = 2.2f / 96.0f;      // ≈ 0.0229167 per fine tower

  const unsigned vb = Verbosity();

  auto printHdr = [&](const char* tag)
  {
    if (vb > 3)
      std::cout << ANSI_BOLD << ANSI_CYAN
                << "[η-X-MAP] " << tag << ANSI_RESET << '\n';
  };

  printHdr("ENTER");

  /* ================================================================
   * 0)  Parameter sanity
   * ============================================================ */
  if (block_eta_bin < 0 || block_eta_bin >= kNCoarseBlocks)
  {
    if (vb)
      std::cerr << ANSI_RED
                << "[η-X-MAP]  ERROR: coarse-η index (" << block_eta_bin
                << ") outside [0," << kNCoarseBlocks-1 << "] – returning NaN\n"
                << ANSI_RESET;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (!std::isfinite(localEta))
  {
    if (vb)
      std::cerr << ANSI_RED
                << "[η-X-MAP]  ERROR: localEta is NaN/Inf – returning NaN\n"
                << ANSI_RESET;
    return std::numeric_limits<float>::quiet_NaN();
  }

  /* ================================================================
   * 1)  Single symmetric fold of the *local* coordinate
   * ============================================================ */
  if (vb > 4)
    std::cout << "  » raw localEta = " << localEta << '\n';

  if (localEta <= -0.5f || localEta > 1.5f)
  {
    const float before = localEta;
    localEta = std::fmod(localEta + 2.0f, 2.0f);          // (0 … 2]
    if (localEta > 1.5f) localEta -= 2.0f;                // (-0.5 … +1.5]

    if (vb > 3)
      std::cout << ANSI_YELLOW
                << "  • one-fold applied: " << before << " → " << localEta
                << ANSI_RESET << '\n';
  }

  if (localEta <= -0.5f || localEta > 1.5f)
  {
    if (vb)
      std::cerr << ANSI_RED
                << "[η-X-MAP]  ERROR: localEta still outside (-0.5,1.5] after fold – NaN\n"
                << ANSI_RESET;
    return std::numeric_limits<float>::quiet_NaN();
  }

  /* ================================================================
   * 2)  Determine the **fine** tower index (integer 0…95)
   * ============================================================ */
  const int iyFine =
        block_eta_bin * kFinePerBlock +
        static_cast<int>(std::floor(localEta + 0.5f));     // nearest fine tower

  if (vb > 4)
    std::cout << "  » iyFine (tower row) = " << iyFine << '\n';

  if (iyFine < 0 || iyFine >= kNFineEtaBins)
  {
    if (vb)
      std::cerr << ANSI_RED
                << "[η-X-MAP]  ERROR: computed iyFine=" << iyFine
                << " outside [0," << kNFineEtaBins-1 << "] – NaN\n"
                << ANSI_RESET;
    return std::numeric_limits<float>::quiet_NaN();
  }

  /* ================================================================
   * 3)  Try to get tower-centre geometry (detailed mode)
   * ============================================================ */
  const RawTowerDefs::keytype key =
        RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, iyFine, /*iphi=*/0);

  const RawTowerGeom* tg = m_geometry ? m_geometry->get_tower_geometry(key) : nullptr;

  if (!tg)
  {
    if (vb > 2)
      std::cout << ANSI_RED
                << "  • geometry not available – falling back to linear η grid\n"
                << ANSI_RESET;
    const float fullEtaIndex =
          static_cast<float>(block_eta_bin) * kFinePerBlock + localEta + 0.5f;
    const float etaFallback = kEtaMin + fullEtaIndex * kDEtaPerFine;

    if (vb > 2)
      std::cout << "  » fallback η  = " << etaFallback << '\n';

    printHdr("EXIT (fallback)");
    return etaFallback;
  }

  /* ================================================================
   * 4)  Projective-aware interpolation *inside* the physical tower
   * ============================================================ */
  const double rC   = std::hypot(tg->get_center_x(), tg->get_center_y());
  const double zC   = tg->get_center_z();
  const double etaC = std::asinh(zC / rC);

  if (vb > 4)
    std::cout << "  » tower-centre ηC = " << etaC << '\n';

  /* signed offset: (+) towards higher iy, (-) towards lower iy        */
  const double offsetFine =
        localEta - ((localEta < 0.5f) ? 0.0 : 1.0);   // now in (-0.5, +0.5]

  /* —— measure actual η-pitch from neighbouring rows ———————— */
  const auto fetchEta = [&](int iy) -> double
  {
    const auto* t = m_geometry->get_tower_geometry(
          RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, iy, 0));
    return t ? std::asinh(t->get_center_z() /
                          std::hypot(t->get_center_x(), t->get_center_y()))
             : std::numeric_limits<double>::quiet_NaN();
  };

  double etaUp = fetchEta(std::min(iyFine + 1, kNFineEtaBins - 1));
  double etaDn = fetchEta(std::max(iyFine - 1, 0));

  double etaPitch = (std::isfinite(etaUp) && std::isfinite(etaDn))
                      ? 0.5 * (etaUp - etaDn)
                      : kDEtaPerFine;

  if (vb > 4)
    std::cout << "  » η-pitch = " << etaPitch << '\n';

  /* ================================================================
   * 5)  Final global η
   * ============================================================ */
  const double globalEta = etaC + offsetFine * etaPitch;

  if (vb > 2)
    std::cout << ANSI_GREEN
              << "  » global η = " << globalEta << ANSI_RESET << '\n';

  printHdr("EXIT (OK)");
  return static_cast<float>(globalEta);
}
