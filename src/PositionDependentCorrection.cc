#include "PositionDependentCorrection.h"

// G4Hits includes
#include <fstream>
#include <TLorentzVector.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <TKey.h>
// G4Cells includes
#include <array>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <regex>
// Tower includes
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

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// MBD
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
#include <utility>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include"TRandom3.h"
#include <algorithm>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree
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

R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)

using namespace std;

PositionDependentCorrection::PositionDependentCorrection(const std::string &name,
                                                         const std::string &filename)
  : SubsysReco(name)
  , detector("HCALIN")
  , outfilename(filename)
  , g4hitntuple(nullptr)
  , g4cellntuple(nullptr)
  , towerntuple(nullptr)
  , clusterntuple(nullptr)
{
  _eventcounter = 0;
  _vz = 10.0;;
}

PositionDependentCorrection::~PositionDependentCorrection()
{
  delete hm;

  if (g4hitntuple)   delete g4hitntuple;
  if (g4cellntuple)  delete g4cellntuple;
  if (towerntuple)   delete towerntuple;
  if (clusterntuple) delete clusterntuple;
}

int PositionDependentCorrection::Init(PHCompositeNode*)
{
  // ------------------------------------------------------------------------------------
  // Enhanced Debug/Checks for initialization
  // ------------------------------------------------------------------------------------
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] PositionDependentCorrection::Init() called. "
              << "Preparing to initialize histograms and managers." << std::endl;
  }

  // ------------------------------------------------------------------------------------
  // Create a Fun4AllHistoManager to handle and keep track of all histograms created
  // in this module. We check if it is allocated properly.
  // ------------------------------------------------------------------------------------
  hm = new Fun4AllHistoManager(Name());
  if (!hm)
  {
    if (Verbosity() > 0)
    {
      std::cerr << "[ERROR] Fun4AllHistoManager allocation failed (hm is null). "
                << "Cannot proceed with initialization." << std::endl;
    }
    // Return an error code if needed
    return -1;
  }

  // ------------------------------------------------------------------------------------
  // Open the output ROOT file, which will store the histograms and other results
  // when the job completes. We also ensure it's opened correctly.
  // ------------------------------------------------------------------------------------
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  std::cout << "[INFO-OUTFILENAME] Writing histograms to: " << outfilename << std::endl;
  if (!outfile || outfile->IsZombie())
  {
    if (Verbosity() > 0)
    {
      std::cerr << "[ERROR] Could not open output file '" << outfilename
                << "' for writing (outfile is null or zombie)." << std::endl;
    }
    // Return an error code if needed
    return -1;
  }

  // ------------------------------------------------------------------------------------
  // Initialize the TriggerAnalyzer pointer. Ensuring we have a valid pointer for
  // any further use in the module.
  // ------------------------------------------------------------------------------------
  trigAna = new TriggerAnalyzer();
  if (!trigAna)
  {
    if (Verbosity() > 0)
    {
      std::cerr << "[ERROR] TriggerAnalyzer pointer is null. "
                << "Cannot proceed with certain functionalities." << std::endl;
    }
    // We continue, but this may lead to limited functionality later.
  }

  // ===========================
  //  CORRELATION PLOTS
  // ===========================
  //
  // h_mass_eta_lt, h_mass_eta_lt_rw:
  //   Each is a 2D histogram that plots the pi0 candidate mass on the x-axis and the tower-eta index
  //   (lead tower in the cluster) on the y-axis.
  //   - h_mass_eta_lt:    unweighted or "raw" correlation of (mass vs. lead tower eta).
  //   - h_mass_eta_lt_rw: the same correlation but with pT-based reweighting applied.
  //   They are used to study how the measured diphoton mass depends on the calorimeter's eta index,
  //   which is crucial for calibrations and position‐dependent corrections.
  h_mass_eta_lt = new TH2F("h_mass_eta_lt", "", 50, 0, 0.5, 96, 0, 96);
  h_mass_eta_lt_rw = new TH2F("h_mass_eta_lt_rw", "", 50, 0, 0.5, 96, 0, 96);

  // h_pt_eta, h_pt_eta_rw:
  //   Two 2D histograms plotting cluster pT on the x-axis and the tower-eta index on the y-axis.
  //   - h_pt_eta:    raw distribution of (pT vs. eta index).
  //   - h_pt_eta_rw: the same distribution but with a reweighting factor.
  //   This helps in checking how often certain pT ranges occur in different eta indices,
  //   which can be related to position corrections or acceptance effects.
  h_pt_eta = new TH2F("h_pt_eta", "", 100, 0, 10, 96, 0, 96);
  h_pt_eta_rw = new TH2F("h_pt_eta_rw", "", 100, 0, 10, 96, 0, 96);

  // h_cemc_etaphi:
  //   A 2D occupancy histogram of tower coordinates in the CEMC (eta index vs. phi index).
  //   Useful for checking the spatial distribution of towers above threshold
  //   (which can highlight dead/hot areas, or show position-dependent coverage).
  h_cemc_etaphi = new TH2F("h_cemc_etaphi", "", 96, 0, 96, 256, 0, 256);

  // ===========================
  //  1D DISTRIBUTIONS
  // ===========================
  //
  // h_InvMass, h_InvMass_w, h_InvMassMix:
  //   1D histograms for the invariant mass of two-cluster systems (candidates for pi0).
  //   - h_InvMass:      unweighted mass distribution.
  //   - h_InvMass_w:    weighted mass distribution, typically used if we apply pT reweighting.
  //   - h_InvMassMix:   used for mixed-event background or combinatorial background studies
  //                     (depending on your method, e.g. to subtract random pair backgrounds).
  //   These are central to pi0 calibration and resolution checks as a function of position.
  h_InvMass = new TH1F("h_InvMass", "Invariant Mass", 500, 0, 1.0);
  h_InvMass_w = new TH1F("h_InvMass_w", "Invariant Mass", 500, 0, 1.0);
  h_InvMassMix = new TH1F("h_InvMassMix", "Invariant Mass", 120, 0, 1.2);

  // h_tower_e:
  //   1D histogram storing the energy of individual towers. Used to get the overall tower energy
  //   distribution, to see how often towers are firing at certain energy, and for QA.
  h_tower_e = new TH1F("h_tower_e","",1000,-1,5);

  // ===========================
  //  CLUSTER QA
  // ===========================
  //
  // h_etaphi_clus:
  //   2D histogram of cluster pseudorapidity vs. cluster phi (in radians).
  //   Good for visualizing cluster distributions in the detector (e.g. identifying hot or dead regions).
  h_etaphi_clus = new TH2F("h_etaphi_clus", "", 140, -1.2, 1.2, 64, -1 * TMath::Pi(), TMath::Pi());

  // h_clusE:
  //   1D histogram of the cluster total energy. Helps examine the energy spectrum of reconstructed clusters.
  h_clusE = new TH1F("h_clusE", "", 100, 0, 10);

  // h_clusE_nTow:
  //   2D histogram correlating the cluster energy (x-axis) with the number of towers in that cluster (y-axis).
  //   Useful QA to see how many towers typically form a cluster of a given energy and for checking clustering performance.
  h_clusE_nTow = new TH2F("h_clusE_nTow","",20,0,20,50,0,50);

  // h_emcal_e_eta:
  //   1D histogram that accumulates tower energy sums in each EMCal eta index.
  //   Good for verifying uniformity or identifying large-scale non-uniformities across eta.
  h_emcal_e_eta = new TH1F("h_emcal_e_eta", "", 96, 0, 96);

  // h_pt1, h_pt2:
  //   1D histograms for the transverse momentum of the first and second photons (clusters) in a pair,
  //   typically used when reconstructing pi0s. They help track the leading vs. subleading cluster pT distribution.
  h_pt1 = new TH1F("h_pt1", "", 100, 0, 5);
  h_pt2 = new TH1F("h_pt2", "", 100, 0, 5);

  // h_nclusters:
  //   1D histogram counting the total number of clusters in each event.
  //   Good for monitoring occupancy and changes in cluster multiplicity with event conditions.
  h_nclusters = new TH1F("h_nclusters", "", 100, 0, 100);

  // ===========================
  //  TRUTH / MATCHING QA
  // ===========================
  //
  // h_truth_eta, h_truth_e, h_truth_pt:
  //   1D histograms for truth-level primary particles (especially photons):
  //   - h_truth_eta: distribution in pseudorapidity,
  //   - h_truth_e:   distribution in energy,
  //   - h_truth_pt:  distribution in transverse momentum.
  //   All used to compare generator-level (or G4-level) kinematics vs. reconstructed data.
  h_truth_eta = new TH1F("h_truth_eta", "", 100, -1.2, 1.2);
  h_truth_e = new TH1F("h_truth_e", "", 100, 0, 10);
  h_truth_pt = new TH1F("h_truth_pt", "", 100, 0, 10);

  // h_delR_recTrth:
  //   1D histogram of the ∆R separation between a reconstructed cluster and the nearest truth photon.
  //   This helps in studying matching efficiency vs. ∆R and finalizing the matching cut.
  h_delR_recTrth = new TH1F("h_delR_recTrth", "", 500, 0, 5);

  // h_matched_res:
  //   2D histogram plotting the ratio of reconstructed E to truth E (x-axis) vs. the cluster pseudorapidity (y-axis).
  //   This is for seeing how the energy resolution or scaling depends on eta, once a cluster is matched to a truth photon.
  h_matched_res = new TH2F("h_matched_res","",100,0,1.5,20,-1,1);

  // h_res_e:
  //   2D histogram of the reconstructed/truth energy ratio (x-axis) vs. reconstructed cluster energy (y-axis).
  //   Good for analyzing resolution or scale shifts as a function of cluster energy.
  h_res_e = new TH2F("h_res_e","",100,0,1.5,20,0,20);

  // h_res_e_phi, h_res_e_eta:
  //   3D histograms that extend the above resolution studies to dependence on phi and eta indices.
  //   - h_res_e_phi: ratio vs. cluster energy, vs. tower phi bin.
  //   - h_res_e_eta: ratio vs. cluster energy, vs. tower eta bin.
  //   These are key to diagnosing position-dependent energy calibration (which ties in with cluster position corrections).
  h_res_e_phi = new TH3F("h_res_e_phi","",100,0,1.5,10,0,20,256,0,256);
  h_res_e_eta = new TH3F("h_res_e_eta","",300,0,1.5,40,0,20,96,0,96);

  // h_m_pt_eta, h_m_ptTr_eta, h_m_ptTr_eta_trKin:
  //   3D histograms tracking mass, pT, and eta simultaneously, used for pi0 calibration:
  //   - h_m_pt_eta:        (mass vs. cluster pT vs. tower eta).
  //   - h_m_ptTr_eta:      (mass vs. truth photon E vs. tower eta).
  //   - h_m_ptTr_eta_trKin:(mass vs. truth-based pi0 candidate kinematics).
  //   These help isolate how mass reconstruction depends on energy and position in the calorimeter.
  h_m_pt_eta = new TH3F("h_m_pt_eta","",70,0,0.7,10,0,10,96,0,96);
  h_m_ptTr_eta= new TH3F("h_m_ptTr_eta","",70,0,0.7,10,0,10,96,0,96);
  h_m_ptTr_eta_trKin = new TH3F("h_m_ptTr_eta_trKin","",70,0,0.7,10,0,10,96,0,96);

  // h_res:
  //   1D histogram capturing a generic resolution measure, typically the ratio E_reco/E_truth for matched clusters.
  //   Summarizes the distribution of how well cluster energy matches truth.
  h_res = new TH1F("h_res", "", 50, 0, 1.5);

  // h_delEta_e_eta, h_delR_e_eta, h_delPhi_e_eta, h_delPhi_e_phi:
  //   3D histograms capturing the differences in eta, R, phi between reco cluster and truth photon:
  //   - h_delEta_e_eta: (∆eta) vs. truth energy vs. tower eta.
  //   - h_delR_e_eta:   (∆R)   vs. truth energy vs. tower eta.
  //   - h_delPhi_e_eta: (∆phi) vs. truth energy vs. tower eta.
  //   - h_delPhi_e_phi: (∆phi) vs. truth energy vs. tower phi.
  //   All are used to analyze and correct position shifts or smearing as a function of location in the detector.
  h_delEta_e_eta = new TH3F("h_delEta_e_eta","",100,-0.1,0.1,10,0,20,96,0,96);
  h_delR_e_eta = new TH3F("h_delR_e_eta","",100,-0.1,0.1,10,0,20,96,0,96);
  h_delPhi_e_eta = new TH3F("h_delPhi_e_eta","",100,-0.3,0.3,20,0,20,96,0,96);
  h_delPhi_e_phi = new TH3F("h_delPhi_e_phi","",100,-0.1,0.1,20,0,20,256,0,256);

  // pr_eta_shower, pr_phi_shower:
  //   TProfile histograms that can store the average shower response or correction factor vs. tower eta or phi.
  //   - pr_eta_shower: profiles the detector response (or residual) as a function of the tower eta index.
  //   - pr_phi_shower: same, but as a function of the tower phi index.
  //   Profiles are often used in calibration contexts because they store mean values in each bin.
  pr_eta_shower = new TProfile("pr_eta_shower","",96,-48.5,47.5, -1,1.5);
  pr_phi_shower = new TProfile("pr_phi_shower","",256,-128.5,127.5, -1,1.5);

  // h_vert_xy:
  //   2D histogram collecting the x–y positions of the global event vertex. Useful for verifying
  //   the distribution of vertex positions in the plane, which can affect angles/incidences.
  h_vert_xy = new TH2F("h_vert_xy","",500,-120,120,500,-120,120);
    
    // 1D histograms for truth and reco vertex Z
  h_truth_vz = new TH1F("h_truth_vz","Truth Vertex Z;z_{truth} (cm);Counts",200,-100,100);
  hm->registerHisto(h_truth_vz);

  h_reco_vz = new TH1F("h_reco_vz","Reco Vertex Z;z_{reco} (cm);Counts",200,-100,100);
  hm->registerHisto(h_reco_vz);

  // 2D histogram for truth vs reco vertex Z
  h2_truthReco_vz = new TH2F("h2_truthReco_vz","Truth vs Reco Vertex Z;z_{truth} (cm);z_{reco} (cm)",
                               200,-100,100, 200,-100,100);
  hm->registerHisto(h2_truthReco_vz);

  // h_truthE:
  //   1D histogram of the truth-level energies of photons or relevant particles.
  //   With 10,000 bins from 0 to 30, it gives a fine resolution of the truth E distribution.
  h_truthE = new TH1F("h_truthE","",10000,0,30);

  // ===========================
  //  POSITION-DEPENDENT STUFF
  // ===========================
  //
  // The following histograms directly relate to the block (2×2) coordinates within a tower cluster.
  // We want to see how the invariant mass or energy resolution changes as a function of block coordinate,
  // which is crucial for implementing position-dependent corrections inside each tower block.

  // h_mass_block_pt[ie][ip]:
  //   A 2D histogram for each block coordinate bin (ie, ip). The x-axis is the di-photon invariant mass,
  //   and the y-axis is the cluster pT. We want to see how mass depends on local block position
  //   to refine calibrations that vary with position in a tower block.
  // h_res_block_E[ie][ip]:
  //   Similarly, for each block coordinate bin, it records the ratio of reconstructed/truth energy
  //   (or some resolution measure) vs. cluster energy.
  //   This helps see if different positions in a tower lead to different calibration biases.
  for (int ie=0; ie<NBinsBlock; ie++){
    for (int ip=0; ip<NBinsBlock; ip++){
      h_mass_block_pt[ie][ip] = new TH2F(Form("h_mass_block_%d_%d_pt",ie,ip),"",100,0,1,5,0,10);
      h_res_block_E[ie][ip]   = new TH2F(Form("h_res_block_%d_%d_E",ie,ip),"",120,0,1.2,5,0,10);
    }
  }



    // 1) Define edges for X (uniform -0.5..1.5, 14 bins => 15 edges)
    Double_t xEdges[15];
    {
      const double xMin = -0.5;
      const double xMax =  1.5;
      const double step = (xMax - xMin) / 14.0; // =2.0/14=0.142857...
      for(int i=0; i<=14; i++)
        xEdges[i] = xMin + i*step;
    }

    // 2) Define edges for Y (uniform -0.5..1.5, 14 bins => 15 edges)
    Double_t yEdges[15];
    {
      const double yMin = -0.5;
      const double yMax =  1.5;
      const double step = (yMax - yMin) / 14.0; // same as above
      for(int i=0; i<=14; i++)
        yEdges[i] = yMin + i*step;
    }

    // 3) Z bin edges: 2..4..6..8..10..12..15..20..30
    //    You already have this in your eEdge[] array; ensure it's a Double_t[].
    static constexpr Double_t eEdges[9] = {2,4,6,8,10,12,15,20,30};

    // 4) Now create the “uncorrected” TH3F with array-of-edges constructor:
    h3_cluster_block_cord_E = new TH3F(
        "h2_cluster_block_cord_E",                // name
        "Uncorrected local block coords vs. E",   // title
        14, xEdges,                               // X: 14 bins
        14, yEdges,                               // Y: 14 bins
        8,  eEdges                                // Z: 8 bins, from eEdges[]
    );

    // 5) Create the “corrected” TH3F likewise:
    h3_cluster_block_cord_E_corrected = new TH3F(
        "h2_cluster_block_cord_E_corrected",
        "Corrected local block coords vs. E",
        14, xEdges,
        14, yEdges,
        8,  eEdges
    );
  // h_block_bin:
  //   1D histogram that may be used to look at the block coordinate axis in discrete steps (0–1, or subdivided further).
  //   Potentially used to label or count how many clusters end up in each sub-cell bin across the 2×2 block space.
  h_block_bin = new TH1F("h_block_bin","",14,-0.5,1.5);

  // For measuring raw phi resolution
    // (A) Create the raw Δφ histograms for each pT bin
  for (int i = 0; i < N_Ebins; i++)
    {
      float eLo = eEdge[i];
      float eHi = eEdge[i+1];
      h_phi_diff_raw_E[i] = new TH1F(
        Form("h_phi_diff_raw_%.0f_%.0f", eLo, eHi),  // histogram name
        Form("#Delta#phi raw, %.0f < pT < %.0f", eLo, eHi),
        200, -0.1, 0.1
      );
        
      h_eta_diff_raw_E[i] = new TH1F(
          Form("h_eta_diff_raw_%.1f_%.1f", eLo, eHi),
          Form("#Delta#eta raw, %.1f < E < %.1f", eLo, eHi),
          200, -0.1, 0.1
      );
  }
    
    // (B) Create the corrected Δφ histograms for each pT bin
    for (int i = 0; i < N_Ebins; i++)
    {
      // Use your already-declared static member array eEdge[]:
      float eLo = eEdge[i];
      float eHi = eEdge[i+1];

      // Now you have a properly named eHi, so you can use it
      h_phi_diff_corrected_E[i] = new TH1F(
        Form("h_phi_diff_corr_%.1f_%.1f", eLo, eHi),
        Form("#Delta#phi corrected, %.1f < pT < %.1f", eLo, eHi),
        200, -0.1, 0.1
      );
        
      h_eta_diff_corrected_E[i] = new TH1F(
          Form("h_eta_diff_corr_%.1f_%.1f", eLo, eHi),
          Form("#Delta#eta corrected, %.1f < E < %.1f", eLo, eHi),
          200, -0.1, 0.1
      );
    }

    

  // Maybe a TProfile to see how #Delta#phi depends on local coordinate or energy
  pr_phi_vs_blockcoord = new TProfile("pr_phi_vs_blockcoord","",14,-0.5,1.5, -0.2,0.2);
    

    // GOOD: use eEdge[i] and eEdge[i+1]
    for (int i = 0; i < N_Ebins; i++)
    {
      // e.g. eEdge[] = {2.0, 3.0, 5.0, 8.0, 12.0}
      float eLo = eEdge[i];
      float eHi = eEdge[i+1];

      h_localPhi_corrected[i] = new TH1F(
        Form("h_localPhi_corrected_%.1f_%.1f", eLo, eHi),
        Form("Corrected local #phi, %.1f < E < %.1f", eLo, eHi),
        50, -0.5, 0.5
      );
        
     h_localEta_corrected[i] = new TH1F(
          Form("h_localEta_corrected_%.1f_%.1f", eLo, eHi),
          Form("Corrected local #eta, %.1f < E < %.1f", eLo, eHi),
          50, -0.5, 0.5
        );
    }
 
    // ──────────────────────────────────────────────────────────────────────────────
    //  Ash-b   and   Log-w0  trial grids  (energy-binned version, PHENIX-like)
    // ──────────────────────────────────────────────────────────────────────────────
    {
      // Geometry constant – cell size at the EMCal front face (sPHENIX barrel)
      constexpr double cellSize = 5.55;     // [cm]

      // Histogram window for Δx  (= xreco – xtrue)
      constexpr double DX_MAX  = 12.0;      // [cm]  covers the full observed range
      constexpr double BIN_W   = 0.10;      // [cm]  keep 0.1-cm resolution
      constexpr int    NBINS   = int( 2*DX_MAX / BIN_W + 0.5 );   // → 240 bins

      // ────────────────────────────────────────────────────────────────────────
      // 1) Ash scan – b parameter in *cm*  (0.00 … 1.60 cm, 0.05-cm grid)
      //    NOTE: the list here ***MUST*** match the list probed later in doAshScan.
      // ────────────────────────────────────────────────────────────────────────
      m_bScan.clear();
      for (double b_cm = 0.00; b_cm <= 1.60 + 1e-9; b_cm += 0.05)   // ← start at 0.00
      {
        const double b_cell = b_cm / cellSize;                      // dimensionless
        m_bScan.push_back( std::round(b_cell*10000.)/10000. );      // 4-dec precision
      }

      // ────────────────────────────────────────────────────────────────────────
      // 2) Log-weight scan – w0 parameter (1.5 … 7.0 in 0.10 steps)
      // ────────────────────────────────────────────────────────────────────────
      m_w0Scan.clear();
      for (double w0 = 1.5; w0 <= 7.0 + 1e-9; w0 += 0.10)
        m_w0Scan.push_back( std::round(w0*100.)/100. );             // 2-dec precision

      // ────────────────────────────────────────────────────────────────────────
      // 3) Declare the (E-slice dependent) histograms  **ONCE**
      // ────────────────────────────────────────────────────────────────────────
      if (!alreadyDeclaredHistograms)
      {
        for (int iE = 0; iE < N_Ebins; ++iE)          // ← energy-bin loop
        {
          // 3a)  Ash  (b scan)
          for (double bVal : m_bScan)
          {
            const TString hName = Form("h_dx_ash_b%.4f_E%d", bVal, iE);
            auto* h = new TH1F(hName,
                               ";x_{reco}-x_{true}  [cm];Counts / 0.1 cm",
                               NBINS, -DX_MAX, +DX_MAX);             // **expanded axis**
            hm->registerHisto(h);
          }

          // 3b)  Log  (w0 scan)
          for (double w0 : m_w0Scan)
          {
            const TString hName = Form("h_dx_log_w0%.2f_E%d", w0, iE);
            auto* h = new TH1F(hName,
                               ";x_{reco}-x_{true}  [cm];Counts / 0.1 cm",
                               NBINS, -DX_MAX, +DX_MAX);             // **expanded axis**
            hm->registerHisto(h);
          }
        }
        alreadyDeclaredHistograms = true;
      }
    }



  rnd = new TRandom3();
    //--------------------------------------------------------------------
    //  b-value reader   (general-N_Ebins version)
    //--------------------------------------------------------------------
    {
      /*--------------------------------------------------------------*/
      /* 0)  First check whether the user even asked for the fit.      */
      /*--------------------------------------------------------------*/
      if (!isFitDoneForPhi && !isFitDoneForEta)
      {
        std::cout << "[INFO]   isFitDoneForPhi && isFitDoneForEta are both false ⇒ "
                     "skip reading bValues.txt\n";
        /* keep the two flags false and return */
      }
      else
      {
        const std::string bFilePath =
            "/sphenix/u/patsfan753/scratch/PDCrun24pp/bParameters/bValues.txt";

        std::ifstream bfile(bFilePath);
        if (!bfile.is_open())
        {
          std::cout << "[WARN]  bValues.txt NOT found at " << bFilePath
                    << "  ➜  disabling φ- and η-corrections.\n";
          isFitDoneForPhi = false;
          isFitDoneForEta = false;
        }
        else
        {
          std::cout << "[INFO]  Reading b-values from " << bFilePath << '\n';

          /*----------------------------------------------------------*/
          /* 1)  Local book-keeping                                   */
          /*----------------------------------------------------------*/
          std::vector<bool> gotBinPhi(N_Ebins,false);
          std::vector<bool> gotBinEta(N_Ebins,false);

          /*  PHI [2,4) 0.123  or  ETA [10,12) 0.456  */
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
            const float  eLo = std::stof(m[2]);
            const float  eHi = std::stof(m[3]);
            const float  bVal= std::stof(m[4]);

            /* find matching energy slice */
            for (int i=0;i<N_Ebins;++i)
            {
              if (std::fabs(eLo-expectedLo[i])<1e-3 &&
                  std::fabs(eHi-expectedHi[i])<1e-3)
              {
                if (dim=="PHI")
                {
                  m_bValsPhi[i] = bVal;
                  gotBinPhi[i]  = true;
                } else {
                  m_bValsEta[i] = bVal;
                  gotBinEta[i]  = true;
                }
                break;
              }
            }
          } // while getline
          bfile.close();

          /*----------------------------------------------------------*/
          /* 2)  Decide whether we have a complete set                */
          /*----------------------------------------------------------*/
          const bool allPhi = std::all_of(gotBinPhi.begin(),gotBinPhi.end(),
                                          [](bool x){return x;});
          const bool allEta = std::all_of(gotBinEta.begin(),gotBinEta.end(),
                                          [](bool x){return x;});

          isFitDoneForPhi &= allPhi;   // keep previous request, but veto if incomplete
          isFitDoneForEta &= allEta;

          /*----------------------------------------------------------*/
          /* 3)  Print a slice-by-slice summary                       */
          /*----------------------------------------------------------*/
          auto printTable = [&](const char* tag,
                                const std::vector<bool>& got,
                                const float*  bArr,
                                bool enabled)
          {
            std::cout << "[INFO]  " << tag << (enabled ? "  (enabled)\n"
                                                       : "  (DISABLED — missing slices)\n");
            for (int i=0;i<N_Ebins;++i)
            {
              std::cout << "        [" << expectedLo[i] << ',' << expectedHi[i]
                        << ")  : ";
              if (got[i])
                std::cout << bArr[i] << '\n';
              else
                std::cout << "-- missing --\n";
            }
          };

          printTable("PHI b-values",gotBinPhi,m_bValsPhi,isFitDoneForPhi);
          printTable("ETA b-values",gotBinEta,m_bValsEta,isFitDoneForEta);
        } // file opened OK
      }   // end “user asked for correction”
    }

    

  // ------------------------------------------------------------------------------------
  // Final check: if we reach here, we have presumably allocated all histograms successfully.
  // Return success to indicate the initialization is complete.
  // ------------------------------------------------------------------------------------
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] PositionDependentCorrection::Init() completed successfully. "
              << "All histograms and objects have been allocated." << std::endl;
  }

  return 0;
}



int PositionDependentCorrection::process_event(PHCompositeNode* topNode)
{
    // Print a debug statement before incrementing _eventcounter
    std::cout << ANSI_BOLD
    << "\n---------  processing event " << _eventcounter << "  ---------\n"
    << ANSI_RESET << std::endl;
    
    // (optional) keep Verbosity-gated details, or delete them:
    if (Verbosity() > 0)
    {
        std::cout << "[DEBUG] PositionDependentCorrection::process_event() called. "
        << "Current event counter = " << _eventcounter << std::endl;
    }
    
    // bump the counter *after* the headline if you want it to start at 0
    ++_eventcounter;
    
    // ----------------------------------------------------------------
    //  normal processing
    // ----------------------------------------------------------------
    process_towers(topNode);
    return Fun4AllReturnCodes::EVENT_OK;
}

float PositionDependentCorrection::retrieveVertexZ(PHCompositeNode* topNode)
{
  float vtx_z = 0;

  // Print debug statement about retrieveVertexZ logic
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] PositionDependentCorrection::retrieveVertexZ() called. "
              << "Checking if getVtx is true..." << std::endl;
  }

  if (getVtx)
  {
    // Attempt to locate the GlobalVertexMap node
    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (!vertexmap)
    {
      // If it's missing, print a warning
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG] PositionDependentCorrection::retrieveVertexZ() - "
                  << "GlobalVertexMap node is missing!" << std::endl;
      }
    }

    // If the vertexmap exists and is not empty, retrieve the first vertex
    if (vertexmap && !vertexmap->empty())
    {
      GlobalVertex* vtx = vertexmap->begin()->second;
      if (vtx)
      {
        vtx_z = vtx->get_z();
        if (std::isfinite(vtx_z)) h_reco_vz->Fill(vtx_z);        // finally fill
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG] Retrieved vertex Z position: "
                    << vtx_z << std::endl;
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG] GlobalVertex pointer was null, cannot read z position." << std::endl;
        }
      }
    }
    else if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] Either GlobalVertexMap is null or empty; no vertex to retrieve." << std::endl;
    }
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] getVtx is false; not retrieving vertex from GlobalVertexMap." << std::endl;
    }
  }

  // Return the found or defaulted z position
  return vtx_z;
}


// ----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////
// fillTowerInfo(...)
// ------------------
// This method loops over the calorimeter towers in the node tree (named "TOWERINFO_CALIB_CEMC")
// and fills various histograms and QA metrics. It updates the total tower energy in 'tower_tot_e',
// marks any towers flagged as "bad" but having significant energy, and records occupancy above a
// threshold 'emcal_hit_threshold'. It also feeds histograms that track tower energy distributions.
///////////////////////////////////////////////////////////////////////////////
void PositionDependentCorrection::fillTowerInfo(
    PHCompositeNode* topNode,
    float emcal_hit_threshold,
    float& tower_tot_e,
    std::vector<float>& ht_eta,
    std::vector<float>& ht_phi)
{
  // Attempt to retrieve a TowerInfoContainer named "TOWERINFO_CALIB_CEMC" from the node tree.
  // If missing or null, we won't fill these histograms.
  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  if (towers)
  {
    // 'size' is the total number of tower channels that exist in this container.
    int size = towers->size();

    // Loop over each tower channel index.
    for (int channel = 0; channel < size; channel++)
    {
      // Get the TowerInfo object for this channel.
      // TowerInfo contains the measured (or calibrated) energy plus flags for tower status.
      TowerInfo* tower = towers->get_tower_at_channel(channel);

      // offlineenergy is the recorded or reconstructed energy of this tower.
      float offlineenergy = tower->get_energy();

      // Each tower is associated with an encoded key that can be decoded
      // into integer eta/phi indices that label where it sits in the calorimeter map.
      unsigned int towerkey = towers->encode_key(channel);
      int ieta = towers->getTowerEtaBin(towerkey);
      int iphi = towers->getTowerPhiBin(towerkey);

      // Determine if the tower is flagged as good or bad
      // (bad if tower->get_isBadChi2() was set, meaning suspicious or corrupt data).
      bool isGood = !(tower->get_isBadChi2());

      // Add this tower's energy to the running total 'tower_tot_e'.
      // This can be used later for overall QA or event-level checks.
      tower_tot_e += offlineenergy;

      // Fill a 1D histogram of tower energies to observe the global energy distribution across towers.
      h_tower_e->Fill(offlineenergy);

      // If the tower is flagged as 'bad' (isGood == false) but still has >0.2 GeV energy,
      // record its ieta and iphi. This can be useful for diagnosing hot or misbehaving towers.
      if (!isGood && offlineenergy > 0.2)
      {
        ht_eta.push_back(ieta);
        ht_phi.push_back(iphi);
      }

      // If the tower is good, fill a 2D histogram (eta index vs energy),
      // to visualize how energy is distributed over eta for good towers.
      if (isGood)
      {
        h_emcal_e_eta->Fill(ieta, offlineenergy);
      }

      // If the tower's energy is above a certain threshold (emcal_hit_threshold),
      // fill a 2D occupancy histogram (ieta vs iphi). This is a common QA check
      // to see where active towers are in the calorimeter.
      if (offlineenergy > emcal_hit_threshold)
      {
        h_cemc_etaphi->Fill(ieta, iphi);
      }
    }
  }
}

void PositionDependentCorrection::fillTruthInfo(
    PHCompositeNode* topNode,
    float& vtx_z,
    std::vector<TLorentzVector>& truth_photons,
    std::vector<TLorentzVector>& truth_meson_photons)
{
  // 'wieght' and 'weight' are separate float variables
  // which can be used for histogram weighting (user-defined).
  float wieght = 1; // (unchanged logic, original variable name)
  float weight = 1; // (unchanged logic, original variable name)

  // -------------------------------------------------------------------------
  // Verbosity block: introduction
  // -------------------------------------------------------------------------
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

  // If Verbosity is enabled, let the user know we've got a valid pointer
  if (Verbosity() > 0)
  {
    std::cout << "  --> Successfully retrieved 'G4TruthInfo' container.\n"
              << "  --> Now analyzing primary particles..." << std::endl;
  }

  // -------------------------------------------------------------------------
  // 1) Primary particles
  // -------------------------------------------------------------------------
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

    // pT-based weight = pT * exp(-3 * pT)
    weight = myVector.Pt() * TMath::Exp(-3 * myVector.Pt());
    h_truth_pt->Fill(myVector.Pt(), weight);

    // Debug print: PID, E, pT, eta, phi
    if (debug)
    {
      std::cout << ANSI_BOLD << "[Debug]" << ANSI_RESET
                << " Primary pid=" << truth->get_pid()
                << "   E=" << energy
                << "  pt=" << myVector.Pt()
                << "  eta=" << myVector.Eta()
                << "  phi=" << myVector.Phi() << std::endl;
    }

    // Store in truth_photons
    truth_photons.push_back(myVector);
  }

  // If Verbosity is enabled, indicate we've finished analyzing primary particles
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

  // If Verbosity is enabled, indicate we've finished secondary analysis
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

  // -------------------------------------------------------------------------
  // Final summary (Verbosity)
  // -------------------------------------------------------------------------
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



float PositionDependentCorrection::doPhiBlockCorr(float localPhi, float bphi)
{
  // If Verbosity is enabled, print out a brief start message and function parameters
  if (Verbosity() > 0)
  {
    std::cout << "PositionDependentCorrection::doPhiBlockCorr - START" << std::endl;
    std::cout << "  --> localPhi = " << localPhi << ", bphi = " << bphi << std::endl;
    std::cout << "  --> Applying the hyperbolic formula shift in [-0.5, +0.5] domain." << std::endl;
  }

  float Xcg_phi = (localPhi < 0.5) ? localPhi : (localPhi - 1.0);

  TF1 finv("finv", "[0]*TMath::ASinH(2*x*sinh(1/(2*[0])))", -0.5, 0.5);
  finv.SetParameter(0, bphi);
  float t_phi = finv.Eval(Xcg_phi);

  // shift back to ~[0,1]
  float corrected = (localPhi < 0.5) ? t_phi : (t_phi + 1.0);

  if (Verbosity() > 0)
  {
    std::cout << "  --> Intermediate Xcg_phi = " << Xcg_phi << ", t_phi = " << t_phi << std::endl;
    std::cout << "  --> Corrected block phi = " << corrected << std::endl;
    std::cout << "PositionDependentCorrection::doPhiBlockCorr - END" << std::endl << std::endl;
  }

  return corrected;
}

// Same logic for localEta => doEtaBlockCorr(...):
float PositionDependentCorrection::doEtaBlockCorr(float localEta, float bEta)
{
  // SHIFT from [0..1] => [-0.5..+0.5]
  float Xcg_eta = (localEta < 0.5f) ? localEta : (localEta - 1.f);

  // asinh transform with parameter bEta
  TF1 finv("finv_eta", "[0]*TMath::ASinH(2*x*sinh(1/(2*[0])))", -0.5, 0.5);
  finv.SetParameter(0, bEta);
  float t_eta = finv.Eval(Xcg_eta);

  // SHIFT back to [0..1]
  float corrected = (localEta < 0.5f) ? t_eta : (t_eta + 1.f);
  return corrected;
}

// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------

/* ========================================================================
 *  Convert local φ within a 2×2 block  →  absolute calorimeter φ  (radians)
 *  ------------------------------------------------------------------------
 *  • block_phi_bin  : 0 … 127   (coarse 2-tower bin)
 *  • localPhi       : [0,1)     (fractional coordinate in that block)
 *  • result         : φ ∈ (−π,π]
 * ======================================================================== */
float PositionDependentCorrection::convertBlockToGlobalPhi(int   block_phi_bin,
                                                           float localPhi)
{
  constexpr int   kNCoarseBlocks = 128;   // 256 fine bins / 2 towers per block
  constexpr int   kNTowerBins    = 256;   // full EMCal division in φ
  constexpr float kRadPerBin     = 2.f * M_PI / kNTowerBins;

  /*--------------------------------------------------
   * 0)  verbose header
   *--------------------------------------------------*/
  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD << "[convertBlockToGlobalPhi] ENTER" << ANSI_RESET
              << "\n  input  block_phi_bin = " << block_phi_bin
              << "  |  localPhi = " << localPhi << std::endl;
  }

  /*--------------------------------------------------
   * 1)  basic sanity checks
   *--------------------------------------------------*/
  if (Verbosity() > 0)
  {
    if (block_phi_bin < 0 || block_phi_bin >= kNCoarseBlocks)
      std::cout << ANSI_YELLOW << "  [WARN] block_phi_bin outside [0,"
                << kNCoarseBlocks-1 << "] !" << ANSI_RESET << std::endl;

    if (localPhi < 0.f || localPhi >= 1.f)
      std::cout << ANSI_YELLOW << "  [WARN] localPhi outside [0,1) !"
                << ANSI_RESET << std::endl;
  }

  /*--------------------------------------------------
   * 2)  compute global tower index
   *--------------------------------------------------*/
  const float coarseTowerPhi = block_phi_bin * 2.f;   // 2 towers / coarse block
  const float fullPhiIndex   = coarseTowerPhi + localPhi * 2.f;  // 0 … 256

  /*--------------------------------------------------
   * 3)  convert to radians and wrap to (−π,π]
   *--------------------------------------------------*/
  float globalPhi = fullPhiIndex * kRadPerBin;        // 0 … 2π
  globalPhi       = TVector2::Phi_mpi_pi(globalPhi);  // −π … +π

  /*--------------------------------------------------
   * 4)  verbose footer
   *--------------------------------------------------*/
  if (Verbosity() > 0)
  {
    std::cout << "  coarseTowerPhi = " << coarseTowerPhi
              << "  |  fullPhiIndex = " << fullPhiIndex << '\n'
              << "  kRadPerBin      = " << kRadPerBin    << '\n'
              << "  → global φ(rad) = " << globalPhi     << '\n'
              << "[convertBlockToGlobalPhi] EXIT\n" << std::endl;
  }

  return globalPhi;
}


// ----------------------------------------------------------------------------
float PositionDependentCorrection::convertBlockToGlobalEta(int block_eta_bin,
                                                           float localEta)
{
  // Each 2×2 “block” spans two towers in η.
  const float coarseTowerEta  = block_eta_bin * 2.0f;          // tower index (0…95)
  const float fullEtaIndex    = coarseTowerEta + localEta*2.0f; // fractional index

  // Geometry constants for the sPHENIX EMCal barrel
  constexpr float ETA_MIN      = -1.1f;          // lower edge of first tower
  constexpr float TOWER_Deltaη = 2.2f / 96.0f;   // ≈ 0.0229167

  // Map *tower index* → centre of that tower → pseudorapidity
  const float globalEta = ETA_MIN + (fullEtaIndex + 0.5f)*TOWER_Deltaη;

  if (Verbosity() > 0)
  {
    std::cout << "convertBlockToGlobalEta: blockEtaBin=" << block_eta_bin
              << "  localEta=" << localEta
              << "  ⇒ fullEtaIndex=" << fullEtaIndex
              << "  ⇒ η=" << globalEta << std::endl;
  }
  return globalEta;
}



////////////////////////////////////////////////////////////////////////
// doAshShift(...)
//    Minimal "ash" formula that parallels doPhiBlockCorr but is simpler.
//    Takes a local-block φ in [0,1], applies b*asinh, returns ∈ [0,1].
////////////////////////////////////////////////////////////////////////
float PositionDependentCorrection::doAshShift(float localPhi, float bVal)
{
  // 1) Shift localPhi (0..1) -> Xcg in (-0.5..+0.5)
  float Xcg = (localPhi < 0.5f) ? localPhi : (localPhi - 1.0f);

  // 2) asinh transform (with no forced clamp to ±0.5):
  //    Xash = bVal * asinh( 2 * Xcg * sinh(1/(2*bVal)) )
  float val = bVal * asinh( 2.f * Xcg * sinh(1.f/(2.f*bVal)) );

  // 3) Shift back to [0,1] if needed
  if (localPhi >= 0.5f) val += 1.0f;

  // Return the final coordinate
  return val;
}


////////////////////////////////////////////////////////////////////////
// doLogWeightCoord(...)
//    Computes local ϕ in [0,1] via the "log‐weight" formula:
//      w_i = max(0, w0 + ln(E_i / sumE)),
//    and an energy-weighted center in tower phi indices (with wrap fix).
////////////////////////////////////////////////////////////////////////
float PositionDependentCorrection::doLogWeightCoord(const std::vector<int>& towerphis,
                                                    const std::vector<float>& towerenergies,
                                                    float w0)
{
  if(towerphis.empty() || towerenergies.empty()) return 0.0f;

  // 1) sum total tower energy
  double sumE=0;
  for(size_t i=0; i<towerenergies.size(); ++i) sumE += towerenergies[i];
  if(sumE < 1e-9) return 0.0f; // no energy => localPhi=0

  // 2) compute weights w_i, handle phi wrap
  int nphibin = 256;    // for sPHENIX EMC
  double sumW=0, sumWPhi=0;
  auto itMax  = std::max_element(towerenergies.begin(), towerenergies.end());
  int refPhi  = towerphis[ std::distance(towerenergies.begin(), itMax) ];

  for(size_t i=0; i<towerenergies.size(); i++)
  {
    double frac = towerenergies[i] / sumE;
    double w    = w0 + log(frac);
    if(w < 0) w=0; // clamp to zero

    // wrap fix
    int phibin = towerphis[i];
    if(phibin - refPhi < -nphibin/2)  phibin += nphibin;
    else if(phibin - refPhi >  nphibin/2)  phibin -= nphibin;

    sumW    += w;
    sumWPhi += w * phibin;
  }
  if(sumW < 1e-9) return 0.0f;

  // 3) average phi bin => local phi
  double avgphi = sumWPhi / sumW;
  // keep in [0..256)
  while(avgphi < 0)        avgphi += nphibin;
  while(avgphi >= nphibin) avgphi -= nphibin;

  // map tower index to [0,2) => then to [-0.5..+1.5], etc.
  double localPhi = fmod(avgphi + 0.5, 2.0) - 0.5;
  // clamp to [0..1)
  if(localPhi < 0.0)  localPhi += 1.0;
  if(localPhi >=1.0)  localPhi -= 1.0;

  return float(localPhi);
}


//--------------------------------------------------------------------------
//  Fill Δx  (Ash-b & Log-w0 scans) – verbose version
//--------------------------------------------------------------------------
void PositionDependentCorrection::fillAshLogDx(
    RawCluster*                     cluster,
    const TLorentzVector&           recoPhoton,
    const TLorentzVector&           truthPhoton,
    const std::pair<float,float>&   blockCord,   // local (η,φ)
    int                             blockPhiBin,
    const std::vector<int>&         tower_phis,
    const std::vector<float>&       tower_energies )
{
  //======================================================================
  // 0) Energy slice
  //======================================================================
  const float eReco = recoPhoton.E();
  int iEslice = -1;
  for (int i = 0; i < N_Ebins; ++i)
    if (eReco >= eEdge[i] && eReco < eEdge[i+1]) { iEslice = i; break; }

  if (iEslice < 0) {
    if (Verbosity() > 0)
      std::cout << "[fillAshLogDx]  E=" << eReco
                << " GeV is outside [" << eEdge[0] << ',' << eEdge[N_Ebins]
                << ") – skip\n";
    return;
  }
  if (Verbosity() > 0)
    std::cout << "[fillAshLogDx]  slice " << iEslice
              << "  (" << eEdge[iEslice] << "–" << eEdge[iEslice+1]
              << " GeV)  recoE=" << eReco << '\n';

  //======================================================================
  // 1) Guard checks & lead-tower geometry
  //======================================================================
  if (!cluster || !m_bemcRec || !m_geometry) return;

  float eMax = -1.f;  RawCluster::TowerConstIterator best{};
  for (auto it = cluster->get_towers().first;
            it != cluster->get_towers().second; ++it)
    if (it->second > eMax) { eMax = it->second; best = it; }
  if (eMax < 1e-6) return;

  RawTowerGeom* leadGeo = m_geometry->get_tower_geometry(best->first);
  if (!leadGeo) return;

  const double rFront = leadGeo->get_center_radius();
  const double zFront = leadGeo->get_center_z();

  //======================================================================
  // 2) Truth photon → shower-depth position
  //======================================================================
  TVector3 pT(truthPhoton.Px(), truthPhoton.Py(), truthPhoton.Pz());
  const double tScale = rFront / pT.Pt();
  float xTsd,yTsd,zTsd;
  m_bemcRec->CorrectShowerDepth(truthPhoton.E(),
                                pT.x()*tScale, pT.y()*tScale, pT.z()*tScale,
                                xTsd,yTsd,zTsd);
  const double xTrue = std::hypot(xTsd,yTsd)*std::atan2(yTsd,xTsd);

  if (Verbosity() > 1)
    std::cout << "  truth (SD):  x=" << xTsd << "  y=" << yTsd
              << "  → xTrue=" << xTrue << '\n';

  //======================================================================
  // 3) Ash-b scan
  //======================================================================
  int nAshOK = 0, nAshMiss = 0;
  if (blockCord.second >= 0.f && blockCord.second <  1.f)
  {
    for (double bCell : m_bScan)
    {
      const double bVal = std::round(bCell*10000.)/10000.;
      const float  localAsh = doAshShift(blockCord.second, bVal);
      const float  globalPhi = (blockPhiBin*2 + localAsh*2.f)*(2.*M_PI/256.);

      const float xA = rFront*std::cos(globalPhi);
      const float yA = rFront*std::sin(globalPhi);

      float xRsd,yRsd,zRsd;
      m_bemcRec->CorrectShowerDepth(eReco, xA,yA,zFront, xRsd,yRsd,zRsd);
      const double xReco = std::hypot(xRsd,yRsd)*std::atan2(yRsd,xRsd);

      const TString hN = Form("h_dx_ash_b%.4f_E%d", bVal, iEslice);
      if (auto* h = dynamic_cast<TH1F*>( hm->getHisto( hN.Data() ) )) {
        h->Fill(xReco - xTrue);
        ++nAshOK;
        if (Verbosity() > 2)
          std::cout << "    [Ash] b=" << bVal
                    << "  Δx=" << xReco-xTrue << "  ► filled\n";
      } else {
        ++nAshMiss;
        if (Verbosity() > 2)
          std::cout << "    [Ash] b=" << bVal
                    << "  ► missing histo " << hN << '\n';
      }
    }
  }

  //======================================================================
  // 4) Log-w0 scan
  //======================================================================
  int nLogOK = 0, nLogMiss = 0;
  if (!tower_phis.empty() && !tower_energies.empty())
  {
    for (double w0 : m_w0Scan)
    {
      const double w0Val = std::round(w0*100.)/100.;
      const float  localLog = doLogWeightCoord(tower_phis,tower_energies,w0Val);
      const float  globalPhi = (blockPhiBin*2 + localLog*2.f)*(2.*M_PI/256.);

      const float xA = rFront*std::cos(globalPhi);
      const float yA = rFront*std::sin(globalPhi);

      float xRsd,yRsd,zRsd;
      m_bemcRec->CorrectShowerDepth(eReco, xA,yA,zFront, xRsd,yRsd,zRsd);
      const double xReco = std::hypot(xRsd,yRsd)*std::atan2(yRsd,xRsd);

        // --- Log section ---------------------------------------------------------
      const TString hN = Form("h_dx_log_w0%.2f_E%d", w0Val, iEslice);
      if (auto* h = dynamic_cast<TH1F*>( hm->getHisto( hN.Data() ) )) {
        h->Fill(xReco - xTrue);
        ++nLogOK;
        if (Verbosity() > 2)
          std::cout << "    [Log] w0=" << w0Val
                    << "  Δx=" << xReco-xTrue << "  ► filled\n";
      } else {
        ++nLogMiss;
        if (Verbosity() > 2)
          std::cout << "    [Log] w0=" << w0Val
                    << "  ► missing histo " << hN << '\n';
      }
    }
  }

  //======================================================================
  // 5) Summary for this cluster
  //======================================================================
  if (Verbosity() > 0)
    std::cout << "  ➜ Ash  filled/miss: " << nAshOK << '/' << nAshMiss
              << "   |   Log  filled/miss: " << nLogOK << '/' << nLogMiss
              << "\n";
}





void PositionDependentCorrection::fillDPhiRawAndCorrected(
    const TLorentzVector&  recoPhoton,
    const TLorentzVector&  truthPhoton,
    const std::pair<float,float>& blockCord,   // local (η,φ) ∈ [0,1]²
    int                     blockPhiBin,       // coarse φ-block 0…127
    float                   delPhi             // raw Δφ (reco – truth)
)
{
  /*--------------------------------------------------
   * 0) entrance summary
   *--------------------------------------------------*/
  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD << "[fillDPhiRawAndCorrected] ENTER" << ANSI_RESET
              << "\n  reco:(E,η,φ)=("  << recoPhoton.E()  << ','
                                        << recoPhoton.Eta() << ','
                                        << recoPhoton.Phi() << ')'
              << "\n  true:(E,η,φ)=("  << truthPhoton.E() << ','
                                        << truthPhoton.Eta() << ','
                                        << truthPhoton.Phi() << ')'
              << "\n  raw blockCord=(ηloc=" << blockCord.first
              << ", φloc="               << blockCord.second << ')'
              << "  blockPhiBin="        << blockPhiBin
              << "\n  input Δφ(raw)= "   << delPhi
              << std::endl;
  }

  /*------------------------------------------------------------
   * 1) energy-slice selection
   *----------------------------------------------------------*/
  const float thisE = recoPhoton.E();        // use full energy, not pT
  int iEbin = -1;
  for (int i = 0; i < N_Ebins; ++i)
    if (thisE >= eEdge[i] && thisE < eEdge[i+1]) { iEbin = i; break; }

  if (iEbin < 0)
  {
    if (Verbosity() > 0)
      std::cout << ANSI_YELLOW
                << "  [WARN] E=" << thisE
                << " GeV is outside histogram range [" << eEdge[0]
                << ',' << eEdge[N_Ebins] << ") – skipping Δφ fill\n"
                << ANSI_RESET;
    return;
  }
  else if (Verbosity() > 0)
  {
    std::cout << "  Energy slice  = " << iEbin
              << "  (" << eEdge[iEbin] << "–" << eEdge[iEbin+1] << " GeV)"
              << std::endl;
  }

  /*------------------------------------------------------------
   * 2) Δφ in (−π,π]
   *----------------------------------------------------------*/
  delPhi = TVector2::Phi_mpi_pi(delPhi);

  /*------------------------------------------------------------
   * 3) RAW histogram fills
   *----------------------------------------------------------*/
  if (h_phi_diff_raw_E[iEbin])
      h_phi_diff_raw_E[iEbin]->Fill(delPhi);
  else if (Verbosity() > 0)
      std::cout << ANSI_YELLOW << "  [WARN] h_phi_diff_raw_E[" << iEbin
                << "] pointer is null – RAW Δφ not filled\n" << ANSI_RESET;

  if (pr_phi_vs_blockcoord)
      pr_phi_vs_blockcoord->Fill(blockCord.second, delPhi);
  else if (Verbosity() > 0)
      std::cout << ANSI_YELLOW << "  [WARN] pr_phi_vs_blockcoord pointer is null\n"
                << ANSI_RESET;

  /*------------------------------------------------------------
   * 4)  Ash-b φ-correction
   *----------------------------------------------------------*/
  if (!isFitDoneForPhi)
  {
    if (Verbosity() > 0)
      std::cout << "  (No φ-fit available ⇒ only RAW filled)\n";
    return;
  }

  float bPhiUsed = m_bValsPhi[iEbin];        // cm if <1, otherwise block units
  constexpr float kCellSize = 5.55f;         // cm

  if (bPhiUsed > 1.f)  bPhiUsed /= kCellSize;  // convert from cm to “b”
  if (Verbosity() > 0)
    std::cout << "  Using bΦ = " << bPhiUsed << '\n';

  if (bPhiUsed <= 1e-9f)      /* nothing to correct */
  {
    if (Verbosity() > 0)
      std::cout << "  bΦ≈0  ⇒ skipping correction\n";
    return;
  }

  /* 4.1  local φ within block */
  float corrBlockPhi = doPhiBlockCorr(blockCord.second, bPhiUsed);

  /* wrap into [0,1) for safety */
  if (corrBlockPhi < 0.f)  corrBlockPhi += 1.f;
  if (corrBlockPhi >= 1.f) corrBlockPhi -= 1.f;

  if (Verbosity() > 0)
  {
    std::cout << "    raw  φloc=" << blockCord.second
              << "  →  corrected φloc=" << corrBlockPhi << '\n';
    if (corrBlockPhi < 0.f || corrBlockPhi >= 1.f)
      std::cout << ANSI_YELLOW
                << "    [WARN] corrected φloc fell outside [0,1) after wrap!\n"
                << ANSI_RESET;
  }

  if (h_localPhi_corrected[iEbin])
      h_localPhi_corrected[iEbin]->Fill(corrBlockPhi);
  else if (Verbosity() > 0)
      std::cout << ANSI_YELLOW << "  [WARN] h_localPhi_corrected[" << iEbin
                << "] pointer is null – local φ not filled\n" << ANSI_RESET;

  /* guard for illegal block index */
  if (blockPhiBin < 0 || blockPhiBin >= 128)
  {
    if (Verbosity() > 0)
      std::cout << ANSI_YELLOW
                << "  [WARN] blockPhiBin=" << blockPhiBin
                << " is outside [0,127] – global φ conversion skipped\n"
                << ANSI_RESET;
    return;
  }

  /* 4.2  convert to calorimeter global φ */
  float corrPhi = convertBlockToGlobalPhi(blockPhiBin, corrBlockPhi);

  if (Verbosity() > 0)
    std::cout << "    block(" << blockPhiBin << ") + φloc=" << corrBlockPhi
              << "  →  global φ=" << corrPhi << '\n';

  /* 4.3  corrected Δφ */
  float delPhiCorr = TVector2::Phi_mpi_pi(corrPhi - truthPhoton.Phi());

  if (h_phi_diff_corrected_E[iEbin])
      h_phi_diff_corrected_E[iEbin]->Fill(delPhiCorr);
  else if (Verbosity() > 0)
      std::cout << ANSI_YELLOW << "  [WARN] h_phi_diff_corrected_E[" << iEbin
                << "] pointer is null – CORR Δφ not filled\n" << ANSI_RESET;

  if (Verbosity() > 0)
  {
    std::cout << "    Δφ(corrected)= " << delPhiCorr
              << "\n[fillDPhiRawAndCorrected] EXIT\n";
  }
}




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
    float weight)
{
  // -----------------------------------------------------------------
  // 0) Basic information about the number of clusters & tower sum
  // -----------------------------------------------------------------
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

  // Fill #clusters in event
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

  // Retrieve all clusters
  RawClusterContainer::ConstRange clusterRange = clusterContainer->getClusters();
  RawClusterContainer::ConstIterator cIt1, cIt2;

  // For optional pT smearing
  float smear = 0.00;

  // Will hold booleans to see if first/second matched to meson photons
  bool match1 = false;
  bool match2 = false;

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

    // Recompute cluster 4-vector given event vertex
    CLHEP::Hep3Vector vertex(0,0,vtx_z);
      // ------------------------------------------------------------
      // Print debug: what's our input vertex z?
      // ------------------------------------------------------------
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG] finalClusterLoop: Using vtx_z=" << vtx_z
                  << " => vertex(0,0," << vtx_z << ")\n";
      }

    CLHEP::Hep3Vector E_vec_1 = RawClusterUtility::GetEVec(*clus1, vertex);
      
      // ------------------------------------------------------------
      // Additional check to see if the returned 3-vector is valid
      // ------------------------------------------------------------
      if (std::isnan(E_vec_1.mag()))
      {
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG] finalClusterLoop: WARNING - E_vec_1.mag() is NaN!\n"
                    << " -> cluster ID=" << clus1->get_id()
                    << "  vtx_z=" << vtx_z
                    << std::endl;
        }
        // Optionally continue or skip here, depending on what you prefer
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
        // Potentially skip or keep going
        continue;
      }

    float clusE   = E_vec_1.mag();
    float clusEta = E_vec_1.pseudoRapidity();
    float clusPhi = E_vec_1.phi();
    float clusPt  = E_vec_1.perp();
    float clusChi2= clus1->get_chi2();
      
      // ------------------------------------------------------------
      // Print debug: check cluster 4-vector values
      // ------------------------------------------------------------
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG] finalClusterLoop: cluster ID=" << clus1->get_id()
                  << " => E="   << clusE
                  << "  Pt="    << clusPt
                  << "  Eta="   << clusEta
                  << "  Phi="   << clusPhi
                  << "  Chi2="  << clusChi2
                  << std::endl;
      }

    // If we want to smear pT:
    clusPt *= rnd->Gaus(1, smear);

    // For QA, we note the lead tower indices
    int lt_eta = clus1->get_lead_tower().first;
    int lt_phi = clus1->get_lead_tower().second;

    // -----------------------------------------------------------------
    // (A) Basic cluster-level cuts
    // -----------------------------------------------------------------
    // 1) Minimum cluster E
    if (clusE < 0.1f)
    {
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG]  => Skipping cluster with E=" << clusE
                  << " (<0.1). \n";
      }
      continue;
    }
    // 2) Large chi2
    if (clusChi2 > 10)
    {
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG]  => Skipping cluster with chi2=" << clusChi2
                  << " (>10000). \n";
      }
      continue;
    }
    // 3) If lead tower index is out of range
    if (lt_eta > 95)
    {
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG]  => Skipping cluster lead_tower eta="
                  << lt_eta << " (>95). \n";
      }
      continue;
    }

    // Fill cluster‐energy QA
    h_clusE->Fill(clusE);

    // Gather the towers that make up this cluster (for block coords, etc.)
    RawCluster::TowerConstRange towerCR = clus1->get_towers();
    std::vector<int> towerEtas, towerPhis;
    std::vector<float> towerEs;

    int nTow=0;
    for (auto tIter = towerCR.first; tIter != towerCR.second; ++tIter)
    {
      nTow++;
      float twE = tIter->second;

      // geometry for bin indices
      RawTowerDefs::keytype tk = tIter->first;
      int iEta = m_geometry->get_tower_geometry(tk)->get_bineta();
      int iPhi = m_geometry->get_tower_geometry(tk)->get_binphi();

      towerEtas.push_back(iEta);
      towerPhis.push_back(iPhi);
      towerEs.push_back(twE);
    }

    // Fill cluster-size correlation
    h_clusE_nTow->Fill(clusE, nTow);

    // Build a TLorentzVector for cluster #1
    TLorentzVector photon1;
    photon1.SetPtEtaPhiE(clusPt, clusEta, clusPhi, clusE);

    // -----------------------------------------------------------------
    // (B) Check pT cut for the first photon
    // -----------------------------------------------------------------
    if (clusPt < pt1ClusCut || clusPt > ptMaxCut)
    {
      if (Verbosity() > 0)
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
      /* 1) energy–weighted tower indices ................................ */
      const float avgEta = getAvgEta (towerEtas, towerEs);      // 0 … 95.999
      const float avgPhi = getAvgPhi (towerPhis, towerEs);      // 0 … 255.999

      /* 2) coarse 2×2-block indices – now truly 0 … 47 / 0 … 127 */
      const int   block_eta_bin = int(std::floor(avgEta)) / 2;
      const int   block_phi_bin = int(std::floor(avgPhi)) / 2;

      /* 3) fine (local) coordinates inside that block, ∈ [0,1)²          */
      const std::pair<float,float> blockCord =
              getBlockCord(towerEtas, towerPhis, towerEs);        // (ηloc, φloc)

      /* 4) find energy slice ........................................... */
      int iEbin = -1;
      for (int i = 0; i < N_Ebins; ++i)
        if (clusE >= eEdge[i] && clusE < eEdge[i+1]) { iEbin = i; break; }

      /* optional detailed print ........................................ */
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG] finalClusterLoop →"
                  << "  blockCord=(" << blockCord.first << ',' << blockCord.second << ")"
                  << "  |  coarse (η,φ)=(" << block_eta_bin << ',' << block_phi_bin << ")"
                  << "  |  iEbin=" << iEbin
                  << "  |  E=" << clusE << "  pT=" << clusPt
                  << "  |  #towers=" << towerEs.size()
                  << '\n';
      }

      // ---------------------------------------------
      // 1)  RAW local-block coordinates  (no b-shift)
      // ---------------------------------------------
      {
        constexpr double kFillW = 1.0;      // weight per fill
        h3_cluster_block_cord_E->Fill(blockCord.first,   // η-coord   (x-axis)
                                      blockCord.second,  // φ-coord   (y-axis)
                                      clusE,             // energy    (z-axis)
                                      kFillW);

        if (Verbosity() > 0)
        {
          const int bx = h3_cluster_block_cord_E->GetXaxis()->FindBin(blockCord.first);
          const int by = h3_cluster_block_cord_E->GetYaxis()->FindBin(blockCord.second);
          const int bz = h3_cluster_block_cord_E->GetZaxis()->FindBin(clusE);

          std::cout << "[RAW-FILL]  ηloc="  << blockCord.first
                    << "  φloc="           << blockCord.second
                    << "  E="              << clusE
                    << "  (binX,Y,Z = "    << bx << "," << by << "," << bz << ")  "
                    << "→ new content = "
                    << h3_cluster_block_cord_E->GetBinContent(bx,by,bz)
                    << '\n';
        }
      }

      // -------------------------------------------------------------
      // 2)  CORRECTED local-block coordinates  (η and/or φ b-shift)
      //     – always fill so raw & corrected have the same statistics
      // -------------------------------------------------------------
      {
        /* ---- original raw local coordinates ---------------- */
        const float rawEta = blockCord.first;    // [0,1]
        const float rawPhi = blockCord.second;   // [0,1]

        /* ---- start from raw; apply shifts only if we have a fit --- */
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

        /* ---- always fill the corrected histogram ------------------ */
        constexpr double kFillW = 1.0;
        h3_cluster_block_cord_E_corrected->Fill(corrEta, corrPhi, clusE, kFillW);

        /* ---- optional diagnostics (verbosity-gated) ---------------- */
        if (Verbosity() > 0)
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
      }   // end corrected-block scope


    // QA: fill pT vs leadTower, fill cluster-level eta/phi
    h_pt_eta->Fill(clusPt, lt_eta);
    h_pt_eta_rw->Fill(clusPt, lt_eta, weight);
    h_etaphi_clus->Fill(clusEta, clusPhi);

    // -----------------------------------------------------------------
    // (D) Compare with TRUTH photons
    // -----------------------------------------------------------------
    for (auto & trPhoton : truth_photons)
    {
      float dR   = photon1.DeltaR(trPhoton);
      float ratio= (photon1.E() / trPhoton.E());
      float dPhi = photon1.Phi() - trPhoton.Phi();

      // keep dPhi in (-pi..+pi)
      if (dPhi > TMath::Pi())   dPhi -= 2.f*TMath::Pi();
      if (dPhi < -TMath::Pi())  dPhi += 2.f*TMath::Pi();

      // fill a distribution of dR
      h_delR_recTrth->Fill(dR);

      // check if matched
      if (dR < 0.03 && ratio < 1.3 && ratio > 0.3)
      {
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG] => cluster1 E=" << photon1.E()
                    << " matched to truthE=" << trPhoton.E()
                    << " dPhi=" << dPhi << " ratio=" << ratio << "\n";
        }
        // fill resolution hists
        h_matched_res->Fill(ratio, photon1.Eta());
        h_res_e->Fill(ratio, photon1.E());
        h_res->Fill(ratio);

        int iLTeta=lt_eta, iLTphi=lt_phi;

        h_res_e_eta->Fill(ratio, trPhoton.E(), iLTeta);
        h_res_e_phi->Fill(ratio, trPhoton.E(), iLTphi);

        h_delEta_e_eta->Fill(photon1.Eta()-trPhoton.Eta(), trPhoton.E(), iLTeta);
        h_delR_e_eta->Fill(dR, trPhoton.E(), iLTeta);
        h_delPhi_e_eta->Fill(dPhi, trPhoton.E(), iLTeta);
        h_delPhi_e_phi->Fill(dPhi, trPhoton.E(), iLTphi);

        h_truthE->Fill(trPhoton.E());

        // fill the ash/log-hist dx
        fillAshLogDx(clus1, photon1, trPhoton, blockCord,
                     block_phi_bin, towerPhis, towerEs);

        fillDPhiRawAndCorrected(photon1, trPhoton,
                                blockCord, block_phi_bin, dPhi);

        // position-dependent resolution
        if (block_eta_bin >=0 && block_eta_bin < NBinsBlock &&
            block_phi_bin >=0 && block_phi_bin < NBinsBlock)
        {
          h_res_block_E[block_eta_bin][block_phi_bin]->Fill(ratio, trPhoton.E());
        }
      }
    } // end loop truth_photons

    // -----------------------------------------------------------------
    // (E) Check if photon1 matches a meson‐decay photon
    // -----------------------------------------------------------------
    match1 = false;
    TLorentzVector ph1_trEtaPhi(0,0,0,0);
    for (auto & trPhot : truth_meson_photons)
    {
      float dR2   = photon1.DeltaR(trPhot);
      float ratio2= (photon1.E() / trPhot.E());
      if (dR2 < 0.03 && ratio2>0.7 && ratio2<1.5)
      {
        ph1_trEtaPhi.SetPtEtaPhiE( clusE / TMath::CosH(trPhot.Eta()),
                                   trPhot.Eta(),
                                   trPhot.Phi(),
                                   clusE );
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG] => meson-decay match found, E="
                    << ph1_trEtaPhi.E() << "  eta="
                    << ph1_trEtaPhi.Eta() << "\n";
        }
        match1 = true;
        break;
      }
    }

    // -----------------------------------------------------------------
    // 2) Inner loop over second cluster => form pi0
    // -----------------------------------------------------------------
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] => Starting INNER loop for cluster pairs.\n";
    }

    RawClusterContainer::ConstRange cRange2 = clusterContainer->getClusters();
    for (cIt2 = cRange2.first; cIt2 != cRange2.second; ++cIt2)
    {
      // avoid pairing with itself
      if (cIt2 == cIt1) continue;

      RawCluster* clus2 = cIt2->second;
      if (!clus2) continue;

      CLHEP::Hep3Vector E_vec_2 = RawClusterUtility::GetEVec(*clus2, vertex);

      float clus2E     = E_vec_2.mag();
      float clus2Eta   = E_vec_2.pseudoRapidity();
      float clus2Phi   = E_vec_2.phi();
      float clus2Pt    = E_vec_2.perp();
      float clus2Chi2  = clus2->get_chi2();

      // skip if pT < pt2ClusCut or pT>ptMaxCut
      if (clus2Pt < pt2ClusCut || clus2Pt > ptMaxCut)
      {
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG]    => cluster2 pT=" << clus2Pt
                    << " fails cut (<" << pt2ClusCut << " or >" << ptMaxCut
                    << "). \n";
        }
        continue;
      }
      if (clus2Chi2 > 10000)
      {
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG]    => cluster2 Chi2=" << clus2Chi2
                    << " (>10000), skipping.\n";
        }
        continue;
      }

      // check alpha
      float alpha = std::fabs(clusE - clus2E)/(clusE + clus2E);
      if (alpha>maxAlpha)
      {
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG]    => alpha=" << alpha
                    << " > maxAlpha=" << maxAlpha
                    << ", skipping.\n";
        }
        continue;
      }

      // Build TLorentzVector for cluster2
      TLorentzVector photon2;
      photon2.SetPtEtaPhiE(clus2Pt, clus2Eta, clus2Phi, clus2E);

      // optional check for meson-decay match #2
      match2 = false;
      TLorentzVector ph2_trEtaPhi(0,0,0,0);

      for (auto & trPhot : truth_meson_photons)
      {
        float dR3   = photon2.DeltaR(trPhot);
        float ratio3= (photon2.E()/ trPhot.E());
        if (dR3<0.02 && ratio3>0.7 && ratio3<1.5)
        {
          ph2_trEtaPhi.SetPtEtaPhiE( clus2E / TMath::CosH(trPhot.Eta()),
                                     trPhot.Eta(),
                                     trPhot.Phi(),
                                     clus2E);
          if (Verbosity() > 0)
          {
            std::cout << "[DEBUG]    => cluster2 matched meson-decay E="
                      << ph2_trEtaPhi.E() << "\n";
          }
          if (match1) match2 = true;
        }
      }

      // form pi0 4-vector
      TLorentzVector pi0_trKin = ph1_trEtaPhi + ph2_trEtaPhi; // truth-based
      TLorentzVector pi0       = photon1 + photon2;           // reco-based

      // check pi0 pT
      if (pi0.Pt()< pi0ptcut)
      {
        if (Verbosity() > 0)
        {
          std::cout << "[DEBUG]    => pi0 pT=" << pi0.Pt()
                    << " < pi0ptcut=" << pi0ptcut
                    << ", skipping.\n";
        }
        continue;
      }

      // If we get here => fill final histos
      h_pt1->Fill(photon1.Pt());
      h_pt2->Fill(photon2.Pt());

      h_InvMass->Fill(pi0.M());
      h_InvMass_w->Fill(pi0.M(), weight);

      // fill mass vs lead tower index
      h_mass_eta_lt->Fill(pi0.M(), lt_eta);
      h_mass_eta_lt_rw->Fill(pi0.M(), lt_eta, weight);

      // fill 3D distribution
      h_m_pt_eta->Fill(pi0.M(), pi0.E(), lt_eta);

      // fill block histogram if in valid bin
      if (block_eta_bin>=0 && block_eta_bin<NBinsBlock &&
          block_phi_bin>=0 && block_phi_bin<NBinsBlock)
      {
        h_mass_block_pt[block_eta_bin][block_phi_bin]->Fill(pi0.M(), pi0.E());
      }

      // if both matched => fill truth-based meson hist
      if (match2 && pi0_trKin.M()>0.001)
      {
        h_m_ptTr_eta->Fill(pi0.M(), truth_meson_photons.at(0).E(), lt_eta);
        h_m_ptTr_eta_trKin->Fill(pi0_trKin.M(), truth_meson_photons.at(0).E(), lt_eta);
      }

    } // end cluster2 loop

  } // end cluster1 loop

  // -----------------------------------------------------------------
  // Final summary
  // -----------------------------------------------------------------
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


int PositionDependentCorrection::process_towers(PHCompositeNode* topNode)
{
  // ------------------------------------------------------------------------
  // 0) Intro message
  // ------------------------------------------------------------------------
  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD << ANSI_CYAN
              << "[process_towers] START" << ANSI_RESET
              << "  => event counter: " << _eventcounter
              << std::endl;
  }

  // Print every 1000 events (for progress tracking)
  if ((_eventcounter % 1000) == 0)
  {
    // color highlight to see it better in the log
    std::cout << ANSI_BOLD << ANSI_YELLOW
              << "[Info] Processing event " << _eventcounter
              << ANSI_RESET << std::endl;
  }

  // Example line (unused):
  // float emcaldownscale = 1000000 / 800;

  float emcal_hit_threshold = 0.5;  // GeV

  if (debug)
  {
    std::cout << ANSI_BOLD << "[DEBUG] " << ANSI_RESET
              << "-----------------------------------" << std::endl;
  }

  // Define cuts
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

  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD
              << "[process_towers] Finished fillTowerInfo, total tower energy="
              << ANSI_RESET << tower_tot_e
              << ", #tower etas=" << ht_eta.size()
              << ", #tower phis=" << ht_phi.size() << std::endl;
  }

  // For later weighting (hard-coded = 1 in this snippet)
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

  // If vertex out of range => skip
  if (std::fabs(vtx_z) > _vz)
  {
    if (Verbosity() > 0)
    {
      std::cout << ANSI_BOLD << ANSI_YELLOW
                << "[process_towers] => Vertex Z out of range (|"
                << vtx_z << "| > " << _vz
                << ") => Skipping event..."
                << ANSI_RESET << std::endl
                << ANSI_BOLD << "[process_towers] END\n" << ANSI_RESET
                << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // ------------------------------------------------------------------------
  // 4) Check geometry node
  // ------------------------------------------------------------------------
  m_geometry = checkTowerGeometry(topNode);
  if (Verbosity() > 0 && m_geometry)
  {
    std::cout << ANSI_BOLD
              << "[process_towers] => Tower geometry retrieved successfully."
              << ANSI_RESET << std::endl;
  }

  // ------------------------------------------------------------------------
  // 5) Retrieve cluster container
  // ------------------------------------------------------------------------
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

  // ------------------------------------------------------------------------
  // 6) Cluster counting
  // ------------------------------------------------------------------------
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

  // ------------------------------------------------------------------------
  // 7) Final cluster-level loop (pi0 reconstruction, etc.)
  // ------------------------------------------------------------------------
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
                   weight);

  // Clear vectors for tower bin indices
  ht_phi.clear();
  ht_eta.clear();

  if (Verbosity() > 0)
  {
    std::cout << "[process_towers] => Completed finalClusterLoop, "
              << "clearing tower vectors now." << std::endl;
  }

  // ------------------------------------------------------------------------
  // 8) Return success code + end message
  // ------------------------------------------------------------------------
  if (Verbosity() > 0)
  {
    std::cout << ANSI_BOLD << "[process_towers] END" << ANSI_RESET
              << "\n" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


int PositionDependentCorrection::End(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] PositionDependentCorrection::End() - Entering End routine." << std::endl;
  }

  // Check if 'outfile' is valid before writing
  if (outfile)
  {
    // Move into the file directory
    outfile->cd();
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] Changed to output file directory. Now writing histograms..." << std::endl;
    }

    // Write all objects to the file
    outfile->Write();
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] Successfully wrote histograms to the output file: "
                << outfilename << std::endl;
      std::cout << "[DEBUG] Closing output file..." << std::endl;
    }
      
      if (Verbosity() > 0)
      {
        std::cout << "\n[INFO] Listing all histograms in the output file:\n";
        std::cout << std::left
                  << std::setw(30) << "Histogram Name"
                  << std::setw(10) << "Class"
                  << std::setw(12) << "#Entries"
                  << std::endl;
        std::cout << "--------------------------------------------------------------\n";

        // Change directory to outfile if not already
        outfile->cd();

        // Grab the list of all objects in this directory
        TList* listKeys = gDirectory->GetListOfKeys();
        if (listKeys)
        {
          // Iterate over each key
          TIter nextKey(listKeys);
          TKey* key = nullptr;
          while ((key = (TKey*)nextKey()))
          {
            // Get the class name and the object name
            const char* className = key->GetClassName();
            const char* objName   = key->GetName();

            // Read the object from file
            TObject* obj = key->ReadObj();
            if (!obj) continue;

            // If it’s a histogram (TH1‐derived), print out entries
            TH1* h1 = dynamic_cast<TH1*>(obj);
            if (h1)
            {
              // Print name, class, and entries
              std::cout << std::left
                        << std::setw(30) << objName
                        << std::setw(10) << className
                        << std::setw(12) << (Long64_t)h1->GetEntries()
                        << std::endl;
            }
            // optional else: if you want to handle TH2, TH3, etc. separately
          }
        }
        std::cout << "[INFO] Finished listing histograms.\n" << std::endl;
      }


    // Close the file
    outfile->Close();
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] Output file closed. Deleting the outfile pointer now." << std::endl;
    }

    // Delete the file pointer
    delete outfile;
    outfile = nullptr;
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cerr << "[ERROR] 'outfile' pointer is null. No histograms could be written!" << std::endl;
    }
  }

//  // Dump histograms via the HistoManager if valid
//  if (hm)
//  {
//    if (Verbosity() > 0)
//    {
//      std::cout << "[DEBUG] Dumping histograms to file '"
//                << outfilename << "' with mode 'UPDATE'." << std::endl;
//    }
//    hm->dumpHistos(outfilename, "UPDATE");
//
//    if (Verbosity() > 0)
//    {
//      std::cout << "[DEBUG] Histograms successfully dumped to file: " << outfilename << std::endl;
//    }
//  }
//  else
//  {
//    if (Verbosity() > 0)
//    {
//      std::cerr << "[ERROR] HistoManager (hm) pointer is null. Cannot dump histograms!" << std::endl;
//    }
//  }

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] PositionDependentCorrection::End() - Routine completed successfully." << std::endl;
  }

  return 0;
}

float PositionDependentCorrection::getWeight(int ieta, float pt)
{
  // Optional debug print
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] getWeight() called with ieta=" << ieta << ", pt=" << pt << std::endl;
  }

  // Maintain exact functionality
  if (ieta < 0 || ieta > 95)
  {
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] ieta is out of valid range [0..95]. Returning 0." << std::endl;
    }
    return 0;
  }

  float val = h_pt_rw[ieta]->GetBinContent(h_pt_rw[ieta]->FindBin(pt));
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Bin content for pt=" << pt
              << " in ieta=" << ieta
              << " is " << val << std::endl;
  }

  if (val == 0)
  {
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] Bin content is 0, returning 0 as weight." << std::endl;
    }
    return 0;
  }

  float result = 1 / val;
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] getWeight(): returning weight = " << result << std::endl;
  }
  return result;
}


// ============================================================================
//  getBlockCord – local (η,φ) of a cluster inside its 2×2 super-cell
//  ----------------------------------------------------------------------------
//  INPUT  towerEtas   int   vector, fine-tower η rows      (0 … 95)
//         towerPhis   int   vector, fine-tower φ columns   (0 … 255)
//         towerEs     float vector, tower energies (GeV)
//
//  OUTPUT std::pair<float,float>   (ηloc , φloc) in detector *tower units*
//
//  Range   −0.5 … +1.5  on each axis   ──> matches the X/Y edges you booked
// ============================================================================
std::pair<float,float>
PositionDependentCorrection::getBlockCord(const std::vector<int>   &towerEtas,
                                          const std::vector<int>   &towerPhis,
                                          const std::vector<float> &towerEs)
{
  // -------------------------------------------------------------------------
  // 0) Verbose banner
  // -------------------------------------------------------------------------
  if (Verbosity() > 0)
  {
    std::cout << ANSI_CYAN << "[getBlockCord] START  nTowers = "
              << towerEs.size() << ANSI_RESET << '\n';
  }

  // -------------------------------------------------------------------------
  // 1) Energy–weighted centres of gravity in *raw tower space*
  // -------------------------------------------------------------------------
  const float etaCoG = getAvgEta(towerEtas , towerEs);   //  0.0 … 95.0
  const float phiCoG = getAvgPhi(towerPhis , towerEs);   //  0.0 … 255.999

  if (Verbosity() > 0)
  {
    std::cout << "  ⟨η⟩_tower = " << etaCoG
              << "   ⟨φ⟩_tower = " << phiCoG << '\n';
  }

  // -------------------------------------------------------------------------
  // 2) Which 2×2 super-cell owns that point?
  //    Each super-cell spans *two* fine towers.
  // -------------------------------------------------------------------------
  const int iBlkEta = static_cast<int>(std::floor(etaCoG)) / 2;   // 0 … 47
  const int iBlkPhi = static_cast<int>(std::floor(phiCoG)) / 2;   // 0 … 127

  if (Verbosity() > 0)
    std::cout << "  → super-cell index (η,φ) = (" << iBlkEta
              << ',' << iBlkPhi << ")\n";

  // -------------------------------------------------------------------------
  // 3) Local coordinates inside that super-cell
  //    (subtract the block origin = 2·index)
  // -------------------------------------------------------------------------
  float locEta = etaCoG - 2.f * iBlkEta;   // −0.5 … +1.5
  float locPhi = phiCoG - 2.f * iBlkPhi;   // −0.5 … +1.5

  // Rare protection against pathological inputs: wrap η into the legal band.
  if (locEta >  1.5f) locEta -= 2.f;
  if (locEta < -0.5f) locEta += 2.f;
    
  if (locPhi >  1.5f) locPhi -= 2.f;
  if (locPhi < -0.5f) locPhi += 2.f;
    

  if (Verbosity() > 0)
  {
    std::cout << "  local(η,φ) inside cell = ("
              << locEta << ',' << locPhi << ")  [−0.5 … +1.5]\n";
  }

  // -------------------------------------------------------------------------
  // 4) Done
  // -------------------------------------------------------------------------
  if (Verbosity() > 0)
  {
    std::cout << ANSI_GREEN << "[getBlockCord] RETURN  "
              << "(ηloc,φloc)=(" << locEta << ',' << locPhi << ')'
              << ANSI_RESET << "\n\n";
  }

  return { locEta , locPhi };
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

  // -------------------------------------------------------------------------
  // 0) Sanity checks
  // -------------------------------------------------------------------------
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


/////////////////////////////////////////////////////////////////////////////////////
// getAvgPhi(...)
// --------------
// Similar to getAvgEta(...), but for phi indices. Here, we must handle wrapping issues
// because phi can wrap around from 255 back to 0 in tower space. The code ensures we
// don't inadvertently shift the average across the phi=0 boundary incorrectly.
/////////////////////////////////////////////////////////////////////////////////////

float PositionDependentCorrection::getAvgPhi(
    const std::vector<int> &towerphis,
    const std::vector<float> &towerenergies)
{
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] getAvgPhi() called with towerphis.size()="
                << towerphis.size() << std::endl;
    }

    // The EMCal has 256 phi bins, so we define nphibin = 256 for shifting logic.
    int nphibin = 256;
    float phimult = 0;
    float phisum  = 0;

    // We iterate over each tower in the cluster, adjusting its phi index if needed
    // so it doesn't jump artificially across the 0/256 boundary.
    for (UInt_t j = 0; j < towerenergies.size(); j++)
    {
        int phibin = towerphis[j];

        // If the difference between 'phibin' and the first tower's phi bin is < -128,
        // it means we've likely wrapped around the boundary, so we add 256.
        if (phibin - towerphis[0] < -nphibin / 2.0)
        {
            phibin += nphibin;
            if (Verbosity() > 0)
            {
              std::cout << "[DEBUG] phi wrap-around (low side): adjusted phibin => "
                        << phibin << std::endl;
            }
        }
        // Conversely, if the difference is > +128, we subtract 256 to wrap around
        // the other way.
        else if (phibin - towerphis[0] > +nphibin / 2.0)
        {
            phibin -= nphibin;
            if (Verbosity() > 0)
            {
              std::cout << "[DEBUG] phi wrap-around (high side): adjusted phibin => "
                        << phibin << std::endl;
            }
        }

        // Double-check that after shifting, we are within ±128 bins of the first tower.
        assert(std::abs(phibin - towerphis[0]) <= nphibin / 2.0);

        // Accumulate energy-weighted phi-bin index
        float energymult = towerenergies[j] * phibin;
        phimult += energymult;
        phisum  += towerenergies[j];
    }

    // Compute the average phi bin from the weighted sums
    float avgphi = phimult / phisum;

    // If the result is negative, shift it back into [0, 256) by adding nphibin
    if (avgphi < 0)
    {
       avgphi += nphibin;
       if (Verbosity() > 0)
       {
         std::cout << "[DEBUG] avgphi was negative; shifted to "
                   << avgphi << " by adding nphibin=" << nphibin << std::endl;
       }
    }

    // Also keep it modulo 256
    avgphi = fmod(avgphi, nphibin);

    // If the result is ≥ 255.5, we shift it back by 256 to ensure it's < 256
    if (avgphi >= 255.5)
    {
       avgphi -= nphibin;
       if (Verbosity() > 0)
       {
         std::cout << "[DEBUG] avgphi >= 255.5; shifted to "
                   << avgphi << " by subtracting nphibin=" << nphibin << std::endl;
       }
    }

    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] getAvgPhi() returning final adjusted phi="
                << avgphi << std::endl;
    }

    // Return the final adjusted phi coordinate in tower-bin space.
    return avgphi;
}
