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

R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)

using namespace std;

const double PositionDependentCorrection::ptEdge[PositionDependentCorrection::N_PT+1] =
{
    0.5, 1.0, 2.0, 4.0, 8.0
};

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
  _vz = 60.0;
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

  // h3_cluster_block_cord_pt:
  //   A 3D histogram that stores (x: local eta-block coordinate, y: local phi-block coordinate, z: cluster energy).
  //   This is used to map out how often clusters occur at each sub-tower position (0–1 in each dimension),
  //   which is the key to understanding the shower position distribution inside towers.
  h3_cluster_block_cord_pt = new TH3F("h2_cluster_block_cord_pt","",14,-0.5,1.5,14,-0.5,1.5,5,0,10);

  // h_block_bin:
  //   1D histogram that may be used to look at the block coordinate axis in discrete steps (0–1, or subdivided further).
  //   Potentially used to label or count how many clusters end up in each sub-cell bin across the 2×2 block space.
  h_block_bin = new TH1F("h_block_bin","",14,-0.5,1.5);

  // For measuring raw phi resolution
    // (A) Create the raw Δφ histograms for each pT bin
  for (int i = 0; i < NPTBINS; i++)
    {
      float ptLo = ptEdge[i];
      float ptHi = ptEdge[i+1];
      h_phi_diff_raw_pt[i] = new TH1F(
        Form("h_phi_diff_raw_%.0f_%.0f", ptLo, ptHi),  // histogram name
        Form("#Delta#phi raw, %.0f < pT < %.0f", ptLo, ptHi),
        200, -0.1, 0.1
      );
  }
    
    // (B) Create the corrected Δφ histograms for each pT bin
    for (int i = 0; i < NPTBINS; i++)
    {
      // Use your already-declared static member array ptEdge[]:
      float ptLo = ptEdge[i];
      float ptHi = ptEdge[i+1];

      // Now you have a properly named ptHi, so you can use it
      h_phi_diff_corrected_pt[i] = new TH1F(
        Form("h_phi_diff_corr_%.1f_%.1f", ptLo, ptHi),
        Form("#Delta#phi corrected, %.1f < pT < %.1f", ptLo, ptHi),
        200, -0.1, 0.1
      );
    }

    

  // Maybe a TProfile to see how #Delta#phi depends on local coordinate or energy
  pr_phi_vs_blockcoord = new TProfile("pr_phi_vs_blockcoord","",14,-0.5,1.5, -0.2,0.2);
    

    // GOOD: use ptEdge[i] and ptEdge[i+1]
    for (int i = 0; i < NPTBINS; i++)
    {
      // e.g. ptEdge[] = {2.0, 3.0, 5.0, 8.0, 12.0}
      float ptLo = ptEdge[i];
      float ptHi = ptEdge[i+1];

      h_localPhi_corrected[i] = new TH1F(
        Form("h_localPhi_corrected_%.1f_%.1f", ptLo, ptHi),
        Form("Corrected local #phi, %.1f < pT < %.1f", ptLo, ptHi),
        50, -0.5, 0.5
      );
    }
 

    {
        // Example: revised parameter scans with a wider range:
        const double bMin   = 0.10;
        const double bMax   = 1.00;
        const double bStep  = 0.05;

        const double w0Min  = 2.00;
        const double w0Max  = 6.00;
        const double w0Step = 0.10;

        // Clear old vectors
        m_bScan.clear();
        m_w0Scan.clear();

        // 1) Fill b‐scan values
        {
          int nb = static_cast<int>( std::round( (bMax - bMin)/bStep ) );
          for (int ib = 0; ib <= nb; ++ib)
          {
            double val = bMin + ib*bStep;
            // Snap to 3 decimals
            val = std::floor(val*1000. + 0.5)/1000.;
            m_bScan.push_back(val);
          }
        }

        // 2) Fill w0‐scan values
        {
          int nw = static_cast<int>( std::round( (w0Max - w0Min)/w0Step ) );
          for (int iw = 0; iw <= nw; ++iw)
          {
            double val = w0Min + iw*w0Step;
            // Snap to 2 decimals
            val = std::floor(val*100. + 0.5)/100.;
            m_w0Scan.push_back(val);
          }
        }

        const int N_B = m_bScan.size();
        const int N_W = m_w0Scan.size();

        // Make sure we only declare histograms once
        if (!alreadyDeclaredHistograms)
        {
          for (int ipt = 0; ipt < N_PT; ++ipt)
          {
            // ---------- Ash approach ----------
            for (int ib = 0; ib < N_B; ++ib)
            {
              double bVal = m_bScan[ib];

              // name example: "h_dx_ash_b0.300_pt0"
              TString hname = Form("h_dx_ash_b%.3f_pt%d", bVal, ipt);
              TH1F* hAsh = new TH1F(hname,
                                    Form("%s; x_{reco}-x_{true}; counts", hname.Data()),
                                    80, -4, 4 );
              hm->registerHisto(hAsh);
            }

            // ---------- Log‐weight approach ----------
            for (int iw = 0; iw < N_W; ++iw)
            {
              double w0Val = m_w0Scan[iw];

              // name example: "h_dx_log_w05.70_pt2"
              TString hname = Form("h_dx_log_w0%.2f_pt%d", w0Val, ipt);
              TH1F* hLog = new TH1F(hname,
                                    Form("%s; x_{reco}-x_{true}; counts", hname.Data()),
                                    80, -4, 4 );
              hm->registerHisto(hLog);
            }
          }
          alreadyDeclaredHistograms = true;
        }
    }


  rnd = new TRandom3();
    
    // Try reading bValues.txt
    {
      const std::string bFilePath = "/sphenix/u/patsfan753/scratch/PDCrun24pp/bParameters/bValues.txt";
      std::ifstream bfile(bFilePath.c_str());
      if (!bfile.is_open())
      {
        // Not found => no correction
        std::cout << "[INFO] bValues.txt NOT found at: " << bFilePath
                  << " => will NOT apply phi correction." << std::endl;
        isFitDoneForB = false;
      }
      else
      {
        std::cout << "[INFO] Found " << bFilePath << "; reading b-values for correction." << std::endl;

        // The 4 pT ranges you expect
        const float expectedLo[4] = {2.0, 3.0, 5.0, 8.0};
        const float expectedHi[4] = {3.0, 5.0, 8.0, 12.0};

        // Keep track of which bins we have found
        bool gotBin[4] = {false, false, false, false};

        // Prepare to read lines
        std::string line;
        while (std::getline(bfile, line))
        {
          if (line.empty() || line[0] == '#') continue; // skip empty or comment lines

          std::istringstream iss(line);
          float ptLo, ptHi, bVal;
          if (!(iss >> ptLo >> ptHi >> bVal))
          {
            // parse fail => skip
            continue;
          }

          // Attempt to match each expected bin
          for (int i = 0; i < 4; i++)
          {
            // If exact equality is safe, do (ptLo == expectedLo[i]) etc.
            // Or if worried about floating precision, do something like:
            auto closeEnough = [](float a, float b) {
              return std::fabs(a - b) < 1e-5;
            };

            // Check if this line’s pT range matches the i-th bin
            if ( closeEnough(ptLo, expectedLo[i]) &&
                 closeEnough(ptHi, expectedHi[i]) )
            {
              m_bVals[i] = bVal;    // store bVal
              gotBin[i]   = true;   // mark found
              break; // Stop checking the other bins
            }
          }
        } // end while reading lines

        bfile.close();

        // Now see if we got *all* 4 bins
        bool allFound = (gotBin[0] && gotBin[1] && gotBin[2] && gotBin[3]);
        if (allFound)
        {
          isFitDoneForB = true;
          std::cout << "[INFO] Successfully parsed b-values for all 4 bins:\n"
                    << "       [2,3] : " << m_bVals[0] << "\n"
                    << "       [3,5] : " << m_bVals[1] << "\n"
                    << "       [5,8] : " << m_bVals[2] << "\n"
                    << "       [8,12]: " << m_bVals[3] << "\n";
        }
        else
        {
          // Not enough lines matched => no correction
          isFitDoneForB = false;
          std::cout << "[WARNING] bValues.txt found but did NOT have the 4 needed lines.\n"
                    << "          => no phi correction.\n";
        }
      }
    } // end reading bValues.txt

    

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

  // localPhi is the block coordinate in [0,1].
  // Shift it to [-0.5,+0.5], apply the hyperbolic formula,
  // shift back. It's basically the same logic as getBlockCordCorr but for phi only.

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

float PositionDependentCorrection::convertBlockToGlobalPhi(int block_phi_bin, float localPhi)
{
  if (Verbosity() > 0)
  {
    std::cout << "PositionDependentCorrection::convertBlockToGlobalPhi - START" << std::endl;
    std::cout << "  --> block_phi_bin: " << block_phi_bin
              << ", localPhi: " << localPhi << std::endl;
    std::cout << "  --> Computing coarseTowerPhi = block_phi_bin * 2.0, etc..." << std::endl;
  }

  // block_phi_bin is the integer bin (0..255 for CEMC).
  // localPhi is in [0,1].
  // so the "global" tower phi index is block_phi_bin*2 + localPhi*(2).
  // Then map to [-pi, +pi].

  float coarseTowerPhi = block_phi_bin * 2.0; // each block is 2 wide
  float fullPhiIndex   = coarseTowerPhi + localPhi*2.0;

  // Now convert that “phi index” to actual radians:
  // typical mapping: index=0 => phi=0, index=256 => phi=2π
  // so 1 index = 2π/256 = ~0.0245 rad
  float radPerBin = (2.0 * M_PI) / 256.0;
  float globalPhi = fullPhiIndex * radPerBin;

  // shift into [-π, π] if needed:
  globalPhi = TVector2::Phi_mpi_pi(globalPhi);

  if (Verbosity() > 0)
  {
    std::cout << "  --> fullPhiIndex = " << fullPhiIndex << ", radPerBin = " << radPerBin << std::endl;
    std::cout << "  --> Final globalPhi (mapped into [-π, π]) = " << globalPhi << std::endl;
    std::cout << "PositionDependentCorrection::convertBlockToGlobalPhi - END" << std::endl << std::endl;
  }

  return globalPhi;
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


void PositionDependentCorrection::fillAshLogDx(
    RawCluster* cluster,               // pass in the cluster we are analyzing
    const TLorentzVector& recoPhoton,
    const TLorentzVector& truthPhoton,
    const std::pair<float,float>& blockCord, // local block coords in [0,1]
    int blockPhiBin,
    const std::vector<int>& tower_phis,
    const std::vector<float>& tower_energies
)
{
  //====================================================================
  // 0) Identify which pT slice the reco photon is in [2..12]
  //====================================================================
  float ptReco = recoPhoton.Pt();
  int iPtSlice = -1;

  for (int i = 0; i < N_PT; ++i)
  {
    if (ptReco >= ptEdge[i] && ptReco < ptEdge[i+1])
    {
      iPtSlice = i;
      break;
    }
  }
  if (iPtSlice < 0)
  {
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] fillAshLogDx => Reco pT=" << ptReco
                << " is not in [2..12], skipping.\n";
    }
    return;
  }

  //====================================================================
  // 1) Retrieve the cluster’s “lead tower” geometry by scanning
  //====================================================================
  // We require a valid cluster & BEmcRec pointer for further steps.
  if (!cluster)
  {
    if (Verbosity() > 0)
    {
      std::cout << "[WARN] fillAshLogDx => 'cluster' is NULL => skipping fill.\n";
    }
    return;
  }
  if (!m_bemcRec)
  {
    if (Verbosity() > 0)
    {
      std::cout << "[WARN] fillAshLogDx => m_bemcRec is NULL => cannot do shower-depth correction.\n"
                << "       => skipping fill.\n";
    }
    return;
  }

  // Loop over all towers in this cluster to find the one with the largest energy.
  RawCluster::TowerConstRange towers = cluster->get_towers();
  RawCluster::TowerConstIterator bestTowIter = towers.first;
  float maxE = -1.f;
  for (auto it = towers.first; it != towers.second; ++it)
  {
    float towerE = it->second;  // energy contributed by this tower
    if (towerE > maxE)
    {
      maxE        = towerE;
      bestTowIter = it;
    }
  }

  // If no tower has positive energy, skip
  if (maxE < 1e-9)
  {
    if (Verbosity() > 0)
    {
      std::cout << "[WARN] fillAshLogDx => no positive-energy tower in cluster => skipping.\n";
    }
    return;
  }

  // Retrieve this tower’s geometry
  RawTowerDefs::keytype bestTowerKey = bestTowIter->first;
  RawTowerGeom* leadTowerGeom = m_geometry->get_tower_geometry(bestTowerKey);
  if (!leadTowerGeom)
  {
    if (Verbosity() > 0)
    {
      std::cout << "[ERROR] fillAshLogDx => tower geom is NULL for key="
                << bestTowerKey << " => skipping.\n";
    }
    return;
  }

  // Extract the integer bin indices for optional QA/logging
  int leadEta = leadTowerGeom->get_bineta();
  int leadPhi = leadTowerGeom->get_binphi();

  // Approx. front-face radius
  double rFront = leadTowerGeom->get_center_radius();

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] fillAshLogDx => found lead tower with binEta=" << leadEta
              << ", binPhi=" << leadPhi << ", E=" << maxE << "\n"
              << "                cluster pT=" << ptReco
              << ", iPtSlice=" << iPtSlice << "\n"
              << "                rFront ~ " << rFront << std::endl;
  }

  //====================================================================
  // 2) Project the TRUTH photon to that same radius => apply CorrectShowerDepth
  //====================================================================
  float pxTruth = truthPhoton.Px();
  float pyTruth = truthPhoton.Py();
  float pzTruth = truthPhoton.Pz();
  double ptTruth = std::sqrt(pxTruth*pxTruth + pyTruth*pyTruth);
  if (ptTruth < 1e-9)
  {
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] fillAshLogDx => truth photon pT ~ 0 => skipping fill.\n";
    }
    return;
  }

  // Solve for how far we need to project in x,y to hit radius = rFront
  double tTruth = rFront / ptTruth;
  double xAtrue = pxTruth * tTruth;
  double yAtrue = pyTruth * tTruth;
  double zAtrue = pzTruth * tTruth;  // If there’s vertex offset, you could add it here

  // Depth-correct for the truth photon’s energy
  float xT=0.f, yT=0.f, zT=0.f;
  m_bemcRec->CorrectShowerDepth(truthPhoton.E(),
                                xAtrue, yAtrue, zAtrue,
                                xT,      yT,     zT);

  // Convert to “arc-length”: xTrue = rTrue * phiTrue
  double rTrue   = std::sqrt(xT*xT + yT*yT);
  double phiTrue = std::atan2(yT, xT);
  double xTrue   = rTrue * phiTrue;

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] fillAshLogDx => truth E=" << truthPhoton.E()
              << ", front-face coords (" << xAtrue << "," << yAtrue << "," << zAtrue << ")\n"
              << "                       => after depth corr => ("
              << xT << "," << yT << "," << zT << ")\n"
              << "                       => xTrue = " << xTrue << "\n";
  }

  //====================================================================
  // 3) Ash formula: local block φ => global => shower-depth => fill
  //====================================================================
  if (blockCord.second < 0.f || blockCord.second > 1.f)
  {
    if (Verbosity() > 0)
    {
      std::cout << "[WARN] fillAshLogDx => local block φ=" << blockCord.second
                << " is not in [0,1], skipping Ash loop.\n";
    }
  }
  else
  {
    for (double bVal : m_bScan)
    {
      float localAsh = doAshShift(blockCord.second, bVal);

      // Convert that local-φ to a global tower φ index => multiply by 2 => total
      float coarsePhiInd = blockPhiBin*2 + localAsh*2.f;
      float globalPhi    = coarsePhiInd * (2.f*M_PI / 256.f);

      float xA = rFront * std::cos(globalPhi);
      float yA = rFront * std::sin(globalPhi);

      // Depth-correct for the *reco* cluster’s energy
      float xC=0.f, yC=0.f, zC=0.f;
      m_bemcRec->CorrectShowerDepth(recoPhoton.E(),
                                    xA, yA, /*zA=*/0.f,
                                    xC, yC, zC);

      double rReco   = std::sqrt(xC*xC + yC*yC);
      double phiReco = std::atan2(yC, xC);
      double xReco   = rReco * phiReco;

      // Fill “h_dx_ash_bXX_ptYY” if it exists
      TString hname = Form("h_dx_ash_b%.3f_pt%d", bVal, iPtSlice);
      TH1F* hAsh = dynamic_cast<TH1F*>(hm->getHisto(hname.Data()));
      if (!hAsh)
      {
        if (Verbosity() > 0)
        {
          std::cout << "[WARN] fillAshLogDx => histogram " << hname
                    << " not found => no fill.\n";
        }
      }
      else
      {
        double dx = xReco - xTrue;
        hAsh->Fill(dx);

        if (Verbosity() > 1)
        {
          std::cout << "[DBX] fillAshLogDx => bVal=" << bVal
                    << ", localAsh=" << localAsh
                    << ", globalPhi=" << globalPhi
                    << ", xReco=" << xReco
                    << " => dx=" << dx << "\n";
        }
      }
    }
  }

  //====================================================================
  // 4) Log-weight formula: localLog => global => shower-depth => fill
  //====================================================================
  if (tower_phis.empty() || tower_energies.empty())
  {
    if (Verbosity() > 0)
    {
      std::cout << "[WARN] fillAshLogDx => tower_phis or tower_energies empty => skipping log-weight loop.\n";
    }
  }
  else
  {
    for (double w0Val : m_w0Scan)
    {
      float localLog = doLogWeightCoord(tower_phis, tower_energies, w0Val);

      float coarsePhiInd = blockPhiBin*2 + localLog*2.f;
      float globalPhi    = coarsePhiInd * (2.f*M_PI / 256.f);

      float xA = rFront * std::cos(globalPhi);
      float yA = rFront * std::sin(globalPhi);

      float xC=0.f, yC=0.f, zC=0.f;
      m_bemcRec->CorrectShowerDepth(recoPhoton.E(),
                                    xA, yA, /*zA=*/0.f,
                                    xC, yC, zC);

      double rReco   = std::sqrt(xC*xC + yC*yC);
      double phiReco = std::atan2(yC, xC);
      double xReco   = rReco * phiReco;

      TString hname2 = Form("h_dx_log_w0%.2f_pt%d", w0Val, iPtSlice);
      TH1F* hLog = dynamic_cast<TH1F*>(hm->getHisto(hname2.Data()));
      if (!hLog)
      {
        if (Verbosity() > 0)
        {
          std::cout << "[WARN] fillAshLogDx => histogram " << hname2
                    << " not found => no fill.\n";
        }
      }
      else
      {
        double dx = xReco - xTrue;
        hLog->Fill(dx);

        if (Verbosity() > 1)
        {
          std::cout << "[DBX] fillAshLogDx => w0Val=" << w0Val
                    << ", localLog=" << localLog
                    << ", globalPhi=" << globalPhi
                    << ", xReco=" << xReco
                    << " => dx=" << dx << "\n";
        }
      }
    }
  }

  //====================================================================
  // Summary
  //====================================================================
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] fillAshLogDx => DONE for pT slice=" << iPtSlice
              << ", cluster pT=" << ptReco
              << "\n                 lead tower (binEta=" << leadEta
              << ", binPhi=" << leadPhi
              << "), possibly filled Ash & Log histograms.\n";
  }
}


// ----------------------------------------------------------------------
//  fillDPhiRawAndCorrected
//    * Extracted from the "raw Δφ" fill and "corrected Δφ" fill
//    * Maintains pT-binning [2..8] for iPtBin
// ----------------------------------------------------------------------
void PositionDependentCorrection::fillDPhiRawAndCorrected(
    const TLorentzVector& recoPhoton,
    const TLorentzVector& truthPhoton,
    const std::pair<float,float>& blockCord,
    int blockPhiBin,
    float delPhi
)
{
  // 1) Determine which of the 4 pT bins we’re in
  float thisPt = recoPhoton.Pt();  // or clus_pt

  int iPtBin = -1;
  for (int i = 0; i < NPTBINS; ++i)  // NPTBINS == 4
  {
    if (thisPt >= ptEdge[i] && thisPt < ptEdge[i+1])
    {
      iPtBin = i;
      break;
    }
  }
  // If the cluster’s pT is not in [2..12], skip
  if (iPtBin < 0)
  {
    return;
  }

  // 2) Compute the uncorrected ∆φ in the range (-π..+π)
  delPhi = TVector2::Phi_mpi_pi(delPhi);

  // Optional debug
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] fillDPhiRawAndCorrected:\n"
              << "   * cluster E=" << recoPhoton.E()
              << "   * truth E="   << truthPhoton.E()
              << "   * ∆R="        << recoPhoton.DeltaR(truthPhoton)
              << "   * ∆φ(raw)="   << delPhi
              << std::endl;
  }

  // 3) Fill uncorrected histogram
  if (h_phi_diff_raw_pt[iPtBin])
  {
    h_phi_diff_raw_pt[iPtBin]->Fill(delPhi);
  }
  else if (Verbosity() > 0)
  {
    std::cout << "[WARNING] h_phi_diff_raw_pt[" << iPtBin
              << "] is null => no fill!" << std::endl;
  }

  // 4) (Optional) fill TProfile vs. local block φ coordinate
  if (pr_phi_vs_blockcoord)
  {
    pr_phi_vs_blockcoord->Fill(blockCord.second, delPhi);
  }

  // 5) If we have read in b-values, apply the correction
  if (isFitDoneForB)
  {
    float bPhiUsed = m_bVals[iPtBin];  // pick the b-value for this pT bin
    if (bPhiUsed > 1e-9)
    {
      // (a) Correct the local block φ coordinate
      float correctedBlockPhi = doPhiBlockCorr(blockCord.second, bPhiUsed);

      // Keep track of how that local φ changes (optional QA)
      if (h_localPhi_corrected[iPtBin])
      {
        h_localPhi_corrected[iPtBin]->Fill(correctedBlockPhi);
      }

      // (b) Convert that corrected local coordinate to global φ
      float correctedPhi = convertBlockToGlobalPhi(blockPhiBin, correctedBlockPhi);

      // (c) Recompute ∆φ in (-π..+π)
      float delPhiCorr = TVector2::Phi_mpi_pi(correctedPhi - truthPhoton.Phi());

      // (d) Fill corrected histogram
      if (h_phi_diff_corrected_pt[iPtBin])
      {
        h_phi_diff_corrected_pt[iPtBin]->Fill(delPhiCorr);
      }
      else if (Verbosity() > 0)
      {
        std::cout << "[WARNING] h_phi_diff_corrected_pt[" << iPtBin
                  << "] is null => no fill!" << std::endl;
      }
    }
    else if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] bPhiUsed=0 => skipping local φ correction.\n";
    }
  } // end if isFitDoneForB
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
  // 0) Trigger decoding and check
  // -----------------------------------------------------------------
//  if (!trigAna)
//  {
//    // Original error message:
//    std::cerr << "[ERROR] No TriggerAnalyzer pointer!\n"
//              << "Cannot determine triggers -> skipping event.\n";
//
//    // Additional verbosity-based info:
//    if (Verbosity() > 0)
//    {
//      std::cout << "[Verbosity] finalClusterLoop: 'trigAna' is nullptr. This prevents decoding triggers." << std::endl;
//      std::cout << "  --> Returning early from finalClusterLoop." << std::endl;
//    }
//    return; // or return Fun4AllReturnCodes::ABORTEVENT; if you prefer.
//  }

//  // Optionally print debug info (note: existing logic checks Verbosity, not the 'debug' variable here):
//  if (Verbosity() > 0)
//  {
//    std::cout << "[Verbosity] finalClusterLoop: About to call trigAna->decodeTriggers() on topNode = "
//              << topNode << std::endl;
//  }
//    if (!isSimulation)  // or however you detect data vs sim
//    {
//       trigAna->decodeTriggers(topNode);
//    }


//  if (Verbosity() > 0)
//  {
//    std::cout << "[Verbosity] finalClusterLoop: Building list of fired triggers from triggerNameMap..." << std::endl;
//  }
//
//  // Build a vector of short-name triggers that fired
//  std::vector<std::string> activeTriggerNames;
//  activeTriggerNames.reserve(triggerNameMap.size());
//
//  for (const auto &kv : triggerNameMap)
//  {
//    const std::string &dbTriggerName   = kv.first;   // e.g. "MBD N&S >= 1"
//    const std::string &histFriendlyStr = kv.second;  // e.g. "MBD_NandS_geq_1"
//
//    if (trigAna->didTriggerFire(dbTriggerName))
//    {
//      activeTriggerNames.push_back(histFriendlyStr);
//
//      if (Verbosity() > 0)
//      {
//        std::cout << "[Verbosity]  -> Trigger fired: \"" << dbTriggerName
//                  << "\" => short name \"" << histFriendlyStr << "\"\n";
//      }
//    }
//  }

//  // If none of our triggers fired, skip the hist-filling
//  if (activeTriggerNames.empty())
//  {
//    if (Verbosity() > 0)
//    {
//      std::cout << "[Verbosity] finalClusterLoop: No triggers from triggerNameMap are active. Skipping event." << std::endl;
//    }
//    return;
//  }
    // -----------------------------------------------------------------
    // 1) Basic cluster counting and QA
    // -----------------------------------------------------------------
    if (debug)
    {
      // Print tower totals in bold green for easy scanning
      std::cout << ANSI_BOLD << ANSI_GREEN
                << "[DEBUG] tower tot E=" << tower_tot_e
                << "    nClusContainer=" << nClusCount
                << ANSI_RESET << std::endl;
    }

    // Fill a histogram with the count of clusters in this event
    h_nclusters->Fill(nClusCount);

    // If the number of clusters is too large, skip
    if (nClusCount > max_nClusCount)
    {
      if (Verbosity() > 0)
      {
        std::cout << ANSI_BOLD << ANSI_YELLOW
                  << "[Verbosity] finalClusterLoop: nClusCount (" << nClusCount
                  << ") exceeds max_nClusCount (" << max_nClusCount
                  << "). Skipping event."
                  << ANSI_RESET << std::endl;
      }
      return;
    }

    // Retrieve cluster range
    RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
    RawClusterContainer::ConstIterator clusterIter, clusterIter2;

    // optional pT smearing
    float smear = 0.00;

    // Booleans for tracking matches to meson photons
    bool match1 = false;
    bool match2 = false;

    // Optionally store kinematics in these vectors (unused, but retained)
    std::vector<float> save_pt;
    std::vector<float> save_eta;
    std::vector<float> save_phi;
    std::vector<float> save_e;

    if (Verbosity() > 0)
    {
      std::cout << ANSI_BOLD
                << "[Verbosity] finalClusterLoop: Beginning outer loop over clusters."
                << ANSI_RESET << std::endl;
    }

    // -----------------------------------------------------------------
    // 2) Outer loop over clusters
    // -----------------------------------------------------------------
    for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; ++clusterIter)
    {
      RawCluster* recoCluster = clusterIter->second;

      // Recompute vector with event vertex
      CLHEP::Hep3Vector vertex(0, 0, vtx_z);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetEVec(*recoCluster, vertex);

      float clusE      = E_vec_cluster.mag();
      float clus_eta   = E_vec_cluster.pseudoRapidity();
      float clus_phi   = E_vec_cluster.phi();
      float clus_pt    = E_vec_cluster.perp();
      float clus_chisq = recoCluster->get_chi2();

      // optional pT smearing (default = 0)
      clus_pt *= rnd->Gaus(1, smear);

      // lead tower (eta, phi) bins from the cluster
      int lt_eta = recoCluster->get_lead_tower().first;
      int lt_phi = recoCluster->get_lead_tower().second;

      // Basic cluster cuts
      if (clusE < 0.1)         continue;
      if (clus_chisq > 10000)  continue;

      // Fill cluster‐energy QA
      h_clusE->Fill(clusE);

      // Gather the towers in this cluster
      RawCluster::TowerConstRange towerCR = recoCluster->get_towers();
      int nTow = 0;
      std::vector<int> tower_etas, tower_phis;
      std::vector<float> tower_energies;

      for (auto toweriter = towerCR.first; toweriter != towerCR.second; ++toweriter)
      {
        nTow++;
        tower_etas.push_back( m_geometry->get_tower_geometry(toweriter->first)->get_bineta() );
        tower_phis.push_back( m_geometry->get_tower_geometry(toweriter->first)->get_binphi() );
        tower_energies.push_back( toweriter->second );
      }

      // Compute local (block) coordinates for the cluster
      std::pair<float,float> blockCord = getBlockCord(tower_etas, tower_phis, tower_energies);

      // Convert local coordinates to discrete block indices
      int block_eta_bin = h_block_bin->FindBin(blockCord.first)  - 1;
      int block_phi_bin = h_block_bin->FindBin(blockCord.second) - 1;

      // Fill 3D histogram: local coord (x,y) vs. cluster energy (z)
      h3_cluster_block_cord_pt->Fill(blockCord.first, blockCord.second, clusE);

      // Fill cluster‐size correlation
      h_clusE_nTow->Fill(clusE, nTow);

      // lead tower range check
      if (lt_eta > 95) continue;

      // Fill pT vs. lead tower
      h_pt_eta->Fill(clus_pt, lt_eta);
      h_pt_eta_rw->Fill(clus_pt, lt_eta, weight);

      // Fill cluster-level eta/phi QA
      h_etaphi_clus->Fill(clus_eta, clus_phi);

      // Build TLorentzVector for the candidate cluster
      TLorentzVector photon1;
      photon1.SetPtEtaPhiE(clus_pt, clus_eta, clus_phi, clusE);

      // ----------------------------------------------------
      // Compare with truth photons (NOT necessarily pi0 decays)
      // ----------------------------------------------------
      for (auto &tr_phot : truth_photons)
      {
        float delR   = photon1.DeltaR(tr_phot);
        float res    = (photon1.E() / tr_phot.E());
        float delPhi = photon1.Phi() - tr_phot.Phi();

        // Keep delPhi in (-pi, +pi)
        if (delPhi > TMath::TwoPi())  delPhi -= TMath::TwoPi();
        if (delPhi < -TMath::TwoPi()) delPhi += TMath::TwoPi();

        // Fill global distribution
        h_delR_recTrth->Fill(delR);

        // Check if cluster is matched
        if (delR < 0.03 && res < 1.3 && res > 0.3)
        {
          if (debug)
          {
            std::cout << ANSI_BOLD << "[DEBUG] match clusE=" << photon1.E()
                      << "  truthE=" << tr_phot.E()
                      << "  delPhi=" << delPhi
                      << ANSI_RESET << std::endl;
          }
          // Fill resolution histograms
          h_matched_res->Fill(res, photon1.Eta());
          h_res_e->Fill(res, photon1.E());
          h_res->Fill(res);

          h_res_e_eta->Fill(res, tr_phot.E(), lt_eta);
          h_res_e_phi->Fill(res, tr_phot.E(), lt_phi);

          h_delEta_e_eta->Fill(photon1.Eta() - tr_phot.Eta(), tr_phot.E(), lt_eta);
          h_delR_e_eta->Fill(delR, tr_phot.E(), lt_eta);
          h_delPhi_e_eta->Fill(delPhi, tr_phot.E(), lt_eta);
          h_delPhi_e_phi->Fill(delPhi, tr_phot.E(), lt_phi);

          h_truthE->Fill(tr_phot.E());
            

            // snippet inside finalClusterLoop, after matching a photon:
            fillAshLogDx(
                recoCluster,          // <--- first argument is RawCluster*
                photon1,              // recoPhoton
                tr_phot,              // truthPhoton
                blockCord,
                block_phi_bin,
                tower_phis,
                tower_energies
            );

          fillDPhiRawAndCorrected(photon1, tr_phot, blockCord, block_phi_bin, delPhi);

            
          // position‐dependent resolution
          if (block_eta_bin >= 0 && block_eta_bin < NBinsBlock &&
              block_phi_bin >= 0 && block_phi_bin < NBinsBlock)
          {
            h_res_block_E[block_eta_bin][block_phi_bin]->Fill(res, tr_phot.E());
          }
        }
      } // end loop over truth_photons

      // ----------------------------------------------------
      // Check if photon1 matches a meson‐decay photon
      // ----------------------------------------------------
      TLorentzVector ph1_trEtaPhi(0,0,0,0);
      match1 = false;

      for (auto &tr_phot : truth_meson_photons)
      {
        float delR = photon1.DeltaR(tr_phot);
        float res  = photon1.E() / tr_phot.E();

        // If matched in ΔR and ratio, store direction from truth, energy from cluster
        if (delR < 0.03 && res > 0.7 && res < 1.5)
        {
          ph1_trEtaPhi.SetPtEtaPhiE(clusE / TMath::CosH(tr_phot.Eta()),
                                    tr_phot.Eta(),
                                    tr_phot.Phi(),
                                    clusE);

          if (debug)
          {
            std::cout << ANSI_BOLD
                      << "[DEBUG] meson-decay match  eta=" << ph1_trEtaPhi.Eta()
                      << " E="   << ph1_trEtaPhi.E()
                      << ANSI_RESET << std::endl;
          }
          match1 = true;
          break; // found a match
        }
      }

      // pT cut for the first photon candidate
      if (clus_pt < pt1ClusCut || clus_pt > ptMaxCut) continue;

      // -----------------------------------------------------------------
      // 3) Inner loop over a second cluster to form pi0
      // -----------------------------------------------------------------
      for (clusterIter2 = clusterEnd.first; clusterIter2 != clusterEnd.second; ++clusterIter2)
      {
        // Avoid pairing the same cluster with itself
        if (clusterIter2 == clusterIter) continue;

        RawCluster* recoCluster2 = clusterIter2->second;
        CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetEVec(*recoCluster2, vertex);

        float clus2E      = E_vec_cluster2.mag();
        float clus2_eta   = E_vec_cluster2.pseudoRapidity();
        float clus2_phi   = E_vec_cluster2.phi();
        float clus2_pt    = E_vec_cluster2.perp();
        float clus2_chisq = recoCluster2->get_chi2();

        if (clus2_pt < pt2ClusCut || clus2_pt > ptMaxCut)    continue;
        if (clus2_chisq > 10000)                            continue;

        TLorentzVector photon2;
        photon2.SetPtEtaPhiE(clus2_pt, clus2_eta, clus2_phi, clus2E);

        // alpha cut
        float alpha = std::fabs(clusE - clus2E) / (clusE + clus2E);
        if (alpha > maxAlpha) continue;

        // check if cluster2 also matches meson photon
        TLorentzVector ph2_trEtaPhi(0,0,0,0);
        match2 = false;

        for (auto &tr_phot : truth_meson_photons)
        {
          float delR = photon2.DeltaR(tr_phot);
          float res  = photon2.E() / tr_phot.E();

          if (delR < 0.02 && res > 0.7 && res < 1.5)
          {
            ph2_trEtaPhi.SetPtEtaPhiE(clus2E / TMath::CosH(tr_phot.Eta()),
                                      tr_phot.Eta(),
                                      tr_phot.Phi(),
                                      clus2E);

            if (debug)
            {
              std::cout << ANSI_BOLD
                        << "[DEBUG] meson-decay match  eta=" << ph2_trEtaPhi.Eta()
                        << " E="   << ph2_trEtaPhi.E()
                        << ANSI_RESET << std::endl;
            }
            if (match1) match2 = true;
          }
        }

        // build pi0 4‐vectors
        TLorentzVector pi0_trKin = ph1_trEtaPhi + ph2_trEtaPhi; // truth-based directions
        TLorentzVector pi0       = photon1 + photon2;           // reco-based

        // check pi0 pT
        if (pi0.Pt() < pi0ptcut) continue;

        // fill single cluster pT
        h_pt1->Fill(photon1.Pt());
        h_pt2->Fill(photon2.Pt());

        // fill diphoton mass distributions
        h_InvMass->Fill(pi0.M());
        h_InvMass_w->Fill(pi0.M(), weight);

        // fill mass vs lead tower eta
        h_mass_eta_lt->Fill(pi0.M(), lt_eta);
        h_mass_eta_lt_rw->Fill(pi0.M(), lt_eta, weight);

        // fill 3D distribution (mass, pi0 E, tower eta)
        h_m_pt_eta->Fill(pi0.M(), pi0.E(), lt_eta);

        // fill block histogram if valid
        if (block_eta_bin >= 0 && block_eta_bin < NBinsBlock &&
            block_phi_bin >= 0 && block_phi_bin < NBinsBlock)
        {
          h_mass_block_pt[block_eta_bin][block_phi_bin]->Fill(pi0.M(), pi0.E());
        }

        // If both clusters matched meson photons => fill truth-based hist
        if (match2 && pi0_trKin.M() > 0.001)
        {
          // you might choose the first meson photon or apply more logic
          h_m_ptTr_eta->Fill(pi0.M(), truth_meson_photons.at(0).E(), lt_eta);
          h_m_ptTr_eta_trKin->Fill(pi0_trKin.M(), truth_meson_photons.at(0).E(), lt_eta);
        }

      } // end inner loop (clusterIter2)

    } // end outer loop (clusterIter)

    // Optionally print a final summary if Verbosity is set:
    if (Verbosity() > 0)
    {
      // Simple tabular summary:
      // Example of an ASCII table row: we can do more columns if you like
      // or just keep it minimal:
      std::cout << ANSI_BOLD << "\n[Verbosity] finalClusterLoop: Completed analysis of clusters.\n"
                << "    -> # of clusters that passed the QA: " << nClusCount << "\n"
                << "    -> Weighted histograms filled with weight = " << weight << "\n"
                << "[Verbosity] finalClusterLoop: End of function.\n"
                << ANSI_RESET << std::endl;
    }

} // end finalClusterLoop


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
  float clus_chisq_cut = 10000;
  float nClus_ptCut    = 0.5;
  int   max_nClusCount = 3000000;

  // ------------------------------------------------------------------------
  // 1) Retrieve vertex z-position
  // ------------------------------------------------------------------------
  float vtx_z = retrieveVertexZ(topNode);
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
  float ptMaxCut   = 100.f;
  float pt1ClusCut = 0.9f;
  float pt2ClusCut = 0.9f;
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


TF1* PositionDependentCorrection::fitHistogram(TH1* h)
{
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] fitHistogram() called. Setting up the fit function parameters now." << std::endl;
  }

  // Maintain exact functionality
  TF1* fitFunc = new TF1("fitFunc",
                         "[0]/[2]/2.5*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x^2 + [6]*x^3",
                         h->GetXaxis()->GetXmin(),
                         h->GetXaxis()->GetXmax());

  fitFunc->SetParameter(0, 5);
  fitFunc->SetParameter(1, target_pi0_mass);
  fitFunc->SetParameter(2, 0.01);
  fitFunc->SetParameter(3, 0.0);
  fitFunc->SetParameter(4, 0.0);
  fitFunc->SetParameter(5, 0.0);
  fitFunc->SetParameter(6, 100);

  fitFunc->SetParLimits(0, 0, 10);
  fitFunc->SetParLimits(1, 0.113, 0.25);
  fitFunc->SetParLimits(2, 0.01, 0.04);
  fitFunc->SetParLimits(3, -2, 1);
  fitFunc->SetParLimits(4, 0, 40);
  fitFunc->SetParLimits(5, -150, 50);
  fitFunc->SetParLimits(6, 0, 200);

  fitFunc->SetRange(0.05, 0.7);

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Performing fit on histogram '" << (h ? h->GetName() : "null")
              << "' using option 'QN'." << std::endl;
  }

  // Perform the fit
  h->Fit("fitFunc", "QN");

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] fitHistogram() completed the fit. Returning TF1 pointer." << std::endl;
  }
  return fitFunc;
}


void PositionDependentCorrection::fitEtaSlices(const std::string& infile,
                                              const std::string& fitOutFile,
                                              const std::string& cdbFile)
{
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] fitEtaSlices() called with infile='" << infile
              << "', fitOutFile='" << fitOutFile
              << "', cdbFile='" << cdbFile << "'." << std::endl;
  }

  // Maintain exact functionality
  TFile* fin = new TFile(infile.c_str());
  if (!fin)
  {
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] fitEtaSlices(): input TFile pointer is null! Exiting." << std::endl;
    }
    std::cout << "PositionDependentCorrection::fitEtaSlices null fin" << std::endl;
    exit(1);
  }

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Successfully opened input file: " << infile << std::endl;
    std::cout << "[DEBUG] Starting to retrieve histograms for fit..." << std::endl;
  }

  std::cout << "getting hists" << std::endl; // preserving original line

  TH1F* h_peak_eta  = new TH1F("h_peak_eta",  "", 96, 0, 96);
  TH1F* h_sigma_eta = new TH1F("h_sigma_eta", "", 96, 0, 96);
  TH1F* h_p3_eta    = new TH1F("h_p3_eta",    "", 96, 0, 96);
  TH1F* h_p4_eta    = new TH1F("h_p4_eta",    "", 96, 0, 96);
  TH1F* h_p5_eta    = new TH1F("h_p5_eta",    "", 96, 0, 96);
  TH1F* h_p6_eta    = new TH1F("h_p6_eta",    "", 96, 0, 96);
  TH1F* h_p0_eta    = new TH1F("h_p0_eta",    "", 96, 0, 96);

  TH1F* h_M_eta[96];
  for (int i = 0; i < 96; i++)
  {
    h_M_eta[i] = (TH1F*) fin->Get(Form("h_mass_eta_lt_rw%d", i));
    if (h_M_eta[i])
    {
      h_M_eta[i]->Scale(1. / h_M_eta[i]->Integral(), "width");
    }
    else
    {
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG] fitEtaSlices(): histogram 'h_mass_eta_lt_rw" << i
                  << "' does not exist in file. This may cause issues but continuing." << std::endl;
      }
    }
  }

  TF1* fitFunOut[96];
  for (int i = 0; i < 96; i++)
  {
    if (!h_M_eta[i])
    {
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG] PositionDependentCorrection::fitEtaSlices null hist 'h_mass_eta_lt_rw"
                  << i << "'. Continuing." << std::endl;
      }
      std::cout << "PositionDependentCorrection::fitEtaSlices null hist" << std::endl;
      continue;
    }

    fitFunOut[i] = fitHistogram(h_M_eta[i]);
    fitFunOut[i]->SetName(Form("f_pi0_eta%d", i));

    float mass_val_out = fitFunOut[i]->GetParameter(1);
    float mass_err_out = fitFunOut[i]->GetParError(1);
    h_peak_eta->SetBinContent(i + 1, mass_val_out);

    if (isnan(h_M_eta[i]->GetEntries()))
    {
      h_peak_eta->SetBinError(i + 1, 0);
      continue;
    }
    h_peak_eta->SetBinError(i + 1, mass_err_out);

    h_sigma_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(2));
    h_sigma_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(2));

    h_p3_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(3));
    h_p3_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(3));

    h_p4_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(4));
    h_p4_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(4));

    h_p5_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(5));
    h_p5_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(5));

    h_p6_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(6));
    h_p6_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(6));

    h_p0_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(0));
    h_p0_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(0));
  }

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Initial fit loops completed. Proceeding with CDBTTree usage..." << std::endl;
  }

  CDBTTree* cdbttree1 = new CDBTTree(cdbFile.c_str());
  CDBTTree* cdbttree2 = new CDBTTree(cdbFile.c_str());
  std::string m_fieldname = "Femc_datadriven_qm1_correction";

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Updating CDBTTree values with pi0 mass corrections." << std::endl;
  }

  for (int i = 0; i < 96; i++)
  {
    for (int j = 0; j < 256; j++)
    {
      float correction = target_pi0_mass / h_peak_eta->GetBinContent(i + 1);
      unsigned int key = TowerInfoDefs::encode_emcal(i, j);
      float val1 = cdbttree1->GetFloatValue(key, m_fieldname);
      cdbttree2->SetFloatValue(key, m_fieldname, val1 * correction);
    }
  }

  cdbttree2->Commit();
  cdbttree2->WriteCDBTTree();
  delete cdbttree2;
  delete cdbttree1;

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Committed updated corrections to the CDBTTree." << std::endl;
    std::cout << "[DEBUG] Creating output fit file: " << fitOutFile << std::endl;
  }

  TFile* fit_out = new TFile(fitOutFile.c_str(), "recreate");
  fit_out->cd();

  for (auto& i : h_M_eta)
  {
    if (i)
    {
      i->Write();
      delete i;
    }
  }

  for (auto& i : fitFunOut)
  {
    if (i) i->Write();
    delete i;
  }

  h_p3_eta->Write();
  h_p4_eta->Write();
  h_p5_eta->Write();
  h_p6_eta->Write();
  h_p0_eta->Write();
  h_sigma_eta->Write();
  h_peak_eta->Write();

  fin->Close();
  delete fin;

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Fit results and histograms written to " << fitOutFile << std::endl;
  }

  std::cout << "finish fitting suc" << std::endl; // original line maintained

  return;
}


std::pair<float,float> PositionDependentCorrection::getBlockCordCorr(
    std::vector<int> etas,
    std::vector<int> phis,
    std::vector<float> Es)
{
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] getBlockCordCorr() called. Computing raw local block coordinate." << std::endl;
  }

  // Maintain exact functionality
  std::pair<float,float> raw_cord = getBlockCord(etas, phis, Es);

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Raw block coordinate: (" << raw_cord.first
              << ", " << raw_cord.second << ")" << std::endl;
  }

  float Xcg_eta = (raw_cord.first  < 0.5) ? raw_cord.first  : (raw_cord.first  - 1.0);
  float Xcg_phi = (raw_cord.second < 0.5) ? raw_cord.second : (raw_cord.second - 1.0);

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Adjusted local coordinates for asinh domain: Xcg_eta="
              << Xcg_eta << ", Xcg_phi=" << Xcg_phi << std::endl;
  }

  float b_eta = 0.155;

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Using b_eta=" << b_eta << " and b_phi=" << b_phi
              << " for correction." << std::endl;
  }

  TF1 finv("finv", "[0]*TMath::ASinH(2*x*sinh(1/(2*[0])))", -0.5, 0.5);

  // Evaluate for eta direction
  finv.SetParameter(0, b_eta);
  float t_eta = finv.Eval(Xcg_eta);

  // Evaluate for phi direction
  finv.SetParameter(0, b_phi);
  float t_phi = finv.Eval(Xcg_phi);

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Corrected local coordinates: t_eta="
              << t_eta << ", t_phi=" << t_phi << std::endl;
  }

  std::pair<float,float> res;
  res.first  = (raw_cord.first  < 0.5) ? t_eta : (t_eta + 1.0);
  res.second = (raw_cord.second < 0.5) ? t_phi : (t_phi + 1.0);

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Final block coordinate after shift: ("
              << res.first << ", " << res.second << ")" << std::endl;
    std::cout << "[DEBUG] getBlockCordCorr() returning corrected coordinates." << std::endl;
  }

  return res;
}

/////////////////////////////////////////////////////////////////////////////////////
// getBlockCord(...)
// -----------------
// This method calculates the "raw" local block coordinate (eta, phi) within a 2×2
// cluster block based on tower indices. It does this by computing a "weighted average"
// in tower (eta, phi) space and then mapping that to the 0–2 range in each dimension.
// For example, if the cluster is near the middle of one tower, the block coordinate
// might be near (0.5, 0.5). If it's near the boundary of a 2×2 set, it might be
// (1.9, 0.1), etc. Then we shift that range into ~[0,1] by subtracting the block indices.
/////////////////////////////////////////////////////////////////////////////////////

std::pair<float,float> PositionDependentCorrection::getBlockCord(
    std::vector<int> etas,
    std::vector<int> phis,
    std::vector<float> Es)
{
  // Maintain original functionality, adding verbosity prints:
  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] getBlockCord() called. Computing average tower indices..." << std::endl;
  }

  // Get the average tower index in eta and phi directions
  // weighted by tower energies (like an energy-weighted center of gravity).
  float avgEta = getAvgEta(etas, Es);
  float avgPhi = getAvgPhi(phis, Es);

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] Calculated avgEta=" << avgEta
              << ", avgPhi=" << avgPhi << std::endl;
    std::cout << "[DEBUG] Converting to integer bins, dividing by 2 for block indexing..." << std::endl;
  }

  // Convert floating tower indices into integer bins, then find the block index
  // by dividing by 2. This concept lumps towers into 2×2 blocks in the geometry.
  int avgIPhi    = std::floor(avgPhi);
  int avgIEta    = std::floor(avgEta);
  int iblockphi  = avgIPhi  / 2;
  int iblocketa  = avgIEta  / 2;

  // Subtract the block indices times 2 so that we get an offset
  // within that 2×2 block (range 0 to 2).
  float interBlockPhi = avgPhi - iblockphi * 2;
  float interBLockEta = avgEta - iblocketa * 2;

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] interBlockPhi=" << interBlockPhi
              << ", interBLockEta=" << interBLockEta << std::endl;
    std::cout << "[DEBUG] Checking for local eta coordinate above 1.5..." << std::endl;
  }

  // If the local eta coordinate ended up above 1.5, shift it by -2.
  // This is typically to keep it in [0,2) range rather than [2, something bigger].
  if (interBLockEta > 1.5)
  {
    interBLockEta -= 2;
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] interBLockEta was >1.5, shifted by 2 => "
                << interBLockEta << std::endl;
    }
  }

  // Return the final "blockCord" pair, which is the coordinate inside that 2×2 block.
  // Typically, we interpret interBLockEta and interBlockPhi in the range [0,1].
  std::pair<float,float> blockCord = {interBLockEta, interBlockPhi};

  if (Verbosity() > 0)
  {
    std::cout << "[DEBUG] getBlockCord() returning block coordinates: ("
              << blockCord.first << ", " << blockCord.second << ")" << std::endl;
  }
  return blockCord;
}


/////////////////////////////////////////////////////////////////////////////////////
// getAvgEta(...)
// --------------
// Calculates the energy-weighted average of the tower 'eta indices' for all towers in
// the cluster. The index itself is an integer tower bin, but we treat it as a float
// when we do the weighting. This effectively returns the cluster's position in tower
// eta-bin space.
/////////////////////////////////////////////////////////////////////////////////////

float PositionDependentCorrection::getAvgEta(
    const std::vector<int> &toweretas,
    const std::vector<float> &towerenergies)
{
    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] getAvgEta() called with toweretas.size()="
                << toweretas.size() << std::endl;
    }

    float etamult = 0;
    float etasum  = 0;

    // If only one tower is present, simply return its integer index
    // without weighting.
    if (toweretas.size() == 1)
    {
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG] Only one tower in cluster => returning that tower's eta index: "
                  << toweretas[0] << std::endl;
      }
      return toweretas[0];
    }

    // Otherwise, sum up eta_index × energy for each tower, and also sum energies.
    for (UInt_t j = 0; j < towerenergies.size(); j++)
    {
        float energymult = towerenergies[j] * toweretas[j];
        etamult += energymult;
        etasum  += towerenergies[j];
    }

    // The ratio gives the "average" tower-eta bin for this cluster.
    float result = etamult / etasum;

    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] getAvgEta() returning "
                << result << " as the energy-weighted average eta-index."
                << std::endl;
    }
    return result;
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

    // Additional offset: we then shift by +0.5 and apply modulo 2, then subtract 0.5.
    // This is a final step to keep the average phi in some narrower block coordinate range
    // for local tower indexing or cluster block identification. The details of why 0.5 is used
    // can be geometry-specific, ensuring the center of a tower is at .5, for instance.
    avgphi = fmod(avgphi + 0.5, 2) - 0.5;

    if (Verbosity() > 0)
    {
      std::cout << "[DEBUG] getAvgPhi() returning final adjusted phi="
                << avgphi << std::endl;
    }

    // Return the final adjusted phi coordinate in tower-bin space.
    return avgphi;
}
