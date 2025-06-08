#include "PositionDependentCorrection.h"

// G4Hits includes
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
#include <calobase/RawTowerDefs.h>

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

    // ------------------------------------------------------------------
    // (A)  Optional one-time load of surveyed tower geometry
    // ------------------------------------------------------------------
    if (m_useSurveyGeometry)
    {
      // ────────────────────────────────────────────────────────────────
      // 0)  Banner
      // ────────────────────────────────────────────────────────────────
      if (Verbosity() > 0)
        std::cout << ANSI_BOLD << ANSI_CYAN
                  << "[PDC]  Loading CDB CALO_TOWER_GEOMETRY  "
                  << "(tag = " << m_cdbTag << ",  ts = " << m_timeStamp << ")"
                  << ANSI_RESET << '\n';

      /* make CDB happy */
      recoConsts* rc = recoConsts::instance();
      rc->set_StringFlag("CDB_GLOBALTAG", m_cdbTag);
      rc->set_uint64Flag("TIMESTAMP",     m_timeStamp);

      // ────────────────────────────────────────────────────────────────
      // 1)  Download + open the tree
      // ────────────────────────────────────────────────────────────────
      const std::string url =
          CDBInterface::instance()->getUrl("CALO_TOWER_GEOMETRY");

      std::unique_ptr<CDBTTree> geomTree(new CDBTTree(url));
      geomTree->LoadCalibrations();                    // <— fills the tree

      if (Verbosity() > 0)
        std::cout << ANSI_GREEN
                  << "      … geometry loaded from\n      " << url
                  << ANSI_RESET << "\n";

      // ────────────────────────────────────────────────────────────────
      // 2)  Extract the rigid φ-offset and print a concise table
      // ────────────────────────────────────────────────────────────────
      if (Verbosity() > 0)
      {
        constexpr int kNPhi = 256;                     // full fine bins
        std::cout << ANSI_BOLD << "[φ-survey summary]" << ANSI_RESET << '\n';

        /* (a) shift = ideal – surveyed of bin 0 */
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

        /* (b) small header */
        std::cout << ANSI_CYAN
                  << "      bin   φ_first [rad]   φ_second [rad]\n"
                  << "      ─────────────────────────────────────\n"
                  << ANSI_RESET;

        /* (c) print either all bins (verbosity>1) or the first few */
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
    
    
  // ────────────────────────────────────────────────────────────────────────
  // χ² acceptance maps (tower η × φ)
  // ────────────────────────────────────────────────────────────────────────
  h2_chi2_tot_etaPhi = new TH2F(
        "h2_chi2_tot_etaPhi",
        "Clusters BEFORE #chi^{2} cut;Tower #eta index;Tower #varphi index;Counts",
        96, 0, 96,        // η  (same binning you already use)
        256, 0, 256);     // φ

  h2_chi2_rej_etaPhi = new TH2F(
        "h2_chi2_rej_etaPhi",
        "Clusters REJECTED by #chi^{2} cut;Tower #eta index;Tower #varphi index;Counts",
        96, 0, 96,
        256, 0, 256);

  /* profile: fill value = 1 if the cluster PASSES the χ² cut, 0 otherwise.
     * The bin mean therefore equals the pass-fraction; 1−mean is the reject-fraction
  */
  p_chi2_pass_etaPhi = new TProfile2D(
        "p_chi2_pass_etaPhi",
        "Pass fraction after #chi^{2} cut;Tower #eta index;Tower #varphi index;⟨pass⟩",
        96, 0, 96,
        256, 0, 256,      // low/high z are ignored in profile constructor
        -0.1, 1.1);       // y-range just aesthetics

  hm->registerHisto(h2_chi2_tot_etaPhi);
  hm->registerHisto(h2_chi2_rej_etaPhi);
  hm->registerHisto(p_chi2_pass_etaPhi);


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
    //  Ash-b   and   Log-w0  trial grids  (energy-binned)
    //  ---------------------------------------------------------------------------
    //  • ***b is now expressed in pure “tower-width units”***
    //      (1.0 = exactly one 5.55-cm tower pitch).
    //  • Δx is still measured in centimetres, so the histogram axis code is
    //    unchanged.
    // ──────────────────────────────────────────────────────────────────────────────
    {
      // Histogram window for Δx  (= xreco − xtrue)
      constexpr double DX_MAX = 12.0;        // [cm]
      constexpr double BIN_W  = 0.10;        // [cm]
      constexpr int    NBINS  = int( 2*DX_MAX / BIN_W + 0.5 );   // → 240 bins

      // 1)  Ash scan  –  b grid in **tower units**  (0.01 … 0.50, step 0.01)
      m_bScan.clear();
      for (double b = 0.01; b <= 0.50 + 1e-9; b += 0.01)   // skip 0.00
          m_bScan.push_back( std::round(b * 10000.) / 10000. );   // 4-dec precision

      // ────────────────────────────────────────────────────────────────────────
      // 2)  Log-weight scan  –  w0 grid (1.5 … 7.0, step 0.10)
      // ────────────────────────────────────────────────────────────────────────
      m_w0Scan.clear();
      for (double w0 = 1.5; w0 <= 7.0 + 1e-9; w0 += 0.10)
        m_w0Scan.push_back( std::round(w0*100.)/100. );       // 2-dec precision

      // ────────────────────────────────────────────────────────────────────────
      // 3)  Histogram booking (done **once**)
      // ────────────────────────────────────────────────────────────────────────
      if (!alreadyDeclaredHistograms)
      {
        for (int iE = 0; iE < N_Ebins; ++iE)            // energy-slice loop
        {
          // 3a)  Ash-b histograms
          for (double bVal : m_bScan)
          {
            const TString hName = Form("h_dx_ash_b%.4f_E%d", bVal, iE);
            auto* h = new TH1F( hName,
                                ";x_{reco}-x_{true}  [cm];Counts / 0.1 cm",
                                NBINS, -DX_MAX, +DX_MAX );
            hm->registerHisto(h);
          }

          // 3b)  Log-w0 histograms
          for (double w0 : m_w0Scan)
          {
            const TString hName = Form("h_dx_log_w0%.2f_E%d", w0, iE);
            auto* h = new TH1F( hName,
                                ";x_{reco}-x_{true}  [cm];Counts / 0.1 cm",
                                NBINS, -DX_MAX, +DX_MAX );
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

//═══════════════════════════════════════════════════════════════════════
//  retrieveVertexZ()
//───────────────────────────────────────────────────────────────────────
//  • Returns the best-guess event-vertex z (cm).
//  • If `getVtx` is false we immediately fall back to 0.
//  • Every failure point prints a colour-coded message when
//    Verbosity() > 0.
//═══════════════════════════════════════════════════════════════════════
float
PositionDependentCorrection::retrieveVertexZ(PHCompositeNode* topNode)
{
  /*-----------------------------------------------------------------*/
  /* 0) function entry                                               */
  /*-----------------------------------------------------------------*/
  if (Verbosity() > 0)
  {
    std::cout << '\n'
              << ANSI_BOLD << ANSI_CYAN << "[retrieveVertexZ] ENTER"
              << ANSI_RESET << '\n'
              << "  getVtx flag ........... ";
    if (getVtx)
      std::cout << ANSI_GREEN << "true";
    else
      std::cout << ANSI_YELLOW << "false";
    std::cout << ANSI_RESET << '\n';
  }

  /* default value if vertex cannot be obtained */
  float vtx_z = 0.f;

  /*-----------------------------------------------------------------*/
  /* 1) short-circuit if user disabled vertex reading                */
  /*-----------------------------------------------------------------*/
  if (!getVtx)
  {
    if (Verbosity() > 0)
      std::cout << ANSI_YELLOW
                << "  getVtx = false  →  using vtx_z = 0\n"
                << ANSI_RESET;
    return vtx_z;
  }

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



// --------------------------------------------------------------------------
//  PositionDependentCorrection
//  Inverse Ash-b correction  (measured  →  true)               2025-06-05
// --------------------------------------------------------------------------
//  Prerequisites:  #include <cmath>     // std::sinh, std::asinh, std::fmod
// --------------------------------------------------------------------------

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

  /* —————————————————— exit —————————————————— */
  if (Verbosity() > 0)
    std::cout << ANSI_GREEN << "[doPhiBlockCorr] EXIT  → "
              << corrected << " (rad-fraction)"
              << ANSI_RESET << "\n\n";

  return corrected;
}

/*==========================================================================*
 *  η-direction  (carbon-copy of φ version)                                 *
 *==========================================================================*/
float PositionDependentCorrection::doEtaBlockCorr(float localEta, float b)
{
  if (Verbosity() > 0)
    std::cout << ANSI_BOLD << "[doEtaBlockCorr] ENTER  "
              << "localEta = " << localEta << " ,  b = " << b
              << ANSI_RESET << '\n';

  if (std::fabs(b) < 1e-9f) {
    if (Verbosity() > 0) std::cout << "  b≈0 → return unchanged\n";
    return localEta;
  }

  if (localEta <= -0.5f || localEta > 1.5f) {
    if (Verbosity() > 1) std::cout << "  fold: " << localEta << " → ";
    localEta = std::fmod(localEta + 2.f, 2.f);
    if (Verbosity() > 1) std::cout << localEta << '\n';
  }

  const float Xmeas = (localEta < 0.5f) ? localEta : localEta - 1.f;
  if (Verbosity() > 1) std::cout << "  Xmeas = " << Xmeas << '\n';

  const double S = std::sinh(1.0 / (2.0 * b));
  if (Verbosity() > 1) std::cout << "  S = sinh(1/2b) = " << S << '\n';

  const float Xtrue = static_cast<float>( b * std::asinh( 2.0 * S * Xmeas ) );
  if (Verbosity() > 1) std::cout << "  Xtrue = " << Xtrue << '\n';

  float corrected = (localEta < 0.5f) ? Xtrue : Xtrue + 1.f;
  if (Verbosity() > 1)
    std::cout << "  corrected (pre-wrap) = " << corrected << '\n';

  if (corrected <= -0.5f || corrected > 1.5f) {
    corrected = std::fmod(corrected + 2.f, 2.f);
    if (corrected > 1.5f) corrected -= 2.f;
    if (Verbosity() > 1)
      std::cout << "  corrected (post-wrap) = " << corrected << '\n';
  }

  if (Verbosity() > 0)
    std::cout << ANSI_GREEN << "[doEtaBlockCorr] EXIT  → "
              << corrected << ANSI_RESET << "\n\n";

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

// ============================================================================
//  convertBlockToGlobalPhi – (blk,loc) → absolute EMCAL φ in (−π , +π]
// ----------------------------------------------------------------------------
//  • Accepts any real local coordinate; performs **exactly one** fold into
//    the canonical domain (−0.5 … +1.5] without ever touching the coarse index.
//  • Returns NaN if the coarse index is illegal or if the folded local
//    coordinate is still out of bounds (should never happen).
// ============================================================================

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
  const float fullPhiIndex =
        static_cast<float>(block_phi_bin) * kFinePerBlock +  // coarse origin
        localPhi +                                           // intra-block
        0.5f;                                                // centre of fine bin

  if (Verbosity() > 1)
    std::cout << "    ▸ fullPhiIndex = " << fullPhiIndex
              << "   (coarse*2 + local + 0.5)\n";

  /* ─── linear → angular & wrap to (−π , +π] ──────────────────────────── */
  float globalPhi = fullPhiIndex * kRadPerBin;  // ideal grid
//  if (m_hasOffset) globalPhi += m_phi0Offset;   // shift to real detector
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


//============================================================================
//  φ at shower depth
//----------------------------------------------------------------------------
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

//============================================================================
//  x at shower depth (cm)
//----------------------------------------------------------------------------
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



//////////////////////////////////////////////////////////////////////
// Forward distortion  ( true  →  measured )
//////////////////////////////////////////////////////////////////////
float PositionDependentCorrection::doAshShift(float localPhi, float b)
{
  /* 0) no distortion when b ≈ 0 ------------------------------------ */
  if (std::fabs(b) < 1e-9f) return localPhi;

  /* 1) one fold into (-0.5 … +1.5] --------------------------------- */
  if (localPhi <= -0.5f || localPhi > 1.5f)
      localPhi = std::fmod(localPhi + 2.f, 2.f);

  /* 2) map to centre-of-gravity X ∈ (-0.5 … +0.5] ------------------ */
  const float Xcg = (localPhi < 0.5f) ? localPhi : localPhi - 1.f;

  /* 3) forward Ash-b ---------------------------------------------- */
  const double S  = std::sinh(1.0 / (2.0 * b));                // sinh(1/2b)
  float t        = static_cast<float>( b * std::asinh( 2.0 * Xcg * S ) );

  /* 4) shift back to (-0.5 … +1.5] --------------------------------- */
  if (localPhi >= 0.5f) t += 1.f;

  /* 5) final safety fold (numerical guard) ------------------------- */
  if (t <= -0.5f || t > 1.5f)
  {
    t = std::fmod(t + 2.f, 2.f);
    if (t > 1.5f) t -= 2.f;
  }
  return t;
}

////////////////////////////////////////////////////////////////////////
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
//
//   NOTE:
//     Verbosity-gated printouts (Verbosity()>0) give a concise trace of
//     the internal steps without flooding the log.
//
//   2025-06-03  --  unified with Ash implementation
////////////////////////////////////////////////////////////////////////
float
PositionDependentCorrection::doLogWeightCoord( const std::vector<int>&   towerPhi,
                                               const std::vector<float>& towerE,
                                               float                    w0 )
{
  /*──────────────────────────────
   * 0) Sanity & early exits
   *──────────────────────────────*/
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
 *  fillAshLogDx(...)  –  centimetre-accurate version
 ************************************************************************/
void PositionDependentCorrection::fillAshLogDx(
        RawCluster*                     cluster,
        const TLorentzVector&           recoPhoton,
        const TLorentzVector&           truthPhoton,
        const std::pair<float,float>&   blockCord,   // (ηloc , φloc)
        int                             blockPhiBin, // 0 … 127
        const std::vector<int>&         tower_phis,
        const std::vector<float>&       tower_energies )
{
  /* 0) energy slice --------------------------------------------------- */
  const float eReco = recoPhoton.E();
  int iEslice = -1;
  for (int i = 0; i < N_Ebins; ++i)
    if (eReco >= eEdge[i] && eReco < eEdge[i+1]) { iEslice = i; break; }
  if (iEslice < 0) return;

  /* 1) guards & geometry --------------------------------------------- */
  if (!cluster || !m_geometry || !m_bemcRec) return;

  float eMax = -1.f;  RawCluster::TowerConstIterator best{};
  for (auto it = cluster->get_towers().first;
            it != cluster->get_towers().second; ++it)
    if (it->second > eMax) { eMax = it->second; best = it; }
  if (eMax < 1e-6) return;
    

  auto* leadGeo = m_geometry->get_tower_geometry(best->first);
  if (!leadGeo) return;

  const double rFront = leadGeo->get_center_radius();
  const double zFront = leadGeo->get_center_z();
  int ixLead = leadGeo->get_binphi();   // 0 … 255   (column / φ)
  int iyLead = leadGeo->get_bineta();   // 0 … 95    (row / η)
    
    
  /* 2) truth reference ----------------------------------------------- */
    const float  phiSDtruth = phiAtShowerDepth(truthPhoton.E(),
                                               rFront, zFront,
                                               truthPhoton.Phi(),
                                               ixLead, iyLead);
    
    const double xTrue = xAtShowerDepth(truthPhoton.E(),
                                        rFront, zFront,
                                        truthPhoton.Phi(),
                                        ixLead, iyLead);

  if (Verbosity() > 1)
    std::cout << ANSI_CYAN
              << "[fillAshLogDx] slice=" << iEslice
              << "  E=" << eReco << " GeV"
              << "  φSDtruth=" << phiSDtruth
              << "  → xTrue="  << xTrue  << " cm"
              << ANSI_RESET << '\n';

  /* 3) local helper --------------------------------------------------- */
  auto tryFill = [&](float        localPhi,
                     const TString& hName,
                     int&         okCtr,
                     int&         missCtr)
  {
    if (localPhi < 0.f || localPhi >= 1.f) return;

    const float  phiFront = convertBlockToGlobalPhi(blockPhiBin, localPhi);
    const double xReco = xAtShowerDepth(eReco, rFront, zFront,
                                          phiFront, ixLead, iyLead);

    if (auto* h = dynamic_cast<TH1F*>(hm->getHisto(hName.Data())))
    {
      h->Fill(xReco - xTrue);
      ++okCtr;
      if (Verbosity() > 2)
        std::cout << "    " << hName
                  << "  φloc=" << localPhi
                  << "  Δx="   << xReco - xTrue
                  << "  ► filled\n";
    }
    else
    {
      ++missCtr;
      if (Verbosity() > 2)
        std::cout << "    " << hName << "  ► missing\n";
    }
  };

  /* 4) Ash-b scan ----------------------------------------------------- */
  int nAshOK = 0, nAshMiss = 0;
  if (blockCord.second >= 0.f && blockCord.second < 1.f)
    for (double b : m_bScan)
    {
      const double bVal   = std::round(b * 10000.) / 10000.;
      const float  phiLoc = doAshShift(blockCord.second, bVal);
      tryFill(phiLoc,
              Form("h_dx_ash_b%.4f_E%d", bVal, iEslice),
              nAshOK, nAshMiss);
    }

  /* 5) Log-w0 scan ---------------------------------------------------- */
  int nLogOK = 0, nLogMiss = 0;
  if (!tower_phis.empty() && !tower_energies.empty())
    for (double w0 : m_w0Scan)
    {
      const double w0Val  = std::round(w0 * 100.) / 100.;
      const float  phiLoc = doLogWeightCoord(tower_phis, tower_energies,
                                             w0Val);
      tryFill(phiLoc,
              Form("h_dx_log_w0%.2f_E%d", w0Val, iEslice),
              nLogOK, nLogMiss);
    }

  /* 6) summary -------------------------------------------------------- */
  if (Verbosity() > 0)
    std::cout << "  ➜ Ash filled/miss = "
              << nAshOK << '/' << nAshMiss
              << "   |  Log filled/miss = "
              << nLogOK << '/' << nLogMiss << '\n';
}



// ==========================================================================
//  fillDPhiRawAndCorrected
// --------------------------------------------------------------------------
//  Computes Δφ between the reconstructed photon and its truth partner
//  – first with *raw* local-block coordinates,
//  – then with the Ash-b–corrected coordinates (if a fit exists).
//
//  ► All φ’s are propagated to **shower depth** so they can be compared
//    1-to-1 with the Δx study.
//
//  Verbosity levels
//  ----------------
//    0 : completely silent (fast production)
//    1 : one-line banner, per-cluster RAW / CORR print
//    2 : geometry details, tower radii, truth φ, warnings
//    3 : everything above  +  individual diagnostic lines
//
//  QA snapshots are printed after               N = 2, 5, 10  clusters,
//  then every 50th cluster afterwards.
// ==========================================================================
void PositionDependentCorrection::fillDPhiRawAndCorrected(
        RawCluster*                  cluster,
        const TLorentzVector&        recoPhoton,
        const TLorentzVector&        truthPhoton,
        const std::pair<float,float>& blkCoord, // <-- renamed
        int                          blockPhiBin,     // 0 … 127  (coarse)
        float /*rawDelPhi – ignored*/ )
{
  /* ────────────────────────  static per-job counters  ───────────────────── */
  static std::uint64_t nSeen      = 0;     // clusters processed
  static std::uint64_t nCorrUsed  = 0;     // clusters where b-fit existed
  static double        sumAbsRaw  = 0.;    // Σ|Δφ_raw|
  static double        sumAbsCorr = 0.;    // Σ|Δφ_corr|

  const bool v1 = Verbosity() > 0;
  const bool v2 = Verbosity() > 1;
  /* ────────────────────────   0) banner & fast sanity   ─────────────────── */
  if (v1)
  {
    std::cout << '\n' << ANSI_BOLD << ANSI_CYAN
              << "[fillDPhiRawAndCorrected] cluster #" << nSeen + 1
              << ANSI_RESET << "\n"
              << "    • blockCoord = (ηloc=" << blkCoord.first
              << " , φloc=" << blkCoord.second << ")\n"
              << "    • blockPhiBin (coarse) = " << blockPhiBin << '\n'
              << "    • E(reco)  = " << recoPhoton.E()
              << " GeV   |   E(truth) = " << truthPhoton.E() << " GeV\n";
  }

  if (!cluster || !m_geometry || !m_bemcRec)
  {
    if (v1) std::cout << ANSI_RED
      << "  ✘ Missing cluster or geometry pointers – SKIP\n"
      << ANSI_RESET;
    return;
  }

  /* ────────────────────────   1) energy slice id   ──────────────────────── */
  const float eReco = recoPhoton.E();
  int iSlice = -1;
  for (int i = 0; i < N_Ebins; ++i)
    if (eReco >= eEdge[i] && eReco < eEdge[i+1]) { iSlice = i; break; }

  if (iSlice < 0)
  {
    if (v1) std::cout << ANSI_YELLOW
      << "  → reco-E outside slice table – SKIP\n" << ANSI_RESET;
    return;
  }
  if (v2)
    std::cout << "    • slice = " << iSlice
              << "  [" << eEdge[iSlice] << " – " << eEdge[iSlice+1] << ")\n";

  /* ────────────────────────   2) lead-tower geometry   ──────────────────── */
  RawCluster::TowerConstIterator lead{};
  float eLead = -1.f;
  for (auto it = cluster->get_towers().first;
            it != cluster->get_towers().second; ++it)
    if (it->second > eLead) { eLead = it->second; lead = it; }

  if (eLead < 1e-6) { if (v1) std::cout << "  ΣE≈0 – SKIP\n"; return; }
    
  const RawTowerGeom* geo = m_geometry->get_tower_geometry( lead->first );
  if (!geo) { if (v1) std::cout << "  ✘ lead-tower geom missing – SKIP\n"; return; }

  const double rFront = geo->get_center_radius();
  const double zFront = geo->get_center_z();
    
  int ixLead = geo->get_binphi();   // 0 … 255   (column / φ)
  int iyLead = geo->get_bineta();   // 0 … 95    (row / η)
    
  if (v2) std::cout << "    • rFront = " << rFront << " cm , zFront = "
                    << zFront  << " cm\n";

  /* ────────────────────────   3) truth φ at shower depth  ──────────────── */
    const float phiSDtruth =
      phiAtShowerDepth(truthPhoton.E(), rFront, zFront,
                       truthPhoton.Phi(),        // unchanged
                       ixLead, iyLead);          // NEW args


  /* ────────────────────────   4) RAW variant   ──────────────────────────── */
  const float rawPhiFront  =
    convertBlockToGlobalPhi( blockPhiBin, blkCoord.second );
    
    const float phiSDrecoRaw =
      phiAtShowerDepth(eReco, rFront, zFront,
                       rawPhiFront, ixLead, iyLead);
    
  const float dPhiRaw =
      TVector2::Phi_mpi_pi( phiSDrecoRaw - phiSDtruth );

  if (h_phi_diff_raw_E[iSlice]) h_phi_diff_raw_E[iSlice]->Fill( dPhiRaw );
  if (pr_phi_vs_blockcoord) pr_phi_vs_blockcoord->Fill(blkCoord.second, dPhiRaw);


  if (v1) std::cout << ANSI_YELLOW
    << "    RAW  φFront=" << rawPhiFront
    << "  φSD(reco)=" << phiSDrecoRaw
    << "  Δφ_raw=" << dPhiRaw << ANSI_RESET << '\n';

  /* ────────────────────────   5) corrected variant (if fit)  ───────────── */
  ++nSeen;
  sumAbsRaw += std::fabs(dPhiRaw);

  if (!isFitDoneForPhi || m_bValsPhi[iSlice] <= 0.f)
  {
    if (v2) std::cout << "    (no φ-fit for this slice)\n";
  }
  else
  {
    const float locCorr =
        doPhiBlockCorr( blkCoord.second, m_bValsPhi[iSlice] );

    const float corrPhiFront =
        convertBlockToGlobalPhi( blockPhiBin, locCorr );
      
      const float phiSDrecoCorr =
        phiAtShowerDepth(eReco, rFront, zFront,
                         corrPhiFront, ixLead, iyLead);
      
    const float dPhiCorr =
        TVector2::Phi_mpi_pi( phiSDrecoCorr - phiSDtruth );

    if (h_phi_diff_corrected_E[iSlice])
      h_phi_diff_corrected_E[iSlice]->Fill( dPhiCorr );
    if (h_localPhi_corrected[iSlice])
      h_localPhi_corrected[iSlice]->Fill( locCorr );

    ++nCorrUsed;
    sumAbsCorr += std::fabs(dPhiCorr);

    if (v1) std::cout << ANSI_GREEN
      << "    CORR b=" << m_bValsPhi[iSlice]
      << "  φloc_corr=" << locCorr
      << "  Δφ_corr=" << dPhiCorr << ANSI_RESET << '\n';
  }

  /* ────────────────────────   6) QA snapshots   ─────────────────────────── */
  auto printSnapshot = [&](const char* tag)
  {
    const double meanRaw = sumAbsRaw  / static_cast<double>(nSeen);
    const double meanCor = nCorrUsed
            ? sumAbsCorr / static_cast<double>(nCorrUsed) : 0.;

    std::cout << ANSI_BOLD
              << "\n──────────  φ-QA  (" << tag << ")  ──────────\n"
              << ANSI_RESET
              << "    clusters processed  : " << nSeen << '\n'
              << "    ⟨|Δφ_raw|⟩         : " << meanRaw << " rad\n";
      if (nCorrUsed)
        std::cout << "    ⟨|Δφ_corr|⟩        : " << meanCor << " rad  ("
                  << (meanCor < meanRaw ? ANSI_GREEN : ANSI_RED)     // colour
                  << (meanCor < meanRaw ? "improvement" : "worse")   // word
                  << ANSI_RESET << ")\n";   // ← add this “)\n” and keep the semicolon
      else
        std::cout << "    (no corrected clusters yet)\n";
    std::cout << "───────────────────────────────────────────────\n\n";
  };

  if (nSeen==2 || nSeen==5 || nSeen==10 || nSeen%50==0) printSnapshot("checkpoint");
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
      
    // ---- χ² QA (denominator) ----
    h2_chi2_tot_etaPhi->Fill(lt_eta, lt_phi);

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
      constexpr float kChi2Cut = 10.0f;   // keep the magic number centralised

      if (clusChi2 > kChi2Cut)
      {
        // ---- χ² QA (numerator & profile) ----
        h2_chi2_rej_etaPhi->Fill(lt_eta, lt_phi);       // rejected count
        p_chi2_pass_etaPhi ->Fill(lt_eta, lt_phi, 0.0); // pass-flag = 0
        if (Verbosity() > 0)
          std::cout << "[χ²-CUT]  cluster Chi2=" << clusChi2
                    << "  (lead η,φ)=" << lt_eta << ',' << lt_phi
                    << "  REJECTED\n";
        continue;                                       // skip the cluster
      }
      else
      {
        // cluster is accepted
        p_chi2_pass_etaPhi->Fill(lt_eta, lt_phi, 1.0);  // pass-flag = 1
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

      int   blkPhiCoarse = -1;             // will be filled by getBlockCord
      auto  blkCoord     = getBlockCord(towerEtas,
                                        towerPhis,
                                        towerEs,
                                        blkPhiCoarse);   // ← new 4-th argument

      /* 4) find energy slice ........................................... */
      int iEbin = -1;
      for (int i = 0; i < N_Ebins; ++i)
        if (clusE >= eEdge[i] && clusE < eEdge[i+1]) { iEbin = i; break; }

      /* optional detailed print ........................................ */
      if (Verbosity() > 0)
      {
        std::cout << "[DEBUG] finalClusterLoop →"
                  << "  blockCord=(" << blkCoord.first << ',' << blkCoord.second << ")"
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
        h3_cluster_block_cord_E->Fill(blkCoord.first,   // η-coord   (x-axis)
                                      blkCoord.second,  // φ-coord   (y-axis)
                                      clusE,             // energy    (z-axis)
                                      kFillW);

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
      }

      // -------------------------------------------------------------
      // 2)  CORRECTED local-block coordinates  (η and/or φ b-shift)
      //     – always fill so raw & corrected have the same statistics
      // -------------------------------------------------------------
      {
        /* ---- original raw local coordinates ---------------- */
        const float rawEta = blkCoord.first;    // [0,1]
        const float rawPhi = blkCoord.second;   // [0,1]

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
      for (auto& trPhoton : truth_photons)
      {
          /* ----------------------------------------------------------
           * 1) fast pre-selection — skip non-matches immediately
           * -------------------------------------------------------- */
          const float dR    = photon1.DeltaR(trPhoton);          // ΔR(reco,truth)
          const float ratio = photon1.E() / trPhoton.E();        // Ereco / Etruth

          if (dR   > 0.03f)                    continue;  // outside match cone
          if (ratio < 0.30f || ratio > 1.30f)  continue;  // wrong energy scale

          /* ----------------------------------------------------------
           * 2) this is the matched truth photon — safe to proceed
           * -------------------------------------------------------- */
          const float dPhi = TVector2::Phi_mpi_pi(
                               photon1.Phi() - trPhoton.Phi() ); // branch-safe

          /* QA that should **only** use the matched pair             */
          h_delR_recTrth->Fill(dR);

          if (Verbosity() > 0)
          {
              std::cout << "[DEBUG] => cluster1 E=" << photon1.E()
                        << " matched to truthE="    << trPhoton.E()
                        << "  dR="   << dR
                        << "  dPhi=" << dPhi
                        << "  ratio="<< ratio << '\n';
          }

          /* ----------------------------------------------------------
           * 3) original “matched” content — unchanged
           * -------------------------------------------------------- */
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

          /* position-dependent studies */
          fillAshLogDx(clus1, photon1, trPhoton,
                       blkCoord, blkPhiCoarse,
                       towerPhis, towerEs);

          fillDPhiRawAndCorrected(clus1, photon1, trPhoton,
                                  blkCoord, blkPhiCoarse, dPhi);

          if (block_eta_bin >= 0 && block_eta_bin < NBinsBlock &&
              block_phi_bin >= 0 && block_phi_bin < NBinsBlock)
          {
              h_res_block_E[block_eta_bin][block_phi_bin]
                  ->Fill(ratio, trPhoton.E());
          }

          /* optional: break after first successful match to save CPU */
          // break;
    }   // end loop truth_photons

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

    // ────────────────────────────────────────────────────────────────────────────────
    // 4) Retrieve calorimeter geometry *once* and measure the rigid φ-offset between
    //    the ideal tower grid and the real detector.
    // ────────────────────────────────────────────────────────────────────────────────
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

        const int   iphi = tg->get_binphi();              // 0 … nPhiBins-1
        const float phiIdeal = (iphi + 0.5F)*kRadPerBin;   // <<< NO “–π” here
        const float delta    = TVector2::Phi_mpi_pi( tg->get_phi() - phiIdeal );
        sumSin += std::sin(delta);
        sumCos += std::cos(delta);
        ++nTowers;

        /* optional: print a few towers for sanity */
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

      /* ------------------------------------------------------------------
       * 4. Final banner
       * ----------------------------------------------------------------*/
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
                                          int                      &blockPhiBinOut)
{
  [[maybe_unused]] constexpr int kNEtaFine        =  96;   // mapping only
  [[maybe_unused]] constexpr int kNPhiFine        = 256;
  constexpr int      kFinePerBlock   = 2;                 // 2×2 super-cell
  constexpr int      kNCoarseBlocks  = kNPhiFine / 2;     // 128

  const std::size_t nT = towerEs.size();

  if (Verbosity() > 0)
    std::cout << ANSI_CYAN << "[getBlockCord] ENTER  nTowers = "
              << nT << ANSI_RESET << '\n';

  /* ─── sanity ─────────────────────────────────────────────────────────── */
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

  /* ─── 4) single symmetric fold helper ───────────────────────────────── */
  auto foldOnce = [&](float &loc, int &coarse, const char *tag)
  {
    if (loc <= -0.5f || loc > 1.5f)
    {
      const float before = loc;
      loc = std::fmod(loc + 2.0f, 2.0f);           // (0 … 2)
      if (loc > 1.5f) { loc -= 2.0f; ++coarse; }   // shift block by +1

      if (coarse == kNCoarseBlocks) coarse = 0;    // φ wrap-around

      if (Verbosity() > 0)
        std::cout << ANSI_YELLOW
                  << "    • " << tag
                  << " folded: " << before << " → " << loc
                  << "  |  blk+1 → " << coarse
                  << ANSI_RESET << '\n';
    }
  };

  foldOnce(locEta, blkEta, "η");   // kept for completeness
  foldOnce(locPhi, blkPhi, "φ");

  /* ─── 5) propagate the *updated* block index and return ─────────────── */
  blockPhiBinOut = blkPhi;

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

// ============================================================================
//  getAvgPhi – energy-weighted φ–centroid  (0 ≤ ⟨φ⟩ < 256)
// ----------------------------------------------------------------------------
//  • Each tower is specified by its *fine* φ-bin (0 … 255) and energy (GeV).
//  • The routine performs **exactly one** wrap relative to the max-energy
//    reference tower so that every adjusted tower lies within ±128 bins
//    of the reference.  No further modulo arithmetic is applied later.
//  • If the caller supplies ill-formed input (size mismatch, ΣE≈0, …) the
//    function degrades gracefully and returns 0 (with an optional warning).
// ============================================================================

float
PositionDependentCorrection::getAvgPhi(const std::vector<int>   &towerPhis,
                                       const std::vector<float> &towerEs)
{
  constexpr int   kNPhiFine  = 256;                // total fine bins
  constexpr float kHalfSpan  = kNPhiFine / 2.0f;   // 128.0
  const std::size_t nT       = towerPhis.size();

  /* ─── sanity ─────────────────────────────────────────────────────────── */
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

  /* ─── 4) diagnostics ─────────────────────────────────────────────────── */
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
