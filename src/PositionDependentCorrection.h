//  ════════════════════════════════════════════════════════════════════
//  PositionDependentCorrection.h
//
//  – Public interface and implementation details are **unchanged**.
//  – Ready to replace the original header and compile without edits.
//  ════════════════════════════════════════════════════════════════════
#pragma once
#ifndef PositionDependentCorrection_H__
#define PositionDependentCorrection_H__

/* =================================================================== *
 *  1.  Standard / sPHENIX includes                                     *
 * =================================================================== */
#include <functional>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <CLHEP/Vector/ThreeVector.h>   // Hep3Vector in inline helpers
#include <fun4all/SubsysReco.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <TString.h>
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/BEmcRec.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/BEmcRecCEMC.h"


/* =================================================================== *
 *  2.  Forward declarations (ROOT & sPHENIX)                           *
 * =================================================================== */
class Fun4AllHistoManager;
class PHCompositeNode;
class RawTowerGeomContainer;
class RawClusterContainer;
class BEmcRecCEMC;
class RawCluster;

/* ROOT */
class TFile;  class TNtuple; class TTree;
class TH1;    class TH1F;    class TH2;    class TH2F; class TH3;    class TH3F;
class TF1;    class TProfile;class TProfile2D;
class TLorentzVector; class TRandom3;

/* =================================================================== *
 *  3.  PositionDependentCorrection class                               *
 * =================================================================== */
class PositionDependentCorrection : public SubsysReco
{
 public:
  inline static std::atomic<std::uint64_t> s_verbosityLevel{0};
  static uint64_t getVerbosityLevel() { return s_verbosityLevel.load(); }
  /* ----------------------------------------------------------------
   *  Constructors, destructors, Fun4All hooks
   * -------------------------------------------------------------- */
  explicit PositionDependentCorrection(const std::string& name  = "PositionDependentCorrection",
                                       const std::string& fname = "MyNtuple.root");
  ~PositionDependentCorrection() override;

  int Init          (PHCompositeNode*) override;
  int process_event (PHCompositeNode*) override;
  int End           (PHCompositeNode*) override;

  /* ----------------------------------------------------------------
   *  Utility setters
   * -------------------------------------------------------------- */
  void setIsSimulation(bool sim)              { m_isSimulation = sim; }
  void Detector(const std::string& n)         { detector = n; }
  void set_timing_cut_width(int t)            { _range   = t; }
  void set_vertex_cut(float v)                { _vz      = v; }
  void apply_vertex_cut(bool on=true)         { m_vtxCut = on; }
  void setMassFitsDone(bool done) { m_massFitsDone = done; }
  void UseSignedVz(bool v = true) { m_useSignedVz = v; }
  void UseSurveyGeometry(bool v=true)         { m_useSurveyGeometry = v; }
    
  // First-pass switch: when true, only uncorrected TH3 views are booked/filled
  void setFirstPassBvaluesOnly(bool v=true)   { m_firstPassBvaluesOnly = v; }

  enum class EBinningMode { kRange , kDiscrete };
  void setBinningMode(EBinningMode m)         { m_binningMode = m; }
  EBinningMode getBinningMode() const         { return m_binningMode; }

  void setBEmcRec(BEmcRecCEMC* p)             { m_bemcRec = p; }
  void setTightVzCut(float cm)                { m_vzTightCut = std::fabs(cm); }

  /* ----------------------------------------------------------------
   *  Stand‑alone helpers that macros might call
   * -------------------------------------------------------------- */
  TF1*  fitHistogram(TH1* h);
  void  fitEtaSlices(const std::string& infile,
                     const std::string& fitOutFile,
                     const std::string& cdbFile);

  /* ----------------------------------------------------------------
   *  “process_…” helpers callable from macros (signatures unchanged)
   * -------------------------------------------------------------- */
  int process_g4hits   (PHCompositeNode*);
  int process_g4cells  (PHCompositeNode*);
  int process_towers   (PHCompositeNode*);
  int process_clusters (PHCompositeNode*);

  /* ****************************************************************
   *  PROTECTED section  – internal algorithms & data                *
   * ****************************************************************/
 protected:
    /* --------------------  vertex‑Z binning ------------------------
     *
     *  – m_vzTightCut      : 10 cm “physics” cut
     *  – vzEdge            : |z| table extended to 150 cm
     *  – getVzSlice()      : absolute‑|z| index          0 … N_VzBins‑1
     *  – getVzSliceSigned(): signed   index (+z first)   0 … 2*N_VzBins‑1
     * -------------------------------------------------------------- */
  bool skipVertexDep {true};
  float m_vzTightCut { 10.f };                 ///< |z| ≤ 10 cm

  static constexpr std::array<float,19> vzEdge = {
         0.f,  5.f, 10.f, 15.f, 20.f, 25.f, 30.f,
        40.f, 50.f, 60.f, 70.f, 80.f, 90.f,
       100.f,110.f,120.f,130.f,140.f,150.f
    };
  static constexpr int   N_VzBins       = vzEdge.size() - 1;
  static constexpr int   N_VzBinsSigned = 2 * N_VzBins;

  float m_vzSliceMax { vzEdge.back() };        ///< == 150 cm

  inline int getVzSlice(float absVz) const     /* 0 … N_VzBins‑1  (|z|) */
    {
      for (int i = 0; i < N_VzBins; ++i)
        if (absVz >= vzEdge[i] && absVz < vzEdge[i + 1]) return i;
      return -1;
    }

  inline int getVzSliceSigned(float vtxZ) const/* 0 … 2*N_VzBins‑1 (±z) */
    {
      const float absVz = std::fabs(vtxZ);
      for (int i = 0; i < N_VzBins; ++i)
        if (absVz >= vzEdge[i] && absVz < vzEdge[i + 1])
          return (vtxZ >= 0.f) ? i : i + N_VzBins;
      return -1;
    }

  /* user can switch signed histograms on/off at runtime */
  bool m_useSignedVz {false};
    
  /* --------------------  major functional helpers ---------------- */
  float retrieveVertexZ(PHCompositeNode*);

  void  fillTowerInfo(PHCompositeNode*, float thr,
                      float& tower_tot_e,
                      std::vector<float>& ht_eta,
                      std::vector<float>& ht_phi);

  void  fillTruthInfo(PHCompositeNode*, float& vtx_z,
                      std::vector<TLorentzVector>& truth_photons,
                      std::vector<TLorentzVector>& truth_meson_photons);

  RawTowerGeomContainer* checkTowerGeometry      (PHCompositeNode*);
  RawClusterContainer*   retrieveClusterContainer(PHCompositeNode*);

  int   countClusters(RawClusterContainer*, float vtx_z,
                      float ptCut,float chi2Cut,int& nCont);

  /* Block corrections, projections, centroids */
  float doPhiBlockCorr(float localPhi,float bphi);
  float doEtaBlockCorr(float localEta,float bEta,float dEtaTower=1.f);

  float convertBlockToGlobalPhi(int blkPhi,float localPhi);
  float convertBlockToGlobalEta(int blkEta,float localEta);

  float  phiAtShowerDepth(float energy,double rF,double zF,
                          float phiF,int ix,int iy)                    const;
  float  etaAtShowerDepth(float energy,double rF,double zF,
                          float phiF,int ix,int iy,float vtx_z)       const;
  double xAtShowerDepth  (float energy,double rF,double zF,
                          float phiF,int ix,int iy)                   const;

  float doAshShift        (float localPhi,float bVal);
  float doLogWeightCoord  (const std::vector<int>& towerPhi,
                           const std::vector<float>& towerE,float w0);
    
  void  bookCommonHistograms      (const std::function<std::string(int)>& makeLabel);
  void  bookSimulationHistograms  (const std::function<std::string(int)>& makeLabel);
  bool  loadBValues               (const std::string& bFilePath);
  bool  loadMassWindowTable       (const std::string& path);

  /* Driving loops */
  static constexpr int    N_Ebins = 8;
  static constexpr double eEdge[N_Ebins+1] = {2,4,6,8,10,12,15,20,30};
    
  void finalClusterLoop(PHCompositeNode*,RawClusterContainer*,float vtx_z,
                        const std::vector<TLorentzVector>& truth_photons,
                        const std::vector<TLorentzVector>& truth_meson_photons,
                        float tower_tot_e,float max_nClusCount,int nClusCount,
                        float maxAlpha,float ptMaxCut,float pt1Cut,float pt2Cut,
                        float pi0ptcut,float weight,bool fillGlobal);

  void fillAshLogDx(RawCluster*,const TLorentzVector& recoPhoton,
                    const TLorentzVector& truthPhoton,
                    const std::pair<float,float>& blkCord,int blkPhiBin,
                    const std::vector<int>& tower_phis,
                    const std::vector<float>& tower_energies);

  void fillDPhiAllVariants(RawCluster*,const TLorentzVector&,
                             const TLorentzVector&,const std::pair<float,float>&,
                             int blkPhiCoarse,float vtxZ,      /* ← named now */
                             TH1F* hRAW[N_Ebins],TH1F* hCP[N_Ebins],
                             bool fillGlobal);

  void fillDEtaAllVariants(RawCluster*,const TLorentzVector&,
                           const TLorentzVector&,const std::pair<float,float>&,
                           int blkEtaCoarse,int blkPhiCoarse,float vtx_z,
                           TH1F* hRAW[N_Ebins],TH1F* hCP[N_Ebins],
                           bool fillGlobal);

  int   getEnergySlice(float E) const;

  float getAvgEta(const std::vector<int>&,  const std::vector<float>&);
  float getAvgPhi(const std::vector<int>&,  const std::vector<float>&);

  std::pair<float,float>
        getBlockCord(const std::vector<int>& towerEtas,
                     const std::vector<int>& towerPhis,
                     const std::vector<float>& towerEs,
                     int& blkPhiOut,int& blkEtaOut);

  void  processClusterPairs(RawClusterContainer*,
                            RawClusterContainer::ConstIterator,
                            const CLHEP::Hep3Vector& vertex,
                            const TLorentzVector& photon1,float clusE,
                            int lt_eta,int lt_phi,
                            int blkEtaCoarse,int blkPhiCoarse,
                            const std::pair<float,float>& blkCoord,
                            float maxAlpha,float ptMaxCut,float pt2Cut,
                            float pi0ptcut,float weight,bool match1,
                            const TLorentzVector& ph1_trEtaPhi,
                            const std::vector<TLorentzVector>& truth_meson_photons,
                            bool isSimulation);

  /* --------------------  compile‑time constants ------------------- */
  static constexpr std::array<float,N_Ebins> expectedLo = []{
    std::array<float,N_Ebins> a{};
    for(int i=0;i<N_Ebins;++i) a[i]=static_cast<float>(eEdge[i]);
    return a;}();
  static constexpr std::array<float,N_Ebins> expectedHi = []{
    std::array<float,N_Ebins> a{};
    for(int i=0;i<N_Ebins;++i) a[i]=static_cast<float>(eEdge[i+1]);
    return a;}();

  static constexpr double kTolFactor = 1.5;
  static constexpr double kTolMinGeV = 0.20;
  static constexpr float  kDiscreteTol = 0.25f;

  inline std::string sliceTag(int i) const
  { return (m_binningMode==EBinningMode::kRange)
           ? Form("%.0f_%.0f",eEdge[i],eEdge[i+1])
           : Form("E%.0f",     eEdge[i]); }

  /* --------------------  data members ----------------------------- */
  /* steering / run‑environment */
  bool              m_isSimulation{false};
  bool              m_useSurveyGeometry{false};
  bool  m_massFitsDone{false};
  EBinningMode      m_binningMode{EBinningMode::kRange};
  bool              m_vtxCut{true};
  bool              alreadyDeclaredHistograms{false};
  bool              m_firstPassBvaluesOnly{false};   ///< fast path: only book/fill uncorrected TH3s


  std::string       detector;
  std::string       outfilename;

  /* geometry / helpers */
  RawTowerGeomContainer* m_geometry{nullptr};
  BEmcRecCEMC*           m_bemcRec{nullptr};
  TriggerAnalyzer*       trigAna{nullptr};

  /* survey tilt */
  float m_phi0Offset{0.f};
  bool  m_hasOffset{false};

  /* CDB */
  uint64_t   m_timeStamp{0};
  std::string m_cdbTag{"MDC2"};

    /* Ash / log scans */
    std::vector<double> m_bScan;
    std::vector<double> m_w0Scan;

    // Default (legacy) b-tables — unchanged behaviour
    float m_bValsPhi[N_Ebins]{};
    float m_bValsEta[N_Ebins]{};

    // Per-η-view b-tables (0: fullEta, 1: etaCore, 2: etaMid, 3: etaEdge)
    static constexpr int kNEtaViews = 4;
    std::array<std::array<float, N_Ebins>, kNEtaViews> m_bValsPhi_view{};
    std::array<std::array<float, N_Ebins>, kNEtaViews> m_bValsEta_view{};
    std::array<bool, kNEtaViews> m_bPhiReady_view{{false,false,false,false}};
    std::array<bool, kNEtaViews> m_bEtaReady_view{{false,false,false,false}};

  /* random */
  TRandom3* rnd{nullptr};

  /* global counters (unchanged) */
  std::atomic<std::uint64_t> m_nWinRAW{0},   m_nWinCP{0},   m_nWinBCorr{0};
  std::atomic<std::uint64_t> m_nWinRAW_Eta{0},m_nWinCP_Eta{0},m_nWinBCorr_Eta{0};
    std::atomic<std::uint64_t> m_phiWinCLUSraw{0}, m_phiWinCLUScp{0},
                               m_phiWinCLUSbcorr{0}, m_phiWinPDCraw{0},
                               m_phiWinPDCcorr{0},  m_phiWinCLUScpEA{0};
    std::atomic<std::uint64_t> m_etaWinCLUSraw{0}, m_etaWinCLUScp{0},
                               m_etaWinCLUSbcorr{0}, m_etaWinPDCraw{0},
                               m_etaWinPDCcorr{0},  m_etaWinCLUScpEA{0};
    
    // PDC side-peak tallies (± one-tower width)
    std::atomic<std::uint64_t> m_pdcRawPos1Tw{0},  m_pdcRawNeg1Tw{0};
    std::atomic<std::uint64_t> m_pdcCorrPos1Tw{0}, m_pdcCorrNeg1Tw{0};

    std::vector<float> m_pdcRawPos1Tw_dphi,  m_pdcRawPos1Tw_eta,   m_pdcRawPos1Tw_vz,   m_pdcRawPos1Tw_E;
    std::vector<float> m_pdcRawNeg1Tw_dphi,  m_pdcRawNeg1Tw_eta,   m_pdcRawNeg1Tw_vz,   m_pdcRawNeg1Tw_E;
    std::vector<float> m_pdcCorrPos1Tw_dphi, m_pdcCorrPos1Tw_eta,  m_pdcCorrPos1Tw_vz,  m_pdcCorrPos1Tw_E;
    std::vector<float> m_pdcCorrNeg1Tw_dphi, m_pdcCorrNeg1Tw_eta,  m_pdcCorrNeg1Tw_vz,  m_pdcCorrNeg1Tw_E;


    // CLUS (clusterizer) side-peak tallies (± one-tower width) for RAW and CP
    std::atomic<std::uint64_t> m_clusRawPos1Tw{0},  m_clusRawNeg1Tw{0};
    std::atomic<std::uint64_t> m_clusCPPos1Tw{0},   m_clusCPNeg1Tw{0};

    std::vector<float> m_clusRawPos1Tw_dphi, m_clusRawPos1Tw_eta,  m_clusRawPos1Tw_vz,  m_clusRawPos1Tw_E;
    std::vector<float> m_clusRawNeg1Tw_dphi, m_clusRawNeg1Tw_eta,  m_clusRawNeg1Tw_vz,  m_clusRawNeg1Tw_E;
    std::vector<float> m_clusCPPos1Tw_dphi,  m_clusCPPos1Tw_eta,   m_clusCPPos1Tw_vz,   m_clusCPPos1Tw_E;
    std::vector<float> m_clusCPNeg1Tw_dphi,  m_clusCPNeg1Tw_eta,   m_clusCPNeg1Tw_vz,   m_clusCPNeg1Tw_E;

    
  /* Δη out‑of‑window counters (|Δη| > 0.04) */
  std::atomic<std::uint64_t> m_etaOutCLUSraw   {0};
  std::atomic<std::uint64_t> m_etaOutCLUScp    {0};
  std::atomic<std::uint64_t> m_etaOutCLUSbcorr {0};
  std::atomic<std::uint64_t> m_etaOutPDCraw    {0};
  std::atomic<std::uint64_t> m_etaOutPDCcorr   {0};
    
    std::array<std::array<std::uint64_t, N_Ebins>, 6> m_phiWinByE{{}};
    std::array<std::array<std::uint64_t, N_Ebins>, 6> m_etaWinByE{{}};

  std::atomic<std::uint64_t> g_nanPhi{0}, g_nanEta{0}, g_nCorrPhiRight{0};

  /* ROOT output */
  Fun4AllHistoManager* hm{nullptr};
  TFile*               outfile{nullptr};

  /* --------------------  HISTOGRAM DECLARATIONS ------------------- */
    
  TH1F* h_phi_diff_raw_E        [N_Ebins]{};
  TH1F* h_phi_diff_cpRaw_E      [N_Ebins]{};
  TH1F* h_phi_diff_cpCorr_E     [N_Ebins]{};
  TH1F* h_phi_diff_cpCorrEA_E   [N_Ebins]{};
  TH1F* h_phi_diff_cpBcorr_E    [N_Ebins]{};
  TH1F* h_phi_diff_corrected_E  [N_Ebins]{};
    
  TH2F* h2_phi_diff_vsEta_RAW_E     [N_Ebins]{};
  TH2F* h2_phi_diff_vsEta_CP_E      [N_Ebins]{};
  TH2F* h2_phi_diff_vsEta_BCORR_E   [N_Ebins]{};
  TH2F* h2_phi_diff_vsEta_PDCraw_E  [N_Ebins]{};
  TH2F* h2_phi_diff_vsEta_PDCcorr_E [N_Ebins]{};
    
  /* |z|‑sliced and signed‑z φ histograms */
  TH1F* h_phi_diff_raw_E_vz        [N_Ebins][N_VzBins]{};
  TH1F* h_phi_diff_corrected_E_vz  [N_Ebins][N_VzBins]{};
  TH1F* h_phi_diff_cpRaw_E_vz      [N_Ebins][N_VzBins]{};
  TH1F* h_phi_diff_cpCorr_E_vz     [N_Ebins][N_VzBins]{};
  TH1F* h_phi_diff_cpBcorr_E_vz    [N_Ebins][N_VzBins]{};

  TH1F* h_phi_diff_raw_E_vzsgn        [N_Ebins][N_VzBinsSigned]{};
  TH1F* h_phi_diff_corrected_E_vzsgn  [N_Ebins][N_VzBinsSigned]{};
  TH1F* h_phi_diff_cpRaw_E_vzsgn      [N_Ebins][N_VzBinsSigned]{};
  TH1F* h_phi_diff_cpCorr_E_vzsgn     [N_Ebins][N_VzBinsSigned]{};
  TH1F* h_phi_diff_cpBcorr_E_vzsgn    [N_Ebins][N_VzBinsSigned]{};

  TH1F* h_eta_diff_raw_E        [N_Ebins]{};
  TH1F* h_eta_diff_cpRaw_E      [N_Ebins]{};
  TH1F* h_eta_diff_cpCorr_E     [N_Ebins]{};
  TH1F* h_eta_diff_cpCorrEA_E   [N_Ebins]{};
  TH1F* h_eta_diff_cpBcorr_E    [N_Ebins]{};
  TH1F* h_eta_diff_corrected_E  [N_Ebins]{};

  TH1F* h_eta_diff_raw_E_vz        [N_Ebins][N_VzBins]{};
  TH1F* h_eta_diff_corrected_E_vz  [N_Ebins][N_VzBins]{};
  TH1F* h_eta_diff_cpRaw_E_vz      [N_Ebins][N_VzBins]{};
  TH1F* h_eta_diff_cpCorr_E_vz     [N_Ebins][N_VzBins]{};
  TH1F* h_eta_diff_cpBcorr_E_vz    [N_Ebins][N_VzBins]{};

  TH1F* h_eta_diff_raw_E_vzsgn        [N_Ebins][N_VzBinsSigned]{};
  TH1F* h_eta_diff_corrected_E_vzsgn  [N_Ebins][N_VzBinsSigned]{};
  TH1F* h_eta_diff_cpRaw_E_vzsgn      [N_Ebins][N_VzBinsSigned]{};
  TH1F* h_eta_diff_cpCorr_E_vzsgn     [N_Ebins][N_VzBinsSigned]{};
  TH1F* h_eta_diff_cpBcorr_E_vzsgn    [N_Ebins][N_VzBinsSigned]{};
    
  TH1* h_tower_e{nullptr};
  TH1* h_clusE{nullptr};
  TH1* h_pt1{nullptr};
  TH1* h_pt2{nullptr};
  TH1* h_nclusters{nullptr};
  TH2* h_mass_eta_lt{nullptr};
  TH2* h_mass_eta_lt_rw{nullptr};
  TH2* h_pt_eta{nullptr};
  TH2* h_pt_eta_rw{nullptr};
  TH1* h_emcal_e_eta{nullptr};
  TH1* h_truth_eta{nullptr};
  TH1* h_truth_e{nullptr};
  TH1F* h_truth_vz{nullptr};
  TH1F* h_reco_vz{nullptr};
  TH2F* h2_truthReco_vz{nullptr};
  TH1* h_pt_rw[96]{};
  TFile* frw{nullptr};
  TH1* h_matched_res{nullptr};
  TH1* h_res_e{nullptr};
  TH3* h_res_e_phi{nullptr};
  TH3* h_res_e_eta{nullptr};
  TH3* h_res_e_eta_pdc{nullptr};
  TH1* h_res{nullptr};
  TH3* h_m_pt_eta{nullptr};
  TH3* h_m_ptTr_eta{nullptr};
  TH3* h_m_ptTr_eta_trKin{nullptr};
  TH3* h_delPhi_e_eta{nullptr};
  TH3* h_delEta_e_eta{nullptr};
  TH3* h_delR_e_eta{nullptr};
  TH3* h_delPhi_e_phi{nullptr};
  TProfile* pr_eta_shower{nullptr};
  TProfile* pr_phi_shower{nullptr};
  TH2* h_vert_xy{nullptr};
  TH1* h_truthE{nullptr};

  TH2F* h2_chi2_tot_etaPhi{nullptr};
  TH2F* h2_chi2_rej_etaPhi{nullptr};
  TProfile2D* p_chi2_pass_etaPhi{nullptr};

  static const int NBinsBlock = 14;
  TH2* h_mass_block_pt[NBinsBlock][NBinsBlock]{};
  TH2* h_res_block_E [NBinsBlock][NBinsBlock]{};

    // Legacy default uncorrected TH3 and corrected TH3 (unchanged)
    TH3* h3_cluster_block_cord_E{nullptr};
    TH3* h3_cluster_block_cord_E_corrected{nullptr};

    // New: uncorrected η-sliced variants (filled in parallel with the default)
    TH3* h3_cluster_block_cord_E_full{nullptr};        // full |eta| <= 1.10
    TH3* h3_cluster_block_cord_E_etaCore{nullptr};     // |eta| <= 0.20
    TH3* h3_cluster_block_cord_E_etaMid{nullptr};      // 0.20 < |eta| <= 0.70
    TH3* h3_cluster_block_cord_E_etaEdge{nullptr};     // 0.70 < |eta| <= 1.10

    // New: corrected counterparts per η-slice (second pass)
    TH3* h3_cluster_block_cord_E_full_corr{nullptr};
    TH3* h3_cluster_block_cord_E_etaCore_corr{nullptr};
    TH3* h3_cluster_block_cord_E_etaMid_corr{nullptr};
    TH3* h3_cluster_block_cord_E_etaEdge_corr{nullptr};
    
//    TH3* h3_cluster_block_cord_E_full{nullptr};        // full |eta| <= 1.10
//    TH3* h3_cluster_block_cord_E_etaCore{nullptr};     // |eta| <= 0.20
//    TH3* h3_cluster_block_cord_E_etaMid{nullptr};      // 0.20 < |eta| <= 0.70
//    TH3* h3_cluster_block_cord_E_etaEdge{nullptr};     // 0.70 < |eta| <= 1.10
//    TH3* h3_cluster_block_cord_E_corrected{nullptr};
    
    
  TH1* h_block_phi{nullptr};
  TH1* h_block_eta{nullptr};
  TH1* h_clus_E_size{nullptr};
  TH1* h_block_bin{nullptr};

  TH2* h_emcal_mbd_correlation{nullptr};
  TH2* h_ohcal_mbd_correlation{nullptr};
  TH2* h_ihcal_mbd_correlation{nullptr};
  TH2* h_emcal_hcal_correlation{nullptr};
  TH2* h_emcal_zdc_correlation{nullptr};
  TH2* h_cemc_etaphi{nullptr};
  TH2* h_hcalin_etaphi{nullptr};
  TH2* h_hcalout_etaphi{nullptr};
  TH2* h_cemc_etaphi_wQA{nullptr};
  TH2* h_hcalin_etaphi_wQA{nullptr};
  TH2* h_hcalout_etaphi_wQA{nullptr};
    
  TH1F*     h_InvMass{nullptr};
  TH1F*     h_InvMass_w{nullptr};
  TH1F*     h_InvMassMix{nullptr};

  TH2F*     h_mE_raw{nullptr};
  TH2F*     h_mE_corr{nullptr};
  TH3F*     h_m_blk_raw{nullptr};
  TH3F*     h_m_blk_corr{nullptr};

  TH2F*     h_clusE_nTow{nullptr};
  TProfile* pr_phi_vs_blockcoord{nullptr};

  TH1F*     h_localPhi_corrected[N_Ebins]{};
  TH1F*     h_localEta_corrected[N_Ebins]{};

  TH1* h_totalzdc_e{nullptr};
  TH1* h_delR_recTrth{nullptr};
  TProfile2D* h_cemc_etaphi_time{nullptr};
  TProfile2D* h_hcalin_etaphi_time{nullptr};
  TProfile2D* h_hcalout_etaphi_time{nullptr};
  TProfile2D* h_cemc_etaphi_badChi2{nullptr};
  TProfile2D* h_hcalin_etaphi_badChi2{nullptr};
  TProfile2D* h_hcalout_etaphi_badChi2{nullptr};

  TH1* hzdctime{nullptr};        TH1* hmbdtime{nullptr};
  TH1* hemcaltime{nullptr};      TH1* hihcaltime{nullptr};
  TH1* hohcaltime{nullptr};

  TH1* hzdctime_cut{nullptr};    TH1* hmbdtime_cut{nullptr};
  TH1* hemcaltime_cut{nullptr};  TH1* hihcaltime_cut{nullptr};
  TH1* hohcaltime_cut{nullptr};

  TH1* hvtx_z_raw{nullptr};
  TH1* hvtx_z_cut{nullptr};
  TH1* hzdcSouthraw{nullptr};
  TH1* hzdcNorthraw{nullptr};
  TH1* hzdcSouthcalib{nullptr};
  TH1* hzdcNorthcalib{nullptr};
  TH2* h_etaphi_clus{nullptr};

  TNtuple* g4hitntuple{nullptr};
  TNtuple* g4cellntuple{nullptr};
  TTree*   towerntuple{nullptr};
  TNtuple* clusterntuple{nullptr};

  /* vectors cached for ntuples */
  std::vector<float> m_energy;
  std::vector<int>   m_etabin;
  std::vector<int>   m_phibin;
  std::vector<int>   m_time;
  std::vector<float> m_hcalin_energy;
  std::vector<int>   m_hcalin_etabin;
  std::vector<int>   m_hcalin_phibin;
  std::vector<int>   m_hcalin_time;
  std::vector<float> m_hcalout_energy;
  std::vector<int>   m_hcalout_etabin;
  std::vector<int>   m_hcalout_phibin;
  std::vector<int>   m_hcalout_time;
  std::vector<float> m_zdc_energy;
  std::vector<int>   m_zdc_index;
  std::vector<int>   m_zdc_side;
  std::vector<float> m_bbc_energy;
  std::vector<int>   m_bbc_type;
  std::vector<int>   m_bbc_side;

  /* misc run bookkeeping */
  int   _eventcounter{0};
  int   _range{1};
  float _vz{0.f};
  bool  dynMaskClus{false};
  bool  debug{false};

  /* trigger labels */
  std::map<std::string,std::string> triggerNameMap {
    {"MBD N&S >= 1","MBD_NandS_geq_1"} };
  std::map<int,std::string>* activeTriggerNameMap{nullptr};

  /* π0 mass window */
  bool  isFitDoneForPhi{true};
  bool  isFitDoneForEta{true};
  float target_pi0_mass{0.145};

  struct MassWindow { float mu{0.f}; float sigma{0.f}; };
  MassWindow m_winRaw [N_Ebins]{};
  MassWindow m_winCorr[N_Ebins]{};

  /* fine‑block tallies */
  std::array<std::atomic<std::uint64_t>,4> m_blkLocCount {{0,0,0,0}};
};




/* =================================================================== *
 *  4.  Inline helpers (unchanged algorithms; whitespace only)          *
 * =================================================================== */
namespace PDC_detail
{
  inline uint64_t vbLevel()
  { return PositionDependentCorrection::getVerbosityLevel(); }

    inline float cg2GlobalPhi(BEmcRecCEMC* rec,
                              float eReco,float xT,float yT)
    {
      float gx = 0.f, gy = 0.f, gz = 0.f;
      rec->Tower2Global(eReco, xT, yT, gx, gy, gz);
      return std::atan2(gy, gx);
    }

    inline float cg2ShowerEta(BEmcRecCEMC* rec,
                              float eReco,float xT,float yT,float vtxZ)
    {
      float gx = 0.f, gy = 0.f, gz = 0.f;
      rec->Tower2Global(eReco, xT, yT, gx, gy, gz);
      return std::asinh((gz - vtxZ) / std::hypot(gx, gy));
    }


    inline float front2ShowerPhi(BEmcRecCEMC* rec,float eReco,
                                 double rF,double zF,float phiF,
                                 int ix,int iy)
    {
      // Transport the chosen front-face point to shower depth
      float xSD,ySD,zSD;
      rec->CorrectShowerDepth(ix,iy,eReco,
                              static_cast<float>(rF*std::cos(phiF)),
                              static_cast<float>(rF*std::sin(phiF)),
                              static_cast<float>(zF),
                              xSD,ySD,zSD);
      float phiSD = std::atan2(ySD,xSD);

      // Owner-branch guard: fold around the owner tower center at shower depth
      TowerGeom g{};
      if (rec->GetTowerGeometry(ix, iy, g))
      {
        float xcC=0.f, ycC=0.f, zcC=0.f;
        rec->CorrectShowerDepth(ix, iy, eReco,
                                g.Xcenter, g.Ycenter, g.Zcenter,
                                xcC, ycC, zcC);
        const float phiC = std::atan2(ycC, xcC);
        const float tw   = 2.f * static_cast<float>(M_PI)
                           / static_cast<float>(rec->GetNx());

        float d = static_cast<float>(
                    std::remainder(static_cast<double>(phiSD - phiC),
                                   static_cast<double>(tw))
                  );
        if (d <= -0.5f*tw) d += tw;
        if (d >   0.5f*tw) d -= tw;

        phiSD = phiC + d;
      }
      return phiSD;
    }


  inline float front2ShowerEta(BEmcRecCEMC* rec,
                                 float        eReco,     // cluster E [GeV]
                                 double       /*unused*/,// legacy rFront
                                 int          ix,        // tower φ index
                                 int          iy,        // tower η index
                                 float        etaF,      // desired η at front
                                 float        /*phiF*/,
                                 float        vtxZ)      // vertex-z [cm]
    {
        // --- geometry lookup -------------------------------------------------
        TowerGeom g{};
        if (!rec || !rec->GetTowerGeometry(ix, iy, g))
            return std::numeric_limits<float>::quiet_NaN();

        const double rC = std::hypot(g.Xcenter, g.Ycenter);
        const double xF = g.Xcenter;
        const double yF = g.Ycenter;
        const double zF = rC * std::sinh(etaF);           // z on front face that yields ηF

        // --- one‑shot transport to shower depth ------------------------------
        float xSD{}, ySD{}, zSD{};
        rec->CorrectShowerDepth(ix, iy, eReco,
                                static_cast<float>(xF),
                                static_cast<float>(yF),
                                static_cast<float>(zF),
                                xSD, ySD, zSD);

        // --- finished: global η after depth spline + saw‑tooth ---------------
        return std::asinh((zSD - vtxZ) /
                          std::hypot(xSD, ySD));
    }


    /* ---------- cluster centre of gravity (Momenta-equivalent) ------ */
    inline bool clusterCentreOfGravity(const RawCluster* cluster,
                                       RawTowerGeomContainer* geo,
                                       const BEmcRecCEMC* rec,
                                       float& x, float& y)
    {
        if (!cluster || !geo || !rec)
        {
            x = y = 0.0f;
            return false;
        }

        // 1) Build the EmcModule hit list exactly like the clusterizer expects
        std::vector<EmcModule> hitlist;
        hitlist.reserve(std::distance(cluster->get_towers().first,
                                      cluster->get_towers().second));

        const int Nx = rec->GetNx();    // width of the tower grid used for linear indexing

        RawCluster::TowerConstRange towers = cluster->get_towers();
        for (auto it = towers.first; it != towers.second; ++it)
        {
            const auto tk   = it->first;
            const double e  = it->second;
            const auto tg   = geo->get_tower_geometry(tk);
            if (!tg) continue;

            const int ix = tg->get_binphi();
            const int iy = tg->get_bineta();

            EmcModule m;
            m.ich = iy * Nx + ix;               // linear index: ich = iy*fNx + ix
            m.amp = static_cast<float>(e);      // tower energy
            m.tof = 0.0f;                       // unused here
            hitlist.push_back(m);
        }

        if (hitlist.empty())
        {
            x = y = 0.0f;
            return false;
        }

        // 2) Call the clusterizer's CoG routine with its threshold and φ-wrapping logic
        float E = 0.0f, px = 0.0f, py = 0.0f, pxx = 0.0f, pyy = 0.0f, pyx = 0.0f;

        // If GetTowerThreshold() is non-const in your build, the const_cast below keeps this helper drop-in:
        const float thr = const_cast<BEmcRecCEMC*>(rec)->GetTowerThreshold();

        rec->Momenta(&hitlist, E, px, py, pxx, pyy, pyx, thr);
        if (E <= 0.0f)
        {
            x = y = 0.0f;
            return false;
        }

        // Momenta already adds back the max-tower indices and wraps φ to [-0.5, Nx-0.5)
        x = px;   // φ-like tower coordinate (xC) in tower units
        y = py;   // η-like tower coordinate (yC) in tower units
        return true;
    }


  /* ---------- light‑weight records for Δη/Δφ bookkeeping ---------- */
  struct PhiRec { const char* tag; float loc; float phi; float d{0.f}; };
  struct EtaRec { const char* tag; float loc; float eta; float d{0.f}; };
} // namespace PDC_detail

#endif  // PositionDependentCorrection_H__
