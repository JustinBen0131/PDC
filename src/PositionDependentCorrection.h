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
#include <TLorentzVector.h>
#include <cstdint>
#include <atomic>

#include <CLHEP/Vector/ThreeVector.h>   // Hep3Vector in inline helpers
#include <fun4all/SubsysReco.h>
#include <calotrigger/TriggerAnalyzer.h>
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawCluster.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawClusterv1.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawClusterContainer.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawTowerDefs.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawTowerGeom.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawTowerGeomContainer.h"
#include <TString.h>
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloReco/BEmcRec.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloReco/BEmcRecCEMC.h"


/* =================================================================== *
 *  2.  Forward declarations (ROOT & sPHENIX)                           *
 * =================================================================== */
class Fun4AllHistoManager;
class PHCompositeNode;
class RawTowerGeomContainer;
class RawClusterContainer;
class BEmcRecCEMC;
class RawCluster;
class RawClusterv1;

/* ROOT */
class TFile;  class TNtuple; class TTree;
class TH1;    class TH1F;    class TH2;    class TH2F; class TH3;    class TH3F;
class TF1;    class TProfile;class TProfile2D;
class TRandom3;

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

    // NEW: single-photon steering (truth-only projections, skip pair loop)
    void setIsSinglePhoton(bool on = true)      { m_isSinglePhoton = on; }

    void set_vertex_cut(float v)                { _vz      = v; }
    void apply_vertex_cut(bool on=true)         { m_vtxCut = on; }
    void setMassFitsDone(bool done)             { m_massFitsDone = done; }
    void UseSignedVz(bool v = true)             { m_useSignedVz = v; }
    void UseSurveyGeometry(bool v=true)         { m_useSurveyGeometry = v; }
      
    // First-pass switch: when true, only uncorrected TH3 views are booked/filled
    void setFirstPassBvaluesOnly(bool v=true)   { m_firstPassBvaluesOnly = v; }

    enum class EBinningMode { kRange , kDiscrete };
    void setBinningMode(EBinningMode m)         { m_binningMode = m; }
    EBinningMode getBinningMode() const         { return m_binningMode; }

    void setBEmcRec(BEmcRecCEMC* p)             { m_bemcRec = p; }
    void setTightVzCut(float cm)                { m_vzTightCut = std::fabs(cm); }

    // --- Quality gates steering ---
    enum class AnchorQualityMode { kNone, kChi2, kProb, kBoth };

    // Choose which anchor-cluster gate to use
    void setAnchorQualityMode(AnchorQualityMode m) { m_anchorQMode = m; }
    void setAnchorChi2Cut(float v)                 { m_anchorChi2Cut = v; }   // default 10
    void setAnchorProbMin(float v)                 { m_anchorProbMin  = v; }   // default 0.05

    // Enable/disable partner (second cluster) probability cut and set its threshold
    void enableSecondClusterProbCut(bool on = true){ m_enableC2ProbCut = on; } // default ON
    void setSecondClusterProbMin(float v)          { m_pairProbMin     = v; }  // default 0.05

    // Enable/disable partner (second cluster) χ² cut and set its threshold
    void enableSecondClusterChi2Cut(bool on = true) { m_enableC2Chi2Cut = on; } // default OFF
    void setSecondClusterChi2Max(float v)           { m_pairChi2Max     = v; }  // default 10.0

    // Optional convenience: expose anchor χ² (useful to tie partner thr to anchor)
    float getAnchorChi2Cut() const                  { return m_anchorChi2Cut; }

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

    int  m_print_first_N_clusters = 12;   // how many clusters to print nicely per job
    bool m_printed_prop_header    = false;
    TH2F* h_dx_prop_vsE = nullptr;        // (optional) quick QA
    TH2F* h_dy_prop_vsE = nullptr;
    
    AnchorQualityMode m_anchorQMode{AnchorQualityMode::kChi2}; // default: χ² gate
    float m_anchorChi2Cut{10.0f};   // χ² threshold when χ² gate is enabled
    float m_anchorProbMin{0.05f};   // probability threshold for anchor (used if kProb/kBoth)

    bool  m_enableC2ProbCut{false}; // default: partner probability cut OFF
    float m_pairProbMin{0.05f};     // threshold if partner-prob cut is enabled

    bool  m_enableC2Chi2Cut{false}; // default: partner χ² cut OFF
    float m_pairChi2Max{10.0f};     // threshold if partner-χ² cut is enabled
    
    
  static constexpr std::array<float,10> vzEdge = {
         0.f,  5.f, 10.f, 15.f, 20.f, 25.f, 30.f,
        40.f, 50.f, 60.f
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

  void  fillTruthInfo(PHCompositeNode*, float& truth_vz,
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

  float doAshShift(float xTower, float b);
  float doLogWeightCoord(const std::vector<int>& towerPhi,
                           const std::vector<float>& towerE,
                           float w0);
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

    void fillAshLogDx(RawCluster*, const TLorentzVector& recoPhoton,
                      const TLorentzVector& truthPhoton,
                      float vtxZ,
                      const std::pair<float,float>& blkCord, int blkPhiBin,
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
    
    // CP vs CP(EA) tie tolerances (units: radians for φ, unitless for η)
    static constexpr float  kTieTolPhi = 1e-8f;
    static constexpr float  kTieTolEta = 1e-8f;

  inline std::string sliceTag(int i) const
  { return (m_binningMode==EBinningMode::kRange)
           ? Form("%.0f_%.0f",eEdge[i],eEdge[i+1])
           : Form("E%.0f",     eEdge[i]); }

  /* --------------------  data members ----------------------------- */
  /* steering / run‑environment */
    bool              m_isSimulation{false};
    bool              m_useSurveyGeometry{false};
    bool              m_massFitsDone{false};
    EBinningMode      m_binningMode{EBinningMode::kRange};

    // NEW: when true, use truth vertex everywhere and skip pair loop
    bool              m_isSinglePhoton{false};

    bool              m_vtxCut{true};
    bool              alreadyDeclaredHistograms{false};
    bool              m_firstPassBvaluesOnly{false};   ///< fast path: only book/fill uncorrected TH3s
    
    struct Pi0CutCounters {
        // Event accounting
      std::uint64_t ev_total         = 0;  // all events entering process_towers
      std::uint64_t ev_processed     = 0;  // events that pass |vz| acceptance
      std::uint64_t ev_vz_skipped    = 0;  // events rejected by |vz| > m_vzSliceMax
      std::uint64_t ev_trig_skipped  = 0;  // events skipped because required data trigger not fired
      std::uint64_t ev_nClusTooLarge = 0;  // events rejected by nClusCount > max_nClusCount

      std::uint64_t c1_prob         = 0; // anchor rejected by probability
      std::uint64_t pairs_prob2     = 0; // partner rejected by probability
      // Cluster‑1 (anchor) screening
      std::uint64_t c1_seen         = 0;   // clusters iterated in OUTER loop
      std::uint64_t c1_E_nan        = 0;   // E_vec_1.mag() NaN
      std::uint64_t c1_E_zero       = 0;   // E_vec_1.mag() ~ 0
      std::uint64_t c1_E_below      = 0;   // clusE < 0.1
      std::uint64_t c1_chi2         = 0;   // χ² (anchor) cut
      std::uint64_t c1_ltEta_oob    = 0;   // lead tower eta > 95
      std::uint64_t c1_pt1_outside  = 0;   // pt1 outside [pt1ClusCut, ptMaxCut]
      std::uint64_t c1_pass         = 0;   // clusters that reached pair loop

        // Pair‑level screening (after a cluster1 passed)
      std::uint64_t pairs_seen        = 0; // all pair trials (including self/null)
      std::uint64_t pairs_self_or_null= 0; // self pairing or null c2
      std::uint64_t pairs_pt2         = 0; // pt2 outside [pt2ClusCut, ptMaxCut]
      std::uint64_t pairs_chi2        = 0; // χ²(c2) pathological guard
      std::uint64_t pairs_chi2_gate   = 0; // χ²(c2) failed configured partner-χ² gate
      std::uint64_t pairs_alpha       = 0; // asymmetry > maxAlpha
      std::uint64_t pairs_mass_nan    = 0; // π0 mass NaN
      std::uint64_t pairs_pi0pt       = 0; // π0 pt cut (simulation mode only)
      std::uint64_t pairs_outOfSlice  = 0; // eLead outside energy slice table (info)
      std::uint64_t pairs_final       = 0; // pairs that reached “legacy” mass fills

      // Side fills (diagnostics)
      std::uint64_t pass1_fills       = 0; // h_mE_raw/h_mE_corr fills
      std::uint64_t block_raw_fills   = 0; // block raw map fills (μ±3σ)
      std::uint64_t block_corr_fills  = 0; // block corr map fills (μ±3σ)
      std::uint64_t truth_pairs_ok    = 0; // both reco photons matched to same truth π0
      std::uint64_t variant_mass_fills= 0; // per‑variant π0 mass histogram fills
    };
    static Pi0CutCounters s_cuts;

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
    
  /* global counters (unchanged + ties added) */
  std::atomic<std::uint64_t> m_nWinRAW{0},   m_nWinCP{0},   m_nWinBCorr{0};
  std::atomic<std::uint64_t> m_nWinRAW_Eta{0},m_nWinCP_Eta{0},m_nWinBCorr_Eta{0};
  std::atomic<std::uint64_t> m_phiWinCLUSraw{0}, m_phiWinCLUScp{0},
                                 m_phiWinCLUSbcorr{0}, m_phiWinPDCraw{0},
                                 m_phiWinPDCcorr{0},  m_phiWinCLUScpEA{0};
  std::atomic<std::uint64_t> m_etaWinCLUSraw{0}, m_etaWinCLUScp{0},
                                 m_etaWinCLUSbcorr{0}, m_etaWinPDCraw{0},
                                 m_etaWinPDCcorr{0},  m_etaWinCLUScpEA{0};


  std::atomic<std::uint64_t> m_phiTieCLUScpEA{0};
  std::atomic<std::uint64_t> m_etaTieCLUScpEA{0};


  std::atomic<std::uint64_t> m_phiNoChange_CP_vs_RAW{0};
  std::atomic<std::uint64_t> m_phiNoChange_EA_vs_RAW{0};
  std::atomic<std::uint64_t> m_etaNoChange_CP_vs_RAW{0};
  std::atomic<std::uint64_t> m_etaNoChange_EA_vs_RAW{0};

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


  std::array<std::uint64_t, N_Ebins> m_phiTieByE_CP_EA{};
  std::array<std::uint64_t, N_Ebins> m_etaTieByE_CP_EA{};


  std::array<std::uint64_t, N_Ebins> m_phiNoChangeByE_CP_vs_RAW{};
  std::array<std::uint64_t, N_Ebins> m_phiNoChangeByE_EA_vs_RAW{};
  std::array<std::uint64_t, N_Ebins> m_etaNoChangeByE_CP_vs_RAW{};
  std::array<std::uint64_t, N_Ebins> m_etaNoChangeByE_EA_vs_RAW{};

  std::atomic<std::uint64_t> g_nanPhi{0}, g_nanEta{0}, g_nCorrPhiRight{0};

  /* ROOT output */
  Fun4AllHistoManager* hm{nullptr};
  TFile*               outfile{nullptr};

  /* --------------------  HISTOGRAM DECLARATIONS ------------------- */
    
  // Nine reconstruction variants for π0 studies
  enum class VarPi0 {
      CLUS_RAW=0, CLUS_CP=1,
      EA_FIT_zDEP=2, EA_FIT_ETADEP=3, EA_FIT_EONLY=4, EA_FIT_EONLY_INCIDENT=5,
      EA_FIT_ZVTXETADEP=6, PDC_RAW=7, PDC_CORR=8,
      NVAR=9
  };

  // Truth π0 photons with their mother info (filled in fillTruthInfo)
  struct TruthPhoton { TLorentzVector p4; int mother_id{0}; int mother_pid{0}; };
  std::vector<TruthPhoton> m_truth_pi0_photons;

  // Per-variant, per-slice π0 mass histograms (8 variants × N_Ebins)
  TH1F* h_m_pi0_var[static_cast<int>(VarPi0::NVAR)][N_Ebins]{};

  // NEW: Per-η-view, per-variant, per-slice π0 mass histograms
  // view index: 0=fullEta (|η|≤1.10), 1=etaCore (|η|≤0.20),
  //             2=etaMid (0.20<|η|≤0.70), 3=etaEdge (0.70<|η|≤1.10)
  TH1F* h_m_pi0_var_eta[kNEtaViews][static_cast<int>(VarPi0::NVAR)][N_Ebins]{};

  TH1F* h_phi_diff_raw_E        [N_Ebins]{};
  TH1F* h_phi_diff_cpRaw_E      [N_Ebins]{};
  TH1F* h_phi_diff_cpCorr_E     [N_Ebins]{};
  TH1F* h_phi_diff_cpCorrEA_E   [N_Ebins]{};
  TH1F* h_phi_diff_cpBcorr_E    [N_Ebins]{};
  TH1F* h_phi_diff_corrected_E  [N_Ebins]{};
    
  TH1F* h_phi_diff_cpCorrEA_fitEnergyOnly_AndIncidentAngle_E[N_Ebins] = {nullptr};
  TH1F* h_phi_diff_cpCorrEA_fitZVTXEtaDep_E[N_Ebins]                  = {nullptr};
  TH1F* h_eta_diff_cpCorrEA_fitEnergyOnly_AndIncidentAngle_E[N_Ebins] = {nullptr};
  TH1F* h_eta_diff_cpCorrEA_fitZVTXEtaDep_E[N_Ebins]                  = {nullptr};
    
    
  // --- NEW: per-slice (E) 2D maps for CP(EA) ---
  // x = b_phi  , y = Δφ(folded); one per energy bin
  TH2F* h2_phi_diffEA_vs_bphi_E [N_Ebins]{};
  // x = b_eta  , y = Δη        ; one per energy bin
  TH2F* h2_eta_diffEA_vs_beta_E [N_Ebins]{};

  // --- EA residuals (φ) – four uniquely-named variants
  TH1F* h_phi_diff_cpCorrEA_fitZVTXDep_E                          [N_Ebins]{};
  TH1F* h_phi_diff_cpCorrEA_fitEtaDep_E                     [N_Ebins]{};
  TH1F* h_phi_diff_cpCorrEA_fitEnergyOnly_E                 [N_Ebins]{};
  TH1F* h_phi_diff_cpCorrEA_fitPhiEnergy_etaEtaDep_E        [N_Ebins]{};

  // --- EA residuals (η) – four uniquely-named variants
  TH1F* h_eta_diff_cpCorrEA_fitZVTXDep_E                          [N_Ebins]{};
  TH1F* h_eta_diff_cpCorrEA_fitEtaDep_E                     [N_Ebins]{};
  TH1F* h_eta_diff_cpCorrEA_fitEnergyOnly_E                 [N_Ebins]{};
  TH1F* h_eta_diff_cpCorrEA_fitPhiEnergy_etaEtaDep_E        [N_Ebins]{};

  // ---------- NEW: |z| ≤ 60 cm gated duplicates (suffix: _0_60vz) ----------
  // Δφ (cluster family + scratch)
  TH1F* h_phi_diff_cpRaw_E_0_60vz                           [N_Ebins]{};
  TH1F* h_phi_diff_cpCorr_E_0_60vz                          [N_Ebins]{};
  TH1F* h_phi_diff_cpCorrEA_fitZVTXDep_E_0_60vz                   [N_Ebins]{};
  TH1F* h_phi_diff_cpCorrEA_fitEtaDep_E_0_60vz              [N_Ebins]{};
  TH1F* h_phi_diff_cpCorrEA_fitEnergyOnly_E_0_60vz          [N_Ebins]{};
  TH1F* h_phi_diff_cpCorrEA_fitPhiEnergy_etaEtaDep_E_0_60vz [N_Ebins]{};
  TH1F* h_phi_diff_cpBcorr_E_0_60vz                         [N_Ebins]{};
  TH1F* h_phi_diff_raw_E_0_60vz                             [N_Ebins]{};
  TH1F* h_phi_diff_corrected_E_0_60vz                       [N_Ebins]{};

  // Δη (cluster family + scratch)
  TH1F* h_eta_diff_cpRaw_E_0_60vz                           [N_Ebins]{};
  TH1F* h_eta_diff_cpCorr_E_0_60vz                          [N_Ebins]{};
  TH1F* h_eta_diff_cpCorrEA_fitZVTXDep_E_0_60vz                   [N_Ebins]{};
  TH1F* h_eta_diff_cpCorrEA_fitEtaDep_E_0_60vz              [N_Ebins]{};
  TH1F* h_eta_diff_cpCorrEA_fitEnergyOnly_E_0_60vz          [N_Ebins]{};
  TH1F* h_eta_diff_cpCorrEA_fitPhiEnergy_etaEtaDep_E_0_60vz [N_Ebins]{};
  TH1F* h_eta_diff_cpBcorr_E_0_60vz                         [N_Ebins]{};
  TH1F* h_eta_diff_raw_E_0_60vz                             [N_Ebins]{};
  TH1F* h_eta_diff_corrected_E_0_60vz                       [N_Ebins]{};
    

  // --- Ash scan profiles: <(Δφ)^2> / <(Δη)^2> vs b (tower-space driven)
  TProfile* p_phi_rms2_vs_b_E [N_Ebins]{};  // x: b, y: <(Δφ)^2>
  TProfile* p_eta_rms2_vs_b_E [N_Ebins]{};  // x: b, y: <(Δη)^2>

  // --- Log scan profiles: <(Δφ)^2> / <(Δη)^2> vs w0 (tower-space driven)
  TProfile* p_phi_rms2_vs_w0_E[N_Ebins]{};
  TProfile* p_eta_rms2_vs_w0_E[N_Ebins]{};

  // --- Agreement diagnostic in tower units: |CP(b) - CP(EA)|
  TProfile* p_abs_dloc_phi_CPea_vs_b_E [N_Ebins]{};
  TProfile* p_abs_dloc_eta_CPea_vs_b_E [N_Ebins]{};
    
  // ----- 6-way corrected winners (Δφ): RAW, CP, EA_geom, EA_fitEta, EA_fitE, EA_mix
  std::uint64_t m_phi6WayWin_RAW{0}, m_phi6WayWin_CP{0},
                  m_phi6WayWin_EA_geom{0}, m_phi6WayWin_EA_fitEta{0},
                  m_phi6WayWin_EA_fitE{0},  m_phi6WayWin_EA_mix{0};
  std::uint64_t m_phi6WayWinByE_RAW[N_Ebins]{}, m_phi6WayWinByE_CP[N_Ebins]{},
                  m_phi6WayWinByE_EA_geom[N_Ebins]{}, m_phi6WayWinByE_EA_fitEta[N_Ebins]{},
                  m_phi6WayWinByE_EA_fitE[N_Ebins]{},  m_phi6WayWinByE_EA_mix[N_Ebins]{};

  // ----- 6-way corrected winners (Δη): RAW, CP, EA_geom, EA_fitEta, EA_fitE, EA_mix
  std::uint64_t m_eta6WayWin_RAW{0}, m_eta6WayWin_CP{0},
                  m_eta6WayWin_EA_geom{0}, m_eta6WayWin_EA_fitEta{0},
                  m_eta6WayWin_EA_fitE{0},  m_eta6WayWin_EA_mix{0};
  std::uint64_t m_eta6WayWinByE_RAW[N_Ebins]{}, m_eta6WayWinByE_CP[N_Ebins]{},
                  m_eta6WayWinByE_EA_geom[N_Ebins]{}, m_eta6WayWinByE_EA_fitEta[N_Ebins]{},
                  m_eta6WayWinByE_EA_fitE[N_Ebins]{},  m_eta6WayWinByE_EA_mix[N_Ebins]{};

    
    /* --------------------  NEW: |z_vtx| slice TH3s (uncorrected)  -------------------- */
    /* COARSE bins for 0–150 cm: 0–10, 10–20, 20–30, 30–45, 45–60*/
    static constexpr int   N_VzH3Bins = 5;
    static constexpr float vzH3Edges[N_VzH3Bins + 1] = { 0.f, 10.f, 20.f, 30.f, 45.f, 60.f};

    /* FINE bins for 0–10 cm: 0–2, 2–4, 4–6, 6–8, 8–10 cm (unchanged) */
    static constexpr int   N_VzH3FineBins = 5;
    static constexpr float vzH3FineEdges[N_VzH3FineBins + 1] = { 0.f, 2.f, 4.f, 6.f, 8.f, 10.f };

    /* Helper: map |z| to 0..N_VzH3Bins-1 (returns -1 if outside 0–150) */
    inline int getVzH3Slice(float absVz) const
    {
      for (int i = 0; i < N_VzH3Bins; ++i)
        if (absVz >= vzH3Edges[i] && absVz < vzH3Edges[i + 1]) return i;
      return -1;
    }

    /* Helper (fine): map |z| to 0..N_VzH3FineBins-1 (returns -1 if outside 0–10) */
    inline int getVzH3FineSlice(float absVz) const
    {
      for (int i = 0; i < N_VzH3FineBins; ++i)
        if (absVz >= vzH3FineEdges[i] && absVz < vzH3FineEdges[i + 1]) return i;
      return -1;
    }

    /* Tags like "z00to10" or "z00to02" */
    inline std::string vzH3RangeTag(int i) const
    {
      return Form("z%02dto%02d", (int)std::lround(vzH3Edges[i]),
                                  (int)std::lround(vzH3Edges[i+1]));
    }
    inline std::string vzH3FineRangeTag(int i) const
    {
      return Form("z%02dto%02d", (int)std::lround(vzH3FineEdges[i]),
                                  (int)std::lround(vzH3FineEdges[i+1]));
    }


    /* Hist arrays: one uncorrected TH3 per |z| bin (coarse and fine) */
    TH3F* h3_blockCoord_E_vz      [N_VzH3Bins]     {};
    TH3F* h3_blockCoord_E_vz_fine [N_VzH3FineBins]{};
    
    /* NEW: Cross-slice variants — per-η-view × per-|z| (coarse and fine) */
    /* view index: 0=fullEta (|η|≤1.10), 1=etaCore (|η|≤0.20),
       2=etaMid (0.20<|η|≤0.70), 3=etaEdge (0.70<|η|≤1.10) */
    TH3F* h3_blockCoord_E_eta_vz      [kNEtaViews][N_VzH3Bins]{};
    TH3F* h3_blockCoord_E_eta_vz_fine [kNEtaViews][N_VzH3FineBins]{};

    /* NEW: per-|z_vtx| b-tables and readiness flags (coarse and fine) */
    std::array<std::array<float, N_Ebins>, N_VzH3Bins>       m_bValsPhi_vz{};
    std::array<std::array<float, N_Ebins>, N_VzH3Bins>       m_bValsEta_vz{};
    std::array<std::array<float, N_Ebins>, N_VzH3FineBins>   m_bValsPhi_vz_fine{};
    std::array<std::array<float, N_Ebins>, N_VzH3FineBins>   m_bValsEta_vz_fine{};
    std::array<bool, N_VzH3Bins>       m_bPhiReady_vz{{}};
    std::array<bool, N_VzH3Bins>       m_bEtaReady_vz{{}};
    std::array<bool, N_VzH3FineBins>   m_bPhiReady_vz_fine{{}};
    std::array<bool, N_VzH3FineBins>   m_bEtaReady_vz_fine{{}};

    /* NEW: corrected counterparts per |z| slice */
    TH3F* h3_blockCoord_E_vz_corr      [N_VzH3Bins]     {};
    TH3F* h3_blockCoord_E_vz_fine_corr [N_VzH3FineBins]{};

    /* NEW: corrected counterparts for cross-slice η-view × |z| (coarse & fine) */
    TH3F* h3_blockCoord_E_eta_vz_corr      [kNEtaViews][N_VzH3Bins]     {};
    TH3F* h3_blockCoord_E_eta_vz_fine_corr [kNEtaViews][N_VzH3FineBins] {};

    /* NEW: b-tables for η×z (coarse & fine) */
    std::array<std::array<std::array<float, N_Ebins>, N_VzH3Bins>,     kNEtaViews> m_bValsPhi_eta_vz{};
    std::array<std::array<std::array<float, N_Ebins>, N_VzH3Bins>,     kNEtaViews> m_bValsEta_eta_vz{};
    std::array<std::array<std::array<float, N_Ebins>, N_VzH3FineBins>, kNEtaViews> m_bValsPhi_eta_vz_fine{};
    std::array<std::array<std::array<float, N_Ebins>, N_VzH3FineBins>, kNEtaViews> m_bValsEta_eta_vz_fine{};

    std::array<std::array<bool, N_VzH3Bins>,     kNEtaViews> m_bPhiReady_eta_vz{};
    std::array<std::array<bool, N_VzH3Bins>,     kNEtaViews> m_bEtaReady_eta_vz{};
    std::array<std::array<bool, N_VzH3FineBins>, kNEtaViews> m_bPhiReady_eta_vz_fine{};
    std::array<std::array<bool, N_VzH3FineBins>, kNEtaViews> m_bEtaReady_eta_vz_fine{};

    // Global QA (REPLACED): signed incidence vs z_vtx and η_SD
    TH3F* h3_alphaPhi_vsVz_vsEta {nullptr};
    TH3F* h3_alphaEta_vsVz_vsEta {nullptr};
    
    TH3F* h3_alphaPhi_vsVz_vsEtaDet {nullptr};
    TH3F* h3_alphaEta_vsVz_vsEtaDet {nullptr};
    // ---------- Step (3): incidence‑aware per‑energy maps ----------
    // α binning: 0 … 0.30 rad in 180 bins (≈1.667 mrad/bin)
    // X^{meas} binning: −0.5 … +0.5 in 200 bins (0.005 tower/bin)
    static constexpr int   NAlphaBins = 180;
    static constexpr float AlphaMin   = 0.0f;
    static constexpr float AlphaMax   = 0.30f;
    static constexpr int   NXmeasBins = 200;
    static constexpr float XmeasMin   = -0.5f;
    static constexpr float XmeasMax   = +0.5f;

    // Per‑E slice, φ‑view: Xφ(meas) vs αφ  and QA spectra/profiles
    TH2F*     h2_XmeasPhi_vsAlpha_E[N_Ebins] {};
    TH1F*     h_alphaPhi_E          [N_Ebins] {};
    
    TH1F*    h_alphaPhi_sgn_mean_fold_core_E[N_Ebins]{};
    TH1F*    h_alphaPhi_sgn_mean_fold_mid_E [N_Ebins]{};
    TH1F*    h_alphaPhi_sgn_mean_fold_edge_E[N_Ebins]{};
    TH1F*    h_alphaPhi_sgn_mean_fold_full_E[N_Ebins]{};
    
    TProfile* p_secAlpha_phi_E      [N_Ebins] {};

    // Per‑E slice, η‑view: Xη(meas) vs αη  and QA spectra/profiles
    TH2F*     h2_XmeasEta_vsAlpha_E[N_Ebins] {};
    TH1F*     h_alphaEta_E          [N_Ebins] {};
    TProfile* p_secAlpha_eta_E      [N_Ebins] {};
    
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
 
