#ifndef PositionDependentCorrection_H__
#define PositionDependentCorrection_H__
#include <array>
#include <fun4all/SubsysReco.h>
#include <string>
#include <vector>
#include <map>
#include <calotrigger/TriggerAnalyzer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>   // ← gives access to ConstIterator
#include <cmath>
#include <TString.h>
// CLHEP types appear in inline code inside the header
#include <CLHEP/Vector/ThreeVector.h>   // Hep3Vector
#include <functional>
#include <atomic>
#include <TH3F.h>
#include "PDCGeo.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/BEmcRec.h"
#include "PDCGeo.h"

class Fun4AllHistoManager;
class PHCompositeNode;
class RawTowerGeomContainer;
class RawClusterContainer;
class BEmcRecCEMC;
class TFile;
class TNtuple;
class TTree;
class TH1;
class TH1F;
class TH3;
class TH3F;
class TH2;
class TH2F;
class TF1;
class TProfile2D;
class TProfile;
class TLorentzVector;
class TRandom3;
class RawCluster;
class PositionDependentCorrection : public SubsysReco

{
 public:
  // --------------------------------------------------------------------
  // Constructor and destructor
  // --------------------------------------------------------------------
  PositionDependentCorrection(const std::string &name = "PositionDependentCorrection",
                              const std::string &fname = "MyNtuple.root");
  virtual ~PositionDependentCorrection();
  void bookCommonHistograms(const std::function<std::string(int)>& makeLabel);
  void bookSimulationHistograms(const std::function<std::string(int)>& makeLabel);
  int Init(PHCompositeNode *);
  int process_event(PHCompositeNode *);
  int End(PHCompositeNode *);

  void setIsSimulation(bool flag)   { m_isSimulation = flag; }
  bool isSimulation()         const { return m_isSimulation; }
  void setMassFitsDone(bool flag = true) { m_massFitsDone = flag; }
    
  int process_g4hits(PHCompositeNode *);
  int process_g4cells(PHCompositeNode *);
  int process_towers(PHCompositeNode *);
  int process_clusters(PHCompositeNode *);

  // --------------------------------------------------------------------
  // Utility setters
  // --------------------------------------------------------------------
  void Detector(const std::string &name) { detector = name; }
  void set_timing_cut_width(const int &t) { _range = t; }
  void set_vertex_cut(const float &v)     { _vz = v; }
  void apply_vertex_cut(bool Vtx_cut)     { m_vtxCut = Vtx_cut; }

  TF1* fitHistogram(TH1* h);
  void fitEtaSlices(const std::string& infile,
                    const std::string& fitOutFile,
                    const std::string& cdbFile);
    
  void setBEmcRec(BEmcRecCEMC* bemcptr)
  {
      m_bemcRec = bemcptr;
  }
  bool alreadyDeclaredHistograms = false;
    
  void UseSurveyGeometry(bool v = true) { m_useSurveyGeometry = v; }
    
  enum class EBinningMode { kRange,
        kDiscrete
  };
    
  void setBinningMode(EBinningMode mode) { m_binningMode = mode; }
  EBinningMode getBinningMode() const    { return m_binningMode; }
  static uint64_t getVerbosityLevel() { return s_verbosityLevel.load(); }

  static inline float doPhiBlockCorr(float u, float b)
  { return PDC::Geo::phi::undoAsh(u, b); }

  static inline float doEtaBlockCorr(float u, float b, float d = 1.f)
  { return PDC::Geo::eta::undoAsh(u, b, d); }

  inline float convertBlockToGlobalPhi(int blk, float loc) const
  {
      /* honour run‑time offset if it was measured */
      struct Off { const PositionDependentCorrection* self; float operator()() const { return self->m_hasOffset ? self->m_phi0Offset : 0.F; } };
      return PDC::Geo::phi::blockToGlobal(blk, loc, Off{this});
  }

  inline float convertBlockToGlobalEta(int blk, float loc) const
  {
      struct Off { const PositionDependentCorrection* self; float operator()() const { return self->m_hasEtaOffset ? self->m_eta0Offset : 0.F; } };
      return PDC::Geo::eta::blockToGlobal(blk, loc, Off{this});
  }

  /* --------------------------------------------------------------------
   * 2) fine‑index helpers (unchanged public signature)
   * -------------------------------------------------------------------- */
  static inline int finePhiIndex(int blk, float loc)
  { return PDC::Geo::finePhiIdx(blk, loc); }

  static inline int fineEtaIndex(int blk, float loc)
  { return PDC::Geo::fineEtaIdx(blk, loc); }
    
 protected:
  /*! Current verbosity that static helpers can consult.                     *
    *  It is updated for every instance in the constructor.                   */
  static std::atomic<uint64_t> s_verbosityLevel;
  /** Measure the rigid φ–offset (barrel tilt) once per job.
    *
    *  – Uses `m_geometry`, fills `m_phi0Offset`, sets `m_hasOffset`.
    *  – Returns *true* on success, *false* if the geometry container is empty.
  */
  bool computeRigidPhiOffset();
  float m_eta0Offset {0.f};
  bool  m_hasEtaOffset {false};
  bool  computeRigidEtaOffset();
    
  // --------------------------------------------------------------------
  // 2) Method prototypes that reference RawClusterContainer, etc.
  // --------------------------------------------------------------------
  float retrieveVertexZ(PHCompositeNode* topNode);

  void fillTowerInfo(PHCompositeNode* topNode,
                     float emcal_hit_threshold,
                     float &tower_tot_e,
                     std::vector<float> &ht_eta,
                     std::vector<float> &ht_phi);

  void fillTruthInfo(PHCompositeNode* topNode,
                     float &vtx_z,
                     std::vector<TLorentzVector> &truth_photons,
                     std::vector<TLorentzVector> &truth_meson_photons);
  BEmcRecCEMC* m_bemcRec = nullptr;
  //  <-- We forward-declared RawTowerGeomContainer & RawClusterContainer above
  RawTowerGeomContainer* checkTowerGeometry(PHCompositeNode* topNode);
  RawClusterContainer* retrieveClusterContainer(PHCompositeNode* topNode);

  int countClusters(RawClusterContainer* clusterContainer,
                    float vtx_z,
                    float nClus_ptCut,
                    float clus_chisq_cut,
                    int &nClusContainer);

  // ───────────────────────────────────────────────────────────────
  // Forward Ash distortion  ( true  →  measured )
  //
  //   PDC::Geo::phi::undoAsh() converts   measured → true
  //   therefore,  true → measured  is obtained with  (‑b).
  //   The sign flip is the only difference.
  //
  //   • keeps all folding & edge fixes identical to the library
  //   • zero duplicate maintenance burden
  // ───────────────────────────────────────────────────────────────
  inline float
  doAshShift(float localPhi, float b)  // header unchanged
  {
        // direct call into the constexpr helper – *no* extra branches
        return PDC::Geo::phi::undoAsh(localPhi, -b);
  }

  float doLogWeightCoord(const std::vector<int>& towerphis,
                           const std::vector<float>& towerenergies,
                           float w0);
    

  void fillBlockCoordinateHistograms(const std::pair<float,float>& blkCoord,
                                       int   blkEtaCoarse,
                                       int   blkPhiCoarse,
                                       float clusE,
                                       int   iEbin,
                                       std::size_t nTowers);

  void processSimulationTruthMatches(
        RawCluster*                     clus1,
        const TLorentzVector&           photon1,
        int                             lt_eta,
        int                             lt_phi,
        const std::pair<float,float>&   blkCoord,
        int                             blkEtaCoarse,
        int                             blkPhiCoarse,
        const std::vector<int>&         towerPhis,
        const std::vector<float>&       towerEs,
        float                           vtx_z,
        const std::vector<TLorentzVector>& truth_photons,
        const std::vector<TLorentzVector>& truth_meson_photons,
        bool&                           match1,
        TLorentzVector&                 ph1_trEtaPhi);
    
    
 void processClusterPairs(
          RawClusterContainer*               clusterContainer,
          RawClusterContainer::ConstIterator cIt1,
          const CLHEP::Hep3Vector&           vertex,
          const TLorentzVector&              photon1,
          float                              clusE,
          int                                lt_eta,
          int                                lt_phi,
          int                                blkEtaCoarse,
          int                                blkPhiCoarse,
          const std::pair<float,float>&      blkCoord,
          float                              maxAlpha,
          float                              ptMaxCut,
          float                              pt2ClusCut,
          float                              pi0ptcut,
          float                              weight,
          bool                               match1,
          const TLorentzVector&              ph1_trEtaPhi,
          const std::vector<TLorentzVector>& truth_meson_photons,
          bool                               isSimulation);              // << NEW
    
  void finalClusterLoop(PHCompositeNode* topNode,
                        RawClusterContainer* clusterContainer,
                        float vtx_z,
                        const std::vector<TLorentzVector> &truth_photons,
                        const std::vector<TLorentzVector> &truth_meson_photons,
                        float tower_tot_e,
                        float max_nClusCount,
                        int nClusCount,
                        float maxAlpha,
                        float ptMaxCut,
                        float pt1ClusCut,
                        float pt2ClusCut,
                        float pi0ptcut,
                        float weight);

  void fillAshLogDx(
        RawCluster* cluster,
        const TLorentzVector& recoPhoton,
        const TLorentzVector& truthPhoton,
        const std::pair<float, float>& blockCord,
        int blockPhiBin,
        const std::vector<int>& tower_phis,
        const std::vector<float>& tower_energies
  );
    

    
  int  getEnergySlice(float E) const;
    
  // --------------------------------------------------------------------
  // 3) Data members
  // --------------------------------------------------------------------
  std::string detector;
  /* ----------------------------------------------------------------
   *   –  m_vzTightCut  … “physics” histograms (default 10 cm)
   *   –  m_vzSliceMax  … upper edge of the vzEdge table   (30 cm)
   * -------------------------------------------------------------- */
  float m_vzTightCut {10.f};
  float m_vzSliceMax {30.f};
  /// steer the tight cut from a macro if desired
  void setTightVzCut(float cm) { m_vzTightCut = std::fabs(cm); }
    
  static constexpr std::array<float,7> vzEdge = { 0, 5, 10, 15, 20, 25, 30};
  static constexpr int N_VzBins = vzEdge.size() - 1;

  std::string outfilename;
  bool m_isSimulation = false;      ///< true = MC, false = real data
  int Getpeaktime(TH1 *h);
  Fun4AllHistoManager *hm = nullptr;
  TFile *outfile = nullptr;
  bool  isFitDoneForPhi = true;
  bool  isFitDoneForEta = true;
    
  bool m_useSurveyGeometry {false};   ///< load φ-tilt from CDB?
  uint64_t m_timeStamp      {0};      ///< can be set by macro if desired
  std::string m_cdbTag      {"MDC2"};
    
  EBinningMode m_binningMode { EBinningMode::kRange };
  static constexpr float kDiscreteTol = 0.25f;

  static constexpr double kTolFactor = 1.5;   ///< Nσ
  static constexpr double kTolMinGeV = 0.20;  ///< lower bound [GeV]
    
  inline std::string sliceTag(int i) const             // extra scope removed
  {
      return (m_binningMode == EBinningMode::kRange)
             ? Form("%.0f_%.0f", eEdge[i], eEdge[i+1])
             : Form("E%.0f",      eEdge[i]);
  }
    
  static constexpr int    N_Ebins = 8;          // 9 edges → 8 slices

  // master edges [GeV]
  static constexpr double eEdge[N_Ebins + 1] =
      { 2, 4, 6, 8, 10, 12, 15, 20, 30 };

  static constexpr std::array<float, N_Ebins> expectedLo =
  []{
      std::array<float, N_Ebins> lo{};
      for (int i = 0; i < N_Ebins; ++i) lo[i] = static_cast<float>(eEdge[i]);
      return lo;
  }();

  static constexpr std::array<float, N_Ebins> expectedHi =
  []{
      std::array<float, N_Ebins> hi{};
      for (int i = 0; i < N_Ebins; ++i) hi[i] = static_cast<float>(eEdge[i+1]);
      return hi;
  }();
    
  /* ---------- π0‑mass‑window (pass‑2 helper) ------------------------------ */
  static constexpr int kMaxEBins = N_Ebins;
  struct MassWindow { float mu{0.f}; float sigma{0.f}; };

  bool        m_massFitsDone {false};                // external flag
  MassWindow  m_winRaw [kMaxEBins] {};               // RAW  (unused for now)
  MassWindow  m_winCorr[kMaxEBins] {};               // CORR (unused for now)

  /* ---------------- helper loaders --------------------------------------- */
  bool loadBValues          (const std::string& path);   // b‑parameters
  bool loadMassWindowTable  (const std::string& path);   // μ/σ table
    
  /* 5-way residual helpers (Δφ / Δη) -- NEW */
  void fillDPhiAllVariants( RawCluster*                cluster,
                              const TLorentzVector&      recoPhoton,
                              const TLorentzVector&      truthPhoton,
                              const std::pair<float,float>& blkCoord,
                              int                        blockPhiBin,
                              float                      vtx_z,
                              TH1F*                      cpRawHistArr [N_Ebins],
                              TH1F*                      cpCorrHistArr[N_Ebins] );

  void fillDEtaAllVariants( RawCluster*                cluster,
                              const TLorentzVector&      recoPhoton,
                              const TLorentzVector&      truthPhoton,
                              const std::pair<float,float>& blkCoord,
                              int                        blockEtaBin,
                              int                        blockPhiBin,
                              float                      vtx_z,
                              TH1F*                      cpRawHistArr [N_Ebins],
                              TH1F*                      cpCorrHistArr[N_Ebins],
                              bool                       fillGlobal = true );


  std::array<std::atomic<std::uint64_t>,4> m_blkLocCount {{0,0,0,0}};
  float m_bValsPhi[N_Ebins]{};
  float m_bValsEta[N_Ebins]{};
  std::vector<double> m_bScan;
  std::vector<double> m_w0Scan;
  TriggerAnalyzer* trigAna{nullptr};
  TH1F* h_localPhi_corrected[N_Ebins];
  TH1F* h_localEta_corrected[N_Ebins];
  TH1F* h_phi_diff_raw_E[N_Ebins];
  TH1F* h_phi_diff_corrected_E[N_Ebins];
  TH1F* h_eta_diff_raw_E[N_Ebins];
  TH1F* h_eta_diff_corrected_E[N_Ebins];
    
  TH1F* h_phi_diff_cpRaw_E     [N_Ebins]{};
  TH1F* h_phi_diff_cpCorr_E    [N_Ebins]{};
  TH1F* h_phi_diff_cpBcorr_E [N_Ebins] {};
    
  TH1F* h_eta_diff_cpRaw_E     [N_Ebins]{};
  TH1F* h_eta_diff_cpCorr_E    [N_Ebins]{};
  TH1F* h_eta_diff_cpBcorr_E   [N_Ebins] {};
    
  TH1F* h_eta_diff_raw_E_vz     [N_Ebins][N_VzBins] {};
  TH1F* h_eta_diff_corrected_E_vz[N_Ebins][N_VzBins]{};

  TH1F* h_eta_diff_cpRaw_E_vz   [N_Ebins][N_VzBins]{};
  TH1F* h_eta_diff_cpCorr_E_vz  [N_Ebins][N_VzBins]{};
  TH1F* h_eta_diff_cpBcorr_E_vz [N_Ebins][N_VzBins]{};
    
  // –– NEW –– global tallies (initialised to zero in ctor)
  mutable std::atomic<std::uint64_t> m_nWinRAW   {0};
  mutable std::atomic<std::uint64_t> m_nWinCP    {0};
  mutable std::atomic<std::uint64_t> m_nWinBCorr {0};
    
  mutable std::atomic<std::uint64_t> m_nWinRAW_Eta   {0};
  mutable std::atomic<std::uint64_t> m_nWinCP_Eta    {0};
  mutable std::atomic<std::uint64_t> m_nWinBCorr_Eta {0};
    
 /* 5‑way Δφ win‑counters */
  std::atomic<std::uint64_t> m_phiWinCLUSraw   {0};
  std::atomic<std::uint64_t> m_phiWinCLUScp    {0};
  std::atomic<std::uint64_t> m_phiWinCLUSbcorr {0};
  std::atomic<std::uint64_t> m_phiWinPDCraw    {0};
  std::atomic<std::uint64_t> m_phiWinPDCcorr   {0};

  std::atomic<std::uint64_t> m_etaWinCLUSraw   {0};
  std::atomic<std::uint64_t> m_etaWinCLUScp    {0};
  std::atomic<std::uint64_t> m_etaWinCLUSbcorr {0};
  std::atomic<std::uint64_t> m_etaWinPDCraw    {0};
  std::atomic<std::uint64_t> m_etaWinPDCcorr   {0};
    
    
  /* optional numerical‑anomaly monitors */
  std::atomic<std::uint64_t> g_nanPhi {0};
  std::atomic<std::uint64_t> g_nanEta {0};
  std::atomic<std::uint64_t> g_nCorrPhiRight {0};  
    
  TProfile* pr_phi_vs_blockcoord = nullptr;
  TH2* h_emcal_mbd_correlation = nullptr;
  TH2* h_ohcal_mbd_correlation = nullptr;
  TH2* h_ihcal_mbd_correlation = nullptr;
  TH2* h_emcal_hcal_correlation = nullptr;
  TH2* h_emcal_zdc_correlation = nullptr;
  TH2* h_clusE_nTow{nullptr};
  TH1* h_InvMass = nullptr;
  TH1* h_InvMass_w = nullptr;
  TH1* h_InvMassMix = nullptr;
  TH2* h_cemc_etaphi = nullptr;
  TH2* h_hcalin_etaphi = nullptr;
  TH2* h_hcalout_etaphi = nullptr;
  TH2* h_cemc_etaphi_wQA = nullptr;
  TH2* h_hcalin_etaphi_wQA = nullptr;
  TH2* h_hcalout_etaphi_wQA = nullptr;
  TH1* h_totalzdc_e;
  TH1* h_delR_recTrth = nullptr;
  TProfile2D* h_cemc_etaphi_time = nullptr;
  TProfile2D* h_hcalin_etaphi_time = nullptr;
  TProfile2D* h_hcalout_etaphi_time = nullptr;
  TProfile2D* h_cemc_etaphi_badChi2 = nullptr;
  TProfile2D* h_hcalin_etaphi_badChi2 = nullptr;
  TProfile2D* h_hcalout_etaphi_badChi2 = nullptr;
  TH1* hzdctime;
  TH1* hmbdtime;
  TH1* hemcaltime;
  TH1* hihcaltime;
  TH1* hohcaltime;
  TH1* hzdctime_cut;
  TH1* hmbdtime_cut;
  TH1* hemcaltime_cut;
  TH1* hihcaltime_cut;
  TH3* h3_cluster_block_cord_E_corrected;
  TH1* hohcaltime_cut;
  TH1* h_tower_e;
  TH1* hvtx_z_raw;
  TH1* hvtx_z_cut;
  TH1* hzdcSouthraw;
  TH1* hzdcNorthraw;
  TH1* hzdcSouthcalib;
  TH1* hzdcNorthcalib;
  TH1* h_clusE;
  TH2* h_etaphi_clus;
  TNtuple *g4hitntuple = nullptr;
  TNtuple *g4cellntuple = nullptr;
  TTree *towerntuple = nullptr;
  TNtuple *clusterntuple = nullptr;
  std::vector<float> m_energy;
  std::vector<int> m_etabin;
  std::vector<int> m_phibin;
  std::vector<int> m_time;
  std::vector<float> m_hcalin_energy;
  std::vector<int> m_hcalin_etabin;
  std::vector<int> m_hcalin_phibin;
  std::vector<int> m_hcalin_time;
  std::vector<float> m_hcalout_energy;
  std::vector<int> m_hcalout_etabin;
  std::vector<int> m_hcalout_phibin;
  std::vector<int> m_hcalout_time;
  std::vector<float> m_zdc_energy;
  std::vector<int> m_zdc_index;
  std::vector<int> m_zdc_side;
  std::vector<float> m_bbc_energy;
  std::vector<int> m_bbc_type;
  std::vector<int> m_bbc_side;
  RawTowerGeomContainer* m_geometry = nullptr;

  int _eventcounter;
  int _range = 1;
  bool m_vtxCut = true;
  bool dynMaskClus = false;
  float _vz = 0;
  bool debug = false;
    
  float m_phi0Offset  = 0.f;   ///< measured once at run time (radians)
  bool  m_hasOffset   = false; ///< guard so we do the measurement only once
  TH1* h_pt1;
  TH1* h_pt2;
  TH1* h_nclusters;
  TH2* h_mass_eta_lt;
  TH2* h_mass_eta_lt_rw;
  TH2* h_pt_eta;
  TH2* h_pt_eta_rw;
  TH1* h_emcal_e_eta;
  TH1* h_truth_eta;
  TH1* h_truth_e;
  TH1F* h_truth_vz       {nullptr};   // truth vz   distribution
  TH1F* h_reco_vz        {nullptr};   // reconstructed vz distribution
  TH2F* h2_truthReco_vz  {nullptr};   // 2-D truth vs reco correlation
  TH1* h_pt_rw[96];
    
  TFile* frw;
  TH1* h_matched_res;
  TH1* h_res_e;
  TH3* h_res_e_phi;
  TH3* h_res_e_eta;
  TH3* h_res_e_eta_pdc;
  TH1* h_res;
  TH3* h_m_pt_eta;
  TH3* h_m_ptTr_eta;
  TH3* h_m_ptTr_eta_trKin;
  TH3* h_delPhi_e_eta;
  TH3* h_delEta_e_eta;
  TH3* h_delR_e_eta;
  TH3* h_delPhi_e_phi;
  TH2F* h_mE_raw  {nullptr};
  TH2F* h_mE_corr {nullptr};
  TH3F* h_m_blk_raw  {nullptr};
  TH3F* h_m_blk_corr {nullptr};
  TProfile* pr_eta_shower;
  TProfile* pr_phi_shower;
  TH2* h_vert_xy;
  TH1* h_truthE;

  TH2F*       h2_chi2_tot_etaPhi   {nullptr};
  TH2F*       h2_chi2_rej_etaPhi   {nullptr};
  TProfile2D* p_chi2_pass_etaPhi  {nullptr};

  static const int NBinsBlock = 14;
  TH2* h_mass_block_pt[NBinsBlock][NBinsBlock];
  TH2* h_res_block_E[NBinsBlock][NBinsBlock];
  TH3* h3_cluster_block_cord_E;
  TH1* h_block_phi;
  TH1* h_block_eta;
  TH1* h_clus_E_size;
  TH1* h_block_bin;

    
  inline int getVzSlice(float vz) const
  {
        for (int i = 0; i < N_VzBins; ++i)
            if (vz >= vzEdge[i] && vz < vzEdge[i + 1]) return i;
        return -1;
  }

  std::map<std::string, std::string> triggerNameMap = {
        {"MBD N&S >= 1",          "MBD_NandS_geq_1"}
//        {"Photon 3 GeV + MBD NS >= 1","Photon_3_GeV_plus_MBD_NS_geq_1"},
//        {"Photon 4 GeV + MBD NS >= 1","Photon_4_GeV_plus_MBD_NS_geq_1"},
//        {"Photon 5 GeV + MBD NS >= 1","Photon_5_GeV_plus_MBD_NS_geq_1"}
  };
  std::map<int, std::string>* activeTriggerNameMap = nullptr;

  float target_pi0_mass = 0.145;
  TRandom3* rnd;
};

#endif  // PositionDependentCorrection_H__
