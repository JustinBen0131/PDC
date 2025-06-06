#ifndef PositionDependentCorrection_H__
#define PositionDependentCorrection_H__

#include <array>
#include <fun4all/SubsysReco.h>  // we do need this
#include <string>                // for std::string
#include <vector>                // for std::vector
#include <map>                   // for std::map
#include <calotrigger/TriggerAnalyzer.h>
#include <calobase/RawCluster.h>

// --------------------------------------------------------------------
// 1) Forward declarations for classes we only use as pointers/references
// --------------------------------------------------------------------
class Fun4AllHistoManager;
class PHCompositeNode;
class RawTowerGeomContainer;   // <-- ADDED
class RawClusterContainer;     // <-- ADDED
class BEmcRecCEMC;

class TFile;
class TNtuple;
class TTree;
class TH1;
class TH1F;
class TH3;
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

  // --------------------------------------------------------------------
  // Fun4All methods
  // --------------------------------------------------------------------
  int Init(PHCompositeNode *);
  int process_event(PHCompositeNode *);
  int End(PHCompositeNode *);

  void setIsSimulation(bool sim) { isSimulation = sim; }
    
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

  // --------------------------------------------------------------------
  // Additional user methods
  // --------------------------------------------------------------------
  float getWeight(int ieta, float pt);

  TF1* fitHistogram(TH1* h);
  void fitEtaSlices(const std::string& infile,
                    const std::string& fitOutFile,
                    const std::string& cdbFile);
    
  void setBEmcRec(BEmcRecCEMC* bemcptr)
  {
      m_bemcRec = bemcptr;
  }
  bool alreadyDeclaredHistograms = false;


 protected:
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
    
  float doPhiBlockCorr(float localPhi, float bphi);
  float doEtaBlockCorr(float localEta, float bEta);
    
  float convertBlockToGlobalPhi(int block_phi_bin, float localPhi);
  float convertBlockToGlobalEta(int block_eta_bin, float localEta);
    
    // PositionDependentCorrection.h
    //----------------------------------------------------
    float  phiAtShowerDepth( float  energy,
                             double rFront,
                             double zFront,
                             float  phiFront,
                             int    ix,          ///< lead-tower fine φ-index
                             int    iy ) const;  ///< lead-tower fine η-index

    double xAtShowerDepth ( float  energy,
                             double rFront,
                             double zFront,
                             float  phiFront,
                             int    ix,
                             int    iy ) const;
    //----------------------------------------------------
    

  float doAshShift(float localPhi, float bVal);
  float doLogWeightCoord(const std::vector<int>& towerphis,
                           const std::vector<float>& towerenergies,
                           float w0);

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

    // ────────────────────────────────────────────────────────────
    //  Public interface
    // ────────────────────────────────────────────────────────────
    void fillDPhiRawAndCorrected( RawCluster*            cluster,
                                  const TLorentzVector&  recoPhoton,
                                  const TLorentzVector&  truthPhoton,
                                  const std::pair<float,float>& blkCoord,
                                  int   blockPhiBin,
                                  float rawDelPhi /* unused – kept for call-site compatibility */ );

  void fillDEtaRawAndCorrected(
        const TLorentzVector& recoPhoton,
        const TLorentzVector& truthPhoton,
        const std::pair<float,float>& blockCord, // local block coords in [0,1]
        int blockEtaBin,
        float delEta
  );

  // --------------------------------------------------------------------
  // 3) Data members
  // --------------------------------------------------------------------
  std::string detector;
  std::string outfilename;
  int Getpeaktime(TH1 *h);
  Fun4AllHistoManager *hm = nullptr;
  TFile *outfile = nullptr;
  bool  isFitDoneForPhi = true;
  bool  isFitDoneForEta = true;
    //--------------------------------------------------------------------
    //  Energy–slice constants  (barrel EMCAL photons)
    //--------------------------------------------------------------------
    static constexpr int    N_Ebins = 8;          // 9 edges → 8 slices

    // master edges [GeV]
    static constexpr double eEdge[N_Ebins + 1] =
        { 2, 4, 6, 8, 10, 12, 15, 20, 30 };

    // derive {E_low} and {E_high} once, at compile-time, via a lambda
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

    // everything that used to be length-4 must now be length-8
    float m_bValsPhi[N_Ebins]{};
    float m_bValsEta[N_Ebins]{};
    
  std::vector<double> m_bScan;
  std::vector<double> m_w0Scan;
  TriggerAnalyzer* trigAna{nullptr};
  TH1F* h_localPhi_corrected[N_Ebins];
  TH1F* h_localEta_corrected[N_Ebins];
  TH1F* h_phi_diff_raw_E[N_Ebins];         // raw Δφ for each bin
  TH1F* h_phi_diff_corrected_E[N_Ebins];   // corrected Δφ for each bin
  TH1F* h_eta_diff_raw_E[N_Ebins];
  TH1F* h_eta_diff_corrected_E[N_Ebins];
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
  bool getVtx = true;
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
  TH1* h_truth_pt;
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
  TProfile* pr_eta_shower;
  TProfile* pr_phi_shower;
  TH2* h_vert_xy;
  TH1* h_truthE;
    
    //  --- χ² QA ------------------------------------------------------------
    //  Denominator = all clusters (before χ² cut)
    //  Numerator   = clusters rejected by χ² cut
    //  Profile     = ⟨pass flag⟩  →   1 = kept    0 = rejected
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
  bool isSimulation = false;

  float getAvgEta(const std::vector<int> &toweretas,
                  const std::vector<float> &towerenergies);
  float getAvgPhi(const std::vector<int> &towerphis,
                  const std::vector<float> &towerenergies);

    std::pair<float,float> getBlockCord(const std::vector<int>&   towerEtas,
                                        const std::vector<int>&   towerPhis,
                                        const std::vector<float>& towerEs,
                                        int&                      blkPhiOut);
  std::map<std::string, std::string> triggerNameMap = {
    {"MBD N&S >= 1", "MBD_NandS_geq_1"}
  };
  std::map<int, std::string>* activeTriggerNameMap = nullptr;

  float target_pi0_mass = 0.145;
  TRandom3* rnd;
};

#endif  // PositionDependentCorrection_H__
