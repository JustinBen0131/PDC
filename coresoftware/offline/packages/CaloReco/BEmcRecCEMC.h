#ifndef CALORECO_BEMCRECCEMC_H
#define CALORECO_BEMCRECCEMC_H

#include "BEmcRec.h"

#include <string>
#include <vector>

// class BEmcProfile;
class EmcModule;

class BEmcRecCEMC : public BEmcRec
{
 public:
  BEmcRecCEMC();
  ~BEmcRecCEMC() override = default;
  void CorrectEnergy(float energy, float x, float y, float &ecorr) override;
  void CorrectECore(float ecore, float x, float y, float &ecorecorr) override;
  void CorrectPosition(float energy, float x, float y, float &xcorr, float &ycorr) override;

  // EA (baseline, already present)
  void CorrectPositionEnergyAware(float Energy, float x, float y,
                                    float& xc, float& yc,
                                    float* out_bphi = nullptr,
                                    float* out_beta = nullptr);
    
  void CorrectPositionEnergyAwareZVTXEtaAndEnergyDep(float Energy, float vtxZ_cm,
                                                       float x, float y,
                                                       float& xc, float& yc);

  void CorrectPositionEnergyAwareZVTXAndEnergyDep(float Energy, float vtxZ_cm, float x, float y,
                                                    float& xc, float& yc);

  void CorrectPositionEnergyAwareEtaAndEnergyDep(float Energy, float x, float y,
                                                   float& xc, float& yc);

  void CorrectPositionEnergyAwareEnergyDepAndIncidentAngle(float Energy, float x, float y,
                                                             float& xc, float& yc);
    
  void CorrectPositionEnergyAwareEnergyDepOnly(float Energy, float x, float y,
                                                 float& xc, float& yc);

  void CorrectShowerDepth(int ix, int iy, float energy, float x, float y, float z, float &xc, float &yc, float &zc) override;


  // ---- PUBLIC (add near other getters) ----
  float lastSignedAlphaPhi() const { return m_lastAlphaPhiSigned; }
  float lastSignedAlphaEta()  const { return m_lastAlphaEtaSigned; }


  void LoadProfile(const std::string &fname) override;
  //  float GetProb(std::vector<EmcModule> HitList, float e, float xg, float yg, float zg, float &chi2, int &ndf) override;
    
  void GetImpactThetaPhi(float xg, float yg, float zg, float &theta, float &phi) override;
    
    struct PhiTiltAB { double a{8.654924e-04}; double b{8.490399e-04}; };

    // Unified φ(E) API — used by CorrectShowerDepth and ComputeIncidenceSD
    inline void SetPhiTiltCoeffs(double a, double b) {  // backward-compat: set all three to same
      m_phiTiltAB_core = {a,b}; m_phiTiltAB_mid = {a,b}; m_phiTiltAB_edge = {a,b};
    }
    inline void SetPhiTiltCoeffsByEta(double a_core,double b_core,
                                      double a_mid ,double b_mid ,
                                      double a_edge,double b_edge) {
      m_phiTiltAB_core = {a_core,b_core};
      m_phiTiltAB_mid  = {a_mid ,b_mid };
      m_phiTiltAB_edge = {a_edge,b_edge};
    }
    inline void EnablePhiTilt(bool on) { m_enablePhiTilt = on; }

    // ---------- NEW: constant φ offsets (rad) per |η| band ----------
    inline void SetPhiTiltOffsets(double phi0_all) {
      m_phi0_core = m_phi0_mid = m_phi0_edge = phi0_all;
    }
    inline void SetPhiTiltOffsetsByEta(double phi0_core, double phi0_mid, double phi0_edge) {
      m_phi0_core = phi0_core; m_phi0_mid = phi0_mid; m_phi0_edge = phi0_edge;
    }

   private:
    // Per-|η| band coefficients for 〈αφ^fold〉(E) = a − b ln E  [radians]
    // |η| ≤ 0.20       → core
    // 0.20 < |η| ≤ 0.70 → mid
    // 0.70 < |η| ≤ 1.10 → edge
    PhiTiltAB m_phiTiltAB_core{ 0.14940,  -0.000238 };
    PhiTiltAB m_phiTiltAB_mid { 0.15177,  -0.000401 };
    PhiTiltAB m_phiTiltAB_edge{ 0.15244,  -0.000877 };

    bool m_enablePhiTilt{true};

    // NEW: constant φ-offsets (rad) per |η| band (defaults = 0)
    double m_phi0_core{0.0};
    double m_phi0_mid {0.0};
    double m_phi0_edge{0.0};

    // Band selector:  core (|η|≤0.20), mid (0.20<|η|≤0.70), edge (0.70<|η|≤1.10)
    inline const PhiTiltAB& phi_ab_for_absEta(double absEta) const {
      if (!std::isfinite(absEta)) return m_phiTiltAB_mid;      // conservative default
      if (absEta <= 0.20) return m_phiTiltAB_core;
      if (absEta <= 0.70) return m_phiTiltAB_mid;
      return m_phiTiltAB_edge;
    }
    inline double phi0_for_absEta(double absEta) const {
      if (!std::isfinite(absEta)) return m_phi0_mid;      // conservative default
      if (absEta <= 0.20) return m_phi0_core;
      if (absEta <= 0.70) return m_phi0_mid;
      return m_phi0_edge;
    }

    // φ(E) pre-rotation built directly from the measured drift.
    // Fit convention on the canvas: ⟨αφ^fold⟩(E) = a − b ln E  [rad],
    // with b < 0 in all three |η| bands. To cancel that drift we rotate by
    // φ(E) = φ0(|η|) + b ln(E/E0), where φ0 removes the constant offset.
    // Note: we intentionally do *not* derive the slope from (B_depth/R)*tan(a) here.
    // φ(E) pre-rotation from your *measured* drift:
    //  <αφ^fold>(E) = a - b ln E  ⇒  pre-rotate by  φ(E) = φ0(|η|) + b ln(E/E0).
    inline double phi_tilt(double E, double absEta, double /*Rcm*/) const {
      if (!m_enablePhiTilt) return 0.0;

      const auto& ab    = phi_ab_for_absEta(absEta);   // {a, b} from your fit (rad, rad/lnGeV)
      const double phi0 = phi0_for_absEta(absEta);     // constant offset (rad), one per |η| band
      const double e    = (E > 0.1) ? static_cast<double>(E) : 0.1;
      constexpr double E0_ref = 3.0;                   // pivot used on your plots

      double phi = phi0 + ab.b * std::log(e / E0_ref); // small (mrad) rotation
      if (phi >  0.01) phi =  0.01;
      if (phi < -0.01) phi = -0.01;
      return phi;
    }

    // keep angle caches
    float m_lastAlphaPhi { std::numeric_limits<float>::quiet_NaN() };
    float m_lastAlphaEta { std::numeric_limits<float>::quiet_NaN() };

    // signed caches: initialize to NaN so missing incidence does not masquerade as α=0
    float m_lastAlphaPhiSigned { std::numeric_limits<float>::quiet_NaN() };
    float m_lastAlphaEtaSigned { std::numeric_limits<float>::quiet_NaN() };

  // Average tower angle with respect to the transverse plane for each bin in eta
  // Only used when the detailed RawTowerGeom objects are not available
  double angles[96] = {2.382132067280038,2.3708136440088796,2.355383841518272,2.3434466029558525,2.32455480538339,2.3135158793606005,
                       2.2999027706894153,2.289397309387001,2.2727994372634956,2.259380682790801,2.2392521417328277,2.226335700118477,
                       2.2095791811782735,2.1935361258633024,2.1754991989110875,2.160931303742618,2.142623616095435,2.1259546238847804,
                       2.1054808114318946,2.0900705081463244,2.0670824189635053,2.050851798834259,2.026500612990249,2.005259135338716,
                       1.9868346921870874,1.966226760986188,1.9455349050274326,1.924706496569077,1.9017832210973762,1.8805537731028699,
                       1.8568682630187272,1.8337068107715466,1.8112862676694754,1.7884337950298892,1.76612732335177,1.7439109866141163,
                       1.7137229887165837,1.6887545076082588,1.6601011126718555,1.633546284539203,1.6126575096324751,1.594318479859694,
                       1.5699679790263024,1.5683722241543814,1.570109507416994,1.5686069317665148,1.5706670797900908,1.5655401209498148,
                       1.5760525326399784,1.5709255737997023,1.5729857218232786,1.5714831461727994,1.573220429435412,1.5716246745634908,
                       1.5472741737300992,1.528935143957318,1.50804636905059,1.4814915409179379,1.4528381459815345,1.4278696648732097,
                       1.3976816669756769,1.375465330238023,1.353158858559904,1.3303063859203177,1.3078858428182465,1.284724390571066,
                       1.2610388804869233,1.2398094324924172,1.2168861570207163,1.1960577485623607,1.1753658926036052,1.154757961402706,
                       1.1363335182510774,1.1150920405995441,1.0907408547555344,1.074510234626288,1.051522145443469,1.0361118421578988,
                       1.015638029705013,0.9989690374943581,0.9806613498471752,0.9660934546787058,0.948056527726491,0.9320134724115198,
                       0.9152569534713161,0.9023405118569655,0.8822119707989924,0.8687932163262975,0.852195344202792,0.8416898829003777,
                       0.8280767742291925,0.817037848206403,0.798146050633941,0.7862088120715213,0.7707790095809137,0.7594605863097555};

  double slopes_[96] = {0, 0.00023905, 4.90448e-05, 0.00027207, 4.62457e-05, 0.000252483, 4.07806e-05, 0.000248309, 3.76606e-05, 0.000257749, 3.96662e-05, 0.000268067, 3.79489e-05, 0.000258554, 4.64794e-05, 0.000252937, 4.00813e-05, 0.000256408, 5.76557e-05, 0.000254708, 6.78359e-05, 0.00025745, 9.10421e-05, 0.000249737, 9.98244e-05, 0.000248255, 0.000129006, 0.000248438, 0.000137128, 0.000252377, 0.000162495, 0.000255381, 0.000178851, 0.000255713, 0.000196615, 0.000257513, 0.000218168, 0.000260783, 0.00023512, 0.000266362, 0.000253771, 0.000258489, 0.000261129, 0.000237658, 0.000229224, 0.000219019, 0.0002103555, 0.000194438, 0.000191487, 0.000169857, 0.000163719, 0.000145249, 0.000140684, 0.00011482, 0.000114154, 0.000100124, 0.00011904, 7.79723e-05, 0.000107485, 5.65278e-05, 9.40123e-05, 3.13618e-05, 8.63381e-05, 6.76016e-06, 7.83596e-05, -1.41566e-05, 7.57772e-05, -3.53334e-05, 6.71255e-05, -5.46557e-05, 6.50159e-05, -7.18416e-05, 6.88217e-05, -8.98952e-05, 6.77945e-05, -0.000104923, 6.87795e-05, -0.000117168, 7.67749e-05, -0.000130278, 8.43921e-05, -0.000137269, 8.63373e-05, -0.000142283, 8.78225e-05, -0.000145142, 9.22318e-05, -0.000121506, 0.000104699, -0.000114551, 0.00010633, -0.0001077, 9.84409e-05, -0.000109789, 8.14774e-05, -8.70686e-05};

  double factor_[8] = {-0.07, 0.03, -0.07, 0.20, 0.15, -0.05, 0.00, -0.07}; // Correction factors for the phi bias within a sector of 8 towers (phi direction)

  bool ComputeIncidenceSD(float E, float x, float y,
                            float& cos_a_phi, float& cos_a_eta,
                            float& a_phi_sgn, float& a_eta_sgn);
};



#endif
