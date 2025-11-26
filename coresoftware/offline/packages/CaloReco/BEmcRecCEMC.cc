#include "BEmcRecCEMC.h"

#include "BEmcCluster.h"
#include "BEmcProfile.h"
#include "TVector3.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <limits>   // for quiet_NaN()
#include <fstream>  // for reading bFit_master.txt
#include <map>
#include <tuple>
#include <algorithm>
#include <cstdlib>  // getenv
#include <string>

namespace {
// Translation-unit local verbosity flag. Flip to true to see logs everywhere.
static bool g_verbose_bfit = false;

// Convenience: behaves like std::cerr only when verbose is on.
#define VLOG() if (!::g_verbose_bfit) {} else std::cerr
} // namespace

namespace {  // ======================= VERBOSE b-value infra =======================

// Key for lookup: (variantName, etaOrPhi, etaRange, zRange)
struct BKey {
  std::string variant, eop, etaRange, zRange;
  bool operator<(const BKey& o) const {
    return std::tie(variant, eop, etaRange, zRange)
         < std::tie(o.variant, o.eop, o.etaRange, o.zRange);
  }
};

struct Fit {
  double m  { std::numeric_limits<double>::quiet_NaN() };
  double b0 { std::numeric_limits<double>::quiet_NaN() };
};

static inline std::string key_str(const BKey& k) {
  return "{variant=" + k.variant + ", eop=" + k.eop +
         ", eta=" + k.etaRange + ", z=" + k.zRange + "}";
}

class BFitDB {
 public:
  static BFitDB& instance() { static BFitDB db; return db; }

  // Returns nullptr if not found
  const Fit* get(const std::string& variant,
                 const std::string& eop,
                 const std::string& etaRange,
                 const std::string& zRange) const {
    BKey k{variant, eop, etaRange, zRange};
    auto it = map_.find(k);
    VLOG() << "[BFitDB::get] lookup " << key_str(k) << " -> "
           << (it == map_.end() ? "MISS" :
               ("HIT(m=" + std::to_string(it->second.m) +
                ", b0=" + std::to_string(it->second.b0) + ")"))
           << "\n";
    return (it == map_.end()) ? nullptr : &it->second;
  }

 private:
  std::map<BKey, Fit> map_;

  BFitDB() { load(); }

  void load() {
    VLOG() << "[BFitDB::load] begin\n";
    const char* env = std::getenv("BEMC_BFIT_MASTER");
    std::string path = env ? std::string(env) : std::string("/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloReco/bFit_master.txt");
    VLOG() << "[BFitDB::load] env BEMC_BFIT_MASTER=" << (env ? env : "<unset>") << "\n";
    VLOG() << "[BFitDB::load] opening file: " << path << "\n";

    std::ifstream fin(path);
    if (!fin) {
      VLOG() << "[BFitDB::load] open FAILED: " << path << "\n";
      std::cerr << "[BEmcRecCEMC] WARNING: cannot open bFit master file: " << path << "\n";
      return;
    }
    VLOG() << "[BFitDB::load] open OK\n";

    // Try to skip an optional header
    {
      std::string header;
      std::streampos pos = fin.tellg();
      if (std::getline(fin, header)) {
        if (header.find("variantName") == std::string::npos) {
          fin.seekg(pos);
          VLOG() << "[BFitDB::load] no header line detected (rewind)\n";
        } else {
          VLOG() << "[BFitDB::load] header: " << header << "\n";
        }
      } else {
        VLOG() << "[BFitDB::load] file empty; aborting\n";
        return;
      }
    }

    std::string variant, eop, etaRange, zRange;
    double m = std::numeric_limits<double>::quiet_NaN();
    double b0 = std::numeric_limits<double>::quiet_NaN();

    size_t row = 0, inserted = 0, overwritten = 0;
    while (fin >> variant >> eop >> etaRange >> zRange >> m >> b0) {
      ++row;
      BKey k{variant, eop, etaRange, zRange};
      auto pre = map_.find(k);
      if (pre != map_.end()) {
        ++overwritten;
        VLOG() << "[BFitDB::load] duplicate key at data row " << row
               << " " << key_str(k) << " -> overwriting old {m="
               << pre->second.m << ", b0=" << pre->second.b0 << "} with {m="
               << m << ", b0=" << b0 << "}\n";
      } else {
        VLOG() << "[BFitDB::load] row " << row << " insert "
               << key_str(k) << " m=" << m << " b0=" << b0 << "\n";
      }
      map_[k] = Fit{m, b0};
      ++inserted;
    }

    VLOG() << "[BFitDB::load] done: rows_read=" << inserted
           << ", unique_keys=" << map_.size()
           << ", overwritten=" << overwritten << "\n";
  }
};

// No clamping — pass-through (leave only a non-finite guard)
inline float clamp_b(float b) {
  if (!std::isfinite(b)) {
    VLOG() << "[clamp_b] non-finite b -> fallback 0.15\n";
    return 0.15f;         // your existing fallback policy
  }
  return b;               // return the fitted/derived value as-is
}

// Build b(E) from fit (b0 + m ln(E/E0)); fallback to 0.15 if missing
inline float b_from(const Fit* fp, float Energy, float E0 = 3.0f) {
  if (!fp) {
    VLOG() << "[b_from] NO FIT available → fallback b=0.15 (E=" << Energy
           << ", E0=" << E0 << ")\n";
    return 0.15f;
  }
  const float lnE = std::log(std::max(Energy, 1e-6f) / E0);
  const float raw = static_cast<float>(fp->b0 + fp->m * lnE);
  const float clamped = clamp_b(raw);
  VLOG() << "[b_from] using {m=" << fp->m << ", b0=" << fp->b0
         << "}, E=" << Energy << ", ln(E/E0)=" << lnE
         << " → b_raw=" << raw << " → b=" << clamped << "\n";
  return clamped;
}

// Map |η| to label
inline std::string etaLabelFromAbsEta(float absEta) {
  std::string lab;
  if      (absEta <= 0.20f) lab = "etaCore";
  else if (absEta <= 0.70f) lab = "etaMid";
  else if (absEta <  1.10f) lab = "etaEdge";
  else                      lab = "fullEta";
  VLOG() << "[etaLabelFromAbsEta] |eta|=" << absEta << " → " << lab << "\n";
  return lab;
}

// Map |z_vtx| to the preferred *fine* label for 0–10 and *coarse* label beyond
inline std::string zLabelFromAbsZ(float z) {
  std::string lab;
  if (z < 10.f) {
    if      (z < 2.f) lab = "z00to02";
    else if (z < 4.f) lab = "z02to04";
    else if (z < 6.f) lab = "z04to06";
    else if (z < 8.f) lab = "z06to08";
    else              lab = "z08to10";
  } else if (z < 20.f) lab = "z10to20";
    else if (z < 30.f) lab = "z20to30";
    else if (z < 45.f) lab = "z30to45";
    else               lab = "z45to60"; // clamp ≥60 to last bin
  VLOG() << "[zLabelFromAbsZ] |z_vtx|=" << z << " → " << lab << "\n";
  return lab;
}

// ======================= δ(E,α) (Delta) loader ======================
struct DeltaFit {
  double c0{std::numeric_limits<double>::quiet_NaN()};
  double m {std::numeric_limits<double>::quiet_NaN()};
  double E0{3.0};
  bool    valid{false};
};

class DeltaDB {
public:
  static DeltaDB& instance() { static DeltaDB d; return d; }

  // tag must be "phi" or "eta"; returns nullptr if not loaded
  const DeltaFit* get(const char* tag) const {
    if (std::string(tag)=="phi" && phi_.valid) return &phi_;
    if (std::string(tag)=="eta" && eta_.valid) return &eta_;
    return nullptr;
  }

private:
  DeltaFit phi_{}, eta_{};

  static bool LoadOne(const std::string& path, const char* expectTag, DeltaFit& out) {
    std::ifstream fin(path);
    if (!fin) return false;
    std::string hdr, tag, keyE0;
    double c0=std::numeric_limits<double>::quiet_NaN();
    double m =std::numeric_limits<double>::quiet_NaN();
    double E0=3.0;
    // expected line format:
    //   deltaFit  <tag>  <c0>  <m>  E0 <E0>
    if (!(fin >> hdr >> tag >> c0 >> m >> keyE0 >> E0)) return false;
    if (hdr != "deltaFit") return false;
    if (std::string(expectTag)!=tag)   return false;
    out.c0 = c0; out.m = m; out.E0 = E0; out.valid = (std::isfinite(c0) && std::isfinite(m) && std::isfinite(E0));
    return out.valid;
  }

  // construct and load both files; look for env overrides first
  DeltaDB() {
    // Derive default directory from BFitDB env if present
    std::string baseDir;
    if (const char* bfit = std::getenv("BEMC_BFIT_MASTER")) {
      std::string p(bfit);
      auto pos = p.find_last_of("/\\");
      baseDir = (pos==std::string::npos) ? std::string(".") : p.substr(0, pos);
    } else {
      // fall back to your existing default directory; adjust if needed
      baseDir = "/sphenix/u/patsfan753/scratch/PDCru n24pp/coresoftware/offline/packages/CaloReco";
    }

    std::string phiPath = (std::getenv("BEMC_DELTA_PHI") ?
                           std::string(std::getenv("BEMC_DELTA_PHI")) :
                           (baseDir + "/deltaFit_master_phi.txt"));
    std::string etaPath = (std::getenv("BEMC_DELTA_ETA") ?
                           std::string(std::getenv("BEMC_DELTA_ETA")) :
                           (baseDir + "/deltaFit_master_eta.txt"));

    bool okP = LoadOne(phiPath, "phi", phi_);
    bool okE = LoadOne(etaPath, "eta", eta_);

      VLOG() << "[DeltaDB] load phi(" << (okP ? "OK " : "MISS ") << ") from " << phiPath << "\n";
      VLOG() << "[DeltaDB] load  eta(" << (okE ? "OK " : "MISS ") << ") from " << etaPath << "\n";
      if (!okP) { VLOG() << "[DeltaDB] WARNING: φ delta fit missing, proceeding with δφ=0\n"; }
      if (!okE) { VLOG() << "[DeltaDB] WARNING: η delta fit missing, proceeding with δη=0\n"; }

  }
};

// convenience: evaluate δ_a(E,α) = (c0 + m ln(E/E0)) * sin α
inline double eval_delta(const DeltaFit* df, double E, double alpha) {
  if (!df || !df->valid) return 0.0;
  const double logR = std::log(std::max(1e-6, E / df->E0));
  const double c1   = df->c0 + df->m * logR;
  return c1 * std::sin(alpha);
}

} // namespace


BEmcRecCEMC::BEmcRecCEMC()
//  : _emcprof(nullptr)
{
  Name("BEmcRecCEMC");
  SetCylindricalGeometry();
}

void BEmcRecCEMC::LoadProfile(const std::string& fname)
{
  //  std::cout << "Infor from BEmcRecCEMC::LoadProfile(): no external file used for shower profile evaluation in CEMC" << std::endl;
  _emcprof = new BEmcProfile(fname);
}

void BEmcRecCEMC::GetImpactThetaPhi(float xg, float yg, float zg, float& theta, float& phi)
{
  theta = 0;
  phi = 0;

  //  float theta = atan(sqrt(xg*xg + yg*yg)/fabs(zg-fVz));
  float rg = std::sqrt((xg * xg) + (yg * yg));
  float theta_twr;
  if (std::fabs(zg) <= 15)
  {
    theta_twr = 0;
  }
  else if (zg > 15)
  {
    theta_twr = std::atan2(zg - 15, rg);
  }
  else
  {
    theta_twr = std::atan2(zg + 15, rg);
  }
  float theta_tr = std::atan2(zg - fVz, rg);
  theta = std::fabs(theta_tr - theta_twr);
  //  phi = atan2(yg,xg);
}

/*
float BEmcRecCEMC::GetProb(vector<EmcModule> HitList, float ecl, float xg, float yg, float zg, float& chi2, int& ndf)
{
  chi2 = 0;
  ndf = 0;
  float prob = -1;

  //  float theta = atan(sqrt(xg*xg + yg*yg)/fabs(zg-fVz));
  float rg = sqrt(xg * xg + yg * yg);
  float theta_twr;
  if (fabs(zg) <= 15)
    theta_twr = 0;
  else if (zg > 15)
    theta_twr = atan2(zg - 15, rg);
  else
    theta_twr = atan2(zg + 15, rg);
  float theta_tr = atan2(zg - fVz, rg);
  float theta = fabs(theta_tr - theta_twr);

  float phi = atan2(yg, xg);
  if (_emcprof != nullptr) prob = _emcprof->GetProb(&HitList, fNx, ecl, theta, phi);

  return prob;
}
*/
/*
float BEmcRecCEMC::GetProb(vector<EmcModule> HitList, float et, float xg, float yg, float zg, float& chi2, int& ndf)
// et, xg, yg, zg not used here
{
  const float thresh = 0.01;
  const int DXY = 3;  // 2 is for 5x5 matrix; 3 for 7x7 matrix
  const int Nmax = 1000;
  float ee[Nmax];
  int iyy[Nmax];
  int izz[Nmax];

  int ich;
  vector<EmcModule>::iterator ph = HitList.begin();

  chi2 = 0;
  ndf = 0;

  int nn = 0;

  while (ph != HitList.end())
  {
    ee[nn] = ph->amp;
    if (ee[nn] > thresh)
    {
      ich = ph->ich;
      izz[nn] = ich % fNx;
      iyy[nn] = ich / fNx;
      nn++;
      if (nn >= Nmax)
      {
      std::cout << "BEmcRec::GetProb: Cluster size is too big. Skipping the rest of the towers" << std::endl;
        break;
      }
    }  // if( ee[nn]
    ++ph;
  }  // while( ph

  if (nn <= 0) return -1;

  int iy0 = -1, iz0 = -1;
  float emax = 0;

  for (int i = 0; i < nn; i++)
  {
    if (ee[i] > emax)
    {
      emax = ee[i];
      iy0 = iyy[i];
      iz0 = izz[i];
    }
  }

  if (emax <= 0) return -1;

  int id;
  float etot = 0;
  float sz = 0;
  float sy = 0;

  for (int idz = -DXY; idz <= DXY; idz++)
  {
    for (int idy = -DXY; idy <= DXY; idy++)
    {
      id = GetTowerID(iy0 + idy, iz0 + idz, nn, iyy, izz, ee);
      if (id >= 0)
      {
        etot += ee[id];
        sz += ee[id] * (iz0 + idz);
        sy += ee[id] * (iy0 + idy);
      }
    }
  }
  float zcg = sz / etot;  // Here cg allowed to be out of range
  float ycg = sy / etot;
  int iz0cg = int(zcg + 0.5);
  int iy0cg = int(ycg + 0.5);
  float ddz = fabs(zcg - iz0cg);
  float ddy = fabs(ycg - iy0cg);

  int isz = 1;
  if (zcg - iz0cg < 0) isz = -1;
  int isy = 1;
  if (ycg - iy0cg < 0) isy = -1;

  // 4 central towers: 43
  //                   12
  // Tower 1 - central one
  float e1, e2, e3, e4;
  e1 = e2 = e3 = e4 = 0;
  id = GetTowerID(iy0cg, iz0cg, nn, iyy, izz, ee);
  if (id >= 0) e1 = ee[id];
  id = GetTowerID(iy0cg, iz0cg + isz, nn, iyy, izz, ee);
  if (id >= 0) e2 = ee[id];
  id = GetTowerID(iy0cg + isy, iz0cg + isz, nn, iyy, izz, ee);
  if (id >= 0) e3 = ee[id];
  id = GetTowerID(iy0cg + isy, iz0cg, nn, iyy, izz, ee);
  if (id >= 0) e4 = ee[id];

  float e1t = (e1 + e2 + e3 + e4) / etot;
  float e2t = (e1 + e2 - e3 - e4) / etot;
  float e3t = (e1 - e2 - e3 + e4) / etot;
  float e4t = (e3) / etot;
  //  float e5t = (e2+e4)/etot;

  float rr = sqrt((0.5 - ddz) * (0.5 - ddz) + (0.5 - ddy) * (0.5 - ddy));

  float c1, c2, c11;

  float logE = log(etot);

  // e1 energy is the most effective for PID if properly tuned !
  // Discrimination power is very sensitive to paramter c1: the bigger it is
  // the better discrimination;
  c1 = 0.95;
  c2 = 0.0066364 * logE + 0.00466667;
  if (c2 < 0) c2 = 0;
  float e1p = c1 - c2 * rr * rr;
  c1 = 0.034 - 0.01523 * logE + 0.0029 * logE * logE;
  float err1 = c1;

  // For e2
  c1 = 0.00844086 + 0.00645359 * logE - 0.00119381 * logE * logE;
  if (etot > 15) c1 = 0.00844086 + 0.00645359 * log(15.) - 0.00119381 * log(15.) * log(15.);  // Const at etot>15GeV
  if (c1 < 0) c1 = 0;
  c2 = 3.9;                                                      // Fixed
  float e2p = sqrt(c1 + 0.25 * c2) - sqrt(c1 + c2 * ddy * ddy);  // =0 at ddy=0.5

  c1 = 0.0212333 + 0.0420473 / etot;
  c2 = 0.090;  // Fixed
  float err2 = c1 + c2 * ddy;
  if (ddy > 0.3) err2 = c1 + c2 * 0.3;  // Const at ddy>0.3

  // For e3
  c1 = 0.0107857 + 0.0056801 * logE - 0.000892016 * logE * logE;
  if (etot > 15) c1 = 0.0107857 + 0.0056801 * log(15.) - 0.000892016 * log(15.) * log(15.);  // Const at etot>15GeV
  if (c1 < 0) c1 = 0;
  c2 = 3.9;                                                      // Fixed
  float e3p = sqrt(c1 + 0.25 * c2) - sqrt(c1 + c2 * ddz * ddz);  // =0 at ddz=0.5

  //  c1 = 0.0200 + 0.042/etot;
  c1 = 0.0167 + 0.058 / etot;
  c2 = 0.090;  // Fixed
  float err3 = c1 + c2 * ddz;
  if (ddz > 0.3) err3 = c1 + c2 * 0.3;  // Const at ddz>0.3

  // For e4
  float e4p = 0.25 - 0.668 * rr + 0.460 * rr * rr;
  c11 = 0.171958 + 0.0142421 * logE - 0.00214827 * logE * logE;
  //  c11 = 0.171085 + 0.0156215*logE - -0.0025809*logE*logE;
  float err4 = 0.102 - 1.43 * c11 * rr + c11 * rr * rr;  // Min is set to x=1.43/2.
  err4 *= 1.1;

  chi2 = 0.;
  chi2 += (e1p - e1t) * (e1p - e1t) / err1 / err1;
  chi2 += (e2p - e2t) * (e2p - e2t) / err2 / err2;
  chi2 += (e3p - e3t) * (e3p - e3t) / err3 / err3;
  chi2 += (e4p - e4t) * (e4p - e4t) / err4 / err4;
  ndf = 4;

  //  chi2 /= 1.1;
  float prob = TMath::Prob(chi2, ndf);

  return prob;
}
*/


void BEmcRecCEMC::SetPhiTiltVariant(ETiltVariant v)
{
  // (a,b) in radians for φ(E) = a − b·lnE.  If you later calibrate EA flavours,
  // update the four CLUS_CP_EA_* entries below.
  static const std::unordered_map<ETiltVariant,std::pair<double,double>> lut = {
      {ETiltVariant::CLUS_RAW, { 2.450797e-03, 8.010320e-04 }},  // no corr, cluster
      {ETiltVariant::CLUS_CP, { 1.738836e-03, 7.853763e-04 }},  // CorrectPosition, cluster
      {ETiltVariant::CLUS_CP_EA_FIT_ZDEP_ETADEP, { 1.886333e-03, 8.111119e-04 }},  // CorrectPosition(EA |z|+|#eta|+E), cluster
      {ETiltVariant::CLUS_CP_EA_FIT_ZDEP, { 1.880064e-03, 8.097306e-04 }},  // CorrectPosition(EA |z|+E), cluster
      {ETiltVariant::CLUS_CP_EA_FIT_ETADEP, { 1.886620e-03, 8.119962e-04 }},  // CorrectPosition(EA |#eta|+E), cluster
      {ETiltVariant::CLUS_CP_EA_FIT_EandIncident, { 1.956265e-03, 8.082764e-04 }},  // CorrectPosition(EA E-only + incident-angle), cluster
      {ETiltVariant::CLUS_CP_EA_FIT_EONLY, { 1.879850e-03, 8.095643e-04 }},  // CorrectPosition(EA E-only), cluster
      {ETiltVariant::CLUS_BCORR, { 1.924901e-03, 8.042676e-04 }},  // b(E) corr, cluster
      {ETiltVariant::PDC_RAW, { 2.542701e-03, 8.062445e-04 }},  // no corr, scratch
      {ETiltVariant::PDC_CORR, { 2.085638e-03, 8.286254e-04 }},  // b(E) corr, scratch

      {ETiltVariant::DEFAULT,               { 8.654924e-04, 8.490399e-04}}   // legacy numbers
  };

  auto it = lut.find(v);
  m_abTiltCurrent = (it != lut.end()) ? it->second : lut.at(ETiltVariant::DEFAULT);
  m_tiltVariant   = v;
}





void BEmcRecCEMC::CorrectShowerDepth(int ix, int iy, float E, float xA, float yA, float zA, float& xC, float& yC, float& zC)
{
  xC = xA;
  yC = yA;
  zC = zA;

  float logE = log(0.1);
  if (E > 0.1)
  {
    logE = std::log(E);
  }

  /* -------------------------------------------------------------- *
    *  Azimuthal-tilt correction  (2025-06-25 calibration)
    *  Internal toggle: set doPhiTilt=false to skip rotation
    *
    *        φ(E) =  a  –  b · ln E      with  E in GeV
    *
    * -------------------------------------------------------------- */
    bool doPhiTilt = false; // <<< set to false internally to disable
    if (doPhiTilt)
    {
        /* variant–dependent tilt constants (set via SetPhiTiltVariant) */
        const double aTilt = m_abTiltCurrent.first;   // rad
        const double bTilt = m_abTiltCurrent.second;  // rad

        const double phi = aTilt - bTilt * static_cast<double>(logE);  // use logE already defined
        const double  c  = std::cos(phi);
        const double  s  = std::sin(phi);

        /* rotate the uncorrected impact point (xA, yA) */
        xC = xA * c - yA * s;
        yC = xA * s + yA * c;
    }
    else
    {
        // skip tilt correction entirely
        xC = xA;
        yC = yA;
  }
  /* -------------------------------------------------------------- */

  // Correction in z
  float rA = std::sqrt((xA * xA) + (yA * yA));
  float theta_twr;
  if (m_UseDetailedGeometry)
  {
    // Read the angle right from the detailed RawTowerGeom objects
    TowerGeom geom0;
    GetTowerGeometry(ix, iy, geom0);
    theta_twr = M_PI / 2 + geom0.rotX;
  }
  else
  {
    // Use the approximate default tower angles array.
    theta_twr = M_PI / 2 - angles[iy];
  }

  float theta_tr = std::atan2(zA - fVz, rA);

  // Shower CG in long. direction
  // Different tuning for the approximate and detailed geometry
  float L = 0;
  if (m_UseDetailedGeometry)
  {
    L = -2.67787 + 0.924138 * logE;
  }
  else
  {
    L = -1.79968 + 0.837322 * logE;
  }
    float dz = L * std::sin(theta_tr - theta_twr) / std::cos(theta_twr);


  if (!m_UseDetailedGeometry)
  {
    // Not strictly speaking a "shower depth correction" but rather a baseline correction
    // The factor 0.10 accounts for the fact that the approximate geometry
    // is projected at 93.5cm, which is roughly 90 % of the actual average EMCal radius (~ 105 cm)
    dz -= fVz * 0.10;
  }

  zC = zA - dz;
//    std::cout << "[CSD] depth shift: zA=" << zA << "  ->  zC=" << zC << "\n";
//    if (!std::isfinite(phi) || !std::isfinite(xC) || !std::isfinite(yC) || !std::isfinite(zC)) {
//      std::cout << "[CSD][WARN] non-finite output (phi/xC/yC/zC)\n";
//    }
    
  // zvtx-dependent pseudorapidity-correction
  // At large zvtx (> 20 cm), the sawtooth shape of the EMCal is no longer projective wrt collision point,
  // it introduces a bias which can be approximately corrected with the linear formula below
  if (m_UseDetailedGeometry && std::abs(fVz) > 20)
  {
    float eta = std::asinh((zC - fVz) / rA);
    if (fVz > 0)
    {
      eta -= slopes_[iy] * (fVz - 20);
    }
    else
    {
      eta -= slopes_[fNy - iy - 1] * (fVz + 20);
    }
    zC = fVz + rA * std::sinh(eta);
  }

  return;
}

void BEmcRecCEMC::CorrectEnergy(float Energy, float /*x*/, float /*y*/,
                                float& Ecorr)
{
  // Corrects the EM Shower Energy for attenuation in fibers and
  // long energy leakage
  //
  // (x,y) - impact position (cm) in Sector frame
  /*
  float sinT;
  float att, leak, corr;
  const float leakPar = 0.0033; // parameter from fit
  const float attPar = 120; // Attenuation in module (cm)
  const float X0 = 2; // radiation length (cm)

  *Ecorr = Energy;
  if( Energy < 0.01 ) return;

  GetImpactAngle(x, y, &sinT); // sinT not used so far
  leak = 2-sqrt(1+leakPar*log(1+Energy)*log(1+Energy));
  att = exp(log(Energy)*X0/attPar);
  corr = leak*att;
  *Ecorr = Energy/corr;
  */
  Ecorr = Energy;
}

void BEmcRecCEMC::CorrectECore(float Ecore, float /*x*/, float /*y*/, float& Ecorr)
{
  // Corrects the EM Shower Core Energy for attenuation in fibers,
  // long energy leakage and angle dependance
  //
  // (x,y) - impact position (cm) in Sector frame

  const float c0 = 0.950;  // For no threshold
  Ecorr = Ecore / c0;
}

void BEmcRecCEMC::CorrectPosition(float Energy, float x, float y,
                                  float& xc, float& yc)
{
  // Corrects the Shower Center of Gravity for the systematic shift due to
  // limited tower size
  //
  // Everything here is in tower units.
  // (x,y) - CG position, (xc,yc) - corrected position

  float bx;
  float by;
  float x0;
  float y0;
  int ix0;
  int iy0;

  if (!m_UseCorrectPosition)
  {
    xc = x;
    yc = y;
    return;
  }

  if (Energy < 0.01)
  {
    return;
  }

  bx = 0.15;
  by = 0.15;

  x0 = x;
  ix0 = EmcCluster::lowint(x0 + 0.5);

  if (std::abs(x0 - ix0) <= 0.5)
  {
    x0 = ix0 + bx * asinh(2. * (x0 - ix0) * sinh(0.5 / bx));
  }
  else
  {
    x0 = x;
    std::cout << "????? Something wrong in BEmcRecCEMC::CorrectPosition: x = "
              << x << " dx = " << x0 - ix0 << std::endl;
  }

  // Correct for phi bias within module of 8 towers
// NOLINTNEXTLINE(bugprone-incorrect-roundings)
  int ix8 = int(x + 0.5) / 8; // that is hokey - suggest lroundf(x)
  float x8 = x + 0.5 - (ix8 * 8) - 4;  // from -4 to +4
  float dx = 0;
  if (m_UseDetailedGeometry)
  {
    // Don't know why there is a different factor for each tower of the sector
    // Just tuned from MC
// NOLINTNEXTLINE(bugprone-incorrect-roundings)
    int local_ix8 = int(x+0.5) - ix8 * 8; // that is hokey - suggest lroundf(x)
    dx = factor_[local_ix8] * x8 / 4.;
  }
  else
  {
    dx = 0.10 * x8 / 4.;
    if (std::fabs(x8) > 3.3)
    {
      dx = 0;  // Don't correct near the module edge
    }
  }

  xc = x0 - dx;
  while (xc < -0.5)
  {
    xc += float(fNx);
  }
  while (xc >= fNx - 0.5)
  {
    xc -= float(fNx);
  }

  y0 = y;
  iy0 = EmcCluster::lowint(y0 + 0.5);

  if (std::abs(y0 - iy0) <= 0.5)
  {
    y0 = iy0 + by * std::asinh(2. * (y0 - iy0) * sinh(0.5 / by));
  }
  else
  {
    y0 = y;
    std::cout << "????? Something wrong in BEmcRecCEMC::CorrectPosition: y = "
              << y << "dy = " << y0 - iy0 << std::endl;
  }
  yc = y0;
}



void BEmcRecCEMC::CorrectPositionEnergyAwareZVTXEtaAndEnergyDep(float Energy, float vtxZ_cm,
                                                                float x, float y,
                                                                float& xc, float& yc)
{
  // Legacy gating
  if (!m_UseCorrectPosition) { xc = x; yc = y; return; }
  if (!std::isfinite(Energy) || !std::isfinite(x) || !std::isfinite(y) || Energy < 0.01f)
  { xc = x; yc = y; return; }

  // Determine |η| slice from tower geometry
  const int ix0 = EmcCluster::lowint(x + 0.5f);
  const int iy0 = EmcCluster::lowint(y + 0.5f);

  TowerGeom g{};
  bool haveGeom = GetTowerGeometry(ix0, iy0, g);

  float absEta = 0.f;
  if (haveGeom) {
    const float r = std::hypot(g.Xcenter, g.Ycenter);
    absEta = std::fabs(std::asinh((r > 0.f ? g.Zcenter / r : 0.f)));
  }
  std::string etaLbl = haveGeom ? etaLabelFromAbsEta(absEta) : std::string("fullEta");

  // Map |z_vtx| to label (fine for 0–10, coarse otherwise; clamp ≥60 to 45–60)
  const float absZ = std::fabs(vtxZ_cm);
  const std::string zLbl = zLabelFromAbsZ(absZ);

  // Fetch per-(η,z) fits: variant zAndEtaAndEnergyDep
  constexpr float E0 = 3.0f;
  const auto* fitPhi = BFitDB::instance().get("zAndEtaAndEnergyDep", "phi", etaLbl, zLbl);
  const auto* fitEta = BFitDB::instance().get("zAndEtaAndEnergyDep", "eta", etaLbl, zLbl);

  const float bx = b_from(fitPhi, Energy, E0);
  const float by = b_from(fitEta, Energy, E0);

  // ----- identical inverse-asinh + ripple/wrap as elsewhere -----
  float x_corr = x;
  {
    if (std::fabs(x - ix0) <= 0.5f) {
      const float Sx = std::sinh(0.5f / bx);
      x_corr = ix0 + bx * std::asinh(2.f * (x - ix0) * Sx);
    }

    // module-of-8 ripple
    int   ix8 = int(x + 0.5f) / 8;                      // NOLINT(bugprone-incorrect-roundings)
    float x8  = x + 0.5f - (ix8 * 8) - 4.f;             // −4 … +4
    float dx  = 0.f;
    if (m_UseDetailedGeometry) {
      int local_ix8 = int(x + 0.5f) - ix8 * 8;          // NOLINT(bugprone-incorrect-roundings)
      dx = static_cast<float>(factor_[local_ix8]) * (x8 / 4.f);
    } else {
      dx = (std::fabs(x8) > 3.3f) ? 0.f : 0.10f * (x8 / 4.f);
    }
    x_corr -= dx;

    while (x_corr < -0.5f)              x_corr += float(fNx);
    while (x_corr >= float(fNx) - 0.5f) x_corr -= float(fNx);
  }
  xc = x_corr;

  float y_corr = y;
  {
    if (std::fabs(y - iy0) <= 0.5f) {
      const float Sy = std::sinh(0.5f / by);
      y_corr = iy0 + by * std::asinh(2.f * (y - iy0) * Sy);
    }
  }
  yc = y_corr;
}



/*
 energy and z dependent -- z < 60 tight cut in CODE
 */

void BEmcRecCEMC::CorrectPositionEnergyAwareZVTXAndEnergyDep(float Energy, float vtxZ_cm,
                                                             float x, float y,
                                                             float& xc, float& yc)

{
  // Legacy gating (match your other correctors)
  if (!m_UseCorrectPosition) { xc = x; yc = y; return; }
  if (!std::isfinite(Energy) || !std::isfinite(x) || !std::isfinite(y) || Energy < 0.01f)
  { xc = x; yc = y; return; }

  // Map |z_vtx| to a z-slice index
  const float absZ = std::fabs(vtxZ_cm);


  // Use per-|z| fits from the master file (zAndEnergyDep, originalEta)
  constexpr float E0 = 3.0f;

  const std::string zLbl = zLabelFromAbsZ(absZ);
  const auto* fitPhi = BFitDB::instance().get("zAndEnergyDep", "phi", "originalEta", zLbl);
  const auto* fitEta = BFitDB::instance().get("zAndEnergyDep", "eta", "originalEta", zLbl);

  const float bx = b_from(fitPhi, Energy, E0);  // φ-direction b(E; z)
  const float by = b_from(fitEta, Energy, E0);  // η-direction b(E; z)

  // ----- identical inverse-asinh + ripple/wrap as your other correctors -----
  const int ix0 = EmcCluster::lowint(x + 0.5f);
  const int iy0 = EmcCluster::lowint(y + 0.5f);

  // φ
  float x_corr = x;
  if (std::fabs(x - ix0) <= 0.5f)
  {
    const float Sx = std::sinh(0.5f / bx);
    x_corr = ix0 + bx * std::asinh(2.f * (x - ix0) * Sx);
  }

  // module-of-8 ripple (same as legacy CorrectPosition)
  int   ix8 = int(x + 0.5f) / 8;                       // NOLINT(bugprone-incorrect-roundings)
  float x8  = x + 0.5f - (ix8 * 8) - 4.f;               // −4 … +4
  float dx  = 0.f;
  if (m_UseDetailedGeometry)
  {
    int local_ix8 = int(x + 0.5f) - ix8 * 8;            // 0..7   // NOLINT(bugprone-incorrect-roundings)
    dx = static_cast<float>(factor_[local_ix8]) * (x8 / 4.f);
  }
  else
  {
    dx = (std::fabs(x8) > 3.3f) ? 0.f : 0.10f * (x8 / 4.f);
  }
  x_corr -= dx;

  // Wrap φ tower coordinate to [−0.5, Nx−0.5)
  while (x_corr < -0.5f)              x_corr += float(fNx);
  while (x_corr >= float(fNx) - 0.5f) x_corr -= float(fNx);

  // η
  float y_corr = y;
  if (std::fabs(y - iy0) <= 0.5f)
  {
    const float Sy = std::sinh(0.5f / by);
    y_corr = iy0 + by * std::asinh(2.f * (y - iy0) * Sy);
  }

  // Output
  xc = x_corr;
  yc = y_corr;
}


/*
 energy and |eta| dependent -- z < 60 tight cut in CODE
 */

void BEmcRecCEMC::CorrectPositionEnergyAwareEtaAndEnergyDep(float Energy, float x, float y,
                                             float& xc, float& yc)
{
  // Legacy gating
  if (!m_UseCorrectPosition) { xc = x; yc = y; return; }
  if (!std::isfinite(Energy) || !std::isfinite(x) || !std::isfinite(y) || Energy < 0.01f)
  { xc = x; yc = y; return; }

  // ---- determine |η| slice from tower geometry (fallback: no-η-dep) ----
  const int ix0 = EmcCluster::lowint(x + 0.5f);
  const int iy0 = EmcCluster::lowint(y + 0.5f);

  TowerGeom g{};
  bool haveGeom = GetTowerGeometry(ix0, iy0, g);

  float absEta = 0.f;
  if (haveGeom)
  {
    const float r = std::hypot(g.Xcenter, g.Ycenter);
    absEta = std::fabs(std::asinh((r > 0.f ? g.Zcenter / r : 0.f)));
  }

  // |η| bins: 0: [0,0.20], 1: (0.20,0.70], 2: (0.70,1.10], 3: no-η-dep (fallback)
  int iEtaBin = 3;
  if (haveGeom)
  {
    if      (absEta <= 0.20f) iEtaBin = 0;
    else if (absEta <= 0.70f) iEtaBin = 1;
    else                      iEtaBin = 2; // up to ~1.10
  }

    // Use fits per η-range from the master file (etaAndEnergyDep, originalZRange)
    constexpr float E0 = 3.0f;

    std::string etaRangeLbl = "fullEta";
    if      (iEtaBin == 0) etaRangeLbl = "etaCore";
    else if (iEtaBin == 1) etaRangeLbl = "etaMid";
    else if (iEtaBin == 2) etaRangeLbl = "etaEdge";

    const auto* fitPhi = BFitDB::instance().get("etaAndEnergyDep", "phi", etaRangeLbl, "originalZRange");
    const auto* fitEta = BFitDB::instance().get("etaAndEnergyDep", "eta", etaRangeLbl, "originalZRange");

    const float bx = b_from(fitPhi, Energy, E0);  // φ-direction b(E)
    const float by = b_from(fitEta, Energy, E0);  // η-direction b(E)

  // ---------------- φ inverse-asinh (legacy form) ----------------
  float x_corr = x;
  if (std::fabs(x - ix0) <= 0.5f)
  {
    const float Sx = std::sinh(0.5f / bx);
    x_corr = ix0 + bx * std::asinh(2.f * (x - ix0) * Sx);
  }

  // module-of-8 ripple (identical to legacy CorrectPosition)
  // NOLINTNEXTLINE(bugprone-incorrect-roundings)
  int   ix8 = int(x + 0.5f) / 8;
  float x8  = x + 0.5f - (ix8 * 8) - 4.f;   // −4 … +4
  float dx  = 0.f;
  if (m_UseDetailedGeometry)
  {
    // NOLINTNEXTLINE(bugprone-incorrect-roundings)
    int local_ix8 = int(x + 0.5f) - ix8 * 8;
    dx = static_cast<float>(factor_[local_ix8]) * (x8 / 4.f);
  }
  else
  {
    dx = (std::fabs(x8) > 3.3f) ? 0.f : 0.10f * (x8 / 4.f);
  }
  x_corr -= dx;

  // wrap φ tower coordinate to [-0.5, Nx-0.5)
  while (x_corr < -0.5f)              x_corr += float(fNx);
  while (x_corr >= float(fNx) - 0.5f) x_corr -= float(fNx);
  xc = x_corr;

  // ---------------- η inverse-asinh (legacy form) ----------------
  float y_corr = y;
  if (std::fabs(y - iy0) <= 0.5f)
  {
    const float Sy = std::sinh(0.5f / by);
    y_corr = iy0 + by * std::asinh(2.f * (y - iy0) * Sy);
  }
  yc = y_corr;
}













bool BEmcRecCEMC::ComputeIncidenceSD(float E, float x, float y,
                                     float& cos_a_phi, float& cos_a_eta,
                                     float& a_phi_sgn, float& a_eta_sgn)
{
  (void)E;  // Front-face incidence does not depend on energy

  const bool dbg = (m_incDbgLevel > 0);

  // ------------------------------------------------------------------
  // Incidence mode (runtime-selectable):
  //
  //   FACE : use detailed front-face tangents only  (legacy behaviour)
  //   MECH : use mechanical rotations rotX/Y/Z only
  //   BOTH : return MECH result (if available) but ALSO compute and
  //          print both FACE and MECH in a compact tabulated line.
  //
  // Mode is set via SetIncidenceMode(...) on this BEmcRecCEMC object.
  // ------------------------------------------------------------------
  const EIncidenceMode kIncidenceMode = m_incMode;

  // "Primary" frame is what defines the returned values.
  const bool primaryIsMech =
      (kIncidenceMode == EIncidenceMode::MECH ||
       kIncidenceMode == EIncidenceMode::BOTH);

  auto dump_and_maybe_stop =
    [&](const char* tag,
        const TVector3& C, const TVector3& F,
        const TVector3& V, const TVector3& un,
        const TVector3& uphi, const TVector3& ueta,
        const TVector3& pF,
        double pn, double pph, double pet,
        double aphi, double aeta,
        double cphi, double ceta)
  {
    if (!dbg) return;
    std::cout << "\n=== Incidence Trace ("<< tag <<") ===\n"
              << "C=("<<C.X()<<","<<C.Y()<<","<<C.Z()<<")  "
              << "F=("<<F.X()<<","<<F.Y()<<","<<F.Z()<<")  "
              << "V=("<<V.X()<<","<<V.Y()<<","<<V.Z()<<")\n"
              << "un =("<<un.X()<<","<<un.Y()<<","<<un.Z()<<")\n"
              << "uphi=("<<uphi.X()<<","<<uphi.Y()<<","<<uphi.Z()<<")  "
              << "ueta=("<<ueta.X()<<","<<ueta.Y()<<","<<ueta.Z()<<")\n"
              << "pF =("<<pF.X()<<","<<pF.Y()<<","<<pF.Z()<<")  |pF|="<<pF.Mag()<<"\n"
              << "pn="<<pn<<"  pφ="<<pph<<"  pη="<<pet<<"\n"
              << "αφ="<<aphi<<"  αη="<<aeta
              << "  cosφ="<<cphi<<"  cosη="<<ceta<<"\n"
              << "noTiltSandbox="<<(m_incNoTiltSandbox?"true":"false")
              << "  primaryIsMech="<<(primaryIsMech?"true":"false")<<"\n"
              << "======================================\n";
    if (m_incHardStop >= 0 && m_incDbgLevel >= m_incHardStop)
    {
      throw std::runtime_error("Incidence debug hard-stop");
    }
  };

  // ────────────────────────────────────────────────────────────────────────────
  // Step 1) Owner tower in "tower units"
  // ────────────────────────────────────────────────────────────────────────────
  const int ix = EmcCluster::lowint(x + 0.5f);
  const int iy = EmcCluster::lowint(y + 0.5f);

  TowerGeom g{};
  if (!GetTowerGeometry(ix, iy, g)) return false; // need detailed geometry

  // Sub-cell offsets δ ∈ (−0.5, +0.5]
  const float dx = x - ix;
  const float dy = y - iy;

  // ────────────────────────────────────────────────────────────────────────────
  // Step 2) FRONT-FACE point F in cm:
  //         F = C + δx * eφ + δy * eη
  //         (dX[], dY[], dZ[] are one-pitch tangents from detailed geometry)
  // ────────────────────────────────────────────────────────────────────────────
  const TVector3 C(g.Xcenter, g.Ycenter, g.Zcenter);
  const TVector3 ephi_face(g.dX[0], g.dY[0], g.dZ[0]);
  const TVector3 eeta_face(g.dX[1], g.dY[1], g.dZ[1]);
  if (ephi_face.Mag() < 1e-9 || eeta_face.Mag() < 1e-9) return false;

  TVector3 ephi = ephi_face;
  TVector3 eeta = eeta_face;

  // Optional no-tilt sandbox: cylinder-aligned φ/η but keep pitch magnitudes
  if (m_incNoTiltSandbox)
  {
    const TVector3 rhat = TVector3(C.X(), C.Y(), 0.0).Unit();
    const TVector3 zhat(0., 0., 1.);
    const TVector3 uphi_cyl = (zhat.Cross(rhat)).Unit();      // φ̂
    const TVector3 ueta_cyl = (rhat.Cross(uphi_cyl)).Unit();  // η̂ = r̂×φ̂
    ephi = uphi_cyl * ephi_face.Mag();
    eeta = ueta_cyl * eeta_face.Mag();
  }

  const TVector3 F = C + dx * ephi + dy * eeta;

  // ────────────────────────────────────────────────────────────────────────────
  // Step 3) Unit ray from vertex to the FRONT FACE
  // ────────────────────────────────────────────────────────────────────────────
  const TVector3 V(0., 0., fVz);
  TVector3 pF = F - V;
  const double pFmag = pF.Mag();
  if (!(pFmag > 0.0)) return false;
  pF *= (1.0 / pFmag);

  // Helper: from a basis {un,uphi,ueta} and ray pF, compute projections,
  // signed incidence angles, and foreshortening cosines.
  auto basis_to_incidence =
    [&](const TVector3& un_b, const TVector3& uphi_b, const TVector3& ueta_b,
        double& pn_b, double& pph_b, double& pet_b,
        double& aphi_b, double& aeta_b,
        double& cosphi_b, double& coseta_b) -> bool
  {
    pn_b  = pF.Dot(un_b);
    pph_b = pF.Dot(uphi_b);
    pet_b = pF.Dot(ueta_b);

    const double apn_b      = std::max(1e-12, std::fabs(pn_b));
    const double aphi_mag_b = std::atan2(std::fabs(pph_b), apn_b);
    const double aeta_mag_b = std::atan2(std::fabs(pet_b), apn_b);

    const double sgn_phi_b = (pph_b >= 0.0) ? +1.0 : -1.0;
    const double sgn_eta_b = (pet_b >= 0.0) ? +1.0 : -1.0;

    aphi_b = sgn_phi_b * aphi_mag_b;
    aeta_b = sgn_eta_b * aeta_mag_b;

    cosphi_b = std::clamp(apn_b / std::sqrt(apn_b*apn_b + pph_b*pph_b),
                          1.0e-6, 1.0);
    coseta_b = std::clamp(apn_b / std::sqrt(apn_b*apn_b + pet_b*pet_b),
                          1.0e-6, 1.0);

    return (std::isfinite(aphi_b)  && std::isfinite(aeta_b) &&
            std::isfinite(cosphi_b) && std::isfinite(coseta_b));
  };

  // ────────────────────────────────────────────────────────────────────────────
  // Step 4a) FACE frame from (possibly sandboxed) front-face tangents
  // ────────────────────────────────────────────────────────────────────────────
  TVector3 un_face, uphi_face, ueta_face;
  bool have_face = false;
  {
    TVector3 n_face = ephi.Cross(eeta);
    const double nmag = n_face.Mag();
    if (nmag > 0.0)
    {
      n_face *= (1.0 / nmag);

      // Outward normal (away from IP): choose sign by n · C
      if (n_face.Dot(C) < 0.0) n_face = -n_face;
      un_face = n_face;

      // φ̂: start from eφ and remove any normal component; normalize
      uphi_face = ephi - un_face * ephi.Dot(un_face);
      const double uphim = uphi_face.Mag();
      if (uphim > 0.0)
      {
        uphi_face *= (1.0 / uphim);

        // η̂: start from eη, orthogonalize to {un,uphi}; normalize
        ueta_face = eeta
                  - un_face   * eeta.Dot(un_face)
                  - uphi_face * eeta.Dot(uphi_face);
        const double uetam = ueta_face.Mag();
        if (uetam > 0.0)
        {
          ueta_face *= (1.0 / uetam);

          // Ensure right-handed triad (n̂ × φ̂ → +η̂)
          if (un_face.Cross(uphi_face).Dot(ueta_face) < 0.0)
            ueta_face = -ueta_face;

          have_face = true;
        }
      }
    }
  }

  // ────────────────────────────────────────────────────────────────────────────
  // Step 4b) MECHANICAL frame from RawTowerGeomv5 rotations
  // ────────────────────────────────────────────────────────────────────────────
  TVector3 un_mech, uphi_mech, ueta_mech;
  bool have_mech = false;
  {
    const double rx = g.rotX;
    const double ry = g.rotY;
    const double rz = g.rotZ;

    const double cx = std::cos(rx), sx = std::sin(rx);
    const double cy = std::cos(ry), sy = std::sin(ry);
    const double cz = std::cos(rz), sz = std::sin(rz);

    auto applyRot = [&](const TVector3& v) -> TVector3
    {
      // Apply Rz * Ry * Rx to local vector v (column-vector convention).
      double vx = v.X(), vy = v.Y(), vz = v.Z();

      // Rx
      const double vx1 = vx;
      const double vy1 =  cx*vy - sx*vz;
      const double vz1 =  sx*vy + cx*vz;

      // Ry
      const double vx2 =  cy*vx1 + sy*vz1;
      const double vy2 =  vy1;
      const double vz2 = -sy*vx1 + cy*vz1;

      // Rz
      const double vx3 =  cz*vx2 - sz*vy2;
      const double vy3 =  sz*vx2 + cz*vy2;
      const double vz3 =  vz2;

      return TVector3(vx3, vy3, vz3);
    };

    // Local axes: ẑ = bar axis, x̂/ŷ = transverse directions
    TVector3 u_axis = applyRot(TVector3(0.,0.,1.)); // tower "z" axis
    TVector3 u_x    = applyRot(TVector3(1.,0.,0.)); // tower "x" axis
    TVector3 u_y    = applyRot(TVector3(0.,1.,0.)); // tower "y" axis

    // Normalize axis and enforce outward direction
    const double amag = u_axis.Mag();
    if (amag > 0.0)
    {
      u_axis *= (1.0 / amag);
      if (u_axis.Dot(C) < 0.0) u_axis = -u_axis;

      un_mech = u_axis;

      // φ̂: project rotated local x̂ into plane ⟂ axis; normalize
      uphi_mech = u_x - un_mech * u_x.Dot(un_mech);
      const double uphim = uphi_mech.Mag();
      if (uphim > 0.0)
      {
        uphi_mech *= (1.0 / uphim);

        // η̂: start from rotated local ŷ, orthogonalize to {un,uphi}; normalize
        ueta_mech = u_y
                  - un_mech   * u_y.Dot(un_mech)
                  - uphi_mech * u_y.Dot(uphi_mech);
        const double uetam = ueta_mech.Mag();
        if (uetam > 0.0)
        {
          ueta_mech *= (1.0 / uetam);

          // Ensure right-handed triad
          if (un_mech.Cross(uphi_mech).Dot(ueta_mech) < 0.0)
            ueta_mech = -ueta_mech;

          // Align mechanical φ-axis with measured face φ tangent
          TVector3 uphi_ref = ephi_face;
          const double refMag = uphi_ref.Mag();
          if (refMag > 0.0)
          {
            uphi_ref *= (1.0 / refMag);
            if (uphi_mech.Dot(uphi_ref) < 0.0)
            {
              uphi_mech = -uphi_mech;
              ueta_mech = -ueta_mech;
            }
          }

          have_mech = true;
        }
      }
    }
  }

  // ────────────────────────────────────────────────────────────────────────────
  // Step 5) Choose primary frame and compute incidence
  // ────────────────────────────────────────────────────────────────────────────
  TVector3 un, uphi, ueta;
  bool have_primary = false;

  switch (kIncidenceMode)
  {
    case EIncidenceMode::FACE:
      if (!have_face) return false;
      un = un_face; uphi = uphi_face; ueta = ueta_face;
      have_primary = true;
      break;

    case EIncidenceMode::MECH:
      if (!have_mech) return false;
      un = un_mech; uphi = uphi_mech; ueta = ueta_mech;
      have_primary = true;
      break;

    case EIncidenceMode::BOTH:
      // Prefer mechanical as the "physics" frame; fall back to FACE if needed.
      if (have_mech)
      {
        un = un_mech; uphi = uphi_mech; ueta = ueta_mech;
        have_primary = true;
      }
      else if (have_face)
      {
        un = un_face; uphi = uphi_face; ueta = ueta_face;
        have_primary = true;
      }
      break;
  }

  if (!have_primary) return false;

  double pn  = 0.0, pph  = 0.0, pet  = 0.0;
  double aphi = 0.0, aeta = 0.0;
  double cphi = 1.0, ceta = 1.0;

  if (!basis_to_incidence(un, uphi, ueta,
                          pn, pph, pet,
                          aphi, aeta, cphi, ceta))
  {
    return false;
  }

  a_phi_sgn = static_cast<float>(aphi);
  a_eta_sgn = static_cast<float>(aeta);
  cos_a_phi = static_cast<float>(cphi);
  cos_a_eta = static_cast<float>(ceta);

  // Cache for QA
  m_lastAlphaPhiSigned = a_phi_sgn;
  m_lastAlphaEtaSigned = a_eta_sgn;

  const bool ok = (std::isfinite(cos_a_phi) && std::isfinite(cos_a_eta) &&
                   std::isfinite(a_phi_sgn) && std::isfinite(a_eta_sgn));

    // ────────────────────────────────────────────────────────────────────────────
    // Step 6) Debug: tabulated FACE vs MECH (or FACE only) + sandbox QA
    //
    //  • QA bookkeeping (m_qas_*) runs regardless of m_incDbgLevel
    //    so PrintIncidenceSandboxQASummary() works even with debug off.
    //  • Per-call CSV ([INC] lines) are only printed when dbg==true.
    // ────────────────────────────────────────────────────────────────────────────
    {
      // Incidence in FACE and MECH frames (if available)
      double pn_face   = 0.0, pph_face   = 0.0, pet_face   = 0.0;
      double aphi_face = 0.0, aeta_face  = 0.0;
      double cphi_face = 1.0, ceta_face  = 1.0;

      double pn_mech   = 0.0, pph_mech   = 0.0, pet_mech   = 0.0;
      double aphi_mech = 0.0, aeta_mech  = 0.0;
      double cphi_mech = 1.0, ceta_mech  = 1.0;

      if (have_face)
      {
        bool ok_face_dbg = basis_to_incidence(
            un_face, uphi_face, ueta_face,
            pn_face, pph_face, pet_face,
            aphi_face, aeta_face, cphi_face, ceta_face);
        (void) ok_face_dbg;
      }

      if (have_mech)
      {
        bool ok_mech_dbg = basis_to_incidence(
            un_mech, uphi_mech, ueta_mech,
            pn_mech, pph_mech, pet_mech,
            aphi_mech, aeta_mech, cphi_mech, ceta_mech);
        (void) ok_mech_dbg;
      }

      // Geometry helpers for digestible output
      const double Rc = std::hypot(C.X(), C.Y());

      // eta_det : tower-center η as seen from IP (z=0)
      const double eta_det = (Rc > 0.0) ? std::asinh(C.Z() / Rc)         : 0.0;
      // eta_SD  : tower-center η as seen from the chosen vertex z = fVz
      const double eta_SD  = (Rc > 0.0) ? std::asinh((C.Z() - fVz) / Rc) : 0.0;

      const double phi_det = std::atan2(C.Y(), C.X());    // tower center φ (detector)
      const double phi_ray = std::atan2(pF.Y(), pF.X());  // ray direction φ from vertex

      // Differences (ray–detector). In the Step-1 sandbox we only require dphi≈0.
      const double dphi = phi_ray - phi_det;
      const double deta = eta_SD  - eta_det;  // kept for information; not a Step-1 QA cut

        // Magnitudes of the FACE / MECH incidence vectors in (φ,η) plane
        const double amag_face = std::sqrt(aphi_face*aphi_face + aeta_face*aeta_face);
        const double amag_mech = std::sqrt(aphi_mech*aphi_mech + aeta_mech*aeta_mech);

        // FACE–MECH incidence differences in (φ,η) plane.
        // These are only physically meaningful when both frames are available
        // (have_face && have_mech), but we define them for all calls for simplicity.
        double dalpha_phi = 0.0;
        double dalpha_eta = 0.0;
        double dalpha_mag = 0.0;
        if (have_face && have_mech)
        {
          dalpha_phi = aphi_mech - aphi_face;
          dalpha_eta = aeta_mech - aeta_face;
          dalpha_mag = std::sqrt(dalpha_phi*dalpha_phi + dalpha_eta*dalpha_eta);
        }

        // Simple tagging / "color coding" via textual labels:
        //  - FACE≈0     : FACE angles essentially zero (perfect cylinder sandbox)
        //  - FACE≈MECH  : FACE and MECH magnitudes agree within ~0.3°
        //  - TILT       : genuinely tilted (mechanical incidence differs from FACE)
        std::string tag = "TILT";
        const double tolZero  = 1.0e-3;  // ~0.057° in radians
        const double tolMatch = 5.0e-3;  // ~0.29° in radians

      if (std::fabs(aphi_face) < tolZero && std::fabs(aeta_face) < tolZero)
      {
        tag = "FACE≈0";
      }
      else if (std::fabs(amag_face - amag_mech) < tolMatch)
      {
        tag = "FACE≈MECH";
      }

      // --- Step-1 sandbox QA bookkeeping: FACE-only, no-tilt sandbox ---
      // Runs irrespective of debug level; PASS/FAIL is decided later.
      if (m_incNoTiltSandbox &&
          kIncidenceMode == EIncidenceMode::FACE)
      {
        ++m_qas_nCalls;
        m_qas_maxAbsDphi    = std::max(m_qas_maxAbsDphi,    std::fabs(dphi));
        m_qas_maxAbsAlphaPh = std::max(m_qas_maxAbsAlphaPh, std::fabs(aphi_face));
        m_qas_maxAbsAlphaEt = std::max(m_qas_maxAbsAlphaEt, std::fabs(aeta_face));

        if (tag == "FACE≈0")
        {
          ++m_qas_nFaceZero;
        }
      }

      // Only print the CSV-style [INC] lines when debug is enabled
      if (dbg)
      {
        // Print header once (CSV-style; easy to grep / load in Python/ROOT)
        static bool s_headerPrinted = false;
          if (!s_headerPrinted)
          {
            std::cout << "# Columns ([INC] lines)\n"
                      << "#  vtxZ_cm, ix, iy, "
                      << "eta_det_IP, phi_det, eta_SD_vtx, phi_ray, "
                      << "dphi, deta, "
                      << "rotX_mrad, rotY_mrad, rotZ_mrad, "
                      << "a_phi_face, a_eta_face, |a_face|, cos_phi_face, cos_eta_face, "
                      << "a_phi_mech, a_eta_mech, |a_mech|, cos_phi_mech, cos_eta_mech, "
                      << "dalpha_phi, dalpha_eta, |dalpha|, "
                      << "tag\n";
            s_headerPrinted = true;
          }


          std::cout << "[INC] "
                    << fVz             << ", "
                    << ix              << ", "
                    << iy              << ", "
                    << eta_det         << ", "
                    << phi_det         << ", "
                    << eta_SD          << ", "
                    << phi_ray         << ", "
                    << dphi            << ", "
                    << deta            << ", "
                    << g.rotX * 1.0e3  << ", "
                    << g.rotY * 1.0e3  << ", "
                    << g.rotZ * 1.0e3  << ", "
                    << aphi_face       << ", "
                    << aeta_face       << ", "
                    << amag_face       << ", "
                    << cphi_face       << ", "
                    << ceta_face       << ", "
                    << aphi_mech       << ", "
                    << aeta_mech       << ", "
                    << amag_mech       << ", "
                    << cphi_mech       << ", "
                    << ceta_mech       << ", "
                    << dalpha_phi      << ", "
                    << dalpha_eta      << ", "
                    << dalpha_mag      << ", "
                    << tag             << "\n";
      }
    }

    // Keep the very detailed vector dump only for high debug levels
    if (dbg && m_incDbgLevel >= 10)
    {
      dump_and_maybe_stop("post-compute", C, F, V, un, uphi, ueta, pF,
                          pn, pph, pet, a_phi_sgn, a_eta_sgn, cos_a_phi, cos_a_eta);
    }

    return ok;
}

// ----------------------------------------------------------------------
// Step-1 sandbox QA summary: FACE + no-tilt test
//   • On PASS: print a single compact line.
//   • On FAIL: print full multi-line details (as before).
// ----------------------------------------------------------------------
void BEmcRecCEMC::PrintIncidenceSandboxQASummary(double tolDphi,
                                                 double tolAlpha) const
{
  const bool pass =
      (m_qas_nCalls > 0 &&
       m_qas_nFaceZero == m_qas_nCalls &&
       m_qas_maxAbsDphi    < tolDphi &&
       m_qas_maxAbsAlphaPh < tolAlpha &&
       m_qas_maxAbsAlphaEt < tolAlpha);

  if (pass)
  {
    std::cout << "\n[QA_STEP1] STATUS=PASS (FACE,noTilt)"
              << "  nCalls="       << m_qas_nCalls
              << "  max|dphi|="    << m_qas_maxAbsDphi    << " rad"
              << "  max|αφ_FACE|=" << m_qas_maxAbsAlphaPh << " rad"
              << "  max|αη_FACE|=" << m_qas_maxAbsAlphaEt << " rad"
              << "\n\n";
    return;
  }

  std::string status = "FAIL";

  std::cout << "\n[QA_STEP1] incidence sandbox (FACE,noTilt) summary\n"
            << "  nCalls       = " << m_qas_nCalls        << "\n"
            << "  nFACE≈0      = " << m_qas_nFaceZero     << "\n"
            << "  max|dphi|    = " << m_qas_maxAbsDphi    << " rad\n"
            << "  max|αφ_FACE| = " << m_qas_maxAbsAlphaPh << " rad\n"
            << "  max|αη_FACE| = " << m_qas_maxAbsAlphaEt << " rad\n"
            << "  STATUS       = " << status << "\n\n";
}




void BEmcRecCEMC::CorrectPositionEnergyAwareEnergyDepAndIncidentAngle(
    float Energy, float x, float y, float& xc, float& yc)
{
  // ---- gating ---------------------------------------------------------
  if (!m_UseCorrectPosition ||
      !std::isfinite(Energy) || !std::isfinite(x) || !std::isfinite(y) ||
      Energy < 0.01f)
  { xc = x; yc = y; return; }

  // ---- baseline b(E) at normal incidence (from master fits) ----------
  constexpr float E0 = 3.0f;
  const auto* fitPhi = BFitDB::instance().get("energyDepOnly","phi","originalEta","originalZRange");
  const auto* fitPhi = BFitDB::instance().get("energyDepOnly","phi","originalEta","originalZRange");
  const auto* fitEta = BFitDB::instance().get("energyDepOnly","eta","originalEta","originalZRange");
  const float bphi_E = b_from(fitPhi, Energy, E0);
  const float beta_E = b_from(fitEta, Energy, E0);

  VLOG() << "[EA+α] ENTER  E=" << Energy << "  (x,y)=(" << x << "," << y << ")"
         << "  bphi_E=" << bphi_E << "  beta_E=" << beta_E << "\n";

  // defaults (no incidence): keep baseline b(E), zero cached angles
  float bx = bphi_E, by = beta_E;
  m_lastAlphaPhi = 0.f; m_lastAlphaEta = 0.f;
  m_lastAlphaPhiSigned = 0.f; m_lastAlphaEtaSigned = 0.f;

  // ---- modular SD-incidence (seam-safe, detailed geometry) -----------
  {
    float cos_a_phi = 1.f, cos_a_eta = 1.f, a_phi_sgn = 0.f, a_eta_sgn = 0.f;

    if (ComputeIncidenceSD(Energy, x, y, cos_a_phi, cos_a_eta, a_phi_sgn, a_eta_sgn))
    {
      // cache angles (unsigned and signed)
      const float cphi = std::clamp(cos_a_phi, 1e-6f, 1.0f);
      const float ceta = std::clamp(cos_a_eta, 1e-6f, 1.0f);
      m_lastAlphaPhi       = std::acos(cphi);
      m_lastAlphaEta       = std::acos(ceta);
      m_lastAlphaPhiSigned = a_phi_sgn;
      m_lastAlphaEtaSigned = a_eta_sgn;

      // foreshortening: b_eff = b(E) / cos α
      bx = clamp_b(bphi_E / cphi);
      by = clamp_b(beta_E / ceta);

      VLOG() << "        [EA+α] cos(a_phi)=" << cphi << "  cos(a_eta)=" << ceta
             << "  → a_phi=" << m_lastAlphaPhi << "  a_eta=" << m_lastAlphaEta
             << "  → bx=" << bx << "  by=" << by << "\n";
    }
    else
    {
      VLOG() << "[EA+α] incidence unavailable → using baseline b(E)\n";
    }
  }

  VLOG() << "[EA+α] EXIT  bphi_E=" << bphi_E << "  beta_E=" << beta_E
         << "  → bx=" << bx << "  by=" << by << "\n";

    // ============================== φ ===================================
    // Apply signed δ_φ(E,αϕ) BEFORE the inverse mapping (guarded by a tiny DeltaDB)
    float x_pre = x;
  #ifdef PDC_USE_DELTA_RUNTIME
    if (DeltaDB::instance().available()) {
      const float c1_phi   = DeltaDB::instance().c1_phi(Energy);    // c0_phi + m_phi*ln(E/E0)
      const float delta_phi= c1_phi * std::sin(m_lastAlphaPhiSigned);
      x_pre = x - delta_phi;
      VLOG() << "        [EA+α] δφ(E,α)=" << delta_phi << "  → x_pre=" << x_pre << "\n";
    }
  #endif

    float x_corr = x_pre;
    {
      const int   ix0_local = EmcCluster::lowint(x_pre + 0.5f);
      const float X         = x_pre - ix0_local;            // in (−0.5, +0.5]
      if (std::fabs(X) <= 0.5f) {
        const float Sx = std::sinh(0.5f / bx);
        x_corr = ix0_local + bx * std::asinh(2.f * X * Sx);
      }

      // module-of-8 ripple (unchanged), evaluated at x_pre
      int   ix8 = int(x_pre + 0.5f) / 8;                     // NOLINT
      float x8  = x_pre + 0.5f - (ix8 * 8) - 4.f;            // −4 … +4
      float dx  = 0.f;
      if (m_UseDetailedGeometry) {
        int local_ix8 = int(x_pre + 0.5f) - ix8 * 8;         // 0..7   // NOLINT
        dx = static_cast<float>(factor_[local_ix8]) * (x8 / 4.f);
      } else {
        dx = (std::fabs(x8) > 3.3f) ? 0.f : 0.10f * (x8 / 4.f);
      }
      x_corr -= dx;

      // wrap φ tower coordinate to [−0.5, Nx−0.5)
      while (x_corr < -0.5f)              x_corr += float(fNx);
      while (x_corr >= float(fNx) - 0.5f) x_corr -= float(fNx);
    }
    xc = x_corr;


    // ============================== η ===================================
    // Apply signed δ_η(E,αη) BEFORE the inverse mapping
    float y_pre = y;
  #ifdef PDC_USE_DELTA_RUNTIME
    if (DeltaDB::instance().available()) {
      const float c1_eta   = DeltaDB::instance().c1_eta(Energy);     // c0_eta + m_eta*ln(E/E0)
      const float delta_eta= c1_eta * std::sin(m_lastAlphaEtaSigned);
      y_pre = y - delta_eta;
      VLOG() << "        [EA+α] δη(E,α)=" << delta_eta << "  → y_pre=" << y_pre << "\n";
    }
  #endif

    float y_corr = y_pre;
    {
      const int   iy0_local = EmcCluster::lowint(y_pre + 0.5f);
      const float Y         = y_pre - iy0_local;            // in (−0.5, +0.5]
      if (std::fabs(Y) <= 0.5f) {
        const float Sy = std::sinh(0.5f / by);
        y_corr = iy0_local + by * std::asinh(2.f * Y * Sy);
      }
    }
    yc = y_corr;
}





void BEmcRecCEMC::CorrectPositionEnergyAwareEnergyDepOnly(float Energy, float x, float y,
                                             float& xc, float& yc)
{
  // ---- legacy gating -------------------------------------------------
  if (!m_UseCorrectPosition) { xc = x; yc = y; return; }
  if (!std::isfinite(Energy) || !std::isfinite(x) || !std::isfinite(y) || Energy < 0.01f)
  { xc = x; yc = y; return; }

    constexpr float E0 = 3.0f;

    // energyDepOnly, originalEta × originalZRange
    const auto* fitPhi = BFitDB::instance().get("energyDepOnly", "phi", "originalEta", "originalZRange");
    const auto* fitEta = BFitDB::instance().get("energyDepOnly", "eta", "originalEta", "originalZRange");

    const float bx = b_from(fitPhi, Energy, E0);  // φ-direction b(E)
    const float by = b_from(fitEta, Energy, E0);  // η-direction b(E)
    

  // ============================== φ ===================================
  float x_corr = x;
  {
    const int ix0 = EmcCluster::lowint(x + 0.5f);
    const float X = x - ix0;                       // in (−0.5, +0.5]
    if (std::fabs(X) <= 0.5f)
    {
      const float Sx = std::sinh(0.5f / bx);
      x_corr = ix0 + bx * std::asinh(2.f * X * Sx);
    }

    // module-of-8 ripple (identical to legacy CorrectPosition)
    // NOLINTNEXTLINE(bugprone-incorrect-roundings)
    int   ix8 = int(x + 0.5f) / 8;
    float x8  = x + 0.5f - (ix8 * 8) - 4.f;        // −4 … +4
    float dx  = 0.f;
    if (m_UseDetailedGeometry)
    {
      // NOLINTNEXTLINE(bugprone-incorrect-roundings)
      int local_ix8 = int(x + 0.5f) - ix8 * 8;     // 0..7
      dx = static_cast<float>(factor_[local_ix8]) * (x8 / 4.f);
    }
    else
    {
      dx = (std::fabs(x8) > 3.3f) ? 0.f : 0.10f * (x8 / 4.f);
    }
    x_corr -= dx;

    // wrap φ tower coordinate to [−0.5, Nx−0.5)
    while (x_corr < -0.5f)              x_corr += float(fNx);
    while (x_corr >= float(fNx) - 0.5f) x_corr -= float(fNx);
  }
  xc = x_corr;

  // ============================== η ===================================
  float y_corr = y;
  {
    const int iy0 = EmcCluster::lowint(y + 0.5f);
    const float Y = y - iy0;                      // in (−0.5, +0.5]
    if (std::fabs(Y) <= 0.5f)
    {
      const float Sy = std::sinh(0.5f / by);
      y_corr = iy0 + by * std::asinh(2.f * Y * Sy);
    }
  }
  yc = y_corr;
}
