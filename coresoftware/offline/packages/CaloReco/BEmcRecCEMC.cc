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
    std::string path = env ? std::string(env) : std::string("/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/bFit_master.txt");
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
      {ETiltVariant::CLUS_RAW, { 2.641645e-03, 8.117525e-04 }},  // no corr, cluster
      {ETiltVariant::CLUS_CP, { 1.925290e-03, 7.915014e-04 }},  // CorrectPosition, cluster
      {ETiltVariant::CLUS_CP_EA_FIT_ZDEP_ETADEP, { 2.091581e-03, 8.209508e-04 }},  // CorrectPosition(EA |z|+|#eta|+E), cluster
      {ETiltVariant::CLUS_CP_EA_FIT_ZDEP, { 2.085659e-03, 8.202723e-04 }},  // CorrectPosition(EA |z|+E), cluster
      {ETiltVariant::CLUS_CP_EA_FIT_ETADEP, { 2.091938e-03, 8.213024e-04 }},  // CorrectPosition(EA |#eta|+E), cluster
      {ETiltVariant::CLUS_CP_EA_FIT_EandIncident, { 2.162642e-03, 8.224042e-04 }},  // CorrectPosition(EA E-only + incident-angle), cluster
      {ETiltVariant::CLUS_CP_EA_FIT_EONLY, { 2.086863e-03, 8.207850e-04 }},  // CorrectPosition(EA E-only), cluster
      {ETiltVariant::CLUS_BCORR, { 2.121298e-03, 8.136292e-04 }},  // b(E) corr, cluster
      {ETiltVariant::PDC_RAW, { 2.739621e-03, 8.204782e-04 }},  // no corr, scratch
      {ETiltVariant::PDC_CORR, { 2.244868e-03, 8.255293e-04 }},  // b(E) corr, scratch

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



void BEmcRecCEMC::CorrectPositionEnergyAwareEnergyDepAndIncidentAngle(float Energy, float x, float y,
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
  // base (angle–zero) scales from the master-file fits
  const float bphi_E = b_from(fitPhi, Energy, E0);  // φ baseline
  const float beta_E = b_from(fitEta,  Energy, E0); // η baseline
  VLOG() << "[EA+α] ENTER  E=" << Energy << "  (x,y)=(" << x << "," << y << ")"
         << "  bphi_E=" << bphi_E << "  beta_E=" << beta_E << "\n";
  // --- geometry-only incidence from detailed tower geometry -----------------------
  // use the same integer anchor as your local tower coordinates (no extra lookup later)
  const int ix_geom = EmcCluster::lowint(x + 0.5f);
  const int iy_geom = EmcCluster::lowint(y + 0.5f);
  // fetch per-tower geometry already filled by SetTowerGeometry(...)
  TowerGeom g{};
  bool haveGeom = GetTowerGeometry(ix_geom, iy_geom, g);   // your existing helper
  float bx = bphi_E;   // defaults (if geometry missing)
  float by = beta_E;
  if (!haveGeom)
  {
    VLOG() << "[EA+α] geometry MISSING for (ix,iy)=(" << ix_geom << "," << iy_geom
           << ")  → using bphi_E/beta_E baselines only\n";
  }
  else
  {
    VLOG() << "[EA+α] geometry OK  (ix,iy)=(" << ix_geom << "," << iy_geom << ")\n"
           << "        C=(" << g.Xcenter << "," << g.Ycenter << "," << g.Zcenter << ")\n"
           << "        dPhi=(" << g.dX[0] << "," << g.dY[0] << "," << g.dZ[0] << ")"
           << "  dEta=("  << g.dX[1] << "," << g.dY[1] << "," << g.dZ[1] << ")\n";
    // --- shower-depth point and ray from actual vertex (energy & local-offset aware) ---
    TVector3 C(g.Xcenter, g.Ycenter, g.Zcenter);
    TVector3 V(0., 0., fVz);
    float x_sd = 0.f, y_sd = 0.f, z_sd = 0.f;
    Tower2Global(Energy, x, y, x_sd, y_sd, z_sd);   // maps (x,y,E) → shower-depth (cm)
    TVector3 H(x_sd, y_sd, z_sd);
    // Ray to shower depth; fall back to tower-center ray only if degenerate
    TVector3 p = (H - V);
    if (p.Mag2() <= 0.0) p = (C - V);
    p = p.Unit();
    // --- raw face-edge directions (provided by detailed geometry) ---
    TVector3 ephi_raw(g.dX[0], g.dY[0], g.dZ[0]);
    TVector3 eeta_raw(g.dX[1], g.dY[1], g.dZ[1]);
    if (ephi_raw.Mag2() < 1e-24 || eeta_raw.Mag2() < 1e-24)
    {
      VLOG() << "[EA+α] degenerate face tangents → skip angle scaling this event\n";
      m_lastAlphaPhi = 0.f;
      m_lastAlphaEta = 0.f;
    }
    else
    {
      // --- mechanical tower axis from face normals (right-handed) ---
      TVector3 n = (ephi_raw.Cross(eeta_raw)).Unit();         // axis from faces (not -C)
      if (n.Dot(-C) < 0.0) n = -n;                            // orient toward IP
      VLOG() << "        n=(" << n.X() << "," << n.Y() << "," << n.Z() << ")"
             << "  p=(" << p.X() << "," << p.Y() << "," << p.Z() << ")\n";
      // --- BUILD φ/η TANGENTS AND NORMAL AT SHOWER DEPTH (projective, depth-aware) ---
      float x_phip=0.f,y_phip=0.f,z_phip=0.f, x_phim=0.f,y_phim=0.f,z_phim=0.f;
      float x_etap=0.f,y_etap=0.f,z_etap=0.f, x_etam=0.f,y_etam=0.f,z_etam=0.f;
      // step ±0.5 tower in φ at the SAME (x,y) anchor, map both to shower depth
      Tower2Global(Energy, x + 0.5f, y, x_phip, y_phip, z_phip);
      Tower2Global(Energy, x - 0.5f, y, x_phim, y_phim, z_phim);
      TVector3 ephi_d(x_phip - x_phim, y_phip - y_phim, z_phip - z_phim);
      // step ±0.5 tower in η at the SAME (x,y) anchor, map both to shower depth
      Tower2Global(Energy, x, y + 0.5f, x_etap, y_etap, z_etap);
      Tower2Global(Energy, x, y - 0.5f, x_etam, y_etam, z_etam);
      TVector3 eeta_d(x_etap - x_etam, y_etap - y_etam, z_etap - z_etam);
      const double Lphi = ephi_d.Mag();
      const double Leta = eeta_d.Mag();
      if (Lphi < 1e-12 || Leta < 1e-12)
      {
        VLOG() << "[EA+α] depth tangents too small → skip angle scaling\n";
        m_lastAlphaPhi = 0.f;
        m_lastAlphaEta = 0.f;
      }
      else
      {
        // ------------------ NEW: orthonormalize the depth tangents ------------------
        // unit φ tangent at depth
        TVector3 uphi = ephi_d * (1.0 / Lphi);
        // start from unit η tangent and remove any φ component (Gram–Schmidt in-plane)
        TVector3 ueta_raw = eeta_d * (1.0 / Leta);
        TVector3 ueta_t   = ueta_raw - (ueta_raw.Dot(uphi)) * uphi;
        double   Leta_t   = ueta_t.Mag();
        TVector3 nd;    // depth-surface normal (orthogonal to both uphi and ueta)
        TVector3 ueta;  // final orthonormal η direction at depth
        if (Leta_t < 1e-12)
        {
          // Rare near-collinearity: build a stable in-plane orthobasis using n as helper
          nd = (uphi.Cross(n)).Unit();
          if (nd.Mag2() < 1e-24)
          {
            // pathological fallback: use a fixed axis to break degeneracy (won't be hit in practice)
            nd = (uphi.Cross(TVector3(0.,0.,1.))).Unit();
          }
          ueta = (nd.Cross(uphi)).Unit();
        }
        else
        {
          ueta = ueta_t * (1.0 / Leta_t);
          nd   = (uphi.Cross(ueta)).Unit();
        }
        // orient normal toward the actual shower ray
        if (nd.Dot(p) < 0.0) nd = -nd;
        // enforce right-handed orthonormal triad exactly
        ueta = (nd.Cross(uphi)).Unit();
        // decompose ray in the orthonormal basis {nd, uphi, ueta} at depth
        const double pn   = std::abs(p.Dot(nd));
        const double pphi = std::abs(p.Dot(uphi));
        const double peta = std::abs(p.Dot(ueta));
        VLOG() << "        [depth⊥] unit-closure=" << (pn*pn + pphi*pphi + peta*peta)
               << "  pn=" << pn << "  pphi=" << pphi << "  peta=" << peta << "\n";
        // incidence angles: cos α_dir = |p·nd| / ||proj_{span{nd,dir}}(p)||  (at depth, orthonormal)
        const double cos_a_phi = std::max(1e-6, pn / std::sqrt(pn*pn + pphi*pphi));
        const double cos_a_eta = std::max(1e-6, pn / std::sqrt(pn*pn + peta*peta));
        m_lastAlphaPhi = static_cast<float>(std::acos(std::min(1.0, cos_a_phi)));
        m_lastAlphaEta = static_cast<float>(std::acos(std::min(1.0, cos_a_eta)));
        VLOG() << "        [depth⊥] cos(a_phi)=" << cos_a_phi << "  cos(a_eta)=" << cos_a_eta
               << "  → a_phi=" << m_lastAlphaPhi << " rad"
               << "  a_eta=" << m_lastAlphaEta << " rad\n";
        // Pure geometric foreshortening at depth (no tunables, no extra modeling)
        bx = clamp_b( bphi_E / static_cast<float>(cos_a_phi) );
        by = clamp_b(  beta_E / static_cast<float>(cos_a_eta) );
        VLOG() << "        [depth⊥] b_eff: bx=" << bx << "  by=" << by << "\n";
        // ---------------- end NEW orthonormalization refinement ---------------------
      }
    }
  }
  VLOG() << "[EA+α] EXIT  bphi_E=" << bphi_E << "  beta_E=" << beta_E
         << "  → bx=" << bx << "  by=" << by << "\n";
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
