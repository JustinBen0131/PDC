#include "BEmcRecCEMC.h"

#include "BEmcCluster.h"
#include "BEmcProfile.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <limits>  // for quiet_NaN()


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
  static const std::unordered_map<ETiltVariant,std::pair<double,double>> lut = {
      {ETiltVariant::CLUS_RAW,    { 2.687127e-03, 8.367887e-04}},  // no corr, cluster
      {ETiltVariant::CLUS_CP,     { 2.001589e-03, 8.143681e-04}},  // CorrectPosition, cluster
      {ETiltVariant::CLUS_CP_EA,  { 2.001589e-03, 8.143681e-04}},  // Energy/angle-aware CP (start equal to CP)
      {ETiltVariant::CLUS_BCORR,  { 2.201905e-03, 8.409767e-04}},  // b(E) corr, cluster
      {ETiltVariant::PDC_RAW,     { 1.236497e-03, 7.982602e-04}},  // no corr, scratch
      {ETiltVariant::PDC_CORR,    { 8.076897e-04, 8.185800e-04}},  // b(E) corr, scratch
      {ETiltVariant::DEFAULT,     { 8.654924e-04, 8.490399e-04}}   // legacy numbers
  };

  auto it = lut.find(v);
  m_abTiltCurrent = (it!=lut.end()) ? it->second : lut.at(ETiltVariant::DEFAULT);
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


void BEmcRecCEMC::CorrectPositionEnergyAware(float Energy, float x, float y,
                                             float& xc, float& yc,
                                             float* out_bphi, float* out_beta)
{
  // initialize optional outputs (so caller never sees stale data)
  if (out_bphi) *out_bphi = std::numeric_limits<float>::quiet_NaN();
  if (out_beta) *out_beta = std::numeric_limits<float>::quiet_NaN();

  // -------- pass-through guards
  xc = x; yc = y;
  if (!std::isfinite(Energy) || !std::isfinite(x) || !std::isfinite(y)) return;
  if (Energy < 0.01f) return;

  // -------- small helpers
  auto clampf = [](float v, float lo, float hi){ return (v < lo) ? lo : (v > hi) ? hi : v; };
  auto norm3  = [&](float v[3]) { float n = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); if (n<1e-12f) n=1e-12f; v[0]/=n; v[1]/=n; v[2]/=n; };
  auto proj_perp_len = [&](const float aHat[3], float vx, float vy, float vz)->float {
    const float k = vx*aHat[0] + vy*aHat[1] + vz*aHat[2];
    const float px = vx - k*aHat[0];
    const float py = vy - k*aHat[1];
    const float pz = vz - k*aHat[2];
    return std::sqrt(px*px + py*py + pz*pz);
  };
  auto wrap_phi = [&](float v)->float {
    if (!(fNx > 0) || !std::isfinite(v)) return v;
    const float lo = -0.5f, hi = float(fNx) - 0.5f, span = float(fNx);
    float t = std::fmod(v - lo, span); if (t < 0.f) t += span;
    float r = lo + t;
    // open-right guard: map exact top back to lo
    if (std::fabs(r - hi) < 1e-7f) r = lo;
    return r;
  };
  auto dot3 = [&](const float a[3], const float b[3]){ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; };

  // -------- indices & geometry
  const int ixHome = EmcCluster::lowint(x + 0.5f);
  const int iyHome = EmcCluster::lowint(y + 0.5f);

  TowerGeom g{};
  if (!GetTowerGeometry(ixHome, iyHome, g)) {
    // conservative fallback: legacy constant b
    const float b = 0.15f;

    // report b-values in fallback
    if (out_bphi) *out_bphi = b;
    if (out_beta) *out_beta = b;

    // φ inverse-asinh (no zero-shift)
    float x0 = x; const int ix0 = EmcCluster::lowint(x0 + 0.5f);
    if (std::fabs(x0 - ix0) <= 0.5f) {
      const float Sx = std::sinh(0.5f / b);
      x0 = ix0 + b * std::asinh( 2.f * (x0 - ix0) * Sx );
    }

    // ripple
    int   ix8 = int(x + 0.5f) / 8;
    float x8  = x + 0.5f - (ix8 * 8) - 4.f;
    float dx  = m_UseDetailedGeometry
              ? static_cast<float>(factor_[int(x+0.5f)-ix8*8]) * (x8/4.f)
              : (std::fabs(x8) > 3.3f ? 0.f : 0.10f * (x8/4.f));
    xc = wrap_phi(x0 - dx);

    // η inverse-asinh (no zero-shift)
    float y0 = y; const int iy0 = EmcCluster::lowint(y0 + 0.5f);
    if (std::fabs(y0 - iy0) <= 0.5f) {
      const float Sy = std::sinh(0.5f / b);
      y0 = iy0 + b * std::asinh( 2.f * (y0 - iy0) * Sy );
    }
    yc = y0;
    return;
  }

  // ============================================================
  // Local basis & shower axis at depth
  // ============================================================
  float ephi[3] = { g.dX[0], g.dY[0], g.dZ[0] };
  float eeta[3] = { g.dX[1], g.dY[1], g.dZ[1] };
  norm3(ephi); norm3(eeta);
  float erad[3] = { ephi[1]*eeta[2] - ephi[2]*eeta[1],
                    ephi[2]*eeta[0] - ephi[0]*eeta[2],
                    ephi[0]*eeta[1] - ephi[1]*eeta[0] };
  norm3(erad);

  float r0E[3];
  CorrectShowerDepth(ixHome, iyHome, Energy, g.Xcenter, g.Ycenter, g.Zcenter,
                     r0E[0], r0E[1], r0E[2]);

  float aTmp[3] = { r0E[0], r0E[1], r0E[2] - fVz };
  norm3(aTmp);
  const float aHat[3] = { aTmp[0], aTmp[1], aTmp[2] };

  // ---- small high-E depth reference (shared by φ & η tempering)
  constexpr float E0_HE_REF = 10.0f;
  float r0E_ref[3];
  CorrectShowerDepth(ixHome, iyHome, E0_HE_REF, g.Xcenter, g.Ycenter, g.Zcenter,
                     r0E_ref[0], r0E_ref[1], r0E_ref[2]);
  float aTmpRef[3] = { r0E_ref[0], r0E_ref[1], r0E_ref[2] - fVz };
  norm3(aTmpRef);
  const float aHatRef[3] = { aTmpRef[0], aTmpRef[1], aTmpRef[2] };

  const float rF[3]  = { g.Xcenter, g.Ycenter, g.Zcenter };
  const float drE[3]   = { r0E[0]-rF[0],    r0E[1]-rF[1],    r0E[2]-rF[2] };
  const float drRef[3] = { r0E_ref[0]-rF[0], r0E_ref[1]-rF[1], r0E_ref[2]-rF[2] };
  const float depthE   = std::fabs(dot3(drE,   aHat));
  const float depthRef = std::fabs(dot3(drRef, aHatRef));

  // Neighbor indices
  const int ixP   = (fNx > 0) ? ((ixHome + 1 + fNx) % fNx) : ixHome; // +φ
  const int ixM   = (fNx > 0) ? ((ixHome - 1 + fNx) % fNx) : ixHome; // −φ  (PATCH)
  const int iyP   = (iyHome + 1 < fNy) ? (iyHome + 1) : iyHome;       // +η
  const int iyDn  = (iyHome > 0) ? (iyHome - 1) : iyHome;             // −η

  // ============================================================
  // φ side — high‑E refinement (kept) + PATCH: ±φ symmetry at high E
  // ============================================================
  constexpr float b0 = 0.15f;

  const float dPhF[3] = { g.dX[0], g.dY[0], g.dZ[0] }; // φ step on the face

  // transport a *face-advanced* point (home + dPhF) to depth at +φ index
  const float cPhP[3] = { g.Xcenter + dPhF[0], g.Ycenter + dPhF[1], g.Zcenter + dPhF[2] };
  float rPhE_p[3];
  CorrectShowerDepth(ixP, iyHome, Energy, cPhP[0], cPhP[1], cPhP[2],
                     rPhE_p[0], rPhE_p[1], rPhE_p[2]);

  // PATCH: mirror to −φ using a face-retreated point (home − dPhF)
  const float cPhM[3] = { g.Xcenter - dPhF[0], g.Ycenter - dPhF[1], g.Zcenter - dPhF[2] };
  float rPhE_m[3];
  CorrectShowerDepth(ixM, iyHome, Energy, cPhM[0], cPhM[1], cPhM[2],
                     rPhE_m[0], rPhE_m[1], rPhE_m[2]);

  // widths: ⟂ aHat (baseline)
  const float wPh_front_perp = proj_perp_len(aHat, dPhF[0], dPhF[1], dPhF[2]);
  const float dPhE_p[3]      = { rPhE_p[0]-r0E[0], rPhE_p[1]-r0E[1], rPhE_p[2]-r0E[2] };
  const float dPhE_m[3]      = { rPhE_m[0]-r0E[0], rPhE_m[1]-r0E[1], rPhE_m[2]-r0E[2] };

  const float wPh_depth_perp_p = proj_perp_len(aHat, dPhE_p[0], dPhE_p[1], dPhE_p[2]);
  const float wPh_depth_perp_m = proj_perp_len(aHat, dPhE_m[0], dPhE_m[1], dPhE_m[2]);

  // along-eφ widths
  const float wPh_front_along = std::fabs(dot3(dPhF, ephi));
  const float wPh_depth_along_p = std::fabs(dot3(dPhE_p, ephi));
  const float wPh_depth_along_m = std::fabs(dot3(dPhE_m, ephi));

  // single-sided sPhi (current baseline)
  const float sPhi_perp_single  = (wPh_depth_perp_p  > 1e-6f) ? (wPh_front_perp  / wPh_depth_perp_p ) : 1.0f;
  const float sPhi_along_single = (wPh_depth_along_p > 1e-6f) ? (wPh_front_along / wPh_depth_along_p) : 1.0f;

  // PATCH: symmetric ±φ widths (arithmetic mean), gated to high E
  float wPh_depth_perp_sym  = 0.f, wPh_depth_along_sym = 0.f;
  if (wPh_depth_perp_p > 0.f && wPh_depth_perp_m > 0.f)
    wPh_depth_perp_sym = 0.5f * (wPh_depth_perp_p + wPh_depth_perp_m);
  else
    wPh_depth_perp_sym = (wPh_depth_perp_p > 0.f ? wPh_depth_perp_p :
                          (wPh_depth_perp_m > 0.f ? wPh_depth_perp_m : 0.f));

  if (wPh_depth_along_p > 0.f && wPh_depth_along_m > 0.f)
    wPh_depth_along_sym = 0.5f * (wPh_depth_along_p + wPh_depth_along_m);
  else
    wPh_depth_along_sym = (wPh_depth_along_p > 0.f ? wPh_depth_along_p :
                           (wPh_depth_along_m > 0.f ? wPh_depth_along_m : 0.f));

  const float sPhi_perp_sym  = (wPh_depth_perp_sym  > 1e-6f) ? (wPh_front_perp  / wPh_depth_perp_sym ) : 1.0f;
  const float sPhi_along_sym = (wPh_depth_along_sym > 1e-6f) ? (wPh_front_along / wPh_depth_along_sym) : 1.0f;

  // high‑E weight for φ blending / zero-shift / temper
  constexpr float E_PHI_HE_LO = 12.0f;
  constexpr float E_PHI_HE_HI = 25.0f;
  float wHE_phi = 0.0f;
  if (Energy > E_PHI_HE_LO) {
    const float t = (Energy - E_PHI_HE_LO) / (E_PHI_HE_HI - E_PHI_HE_LO);
    wHE_phi = (t < 0.f) ? 0.f : (t > 1.f) ? 1.f : t;
  }

  // PATCH: gently fade from single-sided → symmetric with high‑E
  const float wSYM_phi = 0.60f * wHE_phi; // modest weight
  const float sPhi_perp_eff  = (1.0f - wSYM_phi) * sPhi_perp_single  + wSYM_phi * sPhi_perp_sym;
  const float sPhi_along_eff = (1.0f - wSYM_phi) * sPhi_along_single + wSYM_phi * sPhi_along_sym;

  // your existing ⟂ vs along φ blend (kept)
  const float sPhi_blend = (1.0f - wHE_phi) * sPhi_perp_eff + wHE_phi * sPhi_along_eff;
  float bx = b0 * sPhi_blend;

  // --- tiny high‑E depth tempering (kept)
  constexpr float E_TEMPER_LO = 18.0f;
  constexpr float E_TEMPER_HI = 30.0f;
  float wHE_temper = 0.0f;
  if (Energy > E_TEMPER_LO) {
    const float t = (Energy - E_TEMPER_LO) / (E_TEMPER_HI - E_TEMPER_LO);
    wHE_temper = (t < 0.f) ? 0.f : (t > 1.f) ? 1.f : t;
  }
  if (depthE > 1e-3f && depthRef > 1e-3f) {
    const float ratio = clampf(depthRef / depthE, 0.85f, 1.15f);
    const float LAMBDA_HE = 0.12f; // very gentle
    const float boost = std::pow(ratio, LAMBDA_HE * wHE_temper);
    bx *= boost;
  }

  // clamp: baseline [0.14,0.20], tightened a hair at very high E
  const float bx_lo = (1.0f - wHE_phi)*0.14f + wHE_phi*0.135f; // → 0.135 at top
  const float bx_hi = (1.0f - wHE_phi)*0.20f + wHE_phi*0.190f; // → 0.190 at top
  bx = clampf(bx, bx_lo, bx_hi);

  // report bφ now (final)
  if (out_bphi) *out_bphi = bx;

  // tiny tilt-based zero-shift enabled only at high E (kept)
  float xZero = 0.f;
  if (wHE_phi > 0.f) {
    const float a_phi = dot3(aHat, ephi);
    const float a_rad = dot3(aHat, erad);
    const float den   = std::sqrt(a_phi*a_phi + a_rad*a_rad);
    const float sinTx = (den > 1e-9f) ? (a_phi / den) : 0.f;
    const float sin2  = sinTx * sinTx;
    const float amp   = clampf(0.9f * (wPh_front_along / (wPh_depth_along_p + 1e-6f)), 0.75f, 1.25f);
    float xZero_raw   = amp * ( (sinTx > 0.f) ? (-0.20f * sinTx - 0.60f * sin2)
                                              : (-0.20f * sinTx + 0.60f * sin2) );
    xZero_raw = clampf(xZero_raw, -0.30f, 0.30f);
    const int ix0 = EmcCluster::lowint(x + 0.5f);
    const float xShifted = x + xZero_raw;
    const float over = std::fabs(xShifted - ix0) - 0.499f;
    if (over > 0.f) xZero_raw -= std::copysign(over, xZero_raw);
    xZero = wHE_phi * xZero_raw;
  }

    // φ inverse with asymmetric window at shower depth (±φ), CORRECTLY centered.
        // Measure in owner cell: X = x - ix0. Apply zero-shift inside the inverse:
        // Xeff = X - xZero. Solve for t from Xeff; then add the shift back: ix0 + xZero + t.
        float x_corr = x; // start from measured x (no pre-shift)
        {
          const int   ix0  = EmcCluster::lowint(x + 0.5f);   // owner cell (kept fixed by earlier clamp)
          const float X    = x - ix0;                        // in (-0.5, +0.5]
          const float Xeff = X - xZero;                      // zero-shifted coord to invert

          if (std::fabs(X) <= 0.5f)
          {
            // symmetric closed-form (initializer)
            const float Sx    = std::sinh(0.5f / bx);
            const float t_sym = bx * std::asinh( 2.f * Xeff * Sx );

            auto safe_ratio = [&](float num, float den)->float { return (den > 1e-6f) ? (num/den) : 1.0f; };

              // side-resolved half-widths at depth and energy-gated blend
              const float sR_perp  = safe_ratio(wPh_front_perp , wPh_depth_perp_p );
              const float sL_perp  = safe_ratio(wPh_front_perp , wPh_depth_perp_m );
              const float sR_along = safe_ratio(wPh_front_along, wPh_depth_along_p);
              const float sL_along = safe_ratio(wPh_front_along, wPh_depth_along_m);

              const float sR_eff = (1.0f - wHE_phi) * sR_perp  + wHE_phi * sR_along;
              const float sL_eff = (1.0f - wHE_phi) * sL_perp  + wHE_phi * sL_along;

              float wR = clampf(0.5f * sR_eff, 0.35f, 0.80f);
              float wL = clampf(0.5f * sL_eff, 0.35f, 0.80f);

            // forward map for asymmetric window [-wL, +wR] in tower units
            auto X_from_t_asym = [&](float t)->float {
              const float eR = std::exp(-(wR - t)/bx);
              const float eL = std::exp(-(wL + t)/bx);
              const float I0 = 1.f - 0.5f*(eR + eL);
              const float I1 = 0.5f*((wL + t)*eL - (wR - t)*eR) + 0.5f*bx*(eL - eR);
              return t + I1 / std::max(I0, 1e-6f);
            };
            auto dXdt_num = [&](float t)->float {
              const float h = 1e-3f;
              return (X_from_t_asym(t + h) - X_from_t_asym(t - h)) / (2.f*h);
            };

            float t = t_sym;
            const float wASY = 0.70f * wHE_phi;
            if (wASY > 0.f)
            {
              for (int it = 0; it < 2; ++it)
              {
                const float F  = X_from_t_asym(t) - Xeff;   // solve for Xeff (not X+shift)
                const float dF = dXdt_num(t);
                if (std::fabs(dF) < 1e-4f) break;
                t -= F / dF;
                t  = clampf(t, -wL + 1e-4f, +wR - 1e-4f);
              }
              t = (1.f - wASY) * t_sym + wASY * t;
            }

            // shift back to absolute coordinate after inversion
            x_corr = ix0 + xZero + t;
          }
        }

    
  // φ ripple (module-of-8) + wrap (kept, with high‑E de‑emphasis)
  float dx = 0.f;
  {
    const int   ix8 = int(x + 0.5f) / 8;
    const float x8  = x + 0.5f - (ix8 * 8) - 4.f;
    if (m_UseDetailedGeometry) {
      const int local_ix8 = int(x + 0.5f) - ix8 * 8; // nominally 0..7
      const int FACTOR_LEN = int(sizeof(factor_) / sizeof(factor_[0]));
      if (local_ix8 >= 0 && local_ix8 < FACTOR_LEN)
        dx = static_cast<float>(factor_[local_ix8]) * (x8 / 4.f);
      else
        dx = 0.f;
    } else {
      dx = (std::fabs(x8) > 3.3f) ? 0.f : 0.10f * (x8 / 4.f);
    }
    dx *= (1.0f - 0.20f * wHE_phi); // de‑emphasize ripple at high E
    xc = wrap_phi(x_corr - dx);
  }

  // ============================================================
  // η side — kept + PATCH: power‑mean for ±η at high E
  // ============================================================
    const float dEtF[3] = { g.dX[1], g.dY[1], g.dZ[1] };


  const float wEt_front_perp = proj_perp_len(aHat, dEtF[0], dEtF[1], dEtF[2]);

  // +η and −η face points
  const float cEtP[3] = { g.Xcenter + dEtF[0], g.Ycenter + dEtF[1], g.Zcenter + dEtF[2] };
  const float cEtM[3] = { g.Xcenter - dEtF[0], g.Ycenter - dEtF[1], g.Zcenter - dEtF[2] };

  float rEtUp[3], rEtDn[3];
  CorrectShowerDepth(ixHome, iyP,  Energy, cEtP[0], cEtP[1], cEtP[2], rEtUp[0], rEtUp[1], rEtUp[2]);
  CorrectShowerDepth(ixHome, iyDn, Energy, cEtM[0], cEtM[1], cEtM[2], rEtDn[0], rEtDn[1], rEtDn[2]);

  const float dEtUpVec[3] = { rEtUp[0]-r0E[0], rEtUp[1]-r0E[1], rEtUp[2]-r0E[2] };
  const float dEtDnVec[3] = { rEtDn[0]-r0E[0], rEtDn[1]-r0E[1], rEtDn[2]-r0E[2] };

  const float wEt_depth_perpU = proj_perp_len(aHat, dEtUpVec[0], dEtUpVec[1], dEtUpVec[2]);
  const float wEt_depth_perpD = proj_perp_len(aHat, dEtDnVec[0], dEtDnVec[1], dEtDnVec[2]);

  // arithmetic mean (baseline)
  float wEt_depth_perp_sym = 0.f;
  if (wEt_depth_perpU > 0.f && wEt_depth_perpD > 0.f)
    wEt_depth_perp_sym = 0.5f * (wEt_depth_perpU + wEt_depth_perpD);
  else
    wEt_depth_perp_sym = (wEt_depth_perpU > 0.f ? wEt_depth_perpU :
                          (wEt_depth_perpD > 0.f ? wEt_depth_perpD : 0.f));

  // PATCH: power‑mean (p=0.5) to de‑emphasize outliers at high E
  float wHE_eta = 0.0f;
  constexpr float E_ETA_AVG_LO = 12.0f;
  constexpr float E_ETA_AVG_HI = 25.0f;
  if (Energy > E_ETA_AVG_LO) {
    const float t = (Energy - E_ETA_AVG_LO) / (E_ETA_AVG_HI - E_ETA_AVG_LO);
    wHE_eta = (t < 0.f) ? 0.f : (t > 1.f) ? 1.f : t;
  }
  if (wEt_depth_perpU > 0.f && wEt_depth_perpD > 0.f) {
    const float p = 0.5f;
    const float pm = 0.5f * (std::pow(wEt_depth_perpU, p) + std::pow(wEt_depth_perpD, p));
    const float w_pm = std::pow(pm, 1.0f/p);
    const float beta_pm = 0.35f * wHE_eta;  // modest influence
    wEt_depth_perp_sym = (1.0f - beta_pm) * wEt_depth_perp_sym + beta_pm * w_pm;
  }

  const float sEta_plus  = (wEt_depth_perpU   > 1e-6f) ? (wEt_front_perp / wEt_depth_perpU)   : 1.0f;
  const float sEta_sym   = (wEt_depth_perp_sym > 1e-6f) ? (wEt_front_perp / wEt_depth_perp_sym) : 1.0f;

  // main η scale: lowE→ +η only, highE→ symmetric ±η (kept)
  const float sEta_blend = (1.0f - wHE_eta) * sEta_plus + wHE_eta * sEta_sym;

  // tiny along‑eη component at high E (kept)
  const float wEt_front_along = std::fabs(dot3(dEtF, eeta));
  float wEt_depth_alongU = std::fabs(dot3(dEtUpVec, eeta));
  float wEt_depth_alongD = std::fabs(dot3(dEtDnVec, eeta));

  float wEt_depth_along_sym = 0.f;
  if (wEt_depth_alongU > 0.f && wEt_depth_alongD > 0.f) {
    // PATCH: power‑mean for along as well (same p), small weight
    float w_al_sym = 0.5f * (wEt_depth_alongU + wEt_depth_alongD);
    const float p = 0.5f;
    const float pm = 0.5f * (std::pow(wEt_depth_alongU, p) + std::pow(wEt_depth_alongD, p));
    const float w_pm = std::pow(pm, 1.0f/p);
    const float beta_pm_al = 0.25f * wHE_eta; // slightly smaller than ⟂ case
    wEt_depth_along_sym = (1.0f - beta_pm_al) * w_al_sym + beta_pm_al * w_pm;
  } else {
    wEt_depth_along_sym = (wEt_depth_alongU > 0.f ? wEt_depth_alongU :
                           (wEt_depth_alongD > 0.f ? wEt_depth_alongD : 0.f));
  }

  const float sEta_along = (wEt_depth_along_sym > 1e-6f) ? (wEt_front_along / wEt_depth_along_sym) : 1.0f;

  // tiny push toward along‑eη at high E (kept)
  const float alpha_eta_along = 0.15f * wHE_eta;
  const float sEta_final = (1.0f - alpha_eta_along) * sEta_blend + alpha_eta_along * sEta_along;

  float by = b0 * sEta_final;

  // tiny high‑E depth tempering for η (kept)
  if (depthE > 1e-3f && depthRef > 1e-3f) {
    const float ratio = clampf(depthRef / depthE, 0.85f, 1.15f);
    const float LAMBDA_HE = 0.10f;
    const float boost = std::pow(ratio, LAMBDA_HE * wHE_temper);
    by *= boost;
  }

  // Existing high‑E clamp tightening (kept)
  constexpr float E_ETA_TIGHT_LO = 10.0f;
  constexpr float E_ETA_TIGHT_HI = 25.0f;
  float wHigh = 0.0f;
  if (Energy <= E_ETA_TIGHT_LO)      wHigh = 0.0f;
  else if (Energy >= E_ETA_TIGHT_HI) wHigh = 1.0f;
  else                                wHigh = (Energy - E_ETA_TIGHT_LO) / (E_ETA_TIGHT_HI - E_ETA_TIGHT_LO);

  const float by_lo = (1.0f - wHigh)*0.11f + wHigh*0.12f;
  const float by_hi = (1.0f - wHigh)*0.30f + wHigh*0.25f;
  by = clampf(by, by_lo, by_hi);

  // report bη now (final)
  if (out_beta) *out_beta = by;

    // η inverse with asymmetric window at shower depth (±η), gently enabled at high E
    float y_corr = y;
    {
      const int   iy0   = EmcCluster::lowint(y_corr + 0.5f);
      const float Ymeas = y_corr - iy0; // in (-0.5, +0.5]

      if (std::fabs(Ymeas) <= 0.5f)
      {
        // symmetric closed-form (fallback & initializer)
        const float Sy    = std::sinh(0.5f / by);
        const float t_sym = by * std::asinh( 2.f * Ymeas * Sy );

        auto safe_ratio = [&](float num, float den)->float { return (den > 1e-6f) ? (num/den) : 1.0f; };

        // ⟂ to shower axis: +η(Up) and −η(Down) depths
        const float sU_perp = safe_ratio(wEt_front_perp, wEt_depth_perpU);
        const float sD_perp = safe_ratio(wEt_front_perp, wEt_depth_perpD);

        // along eη (tiny weight at high E only)
        const float sU_al   = safe_ratio(wEt_front_along, wEt_depth_alongU);
        const float sD_al   = safe_ratio(wEt_front_along, wEt_depth_alongD);

        const float sU_eff = (1.f - alpha_eta_along) * sU_perp + alpha_eta_along * sU_al;
        const float sD_eff = (1.f - alpha_eta_along) * sD_perp + alpha_eta_along * sD_al;

        float wUp = 0.5f * sU_eff;
        float wDn = 0.5f * sD_eff;

        wUp = clampf(wUp, 0.35f, 0.80f);
        wDn = clampf(wDn, 0.35f, 0.80f);

        // exact asymmetric map (same structure as φ, with Up/Down in place of R/L)
        auto Y_from_t_asym = [&](float t)->float {
          const float eU = std::exp(-(wUp - t)/by);
          const float eD = std::exp(-(wDn + t)/by);
          const float I0 = 1.f - 0.5f*(eU + eD);
          const float I1 = 0.5f*((wDn + t)*eD - (wUp - t)*eU) + 0.5f*by*(eD - eU);
          return t + I1 / std::max(I0, 1e-6f);
        };
        auto dYdt_num = [&](float t)->float {
          const float h = 1e-3f;
          return (Y_from_t_asym(t + h) - Y_from_t_asym(t - h)) / (2.f*h);
        };

        float t = t_sym;
        const float wASY = 0.75f * wHE_eta; // slightly stronger gate for η
        if (wASY > 0.f)
        {
          for (int it = 0; it < 2; ++it)
          {
            const float F  = Y_from_t_asym(t) - Ymeas;
            const float dF = dYdt_num(t);
            if (std::fabs(dF) < 1e-4f) break;
            t -= F / dF;
            t  = clampf(t, -wDn + 1e-4f, +wUp - 1e-4f);
          }
          t = (1.f - wASY) * t_sym + wASY * t;
        }

        y_corr = iy0 + t;
      }
    }
    yc = y_corr;
}







//void BEmcRecCEMC::CorrectPositionEnergyAware(float Energy, float x, float y,
//                                             float& xc, float& yc)
//{
//  // Legacy pass-through / gating
//  if (!m_UseCorrectPosition)
//  { xc = x; yc = y; return; }
//
//  if (!std::isfinite(Energy) || !std::isfinite(x) || !std::isfinite(y) || Energy < 0.01f)
//  { xc = x; yc = y; return; }
//
//  // ---------------- MC‑tuned energy laws (no |η| dep) ----------------
//  constexpr float E0    = 3.0f;        // reference energy for the log
//  constexpr float b0phi = 0.183339f;   // from your "no |η|-dep" φ fit
//  constexpr float mphi  = -0.007931f;
//
//  constexpr float b0eta = 0.191296f;   // from your "no |η|-dep" η fit
//  constexpr float meta  = -0.004557f;
//
//  const float lnE = std::log(std::max(Energy, 1e-6f) / E0);
//
//  auto clamp_b = [](float b)
//  {
//    // light numerical guard; keeps inverse‑asinh well‑behaved
//    return (b < 0.10f) ? 0.10f : (b > 0.30f ? 0.30f : b);
//  };
//
//  const float bx = clamp_b(b0phi + mphi * lnE);   // φ‑direction b
//  const float by = clamp_b(b0eta + meta * lnE);   // η‑direction b
//
//  // ---------------- φ inverse-asinh (identical to legacy) -------------
//  float x0  = x;
//  int   ix0 = EmcCluster::lowint(x0 + 0.5f);
//  if (std::fabs(x0 - ix0) <= 0.5f)
//  {
//    x0 = ix0 + bx * std::asinh(2.f * (x0 - ix0) * std::sinh(0.5f / bx));
//  }
//  else
//  {
//    x0 = x;
//#ifndef NDEBUG
//    std::cout << "????? BEmcRecCEMC::CorrectPositionEnergyAware: x guard triggered; "
//                 "x=" << x << "  dx=" << (x0 - ix0) << std::endl;
//#endif
//  }
//
//  // ---- module-of-8 ripple (legacy)
//  // NOLINTNEXTLINE(bugprone-incorrect-roundings)
//  int   ix8 = int(x + 0.5f) / 8;
//  float x8  = x + 0.5f - (ix8 * 8) - 4.f;  // −4 … +4
//  float dx  = 0.f;
//  if (m_UseDetailedGeometry)
//  {
//    // NOLINTNEXTLINE(bugprone-incorrect-roundings)
//    int local_ix8 = int(x + 0.5f) - ix8 * 8;
//    dx = static_cast<float>(factor_[local_ix8]) * (x8 / 4.f);
//  }
//  else
//  {
//    dx = 0.10f * (x8 / 4.f);
//    if (std::fabs(x8) > 3.3f) dx = 0.f;  // Don’t correct near module edge
//  }
//
//  // ---- compose φ and wrap
//  xc = x0 - dx;
//  while (xc < -0.5f)       { xc += float(fNx); }
//  while (xc >= fNx - 0.5f) { xc -= float(fNx); }
//
//  // ---------------- η inverse-asinh (identical to legacy) -------------
//  float y0  = y;
//  int   iy0 = EmcCluster::lowint(y0 + 0.5f);
//  if (std::fabs(y0 - iy0) <= 0.5f)
//  {
//    y0 = iy0 + by * std::asinh(2.f * (y0 - iy0) * std::sinh(0.5f / by));
//  }
//  else
//  {
//    y0 = y;
//#ifndef NDEBUG
//    std::cout << "????? BEmcRecCEMC::CorrectPositionEnergyAware: y guard triggered; "
//                 "y=" << y << "  dy=" << (y0 - iy0) << std::endl;
//#endif
//  }
//  yc = y0;
//}
