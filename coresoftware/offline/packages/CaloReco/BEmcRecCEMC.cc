#include "BEmcRecCEMC.h"

#include "BEmcCluster.h"
#include "BEmcProfile.h"
#include "TVector3.h"

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
  // (a,b) in radians for φ(E) = a − b·lnE.  If you later calibrate EA flavours,
  // update the four CLUS_CP_EA_* entries below.
  static const std::unordered_map<ETiltVariant,std::pair<double,double>> lut = {
      {ETiltVariant::CLUS_RAW, { 2.641645e-03, 8.117525e-04 }},  // no corr, cluster
      {ETiltVariant::CLUS_CP, { 1.925290e-03, 7.915014e-04 }},  // CorrectPosition, cluster
      {ETiltVariant::CLUS_CP_EA_FIT_ZDEP, { 2.084768e-03, 8.199622e-04 }},  // CorrectPosition(EA |z| + E), cluster
      {ETiltVariant::CLUS_CP_EA_FIT_ETADEP, { 2.090636e-03, 8.208976e-04 }},  // CorrectPosition(EA |#eta|+E), cluster
      {ETiltVariant::CLUS_CP_EA_FIT_EONLY, { 2.085796e-03, 8.203957e-04 }},  // CorrectPosition(EA E-only), cluster
      {ETiltVariant::CLUS_CP_EA_MIX, { 2.101579e-03, 8.261614e-04 }},  // CorrectPosition(EA #varphi:E-only, #eta:|#eta|+E), cluster
      {ETiltVariant::CLUS_BCORR, { 2.119625e-03, 8.130756e-04 }},  // b(E) corr, cluster
      {ETiltVariant::PDC_RAW, { 2.739621e-03, 8.204782e-04 }},  // no corr, scratch
      {ETiltVariant::PDC_CORR, { 2.244311e-03, 8.253109e-04 }},  // b(E) corr, scratch


      
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
    bool doPhiTilt = true; // <<< set to false internally to disable
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

    // Bin order (9 bins total):
    //   0:[0,2)  1:[2,4)  2:[4,6)  3:[6,8)  4:[8,10)  5:[10,20)  6:[20,30)  7:[30,45)  8:[45,60] and above
    auto zIndex = [](float z_cm)->int {
      // |z| is non-negative; return a valid bin for all inputs (≥60 → last bin)
      if (z_cm <  2.f)  return 0;
      if (z_cm <  4.f)  return 1;
      if (z_cm <  6.f)  return 2;
      if (z_cm <  8.f)  return 3;
      if (z_cm < 10.f)  return 4;
      if (z_cm < 20.f)  return 5;
      if (z_cm < 30.f)  return 6;
      if (z_cm < 45.f)  return 7;
      // [45,60) → 8, and anything ≥60 cm also → 8 (clamp high)
      return 8;
    };

    int iz = zIndex(absZ);

  // ---- per-|z| bin log-fit coefficients from your PDC-RAW tables (E0 = 3 GeV) ----
  // φ: bφ(E,|z|) = b0φ + mφ * ln(E/E0)
  // η: bη(E,|z|) = b0η + mη * ln(E/E0)
  //
  // Fine 0–10 cm (2-cm slices) from bValuesPhiOverlay_zFine / bValuesEtaOverlay_zFine:
  //   z00to02:  b0φ=0.181614, mφ=-0.006720   |  b0η=0.187871, mη=-0.002202
  //   z02to04:  b0φ=0.181808, mφ=-0.007140   |  b0η=0.187243, mη=-0.002529
  //   z04to06:  b0φ=0.183499, mφ=-0.007598   |  b0η=0.186277, mη=-0.001018
  //   z06to08:  b0φ=0.177695, mφ=-0.004814   |  b0η=0.189039, mη=-0.002569
  //   z08to10:  b0φ=0.179001, mφ=-0.005619   |  b0η=0.190574, mη=-0.004268
  //
  // Coarse 10–60 cm from bValuesPhiOverlay_zOnly / bValuesEtaOverlay_zOnly:
  //   z10to20:  b0φ=0.181252, mφ=-0.006485   |  b0η=0.192266, mη=-0.002334
  //   z20to30:  b0φ=0.179671, mφ=-0.006219   |  b0η=0.207996, mη=-0.006129
  //   z30to45:  b0φ=0.178555, mφ=-0.006275   |  b0η=0.235485, mη=-0.010342
  //   z45to60:  b0φ=0.176307, mφ=-0.006279   |  b0η=0.322516, mη=-0.002704
  //
  // For 0–10 coarse in that table, you are already using the fine bins above.

  static constexpr float E0 = 3.0f;

  // 9-bin arrays in the order listed above
  static constexpr float B0_PHI_Z[9] = {
    0.181614f, 0.181808f, 0.183499f, 0.177695f, 0.179001f,  // 0..4 : 0–10 cm fine
    0.181252f, 0.179671f, 0.178555f, 0.176307f               // 5..8 : 10–60 cm coarse
  };
  static constexpr float M_PHI_Z[9] = {
   -0.006720f,-0.007140f,-0.007598f,-0.004814f,-0.005619f,  // 0..4
   -0.006485f,-0.006219f,-0.006275f,-0.006279f               // 5..8
  };

  static constexpr float B0_ETA_Z[9] = {
    0.187871f, 0.187243f, 0.186277f, 0.189039f, 0.190574f,  // 0..4
    0.192266f, 0.207996f, 0.235485f, 0.322516f               // 5..8
  };
  static constexpr float M_ETA_Z[9] = {
   -0.002202f,-0.002529f,-0.001018f,-0.002569f,-0.004268f,  // 0..4
   -0.002334f,-0.006129f,-0.010342f,-0.002704f               // 5..8
  };

  // Safety: clamp iz in [0..8]
  if (iz < 0) iz = 0;
  if (iz > 8) iz = 8;

  // Build b(E) for φ/η at this z-slice
  auto clamp_b = [](float b){ return (b < 0.10f) ? 0.10f : (b > 0.30f ? 0.30f : b); };
  const float lnE = std::log(std::max(Energy, 1e-6f) / E0);

  const float bx = clamp_b(B0_PHI_Z[iz] + M_PHI_Z[iz] * lnE); // φ-direction b(E; z)
  const float by = clamp_b(B0_ETA_Z[iz] + M_ETA_Z[iz] * lnE); // η-direction b(E; z)

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
 WITH ETA DEPENDENCY AND ENERGY DEPENDENCY FROM PDC-RAW b VALUE FITS
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

    // ---- per-bin log-fit coefficients from PDC-RAW (E0 = 3 GeV) ----
    // Indexing convention:
    //   [0] |η| ≤ 0.20        (etaCore)
    //   [1] 0.20 < |η| ≤ 0.70 (etaMid)
    //   [2] 0.70 < |η| ≤ 1.10 (etaEdge)
    //   [3] fallback (originalEta; no-η-dep)
    constexpr float E0 = 3.0f;

    // φ: bφ(E,|η|) = b0φ + mφ * ln(E/E0)  (from bValuesPhiOverlay)
    static constexpr float B0_PHI[4] = {
      0.185809f,  // etaCore
      0.182105f,  // etaMid
      0.179244f,  // etaEdge
      0.180775f   // originalEta (fallback)
    };
    static constexpr float M_PHI[4] = {
     -0.006405f, // etaCore
     -0.006936f, // etaMid
     -0.006973f, // etaEdge
     -0.006402f  // originalEta (fallback)
    };

    // η: bη(E,|η|) = b0η + mη * ln(E/E0)  (from bValuesEtaOverlay)
    static constexpr float B0_ETA[4] = {
      0.177320f,  // etaCore
      0.194483f,  // etaMid
      0.196258f,  // etaEdge
      0.188228f   // originalEta (fallback)
    };
    static constexpr float M_ETA[4] = {
     -0.006088f, // etaCore
     -0.005003f, // etaMid
      0.002431f, // etaEdge (NOTE: positive slope)
     -0.002528f  // originalEta (fallback)
    };


  auto clamp_b = [](float b){ return (b < 0.10f) ? 0.10f : (b > 0.30f ? 0.30f : b); };
  const float lnE = std::log(std::max(Energy, 1e-6f) / E0);

  const float bx = clamp_b(B0_PHI[iEtaBin] + M_PHI[iEtaBin] * lnE);
  const float by = clamp_b(B0_ETA[iEtaBin] + M_ETA[iEtaBin] * lnE);

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



/*
 ENERGY ONLY NO ETA DEP FITS FROM PDC-RAW INFO
 */

void BEmcRecCEMC::CorrectPositionEnergyAwareEnergyDepOnly(float Energy, float x, float y,
                                             float& xc, float& yc)
{
  // ---- legacy gating -------------------------------------------------
  if (!m_UseCorrectPosition) { xc = x; yc = y; return; }
  if (!std::isfinite(Energy) || !std::isfinite(x) || !std::isfinite(y) || Energy < 0.01f)
  { xc = x; yc = y; return; }

    // ---- energy-only (no |η| dependence) log-fit laws (from PDC-RAW ‘originalEta’) ----
    // b(E) = b0 + m * ln(E/E0)  with  E0 = 3 GeV
    //  - PHI  (originalEta): b0 = 0.180775 , m = -0.006402
    //  - ETA  (originalEta): b0 = 0.188228 , m = -0.002528
    constexpr float E0     = 3.0f;
    constexpr float B0_PHI = 0.180775f;   // originalEta (PHI)
    constexpr float M_PHI  = -0.006402f;  // originalEta (PHI)
    constexpr float B0_ETA = 0.188228f;   // originalEta (ETA)
    constexpr float M_ETA  = -0.002528f;  // originalEta (ETA)

    const float lnE = std::log(std::max(Energy, 1e-6f) / E0);

    auto clamp_b = [](float b){ return (b < 0.10f) ? 0.10f : (b > 0.30f ? 0.30f : b); };

    // base (angle–zero) scales from energy-only laws
    const float bphi_E = clamp_b(B0_PHI + M_PHI * lnE);   // φ baseline
    const float beta_E = clamp_b(B0_ETA + M_ETA * lnE);   // η baseline

    // --- geometry-only incidence from detailed tower geometry -----------------------
    // use the same integer anchor as your local tower coordinates (no extra lookup later)
    const int ix_geom = EmcCluster::lowint(x + 0.5f);
    const int iy_geom = EmcCluster::lowint(y + 0.5f);

    // fetch per-tower geometry already filled by SetTowerGeometry(...)
    TowerGeom g{};
    bool haveGeom = GetTowerGeometry(ix_geom, iy_geom, g);   // your existing helper

    float bx = bphi_E;   // defaults (if geometry missing)
    float by = beta_E;

    if (haveGeom)
    {
      // local face tangents and depth (fiber) axis
      TVector3 ephi(g.dX[0], g.dY[0], g.dZ[0]); ephi = ephi.Unit();
      TVector3 eeta(g.dX[1], g.dY[1], g.dZ[1]); eeta = eeta.Unit();
      TVector3 n = (ephi.Cross(eeta)).Unit();

      // orient axis toward IP (projective SPACAL)
      TVector3 C(g.Xcenter, g.Ycenter, g.Zcenter);
      if (n.Dot(-C) < 0.0) n = -n;

      // ray from actual vertex to tower center (adequate for incidence)
      TVector3 V(0., 0., fVz);
      TVector3 p = (C - V).Unit();

      // component angles via atan2; use cosine for the geometric foreshortening
      const double pn   = p.Dot(n);
      const double pphi = p.Dot(ephi);
      const double peta = p.Dot(eeta);

      // cos(α_dir) = (p·n) / sqrt( (p·n)^2 + (p·e_dir)^2 )
      const double cos_alpha_phi = std::max(1e-6, std::abs(pn) / std::sqrt(pn*pn + pphi*pphi));
      const double cos_alpha_eta = std::max(1e-6, std::abs(pn) / std::sqrt(pn*pn + peta*peta));

      // record incident angles (radians) for downstream QA fills
      m_lastAlphaPhi = static_cast<float>(std::acos(std::min(1.0, cos_alpha_phi)));
      m_lastAlphaEta = static_cast<float>(std::acos(std::min(1.0, cos_alpha_eta)));

      // geometry-only correction: b_dir^eff = b_dir(E) * sec(α_dir) = b_dir(E) / cos(α_dir)
      bx = clamp_b( bphi_E / static_cast<float>(cos_alpha_phi) );
      by = clamp_b( beta_E / static_cast<float>(cos_alpha_eta) );
    }
    // -------------------------------------------------------------------------------
    // bx/by are now the angle-aware scales used by the inverse-asinh below


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



/*
 energy only in phi no eta dependence and eta dependence in eta
 */
void BEmcRecCEMC::CorrectPositionEnergyAwareEtaDepOnlyForEtaEnergyForPhi(float Energy, float x, float y,
                                             float& xc, float& yc)
{
  // ---- legacy gating -------------------------------------------------
  if (!m_UseCorrectPosition) { xc = x; yc = y; return; }
  if (!std::isfinite(Energy) || !std::isfinite(x) || !std::isfinite(y) || Energy < 0.01f)
  { xc = x; yc = y; return; }

  // ---- determine |η| slice from tower geometry (fallback → no-η-dep) ----
  const int ix0 = EmcCluster::lowint(x + 0.5f);
  const int iy0 = EmcCluster::lowint(y + 0.5f);

  TowerGeom g{};
  const bool haveGeom = GetTowerGeometry(ix0, iy0, g);

  float absEta = 0.f;
  if (haveGeom)
  {
    const float r = std::hypot(g.Xcenter, g.Ycenter);
    absEta = std::fabs((r > 0.f) ? std::asinh(g.Zcenter / r) : 0.f);
  }

  // |η| bins for η-side tuning:
  // 0: |η| ≤ 0.20, 1: (0.20,0.70], 2: (0.70,1.10], 3: fallback (no-η-dep law)
  int iEtaBin = 3;
  if (haveGeom)
  {
    if      (absEta <= 0.20f) iEtaBin = 0;
    else if (absEta <= 0.70f) iEtaBin = 1;
    else                      iEtaBin = 2;
  }

  // ---- log-fit coefficients from MC (E0 = 3 GeV) -------------------------
  constexpr float E0 = 3.0f;

    // φ uses the **no-|η|-dep** category ONLY
    constexpr float B0_PHI_NOETA = 0.183330f;
    constexpr float M_PHI_NOETA  = -0.007932f;

    // η uses the **|η|-dependent** categories (0,1,2), 3=fallback(no-η-dep)
    static constexpr float B0_ETA[4] = { 0.178946f, 0.196326f, 0.200883f, 0.191289f };
    static constexpr float M_ETA [4] = {-0.007106f,-0.006145f,-0.001358f,-0.004542f };


  const float lnE = std::log(std::max(Energy, 1e-6f) / E0);
  auto clamp_b = [](float b){ return (b < 0.10f) ? 0.10f : (b > 0.30f ? 0.30f : b); };

  const float bx = clamp_b(B0_PHI_NOETA + M_PHI_NOETA * lnE);               // φ: no-|η|-dep law
  const float by = clamp_b(B0_ETA[iEtaBin] + M_ETA[iEtaBin] * lnE);         // η: |η|-dep law

  // ============================== φ ===================================
  float x_corr = x;
  {
    const float X = x - ix0;                         // (-0.5, +0.5]
    if (std::fabs(X) <= 0.5f)
    {
      const float Sx = std::sinh(0.5f / bx);
      x_corr = ix0 + bx * std::asinh(2.f * X * Sx);
    }

    // module-of-8 ripple (identical to legacy CorrectPosition)
    // NOLINTNEXTLINE(bugprone-incorrect-roundings)
    int   ix8 = int(x + 0.5f) / 8;
    float x8  = x + 0.5f - (ix8 * 8) - 4.f;          // −4 … +4
    float dx  = 0.f;
    if (m_UseDetailedGeometry)
    {
      // NOLINTNEXTLINE(bugprone-incorrect-roundings)
      int local_ix8 = int(x + 0.5f) - ix8 * 8;       // 0..7
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
    const float Y = y - iy0;                           // (-0.5, +0.5]
    if (std::fabs(Y) <= 0.5f)
    {
      const float Sy = std::sinh(0.5f / by);
      y_corr = iy0 + by * std::asinh(2.f * Y * Sy);
    }
  }
  yc = y_corr;
}
