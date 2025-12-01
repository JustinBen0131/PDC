#pragma once
#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wvirtual-function-default"
#endif
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

////////////////////////////////////////////////////////////
// Standard includes
////////////////////////////////////////////////////////////
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>

#include <phool/recoConsts.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <ffamodules/CDBInterface.h>

#include <g4main/PHG4ParticleGenerator.h>      // particle gun (API: set_eta/phi/z/mom)
#include <TVector2.h>                          // phi wrapping helper
#include <TVector3.h>
#include <TCanvas.h>                           // ROOT canvas for PNG plots
#include <TGraphErrors.h>                      // ROOT graphs with error bars
#include <TStyle.h>                            // ROOT global style
#include <TSystem.h>                           // ROOT system utilities (mkdir)
#include <TLine.h>                             // ROOT line primitive (vertical dividers)
#include <TLatex.h>                            // ROOT text (dummy module labels)
#include <TLegend.h>                           // ROOT legend for block-color mapping
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

////////////////////////////////////////////////////////////
// Local / sPHENIX packages
////////////////////////////////////////////////////////////
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloReco/CaloGeomMapping.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloReco/BEmcRecCEMC.h"

#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawTowerGeom.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawTowerGeomv5.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawTowerGeomContainer.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawTowerDefs.h"

#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic pop
#endif

// --- Load only what we need; local libs last so they override core versions ---
R__LOAD_LIBRARY(libcdbobjects)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libffaobjects.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)
R__LOAD_LIBRARY(libcalibCaloEmc_pi0.so)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libcalo_io.so)    // CaloBase (local)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libcalo_reco.so)  // CaloReco (local)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libPDC.so)        // PDC (local)

#endif // ROOT>=6

////////////////////////////////////////////////////////////
// Main macro (no outputs; terminal prints only)
////////////////////////////////////////////////////////////
void Fun4all_ParticleGunTest(
    // geometry / conditions
    int         runNumber           = 24,          // TIMESTAMP / RUNNUMBER for CDB
    // scan controls
    const std::string& scanMode     = "ring",      // "ring" | "philine" | "etaZvGrid" | "none"
    int         ietaRing            = 48,          // used by scanMode=="ring" and "etaZvGrid"
    double      phiTarget           = 0.0,         // used by scanMode=="philine" (radians)
    // incidence controls
    double      EGeV                = 8.0,         // test energy used in ComputeIncidenceSD
    bool        sandboxNoTilt       = true,        // zero-incidence anchor mode for incidence
    int         incidenceDbgLevel   = 5,           // >=5: structured per-call debug from BEmcRecCEMC
    int         incidenceHardStopLevel = -1,       // -1=no stop; >=0 with dbg>=this → throw
    // execution controls
    bool        actuallyFireGeant   = false        // if true: run 1 event per tower using the gun
)
{
    std::cout << "============================================================\n"
              << " Fun4all_ParticleGunTest - geometry / incidence diagnostic\n"
              << "------------------------------------------------------------\n"
              << "  runNumber          = " << runNumber          << "\n"
              << "  scanMode           = " << scanMode           << "\n"
              << "  ietaRing           = " << ietaRing           << "\n"
              << "  phiTarget [rad]    = " << phiTarget          << "\n"
              << "  EGeV               = " << EGeV               << "\n"
              << "  incidenceDbgLevel  = " << incidenceDbgLevel  << "\n"
              << "  incidenceHardStop  = " << incidenceHardStopLevel << "\n"
              << "  actuallyFireGeant  = " << (actuallyFireGeant ? "true" : "false") << "\n"
              << "  NOTE: this macro does *not* read DSTs or loop over Run-24 events.\n"
              << "        Geometry is built once via CaloGeomMapping::Init/InitRun.\n"
              << "============================================================\n";


  // ────────────────────────────────────────────────────────────────────────────
  // 1) Fun4All + Conditions DB (no input managers, no event loop)
  // ────────────────────────────────────────────────────────────────────────────
  Fun4AllServer* se = Fun4AllServer::instance();
  if (!se)
  {
    std::cerr << "[FATAL] Fun4AllServer::instance() returned nullptr.\n";
    return;
  }

  const int fun4allVerbosity = 1;   // 0=silent, 1=summary, 2+=chatty
  se->Verbosity(fun4allVerbosity);
  std::cout << "[INFO] Fun4AllServer created, verbosity=" << fun4allVerbosity << "\n";

  auto* rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  rc->set_uint64Flag("TIMESTAMP",     runNumber);
  rc->set_IntFlag   ("RUNNUMBER",     runNumber);
  std::cout << "[INFO] recoConsts set: GLOBALTAG=MDC2, TIMESTAMP=" << runNumber
            << ", RUNNUMBER=" << runNumber << "\n";

  // Sanity check: EMCal calibration CDB URL
  {
    const std::string url = CDBInterface::instance()->getUrl("EMCTOWERCALIB");
    std::cout << "[CDB] EMCTOWERCALIB URL = " << url << "\n";
  }

  // ────────────────────────────────────────────────────────────────────────────
  // 2) Build CEMC detailed geometry *without* running any events
  //    We call CaloGeomMapping::Init/InitRun directly, not Fun4AllServer::run.
  // ────────────────────────────────────────────────────────────────────────────
  auto* geomMap = new CaloGeomMapping("CEMC_GeomFiller");
  geomMap->set_detector_name("CEMC");
  geomMap->set_UseDetailedGeometry(true);   // 8-vertex blocks (RawTowerGeomv5)
  const int geomVerbosity = 0;              // 0=quiet, 1=summary, 2+=per-tower
  geomMap->Verbosity(geomVerbosity);

  // Register subsystem (mostly for consistency; we won't use se->run here)
  se->registerSubsystem(geomMap);
  std::cout << "[INIT] Registered CaloGeomMapping(\"CEMC_GeomFiller\") as a subsystem.\n";

  PHCompositeNode* topNode = se->topNode();
  if (!topNode)
  {
    std::cerr << "[FATAL] topNode is null; Fun4AllServer node tree not initialized.\n";
    return;
  }

  std::cout << "[INIT] Manually calling CaloGeomMapping::Init/InitRun (no event loop)...\n";
  int ierr_init = geomMap->Init(topNode);
  int ierr_run  = geomMap->InitRun(topNode);
  std::cout << "[INIT] CaloGeomMapping::Init returned "  << ierr_init
            << ", InitRun returned " << ierr_run << "\n";

  if (ierr_init != 0 || ierr_run != 0)
  {
    std::cerr << "[WARN] CaloGeomMapping init returned non-zero codes; geometry may be incomplete.\n";
  }

  // ────────────────────────────────────────────────────────────────────────────
  // 3) Fetch geometry container from the node tree
  // ────────────────────────────────────────────────────────────────────────────
  auto* geo =
    findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC_DETAILED");
  if (!geo)
  {
    std::cout << "[WARN] TOWERGEOM_CEMC_DETAILED not found, trying TOWERGEOM_CEMC instead.\n";
    geo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  }
  if (!geo)
  {
    std::cerr << "[FATAL] CEMC geometry (detailed or default) not found on node tree.\n";
    return;
  }

  const int nPhi = geo->get_phibins();
  const int nEta = geo->get_etabins();
  std::cout << "[GEOM] CEMC geometry ready: phibins=" << nPhi
            << "  etabins=" << nEta << "\n";

    // ────────────────────────────────────────────────────────────────────────────
    // 4) Configure BEmcRecCEMC with detailed tower geometry (one pass)
    // ────────────────────────────────────────────────────────────────────────────
    auto* bemc = new BEmcRecCEMC();
    bemc->SetCylindricalGeometry();
    bemc->set_UseDetailedGeometry(true);
    bemc->SetTowerThreshold(0.030f); // keep consistent with your cluster builder

    // Flag to track whether any mechanical-incidence calls were actually made
    // in this macro (used to decide whether to print the MECH_QA summary).
    bool incidenceEverCalled = false;

    // Mechanical-incidence QA: reset accumulators before the scan.
    // All scan modes below (ring / phiScanFromCenter / philine / etaZvGrid)
    // call CalculateMechIncidence(...); ResetMechIncidenceQA() clears the
    // static per-process counters that are summarised at the end by
    // PrintMechIncidenceQASummary().
    bemc->ResetMechIncidenceQA();



  // Copy dims then feed towers; note: ix=phi, iy=eta convention in BEmcRecCEMC
  bemc->SetDim(nPhi, nEta);

  std::cout << "[GEOM] Loading tower geometries into BEmcRecCEMC...\n";
  int nLoaded = 0;
  for (int ieta = 0; ieta < nEta; ++ieta)
  {
    for (int iphi = 0; iphi < nPhi; ++iphi)
    {
      const RawTowerDefs::keytype key =
        RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi);
      if (auto* tg = geo->get_tower_geometry(key))
      {
        bemc->SetTowerGeometry(iphi, ieta, *tg);
        ++nLoaded;
      }
    }
  }
  std::cout << "[GEOM] Loaded " << nLoaded << " towers into BEmcRecCEMC.\n";

  bemc->CompleteTowerGeometry();   // derive per-face tangents once
  std::cout << "[GEOM] CompleteTowerGeometry() done (face tangents cached).\n";

    // Incidence tracer controls (local to BEmcRecCEMC)
    // In the MECH-only implementation, we no longer support a "no-tilt sandbox"
    // switch inside BEmcRecCEMC. We keep sandboxNoTilt only as a label in the
    // macro output, but incidence is always computed from the encoded rotX/Y/Z
    // plus detailed geometry.
    bemc->SetIncidenceDebugLevel(incidenceDbgLevel);
    bemc->SetIncidenceHardStopLevel(incidenceHardStopLevel);

    std::cout << "[INC] MECH-only incidence\n"
              << "      IncidenceDebugLevel    = " << incidenceDbgLevel      << "\n"
              << "      IncidenceHardStopLevel = " << incidenceHardStopLevel << "\n";


  // ────────────────────────────────────────────────────────────────────────────
  // 5) Minimal gun (used only if actuallyFireGeant==true)
  //    PHG4ParticleGenerator supports set_eta/phi/mom/vtx.
  //    *** IMPORTANT: we never call se->run() here unless actuallyFireGeant=true,
  //    so no events are processed by default. ***
  // ────────────────────────────────────────────────────────────────────────────
  auto* gun = new PHG4ParticleGenerator("GUN");
  gun->set_name("gamma");
  gun->set_mom_range(EGeV, EGeV);
  gun->set_eta_range(0.0, 0.0);              // radial shot → αφ≈0 for a perfectly projective tower

  se->registerSubsystem(gun);
  std::cout << "[GUN] PHG4ParticleGenerator registered (name=\"gamma\", E=" << EGeV << " GeV).\n"
            << "      No Fun4All event loop will run unless actuallyFireGeant=true.\n";

  // ────────────────────────────────────────────────────────────────────────────
  // 6) Helper: aim/compute (and optionally "fire") at a specific tower center
  // ────────────────────────────────────────────────────────────────────────────
  auto compute_or_fire_at = [&](int ieta, int iphi)
  {
      const RawTowerDefs::keytype key =
          RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi);
      RawTowerGeom* tg = geo->get_tower_geometry(key);
      if (!tg)
      {
          std::cout << "[WARN] Missing geometry for (ieta="<<ieta<<", iphi="<<iphi<<") - skipped.\n";
          return;
      }

      // ------------------------------------------------------------------
      // Tower-center geometry in detector coordinates
      //   Cx, Cy, Cz : tower center in cm
      //   Rc        : cylindrical radius = sqrt(Cx^2 + Cy^2)
      //   phiC      : tower-center azimuth
      //
      // Detector pseudorapidity (eta_det_IP) is defined as seen from the IP:
      //   eta_det_IP = asinh( Cz / Rc )
      //
      // Source-detector pseudorapidity (eta_SD) uses the chosen vertex z_vtx:
      //   eta_SD = asinh( (Cz - z_vtx) / Rc )
      // For this QA we set z_vtx = 0 ⇒ eta_SD == eta_det_IP.
      // ------------------------------------------------------------------
      const double Cx  = tg->get_center_int_x();
      const double Cy  = tg->get_center_int_y();
      const double Cz  = tg->get_center_z();

      const double Rc  = std::hypot(Cx, Cy);            // detector radius of tower center
      const double phiC = std::atan2(Cy, Cx);           // detector φ of tower center

      const double z_vtx = 0.0;                         // vertex at IP for this QA
      const double eta_det_IP = (Rc > 0.0) ? std::asinh(Cz / Rc)           : 0.0;
      const double eta_SD     = (Rc > 0.0) ? std::asinh((Cz - z_vtx) / Rc) : 0.0;

      // Tell BEmcRecCEMC to use this vertex z for its ray construction
      bemc->SetVertexZ(static_cast<float>(z_vtx));

      if (actuallyFireGeant)
      {
          std::cout << "[GUN] Firing one gamma through tower center "
                    << "(ieta="<<ieta<<", iphi="<<iphi<<") "
                    << "at z_vtx="<<z_vtx<<" phi="<<phiC<<"...\n";
          gun->set_vtx(0.0, 0.0, z_vtx);
          gun->set_phi_range(phiC, phiC);
          // Fun4AllServer::instance()->run(1);  // optional if you really want to simulate
      }

      // Compute front-face mechanical incidence at the tower center (x=iphi, y=ieta)
      float cosphi = 0.f, coseta = 0.f, aphis = 0.f, aetas = 0.f;
      incidenceEverCalled = true;
      const bool ok = bemc->CalculateMechIncidence(
              static_cast<float>(EGeV),
              static_cast<float>(iphi), static_cast<float>(ieta),
              cosphi, coseta, aphis, aetas);


          // Per-tower [SCAN] output only when debug is enabled
          if (incidenceDbgLevel > 0)
          {
            // ANSI helpers (local to this lambda)
            const char* ANSI_BOLD   = "\033[1m";
            const char* ANSI_GREEN  = "\033[32m";
            const char* ANSI_YELLOW = "\033[33m";
            const char* ANSI_RED    = "\033[31m";
            const char* ANSI_RESET  = "\033[0m";

            // Mechanical tilts (if detailed geometry is v5)
            double rotX = 0.0, rotY = 0.0, rotZ = 0.0;
            if (auto* g5 = dynamic_cast<RawTowerGeomv5*>(tg))
            {
              rotX = g5->get_rotx();
              rotY = g5->get_roty();
              rotZ = g5->get_rotz();
            }

            const double alpha_mag = std::sqrt(aphis*aphis + aetas*aetas);
            const double tol_alpha = 5.0e-3; // same ~0.29° scale used in MECH QA

            const char* statusColor = ok ? ANSI_GREEN : ANSI_RED;
            const char* statusText  = ok ? "OK"       : "FAIL";

            const char* tiltColor = (alpha_mag < tol_alpha)
                                    ? ANSI_GREEN
                                    : (alpha_mag < 3.0*tol_alpha ? ANSI_YELLOW : ANSI_RED);
            const char* tiltTag   = (alpha_mag < tol_alpha)
                                    ? "within_tol"
                                    : (alpha_mag < 3.0*tol_alpha ? "moderate_tilt" : "large_tilt");

            std::cout << "\n"
                      << ANSI_BOLD << "[SCAN] tower (ieta=" << ieta
                      << ", iphi=" << iphi << ") "
                      << statusColor << statusText << ANSI_RESET << "\n";

            std::cout << "  GEOM:  C=(" << Cx << ", " << Cy << ", " << Cz << ")"
                      << "  R_cyl=" << Rc
                      << "  phi_center=" << phiC
                      << "  eta_det_IP=" << eta_det_IP
                      << "  eta_SD=" << eta_SD
                      << "  z_vtx=" << z_vtx << "\n";

            std::cout << "  ROT :  rotX=" << rotX*1.0e3 << " mrad"
                      << "  rotY=" << rotY*1.0e3 << " mrad"
                      << "  rotZ=" << rotZ*1.0e3 << " mrad\n";

            std::cout << "  INC :  a_phi=" << aphis
                      << "  a_eta=" << aetas
                      << "  |a|=" << alpha_mag
                      << "  cos_phi=" << cosphi
                      << "  cos_eta=" << coseta << "\n";

            std::cout << "  TAGS:  " << tiltColor << tiltTag << ANSI_RESET
                        << "\n";

          }

      };


    // Helper: find iphi whose tower center φ is closest to a target φ
    auto pick_iphi_closest_to = [&](int ieta, double targetPhi)->int
    {
      int best_iphi = 0;
      double best_d = 1e9;
      for (int iphi = 0; iphi < nPhi; ++iphi)
      {
        const RawTowerDefs::keytype key =
          RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi);
        RawTowerGeom* tg = geo->get_tower_geometry(key);
        if (!tg) continue;
        const double phiC = std::atan2(tg->get_center_int_y(), tg->get_center_int_x());
        const double d = std::fabs(TVector2::Phi_mpi_pi(phiC - targetPhi));
        if (d < best_d)
        {
          best_d   = d;
          best_iphi = iphi;
        }
      }
      return best_iphi;
    };

    // ────────────────────────────────────────────────────────────────────────────
    // 7) Scan logic (no clusterizer, no PDC, no DSTs, no event loop)
    //     + helper to map (ieta, iphi) → (sector, module, block, tower-in-block)
    // ────────────────────────────────────────────────────────────────────────────
    struct TowerHardwareIndex
    {
      int sector;        // 0..63 global sector ID
      bool isNorth;      // true = North (ieta >= nEta/2), false = South
      int phiSectorIdx;  // 0..31 sector index within an endcap
      int ietaLocal;     // 0..47 within hemisphere
      int iphiLocal;     // 0..7 within sector
      int module;        // 0..23 module index along η within sector
      int etaInModule;   // 0 or 1 (tower position inside module along η)
      int blockPhi;      // 0..3 block index in φ within module
      int phiInBlock;    // 0 or 1 (tower position inside block along φ)
      int blockInSector; // 0..95 block index within sector
      int towerInBlock;  // 0..3 tower index within 2×2 block
    };

    auto mapTowerToHardwareIndex = [&](int ieta, int iphi) -> TowerHardwareIndex
    {
      TowerHardwareIndex idx{};

      // These constants must match the EMCal segmentation.
      const int nSectorsTotal      = 64;                 // 2 (±η) × 32 (φ)
      const int nPhiSectors        = nSectorsTotal / 2;  // 32 φ-sectors per endcap
      const int phiTowersPerSector = nPhi / nPhiSectors; // 8 φ-towers per sector
      const int etaSplitIndex      = nEta / 2;           // 48: North/South boundary

      // Hemisphere and sector index
      idx.isNorth      = (ieta >= etaSplitIndex);
      idx.phiSectorIdx = iphi / phiTowersPerSector;                     // 0..31
      idx.sector       = (idx.isNorth ? 0 : nPhiSectors) + idx.phiSectorIdx; // 0..63

      // Local coordinates within hemisphere / sector
      idx.ietaLocal = idx.isNorth ? (ieta - etaSplitIndex) : ieta;      // 0..47
      idx.iphiLocal = iphi % phiTowersPerSector;                        // 0..7

      // Modules along η: 24 slices × 2 towers each
      idx.module      = idx.ietaLocal / 2;  // 0..23
      idx.etaInModule = idx.ietaLocal % 2;  // 0 or 1

      // Blocks inside a module along φ: 4 blocks × 2 towers each
      idx.blockPhi   = idx.iphiLocal / 2;   // 0..3
      idx.phiInBlock = idx.iphiLocal % 2;   // 0 or 1

      // Convenience rolled-up indices
      idx.blockInSector = idx.module * 4 + idx.blockPhi;      // 0..95
      idx.towerInBlock  = idx.etaInModule * 2 + idx.phiInBlock; // 0..3

      return idx;
    };

      if (scanMode == "ring")
      {
        const int ieta = (ietaRing >= 0 && ietaRing < nEta)
                           ? ietaRing
                           : nEta/2;
        std::cout << "[MODE] ring - scanning all iphi at fixed ieta=" << ieta << "\n";
        for (int iphi = 0; iphi < nPhi; ++iphi)
          compute_or_fire_at(ieta, iphi);
      }
      else if (scanMode == "phiScanFromCenter")
      {
        // New mode:
        //   • auto-pick the η-ring whose tower center η is closest to 0 (within |η|<1.1),
        //   • pick iphi_zero so that tower center φ is closest to 0 rad,
        //   • scan φIndex = 0..nPhi-1 in order, wrapping iphi = (iphi_zero + φIndex) % nPhi,
        //   • for each "event", print a bold red header with:
        //        event N ---- #eta_center≈0 (ieta=...)  #phiTower=... (phiIndex=...)
        //        phi_center, phi_nominal(φIndex), and dphi_index = phi_center - phi_nominal
        //   • then print source vertex and a compact incidence summary.
        const char* ANSI_RED_BOLD = "\033[1;31m";
        const char* ANSI_RESET    = "\033[0m";

        // Helper: π for nominal φIndex → φ mapping
        const double PI    = std::acos(-1.0);
        const double TWO_PI = 2.0 * PI;

        // 1) Find central η ring in GLOBAL coordinates (η_det at IP closest to 0)
        int    ieta_central = -1;
        double bestAbsEta   = 1.0e9;

        for (int ieta = 0; ieta < nEta; ++ieta)
        {
          // probe tower at φ≈0 for this η to estimate η_det_IP
          const int iphi_probe = pick_iphi_closest_to(ieta, 0.0);
          const RawTowerDefs::keytype key_probe =
            RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi_probe);
          RawTowerGeom* tg_probe = geo->get_tower_geometry(key_probe);
          if (!tg_probe) continue;

          const double Cx  = tg_probe->get_center_int_x();
          const double Cy  = tg_probe->get_center_int_y();
          const double Cz  = tg_probe->get_center_int_z();
          const double Rc  = std::hypot(Cx, Cy);
          if (Rc <= 0.0) continue;

          const double eta_det = std::asinh(Cz / Rc);
          if (std::fabs(eta_det) > 1.1) continue;  // stay within EMCal acceptance

          const double aeta = std::fabs(eta_det);
          if (aeta < bestAbsEta)
          {
            bestAbsEta   = aeta;
            ieta_central = ieta;
          }
        }

        if (ieta_central < 0)
        {
          // Fallback: use user-specified ring or barrel center
          ieta_central = (ietaRing >= 0 && ietaRing < nEta)
                           ? ietaRing
                           : nEta/2;
          std::cout << "[phiScanFromCenter] WARNING: auto-selection of η≈0 failed; "
                    << "falling back to ieta=" << ieta_central << "\n";
        }

        // 2) On that η ring, find iphi whose tower-center φ is closest to 0 rad
        const int iphi_zero = pick_iphi_closest_to(ieta_central, 0.0);

        // Print anchor summary for this φ-scan
        {
          const RawTowerDefs::keytype key_ref =
            RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta_central, iphi_zero);
          RawTowerGeom* tg_ref = geo->get_tower_geometry(key_ref);
          if (tg_ref)
          {
            const double Cx0  = tg_ref->get_center_int_x();
            const double Cy0  = tg_ref->get_center_int_y();
            const double Cz0  = tg_ref->get_center_int_z();
            const double Rc0  = std::hypot(Cx0, Cy0);
            const double phi0 = std::atan2(Cy0, Cx0);
            const double eta0 = (Rc0 > 0.0) ? std::asinh(Cz0 / Rc0) : 0.0;

            std::cout << ANSI_RED_BOLD
                      << "\n[phiScanFromCenter] central (η,φ) anchor selection\n"
                      << ANSI_RESET;
            std::cout << "  chosen ieta_central = " << ieta_central
                      << "  (|eta_det_IP| ≈ " << std::fabs(eta0) << ")\n";
            std::cout << "  iphi_zero (phiIndex=0) = " << iphi_zero
                      << "  phi_center≈" << phi0 << " rad\n";
            std::cout << "  C_ref=(" << Cx0 << "," << Cy0 << "," << Cz0 << ")  "
                      << "R_cyl≈" << Rc0 << " cm\n\n";
          }
        }

        // 3) Vertex is at the IP, z = 0, for this diagnostic
        const double z_vtx = 0.0;
        bemc->SetVertexZ(static_cast<float>(z_vtx));

        const double dphi_step = TWO_PI / static_cast<double>(nPhi);

        std::cout << "[MODE] phiScanFromCenter - scanning φ around global (η≈0, φ≈0)"
                  << " with nPhi=" << nPhi << " towers\n";

        // 4) Scan φIndex=0…nPhi−1, wrapping iphi around the barrel
        for (int phiIndex = 0; phiIndex < nPhi; ++phiIndex)
        {
          const int iphi = (iphi_zero + phiIndex) % nPhi;

          const RawTowerDefs::keytype key =
            RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta_central, iphi);
          RawTowerGeom* tg = geo->get_tower_geometry(key);
          if (!tg)
          {
            std::cout << "[phiScanFromCenter] WARNING: missing geometry for (ieta="
                      << ieta_central << ", iphi=" << iphi << ") - skipped.\n";
            continue;
          }

          // Tower center geometry
          const double Cx  = tg->get_center_int_x();
          const double Cy  = tg->get_center_int_y();
          const double Cz  = tg->get_center_int_z();
          const double Rc  = std::hypot(Cx, Cy);
          const double phi_center = std::atan2(Cy, Cx);
          const double eta_det_IP = (Rc > 0.0) ? std::asinh(Cz / Rc) : 0.0;

          // Nominal "index-based" φ for this phiIndex and difference to true φ_center
          const double phi_nominal = TVector2::Phi_mpi_pi(phiIndex * dphi_step);
          const double dphi_index  = TVector2::Phi_mpi_pi(phi_center - phi_nominal);

          // Optionally: actually fire a Geant photon from IP towards this tower center
          if (actuallyFireGeant)
          {
            gun->set_vtx(0.0, 0.0, z_vtx);
            gun->set_phi_range(phi_center, phi_center);
            // Fun4AllServer::instance()->run(1);  // enable if you really want to simulate
          }

            // Compute front-face mechanical incidence using BEmcRecCEMC
            // (tower coords x=iphi, y=ieta_central)
            float cosphi = 0.f, coseta = 0.f, aphis = 0.f, aetas = 0.f;
            incidenceEverCalled = true;
            const bool ok_inc = bemc->CalculateMechIncidence(
                static_cast<float>(EGeV),
                static_cast<float>(iphi), static_cast<float>(ieta_central),
                cosphi, coseta, aphis, aetas);



            // ------------------------------------------------------------------
            // Pretty, human-readable per-"event" block (ANSI colored, tabulated)
            // ------------------------------------------------------------------
            const char* ANSI_GREEN  = "\033[32m";
            const char* ANSI_YELLOW = "\033[33m";

            // Mechanical tilts (in mrad) if available
            double rotX = 0.0, rotY = 0.0, rotZ = 0.0;
            if (auto* g5 = dynamic_cast<RawTowerGeomv5*>(tg))
            {
              rotX = g5->get_rotx();
              rotY = g5->get_roty();
              rotZ = g5->get_rotz();
            }

            const double alpha_mag = std::sqrt(aphis*aphis + aetas*aetas);
            const double tol_alpha = 5.0e-3;  // same ~0.29° scale as MECH QA

            const char* statusColor = ok_inc ? ANSI_GREEN : ANSI_RED_BOLD;
            const char* statusText  = ok_inc ? "OK"       : "FAIL";

            const char* tiltColor = (alpha_mag < tol_alpha)
                                    ? ANSI_GREEN
                                    : (alpha_mag < 3.0*tol_alpha ? ANSI_YELLOW : ANSI_RED_BOLD);
            const char* tiltTag   = (alpha_mag < tol_alpha)
                                    ? "within_tol"
                                    : (alpha_mag < 3.0*tol_alpha ? "moderate_tilt" : "large_tilt");

            std::cout << ANSI_RED_BOLD
                      << "event " << phiIndex
                      << "  tower(ieta=" << ieta_central
                      << ", iphi=" << iphi << ")  "
                      << "phiIndex=" << phiIndex << "  "
                      << statusColor << statusText << ANSI_RESET << "\n";

            std::cout << "  ANGLES:  phi_center=" << phi_center
                      << "  phi_nominal=" << phi_nominal
                      << "  dphi(center-nom)=" << dphi_index << " rad"
                      << "  eta_det_IP=" << eta_det_IP << "\n";

            std::cout << "  GEOM  :  C=(" << Cx << ", " << Cy << ", " << Cz << ")"
                      << "  R_cyl=" << Rc
                      << "  z_vtx=" << z_vtx << "  (src_vertex=(0,0," << z_vtx << "))\n";

            std::cout << "  ROT   :  rotX=" << rotX*1.0e3 << " mrad"
                      << "  rotY=" << rotY*1.0e3 << " mrad"
                      << "  rotZ=" << rotZ*1.0e3 << " mrad\n";

            std::cout << "  INC   :  a_phi=" << aphis
                      << "  a_eta=" << aetas
                      << "  |a|=" << alpha_mag
                      << "  cos_phi=" << cosphi
                      << "  cos_eta=" << coseta << "\n";

            std::cout << "  TAGS  :  " << tiltColor << tiltTag << ANSI_RESET
                      << "\n\n";

        }
      }
      else if (scanMode == "philine")
      {
        std::cout << "[MODE] philine - scanning all ieta along φ≈" << phiTarget << " rad\n";
        const double target = phiTarget;
        for (int ieta = 0; ieta < nEta; ++ieta)
        {
          const int best_iphi = pick_iphi_closest_to(ieta, target);
          compute_or_fire_at(ieta, best_iphi);
        }
      }
      else if (scanMode == "etaZvGrid")
      {
        std::cout << "[MODE] etaZvGrid - scanning a fixed φ column vs (ieta, z_vtx)\n";

        const double target = phiTarget;
        const int ieta_ref  = std::clamp(ietaRing, 0, nEta - 1);
        const int iphi_fixed = pick_iphi_closest_to(ieta_ref, target);

        // Reference tower info
        {
          const RawTowerDefs::keytype key_ref =
            RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta_ref, iphi_fixed);
          RawTowerGeom* tg_ref = geo->get_tower_geometry(key_ref);
          if (tg_ref)
          {
            const double Cx  = tg_ref->get_center_int_x();
            const double Cy  = tg_ref->get_center_int_y();
            const double Cz  = tg_ref->get_center_int_z();
            const double Rc  = std::hypot(Cx, Cy);
            const double phi = std::atan2(Cy, Cx);
            const double eta_det_ref = (Rc > 0.0) ? std::asinh(Cz / Rc) : 0.0;

            std::cout << "# [etaZvGrid] fixed-phi column scan\n"
                      << "#  reference tower: (ieta_ref="<<ieta_ref
                      << ", iphi_fixed="<<iphi_fixed<<")\n"
                      << "#  C_ref=("<<Cx<<","<<Cy<<","<<Cz<<")  R_ref="<<Rc
                      << "  phi_center_ref="<<phi
                      << "  eta_det_ref="<<eta_det_ref << "\n";
          }
        }

        std::cout << "# Columns (SCAN):\n"
                  << "#  z_vtx_cm, ieta, eta_det, eta_SD, "
                  << "alpha_phi, alpha_eta, cos_phi, cos_eta, "
                  << "phiTilt_mrad, etaTilt_mrad, rotZ_mrad\n";

        const bool   scanZvertex = false;         // set true if you want z=0 and z=20
        const double zVerts[2]   = {0.0, 20.0};
        const int    nZ          = scanZvertex ? 2 : 1;

        for (int ieta = 0; ieta < nEta; ++ieta)
        {
          const RawTowerDefs::keytype key =
            RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi_fixed);
          RawTowerGeom* tg = geo->get_tower_geometry(key);
          if (!tg) continue;

          const double Cx  = tg->get_center_int_x();
          const double Cy  = tg->get_center_int_y();
          const double Cz  = tg->get_center_int_z();
          const double Rc  = std::hypot(Cx, Cy);

          const double eta_det = (Rc > 0.0) ? std::asinh(Cz / Rc) : 0.0;
          if (eta_det < -1.1 || eta_det > 1.1) continue;  // acceptance cut

          const float x_eval = static_cast<float>(iphi_fixed); // fixed φ column
          const float y_eval = static_cast<float>(ieta);       // this η row

          for (int iz = 0; iz < nZ; ++iz)
          {
            const double z_vtx = scanZvertex ? zVerts[iz] : 0.0;

              const double eta_SD = (Rc > 0.0) ? std::asinh((Cz - z_vtx) / Rc) : 0.0;

              // Geometric φ/η tilts of the tower axis w.r.t. the radial direction.
              double phiTilt_mrad = 0.0;
              double etaTilt_mrad = 0.0;
              double rotZ         = 0.0;

              if (auto* g5 = dynamic_cast<RawTowerGeomv5*>(tg))
              {
                // Tower axis from inner → outer face (full 3D).
                TVector3 Cint(g5->get_center_int_x(),
                              g5->get_center_int_y(),
                              g5->get_center_int_z());
                TVector3 Cext(g5->get_center_ext_x(),
                              g5->get_center_ext_y(),
                              g5->get_center_ext_z());
                TVector3 axis = Cext - Cint;

                // Radial direction in transverse plane (IP → inner-face center).
                TVector3 rT(Cint.X(), Cint.Y(), 0.0);
                TVector3 aT(axis.X(), axis.Y(), 0.0);

                if (rT.Mag2() > 0.0 && aT.Mag2() > 0.0 && axis.Mag2() > 0.0)
                {
                  TVector3 e_rT = rT.Unit();
                  TVector3 e_aT = aT.Unit();

                  const double phi_radial = std::atan2(e_rT.Y(), e_rT.X());
                  const double phi_axis   = std::atan2(e_aT.Y(), e_aT.X());
                  const double dphi       = TVector2::Phi_mpi_pi(phi_axis - phi_radial);

                  // φ tilt of tower axis relative to radial (in mrad).
                  phiTilt_mrad = dphi * 1.0e3;

                  // η tilt: difference in polar angle (θ) between axis and radial (r-z plane).
                  TVector3 e_r = Cint.Unit();
                  TVector3 e_a = axis.Unit();

                  const double theta_r = e_r.Theta();  // angle from +z
                  const double theta_a = e_a.Theta();
                  const double dtheta  = theta_a - theta_r;

                  etaTilt_mrad = dtheta * 1.0e3;
                }

                // Keep rotZ from the geometry for debugging / information.
                rotZ = g5->get_rotz();
              }

              bemc->SetVertexZ(static_cast<float>(z_vtx));

              if (actuallyFireGeant)
              {
                const double phiCenter = tg->get_phi();
                gun->set_vtx(0.0, 0.0, z_vtx);
                gun->set_phi_range(phiCenter, phiCenter);
                // Fun4AllServer::instance()->run(1);  // optional if you really want to simulate
              }
              float cosphi = 0.f, coseta = 0.f, aphis = 0.f, aetas = 0.f;
              incidenceEverCalled = true;
              const bool ok = bemc->CalculateMechIncidence(
                  static_cast<float>(EGeV), x_eval, y_eval,
                  cosphi, coseta, aphis, aetas);


              (void)ok;
              std::cout << std::fixed
                        << z_vtx        << ", "
                        << ieta         << ", "
                        << eta_det      << ", "
                        << eta_SD       << ", "
                        << aphis        << ", "
                        << aetas        << ", "
                        << cosphi       << ", "
                        << coseta       << ", "
                        << phiTilt_mrad << ", "
                        << etaTilt_mrad << ", "
                        << rotZ * 1.0e3 << "  # rotZ_mrad\n";

          }
        }
      }
      else if (scanMode == "dissectPhiSectors")
        {
          // Purely geometrical sector-by-sector φ-tilt scan.
          // Uses RawTowerGeomv5 detailed geometry to compute, for each tower:
          //   r_T  = (C_x, C_y, 0)
          //   a_T  = (C_ext,x - C_int,x, C_ext,y - C_int,y, 0)
          //   phi_radial = atan2(r_Ty, r_Tx)
          //   phi_axis   = atan2(a_Ty, a_Tx)
          //   δφ_tower   = wrap(phi_axis - phi_radial)
          //
          // For each sector s (8 φ-towers × 48 η-rows) we also compute:
          //   r_T^(s) = Σ r_T,  a_T^(s) = Σ a_T
          //   φ_radial^(s), φ_axis^(s), Δφ_sec(s) = wrap(φ_axis^(s) - φ_radial^(s))
          //   RMS residual of (δφ_tower - Δφ_sec).
          //
          // Output:
          //   • Per-sector tower table with r_T, a_T, φ_radial, φ_axis, δφ_tower.
          //   • Per-sector compact summary.
          //   • Final sector-by-sector summary + neighbour Δφ_radial checks.

          std::cout << "[MODE] dissectPhiSectors - tower and sector φ-tilts from RawTowerGeomv5\n";

          const char* ANSI_RED_BOLD   = "\033[1;31m";
          const char* ANSI_GREEN      = "\033[32m";
          const char* ANSI_YELLOW     = "\033[33m";
          const char* ANSI_CYAN_BOLD  = "\033[1;36m";
          const char* ANSI_RESET      = "\033[0m";

          const double RAD2DEG = 180.0 / std::acos(-1.0);

          // Geometry-based mapping:
          //   • 96 towers in η, split at ieta = 48 (nEta/2) into North and South halves.
          //   • 256 towers in φ, 32 sectors, 8 towers per sector in φ.
          //   • sector index 0-31  : North (η ∈ [48,95])
          //   • sector index 32-63 : South (η ∈ [0,47])
          const int nSectorsTotal      = 64;                 // 2 (±η) × 32 (φ)
          const int nPhiSectors        = nSectorsTotal / 2;  // 32
          const int phiTowersPerSector = nPhi / nPhiSectors; // 8 for 256 φ towers
          const int etaSplitIndex      = nEta / 2;           // 48 for 96 η towers

          // Safety bound for towers in a sector: 48 η-rows × 8 φ-towers = 384.
          static const int kMaxTowersPerSector = 512;

          // Arrays for final sector-by-sector summary.
          double sec_phiRad[64]      = {0.0};  // φ_radial^(s)
          double sec_phiAxis[64]     = {0.0};  // φ_axis^(s)
          double sec_deltaSec[64]    = {0.0};  // Δφ_sec(s) [rad]
          double sec_rmsTower[64]    = {0.0};  // RMS(δφ_tower) [mrad]
          double sec_rmsResid[64]    = {0.0};  // RMS(δφ_tower - Δφ_sec) [mrad]
          int    sec_nTowers[64]     = {0};
          bool   sec_valid[64]       = {false};

          // Design φ-tilt and tolerances used only for coloring in final summary.
          const double designTilt_rad = 0.0876;   // 87.6 mrad from design slide
          const double tightTol_rad   = 0.005;    // ±5 mrad window → green
          const double looseTol_rad   = 0.010;    // ±10 mrad window → yellow

          // Per-sector loop.
          for (int sector = 0; sector < nSectorsTotal; ++sector)
          {
            const bool isNorth      = (sector < nPhiSectors);       // 0-31: η ≥ 48
            const int  phiSectorIdx = sector % nPhiSectors;         // 0-31
            const int  iphi_min     = phiSectorIdx * phiTowersPerSector;
            const int  iphi_max     = iphi_min + phiTowersPerSector - 1;
            const int  ieta_min     = isNorth ? etaSplitIndex : 0;
            const int  ieta_max     = isNorth ? (nEta - 1)    : (etaSplitIndex - 1);

            // Per-sector accumulators for rigid-body directions and tower δφ stats.
            TVector3 sum_rT(0., 0., 0.);
            TVector3 sum_aT(0., 0., 0.);
            double   sum_delta    = 0.0; // Σ δφ_tower (rad)
            double   sum_delta2   = 0.0; // Σ δφ_tower^2 (rad^2)
            double   sum_phiRad   = 0.0; // Σ φ_radial (rad)
            int      nTowers      = 0;

            // Per-tower storage (for printing and residual RMS).
            int    t_ieta[kMaxTowersPerSector];
            int    t_iphi[kMaxTowersPerSector];
            int    t_phiLocal[kMaxTowersPerSector];
            double t_rTx[kMaxTowersPerSector];
            double t_rTy[kMaxTowersPerSector];
            double t_aTx[kMaxTowersPerSector];
            double t_aTy[kMaxTowersPerSector];
            double t_phiRad[kMaxTowersPerSector];
            double t_phiAxis[kMaxTowersPerSector];
            double t_delta[kMaxTowersPerSector];  // δφ_tower (rad)

              // Fill tower-level information for this sector.
              for (int ieta = ieta_min; ieta <= ieta_max; ++ieta)
              {
                for (int iphi = iphi_min; iphi <= iphi_max; ++iphi)
                {
                  if (nTowers >= kMaxTowersPerSector) continue; // safety

                  const RawTowerDefs::keytype key =
                    RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi);
                  RawTowerGeom* tg = geo->get_tower_geometry(key);
                  if (!tg) continue;

                  auto* g5 = dynamic_cast<RawTowerGeomv5*>(tg);
                  if (!g5) continue;

  #ifndef NDEBUG
                  // Cross-check that the helper mapping agrees with the current sector loop.
                  TowerHardwareIndex idx = mapTowerToHardwareIndex(ieta, iphi);
                  assert(idx.sector == sector);
  #endif

                  // Front / back face centers and bulk center in global frame.
                  TVector3 Cint(g5->get_center_int_x(),
                                g5->get_center_int_y(),
                                g5->get_center_int_z());
                  TVector3 Cext(g5->get_center_ext_x(),
                                g5->get_center_ext_y(),
                                g5->get_center_ext_z());
                  TVector3 axis = Cext - Cint;   // tower axis (3D)

                  if (axis.Mag2() == 0.) continue;

                  TVector3 C(g5->get_center_x(),
                             g5->get_center_y(),
                             g5->get_center_z());

                  // Radial direction in transverse plane and axis projection.
                  TVector3 rT(C.X(), C.Y(), 0.0);
                  TVector3 aT(axis.X(), axis.Y(), 0.0);

                  if (rT.Mag2() == 0. || aT.Mag2() == 0.) continue;

                  TVector3 e_rT = rT.Unit();
                  TVector3 e_aT = aT.Unit();

                  const double phi_radial = std::atan2(e_rT.Y(), e_rT.X());
                  const double phi_axis   = std::atan2(e_aT.Y(), e_aT.X());
                  const double delta_phi  = TVector2::Phi_mpi_pi(phi_axis - phi_radial); // δφ_tower

                  // Record tower row.
                  const int phiLocal = iphi - iphi_min; // 0..7 within the sector

                  t_ieta[nTowers]     = ieta;
                  t_iphi[nTowers]     = iphi;
                  t_phiLocal[nTowers] = phiLocal;
                  t_rTx[nTowers]      = rT.X();
                  t_rTy[nTowers]      = rT.Y();
                  t_aTx[nTowers]      = aT.X();
                  t_aTy[nTowers]      = aT.Y();
                  t_phiRad[nTowers]   = phi_radial;
                  t_phiAxis[nTowers]  = phi_axis;
                  t_delta[nTowers]    = delta_phi;

                  // Accumulate sector sums.
                  sum_rT     += rT;
                  sum_aT     += aT;
                  sum_delta  += delta_phi;
                  sum_delta2 += delta_phi * delta_phi;
                  sum_phiRad += phi_radial;
                  ++nTowers;
                }
              }

            if (nTowers == 0)
            {
              std::cout << "\n" << ANSI_RED_BOLD
                        << "[dissectPhiSectors] Sector " << sector
                        << " (" << (isNorth ? "N" : "S")
                        << ", phiIdx=" << phiSectorIdx
                        << ") has no valid towers."
                        << ANSI_RESET << "\n";
              continue;
            }

            // Sector rigid-body radial / axis directions.
            TVector3 eR = sum_rT.Unit();
            TVector3 eA = sum_aT.Unit();

            const double phiRad_sector  = std::atan2(eR.Y(), eR.X());
            const double phiAxis_sector = std::atan2(eA.Y(), eA.X());
            const double deltaSec_rad   =
              TVector2::Phi_mpi_pi(phiAxis_sector - phiRad_sector); // Δφ_sec(s)

            // Tower-level statistics.
            const double meanDelta = sum_delta / static_cast<double>(nTowers);
            double varDelta        = sum_delta2 / static_cast<double>(nTowers) - meanDelta * meanDelta;
            if (varDelta < 0.) varDelta = 0.;
            const double rmsDelta_rad   = std::sqrt(varDelta);          // RMS(δφ_tower) in rad
            const double rmsDelta_mrad  = rmsDelta_rad * 1.0e3;         // in mrad

            // Residuals δφ_resid = δφ_tower - Δφ_sec(s), RMS in mrad.
            double sumResid2_mrad = 0.0;
            for (int it = 0; it < nTowers; ++it)
            {
              const double resid_rad  = t_delta[it] - deltaSec_rad;
              const double resid_mrad = resid_rad * 1.0e3;
              sumResid2_mrad += resid_mrad * resid_mrad;
            }
            const double rmsResid_mrad =
              std::sqrt(sumResid2_mrad / static_cast<double>(nTowers));

            // Store for global summary.
            sec_phiRad[sector]   = phiRad_sector;
            sec_phiAxis[sector]  = phiAxis_sector;
            sec_deltaSec[sector] = deltaSec_rad;
            sec_rmsTower[sector] = rmsDelta_mrad;
            sec_rmsResid[sector] = rmsResid_mrad;
            sec_nTowers[sector]  = nTowers;
            sec_valid[sector]    = true;

            // -------------------------------
            // Per-sector tower-level table
            // -------------------------------
            std::cout << "\n" << ANSI_CYAN_BOLD
                      << "[dissectPhiSectors] Sector " << sector
                      << " (" << (isNorth ? "N" : "S")
                      << ", phiIdx=" << phiSectorIdx
                      << ")  ieta=[" << ieta_min << "," << ieta_max
                      << "], iphi=[" << iphi_min << "," << iphi_max << "]  Ntow=" << nTowers
                      << ANSI_RESET << "\n";

            std::cout << "  sec hemi phiIdx  ieta iphi  φLoc"
                      << "        rTx        rTy        aTx        aTy"
                      << "    φ_radial      φ_axis   δφ_tower[mrad]\n";
            std::cout << "  ----------------------------------------------------------------------------------------\n";

            for (int it = 0; it < nTowers; ++it)
            {
              const double delta_mrad = t_delta[it] * 1.0e3;

              std::cout << "  "
                        << sector
                        << "   " << (isNorth ? "N" : "S")
                        << "      " << phiSectorIdx
                        << "   " << t_ieta[it]
                        << "   " << t_iphi[it]
                        << "    " << t_phiLocal[it]
                        << "   " << t_rTx[it]
                        << "   " << t_rTy[it]
                        << "   " << t_aTx[it]
                        << "   " << t_aTy[it]
                        << "   " << t_phiRad[it]
                        << "   " << t_phiAxis[it]
                        << "   " << delta_mrad
                        << "\n";
            }

            // -------------------------------
            // Per-sector compact summary
            // -------------------------------
            const double deltaSec_mrad = deltaSec_rad * 1.0e3;
            const double diffToDesign  = deltaSec_rad - designTilt_rad;
            const double absDiff       = std::fabs(diffToDesign);

            const char* colorTilt =
              (absDiff < tightTol_rad) ? ANSI_GREEN :
              (absDiff < looseTol_rad) ? ANSI_YELLOW : ANSI_RED_BOLD;

            std::cout << "  ----------------------------------------------------------------------------------------\n";
            std::cout << "  Sector summary:\n";
            std::cout << "    φ_radial^(s) = " << phiRad_sector
                      << " rad  (" << phiRad_sector * RAD2DEG << " deg)\n";
            std::cout << "    φ_axis^(s)   = " << phiAxis_sector
                      << " rad  (" << phiAxis_sector * RAD2DEG << " deg)\n";
            std::cout << "    Δφ_sec(s)    = " << colorTilt << deltaSec_mrad << ANSI_RESET
                      << " mrad   [design ≈ 87.6 mrad]\n";
            std::cout << "    RMS(δφ_tower)        = " << rmsDelta_mrad  << " mrad\n";
            std::cout << "    RMS(δφ_tower-Δφ_sec) = " << rmsResid_mrad << " mrad\n";
          } // end sector loop

            // ------------------------------------------------------------------
            // Global sector-by-sector summary (North and South)
            // Pretty, column-aligned ANSI table + global stats.
            // ------------------------------------------------------------------
            std::cout << "\n[dissectPhiSectors] Sector-by-sector summary (rigid-body yaw & residuals)\n";

            // Column widths for a clean, linearized layout.
            const int wSec       = 4;   // sector index
            const int wHemi      = 4;   // N / S
            const int wPhiIdx    = 6;   // 0..31
            const int wNtow      = 6;   // number of towers
            const int wPhiRadDeg = 16;  // φ_radial[deg]
            const int wPhiAxisDeg= 16;  // φ_axis[deg]
            const int wDphiSec   = 14;  // Δφ_sec[mrad]
            const int wRmsTower  = 22;  // RMS(δφ_tower)[mrad]
            const int wRmsResid  = 22;  // RMS(δφ_resid)[mrad]

            const char* ANSI_HDR      = "\033[1;36m"; // bold cyan header
            const char* ANSI_HDR_LINE = "\033[36m";   // cyan underline
            const char* ANSI_SUM_HDR  = "\033[1;35m"; // bold magenta for global summary

            // Header row
            std::cout << "  "
                      << ANSI_HDR
                      << std::left  << std::setw(wSec)       << "sec"   << " "
                      << std::left  << std::setw(wHemi)      << "hemi"  << " "
                      << std::left  << std::setw(wPhiIdx)    << "phiIdx"<< " "
                      << std::right << std::setw(wNtow)      << "Ntow"  << " "
                      << std::right << std::setw(wPhiRadDeg) << "φ_radial[deg]"        << " "
                      << std::right << std::setw(wPhiAxisDeg)<< "φ_axis[deg]"          << " "
                      << std::right << std::setw(wDphiSec)   << "Δφ_sec[mrad]"         << " "
                      << std::right << std::setw(wRmsTower)  << "RMS(δφ_tower)[mrad]"  << " "
                      << std::right << std::setw(wRmsResid)  << "RMS(δφ_resid)[mrad]"
                      << ANSI_RESET << "\n";

            const int totalWidth =
                2 + wSec + 1 + wHemi + 1 + wPhiIdx + 1 + wNtow + 1 +
                wPhiRadDeg + 1 + wPhiAxisDeg + 1 + wDphiSec + 1 +
                wRmsTower + 1 + wRmsResid;

            std::cout << "  "
                      << ANSI_HDR_LINE << std::string(totalWidth, '-') << ANSI_RESET
                      << "\n";

            std::cout << std::fixed;

            // Global accumulators for sector-level stats
            double globalSumTilt      = 0.0;   // Σ Δφ_sec (rad)
            double globalSumTilt2     = 0.0;   // Σ Δφ_sec^2 (rad^2)
            double globalMinTilt_rad  =  1.0e9;
            double globalMaxTilt_rad  = -1.0e9;
            double globalSumRmsTower  = 0.0;   // Σ RMS(δφ_tower) (mrad)
            double globalSumRmsResid  = 0.0;   // Σ RMS(δφ_tower-Δφ_sec) (mrad)
            int    globalTiltN        = 0;     // number of valid sectors

            for (int sector = 0; sector < nSectorsTotal; ++sector)
            {
              if (!sec_valid[sector]) continue;

              const bool  isNorth      = (sector < nPhiSectors);
              const int   phiSectorIdx = sector % nPhiSectors;
              const char* hemi         = isNorth ? "N" : "S";

              const double phiRad_deg    = sec_phiRad[sector]  * RAD2DEG;
              const double phiAxis_deg   = sec_phiAxis[sector] * RAD2DEG;
              const double deltaSec_rad  = sec_deltaSec[sector];
              const double deltaSec_mrad = deltaSec_rad * 1.0e3;

              const double diffToDesign = deltaSec_rad - designTilt_rad;
              const double absDiff      = std::fabs(diffToDesign);
              const char*  colorTilt    =
                (absDiff < tightTol_rad) ? ANSI_GREEN :
                (absDiff < looseTol_rad) ? ANSI_YELLOW : ANSI_RED_BOLD;

              // Accumulate global stats
              globalSumTilt   += deltaSec_rad;
              globalSumTilt2  += deltaSec_rad * deltaSec_rad;
              if (deltaSec_rad < globalMinTilt_rad) globalMinTilt_rad = deltaSec_rad;
              if (deltaSec_rad > globalMaxTilt_rad) globalMaxTilt_rad = deltaSec_rad;
              globalSumRmsTower += sec_rmsTower[sector];
              globalSumRmsResid += sec_rmsResid[sector];
              ++globalTiltN;

              std::cout << "  "
                        << std::right << std::setw(wSec)    << sector << " ";

              // hemi column (optionally tinted: green for N, yellow for S)
              const char* hemiColor = isNorth ? ANSI_GREEN : ANSI_YELLOW;
              std::cout << hemiColor
                        << std::left  << std::setw(wHemi)   << hemi
                        << ANSI_RESET << " ";

              std::cout << std::right << std::setw(wPhiIdx) << phiSectorIdx << " "
                        << std::right << std::setw(wNtow)   << sec_nTowers[sector] << " "
                        << std::right << std::setw(wPhiRadDeg)
                        << std::setprecision(4) << phiRad_deg  << " "
                        << std::right << std::setw(wPhiAxisDeg)
                        << std::setprecision(4) << phiAxis_deg << " ";

              // Δφ_sec column, colour-coded vs design
              std::cout << colorTilt
                        << std::right << std::setw(wDphiSec)
                        << std::setprecision(3) << deltaSec_mrad
                        << ANSI_RESET << " ";

              // RMS columns
              std::cout << std::right << std::setw(wRmsTower)
                        << std::setprecision(4) << sec_rmsTower[sector] << " "
                        << std::right << std::setw(wRmsResid)
                        << std::setprecision(4) << sec_rmsResid[sector]
                        << "\n";
            }

            if (globalTiltN > 0)
            {
              const double meanTilt_rad   = globalSumTilt / static_cast<double>(globalTiltN);
              const double meanTilt_mrad  = meanTilt_rad * 1.0e3;

              double varTilt_rad = globalSumTilt2 / static_cast<double>(globalTiltN)
                                 - meanTilt_rad * meanTilt_rad;
              if (varTilt_rad < 0.0) varTilt_rad = 0.0;
              const double rmsSpread_mrad = std::sqrt(varTilt_rad) * 1.0e3;

              const double meanRmsTower_mrad = globalSumRmsTower / static_cast<double>(globalTiltN);
              const double meanRmsResid_mrad = globalSumRmsResid / static_cast<double>(globalTiltN);
              const double minTilt_mrad      = globalMinTilt_rad * 1.0e3;
              const double maxTilt_mrad      = globalMaxTilt_rad * 1.0e3;

              std::cout << "\n" << ANSI_SUM_HDR
                        << "[dissectPhiSectors] Global sector tilt summary ("
                        << globalTiltN << " sectors)\n"
                        << ANSI_RESET;

              std::cout << "  <Δφ_sec(axis-radial)>        = "
                        << std::fixed << std::setprecision(3) << meanTilt_mrad
                        << " mrad  (" << std::setprecision(6) << meanTilt_rad << " rad)\n";

              std::cout << "  RMS[Δφ_sec - global mean]    = "
                        << std::setprecision(3) << rmsSpread_mrad << " mrad\n";

              std::cout << "  min/max Δφ_sec               = "
                        << std::setprecision(3) << minTilt_mrad
                        << " / " << maxTilt_mrad << " mrad\n";

              std::cout << "  <RMS(δφ_tower)>              = "
                        << std::setprecision(4) << meanRmsTower_mrad << " mrad\n";

              std::cout << "  <RMS(δφ_tower-Δφ_sec)>       = "
                        << std::setprecision(4) << meanRmsResid_mrad << " mrad\n";

                // ----------------------------------------------------------------
                // Build ROOT graphs *without* error bars:
                //   1) |Δφ_sec|[mrad] vs sector index (0..63)
                //   2) RMS(δφ_resid)[mrad] vs sector index (0..63)
                //
                // All ey-errors are explicitly set to zero so no error bars are
                // drawn; we just use black circle markers. We also add a bit of
                // headroom on the y-axis so the labels do not clash with tick
                // labels and the cloud of points has visible margins.
                // ----------------------------------------------------------------
                const int maxPoints = nSectorsTotal;

                double xSec[64];
                double yAbsDelta_mrad[64];
                double exSec[64];
                double eyAbsDelta_mrad[64];

                double yRmsResid_mrad[64];
                double exSecRms[64];
                double eyRmsResid_mrad[64];

                // Track min/max for y-axis padding
                double minAbsDelta =  1.0e9;
                double maxAbsDelta = -1.0e9;
                double minRmsResid =  1.0e9;
                double maxRmsResid = -1.0e9;

                int nPoints = 0;
                for (int sector = 0; sector < nSectorsTotal; ++sector)
                {
                  if (!sec_valid[sector]) continue;

                  const int ntow = sec_nTowers[sector];
                  if (ntow <= 0) continue;

                  const double deltaSec_mrad    = sec_deltaSec[sector] * 1.0e3;
                  const double absDeltaSec_mrad = std::fabs(deltaSec_mrad);
                  const double rmsResid         = sec_rmsResid[sector];

                  xSec[nPoints]            = static_cast<double>(sector);

                  // y-values for |Δφ_sec|, no error bars
                  yAbsDelta_mrad[nPoints]  = absDeltaSec_mrad;
                  exSec[nPoints]           = 0.0;
                  eyAbsDelta_mrad[nPoints] = 0.0;

                  // y-values for RMS(δφ_resid), no error bars
                  yRmsResid_mrad[nPoints]  = rmsResid;
                  exSecRms[nPoints]        = 0.0;
                  eyRmsResid_mrad[nPoints] = 0.0;

                  // Track min/max for padding
                  if (absDeltaSec_mrad < minAbsDelta) minAbsDelta = absDeltaSec_mrad;
                  if (absDeltaSec_mrad > maxAbsDelta) maxAbsDelta = absDeltaSec_mrad;
                  if (rmsResid         < minRmsResid) minRmsResid = rmsResid;
                  if (rmsResid         > maxRmsResid) maxRmsResid = rmsResid;

                  ++nPoints;
                }

                if (nPoints > 0)
                {
                  // Disable stats box on plots
                  gStyle->SetOptStat(0);

                  // ----------------------------
                  // Graph 1: |Δφ_sec| vs sector
                  // ----------------------------
                  TGraphErrors* gDeltaSec = new TGraphErrors(
                      nPoints, xSec, yAbsDelta_mrad, exSec, eyAbsDelta_mrad);
                  gDeltaSec->SetTitle("|#Delta#phi_{sec}| vs sector;sector index;|#Delta#phi_{sec}| [mrad]");
                  gDeltaSec->SetMarkerStyle(20);   // filled circle
                  gDeltaSec->SetMarkerSize(1.1);
                  gDeltaSec->SetMarkerColor(kBlack);
                  gDeltaSec->SetLineColor(kBlack);
                  gDeltaSec->SetLineWidth(1);

                  TCanvas* cDelta = new TCanvas(
                      "cDeltaSec", "|Delta phi_sec| vs sector", 1200, 600);
                  cDelta->SetGridx();
                  cDelta->SetGridy();
                  cDelta->SetMargin(0.14, 0.04, 0.14, 0.06); // extra left/bottom margin

                  gDeltaSec->Draw("AP");
                  gDeltaSec->GetXaxis()->SetLimits(-0.5, nSectorsTotal - 0.5);
                  gDeltaSec->GetXaxis()->SetNdivisions(512);
                  gDeltaSec->GetXaxis()->SetTitleOffset(1.1);
                  gDeltaSec->GetYaxis()->SetTitleOffset(1.6); // push label away from ticks

                  // Add some vertical padding so points are not crammed
                  double spanAbs = maxAbsDelta - minAbsDelta;
                  if (spanAbs <= 0.0) spanAbs = std::max(1.0e-3, 0.001 * maxAbsDelta);
                  double padAbs  = 0.25 * spanAbs; // 25% padding top/bottom
                  gDeltaSec->GetYaxis()->SetRangeUser(minAbsDelta - padAbs,
                                                      maxAbsDelta + padAbs);

                  cDelta->SaveAs("dissectPhiSectors_DeltaPhiSec_vs_sector.png");
                  std::cout << "[dissectPhiSectors] Saved PNG: dissectPhiSectors_DeltaPhiSec_vs_sector.png\n";

                  // ----------------------------------------
                  // Graph 2: RMS(δφ_resid) vs sector
                  // ----------------------------------------
                  TGraphErrors* gRmsResid = new TGraphErrors(
                      nPoints, xSec, yRmsResid_mrad, exSecRms, eyRmsResid_mrad);
                  gRmsResid->SetTitle("RMS(#delta#phi_{resid}) vs sector;sector index;RMS(#delta#phi_{resid}) [mrad]");
                  gRmsResid->SetMarkerStyle(20);   // filled circle
                  gRmsResid->SetMarkerSize(1.1);
                  gRmsResid->SetMarkerColor(kBlack);
                  gRmsResid->SetLineColor(kBlack);
                  gRmsResid->SetLineWidth(1);

                  TCanvas* cRms = new TCanvas(
                      "cRmsResid", "RMS delta phi_resid vs sector", 1200, 600);
                  cRms->SetGridx();
                  cRms->SetGridy();
                  cRms->SetMargin(0.14, 0.04, 0.14, 0.06);

                  gRmsResid->Draw("AP");
                  gRmsResid->GetXaxis()->SetLimits(-0.5, nSectorsTotal - 0.5);
                  gRmsResid->GetXaxis()->SetNdivisions(512);
                  gRmsResid->GetXaxis()->SetTitleOffset(1.1);
                  gRmsResid->GetYaxis()->SetTitleOffset(1.6);

                  double spanRms = maxRmsResid - minRmsResid;
                  if (spanRms <= 0.0) spanRms = std::max(1.0e-3, 0.001 * maxRmsResid);
                  double padRms  = 0.25 * spanRms;
                  gRmsResid->GetYaxis()->SetRangeUser(minRmsResid - padRms,
                                                      maxRmsResid + padRms);

                  cRms->SaveAs("dissectPhiSectors_RMSResid_vs_sector.png");
                  std::cout << "[dissectPhiSectors] Saved PNG: dissectPhiSectors_RMSResid_vs_sector.png\n";
                }

            }


            // ------------------------------------------------------------------
            // Neighbour-sector Δφ_radial checks (~11.25° step between φ wedges)
            // ------------------------------------------------------------------
            const double stepExpected =
              2.0 * std::acos(-1.0) / static_cast<double>(nPhiSectors); // 2π / 32

            // We will build two arrays:
            //   • dphiNorth[i] = Δφ_radial between sector i and i+1  (i = 0…30)
            //   • dphiSouth[i] = Δφ_radial between sector 32+i and 33+i
            const int nSteps = nPhiSectors - 1; // 31 neighbour gaps per hemisphere

            double dphiNorth[32] = {0.0};
            double dphiSouth[32] = {0.0};
            bool   haveNorth[32] = {false};
            bool   haveSouth[32] = {false};

            for (int i = 0; i < nSteps; ++i)
            {
              // North: 0→1, 1→2, …, 30→31
              if (sec_valid[i] && sec_valid[i + 1])
              {
                dphiNorth[i] = TVector2::Phi_mpi_pi(sec_phiRad[i + 1] - sec_phiRad[i]);
                haveNorth[i] = true;
              }

              // South: 32→33, 33→34, …, 62→63
              const int sA = nPhiSectors + i;   // 32+i
              const int sB = sA + 1;           // 33+i
              if (sec_valid[sA] && sec_valid[sB])
              {
                dphiSouth[i] = TVector2::Phi_mpi_pi(sec_phiRad[sB] - sec_phiRad[sA]);
                haveSouth[i] = true;
              }
            }

            // Column widths: include extra columns showing the two φ values
            // that go into each neighbour Δφ_radial.
            // Make the "s→s+1" and φ columns wider so the numeric columns
            // are visually separated and not bunched up.
            const int wPair   = 16;  // "North s→s+1" / "South s→s+1"
            const int wPhiVal = 16;  // each φ_s[rad] / φ_{s+1}[rad] column
            const int wRad    = 18;  // Δφ_radial[rad]
            const int wDeg    = 12;  // Δφ[deg]

            std::cout << "\n[dissectPhiSectors] Neighbour-sector Δφ_radial checks (North vs South)\n";

            // Header row - one block for North, one for South.
            std::cout << "  "
                      << std::left << std::setw(wPair)   << "North s→s+1"
                      << std::left << std::setw(wPhiVal) << "φ_s[rad]"
                      << std::left << std::setw(wPhiVal) << "φ_{s+1}[rad]"
                      << std::left << std::setw(wRad)    << "Δφ_radial[rad]"
                      << std::left << std::setw(wDeg)    << "Δφ[deg]"
                      << "   "
                      << std::left << std::setw(wPair)   << "South s→s+1"
                      << std::left << std::setw(wPhiVal) << "φ_s[rad]"
                      << std::left << std::setw(wPhiVal) << "φ_{s+1}[rad]"
                      << std::left << std::setw(wRad)    << "Δφ_radial[rad]"
                      << std::left << std::setw(wDeg)    << "Δφ[deg]"
                      << "\n";

            const int halfWidth =
              wPair + wPhiVal + wPhiVal + wRad + wDeg + 3;  // +3 for the "   "
            const int totalWidthPairs = 2 * halfWidth;

            std::cout << "  " << std::string(totalWidthPairs, '-') << "\n";

            std::cout << std::fixed;

            // Data rows
            for (int i = 0; i < nSteps; ++i)
            {
              std::cout << "  ";

              // ---------- Left half: North pair i→i+1 ----------
              std::string northPair = haveNorth[i]
                ? (std::to_string(i) + "→" + std::to_string(i + 1))
                : "";

              if (!northPair.empty())
              {
                std::cout << ANSI_RED_BOLD
                          << std::left << std::setw(wPair) << northPair
                          << ANSI_RESET;
              }
              else
              {
                std::cout << std::left << std::setw(wPair) << " ";
              }

              if (haveNorth[i])
              {
                const double phiA = sec_phiRad[i];
                const double phiB = sec_phiRad[i + 1];

                std::cout << std::left << std::setw(wPhiVal)
                          << std::setprecision(6) << phiA
                          << std::left << std::setw(wPhiVal)
                          << std::setprecision(6) << phiB
                          << std::left << std::setw(wRad)
                          << std::setprecision(6) << dphiNorth[i]
                          << std::left << std::setw(wDeg)
                          << std::setprecision(2) << dphiNorth[i] * RAD2DEG;
              }
              else
              {
                std::cout << std::left << std::setw(wPhiVal) << " "
                          << std::left << std::setw(wPhiVal) << " "
                          << std::left << std::setw(wRad)    << " "
                          << std::left << std::setw(wDeg)    << " ";
              }

              // Spacer between North and South halves
              std::cout << "   ";

              // ---------- Right half: South pair (32+i)→(33+i) ----------
              const int sA = nPhiSectors + i;
              const int sB = sA + 1;

              std::string southPair = haveSouth[i]
                ? (std::to_string(sA) + "→" + std::to_string(sB))
                : "";

              if (!southPair.empty())
              {
                std::cout << ANSI_RED_BOLD
                          << std::left << std::setw(wPair) << southPair
                          << ANSI_RESET;
              }
              else
              {
                std::cout << std::left << std::setw(wPair) << " ";
              }

              if (haveSouth[i])
              {
                const double phiA = sec_phiRad[sA];
                const double phiB = sec_phiRad[sB];

                std::cout << std::left << std::setw(wPhiVal)
                          << std::setprecision(6) << phiA
                          << std::left << std::setw(wPhiVal)
                          << std::setprecision(6) << phiB
                          << std::left << std::setw(wRad)
                          << std::setprecision(6) << dphiSouth[i]
                          << std::left << std::setw(wDeg)
                          << std::setprecision(2) << dphiSouth[i] * RAD2DEG;
              }
              else
              {
                std::cout << std::left << std::setw(wPhiVal) << " "
                          << std::left << std::setw(wPhiVal) << " "
                          << std::left << std::setw(wRad)    << " "
                          << std::left << std::setw(wDeg)    << " ";
              }

              std::cout << "\n";
            }

            std::cout << "\n  Expected step ≈ " << stepExpected
                      << " rad ≈ " << stepExpected * RAD2DEG << " deg\n";

            std::cout << "\n[dissectPhiSectors] Expectation checks:\n"
                      << "  • Within each sector, Δφ_sec(s) ≈ +0.0876 rad (~87.6 mrad) if the encoded\n"
                      << "    geometry matches the design (rigid yaw of each 2×2 wedge).\n"
                      << "  • Between neighbouring φ sectors, Δφ_radial ≈ 2π/32 ≈ "
                      << stepExpected << " rad (~" << stepExpected * RAD2DEG << "°)\n"
                      << "    for both North (0-31) and South (32-63) sector centers.\n";
          }

      else if (scanMode == "dissectEtaBlocks")
      {
        std::cout << "[MODE] dissectEtaBlocks - η-tilts of towers, 2×2 blocks, modules and sectors\n";

        const char* ANSI_RED_BOLD      = "\033[1;31m";
        const char* ANSI_GREEN         = "\033[32m";
        const char* ANSI_YELLOW        = "\033[33m";
        const char* ANSI_CYAN_BOLD     = "\033[1;36m";
        const char* ANSI_MAGENTA_BOLD  = "\033[1;35m";
        const char* ANSI_RESET         = "\033[0m";

        enum
        {
          kSectorsTotal     = 64,                  // 2 (±η) × 32 (φ)
          kPhiSectors       = kSectorsTotal / 2,   // 32 φ-sectors per endcap
          kModulesPerSector = 24,                  // 24 η-modules per sector
          kBlocksPerModule  = 4,                   // 4 φ-blocks per module
          kTowersPerBlock   = 4                    // 2×2 towers per block
        };

        const int phiTowersPerSector = nPhi / kPhiSectors; // expected 8
        const int etaSplitIndex      = nEta / 2;           // expected 48

        assert(phiTowersPerSector * kPhiSectors == nPhi);
        assert(2 * kModulesPerSector == etaSplitIndex);
        assert(2 * kBlocksPerModule  == phiTowersPerSector);
        assert(kTowersPerBlock == 4);

          // Per-sector summary
          double sectorTilt_rad      [kSectorsTotal] = {0.0};
          double sectorRmsBlocks_rad [kSectorsTotal] = {0.0};
          double sectorRmsModules_rad[kSectorsTotal] = {0.0};
          double sectorRmsTowers_rad [kSectorsTotal] = {0.0};
          int    sectorNTowers       [kSectorsTotal] = {0};
          bool   sectorValid         [kSectorsTotal] = {false};
          bool   sectorIsNorth       [kSectorsTotal] = {false};
          int    sectorPhiIdx        [kSectorsTotal] = {0};

          // Global sector-level accumulators
          double globalSumTilt_rad        = 0.0;
          double globalSumTilt2_rad       = 0.0;
          double globalSumRmsBlocks_mrad  = 0.0;
          double globalSumRmsModules_mrad = 0.0;
          double globalSumRmsTowers_mrad  = 0.0;
          double globalSumAbsTilt_mrad    = 0.0;
          double northSumTilt_rad         = 0.0;
          double southSumTilt_rad         = 0.0;
          int    globalNSectors           = 0;
          int    northNSectors            = 0;
          int    southNSectors            = 0;

          // Global per-module and per-block accumulators (across sectors), split by hemisphere.
          // These are used to build hemisphere-averaged module and block summaries.
          //
          // North hemisphere (sectors 0..31)
          double moduleTiltSumNorth_rad        [kModulesPerSector]            = {0.0};
          double moduleTiltSumNorth2_rad       [kModulesPerSector]            = {0.0};
          int    moduleTiltNorthNSectors       [kModulesPerSector]            = {0};

          double moduleRmsBlocksSumNorth_mrad  [kModulesPerSector]            = {0.0};
          double moduleRmsBlocksSumNorth2_mrad [kModulesPerSector]            = {0.0};
          int    moduleRmsBlocksNorthNSectors  [kModulesPerSector]            = {0};

          // Per-block (within module) η-tilt accumulators for North hemisphere.
          double blockTiltSumNorth_rad         [kModulesPerSector][kBlocksPerModule] = {{0.0}};
          double blockTiltSumNorth2_rad        [kModulesPerSector][kBlocksPerModule] = {{0.0}};
          int    blockTiltNorthNSectors        [kModulesPerSector][kBlocksPerModule] = {{0}};

          // South hemisphere (sectors 32..63)
          double moduleTiltSumSouth_rad        [kModulesPerSector]            = {0.0};
          double moduleTiltSumSouth2_rad       [kModulesPerSector]            = {0.0};
          int    moduleTiltSouthNSectors       [kModulesPerSector]            = {0};

          double moduleRmsBlocksSumSouth_mrad  [kModulesPerSector]            = {0.0};
          double moduleRmsBlocksSumSouth2_mrad [kModulesPerSector]            = {0.0};
          int    moduleRmsBlocksSouthNSectors  [kModulesPerSector]            = {0};

          // Per-block (within module) η-tilt accumulators for South hemisphere.
          double blockTiltSumSouth_rad         [kModulesPerSector][kBlocksPerModule] = {{0.0}};
          double blockTiltSumSouth2_rad        [kModulesPerSector][kBlocksPerModule] = {{0.0}};
          int    blockTiltSouthNSectors        [kModulesPerSector][kBlocksPerModule] = {{0}};

        const int wSec   = 4;
        const int wMod   = 4;
        const int wBlk   = 4;
        const int wTinB  = 5;
        const int wIeta  = 5;
        const int wIphi  = 5;
        const int wTh    = 16;
        const int wDelta = 18;

        for (int sector = 0; sector < kSectorsTotal; ++sector)
          {
            const bool isNorth       = (sector < kPhiSectors);
            const int  phiSectorIdx  = sector % kPhiSectors;
            const int  iphi_min      = phiSectorIdx * phiTowersPerSector;
            const int  iphi_max      = iphi_min + phiTowersPerSector - 1;
            const int  ieta_min      = isNorth ? etaSplitIndex : 0;
            const int  ieta_max      = isNorth ? (nEta - 1)    : (etaSplitIndex - 1);
            const int  sectorDisplay = sector + 1;

            sectorIsNorth[sector] = isNorth;
            sectorPhiIdx[sector]  = phiSectorIdx;

            std::cout << "\n"
                      << ANSI_CYAN_BOLD
                      << "======================================================================\n"
                      << "[dissectEtaBlocks] Sector " << sectorDisplay << " / " << kSectorsTotal
                      << "  (ID=" << sector << ", hemi=" << (isNorth ? "N" : "S")
                      << ", phiIdx=" << phiSectorIdx << ")\n"
                      << "  η-range: [" << ieta_min << "," << ieta_max << "]"
                      << "  ,  φ-range: [" << iphi_min << "," << iphi_max << "]\n"
                      << "======================================================================"
                      << ANSI_RESET << "\n\n";

          // Per-sector storage
          double towerTilt_rad      [kModulesPerSector][kBlocksPerModule][kTowersPerBlock];
          bool   towerValid         [kModulesPerSector][kBlocksPerModule][kTowersPerBlock];
          double blockTilt_rad      [kModulesPerSector][kBlocksPerModule];
          double blockRms_rad       [kModulesPerSector][kBlocksPerModule];
          bool   blockValid         [kModulesPerSector][kBlocksPerModule];
          double moduleTilt_rad     [kModulesPerSector];
          double moduleRmsBlocks_rad[kModulesPerSector];
          double moduleRmsTowers_rad[kModulesPerSector];
          bool   moduleValid        [kModulesPerSector];

          for (int m = 0; m < kModulesPerSector; ++m)
          {
            moduleTilt_rad[m]        = 0.0;
            moduleRmsBlocks_rad[m]   = 0.0;
            moduleRmsTowers_rad[m]   = 0.0;
            moduleValid[m]           = false;
            for (int b = 0; b < kBlocksPerModule; ++b)
            {
              blockTilt_rad[m][b] = 0.0;
              blockRms_rad[m][b]  = 0.0;
              blockValid[m][b]    = false;
              for (int t = 0; t < kTowersPerBlock; ++t)
              {
                towerTilt_rad[m][b][t] = 0.0;
                towerValid[m][b][t]    = false;
              }
            }
          }

          double sumBlockTilt_sector = 0.0;
          int    nBlocks_sector      = 0;
          int    nTowers_sector      = 0;

          // ───────────────────────────────────────────────────────────────
          // First pass: for each module and block, compute δθ_tower,
          // Δθ_block, RMS_block and print 2×2 block-level tables.
          // δθ_tower = θ_axis - θ_radial in the r-z plane.
          // ───────────────────────────────────────────────────────────────
          for (int module = 0; module < kModulesPerSector; ++module)
              {
                const int ietaLocal0    = 2 * module;
                const int ietaLocal1    = ietaLocal0 + 1;
                const int ieta0         = isNorth ? (etaSplitIndex + ietaLocal0) : ietaLocal0;
                const int ieta1         = isNorth ? (etaSplitIndex + ietaLocal1) : ietaLocal1;
                const int moduleDisplay = module + 1;

                std::cout << "\n" << ANSI_YELLOW
                          << "  [Sector " << sector << " (" << (isNorth ? "N" : "S")
                          << ", phiIdx=" << phiSectorIdx << ")] "
                          << "Module " << moduleDisplay << " / " << kModulesPerSector
                          << "  (ietaLocal=" << ietaLocal0 << "," << ietaLocal1
                          << " → global ieta=" << ieta0 << "," << ieta1 << ")"
                          << ANSI_RESET << "\n";

            for (int blockPhi = 0; blockPhi < kBlocksPerModule; ++blockPhi)
              {
              const int blockInSector = module * kBlocksPerModule + blockPhi;
              const int blockDisplay  = blockPhi + 1;


              double thetaRadial[kTowersPerBlock] = {0.0, 0.0, 0.0, 0.0};
              double thetaAxis  [kTowersPerBlock] = {0.0, 0.0, 0.0, 0.0};
              double deltaTheta [kTowersPerBlock] = {0.0, 0.0, 0.0, 0.0};
              int    gIeta      [kTowersPerBlock] = {0, 0, 0, 0};
              int    gIphi      [kTowersPerBlock] = {0, 0, 0, 0};
              bool   hasTower   [kTowersPerBlock] = {false, false, false, false};

              // 4 towers in the 2×2 block
              for (int etaInModule = 0; etaInModule < 2; ++etaInModule)
              {
                for (int phiInBlock = 0; phiInBlock < 2; ++phiInBlock)
                {
                  const int towerInBlock = etaInModule * 2 + phiInBlock;

                  const int ietaLocal = 2 * module + etaInModule;
                  const int iphiLocal = 2 * blockPhi + phiInBlock;

                  const int ieta = isNorth ? (etaSplitIndex + ietaLocal) : ietaLocal;
                  const int iphi = iphi_min + iphiLocal;

                  if (ieta < ieta_min || ieta > ieta_max ||
                      iphi < iphi_min || iphi > iphi_max)
                  {
                    continue;
                  }

                  const RawTowerDefs::keytype key =
                    RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi);
                  RawTowerGeom* tg = geo->get_tower_geometry(key);
                  auto* g5        = dynamic_cast<RawTowerGeomv5*>(tg);
                  if (!g5) continue;

                  TVector3 Cint(g5->get_center_int_x(),
                                g5->get_center_int_y(),
                                g5->get_center_int_z());
                  TVector3 Cext(g5->get_center_ext_x(),
                                g5->get_center_ext_y(),
                                g5->get_center_ext_z());
                  TVector3 axis = Cext - Cint;

                  if (Cint.Mag2() == 0.0 || axis.Mag2() == 0.0) continue;

                  TVector3 e_r = Cint.Unit();  // radial line IP→tower center
                  TVector3 e_a = axis.Unit();  // tower axis direction

                  const double theta_radial = e_r.Theta();                // angle in r-z plane
                  const double theta_axis   = e_a.Theta();
                  const double dtheta       = theta_axis - theta_radial;  // δθ_tower [rad]

                    thetaRadial[towerInBlock] = theta_radial;
                    thetaAxis  [towerInBlock] = theta_axis;
                    deltaTheta [towerInBlock] = dtheta;
                    gIeta      [towerInBlock] = ieta;
                    gIphi      [towerInBlock] = iphi;
                    hasTower   [towerInBlock] = true;

  #ifndef NDEBUG
                    // Cross-check logical → geometrical mapping for this tower.
                    // This asserts that the (sector, module, blockPhi, towerInBlock)
                    // bookkeeping matches the generic mapTowerToHardwareIndex() helper,
                    // and that we never fill the same logical tower twice.
                    {
                      TowerHardwareIndex hw = mapTowerToHardwareIndex(ieta, iphi);
                      assert(hw.sector       == sector);
                      assert(hw.module       == module);
                      assert(hw.blockPhi     == blockPhi);
                      assert(hw.towerInBlock == towerInBlock);
                      assert(!towerValid[module][blockPhi][towerInBlock]);
                    }
  #endif

                    towerTilt_rad [module][blockPhi][towerInBlock] = dtheta;
                    towerValid    [module][blockPhi][towerInBlock] = true;
                }
              }

              int    nValidInBlock = 0;
              double sumDeltaBlock = 0.0;
              for (int t = 0; t < kTowersPerBlock; ++t)
              {
                if (!hasTower[t]) continue;
                ++nValidInBlock;
                sumDeltaBlock += deltaTheta[t];
              }

              if (nValidInBlock == 0) continue;

              const double deltaBlock = sumDeltaBlock / static_cast<double>(nValidInBlock);

              double sumResid2 = 0.0;
              for (int t = 0; t < kTowersPerBlock; ++t)
              {
                if (!hasTower[t]) continue;
                const double resid = deltaTheta[t] - deltaBlock;
                sumResid2 += resid * resid;
              }
              const double rmsBlock =
                std::sqrt(sumResid2 / static_cast<double>(nValidInBlock));

              blockTilt_rad[module][blockPhi] = deltaBlock;
              blockRms_rad [module][blockPhi] = rmsBlock;
              blockValid   [module][blockPhi] = true;

              sumBlockTilt_sector += deltaBlock;
              ++nBlocks_sector;
              nTowers_sector += nValidInBlock;

              const int totalWidth =
                    wSec + wMod + wBlk + wTinB + wIeta + wIphi +
                    2 * wTh + wDelta + 8;

              std::cout << "\n    Block " << blockDisplay << " / " << kBlocksPerModule
                            << "  (module=" << module
                            << ", blockPhi=" << blockPhi
                            << ", blockInSector=" << blockInSector << ")\n";

              std::cout << "      "
                        << std::left << std::setw(wSec)   << "sec"   << " "
                        << std::setw(wMod)   << "mod"   << " "
                        << std::setw(wBlk)   << "blk"   << " "
                        << std::setw(wTinB)  << "tInB"  << " "
                        << std::setw(wIeta)  << "ieta"  << " "
                        << std::setw(wIphi)  << "iphi"  << " "
                        << std::setw(wTh)    << "theta_radial[rad]" << " "
                        << std::setw(wTh)    << "theta_axis[rad]"   << " "
                        << std::setw(wDelta) << "delta_theta[mrad]"
                        << "\n";
              std::cout << "      "
                        << std::string(totalWidth, '-')
                        << "\n";

              for (int t = 0; t < kTowersPerBlock; ++t)
              {
                if (!hasTower[t]) continue;
                const double delta_mrad = deltaTheta[t] * 1.0e3;

                std::cout << "      "
                          << std::right << std::setw(wSec)  << sector << " "
                          << std::setw(wMod)  << module << " "
                          << std::setw(wBlk)  << blockPhi << " "
                          << std::setw(wTinB) << t << " "
                          << std::setw(wIeta) << gIeta[t] << " "
                          << std::setw(wIphi) << gIphi[t] << " "
                          << std::setw(wTh)   << std::setprecision(6) << thetaRadial[t] << " "
                          << std::setw(wTh)   << std::setprecision(6) << thetaAxis[t]   << " "
                          << std::setw(wDelta)<< std::setprecision(3) << delta_mrad
                          << "\n";
              }

              std::cout << "      --> Δθ_block = "
                        << std::setprecision(3) << deltaBlock * 1.0e3
                        << " mrad,  RMS_block = "
                        << std::setprecision(3) << rmsBlock * 1.0e3
                        << " mrad\n";
            } // end blocks in this module

            // Module-level rigid-body tilt and RMS’s
            double sumBlockTilt_module = 0.0;
            int    nBlocks_module      = 0;
            for (int b = 0; b < kBlocksPerModule; ++b)
            {
              if (!blockValid[module][b]) continue;
              sumBlockTilt_module += blockTilt_rad[module][b];
              ++nBlocks_module;
            }

            if (nBlocks_module == 0) continue;

            const double theta_module =
              sumBlockTilt_module / static_cast<double>(nBlocks_module);
            moduleTilt_rad[module] = theta_module;
            moduleValid   [module] = true;

            double sumBlockResid2_m = 0.0;
            for (int b = 0; b < kBlocksPerModule; ++b)
            {
              if (!blockValid[module][b]) continue;
              const double resid = blockTilt_rad[module][b] - theta_module;
              sumBlockResid2_m += resid * resid;
            }
            moduleRmsBlocks_rad[module] =
              std::sqrt(sumBlockResid2_m / static_cast<double>(nBlocks_module));

            double sumTowerResid2_m = 0.0;
            int    nTowers_module   = 0;
            for (int b = 0; b < kBlocksPerModule; ++b)
            {
              for (int t = 0; t < kTowersPerBlock; ++t)
              {
                if (!towerValid[module][b][t]) continue;
                const double resid = towerTilt_rad[module][b][t] - theta_module;
                sumTowerResid2_m += resid * resid;
                ++nTowers_module;
              }
            }
            if (nTowers_module > 0)
            {
              moduleRmsTowers_rad[module] =
                std::sqrt(sumTowerResid2_m / static_cast<double>(nTowers_module));
            }

                  std::cout << ANSI_MAGENTA_BOLD
                            << "  Module summary [sector=" << sector
                            << ", module " << (module + 1) << " / " << kModulesPerSector << "]: "
                            << "Θ_module = " << std::fixed << std::setprecision(3)
                            << theta_module * 1.0e3 << " mrad,  "
                            << "RMS_blocks|m = " << moduleRmsBlocks_rad[module] * 1.0e3
                            << " mrad,  "
                            << "RMS_towers|m = " << moduleRmsTowers_rad[module] * 1.0e3
                            << " mrad"
                            << ANSI_RESET << "\n";
          } // end loop over modules in this sector

          if (nBlocks_sector == 0)
          {
            std::cout << ANSI_RED_BOLD
                      << "  [dissectEtaBlocks] Sector " << sector
                      << " has no valid blocks, skipping."
                      << ANSI_RESET << "\n";
            continue;
          }

          // Sector-level tilt: average block tilt across all 96 blocks
          const double phi_sector =
            sumBlockTilt_sector / static_cast<double>(nBlocks_sector);
          sectorTilt_rad[sector] = phi_sector;
          sectorNTowers [sector] = nTowers_sector;

          // RMS of block tilts in the sector
          double sumBlockResid2_s = 0.0;
          for (int m = 0; m < kModulesPerSector; ++m)
          {
            for (int b = 0; b < kBlocksPerModule; ++b)
            {
              if (!blockValid[m][b]) continue;
              const double resid = blockTilt_rad[m][b] - phi_sector;
              sumBlockResid2_s += resid * resid;
            }
          }
          sectorRmsBlocks_rad[sector] =
            std::sqrt(sumBlockResid2_s / static_cast<double>(nBlocks_sector));

          // RMS of module tilts in the sector
          double sumModuleResid2_s = 0.0;
          int    nModules_sector   = 0;
          for (int m = 0; m < kModulesPerSector; ++m)
          {
            if (!moduleValid[m]) continue;
            const double resid = moduleTilt_rad[m] - phi_sector;
            sumModuleResid2_s += resid * resid;
            ++nModules_sector;
          }
          if (nModules_sector > 0)
          {
            sectorRmsModules_rad[sector] =
              std::sqrt(sumModuleResid2_s / static_cast<double>(nModules_sector));
          }

              // RMS of tower tilts in the sector
              double sumTowerResid2_s = 0.0;
              int    nTowersUsed_s    = 0;
              for (int m = 0; m < kModulesPerSector; ++m)
              {
                for (int b = 0; b < kBlocksPerModule; ++b)
                {
                  for (int t = 0; t < kTowersPerBlock; ++t)
                  {
                    if (!towerValid[m][b][t]) continue;
                    const double resid = towerTilt_rad[m][b][t] - phi_sector;
                    sumTowerResid2_s += resid * resid;
                    ++nTowersUsed_s;
                  }
                }
              }
              if (nTowersUsed_s > 0)
              {
                sectorRmsTowers_rad[sector] =
                  std::sqrt(sumTowerResid2_s / static_cast<double>(nTowersUsed_s));
              }

                  sectorValid[sector] = true;

                  std::cout << ANSI_GREEN
                            << "  Sector summary [sector " << (sector + 1) << " / " << kSectorsTotal
                            << " (" << (isNorth ? "N" : "S")
                            << ", phiIdx=" << phiSectorIdx << ")]: "
                            << "Φ_sector = " << std::fixed << std::setprecision(3)
                            << phi_sector * 1.0e3 << " mrad,  "
                            << "RMS_blocks|s = " << sectorRmsBlocks_rad[sector] * 1.0e3
                            << " mrad,  "
                            << "RMS_modules|s = " << sectorRmsModules_rad[sector] * 1.0e3
                            << " mrad,  "
                            << "RMS_towers|s = " << sectorRmsTowers_rad[sector] * 1.0e3
                            << " mrad"
                            << ANSI_RESET << "\n\n";

              // ------------------------------------------------------------------
              // Enhanced: per-sector plot of block tilts vs module index
              //
              // X-axis layout:
              //   * Module m (0..23) occupies the bin [m, m+1]
              //   * Bin center at x = m + 0.5
              //   * Four block markers for that module are drawn inside that bin
              //     with small offsets around the center.
              //
              // Visual behaviour:
              //   * Vertical dashed lines at every integer x = 0..24 (bin edges)
              //   * Tick marks have no numeric labels
              //   * Dummy module labels "1..24" drawn at x = m+0.5 between the lines
              //   * All markers for module m live strictly between its two bin edges
              // ------------------------------------------------------------------
              {
                // One array per block (0..3), up to one point per module
                double xBlock[kBlocksPerModule][kModulesPerSector];
                double yBlock[kBlocksPerModule][kModulesPerSector];
                int    nBlock[kBlocksPerModule] = {0, 0, 0, 0};

                double yMin = 0.0;
                double yMax = 0.0;
                bool   firstPoint = true;

                // Horizontal offset inside each module bin around the center x = m + 0.5
                const double dx = 0.20;

                // Fill per-block arrays and track global y-range
                for (int m = 0; m < kModulesPerSector; ++m)
                {
                  const double moduleCenterX = static_cast<double>(m) + 0.5; // bin center: [m,m+1]

                  for (int b = 0; b < kBlocksPerModule; ++b)
                  {
                    if (!blockValid[m][b]) continue;

                    // blockPhi = 0..3 → four offsets inside the bin [m,m+1]
                    const double x = moduleCenterX + (static_cast<double>(b) - 1.5) * dx;
                    const double tilt_mrad = blockTilt_rad[m][b] * 1.0e3;

                    if (nBlock[b] >= kModulesPerSector) continue; // safety guard

                    xBlock[b][nBlock[b]] = x;
                    yBlock[b][nBlock[b]] = tilt_mrad;
                    ++nBlock[b];

                    if (firstPoint)
                    {
                      yMin = yMax = tilt_mrad;
                      firstPoint = false;
                    }
                    else
                    {
                      if (tilt_mrad < yMin) yMin = tilt_mrad;
                      if (tilt_mrad > yMax) yMax = tilt_mrad;
                    }
                  }
                }

                int totalPoints = 0;
                for (int b = 0; b < kBlocksPerModule; ++b)
                {
                  totalPoints += nBlock[b];
                }

                if (totalPoints > 0)
                {
                  // Ensure output directory exists (idempotent)
                  if (gSystem)
                  {
                    gSystem->mkdir("sector_module_summariesPNGs", true);
                  }

                    std::string cname  = "cSectorModuleTilt_" + std::to_string(sector);
                    std::string ctitle = "Sector " + std::to_string(sector);

                      // Wider canvas with margins tuned so the title has room at the top
                      TCanvas* cSec = new TCanvas(cname.c_str(), ctitle.c_str(), 1600, 450);
                      // SetMargin(left, right, bottom, top) in NDC
                      cSec->SetMargin(0.06, 0.02, 0.14, 0.08);
                      cSec->SetGrid(0, 0); // no grid lines

                  // Choose first non-empty block as the frame graph (for axes)
                  int frameBlock = -1;
                  for (int b = 0; b < kBlocksPerModule; ++b)
                  {
                    if (nBlock[b] > 0)
                    {
                      frameBlock = b;
                      break;
                    }
                  }

                  if (frameBlock >= 0)
                  {
                    TGraph* gFrame = new TGraph(nBlock[frameBlock],
                                                xBlock[frameBlock],
                                                yBlock[frameBlock]);

                      std::string gtitle = "Sector " + std::to_string(sector)
                                           + " -- Per Block (within each module) #Delta #theta_{block} (m, b) [mrad]"
                                           + ";Module Index;#Delta #theta_{block} (m, b) [mrad]";
                      gFrame->SetTitle(gtitle.c_str());
                      gFrame->SetMarkerStyle(20);   // filled circle
                      
                    gFrame->SetMarkerSize(1.0);
                    gFrame->SetMarkerColor(kBlue); // temporary, will overdraw with colored graphs
                    gFrame->SetLineColor(0);       // no connecting line

                    gFrame->Draw("AP");

                      // X axis: modules 1..24 as bins; range [0,24], edges at integer x
                      const double xMin = 0.0;
                      const double xMax = static_cast<double>(kModulesPerSector);

                      gFrame->GetXaxis()->SetLimits(xMin, xMax);
                      gFrame->GetXaxis()->SetNdivisions(kModulesPerSector, 0, 0);
                      gFrame->GetXaxis()->SetTickLength(0.02);
                      gFrame->GetXaxis()->CenterTitle(false);

                      // X-axis title and styling
                      gFrame->GetXaxis()->SetTitle("Module Index");
                      gFrame->GetXaxis()->SetTitleSize(0.045);   // readable size
                      gFrame->GetXaxis()->SetTitleOffset(1.0);   // distance from axis
                      // Hide numeric tick labels; use dummy module labels instead
                      gFrame->GetXaxis()->SetLabelSize(0.0);

                      // Y axis: range with modest padding so data fills more of the frame
                      double span = yMax - yMin;
                      double pad  = (span > 0.0) ? 0.12 * span : 0.5;

                      const double yLow  = yMin - pad;
                      const double yHigh = yMax + pad;

                      // Use TLatex-style label for Δη block tilt
                      gFrame->GetYaxis()->SetTitle("#Delta #theta_{block} (m, b) [mrad]");
                      gFrame->GetYaxis()->CenterTitle(true);
                      gFrame->GetYaxis()->SetTitleSize(0.050);   // larger y-axis title font
                      gFrame->GetYaxis()->SetTitleOffset(0.95);  // slight spacing from axis
                      gFrame->GetYaxis()->SetLabelSize(0.035);   // slightly larger tick labels
                      gFrame->GetYaxis()->SetRangeUser(yLow, yHigh);

                      // Block color scheme:
                      //   block 0 (block 1): kBlue
                      //   block 1 (block 2): kRed
                      //   block 2 (block 3): kGreen+2
                      //   block 3 (block 4): kBlack
                      int blockColor[kBlocksPerModule] = { kBlue, kRed, kGreen + 2, kBlack };

                      // Vertical dashed lines at every integer x = 0..24 (bin edges)
                      for (int edge = 0; edge <= kModulesPerSector; ++edge)
                      {
                        const double xEdge = static_cast<double>(edge);
                        TLine* l = new TLine(xEdge, yLow, xEdge, yHigh);
                        l->SetLineStyle(2);        // dashed
                        l->SetLineWidth(1);
                        l->SetLineColor(kGray + 1);
                        l->Draw("SAME");
                      }

                      // Draw all blocks on top of the frame with their colors
                      // and keep pointers so we can build a legend.
                      TGraph* gBlock[kBlocksPerModule] = { nullptr, nullptr, nullptr, nullptr };

                      for (int b = 0; b < kBlocksPerModule; ++b)
                      {
                        if (nBlock[b] <= 0) continue;

                        TGraph* gB = new TGraph(nBlock[b], xBlock[b], yBlock[b]);
                        gB->SetMarkerStyle(20);
                        gB->SetMarkerSize(1.1);
                        gB->SetMarkerColor(blockColor[b]);
                        gB->SetLineColor(0); // markers only
                        gB->Draw("P SAME");

                        gBlock[b] = gB;
                      }

                      // Add legend mapping colors to blocks within each module
                      {
                        // Place legend depending on hemisphere:
                        //   • North sectors (0..31, isNorth=true)  → bottom-right
                        //   • South sectors (32..63, isNorth=false) → top-left
                        double x1, y1, x2, y2;
                        if (isNorth)
                        {
                          // Same as previous behaviour: lower-right corner
                          x1 = 0.78; y1 = 0.30; x2 = 0.97; y2 = 0.52;
                        }
                        else
                        {
                          // Move legend to upper-left corner for South sectors
                          x1 = 0.10; y1 = 0.66; x2 = 0.30; y2 = 0.88;
                        }

                        TLegend* leg = new TLegend(x1, y1, x2, y2); // NDC coordinates
                        leg->SetBorderSize(0);
                        leg->SetFillStyle(0);
                        leg->SetTextSize(0.04);
                        leg->SetTextAlign(12); // left-justified text inside legend box

                        // Left-aligned, bold header above the marker labels
                        leg->SetHeader("#bf{Blocks within Modules}", "L");

                        for (int b = 0; b < kBlocksPerModule; ++b)
                        {
                          if (!gBlock[b]) continue;
                          leg->AddEntry(gBlock[b], Form("Block %d", b + 1), "P");
                        }

                        leg->Draw();
                      }


                      // Dummy module labels centered in each bin between dashed lines
                      TLatex latex;
                      latex.SetTextAlign(22);  // center horizontally and vertically
                      latex.SetTextSize(0.030);

                      // Place labels a bit above the bottom of the plotting region
                      const double yLabel = yLow + 0.05 * (yHigh - yLow);

                      for (int m = 0; m < kModulesPerSector; ++m)
                      {
                        const double xCenter = static_cast<double>(m) + 0.5;
                        latex.DrawLatex(xCenter, yLabel, Form("%d", m + 1));
                      }

                      // Save PNG for this sector
                      std::string outName = "sector_module_summariesPNGs/sector_"
                                            + std::to_string(sector) + ".png";
                      cSec->SaveAs(outName.c_str());
                  }

                  delete cSec;
                }
              }

              // Accumulate sector-level statistics
              globalSumTilt_rad        += phi_sector;
              globalSumTilt2_rad       += phi_sector * phi_sector;
              globalSumRmsBlocks_mrad  += sectorRmsBlocks_rad [sector] * 1.0e3;
              globalSumRmsModules_mrad += sectorRmsModules_rad[sector] * 1.0e3;
              globalSumRmsTowers_mrad  += sectorRmsTowers_rad[sector] * 1.0e3;
              globalSumAbsTilt_mrad    += std::fabs(phi_sector * 1.0e3);

              if (isNorth)
              {
                northSumTilt_rad += phi_sector;
                ++northNSectors;
              }
              else
              {
                southSumTilt_rad += phi_sector;
                ++southNSectors;
              }

              ++globalNSectors;

              // Accumulate per-module and per-block statistics across sectors, split by hemisphere:
              //   - Φ_module (mean η-tilt per module in this sector)
              //   - RMS of block tilts within the module (in mrad)
              //   - Φ_block (mean η-tilt per 2×2 block inside each module)
              for (int m = 0; m < kModulesPerSector; ++m)
              {
                if (!moduleValid[m]) continue;

                const double tilt_m_rad      = moduleTilt_rad[m];
                const double rmsBlock_m_mrad = moduleRmsBlocks_rad[m] * 1.0e3;

                if (isNorth)
                {
                  // North hemisphere (sectors 0..31): module-level stats
                  moduleTiltSumNorth_rad[m]  += tilt_m_rad;
                  moduleTiltSumNorth2_rad[m] += tilt_m_rad * tilt_m_rad;
                  moduleTiltNorthNSectors[m]++;

                  moduleRmsBlocksSumNorth_mrad [m] += rmsBlock_m_mrad;
                  moduleRmsBlocksSumNorth2_mrad[m] += rmsBlock_m_mrad * rmsBlock_m_mrad;
                  moduleRmsBlocksNorthNSectors [m]++;
                }
                else
                {
                  // South hemisphere (sectors 32..63): module-level stats
                  moduleTiltSumSouth_rad[m]  += tilt_m_rad;
                  moduleTiltSumSouth2_rad[m] += tilt_m_rad * tilt_m_rad;
                  moduleTiltSouthNSectors[m]++;

                  moduleRmsBlocksSumSouth_mrad [m] += rmsBlock_m_mrad;
                  moduleRmsBlocksSumSouth2_mrad[m] += rmsBlock_m_mrad * rmsBlock_m_mrad;
                  moduleRmsBlocksSouthNSectors [m]++;
                }

                // Per-block stats: loop over the 4 φ-blocks in this module.
                for (int b = 0; b < kBlocksPerModule; ++b)
                {
                  if (!blockValid[m][b]) continue;

                  const double tilt_block_rad = blockTilt_rad[m][b];

                  if (isNorth)
                  {
                    blockTiltSumNorth_rad [m][b] += tilt_block_rad;
                    blockTiltSumNorth2_rad[m][b] += tilt_block_rad * tilt_block_rad;
                    blockTiltNorthNSectors [m][b]++;
                  }
                  else
                  {
                    blockTiltSumSouth_rad [m][b] += tilt_block_rad;
                    blockTiltSumSouth2_rad[m][b] += tilt_block_rad * tilt_block_rad;
                    blockTiltSouthNSectors [m][b]++;
                  }
                }
              }

          } // end loop over sectors

          // Global summary across all sectors
          if (globalNSectors > 0)
          {
              const double meanTilt_rad =
                globalSumTilt_rad / static_cast<double>(globalNSectors);
              double varTilt_rad =
                globalSumTilt2_rad / static_cast<double>(globalNSectors)
                - meanTilt_rad * meanTilt_rad;
              if (varTilt_rad < 0.0) varTilt_rad = 0.0;
              const double rmsTiltSpread_rad = std::sqrt(varTilt_rad);

              const double meanTilt_mrad      = meanTilt_rad * 1.0e3;
              const double rmsTiltSpread_mrad = rmsTiltSpread_rad * 1.0e3;
              const double meanRmsBlocks_mrad =
                globalSumRmsBlocks_mrad  / static_cast<double>(globalNSectors);
              const double meanRmsModules_mrad =
                globalSumRmsModules_mrad / static_cast<double>(globalNSectors);
              const double meanRmsTowers_mrad =
                globalSumRmsTowers_mrad  / static_cast<double>(globalNSectors);
              const double meanAbsTilt_mrad =
                globalSumAbsTilt_mrad    / static_cast<double>(globalNSectors);
              const double meanNorthTilt_mrad =
                (northNSectors > 0)
                ? (northSumTilt_rad / static_cast<double>(northNSectors)) * 1.0e3
                : 0.0;
              const double meanSouthTilt_mrad =
                (southNSectors > 0)
                ? (southSumTilt_rad / static_cast<double>(southNSectors)) * 1.0e3
                : 0.0;

              std::cout << "\n" << ANSI_MAGENTA_BOLD
                        << "[dissectEtaBlocks] Global η-tilt summary over "
                        << globalNSectors << " sectors\n"
                        << ANSI_RESET;

              std::cout << "  Signed <Φ_sector>            = "
                        << std::fixed << std::setprecision(3)
                        << meanTilt_mrad << " mrad  ("
                        << std::setprecision(6) << meanTilt_rad << " rad)\n";
              std::cout << "  <|Φ_sector|>                 = "
                        << std::setprecision(3) << meanAbsTilt_mrad << " mrad\n";
              std::cout << "  North  signed <Φ_sector>     = "
                        << std::setprecision(3) << meanNorthTilt_mrad << " mrad  over "
                        << northNSectors << " sectors\n";
              std::cout << "  South  signed <Φ_sector>     = "
                        << std::setprecision(3) << meanSouthTilt_mrad << " mrad  over "
                        << southNSectors << " sectors\n";
              std::cout << "  RMS[Φ_sector - global mean]  = "
                        << std::setprecision(3)
                        << rmsTiltSpread_mrad << " mrad\n";
              std::cout << "  <RMS_blocks|s>               = "
                        << std::setprecision(3)
                        << meanRmsBlocks_mrad << " mrad\n";
              std::cout << "  <RMS_modules|s>              = "
                        << std::setprecision(3)
                        << meanRmsModules_mrad << " mrad\n";
              std::cout << "  <RMS_towers|s>               = "
                        << std::setprecision(3)
                        << meanRmsTowers_mrad << " mrad\n";

              // ------------------------------------------------------------------
              // Global per-module summary over sectors, split by hemisphere:
              //   For each module m (0..23):
              //     <Φ_module>[mrad] ± σ_sector(Φ_module)[mrad]
              //     <RMS_block|mod>[mrad] ± σ_sector(RMS_block|mod)[mrad]
              //     and per-block means <Φ_block_b>[mrad],  b = 1..4.
              // ------------------------------------------------------------------
              // Compact column widths for a single "value ± error" per column
              const int wModIdx = 4;
              const int wPhiCol = 24;  // Φ_module ± σ
              const int wRmsCol = 24;  // RMS_block|mod ± σ
              const int wBlkCol = 22;  // Φ_blockX ± σ

              // Total width for separator line (approximate but generous)
              const int wTotal =
                  2 + wModIdx + 1 + wPhiCol + 1 + wRmsCol + 1 +
                  4 * wBlkCol + 3;

              // Helper lambda to print one table (North or South), including block-averaged tilts.
              // Now prints "value ± σ" and colors block columns:
              //   block 1 → dark blue, block 2 → dark red, block 3 → dark green, block 4 → default.
              auto printModuleSummary = [&](const char* label,
                                            const double* tiltSum_rad,
                                            const double* tiltSum2_rad,
                                            const int*    tiltNSectors,
                                            const double* rmsSum_mrad,
                                            const double* rmsSum2_mrad,
                                            const int*    rmsNSectors,
                                            double        blockTiltSum_rad[][kBlocksPerModule],
                                            double        blockTiltSum2_rad[][kBlocksPerModule],
                                            int           blockTiltNSectors[][kBlocksPerModule])
              {
                // ANSI colors for block columns
                const char* ANSI_DARK_BLUE  = "\033[34m";
                const char* ANSI_DARK_RED   = "\033[31m";
                const char* ANSI_DARK_GREEN = "\033[32m";

                std::cout << "\n" << ANSI_CYAN_BOLD
                          << "  Module-wise global summary over sectors (" << label << ")\n"
                          << ANSI_RESET;

                // Header row (compact labels with ±σ baked in)
                std::cout << "  "
                          << std::left  << std::setw(wModIdx) << "mod" << " "
                          << std::right << std::setw(wPhiCol) << "Φ_mod±σ [mrad]" << " "
                          << std::right << std::setw(wRmsCol) << "RMS_blk±σ [mrad]" << " "
                          << std::right << std::setw(wBlkCol) << "Φ_b1±σ [mrad]" << " "
                          << std::right << std::setw(wBlkCol) << "Φ_b2±σ [mrad]" << " "
                          << std::right << std::setw(wBlkCol) << "Φ_b3±σ [mrad]" << " "
                          << std::right << std::setw(wBlkCol) << "Φ_b4±σ [mrad]"
                          << "\n";

                // Separator line
                std::cout << "  " << std::string(wTotal, '-') << "\n";

                // Data rows
                for (int m = 0; m < kModulesPerSector; ++m)
                {
                  const int nTilt = tiltNSectors[m];
                  const int nRms  = rmsNSectors[m];

                  if (nTilt <= 0 && nRms <= 0) continue;

                  double meanTilt_mrad      = 0.0;
                  double sigmaTilt_mrad     = 0.0;
                  double meanRmsBlock_mrad  = 0.0;
                  double sigmaRmsBlock_mrad = 0.0;

                  // Module-level Φ_module mean and σ across sectors
                  if (nTilt > 0)
                  {
                    const double meanTilt_rad_mod =
                      tiltSum_rad[m] / static_cast<double>(nTilt);
                    double varTilt_rad_mod =
                      tiltSum2_rad[m] / static_cast<double>(nTilt)
                      - meanTilt_rad_mod * meanTilt_rad_mod;
                    if (varTilt_rad_mod < 0.0) varTilt_rad_mod = 0.0;

                    meanTilt_mrad  = meanTilt_rad_mod * 1.0e3;
                    sigmaTilt_mrad = std::sqrt(varTilt_rad_mod) * 1.0e3;
                  }

                  // Module-level RMS_block|mod mean and σ across sectors
                  if (nRms > 0)
                  {
                    meanRmsBlock_mrad =
                      rmsSum_mrad[m] / static_cast<double>(nRms);
                    double varRms_mrad_mod =
                      rmsSum2_mrad[m] / static_cast<double>(nRms)
                      - meanRmsBlock_mrad * meanRmsBlock_mrad;
                    if (varRms_mrad_mod < 0.0) varRms_mrad_mod = 0.0;

                    sigmaRmsBlock_mrad = std::sqrt(varRms_mrad_mod);
                  }

                  // Format "value ± σ" strings for the module-level columns
                  std::string phiColStr =
                    (nTilt > 0)
                    ? std::string(Form("%.3f ± %.3f", meanTilt_mrad, sigmaTilt_mrad))
                    : std::string("-");
                  std::string rmsColStr =
                    (nRms > 0)
                    ? std::string(Form("%.3f ± %.3f", meanRmsBlock_mrad, sigmaRmsBlock_mrad))
                    : std::string("-");

                  // Per-block hemisphere-averaged tilts with errors for this module.
                  std::string blkColStr[kBlocksPerModule];
                  for (int b = 0; b < kBlocksPerModule; ++b)
                  {
                    const int nB = blockTiltNSectors[m][b];
                    if (nB > 0)
                    {
                      const double meanTilt_block_rad =
                        blockTiltSum_rad[m][b] / static_cast<double>(nB);
                      double varTilt_block_rad =
                        blockTiltSum2_rad[m][b] / static_cast<double>(nB)
                        - meanTilt_block_rad * meanTilt_block_rad;
                      if (varTilt_block_rad < 0.0) varTilt_block_rad = 0.0;

                      const double meanTilt_block_mrad  = meanTilt_block_rad * 1.0e3;
                      const double sigmaTilt_block_mrad = std::sqrt(varTilt_block_rad) * 1.0e3;

                      blkColStr[b] = std::string(
                        Form("%.3f ± %.3f", meanTilt_block_mrad, sigmaTilt_block_mrad));
                    }
                    else
                    {
                      blkColStr[b] = "-";
                    }
                  }

                  // Print the row: module index, Φ_mod±σ, RMS_blk±σ, then colored blocks
                  std::cout << "  "
                            << std::right << std::setw(wModIdx) << (m + 1) << " "
                            << std::right << std::setw(wPhiCol) << phiColStr << " "
                            << std::right << std::setw(wRmsCol) << rmsColStr << " ";

                  // Block 1 → dark blue
                  std::cout << ANSI_DARK_BLUE
                            << std::right << std::setw(wBlkCol) << blkColStr[0]
                            << ANSI_RESET << " ";

                  // Block 2 → dark red
                  std::cout << ANSI_DARK_RED
                            << std::right << std::setw(wBlkCol) << blkColStr[1]
                            << ANSI_RESET << " ";

                  // Block 3 → dark green
                  std::cout << ANSI_DARK_GREEN
                            << std::right << std::setw(wBlkCol) << blkColStr[2]
                            << ANSI_RESET << " ";

                  // Block 4 → default color
                  std::cout << std::right << std::setw(wBlkCol) << blkColStr[3]
                            << "\n";
                }
              };

              // Print separate tables for North (0..31) and South (32..63)
              printModuleSummary(
                  "North hemisphere (sectors 0-31)",
                  moduleTiltSumNorth_rad,
                  moduleTiltSumNorth2_rad,
                  moduleTiltNorthNSectors,
                  moduleRmsBlocksSumNorth_mrad,
                  moduleRmsBlocksSumNorth2_mrad,
                  moduleRmsBlocksNorthNSectors,
                  blockTiltSumNorth_rad,
                  blockTiltSumNorth2_rad,
                  blockTiltNorthNSectors);

              printModuleSummary(
                  "South hemisphere (sectors 32-63)",
                  moduleTiltSumSouth_rad,
                  moduleTiltSumSouth2_rad,
                  moduleTiltSouthNSectors,
                  moduleRmsBlocksSumSouth_mrad,
                  moduleRmsBlocksSumSouth2_mrad,
                  moduleRmsBlocksSouthNSectors,
                  blockTiltSumSouth_rad,
                  blockTiltSumSouth2_rad,
                  blockTiltSouthNSectors);

              // ------------------------------------------------------------------
              // Hemisphere-averaged block-tilt summary plots with error bars
              //   - One PNG for North sectors 0-31
              //   - One PNG for South sectors 32-63
              //
              // For each module m and block b=0..3, we plot:
              //   <Φ_block>[mrad] with vertical error bar σ_sector(Φ_block)[mrad].
              // Layout (axes, dashed vertical lines, legend) matches the per-sector plots.
              // ------------------------------------------------------------------
              auto makeBlockSummaryGraph = [&](const char* hemiLabel,
                                               const char* fileLabel,
                                               double      tiltSum_rad[][kBlocksPerModule],
                                               double      tiltSum2_rad[][kBlocksPerModule],
                                               int         tiltNSectors[][kBlocksPerModule])
              {
                const int nMod   = kModulesPerSector;
                const int nBlock = kBlocksPerModule;

                // One TGraphErrors worth of arrays per block.
                double xBlock [kBlocksPerModule][kModulesPerSector];
                double yBlock [kBlocksPerModule][kModulesPerSector];
                double exBlock[kBlocksPerModule][kModulesPerSector];
                double eyBlock[kBlocksPerModule][kModulesPerSector];
                int    nPoints[kBlocksPerModule] = {0, 0, 0, 0};

                double yMin = 0.0;
                double yMax = 0.0;
                bool   firstPoint = true;

                const double dx = 0.20;  // same intra-module offsets as per-sector plots

                // Build per-block mean and RMS across sectors.
                for (int m = 0; m < nMod; ++m)
                {
                  const double moduleCenterX = static_cast<double>(m) + 0.5;

                  for (int b = 0; b < nBlock; ++b)
                  {
                    const int nTilt = tiltNSectors[m][b];
                    if (nTilt <= 0) continue;

                    const double meanTilt_rad =
                      tiltSum_rad[m][b] / static_cast<double>(nTilt);
                    double varTilt_rad =
                      tiltSum2_rad[m][b] / static_cast<double>(nTilt)
                      - meanTilt_rad * meanTilt_rad;
                    if (varTilt_rad < 0.0) varTilt_rad = 0.0;

                    const double meanTilt_mrad  = meanTilt_rad * 1.0e3;
                    const double sigmaTilt_mrad = std::sqrt(varTilt_rad) * 1.0e3;

                    const int idx = nPoints[b];
                    const double x = moduleCenterX + (static_cast<double>(b) - 1.5) * dx;

                    xBlock [b][idx] = x;
                    yBlock [b][idx] = meanTilt_mrad;
                    exBlock[b][idx] = 0.0;
                    eyBlock[b][idx] = sigmaTilt_mrad;
                    nPoints[b]++;

                    if (firstPoint)
                    {
                      yMin = yMax = meanTilt_mrad;
                      firstPoint  = false;
                    }
                    else
                    {
                      if (meanTilt_mrad < yMin) yMin = meanTilt_mrad;
                      if (meanTilt_mrad > yMax) yMax = meanTilt_mrad;
                    }
                  }
                }

                int totalPoints = 0;
                for (int b = 0; b < nBlock; ++b)
                {
                  totalPoints += nPoints[b];
                }
                if (totalPoints == 0) return;

                if (gSystem)
                {
                  gSystem->mkdir("sector_module_summariesPNGs", true);
                }

                std::string cname  = std::string("cSummaryBlockTilt_") + fileLabel;
                std::string ctitle = std::string("Mean #theta tilt per block (")
                                   + hemiLabel + ")";

                TCanvas* cSum = new TCanvas(cname.c_str(), ctitle.c_str(), 1600, 450);
                cSum->SetMargin(0.06, 0.02, 0.14, 0.08);
                cSum->SetGrid(0, 0);

                // Choose first block with data as frame for axes.
                int frameBlock = -1;
                for (int b = 0; b < nBlock; ++b)
                {
                  if (nPoints[b] > 0)
                  {
                    frameBlock = b;
                    break;
                  }
                }
                if (frameBlock < 0) return;

                  TGraphErrors* gFrame = new TGraphErrors(
                      nPoints[frameBlock],
                      xBlock [frameBlock],
                      yBlock [frameBlock],
                      exBlock[frameBlock],
                      eyBlock[frameBlock]);

                  std::string gTitle = std::string("Mean #Delta #theta_{block} (m, b) [mrad] per block (")
                                     + hemiLabel + ");Module Index;#Delta #theta_{block} (m, b) [mrad]";
                  gFrame->SetTitle(gTitle.c_str());
                  gFrame->SetMarkerStyle(20);
                  gFrame->SetMarkerSize(1.0);
                  gFrame->SetMarkerColor(kBlack);
                  gFrame->SetLineColor(0);  // markers + error bars only
                  gFrame->Draw("AP");

                  // X axis (modules 1..24 as bins)
                  const double xMin = 0.0;
                  const double xMax = static_cast<double>(kModulesPerSector);

                  gFrame->GetXaxis()->SetLimits(xMin, xMax);
                  gFrame->GetXaxis()->SetNdivisions(kModulesPerSector, 0, 0);
                  gFrame->GetXaxis()->SetTickLength(0.02);
                  gFrame->GetXaxis()->CenterTitle(false);
                  gFrame->GetXaxis()->SetTitle("Module Index");
                  gFrame->GetXaxis()->SetTitleSize(0.045);
                  gFrame->GetXaxis()->SetTitleOffset(1.0);
                  gFrame->GetXaxis()->SetLabelSize(0.0); // hide numeric labels

                  // Y axis range with padding
                  double span = yMax - yMin;
                  double pad  = (span > 0.0) ? 0.15 * span : 1.0;
                  const double yLow  = yMin - pad;
                  const double yHigh = yMax + pad;

                  gFrame->GetYaxis()->SetTitle("#Delta #theta_{block} (m, b) [mrad]");
                  gFrame->GetYaxis()->CenterTitle(true);
                  gFrame->GetYaxis()->SetTitleSize(0.050);
                  gFrame->GetYaxis()->SetTitleOffset(0.95);
                  gFrame->GetYaxis()->SetLabelSize(0.035);
                  gFrame->GetYaxis()->SetRangeUser(yLow, yHigh);


                // Vertical dashed lines at module edges
                for (int edge = 0; edge <= kModulesPerSector; ++edge)
                {
                  const double xEdge = static_cast<double>(edge);
                  TLine* l = new TLine(xEdge, yLow, xEdge, yHigh);
                  l->SetLineStyle(2);
                  l->SetLineWidth(1);
                  l->SetLineColor(kGray + 1);
                  l->Draw("SAME");
                }

                // Draw per-block graphs on top of the frame.
                int blockColor[kBlocksPerModule] = { kBlue, kRed, kGreen + 2, kBlack };
                TGraphErrors* gBlock[kBlocksPerModule] = { nullptr, nullptr, nullptr, nullptr };

                for (int b = 0; b < nBlock; ++b)
                {
                  if (nPoints[b] <= 0) continue;

                  TGraphErrors* gB = new TGraphErrors(
                      nPoints[b],
                      xBlock [b],
                      yBlock [b],
                      exBlock[b],
                      eyBlock[b]);

                  gB->SetMarkerStyle(20);
                  gB->SetMarkerSize(1.1);
                  gB->SetMarkerColor(blockColor[b]);
                  gB->SetLineColor(0); // markers + error bars only
                  gB->Draw("P SAME");

                  gBlock[b] = gB;
                }

                // Legend: same style/placement as per-sector plots
                double x1, y1, x2, y2;
                if (std::string(fileLabel) == "sectors_0_31")
                {
                  // North → bottom-right
                  x1 = 0.78; y1 = 0.30; x2 = 0.97; y2 = 0.52;
                }
                else
                {
                  // South → top-left
                  x1 = 0.10; y1 = 0.66; x2 = 0.30; y2 = 0.88;
                }

                TLegend* leg = new TLegend(x1, y1, x2, y2);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.040);
                leg->SetTextAlign(12);
                leg->SetHeader("#bf{Blocks within Modules}", "L");

                for (int b = 0; b < nBlock; ++b)
                {
                  if (!gBlock[b]) continue;
                  leg->AddEntry(gBlock[b], Form("Block %d", b + 1), "P");
                }
                leg->Draw();

                // Dummy module labels 1..24 centred in each bin
                TLatex latex;
                latex.SetTextAlign(22);
                latex.SetTextSize(0.030);
                const double yLabel = yLow + 0.05 * (yHigh - yLow);

                for (int m = 0; m < kModulesPerSector; ++m)
                {
                  const double xCenter = static_cast<double>(m) + 0.5;
                  latex.DrawLatex(xCenter, yLabel, Form("%d", m + 1));
                }

                std::string outName = std::string("sector_module_summariesPNGs/summary_sector_")
                                      + fileLabel + ".png";
                cSum->SaveAs(outName.c_str());
              };

              // North-hemisphere block averages: sectors 0-31
              makeBlockSummaryGraph(
                  "North hemisphere (sectors 0-31)",
                  "sectors_0_31",
                  blockTiltSumNorth_rad,
                  blockTiltSumNorth2_rad,
                  blockTiltNorthNSectors);

              // South-hemisphere block averages: sectors 32-63
              makeBlockSummaryGraph(
                  "South hemisphere (sectors 32-63)",
                  "sectors_32_63",
                  blockTiltSumSouth_rad,
                  blockTiltSumSouth2_rad,
                  blockTiltSouthNSectors);

              // Existing per-sector table
              const int wSecT = 4;
              const int wHemi = 4;
              const int wPhi  = 6;
              const int wTilt = 18;
              const int wRmsB = 18;
              const int wRmsM = 18;
              const int wRmsT = 18;


              std::cout << "\n" << ANSI_CYAN_BOLD
                        << "  sec hemi phiIdx  Φ_sector[mrad]"
                        << "  RMS_blocks|s[mrad]  RMS_modules|s[mrad]"
                        << "  RMS_towers|s[mrad]"
                        << ANSI_RESET << "\n";

              for (int sector = 0; sector < kSectorsTotal; ++sector)
              {
                if (!sectorValid[sector]) continue;
                const char* hemi = sectorIsNorth[sector] ? "N" : "S";

                const double tilt_mrad_sec = sectorTilt_rad[sector] * 1.0e3;
                const char*  tiltColor     = (tilt_mrad_sec >= 0.0) ? ANSI_GREEN : ANSI_RED_BOLD;

                std::cout << "  "
                          << std::right << std::setw(wSecT) << sector << " "
                          << std::left  << std::setw(wHemi) << hemi   << " "
                          << std::right << std::setw(wPhi)  << sectorPhiIdx[sector] << " ";

                std::cout << tiltColor
                          << std::right << std::setw(wTilt)
                          << std::setprecision(3) << tilt_mrad_sec
                          << ANSI_RESET << " ";

                std::cout << std::right << std::setw(wRmsB)
                          << std::setprecision(3) << sectorRmsBlocks_rad[sector] * 1.0e3 << " "
                          << std::right << std::setw(wRmsM)
                          << std::setprecision(3) << sectorRmsModules_rad[sector] * 1.0e3 << " "
                          << std::right << std::setw(wRmsT)
                          << std::setprecision(3) << sectorRmsTowers_rad[sector] * 1.0e3
                          << "\n";
              }
           }
        }
      else
      {
        std::cout << "[INFO] scanMode='" << scanMode
                  << "' → nothing to do (use 'ring', 'phiScanFromCenter', "
                  << "'philine', 'etaZvGrid', 'dissectPhiSectors', or 'dissectEtaBlocks').\n";
    }
    // After scan: print mechanical-incidence QA summary (if any incidence calls
    // were made).  The tolerances (in radians) are the acceptable maxima for
    // |αφ_mech| and |αη_mech| over the whole scan.
    const double tolAlphaPhi = 5.0e-3;  // ≈ 0.29°
    const double tolAlphaEta = 5.0e-3;  // ≈ 0.29°

    if (incidenceEverCalled)
    {
      bemc->PrintMechIncidenceQASummary(tolAlphaPhi, tolAlphaEta);
    }
    else
    {
      std::cout << "[MECH_QA] No mechanical-incidence calls were made in this scan; "
                << "skipping incidence summary.\n";
    }



  // ────────────────────────────────────────────────────────────────────────────
  // 8) Wrap up - no event loop was run in this macro
  // ────────────────────────────────────────────────────────────────────────────
  se->End();
  std::cout << "[INFO] Fun4all_ParticleGunTest finished (no Run-24 events were processed).\n";
  std::cout << "============================================================\n";
  return;
}
