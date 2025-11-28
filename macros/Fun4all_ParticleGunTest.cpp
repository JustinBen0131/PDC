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
#include <cassert>
#include <cmath>
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
  // ────────────────────────────────────────────────────────────────────────────
  // 0) High-level run summary (no events are processed unless *you* call se->run)
  // ────────────────────────────────────────────────────────────────────────────
  std::cout << "============================================================\n"
            << " Fun4all_ParticleGunTest – geometry / incidence diagnostic\n"
            << "------------------------------------------------------------\n"
            << "  runNumber          = " << runNumber          << "\n"
            << "  scanMode           = " << scanMode           << "\n"
            << "  ietaRing           = " << ietaRing           << "\n"
            << "  phiTarget [rad]    = " << phiTarget          << "\n"
            << "  EGeV               = " << EGeV               << "\n"
            << "  sandboxNoTilt      = " << (sandboxNoTilt ? "true" : "false") << "\n"
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

    // Incidence mode selection:
    //  - sandboxNoTilt == true  → QA STEP 1: ideal-cylinder sandbox, FACE-only incidence.
    //  - sandboxNoTilt == false → QA STEP 2: real tilted geometry, BOTH (FACE + MECH).
    if (sandboxNoTilt)
    {
      bemc->SetIncidenceMode(BEmcRecCEMC::EIncidenceMode::FACE);
    }
    else
    {
      bemc->SetIncidenceMode(BEmcRecCEMC::EIncidenceMode::BOTH);
    }

    // Reset QA STEP 2 (FACE vs MECH) accumulators before the scan.
    bemc->ResetIncidenceFaceMechQA();



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
  bemc->SetIncidenceNoTiltSandbox(sandboxNoTilt);
  bemc->SetIncidenceDebugLevel(incidenceDbgLevel);
  bemc->SetIncidenceHardStopLevel(incidenceHardStopLevel);
  std::cout << "[INC] IncidenceNoTiltSandbox=" << (sandboxNoTilt ? "true" : "false")
            << "  IncidenceDebugLevel=" << incidenceDbgLevel
            << "  IncidenceHardStopLevel=" << incidenceHardStopLevel << "\n";

  // ────────────────────────────────────────────────────────────────────────────
  // 5) Minimal gun (used only if actuallyFireGeant==true)
  //    PHG4ParticleGenerator supports set_eta/phi/mom/vtx.
  //    *** IMPORTANT: we never call se->run() here unless actuallyFireGeant=true,
  //    so no events are processed by default. ***
  // ────────────────────────────────────────────────────────────────────────────
  auto* gun = new PHG4ParticleGenerator("GUN");
  gun->set_name("gamma");
  gun->set_mom_range(EGeV, EGeV);
  gun->set_eta_range(0.0, 0.0);              // radial shot → αφ≈0 in sandbox
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
          std::cout << "[WARN] Missing geometry for (ieta="<<ieta<<", iphi="<<iphi<<") – skipped.\n";
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
        // Source–detector pseudorapidity (eta_SD) uses the chosen vertex z_vtx:
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

        // Compute front-face incidence at the tower center (x=iphi, y=ieta)
        float cosphi = 0.f, coseta = 0.f, aphis = 0.f, aetas = 0.f;
        const bool ok = bemc->ComputeIncidenceSD(
            static_cast<float>(EGeV),
            static_cast<float>(iphi), static_cast<float>(ieta),
            cosphi, coseta, aphis, aetas);

        // Per-tower [SCAN] output only when debug is enabled
        if (incidenceDbgLevel > 0)
        {
          std::cout << "[SCAN]"
                    << " ok="        << ok
                    << "  ieta="     << ieta            // tower index in η (0…nEta-1)
                    << "  iphi="     << iphi            // tower index in φ (0…nPhi-1)
                    << "  Cx="       << Cx              // tower-center x (cm)
                    << "  Cy="       << Cy              // tower-center y (cm)
                    << "  Cz="       << Cz              // tower-center z (cm)
                    << "  R_cyl="    << Rc              // sqrt(Cx^2 + Cy^2), used in η
                    << "  z_vtx="    << z_vtx          // vertex z used in incidence calc
                    << "  eta_det_IP=" << eta_det_IP   // asinh(Cz / R_cyl)
                    << "  eta_SD="   << eta_SD         // asinh((Cz - z_vtx) / R_cyl)
                    << "  phi_center="<< phiC          // detector φ of tower center
                    << "  a_phi="    << aphis          // signed incidence in φ from BEmcRecCEMC
                    << "  a_eta="    << aetas          // signed incidence in η from BEmcRecCEMC
                    << "  cos_phi="  << cosphi         // foreshortening cos(αφ)
                    << "  cos_eta="  << coseta         // foreshortening cos(αη)
                    << "  sandbox="  << (sandboxNoTilt ? "1" : "0")
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
    // ────────────────────────────────────────────────────────────────────────────
      if (scanMode == "ring")
      {
        const int ieta = (ietaRing >= 0 && ietaRing < nEta)
                           ? ietaRing
                           : nEta/2;
        std::cout << "[MODE] ring – scanning all iphi at fixed ieta=" << ieta << "\n";
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

        std::cout << "[MODE] phiScanFromCenter – scanning φ around global (η≈0, φ≈0)"
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
                      << ieta_central << ", iphi=" << iphi << ") – skipped.\n";
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

          // Compute front-face incidence using BEmcRecCEMC (tower coords x=iphi, y=ieta_central)
          float cosphi = 0.f, coseta = 0.f, aphis = 0.f, aetas = 0.f;
          const bool ok_inc = bemc->ComputeIncidenceSD(
          static_cast<float>(EGeV),
          static_cast<float>(iphi), static_cast<float>(ieta_central),
          cosphi, coseta, aphis, aetas);

          // ------------------------------------------------------------------
          // Pretty, human-readable per-"event" block (ANSI bold red header)
          // ------------------------------------------------------------------
          std::cout << ANSI_RED_BOLD
                      << "event " << phiIndex
                      << " ---- etaTower=" << ieta_central
                      << "  etaGlobal=" << eta_det_IP
                      << "  phiTower=" << iphi
                      << "  phiGlobal=" << phi_center
                      << "\n"
                      << "         phi_nominal_fromIndex=" << phi_nominal
                      << "  dphi(center - nominal)=" << dphi_index << " rad"
                      << "\n"
                      << ANSI_RESET;

          // Source vertex (always at IP for this diagnostic)
          std::cout << "  src_vertex: (x=0, y=0, z=" << z_vtx << ")\n";

          // Tower center + η information
          std::cout << "  tower_center: "
                       << "(Cx=" << Cx << ", Cy=" << Cy << ", Cz=" << Cz << ") "
                       << "R_cyl=" << Rc << "  eta_det_IP=" << eta_det_IP << "\n";


          // Source vertex (always at IP for this diagnostic)
          std::cout << "  src_vertex: (x=0, y=0, z=" << z_vtx << ")\n";

          // Tower center + η information
          std::cout << "  tower_center: "
                    << "(Cx=" << Cx << ", Cy=" << Cy << ", Cz=" << Cz << ") "
                    << "R_cyl=" << Rc << "  eta_det_IP=" << eta_det_IP << "\n";

          // Incidence summary (ComputeIncidenceSD also prints [INC] lines if dbg>0)
          std::cout << "  incidence: ok=" << ok_inc
                    << "  a_phi=" << aphis
                    << "  a_eta=" << aetas
                    << "  cos_phi=" << cosphi
                    << "  cos_eta=" << coseta
                    << "  sandbox=" << (sandboxNoTilt ? "1" : "0")
                    << "\n\n";
        }
      }
      else if (scanMode == "philine")
      {
        std::cout << "[MODE] philine – scanning all ieta along φ≈" << phiTarget << " rad\n";
        const double target = phiTarget;
        for (int ieta = 0; ieta < nEta; ++ieta)
        {
          const int best_iphi = pick_iphi_closest_to(ieta, target);
          compute_or_fire_at(ieta, best_iphi);
        }
      }
      else if (scanMode == "etaZvGrid")
      {
        std::cout << "[MODE] etaZvGrid – scanning a fixed φ column vs (ieta, z_vtx)\n";

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

        std::cout << "# sandbox="<<(sandboxNoTilt ? "1" : "0")<<"\n";
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

            double rotX = 0.0, rotY = 0.0, rotZ = 0.0;
            if (auto* g5 = dynamic_cast<RawTowerGeomv5*>(tg))
            {
              rotX = g5->get_rotx();
              rotY = g5->get_roty();
              rotZ = g5->get_rotz();
            }
            const double phiTilt_mrad = rotX * 1.0e3;
            const double etaTilt_mrad = rotY * 1.0e3;

            bemc->SetVertexZ(static_cast<float>(z_vtx));

            if (actuallyFireGeant)
            {
              const double phiCenter = tg->get_phi();
              gun->set_vtx(0.0, 0.0, z_vtx);
              gun->set_phi_range(phiCenter, phiCenter);
              // Fun4AllServer::instance()->run(1);  // optional if you really want to simulate
            }

            float cosphi = 0.f, coseta = 0.f, aphis = 0.f, aetas = 0.f;
            const bool ok = bemc->ComputeIncidenceSD(
                static_cast<float>(EGeV), x_eval, y_eval,
                cosphi, coseta, aphis, aetas);

            (void)ok;
            std::cout << std::fixed
                      << z_vtx         << ", "
                      << ieta          << ", "
                      << eta_det       << ", "
                      << eta_SD        << ", "
                      << aphis         << ", "
                      << aetas         << ", "
                      << cosphi        << ", "
                      << coseta        << ", "
                      << phiTilt_mrad  << ", "
                      << etaTilt_mrad  << ", "
                      << rotZ * 1.0e3  << "  # rotZ_mrad\n";
          }
        }
      }
      else
      {
        std::cout << "[INFO] scanMode='"<<scanMode
                  << "' → nothing to do (use 'ring', 'phiScanFromCenter', 'philine', or 'etaZvGrid').\n";
    }


  // After scan: print Step-1 sandbox QA summary (if any incidence calls were made)
  // For the ideal-cylinder test we expect STATUS=PASS and nFACE≈0 == nCalls.
    // After scan: print QA summary appropriate to this run.
    //
    //  - sandboxNoTilt == true  → QA STEP 1 (ideal-cylinder sandbox, FACE,noTilt)
    //  - sandboxNoTilt == false → QA STEP 2 (real geometry, FACE vs MECH)
    if (sandboxNoTilt)
    {
      bemc->PrintIncidenceSandboxQASummary(1.0e-3, 1.0e-3);
    }
    else
    {
      // tol|Δa| is the "ok" scale for FACE vs MECH disagreement, ~0.3° here.
      bemc->PrintIncidenceFaceMechQASummary(5.0e-3);
    }

  // ────────────────────────────────────────────────────────────────────────────
  // 8) Wrap up – no event loop was run in this macro
  // ────────────────────────────────────────────────────────────────────────────
  se->End();
  std::cout << "[INFO] Fun4all_ParticleGunTest finished (no Run-24 events were processed).\n";
  std::cout << "============================================================\n";
  return;
}
