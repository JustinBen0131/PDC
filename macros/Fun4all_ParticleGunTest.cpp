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

#include <g4main/PHG4ParticleGenerator.h>      // particle gun (API: set_eta/phi/z/mom) [docs cited]
#include <TVector2.h>                           // phi wrapping helper
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
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawTowerGeomContainer.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/coresoftware/offline/packages/CaloBase/RawTowerDefs.h"

#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic pop
#endif

// --- Load only what we need; no PDC, no eval libs ---
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffaobjects.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcdbobjects)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libcalo_io.so)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libcalo_reco.so)
// Optional (only if actually firing Geant events; safe to keep loaded)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)

#endif // ROOT>=6

////////////////////////////////////////////////////////////
// Main macro (no outputs; terminal prints only)
////////////////////////////////////////////////////////////
void Fun4all_ParticleGunTest(
    // geometry / conditions
    int runNumber                 = 24,          // TIMESTAMP/RUNNUMBER for CDB
    // scan controls
    const std::string& scanMode   = "ring",      // "ring" | "philine" | "none"
    int ietaRing                  = 48,          // used by scanMode=="ring"
    double phiTarget              = 0.0,         // used by scanMode=="philine" (radians)
    // incidence controls
    double EGeV                   = 8.0,         // test energy used in ComputeIncidenceSD (front-face is E‑independent)
    bool sandboxNoTilt            = true,        // zero-incidence anchor mode
    int  incidenceDbgLevel        = 5,           // >=5: structured dump per call
    int  incidenceHardStopLevel   = -1,          // -1=no stop; >=0 with dbg>=this → throw after dump
    // execution controls
    bool actuallyFireGeant        = false        // if true: run 1 event per tower using the gun
)
{
  // ────────────────────────────────────────────────────────────────────────────
  // 0) Fun4All + Conditions DB
  // ────────────────────────────────────────────────────────────────────────────
  Fun4AllServer* se = Fun4AllServer::instance();
  if (!se) { std::cerr << "F4A server null\n"; return; }
  se->Verbosity(0);

  auto* rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  rc->set_uint64Flag("TIMESTAMP", runNumber);
  rc->set_IntFlag   ("RUNNUMBER", runNumber);

  // prove CDB key is reachable (debug; harmless)
  (void) CDBInterface::instance()->getUrl("EMCTOWERCALIB");

  // ────────────────────────────────────────────────────────────────────────────
  // 1) Put CEMC detailed geometry on the node tree
  // ────────────────────────────────────────────────────────────────────────────
  auto* geomMap = new CaloGeomMapping("CEMC_GeomFiller");
  geomMap->set_detector_name("CEMC");
  geomMap->set_UseDetailedGeometry(true);   // we want the 8‑vertex blocks
  geomMap->Verbosity(0);
  se->registerSubsystem(geomMap);

  // Force InitRun so the geometry nodes are materialized
  // (SubsysReco::InitRun happens on the first run() call)
  se->run(0);   // 0 events → triggers InitRun without processing data

  // ────────────────────────────────────────────────────────────────────────────
  // 2) Fetch geometry from the node tree
  // ────────────────────────────────────────────────────────────────────────────
  PHCompositeNode* topNode = Fun4AllServer::instance()->topNode();
  auto* geo =
    findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC_DETAILED");
  if (!geo)
    geo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if (!geo) {
    std::cerr << "FATAL: CEMC geometry (detailed or default) not found on node tree.\n";
    return;
  }

  // ────────────────────────────────────────────────────────────────────────────
  // 3) Configure BEmcRecCEMC with detailed tower geometry (one pass)
  // ────────────────────────────────────────────────────────────────────────────
  auto* bemc = new BEmcRecCEMC();
  bemc->SetCylindricalGeometry();
  bemc->set_UseDetailedGeometry(true);
  bemc->SetTowerThreshold(0.030f);           // keep consistent with your cluster builder
  bemc->EnablePhiTilt(false);                // explicit: no φ(E) pre-rotation here

  // copy dims then feed towers; note: ix=phi, iy=eta convention
  bemc->SetDim(geo->get_phibins(), geo->get_etabins());
  for (int ieta = 0; ieta < geo->get_etabins(); ++ieta)
  {
    for (int iphi = 0; iphi < geo->get_phibins(); ++iphi)
    {
      const RawTowerDefs::keytype key =
        RawTowerDefs::encode_towerid( geo->get_calorimeter_id(), ieta, iphi ); // CEMC
      if (auto* tg = geo->get_tower_geometry(key))
      {
        bemc->SetTowerGeometry(iphi, ieta, *tg);
      }
    }
  }
  bemc->CompleteTowerGeometry();   // derive per-face tangents once

  // Incidence tracer controls (the tiny toggles you added)
  bemc->SetIncidenceNoTiltSandbox(sandboxNoTilt);
  bemc->SetIncidenceDebugLevel(incidenceDbgLevel);
  bemc->SetIncidenceHardStopLevel(incidenceHardStopLevel);

  // ────────────────────────────────────────────────────────────────────────────
  // 4) Minimal gun (used only if actuallyFireGeant==true)
  //    (PHG4ParticleGenerator supports set_eta/phi/mom/vtx setters)
  // ────────────────────────────────────────────────────────────────────────────
  auto* gun = new PHG4ParticleGenerator("GUN");
  gun->set_name("gamma");
  gun->set_mom_range(EGeV, EGeV);
  gun->set_eta_range(0.0, 0.0);           // radial shot → αφ≈0 in sandbox
  se->registerSubsystem(gun);             // harmless if we don't run events

    // ────────────────────────────────────────────────────────────────────────────
    // 5) Helper: aim/compute (and optionally "fire") at a specific tower center
    // ────────────────────────────────────────────────────────────────────────────
    auto compute_or_fire_at = [&](int ieta, int iphi)
    {
      const RawTowerDefs::keytype key =
        RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi);
      RawTowerGeom* tg = geo->get_tower_geometry(key);
      if (!tg) return;

        const double zc   = tg->get_center_int_z();                 // front-face center z
        const double phic = std::atan2(tg->get_center_int_y(),      // φ from front-face center
                                       tg->get_center_int_x());

      // Set "event" vertex (beamline) at this tower's face z for the incidence calc
      bemc->SetVertexZ(static_cast<float>(zc));

      if (actuallyFireGeant)
      {
        // Fire exactly one gamma through this tower's center
        gun->set_vtx(0.0, 0.0, zc);
        gun->set_phi_range(phic, phic);
        Fun4AllServer::instance()->run(1);       // one event (requires a G4 stack if you later add it)
      }

      // Compute front-face incidence at the tower center (x=iphi, y=ieta)
      float cosphi=0, coseta=0, aphis=0, aetas=0;
      const bool ok = bemc->ComputeIncidenceSD(
          static_cast<float>(EGeV),
          static_cast<float>(iphi), static_cast<float>(ieta),
          cosphi, coseta, aphis, aetas);

      // One-line summary (your structured per-call tracer will also dump if dbg>0)
      std::cout << "[SCAN] ok="<<ok
                << "  ieta="<<ieta<<" iphi="<<iphi
                << "  zc="<<zc<<"  phi="<<phic
                << "  a_phi="<<aphis<<"  a_eta="<<aetas
                << "  cos_phi="<<cosphi<<"  cos_eta="<<coseta
                << "  sandbox="<<(sandboxNoTilt?"1":"0")
                << std::endl;
    };

    // Small helper: find iphi whose tower center φ is closest to target
    auto pick_iphi_closest_to = [&](int ieta, double targetPhi)->int
    {
      int best_iphi = 0;
      double best_d = 1e9;
      for (int iphi = 0; iphi < geo->get_phibins(); ++iphi)
      {
        const RawTowerDefs::keytype key =
          RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi);
        RawTowerGeom* tg = geo->get_tower_geometry(key);
        if (!tg) continue;
        const double phiC = std::atan2(tg->get_center_int_y(), tg->get_center_int_x()); // front-face φ
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
    // 6) Scan logic (no clusterizer, no PDC, no files)
    // ────────────────────────────────────────────────────────────────────────────
    if (scanMode == "ring")
    {
      const int ieta = (ietaRing >= 0 && ietaRing < geo->get_etabins())
                         ? ietaRing : geo->get_etabins()/2;
      for (int iphi = 0; iphi < geo->get_phibins(); ++iphi)
        compute_or_fire_at(ieta, iphi);
    }
    else if (scanMode == "philine")
    {
      const double target = phiTarget;  // radians
      for (int ieta = 0; ieta < geo->get_etabins(); ++ieta)
      {
        const int best_iphi = pick_iphi_closest_to(ieta, target);
        compute_or_fire_at(ieta, best_iphi);
      }
    }
    else if (scanMode == "etaZvGrid")
    {
      // ------------------------------------------------------------------
      // Eta vs vertex-z table along a single φ column:
      //  * choose one fixed iphi column using reference row ietaRing (closest to phiTarget)
      //  * for each ieta on that column:
      //      - compute detector η from tower front-face center (geometry only)
      //      - for z_vtx = 0 and 20 cm, compute shower-depth η and incidence
      //  * print: z_vtx_cm, ieta, eta_det, eta_SD, alpha_phi, alpha_eta, cos_phi, cos_eta
      // ------------------------------------------------------------------

      const double target = phiTarget;               // radians
      const int nEta = geo->get_etabins();
      const int ieta_ref = std::clamp(ietaRing, 0, nEta - 1);
      const int iphi_fixed = pick_iphi_closest_to(ieta_ref, target);

      // Header context from the reference tower
      {
        const RawTowerDefs::keytype key_ref =
          RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta_ref, iphi_fixed);
        RawTowerGeom* tg_ref = geo->get_tower_geometry(key_ref);
        if (tg_ref)
        {
            const double Cx  = tg_ref->get_center_int_x();   // front-face center
            const double Cy  = tg_ref->get_center_int_y();
            const double Cz  = tg_ref->get_center_int_z();
            const double Rc  = std::hypot(Cx, Cy);
            const double phi = std::atan2(Cy, Cx);           // front-face φ
            const double eta_det_ref = (Rc > 0.0) ? std::asinh(Cz / Rc) : 0.0;

          std::cout << "# [etaZvGrid] fixed-phi column scan\n"
                    << "#  reference tower: (ieta_ref="<<ieta_ref<<", iphi_fixed="<<iphi_fixed<<")\n"
                    << "#  C_ref=("<<Cx<<","<<Cy<<","<<Cz<<")  R_ref="<<Rc
                    << "  phi_center_ref="<<phi
                    << "  eta_det_ref="<<eta_det_ref << "\n";
        }
      }

      std::cout << "# sandbox="<<(sandboxNoTilt ? 1 : 0)<<"\n";
      std::cout << "# Columns:\n"
                << "#  z_vtx_cm, ieta, eta_det, eta_SD, alpha_phi, alpha_eta, cos_phi, cos_eta\n";

      // Two vertices: red (0 cm), blue (20 cm)
      const double zVertices[2] = {0.0, 20.0};

      for (int ieta = 0; ieta < nEta; ++ieta)
      {
        // Tower at (ieta, iphi_fixed)
        const RawTowerDefs::keytype key =
          RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi_fixed);
        RawTowerGeom* tg = geo->get_tower_geometry(key);
        if (!tg) continue;

          const double Cx  = tg->get_center_int_x();   // front-face center
          const double Cy  = tg->get_center_int_y();
          const double Cz  = tg->get_center_int_z();
          const double Rc  = std::hypot(Cx, Cy);

          // detector η from front-face geometry (vertex-independent)
          const double eta_det = (Rc > 0.0) ? std::asinh(Cz / Rc) : 0.0;

        // Keep only towers within detector acceptance [-1.1, 1.1]
        if (eta_det < -1.1 || eta_det > 1.1) continue;

        const float x_eval = static_cast<float>(iphi_fixed); // fixed φ column
        const float y_eval = static_cast<float>(ieta);       // this η row

        for (double z_vtx : zVertices)
        {
          // Shower-depth η: ray from (0,0,z_vtx) to the tower front-face center
          const double eta_SD = (Rc > 0.0) ? std::asinh((Cz - z_vtx) / Rc) : 0.0;

          // Set event vertex z for ComputeIncidenceSD
          bemc->SetVertexZ(static_cast<float>(z_vtx));

          if (actuallyFireGeant)
          {
            const double phiCenter = tg->get_phi();
            gun->set_vtx(0.0, 0.0, z_vtx);
            gun->set_phi_range(phiCenter, phiCenter);
            Fun4AllServer::instance()->run(1);
          }

          float cosphi = 0.f, coseta = 0.f, aphis = 0.f, aetas = 0.f;
          const bool ok = bemc->ComputeIncidenceSD(
              static_cast<float>(EGeV), x_eval, y_eval,
              cosphi, coseta, aphis, aetas);

          // Compact row: same hit point for both z_vtx values at this (ieta, iphi_fixed)
          (void)ok; // printed values are finite-by-construction when ok==true
          std::cout << std::fixed
                    << z_vtx   << ", "
                    << ieta    << ", "
                    << eta_det << ", "
                    << eta_SD  << ", "
                    << aphis   << ", "
                    << aetas   << ", "
                    << cosphi  << ", "
                    << coseta  << "\n";
        }
      }
    }
    else
    {
      std::cout << "[INFO] scanMode='"<<scanMode
                << "' → nothing to do (use 'ring', 'philine', or 'etaZvGrid').\n";
    }

    // Wrap up (no outputs)
    se->End();
    std::cout << "[INFO] ParticleGunTest finished.\n";
    return;
  }

