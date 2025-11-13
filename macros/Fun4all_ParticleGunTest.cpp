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

      const double zc   = tg->get_center_z();   // face-center z
      const double phic = tg->get_phi();        // tower-center φ

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
        const double d = std::fabs(TVector2::Phi_mpi_pi(tg->get_phi() - targetPhi));
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
      // Rigorous grid:
      //  * pick a single tower (ietaRing, iphi*) along a φ line (phiTarget)
      //  * for that tower, scan η within the block: y = ietaRing + yOff
      //  * for each yOff, scan z_vtx = Zc + zRel
      //  * use detector-η from geometry (eta_det) as label (no vertex dependence)
      // ------------------------------------------------------------------
      const int ieta = std::clamp(ietaRing, 0, geo->get_etabins() - 1);
      const double target = phiTarget;  // radians

      const int iphi = pick_iphi_closest_to(ieta, target);
      const RawTowerDefs::keytype key =
        RawTowerDefs::encode_towerid(geo->get_calorimeter_id(), ieta, iphi);
      RawTowerGeom* tg = geo->get_tower_geometry(key);
      if (!tg)
      {
        std::cerr << "[etaZvGrid] ERROR: no geometry for (ieta,iphi)=("
                  << ieta << "," << iphi << ")\n";
      }
      else
      {
        const double Cx  = tg->get_center_x();
        const double Cy  = tg->get_center_y();
        const double Cz  = tg->get_center_z();
        const double Rc  = std::hypot(Cx, Cy);
        const double phi = tg->get_phi();
        const double eta_det = (Rc > 0.0) ? std::asinh(Cz / Rc) : 0.0;

        std::cout << "# [etaZvGrid] tower (ieta="<<ieta<<", iphi="<<iphi<<")\n"
                  << "#  C=("<<Cx<<","<<Cy<<","<<Cz<<")  R="<<Rc
                  << "  phi_center="<<phi<<"  eta_det(center)="<<eta_det
                  << "  sandbox="<<(sandboxNoTilt?"1":"0")<<"\n";

        std::cout << "# yOff, zRel, z_vtx, eta_det(label), "
                  << "alpha_phi, alpha_eta, cos_phi, cos_eta\n";

        // η offsets inside this tower’s block: y = ieta + yOff, yOff ∈ [-0.5,+0.5]
        const int N_ETA_SAMPLES = 11;
        const int N_Z_SAMPLES   = 13;
        const double yOff_min = -0.5;
        const double yOff_max = +0.5;
        const double yOff_step = (N_ETA_SAMPLES > 1)
                               ? (yOff_max - yOff_min) / double(N_ETA_SAMPLES - 1)
                               : 1.0;

        // zRel ∈ [-60,60] cm relative to Cz (so zRel=0 ⇒ z_vtx=Cz is the anchor)
        const double zRel_min = -60.0;
        const double zRel_max =  60.0;
        const double zRel_step = (N_Z_SAMPLES > 1)
                               ? (zRel_max - zRel_min) / double(N_Z_SAMPLES - 1)
                               : 1.0;

        for (int iyOff = 0; iyOff < N_ETA_SAMPLES; ++iyOff)
        {
          const double yOff = yOff_min + iyOff * yOff_step;
          const float  y_eval = static_cast<float>(ieta) + static_cast<float>(yOff);
          const float  x_eval = static_cast<float>(iphi); // stay at tower φ center

          for (int iz = 0; iz < N_Z_SAMPLES; ++iz)
          {
            const double zRel  = zRel_min + iz * zRel_step;
            const double z_vtx = Cz + zRel;

            // Set “event” vertex z; ComputeIncidenceSD uses fVz
            bemc->SetVertexZ(static_cast<float>(z_vtx));

            // Optionally fire Geant (not needed for pure geometry incidence)
            if (actuallyFireGeant)
            {
              gun->set_vtx(0.0, 0.0, z_vtx);
              gun->set_phi_range(phi, phi);
              Fun4AllServer::instance()->run(1);
            }

            float cosphi=0, coseta=0, aphis=0, aetas=0;
            const bool ok = bemc->ComputeIncidenceSD(
                static_cast<float>(EGeV),
                x_eval, y_eval,
                cosphi, coseta, aphis, aetas);

            // Each row is a clean, organized snapshot; full internal vectors
            // (C,F,V,un,uphi,ueta,pF,pn,pph,pet) still printed by the tracer if dbg>0.
            std::cout << std::fixed
                      << yOff << ", "
                      << zRel << ", "
                      << z_vtx << ", "
                      << eta_det << ", "
                      << aphis << ", "
                      << aetas << ", "
                      << cosphi << ", "
                      << coseta << ", "
                      << (sandboxNoTilt ? 1 : 0) << ", "
                      << (ok ? 1 : 0) << "\n";
          }
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

