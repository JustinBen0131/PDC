#pragma once
#if defined(__GNUC__) && !defined(__clang__)
  // real GCC – it understands the flag
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wvirtual-function-default"
#endif
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

////////////////////////////////////////////////////////////
// Standard includes
////////////////////////////////////////////////////////////
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <calobase/TowerInfoDefs.h>
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <phool/recoConsts.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <caloreco/RawClusterDeadHotMask.h>
#include <caloreco/CaloGeomMapping.h>

#include <caloreco/RawClusterBuilderTemplate.h>   // Use the built-in template
#include <caloreco/RawTowerCalibration.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerStatus.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <caloreco/RawClusterPositionCorrection.h>

#include <g4detectors/PHG4FullProjSpacalCellReco.h>
#include <g4calo/RawTowerBuilder.h>
#include <g4calo/RawTowerDigitizer.h>
#include <calowaveformsim/CaloWaveformSim.h>
#include <phparameter/PHParameterUtils.h>

#include <globalvertex/GlobalVertexReco.h>   // For global vertex, e.g. MBD
#include <globalvertex/GlobalVertex.h>       // if you need the enum VTXTYPE

#include <ffamodules/CDBInterface.h>
#include <ffaobjects/EventHeader.h>
#include <Event/Event.h>
#include <mbd/MbdReco.h>
#include <GlobalVariables.C>
#include "Calo_Calib.C"

////////////////////////////////////////////////////////////
// Possibly your custom code
////////////////////////////////////////////////////////////
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src/PositionDependentCorrection.h"

#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/BEmcRecCEMC.h"
#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic pop
#endif

R__LOAD_LIBRARY(libcdbobjects)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libffaobjects.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)
R__LOAD_LIBRARY(libcalibCaloEmc_pi0.so)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libPDC.so)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libcalo_reco.so)

#endif // if ROOT_VERSION >= 6.00.0

////////////////////////////////////////////////////////////
// The main Fun4All macro
////////////////////////////////////////////////////////////
void Fun4All_PDC(int nevents = 0,
                 const std::string &caloListFile = "dst_calo_cluster.list",
                 const std::string &g4HitsListFile = "g4hits.list",
                 const string &outFileName       = "")
{

  const bool useCustomClusterizer = true;  // <== Flip this to false to disable
    
  // 0) Setup Fun4AllServer
  Fun4AllServer *se = Fun4AllServer::instance();
  if (!se)
  {
    std::cerr << "[ERROR] Could not get Fun4AllServer instance!" << std::endl;
    gSystem->Exit(1);
  }

  int verbosity = 2; // or 100 for extremely verbose
  se->Verbosity(0);

  // 1) Setup recoConsts, enable conditions DB
  recoConsts *rc = recoConsts::instance();
  if (!rc)
  {
    std::cerr << "[ERROR] recoConsts::instance() returned nullptr" << std::endl;
    gSystem->Exit(1);
  }
  // We want the "MDC2" global tag
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  // We'll guess that the "TIMESTAMP" or run # is 21 (or parse from the file)
  rc->set_uint64Flag("TIMESTAMP", 24);
  rc->set_IntFlag("RUNNUMBER", 24);
    
  CDBInterface::instance()->getUrl("EMCTOWERCALIB");

  ////////////////////////////////////////////////////////////
  // 2) Register TWO input managers:
  //    (a) 'in0' for the calo DST list
  //    (b) 'in1' for G4Hits
  ////////////////////////////////////////////////////////////
  if (verbosity > 0)
    {
      std::cout << "[DEBUG] Registering input managers for calo='"
                << caloListFile << "'  and hits='" << g4HitsListFile << "'" << std::endl;
  }


  // (a) for the calo cluster DST
  Fun4AllInputManager* in0 = new Fun4AllDstInputManager("in0");
  if (!in0)
  {
    std::cerr << "[ERROR] Could not create Fun4AllDstInputManager for in0" << std::endl;
    gSystem->Exit(1);
  }
  // The second argument '1' means topNode level or segment, depending on your usage
  in0->AddFile(caloListFile);
  se->registerInputManager(in0);

  // (b) for G4Hits
  Fun4AllInputManager* in1 = new Fun4AllDstInputManager("in1");
  if (!in1)
  {
    std::cerr << "[ERROR] Could not create Fun4AllDstInputManager for in1" << std::endl;
    gSystem->Exit(1);
  }
  in1->AddFile(g4HitsListFile);
  se->registerInputManager(in1);


  GlobalVertexReco *gvertex = new GlobalVertexReco("GlobalVertexReco");
  gvertex->Verbosity(0);
  se->registerSubsystem(gvertex);
    
  BEmcRecCEMC* bemcPtr = new BEmcRecCEMC();
  bemcPtr->SetCylindricalGeometry();        // CEMC is a cylinder

  /* 4a) Let CaloGeomMapping put "TOWERGEOM_CEMC_DETAILED" on the node-tree */
  CaloGeomMapping* geomMap = new CaloGeomMapping("CEMC_GeomFiller");
  geomMap->set_detector_name("CEMC");
  geomMap->set_UseDetailedGeometry(true);   // we want the 8-vertex blocks
  geomMap->Verbosity(0);
  se->registerSubsystem(geomMap);           // register *before* anything that uses it

  /* 4b) AFTER CaloGeomMapping is in place, copy the geometry into BEmcRec */
  {
      PHCompositeNode* topNode = Fun4AllServer::instance()->topNode();

      // --- 1) fetch the container written by CaloGeomMapping
      auto* geo = findNode::getClass<RawTowerGeomContainer>(topNode,
                                                            "TOWERGEOM_CEMC_DETAILED");
      if (!geo)
          geo = findNode::getClass<RawTowerGeomContainer>(topNode,
                                                          "TOWERGEOM_CEMC");   // fallback
      if (!geo)
      {
        std::cerr << "### FATAL: CEMC geometry not found on the node tree – abort\n";
        gSystem->Exit(1);
      }

      // --- 2) configure BEmcRec for detailed geometry
      bemcPtr->set_UseDetailedGeometry(true);
      bemcPtr->SetDim(geo->get_phibins(),      // Nx  (φ-bins)
                      geo->get_etabins());     // Ny  (η-bins)

      // --- 3) feed every tower once (ix = φ, iy = η !)
      for (int ieta = 0; ieta < geo->get_etabins(); ++ieta)
      {
        for (int iphi = 0; iphi < geo->get_phibins(); ++iphi)
        {
          const RawTowerDefs::keytype key =
            RawTowerDefs::encode_towerid( geo->get_calorimeter_id(), ieta, iphi );
          if (auto* tg = geo->get_tower_geometry(key))
          {
            bemcPtr->SetTowerGeometry(iphi, ieta, *tg);
          }
        }
      }
      bemcPtr->CompleteTowerGeometry();   // derive dX/dY/dZ once
  }
    
  ////////////////////////////////////////////////////////////
  // 6) Our PositionDependentCorrection code
  ////////////////////////////////////////////////////////////
  const bool isSimulation =
          !g4HitsListFile.empty() &&
          g4HitsListFile != "g4hits.list" &&       // the default placeholder
          gSystem->AccessPathName( g4HitsListFile.c_str(), kReadPermission ) == kFALSE;

    
  string finalOut = outFileName.empty()
                    ? (isSimulation ? "output_PositionDep_sim.root"
                                    : "output_PositionDep_data.root")
                    : outFileName;
    
    
  if (useCustomClusterizer)
    {
      std::cout << "Building clusters (using custom RawClusterBuilderTemplate)" << std::endl;

      RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
      ClusterBuilder->Detector("CEMC");
      ClusterBuilder->set_threshold_energy(0.030f);  // threshold in GeV

      std::string emc_prof = getenv("CALIBRATIONROOT");
      emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
      ClusterBuilder->LoadProfile(emc_prof);

      ClusterBuilder->set_UseTowerInfo(1);  // Use TowerInfo
      ClusterBuilder->Verbosity(1);

      se->registerSubsystem(ClusterBuilder);
    }
    else
    {
      std::cout << "Skipping custom cluster builder (useCustomClusterizer = false)" << std::endl;
  }

    
  // ========================================================================
  //  Position-dependent π0-correction module
  //  ----------------------------------------------------------------------
  //  How to drive the two most common switches from your macro:
  //
  //      • Simulation / Data
  //          pdc->setIsSimulation(true);   // Monte-Carlo
  //          pdc->setIsSimulation(false);  // Real data
  //
  //      • π0-mass windows (only after you have the pass-1 μ/σ table!)
  //          pdc->setMassFitsDone(true);   // enable mass-window cuts (pass-2)
  //          pdc->setMassFitsDone(false);  // disable  (default / pass-1)
  // ========================================================================

  auto* pdc = new PositionDependentCorrection("PositionDepCorr", finalOut);
  /*------------------------------------------------------------
      1)  Mandatory hooks and global switches
  ------------------------------------------------------------*/
  pdc->setBEmcRec(bemcPtr);          // lead-tower finder from the clusteriser
  pdc->setIsSimulation(isSimulation);
  pdc->UseSurveyGeometry(true);      // load barrel-tilt from CDB (recommended)
  /*------------------------------------------------------------
      2)  π0-mass-window support
    ------------------------------------------------------------*/
  // Keep *false* for the first (pass-1) job that *produces* the μ/σ table.
  // Set *true* for the second (pass-2) job that *consumes* the table and
  // applies the slice-dependent mass cut.
  pdc->setMassFitsDone(false);
  /*------------------------------------------------------------
      4)  Energy-binning mode
  ------------------------------------------------------------*/
  pdc->setBinningMode(PositionDependentCorrection::EBinningMode::kRange);

  pdc->Verbosity(0);                 // 0 = silent → raise for debugging

  // Finally register the module with Fun4All
  se->registerSubsystem(pdc);

  std::cout << "[INFO] Running Fun4All with nevents = " << nevents << std::endl;
  se->run(nevents);

  se->End();


  std::cout << "[INFO] Done. Exiting macro." << std::endl;
  return;
}
