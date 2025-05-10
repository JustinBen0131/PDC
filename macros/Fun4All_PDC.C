#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

////////////////////////////////////////////////////////////
// Standard includes
////////////////////////////////////////////////////////////
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <caloreco/CaloTowerStatus.h>
#include <phool/recoConsts.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <caloreco/RawClusterBuilderTemplate.h>   // Use the built-in template
#include <caloreco/RawTowerCalibration.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/RawClusterPositionCorrection.h>

#include <g4detectors/PHG4FullProjSpacalCellReco.h>
#include <g4calo/RawTowerBuilder.h>
#include <g4calo/RawTowerDigitizer.h>
#include <calowaveformsim/CaloWaveformSim.h>

#include <globalvertex/GlobalVertexReco.h>   // For global vertex, e.g. MBD
#include <globalvertex/GlobalVertex.h>       // if you need the enum VTXTYPE

#include <ffamodules/CDBInterface.h>
#include <ffaobjects/EventHeader.h>
#include <Event/Event.h>
#include <mbd/MbdReco.h>
#include "Calo_Calib.C"

////////////////////////////////////////////////////////////
// Possibly your custom code
////////////////////////////////////////////////////////////
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src/PositionDependentCorrection.h"

R__LOAD_LIBRARY(libcdbobjects)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libffaobjects.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)
R__LOAD_LIBRARY(libcalibCaloEmc_pi0.so)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libPDC.so)

////////////////////////////////////////////////////////////
// The main Fun4All macro
////////////////////////////////////////////////////////////
void Fun4All_PDC(int nevents = 0,
                 const std::string &caloListFile = "dst_calo_cluster.list",
                 const std::string &g4HitsListFile = "g4hits.list",
                 const string &outFileName       = "")
{
  // 0) Setup Fun4AllServer
  Fun4AllServer *se = Fun4AllServer::instance();
  if (!se)
  {
    std::cerr << "[ERROR] Could not get Fun4AllServer instance!" << std::endl;
    gSystem->Exit(1);
  }

  int verbosity = 0; // or 100 for extremely verbose
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
  rc->set_uint64Flag("TIMESTAMP", 21);
  rc->set_IntFlag("RUNNUMBER", 21);

  // Allow the conditions DB
  // (some older code might do 'Enable::CDB = true;' or similar)
  CDBInterface::instance()->Verbosity(0);

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
  in0->AddListFile(caloListFile, 1);
  se->registerInputManager(in0);

  // (b) for G4Hits
  Fun4AllInputManager* in1 = new Fun4AllDstInputManager("in1");
  if (!in1)
  {
    std::cerr << "[ERROR] Could not create Fun4AllDstInputManager for in1" << std::endl;
    gSystem->Exit(1);
  }
  in1->AddListFile(g4HitsListFile, 1);
  se->registerInputManager(in1);

  // If you want, you can unify them into a single manager that merges them,
  // but here we are explicitly using two. The next step might be to handle
  // run # mismatch if your G4Hits has a different run number, etc.

  ////////////////////////////////////////////////////////////
  // 3) Reconstruct global vertex for MBD (for jets)
  ////////////////////////////////////////////////////////////
  std::cout << "[DEBUG] Setting up GlobalVertexReco for MBD-based vertex" << std::endl;
  GlobalVertexReco *gvertex = new GlobalVertexReco("GlobalVertexReco");
  gvertex->Verbosity(0);
  se->registerSubsystem(gvertex);

  ////////////////////////////////////////////////////////////
  // 4) Pipeline: from G4Hits => Calo cells => Raw towers => ...
  //    same as your old steps, but we reference the built-in
  //    RawClusterBuilderTemplate, not a custom MyRawClusterBuilderTemplate.
  ////////////////////////////////////////////////////////////
  bool doG4Reco = true; // presumably we want to convert from G4Hits if present
  if(doG4Reco)
  {
    PHG4FullProjSpacalCellReco *cemc_cells
      = new PHG4FullProjSpacalCellReco("CEMCCellReco");
    cemc_cells->Detector("CEMC");
    cemc_cells->Verbosity(0);
    se->registerSubsystem(cemc_cells);

    RawTowerBuilder *towerBuilder = new RawTowerBuilder("RawTowerBuilderCEMC");
    towerBuilder->Detector("CEMC");
    towerBuilder->set_sim_tower_node_prefix("SIM");
    towerBuilder->Verbosity(0);
    se->registerSubsystem(towerBuilder);

    RawTowerDigitizer *digitizer = new RawTowerDigitizer("EmcRawTowerDigitizer");
    digitizer->Detector("CEMC");
    digitizer->set_digi_algorithm(RawTowerDigitizer::kNo_digitization);
    digitizer->Verbosity(0);
    se->registerSubsystem(digitizer);

    RawTowerCalibration *calib1
      = new RawTowerCalibration("RawTowerCalibrationCEMC");
    calib1->Detector("CEMC");
    // e.g. simple linear sampling fraction
    calib1->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
    calib1->set_calib_const_GeV_ADC(1.0/0.02); // example
    se->registerSubsystem(calib1);

    // Another tower-level calibration => TOWERINFO_CALIB2_CEMC
    CaloTowerCalib *calib2 = new CaloTowerCalib("CaloTowerCalibCEMC");
    calib2->set_detector_type(CaloTowerDefs::CEMC);
    calib2->set_inputNodePrefix("TOWERINFO_CALIB_");
    calib2->set_directURL("/sphenix/u/bseidlitz/work/analysis/EMCal_pi0_Calib_2023/macros/mc_calib_NoDigi.root");
    se->registerSubsystem(calib2);

//    // Finally, cluster building from TOWERINFO_CALIB2_
//    RawClusterBuilderTemplate *clusBuilder
//      = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
//    clusBuilder->Detector("CEMC");
//    clusBuilder->set_threshold_energy(0.07);
//    {
//      char* calroot = getenv("CALIBRATIONROOT");
//      if(!calroot)
//      {
//          std::cerr << "[WARNING] CALIBRATIONROOT not set. EmcProfile path may fail." << std::endl;
//      }
//      std::string emc_prof = std::string(calroot ? calroot : "") + "/EmcProfile/CEMCprof_Thresh30MeV.root";
//      try
//      {
//          clusBuilder->LoadProfile(emc_prof);
//      }
//      catch(const std::exception &ex)
//      {
//          std::cerr << "[EXCEPTION] while loading cluster profile: " << ex.what() << std::endl;
//      }
//    }
//    clusBuilder->set_UseTowerInfo(1);
//
//    se->registerSubsystem(clusBuilder);
   }

  ////////////////////////////////////////////////////////////
  // 5) handle additional calibration
  ////////////////////////////////////////////////////////////
  Process_Calo_Calib();

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
  PositionDependentCorrection *pdc
    = new PositionDependentCorrection("PositionDepCorr", finalOut);
  pdc->setIsSimulation(isSimulation);
  pdc->Verbosity(verbosity);
  se->registerSubsystem(pdc);

  ////////////////////////////////////////////////////////////
  // 7) Optionally set up an output manager if needed
  //    e.g. Fun4AllDstOutputManager *out = ...
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // 8) Print the Fun4All server structure
  ////////////////////////////////////////////////////////////
//  se->Print("ALL");

  std::cout << "[INFO] Running Fun4All with nevents = " << nevents << std::endl;
  se->run(nevents);

  ////////////////////////////////////////////////////////////
  // 9) End and cleanup
  ////////////////////////////////////////////////////////////
  se->End();
//  se->PrintTimer();
  delete se;

  std::cout << "[INFO] Done. Exiting macro." << std::endl;
  gSystem->Exit(0);
}

#endif // if ROOT_VERSION >= 6.00.0
