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

  // Allow the conditions DB
  // (some older code might do 'Enable::CDB = true;' or similar)
//  CDBInterface::instance()->Verbosity(2);
//  CDBInterface::instance()->getUrl("EMCTOWERCALIB");
    
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


    bool doG4Reco = false; // presumably we want to convert from G4Hits if present
    if (doG4Reco)
    {
      // [DEBUG PRINT] Announce that we are entering the doG4Reco block
      std::cout << "[DEBUG] doG4Reco is TRUE => Setting up G4Hits to tower pipeline." << std::endl;

        PHG4FullProjSpacalCellReco *cemc_cells =
            new PHG4FullProjSpacalCellReco("CEMCCellReco");
        cemc_cells->Detector("CEMC");
        cemc_cells->Verbosity(0);
        se->registerSubsystem(cemc_cells);
//
//      // [DEBUG PRINT] Show the environment variable for CALIBRATIONROOT
//      const char* calrootEnv = std::getenv("CALIBRATIONROOT");
//      std::cout << "[DEBUG] CALIBRATIONROOT = "
//                << (calrootEnv ? calrootEnv : "(null)") << std::endl;
//
//      // [DEBUG PRINT] Print out the path we are about to load
//      std::string lightCollectionFile =
//          std::string(calrootEnv ? calrootEnv : "") + "/CEMC/LightCollection/Prototype3Module.xml";
//      std::cout << "[DEBUG] Loading light-collection data from: "
//                << lightCollectionFile << std::endl;
//
//      cemc_cells->get_light_collection_model().load_data_file(
//        lightCollectionFile,
//        "data_grid_light_guide_efficiency",
//        "data_grid_fiber_trans"
//      );
//
//      // [DEBUG PRINT] Confirm we are about to register this subsystem
//      std::cout << "[DEBUG] Registering 'cemc_cells' subsystem => PHG4FullProjSpacalCellReco" << std::endl;
////      se->registerSubsystem(cemc_cells);

      RawTowerBuilder *towerBuilder = new RawTowerBuilder("EmcRawTowerBuilder");
      towerBuilder->Detector("CEMC");
      towerBuilder->set_sim_tower_node_prefix("SIM");
      towerBuilder->Verbosity(0);

      // [DEBUG PRINT] Another note about the next subsystem
      std::cout << "[DEBUG] Registering 'towerBuilder' => RawTowerBuilder" << std::endl;
      se->registerSubsystem(towerBuilder);

//      const double sampling_fraction = 2e-02;
//      const double photoelectron_per_GeV = 500;
//      RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kNo_digitization;
//
//      // [DEBUG PRINT] Show which digitization mode was chosen
//      std::cout << "[DEBUG] TowerDigi algorithm = " << TowerDigi << " (kNo_digitization=0, kSimple_photon_digitization=1, etc.)" << std::endl;
//
//      RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("EmcRawTowerDigitizer");
//      TowerDigitizer->Detector("CEMC");
//      TowerDigitizer->Verbosity(0);
//      TowerDigitizer->set_digi_algorithm(TowerDigi);
//
//      if (TowerDigi == RawTowerDigitizer::kSimple_photon_digitization)
//      {
//        // [DEBUG PRINT] We are in the photon-digitization block
//        std::cout << "[DEBUG] TowerDigi == kSimple_photon_digitization => setting variable pedestals, zero-suppression, etc." << std::endl;
//
//        TowerDigitizer->set_variable_pedestal(true);
//        TowerDigitizer->set_pedstal_central_ADC(0);
//        TowerDigitizer->set_pedstal_width_ADC(0);
//        TowerDigitizer->set_photonelec_ADC(1);
//        TowerDigitizer->set_photonelec_yield_visible_GeV(photoelectron_per_GeV / sampling_fraction);
//        TowerDigitizer->set_variable_zero_suppression(true);
//        TowerDigitizer->set_zero_suppression_ADC(0);
//
//        // [DEBUG PRINT] Print out we are filling from "EMCTOWERCALIB"
//        std::cout << "[DEBUG] TowerDigitizer => PHParameterUtils::FillPHParametersFromCDB(..., \"EMCTOWERCALIB\")" << std::endl;
//        PHParameterUtils::FillPHParametersFromCDB(TowerDigitizer->GetParameters(), "EMCTOWERCALIB");
//        // TowerDigitizer->GetParameters().ReadFromCDB("EMCTOWERCALIB");
//      }
//
//      // [DEBUG PRINT] Next subsystem: TowerDigitizer
//      std::cout << "[DEBUG] Registering 'TowerDigitizer' => RawTowerDigitizer" << std::endl;
//      se->registerSubsystem(TowerDigitizer);

//      RawTowerCalibration *TowerCalibration = new RawTowerCalibration("EmcRawTowerCalibration");
//      TowerCalibration->Detector("CEMC");
//
//      // e.g. simple linear sampling fraction
//      if (TowerDigi == RawTowerDigitizer::kNo_digitization)
//      {
//        std::cout << "[DEBUG] TowerCalibration => set_calib_algorithm(kSimple_linear_calibration)" << std::endl;
//        TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
//        TowerCalibration->set_calib_const_GeV_ADC(1.0 / sampling_fraction);
//      }

//      TowerCalibration->set_pedstal_ADC(0);
//      if (TowerDigi == RawTowerDigitizer::kSimple_photon_digitization)
//      {
//        // [DEBUG PRINT] Show we are in the photon-digitization block for calibration
//        std::cout << "[DEBUG] TowerCalibration => set_calib_algorithm(kTower_by_tower_calibration) and read from EMCTOWERCALIB" << std::endl;
//
//        TowerCalibration->set_calib_algorithm(RawTowerCalibration::kTower_by_tower_calibration);
//        //TowerCalibration->GetCalibrationParameters().ReadFromCDB("EMCTOWERCALIB");
//        PHParameterUtils::FillPHParametersFromCDB(TowerCalibration->GetCalibrationParameters(), "EMCTOWERCALIB");
//        TowerCalibration->set_variable_GeV_ADC(true);
//        TowerCalibration->set_variable_pedestal(true);
//      }
//
//      // [DEBUG PRINT] Register the calibration subsystem
//      std::cout << "[DEBUG] Registering 'TowerCalibration' => RawTowerCalibration" << std::endl;
//      se->registerSubsystem(TowerCalibration);

//      std::string calib_fname = "/sphenix/u/bseidlitz/work/analysis/EMCal_pi0_Calib_2023/macros/mc_calib_NoDigi.root";
//      CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
//      calibEMC->set_detector_type(TowerInfoDefs::CEMC);
//      calibEMC->set_directURL(calib_fname.c_str());
//      calibEMC->set_inputNodePrefix("TOWERINFO_RAW_");
//      calibEMC->set_outputNodePrefix("TOWERINFO_CALIB_");

      // [DEBUG PRINT] We also register the final "CaloTowerCalib" => "CEMCCALIB"
//      std::cout << "[DEBUG] Registering 'calibEMC' => CaloTowerCalib with directURL=" << calib_fname << std::endl;
//      se->registerSubsystem(calibEMC);

      RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
      ClusterBuilder->Detector("CEMC");
      ClusterBuilder->set_threshold_energy(0.070);

      std::string emc_prof = getenv("CALIBRATIONROOT");
      emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
      // [DEBUG PRINT] Show the path for CEMC profile
      std::cout << "[DEBUG] In RawClusterBuilderTemplate => Loading profile: " << emc_prof << std::endl;
      ClusterBuilder->LoadProfile(emc_prof);

      ClusterBuilder->set_UseTowerInfo(1);
      ClusterBuilder->set_UseAltZVertex(3);

      // [DEBUG PRINT] Finally, register the cluster builder
      std::cout << "[DEBUG] Registering 'ClusterBuilder' => RawClusterBuilderTemplate" << std::endl;
      se->registerSubsystem(ClusterBuilder);

      // [DEBUG PRINT] Done setting up doG4Reco pipeline
      std::cout << "[DEBUG] doG4Reco block completed." << std::endl;
    }

//    bool waveformRec = false; // Or true if you want to enable it
//    if (waveformRec)  // user-provided block
//    {
//      TRandom3 randGen;
//      // get seed
//      unsigned int seed = PHRandomSeed();
//      randGen.SetSeed(seed);
//      // a int from 0 to 3259
//      int sequence = randGen.Integer(3260);
//      // pad the name
//      std::ostringstream opedname;
//      opedname << "/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/Pedestal/DSTendofpp/pedestal-54256-0"
//               << std::setw(4) << std::setfill('0') << sequence << ".root";
//
//      std::string pedestalname = opedname.str();
//
//      //pedestalname = "pedestal-00046796.root";
//      Fun4AllInputManager *hitsin = new Fun4AllNoSyncDstInputManager("DST2");
//      hitsin->AddFile(pedestalname);
//      hitsin->Repeat();
//      se->registerInputManager(hitsin);
//
//
//      CaloWaveformSim* caloWaveformSim = new CaloWaveformSim();
//      caloWaveformSim->set_detector_type(CaloTowerDefs::CEMC);
//      caloWaveformSim->set_detector("CEMC");
//      caloWaveformSim->set_nsamples(12);
//      caloWaveformSim->set_pedestalsamples(12);
//      caloWaveformSim->set_timewidth(0.2);
//      caloWaveformSim->set_peakpos(6);
//      caloWaveformSim->set_calibName("cemc_pi0_twrSlope_v1_default");
//      //caloWaveformSim->get_light_collection_model().load_data_file(
//      //string(getenv("CALIBRATIONROOT")) +
//      //string("/CEMC/LightCollection/Prototype3Module.xml"),
//      //"data_grid_light_guide_efficiency", "data_grid_fiber_trans");
//      se->registerSubsystem(caloWaveformSim);
//
//      CaloTowerBuilder * ca2 = new CaloTowerBuilder();
//      ca2->set_detector_type(CaloTowerDefs::CEMC);
//      ca2->set_nsamples(12);
//      ca2->set_dataflag(false);
//      ca2->set_processing_type(CaloWaveformProcessing::TEMPLATE);
//      ca2->set_builder_type(CaloTowerDefs::kWaveformTowerv2);
//      ca2->set_softwarezerosuppression(true, 60);
//      se->registerSubsystem(ca2);
//
//      CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
//      statusEMC->set_detector_type(CaloTowerDefs::CEMC);
//      se->registerSubsystem(statusEMC);
//
//      std::cout << "Calibrating EMCal" << std::endl;
//      CaloTowerCalib *calibEMCchunk = new CaloTowerCalib("CEMCCALIB");
//      calibEMCchunk->set_detector_type(CaloTowerDefs::CEMC);
//      calibEMCchunk->set_outputNodePrefix("TOWERINFO_CALIB_");
//      se->registerSubsystem(calibEMCchunk);
//
//      // string MC_Calib = "calib.root";
//      std::string MC_Calib = CDBInterface::instance()->getUrl("CEMC_MC_RECALIB");
//      std::cout << "MC_Calib " << MC_Calib.c_str() << std::endl;
//      CaloTowerCalib *calibEMC2 = new CaloTowerCalib("CEMCCALIB2");
//      calibEMC2->set_detector_type(CaloTowerDefs::CEMC);
//      calibEMC2->set_directURL(MC_Calib.c_str());
//      calibEMC2->set_inputNodePrefix("TOWERINFO_CALIB_");
//      calibEMC2->set_outputNodePrefix("TOWERINFO_CALIB_");
//      calibEMC2->set_doCalibOnly(true);
//      se->registerSubsystem(calibEMC2);
//
//      RawClusterBuilderTemplate *ClusterBuilderWave = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
//      ClusterBuilderWave->Detector("CEMC");
//      ClusterBuilderWave->set_threshold_energy(0.070);
//      std::string emc_prof2 = getenv("CALIBRATIONROOT");
//      emc_prof2 += "/EmcProfile/CEMCprof_Thresh30MeV.root";
//      ClusterBuilderWave->LoadProfile(emc_prof2);
//      ClusterBuilderWave->set_UseTowerInfo(1);
//      ClusterBuilderWave->set_UseAltZVertex(3);
//      se->registerSubsystem(ClusterBuilderWave);
//    }

  ////////////////////////////////////////////////////////////
  // 5) handle additional calibration
  ////////////////////////////////////////////////////////////
//  Process_Calo_Calib();
    
    
  BEmcRecCEMC *bemcPtr = new BEmcRecCEMC();
  bemcPtr->SetCylindricalGeometry();

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
  pdc->setBEmcRec(bemcPtr);
  pdc->setIsSimulation(isSimulation);
  pdc->UseSurveyGeometry(true);
  // optionally specify a different tag / timestamp
  // pdc->SetCDBTag("MDC3");
  // pdc->SetTimeStamp(runNumber);
    
  pdc->Verbosity(3);
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
  //delete se;


  std::cout << "[INFO] Done. Exiting macro." << std::endl;
  return;
}
