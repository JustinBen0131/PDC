#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

// Standard includes
#include <caloreco/RawClusterPositionCorrection.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <g4detectors/PHG4FullProjSpacalCellReco.h>
#include <g4calo/RawTowerDigitizer.h>
#include <caloreco/RawTowerCalibration.h>
#include <calowaveformsim/CaloWaveformSim.h>
#include <phparameter/PHParameterUtils.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/SubsysReco.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerStatus.h>
#include <phool/PHRandomSeed.h>
#include <phool/PHCompositeNode.h>

// Your custom input manager that ignores run number mismatch:
#include "IgnoreRunNoNoSyncManager.h"  // <--- YOU NEED THIS HEADER

#include <G4_CEmc_Spacal.C>
#include <phool/recoConsts.h>
#include <ffamodules/CDBInterface.h>
#include <GlobalVariables.C>
#include <litecaloeval/LiteCaloEval.h>
#include <g4calo/RawTowerBuilder.h>

// Now we include our new advanced G4HitChecker subsystem:
#include "G4HitChecker.h"

// If you have your custom code:
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src/PositionDependentCorrection.h"

// #include <calotrigger/TriggerRunInfoReco.h> // omit if you don't need real trigger prescales
#include <ffaobjects/EventHeader.h>
#include <calib_emc_pi0/pi0EtaByEta.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>

// Load libraries
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
//R__LOAD_LIBRARY(libcalotrigger.so) // if you need trigger

#endif

// -----------------------------------------------------------------------
// Optional SubsysReco that overwrites run # in EventHeader
// -----------------------------------------------------------------------
#include <Event/Event.h>
#include <phool/getClass.h>

class ForceRunNumberReco : public SubsysReco
{
public:
  ForceRunNumberReco(const int newRunNo)
    : SubsysReco("ForceRunNumberReco"), m_newRunNo(newRunNo) {}

  int process_event(PHCompositeNode* topNode) override
  {
    EventHeader* evtheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
    if(evtheader)
    {
      // Overwrite the run number to e.g. 17
      evtheader->set_RunNumber(m_newRunNo);
    }
    return 0; // EVENT_OK
  }
private:
  int m_newRunNo;
};

// -----------------------------------------------------------------------

void Fun4All_PDC(int nevents = 0,
                 const std::string &fname = "condor/listFiles/gamma_dst_calo_cluster.list",
                 const std::string &fnamehits = "condor/listFiles/gamma_g4hits.list")
{
  std::cout << "[DEBUG] Entering Fun4All_PDC with nevents = " << nevents
            << ", fname = " << fname
            << ", fnamehits = " << fnamehits << std::endl;

  // 1) Setup Fun4AllServer
  Fun4AllServer *se = Fun4AllServer::instance();
  if (!se)
  {
    std::cerr << "[ERROR] Fun4AllServer::instance() returned nullptr. Exiting." << std::endl;
    gSystem->Exit(1);
  }

  // Increase verbosity for debugging
  int verbosity = 100;
  se->Verbosity(verbosity);

  // 2) recoConsts
  recoConsts *rc = recoConsts::instance();
  if (!rc)
  {
    std::cerr << "[ERROR] recoConsts::instance() returned nullptr. Exiting." << std::endl;
    gSystem->Exit(1);
  }

  // 3) Read first line from your DST file list
  std::ifstream file(fname);
  if (!file.is_open())
  {
    std::cerr << "[ERROR] Unable to open file: " << fname << std::endl;
    gSystem->Exit(1);
  }
  std::string first_file;
  if (!std::getline(file, first_file))
  {
    std::cerr << "[ERROR] Could not read first line from file: " << fname << std::endl;
    file.close();
    gSystem->Exit(1);
  }
  file.close();
  std::cout << "[DEBUG] First file in list is: " << first_file << std::endl;

  // 4) Extract run number from the first file
  auto runseg = Fun4AllUtils::GetRunSegment(first_file);
  int runnumber = runseg.first;
  std::cout << "[INFO] run number = " << runnumber << std::endl;

  // 5) Set up conditions DB flags
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  rc->set_IntFlag("RUNNUMBER", runnumber);

  try
  {
    CDBInterface::instance()->Verbosity(1);
    CDBInterface::instance()->getUrl("EMCTOWERCALIB");
  }
  catch (const std::exception &e)
  {
    std::cerr << "[EXCEPTION] CDBInterface usage: " << e.what() << std::endl;
  }

  // 6) Register input manager for DST
  Fun4AllDstInputManager *in = new Fun4AllDstInputManager("DST_TOWERS");
  if(!in)
  {
    std::cerr << "[ERROR] Creating Fun4AllDstInputManager (DST_TOWERS) failed." << std::endl;
    gSystem->Exit(1);
  }
  in->AddListFile(fname, 0);
  se->registerInputManager(in);

  // 7) Check hits file for .root lines
  bool hasHits = false;
  {
    std::ifstream hitsCheck(fnamehits);
    if(hitsCheck.is_open())
    {
      std::string testLine;
      while(std::getline(hitsCheck, testLine))
      {
        if(!testLine.empty() && testLine.find(".root") != std::string::npos)
        {
          hasHits = true;
          break;
        }
      }
      hitsCheck.close();
    }
  }

  // 8) If we do have hits, create an IgnoreRunNoNoSyncManager so we skip mismatch errors
  if(hasHits)
  {
    std::cout << "[DEBUG] Creating IgnoreRunNoNoSyncManager for DST_TOWERS2" << std::endl;
    IgnoreRunNoNoSyncManager *in2 = new IgnoreRunNoNoSyncManager("DST_TOWERS2");
    if(!in2)
    {
      std::cerr << "[ERROR] Failed to allocate IgnoreRunNoNoSyncManager (DST_TOWERS2)." << std::endl;
      gSystem->Exit(1);
    }

    // Optionally force the run number
    ForceRunNumberReco *fixRun = new ForceRunNumberReco(runnumber);
    in2->registerSubsystem(fixRun);

    in2->AddListFile(fnamehits, 0);
    se->registerInputManager(in2);
  }
  else
  {
    std::cout << "[INFO] No valid .root lines in " << fnamehits
              << " => skipping second input manager." << std::endl;
  }

  // 9) Build output histogram file name
  std::string filename = first_file.substr(first_file.find_last_of("/\\") + 1);
  std::string OutFile = Form("OUTHIST_iter_%s", filename.c_str());
  std::cout << "[DEBUG] Output histogram file will be: " << OutFile << std::endl;

  // Decide whether to do "old style" from G4 hits or waveforms
  bool recoFromHits_old = true;
  bool waveformRec = false;
  bool isSim = hasHits; // if we have G4 hits

  // ----------------------------------------------------------------------
  // >>> HERE IS WHERE WE ADD THE G4HitChecker SUBSYSTEM BEFORE ANY RECO <<<
  // ----------------------------------------------------------------------
  {
    // Suppose we allow 50 consecutive zero-hit events before exiting:
    G4HitChecker *hitcheck = new G4HitChecker("CEMC", /* maxZeroAllowed= */ 50);
    hitcheck->Verbosity(5); // can be adjusted from 0..10
    se->registerSubsystem(hitcheck);
  }

  // If we do want old style G4Hits -> TOWER pipeline
  if(recoFromHits_old && isSim)
  {
    // a) CEMC cell reco
    PHG4FullProjSpacalCellReco* cemc_cells = new PHG4FullProjSpacalCellReco("CEMCCYLCELLRECO");
    cemc_cells->Detector("CEMC");
    cemc_cells->Verbosity(verbosity);
    try
    {
      char* calibroot = getenv("CALIBRATIONROOT");
      if(!calibroot)
      {
        std::cerr << "[WARNING] CALIBRATIONROOT not set. Possibly can't load model file." << std::endl;
      }
      std::string model_file = std::string(calibroot ? calibroot : "") + "/CEMC/LightCollection/Prototype3Module.xml";
      cemc_cells->get_light_collection_model().load_data_file(
        model_file, "data_grid_light_guide_efficiency", "data_grid_fiber_trans");
    }
    catch(const std::exception &ex)
    {
      std::cerr << "[EXCEPTION] loading model file: " << ex.what() << std::endl;
    }
    se->registerSubsystem(cemc_cells);

    // b) RawTowerBuilder
    RawTowerBuilder *TowerBuilder = new RawTowerBuilder("EmcRawTowerBuilder");
    TowerBuilder->Detector("CEMC");
    TowerBuilder->set_sim_tower_node_prefix("SIM");
    TowerBuilder->Verbosity(0);
    se->registerSubsystem(TowerBuilder);

    // c) Digitizer
    const double sampling_fraction = 2e-02;
    RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kNo_digitization;
    RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("EmcRawTowerDigitizer");
    TowerDigitizer->Detector("CEMC");
    TowerDigitizer->Verbosity(verbosity);
    TowerDigitizer->set_digi_algorithm(TowerDigi);
    se->registerSubsystem(TowerDigitizer);

    // d) TowerCalibration
    RawTowerCalibration *TowerCalibration = new RawTowerCalibration("EmcRawTowerCalibration");
    TowerCalibration->Detector("CEMC");
    if (TowerDigi == RawTowerDigitizer::kNo_digitization)
    {
      TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
      TowerCalibration->set_calib_const_GeV_ADC(1.0 / sampling_fraction);
    }
    TowerCalibration->set_pedstal_ADC(0);
    se->registerSubsystem(TowerCalibration);

    // e) Another CaloTowerCalib for slopes etc.
    std::string calib_fname = "/sphenix/u/username/path/mc_calib_NoDigi.root";
    CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
    calibEMC->set_detector_type(CaloTowerDefs::CEMC);
    calibEMC->set_directURL(calib_fname.c_str());
    calibEMC->set_inputNodePrefix("TOWERINFO_RAW_");
    calibEMC->set_outputNodePrefix("TOWERINFO_CALIB_");
    se->registerSubsystem(calibEMC);

    // f) RawClusterBuilderTemplate
    RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
    ClusterBuilder->Detector("CEMC");
    ClusterBuilder->set_threshold_energy(0.07);
    {
      // Optionally load shape profile from e.g. $CALIBRATIONROOT/EmcProfile/CEMCprof_Thresh30MeV.root
      char* calroot = getenv("CALIBRATIONROOT");
      if(!calroot)
      {
        std::cerr << "[WARNING] CALIBRATIONROOT not set. EmcProfile path may fail." << std::endl;
      }
      std::string emc_prof = std::string(calroot ? calroot : "") + "/EmcProfile/CEMCprof_Thresh30MeV.root";
      try
      {
        ClusterBuilder->LoadProfile(emc_prof);
      }
      catch(const std::exception &ex)
      {
        std::cerr << "[EXCEPTION] while loading cluster profile: " << ex.what() << std::endl;
      }
    }
    ClusterBuilder->set_UseTowerInfo(1);
    se->registerSubsystem(ClusterBuilder);
  }

  // Possibly waveformsim pipeline
  if(waveformRec)
  {
    // Omitted
  }

  // 10) If you had real triggers, you'd do:
  // TriggerRunInfoReco *triggerruninforeco = new TriggerRunInfoReco();
  // se->registerSubsystem(triggerruninforeco);

  // 11) Register your PositionDependentCorrection module
  PositionDependentCorrection *pdc = new PositionDependentCorrection("PositionDepCorr", OutFile);
  se->registerSubsystem(pdc);

  // 12) Optional debug print of the entire server structure
  se->Print("ALL");
  if (se->topNode()) se->topNode()->print();

  // 13) Run
  std::cout << "[DEBUG] Running Fun4All with nevents = " << nevents << std::endl;
  se->run(nevents);

  // 14) End
  se->End();

  // 15) Print timing
  se->PrintTimer();

  // 16) Cleanup
  delete se;

  std::cout << "[INFO] Done! Exiting." << std::endl;
  gSystem->Exit(0);
}
