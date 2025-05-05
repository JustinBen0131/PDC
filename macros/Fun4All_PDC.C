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

// If you have a custom manager for ignoring run # mismatches
#include "IgnoreRunNoNoSyncManager.h"

#include <G4_CEmc_Spacal.C>
#include <phool/recoConsts.h>
#include <ffamodules/CDBInterface.h>
#include <GlobalVariables.C>
#include <litecaloeval/LiteCaloEval.h>
#include <g4calo/RawTowerBuilder.h>

// Optionally, a G4HitChecker subsystem:
#include "G4HitChecker.h"

// If you have custom code:
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src/PositionDependentCorrection.h"

// #include <calotrigger/TriggerRunInfoReco.h> // if you need real trigger prescales
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

//////////////////////////////////////////////////////////////////////////
// SubsysReco that overwrites run # in EventHeader (optional)
//////////////////////////////////////////////////////////////////////////
#include <Event/Event.h>
#include <phool/getClass.h>

///////////////////////////////////////////////////////////////////
// 1) A simple derived class that "adds" set_towerinfo_node_prefix(...)
//    without changing the real RawClusterBuilderTemplate code.
///////////////////////////////////////////////////////////////////
class MyRawClusterBuilderTemplate : public RawClusterBuilderTemplate
{
public:
  // Inherit RawClusterBuilderTemplate's constructor(s)
  using RawClusterBuilderTemplate::RawClusterBuilderTemplate;

  // Provide a dummy version of set_towerinfo_node_prefix(...)
  void set_towerinfo_node_prefix(const std::string &prefix)
  {
    // Do nothing, or just log a message:
    std::cout
      << "[WARNING] 'set_towerinfo_node_prefix' not truly supported "
      << "by this older version. Ignoring prefix='" << prefix << "'\n";
  }
};

///////////////////////////////////////////////////////////////////
// ForceRunNumberReco: Overwrites run # in EventHeader if desired
///////////////////////////////////////////////////////////////////
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
      // Overwrite the run number to e.g. 17 (or any desired run)
      evtheader->set_RunNumber(m_newRunNo);
    }
    return 0; // EVENT_OK
  }
private:
  int m_newRunNo;
};

///////////////////////////////////////////////////////////////////
// A combined Fun4All macro that can do either old-style G4Hits->towers
// and/or waveformsim pipeline for CEMC.
//
//  - 'recoFromHits_old': if true, do the standard cell reco -> tower -> calibrate
//  - 'waveformRec': if true, do the waveform pipeline
//
// Node name collisions are avoided by carefully setting output prefixes.
//
// The position-dependent correction module is also included at the end.
//
// We now use MyRawClusterBuilderTemplate in place of RawClusterBuilderTemplate
// to preserve calls to set_towerinfo_node_prefix() even though the base class
// lacks that function in your environment.
///////////////////////////////////////////////////////////////////
void Fun4All_PDC(int nevents = 0,
                 const std::string &fname = "condor/listFiles/gamma_dst_calo_cluster.list",
                 const std::string &fnamehits = "condor/listFiles/gamma_g4hits.list")
{
  //---------------------------------------------------------------------
  // 0) Setup the server and read the first input file to get run number
  //---------------------------------------------------------------------
  std::cout << "[DEBUG] Entering Fun4All_PDC with nevents = " << nevents
            << ", fname = " << fname
            << ", fnamehits = " << fnamehits << std::endl;

  Fun4AllServer *se = Fun4AllServer::instance();
  if (!se)
  {
    std::cerr << "[ERROR] Fun4AllServer::instance() returned nullptr. Exiting." << std::endl;
    gSystem->Exit(1);
  }
  int verbosity = 100;
  se->Verbosity(verbosity);

  recoConsts *rc = recoConsts::instance();
  if (!rc)
  {
    std::cerr << "[ERROR] recoConsts::instance() returned nullptr. Exiting." << std::endl;
    gSystem->Exit(1);
  }

  // read the first line from your DST list => extract run number
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

  auto runseg = Fun4AllUtils::GetRunSegment(first_file);
  int runnumber = runseg.first;
  std::cout << "[INFO] run number = " << runnumber << std::endl;

  // set up conditions DB
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

  //---------------------------------------------------------------------
  // 1) Register a DST input manager for "DST_TOWERS"
  //---------------------------------------------------------------------
  Fun4AllDstInputManager *in = new Fun4AllDstInputManager("DST_TOWERS");
  if(!in)
  {
    std::cerr << "[ERROR] Creating Fun4AllDstInputManager (DST_TOWERS) failed." << std::endl;
    gSystem->Exit(1);
  }
  in->AddListFile(fname, 0);
  se->registerInputManager(in);

  //---------------------------------------------------------------------
  // 2) Optionally register a second input manager for "DST_TOWERS2" if we have hits
  //---------------------------------------------------------------------
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

  if(hasHits)
  {
    // We might skip mismatch errors:
    std::cout << "[DEBUG] Creating IgnoreRunNoNoSyncManager for DST_TOWERS2" << std::endl;
    IgnoreRunNoNoSyncManager *in2 = new IgnoreRunNoNoSyncManager("DST_TOWERS2");
    if(!in2)
    {
      std::cerr << "[ERROR] Failed to allocate IgnoreRunNoNoSyncManager (DST_TOWERS2)." << std::endl;
      gSystem->Exit(1);
    }
    // Optionally force run number
    ForceRunNumberReco *fixRun = new ForceRunNumberReco(runnumber);
    in2->registerSubsystem(fixRun);

    in2->AddListFile(fnamehits, 1);
    se->registerInputManager(in2);
  }
  else
  {
    std::cout << "[INFO] No valid .root lines in " << fnamehits
              << " => skipping second input manager." << std::endl;
  }

  //---------------------------------------------------------------------
  // 3) Build output histogram file name
  //---------------------------------------------------------------------
  std::string filename = first_file.substr(first_file.find_last_of("/\\") + 1);
  std::string OutFile = Form("OUTHIST_iter_%s", filename.c_str());
  std::cout << "[DEBUG] Output histogram file will be: " << OutFile << std::endl;

  //---------------------------------------------------------------------
  // 4) Decide pipeline flags
  //---------------------------------------------------------------------
  bool recoFromHits_old = true;   // old style from G4 hits
  bool waveformRec      = false;   // do waveformsim pipeline
  bool isSim            = hasHits; // if we have G4 hits

  //---------------------------------------------------------------------
  // 5) Optionally add G4HitChecker subsystem if desired
  //---------------------------------------------------------------------
  {
    // Suppose we allow 50 consecutive zero-hit events before exiting
    G4HitChecker *hitcheck = new G4HitChecker("CEMC", /*maxZeroAllowed=*/50);
    hitcheck->Verbosity(5);
    se->registerSubsystem(hitcheck);
  }

  //---------------------------------------------------------------------
  // 6) "Old style" G4Hits -> tower pipeline
  //---------------------------------------------------------------------
  if(recoFromHits_old && isSim)
  {
    // (a) CEMC cell reco
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

    // (b) RawTowerBuilder
    RawTowerBuilder *TowerBuilder = new RawTowerBuilder("EmcRawTowerBuilder");
    TowerBuilder->Detector("CEMC");
    TowerBuilder->set_sim_tower_node_prefix("SIM");
    TowerBuilder->Verbosity(0);
    se->registerSubsystem(TowerBuilder);

    // (c) Digitizer
    const double sampling_fraction = 2.0e-02;
    RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kNo_digitization;
    RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("EmcRawTowerDigitizer");
    TowerDigitizer->Detector("CEMC");
    TowerDigitizer->Verbosity(verbosity);
    TowerDigitizer->set_digi_algorithm(TowerDigi);
    se->registerSubsystem(TowerDigitizer);

    // (d) TowerCalibration => writes TOWERINFO_CALIB_CEMC by default
    RawTowerCalibration *TowerCalibration = new RawTowerCalibration("EmcRawTowerCalibration");
    TowerCalibration->Detector("CEMC");
    if (TowerDigi == RawTowerDigitizer::kNo_digitization)
    {
      TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
      TowerCalibration->set_calib_const_GeV_ADC(1.0 / sampling_fraction);
    }
    TowerCalibration->set_pedstal_ADC(0);
    se->registerSubsystem(TowerCalibration);

    // (e) Another CaloTowerCalib => rename its output to avoid collisions
    // We rename to produce TOWERINFO_CALIB2_CEMC
    std::string calib_fname = "/sphenix/u/bseidlitz/work/analysis/EMCal_pi0_Calib_2023/macros/mc_calib_NoDigi.root";
    CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
    calibEMC->set_detector_type(CaloTowerDefs::CEMC);
    calibEMC->set_inputNodePrefix("TOWERINFO_CALIB_");
    calibEMC->set_outputNodePrefix("TOWERINFO_CALIB2_");
    calibEMC->set_directURL(calib_fname.c_str());
    se->registerSubsystem(calibEMC);

    // (f) RawCluster building from TOWERINFO_CALIB2_, but using our derived class
    MyRawClusterBuilderTemplate *ClusterBuilder
      = new MyRawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");

    ClusterBuilder->Detector("CEMC");
    ClusterBuilder->set_threshold_energy(0.07);
    {
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

    // If you want to read from TOWERINFO_CALIB2_:
    ClusterBuilder->set_towerinfo_node_prefix("TOWERINFO_CALIB2_");

    se->registerSubsystem(ClusterBuilder);
  }

  //---------------------------------------------------------------------
  // 7) Waveform pipeline
  //---------------------------------------------------------------------
  if(waveformRec)
  {
    // Example "CaloWaveformSim" for CEMC
    CaloWaveformSim* wfSim = new CaloWaveformSim("WaveSimCEMC");
    wfSim->set_detector_type(CaloTowerDefs::CEMC);
    wfSim->set_detector("CEMC");
    wfSim->set_nsamples(12);
    wfSim->set_pedestalsamples(12);
    wfSim->set_timewidth(0.2);
    wfSim->set_peakpos(6);

    // If needed, you can specify a slope calibration name
    // e.g. "cemc_pi0_twrSlope_v1_default"
    wfSim->set_calibName("cemc_pi0_twrSlope_v1_default");

    // If you need further modeling:
    // wfSim->get_light_collection_model().load_data_file(...);

    se->registerSubsystem(wfSim);

    // Next a "CaloTowerBuilder" that uses the waveforms to produce tower info
    CaloTowerBuilder *waveBuilder = new CaloTowerBuilder("CaloWaveBuilderCEMC");
    waveBuilder->set_detector_type(CaloTowerDefs::CEMC);
    waveBuilder->set_nsamples(12);
    waveBuilder->set_dataflag(false);
    waveBuilder->set_processing_type(CaloWaveformProcessing::TEMPLATE);
    waveBuilder->set_builder_type(CaloTowerDefs::kWaveformTowerv2);
    waveBuilder->set_softwarezerosuppression(true, 60); // example threshold
    se->registerSubsystem(waveBuilder);

    // We can also add "CaloTowerStatus" to mask certain bad towers if needed
    CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
    statusEMC->set_detector_type(CaloTowerDefs::CEMC);
    se->registerSubsystem(statusEMC);

    // Calibrate the waveform-based towers => produce TOWERINFO_WF_CALIB_CEMC or so
    CaloTowerCalib *wfCalib = new CaloTowerCalib("WaveCalibCEMC");
    wfCalib->set_detector_type(CaloTowerDefs::CEMC);
    wfCalib->set_inputNodePrefix("TOWERINFO_WF_"); // or "TOWERINFO_" if you set above
    wfCalib->set_outputNodePrefix("TOWERINFO_WF_CALIB_");
    // Possibly set a direct URL if you have a slope file
    // wfCalib->set_directURL("/path/to/waveformSlope.root");
    se->registerSubsystem(wfCalib);

    // Optionally do a second calibration pass for waveforms
    // or apply some MC re-calib:
    CaloTowerCalib *wfCalib2 = new CaloTowerCalib("WaveCalibCEMC2");
    wfCalib2->set_detector_type(CaloTowerDefs::CEMC);
    // example get from conditions DB:
    std::string MC_Calib = CDBInterface::instance()->getUrl("CEMC_MC_RECALIB");
    wfCalib2->set_directURL(MC_Calib.c_str());
    wfCalib2->set_inputNodePrefix("TOWERINFO_WF_CALIB_");
    wfCalib2->set_outputNodePrefix("TOWERINFO_WF_CALIB2_");
    // if you only want to do calibration, not re-building waveforms
    wfCalib2->set_doCalibOnly(true);
    se->registerSubsystem(wfCalib2);

    // Then build clusters from waveforms
    // but again, we use our derived MyRawClusterBuilderTemplate
    MyRawClusterBuilderTemplate *ClusterBuilderWF
      = new MyRawClusterBuilderTemplate("WaveRawClusterBuilder");

    ClusterBuilderWF->Detector("CEMC");
    ClusterBuilderWF->set_threshold_energy(0.07);
    {
      char* calroot = getenv("CALIBRATIONROOT");
      if(!calroot)
      {
        std::cerr << "[WARNING] waveformsim: CALIBRATIONROOT not set. " << std::endl;
      }
      std::string emc_prof = std::string(calroot ? calroot : "")
                           + "/EmcProfile/CEMCprof_Thresh30MeV.root";
      ClusterBuilderWF->LoadProfile(emc_prof);
    }
    ClusterBuilderWF->set_UseTowerInfo(1);

    // This will read from "TOWERINFO_WF_CALIB2_"
    ClusterBuilderWF->set_towerinfo_node_prefix("TOWERINFO_WF_CALIB2_");

    se->registerSubsystem(ClusterBuilderWF);
  }

  //---------------------------------------------------------------------
  // 8) Position‐Dependent Correction module
  //---------------------------------------------------------------------
  PositionDependentCorrection *pdc = new PositionDependentCorrection("PositionDepCorr", OutFile);
  pdc->setIsSimulation(isSim);
  se->registerSubsystem(pdc);

  //---------------------------------------------------------------------
  // 9) (Optional) Print entire server structure, then run
  //---------------------------------------------------------------------
  se->Print("ALL");
  if (se->topNode()) se->topNode()->print();

  std::cout << "[DEBUG] Running Fun4All with nevents = " << nevents << std::endl;
  se->run(nevents);

  //---------------------------------------------------------------------
  // 10) End job
  //---------------------------------------------------------------------
  se->End();
  se->PrintTimer();
  delete se;

  std::cout << "[INFO] Done! Exiting." << std::endl;
  gSystem->Exit(0);
}
