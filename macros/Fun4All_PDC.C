#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

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
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <G4_CEmc_Spacal.C>
#include <phool/recoConsts.h>
#include <ffamodules/CDBInterface.h>
#include <GlobalVariables.C>
#include <litecaloeval/LiteCaloEval.h>
#include <g4calo/RawTowerBuilder.h>
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src/PositionDependentCorrection.h"
#include <calotrigger/TriggerRunInfoReco.h>
#include <calib_emc_pi0/pi0EtaByEta.h>

R__LOAD_LIBRARY(libcdbobjects)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)
R__LOAD_LIBRARY(libcalibCaloEmc_pi0.so)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libPDC.so)
R__LOAD_LIBRARY(libcalotrigger.so)

#endif

void Fun4All_PDC(int nevents = 5e2, const std::string &fname = "condor/listFiles/gamma_dst_calo_cluster.list",const std::string &fnamehits = "condor/listFiles/gamma_g4hits.list")
{
  Fun4AllServer *se = Fun4AllServer::instance();
 // se->Verbosity(0);
  int verbosity = 0;
  se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();

  ifstream file(fname);
  string first_file;
  getline(file, first_file);

  //===============
  // conditions DB flags
  //===============
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(first_file);
  int runnumber = runseg.first;
  cout << "run number = " << runnumber << endl;

  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  CDBInterface::instance()->Verbosity(1);
  CDBInterface::instance()->getUrl("EMCTOWERCALIB");

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST_TOWERS");
  in->AddListFile(fname,1);
  se->registerInputManager(in);

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DST_TOWERS2");
  in2->AddListFile(fnamehits,1);
  se->registerInputManager(in2);

  std::string filename = first_file.substr(first_file.find_last_of("/\\") + 1);
  std::string OutFile = Form("OUTHIST_iter_%s",filename.c_str());


  ////////////////////
  // hit reco
  bool recoFromHits_old = true;
  bool waveformRec = false;


  if (recoFromHits_old){
    PHG4FullProjSpacalCellReco *cemc_cells = new PHG4FullProjSpacalCellReco("CEMCCYLCELLRECO");
    cemc_cells->Detector("CEMC");
    cemc_cells->Verbosity(verbosity);
    cemc_cells->get_light_collection_model().load_data_file(
    string(getenv("CALIBRATIONROOT")) + string("/CEMC/LightCollection/Prototype3Module.xml"), "data_grid_light_guide_efficiency", "data_grid_fiber_trans");
    se->registerSubsystem(cemc_cells);

    RawTowerBuilder *TowerBuilder = new RawTowerBuilder("EmcRawTowerBuilder");
    TowerBuilder->Detector("CEMC");
    TowerBuilder->set_sim_tower_node_prefix("SIM");
    TowerBuilder->Verbosity(verbosity);
    se->registerSubsystem(TowerBuilder);

    const double sampling_fraction = 2e-02;
    const double photoelectron_per_GeV = 500;
    //RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kSimple_photon_digitization;
    RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kNo_digitization;

    RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("EmcRawTowerDigitizer");
    TowerDigitizer->Detector("CEMC");
    TowerDigitizer->Verbosity(verbosity);
    TowerDigitizer->set_digi_algorithm(TowerDigi);
    if (TowerDigi == RawTowerDigitizer::kSimple_photon_digitization){
      TowerDigitizer->set_variable_pedestal(true);
      TowerDigitizer->set_pedstal_central_ADC(0);
      TowerDigitizer->set_pedstal_width_ADC(0);
      TowerDigitizer->set_photonelec_ADC(1);
      TowerDigitizer->set_photonelec_yield_visible_GeV(photoelectron_per_GeV / sampling_fraction);
      TowerDigitizer->set_variable_zero_suppression(true);
      TowerDigitizer->set_zero_suppression_ADC(0);
      PHParameterUtils::FillPHParametersFromCDB(TowerDigitizer->GetParameters(),"EMCTOWERCALIB");
      //TowerDigitizer->GetParameters().ReadFromCDB("EMCTOWERCALIB");
    }
    se->registerSubsystem(TowerDigitizer);


    RawTowerCalibration *TowerCalibration = new RawTowerCalibration("EmcRawTowerCalibration");
    TowerCalibration->Detector("CEMC");
    if (TowerDigi == RawTowerDigitizer::kNo_digitization){
      TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
      TowerCalibration->set_calib_const_GeV_ADC(1.0/sampling_fraction);
    }
    TowerCalibration->set_pedstal_ADC(0);
    if (TowerDigi == RawTowerDigitizer::kSimple_photon_digitization){
      TowerCalibration->set_calib_algorithm(RawTowerCalibration::kTower_by_tower_calibration);
      //TowerCalibration->GetCalibrationParameters().ReadFromCDB("EMCTOWERCALIB");
      PHParameterUtils::FillPHParametersFromCDB(TowerCalibration->GetCalibrationParameters(),"EMCTOWERCALIB");
      TowerCalibration->set_variable_GeV_ADC(true);
      TowerCalibration->set_variable_pedestal(true);
    }
    se->registerSubsystem(TowerCalibration);

    string calib_fname = "/sphenix/u/bseidlitz/work/analysis/EMCal_pi0_Calib_2023/macros/mc_calib_NoDigi.root";
    CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
    calibEMC->set_detector_type(CaloTowerDefs::CEMC);
    calibEMC->set_directURL(calib_fname.c_str());
    calibEMC->set_inputNodePrefix("TOWERINFO_RAW_");
    calibEMC->set_outputNodePrefix("TOWERINFO_CALIB_");
    se->registerSubsystem(calibEMC);

    RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
    ClusterBuilder->Detector("CEMC");
    ClusterBuilder->set_threshold_energy(0.070);
    std::string emc_prof = getenv("CALIBRATIONROOT");
    emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
    ClusterBuilder->LoadProfile(emc_prof);
    ClusterBuilder->set_UseTowerInfo(1);
    se->registerSubsystem(ClusterBuilder);
  }

  if (waveformRec){

    TRandom3 randGen;
    //get seed
    unsigned int seed = PHRandomSeed();
    randGen.SetSeed(seed);
    // a int from 0 to 3259
    int sequence = randGen.Integer(3260);
    // pad the name
    std::ostringstream opedname;
    opedname << "/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/Pedestal/DSTendofpp/pedestal-54256-0" << std::setw(4) << std::setfill('0') << sequence << ".root";

    std::string pedestalname = opedname.str();
    
    //pedestalname = "pedestal-00046796.root";
    Fun4AllInputManager *hitsin = new Fun4AllNoSyncDstInputManager("DST2");
    hitsin->AddFile(pedestalname);
    hitsin->Repeat();
    se->registerInputManager(hitsin);
    
    CEMC_Cells();
    CEMC_Towers();

    CaloWaveformSim* caloWaveformSim = new CaloWaveformSim();
    caloWaveformSim->set_detector_type(CaloTowerDefs::CEMC);
    caloWaveformSim->set_detector("CEMC");
    caloWaveformSim->set_nsamples(12);
    caloWaveformSim->set_pedestalsamples(12);
    caloWaveformSim->set_timewidth(0.2);
    caloWaveformSim->set_peakpos(6);
    caloWaveformSim->set_calibName("cemc_pi0_twrSlope_v1_default");
    //caloWaveformSim->get_light_collection_model().load_data_file(
    //string(getenv("CALIBRATIONROOT")) +
    //string("/CEMC/LightCollection/Prototype3Module.xml"),
    //"data_grid_light_guide_efficiency", "data_grid_fiber_trans");
    se->registerSubsystem(caloWaveformSim);

    CaloTowerBuilder * ca2 = new CaloTowerBuilder();
    ca2->set_detector_type(CaloTowerDefs::CEMC);
    ca2->set_nsamples(12);
    ca2->set_dataflag(false);
    ca2->set_processing_type(CaloWaveformProcessing::TEMPLATE);
    ca2->set_builder_type(CaloTowerDefs::kWaveformTowerv2);
    ca2->set_softwarezerosuppression(true, 60);
    se->registerSubsystem(ca2);

    CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
    statusEMC->set_detector_type(CaloTowerDefs::CEMC);
    se->registerSubsystem(statusEMC);

    std::cout << "Calibrating EMCal" << std::endl;
    CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
    calibEMC->set_detector_type(CaloTowerDefs::CEMC);
    calibEMC->set_outputNodePrefix("TOWERINFO_CALIB_");
    se->registerSubsystem(calibEMC);

   // string MC_Calib = "calib.root";
    std::string MC_Calib = CDBInterface::instance()->getUrl("CEMC_MC_RECALIB");
    cout << "MC_Calib " << MC_Calib.c_str() << endl;
    CaloTowerCalib *calibEMC2 = new CaloTowerCalib("CEMCCALIB2");
    calibEMC2->set_detector_type(CaloTowerDefs::CEMC);
    calibEMC2->set_directURL(MC_Calib.c_str());
    calibEMC2->set_inputNodePrefix("TOWERINFO_CALIB_");
    calibEMC2->set_outputNodePrefix("TOWERINFO_CALIB_");
    calibEMC2->set_doCalibOnly(true);
    se->registerSubsystem(calibEMC2);
  
    RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
    ClusterBuilder->Detector("CEMC");
    ClusterBuilder->set_threshold_energy(0.070);
    std::string emc_prof = getenv("CALIBRATIONROOT");
    emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
    ClusterBuilder->LoadProfile(emc_prof);
    ClusterBuilder->set_UseTowerInfo(1);
    se->registerSubsystem(ClusterBuilder);


  }

  ////////////////////////////////////////////////
  // analysis modules
  ////////////////////////////////////////////////

  TriggerRunInfoReco *triggerruninforeco = new TriggerRunInfoReco();
  se->registerSubsystem(triggerruninforeco);
    
  PositionDependentCorrection* pdc = new PositionDependentCorrection("PositionDepCorr", OutFile);
  se->registerSubsystem(pdc);

  se->run(nevents);
  se->End();
  se->PrintTimer();
  delete se;

  std::cout << "All done!" << std::endl;
  gSystem->Exit(0);

}
