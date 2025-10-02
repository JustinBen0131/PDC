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
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fstream>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <phool/recoConsts.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/RawClusterDeadHotMask.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/CaloGeomMapping.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/RawClusterBuilderTemplate.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/RawTowerCalibration.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/CaloTowerCalib.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/CaloTowerStatus.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/BEmcRecCEMC.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_BEMC_clusterizer/CaloTowerBuilder.h"

#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawTowerGeom.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawTowerGeomContainer_Cylinderv1.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/TowerInfoDefs.h"
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src_CaloBASE/RawClusterv2.h"

#include <g4detectors/PHG4FullProjSpacalCellReco.h>
#include <g4calo/RawTowerBuilder.h>
#include <g4calo/RawTowerDigitizer.h>
#include <calowaveformsim/CaloWaveformSim.h>
#include <phparameter/PHParameterUtils.h>
#include <globalvertex/GlobalVertexReco.h>
#include <globalvertex/GlobalVertex.h>
#include <ffamodules/CDBInterface.h>
#include <ffaobjects/EventHeader.h>
#include <Event/Event.h>
#include <mbd/MbdReco.h>
#include <GlobalVariables.C>
#include "/sphenix/u/patsfan753/scratch/PDCrun24pp/src/PositionDependentCorrection.h"

#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic pop
#endif

R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libcalo_io.so)    // <— CaloBase (your local)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libcalo_reco.so)  // <— CaloReco (your local)
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
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libPDC.so)



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
    // --- decide runyear and mode up front ---
    // default assumptions
    int  runNo        = 24;     // run-24 data default
    bool isSimulation = false;

    // 0) explicit overrides from the environment
    if (const char* pm = gSystem->Getenv("PDC_MODE")) {
      std::string m(pm); for (auto& c: m) c = std::tolower(c);
      if (m=="sim" || m=="simulation") isSimulation = true;
      if (m=="data")                    isSimulation = false;
    }
    if (const char* er = gSystem->Getenv("PDC_RUNNUMBER")) {
      if (*er) runNo = std::atoi(er);
    }

    // 1) fallbacks if nothing was exported
    if (!gSystem->Getenv("PDC_MODE")) {
      // Heuristic: readable G4Hits list ⇒ sim
      isSimulation =
          !g4HitsListFile.empty() &&
          g4HitsListFile != "g4hits.list" &&
          gSystem->AccessPathName(g4HitsListFile.c_str(), kReadPermission) == kFALSE;
    }
    if (!gSystem->Getenv("PDC_RUNNUMBER")) {
      // Light sniff for minbias ⇒ run21
      auto toLower = [](std::string s){ for (auto& c : s) c = std::tolower(c); return s; };
      const std::string a = toLower(caloListFile);
      const std::string b = toLower(g4HitsListFile);
      const bool looksMB =
          a.find("minbias") != std::string::npos || a.find("_mb") != std::string::npos ||
          b.find("minbias") != std::string::npos || b.find("_mb") != std::string::npos ||
          a.find("run21")   != std::string::npos || b.find("run21")   != std::string::npos;
      if (looksMB) runNo = 21;
    }

    // Conditions DB flags (standard pattern in sPHENIX)
    rc->set_StringFlag("CDB_GLOBALTAG","MDC2");
    rc->set_uint64Flag("TIMESTAMP",    runNo);
    rc->set_IntFlag   ("RUNNUMBER",    runNo);
    std::cout << "[CFG] Using RUNNUMBER/TIMESTAMP = " << runNo
              << "  | isSimulation=" << (isSimulation?"true":"false") << std::endl;

    // Optional: touch one known key to prove we can reach the CDB (debug)
    CDBInterface::instance()->getUrl("EMCTOWERCALIB");


    if (verbosity > 0)
    {
      std::cout << "[DEBUG] Registering input managers for calo='"
                << caloListFile << "'  and hits='" << g4HitsListFile << "'\n"
                << "[CFG] isSimulation = " << (isSimulation ? "true" : "false") << std::endl;
    }

    auto* in0 = new Fun4AllDstInputManager("in0");
    bool looksList = (caloListFile.size() > 5 && caloListFile.substr(caloListFile.size()-5) == ".list");
    if (!looksList) {
      std::ifstream fin(caloListFile);
      std::string first;
      if (fin && std::getline(fin, first)) {
        // if the first line looks like a ROOT file path, treat it as a list
        looksList = (first.find(".root") != std::string::npos);
      }
    }
    if (looksList) in0->AddListFile(caloListFile);
    else           in0->AddFile(caloListFile);
    se->registerInputManager(in0);

    // (b) for G4Hits (only in simulation)
    if (isSimulation)
    {
      auto* in1 = new Fun4AllDstInputManager("in1");
      if (g4HitsListFile.size() > 5 && g4HitsListFile.substr(g4HitsListFile.size()-5) == ".list")
        in1->AddListFile(g4HitsListFile);
      else
        in1->AddFile(g4HitsListFile);
      se->registerInputManager(in1);
    }
    else
    {
      std::cout << "[INFO] Data run detected: not registering a G4Hits input manager." << std::endl;
    }

  GlobalVertexReco *gvertex = new GlobalVertexReco("GlobalVertexReco");
  gvertex->Verbosity(0);
  se->registerSubsystem(gvertex);
    
    BEmcRecCEMC* bemcPtr = new BEmcRecCEMC();
    bemcPtr->SetCylindricalGeometry();        // CEMC is a cylinder
    bemcPtr->set_UseDetailedGeometry(true);   // match builder ripple/geometry
    bemcPtr->SetTowerThreshold(0.030f);       // match builder Momenta threshold

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
    // isSimulation was decided above; reuse it here.
    string finalOut = outFileName.empty()
                      ? (isSimulation ? "output_PositionDep_sim.root"
                                      : "output_PositionDep_data.root")
                      : outFileName;
    
    
  if (useCustomClusterizer)
    {
      std::cout << "Building clusters (using custom RawClusterBuilderTemplate)" << std::endl;

        RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
        ClusterBuilder->Detector("CEMC");

        // Make builder geometry match PDC bemc: cylinder + detailed
        ClusterBuilder->SetCylindricalGeometry();
        ClusterBuilder->set_UseDetailedGeometry(true);
        ClusterBuilder->WriteClusterV2(true);

        // Use the same tower threshold you want PDC to see (optional but keeps Momenta identical)
        ClusterBuilder->set_threshold_energy(0.030f);

        std::string emc_prof = getenv("CALIBRATIONROOT");
        emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
        ClusterBuilder->LoadProfile(emc_prof);

        // TowerInfo mode → node name becomes CLUSTERINFO_CEMC, v2 → CLUSTERINFO_CEMC_V2
        ClusterBuilder->set_UseTowerInfo(1);

        ClusterBuilder->Verbosity(1);
        se->registerSubsystem(ClusterBuilder);

    }
    else
    {
      std::cout << "Skipping custom cluster builder (useCustomClusterizer = false)" << std::endl;
  }

  auto* pdc = new PositionDependentCorrection("PositionDepCorr", finalOut);
  /*------------------------------------------------------------
      1)  Mandatory hooks and global switches
  ------------------------------------------------------------*/
  pdc->setBEmcRec(bemcPtr);          // lead-tower finder from the clusteriser
  pdc->setIsSimulation(isSimulation);
  pdc->setFirstPassBvaluesOnly(false);
  pdc->UseSurveyGeometry(true);      // load barrel-tilt from CDB (recommended)
  pdc->UseSignedVz(false);
    // Automatically configure single-photon vs MinBias:
    //  - runNo is resolved above from $PDC_RUNNUMBER or by filename sniffing
    //  - default = true (single-photon), MinBias (runNo==21) => false
    const bool pdc_isMinBias = (runNo == 21);
    pdc->setIsSinglePhoton(!pdc_isMinBias);
  /*------------------------------------------------------------
      2)  π0-mass-window support
    ------------------------------------------------------------*/
  // Keep *false* for the first (pass-1) job that *produces* the μ/σ table.
  // Set *true* for the second (pass-2) job that *consumes* the table and
  // applies the slice-dependent mass cut.
  pdc->setMassFitsDone(false);
    

    /*------------------------------------------------------------
        3)  Anchor/partner quality gates (configurable)
        ------------------------------------------------------------
        AnchorQualityMode controls anchor-cluster screening in finalClusterLoop():
          - kNone : no anchor quality gate (accept all anchors).
          - kChi2 : require clusChi2 <= m_anchorChi2Cut (χ² gate).
          - kProb : require clusProb >= m_anchorProbMin  (probability gate).
          - kBoth : require BOTH χ² and probability gates to pass.

        Partner (second-cluster) probability gate in processClusterPairs():
          - enableSecondClusterProbCut(false)  → no partner probability cut (default)
          - enableSecondClusterProbCut(true)   → require c2Prob >= m_pairProbMin.
    ------------------------------------------------------------*/
//    // Default job settings (safe & explicit):
    pdc->setAnchorQualityMode(PositionDependentCorrection::AnchorQualityMode::kChi2);
    pdc->setAnchorChi2Cut(4.0f);
    pdc->enableSecondClusterProbCut(false);
    
    // Partner (cluster‑2) χ² gate – enable and set threshold
    pdc->enableSecondClusterChi2Cut(true);
    pdc->setSecondClusterChi2Max(4.0f);

//    // To switch to probability gates at 0.50, uncomment:
//     pdc->setAnchorQualityMode(PositionDependentCorrection::AnchorQualityMode::kProb);
//     pdc->setAnchorProbMin(0.50f);
//     pdc->enableSecondClusterProbCut(true);
//     pdc->setSecondClusterProbMin(0.50f);
    
  /*------------------------------------------------------------
      4)  Energy-binning mode
  ------------------------------------------------------------*/
  pdc->setBinningMode(PositionDependentCorrection::EBinningMode::kRange);

  pdc->Verbosity(2);                 // 0 = silent → raise for debugging

  // Finally register the module with Fun4All
  se->registerSubsystem(pdc);

  std::cout << "[INFO] Running Fun4All with nevents = " << nevents << std::endl;
  se->run(nevents);

  se->End();


  std::cout << "[INFO] Done. Exiting macro." << std::endl;
  return;
}
