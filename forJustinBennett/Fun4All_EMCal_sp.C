#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>

#include <phool/recoConsts.h>


R__LOAD_LIBRARY(libcdbobjects)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libcalo_io.so)
R__LOAD_LIBRARY(libphool.so);
R__LOAD_LIBRARY(libmbd_io.so)
R__LOAD_LIBRARY(libcentrality_io.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libphg4hit.so)
R__LOAD_LIBRARY(libSubsysReco.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4calo.so)

//#include <caloana/CaloAna_spi0.h>
//R__LOAD_LIBRARY(libcaloana.so)

void Fun4All_EMCal_sp(int nevents = 1, const std::string &fname = "inputdata.txt",const std::string &fnamehits = "inputdatahits.txt")
{

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(10000000);

  // se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();
  ifstream file(fname);
  string first_file;
  getline(file, first_file);
  std::cout << fname << std::endl;
  //===============
  // conditions DB flags
  //===============
  std::cout << first_file << std::endl;
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(first_file);
  int runnumber = runseg.first;
  cout << "run number = " << runnumber << endl;

  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST_TOWERS");
  in->Verbosity(0);
  in->AddListFile(fname);
  se->registerInputManager(in);

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DST_TOWERS2");
  in2->Verbosity(0);
  in2->AddListFile(fnamehits);
  se->registerInputManager(in2);

  se->Print("NODETREE");
  se->run(nevents);
  se->End();
  se->PrintTimer();
  delete se;

  TFile* f_done_signal = new TFile("DONE.root","recreate");
  std::cout << "All done!" << std::endl;
  gSystem->Exit(0);
}

