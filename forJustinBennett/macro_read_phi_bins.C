#include <phool/recoConsts.h>
#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libffamodules.so)

void macro_read_phi_bins()
{
  recoConsts *rc = recoConsts::instance();
  rc->set_uint64Flag("TIMESTAMP", 0); // Actually, you can put anything, it just needs to be defined
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2"); // idem 
  std::string inName=CDBInterface::instance()->getUrl("CALO_TOWER_GEOMETRY");
  // (std::string &) "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CALO_TOWER_GEOMETRY/dd/70/dd7077d8261e05e99941d5fb5f29fc7f_calo_geom_mapping.root"
  CDBTTree * cdbttree = new CDBTTree(inName);
  
  cdbttree->LoadCalibrations();
  for (int i = 0; i < 255; i++)
  {
    float phi_first = cdbttree->GetDoubleValue(i, "cemc_phi_first");
    float phi_second = cdbttree->GetDoubleValue(i, "cemc_phi_second");
    if (i == 0) {
      float shift = M_PI /	2 - phi_first;
      std::cout << "shift = " << shift << "[rad] = " << shift * 180. / M_PI << " [°]" << std::endl;
    }
    std::cout << "phi bin " << i << ": [" << phi_first << ", " << phi_second << "]" << std::endl;
  }
}
