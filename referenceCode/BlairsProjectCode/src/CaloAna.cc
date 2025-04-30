#include "CaloAna.h"

// G4Hits includes
#include <TLorentzVector.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

// Tower includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// MBD
#include <mbd/BbcGeom.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtHit.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNtuple.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TTree.h>

#include <Event/Event.h>
#include <Event/packet.h>
#include <cassert>
#include <sstream>
#include <string>

#include <iostream>
#include <utility>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include"TRandom3.h"

#include <cdbobjects/CDBTTree.h>  // for CDBTTree
#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

#include <g4main/PHG4TruthInfoContainer.h>

R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)

using namespace std;

CaloAna::CaloAna(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , detector("HCALIN")
  , outfilename(filename)
{
  _eventcounter = 0;
}

CaloAna::~CaloAna()
{
  delete hm;
  delete g4hitntuple;
  delete g4cellntuple;
  delete towerntuple;
  delete clusterntuple;
}

int CaloAna::Init(PHCompositeNode*)
{
  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here

  outfile = new TFile(outfilename.c_str(), "RECREATE");

  // correlation plots
  h_mass_eta_lt = new TH2F("h_mass_eta_lt", "", 50, 0, 0.5,96,0,96);
  h_mass_eta_lt_rw = new TH2F("h_mass_eta_lt_rw", "", 50, 0, 0.5,96,0,96);
  h_pt_eta = new TH2F("h_pt_eta", "", 100, 0, 10,96,0,96);
  h_pt_eta_rw = new TH2F("h_pt_eta_rw", "", 100, 0, 10,96,0,96);

  h_cemc_etaphi = new TH2F("h_cemc_etaphi", "", 96, 0, 96, 256, 0, 256);

  // 1D distributions
  h_InvMass = new TH1F("h_InvMass", "Invariant Mass", 500, 0, 1.0);
  h_InvMass_w = new TH1F("h_InvMass_w", "Invariant Mass", 500, 0, 1.0);
  h_InvMassMix = new TH1F("h_InvMassMix", "Invariant Mass", 120, 0, 1.2);
  h_tower_e = new TH1F("h_tower_e","",1000,-1,5);

  // cluster QA
  h_etaphi_clus = new TH2F("h_etaphi_clus", "", 140, -1.2, 1.2, 64, -1 * TMath::Pi(), TMath::Pi());
  h_clusE = new TH1F("h_clusE", "", 100, 0, 10);
  h_clusE_nTow = new TH2F("h_clusE_nTow","",20,0,20,50,0,50);

  h_emcal_e_eta = new TH1F("h_emcal_e_eta", "", 96, 0, 96);

  h_pt1 = new TH1F("h_pt1", "", 100, 0, 5);
  h_pt2 = new TH1F("h_pt2", "", 100, 0, 5);

  h_nclusters = new TH1F("h_nclusters", "", 100, 0, 100);

  h_truth_eta = new TH1F("h_truth_eta", "", 100, -1.2, 1.2);
  h_truth_e = new TH1F("h_truth_e", "", 100, 0, 10);
  h_truth_pt = new TH1F("h_truth_pt", "", 100, 0, 10);
  h_delR_recTrth = new TH1F("h_delR_recTrth", "", 500, 0, 5);
  h_matched_res = new TH2F("h_matched_res","",100,0,1.5,20,-1,1);
  h_res_e = new TH2F("h_res_e","",100,0,1.5,20,0,20);
  h_res_e_phi = new TH3F("h_res_e_phi","",100,0,1.5,10,0,20,256,0,256);
  h_res_e_eta = new TH3F("h_res_e_eta","",300,0,1.5,40,0,20,96,0,96);
  h_m_pt_eta = new TH3F("h_m_pt_eta","",70,0,0.7,10,0,10,96,0,96);
  h_m_ptTr_eta= new TH3F("h_m_ptTr_eta","",70,0,0.7,10,0,10,96,0,96);
  h_m_ptTr_eta_trKin = new TH3F("h_m_ptTr_eta_trKin","",70,0,0.7,10,0,10,96,0,96);
  h_res = new TH1F("h_res", "", 50, 0, 1.5);
  h_delEta_e_eta = new TH3F("h_delEta_e_eta","",100,-0.1,0.1,10,0,20,96,0,96);
  h_delR_e_eta = new TH3F("h_delR_e_eta","",100,-0.1,0.1,10,0,20,96,0,96);
  h_delPhi_e_eta = new TH3F("h_delPhi_e_eta","",100,-0.3,0.3,20,0,20,96,0,96);
  h_delPhi_e_phi = new TH3F("h_delPhi_e_phi","",100,-0.1,0.1,20,0,20,256,0,256);
  pr_eta_shower = new TProfile("pr_eta_shower","",96,-48.5,47.5, -1,1.5);
  pr_phi_shower = new TProfile("pr_phi_shower","",256,-128.5,127.5, -1,1.5);
  h_vert_xy = new TH2F("h_vert_xy","",500,-120,120,500,-120,120);
  h_truthE = new TH1F("h_truthE","",10000,0,30);


  // POSITION DEPENDENT STUFF
  for (int ie=0; ie<NBinsBlock; ie++){
    for (int ip=0; ip<NBinsBlock; ip++){
      h_mass_block_pt[ie][ip] = new TH2F(Form("h_mass_block_%d_%d_pt",ie,ip),"",100,0,1,5,0,10);
      h_res_block_E[ie][ip] = new TH2F(Form("h_res_block_%d_%d_E",ie,ip),"",120,0,1.2,5,0,10);
    }
  }

  h3_cluster_block_cord_pt = new TH3F("h2_cluster_block_cord_pt","",14,-0.5,1.5,14,-0.5,1.5,5,0,10);
  h_block_bin = new TH1F("h_block_bin","",14,-0.5,1.5);


  //////////////////////////
  // pT rewieghting
  frw = new TFile("/sphenix/u/bseidlitz/work/analysis/EMCal_pi0_Calib_2023/macros/rw_pt.root");
  //for(int i=0; i<96; i++) h_pt_rw[i] = (TH1F*) frw->Get(Form("h_pt_eta%d", i));

  rnd = new TRandom3(); 


  return 0;
}

int CaloAna::process_event(PHCompositeNode* topNode)
{
  _eventcounter++;

  process_towers(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_towers(PHCompositeNode* topNode)
{
  if ((_eventcounter % 1000) == 0) std::cout << _eventcounter << std::endl;

  // float emcaldownscale = 1000000 / 800;

  float emcal_hit_threshold = 0.5;  // GeV
  if(debug) std::cout << "-----------------------------------" << std::endl;

  // cuts
  float maxAlpha = 0.6;
  float clus_chisq_cut = 10000;
  float nClus_ptCut = 0.5;
  int max_nClusCount = 3000000;

  //-----------------------get vertex----------------------------------------//

  float vtx_z = 0;
  if (getVtx)
  {
    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (!vertexmap)
    {
      // std::cout << PHWHERE << " Fatal Error - GlobalVertexMap node is missing"<< std::endl;
      std::cout << "CaloAna GlobalVertexMap node is missing" << std::endl;
      // return Fun4AllReturnCodes::ABORTRUN;
    }
    if (vertexmap && !vertexmap->empty())
    {
      GlobalVertex* vtx = vertexmap->begin()->second;
      if (vtx)
      {
        vtx_z = vtx->get_z();
      }
    }
  }



  //////////////////////////////////////////////
  //         towers 
  vector<float> ht_eta;
  vector<float> ht_phi;
  float tower_tot_e = 0;

  // if (!m_vtxCut || abs(vtx_z) > _vz)  return Fun4AllReturnCodes::EVENT_OK;

  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  if (towers)
  {
    int size = towers->size();  // online towers should be the same!
    for (int channel = 0; channel < size; channel++)
    {
      TowerInfo* tower = towers->get_tower_at_channel(channel);
      float offlineenergy = tower->get_energy();
      unsigned int towerkey = towers->encode_key(channel);
      int ieta = towers->getTowerEtaBin(towerkey);
      int iphi = towers->getTowerPhiBin(towerkey);
      bool isGood = !(tower->get_isBadChi2());
      tower_tot_e += offlineenergy;
      h_tower_e->Fill(offlineenergy);
      if (!isGood && offlineenergy > 0.2)
      {
        ht_eta.push_back(ieta);
        ht_phi.push_back(iphi);
      }
      if (isGood) h_emcal_e_eta->Fill(ieta, offlineenergy);
      if (offlineenergy > emcal_hit_threshold)
      {
        h_cemc_etaphi->Fill(ieta, iphi);
      }
    }
  }


  float weight =1;

  /////////////////////////////////////////////////
  //// Truth info
  /////////////////////////////////////////////////
  float wieght = 1;
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  vector <TLorentzVector> truth_photons;
  vector <TLorentzVector> truth_meson_photons;
  if (truthinfo)
  {
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
    {
      // Get truth particle
      const PHG4Particle* truth = iter->second;
      if (!truthinfo->is_primary(truth)) continue;
      TLorentzVector myVector;
      myVector.SetXYZT(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());

      float energy = myVector.E();
      h_truth_eta->Fill(myVector.Eta());
      h_truth_e->Fill(energy, wieght);
      weight = myVector.Pt()*TMath::Exp(-3*myVector.Pt());
      h_truth_pt->Fill(myVector.Pt(),weight);
      if (debug) std::cout << "primary pid=" <<  truth->get_pid()  << "   E=" << energy << "  pt=" << myVector.Pt() << "  eta=" << myVector.Eta() << "  phi=" << myVector.Phi() << std::endl;
       truth_photons.push_back(myVector);
    }
    //if ( truth_photons.at(0).M() > 0.1 && (  truth_photons.at(0).Pt() < 2 || truth_photons.at(0).Pt() > 4 )) return Fun4AllReturnCodes::EVENT_OK; 


    ////////////////
    // secondaries 
    ///////////////
    PHG4TruthInfoContainer::Range second_range = truthinfo->GetSecondaryParticleRange();
    float m_g4 = 0;

    for (PHG4TruthInfoContainer::ConstIterator siter = second_range.first;
          siter != second_range.second; ++siter) {
      if (m_g4 >= 19999) break;
      // Get photons from pi0 decays 
      const PHG4Particle *truth = siter->second;

      if (truth->get_pid() == 22 ) {
        PHG4Particle *parent = truthinfo->GetParticle(truth->get_parent_id());
        if (parent->get_pid() == 111 || parent->get_pid() == 221 ) {
          float phot_pt = sqrt(truth->get_px() * truth->get_px()
                            + truth->get_py() * truth->get_py());
          if  (phot_pt < 0.1) continue;
          //float phot_pz = truth->get_pz();
          float phot_e = truth->get_e();
          float phot_phi = atan2(truth->get_py(), truth->get_px());
          float phot_eta = atanh(truth->get_pz() / sqrt(truth->get_px()*truth->get_px()+truth->get_py()*truth->get_py()+truth->get_pz()*truth->get_pz()));

          TLorentzVector myVector;
          myVector.SetXYZT(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
          truth_meson_photons.push_back(myVector);

          if (debug) std::cout<< "2nd photons  pt=" <<  phot_pt <<   " e=" << phot_e << " phi=" << phot_phi << " eta=" << phot_eta << endl; 
       }
      }
    }
   
    PHG4TruthInfoContainer::VtxRange vtxrange = truthinfo->GetVtxRange();
    int n_vertex = 0;
    for (PHG4TruthInfoContainer::ConstVtxIterator iter = vtxrange.first; iter != vtxrange.second; ++iter) {
       PHG4VtxPoint *vtx = iter->second;
       h_vert_xy->Fill(vtx->get_x(),vtx->get_y());
       if ( vtx->get_id() == 1 ) vtx_z = vtx->get_z();
       if ( vtx->get_id() == 1 && false)
         std::cout << "vx=" <<  vtx->get_x() << "  vy=" << vtx->get_y() << "   vz=" << vtx->get_z() << "  id=" << vtx->get_id() << std::endl;
       n_vertex++;
       if (n_vertex >= 100000) break;
    }
  }

  if (fabs(vtx_z) > _vz) return Fun4AllReturnCodes::EVENT_OK; 

  //////////////////////////////////////////
  // geometry for hot tower/cluster masking
  std::string towergeomnodename = "TOWERGEOM_CEMC";
  RawTowerGeomContainer* m_geometry = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!m_geometry)
  {
    std::cout << Name() << "::"
              << "CreateNodeTree"
              << ": Could not find node " << towergeomnodename << std::endl;
    throw std::runtime_error("failed to find TOWERGEOM node in RawClusterDeadHotMask::CreateNodeTree");
  }


  //////////////////////////////////////////
  // Cluster counting

  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
  //RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
  if (!clusterContainer)
  {
    std::cout << PHWHERE << "funkyCaloStuff::process_event - Fatal Error - CLUSTER_CEMC node is missing. " << std::endl;
    return 0;
  }

  int nClusContainer = 0;

  RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
  RawClusterContainer::ConstIterator clusterIter;
  RawClusterContainer::ConstIterator clusterIter2;
  int nClusCount = 0;
  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
  {
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(0, 0, vtx_z);
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetEVec(*recoCluster, vertex);

    float clus_pt = E_vec_cluster.perp();
    float clus_chisq = recoCluster->get_chi2();

    if (debug && clus_pt > 0.1) std::cout << "clus #" << nClusCount  << "  E=" << E_vec_cluster.mag() << "    eta=" << E_vec_cluster.pseudoRapidity() <<  "  chi2=" << clus_chisq<< std::endl;

    nClusContainer++;
    if (clus_pt < nClus_ptCut) continue;


    if (clus_chisq > clus_chisq_cut) continue;

    nClusCount++;
  }

  if(debug) std::cout << "tower tot E=" << tower_tot_e << "    nClusContainer=" << nClusContainer << std::endl;

  h_nclusters->Fill(nClusCount);

  if (nClusCount > max_nClusCount) return Fun4AllReturnCodes::EVENT_OK;

  float ptMaxCut = 100;

  float pt1ClusCut = 0.9;  // 1.3
  float pt2ClusCut = 0.9;  // 0.7
  float pi0ptcut = 0; //1.22 * (pt1ClusCut + pt2ClusCut);


  
  //////////////////////////////////////////
  // clusters
  ////////////////////////////////////////////

  vector<float> save_pt;
  vector<float> save_eta;
  vector<float> save_phi;
  vector<float> save_e;

  float smear = 0.00;
  bool match1 = false;
  bool match2 = false;

  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
  {
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(0, 0, vtx_z);
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetEVec(*recoCluster, vertex);

    float clusE      = E_vec_cluster.mag();
    float clus_eta   = E_vec_cluster.pseudoRapidity();
    float clus_phi   = E_vec_cluster.phi();
    float clus_chisq = recoCluster->get_chi2();
    float clus_pt    = E_vec_cluster.perp();
    clus_pt *= rnd->Gaus(1,smear);

    int lt_eta = recoCluster->get_lead_tower().first; 
    int lt_phi = recoCluster->get_lead_tower().second;

    h_clusE->Fill(clusE);
    if (clusE < 0.1) continue;

    if (clus_chisq > clus_chisq_cut) continue;

    // loop over the towers in the cluster
    RawCluster::TowerConstRange towerCR = recoCluster->get_towers();
    RawCluster::TowerConstIterator toweriter;
    int nTow = 0;
    std::vector<int> tower_etas, tower_phis;
    std::vector<float> tower_energies;
    for (toweriter = towerCR.first; toweriter != towerCR.second; ++toweriter)
    {
      nTow++;
      tower_etas.push_back( m_geometry->get_tower_geometry(toweriter->first)->get_bineta());
      tower_phis.push_back(m_geometry->get_tower_geometry(toweriter->first)->get_binphi());
      tower_energies.push_back(toweriter->second);
    }

    std::pair<float,float> blockCord = getBlockCord(tower_etas,tower_phis,tower_energies);
    int block_eta_bin = h_block_bin->FindBin(blockCord.first) -1;
    int block_phi_bin = h_block_bin->FindBin(blockCord.second)-1;

    h3_cluster_block_cord_pt->Fill(blockCord.first,blockCord.second,clusE);

    h_clusE_nTow->Fill(clusE,nTow);

    if (lt_eta > 95) continue;


    h_pt_eta->Fill(clus_pt,lt_eta);
    h_pt_eta_rw->Fill(clus_pt,lt_eta,weight);
    h_etaphi_clus->Fill(clus_eta, clus_phi);

    TLorentzVector photon1;
    photon1.SetPtEtaPhiE(clus_pt, clus_eta, clus_phi, clusE);

    for (auto tr_phot : truth_photons){
      float delR = photon1.DeltaR(tr_phot);
      h_delR_recTrth->Fill(delR);
      float res = photon1.E()/tr_phot.E();
      float delPhi = photon1.Phi()-tr_phot.Phi();
      if (delPhi > TMath::TwoPi()) delPhi -= TMath::TwoPi();
      if (delPhi < - TMath::TwoPi()) delPhi += TMath::TwoPi();
      if (delR < 0.03 && res < 1.3 && res > 0.3){
        if(debug) std::cout << "match clusE=" << photon1.E() << "  truthE=" << tr_phot.E() <<  " delPhi=" << delPhi << std::endl;
        h_matched_res->Fill(res,photon1.Eta());
        h_res_e->Fill(res,photon1.E());
        h_res_e_eta->Fill(res,tr_phot.E(),lt_eta);
        h_res_e_phi->Fill(res,tr_phot.E(),lt_phi);
        h_res->Fill(res);
        h_delEta_e_eta->Fill(photon1.Eta()-tr_phot.Eta(),tr_phot.E(),lt_eta);
        h_delR_e_eta->Fill(delR,tr_phot.E(),lt_eta);
        h_delPhi_e_eta->Fill(delPhi,tr_phot.E(),lt_eta);
        h_delPhi_e_phi->Fill(delPhi,tr_phot.E(),lt_phi);
        h_truthE->Fill(tr_phot.E());
        if( block_eta_bin >= 0 && block_eta_bin < NBinsBlock && block_phi_bin >= 0 && block_phi_bin < NBinsBlock){
          h_res_block_E[block_eta_bin][block_phi_bin]->Fill(res,tr_phot.E()); 
        }
      }
          
    }

    TLorentzVector ph1_trEtaPhi;
    ph1_trEtaPhi.SetPtEtaPhiE(0,0, 0, 0);
    for (auto tr_phot : truth_meson_photons){
      float delR = photon1.DeltaR(tr_phot);
      float res = photon1.E()/tr_phot.E();
      if (delR < 0.03 && res < 1.5 && res > 0.7){
        ph1_trEtaPhi.SetPtEtaPhiE(clusE/TMath::CosH(tr_phot.Eta()),tr_phot.Eta(), tr_phot.Phi(),clusE);
        if(debug) std::cout << "match  eta=" << ph1_trEtaPhi.Eta() << " E=" << ph1_trEtaPhi.E() << std::endl;
        match1 = true; 
       break;
      }
    }

    if (clus_pt < pt1ClusCut || clus_pt > ptMaxCut) continue;

    for (clusterIter2 = clusterEnd.first; clusterIter2 != clusterEnd.second; clusterIter2++)
    {
      if (clusterIter == clusterIter2)
      {
        continue;
      }
      RawCluster* recoCluster2 = clusterIter2->second;

      CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetEVec(*recoCluster2, vertex);

      float clus2E = E_vec_cluster2.mag();
      float clus2_eta = E_vec_cluster2.pseudoRapidity();
      float clus2_phi = E_vec_cluster2.phi();
      float clus2_pt = E_vec_cluster2.perp();
      //clus2_pt *= rnd->Gaus(1,smear);
      float clus2_chisq = recoCluster2->get_chi2();

      if (clus2_pt < pt2ClusCut || clus2_pt > ptMaxCut) continue;
      if (clus2_chisq > clus_chisq_cut) continue;


      TLorentzVector photon2;
      photon2.SetPtEtaPhiE(clus2_pt, clus2_eta, clus2_phi, clus2E);

      if (fabs(clusE - clus2E) / (clusE + clus2E) > maxAlpha) continue;


      TLorentzVector ph2_trEtaPhi;
      ph2_trEtaPhi.SetPtEtaPhiE(0,0, 0, 0);
      for (auto tr_phot : truth_meson_photons){
        float delR = photon2.DeltaR(tr_phot);
        float res = photon2.E()/tr_phot.E();
        if (delR < 0.02 && res < 1.5 && res > 0.7){
         ph2_trEtaPhi.SetPtEtaPhiE(clus2E/TMath::CosH(tr_phot.Eta()),tr_phot.Eta(), tr_phot.Phi(),clus2E);
         if(debug) std::cout << "match  eta=" << ph2_trEtaPhi.Eta() << " E=" << ph2_trEtaPhi.E() << std::endl;
         if (match1) match2 = true;
        }
      }
          
      TLorentzVector pi0_trKin = ph1_trEtaPhi + ph2_trEtaPhi;

      TLorentzVector pi0 = photon1 + photon2;
      if (pi0.Pt() < pi0ptcut) continue;

      h_pt1->Fill(photon1.Pt());
      h_pt2->Fill(photon2.Pt());
      h_InvMass->Fill(pi0.M());
      h_InvMass_w->Fill(pi0.M(),weight);
      h_mass_eta_lt->Fill(pi0.M(),lt_eta);
      h_mass_eta_lt_rw->Fill(pi0.M(),lt_eta);
      h_m_pt_eta->Fill(pi0.M(),pi0.E(),lt_eta);
      if( block_eta_bin >= 0 && block_eta_bin < NBinsBlock && block_phi_bin >= 0 && block_phi_bin < NBinsBlock){
          h_mass_block_pt[block_eta_bin][block_phi_bin]->Fill(pi0.M(),pi0.E());
      }
      
      if (match2 && pi0_trKin.M() > 0.001){
        h_m_ptTr_eta->Fill(pi0.M(),truth_meson_photons.at(0).E(),lt_eta);
        h_m_ptTr_eta_trKin->Fill(pi0_trKin.M(),truth_meson_photons.at(0).E(),lt_eta);
        //std::cout << pi0_trKin.M() << std::endl;
      }
      
    }
  }  // clus1 loop


  ht_phi.clear();
  ht_eta.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}


int CaloAna::End(PHCompositeNode* /*topNode*/)
{
  outfile->cd();

  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");
  return 0;
}

float CaloAna::getWeight(int ieta, float pt){
  if (ieta < 0 || ieta > 95) return 0;
  float val = h_pt_rw[ieta]->GetBinContent(h_pt_rw[ieta]->FindBin(pt));
  if (val==0) return 0;
  return 1/val;
}


TF1* CaloAna::fitHistogram(TH1* h)
{
  TF1* fitFunc = new TF1("fitFunc", "[0]/[2]/2.5*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x^2 + [6]*x^3", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

  fitFunc->SetParameter(0, 5);
  fitFunc->SetParameter(1, target_pi0_mass);
  fitFunc->SetParameter(2, 0.01);
  fitFunc->SetParameter(3, 0.0);
  fitFunc->SetParameter(4, 0.0);
  fitFunc->SetParameter(5, 0.0);
  fitFunc->SetParameter(6, 100);

  fitFunc->SetParLimits(0, 0,10);
  fitFunc->SetParLimits(1, 0.113, 0.25);
  fitFunc->SetParLimits(2, 0.01, 0.04);
  fitFunc->SetParLimits(3,-2 ,1 );
  fitFunc->SetParLimits(4,0 ,40 );
  fitFunc->SetParLimits(5, -150,50 );
  fitFunc->SetParLimits(6, 0,200 );

  fitFunc->SetRange(0.05, 0.7);

  // Perform the fit
  h->Fit("fitFunc", "QN");

  return fitFunc;
}


void CaloAna::fitEtaSlices(const std::string& infile, const std::string& fitOutFile, const std::string& cdbFile)
{
  TFile* fin = new TFile(infile.c_str());
  std::cout << "getting hists" << std::endl;
  TH1F* h_peak_eta = new TH1F("h_peak_eta", "", 96, 0, 96);
  TH1F* h_sigma_eta = new TH1F("h_sigma_eta", "", 96, 0, 96);
  TH1F* h_p3_eta = new TH1F("h_p3_eta", "", 96, 0, 96);
  TH1F* h_p4_eta = new TH1F("h_p4_eta", "", 96, 0, 96);
  TH1F* h_p5_eta = new TH1F("h_p5_eta", "", 96, 0, 96);
  TH1F* h_p6_eta = new TH1F("h_p6_eta", "", 96, 0, 96);
  TH1F* h_p0_eta = new TH1F("h_p0_eta", "", 96, 0, 96);
  if (!fin)
  {
    std::cout << "CaloAna::fitEtaSlices null fin" << std::endl;
    exit(1);
  }
  TH1F* h_M_eta[96];
  for (int i = 0; i < 96; i++)
  {
    h_M_eta[i] = (TH1F*) fin->Get(Form("h_mass_eta_lt_rw%d", i));
    h_M_eta[i]->Scale(1./h_M_eta[i]->Integral(),"width");
  }

  TF1* fitFunOut[96];
  for (int i = 0; i < 96; i++)
  {
    if (!h_M_eta[i])
    {
      std::cout << "CaloAna::fitEtaSlices null hist" << std::endl;
    }

    fitFunOut[i] = fitHistogram(h_M_eta[i]);
    fitFunOut[i]->SetName(Form("f_pi0_eta%d",i));
    float mass_val_out = fitFunOut[i]->GetParameter(1);
    float mass_err_out = fitFunOut[i]->GetParError(1);
    h_peak_eta->SetBinContent(i + 1, mass_val_out);
    if (isnan(h_M_eta[i]->GetEntries())){
       h_peak_eta->SetBinError(i + 1, 0);
       continue;
    }
    h_peak_eta->SetBinError(i + 1, mass_err_out);
    h_sigma_eta->SetBinContent(i+1,fitFunOut[i]->GetParameter(2));
    h_sigma_eta->SetBinError(i+1,fitFunOut[i]->GetParError(2));
    h_p3_eta->SetBinContent(i+1,fitFunOut[i]->GetParameter(3));
    h_p3_eta->SetBinError(i+1,fitFunOut[i]->GetParError(3));
    h_p4_eta->SetBinContent(i+1,fitFunOut[i]->GetParameter(4));
    h_p4_eta->SetBinError(i+1,fitFunOut[i]->GetParError(4));
    h_p5_eta->SetBinContent(i+1,fitFunOut[i]->GetParameter(5));
    h_p5_eta->SetBinError(i+1,fitFunOut[i]->GetParError(5));
    h_p6_eta->SetBinContent(i+1,fitFunOut[i]->GetParameter(6));
    h_p6_eta->SetBinError(i+1,fitFunOut[i]->GetParError(6));
    h_p0_eta->SetBinContent(i+1,fitFunOut[i]->GetParameter(0));
    h_p0_eta->SetBinError(i+1,fitFunOut[i]->GetParError(0));
  }

  CDBTTree* cdbttree1 = new CDBTTree(cdbFile.c_str());
  CDBTTree* cdbttree2 = new CDBTTree(cdbFile.c_str());

  std::string m_fieldname = "Femc_datadriven_qm1_correction";

  for (int i = 0; i < 96; i++)
  {
    for (int j = 0; j < 256; j++)
    {
      float correction = target_pi0_mass / h_peak_eta->GetBinContent(i + 1);
      unsigned int key = TowerInfoDefs::encode_emcal(i, j);
      float val1 = cdbttree1->GetFloatValue(key, m_fieldname);
      cdbttree2->SetFloatValue(key, m_fieldname, val1 * correction);
    }
  }

  cdbttree2->Commit();
  cdbttree2->WriteCDBTTree();
  delete cdbttree2;
  delete cdbttree1;

  TFile* fit_out = new TFile(fitOutFile.c_str(), "recreate");
  fit_out->cd();
  for (auto& i : h_M_eta)
  {
    i->Write();
    delete i;
  }
  for (auto& i : fitFunOut)
  {
    i->Write();
    delete i;
  }

  h_p3_eta->Write();
  h_p4_eta->Write();
  h_p5_eta->Write();
  h_p6_eta->Write();
  h_p0_eta->Write();
  h_sigma_eta->Write();
  h_peak_eta->Write();
  fin->Close();

  std::cout << "finish fitting suc" << std::endl;

  return;
}


std::pair<float,float> CaloAna::getBlockCordCorr(std::vector<int> etas, std::vector<int> phis, std::vector<float> Es){
  std::pair<float,float> raw_cord = getBlockCord(etas, phis, Es);
  float Xcg_eta = raw_cord.first < 0.5 ? raw_cord.first : raw_cord.first - 1;
  float Xcg_phi = raw_cord.second < 0.5 ? raw_cord.second : raw_cord.second - 1;

  float b_eta = 0.155;
  float b_phi= 0.11;

  TF1 finv = TF1("finv", "[0]*TMath::ASinH(2*x*sinh(1/(2*[0])))", -0.5, 0.5);

  finv.SetParameter(0,b_eta);
  float t_eta = finv.Eval(Xcg_eta);

  finv.SetParameter(0,b_phi);
  float t_phi = finv.Eval(Xcg_phi);

  std::pair<float,float> res;
  res.first  =  raw_cord.first  < 0.5 ? t_eta : t_eta+1;
  res.second =  raw_cord.second < 0.5 ? t_phi : t_phi+1;

  return res;
}

std::pair<float,float> CaloAna::getBlockCord(std::vector<int> etas, std::vector<int> phis, std::vector<float> Es){
  float avgEta = getAvgEta(etas,Es);
  float avgPhi = getAvgPhi(phis,Es);
  int avgIPhi = std::floor(avgPhi);
  int avgIEta = std::floor(avgEta);
  int iblockphi = avgIPhi/2;
  int iblocketa = avgIEta/2;
  float interBlockPhi = avgPhi - iblockphi*2;
  float interBLockEta = avgEta - iblocketa*2;
  if (interBLockEta > 1.5) interBLockEta -= 2;
  std::pair<float,float> blockCord = {interBLockEta,interBlockPhi};
  //std::cout << "-----------------" << std::endl <<  "cluster size=" << Es.size() << std::endl <<  "avgEtaPhi=" << avgEta << "," << avgPhi << std::endl  << "  iblock=" << blockCord.first << "," << blockCord.second << std::endl;
  return blockCord;
}




float CaloAna::getAvgEta(const std::vector<int> &toweretas, const std::vector<float> &towerenergies) {
    float etamult = 0;
    float etasum = 0;
    if(toweretas.size()==1) return toweretas[0];
    for (UInt_t j = 0; j < towerenergies.size(); j++) {
        float energymult = towerenergies.at(j) * toweretas.at(j);
        etamult += energymult;
        etasum += towerenergies.at(j);
    }
    return etamult / etasum;
}


float CaloAna::getAvgPhi(const std::vector<int> &towerphis, const std::vector<float> &towerenergies) {
    int nphibin=256;
    float phimult = 0;
    float phisum  = 0;

    for (UInt_t j = 0; j < towerenergies.size(); j++) {
        int phibin = towerphis.at(j);
        if (phibin - towerphis.at(0) < -nphibin / 2.0) {
            phibin += nphibin;
        }

        else if (phibin - towerphis.at(0) > +nphibin / 2.0) {
            phibin -= nphibin;
        }
        assert(std::abs(phibin - towerphis.at(0)) <= nphibin / 2.0);
        float energymult = towerenergies.at(j) * phibin;
        phimult += energymult;
        phisum += towerenergies.at(j);
    }

    float avgphi = phimult / phisum;
   if (avgphi < 0) {
       avgphi += nphibin;
   }
   avgphi = fmod(avgphi, nphibin);
   if (avgphi >= 255.5) {
       avgphi -= nphibin;
   }
   avgphi = fmod(avgphi + 0.5, 2) - 0.5;
   return avgphi;
}
