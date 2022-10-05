#define BackgroundSources_cxx
// The class definition in BackgroundSources.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("BackgroundSources.C")
// root> T->Process("BackgroundSources.C","some options")
// root> T->Process("BackgroundSources.C+")
//


#include "BackgroundSources.h"
#include <TMath.h>
#include <cmath>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "utils.C"


using namespace std;


void BackgroundSources::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  cout << " --------------------------" << endl;
  cout << "     Background sources    " << endl;
  cout << " --------------------------" << endl;
}


void BackgroundSources::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  TString outFileName = option;

  if(outFileName.Contains("isMC")){
    isMC = true;
    outFileName.Resize(outFileName.Length()-5);
  }
  else isMC = false;

  // check if outputfile exists
  if(gSystem->AccessPathName(outFileName)){
    my_file = new TFile(outFileName, "RECREATE");  
  }
  else{
    my_file = new TFile(outFileName, "UPDATE");  
  }
  my_file->cd();

  // if MC, get the correct content of the GenPart branches
  if(isMC){
    nGenPart = {fReader, "nGenPart"};
    GenPart_eta = {fReader, "GenPart_eta"};
    GenPart_mass = {fReader, "GenPart_mass"};
    GenPart_phi = {fReader, "GenPart_phi"};
    GenPart_pt = {fReader, "GenPart_pt"};
    GenPart_vx = {fReader, "GenPart_vx"};
    GenPart_vy = {fReader, "GenPart_vy"};
    GenPart_vz = {fReader, "GenPart_vz"};
    GenPart_genPartIdxMother = {fReader, "GenPart_genPartIdxMother"};
    GenPart_pdgId = {fReader, "GenPart_pdgId"};
    GenPart_status = {fReader, "GenPart_status"};
    GenPart_statusFlags = {fReader, "GenPart_statusFlags"};
    Muon_genPartIdx = {fReader, "Muon_genPartIdx"};
    Pileup_nPU = {fReader, "Pileup_nPU"};
    Pileup_nTrueInt = {fReader, "Pileup_nTrueInt"};
  }

  // getting the signal tree ready
  sources = new TTree("sources", "sources");
  
  sources->Branch("hnl_mass", &hnl_mass);
  sources->Branch("hnl_pt", &hnl_pt);
  sources->Branch("hnl_charge", &hnl_charge);
  sources->Branch("trgmu_softid", &trgmu_softid);
  sources->Branch("mu_looseid", &mu_looseid);
  sources->Branch("pi_packedcandhashighpurity", &pi_packedcandhashighpurity);
  sources->Branch("trgmu_charge", &trgmu_charge);
  sources->Branch("mu_charge", &mu_charge);
  sources->Branch("trgmu_mu_mass", &trgmu_mu_mass);
  sources->Branch("trgmu_pi_mass", &trgmu_pi_mass);
  sources->Branch("mu0_pt", &mu0_pt);
  sources->Branch("mu0_eta", &mu0_eta);
  sources->Branch("mu0_dxysig", &mu0_dxysig);
  sources->Branch("mu_pt", &mu_pt);
  sources->Branch("mu_eta", &mu_eta);
  sources->Branch("mu_dxysig", &mu_dxysig);
  sources->Branch("pi_pt", &pi_pt);
  sources->Branch("pi_eta", &pi_eta);
  sources->Branch("pi_dxysig", &pi_dxysig);
  sources->Branch("sv_lxy", &sv_lxy);
  sources->Branch("sv_lxyz", &sv_lxyz);
  sources->Branch("sv_lxysig", &sv_lxysig);
  sources->Branch("mu_dxysig", &mu_dxysig);
  sources->Branch("pi_dxysig", &pi_dxysig);
  sources->Branch("hnl_cos2d", &hnl_cos2d);
  
  sources->Branch("mu0_ismatched", &mu0_ismatched);
  sources->Branch("mu0_isfake", &mu0_isfake);
  sources->Branch("mu0_genmatching_deltar", &mu0_genmatching_deltar);
  sources->Branch("mu0_genmatching_deltaptrel", &mu0_genmatching_deltaptrel);

  sources->Branch("mu0_directBdecay", &mu0_directBdecay);
  sources->Branch("mu0_BtoDdecay", &mu0_BtoDdecay);
  sources->Branch("mu0_BtoDtoKdecay", &mu0_BtoDtoKdecay);
  sources->Branch("mu0_Btotaudecay", &mu0_Btotaudecay);
  sources->Branch("mu0_BtoJPsidecay", &mu0_BtoJPsidecay);
  sources->Branch("mu0_BtoDexciteddecay", &mu0_BtoDexciteddecay);
  sources->Branch("mu0_nonBdecay", &mu0_nonBdecay);
  sources->Branch("mu0_other", &mu0_other);

  sources->Branch("mu_ismatched", &mu_ismatched);
  sources->Branch("mu_isfake", &mu_isfake);
  sources->Branch("mu_genmatching_deltar", &mu_genmatching_deltar);
  sources->Branch("mu_genmatching_deltaptrel", &mu_genmatching_deltaptrel);

  sources->Branch("mu_directBdecay", &mu_directBdecay);
  sources->Branch("mu_BtoDdecay", &mu_BtoDdecay);
  sources->Branch("mu_BtoDtoKdecay", &mu_BtoDtoKdecay);
  sources->Branch("mu_Btotaudecay", &mu_Btotaudecay);
  sources->Branch("mu_BtoJPsidecay", &mu_BtoJPsidecay);
  sources->Branch("mu_BtoDexciteddecay", &mu_BtoDexciteddecay);
  sources->Branch("mu_nonBdecay", &mu_nonBdecay);
  sources->Branch("mu_other", &mu_other);

  sources->Branch("pi_ismatched", &pi_ismatched);
  sources->Branch("pi_isfake", &pi_isfake);
  sources->Branch("pi_genmatching_deltar", &pi_genmatching_deltar);
  sources->Branch("pi_genmatching_deltaptrel", &pi_genmatching_deltaptrel);

  sources->Branch("pi_directBdecay", &pi_directBdecay);
  sources->Branch("pi_BtoDdecay", &pi_BtoDdecay);
  sources->Branch("pi_BtoDtoKdecay", &pi_BtoDtoKdecay);
  sources->Branch("pi_Btotaudecay", &pi_Btotaudecay);
  sources->Branch("pi_BtoJPsidecay", &pi_BtoJPsidecay);
  sources->Branch("pi_BtoDexciteddecay", &pi_BtoDexciteddecay);
  sources->Branch("pi_nonBdecay", &pi_nonBdecay);
  sources->Branch("pi_other", &pi_other);

  sources->Branch("cand_isnotmatched", &cand_isnotmatched);
  sources->Branch("cand_isfullmatched", &cand_isfullmatched);
  sources->Branch("cand_ispartialmatched", &cand_ispartialmatched);
  sources->Branch("cand_mu0mu_samemother", &cand_mu0mu_samemother);
  sources->Branch("cand_mu0pi_samemother", &cand_mu0pi_samemother);
  sources->Branch("cand_mupi_samemother", &cand_mupi_samemother);
  sources->Branch("cand_mumupi_samemother", &cand_mumupi_samemother);
}


std::vector<pair<int, float>> BackgroundSources::getGenCandidateMap(float part_pt, float part_eta, float part_phi, float deltaR_max, float deltaPtRel_max, int nGen) const{
  /*
   * Function that returns the list of possible gen-matched candidates (index + deltaR)
   */

  std::vector<pair<int, float>> cand_genidx_deltaR;

  for(int iGen(0); iGen<nGen; ++iGen){
    if(abs(GenPart_pdgId[iGen])!=13 && abs(GenPart_pdgId[iGen])!=211 && abs(GenPart_pdgId[iGen])!=321) continue;
    float gen_pt = GenPart_pt[iGen];
    float gen_eta = GenPart_eta[iGen];
    float gen_phi = GenPart_phi[iGen];
    float deltaR = reco::deltaR(part_eta, part_phi, gen_eta, gen_phi);
    float deltaPtRel = fabs(part_pt - gen_pt) / gen_pt;
  
    // impose condition on deltaR and deltaPtRel
    if(deltaR < deltaR_max && deltaPtRel < deltaPtRel_max){
      //std::cout << "reco muon pt " << part_pt << " eta " << part_eta << " phi " << part_phi << " ismatched to gen " << iGen << " pt " << gen_pt << " eta " << gen_eta << " phi " << gen_phi << " pdgid " << GenPart_pdgId[iGen] << " deltaR " << deltaR << " deltaPtRel " << deltaPtRel << std::endl;
      pair<int, float> cand_genidx_deltaR_tmp;
      cand_genidx_deltaR_tmp.first = iGen;
      cand_genidx_deltaR_tmp.second = deltaR;
      cand_genidx_deltaR.push_back(cand_genidx_deltaR_tmp);
    }
  }

  // sort gen candidates in increasing deltaR
  sort(cand_genidx_deltaR.begin(), cand_genidx_deltaR.end(), [](const pair<int, float> &pair_i, const pair<int, float> &pair_j){
    return pair_i.second < pair_j.second;
  });

  return cand_genidx_deltaR;
}



std::vector<int> BackgroundSources::getAncestor(int particle_idx, std::vector<int> list) const{
  /*
   * Recursive function to get the decay chain starting from a given particle
   */

  // find the mother
  int mother_idx = GenPart_genPartIdxMother[particle_idx];
  int mother_pdgid = GenPart_pdgId[mother_idx];
  int mother_eta = GenPart_eta[mother_idx];

  //std::cout << "mother is " << mother_pdgid << std::endl;

  if(mother_pdgid == 0 || fabs(mother_pdgid) < 8 || fabs(mother_pdgid) == 21 || fabs(mother_eta) > 1000){
    //TODO add condition on diquarks and excitedBs and oscillations
    std::cout << "--- End of decay chain ---" << std::endl;
    return list;
  }
  else{
    std::cout << "mother is " << mother_pdgid << std::endl;
    vector<int> list_tmp = list;
    list_tmp.push_back(mother_pdgid);
    return getAncestor(mother_idx, list_tmp);
  }
}

Int_t BackgroundSources::getBmesonIdx(int particle_idx, int b_meson_idx) const{
  /*
   * Recursive function to get the decay chain starting from a given particle
   */

  // find the mother
  int mother_idx = GenPart_genPartIdxMother[particle_idx];
  int mother_pdgid = GenPart_pdgId[mother_idx];
  int mother_eta = GenPart_eta[mother_idx];

  if(mother_pdgid == 0 || fabs(mother_pdgid) < 8 || fabs(mother_pdgid) == 21 || fabs(mother_eta) > 1000){
    return b_meson_idx;
  }
  else{
    b_meson_idx = mother_idx;
    return getBmesonIdx(mother_idx, b_meson_idx);
  }
}


Bool_t BackgroundSources::checkCommonAncestor(int mu0_idx, int mu_idx, int pi_idx) const{
  /*
   * Function that checks if particles are originating from common B meson
   */

  int initialisator_mu0 = -1;
  int initialisator_mu = -1;
  int initialisator_pi = -1;

  int b_meson_mu0_idx = getBmesonIdx(mu0_idx, initialisator_mu0);
  int b_meson_mu_idx = getBmesonIdx(mu_idx, initialisator_mu);
  int b_meson_pi_idx = getBmesonIdx(pi_idx, initialisator_pi);

  Bool_t flag = false;
  if(b_meson_mu0_idx == b_meson_mu_idx && b_meson_mu_idx == b_meson_pi_idx){
    flag = true;
  }

  return flag;
}


Bool_t BackgroundSources::checkCommonAncestor(int part1_idx, int part2_idx) const{
  /*
   * Function that checks if particles are originating from common B meson
   */

  int initialisator_part1 = -1;
  int initialisator_part2 = -1;

  int b_meson_part1_idx = getBmesonIdx(part1_idx, initialisator_part1);
  int b_meson_part2_idx = getBmesonIdx(part2_idx, initialisator_part2);

  Bool_t flag = false;
  if(b_meson_part1_idx == b_meson_part2_idx){
    flag = true;
  }

  return flag;
}


Bool_t isPartOf(int pdgid, vector<int> list){
  Bool_t flag = false;
  for(unsigned int i(0); i<list.size(); ++i){
    if(fabs(pdgid) == list[i]){
      flag = true;
      break;
    }
  }
  return flag;
}
  

Bool_t BackgroundSources::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetLocalEntry(entry);
  //if(entry>1939) return false;
  cout << endl << "--- Entry " << entry << " ---" << endl;

  // number of candidates in the event
  UInt_t nCand_sig = *nBToMuMuPi; 

  //   ----- Signal Channel -----  //

  if(nCand_sig > 0){ // at least one candidate per event
    // selecting the candidate as the one having the largest hnl pt
    // - create candIdx - cos2d pairs
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_cos2D);
    // - sort it in decreasing cos2d
    stable_sort(pair_candIdx_desc_cos2d_sig.begin(), pair_candIdx_desc_cos2d_sig.end(), sortcansbydesc);

    // - then privilege OS cand over SS ones
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sig, BToMuMuPi_hnl_charge);
    stable_sort(pair_candIdx_desc_cos2d_sign_sig.begin(), pair_candIdx_desc_cos2d_sign_sig.end(), sortcansbydesc_opp);

    // - for signal, priviledge matched candidates
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_matched_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sign_sig, BToMuMuPi_isMatched);
    stable_sort(pair_candIdx_desc_cos2d_sign_matched_sig.begin(), pair_candIdx_desc_cos2d_sign_matched_sig.end(), sortcansbydesc);

    // - priviledge slimmed over dsa candidates
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_matched_muon_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sign_sig, BToMuMuPi_sel_mu_idx, Muon_isDSAMuon);
    stable_sort(pair_candIdx_desc_cos2d_sign_matched_muon_sig.begin(), pair_candIdx_desc_cos2d_sign_matched_muon_sig.end(), sortcansbydesc_opp);

    // - and select the candidate
    UInt_t selectedCandIdx_sig = pair_candIdx_desc_cos2d_sign_matched_sig[0].first;
  
    hnl_mass = BToMuMuPi_hnl_mass[selectedCandIdx_sig];
    hnl_pt = BToMuMuPi_hnl_pt[selectedCandIdx_sig];
    hnl_charge = BToMuMuPi_hnl_charge[selectedCandIdx_sig];
    trgmu_softid = Muon_softId[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
    mu_looseid = Muon_looseId[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
    if(BToMuMuPi_pi_ispacked[selectedCandIdx_sig] && BToMuMuPi_pi_highPurityFlag[selectedCandIdx_sig]){
      pi_packedcandhashighpurity = 1;
    }
    else{
      pi_packedcandhashighpurity = 0;
    }
    trgmu_charge = Muon_charge[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
    mu_charge = Muon_charge[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
    trgmu_mu_mass = BToMuMuPi_trgmu_mu_mass[selectedCandIdx_sig];
    trgmu_pi_mass = BToMuMuPi_trgmu_pi_mass[selectedCandIdx_sig];
    mu0_pt = BToMuMuPi_trg_mu_pt[selectedCandIdx_sig];
    mu0_eta = BToMuMuPi_trg_mu_eta[selectedCandIdx_sig];
    mu_pt = BToMuMuPi_fit_mu_pt[selectedCandIdx_sig];
    mu_eta = BToMuMuPi_fit_mu_eta[selectedCandIdx_sig];
    pi_pt = BToMuMuPi_fit_pi_pt[selectedCandIdx_sig];
    pi_eta = BToMuMuPi_fit_pi_eta[selectedCandIdx_sig];
    sv_lxy = BToMuMuPi_sv_lxy[selectedCandIdx_sig];
    sv_lxysig = BToMuMuPi_sv_lxy_sig[selectedCandIdx_sig];
    sv_lxyz = BToMuMuPi_sv_lxyz[selectedCandIdx_sig];
    mu0_dxysig = Muon_dxyS[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
    mu_dxysig = Muon_dxyS[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
    pi_dxysig = BToMuMuPi_pi_dxyS[selectedCandIdx_sig];
    hnl_cos2d = BToMuMuPi_hnl_cos2D[selectedCandIdx_sig];
   

    // --- 1. Study of the matching --- //

    // get the number of gen particles
    UInt_t nGen = *nGenPart;

    // those conditions are similar to those used for the signal gen-matching
    float deltaR_max = 0.015;
    float deltaPtRel_max = 0.1;

    int mu0_genmatched_idx = -1;
    int mu_genmatched_idx = -1;
    int pi_genmatched_idx = -1;

    // match the primary muon
    float mu0_pt = BToMuMuPi_trg_mu_pt[selectedCandIdx_sig];
    float mu0_eta = BToMuMuPi_trg_mu_eta[selectedCandIdx_sig];
    float mu0_phi = BToMuMuPi_trg_mu_phi[selectedCandIdx_sig];

    std::vector<pair<int, float>> cand_mu0_genidx_deltaR = getGenCandidateMap(mu0_pt, mu0_eta, mu0_phi, deltaR_max, deltaPtRel_max, nGen);

    if(cand_mu0_genidx_deltaR.size() == 0){
      mu0_ismatched = 0;
      //std::cout << "no match" << std::endl;
    }
    else{
      mu0_ismatched = 1;
      mu0_genmatched_idx = cand_mu0_genidx_deltaR[0].first;
      std::cout << "primary muon is first matched " << mu0_genmatched_idx << std::endl;
      if(abs(GenPart_pdgId[mu0_genmatched_idx])!=13) mu0_isfake = 1; 
      else mu0_isfake = 0;
      mu0_genmatching_deltar = cand_mu0_genidx_deltaR[0].second;
      mu0_genmatching_deltaptrel = fabs(mu0_pt - GenPart_pt[mu0_genmatched_idx]) / GenPart_pt[mu0_genmatched_idx];
    }

    // match the displaced muon
    float mu_pt = BToMuMuPi_fit_mu_pt[selectedCandIdx_sig];
    float mu_eta = BToMuMuPi_fit_mu_eta[selectedCandIdx_sig];
    float mu_phi = BToMuMuPi_fit_mu_phi[selectedCandIdx_sig];

    std::vector<pair<int, float>> cand_mu_genidx_deltaR = getGenCandidateMap(mu_pt, mu_eta, mu_phi, deltaR_max, deltaPtRel_max, nGen);

    if(cand_mu_genidx_deltaR.size() == 0){
      mu_ismatched = 0;
      //std::cout << "no match" << std::endl;
    }
    else{
      mu_ismatched = 1;
      mu_genmatched_idx = cand_mu_genidx_deltaR[0].first;
      std::cout << "muon is first matched " << mu_genmatched_idx << std::endl;
      if(abs(GenPart_pdgId[mu_genmatched_idx])!=13) mu_isfake = 1; 
      else mu_isfake = 0;
      mu_genmatching_deltar = cand_mu_genidx_deltaR[0].second;
      mu_genmatching_deltaptrel = fabs(mu_pt - GenPart_pt[mu_genmatched_idx]) / GenPart_pt[mu_genmatched_idx];
      //std::cout << cand_mu_genidx_deltaR[0].first << " " << cand_mu_genidx_deltaR[0].second << std::endl;
    }


    // match the displaced pion
    float pi_pt = BToMuMuPi_fit_pi_pt[selectedCandIdx_sig];
    float pi_eta = BToMuMuPi_fit_pi_eta[selectedCandIdx_sig];
    float pi_phi = BToMuMuPi_fit_pi_phi[selectedCandIdx_sig];

    std::vector<pair<int, float>> cand_pi_genidx_deltaR = getGenCandidateMap(pi_pt, pi_eta, pi_phi, deltaR_max, deltaPtRel_max, nGen);

    if(cand_pi_genidx_deltaR.size() == 0 || cand_pi_genidx_deltaR[0].first == -1){
      pi_ismatched = 0;
      //std::cout << "no match" << std::endl;
    }
    else{
      pi_ismatched = 1;
      pi_genmatched_idx = cand_pi_genidx_deltaR[0].first;
      std::cout << "pion is first matched " << pi_genmatched_idx << std::endl;
      if(abs(GenPart_pdgId[pi_genmatched_idx])!=211 && abs(GenPart_pdgId[pi_genmatched_idx])!=321) pi_isfake = 1; 
      else pi_isfake = 0;
      pi_genmatching_deltar = cand_pi_genidx_deltaR[0].second;
      pi_genmatching_deltaptrel = fabs(pi_pt - GenPart_pt[pi_genmatched_idx]) / GenPart_pt[pi_genmatched_idx];
      //std::cout << cand_pi_genidx_deltaR[0].first << " " << cand_pi_genidx_deltaR[0].second << std::endl;
    }

    // if gen particle matched to two reco particles, match it to the one that has the smallest deltaR
    // avoid double-counting between mu0 and mu
    if(cand_mu0_genidx_deltaR.size() !=0 && cand_mu_genidx_deltaR.size() != 0 && cand_mu0_genidx_deltaR[0].first == cand_mu_genidx_deltaR[0].first){
      if(cand_mu0_genidx_deltaR[0].second < cand_mu_genidx_deltaR[0].second){
        if(cand_mu_genidx_deltaR.size() > 1){
          mu_ismatched = 1;
          mu_genmatched_idx = cand_mu_genidx_deltaR[1].first;
          if(abs(GenPart_pdgId[mu_genmatched_idx])!=13) mu_isfake = 1; 
          else mu_isfake = 0;
          mu_genmatching_deltar = cand_mu_genidx_deltaR[1].second;
          mu_genmatching_deltaptrel = fabs(mu_pt - GenPart_pt[mu_genmatched_idx]) / GenPart_pt[mu_genmatched_idx];
        }
        else{
          mu_genmatched_idx = -1;
          mu_ismatched = 0;
        }
      }
      else if(cand_mu0_genidx_deltaR[0].second > cand_mu_genidx_deltaR[0].second){
        if(cand_mu0_genidx_deltaR.size() > 1){
          mu0_ismatched = 1;
          mu0_genmatched_idx = cand_mu0_genidx_deltaR[1].first;
          if(abs(GenPart_pdgId[mu0_genmatched_idx])!=13) mu0_isfake = 1; 
          else mu0_isfake = 0;
          mu0_genmatching_deltar = cand_mu0_genidx_deltaR[1].second;
          mu0_genmatching_deltaptrel = fabs(mu0_pt - GenPart_pt[mu0_genmatched_idx]) / GenPart_pt[mu0_genmatched_idx];
        }
        else{
          mu0_genmatched_idx = -1;
          mu0_ismatched = 0;
        }
      }
      //std::cout << "updated matching mu0: " << mu0_genmatched_idx << std::endl;
      //std::cout << "updated matching mu: " << mu_genmatched_idx << std::endl;
    }

    // avoid double-counting between mu0 and pi
    if(cand_mu0_genidx_deltaR.size() !=0 && cand_pi_genidx_deltaR.size() != 0 && cand_mu0_genidx_deltaR[0].first == cand_pi_genidx_deltaR[0].first){
      if(cand_mu0_genidx_deltaR[0].second < cand_pi_genidx_deltaR[0].second){
        if(cand_pi_genidx_deltaR.size() > 1){
          pi_ismatched = 1;
          pi_genmatched_idx = cand_pi_genidx_deltaR[1].first;
          if(abs(GenPart_pdgId[pi_genmatched_idx])!=13) pi_isfake = 1; 
          else pi_isfake = 0;
          pi_genmatching_deltar = cand_pi_genidx_deltaR[1].second;
          pi_genmatching_deltaptrel = fabs(pi_pt - GenPart_pt[pi_genmatched_idx]) / GenPart_pt[pi_genmatched_idx];
        }
        else{
          pi_genmatched_idx = -1;
          pi_ismatched = 0;
        }
      }
      else if(cand_mu0_genidx_deltaR[0].second > cand_pi_genidx_deltaR[0].second){
        if(cand_mu0_genidx_deltaR.size() > 1){
          mu0_ismatched = 1;
          mu0_genmatched_idx = cand_mu0_genidx_deltaR[1].first;
          if(abs(GenPart_pdgId[mu0_genmatched_idx])!=13) mu0_isfake = 1; 
          else mu0_isfake = 0;
          mu0_genmatching_deltar = cand_mu0_genidx_deltaR[1].second;
          mu0_genmatching_deltaptrel = fabs(mu0_pt - GenPart_pt[mu0_genmatched_idx]) / GenPart_pt[mu0_genmatched_idx];
        }
        else{
          mu0_genmatched_idx = -1;
          mu0_ismatched = 0;
        }
      }
      //std::cout << "updated matching mu0: " << mu0_genmatched_idx << std::endl;
      //std::cout << "updated matching pi " << pi_genmatched_idx << std::endl;
    }

    // avoid double-counting between mu and pi
    if(cand_mu_genidx_deltaR.size() !=0 && cand_pi_genidx_deltaR.size() != 0 && cand_mu_genidx_deltaR[0].first == cand_pi_genidx_deltaR[0].first){
      if(cand_mu_genidx_deltaR[0].second < cand_pi_genidx_deltaR[0].second){
        if(cand_pi_genidx_deltaR.size() > 1){
          pi_ismatched = 1;
          pi_genmatched_idx = cand_pi_genidx_deltaR[1].first;
          if(abs(GenPart_pdgId[pi_genmatched_idx])!=13) pi_isfake = 1; 
          else pi_isfake = 0;
          pi_genmatching_deltar = cand_pi_genidx_deltaR[1].second;
          pi_genmatching_deltaptrel = fabs(pi_pt - GenPart_pt[pi_genmatched_idx]) / GenPart_pt[pi_genmatched_idx];
        }
        else{
          pi_genmatched_idx = -1;
          pi_ismatched = 0;
        }
      }
      else if(cand_mu_genidx_deltaR[0].second > cand_pi_genidx_deltaR[0].second){
        if(cand_mu_genidx_deltaR.size() > 1){
          mu_ismatched = 1;
          mu_genmatched_idx = cand_mu_genidx_deltaR[1].first;
          if(abs(GenPart_pdgId[mu_genmatched_idx])!=13) mu_isfake = 1; 
          else mu_isfake = 0;
          mu_genmatching_deltar = cand_mu_genidx_deltaR[1].second;
          mu_genmatching_deltaptrel = fabs(mu_pt - GenPart_pt[mu_genmatched_idx]) / GenPart_pt[mu_genmatched_idx];
        }
        else{
          mu_genmatched_idx = -1;
          mu_ismatched = 0;
        }
      }
      //std::cout << "updated matching mu: " << mu_genmatched_idx << std::endl;
      //std::cout << "updated matching pi " << pi_genmatched_idx << std::endl;
    }

    // --- 2. Study of the sources --- //
      
    std::vector<int> list_decay_mu0;
    std::vector<int> list_decay_mu;
    std::vector<int> list_decay_pi;

    // case where the primary muon is matched
    if(mu0_ismatched){

      // get the ancestor chain
      std::cout << "-- Primary muon --" << std::endl;
      std::vector<int> initialisator_mu0;
      list_decay_mu0 = getAncestor(mu0_genmatched_idx, initialisator_mu0);

      mu0_directBdecay = 0;
      mu0_BtoDdecay = 0;
      mu0_BtoDtoKdecay = 0;
      mu0_Btotaudecay = 0;
      mu0_BtoJPsidecay = 0;
      mu0_BtoDexciteddecay = 0;
      mu0_nonBdecay = 0;
      mu0_other = 0;

      if(list_decay_mu0.size() == 0){
        // remove pathological events
        mu0_ismatched = 0;
      }
      else if(list_decay_mu0.size() == 1 && isPartOf(list_decay_mu0[0], B_mesons)){
        std::cout << "-> direct B decay" << std::endl;
        mu0_directBdecay = 1;
      }
      else if(list_decay_mu0.size() == 2 && isPartOf(list_decay_mu0[0], D_mesons) && isPartOf(list_decay_mu0[1], B_mesons)){
        std::cout << "-> BToD decay" << std::endl;
        mu0_BtoDdecay = 1;
      }
      else if(list_decay_mu0.size() == 3 && isPartOf(list_decay_mu0[0], D_mesons) && isPartOf(list_decay_mu0[1], K_mesons) && isPartOf(list_decay_mu0[2], B_mesons)){
        std::cout << "-> BToDtoK decay" << std::endl;
        mu0_BtoDtoKdecay = 1;
      }
      else if(list_decay_mu0.size() == 2 && fabs(list_decay_mu0[0]) == 15 && isPartOf(list_decay_mu0[1], B_mesons)){
        std::cout << "-> BTotau decay" << std::endl;
        mu0_Btotaudecay = 1;
      }
      else if(list_decay_mu0.size() == 2 && fabs(list_decay_mu0[0]) == 443 && isPartOf(list_decay_mu0[1], B_mesons)){
        //isPartOf(443, list_decay_mu0) || isPartOf(-443, list_decay_mu0)) && ){
        std::cout << "-> BToJpsi decay" << std::endl;
        mu0_BtoJPsidecay = 1;
      }
      else if(list_decay_mu0.size() > 2 && (isPartOf(list_decay_mu0[1], excitedD_mesons) || isPartOf(list_decay_mu0[list_decay_mu0.size()-2], excitedD_mesons)) && isPartOf(list_decay_mu0[list_decay_mu0.size()-1], B_mesons)){
        std::cout << "-> BToDexcited decay" << std::endl;
        mu0_BtoDexciteddecay = 1;
      }
      else if(!isPartOf(list_decay_mu0[list_decay_mu0.size()-1], B_mesons)){
        std::cout << "-> not a B decay" << std::endl;
        mu0_nonBdecay = 1;
      }
      else{
        std::cout << "-> other" << std::endl;
        mu0_other = 1;
      }
    
      std::cout << std::endl;
    }

    // case where the displaced muon is matched
    if(mu_ismatched){

      // get the ancestor chain
      std::cout << "-- Displaced muon --" << std::endl;
      std::vector<int> initialisator_mu;
      list_decay_mu = getAncestor(mu_genmatched_idx, initialisator_mu);

      mu_directBdecay = 0;
      mu_BtoDdecay = 0;
      mu_BtoDtoKdecay = 0;
      mu_Btotaudecay = 0;
      mu_BtoJPsidecay = 0;
      mu_BtoDexciteddecay = 0;
      mu_nonBdecay = 0;
      mu_other = 0;

      if(list_decay_mu.size() == 0){
        // remove pathological events
        mu_ismatched = 0;
      }
      else if(list_decay_mu.size() == 1 && isPartOf(list_decay_mu[0], B_mesons)){
        std::cout << "-> direct B decay" << std::endl;
        mu_directBdecay = 1;
      }
      else if(list_decay_mu.size() == 2 && isPartOf(list_decay_mu[0], D_mesons) && isPartOf(list_decay_mu[1], B_mesons)){
        std::cout << "-> BToD decay" << std::endl;
        mu_BtoDdecay = 1;
      }
      else if(list_decay_mu.size() == 3 && isPartOf(list_decay_mu[0], D_mesons) && isPartOf(list_decay_mu[1], K_mesons) && isPartOf(list_decay_mu[2], B_mesons)){
        std::cout << "-> BToDtoK decay" << std::endl;
        mu_BtoDtoKdecay = 1;
      }
      else if(list_decay_mu.size() == 2 && fabs(list_decay_mu[0]) == 15 && isPartOf(list_decay_mu[1], B_mesons)){
        std::cout << "-> BTotau decay" << std::endl;
        mu_Btotaudecay = 1;
      }
      else if(list_decay_mu.size() == 2 && fabs(list_decay_mu[0]) == 443 && isPartOf(list_decay_mu[1], B_mesons)){
        std::cout << "-> BToJpsi decay" << std::endl;
        mu_BtoJPsidecay = 1;
      }
      else if(list_decay_mu.size() > 2 && (isPartOf(list_decay_mu[1], excitedD_mesons) || isPartOf(list_decay_mu[list_decay_mu.size()-2], excitedD_mesons)) && isPartOf(list_decay_mu[list_decay_mu.size()-1], B_mesons)){
        std::cout << "-> BToDexcited decay" << std::endl;
        mu_BtoDexciteddecay = 1;
      }
      else if(!isPartOf(list_decay_mu[list_decay_mu.size()-1], B_mesons)){
        std::cout << "-> not a B decay" << std::endl;
        mu_nonBdecay = 1;
      }
      else{
        std::cout << "-> other" << std::endl;
        mu_other = 1;
      }
    
      std::cout << std::endl;
    }

    // case where the displaced muon is matched
    if(pi_ismatched){

      // get the ancestor chain
      std::cout << "-- Displaced pion --" << std::endl;
      std::vector<int> initialisator_pi;
      list_decay_pi = getAncestor(pi_genmatched_idx, initialisator_pi);

      pi_directBdecay = 0;
      pi_BtoDdecay = 0;
      pi_BtoDtoKdecay = 0;
      pi_Btotaudecay = 0;
      pi_BtoJPsidecay = 0;
      pi_BtoDexciteddecay = 0;
      pi_nonBdecay = 0;
      pi_other = 0;

      if(list_decay_pi.size() == 0){
        // remove pathological events
        pi_ismatched = 0;
      }
      else if(list_decay_pi.size() == 1 && isPartOf(list_decay_pi[0], B_mesons)){
        std::cout << "-> direct B decay" << std::endl;
        pi_directBdecay = 1;
      }
      else if(list_decay_pi.size() == 2 && isPartOf(list_decay_pi[0], D_mesons) && isPartOf(list_decay_pi[1], B_mesons)){
        std::cout << "-> BToD decay" << std::endl;
        pi_BtoDdecay = 1;
      }
      else if(list_decay_pi.size() == 3 && isPartOf(list_decay_pi[0], D_mesons) && isPartOf(list_decay_pi[1], K_mesons) && isPartOf(list_decay_pi[2], B_mesons)){
        std::cout << "-> BToDtoK decay" << std::endl;
        pi_BtoDtoKdecay = 1;
      }
      else if(list_decay_pi.size() == 2 && fabs(list_decay_pi[0]) == 15 && isPartOf(list_decay_pi[1], B_mesons)){
        std::cout << "-> BTotau decay" << std::endl;
        pi_Btotaudecay = 1;
      }
      else if(list_decay_pi.size() == 2 && fabs(list_decay_pi[0]) == 443 && isPartOf(list_decay_pi[1], B_mesons)){
        std::cout << "-> BToJpsi decay" << std::endl;
        pi_BtoJPsidecay = 1;
      }
      else if(list_decay_pi.size() > 2 && (isPartOf(list_decay_pi[1], excitedD_mesons) || isPartOf(list_decay_pi[list_decay_pi.size()-2], excitedD_mesons)) && isPartOf(list_decay_pi[list_decay_pi.size()-1], B_mesons)){
        std::cout << "-> BToDexcited decay" << std::endl;
        pi_BtoDexciteddecay = 1;
      }
      else if(!isPartOf(list_decay_pi[list_decay_pi.size()-1], B_mesons)){
        std::cout << "-> not a B decay" << std::endl;
        pi_nonBdecay = 1;
      }
      else{
        std::cout << "-> other" << std::endl;
        pi_other = 1;
      }
    
      std::cout << std::endl;
    }

    cand_isnotmatched = 0;
    cand_isfullmatched = 0;
    cand_ispartialmatched = 0;

    if(mu0_ismatched + mu_ismatched + pi_ismatched == 0){
      cand_isnotmatched = 1;
    }
    else if(mu0_ismatched + mu_ismatched + pi_ismatched == 3){
      cand_isfullmatched = 1;
    }
    else{
      cand_ispartialmatched = 1;
    }

    // match the mumupi candidate
    std::cout << "-- mumupi candidate --" << std::endl;

    if(mu0_ismatched && mu_ismatched){
      if(checkCommonAncestor(mu0_genmatched_idx, mu_genmatched_idx)){
        std::cout << "mu0 mu originating from same B meson" << std::endl;
        cand_mu0mu_samemother = 1;
      }
      else{
        cand_mu0mu_samemother = 0;
      }
    }
    
    if(mu0_ismatched && pi_ismatched){
      if(checkCommonAncestor(mu0_genmatched_idx, pi_genmatched_idx)){
        std::cout << "mu0 pi originating from same B meson" << std::endl;
        cand_mu0pi_samemother = 1;
      }
      else{
        cand_mu0pi_samemother = 0;
      }
    }

    if(mu_ismatched && pi_ismatched){
      if(checkCommonAncestor(mu_genmatched_idx, pi_genmatched_idx)){
        std::cout << "mu pi originating from same B meson" << std::endl;
        cand_mupi_samemother = 1;
      }
      else{
        cand_mupi_samemother = 0;
      }

      if(checkCommonAncestor(mu0_genmatched_idx, mu_genmatched_idx, pi_genmatched_idx)){
        std::cout << "mumupi originating from same B meson" << std::endl;
        cand_mumupi_samemother = 1;
      }
      else{
        cand_mumupi_samemother = 0;
      }
    }
  
    sources->Fill();

  }// end at least one candidate in the event
      

  return kTRUE;
}


void BackgroundSources::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}


void BackgroundSources::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  my_file->cd();

  sources->Write("", TObject::kOverwrite);

  my_file->Close();

  cout << "- End Background Sources -" << endl;
}

