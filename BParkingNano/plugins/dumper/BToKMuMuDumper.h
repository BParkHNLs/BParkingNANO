//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 23 15:15:57 2021 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: bparknano_nj89.root
//////////////////////////////////////////////////////////

#ifndef BToKMuMuDumper_h
#define BToKMuMuDumper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class BToKMuMuDumper : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> nBToKMuMu = {fReader, "nBToKMuMu"};
   TTreeReaderArray<Float_t> BToKMuMu_b_iso03 = {fReader, "BToKMuMu_b_iso03"};
   TTreeReaderArray<Float_t> BToKMuMu_b_iso03_close = {fReader, "BToKMuMu_b_iso03_close"};
   TTreeReaderArray<Float_t> BToKMuMu_b_iso04 = {fReader, "BToKMuMu_b_iso04"};
   TTreeReaderArray<Float_t> BToKMuMu_b_iso04_close = {fReader, "BToKMuMu_b_iso04_close"};
   TTreeReaderArray<Float_t> BToKMuMu_cos2D = {fReader, "BToKMuMu_cos2D"};
   TTreeReaderArray<Float_t> BToKMuMu_eta = {fReader, "BToKMuMu_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_cos2D = {fReader, "BToKMuMu_fit_cos2D"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_eta = {fReader, "BToKMuMu_fit_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_k_eta = {fReader, "BToKMuMu_fit_k_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_k_phi = {fReader, "BToKMuMu_fit_k_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_k_pt = {fReader, "BToKMuMu_fit_k_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l1_eta = {fReader, "BToKMuMu_fit_l1_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l1_phi = {fReader, "BToKMuMu_fit_l1_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l1_pt = {fReader, "BToKMuMu_fit_l1_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l2_eta = {fReader, "BToKMuMu_fit_l2_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l2_phi = {fReader, "BToKMuMu_fit_l2_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_l2_pt = {fReader, "BToKMuMu_fit_l2_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_mass = {fReader, "BToKMuMu_fit_mass"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_massErr = {fReader, "BToKMuMu_fit_massErr"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_phi = {fReader, "BToKMuMu_fit_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_fit_pt = {fReader, "BToKMuMu_fit_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_k_iso03 = {fReader, "BToKMuMu_k_iso03"};
   TTreeReaderArray<Float_t> BToKMuMu_k_iso03_close = {fReader, "BToKMuMu_k_iso03_close"};
   TTreeReaderArray<Float_t> BToKMuMu_k_iso04 = {fReader, "BToKMuMu_k_iso04"};
   TTreeReaderArray<Float_t> BToKMuMu_k_iso04_close = {fReader, "BToKMuMu_k_iso04_close"};
   TTreeReaderArray<Float_t> BToKMuMu_l1_iso03 = {fReader, "BToKMuMu_l1_iso03"};
   TTreeReaderArray<Float_t> BToKMuMu_l1_iso03_close = {fReader, "BToKMuMu_l1_iso03_close"};
   TTreeReaderArray<Float_t> BToKMuMu_l1_iso04 = {fReader, "BToKMuMu_l1_iso04"};
   TTreeReaderArray<Float_t> BToKMuMu_l1_iso04_close = {fReader, "BToKMuMu_l1_iso04_close"};
   TTreeReaderArray<Float_t> BToKMuMu_l2_iso03 = {fReader, "BToKMuMu_l2_iso03"};
   TTreeReaderArray<Float_t> BToKMuMu_l2_iso03_close = {fReader, "BToKMuMu_l2_iso03_close"};
   TTreeReaderArray<Float_t> BToKMuMu_l2_iso04 = {fReader, "BToKMuMu_l2_iso04"};
   TTreeReaderArray<Float_t> BToKMuMu_l2_iso04_close = {fReader, "BToKMuMu_l2_iso04_close"};
   TTreeReaderArray<Float_t> BToKMuMu_l_xy = {fReader, "BToKMuMu_l_xy"};
   TTreeReaderArray<Float_t> BToKMuMu_l_xy_unc = {fReader, "BToKMuMu_l_xy_unc"};
   TTreeReaderArray<Float_t> BToKMuMu_mass = {fReader, "BToKMuMu_mass"};
   TTreeReaderArray<Float_t> BToKMuMu_maxDR = {fReader, "BToKMuMu_maxDR"};
   TTreeReaderArray<Float_t> BToKMuMu_minDR = {fReader, "BToKMuMu_minDR"};
   TTreeReaderArray<Float_t> BToKMuMu_mllErr_llfit = {fReader, "BToKMuMu_mllErr_llfit"};
   TTreeReaderArray<Float_t> BToKMuMu_mll_fullfit = {fReader, "BToKMuMu_mll_fullfit"};
   TTreeReaderArray<Float_t> BToKMuMu_ll_sv_prob = {fReader, "BToKMuMu_ll_sv_prob"};
   TTreeReaderArray<Float_t> BToKMuMu_mll_llfit = {fReader, "BToKMuMu_mll_llfit"};
   TTreeReaderArray<Float_t> BToKMuMu_mll_raw = {fReader, "BToKMuMu_mll_raw"};
   TTreeReaderArray<Float_t> BToKMuMu_phi = {fReader, "BToKMuMu_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_pt = {fReader, "BToKMuMu_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_svprob = {fReader, "BToKMuMu_svprob"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_ex = {fReader, "BToKMuMu_vtx_ex"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_ey = {fReader, "BToKMuMu_vtx_ey"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_ez = {fReader, "BToKMuMu_vtx_ez"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_x = {fReader, "BToKMuMu_vtx_x"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_y = {fReader, "BToKMuMu_vtx_y"};
   TTreeReaderArray<Float_t> BToKMuMu_vtx_z = {fReader, "BToKMuMu_vtx_z"};
   TTreeReaderArray<Int_t> BToKMuMu_charge = {fReader, "BToKMuMu_charge"};
   TTreeReaderArray<Int_t> BToKMuMu_isMatched = {fReader, "BToKMuMu_isMatched"};
   TTreeReaderArray<Int_t> BToKMuMu_kIdx = {fReader, "BToKMuMu_kIdx"};
   TTreeReaderArray<Int_t> BToKMuMu_l1Idx = {fReader, "BToKMuMu_l1Idx"};
   TTreeReaderArray<Int_t> BToKMuMu_l2Idx = {fReader, "BToKMuMu_l2Idx"};
   TTreeReaderArray<Float_t> BToKMuMu_k_DCASig = {fReader, "BToKMuMu_k_DCASig"};
   TTreeReaderArray<Float_t> BToKMuMu_k_dxy = {fReader, "BToKMuMu_k_dxy"};
   TTreeReaderArray<Float_t> BToKMuMu_k_dxyS = {fReader, "BToKMuMu_k_dxyS"};
   TTreeReaderArray<Float_t> BToKMuMu_k_dz = {fReader, "BToKMuMu_k_dz"};
   TTreeReaderArray<Float_t> BToKMuMu_k_dzS = {fReader, "BToKMuMu_k_dzS"};
   TTreeReaderArray<Float_t> BToKMuMu_k_eta = {fReader, "BToKMuMu_k_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_k_mass = {fReader, "BToKMuMu_k_mass"};
   TTreeReaderArray<Int_t> BToKMuMu_k_charge = {fReader, "BToKMuMu_k_charge"};
   TTreeReaderArray<Float_t> BToKMuMu_k_phi = {fReader, "BToKMuMu_k_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_k_pt = {fReader, "BToKMuMu_k_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_k_vx = {fReader, "BToKMuMu_k_vx"};
   TTreeReaderArray<Float_t> BToKMuMu_k_vy = {fReader, "BToKMuMu_k_vy"};
   TTreeReaderArray<Float_t> BToKMuMu_k_vz = {fReader, "BToKMuMu_k_vz"};
   TTreeReaderArray<Int_t> BToKMuMu_k_highPurityFlag = {fReader, "BToKMuMu_k_highPurityFlag"};
   TTreeReaderArray<Int_t> BToKMuMu_k_ispacked = {fReader, "BToKMuMu_k_ispacked"};
   TTreeReaderArray<Int_t> BToKMuMu_k_islost = {fReader, "BToKMuMu_k_islost"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_k_genIdx = {fReader, "BToKMuMu_matching_k_genIdx"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_k_motherPdgId = {fReader, "BToKMuMu_matching_k_motherPdgId"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_l1_genIdx = {fReader, "BToKMuMu_matching_l1_genIdx"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_l1_motherPdgId = {fReader, "BToKMuMu_matching_l1_motherPdgId"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_l2_genIdx = {fReader, "BToKMuMu_matching_l2_genIdx"};
   TTreeReaderArray<Int_t> BToKMuMu_matching_l2_motherPdgId = {fReader, "BToKMuMu_matching_l2_motherPdgId"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_b_pt       = {fReader, "BToKMuMu_matched_b_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_b_eta      = {fReader, "BToKMuMu_matched_b_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_b_phi      = {fReader, "BToKMuMu_matched_b_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_b_mass     = {fReader, "BToKMuMu_matched_b_mass"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_l1_pt      = {fReader, "BToKMuMu_matched_l1_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_l1_eta     = {fReader, "BToKMuMu_matched_l1_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_l1_phi     = {fReader, "BToKMuMu_matched_l1_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_l1_mass    = {fReader, "BToKMuMu_matched_l1_mass"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_l2_pt      = {fReader, "BToKMuMu_matched_l2_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_l2_eta     = {fReader, "BToKMuMu_matched_l2_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_l2_phi     = {fReader, "BToKMuMu_matched_l2_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_l2_mass    = {fReader, "BToKMuMu_matched_l2_mass"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_k_pt       = {fReader, "BToKMuMu_matched_k_pt"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_k_eta      = {fReader, "BToKMuMu_matched_k_eta"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_k_phi      = {fReader, "BToKMuMu_matched_k_phi"};
   TTreeReaderArray<Float_t> BToKMuMu_matched_k_mass     = {fReader, "BToKMuMu_matched_k_mass"};
   TTreeReaderArray<Int_t> BToKMuMu_n_k_used = {fReader, "BToKMuMu_n_k_used"};
   TTreeReaderArray<Int_t> BToKMuMu_n_l1_used = {fReader, "BToKMuMu_n_l1_used"};
   TTreeReaderArray<Int_t> BToKMuMu_n_l2_used = {fReader, "BToKMuMu_n_l2_used"};
   TTreeReaderArray<Int_t> BToKMuMu_pdgId = {fReader, "BToKMuMu_pdgId"};
   TTreeReaderValue<UInt_t> nMuon = {fReader, "nMuon"};
   TTreeReaderArray<Int_t> Muon_isSlimmedMuon = {fReader, "Muon_isSlimmedMuon"};
   TTreeReaderArray<Int_t> Muon_isDSAMuon = {fReader, "Muon_isDSAMuon"};
   TTreeReaderArray<Float_t> Muon_caloCompatibility = {fReader, "Muon_caloCompatibility"};
   TTreeReaderArray<Float_t> Muon_dxy = {fReader, "Muon_dxy"};
   TTreeReaderArray<Float_t> Muon_dxyS = {fReader, "Muon_dxyS"};
   TTreeReaderArray<Float_t> Muon_dz = {fReader, "Muon_dz"};
   TTreeReaderArray<Float_t> Muon_dzS = {fReader, "Muon_dzS"};
   TTreeReaderArray<Float_t> Muon_eta = {fReader, "Muon_eta"};
   TTreeReaderArray<Float_t> Muon_globalNormalisedChi2 = {fReader, "Muon_globalNormalisedChi2"};
   //TTreeReaderArray<Float_t> Muon_ip3d = {fReader, "Muon_ip3d"};
   TTreeReaderArray<Float_t> Muon_kinkFinderChi2 = {fReader, "Muon_kinkFinderChi2"};
   TTreeReaderArray<Float_t> Muon_localPositionChi2 = {fReader, "Muon_localPositionChi2"};
   TTreeReaderArray<Float_t> Muon_mass = {fReader, "Muon_mass"};
   TTreeReaderArray<Float_t> Muon_matched_dpt = {fReader, "Muon_matched_dpt"};
   TTreeReaderArray<Float_t> Muon_matched_dr = {fReader, "Muon_matched_dr"};
   TTreeReaderArray<Float_t> Muon_pfiso03Rel_all = {fReader, "Muon_pfiso03Rel_all"};
   //TTreeReaderArray<Float_t> Muon_pfiso03Rel_ch = {fReader, "Muon_pfiso03Rel_ch"};
   //TTreeReaderArray<Float_t> Muon_pfiso03Rel_n = {fReader, "Muon_pfiso03Rel_n"};
   //TTreeReaderArray<Float_t> Muon_pfiso03Rel_pho = {fReader, "Muon_pfiso03Rel_pho"};
   //TTreeReaderArray<Float_t> Muon_pfiso03Rel_pu = {fReader, "Muon_pfiso03Rel_pu"};
   //TTreeReaderArray<Float_t> Muon_pfiso03Rel_trk = {fReader, "Muon_pfiso03Rel_trk"};
   TTreeReaderArray<Float_t> Muon_pfiso03_all = {fReader, "Muon_pfiso03_all"};
   //TTreeReaderArray<Float_t> Muon_pfiso03_ch = {fReader, "Muon_pfiso03_ch"};
   //TTreeReaderArray<Float_t> Muon_pfiso03_n = {fReader, "Muon_pfiso03_n"};
   //TTreeReaderArray<Float_t> Muon_pfiso03_pho = {fReader, "Muon_pfiso03_pho"};
   //TTreeReaderArray<Float_t> Muon_pfiso03_pu = {fReader, "Muon_pfiso03_pu"};
   //TTreeReaderArray<Float_t> Muon_pfiso03_trk = {fReader, "Muon_pfiso03_trk"};
   TTreeReaderArray<Float_t> Muon_pfiso04Rel_all = {fReader, "Muon_pfiso04Rel_all"};
   //TTreeReaderArray<Float_t> Muon_pfiso04Rel_ch = {fReader, "Muon_pfiso04Rel_ch"};
   //TTreeReaderArray<Float_t> Muon_pfiso04Rel_n = {fReader, "Muon_pfiso04Rel_n"};
   //TTreeReaderArray<Float_t> Muon_pfiso04Rel_pho = {fReader, "Muon_pfiso04Rel_pho"};
   //TTreeReaderArray<Float_t> Muon_pfiso04Rel_pu = {fReader, "Muon_pfiso04Rel_pu"};
   TTreeReaderArray<Float_t> Muon_pfiso04_all = {fReader, "Muon_pfiso04_all"};
   //TTreeReaderArray<Float_t> Muon_pfiso04_ch = {fReader, "Muon_pfiso04_ch"};
   //TTreeReaderArray<Float_t> Muon_pfiso04_n = {fReader, "Muon_pfiso04_n"};
   //TTreeReaderArray<Float_t> Muon_pfiso04_pho = {fReader, "Muon_pfiso04_pho"};
   //TTreeReaderArray<Float_t> Muon_pfiso04_pu = {fReader, "Muon_pfiso04_pu"};
   //TTreeReaderArray<Float_t> Muon_pfiso05Rel_trk = {fReader, "Muon_pfiso05Rel_trk"};
   //TTreeReaderArray<Float_t> Muon_pfiso05_trk = {fReader, "Muon_pfiso05_trk"};
   TTreeReaderArray<Float_t> Muon_phi = {fReader, "Muon_phi"};
   TTreeReaderArray<Float_t> Muon_pt = {fReader, "Muon_pt"};
   //TTreeReaderArray<Float_t> Muon_ptErr = {fReader, "Muon_ptErr"};
   TTreeReaderArray<Float_t> Muon_segmentCompatibility = {fReader, "Muon_segmentCompatibility"};
   //TTreeReaderArray<Float_t> Muon_sip3d = {fReader, "Muon_sip3d"};
   TTreeReaderArray<Float_t> Muon_validHitFraction = {fReader, "Muon_validHitFraction"};
   TTreeReaderArray<Float_t> Muon_vx = {fReader, "Muon_vx"};
   TTreeReaderArray<Float_t> Muon_vy = {fReader, "Muon_vy"};
   TTreeReaderArray<Float_t> Muon_vz = {fReader, "Muon_vz"};
   TTreeReaderArray<Int_t> Muon_charge = {fReader, "Muon_charge"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu10p5_IP3p5 = {fReader, "Muon_fired_HLT_Mu10p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu12_IP6 = {fReader, "Muon_fired_HLT_Mu12_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu7_IP4 = {fReader, "Muon_fired_HLT_Mu7_IP4"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP3 = {fReader, "Muon_fired_HLT_Mu8_IP3"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP5 = {fReader, "Muon_fired_HLT_Mu8_IP5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP6 = {fReader, "Muon_fired_HLT_Mu8_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8p5_IP3p5 = {fReader, "Muon_fired_HLT_Mu8p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP4 = {fReader, "Muon_fired_HLT_Mu9_IP4"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP5 = {fReader, "Muon_fired_HLT_Mu9_IP5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP6 = {fReader, "Muon_fired_HLT_Mu9_IP6"};
   TTreeReaderArray<Int_t> Muon_isTriggering = {fReader, "Muon_isTriggering"};
   TTreeReaderArray<Int_t> Muon_isTriggeringBPark = {fReader, "Muon_isTriggeringBPark"};
   TTreeReaderArray<Int_t> Muon_looseId = {fReader, "Muon_looseId"};
   TTreeReaderArray<Int_t> Muon_numberOfPixelLayers = {fReader, "Muon_numberOfPixelLayers"};
   TTreeReaderArray<Int_t> Muon_numberOfStations = {fReader, "Muon_numberOfStations"};
   TTreeReaderArray<Int_t> Muon_numberOfTrackerLayers = {fReader, "Muon_numberOfTrackerLayers"};
   TTreeReaderArray<Int_t> Muon_numberOfValidMuonHits = {fReader, "Muon_numberOfValidMuonHits"};
   TTreeReaderArray<Int_t> Muon_numberOfValidPixelHits = {fReader, "Muon_numberOfValidPixelHits"};
   TTreeReaderArray<Int_t> Muon_pdgId = {fReader, "Muon_pdgId"};
   TTreeReaderArray<Int_t> Muon_trackerHighPurityFlag = {fReader, "Muon_trackerHighPurityFlag"};
   TTreeReaderArray<Int_t> Muon_inTimeMuon = {fReader, "Muon_inTimeMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalMuon = {fReader, "Muon_isGlobalMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalNotTrackerMuon = {fReader, "Muon_isGlobalNotTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalOrTrackerMuon = {fReader, "Muon_isGlobalOrTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isPF = {fReader, "Muon_isPF"};
   TTreeReaderArray<Int_t> Muon_isTrackerMuon = {fReader, "Muon_isTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isTrackerNotGlobalMuon = {fReader, "Muon_isTrackerNotGlobalMuon"};
   TTreeReaderArray<Int_t> Muon_mediumId = {fReader, "Muon_mediumId"};
   TTreeReaderArray<UChar_t> Muon_pfIsoId = {fReader, "Muon_pfIsoId"};
   TTreeReaderArray<Int_t> Muon_softId = {fReader, "Muon_softId"};
   TTreeReaderArray<Int_t> Muon_tightId = {fReader, "Muon_tightId"};
   TTreeReaderArray<UChar_t> Muon_tkIsoId = {fReader, "Muon_tkIsoId"};
   TTreeReaderArray<Int_t> Muon_triggerIdLoose = {fReader, "Muon_triggerIdLoose"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetAll = {fReader, "fixedGridRhoFastjetAll"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentral = {fReader, "fixedGridRhoFastjetCentral"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralCalo = {fReader, "fixedGridRhoFastjetCentralCalo"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralChargedPileUp = {fReader, "fixedGridRhoFastjetCentralChargedPileUp"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralNeutral = {fReader, "fixedGridRhoFastjetCentralNeutral"};
   //TTreeReaderValue<UInt_t> nProbeTracks = {fReader, "nProbeTracks"};
   //TTreeReaderArray<Float_t> ProbeTracks_DCASig = {fReader, "ProbeTracks_DCASig"};
   //TTreeReaderArray<Float_t> ProbeTracks_dxy = {fReader, "ProbeTracks_dxy"};
   //TTreeReaderArray<Float_t> ProbeTracks_dxyS = {fReader, "ProbeTracks_dxyS"};
   //TTreeReaderArray<Float_t> ProbeTracks_dz = {fReader, "ProbeTracks_dz"};
   //TTreeReaderArray<Float_t> ProbeTracks_dzS = {fReader, "ProbeTracks_dzS"};
   //TTreeReaderArray<Float_t> ProbeTracks_eta = {fReader, "ProbeTracks_eta"};
   //TTreeReaderArray<Float_t> ProbeTracks_mass = {fReader, "ProbeTracks_mass"};
   //TTreeReaderArray<Float_t> ProbeTracks_phi = {fReader, "ProbeTracks_phi"};
   //TTreeReaderArray<Float_t> ProbeTracks_pt = {fReader, "ProbeTracks_pt"};
   //TTreeReaderArray<Float_t> ProbeTracks_vx = {fReader, "ProbeTracks_vx"};
   //TTreeReaderArray<Float_t> ProbeTracks_vy = {fReader, "ProbeTracks_vy"};
   //TTreeReaderArray<Float_t> ProbeTracks_vz = {fReader, "ProbeTracks_vz"};
   //TTreeReaderArray<Int_t> ProbeTracks_charge = {fReader, "ProbeTracks_charge"};
   //TTreeReaderArray<Int_t> ProbeTracks_isLostTrk = {fReader, "ProbeTracks_isLostTrk"};
   //TTreeReaderArray<Int_t> ProbeTracks_isPacked = {fReader, "ProbeTracks_isPacked"};
   //TTreeReaderArray<Int_t> ProbeTracks_nValidHits = {fReader, "ProbeTracks_nValidHits"};
   //TTreeReaderArray<Int_t> ProbeTracks_pdgId = {fReader, "ProbeTracks_pdgId"};
   //TTreeReaderArray<Bool_t> ProbeTracks_isMatchedToEle = {fReader, "ProbeTracks_isMatchedToEle"};
   //TTreeReaderArray<Bool_t> ProbeTracks_isMatchedToLooseMuon = {fReader, "ProbeTracks_isMatchedToLooseMuon"};
   //TTreeReaderArray<Bool_t> ProbeTracks_isMatchedToMediumMuon = {fReader, "ProbeTracks_isMatchedToMediumMuon"};
   //TTreeReaderArray<Bool_t> ProbeTracks_isMatchedToMuon = {fReader, "ProbeTracks_isMatchedToMuon"};
   //TTreeReaderArray<Bool_t> ProbeTracks_isMatchedToSoftMuon = {fReader, "ProbeTracks_isMatchedToSoftMuon"};
   TTreeReaderValue<UChar_t> HLT_Mu7_IP4 = {fReader, "HLT_Mu7_IP4"};
   TTreeReaderValue<UChar_t> HLT_Mu8_IP6 = {fReader, "HLT_Mu8_IP6"};
   TTreeReaderValue<UChar_t> HLT_Mu8_IP5 = {fReader, "HLT_Mu8_IP5"};
   TTreeReaderValue<UChar_t> HLT_Mu8_IP3 = {fReader, "HLT_Mu8_IP3"};
   TTreeReaderValue<UChar_t> HLT_Mu8p5_IP3p5 = {fReader, "HLT_Mu8p5_IP3p5"};
   TTreeReaderValue<UChar_t> HLT_Mu9_IP6 = {fReader, "HLT_Mu9_IP6"};
   TTreeReaderValue<UChar_t> HLT_Mu9_IP5 = {fReader, "HLT_Mu9_IP5"};
   TTreeReaderValue<UChar_t> HLT_Mu9_IP4 = {fReader, "HLT_Mu9_IP4"};
   TTreeReaderValue<UChar_t> HLT_Mu10p5_IP3p5 = {fReader, "HLT_Mu10p5_IP3p5"};
   TTreeReaderValue<UChar_t> HLT_Mu12_IP6 = {fReader, "HLT_Mu12_IP6"};
   TTreeReaderValue<UChar_t> L1_SingleMu7er1p5 = {fReader, "L1_SingleMu7er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu8er1p5 = {fReader, "L1_SingleMu8er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu9er1p5 = {fReader, "L1_SingleMu9er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu10er1p5 = {fReader, "L1_SingleMu10er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu12er1p5 = {fReader, "L1_SingleMu12er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu22 = {fReader, "L1_SingleMu22"};
   TTreeReaderValue<UInt_t> nTrigObj = {fReader, "nTrigObj"};
   TTreeReaderArray<Float_t> TrigObj_pt = {fReader, "TrigObj_pt"};
   TTreeReaderArray<Float_t> TrigObj_eta = {fReader, "TrigObj_eta"};
   TTreeReaderArray<Float_t> TrigObj_phi = {fReader, "TrigObj_phi"};
   TTreeReaderArray<Float_t> TrigObj_l1pt = {fReader, "TrigObj_l1pt"};
   TTreeReaderArray<Float_t> TrigObj_l1pt_2 = {fReader, "TrigObj_l1pt_2"};
   TTreeReaderArray<Float_t> TrigObj_l2pt = {fReader, "TrigObj_l2pt"};
   TTreeReaderArray<Int_t> TrigObj_id = {fReader, "TrigObj_id"};
   TTreeReaderArray<Int_t> TrigObj_l1iso = {fReader, "TrigObj_l1iso"};
   TTreeReaderArray<Int_t> TrigObj_l1charge = {fReader, "TrigObj_l1charge"};
   TTreeReaderArray<Int_t> TrigObj_filterBits = {fReader, "TrigObj_filterBits"};
   TTreeReaderValue<UInt_t> nOtherPV = {fReader, "nOtherPV"};
   TTreeReaderArray<Float_t> OtherPV_z = {fReader, "OtherPV_z"};
   TTreeReaderValue<Float_t> PV_ndof = {fReader, "PV_ndof"};
   TTreeReaderValue<Float_t> PV_x = {fReader, "PV_x"};
   TTreeReaderValue<Float_t> PV_y = {fReader, "PV_y"};
   TTreeReaderValue<Float_t> PV_z = {fReader, "PV_z"};
   TTreeReaderValue<Float_t> PV_chi2 = {fReader, "PV_chi2"};
   TTreeReaderValue<Float_t> PV_score = {fReader, "PV_score"};
   TTreeReaderValue<Int_t> PV_npvs = {fReader, "PV_npvs"};
   TTreeReaderValue<Int_t> PV_npvsGood = {fReader, "PV_npvsGood"};
   TTreeReaderValue<UInt_t> nSV = {fReader, "nSV"};
   TTreeReaderArray<Float_t> SV_dlen = {fReader, "SV_dlen"};
   TTreeReaderArray<Float_t> SV_dlenSig = {fReader, "SV_dlenSig"};
   TTreeReaderArray<Float_t> SV_pAngle = {fReader, "SV_pAngle"};
   TTreeReaderArray<Float_t> SV_chi2 = {fReader, "SV_chi2"};
   TTreeReaderArray<Float_t> SV_eta = {fReader, "SV_eta"};
   TTreeReaderArray<Float_t> SV_mass = {fReader, "SV_mass"};
   TTreeReaderArray<Float_t> SV_ndof = {fReader, "SV_ndof"};
   TTreeReaderArray<Float_t> SV_phi = {fReader,"SV_phi"};
   TTreeReaderArray<Float_t> SV_pt = {fReader, "SV_pt"};
   TTreeReaderArray<Float_t> SV_x = {fReader, "SV_x"};
   TTreeReaderArray<Float_t> SV_y = {fReader, "SV_y"};
   TTreeReaderArray<Float_t> SV_z = {fReader, "SV_z"};
   TTreeReaderValue<Float_t> Pileup_nTrueInt = {fReader, "PV_ndof"};

   BToKMuMuDumper(TTree * /*tree*/ =0) { }
   virtual ~BToKMuMuDumper() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   // output file
   TFile* my_file;  

   // some branches will be added only if sample is MC
   Bool_t isMC;

   // some option
   // this option is intentionally hardcoded
   Bool_t do_fillhistograms = false;

   // trees to fill
   TTree* control_tree;
   
   // filling variables
   ULong64_t the_event= -99;
   Int_t the_run = -99;
   Int_t the_lumi = -99;

   Int_t the_pv_npvs = -99;

   Int_t the_hlt_mu7_ip4 = -99;
   Int_t the_hlt_mu8_ip6 = -99;
   Int_t the_hlt_mu8_ip5 = -99;
   Int_t the_hlt_mu8_ip3 = -99;
   Int_t the_hlt_mu8p5_ip3p5 = -99;
   Int_t the_hlt_mu9_ip6 = -99;
   Int_t the_hlt_mu9_ip5 = -99;
   Int_t the_hlt_mu9_ip4 = -99;
   Int_t the_hlt_mu10p5_ip3p5 = -99;
   Int_t the_hlt_mu12_ip6 = -99;

   Float_t the_ctrl_b_pt = -99.;
   Float_t the_ctrl_b_eta = -99.;
   Float_t the_ctrl_b_y = -99.;
   Float_t the_ctrl_b_phi = -99.;
   Float_t the_ctrl_b_mass = -99.;
   Int_t the_ctrl_b_charge = -99;
   Int_t the_ctrl_b_pdgid = -99;
   Float_t the_ctrl_b_cos2d = -99.;
   Float_t the_ctrl_b_iso03 = -99.;
   Float_t the_ctrl_b_iso03_close = -99.;
   Float_t the_ctrl_b_iso04 = -99.;
   Float_t the_ctrl_b_iso04_close = -99.;

   Float_t the_ctrl_k_pt = -99.;
   Float_t the_ctrl_k_eta = -99.;
   Float_t the_ctrl_k_phi = -99.;
   Int_t the_ctrl_k_charge = -99;
   Float_t the_ctrl_k_iso03 = -99.;
   Float_t the_ctrl_k_iso03_close = -99.;
   Float_t the_ctrl_k_iso04 = -99.;
   Float_t the_ctrl_k_iso04_close = -99.;

   Float_t the_ctrl_l1_pt = -99.;
   Float_t the_ctrl_l1_eta = -99.;
   Float_t the_ctrl_l1_phi = -99.;
   Int_t the_ctrl_l1_charge = -99;
   Float_t the_ctrl_l1_dxy = -99.;
   Float_t the_ctrl_l1_dxysig = -99.;
   Float_t the_ctrl_l1_dz = -99.;
   Float_t the_ctrl_l1_dzsig = -99.;
   Float_t the_ctrl_l1_iso03= -99.;
   Float_t the_ctrl_l1_iso03_close= -99.;
   Float_t the_ctrl_l1_iso04= -99.;
   Float_t the_ctrl_l1_iso04_close= -99.;
   Int_t the_ctrl_l1_looseid = -99;
   Int_t the_ctrl_l1_mediumid = -99;
   Int_t the_ctrl_l1_tightid = -99;
   Int_t the_ctrl_l1_softid = -99;
   Int_t the_ctrl_l1_pfisoid = -99;
   Int_t the_ctrl_l1_trkisoid = -99;
   Int_t the_ctrl_l1_triggerlooseid = -99;
   Float_t the_ctrl_l1_istriggering = -99.;
   Int_t the_ctrl_l1_fired_hlt_mu7_ip4 = -99;
   Int_t the_ctrl_l1_fired_hlt_mu8_ip6 = -99;
   Int_t the_ctrl_l1_fired_hlt_mu8_ip5 = -99;
   Int_t the_ctrl_l1_fired_hlt_mu8_ip3 = -99;
   Int_t the_ctrl_l1_fired_hlt_mu8p5_ip3p5 = -99;
   Int_t the_ctrl_l1_fired_hlt_mu9_ip6 = -99;
   Int_t the_ctrl_l1_fired_hlt_mu9_ip5 = -99;
   Int_t the_ctrl_l1_fired_hlt_mu9_ip4 = -99;
   Int_t the_ctrl_l1_fired_hlt_mu10p5_ip3p5 = -99;
   Int_t the_ctrl_l1_fired_hlt_mu12_ip6 = -99;

   Float_t the_ctrl_l2_pt = -99.;
   Float_t the_ctrl_l2_eta = -99.;
   Float_t the_ctrl_l2_phi = -99.;
   Int_t the_ctrl_l2_charge = -99;
   Float_t the_ctrl_l2_dxy = -99.;
   Float_t the_ctrl_l2_dxysig = -99.;
   Float_t the_ctrl_l2_dz = -99.;
   Float_t the_ctrl_l2_dzsig = -99.;
   Float_t the_ctrl_l2_iso03 = -99.;
   Float_t the_ctrl_l2_iso03_close = -99.;
   Float_t the_ctrl_l2_iso04 = -99.;
   Float_t the_ctrl_l2_iso04_close = -99.;
   Int_t the_ctrl_l2_looseid = -99;
   Int_t the_ctrl_l2_mediumid = -99;
   Int_t the_ctrl_l2_tightid = -99;
   Int_t the_ctrl_l2_softid = -99;
   Int_t the_ctrl_l2_pfisoid = -99;
   Int_t the_ctrl_l2_trkisoid = -99;
   Int_t the_ctrl_l2_triggerlooseid = -99;
   Float_t the_ctrl_l2_istriggering = -99.;

   Float_t the_ctrl_dimu_mass = -99.;
   Float_t the_ctrl_dimu_sv_prob = -99.;

   Float_t the_ctrl_maxdr = -99.;
   Float_t the_ctrl_mindr = -99.;

   Float_t the_ctrl_sv_x = -99.;
   Float_t the_ctrl_sv_y = -99.;
   Float_t the_ctrl_sv_z = -99.;
   Float_t the_ctrl_sv_lxy = -99.;
   Float_t the_ctrl_sv_lxysig = -99.;
   Float_t the_ctrl_sv_prob = -99.;

   Int_t the_ctrl_ismatched = -99;
   Float_t the_ctrl_matched_b_pt = -99.;
   Float_t the_ctrl_matched_b_eta = -99.;
   Float_t the_ctrl_matched_b_y = -99.;

   Float_t the_ctrl_weight_hlt_A1 = -99.;
   Float_t the_ctrl_weight_hlt_A1_6 = -99.;

   Float_t the_ctrl_weight_pu_sig_A = -99.;
   Float_t the_ctrl_weight_pu_sig_B = -99.;
   Float_t the_ctrl_weight_pu_sig_C = -99.;
   Float_t the_ctrl_weight_pu_sig_D = -99.;
   Float_t the_ctrl_weight_pu_sig_tot = -99.;

   // histograms
   TH1F* ctrlhist_ncand_perevent;
   TH1F* ctrlhist_ncand_matched_perevent;
   TH1F* ctrlhist_selection_efficiency_bpt_allevents;
   TH1F* ctrlhist_selection_efficiency_bpt_eventswithmultcands;
   TH1F* ctrlhist_selection_efficiency_svprob_allevents;
   TH1F* ctrlhist_selection_efficiency_svprob_eventswithmultcands;
   TH1F* ctrlhist_selection_efficiency_cos2d_allevents;
   TH1F* ctrlhist_selection_efficiency_cos2d_eventswithmultcands;
   TH1F* ctrlhist_selection_efficiency_kpt_allevents;
   TH1F* ctrlhist_selection_efficiency_kpt_eventswithmultcands;
   TH1F* ctrlhist_ncand_wtriggeringmuon;

   ClassDef(BToKMuMuDumper,0);
};

#endif

#ifdef BToKMuMuDumper_cxx
void BToKMuMuDumper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t BToKMuMuDumper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef BToKMuMuDumper_cxx
