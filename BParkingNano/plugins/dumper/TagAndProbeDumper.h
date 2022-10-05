//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jul 24 21:08:29 2021 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: bparknano_tag_and_probe_v1.root
//////////////////////////////////////////////////////////

#ifndef TagAndProbeDumper_h
#define TagAndProbeDumper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class TagAndProbeDumper : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> nJPsiToMuMu = {fReader, "nJPsiToMuMu"};
   TTreeReaderArray<Float_t> JPsiToMuMu_charge = {fReader, "JPsiToMuMu_charge"};
   TTreeReaderArray<Float_t> JPsiToMuMu_cos2D = {fReader, "JPsiToMuMu_cos2D"};
   TTreeReaderArray<Float_t> JPsiToMuMu_deltaR = {fReader, "JPsiToMuMu_deltaR"};
   TTreeReaderArray<Float_t> JPsiToMuMu_eta = {fReader, "JPsiToMuMu_eta"};
   TTreeReaderArray<Float_t> JPsiToMuMu_lep1_eta = {fReader, "JPsiToMuMu_lep1_eta"};
   TTreeReaderArray<Float_t> JPsiToMuMu_lep1_phi = {fReader, "JPsiToMuMu_lep1_phi"};
   TTreeReaderArray<Float_t> JPsiToMuMu_lep1_pt = {fReader, "JPsiToMuMu_lep1_pt"};
   TTreeReaderArray<Float_t> JPsiToMuMu_lep2_eta = {fReader, "JPsiToMuMu_lep2_eta"};
   TTreeReaderArray<Float_t> JPsiToMuMu_lep2_phi = {fReader, "JPsiToMuMu_lep2_phi"};
   TTreeReaderArray<Float_t> JPsiToMuMu_lep2_pt = {fReader, "JPsiToMuMu_lep2_pt"};
   TTreeReaderArray<Float_t> JPsiToMuMu_lep_vzdiff = {fReader, "JPsiToMuMu_lep_vzdiff"};
   TTreeReaderArray<Float_t> JPsiToMuMu_lxy = {fReader, "JPsiToMuMu_lxy"};
   TTreeReaderArray<Float_t> JPsiToMuMu_lxy_sig = {fReader, "JPsiToMuMu_lxy_sig"};
   TTreeReaderArray<Float_t> JPsiToMuMu_mass = {fReader, "JPsiToMuMu_mass"};
   TTreeReaderArray<Float_t> JPsiToMuMu_phi = {fReader, "JPsiToMuMu_phi"};
   TTreeReaderArray<Float_t> JPsiToMuMu_pt = {fReader, "JPsiToMuMu_pt"};
   TTreeReaderArray<Float_t> JPsiToMuMu_sv_chi2 = {fReader, "JPsiToMuMu_sv_chi2"};
   TTreeReaderArray<Float_t> JPsiToMuMu_sv_ndof = {fReader, "JPsiToMuMu_sv_ndof"};
   TTreeReaderArray<Float_t> JPsiToMuMu_sv_prob = {fReader, "JPsiToMuMu_sv_prob"};
   TTreeReaderArray<Int_t> JPsiToMuMu_isMatched = {fReader, "JPsiToMuMu_isMatched"};
   TTreeReaderArray<Int_t> JPsiToMuMu_lep1_idx = {fReader, "JPsiToMuMu_lep1_idx"};
   TTreeReaderArray<Int_t> JPsiToMuMu_lep2_idx = {fReader, "JPsiToMuMu_lep2_idx"};
   TTreeReaderValue<UInt_t> nMuon = {fReader, "nMuon"};
   TTreeReaderArray<Float_t> Muon_caloCompatibility = {fReader, "Muon_caloCompatibility"};
   TTreeReaderArray<Float_t> Muon_dxy = {fReader, "Muon_dxy"};
   TTreeReaderArray<Float_t> Muon_dxy_BS = {fReader, "Muon_dxy_BS"};
   TTreeReaderArray<Float_t> Muon_dxyS = {fReader, "Muon_dxyS"};
   TTreeReaderArray<Float_t> Muon_dxyS_BS = {fReader, "Muon_dxyS_BS"};
   TTreeReaderArray<Float_t> Muon_dxyS_BS_alaRdst = {fReader, "Muon_dxyS_BS_alaRdst"};
   TTreeReaderArray<Float_t> Muon_dz = {fReader, "Muon_dz"};
   TTreeReaderArray<Float_t> Muon_dzS = {fReader, "Muon_dzS"};
   TTreeReaderArray<Float_t> Muon_eta = {fReader, "Muon_eta"};
   TTreeReaderArray<Float_t> Muon_globalNormalisedChi2 = {fReader, "Muon_globalNormalisedChi2"};
   TTreeReaderArray<Float_t> Muon_kinkFinderChi2 = {fReader, "Muon_kinkFinderChi2"};
   TTreeReaderArray<Float_t> Muon_localPositionChi2 = {fReader, "Muon_localPositionChi2"};
   TTreeReaderArray<Float_t> Muon_mass = {fReader, "Muon_mass"};
   TTreeReaderArray<Float_t> Muon_matched_dpt = {fReader, "Muon_matched_dpt"};
   TTreeReaderArray<Float_t> Muon_matched_dr = {fReader, "Muon_matched_dr"};
   TTreeReaderArray<Float_t> Muon_pfiso03Rel_all = {fReader, "Muon_pfiso03Rel_all"};
   TTreeReaderArray<Float_t> Muon_pfiso03_all = {fReader, "Muon_pfiso03_all"};
   TTreeReaderArray<Float_t> Muon_pfiso04Rel_all = {fReader, "Muon_pfiso04Rel_all"};
   TTreeReaderArray<Float_t> Muon_pfiso04_all = {fReader, "Muon_pfiso04_all"};
   TTreeReaderArray<Float_t> Muon_phi = {fReader, "Muon_phi"};
   TTreeReaderArray<Float_t> Muon_pt = {fReader, "Muon_pt"};
   TTreeReaderArray<Float_t> Muon_segmentCompatibility = {fReader, "Muon_segmentCompatibility"};
   TTreeReaderArray<Float_t> Muon_validHitFraction = {fReader, "Muon_validHitFraction"};
   TTreeReaderArray<Float_t> Muon_vx = {fReader, "Muon_vx"};
   TTreeReaderArray<Float_t> Muon_vy = {fReader, "Muon_vy"};
   TTreeReaderArray<Float_t> Muon_vz = {fReader, "Muon_vz"};
   TTreeReaderArray<Int_t> Muon_charge = {fReader, "Muon_charge"};
   TTreeReaderArray<Int_t> Muon_fired_DST_DoubleMu1_noVtx_CaloScouting_v2 = {fReader, "Muon_fired_DST_DoubleMu1_noVtx_CaloScouting_v2"};
   TTreeReaderArray<Int_t> Muon_fired_DST_DoubleMu3_noVtx_CaloScouting_v6 = {fReader, "Muon_fired_DST_DoubleMu3_noVtx_CaloScouting_v6"};
   TTreeReaderArray<Int_t> Muon_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3 = {fReader, "Muon_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13 = {fReader, "Muon_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1 = {fReader, "Muon_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_L2Mu23NoVtx_2Cha_v1 = {fReader, "Muon_fired_HLT_L2Mu23NoVtx_2Cha_v1"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu10p5_IP3p5 = {fReader, "Muon_fired_HLT_Mu10p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu12_IP6 = {fReader, "Muon_fired_HLT_Mu12_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu7_IP4 = {fReader, "Muon_fired_HLT_Mu7_IP4"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu7p5_Track7_Jpsi_v11 = {fReader, "Muon_fired_HLT_Mu7p5_Track7_Jpsi_v11"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP3 = {fReader, "Muon_fired_HLT_Mu8_IP3"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP5 = {fReader, "Muon_fired_HLT_Mu8_IP5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP6 = {fReader, "Muon_fired_HLT_Mu8_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_v1 = {fReader, "Muon_fired_HLT_Mu8_v1"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_v12 = {fReader, "Muon_fired_HLT_Mu8_v12"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8p5_IP3p5 = {fReader, "Muon_fired_HLT_Mu8p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP4 = {fReader, "Muon_fired_HLT_Mu9_IP4"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP5 = {fReader, "Muon_fired_HLT_Mu9_IP5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP6 = {fReader, "Muon_fired_HLT_Mu9_IP6"};
   TTreeReaderArray<Int_t> Muon_inTimeMuon = {fReader, "Muon_inTimeMuon"};
   TTreeReaderArray<Int_t> Muon_isDSAMuon = {fReader, "Muon_isDSAMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalMuon = {fReader, "Muon_isGlobalMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalNotTrackerMuon = {fReader, "Muon_isGlobalNotTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalOrTrackerMuon = {fReader, "Muon_isGlobalOrTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isPF = {fReader, "Muon_isPF"};
   TTreeReaderArray<Int_t> Muon_isSlimmedMuon = {fReader, "Muon_isSlimmedMuon"};
   TTreeReaderArray<Int_t> Muon_isTrackerMuon = {fReader, "Muon_isTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isTrackerNotGlobalMuon = {fReader, "Muon_isTrackerNotGlobalMuon"};
   TTreeReaderArray<Int_t> Muon_isTriggering = {fReader, "Muon_isTriggering"};
   TTreeReaderArray<Int_t> Muon_isTriggeringBPark = {fReader, "Muon_isTriggeringBPark"};
   TTreeReaderArray<Int_t> Muon_looseId = {fReader, "Muon_looseId"};
   TTreeReaderArray<Int_t> Muon_mediumId = {fReader, "Muon_mediumId"};
   TTreeReaderArray<Int_t> Muon_numberOfPixelLayers = {fReader, "Muon_numberOfPixelLayers"};
   TTreeReaderArray<Int_t> Muon_numberOfStations = {fReader, "Muon_numberOfStations"};
   TTreeReaderArray<Int_t> Muon_numberOfTrackerLayers = {fReader, "Muon_numberOfTrackerLayers"};
   TTreeReaderArray<Int_t> Muon_numberOfValidMuonHits = {fReader, "Muon_numberOfValidMuonHits"};
   TTreeReaderArray<Int_t> Muon_numberOfValidPixelHits = {fReader, "Muon_numberOfValidPixelHits"};
   TTreeReaderArray<Int_t> Muon_pdgId = {fReader, "Muon_pdgId"};
   TTreeReaderArray<Int_t> Muon_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2 = {fReader, "Muon_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2"};
   TTreeReaderArray<Int_t> Muon_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6 = {fReader, "Muon_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6"};
   TTreeReaderArray<Int_t> Muon_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3 = {fReader, "Muon_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13 = {fReader, "Muon_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1 = {fReader, "Muon_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_L2Mu23NoVtx_2Cha_v1 = {fReader, "Muon_prescale_HLT_L2Mu23NoVtx_2Cha_v1"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu10p5_IP3p5 = {fReader, "Muon_prescale_HLT_Mu10p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu12_IP6 = {fReader, "Muon_prescale_HLT_Mu12_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu7_IP4 = {fReader, "Muon_prescale_HLT_Mu7_IP4"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu7p5_Track7_Jpsi_v11 = {fReader, "Muon_prescale_HLT_Mu7p5_Track7_Jpsi_v11"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP3 = {fReader, "Muon_prescale_HLT_Mu8_IP3"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP5 = {fReader, "Muon_prescale_HLT_Mu8_IP5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP6 = {fReader, "Muon_prescale_HLT_Mu8_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_v1 = {fReader, "Muon_prescale_HLT_Mu8_v1"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_v12 = {fReader, "Muon_prescale_HLT_Mu8_v12"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8p5_IP3p5 = {fReader, "Muon_prescale_HLT_Mu8p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP4 = {fReader, "Muon_prescale_HLT_Mu9_IP4"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP5 = {fReader, "Muon_prescale_HLT_Mu9_IP5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP6 = {fReader, "Muon_prescale_HLT_Mu9_IP6"};
   TTreeReaderArray<Int_t> Muon_softId = {fReader, "Muon_softId"};
   TTreeReaderArray<Int_t> Muon_tightId = {fReader, "Muon_tightId"};
   TTreeReaderArray<Int_t> Muon_trackerHighPurityFlag = {fReader, "Muon_trackerHighPurityFlag"};
   TTreeReaderArray<Int_t> Muon_triggerIdLoose = {fReader, "Muon_triggerIdLoose"};
   TTreeReaderArray<UChar_t> Muon_pfIsoId = {fReader, "Muon_pfIsoId"};
   TTreeReaderArray<UChar_t> Muon_tkIsoId = {fReader, "Muon_tkIsoId"};
   //TTreeReaderValue<UInt_t> nTriggerMuon = {fReader, "nTriggerMuon"};
   //TTreeReaderArray<Float_t> TriggerMuon_eta = {fReader, "TriggerMuon_eta"};
   //TTreeReaderArray<Float_t> TriggerMuon_mass = {fReader, "TriggerMuon_mass"};
   //TTreeReaderArray<Float_t> TriggerMuon_phi = {fReader, "TriggerMuon_phi"};
   //TTreeReaderArray<Float_t> TriggerMuon_pt = {fReader, "TriggerMuon_pt"};
   //TTreeReaderArray<Float_t> TriggerMuon_vx = {fReader, "TriggerMuon_vx"};
   //TTreeReaderArray<Float_t> TriggerMuon_vy = {fReader, "TriggerMuon_vy"};
   //TTreeReaderArray<Float_t> TriggerMuon_vz = {fReader, "TriggerMuon_vz"};
   //TTreeReaderArray<Int_t> TriggerMuon_charge = {fReader, "TriggerMuon_charge"};
   //TTreeReaderArray<Int_t> TriggerMuon_pdgId = {fReader, "TriggerMuon_pdgId"};
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
   //TTreeReaderArray<Float_t> ProbeTracks_dzTrg = {fReader, "ProbeTracks_dzTrg"};
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
   //TTreeReaderValue<UChar_t> HLT_Mu7p5_L2Mu2_Jpsi_v10 = {fReader, "HLT_Mu7p5_L2Mu2_Jpsi_v10"};
   //TTreeReaderValue<UChar_t> HLT_Mu7p5_Track2_Jpsi_v11 = {fReader, "HLT_Mu7p5_Track2_Jpsi_v11"};
   //TTreeReaderValue<UChar_t> HLT_Mu7p5_Track7_Jpsi_v11 = {fReader, "HLT_Mu7p5_Track7_Jpsi_v11"};
   //TTreeReaderValue<UChar_t> HLT_Mu8_TrkIsoVVL_v12 = {fReader, "HLT_Mu8_TrkIsoVVL_v12"};
   //TTreeReaderValue<UChar_t> HLT_Mu8_v12 = {fReader, "HLT_Mu8_v12"};
   //TTreeReaderValue<UChar_t> HLT_Mu17_v13 = {fReader, "HLT_Mu17_v13"};
   //TTreeReaderValue<UChar_t> HLT_Mu19_v4 = {fReader, "HLT_Mu19_v4"};
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
   TTreeReaderArray<Float_t> SV_phi = {fReader, "SV_phi"};
   TTreeReaderArray<Float_t> SV_pt = {fReader, "SV_pt"};
   TTreeReaderArray<Float_t> SV_x = {fReader, "SV_x"};
   TTreeReaderArray<Float_t> SV_y = {fReader, "SV_y"};
   TTreeReaderArray<Float_t> SV_z = {fReader, "SV_z"};
   TTreeReaderValue<Float_t> Pileup_nTrueInt = {fReader, "PV_ndof"};


   TagAndProbeDumper(TTree * /*tree*/ =0) { }
   virtual ~TagAndProbeDumper() { }
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

   TFile* my_file;  

   Bool_t isMC;

   TF1* double_gauss_1;
   TF1* double_gauss_2;

   TF1* gauss_function_0p01;
   TF1* gauss_function_0p02;
   TF1* gauss_function_0p03;
   TF1* gauss_function_0p04;
   TF1* gauss_function_0p05;
   TF1* gauss_function_0p06;
   TF1* gauss_function_0p07;
   TF1* gauss_function_0p08;
   TF1* gauss_function_0p09;
   TF1* gauss_function_0p1;
   TF1* gauss_function_0p15;
   TF1* gauss_function_0p2;
   TF1* gauss_function_0p25;
   TF1* gauss_function_0p3;
   TF1* gauss_function_0p35;
   TF1* gauss_function_0p4;
   TF1* gauss_function_0p45;
   TF1* gauss_function_0p5;

   TTree* tree;

   Float_t the_pt;
   Float_t the_eta;
   Float_t the_phi;
   Float_t the_mass;
   Float_t the_cos2d;
   Float_t the_deltar;
   Float_t the_lxy;
   Float_t the_lxy_sig;
   Int_t the_ismatched;

   Float_t the_tag_pt;
   Float_t the_tag_eta;
   Float_t the_tag_phi;
   Float_t the_tag_dxy;
   Float_t the_tag_dxy_bs;
   Float_t the_tag_dz;
   Float_t the_tag_dxy_sig;
   Float_t the_tag_dxy_sig_bs;
   Float_t the_tag_dxy_sig_bs_rdst;
   Float_t the_tag_dz_sig;
   Int_t the_tag_isloose;
   Int_t the_tag_issoft;
   Int_t the_tag_ismedium;
   Int_t the_tag_istight;
   Int_t the_tag_fired_HLT_Mu7_IP4;
   Int_t the_tag_fired_HLT_Mu8_IP6;
   Int_t the_tag_fired_HLT_Mu8_IP5;
   Int_t the_tag_fired_HLT_Mu8_IP3;
   Int_t the_tag_fired_HLT_Mu8p5_IP3p5;
   Int_t the_tag_fired_HLT_Mu9_IP6;
   Int_t the_tag_fired_HLT_Mu9_IP5;
   Int_t the_tag_fired_HLT_Mu9_IP4;
   Int_t the_tag_fired_HLT_Mu10p5_IP3p5;
   Int_t the_tag_fired_HLT_Mu12_IP6;
   Int_t the_tag_fired_HLT_Mu8_v1;
   Int_t the_tag_fired_HLT_Mu8_v12;
   Int_t the_tag_fired_HLT_Mu7p5_Track7_Jpsi_v11;
   Int_t the_tag_fired_HLT_L2Mu23NoVtx_2Cha_v1;
   Int_t the_tag_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1;
   Int_t the_tag_fired_DST_DoubleMu1_noVtx_CaloScouting_v2;
   Int_t the_tag_fired_DST_DoubleMu3_noVtx_CaloScouting_v6;
   Int_t the_tag_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3;
   Int_t the_tag_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13;
   Int_t the_tag_prescale_HLT_Mu7_IP4;
   Int_t the_tag_prescale_HLT_Mu8_IP6;
   Int_t the_tag_prescale_HLT_Mu8_IP5;
   Int_t the_tag_prescale_HLT_Mu8_IP3;
   Int_t the_tag_prescale_HLT_Mu8p5_IP3p5;
   Int_t the_tag_prescale_HLT_Mu9_IP6;
   Int_t the_tag_prescale_HLT_Mu9_IP5;
   Int_t the_tag_prescale_HLT_Mu9_IP4;
   Int_t the_tag_prescale_HLT_Mu10p5_IP3p5;
   Int_t the_tag_prescale_HLT_Mu12_IP6;
   Int_t the_tag_prescale_HLT_Mu8_v1;
   Int_t the_tag_prescale_HLT_Mu8_v12;
   Int_t the_tag_prescale_HLT_Mu7p5_Track7_Jpsi_v11;
   Int_t the_tag_prescale_HLT_L2Mu23NoVtx_2Cha_v1;
   Int_t the_tag_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1;
   Int_t the_tag_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2;
   Int_t the_tag_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6;
   Int_t the_tag_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3;
   Int_t the_tag_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13;

   Float_t the_probe_pt;
   Float_t the_probe_eta;
   Float_t the_probe_phi;
   Float_t the_probe_dxy;
   Float_t the_probe_dxy_bs;
   Float_t the_probe_dxy_bs_uncorrected;
   Float_t the_weight_dxy_bs;
   Float_t the_probe_dz;
   Float_t the_probe_dxy_sig;

   //Float_t the_smeared_corr_0p01;
   //Float_t the_smeared_corr_0p02;
   Float_t the_smeared_corr_0p03;
   //Float_t the_smeared_corr_0p04;
   //Float_t the_smeared_corr_0p05;
   //Float_t the_smeared_corr_0p06;
   //Float_t the_smeared_corr_0p07;
   //Float_t the_smeared_corr_0p08;
   //Float_t the_smeared_corr_0p09;
   //Float_t the_smeared_corr_0p1;
   //Float_t the_smeared_corr_0p15;
   //Float_t the_smeared_corr_0p2;
   //Float_t the_smeared_corr_0p25;
   //Float_t the_smeared_corr_0p3;
   //Float_t the_smeared_corr_0p35;
   //Float_t the_smeared_corr_0p4;
   //Float_t the_smeared_corr_0p45;
   //Float_t the_smeared_corr_0p5;

   //Float_t the_double_gauss_corr_1;
   //Float_t the_double_gauss_corr_2;

   Float_t the_probe_dxy_sig_bs_uncorrected;
   Float_t the_probe_dxy_sig_bs_corrected_linearscale;
   Float_t the_weight_dxy_sig_bs;

   Float_t the_probe_dxy_sig_bs;
   Float_t the_probe_dxy_sig_bs_rdst;
   Float_t the_probe_dz_sig;
   Int_t the_probe_isloose;
   Int_t the_probe_issoft;
   Int_t the_probe_ismedium;
   Int_t the_probe_istight;
   Int_t the_probe_fired_HLT_Mu7_IP4;
   Int_t the_probe_fired_HLT_Mu8_IP6;
   Int_t the_probe_fired_HLT_Mu8_IP5;
   Int_t the_probe_fired_HLT_Mu8_IP3;
   Int_t the_probe_fired_HLT_Mu8p5_IP3p5;
   Int_t the_probe_fired_HLT_Mu9_IP6;
   Int_t the_probe_fired_HLT_Mu9_IP5;
   Int_t the_probe_fired_HLT_Mu9_IP4;
   Int_t the_probe_fired_HLT_Mu10p5_IP3p5;
   Int_t the_probe_fired_HLT_Mu12_IP6;
   Int_t the_probe_fired_HLT_Mu8_v1;
   Int_t the_probe_fired_HLT_Mu8_v12;
   Int_t the_probe_fired_HLT_Mu7p5_Track7_Jpsi_v11;
   Int_t the_probe_fired_HLT_L2Mu23NoVtx_2Cha_v1;
   Int_t the_probe_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1;
   Int_t the_probe_fired_DST_DoubleMu1_noVtx_CaloScouting_v2;
   Int_t the_probe_fired_DST_DoubleMu3_noVtx_CaloScouting_v6;
   Int_t the_probe_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3;
   Int_t the_probe_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13;
   Int_t the_probe_prescale_HLT_Mu7_IP4;
   Int_t the_probe_prescale_HLT_Mu8_IP6;
   Int_t the_probe_prescale_HLT_Mu8_IP5;
   Int_t the_probe_prescale_HLT_Mu8_IP3;
   Int_t the_probe_prescale_HLT_Mu8p5_IP3p5;
   Int_t the_probe_prescale_HLT_Mu9_IP6;
   Int_t the_probe_prescale_HLT_Mu9_IP5;
   Int_t the_probe_prescale_HLT_Mu9_IP4;
   Int_t the_probe_prescale_HLT_Mu10p5_IP3p5;
   Int_t the_probe_prescale_HLT_Mu12_IP6;
   Int_t the_probe_prescale_HLT_Mu8_v1;
   Int_t the_probe_prescale_HLT_Mu8_v12;
   Int_t the_probe_prescale_HLT_Mu7p5_Track7_Jpsi_v11;
   Int_t the_probe_prescale_HLT_L2Mu23NoVtx_2Cha_v1;
   Int_t the_probe_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1;
   Int_t the_probe_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2;
   Int_t the_probe_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6;
   Int_t the_probe_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3;
   Int_t the_probe_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13;

   Int_t the_probe_fired_BParkingHLT;

   Int_t the_pv_npvs;

   Float_t the_weight_hlt_A1;
   Float_t the_weight_hlt_A1_6;
   Float_t the_weight_pu_A;
   Float_t the_weight_pu_B;
   Float_t the_weight_pu_C;
   Float_t the_weight_pu_D;
   Float_t the_weight_pu_tot;

   ClassDef(TagAndProbeDumper,0);

};

#endif

#ifdef TagAndProbeDumper_cxx
void TagAndProbeDumper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t TagAndProbeDumper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef TagAndProbeDumper_cxx
