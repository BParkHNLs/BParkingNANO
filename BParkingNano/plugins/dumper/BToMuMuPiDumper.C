#define BToMuMuPiDumper_cxx
// The class definition in BToMuMuPiDumper.h has been generated automatically
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
// root> T->Process("BToMuMuPiDumper.C")
// root> T->Process("BToMuMuPiDumper.C","some options")
// root> T->Process("BToMuMuPiDumper.C+")
//


#include "BToMuMuPiDumper.h"
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


void BToMuMuPiDumper::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  cout << " --------------------------" << endl;
  cout << "       B->MuMuPi Dumper    " << endl;
  cout << " --------------------------" << endl;
}


void BToMuMuPiDumper::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  TString outFileName = option;

  isSignalMC = false;
  isMC = false;

  if(outFileName.Contains("isSignalMC")){
    isMC = true;
    isSignalMC = true;
    outFileName.Resize(outFileName.Length()-11);
  }
  else if(outFileName.Contains("isMC")){
    isMC = true;
    outFileName.Resize(outFileName.Length()-5);
  }

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
  signal_tree = new TTree("signal_tree", "signal_tree");

  signal_tree->Branch("event", &the_event);
  signal_tree->Branch("run", &the_run);
  signal_tree->Branch("lumi", &the_lumi);

  signal_tree->Branch("b_pt", &the_sig_b_pt);
  signal_tree->Branch("b_eta", &the_sig_b_eta);
  signal_tree->Branch("b_phi", &the_sig_b_phi);
  signal_tree->Branch("b_mass", &the_sig_b_mass);
  signal_tree->Branch("b_charge", &the_sig_b_charge);
  signal_tree->Branch("b_pdgid", &the_sig_b_pdgid);

  signal_tree->Branch("hnl_pt", &the_sig_hnl_pt);
  signal_tree->Branch("hnl_eta", &the_sig_hnl_eta);
  signal_tree->Branch("hnl_phi", &the_sig_hnl_phi);
  signal_tree->Branch("hnl_mass", &the_sig_hnl_mass);
  signal_tree->Branch("hnl_charge", &the_sig_hnl_charge);
  signal_tree->Branch("hnl_ct", &the_sig_hnl_ct);
  signal_tree->Branch("hnl_cos2d", &the_sig_hnl_cos2d);
  //signal_tree->Branch("hnl_iso03", &the_sig_hnl_iso03);
  //signal_tree->Branch("hnl_iso03_close", &the_sig_hnl_iso03_close);
  //signal_tree->Branch("hnl_iso03_rel_close", &the_sig_hnl_iso03_rel_close);
  //signal_tree->Branch("hnl_iso04", &the_sig_hnl_iso04);
  //signal_tree->Branch("hnl_iso04_close", &the_sig_hnl_iso04_close);
  //signal_tree->Branch("hnl_iso04_rel_close", &the_sig_hnl_iso04_rel_close);

  signal_tree->Branch("mu0_pt", &the_sig_mu0_pt);
  signal_tree->Branch("mu0_eta", &the_sig_mu0_eta);
  signal_tree->Branch("mu0_phi", &the_sig_mu0_phi);
  signal_tree->Branch("mu0_charge", &the_sig_mu0_charge);
  signal_tree->Branch("mu0_dxy", &the_sig_mu0_dxy);
  signal_tree->Branch("mu0_dxy_bs", &the_sig_mu0_dxy_bs);
  signal_tree->Branch("mu0_dxysig", &the_sig_mu0_dxysig);
  signal_tree->Branch("mu0_dxysig_bs", &the_sig_mu0_dxysig_bs);
  signal_tree->Branch("mu0_dxysig_bs_linearscale", &the_sig_mu0_dxysig_bs_linearscale);
  signal_tree->Branch("mu0_dxysig_bs_scale1p12_smear0p03", &the_sig_mu0_dxysig_bs_scale1p12_smear0p03);
  //signal_tree->Branch("mu0_dxysig_bs_smeared", &the_sig_mu0_dxysig_bs_smeared);
  signal_tree->Branch("mu0_dxysig_bs_rdst", &the_sig_mu0_dxysig_bs_rdst);
  signal_tree->Branch("mu0_dz", &the_sig_mu0_dz);
  signal_tree->Branch("mu0_dzsig", &the_sig_mu0_dzsig);
  signal_tree->Branch("mu0_triggermatching_dr", &the_sig_mu0_triggermatching_dr);
  signal_tree->Branch("mu0_triggermatching_dpt", &the_sig_mu0_triggermatching_dpt);
  //signal_tree->Branch("mu0_ip3d", &the_sig_mu0_ip3d);
  //signal_tree->Branch("mu0_ip3dsig", &the_sig_mu0_ip3dsig);
  signal_tree->Branch("mu0_pfiso03", &the_sig_mu0_pfiso03);
  signal_tree->Branch("mu0_pfiso03_rel", &the_sig_mu0_pfiso03_rel);
  //signal_tree->Branch("mu0_iso03", &the_sig_mu0_iso03);
  //signal_tree->Branch("mu0_iso03_close", &the_sig_mu0_iso03_close);
  //signal_tree->Branch("mu0_iso03_rel_close", &the_sig_mu0_iso03_rel_close);
  signal_tree->Branch("mu0_pfiso04", &the_sig_mu0_pfiso04);
  signal_tree->Branch("mu0_pfiso04_rel", &the_sig_mu0_pfiso04_rel);
  //signal_tree->Branch("mu0_iso04", &the_sig_mu0_iso04);
  //signal_tree->Branch("mu0_iso04_close", &the_sig_mu0_iso04_close);
  //signal_tree->Branch("mu0_iso04_rel_close", &the_sig_mu0_iso04_rel_close);
  signal_tree->Branch("mu0_looseid", &the_sig_mu0_looseid);
  signal_tree->Branch("mu0_mediumid", &the_sig_mu0_mediumid);
  signal_tree->Branch("mu0_tightid", &the_sig_mu0_tightid);
  signal_tree->Branch("mu0_softid", &the_sig_mu0_softid);
  signal_tree->Branch("mu0_pfisoid", &the_sig_mu0_pfisoid);
  signal_tree->Branch("mu0_trkisoid", &the_sig_mu0_trkisoid);
  signal_tree->Branch("mu0_triggerlooseid", &the_sig_mu0_triggerlooseid);
  signal_tree->Branch("mu0_istriggering", &the_sig_mu0_istriggering);
  signal_tree->Branch("mu0_isPF", &the_sig_mu0_isPF);
  signal_tree->Branch("mu0_isglobalmuon", &the_sig_mu0_isglobalmuon);
  signal_tree->Branch("mu0_istrackermuon", &the_sig_mu0_istrackermuon);
  signal_tree->Branch("mu0_isglobalortrackermuon", &the_sig_mu0_isglobalortrackermuon);
  signal_tree->Branch("mu0_isglobalnottrackermuon", &the_sig_mu0_isglobalnottrackermuon);
  signal_tree->Branch("mu0_istrackernotglobalmuon", &the_sig_mu0_istrackernotglobalmuon);
  signal_tree->Branch("mu0_intimemuon", &the_sig_mu0_intimemuon);
  //signal_tree->Branch("mu0_segmentcompatibility", &the_sig_mu0_segmentcompatibility);
  //signal_tree->Branch("mu0_calocompatibility", &the_sig_mu0_calocompatibility);
  //signal_tree->Branch("mu0_validhitfraction", &the_sig_mu0_validhitfraction);
  //signal_tree->Branch("mu0_kinkfinderchi2", &the_sig_mu0_kinkfinderchi2);
  //signal_tree->Branch("mu0_globalnormalisedchi2", &the_sig_mu0_globalnormalisedchi2);
  //signal_tree->Branch("mu0_localpositionchi2", &the_sig_mu0_localpositionchi2);
  //signal_tree->Branch("mu0_trackerhighpurityflag", &the_sig_mu0_trackerhighpurityflag);
  //signal_tree->Branch("mu0_numberofvalidmuonhits", &the_sig_mu0_numberofvalidmuonhits);
  //signal_tree->Branch("mu0_numberofvalidpixelhits", &the_sig_mu0_numberofvalidpixelhits);
  //signal_tree->Branch("mu0_numberoftrackerlayers", &the_sig_mu0_numberoftrackerlayers);
  //signal_tree->Branch("mu0_numberofpixellayers", &the_sig_mu0_numberofpixellayers);
  //signal_tree->Branch("mu0_numberofstations", &the_sig_mu0_numberofstations);
  //signal_tree->Branch("mu0_fired_hlt_mu7_ip4", &the_sig_mu0_fired_hlt_mu7_ip4);
  //signal_tree->Branch("mu0_fired_hlt_mu8_ip3", &the_sig_mu0_fired_hlt_mu8_ip3);
  //signal_tree->Branch("mu0_fired_hlt_mu8_ip5", &the_sig_mu0_fired_hlt_mu8_ip5);
  //signal_tree->Branch("mu0_fired_hlt_mu8_ip6", &the_sig_mu0_fired_hlt_mu8_ip6);
  //signal_tree->Branch("mu0_fired_hlt_mu8p5_ip3p5", &the_sig_mu0_fired_hlt_mu8p5_ip3p5);
  //signal_tree->Branch("mu0_fired_hlt_mu9_ip4", &the_sig_mu0_fired_hlt_mu9_ip4);
  //signal_tree->Branch("mu0_fired_hlt_mu9_ip5", &the_sig_mu0_fired_hlt_mu9_ip5);
  //signal_tree->Branch("mu0_fired_hlt_mu9_ip6", &the_sig_mu0_fired_hlt_mu9_ip6);
  //signal_tree->Branch("mu0_fired_hlt_mu10p5_ip3p5", &the_sig_mu0_fired_hlt_mu10p5_ip3p5);
  //signal_tree->Branch("mu0_fired_hlt_mu12_ip6", &the_sig_mu0_fired_hlt_mu12_ip6);
  //signal_tree->Branch("mu0_prescale_hlt_mu7_ip4", &the_sig_mu0_prescale_hlt_mu7_ip4);
  //signal_tree->Branch("mu0_prescale_hlt_mu8_ip3", &the_sig_mu0_prescale_hlt_mu8_ip3);
  //signal_tree->Branch("mu0_prescale_hlt_mu8_ip5", &the_sig_mu0_prescale_hlt_mu8_ip5);
  //signal_tree->Branch("mu0_prescale_hlt_mu8_ip6", &the_sig_mu0_prescale_hlt_mu8_ip6);
  //signal_tree->Branch("mu0_prescale_hlt_mu8p5_ip3p5", &the_sig_mu0_prescale_hlt_mu8p5_ip3p5);
  //signal_tree->Branch("mu0_prescale_hlt_mu9_ip4", &the_sig_mu0_prescale_hlt_mu9_ip4);
  //signal_tree->Branch("mu0_prescale_hlt_mu9_ip5", &the_sig_mu0_prescale_hlt_mu9_ip5);
  //signal_tree->Branch("mu0_prescale_hlt_mu9_ip6", &the_sig_mu0_prescale_hlt_mu9_ip6);
  //signal_tree->Branch("mu0_prescale_hlt_mu10p5_ip3p5", &the_sig_mu0_prescale_hlt_mu10p5_ip3p5);
  //signal_tree->Branch("mu0_prescale_hlt_mu12_ip6", &the_sig_mu0_prescale_hlt_mu12_ip6);

  signal_tree->Branch("mu0_triggering_pt", &the_sig_mu0_triggering_pt);
  signal_tree->Branch("mu0_triggering_eta", &the_sig_mu0_triggering_eta);
  signal_tree->Branch("mu0_triggering_phi", &the_sig_mu0_triggering_phi);
  signal_tree->Branch("mu0_triggering_dxy", &the_sig_mu0_triggering_dxy);
  signal_tree->Branch("mu0_triggering_dxy_bs", &the_sig_mu0_triggering_dxy_bs);
  signal_tree->Branch("mu0_triggering_dxysig", &the_sig_mu0_triggering_dxysig);
  signal_tree->Branch("mu0_triggering_dxysig_bs", &the_sig_mu0_triggering_dxysig_bs);
  signal_tree->Branch("mu0_triggering_dxysig_bs_linearscale", &the_sig_mu0_triggering_dxysig_bs_linearscale);
  signal_tree->Branch("mu0_triggering_dxysig_bs_scale1p12_smear0p03", &the_sig_mu0_triggering_dxysig_bs_scale1p12_smear0p03);
  signal_tree->Branch("mu0_triggering_dxysig_bs_rdst", &the_sig_mu0_triggering_dxysig_bs_rdst);

  signal_tree->Branch("mu_pt", &the_sig_mu_pt);
  signal_tree->Branch("mu_eta", &the_sig_mu_eta);
  signal_tree->Branch("mu_phi", &the_sig_mu_phi);
  signal_tree->Branch("mu_charge", &the_sig_mu_charge);
  signal_tree->Branch("mu_dxy", &the_sig_mu_dxy);
  signal_tree->Branch("mu_dxy_bs", &the_sig_mu_dxy_bs);
  signal_tree->Branch("mu_dxysig", &the_sig_mu_dxysig);
  signal_tree->Branch("mu_dxysig_bs", &the_sig_mu_dxysig_bs);
  signal_tree->Branch("mu_dxysig_bs_linearscale", &the_sig_mu_dxysig_bs_linearscale);
  signal_tree->Branch("mu_dxysig_bs_scale1p12_smear0p03", &the_sig_mu_dxysig_bs_scale1p12_smear0p03);
  signal_tree->Branch("mu_dxysig_bs_rdst", &the_sig_mu_dxysig_bs_rdst);
  signal_tree->Branch("mu_dz", &the_sig_mu_dz);
  signal_tree->Branch("mu_dzsig", &the_sig_mu_dzsig);
  signal_tree->Branch("mu_triggermatching_dr", &the_sig_mu_triggermatching_dr);
  signal_tree->Branch("mu_triggermatching_dpt", &the_sig_mu_triggermatching_dpt);
  signal_tree->Branch("mu_ismatchedtoslimmedmuon", &the_sig_mu_ismatchedtoslimmedmuon);
  signal_tree->Branch("mu_indexmatchedslimmedmuon", &the_sig_mu_indexmatchedslimmedmuon);
  //signal_tree->Branch("mu_dsatoslimmedmatching_deltar", &the_sig_mu_dsatoslimmedmatching_deltar);
  //signal_tree->Branch("mu_dsatoslimmedmatching_deltaptrel", &the_sig_mu_dsatoslimmedmatching_deltaptrel);
  //signal_tree->Branch("mu_dsatoslimmedmatching_deltadxyrel", &the_sig_mu_dsatoslimmedmatching_deltadxyrel);
  //signal_tree->Branch("mu_dsatoslimmedmatching_deltadzrel", &the_sig_mu_dsatoslimmedmatching_deltadzrel);
  signal_tree->Branch("mu_passdsaid", &the_sig_mu_passdsaid);
  //signal_tree->Branch("mu_ip3d", &the_sig_mu_ip3d);
  //signal_tree->Branch("mu_ip3dsig", &the_sig_mu_ip3dsig);
  signal_tree->Branch("mu_pfiso03", &the_sig_mu_pfiso03);
  signal_tree->Branch("mu_pfiso03_rel", &the_sig_mu_pfiso03_rel);
  //signal_tree->Branch("mu_iso03", &the_sig_mu_iso03);
  //signal_tree->Branch("mu_iso03_close", &the_sig_mu_iso03_close);
  //signal_tree->Branch("mu_iso03_rel_close", &the_sig_mu_iso03_rel_close);
  signal_tree->Branch("mu_pfiso04", &the_sig_mu_pfiso04);
  signal_tree->Branch("mu_pfiso04_rel", &the_sig_mu_pfiso04_rel);
  //signal_tree->Branch("mu_iso04", &the_sig_mu_iso04);
  //signal_tree->Branch("mu_iso04_close", &the_sig_mu_iso04_close);
  //signal_tree->Branch("mu_iso04_rel_close", &the_sig_mu_iso04_rel_close);
  //signal_tree->Branch("mu_isloose", &the_sig_mu_isloose);
  //signal_tree->Branch("mu_ismedium", &the_sig_mu_ismedium);
  //signal_tree->Branch("mu_istight", &the_sig_mu_istight);
  //signal_tree->Branch("mu_issoft", &the_sig_mu_issoft);
  //signal_tree->Branch("mu_istriggering", &the_sig_mu_istriggering);
  signal_tree->Branch("mu_looseid", &the_sig_mu_looseid);
  signal_tree->Branch("mu_mediumid", &the_sig_mu_mediumid);
  signal_tree->Branch("mu_tightid", &the_sig_mu_tightid);
  signal_tree->Branch("mu_softid", &the_sig_mu_softid);
  signal_tree->Branch("mu_pfisoid", &the_sig_mu_pfisoid);
  signal_tree->Branch("mu_trkisoid", &the_sig_mu_trkisoid);
  signal_tree->Branch("mu_triggerlooseid", &the_sig_mu_triggerlooseid);
  signal_tree->Branch("mu_whnlid", &the_sig_mu_whnlid);
  signal_tree->Branch("mu_customisedid", &the_sig_mu_customisedid);
  signal_tree->Branch("mu_istriggering", &the_sig_mu_istriggering);
  signal_tree->Branch("mu_isslimmed", &the_sig_mu_isslimmed);
  signal_tree->Branch("mu_isdsa", &the_sig_mu_isdsa);
  signal_tree->Branch("mu_isPF", &the_sig_mu_isPF);
  signal_tree->Branch("mu_isglobalmuon", &the_sig_mu_isglobalmuon);
  signal_tree->Branch("mu_istrackermuon", &the_sig_mu_istrackermuon);
  //signal_tree->Branch("mu_isglobalortrackermuon", &the_sig_mu_isglobalortrackermuon);
  //signal_tree->Branch("mu_isglobalnottrackermuon", &the_sig_mu_isglobalnottrackermuon);
  //signal_tree->Branch("mu_istrackernotglobalmuon", &the_sig_mu_istrackernotglobalmuon);
  //signal_tree->Branch("mu_intimemuon", &the_sig_mu_intimemuon);
  signal_tree->Branch("mu_segmentcompatibility", &the_sig_mu_segmentcompatibility);
  signal_tree->Branch("mu_calocompatibility", &the_sig_mu_calocompatibility);
  signal_tree->Branch("mu_validhitfraction", &the_sig_mu_validhitfraction);
  signal_tree->Branch("mu_kinkfinderchi2", &the_sig_mu_kinkfinderchi2);
  //signal_tree->Branch("mu_globalnormalisedchi2", &the_sig_mu_globalnormalisedchi2);
  signal_tree->Branch("mu_localpositionchi2", &the_sig_mu_localpositionchi2);
  signal_tree->Branch("mu_trackerhighpurityflag", &the_sig_mu_trackerhighpurityflag);
  //signal_tree->Branch("mu_numberofvalidmuonhits", &the_sig_mu_numberofvalidmuonhits);
  signal_tree->Branch("mu_numberofvalidpixelhits", &the_sig_mu_numberofvalidpixelhits);
  signal_tree->Branch("mu_numberoftrackerlayers", &the_sig_mu_numberoftrackerlayers);
  //signal_tree->Branch("mu_numberofpixellayers", &the_sig_mu_numberofpixellayers);
  signal_tree->Branch("mu_numberofstations", &the_sig_mu_numberofstations);

  signal_tree->Branch("mu_triggering_pt", &the_sig_mu_triggering_pt);
  signal_tree->Branch("mu_triggering_eta", &the_sig_mu_triggering_eta);
  signal_tree->Branch("mu_triggering_phi", &the_sig_mu_triggering_phi);
  signal_tree->Branch("mu_triggering_charge", &the_sig_mu_triggering_charge);
  signal_tree->Branch("mu_triggering_dxy", &the_sig_mu_triggering_dxy);
  signal_tree->Branch("mu_triggering_dxy_bs", &the_sig_mu_triggering_dxy_bs);
  signal_tree->Branch("mu_triggering_dxysig", &the_sig_mu_triggering_dxysig);
  signal_tree->Branch("mu_triggering_dxysig_bs", &the_sig_mu_triggering_dxysig_bs);
  signal_tree->Branch("mu_triggering_dxysig_bs_linearscale", &the_sig_mu_triggering_dxysig_bs_linearscale);
  signal_tree->Branch("mu_triggering_dxysig_bs_scale1p12_smear0p03", &the_sig_mu_triggering_dxysig_bs_scale1p12_smear0p03);
  signal_tree->Branch("mu_triggering_dxysig_bs_rdst", &the_sig_mu_triggering_dxysig_bs_rdst);

  signal_tree->Branch("pi_pt", &the_sig_pi_pt);
  signal_tree->Branch("pi_eta", &the_sig_pi_eta);
  signal_tree->Branch("pi_phi", &the_sig_pi_phi);
  signal_tree->Branch("pi_charge", &the_sig_pi_charge);
  signal_tree->Branch("pi_dcasig", &the_sig_pi_dcasig);
  signal_tree->Branch("pi_dxy", &the_sig_pi_dxy);
  signal_tree->Branch("pi_dz", &the_sig_pi_dz);
  signal_tree->Branch("pi_dxysig", &the_sig_pi_dxysig);
  signal_tree->Branch("pi_dzsig", &the_sig_pi_dzsig);
  signal_tree->Branch("pi_ispacked", &the_sig_pi_ispacked);
  signal_tree->Branch("pi_islost", &the_sig_pi_islost);
  //signal_tree->Branch("pi_mu0_dr", &the_sig_pi_mu0_dr);
  //signal_tree->Branch("pi_ismatchedtomuon", &the_sig_pi_ismatchedtomuon);
  //signal_tree->Branch("pi_chi2", &the_sig_pi_chi2);
  //signal_tree->Branch("pi_normalisedchi2", &the_sig_pi_normalisedChi2);
  //signal_tree->Branch("pi_validfraction", &the_sig_pi_validFraction);
  //signal_tree->Branch("pi_ndof", &the_sig_pi_ndof);
  //signal_tree->Branch("pi_numberofvalidhits", &the_sig_pi_numberOfValidHits);
  //signal_tree->Branch("pi_numberoflosthits", &the_sig_pi_numberOfLostHits);
  //signal_tree->Branch("pi_numberofvalidpixelhits", &the_sig_pi_numberOfValidPixelHits);
  //signal_tree->Branch("pi_numberoftrackerlayers", &the_sig_pi_numberOfTrackerLayers);
  //signal_tree->Branch("pi_numberofpixellayers", &the_sig_pi_numberOfPixelLayers);
  //signal_tree->Branch("pi_qualityindex", &the_sig_pi_qualityIndex);
  signal_tree->Branch("pi_highpurityflag", &the_sig_pi_highPurityFlag);
  signal_tree->Branch("pi_packedcandhashighpurity", &the_sig_pi_packedcandhashighpurity);
  //signal_tree->Branch("pi_matchedtomuon_loose", &the_sig_pi_matchedtomuon_loose);
  //signal_tree->Branch("pi_matchedtomuon_medium", &the_sig_pi_matchedtomuon_medium);
  //signal_tree->Branch("pi_matchedtomuon_tight", &the_sig_pi_matchedtomuon_tight);

  //signal_tree->Branch("dimu_mass", &the_sig_dimu_mass);
  //signal_tree->Branch("dimu_pt", &the_sig_dimu_pt);
  signal_tree->Branch("mu0_mu_mass", &the_sig_mu0_mu_mass);
  signal_tree->Branch("mu0_mu_pt", &the_sig_mu0_mu_pt);
  signal_tree->Branch("mu0_pi_mass", &the_sig_mu0_pi_mass);
  signal_tree->Branch("mu0_pi_pt", &the_sig_mu0_pi_pt);
  signal_tree->Branch("dimu_lxy", &the_sig_dimu_lxy);
  signal_tree->Branch("dimu_lxyz", &the_sig_dimu_lxyz);
  //signal_tree->Branch("dimu_vxdiff", &the_sig_dimu_vxdiff);
  //signal_tree->Branch("dimu_vydiff", &the_sig_dimu_vydiff);
  //signal_tree->Branch("dimu_vzdiff", &the_sig_dimu_vzdiff);

  //signal_tree->Branch("cos_theta_star_pion", &the_sig_cos_theta_star_pion);
  //signal_tree->Branch("cos_theta_star_muon", &the_sig_cos_theta_star_muon);
  //signal_tree->Branch("cos_theta_star_sum", &the_sig_cos_theta_star_sum);

  //signal_tree->Branch("px_diff_hnl_daughters_lab", &the_sig_px_diff_hnl_daughters_lab);
  //signal_tree->Branch("py_diff_hnl_daughters_lab", &the_sig_py_diff_hnl_daughters_lab);
  //signal_tree->Branch("pz_diff_hnl_daughters_lab", &the_sig_pz_diff_hnl_daughters_lab);
  //signal_tree->Branch("energy_diff_prefithnl_daughters_lab", &the_sig_energy_diff_prefithnl_daughters_lab);
  //signal_tree->Branch("px_diff_prefithnl_daughters_lab", &the_sig_px_diff_prefithnl_daughters_lab);
  //signal_tree->Branch("py_diff_prefithnl_daughters_lab", &the_sig_py_diff_prefithnl_daughters_lab);
  //signal_tree->Branch("pz_diff_prefithnl_daughters_lab", &the_sig_pz_diff_prefithnl_daughters_lab);

  signal_tree->Branch("deltar_mu_pi", &the_sig_deltar_mu_pi);
  signal_tree->Branch("deltar_mu0_hnl", &the_sig_deltar_mu0_hnl);
  signal_tree->Branch("deltar_mu0_mu", &the_sig_deltar_mu0_mu);
  signal_tree->Branch("deltar_mu0_pi", &the_sig_deltar_mu0_pi);
  //signal_tree->Branch("deltaeta_mu_pi", &the_sig_deltaeta_mu_pi);
  //signal_tree->Branch("deltaeta_mu0_mu", &the_sig_deltaeta_mu0_mu);
  //signal_tree->Branch("deltaeta_mu0_pi", &the_sig_deltaeta_mu0_pi);
  //signal_tree->Branch("deltaeta_mu0_hnl", &the_sig_deltaeta_mu0_hnl);
  //signal_tree->Branch("deltaphi_mu_pi", &the_sig_deltaphi_mu_pi);
  //signal_tree->Branch("deltaphi_mu0_mu", &the_sig_deltaphi_mu0_mu);
  //signal_tree->Branch("deltaphi_mu0_pi", &the_sig_deltaphi_mu0_pi);
  //signal_tree->Branch("deltaphi_mu0_hnl", &the_sig_deltaphi_mu0_hnl);

  //signal_tree->Branch("deltae_pi_fit_pi", &the_sig_deltapt_pi_fit_pi);
  //signal_tree->Branch("deltae_mu_fit_mu", &the_sig_deltapt_mu_fit_mu);
  //signal_tree->Branch("deltae_hnl_fit_hnl", &the_sig_deltapt_hnl_fit_hnl);
  //signal_tree->Branch("deltapt_pi_fit_pi", &the_sig_deltapt_pi_fit_pi);
  //signal_tree->Branch("deltapt_mu_fit_mu", &the_sig_deltapt_mu_fit_mu);
  //signal_tree->Branch("deltapt_hnl_fit_hnl", &the_sig_deltapt_hnl_fit_hnl);
  //signal_tree->Branch("deltapx_pi_fit_pi", &the_sig_deltapx_pi_fit_pi);
  //signal_tree->Branch("deltapx_mu_fit_mu", &the_sig_deltapx_mu_fit_mu);
  //signal_tree->Branch("deltapx_hnl_fit_hnl", &the_sig_deltapx_hnl_fit_hnl);
  //signal_tree->Branch("deltapy_pi_fit_pi", &the_sig_deltapy_pi_fit_pi);
  //signal_tree->Branch("deltapy_mu_fit_mu", &the_sig_deltapy_mu_fit_mu);
  //signal_tree->Branch("deltapy_hnl_fit_hnl", &the_sig_deltapy_hnl_fit_hnl);
  //signal_tree->Branch("deltapz_pi_fit_pi", &the_sig_deltapz_pi_fit_pi);
  //signal_tree->Branch("deltapz_mu_fit_mu", &the_sig_deltapz_mu_fit_mu);
  //signal_tree->Branch("deltapz_hnl_fit_hnl", &the_sig_deltapz_hnl_fit_hnl);
  //signal_tree->Branch("deltaeta_pi_fit_pi", &the_sig_deltaeta_pi_fit_pi);
  //signal_tree->Branch("deltaeta_mu_fit_mu", &the_sig_deltaeta_mu_fit_mu);
  //signal_tree->Branch("deltaeta_hnl_fit_hnl", &the_sig_deltaeta_hnl_fit_hnl);
  //signal_tree->Branch("deltaphi_pi_fit_pi", &the_sig_deltaphi_pi_fit_pi);
  //signal_tree->Branch("deltaphi_mu_fit_mu", &the_sig_deltaphi_mu_fit_mu);
  //signal_tree->Branch("deltaphi_hnl_fit_hnl", &the_sig_deltaphi_hnl_fit_hnl);

  signal_tree->Branch("sv_chi2", &the_sig_sv_chi2);
  signal_tree->Branch("sv_lxy", &the_sig_sv_lxy);
  //signal_tree->Branch("sv_pv_lxy", &the_sig_sv_pv_lxy);
  //signal_tree->Branch("sv_pv_lxyz", &the_sig_sv_pv_lxyz);
  signal_tree->Branch("sv_lxysig", &the_sig_sv_lxysig);
  signal_tree->Branch("sv_lxyz", &the_sig_sv_lxyz);
  signal_tree->Branch("sv_prob", &the_sig_sv_prob);
  //signal_tree->Branch("sv_x", &the_sig_sv_x);
  //signal_tree->Branch("sv_y", &the_sig_sv_y);
  //signal_tree->Branch("sv_z", &the_sig_sv_z);

  signal_tree->Branch("pi_mu0_vzdiff", &the_sig_pi_mu0_vzdiff);

  signal_tree->Branch("ismatched", &the_sig_ismatched);
  signal_tree->Branch("mu0_ismatched", &the_sig_mu0_ismatched);
  signal_tree->Branch("mu_ismatched", &the_sig_mu_ismatched);
  signal_tree->Branch("pi_ismatched", &the_sig_pi_ismatched);
  signal_tree->Branch("mupi_mass_reco_gen_reldiff", &the_sig_mupi_mass_reco_gen_reldiff);
  signal_tree->Branch("lxy_reco_gen_reldiff", &the_sig_lxy_reco_gen_reldiff);

  if(isMC){
    signal_tree->Branch("weight_mu0_softid", &the_sig_weight_mu0_softid);
    signal_tree->Branch("weight_mu_looseid", &the_sig_weight_mu_looseid);

    signal_tree->Branch("weight_mu0_dxy_bs", &the_sig_weight_mu0_dxy_bs);
    signal_tree->Branch("weight_mu_dxy_bs", &the_sig_weight_mu_dxy_bs);
    signal_tree->Branch("weight_mu0_dxysig_bs", &the_sig_weight_mu0_dxysig_bs);
    signal_tree->Branch("weight_mu_dxysig_bs", &the_sig_weight_mu_dxysig_bs);

    //signal_tree->Branch("weight_hlt_A1", &the_sig_weight_hlt_A1);
    //signal_tree->Branch("weight_hlt_A1_6", &the_sig_weight_hlt_A1_6);
    //signal_tree->Branch("weight_hlt_HLT_Mu9_IP6_A1_6", &the_sig_weight_hlt_HLT_Mu9_IP6_A1_6);
    //signal_tree->Branch("weight_hlt_HLT_Mu9_IP6_A1_6_v2", &the_sig_weight_hlt_HLT_Mu9_IP6_A1_6_v2);
    //signal_tree->Branch("weight_hlt_A1_6_B1", &the_sig_weight_hlt_A1_6_B1);
    //signal_tree->Branch("weight_hlt_D1_pteta_v1", &the_sig_weight_hlt_D1_pteta_v1);
    //signal_tree->Branch("weight_hlt_D1_ptdxysig_v1", &the_sig_weight_hlt_D1_ptdxysig_v1);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu12_IP6_pteta_v1", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu12_IP6_pteta_v1);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6);
    //signal_tree->Branch("weight_hlt_fullBPark_pteta_v1", &the_sig_weight_hlt_fullBPark_pteta_v1);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_pteta", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_pteta);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_ptdxysig", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_ptdxysig);
    //signal_tree->Branch("weight_hlt_fullBPark_ptetadxysig_max5e6", &the_sig_weight_hlt_fullBPark_ptetadxysig_max5e6);
    //signal_tree->Branch("weight_hlt_D1_ptetadxysig_max5e6", &the_sig_weight_hlt_D1_ptetadxysig_max5e6);
    //signal_tree->Branch("weight_hlt_fullA_V06_tag_fired_HLT_Mu9_IP6_pteta_max3e6", &the_sig_weight_hlt_fullA_V06_tag_fired_HLT_Mu9_IP6_pteta_max3e6);
    //signal_tree->Branch("weight_hlt_fullA_tag_fired_anyBParkHLT_pteta_max3e6", &the_sig_weight_hlt_fullA_tag_fired_anyBParkHLT_pteta_max3e6);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_DST_DoubleMu1_pteta_max3e6", &the_sig_weight_hlt_D1_tag_fired_DST_DoubleMu1_pteta_max3e6);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max5e6_v2", &the_sig_weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max5e6_v2);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_anyBParkHLT_pteta_max5e6_v2", &the_sig_weight_hlt_fullBPark_tag_fired_anyBParkHLT_pteta_max5e6_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max3e6_v3", &the_sig_weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max3e6_v3);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2);
    
    signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma);

    signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma);

    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma);

    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma);


    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6_smalltable_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6_smalltable_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_bsrdst_newcat_max3e6_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_bsrdst_newcat_max3e6_smalltable_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_bsrdst_inicat_max3e6_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_bsrdst_inicat_max3e6_smalltable_v2);

    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6_smalltable_v2);
    //signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6_smalltable_v2);

    signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale_smalltable_v2);
    signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03_smalltable_v2);


    signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma);

    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma);
    //signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma);

    signal_tree->Branch("weight_pu_qcd_A", &the_sig_weight_pu_qcd_A);
    signal_tree->Branch("weight_pu_qcd_B", &the_sig_weight_pu_qcd_B);
    signal_tree->Branch("weight_pu_qcd_C", &the_sig_weight_pu_qcd_C);
    signal_tree->Branch("weight_pu_qcd_D", &the_sig_weight_pu_qcd_D);
    signal_tree->Branch("weight_pu_qcd_tot", &the_sig_weight_pu_qcd_tot);
    signal_tree->Branch("weight_pu_sig_A", &the_sig_weight_pu_sig_A);
    signal_tree->Branch("weight_pu_sig_B", &the_sig_weight_pu_sig_B);
    signal_tree->Branch("weight_pu_sig_C", &the_sig_weight_pu_sig_C);
    signal_tree->Branch("weight_pu_sig_D", &the_sig_weight_pu_sig_D);
    signal_tree->Branch("weight_pu_sig_tot", &the_sig_weight_pu_sig_tot);
  }

  if(isMC){
    signal_tree->Branch("gen_mu0_mu_lxy", &the_gen_mu0_mu_lxy);
    signal_tree->Branch("gen_mu0_mu_lxyz", &the_gen_mu0_mu_lxyz);
    signal_tree->Branch("gen_b_pt", &the_gen_b_pt);
    signal_tree->Branch("gen_b_eta", &the_gen_b_eta);
    signal_tree->Branch("gen_b_phi", &the_gen_b_phi);
    signal_tree->Branch("gen_b_mass", &the_gen_b_mass);
    signal_tree->Branch("gen_b_pdgid", &the_gen_b_pdgid);
    signal_tree->Branch("gen_hnl_ct", &the_gen_hnl_ct);
    signal_tree->Branch("gen_hnl_lxy", &the_gen_hnl_lxy);
    signal_tree->Branch("gen_hnl_lxyz", &the_gen_hnl_lxyz);
    signal_tree->Branch("gen_hnl_pt", &the_gen_hnl_pt);
    signal_tree->Branch("gen_hnl_eta", &the_gen_hnl_eta);
    signal_tree->Branch("gen_hnl_phi", &the_gen_hnl_phi);
    signal_tree->Branch("gen_hnl_mass", &the_gen_hnl_mass);
    signal_tree->Branch("gen_hnl_vx", &the_gen_hnl_vx);
    signal_tree->Branch("gen_hnl_vy", &the_gen_hnl_vy);
    signal_tree->Branch("gen_hnl_vz", &the_gen_hnl_vz);
    signal_tree->Branch("gen_mu0_pt", &the_gen_mu0_pt);
    signal_tree->Branch("gen_mu0_eta", &the_gen_mu0_eta);
    signal_tree->Branch("gen_mu0_phi", &the_gen_mu0_phi);
    signal_tree->Branch("gen_mu0_vx", &the_gen_mu0_vx);
    signal_tree->Branch("gen_mu0_vy", &the_gen_mu0_vy);
    signal_tree->Branch("gen_mu0_vz", &the_gen_mu0_vz);
    signal_tree->Branch("gen_mu_pt", &the_gen_mu_pt);
    signal_tree->Branch("gen_mu_eta", &the_gen_mu_eta);
    signal_tree->Branch("gen_mu_phi", &the_gen_mu_phi);
    signal_tree->Branch("gen_mu_vx", &the_gen_mu_vx);
    signal_tree->Branch("gen_mu_vy", &the_gen_mu_vy);
    signal_tree->Branch("gen_mu_vz", &the_gen_mu_vz);
    signal_tree->Branch("gen_pi_pt", &the_gen_pi_pt);
    signal_tree->Branch("gen_pi_eta", &the_gen_pi_eta);
    signal_tree->Branch("gen_pi_phi", &the_gen_pi_phi);
    signal_tree->Branch("gen_pi_vx", &the_gen_pi_vx);
    signal_tree->Branch("gen_pi_vy", &the_gen_pi_vy);
    signal_tree->Branch("gen_pi_vz", &the_gen_pi_vz);
  }

  signal_tree->Branch("pv_npvs", &the_pv_npvs);

  //signal_tree->Branch("hlt_mu7_ip4", &the_hlt_mu7_ip4);
  //signal_tree->Branch("hlt_mu8_ip6", &the_hlt_mu8_ip6);
  //signal_tree->Branch("hlt_mu8_ip5", &the_hlt_mu8_ip5);
  //signal_tree->Branch("hlt_mu8_ip3", &the_hlt_mu8_ip3);
  //signal_tree->Branch("hlt_mu8p5_ip3p5", &the_hlt_mu8p5_ip3p5);
  //signal_tree->Branch("hlt_mu9_ip6", &the_hlt_mu9_ip6);
  //signal_tree->Branch("hlt_mu9_ip5", &the_hlt_mu9_ip5);
  //signal_tree->Branch("hlt_mu9_ip4", &the_hlt_mu9_ip4);
  //signal_tree->Branch("hlt_mu10p5_ip3p5", &the_hlt_mu10p5_ip3p5);
  //signal_tree->Branch("hlt_mu12_ip6", &the_hlt_mu12_ip6);


  //// generate toys for smearing
  ////TF1* gauss_function = new TF1("gauss_function", "gauss",-10, 10);
  //gauss_function = new TF1("gauss_function", "gaus",-10, 10);
  //gauss_function->SetParameters(1, 1.1, 0.1);
  ////gauss_function->SetParameters(1,0,1);
  
  //TODO add flag isMC here and for the random choice
  gauss_function_0p03 = new TF1("gauss_function_0p03", "gaus",-10, 10);
  gauss_function_0p03->SetParameters(1, 1., 0.03);


  // defining histograms
  if(do_fillhistograms){
    my_file->mkdir("signal_channel");
    my_file->cd("signal_channel");

    sighist_ncand_perevent = new TH1F("sighist_ncand_perevent", "sighist_ncand_perevent", 10, 0, 10);
    sighist_ncand_matched_perevent = new TH1F("sighist_ncand_matched_perevent", "sighist_ncand_matched__perevent", 10, 0, 10);

    sighist_selection_efficiency_hnlpt_allevents = new TH1F("sighist_selection_efficiency_hnlpt_allevents", "Efficiency of selection of candidate with largest hnl pT (all events)", 2, 0, 2);
    sighist_selection_efficiency_hnlpt_eventswithmultcands = new TH1F("sighist_selection_efficiency_hnlpt_eventswithmultcands", "Efficiency of selection of candidate with largest hnl pT (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_bpt_allevents = new TH1F("sighist_selection_efficiency_bpt_allevents", "Efficiency of selection of candidate with largest b pT (all events)", 2, 0, 2);
    sighist_selection_efficiency_bpt_eventswithmultcands = new TH1F("sighist_selection_efficiency_bpt_eventswithmultcands", "Efficiency of selection of candidate with largest b pT (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_mu0pt_allevents = new TH1F("sighist_selection_efficiency_mu0pt_allevents", "Efficiency of selection of candidate with largest trigger muon pT (all events)", 2, 0, 2);
    sighist_selection_efficiency_mu0pt_eventswithmultcands = new TH1F("sighist_selection_efficiency_mu0pt_eventswithmultcands", "Efficiency of selection of candidate with largest trigger muon pT (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_pipt_allevents = new TH1F("sighist_selection_efficiency_pipt_allevents", "Efficiency of selection of candidate with largest pi pT (all events)", 2, 0, 2);
    sighist_selection_efficiency_pipt_eventswithmultcands = new TH1F("sighist_selection_efficiency_pipt_eventswithmultcands", "Efficiency of selection of candidate with largest pi pT (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_svprob_allevents = new TH1F("sighist_selection_efficiency_svprob_allevents", "Efficiency of selection of candidate with largest vertex prob. (all events)", 2, 0, 2);
    sighist_selection_efficiency_svprob_eventswithmultcands = new TH1F("sighist_selection_efficiency_svprob_eventswithmultcands", "Efficiency of selection of candidate with largest vertex prob. (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_svchi2_allevents = new TH1F("sighist_selection_efficiency_svchi2_allevents", "Efficiency of selection of candidate with largest vertex chi2 (all events)", 2, 0, 2);
    sighist_selection_efficiency_svchi2_eventswithmultcands = new TH1F("sighist_selection_efficiency_svchi2_eventswithmultcands", "Efficiency of selection of candidate with largest vertex chi2 (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_cos2d_allevents = new TH1F("sighist_selection_efficiency_cos2d_allevents", "Efficiency of selection of candidate with largest vertex cos(angle) (all events)", 2, 0, 2);
    sighist_selection_efficiency_cos2d_eventswithmultcands = new TH1F("sighist_selection_efficiency_cos2d_eventswithmultcands", "Efficiency of selection of candidate with largest vertex cos(angle) (events with multiple candidates)", 2, 0, 2);
    //sighist_selection_efficiency_hnliso4_allevents = new TH1F("sighist_selection_efficiency_hnliso4_allevents", "Efficiency of selection of candidate with smallest hnl isolation (all events)", 2, 0, 2);
    //sighist_selection_efficiency_hnliso4_eventswithmultcands = new TH1F("sighist_selection_efficiency_hnliso4_eventswithmultcands", "Efficiency of selection of candidate with smallest hnl isolation (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_dr_allevents = new TH1F("sighist_selection_efficiency_dr_allevents", "Efficiency of selection of candidate with smallest dR(mu0, hnl) (all events)", 2, 0, 2);
    sighist_selection_efficiency_dr_eventswithmultcands = new TH1F("sighist_selection_efficiency_dr_eventswithmultcands", "Efficiency of selection of candidate with smallest dR(mu0, hnl) (events with multiple candidates)", 2, 0, 2);

    my_file->cd();
  } // end define histograms
}


Bool_t BToMuMuPiDumper::Process(Long64_t entry)
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
  //cout << endl << "--- Entry " << entry << " ---" << endl;

  // for data, we skip the event in case it doesn't pass the lumi mask
  if(!isMC && lumiMask(*run, *luminosityBlock) == false) return false;

  // number of candidates in the event
  UInt_t nCand_sig = *nBToMuMuPi; 

  the_event = *event; 
  the_run = *run;
  the_lumi = *luminosityBlock; 

  the_pv_npvs = *PV_npvs;

  //the_hlt_mu7_ip4 = *HLT_Mu7_IP4;
  //the_hlt_mu8_ip6 = *HLT_Mu8_IP6;
  //the_hlt_mu8_ip5 = *HLT_Mu8_IP5;
  //the_hlt_mu8_ip3 = *HLT_Mu8_IP3;
  //the_hlt_mu8p5_ip3p5 = *HLT_Mu8p5_IP3p5;
  //the_hlt_mu9_ip6 = *HLT_Mu9_IP6;
  //the_hlt_mu9_ip5 = *HLT_Mu9_IP5;
  //the_hlt_mu9_ip4 = *HLT_Mu9_IP4;
  //the_hlt_mu10p5_ip3p5 = *HLT_Mu10p5_IP3p5;
  //the_hlt_mu12_ip6 = *HLT_Mu12_IP6;

  //   ----- Signal Channel -----  //

  if(nCand_sig > 0){ // at least one candidate per event
    // selecting the candidate as the one having the largest hnl pt
    // - create candIdx - cos2d pairs
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_cos2D);
    //std::cout << "initial vector" << std::endl;
    //for(unsigned int i(0); i<pair_candIdx_desc_cos2d_sig.size(); ++i){
    //  std::cout << pair_candIdx_desc_cos2d_sig[i].first << " " << pair_candIdx_desc_cos2d_sig[i].second << " " << BToMuMuPi_hnl_charge[pair_candIdx_desc_cos2d_sig[i].first] << std::endl;
    //}
    // - sort it in decreasing cos2d
    stable_sort(pair_candIdx_desc_cos2d_sig.begin(), pair_candIdx_desc_cos2d_sig.end(), sortcansbydesc);
    //std::cout << std::endl << "after sort in cos2d" << std::endl;
    //for(unsigned int i(0); i<pair_candIdx_desc_cos2d_sig.size(); ++i){
    //  std::cout << pair_candIdx_desc_cos2d_sig[i].first << " " << pair_candIdx_desc_cos2d_sig[i].second << " " << BToMuMuPi_hnl_charge[pair_candIdx_desc_cos2d_sig[i].first] << std::endl;
    //}

    // - then privilege OS cand over SS ones
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sig, BToMuMuPi_hnl_charge);
    stable_sort(pair_candIdx_desc_cos2d_sign_sig.begin(), pair_candIdx_desc_cos2d_sign_sig.end(), sortcansbydesc_opp);
    //std::cout << std::endl << "after sort in sign" << std::endl;
    //for(unsigned int i(0); i<pair_candIdx_desc_cos2d_sign_sig.size(); ++i){
    //  std::cout << pair_candIdx_desc_cos2d_sign_sig[i].first << " " << BToMuMuPi_hnl_cos2D[pair_candIdx_desc_cos2d_sign_sig[i].first] << " " << pair_candIdx_desc_cos2d_sign_sig[i].second << std::endl;
    //}

    // - for signal, priviledge matched candidates
    //vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_matched_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sign_sig, BToMuMuPi_isMatched);
    //stable_sort(pair_candIdx_desc_cos2d_sign_matched_sig.begin(), pair_candIdx_desc_cos2d_sign_matched_sig.end(), sortcansbydesc);

    // - priviledge slimmed over dsa candidates
    //vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_matched_muon_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sign_sig, BToMuMuPi_mu_idx, Muon_isDSAMuon);
    //stable_sort(pair_candIdx_desc_cos2d_sign_matched_muon_sig.begin(), pair_candIdx_desc_cos2d_sign_matched_muon_sig.end(), sortcansbydesc_opp);

    // - and select the candidate
    UInt_t selectedCandIdx_sig = pair_candIdx_desc_cos2d_sign_sig[0].first;
    //std::cout << std::endl << "selected candidate: " << selectedCandIdx_sig << std::endl;

    // for signal MC, only keep the matched events
    if(isSignalMC && BToMuMuPi_isMatched[selectedCandIdx_sig] != 1) return false;

    // keep events in the signal region only //TODO add flag for signal region and control region?
    //if(BToMuMuPi_hnl_charge[selectedCandIdx_sig] != 0) return false;

    // mc corrections
    //float corr_dxy_bs = 1.; //isMC ? 1.21 : 1.; 
    //float corr_dxysig_bs = 1.; //isMC ? 1.10 : 1.; 

    //// smearing
    ////gauss_function->SetParameters(1, corr_dxysig_bs, 0.1);
    //float smearing_dxysig_bs = gauss_function->GetRandom(); 
    ////std::cout << "smearing " << smearing_dxysig_bs << std::endl;
    float smeared_corr_0p03 = gauss_function_0p03->GetRandom(); 

    // fill the signal_tree
    if(BToMuMuPi_mu0_pt[selectedCandIdx_sig] == Muon_pt[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]){ // temporary condition, skip events with faulty indexing
      the_sig_b_pt = BToMuMuPi_pt[selectedCandIdx_sig];
      the_sig_b_eta = BToMuMuPi_eta[selectedCandIdx_sig];
      the_sig_b_phi = BToMuMuPi_phi[selectedCandIdx_sig];
      the_sig_b_mass = BToMuMuPi_mass[selectedCandIdx_sig];
      the_sig_b_charge = BToMuMuPi_charge[selectedCandIdx_sig];
      the_sig_b_pdgid = BToMuMuPi_pdgId[selectedCandIdx_sig];

      the_sig_hnl_pt = BToMuMuPi_hnl_pt[selectedCandIdx_sig];
      the_sig_hnl_eta = BToMuMuPi_hnl_eta[selectedCandIdx_sig];
      the_sig_hnl_phi = BToMuMuPi_hnl_phi[selectedCandIdx_sig];
      the_sig_hnl_mass = BToMuMuPi_hnl_mass[selectedCandIdx_sig];
      the_sig_hnl_charge = BToMuMuPi_hnl_charge[selectedCandIdx_sig];
      the_sig_hnl_ct = BToMuMuPi_hnl_ct[selectedCandIdx_sig];
      the_sig_hnl_cos2d = BToMuMuPi_hnl_cos2D[selectedCandIdx_sig];
      //the_sig_hnl_iso03 = BToMuMuPi_hnl_iso03[selectedCandIdx_sig];
      //the_sig_hnl_iso03_close = BToMuMuPi_hnl_iso03_close[selectedCandIdx_sig];
      //the_sig_hnl_iso03_rel_close = BToMuMuPi_hnl_iso03_rel_close[selectedCandIdx_sig];
      //the_sig_hnl_iso04 = BToMuMuPi_hnl_iso04[selectedCandIdx_sig];
      //the_sig_hnl_iso04_close = BToMuMuPi_hnl_iso04_close[selectedCandIdx_sig];
      //the_sig_hnl_iso04_rel_close = BToMuMuPi_hnl_iso04_rel_close[selectedCandIdx_sig];

      the_sig_mu0_pt = BToMuMuPi_mu0_pt[selectedCandIdx_sig];
      //the_sig_mu0_pt = Muon_pt[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_eta = BToMuMuPi_mu0_eta[selectedCandIdx_sig];
      //the_sig_mu0_eta = Muon_eta[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_phi = BToMuMuPi_mu0_phi[selectedCandIdx_sig];
      //the_sig_mu0_phi = Muon_phi[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_charge = Muon_charge[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_dxy = fabs(Muon_dxy[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]);
      //the_sig_mu0_dxy_bs = corr_dxy_bs * Muon_dxy_BS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_dxy_bs = fabs(Muon_dxy_BS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]);
      the_sig_mu0_dxysig = fabs(Muon_dxyS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]);
      the_sig_mu0_dxysig_bs =  fabs(Muon_dxyS_BS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]);
      the_sig_mu0_dxysig_bs_linearscale =  isMC ? (1.178 - 0.002968 * fabs(Muon_dxyS_BS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]])) * fabs(Muon_dxyS_BS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]) : fabs(Muon_dxyS_BS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]);
      the_sig_mu0_dxysig_bs_scale1p12_smear0p03 = isMC ? 1.12 * smeared_corr_0p03 * fabs(Muon_dxyS_BS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]) : fabs(Muon_dxyS_BS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]); 
      //the_sig_mu0_dxysig_bs_smeared = (corr_dxysig_bs + smearing_dxysig_bs) * Muon_dxyS_BS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_dxysig_bs_rdst = fabs(Muon_dxyS_BS_alaRdst[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]);
      the_sig_mu0_dz = Muon_dz[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_dzsig = Muon_dzS[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_triggermatching_dr = Muon_matched_dr[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_triggermatching_dpt = Muon_matched_dpt[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_ip3d = Muon_ip3d[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_ip3dsig = Muon_sip3d[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_pfiso03 = Muon_pfiso03_all[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_pfiso03_rel = Muon_pfiso03Rel_all[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_iso03 = BToMuMuPi_mu0_iso03[selectedCandIdx_sig];
      //the_sig_mu0_iso03_close = BToMuMuPi_mu0_iso03_close[selectedCandIdx_sig];
      //the_sig_mu0_iso03_rel_close = BToMuMuPi_mu0_iso03_rel_close[selectedCandIdx_sig];
      the_sig_mu0_pfiso04 = Muon_pfiso04_all[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_pfiso04_rel = Muon_pfiso04Rel_all[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_iso04 = BToMuMuPi_mu0_iso04[selectedCandIdx_sig];
      //the_sig_mu0_iso04_close = BToMuMuPi_mu0_iso04_close[selectedCandIdx_sig];
      //the_sig_mu0_iso04_rel_close = BToMuMuPi_mu0_iso04_rel_close[selectedCandIdx_sig];
      the_sig_mu0_looseid = Muon_looseId[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_mediumid = Muon_mediumId[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_tightid = Muon_tightId[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_softid = Muon_softId[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_pfisoid = Muon_pfIsoId[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_trkisoid = Muon_tkIsoId[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_triggerlooseid = Muon_triggerIdLoose[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_istriggering = Muon_isTriggeringBPark[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_isPF = Muon_isPF[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_isglobalmuon = Muon_isGlobalMuon[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_istrackermuon = Muon_isTrackerMuon[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_isglobalortrackermuon = Muon_isGlobalOrTrackerMuon[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_isglobalnottrackermuon = Muon_isGlobalNotTrackerMuon[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_istrackernotglobalmuon = Muon_isTrackerNotGlobalMuon[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      the_sig_mu0_intimemuon = Muon_inTimeMuon[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_segmentcompatibility = Muon_segmentCompatibility[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_calocompatibility = Muon_caloCompatibility[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_validhitfraction = Muon_validHitFraction[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_kinkfinderchi2 = Muon_kinkFinderChi2[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_globalnormalisedchi2 = Muon_globalNormalisedChi2[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_localpositionchi2 = Muon_localPositionChi2[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_trackerhighpurityflag = Muon_trackerHighPurityFlag[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_numberofvalidmuonhits = Muon_numberOfValidMuonHits[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_numberofvalidpixelhits = Muon_numberOfValidPixelHits[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_numberoftrackerlayers = Muon_numberOfTrackerLayers[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_numberofpixellayers = Muon_numberOfPixelLayers[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_numberofstations = Muon_numberOfStations[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_fired_hlt_mu7_ip4 = Muon_fired_HLT_Mu7_IP4[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_fired_hlt_mu8_ip3 = Muon_fired_HLT_Mu8_IP3[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_fired_hlt_mu8_ip5 = Muon_fired_HLT_Mu8_IP5[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_fired_hlt_mu8_ip6 = Muon_fired_HLT_Mu8_IP6[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_fired_hlt_mu8p5_ip3p5 = Muon_fired_HLT_Mu8p5_IP3p5[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_fired_hlt_mu9_ip4 = Muon_fired_HLT_Mu9_IP4[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_fired_hlt_mu9_ip5 = Muon_fired_HLT_Mu9_IP5[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_fired_hlt_mu9_ip6 = Muon_fired_HLT_Mu9_IP6[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_fired_hlt_mu10p5_ip3p5 = Muon_fired_HLT_Mu10p5_IP3p5[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_fired_hlt_mu12_ip6 = Muon_fired_HLT_Mu12_IP6[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_prescale_hlt_mu7_ip4 = Muon_prescale_HLT_Mu7_IP4[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_prescale_hlt_mu8_ip3 = Muon_prescale_HLT_Mu8_IP3[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_prescale_hlt_mu8_ip5 = Muon_prescale_HLT_Mu8_IP5[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_prescale_hlt_mu8_ip6 = Muon_prescale_HLT_Mu8_IP6[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_prescale_hlt_mu8p5_ip3p5 = Muon_prescale_HLT_Mu8p5_IP3p5[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_prescale_hlt_mu9_ip4 = Muon_prescale_HLT_Mu9_IP4[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_prescale_hlt_mu9_ip5 = Muon_prescale_HLT_Mu9_IP5[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_prescale_hlt_mu9_ip6 = Muon_prescale_HLT_Mu9_IP6[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_prescale_hlt_mu10p5_ip3p5 = Muon_prescale_HLT_Mu10p5_IP3p5[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];
      //the_sig_mu0_prescale_hlt_mu12_ip6 = Muon_prescale_HLT_Mu12_IP6[BToMuMuPi_mu0_idx[selectedCandIdx_sig]];

      the_sig_mu0_triggering_pt = -99.;
      the_sig_mu0_triggering_eta = -99.;
      the_sig_mu0_triggering_phi = -99.;
      the_sig_mu0_triggering_dxy = -99.;
      the_sig_mu0_triggering_dxy_bs = -99.;
      the_sig_mu0_triggering_dxysig = -99.;
      the_sig_mu0_triggering_dxysig_bs = -99.;
      the_sig_mu0_triggering_dxysig_bs_rdst = -99.;
      the_sig_mu0_triggering_dxysig_bs_linearscale = -99.;
      the_sig_mu0_triggering_dxysig_bs_scale1p12_smear0p03 = -99.;
      if(the_sig_mu0_istriggering == 1){
        the_sig_mu0_triggering_pt = the_sig_mu0_pt;
        the_sig_mu0_triggering_eta = the_sig_mu0_eta;
        the_sig_mu0_triggering_phi = the_sig_mu0_phi;
        the_sig_mu0_triggering_dxy = the_sig_mu0_dxy;
        the_sig_mu0_triggering_dxy_bs = the_sig_mu0_dxy_bs;
        the_sig_mu0_triggering_dxysig = the_sig_mu0_dxysig;
        //the_sig_mu0_triggering_dxysig_bs = corr_dxysig_bs * the_sig_mu0_dxysig_bs;
        the_sig_mu0_triggering_dxysig_bs = the_sig_mu0_dxysig_bs; // corrected above
        the_sig_mu0_triggering_dxysig_bs_linearscale = the_sig_mu0_dxysig_bs_linearscale; // corrected above
        the_sig_mu0_triggering_dxysig_bs_scale1p12_smear0p03 = the_sig_mu0_dxysig_bs_scale1p12_smear0p03; // corrected above
        the_sig_mu0_triggering_dxysig_bs_rdst = the_sig_mu0_dxysig_bs_rdst;
      }

      the_sig_mu_pt = BToMuMuPi_fit_mu_pt[selectedCandIdx_sig];
      the_sig_mu_eta = BToMuMuPi_fit_mu_eta[selectedCandIdx_sig];
      the_sig_mu_phi = BToMuMuPi_fit_mu_phi[selectedCandIdx_sig]; 
      the_sig_mu_charge = Muon_charge[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_dxy = fabs(Muon_dxy[BToMuMuPi_mu_idx[selectedCandIdx_sig]]);
      the_sig_mu_dxy_bs = fabs(Muon_dxy_BS[BToMuMuPi_mu_idx[selectedCandIdx_sig]]);
      the_sig_mu_dxysig = fabs(Muon_dxyS[BToMuMuPi_mu_idx[selectedCandIdx_sig]]);
      //the_sig_mu_dxysig_bs = corr_dxysig_bs * Muon_dxyS_BS[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_dxysig_bs =  fabs(Muon_dxyS_BS[BToMuMuPi_mu_idx[selectedCandIdx_sig]]);
      the_sig_mu_dxysig_bs_linearscale =  isMC ? (1.178 - 0.002968 * fabs(Muon_dxyS_BS[BToMuMuPi_mu_idx[selectedCandIdx_sig]])) * fabs(Muon_dxyS_BS[BToMuMuPi_mu_idx[selectedCandIdx_sig]]) : fabs(Muon_dxyS_BS[BToMuMuPi_mu_idx[selectedCandIdx_sig]]);
      the_sig_mu_dxysig_bs_scale1p12_smear0p03 = isMC ? 1.12 * smeared_corr_0p03 * fabs(Muon_dxyS_BS[BToMuMuPi_mu_idx[selectedCandIdx_sig]]) : fabs(Muon_dxyS_BS[BToMuMuPi_mu_idx[selectedCandIdx_sig]]); 
      the_sig_mu_dxysig_bs_rdst = fabs(Muon_dxyS_BS_alaRdst[BToMuMuPi_mu_idx[selectedCandIdx_sig]]);
      the_sig_mu_dz = Muon_dz[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_dzsig = Muon_dzS[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_triggermatching_dr = Muon_matched_dr[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_triggermatching_dpt = Muon_matched_dpt[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_ip3d = Muon_ip3d[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_ip3dsig = Muon_sip3d[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_pfiso03 = Muon_pfiso03_all[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_pfiso03_rel = Muon_pfiso03Rel_all[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_iso03 = BToMuMuPi_sel_mu_iso03[selectedCandIdx_sig];
      //the_sig_mu_iso03_close = BToMuMuPi_sel_mu_iso03_close[selectedCandIdx_sig];
      //the_sig_mu_iso03_rel_close = BToMuMuPi_sel_mu_iso03_rel_close[selectedCandIdx_sig];
      the_sig_mu_pfiso04 = Muon_pfiso04_all[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_pfiso04_rel = Muon_pfiso04Rel_all[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_iso04 = BToMuMuPi_sel_mu_iso04[selectedCandIdx_sig];
      //the_sig_mu_iso04_close = BToMuMuPi_sel_mu_iso04_close[selectedCandIdx_sig];
      //the_sig_mu_iso04_rel_close = BToMuMuPi_sel_mu_iso04_rel_close[selectedCandIdx_sig];
      //the_sig_mu_isloose = BToMuMuPi_sel_mu_isLoose[selectedCandIdx_sig];
      //the_sig_mu_ismedium = BToMuMuPi_sel_mu_isMedium[selectedCandIdx_sig];
      //the_sig_mu_istight = BToMuMuPi_sel_mu_isTight[selectedCandIdx_sig];
      //the_sig_mu_issoft = BToMuMuPi_sel_mu_isSoft[selectedCandIdx_sig];
      the_sig_mu_looseid = Muon_looseId[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_mediumid = Muon_mediumId[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_tightid = Muon_tightId[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_softid = Muon_softId[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_pfisoid = Muon_pfIsoId[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_trkisoid = Muon_tkIsoId[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_triggerlooseid = Muon_triggerIdLoose[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_istriggering = Muon_isTriggeringBPark[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_isslimmed = Muon_isSlimmedMuon[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_isdsa = Muon_isDSAMuon[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_isPF = Muon_isPF[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_isglobalmuon = Muon_isGlobalMuon[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_istrackermuon = Muon_isTrackerMuon[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_isglobalortrackermuon = Muon_isGlobalOrTrackerMuon[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_isglobalnottrackermuon = Muon_isGlobalNotTrackerMuon[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_istrackernotglobalmuon = Muon_isTrackerNotGlobalMuon[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_intimemuon = Muon_inTimeMuon[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_segmentcompatibility = Muon_segmentCompatibility[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_calocompatibility = Muon_caloCompatibility[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_validhitfraction = Muon_validHitFraction[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_kinkfinderchi2 = Muon_kinkFinderChi2[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_globalnormalisedchi2 = Muon_globalNormalisedChi2[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_localpositionchi2 = Muon_localPositionChi2[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_trackerhighpurityflag = Muon_trackerHighPurityFlag[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_numberofvalidmuonhits = Muon_numberOfValidMuonHits[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_numberofvalidpixelhits = Muon_numberOfValidPixelHits[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_numberoftrackerlayers = Muon_numberOfTrackerLayers[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_numberofpixellayers = Muon_numberOfPixelLayers[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_numberofstations = Muon_numberOfStations[BToMuMuPi_mu_idx[selectedCandIdx_sig]];
      // add customised muon ids
      Bool_t isgoodGlobalMuon = 0;
      if(the_sig_mu_isglobalmuon==1 && the_sig_mu_localpositionchi2<12 && the_sig_mu_kinkfinderchi2<20) isgoodGlobalMuon = 1;
      if(the_sig_mu_looseid==1 && ((isgoodGlobalMuon==1 && the_sig_mu_segmentcompatibility>0.303) || (isgoodGlobalMuon==0 && the_sig_mu_segmentcompatibility>0.451))){
        the_sig_mu_whnlid = 1;
      }
      else{
        the_sig_mu_whnlid = 0;
      }

      if(the_sig_mu_looseid==1 && the_sig_mu_intimemuon==1 && the_sig_mu_trackerhighpurityflag==1 && ((the_sig_mu_isglobalmuon==1 && the_sig_mu_numberofstations>0 && the_sig_mu_numberoftrackerlayers<18) || (the_sig_mu_isglobalmuon!=1 && the_sig_mu_calocompatibility>0.05 && the_sig_mu_numberoftrackerlayers>6 && the_sig_mu_numberoftrackerlayers<16 && the_sig_mu_numberofvalidpixelhits<6))){
        the_sig_mu_customisedid = 1;
      }
      else{
        the_sig_mu_customisedid = 0;
      }

      // displaced triggering muon
      the_sig_mu_triggering_pt = -99.;
      the_sig_mu_triggering_eta = -99.;
      the_sig_mu_triggering_phi = -99.;
      the_sig_mu_triggering_dxy = -99.;
      the_sig_mu_triggering_dxy_bs = -99.;
      the_sig_mu_triggering_dxysig = -99.;
      the_sig_mu_triggering_dxysig_bs = -99.;
      the_sig_mu_triggering_dxysig_bs_rdst = -99.;
      the_sig_mu_triggering_dxysig_bs_linearscale = -99.;
      the_sig_mu_triggering_dxysig_bs_scale1p12_smear0p03 = -99.;
      if(the_sig_mu_istriggering == 1){
        the_sig_mu_triggering_pt = the_sig_mu_pt;
        the_sig_mu_triggering_eta = the_sig_mu_eta;
        the_sig_mu_triggering_phi = the_sig_mu_phi;
        the_sig_mu_triggering_dxy = the_sig_mu_dxy;
        the_sig_mu_triggering_dxy_bs = the_sig_mu_dxy_bs;
        the_sig_mu_triggering_dxysig = the_sig_mu_dxysig;
        //the_sig_mu_triggering_dxysig_bs = corr_dxysig_bs * the_sig_mu_dxysig_bs;
        the_sig_mu_triggering_dxysig_bs = the_sig_mu_dxysig_bs; // corrected above
        the_sig_mu_triggering_dxysig_bs_linearscale = the_sig_mu_dxysig_bs_linearscale; // corrected above
        the_sig_mu_triggering_dxysig_bs_scale1p12_smear0p03 = the_sig_mu_dxysig_bs_scale1p12_smear0p03; // corrected above
        the_sig_mu_triggering_dxysig_bs_rdst = the_sig_mu_dxysig_bs_rdst;
      }

      the_sig_pi_pt = BToMuMuPi_fit_pi_pt[selectedCandIdx_sig];
      the_sig_pi_eta = BToMuMuPi_fit_pi_eta[selectedCandIdx_sig];
      the_sig_pi_phi = BToMuMuPi_fit_pi_phi[selectedCandIdx_sig]; 
      the_sig_pi_charge = BToMuMuPi_pi_charge[selectedCandIdx_sig];
      the_sig_pi_dcasig = BToMuMuPi_pi_DCASig[selectedCandIdx_sig];
      the_sig_pi_dxy = BToMuMuPi_pi_dxy[selectedCandIdx_sig];
      the_sig_pi_dz = BToMuMuPi_pi_dz[selectedCandIdx_sig];
      the_sig_pi_dxysig = BToMuMuPi_pi_dxyS[selectedCandIdx_sig];
      the_sig_pi_dzsig = BToMuMuPi_pi_dzS[selectedCandIdx_sig];
      the_sig_pi_ispacked = BToMuMuPi_pi_ispacked[selectedCandIdx_sig];
      the_sig_pi_islost = BToMuMuPi_pi_islost[selectedCandIdx_sig];
      //the_sig_pi_mu0_dr = ProbeTracks_drTrg[BToMuMuPi_pi_idx[selectedCandIdx_sig]];
      //the_sig_pi_ismatchedtomuon = ProbeTracks_isMatchedToMuon[BToMuMuPi_pi_idx[selectedCandIdx_sig]];
      //the_sig_pi_chi2 = BToMuMuPi_pi_chi2[selectedCandIdx_sig];
      //the_sig_pi_normalisedChi2 = BToMuMuPi_pi_normalisedChi2[selectedCandIdx_sig];
      //the_sig_pi_validFraction = BToMuMuPi_pi_validFraction[selectedCandIdx_sig];
      //the_sig_pi_ndof = BToMuMuPi_pi_ndof[selectedCandIdx_sig];
      //the_sig_pi_numberOfValidHits = BToMuMuPi_pi_numberOfValidHits[selectedCandIdx_sig];
      //the_sig_pi_numberOfLostHits = BToMuMuPi_pi_numberOfLostHits[selectedCandIdx_sig];
      //the_sig_pi_numberOfValidPixelHits = BToMuMuPi_pi_numberOfValidPixelHits[selectedCandIdx_sig];
      //the_sig_pi_numberOfTrackerLayers = BToMuMuPi_pi_numberOfTrackerLayers[selectedCandIdx_sig];
      //the_sig_pi_numberOfPixelLayers = BToMuMuPi_pi_numberOfPixelLayers[selectedCandIdx_sig];
      //the_sig_pi_qualityIndex = BToMuMuPi_pi_qualityIndex[selectedCandIdx_sig];
      the_sig_pi_highPurityFlag = BToMuMuPi_pi_highPurityFlag[selectedCandIdx_sig];
      if(the_sig_pi_ispacked && the_sig_pi_highPurityFlag){
        the_sig_pi_packedcandhashighpurity = 1;
      }
      else{
        the_sig_pi_packedcandhashighpurity = 0;
      }

      // track to muon matching
      //bool pi_matchedtomuon_loose = 0;
      //bool pi_matchedtomuon_medium = 0;
      //bool pi_matchedtomuon_tight = 0;
      //if(do_tracktomuon_matching){
      //  for(int iMuon(0); iMuon<fabs(*nMuon); ++iMuon){
      //    // do not consider muons in the signal final state
      //    if(iMuon == BToMuMuPi_mu0_idx[selectedCandIdx_sig] || iMuon == BToMuMuPi_mu_idx[selectedCandIdx_sig]) continue;
      //    // compute the deltaR and deltaPtRel between the given track and the muon (unfitted values)
      //    float deltaR_track_muon = reco::deltaR(Muon_eta[iMuon], Muon_phi[iMuon], BToMuMuPi_pi_eta[selectedCandIdx_sig], BToMuMuPi_pi_phi[selectedCandIdx_sig]);
      //    float deltaPtRel = fabs(Muon_pt[iMuon] - BToMuMuPi_pi_pt[selectedCandIdx_sig]) / BToMuMuPi_pi_pt[selectedCandIdx_sig];

      //    // match if any muon fulfills requirements
      //    if(deltaR_track_muon < 0.3 && deltaPtRel < 1.){
      //      pi_matchedtomuon_loose = 1;
      //    }
      //    if(deltaR_track_muon < 0.3 && deltaPtRel < 0.3){
      //      pi_matchedtomuon_medium = 1;
      //    }
      //    if(deltaR_track_muon < 0.1 && deltaPtRel < 0.3){
      //      pi_matchedtomuon_tight = 1;
      //    }
      //  }
      //}
      //the_sig_pi_matchedtomuon_loose = pi_matchedtomuon_loose;
      //the_sig_pi_matchedtomuon_medium = pi_matchedtomuon_medium;
      //the_sig_pi_matchedtomuon_tight = pi_matchedtomuon_tight;

      the_sig_mu0_mu_mass = BToMuMuPi_mu0_mu_mass[selectedCandIdx_sig];
      the_sig_mu0_mu_pt = BToMuMuPi_mu0_mu_pt[selectedCandIdx_sig];
      the_sig_mu0_pi_mass = BToMuMuPi_mu0_pi_mass[selectedCandIdx_sig];
      the_sig_mu0_pi_pt = BToMuMuPi_mu0_pi_pt[selectedCandIdx_sig];

      the_sig_dimu_lxy = BToMuMuPi_dimu_Lxy[selectedCandIdx_sig];
      the_sig_dimu_lxyz = BToMuMuPi_dimu_Lxyz[selectedCandIdx_sig];
      //the_sig_dimu_vxdiff = BToMuMuPi_dimu_vxdiff[selectedCandIdx_sig];
      //the_sig_dimu_vydiff = BToMuMuPi_dimu_vydiff[selectedCandIdx_sig];
      //the_sig_dimu_vzdiff = BToMuMuPi_dimu_vzdiff[selectedCandIdx_sig];

      //the_sig_cos_theta_star_pion = BToMuMuPi_cos_theta_star_pion[selectedCandIdx_sig];
      //the_sig_cos_theta_star_muon = BToMuMuPi_cos_theta_star_muon[selectedCandIdx_sig];
      //the_sig_cos_theta_star_sum = fabs(BToMuMuPi_cos_theta_star_pion[selectedCandIdx_sig] + BToMuMuPi_cos_theta_star_muon[selectedCandIdx_sig]);

      //the_sig_px_diff_hnl_daughters_lab = BToMuMuPi_px_diff_hnl_daughters_lab[selectedCandIdx_sig];
      //the_sig_py_diff_hnl_daughters_lab = BToMuMuPi_py_diff_hnl_daughters_lab[selectedCandIdx_sig];
      //the_sig_pz_diff_hnl_daughters_lab = BToMuMuPi_pz_diff_hnl_daughters_lab[selectedCandIdx_sig];
      //the_sig_energy_diff_prefithnl_daughters_lab = BToMuMuPi_energy_diff_prefithnl_daughters_lab[selectedCandIdx_sig];
      //the_sig_px_diff_prefithnl_daughters_lab = BToMuMuPi_px_diff_prefithnl_daughters_lab[selectedCandIdx_sig];
      //the_sig_py_diff_prefithnl_daughters_lab = BToMuMuPi_py_diff_prefithnl_daughters_lab[selectedCandIdx_sig];
      //the_sig_pz_diff_prefithnl_daughters_lab = BToMuMuPi_pz_diff_prefithnl_daughters_lab[selectedCandIdx_sig];

      the_sig_deltar_mu_pi = BToMuMuPi_dr_mu_pi[selectedCandIdx_sig];
      the_sig_deltar_mu0_hnl = BToMuMuPi_dr_mu0_hnl[selectedCandIdx_sig];
      the_sig_deltar_mu0_mu = BToMuMuPi_dr_mu0_mu[selectedCandIdx_sig]; 
      the_sig_deltar_mu0_pi = BToMuMuPi_dr_mu0_pi[selectedCandIdx_sig]; 
      //the_sig_deltaeta_mu_pi = BToMuMuPi_deta_mu_pi[selectedCandIdx_sig];
      //the_sig_deltaeta_mu0_hnl = BToMuMuPi_deta_mu0_hnl[selectedCandIdx_sig];
      //the_sig_deltaeta_mu0_mu = BToMuMuPi_deta_mu0_mu[selectedCandIdx_sig];
      //the_sig_deltaeta_mu0_pi = BToMuMuPi_deta_mu0_pi[selectedCandIdx_sig];
      //the_sig_deltaphi_mu_pi = BToMuMuPi_dphi_mu_pi[selectedCandIdx_sig];
      //the_sig_deltaphi_mu0_hnl = BToMuMuPi_dphi_mu0_hnl[selectedCandIdx_sig];
      //the_sig_deltaphi_mu0_mu = BToMuMuPi_dphi_mu0_mu[selectedCandIdx_sig];
      //the_sig_deltaphi_mu0_pi = BToMuMuPi_dphi_mu0_pi[selectedCandIdx_sig];

      //the_sig_deltae_pi_fit_pi = BToMuMuPi_de_pi_fit_pi[selectedCandIdx_sig];
      //the_sig_deltae_mu_fit_mu = BToMuMuPi_de_mu_fit_mu[selectedCandIdx_sig];
      //the_sig_deltae_hnl_fit_hnl = BToMuMuPi_de_hnl_fit_hnl[selectedCandIdx_sig];
      //the_sig_deltapt_pi_fit_pi = BToMuMuPi_dpt_pi_fit_pi[selectedCandIdx_sig];
      //the_sig_deltapt_mu_fit_mu = BToMuMuPi_dpt_mu_fit_mu[selectedCandIdx_sig];
      //the_sig_deltapt_hnl_fit_hnl = BToMuMuPi_dpt_hnl_fit_hnl[selectedCandIdx_sig];
      //the_sig_deltapx_pi_fit_pi = BToMuMuPi_dpx_pi_fit_pi[selectedCandIdx_sig];
      //the_sig_deltapx_mu_fit_mu = BToMuMuPi_dpx_mu_fit_mu[selectedCandIdx_sig];
      //the_sig_deltapx_hnl_fit_hnl = BToMuMuPi_dpx_hnl_fit_hnl[selectedCandIdx_sig];
      //the_sig_deltapy_pi_fit_pi = BToMuMuPi_dpy_pi_fit_pi[selectedCandIdx_sig];
      //the_sig_deltapy_mu_fit_mu = BToMuMuPi_dpy_mu_fit_mu[selectedCandIdx_sig];
      //the_sig_deltapy_hnl_fit_hnl = BToMuMuPi_dpy_hnl_fit_hnl[selectedCandIdx_sig];
      //the_sig_deltapz_pi_fit_pi = BToMuMuPi_dpz_pi_fit_pi[selectedCandIdx_sig];
      //the_sig_deltapz_mu_fit_mu = BToMuMuPi_dpz_mu_fit_mu[selectedCandIdx_sig];
      //the_sig_deltapz_hnl_fit_hnl = BToMuMuPi_dpz_hnl_fit_hnl[selectedCandIdx_sig];
      //the_sig_deltaeta_pi_fit_pi = BToMuMuPi_deta_pi_fit_pi[selectedCandIdx_sig];
      //the_sig_deltaeta_mu_fit_mu = BToMuMuPi_deta_mu_fit_mu[selectedCandIdx_sig];
      //the_sig_deltaeta_hnl_fit_hnl = BToMuMuPi_deta_hnl_fit_hnl[selectedCandIdx_sig];
      //the_sig_deltaphi_pi_fit_pi = BToMuMuPi_dphi_pi_fit_pi[selectedCandIdx_sig];
      //the_sig_deltaphi_mu_fit_mu = BToMuMuPi_dphi_mu_fit_mu[selectedCandIdx_sig];
      //the_sig_deltaphi_hnl_fit_hnl = BToMuMuPi_dphi_hnl_fit_hnl[selectedCandIdx_sig];
      ////float deltaphi_mu_pi = fabs(BToMuMuPi_fit_mu_phi[selectedCandIdx_sig] - BToMuMuPi_fit_pi_phi[selectedCandIdx_sig]);
      ////the_sig_deltaphi_mu_pi = deltaphi_mu_pi > M_PI ? deltaphi_mu_pi : deltaphi_mu_pi - 2 * M_PI;

      the_sig_sv_chi2 = BToMuMuPi_sv_chi2[selectedCandIdx_sig];
      the_sig_sv_lxy = BToMuMuPi_sv_lxy[selectedCandIdx_sig];
      the_sig_sv_lxysig = BToMuMuPi_sv_lxy_sig[selectedCandIdx_sig];
      the_sig_sv_lxyz = BToMuMuPi_sv_lxyz[selectedCandIdx_sig];
      the_sig_sv_prob = BToMuMuPi_sv_prob[selectedCandIdx_sig];
      //the_sig_sv_x = BToMuMuPi_sv_x[selectedCandIdx_sig];
      //the_sig_sv_y = BToMuMuPi_sv_y[selectedCandIdx_sig];
      //the_sig_sv_z = BToMuMuPi_sv_z[selectedCandIdx_sig];

      the_sig_ismatched = BToMuMuPi_isMatched[selectedCandIdx_sig];
      the_sig_mu0_ismatched = BToMuMuPi_mu0_isMatched[selectedCandIdx_sig];
      the_sig_mu_ismatched = BToMuMuPi_mu_isMatched[selectedCandIdx_sig];
      the_sig_pi_ismatched = BToMuMuPi_pi_isMatched[selectedCandIdx_sig];
      the_sig_mupi_mass_reco_gen_reldiff = BToMuMuPi_mupi_mass_reco_gen_reldiff[selectedCandIdx_sig];
      the_sig_lxy_reco_gen_reldiff = BToMuMuPi_lxy_reco_gen_reldiff[selectedCandIdx_sig];

      // additionnal displacement quantities
      //float dist_sv_pv_xy = sqrt((BToMuMuPi_sv_x[selectedCandIdx_sig] - *PV_x) * (BToMuMuPi_sv_x[selectedCandIdx_sig] - *PV_x) + (BToMuMuPi_sv_y[selectedCandIdx_sig] - *PV_y) * (BToMuMuPi_sv_y[selectedCandIdx_sig] - *PV_y));
      //float dist_sv_pv_xyz = sqrt((BToMuMuPi_sv_x[selectedCandIdx_sig] - *PV_x) * (BToMuMuPi_sv_x[selectedCandIdx_sig] - *PV_x) + (BToMuMuPi_sv_y[selectedCandIdx_sig] - *PV_y) * (BToMuMuPi_sv_y[selectedCandIdx_sig] - *PV_y) + (BToMuMuPi_sv_z[selectedCandIdx_sig] - *PV_z) * (BToMuMuPi_sv_z[selectedCandIdx_sig] - *PV_z));
      //the_sig_sv_pv_lxy = dist_sv_pv_xy;
      //the_sig_sv_pv_lxyz = dist_sv_pv_xyz;

      the_sig_pi_mu0_vzdiff = BToMuMuPi_pi_mu0_vzdiff[selectedCandIdx_sig];


      // getting the displacement at gen level
      if(isSignalMC){
        UInt_t nGen = *nGenPart;

        float hnl_vx(0.), hnl_vy(0.), hnl_vz(0.);
        float mother_vx(0.), mother_vy(0.), mother_vz(0.);
        float mu0_vx(0.), mu0_vy(0.), mu0_vz(0.);
        float mu_vx(0.), mu_vy(0.), mu_vz(0.);
        int mother_idx(-99), mu0_idx(-99), mu_idx(-99);

        if(BToMuMuPi_isMatched[selectedCandIdx_sig]==1){
          mu0_vx = GenPart_vx[Muon_genPartIdx[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]];
          mu0_vy = GenPart_vy[Muon_genPartIdx[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]];
          mu0_vz = GenPart_vz[Muon_genPartIdx[BToMuMuPi_mu0_idx[selectedCandIdx_sig]]];
          mu_vx = GenPart_vx[Muon_genPartIdx[BToMuMuPi_mu_idx[selectedCandIdx_sig]]];
          mu_vy = GenPart_vy[Muon_genPartIdx[BToMuMuPi_mu_idx[selectedCandIdx_sig]]];
          mu_vz = GenPart_vz[Muon_genPartIdx[BToMuMuPi_mu_idx[selectedCandIdx_sig]]];

          the_gen_mu0_mu_lxy = sqrt((mu0_vx - mu_vx) * (mu0_vx - mu_vx) + (mu0_vy - mu_vy) * (mu0_vy - mu_vy));
          the_gen_mu0_mu_lxyz = sqrt((mu0_vx - mu_vx) * (mu0_vx - mu_vx) + (mu0_vy - mu_vy) * (mu0_vy - mu_vy) + (mu0_vz - mu_vz) * (mu0_vz - mu_vz));
        }
        else{
          the_gen_mu0_mu_lxy = -99.;
          the_gen_mu0_mu_lxyz = -99;
        }

        // gen information (no matching)

        // find idx of gen particles of interest
        int gen_hnl_idx(-99), gen_b_idx(-99), gen_mu0_idx(-99), gen_mu_idx(-99), gen_pi_idx(-99);

        for(unsigned int iGen(0); iGen < nGen; ++iGen){
          if(abs(GenPart_pdgId[iGen])==9900015){
            gen_hnl_idx = iGen;
            gen_b_idx = GenPart_genPartIdxMother[iGen]; 
            break;
          }
        }
        for(unsigned int iGen(0); iGen < nGen; ++iGen){
          if(abs(GenPart_pdgId[iGen])==13 && GenPart_genPartIdxMother[iGen]==gen_b_idx)
            gen_mu0_idx = iGen;
          if((abs(GenPart_pdgId[iGen])==13) && GenPart_genPartIdxMother[iGen]==gen_hnl_idx)
            gen_mu_idx = iGen;
          if(abs(GenPart_pdgId[iGen])==211 && GenPart_genPartIdxMother[iGen]==gen_hnl_idx)
            gen_pi_idx = iGen;
        }

        // get quantities that need more than one object
        if(gen_hnl_idx!=-99 && gen_mu_idx!=-99){

          ROOT::Math::PtEtaPhiMVector hnl_p4(GenPart_pt[gen_hnl_idx], GenPart_eta[gen_hnl_idx], GenPart_phi[gen_hnl_idx], GenPart_mass[gen_hnl_idx]);
          Float_t hnl_betagamma = hnl_p4.Beta() * hnl_p4.Gamma();

          the_gen_hnl_lxyz = get3Ddisp(GenPart_vx[gen_hnl_idx], GenPart_vx[gen_mu_idx],
              GenPart_vy[gen_hnl_idx], GenPart_vy[gen_mu_idx],
              GenPart_vz[gen_hnl_idx], GenPart_vz[gen_mu_idx]);

          the_gen_hnl_lxy  = get2Ddisp(GenPart_vx[gen_hnl_idx], GenPart_vx[gen_mu_idx],
              GenPart_vy[gen_hnl_idx], GenPart_vy[gen_mu_idx]);

          the_gen_hnl_ct = the_gen_hnl_lxyz / hnl_betagamma * 10.; // factor 10 to convert cm to mm
          //std::cout << "HNL pt,eta,phi,m"<< GenPart_pt[gen_hnl_idx] << " " << GenPart_eta[gen_hnl_idx] << " " << GenPart_phi[gen_hnl_idx] << " " << GenPart_mass[gen_hnl_idx] << std::endl;
          //std::cout << "HNL beta,gamma=" <<  hnl_p4.Beta() << " " << hnl_p4.Gamma() << std::endl;
          //std::cout << "HNL Lxy,Lxyz  =" <<  the_gen_hnl_lxy << " " << the_gen_hnl_lxyz << std::endl;
        }

        // set quantities for each object
        if(gen_b_idx!=-99){
          the_gen_b_pt = GenPart_pt[gen_b_idx];
          the_gen_b_eta = GenPart_eta[gen_b_idx];
          the_gen_b_phi = GenPart_phi[gen_b_idx];
          the_gen_b_mass = GenPart_mass[gen_b_idx];
          the_gen_b_pdgid = GenPart_pdgId[gen_b_idx];
        }
        if(gen_hnl_idx!=-99){
          the_gen_hnl_pt = GenPart_pt[gen_hnl_idx];
          the_gen_hnl_eta = GenPart_eta[gen_hnl_idx];
          the_gen_hnl_phi = GenPart_phi[gen_hnl_idx];
          the_gen_hnl_mass = GenPart_mass[gen_hnl_idx];
          the_gen_hnl_vx = GenPart_vx[gen_hnl_idx];
          the_gen_hnl_vy = GenPart_vy[gen_hnl_idx];
          the_gen_hnl_vz = GenPart_vz[gen_hnl_idx];
        }
      }

      if(isMC){
        // trigger scale factor
        //the_sig_weight_hlt_A1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_DST_DoubleMu1_A1_v1/scaleFactor_results_cat_pt_eta_fit.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_A1_6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_DST_DoubleMu1_A1_6_v1/scaleFactor_results_cat_pt_eta_fit.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_HLT_Mu9_IP6_A1_6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_HLT_Mu9_IP6_A1_6/scaleFactor_results_cat_pt_eta_fit.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_HLT_Mu9_IP6_A1_6_v2 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullA_fired_HLT_Mu9_IP6_max5e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;

        //the_sig_weight_hlt_A1_6_B1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_DST_DoubleMu1_A1_6_B1_v1/scaleFactor_results_cat_pt_eta_fit.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;

        //the_sig_weight_hlt_D1_pteta_v1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_anyBParkHLT_pteta_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_D1_ptdxysig_v1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_anyBParkHLT_ptdxysig_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)) : 1.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu12_IP6_pteta_v1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu12_IP6_pteta_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_pteta_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;

        //the_sig_weight_hlt_fullBPark_pteta_v1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_anyBParkHLT_max5e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_pteta = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_pteta_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_ptdxysig = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_ptdxysig_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)) : 1.;

        //the_sig_weight_hlt_fullBPark_ptetadxysig_max5e6 = isMC ? getTriggerScaleFactor_full(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_D1_ptetadxysig_max5e6 = isMC ? getTriggerScaleFactor_D1(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) : 1.;

        //the_sig_weight_hlt_fullA_V06_tag_fired_HLT_Mu9_IP6_pteta_max3e6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullA_V06_tag_fired_HLT_Mu9_IP6_pteta_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_fullA_tag_fired_anyBParkHLT_pteta_max3e6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullA_tag_fired_anyBParkHLT_pteta_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;

        //the_sig_weight_hlt_D1_tag_fired_DST_DoubleMu1_pteta_max3e6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_DST_DoubleMu1_pteta_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;

        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max5e6_v2 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_anyBParkHLT_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_fullBPark_tag_fired_anyBParkHLT_pteta_max5e6_v2 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_anyBParkHLT_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max3e6_v3 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_anyBParkHLT_pteta_max3e6_v3/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) : 1.;

        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2 = isMC ? getTriggerScaleFactor_D1(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) : 1.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2 = isMC ? getTriggerScaleFactor_fullBPark(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) : 1.;

        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = -99.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = -99.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = -99.;
        //if(isMC){
        //  if(the_sig_mu0_istriggering && the_sig_mu_istriggering){
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = getTriggerScaleFactor_D1(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = getTriggerScaleFactor_fullBPark(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        //  }
        //  else if(the_sig_mu0_istriggering && !the_sig_mu_istriggering){
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) * (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta))); 
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) * (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta))); 
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = getTriggerScaleFactor_D1(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * (1-getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta)));
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = getTriggerScaleFactor_fullBPark(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * (1-getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta)));
        //  }
        //  else if(!the_sig_mu0_istriggering && the_sig_mu_istriggering){
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta))) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta))) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = (1-getTriggerScaleFactor_D1(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta))) * getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = (1-getTriggerScaleFactor_fullBPark(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta))) * getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        //  }
        //  else if(!the_sig_mu0_istriggering && !the_sig_mu_istriggering){
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta))) * (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta))); 
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta))) * (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta))); 
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = (1-getTriggerScaleFactor_D1(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta))) * (1-getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta)));
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = (1-getTriggerScaleFactor_fullBPark(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta))) * (1-getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta)));
        //  }
        ////}

        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = -99.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = -99.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = -99.;
        //  if(the_sig_mu0_istriggering && the_sig_mu_istriggering){
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_fullBPark(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        //  }
        //  else if(the_sig_mu0_istriggering && !the_sig_mu_istriggering){
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)); 
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta));
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_eta)); 
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_fullBPark(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta));
        //  }
        //  else if(!the_sig_mu0_istriggering && the_sig_mu_istriggering){
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        //  }
        //  else if(!the_sig_mu0_istriggering && !the_sig_mu_istriggering){
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = 1.;
        //    the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = 1.;
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = 1.;
        //    the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = 1.;
        //  }
        //}

        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma = -99.;

        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2 = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma = -99.;

        the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2 = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma = -99.;

        the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2 = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma = -99.;

        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6_smalltable_v2 = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6_smalltable_v2 = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_bsrdst_newcat_max3e6_smalltable_v2 = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_bsrdst_inicat_max3e6_smalltable_v2 = -99.;
        
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6_smalltable_v2 = -99.;
        //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6_smalltable_v2 = -99.;

        the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale_smalltable_v2 = -99;
        the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03_smalltable_v2 = -99;

        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = -99.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma = -99.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma = -99.;

        the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2 = -99.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma = -99.;
        //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma = -99.;

        if(the_sig_mu0_istriggering && the_sig_mu_istriggering){
          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_plus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_plus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_minus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_minus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 

          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors_plus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors_plus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors_minus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors_minus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs)); 

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor_D1_plus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor_D1_minus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1_bs(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_D1_bs(the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor_D1_bs_plus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_D1_bs(the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor_D1_bs_minus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_D1_bs(the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs), fabs(the_sig_mu_eta));

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs_rdst)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs_rdst)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs_rdst)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs_rdst)); 

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs_rdst)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs_rdst)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs_rdst)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs_rdst)); 
          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale/scale_factors.root", the_sig_mu0_pt, the_sig_mu0_dxysig_bs_linearscale) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale/scale_factors.root", the_sig_mu_pt, the_sig_mu_dxysig_bs_linearscale); 
          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03/scale_factors.root", the_sig_mu0_pt, the_sig_mu0_dxysig_bs_scale1p12_smear0p03) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03/scale_factors.root", the_sig_mu_pt, the_sig_mu_dxysig_bs_scale1p12_smear0p03); 


          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_plus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_plus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_minus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_minus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 

          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_fullBPark(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor_fullBPark_plus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_fullBPark_plus_one_sigma(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor_fullBPark_minus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta)) * getTriggerScaleFactor_fullBPark_minus_one_sigma(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        }

        else if(the_sig_mu0_istriggering && !the_sig_mu_istriggering){
          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_plus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_minus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)); 

          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors_plus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors_minus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs)); 

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor_D1_plus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor_D1_minus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta));

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1_bs(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs), fabs(the_sig_mu0_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor_D1_bs_plus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs), fabs(the_sig_mu0_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor_D1_bs_minus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs), fabs(the_sig_mu0_eta));

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs_rdst)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs_rdst)); 

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs_rdst)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig_bs_rdst)); 

          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale/scale_factors.root", the_sig_mu0_pt, the_sig_mu0_dxysig_bs_linearscale); // abs removed as applied above 
          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03/scale_factors.root", the_sig_mu0_pt, the_sig_mu0_dxysig_bs_scale1p12_smear0p03); 

          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)); 
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_plus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)); 
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_minus_one_sigma.root", the_sig_mu0_pt, fabs(the_sig_mu0_dxysig)); 

          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_fullBPark(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta));
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor_fullBPark_plus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta));
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor_fullBPark_minus_one_sigma(the_sig_mu0_pt, fabs(the_sig_mu0_dxysig), fabs(the_sig_mu0_eta));
        }

        else if(!the_sig_mu0_istriggering && the_sig_mu_istriggering){
          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_plus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_minus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 

          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors_plus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6/scale_factors_minus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs)); 

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor_D1_plus_one_sigma(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor_D1_minus_one_sigma(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1_bs(the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor_D1_bs_plus_one_sigma(the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor_D1_bs_minus_one_sigma(the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs), fabs(the_sig_mu_eta));

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs_rdst)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs_rdst)); 

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs_rdst)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig_bs_rdst)); 

          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale/scale_factors.root", the_sig_mu_pt, the_sig_mu_dxysig_bs_linearscale); 
          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03/scale_factors.root", the_sig_mu_pt, the_sig_mu_dxysig_bs_scale1p12_smear0p03); 


          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_plus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2/scale_factors_minus_one_sigma.root", the_sig_mu_pt, fabs(the_sig_mu_dxysig)); 

          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma = getTriggerScaleFactor_fullBPark_plus_one_sigma(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma = getTriggerScaleFactor_fullBPark_minus_one_sigma(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        }

        else if(!the_sig_mu0_istriggering && !the_sig_mu_istriggering){
          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2 = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma = 1.;

          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2 = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma = 1.;

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma = 1.;

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2 = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_plus_one_sigma = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysigbs_max5e6_v2_smalltable_v2_minus_one_sigma = 1.;

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6_smalltable_v2 = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_inicat_max3e6_smalltable_v2 = 1.;

          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_softid_max3e6_smalltable_v2 = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid_max3e6_smalltable_v2 = 1.;

          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_linearscale_smalltable_v2 = 1.;
          the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_scale1p12_smear0p03_smalltable_v2 = 1.;

          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2 = 1.;
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_plus_one_sigma = 1.;
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2_smalltable_v2_minus_one_sigma = 1.;

          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = 1.;
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_plus_one_sigma = 1.;
          //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2_minus_one_sigma = 1.;
        }

        // pile-up weight
        the_sig_weight_pu_qcd_A = isMC ? getPUWeight("pileup_weight_dataA_mcAutumn18.root", *Pileup_nTrueInt) : 1.;
        the_sig_weight_pu_qcd_B = isMC ? getPUWeight("pileup_weight_dataB_mcAutumn18.root", *Pileup_nTrueInt) : 1.;
        the_sig_weight_pu_qcd_C = isMC ? getPUWeight("pileup_weight_dataC_mcAutumn18.root", *Pileup_nTrueInt) : 1.;
        the_sig_weight_pu_qcd_D = isMC ? getPUWeight("pileup_weight_dataD_mcAutumn18.root", *Pileup_nTrueInt) : 1.;
        the_sig_weight_pu_qcd_tot = isMC ? getPUWeight("pileup_weight_datatot_mcAutumn18.root", *Pileup_nTrueInt) : 1.;

        the_sig_weight_pu_sig_A = isMC ? getPUWeight("pileup_weight_dataA_sigAug21.root", *Pileup_nTrueInt) : 1.;
        the_sig_weight_pu_sig_B = isMC ? getPUWeight("pileup_weight_dataB_sigAug21.root", *Pileup_nTrueInt) : 1.;
        the_sig_weight_pu_sig_C = isMC ? getPUWeight("pileup_weight_dataC_sigAug21.root", *Pileup_nTrueInt) : 1.;
        the_sig_weight_pu_sig_D = isMC ? getPUWeight("pileup_weight_dataD_sigAug21.root", *Pileup_nTrueInt) : 1.;
        the_sig_weight_pu_sig_tot = isMC ? getPUWeight("pileup_weight_datatot_sigAug21.root", *Pileup_nTrueInt) : 1.;

        // lepton scale factor
        the_sig_weight_mu0_softid = isMC ? getLeptonScaleFactor("RunABCD_SF_MuonID_2018.root", "softid", the_sig_mu0_pt, abs(the_sig_mu0_eta)) : 1.;
        the_sig_weight_mu_looseid = isMC ? getLeptonScaleFactor("RunABCD_SF_MuonID_2018.root", "looseid", the_sig_mu_pt, abs(the_sig_mu_eta)) : 1.;

        // mc corrections
        the_sig_weight_mu0_dxy_bs = isMC ? getMCCorrection("mc_weight_probe_dxy_bs.root", the_sig_mu0_dxy_bs, 0.3) : 1.;
        the_sig_weight_mu_dxy_bs = isMC ? getMCCorrection("mc_weight_probe_dxy_bs.root", the_sig_mu_dxy_bs, 0.3) : 1.;
        the_sig_weight_mu0_dxysig_bs = isMC ? getMCCorrection("mc_weight_probe_dxy_sig_bs.root", the_sig_mu0_dxysig_bs, 60) : 1.;
        the_sig_weight_mu_dxysig_bs = isMC ? getMCCorrection("mc_weight_probe_dxy_sig_bs.root", the_sig_mu_dxysig_bs, 60) : 1.;


      } // end isMC
        
      signal_tree->Fill();
    } // end sound index
  }// end at least one candidate in the event


  //   ----- Histograms -----  //

  if(do_fillhistograms){

    // number of matched candidates in the event (only saved if nCand non zero)
    UInt_t nMatchedCand_sig = 0;

    sighist_ncand_perevent->Fill(nCand_sig);

    if(nCand_sig > 0){
      for(unsigned int iCand(0); iCand < nCand_sig; ++iCand){
        //cout << "cand " << iCand << " isMatched: " << BToMuMuPi_isMatched[iCand] << " b pt: " <<  BToMuMuPi_pt[iCand] << " hnl charge " << BToMuMuPi_hnl_charge[iCand] << endl;
        if(BToMuMuPi_isMatched[iCand] == 1) ++nMatchedCand_sig;
      }

      sighist_ncand_matched_perevent->Fill(nMatchedCand_sig);

      vector<pair<int,float>> pair_candIdx_desc_hnlpT_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_pt);
      vector<pair<int,float>> pair_candIdx_desc_bpt_sig = createPairWithDesc(nCand_sig, BToMuMuPi_pt);
      vector<pair<int,float>> pair_candIdx_desc_mu0pt_sig = createPairWithDesc(nCand_sig, BToMuMuPi_mu0_pt);
      vector<pair<int,float>> pair_candIdx_desc_pipt_sig = createPairWithDesc(nCand_sig, BToMuMuPi_fit_pi_pt);
      vector<pair<int,float>> pair_candIdx_desc_svprob_sig = createPairWithDesc(nCand_sig, BToMuMuPi_sv_prob);
      vector<pair<int,float>> pair_candIdx_desc_svchi2_sig = createPairWithDesc(nCand_sig, BToMuMuPi_sv_chi2);
      vector<pair<int,float>> pair_candIdx_desc_cos2D_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_cos2D);
      vector<pair<int,float>> pair_candIdx_desc_dr_sig = createPairWithDesc(nCand_sig, BToMuMuPi_dr_mu0_hnl);
      //vector<pair<int,float>> pair_candIdx_desc_hnliso4_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_iso04_close);

      sort(pair_candIdx_desc_hnlpT_sig.begin(), pair_candIdx_desc_hnlpT_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_bpt_sig.begin(), pair_candIdx_desc_bpt_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_mu0pt_sig.begin(), pair_candIdx_desc_mu0pt_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_pipt_sig.begin(), pair_candIdx_desc_pipt_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_svprob_sig.begin(), pair_candIdx_desc_svprob_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_svchi2_sig.begin(), pair_candIdx_desc_svchi2_sig.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_cos2D_sig.begin(), pair_candIdx_desc_cos2D_sig.end(), sortcansbydesc);
      //sort(pair_candIdx_desc_hnliso4_sig.begin(), pair_candIdx_desc_hnliso4_sig.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_dr_sig.begin(), pair_candIdx_desc_dr_sig.end(), sortcansbydesc_opp);

      // then privilege OS cand over SS ones
      vector<pair<int,float>> pair_candIdx_desc_hnlpT_sig_up = updatePairWithDesc(pair_candIdx_desc_hnlpT_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_bpt_sig_up = updatePairWithDesc(pair_candIdx_desc_bpt_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_mu0pt_sig_up = updatePairWithDesc(pair_candIdx_desc_mu0pt_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_pipt_sig_up = updatePairWithDesc(pair_candIdx_desc_pipt_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_svprob_sig_up = updatePairWithDesc(pair_candIdx_desc_svprob_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_svchi2_sig_up = updatePairWithDesc(pair_candIdx_desc_svchi2_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_cos2D_sig_up = updatePairWithDesc(pair_candIdx_desc_cos2D_sig, BToMuMuPi_hnl_charge);
      //vector<pair<int,float>> pair_candIdx_desc_hnliso4_sig_up = updatePairWithDesc(pair_candIdx_desc_hnliso4_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_dr_sig_up = updatePairWithDesc(pair_candIdx_desc_dr_sig, BToMuMuPi_hnl_charge);

      sort(pair_candIdx_desc_hnlpT_sig_up.begin(), pair_candIdx_desc_hnlpT_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_bpt_sig_up.begin(), pair_candIdx_desc_bpt_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_mu0pt_sig_up.begin(), pair_candIdx_desc_mu0pt_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_pipt_sig_up.begin(), pair_candIdx_desc_pipt_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_svprob_sig_up.begin(), pair_candIdx_desc_svprob_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_svchi2_sig_up.begin(), pair_candIdx_desc_svchi2_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_cos2D_sig_up.begin(), pair_candIdx_desc_cos2D_sig_up.end(), sortcansbydesc_opp);
      //sort(pair_candIdx_desc_hnliso4_sig_up.begin(), pair_candIdx_desc_hnliso4_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_dr_sig_up.begin(), pair_candIdx_desc_dr_sig_up.end(), sortcansbydesc_opp);

      //if(nCand_sig>1){
      //for(unsigned int iPair(0); iPair < pair_candIdx_desc_bpt_sig.size(); ++iPair){
      //  cout << "idx: " << pair_candIdx_desc_bpt_sig[iPair].first << " b pt: " << pair_candIdx_desc_bpt_sig[iPair].second << endl;
      //}
      //for(unsigned int iPair(0); iPair < pair_candIdx_desc_bpt_sig_up.size(); ++iPair){
      //  cout << "idx: " << pair_candIdx_desc_bpt_sig_up[iPair].first << " charge: " << pair_candIdx_desc_bpt_sig_up[iPair].second << endl;
      //}
      //}
      //cout << "selected cand idx: " << selectedCandIdx_sig << " is matched: " << BToMuMuPi_isMatched[selectedCandIdx_sig] << endl;

      // number of matched selected candidate per event (0 or 1)
      UInt_t selEff_hnlpt_sig = BToMuMuPi_isMatched[pair_candIdx_desc_hnlpT_sig_up[0].first];
      UInt_t selEff_bpt_sig = BToMuMuPi_isMatched[pair_candIdx_desc_bpt_sig_up[0].first];
      UInt_t selEff_mu0pt_sig = BToMuMuPi_isMatched[pair_candIdx_desc_mu0pt_sig_up[0].first];
      UInt_t selEff_pipt_sig = BToMuMuPi_isMatched[pair_candIdx_desc_pipt_sig_up[0].first];
      UInt_t selEff_svprob_sig = BToMuMuPi_isMatched[pair_candIdx_desc_svprob_sig_up[0].first];
      UInt_t selEff_svchi2_sig = BToMuMuPi_isMatched[pair_candIdx_desc_svchi2_sig_up[0].first];
      UInt_t selEff_cos2d_sig = BToMuMuPi_isMatched[pair_candIdx_desc_cos2D_sig_up[0].first];
      //UInt_t selEff_hnliso4_sig = BToMuMuPi_isMatched[pair_candIdx_desc_hnliso4_sig_up[0].first];
      UInt_t selEff_dr_sig = BToMuMuPi_isMatched[pair_candIdx_desc_dr_sig_up[0].first];

      if(nMatchedCand_sig != 0){
        sighist_selection_efficiency_hnlpt_allevents->Fill(selEff_hnlpt_sig);
        sighist_selection_efficiency_bpt_allevents->Fill(selEff_bpt_sig);
        sighist_selection_efficiency_mu0pt_allevents->Fill(selEff_mu0pt_sig);
        sighist_selection_efficiency_pipt_allevents->Fill(selEff_pipt_sig);
        sighist_selection_efficiency_svprob_allevents->Fill(selEff_svprob_sig);
        sighist_selection_efficiency_svchi2_allevents->Fill(selEff_svchi2_sig);
        sighist_selection_efficiency_cos2d_allevents->Fill(selEff_cos2d_sig);
        //sighist_selection_efficiency_hnliso4_allevents->Fill(selEff_hnliso4_sig);
        sighist_selection_efficiency_dr_allevents->Fill(selEff_dr_sig);

        if(nCand_sig > 1) sighist_selection_efficiency_hnlpt_eventswithmultcands->Fill(selEff_hnlpt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_bpt_eventswithmultcands->Fill(selEff_bpt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_mu0pt_eventswithmultcands->Fill(selEff_mu0pt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_pipt_eventswithmultcands->Fill(selEff_pipt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_svprob_eventswithmultcands->Fill(selEff_svprob_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_svchi2_eventswithmultcands->Fill(selEff_svchi2_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_cos2d_eventswithmultcands->Fill(selEff_cos2d_sig);
        //if(nCand_sig > 1) sighist_selection_efficiency_hnliso4_eventswithmultcands->Fill(selEff_hnliso4_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_dr_eventswithmultcands->Fill(selEff_dr_sig);
      }
    } // end at least one candidate
  } // end fill histograms

  return kTRUE;
}


void BToMuMuPiDumper::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}


void BToMuMuPiDumper::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  my_file->cd();
  if(do_fillhistograms) my_file->Write();

  signal_tree->Write("", TObject::kOverwrite);

  my_file->Close();

  cout << "- End B->MuMuPi Dumper -" << endl;
}

