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

  signal_tree->Branch("trgmu_pt", &the_sig_trgmu_pt);
  signal_tree->Branch("trgmu_eta", &the_sig_trgmu_eta);
  signal_tree->Branch("trgmu_phi", &the_sig_trgmu_phi);
  signal_tree->Branch("trgmu_charge", &the_sig_trgmu_charge);
  signal_tree->Branch("trgmu_dxy", &the_sig_trgmu_dxy);
  signal_tree->Branch("trgmu_dxysig", &the_sig_trgmu_dxysig);
  signal_tree->Branch("trgmu_dz", &the_sig_trgmu_dz);
  signal_tree->Branch("trgmu_dzsig", &the_sig_trgmu_dzsig);
  signal_tree->Branch("trgmu_triggermatching_dr", &the_sig_trgmu_triggermatching_dr);
  signal_tree->Branch("trgmu_triggermatching_dpt", &the_sig_trgmu_triggermatching_dpt);
  //signal_tree->Branch("trgmu_ip3d", &the_sig_trgmu_ip3d);
  //signal_tree->Branch("trgmu_ip3dsig", &the_sig_trgmu_ip3dsig);
  signal_tree->Branch("trgmu_pfiso03", &the_sig_trgmu_pfiso03);
  signal_tree->Branch("trgmu_pfiso03_rel", &the_sig_trgmu_pfiso03_rel);
  signal_tree->Branch("trgmu_iso03", &the_sig_trgmu_iso03);
  signal_tree->Branch("trgmu_iso03_close", &the_sig_trgmu_iso03_close);
  signal_tree->Branch("trgmu_iso03_rel_close", &the_sig_trgmu_iso03_rel_close);
  signal_tree->Branch("trgmu_pfiso04", &the_sig_trgmu_pfiso04);
  signal_tree->Branch("trgmu_pfiso04_rel", &the_sig_trgmu_pfiso04_rel);
  signal_tree->Branch("trgmu_iso04", &the_sig_trgmu_iso04);
  signal_tree->Branch("trgmu_iso04_close", &the_sig_trgmu_iso04_close);
  signal_tree->Branch("trgmu_iso04_rel_close", &the_sig_trgmu_iso04_rel_close);
  signal_tree->Branch("trgmu_looseid", &the_sig_trgmu_looseid);
  signal_tree->Branch("trgmu_mediumid", &the_sig_trgmu_mediumid);
  signal_tree->Branch("trgmu_tightid", &the_sig_trgmu_tightid);
  signal_tree->Branch("trgmu_softid", &the_sig_trgmu_softid);
  signal_tree->Branch("trgmu_pfisoid", &the_sig_trgmu_pfisoid);
  signal_tree->Branch("trgmu_trkisoid", &the_sig_trgmu_trkisoid);
  signal_tree->Branch("trgmu_triggerlooseid", &the_sig_trgmu_triggerlooseid);
  signal_tree->Branch("trgmu_istriggering", &the_sig_trgmu_istriggering);
  signal_tree->Branch("trgmu_isPF", &the_sig_trgmu_isPF);
  signal_tree->Branch("trgmu_isglobalmuon", &the_sig_trgmu_isglobalmuon);
  signal_tree->Branch("trgmu_istrackermuon", &the_sig_trgmu_istrackermuon);
  signal_tree->Branch("trgmu_isglobalortrackermuon", &the_sig_trgmu_isglobalortrackermuon);
  signal_tree->Branch("trgmu_isglobalnottrackermuon", &the_sig_trgmu_isglobalnottrackermuon);
  signal_tree->Branch("trgmu_istrackernotglobalmuon", &the_sig_trgmu_istrackernotglobalmuon);
  signal_tree->Branch("trgmu_intimemuon", &the_sig_trgmu_intimemuon);
  signal_tree->Branch("trgmu_segmentcompatibility", &the_sig_trgmu_segmentcompatibility);
  signal_tree->Branch("trgmu_calocompatibility", &the_sig_trgmu_calocompatibility);
  signal_tree->Branch("trgmu_validhitfraction", &the_sig_trgmu_validhitfraction);
  signal_tree->Branch("trgmu_kinkfinderchi2", &the_sig_trgmu_kinkfinderchi2);
  signal_tree->Branch("trgmu_globalnormalisedchi2", &the_sig_trgmu_globalnormalisedchi2);
  signal_tree->Branch("trgmu_localpositionchi2", &the_sig_trgmu_localpositionchi2);
  signal_tree->Branch("trgmu_trackerhighpurityflag", &the_sig_trgmu_trackerhighpurityflag);
  signal_tree->Branch("trgmu_numberofvalidmuonhits", &the_sig_trgmu_numberofvalidmuonhits);
  signal_tree->Branch("trgmu_numberofvalidpixelhits", &the_sig_trgmu_numberofvalidpixelhits);
  signal_tree->Branch("trgmu_numberoftrackerlayers", &the_sig_trgmu_numberoftrackerlayers);
  signal_tree->Branch("trgmu_numberofpixellayers", &the_sig_trgmu_numberofpixellayers);
  signal_tree->Branch("trgmu_numberofstations", &the_sig_trgmu_numberofstations);
  signal_tree->Branch("trgmu_fired_hlt_mu7_ip4", &the_sig_trgmu_fired_hlt_mu7_ip4);
  signal_tree->Branch("trgmu_fired_hlt_mu8_ip3", &the_sig_trgmu_fired_hlt_mu8_ip3);
  signal_tree->Branch("trgmu_fired_hlt_mu8_ip5", &the_sig_trgmu_fired_hlt_mu8_ip5);
  signal_tree->Branch("trgmu_fired_hlt_mu8_ip6", &the_sig_trgmu_fired_hlt_mu8_ip6);
  signal_tree->Branch("trgmu_fired_hlt_mu8p5_ip3p5", &the_sig_trgmu_fired_hlt_mu8p5_ip3p5);
  signal_tree->Branch("trgmu_fired_hlt_mu9_ip4", &the_sig_trgmu_fired_hlt_mu9_ip4);
  signal_tree->Branch("trgmu_fired_hlt_mu9_ip5", &the_sig_trgmu_fired_hlt_mu9_ip5);
  signal_tree->Branch("trgmu_fired_hlt_mu9_ip6", &the_sig_trgmu_fired_hlt_mu9_ip6);
  signal_tree->Branch("trgmu_fired_hlt_mu10p5_ip3p5", &the_sig_trgmu_fired_hlt_mu10p5_ip3p5);
  signal_tree->Branch("trgmu_fired_hlt_mu12_ip6", &the_sig_trgmu_fired_hlt_mu12_ip6);
  signal_tree->Branch("trgmu_prescale_hlt_mu7_ip4", &the_sig_trgmu_prescale_hlt_mu7_ip4);
  signal_tree->Branch("trgmu_prescale_hlt_mu8_ip3", &the_sig_trgmu_prescale_hlt_mu8_ip3);
  signal_tree->Branch("trgmu_prescale_hlt_mu8_ip5", &the_sig_trgmu_prescale_hlt_mu8_ip5);
  signal_tree->Branch("trgmu_prescale_hlt_mu8_ip6", &the_sig_trgmu_prescale_hlt_mu8_ip6);
  signal_tree->Branch("trgmu_prescale_hlt_mu8p5_ip3p5", &the_sig_trgmu_prescale_hlt_mu8p5_ip3p5);
  signal_tree->Branch("trgmu_prescale_hlt_mu9_ip4", &the_sig_trgmu_prescale_hlt_mu9_ip4);
  signal_tree->Branch("trgmu_prescale_hlt_mu9_ip5", &the_sig_trgmu_prescale_hlt_mu9_ip5);
  signal_tree->Branch("trgmu_prescale_hlt_mu9_ip6", &the_sig_trgmu_prescale_hlt_mu9_ip6);
  signal_tree->Branch("trgmu_prescale_hlt_mu10p5_ip3p5", &the_sig_trgmu_prescale_hlt_mu10p5_ip3p5);
  signal_tree->Branch("trgmu_prescale_hlt_mu12_ip6", &the_sig_trgmu_prescale_hlt_mu12_ip6);

  signal_tree->Branch("mu_pt", &the_sig_mu_pt);
  signal_tree->Branch("mu_eta", &the_sig_mu_eta);
  signal_tree->Branch("mu_phi", &the_sig_mu_phi);
  signal_tree->Branch("mu_charge", &the_sig_mu_charge);
  signal_tree->Branch("mu_dxy", &the_sig_mu_dxy);
  signal_tree->Branch("mu_dxysig", &the_sig_mu_dxysig);
  signal_tree->Branch("mu_dz", &the_sig_mu_dz);
  signal_tree->Branch("mu_dzsig", &the_sig_mu_dzsig);
  signal_tree->Branch("mu_triggermatching_dr", &the_sig_mu_triggermatching_dr);
  signal_tree->Branch("mu_triggermatching_dpt", &the_sig_mu_triggermatching_dpt);
  signal_tree->Branch("mu_ismatchedtoslimmedmuon", &the_sig_mu_ismatchedtoslimmedmuon);
  signal_tree->Branch("mu_indexmatchedslimmedmuon", &the_sig_mu_indexmatchedslimmedmuon);
  signal_tree->Branch("mu_dsatoslimmedmatching_deltar", &the_sig_mu_dsatoslimmedmatching_deltar);
  signal_tree->Branch("mu_dsatoslimmedmatching_deltaptrel", &the_sig_mu_dsatoslimmedmatching_deltaptrel);
  signal_tree->Branch("mu_dsatoslimmedmatching_deltadxyrel", &the_sig_mu_dsatoslimmedmatching_deltadxyrel);
  signal_tree->Branch("mu_dsatoslimmedmatching_deltadzrel", &the_sig_mu_dsatoslimmedmatching_deltadzrel);
  signal_tree->Branch("mu_passdsaid", &the_sig_mu_passdsaid);
  //signal_tree->Branch("mu_ip3d", &the_sig_mu_ip3d);
  //signal_tree->Branch("mu_ip3dsig", &the_sig_mu_ip3dsig);
  signal_tree->Branch("mu_pfiso03", &the_sig_mu_pfiso03);
  signal_tree->Branch("mu_pfiso03_rel", &the_sig_mu_pfiso03_rel);
  signal_tree->Branch("mu_iso03", &the_sig_mu_iso03);
  signal_tree->Branch("mu_iso03_close", &the_sig_mu_iso03_close);
  signal_tree->Branch("mu_iso03_rel_close", &the_sig_mu_iso03_rel_close);
  signal_tree->Branch("mu_pfiso04", &the_sig_mu_pfiso04);
  signal_tree->Branch("mu_pfiso04_rel", &the_sig_mu_pfiso04_rel);
  signal_tree->Branch("mu_iso04", &the_sig_mu_iso04);
  signal_tree->Branch("mu_iso04_close", &the_sig_mu_iso04_close);
  signal_tree->Branch("mu_iso04_rel_close", &the_sig_mu_iso04_rel_close);
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
  signal_tree->Branch("mu_isglobalortrackermuon", &the_sig_mu_isglobalortrackermuon);
  signal_tree->Branch("mu_isglobalnottrackermuon", &the_sig_mu_isglobalnottrackermuon);
  signal_tree->Branch("mu_istrackernotglobalmuon", &the_sig_mu_istrackernotglobalmuon);
  signal_tree->Branch("mu_intimemuon", &the_sig_mu_intimemuon);
  signal_tree->Branch("mu_segmentcompatibility", &the_sig_mu_segmentcompatibility);
  signal_tree->Branch("mu_calocompatibility", &the_sig_mu_calocompatibility);
  signal_tree->Branch("mu_validhitfraction", &the_sig_mu_validhitfraction);
  signal_tree->Branch("mu_kinkfinderchi2", &the_sig_mu_kinkfinderchi2);
  signal_tree->Branch("mu_globalnormalisedchi2", &the_sig_mu_globalnormalisedchi2);
  signal_tree->Branch("mu_localpositionchi2", &the_sig_mu_localpositionchi2);
  signal_tree->Branch("mu_trackerhighpurityflag", &the_sig_mu_trackerhighpurityflag);
  signal_tree->Branch("mu_numberofvalidmuonhits", &the_sig_mu_numberofvalidmuonhits);
  signal_tree->Branch("mu_numberofvalidpixelhits", &the_sig_mu_numberofvalidpixelhits);
  signal_tree->Branch("mu_numberoftrackerlayers", &the_sig_mu_numberoftrackerlayers);
  signal_tree->Branch("mu_numberofpixellayers", &the_sig_mu_numberofpixellayers);
  signal_tree->Branch("mu_numberofstations", &the_sig_mu_numberofstations);

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
  //signal_tree->Branch("pi_trgmu_dr", &the_sig_pi_trgmu_dr);
  //signal_tree->Branch("pi_ismatchedtomuon", &the_sig_pi_ismatchedtomuon);
  signal_tree->Branch("pi_chi2", &the_sig_pi_chi2);
  signal_tree->Branch("pi_normalisedchi2", &the_sig_pi_normalisedChi2);
  signal_tree->Branch("pi_validfraction", &the_sig_pi_validFraction);
  signal_tree->Branch("pi_ndof", &the_sig_pi_ndof);
  signal_tree->Branch("pi_numberofvalidhits", &the_sig_pi_numberOfValidHits);
  signal_tree->Branch("pi_numberoflosthits", &the_sig_pi_numberOfLostHits);
  signal_tree->Branch("pi_numberofvalidpixelhits", &the_sig_pi_numberOfValidPixelHits);
  signal_tree->Branch("pi_numberoftrackerlayers", &the_sig_pi_numberOfTrackerLayers);
  signal_tree->Branch("pi_numberofpixellayers", &the_sig_pi_numberOfPixelLayers);
  signal_tree->Branch("pi_qualityindex", &the_sig_pi_qualityIndex);
  signal_tree->Branch("pi_highpurityflag", &the_sig_pi_highPurityFlag);
  signal_tree->Branch("pi_packedcandhashighpurity", &the_sig_pi_packedcandhashighpurity);
  signal_tree->Branch("pi_matchedtomuon_loose", &the_sig_pi_matchedtomuon_loose);
  signal_tree->Branch("pi_matchedtomuon_medium", &the_sig_pi_matchedtomuon_medium);
  signal_tree->Branch("pi_matchedtomuon_tight", &the_sig_pi_matchedtomuon_tight);

  //signal_tree->Branch("dimu_mass", &the_sig_dimu_mass);
  //signal_tree->Branch("dimu_pt", &the_sig_dimu_pt);
  signal_tree->Branch("trgmu_mu_mass", &the_sig_trgmu_mu_mass);
  signal_tree->Branch("trgmu_mu_pt", &the_sig_trgmu_mu_pt);
  signal_tree->Branch("trgmu_pi_mass", &the_sig_trgmu_pi_mass);
  signal_tree->Branch("trgmu_pi_pt", &the_sig_trgmu_pi_pt);
  signal_tree->Branch("dimu_lxy", &the_sig_dimu_lxy);
  signal_tree->Branch("dimu_lxyz", &the_sig_dimu_lxyz);
  signal_tree->Branch("dimu_vxdiff", &the_sig_dimu_vxdiff);
  signal_tree->Branch("dimu_vydiff", &the_sig_dimu_vydiff);
  signal_tree->Branch("dimu_vzdiff", &the_sig_dimu_vzdiff);

  signal_tree->Branch("cos_theta_star_pion", &the_sig_cos_theta_star_pion);
  signal_tree->Branch("cos_theta_star_muon", &the_sig_cos_theta_star_muon);
  signal_tree->Branch("cos_theta_star_sum", &the_sig_cos_theta_star_sum);

  signal_tree->Branch("px_diff_hnl_daughters_lab", &the_sig_px_diff_hnl_daughters_lab);
  signal_tree->Branch("py_diff_hnl_daughters_lab", &the_sig_py_diff_hnl_daughters_lab);
  signal_tree->Branch("pz_diff_hnl_daughters_lab", &the_sig_pz_diff_hnl_daughters_lab);
  signal_tree->Branch("energy_diff_prefithnl_daughters_lab", &the_sig_energy_diff_prefithnl_daughters_lab);
  signal_tree->Branch("px_diff_prefithnl_daughters_lab", &the_sig_px_diff_prefithnl_daughters_lab);
  signal_tree->Branch("py_diff_prefithnl_daughters_lab", &the_sig_py_diff_prefithnl_daughters_lab);
  signal_tree->Branch("pz_diff_prefithnl_daughters_lab", &the_sig_pz_diff_prefithnl_daughters_lab);

  signal_tree->Branch("deltar_mu_pi", &the_sig_deltar_mu_pi);
  signal_tree->Branch("deltar_trgmu_hnl", &the_sig_deltar_trgmu_hnl);
  signal_tree->Branch("deltar_trgmu_mu", &the_sig_deltar_trgmu_mu);
  signal_tree->Branch("deltar_trgmu_pi", &the_sig_deltar_trgmu_pi);
  signal_tree->Branch("deltaeta_mu_pi", &the_sig_deltaeta_mu_pi);
  signal_tree->Branch("deltaeta_trgmu_mu", &the_sig_deltaeta_trgmu_mu);
  signal_tree->Branch("deltaeta_trgmu_pi", &the_sig_deltaeta_trgmu_pi);
  signal_tree->Branch("deltaeta_trgmu_hnl", &the_sig_deltaeta_trgmu_hnl);
  signal_tree->Branch("deltaphi_mu_pi", &the_sig_deltaphi_mu_pi);
  signal_tree->Branch("deltaphi_trgmu_mu", &the_sig_deltaphi_trgmu_mu);
  signal_tree->Branch("deltaphi_trgmu_pi", &the_sig_deltaphi_trgmu_pi);
  signal_tree->Branch("deltaphi_trgmu_hnl", &the_sig_deltaphi_trgmu_hnl);

  signal_tree->Branch("deltae_pi_fit_pi", &the_sig_deltapt_pi_fit_pi);
  signal_tree->Branch("deltae_mu_fit_mu", &the_sig_deltapt_mu_fit_mu);
  signal_tree->Branch("deltae_hnl_fit_hnl", &the_sig_deltapt_hnl_fit_hnl);
  signal_tree->Branch("deltapt_pi_fit_pi", &the_sig_deltapt_pi_fit_pi);
  signal_tree->Branch("deltapt_mu_fit_mu", &the_sig_deltapt_mu_fit_mu);
  signal_tree->Branch("deltapt_hnl_fit_hnl", &the_sig_deltapt_hnl_fit_hnl);
  signal_tree->Branch("deltapx_pi_fit_pi", &the_sig_deltapx_pi_fit_pi);
  signal_tree->Branch("deltapx_mu_fit_mu", &the_sig_deltapx_mu_fit_mu);
  signal_tree->Branch("deltapx_hnl_fit_hnl", &the_sig_deltapx_hnl_fit_hnl);
  signal_tree->Branch("deltapy_pi_fit_pi", &the_sig_deltapy_pi_fit_pi);
  signal_tree->Branch("deltapy_mu_fit_mu", &the_sig_deltapy_mu_fit_mu);
  signal_tree->Branch("deltapy_hnl_fit_hnl", &the_sig_deltapy_hnl_fit_hnl);
  signal_tree->Branch("deltapz_pi_fit_pi", &the_sig_deltapz_pi_fit_pi);
  signal_tree->Branch("deltapz_mu_fit_mu", &the_sig_deltapz_mu_fit_mu);
  signal_tree->Branch("deltapz_hnl_fit_hnl", &the_sig_deltapz_hnl_fit_hnl);
  signal_tree->Branch("deltaeta_pi_fit_pi", &the_sig_deltaeta_pi_fit_pi);
  signal_tree->Branch("deltaeta_mu_fit_mu", &the_sig_deltaeta_mu_fit_mu);
  signal_tree->Branch("deltaeta_hnl_fit_hnl", &the_sig_deltaeta_hnl_fit_hnl);
  signal_tree->Branch("deltaphi_pi_fit_pi", &the_sig_deltaphi_pi_fit_pi);
  signal_tree->Branch("deltaphi_mu_fit_mu", &the_sig_deltaphi_mu_fit_mu);
  signal_tree->Branch("deltaphi_hnl_fit_hnl", &the_sig_deltaphi_hnl_fit_hnl);

  signal_tree->Branch("sv_chi2", &the_sig_sv_chi2);
  signal_tree->Branch("sv_lxy", &the_sig_sv_lxy);
  signal_tree->Branch("sv_pv_lxy", &the_sig_sv_pv_lxy);
  signal_tree->Branch("sv_pv_lxyz", &the_sig_sv_pv_lxyz);
  signal_tree->Branch("sv_lxysig", &the_sig_sv_lxysig);
  signal_tree->Branch("sv_lxyz", &the_sig_sv_lxyz);
  signal_tree->Branch("sv_prob", &the_sig_sv_prob);
  signal_tree->Branch("sv_x", &the_sig_sv_x);
  signal_tree->Branch("sv_y", &the_sig_sv_y);
  signal_tree->Branch("sv_z", &the_sig_sv_z);

  signal_tree->Branch("pi_mu_vzdiff", &the_sig_pi_mu_vzdiff);

  signal_tree->Branch("ismatched", &the_sig_ismatched);
  signal_tree->Branch("trgmu_ismatched", &the_sig_trgmu_ismatched);
  signal_tree->Branch("mu_ismatched", &the_sig_mu_ismatched);
  signal_tree->Branch("pi_ismatched", &the_sig_pi_ismatched);
  signal_tree->Branch("mupi_mass_reco_gen_reldiff", &the_sig_mupi_mass_reco_gen_reldiff);
  signal_tree->Branch("lxy_reco_gen_reldiff", &the_sig_lxy_reco_gen_reldiff);

  signal_tree->Branch("weight_hlt_A1", &the_sig_weight_hlt_A1);
  signal_tree->Branch("weight_hlt_A1_6", &the_sig_weight_hlt_A1_6);
  signal_tree->Branch("weight_hlt_HLT_Mu9_IP6_A1_6", &the_sig_weight_hlt_HLT_Mu9_IP6_A1_6);
  signal_tree->Branch("weight_hlt_HLT_Mu9_IP6_A1_6_v2", &the_sig_weight_hlt_HLT_Mu9_IP6_A1_6_v2);
  signal_tree->Branch("weight_hlt_A1_6_B1", &the_sig_weight_hlt_A1_6_B1);
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
  signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2);
  signal_tree->Branch("weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max5e6_v2", &the_sig_weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max5e6_v2);
  signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2);
  signal_tree->Branch("weight_hlt_fullBPark_tag_fired_anyBParkHLT_pteta_max5e6_v2", &the_sig_weight_hlt_fullBPark_tag_fired_anyBParkHLT_pteta_max5e6_v2);
  signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable);
  signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2);
  signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable);
  signal_tree->Branch("weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max3e6_v3", &the_sig_weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max3e6_v3);
  signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2);
  signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable);
  signal_tree->Branch("weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2);
  signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2);
  signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable);
  signal_tree->Branch("weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2", &the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2);
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

  if(isMC){
    signal_tree->Branch("gen_trgmu_mu_lxy", &the_gen_trgmu_mu_lxy);
    signal_tree->Branch("gen_trgmu_mu_lxyz", &the_gen_trgmu_mu_lxyz);
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
    signal_tree->Branch("gen_trgmu_pt", &the_gen_trgmu_pt);
    signal_tree->Branch("gen_trgmu_eta", &the_gen_trgmu_eta);
    signal_tree->Branch("gen_trgmu_phi", &the_gen_trgmu_phi);
    signal_tree->Branch("gen_trgmu_vx", &the_gen_trgmu_vx);
    signal_tree->Branch("gen_trgmu_vy", &the_gen_trgmu_vy);
    signal_tree->Branch("gen_trgmu_vz", &the_gen_trgmu_vz);
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

  signal_tree->Branch("hlt_mu7_ip4", &the_hlt_mu7_ip4);
  signal_tree->Branch("hlt_mu8_ip6", &the_hlt_mu8_ip6);
  signal_tree->Branch("hlt_mu8_ip5", &the_hlt_mu8_ip5);
  signal_tree->Branch("hlt_mu8_ip3", &the_hlt_mu8_ip3);
  signal_tree->Branch("hlt_mu8p5_ip3p5", &the_hlt_mu8p5_ip3p5);
  signal_tree->Branch("hlt_mu9_ip6", &the_hlt_mu9_ip6);
  signal_tree->Branch("hlt_mu9_ip5", &the_hlt_mu9_ip5);
  signal_tree->Branch("hlt_mu9_ip4", &the_hlt_mu9_ip4);
  signal_tree->Branch("hlt_mu10p5_ip3p5", &the_hlt_mu10p5_ip3p5);
  signal_tree->Branch("hlt_mu12_ip6", &the_hlt_mu12_ip6);


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
    sighist_selection_efficiency_trgmupt_allevents = new TH1F("sighist_selection_efficiency_trgmupt_allevents", "Efficiency of selection of candidate with largest trigger muon pT (all events)", 2, 0, 2);
    sighist_selection_efficiency_trgmupt_eventswithmultcands = new TH1F("sighist_selection_efficiency_trgmupt_eventswithmultcands", "Efficiency of selection of candidate with largest trigger muon pT (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_pipt_allevents = new TH1F("sighist_selection_efficiency_pipt_allevents", "Efficiency of selection of candidate with largest pi pT (all events)", 2, 0, 2);
    sighist_selection_efficiency_pipt_eventswithmultcands = new TH1F("sighist_selection_efficiency_pipt_eventswithmultcands", "Efficiency of selection of candidate with largest pi pT (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_svprob_allevents = new TH1F("sighist_selection_efficiency_svprob_allevents", "Efficiency of selection of candidate with largest vertex prob. (all events)", 2, 0, 2);
    sighist_selection_efficiency_svprob_eventswithmultcands = new TH1F("sighist_selection_efficiency_svprob_eventswithmultcands", "Efficiency of selection of candidate with largest vertex prob. (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_svchi2_allevents = new TH1F("sighist_selection_efficiency_svchi2_allevents", "Efficiency of selection of candidate with largest vertex chi2 (all events)", 2, 0, 2);
    sighist_selection_efficiency_svchi2_eventswithmultcands = new TH1F("sighist_selection_efficiency_svchi2_eventswithmultcands", "Efficiency of selection of candidate with largest vertex chi2 (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_cos2d_allevents = new TH1F("sighist_selection_efficiency_cos2d_allevents", "Efficiency of selection of candidate with largest vertex cos(angle) (all events)", 2, 0, 2);
    sighist_selection_efficiency_cos2d_eventswithmultcands = new TH1F("sighist_selection_efficiency_cos2d_eventswithmultcands", "Efficiency of selection of candidate with largest vertex cos(angle) (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_hnliso4_allevents = new TH1F("sighist_selection_efficiency_hnliso4_allevents", "Efficiency of selection of candidate with smallest hnl isolation (all events)", 2, 0, 2);
    sighist_selection_efficiency_hnliso4_eventswithmultcands = new TH1F("sighist_selection_efficiency_hnliso4_eventswithmultcands", "Efficiency of selection of candidate with smallest hnl isolation (events with multiple candidates)", 2, 0, 2);
    sighist_selection_efficiency_dr_allevents = new TH1F("sighist_selection_efficiency_dr_allevents", "Efficiency of selection of candidate with smallest dR(trgmu, hnl) (all events)", 2, 0, 2);
    sighist_selection_efficiency_dr_eventswithmultcands = new TH1F("sighist_selection_efficiency_dr_eventswithmultcands", "Efficiency of selection of candidate with smallest dR(trgmu, hnl) (events with multiple candidates)", 2, 0, 2);

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

  the_hlt_mu7_ip4 = *HLT_Mu7_IP4;
  the_hlt_mu8_ip6 = *HLT_Mu8_IP6;
  the_hlt_mu8_ip5 = *HLT_Mu8_IP5;
  the_hlt_mu8_ip3 = *HLT_Mu8_IP3;
  the_hlt_mu8p5_ip3p5 = *HLT_Mu8p5_IP3p5;
  the_hlt_mu9_ip6 = *HLT_Mu9_IP6;
  the_hlt_mu9_ip5 = *HLT_Mu9_IP5;
  the_hlt_mu9_ip4 = *HLT_Mu9_IP4;
  the_hlt_mu10p5_ip3p5 = *HLT_Mu10p5_IP3p5;
  the_hlt_mu12_ip6 = *HLT_Mu12_IP6;

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

    // fill the signal_tree
    if(BToMuMuPi_trg_mu_pt[selectedCandIdx_sig] == Muon_pt[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]]){ // temporary condition, skip events with faulty indexing
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

      the_sig_trgmu_pt = BToMuMuPi_trg_mu_pt[selectedCandIdx_sig];
      //the_sig_trgmu_pt = Muon_pt[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_eta = BToMuMuPi_trg_mu_eta[selectedCandIdx_sig];
      //the_sig_trgmu_eta = Muon_eta[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_phi = BToMuMuPi_trg_mu_phi[selectedCandIdx_sig];
      //the_sig_trgmu_phi = Muon_phi[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_charge = Muon_charge[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_dxy = Muon_dxy[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_dxysig = Muon_dxyS[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_dz = Muon_dz[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_dzsig = Muon_dzS[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_triggermatching_dr = Muon_matched_dr[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_triggermatching_dpt = Muon_matched_dpt[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      //the_sig_trgmu_ip3d = Muon_ip3d[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      //the_sig_trgmu_ip3dsig = Muon_sip3d[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_pfiso03 = Muon_pfiso03_all[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_pfiso03_rel = Muon_pfiso03Rel_all[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_iso03 = BToMuMuPi_trg_mu_iso03[selectedCandIdx_sig];
      the_sig_trgmu_iso03_close = BToMuMuPi_trg_mu_iso03_close[selectedCandIdx_sig];
      the_sig_trgmu_iso03_rel_close = BToMuMuPi_trg_mu_iso03_rel_close[selectedCandIdx_sig];
      the_sig_trgmu_pfiso04 = Muon_pfiso04_all[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_pfiso04_rel = Muon_pfiso04Rel_all[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_iso04 = BToMuMuPi_trg_mu_iso04[selectedCandIdx_sig];
      the_sig_trgmu_iso04_close = BToMuMuPi_trg_mu_iso04_close[selectedCandIdx_sig];
      the_sig_trgmu_iso04_rel_close = BToMuMuPi_trg_mu_iso04_rel_close[selectedCandIdx_sig];
      the_sig_trgmu_looseid = Muon_looseId[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_mediumid = Muon_mediumId[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_tightid = Muon_tightId[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_softid = Muon_softId[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_pfisoid = Muon_pfIsoId[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_trkisoid = Muon_tkIsoId[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_triggerlooseid = Muon_triggerIdLoose[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_istriggering = Muon_isTriggeringBPark[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_isPF = Muon_isPF[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_isglobalmuon = Muon_isGlobalMuon[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_istrackermuon = Muon_isTrackerMuon[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_isglobalortrackermuon = Muon_isGlobalOrTrackerMuon[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_isglobalnottrackermuon = Muon_isGlobalNotTrackerMuon[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_istrackernotglobalmuon = Muon_isTrackerNotGlobalMuon[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_intimemuon = Muon_inTimeMuon[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_segmentcompatibility = Muon_segmentCompatibility[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_calocompatibility = Muon_caloCompatibility[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_validhitfraction = Muon_validHitFraction[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_kinkfinderchi2 = Muon_kinkFinderChi2[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_globalnormalisedchi2 = Muon_globalNormalisedChi2[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_localpositionchi2 = Muon_localPositionChi2[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_trackerhighpurityflag = Muon_trackerHighPurityFlag[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_numberofvalidmuonhits = Muon_numberOfValidMuonHits[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_numberofvalidpixelhits = Muon_numberOfValidPixelHits[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_numberoftrackerlayers = Muon_numberOfTrackerLayers[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_numberofpixellayers = Muon_numberOfPixelLayers[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_numberofstations = Muon_numberOfStations[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_fired_hlt_mu7_ip4 = Muon_fired_HLT_Mu7_IP4[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_fired_hlt_mu8_ip3 = Muon_fired_HLT_Mu8_IP3[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_fired_hlt_mu8_ip5 = Muon_fired_HLT_Mu8_IP5[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_fired_hlt_mu8_ip6 = Muon_fired_HLT_Mu8_IP6[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_fired_hlt_mu8p5_ip3p5 = Muon_fired_HLT_Mu8p5_IP3p5[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_fired_hlt_mu9_ip4 = Muon_fired_HLT_Mu9_IP4[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_fired_hlt_mu9_ip5 = Muon_fired_HLT_Mu9_IP5[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_fired_hlt_mu9_ip6 = Muon_fired_HLT_Mu9_IP6[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_fired_hlt_mu10p5_ip3p5 = Muon_fired_HLT_Mu10p5_IP3p5[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_fired_hlt_mu12_ip6 = Muon_fired_HLT_Mu12_IP6[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_prescale_hlt_mu7_ip4 = Muon_prescale_HLT_Mu7_IP4[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_prescale_hlt_mu8_ip3 = Muon_prescale_HLT_Mu8_IP3[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_prescale_hlt_mu8_ip5 = Muon_prescale_HLT_Mu8_IP5[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_prescale_hlt_mu8_ip6 = Muon_prescale_HLT_Mu8_IP6[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_prescale_hlt_mu8p5_ip3p5 = Muon_prescale_HLT_Mu8p5_IP3p5[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_prescale_hlt_mu9_ip4 = Muon_prescale_HLT_Mu9_IP4[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_prescale_hlt_mu9_ip5 = Muon_prescale_HLT_Mu9_IP5[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_prescale_hlt_mu9_ip6 = Muon_prescale_HLT_Mu9_IP6[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_prescale_hlt_mu10p5_ip3p5 = Muon_prescale_HLT_Mu10p5_IP3p5[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];
      the_sig_trgmu_prescale_hlt_mu12_ip6 = Muon_prescale_HLT_Mu12_IP6[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]];

      the_sig_mu_pt = BToMuMuPi_fit_mu_pt[selectedCandIdx_sig];
      the_sig_mu_eta = BToMuMuPi_fit_mu_eta[selectedCandIdx_sig];
      the_sig_mu_phi = BToMuMuPi_fit_mu_phi[selectedCandIdx_sig]; 
      the_sig_mu_charge = Muon_charge[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_dxy = Muon_dxy[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_dxysig = Muon_dxyS[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_dz = Muon_dz[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_dzsig = Muon_dzS[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_triggermatching_dr = Muon_matched_dr[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_triggermatching_dpt = Muon_matched_dpt[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_ip3d = Muon_ip3d[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      //the_sig_mu_ip3dsig = Muon_sip3d[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_pfiso03 = Muon_pfiso03_all[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_pfiso03_rel = Muon_pfiso03Rel_all[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_iso03 = BToMuMuPi_sel_mu_iso03[selectedCandIdx_sig];
      the_sig_mu_iso03_close = BToMuMuPi_sel_mu_iso03_close[selectedCandIdx_sig];
      the_sig_mu_iso03_rel_close = BToMuMuPi_sel_mu_iso03_rel_close[selectedCandIdx_sig];
      the_sig_mu_pfiso04 = Muon_pfiso04_all[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_pfiso04_rel = Muon_pfiso04Rel_all[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_iso04 = BToMuMuPi_sel_mu_iso04[selectedCandIdx_sig];
      the_sig_mu_iso04_close = BToMuMuPi_sel_mu_iso04_close[selectedCandIdx_sig];
      the_sig_mu_iso04_rel_close = BToMuMuPi_sel_mu_iso04_rel_close[selectedCandIdx_sig];
      //the_sig_mu_isloose = BToMuMuPi_sel_mu_isLoose[selectedCandIdx_sig];
      //the_sig_mu_ismedium = BToMuMuPi_sel_mu_isMedium[selectedCandIdx_sig];
      //the_sig_mu_istight = BToMuMuPi_sel_mu_isTight[selectedCandIdx_sig];
      //the_sig_mu_issoft = BToMuMuPi_sel_mu_isSoft[selectedCandIdx_sig];
      the_sig_mu_looseid = Muon_looseId[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_mediumid = Muon_mediumId[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_tightid = Muon_tightId[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_softid = Muon_softId[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_pfisoid = Muon_pfIsoId[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_trkisoid = Muon_tkIsoId[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_triggerlooseid = Muon_triggerIdLoose[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_istriggering = Muon_isTriggeringBPark[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_isslimmed = Muon_isSlimmedMuon[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_isdsa = Muon_isDSAMuon[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_isPF = Muon_isPF[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_isglobalmuon = Muon_isGlobalMuon[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_istrackermuon = Muon_isTrackerMuon[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_isglobalortrackermuon = Muon_isGlobalOrTrackerMuon[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_isglobalnottrackermuon = Muon_isGlobalNotTrackerMuon[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_istrackernotglobalmuon = Muon_isTrackerNotGlobalMuon[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_intimemuon = Muon_inTimeMuon[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_segmentcompatibility = Muon_segmentCompatibility[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_calocompatibility = Muon_caloCompatibility[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_validhitfraction = Muon_validHitFraction[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_kinkfinderchi2 = Muon_kinkFinderChi2[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_globalnormalisedchi2 = Muon_globalNormalisedChi2[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_localpositionchi2 = Muon_localPositionChi2[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_trackerhighpurityflag = Muon_trackerHighPurityFlag[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_numberofvalidmuonhits = Muon_numberOfValidMuonHits[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_numberofvalidpixelhits = Muon_numberOfValidPixelHits[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_numberoftrackerlayers = Muon_numberOfTrackerLayers[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_numberofpixellayers = Muon_numberOfPixelLayers[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
      the_sig_mu_numberofstations = Muon_numberOfStations[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]];
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
      //the_sig_pi_trgmu_dr = ProbeTracks_drTrg[BToMuMuPi_pi_idx[selectedCandIdx_sig]];
      //the_sig_pi_ismatchedtomuon = ProbeTracks_isMatchedToMuon[BToMuMuPi_pi_idx[selectedCandIdx_sig]];
      the_sig_pi_chi2 = BToMuMuPi_pi_chi2[selectedCandIdx_sig];
      the_sig_pi_normalisedChi2 = BToMuMuPi_pi_normalisedChi2[selectedCandIdx_sig];
      the_sig_pi_validFraction = BToMuMuPi_pi_validFraction[selectedCandIdx_sig];
      the_sig_pi_ndof = BToMuMuPi_pi_ndof[selectedCandIdx_sig];
      the_sig_pi_numberOfValidHits = BToMuMuPi_pi_numberOfValidHits[selectedCandIdx_sig];
      the_sig_pi_numberOfLostHits = BToMuMuPi_pi_numberOfLostHits[selectedCandIdx_sig];
      the_sig_pi_numberOfValidPixelHits = BToMuMuPi_pi_numberOfValidPixelHits[selectedCandIdx_sig];
      the_sig_pi_numberOfTrackerLayers = BToMuMuPi_pi_numberOfTrackerLayers[selectedCandIdx_sig];
      the_sig_pi_numberOfPixelLayers = BToMuMuPi_pi_numberOfPixelLayers[selectedCandIdx_sig];
      the_sig_pi_qualityIndex = BToMuMuPi_pi_qualityIndex[selectedCandIdx_sig];
      the_sig_pi_highPurityFlag = BToMuMuPi_pi_highPurityFlag[selectedCandIdx_sig];
      if(the_sig_pi_ispacked && the_sig_pi_highPurityFlag){
        the_sig_pi_packedcandhashighpurity = 1;
      }
      else{
        the_sig_pi_packedcandhashighpurity = 0;
      }

      // track to muon matching
      bool pi_matchedtomuon_loose = 0;
      bool pi_matchedtomuon_medium = 0;
      bool pi_matchedtomuon_tight = 0;
      /*
      for(unsigned int iMuon(0); iMuon<*nMuon; ++iMuon){
        // do not consider muons in the signal final state
        if(iMuon == BToMuMuPi_trg_mu_idx[selectedCandIdx_sig] || iMuon == BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]) continue;
        // compute the deltaR and deltaPtRel between the given track and the muon (unfitted values)
        float deltaR_track_muon = reco::deltaR(Muon_eta[iMuon], Muon_phi[iMuon], BToMuMuPi_pi_eta[selectedCandIdx_sig], BToMuMuPi_pi_phi[selectedCandIdx_sig]);
        float deltaPtRel = fabs(Muon_pt[iMuon] - BToMuMuPi_pi_pt[selectedCandIdx_sig]) / BToMuMuPi_pi_pt[selectedCandIdx_sig];

        // match if any muon fulfills requirements
        if(deltaR_track_muon < 0.3 && deltaPtRel < 1.){
          pi_matchedtomuon_loose = 1;
        }
        if(deltaR_track_muon < 0.3 && deltaPtRel < 0.3){
          pi_matchedtomuon_medium = 1;
        }
        if(deltaR_track_muon < 0.1 && deltaPtRel < 0.3){
          pi_matchedtomuon_tight = 1;
        }
      }
      */
      the_sig_pi_matchedtomuon_loose = pi_matchedtomuon_loose;
      the_sig_pi_matchedtomuon_medium = pi_matchedtomuon_medium;
      the_sig_pi_matchedtomuon_tight = pi_matchedtomuon_tight;


      the_sig_trgmu_mu_mass = BToMuMuPi_trgmu_mu_mass[selectedCandIdx_sig];
      the_sig_trgmu_mu_pt = BToMuMuPi_trgmu_mu_pt[selectedCandIdx_sig];
      the_sig_trgmu_pi_mass = BToMuMuPi_trgmu_pi_mass[selectedCandIdx_sig];
      the_sig_trgmu_pi_pt = BToMuMuPi_trgmu_pi_pt[selectedCandIdx_sig];

      the_sig_dimu_lxy = BToMuMuPi_dimu_Lxy[selectedCandIdx_sig];
      the_sig_dimu_lxyz = BToMuMuPi_dimu_Lxyz[selectedCandIdx_sig];
      the_sig_dimu_vxdiff = BToMuMuPi_dimu_vxdiff[selectedCandIdx_sig];
      the_sig_dimu_vydiff = BToMuMuPi_dimu_vydiff[selectedCandIdx_sig];
      the_sig_dimu_vzdiff = BToMuMuPi_dimu_vzdiff[selectedCandIdx_sig];

      the_sig_cos_theta_star_pion = BToMuMuPi_cos_theta_star_pion[selectedCandIdx_sig];
      the_sig_cos_theta_star_muon = BToMuMuPi_cos_theta_star_muon[selectedCandIdx_sig];
      the_sig_cos_theta_star_sum = fabs(BToMuMuPi_cos_theta_star_pion[selectedCandIdx_sig] + BToMuMuPi_cos_theta_star_muon[selectedCandIdx_sig]);

      the_sig_px_diff_hnl_daughters_lab = BToMuMuPi_px_diff_hnl_daughters_lab[selectedCandIdx_sig];
      the_sig_py_diff_hnl_daughters_lab = BToMuMuPi_py_diff_hnl_daughters_lab[selectedCandIdx_sig];
      the_sig_pz_diff_hnl_daughters_lab = BToMuMuPi_pz_diff_hnl_daughters_lab[selectedCandIdx_sig];
      the_sig_energy_diff_prefithnl_daughters_lab = BToMuMuPi_energy_diff_prefithnl_daughters_lab[selectedCandIdx_sig];
      the_sig_px_diff_prefithnl_daughters_lab = BToMuMuPi_px_diff_prefithnl_daughters_lab[selectedCandIdx_sig];
      the_sig_py_diff_prefithnl_daughters_lab = BToMuMuPi_py_diff_prefithnl_daughters_lab[selectedCandIdx_sig];
      the_sig_pz_diff_prefithnl_daughters_lab = BToMuMuPi_pz_diff_prefithnl_daughters_lab[selectedCandIdx_sig];

      the_sig_deltar_mu_pi = BToMuMuPi_dr_mu_pi[selectedCandIdx_sig];
      the_sig_deltar_trgmu_hnl = BToMuMuPi_dr_trgmu_hnl[selectedCandIdx_sig];
      the_sig_deltar_trgmu_mu = BToMuMuPi_dr_trgmu_mu[selectedCandIdx_sig]; 
      the_sig_deltar_trgmu_pi = BToMuMuPi_dr_trgmu_pi[selectedCandIdx_sig]; 
      the_sig_deltaeta_mu_pi = BToMuMuPi_deta_mu_pi[selectedCandIdx_sig];
      the_sig_deltaeta_trgmu_hnl = BToMuMuPi_deta_trgmu_hnl[selectedCandIdx_sig];
      the_sig_deltaeta_trgmu_mu = BToMuMuPi_deta_trgmu_mu[selectedCandIdx_sig];
      the_sig_deltaeta_trgmu_pi = BToMuMuPi_deta_trgmu_pi[selectedCandIdx_sig];
      the_sig_deltaphi_mu_pi = BToMuMuPi_dphi_mu_pi[selectedCandIdx_sig];
      the_sig_deltaphi_trgmu_hnl = BToMuMuPi_dphi_trgmu_hnl[selectedCandIdx_sig];
      the_sig_deltaphi_trgmu_mu = BToMuMuPi_dphi_trgmu_mu[selectedCandIdx_sig];
      the_sig_deltaphi_trgmu_pi = BToMuMuPi_dphi_trgmu_pi[selectedCandIdx_sig];

      the_sig_deltae_pi_fit_pi = BToMuMuPi_de_pi_fit_pi[selectedCandIdx_sig];
      the_sig_deltae_mu_fit_mu = BToMuMuPi_de_mu_fit_mu[selectedCandIdx_sig];
      the_sig_deltae_hnl_fit_hnl = BToMuMuPi_de_hnl_fit_hnl[selectedCandIdx_sig];
      the_sig_deltapt_pi_fit_pi = BToMuMuPi_dpt_pi_fit_pi[selectedCandIdx_sig];
      the_sig_deltapt_mu_fit_mu = BToMuMuPi_dpt_mu_fit_mu[selectedCandIdx_sig];
      the_sig_deltapt_hnl_fit_hnl = BToMuMuPi_dpt_hnl_fit_hnl[selectedCandIdx_sig];
      the_sig_deltapx_pi_fit_pi = BToMuMuPi_dpx_pi_fit_pi[selectedCandIdx_sig];
      the_sig_deltapx_mu_fit_mu = BToMuMuPi_dpx_mu_fit_mu[selectedCandIdx_sig];
      the_sig_deltapx_hnl_fit_hnl = BToMuMuPi_dpx_hnl_fit_hnl[selectedCandIdx_sig];
      the_sig_deltapy_pi_fit_pi = BToMuMuPi_dpy_pi_fit_pi[selectedCandIdx_sig];
      the_sig_deltapy_mu_fit_mu = BToMuMuPi_dpy_mu_fit_mu[selectedCandIdx_sig];
      the_sig_deltapy_hnl_fit_hnl = BToMuMuPi_dpy_hnl_fit_hnl[selectedCandIdx_sig];
      the_sig_deltapz_pi_fit_pi = BToMuMuPi_dpz_pi_fit_pi[selectedCandIdx_sig];
      the_sig_deltapz_mu_fit_mu = BToMuMuPi_dpz_mu_fit_mu[selectedCandIdx_sig];
      the_sig_deltapz_hnl_fit_hnl = BToMuMuPi_dpz_hnl_fit_hnl[selectedCandIdx_sig];
      the_sig_deltaeta_pi_fit_pi = BToMuMuPi_deta_pi_fit_pi[selectedCandIdx_sig];
      the_sig_deltaeta_mu_fit_mu = BToMuMuPi_deta_mu_fit_mu[selectedCandIdx_sig];
      the_sig_deltaeta_hnl_fit_hnl = BToMuMuPi_deta_hnl_fit_hnl[selectedCandIdx_sig];
      the_sig_deltaphi_pi_fit_pi = BToMuMuPi_dphi_pi_fit_pi[selectedCandIdx_sig];
      the_sig_deltaphi_mu_fit_mu = BToMuMuPi_dphi_mu_fit_mu[selectedCandIdx_sig];
      the_sig_deltaphi_hnl_fit_hnl = BToMuMuPi_dphi_hnl_fit_hnl[selectedCandIdx_sig];
      //float deltaphi_mu_pi = fabs(BToMuMuPi_fit_mu_phi[selectedCandIdx_sig] - BToMuMuPi_fit_pi_phi[selectedCandIdx_sig]);
      //the_sig_deltaphi_mu_pi = deltaphi_mu_pi > M_PI ? deltaphi_mu_pi : deltaphi_mu_pi - 2 * M_PI;

      the_sig_sv_chi2 = BToMuMuPi_sv_chi2[selectedCandIdx_sig];
      the_sig_sv_lxy = BToMuMuPi_sv_lxy[selectedCandIdx_sig];
      the_sig_sv_lxysig = BToMuMuPi_sv_lxy_sig[selectedCandIdx_sig];
      the_sig_sv_lxyz = BToMuMuPi_sv_lxyz[selectedCandIdx_sig];
      the_sig_sv_prob = BToMuMuPi_sv_prob[selectedCandIdx_sig];
      the_sig_sv_x = BToMuMuPi_sv_x[selectedCandIdx_sig];
      the_sig_sv_y = BToMuMuPi_sv_y[selectedCandIdx_sig];
      the_sig_sv_z = BToMuMuPi_sv_z[selectedCandIdx_sig];

      the_sig_ismatched = BToMuMuPi_isMatched[selectedCandIdx_sig];
      the_sig_trgmu_ismatched = BToMuMuPi_trg_mu_isMatched[selectedCandIdx_sig];
      the_sig_mu_ismatched = BToMuMuPi_sel_mu_isMatched[selectedCandIdx_sig];
      the_sig_pi_ismatched = BToMuMuPi_pi_isMatched[selectedCandIdx_sig];
      the_sig_mupi_mass_reco_gen_reldiff = BToMuMuPi_mupi_mass_reco_gen_reldiff[selectedCandIdx_sig];
      the_sig_lxy_reco_gen_reldiff = BToMuMuPi_lxy_reco_gen_reldiff[selectedCandIdx_sig];

      // additionnal displacement quantities
      float dist_sv_pv_xy = sqrt((BToMuMuPi_sv_x[selectedCandIdx_sig] - *PV_x) * (BToMuMuPi_sv_x[selectedCandIdx_sig] - *PV_x) + (BToMuMuPi_sv_y[selectedCandIdx_sig] - *PV_y) * (BToMuMuPi_sv_y[selectedCandIdx_sig] - *PV_y));
      float dist_sv_pv_xyz = sqrt((BToMuMuPi_sv_x[selectedCandIdx_sig] - *PV_x) * (BToMuMuPi_sv_x[selectedCandIdx_sig] - *PV_x) + (BToMuMuPi_sv_y[selectedCandIdx_sig] - *PV_y) * (BToMuMuPi_sv_y[selectedCandIdx_sig] - *PV_y) + (BToMuMuPi_sv_z[selectedCandIdx_sig] - *PV_z) * (BToMuMuPi_sv_z[selectedCandIdx_sig] - *PV_z));
      the_sig_sv_pv_lxy = dist_sv_pv_xy;
      the_sig_sv_pv_lxyz = dist_sv_pv_xyz;

      the_sig_pi_mu_vzdiff = BToMuMuPi_pi_mu_vzdiff[selectedCandIdx_sig];


      /*
      // getting the displacement at gen level
      if(isMC){
        UInt_t nGen = *nGenPart;

        float hnl_vx(0.), hnl_vy(0.), hnl_vz(0.);
        float mother_vx(0.), mother_vy(0.), mother_vz(0.);
        float trgmu_vx(0.), trgmu_vy(0.), trgmu_vz(0.);
        float mu_vx(0.), mu_vy(0.), mu_vz(0.);
        int mother_idx(-99), trgmu_idx(-99), mu_idx(-99);

        if(BToMuMuPi_isMatched[selectedCandIdx_sig]==1){
          trgmu_vx = GenPart_vx[Muon_genPartIdx[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]]];
          trgmu_vy = GenPart_vy[Muon_genPartIdx[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]]];
          trgmu_vz = GenPart_vz[Muon_genPartIdx[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]]];
          mu_vx = GenPart_vx[Muon_genPartIdx[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]]];
          mu_vy = GenPart_vy[Muon_genPartIdx[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]]];
          mu_vz = GenPart_vz[Muon_genPartIdx[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]]];

          the_gen_trgmu_mu_lxy = sqrt((trgmu_vx - mu_vx) * (trgmu_vx - mu_vx) + (trgmu_vy - mu_vy) * (trgmu_vy - mu_vy));
          the_gen_trgmu_mu_lxyz = sqrt((trgmu_vx - mu_vx) * (trgmu_vx - mu_vx) + (trgmu_vy - mu_vy) * (trgmu_vy - mu_vy) + (trgmu_vz - mu_vz) * (trgmu_vz - mu_vz));
        }
        else{
          the_gen_trgmu_mu_lxy = -99.;
          the_gen_trgmu_mu_lxyz = -99;
        }

        // gen information (no matching)

        // find idx of gen particles of interest
        int gen_hnl_idx(-99), gen_b_idx(-99), gen_trgmu_idx(-99), gen_mu_idx(-99), gen_pi_idx(-99);

        for(unsigned int iGen(0); iGen < nGen; ++iGen){
          if(abs(GenPart_pdgId[iGen])==9900015){
            gen_hnl_idx = iGen;
            gen_b_idx = GenPart_genPartIdxMother[iGen]; 
            break;
          }
        }
        for(unsigned int iGen(0); iGen < nGen; ++iGen){
          if(abs(GenPart_pdgId[iGen])==13 && GenPart_genPartIdxMother[iGen]==gen_b_idx)
            gen_trgmu_idx = iGen;
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
    */

      // trigger scale factor
      //the_sig_weight_hlt_A1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_DST_DoubleMu1_A1_v1/scaleFactor_results_cat_pt_eta_fit.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_A1_6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_DST_DoubleMu1_A1_6_v1/scaleFactor_results_cat_pt_eta_fit.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_HLT_Mu9_IP6_A1_6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_HLT_Mu9_IP6_A1_6/scaleFactor_results_cat_pt_eta_fit.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_HLT_Mu9_IP6_A1_6_v2 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullA_fired_HLT_Mu9_IP6_max5e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;

      //the_sig_weight_hlt_A1_6_B1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_DST_DoubleMu1_A1_6_B1_v1/scaleFactor_results_cat_pt_eta_fit.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;

      //the_sig_weight_hlt_D1_pteta_v1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_anyBParkHLT_pteta_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_D1_ptdxysig_v1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_anyBParkHLT_ptdxysig_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig)) : 1.;
      //the_sig_weight_hlt_D1_tag_fired_HLT_Mu12_IP6_pteta_v1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu12_IP6_pteta_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_D1_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_pteta_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;

      //the_sig_weight_hlt_fullBPark_pteta_v1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_anyBParkHLT_max5e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_pteta = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_pteta_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_ptdxysig = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu10p5_IP3p5_or_HLT_Mu8_IP3_or_HLT_Mu12_IP6_ptdxysig_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig)) : 1.;

      //the_sig_weight_hlt_fullBPark_ptetadxysig_max5e6 = isMC ? getTriggerScaleFactor_full(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_D1_ptetadxysig_max5e6 = isMC ? getTriggerScaleFactor_D1(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta)) : 1.;

      //the_sig_weight_hlt_fullA_V06_tag_fired_HLT_Mu9_IP6_pteta_max3e6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullA_V06_tag_fired_HLT_Mu9_IP6_pteta_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      //the_sig_weight_hlt_fullA_tag_fired_anyBParkHLT_pteta_max3e6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullA_tag_fired_anyBParkHLT_pteta_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;

      //the_sig_weight_hlt_D1_tag_fired_DST_DoubleMu1_pteta_max3e6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_DST_DoubleMu1_pteta_max3e6/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;

      the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      the_sig_weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max5e6_v2 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_anyBParkHLT_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      the_sig_weight_hlt_fullBPark_tag_fired_anyBParkHLT_pteta_max5e6_v2 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_anyBParkHLT_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;
      the_sig_weight_hlt_D1_tag_fired_anyBParkHLT_pteta_max3e6_v3 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_anyBParkHLT_pteta_max3e6_v3/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) : 1.;

      the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2 = isMC ? getTriggerScaleFactor_D1(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta)) : 1.;
      the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2 = isMC ? getTriggerScaleFactor_fullBPark(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta)) : 1.;

      the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = -99.;
      the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = -99.;
      the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = -99.;
      the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = -99.;
      if(isMC){
        if(the_sig_trgmu_istriggering && the_sig_mu_istriggering){
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = getTriggerScaleFactor_D1(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta)) * getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = getTriggerScaleFactor_fullBPark(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta)) * getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        }
        else if(the_sig_trgmu_istriggering && !the_sig_mu_istriggering){
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) * (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta))); 
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) * (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta))); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = getTriggerScaleFactor_D1(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta)) * (1-getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta)));
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = getTriggerScaleFactor_fullBPark(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta)) * (1-getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta)));
        }
        else if(!the_sig_trgmu_istriggering && the_sig_mu_istriggering){
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta))) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta))) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = (1-getTriggerScaleFactor_D1(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta))) * getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = (1-getTriggerScaleFactor_fullBPark(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta))) * getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        }
        else if(!the_sig_trgmu_istriggering && !the_sig_mu_istriggering){
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta))) * (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta))); 
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable = (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta))) * (1-getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta))); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = (1-getTriggerScaleFactor_D1(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta))) * (1-getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta)));
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable = (1-getTriggerScaleFactor_fullBPark(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta))) * (1-getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta)));
        }
      //}

      the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = -99.;
      the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = -99.;
      the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = -99.;
        if(the_sig_trgmu_istriggering && the_sig_mu_istriggering){
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)) * getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta)) * getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_fullBPark(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta)) * getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        }
        else if(the_sig_trgmu_istriggering && !the_sig_mu_istriggering){
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_trgmu_pt, fabs(the_sig_trgmu_eta)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta));
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_fullBPark(the_sig_trgmu_pt, fabs(the_sig_trgmu_dxysig), fabs(the_sig_trgmu_eta));
        }
        else if(!the_sig_trgmu_istriggering && the_sig_mu_istriggering){
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2/scale_factors.root", the_sig_mu_pt, fabs(the_sig_mu_eta)); 
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_D1(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = getTriggerScaleFactor_fullBPark(the_sig_mu_pt, fabs(the_sig_mu_dxysig), fabs(the_sig_mu_eta));
        }
        else if(!the_sig_trgmu_istriggering && !the_sig_mu_istriggering){
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_pteta_max5e6_v2_smalltable_v2 = 1.;
          //the_sig_weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = 1.;
          the_sig_weight_hlt_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2_smalltable_v2 = 1.;
        }
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
      vector<pair<int,float>> pair_candIdx_desc_trgmupt_sig = createPairWithDesc(nCand_sig, BToMuMuPi_trg_mu_pt);
      vector<pair<int,float>> pair_candIdx_desc_pipt_sig = createPairWithDesc(nCand_sig, BToMuMuPi_fit_pi_pt);
      vector<pair<int,float>> pair_candIdx_desc_svprob_sig = createPairWithDesc(nCand_sig, BToMuMuPi_sv_prob);
      vector<pair<int,float>> pair_candIdx_desc_svchi2_sig = createPairWithDesc(nCand_sig, BToMuMuPi_sv_chi2);
      vector<pair<int,float>> pair_candIdx_desc_cos2D_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_cos2D);
      vector<pair<int,float>> pair_candIdx_desc_dr_sig = createPairWithDesc(nCand_sig, BToMuMuPi_dr_trgmu_hnl);
      vector<pair<int,float>> pair_candIdx_desc_hnliso4_sig = createPairWithDesc(nCand_sig, BToMuMuPi_hnl_iso04_close);

      sort(pair_candIdx_desc_hnlpT_sig.begin(), pair_candIdx_desc_hnlpT_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_bpt_sig.begin(), pair_candIdx_desc_bpt_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_trgmupt_sig.begin(), pair_candIdx_desc_trgmupt_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_pipt_sig.begin(), pair_candIdx_desc_pipt_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_svprob_sig.begin(), pair_candIdx_desc_svprob_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_svchi2_sig.begin(), pair_candIdx_desc_svchi2_sig.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_cos2D_sig.begin(), pair_candIdx_desc_cos2D_sig.end(), sortcansbydesc);
      sort(pair_candIdx_desc_hnliso4_sig.begin(), pair_candIdx_desc_hnliso4_sig.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_dr_sig.begin(), pair_candIdx_desc_dr_sig.end(), sortcansbydesc_opp);

      // then privilege OS cand over SS ones
      vector<pair<int,float>> pair_candIdx_desc_hnlpT_sig_up = updatePairWithDesc(pair_candIdx_desc_hnlpT_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_bpt_sig_up = updatePairWithDesc(pair_candIdx_desc_bpt_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_trgmupt_sig_up = updatePairWithDesc(pair_candIdx_desc_trgmupt_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_pipt_sig_up = updatePairWithDesc(pair_candIdx_desc_pipt_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_svprob_sig_up = updatePairWithDesc(pair_candIdx_desc_svprob_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_svchi2_sig_up = updatePairWithDesc(pair_candIdx_desc_svchi2_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_cos2D_sig_up = updatePairWithDesc(pair_candIdx_desc_cos2D_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_hnliso4_sig_up = updatePairWithDesc(pair_candIdx_desc_hnliso4_sig, BToMuMuPi_hnl_charge);
      vector<pair<int,float>> pair_candIdx_desc_dr_sig_up = updatePairWithDesc(pair_candIdx_desc_dr_sig, BToMuMuPi_hnl_charge);

      sort(pair_candIdx_desc_hnlpT_sig_up.begin(), pair_candIdx_desc_hnlpT_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_bpt_sig_up.begin(), pair_candIdx_desc_bpt_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_trgmupt_sig_up.begin(), pair_candIdx_desc_trgmupt_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_pipt_sig_up.begin(), pair_candIdx_desc_pipt_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_svprob_sig_up.begin(), pair_candIdx_desc_svprob_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_svchi2_sig_up.begin(), pair_candIdx_desc_svchi2_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_cos2D_sig_up.begin(), pair_candIdx_desc_cos2D_sig_up.end(), sortcansbydesc_opp);
      sort(pair_candIdx_desc_hnliso4_sig_up.begin(), pair_candIdx_desc_hnliso4_sig_up.end(), sortcansbydesc_opp);
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
      UInt_t selEff_trgmupt_sig = BToMuMuPi_isMatched[pair_candIdx_desc_trgmupt_sig_up[0].first];
      UInt_t selEff_pipt_sig = BToMuMuPi_isMatched[pair_candIdx_desc_pipt_sig_up[0].first];
      UInt_t selEff_svprob_sig = BToMuMuPi_isMatched[pair_candIdx_desc_svprob_sig_up[0].first];
      UInt_t selEff_svchi2_sig = BToMuMuPi_isMatched[pair_candIdx_desc_svchi2_sig_up[0].first];
      UInt_t selEff_cos2d_sig = BToMuMuPi_isMatched[pair_candIdx_desc_cos2D_sig_up[0].first];
      UInt_t selEff_hnliso4_sig = BToMuMuPi_isMatched[pair_candIdx_desc_hnliso4_sig_up[0].first];
      UInt_t selEff_dr_sig = BToMuMuPi_isMatched[pair_candIdx_desc_dr_sig_up[0].first];

      if(nMatchedCand_sig != 0){
        sighist_selection_efficiency_hnlpt_allevents->Fill(selEff_hnlpt_sig);
        sighist_selection_efficiency_bpt_allevents->Fill(selEff_bpt_sig);
        sighist_selection_efficiency_trgmupt_allevents->Fill(selEff_trgmupt_sig);
        sighist_selection_efficiency_pipt_allevents->Fill(selEff_pipt_sig);
        sighist_selection_efficiency_svprob_allevents->Fill(selEff_svprob_sig);
        sighist_selection_efficiency_svchi2_allevents->Fill(selEff_svchi2_sig);
        sighist_selection_efficiency_cos2d_allevents->Fill(selEff_cos2d_sig);
        sighist_selection_efficiency_hnliso4_allevents->Fill(selEff_hnliso4_sig);
        sighist_selection_efficiency_dr_allevents->Fill(selEff_dr_sig);

        if(nCand_sig > 1) sighist_selection_efficiency_hnlpt_eventswithmultcands->Fill(selEff_hnlpt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_bpt_eventswithmultcands->Fill(selEff_bpt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_trgmupt_eventswithmultcands->Fill(selEff_trgmupt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_pipt_eventswithmultcands->Fill(selEff_pipt_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_svprob_eventswithmultcands->Fill(selEff_svprob_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_svchi2_eventswithmultcands->Fill(selEff_svchi2_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_cos2d_eventswithmultcands->Fill(selEff_cos2d_sig);
        if(nCand_sig > 1) sighist_selection_efficiency_hnliso4_eventswithmultcands->Fill(selEff_hnliso4_sig);
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

