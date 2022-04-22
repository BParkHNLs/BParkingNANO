//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 23 15:15:57 2021 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: bparknano_nj89.root
//////////////////////////////////////////////////////////

#ifndef BToMuMuPiDumper_h
#define BToMuMuPiDumper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class BToMuMuPiDumper : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> nBToMuMuPi = {fReader, "nBToMuMuPi"};
   //TTreeReaderArray<Float_t> BToMuMuPi_dilepton_mass = {fReader, "BToMuMuPi_dilepton_mass"};
   //TTreeReaderArray<Float_t> BToMuMuPi_dilepton_pt = {fReader, "BToMuMuPi_dilepton_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_trgmu_mu_mass = {fReader, "BToMuMuPi_trgmu_mu_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_trgmu_mu_pt = {fReader, "BToMuMuPi_trgmu_mu_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_trgmu_pi_mass = {fReader, "BToMuMuPi_trgmu_pi_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_trgmu_pi_pt = {fReader, "BToMuMuPi_trgmu_pi_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_dimu_Lxy = {fReader, "BToMuMuPi_dimu_Lxy"};
   TTreeReaderArray<Float_t> BToMuMuPi_dimu_Lxyz = {fReader, "BToMuMuPi_dimu_Lxyz"};
   TTreeReaderArray<Float_t> BToMuMuPi_dimu_vxdiff = {fReader, "BToMuMuPi_dimu_vxdiff"};
   TTreeReaderArray<Float_t> BToMuMuPi_dimu_vydiff = {fReader, "BToMuMuPi_dimu_vydiff"};
   TTreeReaderArray<Float_t> BToMuMuPi_dimu_vzdiff = {fReader, "BToMuMuPi_dimu_vzdiff"};
   TTreeReaderArray<Float_t> BToMuMuPi_dr_mu_pi = {fReader, "BToMuMuPi_dr_mu_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dr_trgmu_hnl = {fReader, "BToMuMuPi_dr_trgmu_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_dr_trgmu_mu = {fReader, "BToMuMuPi_dr_trgmu_mu"};
   TTreeReaderArray<Float_t> BToMuMuPi_dr_trgmu_pi = {fReader, "BToMuMuPi_dr_trgmu_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_de_pi_fit_pi = {fReader, "BToMuMuPi_de_pi_fit_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_de_mu_fit_mu = {fReader, "BToMuMuPi_de_mu_fit_mu"};
   TTreeReaderArray<Float_t> BToMuMuPi_de_hnl_fit_hnl = {fReader, "BToMuMuPi_de_hnl_fit_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_deta_mu_pi = {fReader, "BToMuMuPi_deta_mu_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_deta_trgmu_hnl = {fReader, "BToMuMuPi_deta_trgmu_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_deta_trgmu_mu = {fReader, "BToMuMuPi_deta_trgmu_mu"};
   TTreeReaderArray<Float_t> BToMuMuPi_deta_trgmu_pi = {fReader, "BToMuMuPi_deta_trgmu_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dphi_mu_pi = {fReader, "BToMuMuPi_dphi_mu_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dphi_trgmu_hnl = {fReader, "BToMuMuPi_dphi_trgmu_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_dphi_trgmu_mu = {fReader, "BToMuMuPi_dphi_trgmu_mu"};
   TTreeReaderArray<Float_t> BToMuMuPi_dphi_trgmu_pi = {fReader, "BToMuMuPi_dphi_trgmu_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpt_pi_fit_pi = {fReader, "BToMuMuPi_dpt_pi_fit_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpt_mu_fit_mu = {fReader, "BToMuMuPi_dpt_mu_fit_mu"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpt_hnl_fit_hnl = {fReader, "BToMuMuPi_dpt_hnl_fit_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpx_pi_fit_pi = {fReader, "BToMuMuPi_dpx_pi_fit_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpx_mu_fit_mu = {fReader, "BToMuMuPi_dpx_mu_fit_mu"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpx_hnl_fit_hnl = {fReader, "BToMuMuPi_dpx_hnl_fit_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpy_pi_fit_pi = {fReader, "BToMuMuPi_dpy_pi_fit_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpy_mu_fit_mu = {fReader, "BToMuMuPi_dpy_mu_fit_mu"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpy_hnl_fit_hnl = {fReader, "BToMuMuPi_dpy_hnl_fit_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpz_pi_fit_pi = {fReader, "BToMuMuPi_dpz_pi_fit_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpz_mu_fit_mu = {fReader, "BToMuMuPi_dpz_mu_fit_mu"};
   TTreeReaderArray<Float_t> BToMuMuPi_dpz_hnl_fit_hnl = {fReader, "BToMuMuPi_dpz_hnl_fit_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_deta_pi_fit_pi = {fReader, "BToMuMuPi_deta_pi_fit_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_deta_mu_fit_mu = {fReader, "BToMuMuPi_deta_mu_fit_mu"};
   TTreeReaderArray<Float_t> BToMuMuPi_deta_hnl_fit_hnl = {fReader, "BToMuMuPi_deta_hnl_fit_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_dphi_pi_fit_pi = {fReader, "BToMuMuPi_dphi_pi_fit_pi"};
   TTreeReaderArray<Float_t> BToMuMuPi_dphi_mu_fit_mu = {fReader, "BToMuMuPi_dphi_mu_fit_mu"};
   TTreeReaderArray<Float_t> BToMuMuPi_dphi_hnl_fit_hnl = {fReader, "BToMuMuPi_dphi_hnl_fit_hnl"};
   TTreeReaderArray<Float_t> BToMuMuPi_eta = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_mu_eta = {fReader, "BToMuMuPi_fit_mu_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_mu_mass = {fReader, "BToMuMuPi_fit_mu_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_mu_phi = {fReader, "BToMuMuPi_fit_mu_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_mu_pt = {fReader, "BToMuMuPi_fit_mu_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_pi_eta = {fReader, "BToMuMuPi_fit_pi_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_pi_mass = {fReader, "BToMuMuPi_fit_pi_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_pi_phi = {fReader, "BToMuMuPi_fit_pi_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_fit_pi_pt = {fReader, "BToMuMuPi_fit_pi_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_cos2D = {fReader, "BToMuMuPi_hnl_cos2D"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_ct = {fReader, "BToMuMuPi_hnl_ct"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_eta = {fReader, "BToMuMuPi_hnl_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_iso03 = {fReader, "BToMuMuPi_hnl_iso03"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_iso03_close = {fReader, "BToMuMuPi_hnl_iso03_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_iso03_rel_close = {fReader, "BToMuMuPi_hnl_iso03_rel_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_iso04 = {fReader, "BToMuMuPi_hnl_iso04"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_iso04_close = {fReader, "BToMuMuPi_hnl_iso04_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_iso04_rel_close = {fReader, "BToMuMuPi_hnl_iso04_rel_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_mass = {fReader, "BToMuMuPi_hnl_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_masserr = {fReader, "BToMuMuPi_hnl_masserr"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_phi = {fReader, "BToMuMuPi_hnl_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_hnl_pt = {fReader, "BToMuMuPi_hnl_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_mass = {fReader, "BToMuMuPi_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_phi = {fReader, "BToMuMuPi_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_cos_theta_star_pion = {fReader, "BToMuMuPi_cos_theta_star_pion"};
   TTreeReaderArray<Float_t> BToMuMuPi_cos_theta_star_muon = {fReader, "BToMuMuPi_cos_theta_star_muon"};
   TTreeReaderArray<Float_t> BToMuMuPi_px_diff_hnl_daughters_lab = {fReader, "BToMuMuPi_px_diff_hnl_daughters_lab"};
   TTreeReaderArray<Float_t> BToMuMuPi_py_diff_hnl_daughters_lab = {fReader, "BToMuMuPi_py_diff_hnl_daughters_lab"};
   TTreeReaderArray<Float_t> BToMuMuPi_pz_diff_hnl_daughters_lab = {fReader, "BToMuMuPi_pz_diff_hnl_daughters_lab"};
   TTreeReaderArray<Float_t> BToMuMuPi_energy_diff_prefithnl_daughters_lab = {fReader, "BToMuMuPi_energy_diff_prefithnl_daughters_lab"};
   TTreeReaderArray<Float_t> BToMuMuPi_px_diff_prefithnl_daughters_lab = {fReader, "BToMuMuPi_px_diff_prefithnl_daughters_lab"};
   TTreeReaderArray<Float_t> BToMuMuPi_py_diff_prefithnl_daughters_lab = {fReader, "BToMuMuPi_py_diff_prefithnl_daughters_lab"};
   TTreeReaderArray<Float_t> BToMuMuPi_pz_diff_prefithnl_daughters_lab = {fReader, "BToMuMuPi_pz_diff_prefithnl_daughters_lab"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_DCASig = {fReader, "BToMuMuPi_pi_DCASig"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_dxy = {fReader, "BToMuMuPi_pi_dxy"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_dxyS = {fReader, "BToMuMuPi_pi_dxyS"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_dz = {fReader, "BToMuMuPi_pi_dz"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_dzS = {fReader, "BToMuMuPi_pi_dzS"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_eta = {fReader, "BToMuMuPi_pi_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_mass = {fReader, "BToMuMuPi_pi_mass"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_phi = {fReader, "BToMuMuPi_pi_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_pt = {fReader, "BToMuMuPi_pi_pt"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_vx = {fReader, "BToMuMuPi_pi_vx"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_vy = {fReader, "BToMuMuPi_pi_vy"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_vz = {fReader, "BToMuMuPi_pi_vz"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_chi2 = {fReader, "BToMuMuPi_pi_chi2"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_normalisedChi2 = {fReader, "BToMuMuPi_pi_normalisedChi2"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_validFraction = {fReader, "BToMuMuPi_pi_validFraction"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_ndof = {fReader, "BToMuMuPi_pi_ndof"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_numberOfValidHits = {fReader, "BToMuMuPi_pi_numberOfValidHits"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_numberOfLostHits = {fReader, "BToMuMuPi_pi_numberOfLostHits"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_numberOfValidPixelHits = {fReader, "BToMuMuPi_pi_numberOfValidPixelHits"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_numberOfTrackerLayers = {fReader, "BToMuMuPi_pi_numberOfTrackerLayers"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_numberOfPixelLayers = {fReader, "BToMuMuPi_pi_numberOfPixelLayers"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_qualityIndex = {fReader, "BToMuMuPi_pi_qualityIndex"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_highPurityFlag = {fReader, "BToMuMuPi_pi_highPurityFlag"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_charge = {fReader, "BToMuMuPi_pi_charge"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_ispacked = {fReader, "BToMuMuPi_pi_ispacked"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_islost = {fReader, "BToMuMuPi_pi_islost"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_pdgId = {fReader, "BToMuMuPi_pi_pdgid"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_iso03 = {fReader, "BToMuMuPi_pi_iso03"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_iso03_close = {fReader, "BToMuMuPi_pi_iso03_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_iso03_rel_close = {fReader, "BToMuMuPi_pi_iso03_rel_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_iso04 = {fReader, "BToMuMuPi_pi_iso04"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_iso04_close = {fReader, "BToMuMuPi_pi_iso04_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_iso04_rel_close = {fReader, "BToMuMuPi_pi_iso04_rel_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_pi_mu_vzdiff = {fReader, "BToMuMuPi_pi_mu_vzdiff"};
   TTreeReaderArray<Float_t> BToMuMuPi_pt = {fReader, "BToMuMuPi_pt"};
   //TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_dxy = {fReader, "BToMuMuPi_sel_mu_dxy"};
   //TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_dz = {fReader, "BToMuMuPi_sel_mu_dz"};
   //TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_ip3d = {fReader, "BToMuMuPi_sel_mu_ip3d"};
   //TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_isLoose = {fReader, "BToMuMuPi_sel_mu_isLoose"};
   //TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_isMedium = {fReader, "BToMuMuPi_sel_mu_isMedium"};
   //TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_isSoft = {fReader, "BToMuMuPi_sel_mu_isSoft"};
   //TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_isTight = {fReader, "BToMuMuPi_sel_mu_isTight"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso03 = {fReader, "BToMuMuPi_sel_mu_iso03"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso03_close = {fReader, "BToMuMuPi_sel_mu_iso03_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso03_rel_close = {fReader, "BToMuMuPi_sel_mu_iso03_rel_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso04 = {fReader, "BToMuMuPi_sel_mu_iso04"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso04_close = {fReader, "BToMuMuPi_sel_mu_iso04_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso04_rel_close = {fReader, "BToMuMuPi_sel_mu_iso04_rel_close"};
   //TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_sip3d = {fReader, "BToMuMuPi_sel_mu_sip3d"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_chi2 = {fReader, "BToMuMuPi_sv_chi2"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_lxy = {fReader, "BToMuMuPi_sv_lxy"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_lxy_sig = {fReader, "BToMuMuPi_sv_lxy_sig"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_lxye = {fReader, "BToMuMuPi_sv_lxye"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_lxyz = {fReader, "BToMuMuPi_sv_lxyz"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_prob = {fReader, "BToMuMuPi_sv_prob"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_x = {fReader, "BToMuMuPi_sv_x"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_xe = {fReader, "BToMuMuPi_sv_xe"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_y = {fReader, "BToMuMuPi_sv_y"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_ye = {fReader, "BToMuMuPi_sv_ye"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_z = {fReader, "BToMuMuPi_sv_z"};
   TTreeReaderArray<Float_t> BToMuMuPi_sv_ze = {fReader, "BToMuMuPi_sv_ze"};
   //TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_dxy = {fReader, "BToMuMuPi_trg_mu_dxy"};
   //TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_dz = {fReader, "BToMuMuPi_trg_mu_dz"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_eta = {fReader, "BToMuMuPi_trg_mu_eta"};
   //TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_ip3d = {fReader, "BToMuMuPi_trg_mu_ip3d"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso03 = {fReader, "BToMuMuPi_trg_mu_iso03"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso03_close = {fReader, "BToMuMuPi_trg_mu_iso03_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso03_rel_close = {fReader, "BToMuMuPi_trg_mu_iso03_rel_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso04 = {fReader, "BToMuMuPi_trg_mu_iso04"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso04_close = {fReader, "BToMuMuPi_trg_mu_iso04_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso04_rel_close = {fReader, "BToMuMuPi_trg_mu_iso04_rel_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_phi = {fReader, "BToMuMuPi_trg_mu_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_pt = {fReader, "BToMuMuPi_trg_mu_pt"};
   //TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_sip3d = {fReader, "BToMuMuPi_trg_mu_sip3d"};
   TTreeReaderArray<Int_t> BToMuMuPi_charge = {fReader, "BToMuMuPi_charge"};
   TTreeReaderArray<Int_t> BToMuMuPi_hnl_charge = {fReader, "BToMuMuPi_hnl_charge"};
   TTreeReaderArray<Int_t> BToMuMuPi_isMatched = {fReader, "BToMuMuPi_isMatched"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_isMatched = {fReader, "BToMuMuPi_pi_isMatched"};
   TTreeReaderArray<Int_t> BToMuMuPi_sel_mu_isMatched = {fReader, "BToMuMuPi_sel_mu_isMatched"};
   TTreeReaderArray<Int_t> BToMuMuPi_trg_mu_isMatched = {fReader, "BToMuMuPi_trg_mu_isMatched"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_pi_genIdx = {fReader, "BToMuMuPi_matching_pi_genIdx"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_pi_motherPdgId = {fReader, "BToMuMuPi_matching_pi_motherPdgId"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_sel_mu_genIdx = {fReader, "BToMuMuPi_matching_sel_mu_genIdx"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_sel_mu_motherPdgId = {fReader, "BToMuMuPi_matching_sel_mu_motherPdgId"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_trg_mu_genIdx = {fReader, "BToMuMuPi_matching_trg_mu_genIdx"};
   TTreeReaderArray<Int_t> BToMuMuPi_matching_trg_mu_motherPdgId = {fReader, "BToMuMuPi_matching_trg_mu_motherPdgId"};
   TTreeReaderArray<Float_t> BToMuMuPi_mupi_mass_reco_gen_reldiff = {fReader, "BToMuMuPi_mupi_mass_reco_gen_reldiff"};
   TTreeReaderArray<Float_t> BToMuMuPi_lxy_reco_gen_reldiff = {fReader, "BToMuMuPi_lxy_reco_gen_reldiff"};
   TTreeReaderArray<Int_t> BToMuMuPi_pdgId = {fReader, "BToMuMuPi_pdgId"};
   TTreeReaderArray<Int_t> BToMuMuPi_pi_idx = {fReader, "BToMuMuPi_pi_idx"};
   TTreeReaderArray<Int_t> BToMuMuPi_sel_mu_idx = {fReader, "BToMuMuPi_sel_mu_idx"};
   TTreeReaderArray<Int_t> BToMuMuPi_trg_mu_idx = {fReader, "BToMuMuPi_trg_mu_idx"};


   TTreeReaderValue<UInt_t> nMuon = {fReader, "nMuon"};
   TTreeReaderArray<Int_t> Muon_isSlimmedMuon = {fReader, "Muon_isSlimmedMuon"};
   TTreeReaderArray<Int_t> Muon_isDSAMuon = {fReader, "Muon_isDSAMuon"};
   TTreeReaderArray<Int_t> Muon_isMatchedToSlimmedMuon = {fReader, "Muon_isMatchedToSlimmedMuon"};
   TTreeReaderArray<Int_t> Muon_indexMatchedSlimmedMuon = {fReader, "Muon_indexMatchedSlimmedMuon"};
   TTreeReaderArray<Int_t> Muon_passDSAMuonID = {fReader, "Muon_passDSAMuonID"};
   TTreeReaderArray<Float_t> Muon_dsaToSlimmedMatching_deltaR = {fReader, "Muon_dsaToSlimmedMatching_deltaR"};
   TTreeReaderArray<Float_t> Muon_dsaToSlimmedMatching_deltaPtRel = {fReader, "Muon_dsaToSlimmedMatching_deltaPtRel"};
   TTreeReaderArray<Float_t> Muon_dsaToSlimmedMatching_deltadxyRel = {fReader, "Muon_dsaToSlimmedMatching_deltadxyRel"};
   TTreeReaderArray<Float_t> Muon_dsaToSlimmedMatching_deltadzRel = {fReader, "Muon_dsaToSlimmedMatching_deltadzRel"};
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
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu10p5_IP3p5 = {fReader, "Muon_prescale_HLT_Mu10p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu12_IP6 = {fReader, "Muon_prescale_HLT_Mu12_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu7_IP4 = {fReader, "Muon_prescale_HLT_Mu7_IP4"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP3 = {fReader, "Muon_prescale_HLT_Mu8_IP3"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP5 = {fReader, "Muon_prescale_HLT_Mu8_IP5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP6 = {fReader, "Muon_prescale_HLT_Mu8_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8p5_IP3p5 = {fReader, "Muon_prescale_HLT_Mu8p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP4 = {fReader, "Muon_prescale_HLT_Mu9_IP4"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP5 = {fReader, "Muon_prescale_HLT_Mu9_IP5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP6 = {fReader, "Muon_prescale_HLT_Mu9_IP6"};
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
   //TTreeReaderArray<Float_t> ProbeTracks_dzTrg = {fReader, "ProbeTracks_dzTrg"};
   //TTreeReaderArray<Float_t> ProbeTracks_drTrg = {fReader, "ProbeTracks_drTrg"};
   //TTreeReaderArray<Float_t> ProbeTracks_eta = {fReader, "ProbeTracks_eta"};
   //TTreeReaderArray<Float_t> ProbeTracks_mass = {fReader, "ProbeTracks_mass"};
   //TTreeReaderArray<Float_t> ProbeTracks_phi = {fReader, "ProbeTracks_phi"};
   //TTreeReaderArray<Float_t> ProbeTracks_pt = {fReader, "ProbeTracks_pt"};
   //TTreeReaderArray<Float_t> ProbeTracks_vx = {fReader, "ProbeTracks_vx"};
   //TTreeReaderArray<Float_t> ProbeTracks_vy = {fReader, "ProbeTracks_vy"};
   //TTreeReaderArray<Float_t> ProbeTracks_vz = {fReader, "ProbeTracks_vz"};
   //TTreeReaderArray<Float_t> ProbeTracks_chi2 = {fReader, "ProbeTracks_chi2"};
   //TTreeReaderArray<Float_t> ProbeTracks_normalisedChi2 = {fReader, "ProbeTracks_normalisedChi2"};
   //TTreeReaderArray<Float_t> ProbeTracks_validFraction = {fReader, "ProbeTracks_validFraction"};
   //TTreeReaderArray<Int_t> ProbeTracks_ndof = {fReader, "ProbeTracks_ndof"};
   //TTreeReaderArray<Int_t> ProbeTracks_numberOfValidHits = {fReader, "ProbeTracks_numberOfValidHits"};
   //TTreeReaderArray<Int_t> ProbeTracks_numberOfLostHits = {fReader, "ProbeTracks_numberOfLostHits"};
   //TTreeReaderArray<Int_t> ProbeTracks_numberOfValidPixelHits = {fReader, "ProbeTracks_numberOfValidPixelHits"};
   //TTreeReaderArray<Int_t> ProbeTracks_numberOfTrackerLayers = {fReader, "ProbeTracks_numberOfTrackerLayers"};
   //TTreeReaderArray<Int_t> ProbeTracks_numberOfPixelLayers = {fReader, "ProbeTracks_numberOfPixelLayers"};
   //TTreeReaderArray<Int_t> ProbeTracks_qualityIndex = {fReader, "ProbeTracks_qualityIndex"};
   //TTreeReaderArray<Int_t> ProbeTracks_highPurityFlag = {fReader, "ProbeTracks_highPurityFlag"};
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
   // branches intentionally pointed to wrong variable, to avoid crash when running on data
   TTreeReaderValue<UInt_t> nGenPart = {fReader, "nBToMuMuPi"};
   TTreeReaderArray<Float_t> GenPart_eta = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_mass = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_phi = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_pt = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_vx = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_vy = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Float_t> GenPart_vz = {fReader, "BToMuMuPi_eta"};
   TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {fReader, "BToMuMuPi_charge"};
   TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "BToMuMuPi_charge"};
   TTreeReaderArray<Int_t> GenPart_status = {fReader, "BToMuMuPi_charge"};
   TTreeReaderArray<Int_t> GenPart_statusFlags = {fReader, "BToMuMuPi_charge"};
   TTreeReaderArray<Int_t> Muon_genPartIdx = {fReader, "BToMuMuPi_charge"};
   // the correct value will be filled for MC only
   TTreeReaderValue<Int_t> Pileup_nPU = {fReader, "PV_npvs"};
   TTreeReaderValue<Float_t> Pileup_nTrueInt = {fReader, "PV_ndof"};


   BToMuMuPiDumper(TTree * /*tree*/ =0) { }
   virtual ~BToMuMuPiDumper() { }
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
   TTree* signal_tree;
   
   // filling variables
   ULong64_t the_event = -99;
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

   Float_t the_sig_b_pt = -99.;
   Float_t the_sig_b_eta = -99.;
   Float_t the_sig_b_phi = -99.;
   Float_t the_sig_b_mass = -99.;
   Int_t the_sig_b_charge = -99;
   Int_t the_sig_b_pdgid = -99;

   Float_t the_sig_hnl_pt = -99.;
   Float_t the_sig_hnl_eta = -99.;
   Float_t the_sig_hnl_phi = -99.;
   Float_t the_sig_hnl_mass = -99.;
   Int_t the_sig_hnl_charge = -99;
   Float_t the_sig_hnl_ct = -99;
   Float_t the_sig_hnl_cos2d = -99.;
   //Float_t the_sig_hnl_iso03 = -99.;
   //Float_t the_sig_hnl_iso03_close = -99.;
   //Float_t the_sig_hnl_iso03_rel_close = -99.;
   //Float_t the_sig_hnl_iso04 = -99.;
   //Float_t the_sig_hnl_iso04_close = -99.;
   //Float_t the_sig_hnl_iso04_rel_close = -99.;

   Float_t the_sig_trgmu_pt = -99.;
   Float_t the_sig_trgmu_eta = -99.;
   Float_t the_sig_trgmu_phi = -99.;
   Int_t the_sig_trgmu_charge = -99;
   Float_t the_sig_trgmu_dxy = -99.;
   Float_t the_sig_trgmu_dxysig = -99.;
   Float_t the_sig_trgmu_dz = -99.;
   Float_t the_sig_trgmu_dzsig = -99.;
   //Float_t the_sig_trgmu_ip3d = -99.;
   //Float_t the_sig_trgmu_ip3dsig = -99.;
   Float_t the_sig_trgmu_pfiso03 = -99.;
   Float_t the_sig_trgmu_pfiso03_rel = -99.;
   Float_t the_sig_trgmu_iso03 = -99.;
   Float_t the_sig_trgmu_iso03_close = -99.;
   Float_t the_sig_trgmu_iso03_rel_close = -99.;
   Float_t the_sig_trgmu_pfiso04 = -99.;
   Float_t the_sig_trgmu_pfiso04_rel = -99.;
   Float_t the_sig_trgmu_iso04 = -99.;
   Float_t the_sig_trgmu_iso04_close = -99.;
   Float_t the_sig_trgmu_iso04_rel_close = -99.;
   Int_t the_sig_trgmu_looseid = -99;
   Int_t the_sig_trgmu_mediumid = -99;
   Int_t the_sig_trgmu_tightid = -99;
   Int_t the_sig_trgmu_softid = -99;
   Int_t the_sig_trgmu_pfisoid = -99;
   Int_t the_sig_trgmu_trkisoid = -99;
   Int_t the_sig_trgmu_triggerlooseid = -99;
   Int_t the_sig_trgmu_istriggering = -99;
   Int_t the_sig_trgmu_isPF = -99;
   Int_t the_sig_trgmu_isglobalmuon = -99;
   Int_t the_sig_trgmu_istrackermuon = -99;
   Int_t the_sig_trgmu_isglobalortrackermuon = -99;
   Int_t the_sig_trgmu_isglobalnottrackermuon = -99;
   Int_t the_sig_trgmu_istrackernotglobalmuon = -99;
   Int_t the_sig_trgmu_intimemuon = -99;
   Float_t the_sig_trgmu_segmentcompatibility = -99.;
   Float_t the_sig_trgmu_calocompatibility = -99.;
   Float_t the_sig_trgmu_validhitfraction = -99.;
   Float_t the_sig_trgmu_kinkfinderchi2 = -99.;
   Float_t the_sig_trgmu_globalnormalisedchi2 = -99.;
   Float_t the_sig_trgmu_localpositionchi2 = -99.;
   Int_t the_sig_trgmu_trackerhighpurityflag = -99;
   Int_t the_sig_trgmu_numberofvalidmuonhits = -99;
   Int_t the_sig_trgmu_numberofvalidpixelhits = -99;
   Int_t the_sig_trgmu_numberoftrackerlayers = -99;
   Int_t the_sig_trgmu_numberofpixellayers = -99;
   Int_t the_sig_trgmu_numberofstations = -99;
   Int_t the_sig_trgmu_fired_hlt_mu7_ip4 = -99;
   Int_t the_sig_trgmu_fired_hlt_mu8_ip6 = -99;
   Int_t the_sig_trgmu_fired_hlt_mu8_ip5 = -99;
   Int_t the_sig_trgmu_fired_hlt_mu8_ip3 = -99;
   Int_t the_sig_trgmu_fired_hlt_mu8p5_ip3p5 = -99;
   Int_t the_sig_trgmu_fired_hlt_mu9_ip6 = -99;
   Int_t the_sig_trgmu_fired_hlt_mu9_ip5 = -99;
   Int_t the_sig_trgmu_fired_hlt_mu9_ip4 = -99;
   Int_t the_sig_trgmu_fired_hlt_mu10p5_ip3p5 = -99;
   Int_t the_sig_trgmu_fired_hlt_mu12_ip6 = -99;
   Int_t the_sig_trgmu_prescale_hlt_mu7_ip4 = -99;
   Int_t the_sig_trgmu_prescale_hlt_mu8_ip6 = -99;
   Int_t the_sig_trgmu_prescale_hlt_mu8_ip5 = -99;
   Int_t the_sig_trgmu_prescale_hlt_mu8_ip3 = -99;
   Int_t the_sig_trgmu_prescale_hlt_mu8p5_ip3p5 = -99;
   Int_t the_sig_trgmu_prescale_hlt_mu9_ip6 = -99;
   Int_t the_sig_trgmu_prescale_hlt_mu9_ip5 = -99;
   Int_t the_sig_trgmu_prescale_hlt_mu9_ip4 = -99;
   Int_t the_sig_trgmu_prescale_hlt_mu10p5_ip3p5 = -99;
   Int_t the_sig_trgmu_prescale_hlt_mu12_ip6 = -99;

   Float_t the_sig_mu_pt = -99.;
   Float_t the_sig_mu_eta = -99.;
   Float_t the_sig_mu_phi = -99.;
   Int_t the_sig_mu_charge = -99;
   Float_t the_sig_mu_dxy = -99.;
   Float_t the_sig_mu_dxysig = -99.;
   Float_t the_sig_mu_dz = -99.;
   Float_t the_sig_mu_dzsig = -99.;
   Int_t the_sig_mu_ismatchedtoslimmedmuon = -99;
   Int_t the_sig_mu_indexmatchedslimmedmuon = -99;
   Float_t the_sig_mu_dsatoslimmedmatching_deltar = -99;
   Float_t the_sig_mu_dsatoslimmedmatching_deltaptrel = -99;
   Float_t the_sig_mu_dsatoslimmedmatching_deltadxyrel = -99;
   Float_t the_sig_mu_dsatoslimmedmatching_deltadzrel = -99;
   Int_t the_sig_mu_passdsaid = -99;
   //Float_t the_sig_mu_ip3d = -99.;
   //Float_t the_sig_mu_ip3dsig = -99.;
   Float_t the_sig_mu_pfiso03 = -99.;
   Float_t the_sig_mu_pfiso03_rel = -99.;
   Float_t the_sig_mu_iso03 = -99.;
   Float_t the_sig_mu_iso03_close = -99.;
   Float_t the_sig_mu_iso03_rel_close = -99.;
   Float_t the_sig_mu_pfiso04 = -99.;
   Float_t the_sig_mu_pfiso04_rel = -99.;
   Float_t the_sig_mu_iso04 = -99.;
   Float_t the_sig_mu_iso04_close = -99.;
   Float_t the_sig_mu_iso04_rel_close = -99.;
   Int_t the_sig_mu_looseid = -99;
   Int_t the_sig_mu_mediumid = -99;
   Int_t the_sig_mu_tightid = -99;
   Int_t the_sig_mu_softid = -99;
   Int_t the_sig_mu_pfisoid = -99;
   Int_t the_sig_mu_trkisoid = -99;
   Int_t the_sig_mu_triggerlooseid = -99;
   Int_t the_sig_mu_whnlid = -99;
   Int_t the_sig_mu_customisedid = -99;
   Int_t the_sig_mu_istriggering = -99;
   Int_t the_sig_mu_isslimmed = -99;
   Int_t the_sig_mu_isdsa = -99;
   Int_t the_sig_mu_isPF = -99;
   Int_t the_sig_mu_isglobalmuon = -99;
   Int_t the_sig_mu_istrackermuon = -99;
   Int_t the_sig_mu_isglobalortrackermuon = -99;
   Int_t the_sig_mu_isglobalnottrackermuon = -99;
   Int_t the_sig_mu_istrackernotglobalmuon = -99;
   Int_t the_sig_mu_intimemuon = -99;
   Float_t the_sig_mu_segmentcompatibility = -99.;
   Float_t the_sig_mu_calocompatibility = -99.;
   Float_t the_sig_mu_validhitfraction = -99.;
   Float_t the_sig_mu_kinkfinderchi2 = -99.;
   Float_t the_sig_mu_globalnormalisedchi2 = -99.;
   Float_t the_sig_mu_localpositionchi2 = -99.;
   Int_t the_sig_mu_trackerhighpurityflag = -99;
   Int_t the_sig_mu_numberofvalidmuonhits = -99;
   Int_t the_sig_mu_numberofvalidpixelhits = -99;
   Int_t the_sig_mu_numberoftrackerlayers = -99;
   Int_t the_sig_mu_numberofpixellayers = -99;
   Int_t the_sig_mu_numberofstations = -99;

   Float_t the_sig_pi_pt = -99.;
   Float_t the_sig_pi_eta = -99.;
   Float_t the_sig_pi_phi = -99.;
   Int_t the_sig_pi_charge = -99;
   Float_t the_sig_pi_dcasig = -99.;
   Float_t the_sig_pi_dxy = -99.;
   Float_t the_sig_pi_dz = -99.;
   Float_t the_sig_pi_dxysig = -99.;
   Float_t the_sig_pi_dzsig = -99.;
   Int_t the_sig_pi_ispacked = -99;
   Int_t the_sig_pi_islost = -99;
   //Float_t the_sig_pi_trgmu_dr = -99.;
   //Int_t the_sig_pi_ismatchedtomuon = -99.;
   Float_t the_sig_pi_chi2 = -99.;
   Float_t the_sig_pi_normalisedChi2 = -99.;
   Float_t the_sig_pi_validFraction = -99.;
   Int_t the_sig_pi_ndof = -99;
   Int_t the_sig_pi_numberOfValidHits = -99;
   Int_t the_sig_pi_numberOfLostHits = -99;
   Int_t the_sig_pi_numberOfValidPixelHits = -99;
   Int_t the_sig_pi_numberOfTrackerLayers = -99;
   Int_t the_sig_pi_numberOfPixelLayers = -99;
   Int_t the_sig_pi_qualityIndex = -99;
   Int_t the_sig_pi_highPurityFlag = -99;
   Int_t the_sig_pi_packedcandhashighpurity = -99;
   Int_t the_sig_pi_matchedtomuon_loose = -99;
   Int_t the_sig_pi_matchedtomuon_medium = -99;
   Int_t the_sig_pi_matchedtomuon_tight = -99;

   Float_t the_sig_trgmu_mu_mass = -99.;
   Float_t the_sig_trgmu_mu_pt = -99.;
   Float_t the_sig_trgmu_pi_mass = -99.;
   Float_t the_sig_trgmu_pi_pt = -99.;
   //Float_t the_sig_dimu_mass = -99.;
   //Float_t the_sig_dimu_pt = -99.;
   Float_t the_sig_dimu_lxy = -99.;
   Float_t the_sig_dimu_lxyz = -99.;
   Float_t the_sig_dimu_vxdiff = -99.;
   Float_t the_sig_dimu_vydiff = -99.;
   Float_t the_sig_dimu_vzdiff = -99.;

   Float_t the_sig_cos_theta_star_pion = -99.;
   Float_t the_sig_cos_theta_star_muon = -99.;
   Float_t the_sig_cos_theta_star_sum = -99.;

   Float_t the_sig_px_diff_hnl_daughters_lab = -99.;
   Float_t the_sig_py_diff_hnl_daughters_lab = -99.;
   Float_t the_sig_pz_diff_hnl_daughters_lab = -99.;
   Float_t the_sig_energy_diff_prefithnl_daughters_lab = -99.;
   Float_t the_sig_px_diff_prefithnl_daughters_lab = -99.;
   Float_t the_sig_py_diff_prefithnl_daughters_lab = -99.;
   Float_t the_sig_pz_diff_prefithnl_daughters_lab = -99.;

   Float_t the_sig_deltar_mu_pi = -99.;
   Float_t the_sig_deltar_trgmu_hnl = -99.;
   Float_t the_sig_deltar_trgmu_mu = -99.;
   Float_t the_sig_deltar_trgmu_pi = -99.;
   Float_t the_sig_deltaeta_mu_pi = -99.;
   Float_t the_sig_deltaeta_trgmu_hnl = -99.;
   Float_t the_sig_deltaeta_trgmu_mu = -99.;
   Float_t the_sig_deltaeta_trgmu_pi = -99.;
   Float_t the_sig_deltaphi_mu_pi = -99.;
   Float_t the_sig_deltaphi_trgmu_hnl = -99.;
   Float_t the_sig_deltaphi_trgmu_mu = -99.;
   Float_t the_sig_deltaphi_trgmu_pi = -99.;

   Float_t the_sig_deltae_pi_fit_pi = -99.;
   Float_t the_sig_deltae_mu_fit_mu = -99.;
   Float_t the_sig_deltae_hnl_fit_hnl = -99.;
   Float_t the_sig_deltapt_pi_fit_pi = -99.;
   Float_t the_sig_deltapt_mu_fit_mu = -99.;
   Float_t the_sig_deltapt_hnl_fit_hnl = -99.;
   Float_t the_sig_deltapx_pi_fit_pi = -99.;
   Float_t the_sig_deltapx_mu_fit_mu = -99.;
   Float_t the_sig_deltapx_hnl_fit_hnl = -99.;
   Float_t the_sig_deltapy_pi_fit_pi = -99.;
   Float_t the_sig_deltapy_mu_fit_mu = -99.;
   Float_t the_sig_deltapy_hnl_fit_hnl = -99.;
   Float_t the_sig_deltapz_pi_fit_pi = -99.;
   Float_t the_sig_deltapz_mu_fit_mu = -99.;
   Float_t the_sig_deltapz_hnl_fit_hnl = -99.;
   Float_t the_sig_deltaeta_pi_fit_pi = -99.;
   Float_t the_sig_deltaeta_mu_fit_mu = -99.;
   Float_t the_sig_deltaeta_hnl_fit_hnl = -99.;
   Float_t the_sig_deltaphi_pi_fit_pi = -99.;
   Float_t the_sig_deltaphi_mu_fit_mu = -99.;
   Float_t the_sig_deltaphi_hnl_fit_hnl = -99.;

   Float_t the_sig_sv_chi2 = -99.;
   Float_t the_sig_sv_lxy = -99.;
   Float_t the_sig_sv_lxyz = -99.;
   Float_t the_sig_sv_pv_lxy = -99.;
   Float_t the_sig_sv_pv_lxyz = -99.;
   Float_t the_sig_sv_lxysig = -99.;
   Float_t the_sig_sv_prob = -99.;
   Float_t the_sig_sv_x = -99.;
   Float_t the_sig_sv_y = -99.;
   Float_t the_sig_sv_z = -99.;

   Float_t the_sig_pi_mu_vzdiff = -99.;

   Int_t the_sig_ismatched = -99;
   Int_t the_sig_trgmu_ismatched = -99;
   Int_t the_sig_mu_ismatched = -99;
   Int_t the_sig_pi_ismatched = -99;
   Float_t the_sig_mupi_mass_reco_gen_reldiff = -99.;
   Float_t the_sig_lxy_reco_gen_reldiff = -99.;

   // weights
   Float_t the_sig_weight_hlt_A1 = -99.;
   Float_t the_sig_weight_hlt_A1_6 = -99.;
   Float_t the_sig_weight_hlt_HLT_Mu9_IP6_A1_6 = -99.;
   Float_t the_sig_weight_hlt_A1_6_B1 = -99.;
   Float_t the_sig_weight_pu_qcd_A = -99.;
   Float_t the_sig_weight_pu_qcd_B = -99.;
   Float_t the_sig_weight_pu_qcd_C = -99.;
   Float_t the_sig_weight_pu_qcd_D = -99.;
   Float_t the_sig_weight_pu_qcd_tot = -99.;
   Float_t the_sig_weight_pu_sig_A = -99.;
   Float_t the_sig_weight_pu_sig_B = -99.;
   Float_t the_sig_weight_pu_sig_C = -99.;
   Float_t the_sig_weight_pu_sig_D = -99.;
   Float_t the_sig_weight_pu_sig_tot = -99.;
   
   // these two are from matching information
   Float_t the_gen_trgmu_mu_lxy = -99.;
   Float_t the_gen_trgmu_mu_lxyz = -99.;

   // these have the gen information, independently on the matching
   Float_t the_gen_b_pt = -99.;
   Float_t the_gen_b_eta = -99.;
   Float_t the_gen_b_phi = -99.;
   Float_t the_gen_b_mass = -99.;
   Int_t the_gen_b_pdgid = -99;
   Float_t the_gen_hnl_ct = -99.;
   Float_t the_gen_hnl_lxy = -99.;
   Float_t the_gen_hnl_lxyz = -99.;
   Float_t the_gen_hnl_pt = -99.;
   Float_t the_gen_hnl_eta = -99.;
   Float_t the_gen_hnl_phi = -99.;
   Float_t the_gen_hnl_mass = -99.;
   Float_t the_gen_hnl_vx = -99.;
   Float_t the_gen_hnl_vy = -99.;
   Float_t the_gen_hnl_vz = -99.;
   Float_t the_gen_trgmu_pt = -99.;
   Float_t the_gen_trgmu_eta = -99.;
   Float_t the_gen_trgmu_phi = -99.;
   Float_t the_gen_trgmu_vx = -99.;
   Float_t the_gen_trgmu_vy = -99.;
   Float_t the_gen_trgmu_vz = -99.;
   Float_t the_gen_mu_pt = -99.;
   Float_t the_gen_mu_eta = -99.;
   Float_t the_gen_mu_phi = -99.;
   Float_t the_gen_mu_vx = -99.;
   Float_t the_gen_mu_vy = -99.;
   Float_t the_gen_mu_vz = -99.;
   Float_t the_gen_pi_pt = -99.;
   Float_t the_gen_pi_eta = -99.;
   Float_t the_gen_pi_phi = -99.;
   Float_t the_gen_pi_vx = -99.;
   Float_t the_gen_pi_vy = -99.;
   Float_t the_gen_pi_vz = -99.;

   // histograms
   TH1F* sighist_ncand_perevent;
   TH1F* sighist_ncand_matched_perevent;
   TH1F* sighist_selection_efficiency_hnlpt_allevents;
   TH1F* sighist_selection_efficiency_hnlpt_eventswithmultcands;
   TH1F* sighist_selection_efficiency_bpt_allevents;
   TH1F* sighist_selection_efficiency_bpt_eventswithmultcands;
   TH1F* sighist_selection_efficiency_trgmupt_allevents;
   TH1F* sighist_selection_efficiency_trgmupt_eventswithmultcands;
   TH1F* sighist_selection_efficiency_pipt_allevents;
   TH1F* sighist_selection_efficiency_pipt_eventswithmultcands;
   TH1F* sighist_selection_efficiency_svprob_allevents;
   TH1F* sighist_selection_efficiency_svprob_eventswithmultcands;
   TH1F* sighist_selection_efficiency_svchi2_allevents;
   TH1F* sighist_selection_efficiency_svchi2_eventswithmultcands;
   TH1F* sighist_selection_efficiency_cos2d_allevents;
   TH1F* sighist_selection_efficiency_cos2d_eventswithmultcands;
   TH1F* sighist_selection_efficiency_hnliso4_allevents;
   TH1F* sighist_selection_efficiency_hnliso4_eventswithmultcands;
   TH1F* sighist_selection_efficiency_dr_allevents;
   TH1F* sighist_selection_efficiency_dr_eventswithmultcands;

   ClassDef(BToMuMuPiDumper,0);
};

#endif

#ifdef BToMuMuPiDumper_cxx
void BToMuMuPiDumper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t BToMuMuPiDumper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef BToMuMuPiDumper_cxx
