//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 23 15:15:57 2021 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: bparknano_nj89.root
//////////////////////////////////////////////////////////

#ifndef BackgroundSources_h
#define BackgroundSources_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class BackgroundSources : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> nBToMuMuPi = {fReader, "nBToMuMuPi"};
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
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso03 = {fReader, "BToMuMuPi_sel_mu_iso03"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso03_close = {fReader, "BToMuMuPi_sel_mu_iso03_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso03_rel_close = {fReader, "BToMuMuPi_sel_mu_iso03_rel_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso04 = {fReader, "BToMuMuPi_sel_mu_iso04"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso04_close = {fReader, "BToMuMuPi_sel_mu_iso04_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_sel_mu_iso04_rel_close = {fReader, "BToMuMuPi_sel_mu_iso04_rel_close"};
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
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_eta = {fReader, "BToMuMuPi_trg_mu_eta"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso03 = {fReader, "BToMuMuPi_trg_mu_iso03"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso03_close = {fReader, "BToMuMuPi_trg_mu_iso03_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso03_rel_close = {fReader, "BToMuMuPi_trg_mu_iso03_rel_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso04 = {fReader, "BToMuMuPi_trg_mu_iso04"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso04_close = {fReader, "BToMuMuPi_trg_mu_iso04_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_iso04_rel_close = {fReader, "BToMuMuPi_trg_mu_iso04_rel_close"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_phi = {fReader, "BToMuMuPi_trg_mu_phi"};
   TTreeReaderArray<Float_t> BToMuMuPi_trg_mu_pt = {fReader, "BToMuMuPi_trg_mu_pt"};
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


   BackgroundSources(TTree * /*tree*/ =0) { }
   virtual ~BackgroundSources() { }
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
   std::vector<pair<int, float>> getGenCandidateMap(float, float, float, float, float, int) const;
   virtual std::vector<int> getAncestor(int, std::vector<int>) const;
   virtual Int_t getBmesonIdx(int particle_idx, int b_meson_idx) const;
   virtual Bool_t checkCommonAncestor(int, int, int) const;
   virtual Bool_t checkCommonAncestor(int, int) const;

   // output file
   TFile* my_file;  

   // some branches will be added only if sample is MC
   Bool_t isMC;

   // trees to fill
   TTree* sources;
      
   std::vector<int> B_mesons = {511, 521, 531, 541};
   std::vector<int> D_mesons = {411, 421, 431, 441};
   std::vector<int> excitedD_mesons = {10411, 10421, 413, 423, 415, 425, 10431, 433, 435};
   std::vector<int> K_mesons = {311, 321};
   
   Float_t hnl_mass = -99.;
   Float_t hnl_pt = -99.;
   Int_t hnl_charge = -99;
   Int_t trgmu_softid = -99;
   Int_t mu_looseid = -99;
   Int_t pi_packedcandhashighpurity = -99;
   Int_t trgmu_charge = -99;
   Int_t mu_charge = -99;
   Float_t trgmu_mu_mass = -99.;
   Float_t trgmu_pi_mass = -99.;
   Float_t mu0_pt = -99.;
   Float_t mu0_eta = -99.;
   Float_t mu_pt = -99.;
   Float_t mu_eta = -99.;
   Float_t pi_pt = -99.;
   Float_t pi_eta = -99.;
   Float_t sv_lxy = -99.;
   Float_t sv_lxyz = -99.;
   Float_t sv_lxysig = -99.;
   Float_t mu0_dxysig = -99.;
   Float_t mu_dxysig = -99.;
   Float_t pi_dxysig = -99.;
   Float_t hnl_cos2d = -99.;

   Int_t mu0_ismatched = -99;
   Int_t mu0_isfake = -99;
   Float_t mu0_genmatching_deltar = -99.;
   Float_t mu0_genmatching_deltaptrel = -99.;

   Int_t mu0_directBdecay = -99;
   Int_t mu0_BtoDdecay = -99;
   Int_t mu0_BtoDtoKdecay = -99;
   Int_t mu0_Btotaudecay = -99;
   Int_t mu0_BtoJPsidecay = -99;
   Int_t mu0_BtoDexciteddecay = -99;
   Int_t mu0_nonBdecay = -99;
   Int_t mu0_other = -99;

   Int_t mu_ismatched = -99;
   Int_t mu_isfake = -99;
   Float_t mu_genmatching_deltar = -99.;
   Float_t mu_genmatching_deltaptrel = -99.;

   Int_t mu_directBdecay = -99;
   Int_t mu_BtoDdecay = -99;
   Int_t mu_BtoDtoKdecay = -99;
   Int_t mu_Btotaudecay = -99;
   Int_t mu_BtoJPsidecay = -99;
   Int_t mu_BtoDexciteddecay = -99;
   Int_t mu_nonBdecay = -99;
   Int_t mu_other = -99;

   Int_t pi_ismatched = -99;
   Int_t pi_isfake = -99;
   Float_t pi_genmatching_deltar = -99.;
   Float_t pi_genmatching_deltaptrel = -99.;

   Int_t pi_directBdecay = -99;
   Int_t pi_BtoDdecay = -99;
   Int_t pi_BtoDtoKdecay = -99;
   Int_t pi_Btotaudecay = -99;
   Int_t pi_BtoJPsidecay = -99;
   Int_t pi_BtoDexciteddecay = -99;
   Int_t pi_nonBdecay = -99;
   Int_t pi_other = -99;

   Int_t cand_isnotmatched = -99;
   Int_t cand_isfullmatched = -99;
   Int_t cand_ispartialmatched = -99;
   Int_t cand_mu0mu_samemother = -99;
   Int_t cand_mu0pi_samemother = -99;
   Int_t cand_mupi_samemother = -99;
   Int_t cand_mumupi_samemother = -99;
   
   ClassDef(BackgroundSources,0);
};

#endif

#ifdef BackgroundSources_cxx
void BackgroundSources::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t BackgroundSources::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef BackgroundSources_cxx
