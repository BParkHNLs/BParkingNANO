#define TagAndProbeDumper_cxx
// The class definition in TagAndProbeDumper.h has been generated automatically
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
// root> T->Process("TagAndProbeDumper.C")
// root> T->Process("TagAndProbeDumper.C","some options")
// root> T->Process("TagAndProbeDumper.C+")
//


#include "TagAndProbeDumper.h"
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include "utils.C"

void TagAndProbeDumper::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  cout << " --------------------------" << endl;
  cout << "   Tag And Probe Dumper    " << endl;
  cout << " --------------------------" << endl;
}

void TagAndProbeDumper::SlaveBegin(TTree * /*tree*/)
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

  if(isMC){
    Pileup_nTrueInt = {fReader, "Pileup_nTrueInt"};
  }

  tree = new TTree("tree", "tree");
  tree->Branch("pt", &the_pt);
  tree->Branch("eta", &the_eta);
  tree->Branch("phi", &the_phi);
  tree->Branch("mass", &the_mass);
  tree->Branch("cos2d", &the_cos2d);
  tree->Branch("deltar", &the_deltar);
  tree->Branch("lxy", &the_lxy);
  tree->Branch("lxy_sig", &the_lxy_sig);
  tree->Branch("ismatched", &the_ismatched);

  tree->Branch("tag_pt", &the_tag_pt);
  tree->Branch("tag_eta", &the_tag_eta);
  tree->Branch("tag_phi", &the_tag_phi);
  tree->Branch("tag_dxy", &the_tag_dxy);
  tree->Branch("tag_dxy_bs", &the_tag_dxy_bs);
  tree->Branch("tag_dz", &the_tag_dz);
  tree->Branch("tag_dxy_sig", &the_tag_dxy_sig);
  tree->Branch("tag_dxy_sig_bs", &the_tag_dxy_sig_bs);
  tree->Branch("tag_dxy_sig_bs_rdst", &the_tag_dxy_sig_bs_rdst);
  tree->Branch("tag_dz_sig", &the_tag_dz_sig);
  tree->Branch("tag_isloose", &the_tag_isloose);
  tree->Branch("tag_issoft", &the_tag_issoft);
  tree->Branch("tag_ismedium", &the_tag_ismedium);
  tree->Branch("tag_istight", &the_tag_istight);
  tree->Branch("tag_fired_HLT_Mu7_IP4", &the_tag_fired_HLT_Mu7_IP4);
  tree->Branch("tag_fired_HLT_Mu8_IP6", &the_tag_fired_HLT_Mu8_IP6);
  tree->Branch("tag_fired_HLT_Mu8_IP5", &the_tag_fired_HLT_Mu8_IP5);
  tree->Branch("tag_fired_HLT_Mu8_IP3", &the_tag_fired_HLT_Mu8_IP3);
  tree->Branch("tag_fired_HLT_Mu8p5_IP3p5", &the_tag_fired_HLT_Mu8p5_IP3p5);
  tree->Branch("tag_fired_HLT_Mu9_IP6", &the_tag_fired_HLT_Mu9_IP6);
  tree->Branch("tag_fired_HLT_Mu9_IP5", &the_tag_fired_HLT_Mu9_IP5);
  tree->Branch("tag_fired_HLT_Mu9_IP4", &the_tag_fired_HLT_Mu9_IP4);
  tree->Branch("tag_fired_HLT_Mu10p5_IP3p5", &the_tag_fired_HLT_Mu10p5_IP3p5);
  tree->Branch("tag_fired_HLT_Mu12_IP6", &the_tag_fired_HLT_Mu12_IP6);
  tree->Branch("tag_fired_HLT_Mu8_v1", &the_tag_fired_HLT_Mu8_v1);
  tree->Branch("tag_fired_HLT_Mu8_v12", &the_tag_fired_HLT_Mu8_v12);
  tree->Branch("tag_fired_HLT_Mu7p5_Track7_Jpsi_v11", &the_tag_fired_HLT_Mu7p5_Track7_Jpsi_v11);
  tree->Branch("tag_fired_HLT_L2Mu23NoVtx_2Cha_v1", &the_tag_fired_HLT_L2Mu23NoVtx_2Cha_v1);
  tree->Branch("tag_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1", &the_tag_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1);
  tree->Branch("tag_fired_DST_DoubleMu1_noVtx_CaloScouting_v2", &the_tag_fired_DST_DoubleMu1_noVtx_CaloScouting_v2);
  tree->Branch("tag_fired_DST_DoubleMu3_noVtx_CaloScouting_v6", &the_tag_fired_DST_DoubleMu3_noVtx_CaloScouting_v6);
  tree->Branch("tag_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3", &the_tag_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3);
  tree->Branch("tag_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13", &the_tag_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13);
  tree->Branch("tag_prescale_HLT_Mu7_IP4", &the_tag_prescale_HLT_Mu7_IP4);
  tree->Branch("tag_prescale_HLT_Mu8_IP6", &the_tag_prescale_HLT_Mu8_IP6);
  tree->Branch("tag_prescale_HLT_Mu8_IP5", &the_tag_prescale_HLT_Mu8_IP5);
  tree->Branch("tag_prescale_HLT_Mu8_IP3", &the_tag_prescale_HLT_Mu8_IP3);
  tree->Branch("tag_prescale_HLT_Mu8p5_IP3p5", &the_tag_prescale_HLT_Mu8p5_IP3p5);
  tree->Branch("tag_prescale_HLT_Mu9_IP6", &the_tag_prescale_HLT_Mu9_IP6);
  tree->Branch("tag_prescale_HLT_Mu9_IP5", &the_tag_prescale_HLT_Mu9_IP5);
  tree->Branch("tag_prescale_HLT_Mu9_IP4", &the_tag_prescale_HLT_Mu9_IP4);
  tree->Branch("tag_prescale_HLT_Mu10p5_IP3p5", &the_tag_prescale_HLT_Mu10p5_IP3p5);
  tree->Branch("tag_prescale_HLT_Mu12_IP6", &the_tag_prescale_HLT_Mu12_IP6);
  tree->Branch("tag_prescale_HLT_Mu8_v1", &the_tag_prescale_HLT_Mu8_v1);
  tree->Branch("tag_prescale_HLT_Mu8_v12", &the_tag_prescale_HLT_Mu8_v12);
  tree->Branch("tag_prescale_HLT_Mu7p5_Track7_Jpsi_v11", &the_tag_prescale_HLT_Mu7p5_Track7_Jpsi_v11);
  tree->Branch("tag_prescale_HLT_L2Mu23NoVtx_2Cha_v1", &the_tag_prescale_HLT_L2Mu23NoVtx_2Cha_v1);
  tree->Branch("tag_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1", &the_tag_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1);
  tree->Branch("tag_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2", &the_tag_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2);
  tree->Branch("tag_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6", &the_tag_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6);
  tree->Branch("tag_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3", &the_tag_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3);
  tree->Branch("tag_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13", &the_tag_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13);


  tree->Branch("probe_pt", &the_probe_pt);
  tree->Branch("probe_eta", &the_probe_eta);
  tree->Branch("probe_phi", &the_probe_phi);
  tree->Branch("probe_dxy", &the_probe_dxy);
  tree->Branch("probe_dxy_bs", &the_probe_dxy_bs);
  tree->Branch("probe_dxy_bs_uncorrected", &the_probe_dxy_bs_uncorrected);
  tree->Branch("weight_dxy_bs", &the_weight_dxy_bs);
  tree->Branch("probe_dz", &the_probe_dz);
  tree->Branch("probe_dxy_sig", &the_probe_dxy_sig);

  //tree->Branch("smeared_corr_0p01", &the_smeared_corr_0p01);
  //tree->Branch("smeared_corr_0p02", &the_smeared_corr_0p02);
  //tree->Branch("smeared_corr_0p03", &the_smeared_corr_0p03);
  //tree->Branch("smeared_corr_0p04", &the_smeared_corr_0p04);
  //tree->Branch("smeared_corr_0p05", &the_smeared_corr_0p05);
  //tree->Branch("smeared_corr_0p06", &the_smeared_corr_0p06);
  //tree->Branch("smeared_corr_0p07", &the_smeared_corr_0p07);
  //tree->Branch("smeared_corr_0p08", &the_smeared_corr_0p08);
  //tree->Branch("smeared_corr_0p09", &the_smeared_corr_0p09);
  //tree->Branch("smeared_corr_0p1", &the_smeared_corr_0p1);
  //tree->Branch("smeared_corr_0p15", &the_smeared_corr_0p15);
  //tree->Branch("smeared_corr_0p2", &the_smeared_corr_0p2);
  //tree->Branch("smeared_corr_0p3", &the_smeared_corr_0p3);
  //tree->Branch("smeared_corr_0p35", &the_smeared_corr_0p35);
  //tree->Branch("smeared_corr_0p4", &the_smeared_corr_0p4);
  //tree->Branch("smeared_corr_0p45", &the_smeared_corr_0p45);
  //tree->Branch("smeared_corr_0p5", &the_smeared_corr_0p5);

  //tree->Branch("double_gauss_corr_1", &the_double_gauss_corr_1);
  //tree->Branch("double_gauss_corr_2", &the_double_gauss_corr_2);

  tree->Branch("probe_dxy_sig_bs_uncorrected", &the_probe_dxy_sig_bs_uncorrected);
  tree->Branch("probe_dxy_sig_bs", &the_probe_dxy_sig_bs);
  tree->Branch("probe_dxy_sig_bs_corrected_linearscale", &the_probe_dxy_sig_bs_corrected_linearscale);
  tree->Branch("weight_dxy_sig_bs", &the_weight_dxy_sig_bs);

  tree->Branch("probe_dxy_sig_bs_rdst", &the_probe_dxy_sig_bs_rdst);
  tree->Branch("probe_dz_sig", &the_probe_dz_sig);
  tree->Branch("probe_isloose", &the_probe_isloose);
  tree->Branch("probe_issoft", &the_probe_issoft);
  tree->Branch("probe_ismedium", &the_probe_ismedium);
  tree->Branch("probe_istight", &the_probe_istight);
  tree->Branch("probe_fired_HLT_Mu7_IP4", &the_probe_fired_HLT_Mu7_IP4);
  tree->Branch("probe_fired_HLT_Mu8_IP6", &the_probe_fired_HLT_Mu8_IP6);
  tree->Branch("probe_fired_HLT_Mu8_IP5", &the_probe_fired_HLT_Mu8_IP5);
  tree->Branch("probe_fired_HLT_Mu8_IP3", &the_probe_fired_HLT_Mu8_IP3);
  tree->Branch("probe_fired_HLT_Mu8p5_IP3p5", &the_probe_fired_HLT_Mu8p5_IP3p5);
  tree->Branch("probe_fired_HLT_Mu9_IP6", &the_probe_fired_HLT_Mu9_IP6);
  tree->Branch("probe_fired_HLT_Mu9_IP5", &the_probe_fired_HLT_Mu9_IP5);
  tree->Branch("probe_fired_HLT_Mu9_IP4", &the_probe_fired_HLT_Mu9_IP4);
  tree->Branch("probe_fired_HLT_Mu10p5_IP3p5", &the_probe_fired_HLT_Mu10p5_IP3p5);
  tree->Branch("probe_fired_HLT_Mu12_IP6", &the_probe_fired_HLT_Mu12_IP6);
  tree->Branch("probe_fired_HLT_Mu8_v1", &the_probe_fired_HLT_Mu8_v1);
  tree->Branch("probe_fired_HLT_Mu8_v12", &the_probe_fired_HLT_Mu8_v12);
  tree->Branch("probe_fired_HLT_Mu7p5_Track7_Jpsi_v11", &the_probe_fired_HLT_Mu7p5_Track7_Jpsi_v11);
  tree->Branch("probe_fired_HLT_L2Mu23NoVtx_2Cha_v1", &the_probe_fired_HLT_L2Mu23NoVtx_2Cha_v1);
  tree->Branch("probe_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1", &the_probe_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1);
  tree->Branch("probe_fired_DST_DoubleMu1_noVtx_CaloScouting_v2", &the_probe_fired_DST_DoubleMu1_noVtx_CaloScouting_v2);
  tree->Branch("probe_fired_DST_DoubleMu3_noVtx_CaloScouting_v6", &the_probe_fired_DST_DoubleMu3_noVtx_CaloScouting_v6);
  tree->Branch("probe_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3", &the_probe_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3);
  tree->Branch("probe_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13", &the_probe_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13);
  tree->Branch("probe_prescale_HLT_Mu7_IP4", &the_probe_prescale_HLT_Mu7_IP4);
  tree->Branch("probe_prescale_HLT_Mu8_IP6", &the_probe_prescale_HLT_Mu8_IP6);
  tree->Branch("probe_prescale_HLT_Mu8_IP5", &the_probe_prescale_HLT_Mu8_IP5);
  tree->Branch("probe_prescale_HLT_Mu8_IP3", &the_probe_prescale_HLT_Mu8_IP3);
  tree->Branch("probe_prescale_HLT_Mu8p5_IP3p5", &the_probe_prescale_HLT_Mu8p5_IP3p5);
  tree->Branch("probe_prescale_HLT_Mu9_IP6", &the_probe_prescale_HLT_Mu9_IP6);
  tree->Branch("probe_prescale_HLT_Mu9_IP5", &the_probe_prescale_HLT_Mu9_IP5);
  tree->Branch("probe_prescale_HLT_Mu9_IP4", &the_probe_prescale_HLT_Mu9_IP4);
  tree->Branch("probe_prescale_HLT_Mu10p5_IP3p5", &the_probe_prescale_HLT_Mu10p5_IP3p5);
  tree->Branch("probe_prescale_HLT_Mu12_IP6", &the_probe_prescale_HLT_Mu12_IP6);
  tree->Branch("probe_prescale_HLT_Mu8_v1", &the_probe_prescale_HLT_Mu8_v1);
  tree->Branch("probe_prescale_HLT_Mu8_v12", &the_probe_prescale_HLT_Mu8_v12);
  tree->Branch("probe_prescale_HLT_Mu7p5_Track7_Jpsi_v11", &the_probe_prescale_HLT_Mu7p5_Track7_Jpsi_v11);
  tree->Branch("probe_prescale_HLT_L2Mu23NoVtx_2Cha_v1", &the_probe_prescale_HLT_L2Mu23NoVtx_2Cha_v1);
  tree->Branch("probe_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1", &the_probe_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1);
  tree->Branch("probe_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2", &the_probe_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2);
  tree->Branch("probe_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6", &the_probe_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6);
  tree->Branch("probe_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3", &the_probe_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3);
  tree->Branch("probe_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13", &the_probe_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13);
  tree->Branch("probe_fired_BParkingHLT", &the_probe_fired_BParkingHLT);

  tree->Branch("pv_npvs", &the_pv_npvs);

  tree->Branch("weight_hlt_A1", &the_weight_hlt_A1);
  tree->Branch("weight_hlt_A1_6", &the_weight_hlt_A1_6);
  tree->Branch("weight_pu_A", &the_weight_pu_A);
  tree->Branch("weight_pu_B", &the_weight_pu_B);
  tree->Branch("weight_pu_C", &the_weight_pu_C);
  tree->Branch("weight_pu_D", &the_weight_pu_D);
  tree->Branch("weight_pu_tot", &the_weight_pu_tot);

  //double_gauss_1 = new TF1("double_gauss_1","[0]/sqrt(2.*TMath::Pi())/[2]*exp(-(x-[1])*(x-[1])/2./[2]/[2])+(1-[0])/sqrt(2.*TMath::Pi())/[4]*exp(-(x-[3])*(x-[3])/2./[4]/[4])", 0.9, 1.3);
  //double_gauss_1->SetParameter(0, 0.62);
  //double_gauss_1->SetParameter(1, 1.042);
  //double_gauss_1->SetParameter(2, 0.034);
  //double_gauss_1->SetParameter(3, 1.144);
  //double_gauss_1->SetParameter(4, 0.026);

  //double_gauss_2 = new TF1("double_gauss_2","[0]/sqrt(2.*TMath::Pi())/[2]*exp(-(x-[1])*(x-[1])/2./[2]/[2])+(1-[0])/sqrt(2.*TMath::Pi())/[4]*exp(-(x-[3])*(x-[3])/2./[4]/[4])", 0.9, 1.3);
  //double_gauss_2->SetParameter(0, 0.4);
  //double_gauss_2->SetParameter(1, 1.03);
  //double_gauss_2->SetParameter(2, 0.030);
  //double_gauss_2->SetParameter(3, 1.13);
  //double_gauss_2->SetParameter(4, 0.05);

  //gauss_function_0p01 = new TF1("gauss_function_0p01", "gaus",-10, 10);
  //gauss_function_0p01->SetParameters(1, 1., 0.01);
  //gauss_function_0p02 = new TF1("gauss_function_0p02", "gaus",-10, 10);
  //gauss_function_0p02->SetParameters(1, 1., 0.02);
  gauss_function_0p03 = new TF1("gauss_function_0p03", "gaus",-10, 10);
  gauss_function_0p03->SetParameters(1, 1., 0.03);
  //gauss_function_0p04 = new TF1("gauss_function_0p04", "gaus",-10, 10);
  //gauss_function_0p04->SetParameters(1, 1., 0.04);
  //gauss_function_0p05 = new TF1("gauss_function_0p05", "gaus",-10, 10);
  //gauss_function_0p05->SetParameters(1, 1., 0.05);
  //gauss_function_0p06 = new TF1("gauss_function_0p06", "gaus",-10, 10);
  //gauss_function_0p06->SetParameters(1, 1., 0.06);
  //gauss_function_0p07 = new TF1("gauss_function_0p07", "gaus",-10, 10);
  //gauss_function_0p07->SetParameters(1, 1., 0.07);
  //gauss_function_0p08 = new TF1("gauss_function_0p08", "gaus",-10, 10);
  //gauss_function_0p08->SetParameters(1, 1., 0.08);
  //gauss_function_0p09 = new TF1("gauss_function_0p09", "gaus",-10, 10);
  //gauss_function_0p09->SetParameters(1, 1., 0.09);
  //gauss_function_0p1 = new TF1("gauss_function_0p1", "gaus",-10, 10);
  //gauss_function_0p1->SetParameters(1, 1., 0.1);
  //gauss_function_0p15 = new TF1("gauss_function_0p15", "gaus",-10, 10);
  //gauss_function_0p15->SetParameters(1, 1., 0.15);
  //gauss_function_0p2 = new TF1("gauss_function_0p2", "gaus",-10, 10);
  //gauss_function_0p2->SetParameters(1, 1., 0.2);
  //gauss_function_0p25 = new TF1("gauss_function_0p25", "gaus",-10, 10);
  //gauss_function_0p25->SetParameters(1, 1., 0.25);
  //gauss_function_0p3 = new TF1("gauss_function_0p3", "gaus",-10, 10);
  //gauss_function_0p3->SetParameters(1, 1., 0.3);
  //gauss_function_0p35 = new TF1("gauss_function_0p35", "gaus",-10, 10);
  //gauss_function_0p35->SetParameters(1, 1., 0.35);
  //gauss_function_0p4 = new TF1("gauss_function_0p4", "gaus",-10, 10);
  //gauss_function_0p4->SetParameters(1, 1., 0.4);
  //gauss_function_0p45 = new TF1("gauss_function_0p45", "gaus",-10, 10);
  //gauss_function_0p45->SetParameters(1, 1., 0.45);
  //gauss_function_0p5 = new TF1("gauss_function_0p5", "gaus",-10, 10);
  //gauss_function_0p5->SetParameters(1, 1., 0.5);

}

Bool_t TagAndProbeDumper::Process(Long64_t entry)
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
  //if(entry > 10000) return false;
  //std::cout << entry << std::endl;

  // for data, we skip the event in case it doesn't pass the lumi mask
  if(!isMC && lumiMask(*run, *luminosityBlock) == false) return false;

  // initial strategy
  //if(Muon_fired_DST_DoubleMu1_noVtx_CaloScouting_v2[JPsiToMuMu_lep1_idx[0]] != 1) return false;
  //if(Muon_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2[JPsiToMuMu_lep1_idx[0]] != 1) return false;

  //if(Muon_fired_DST_DoubleMu1_noVtx_CaloScouting_v2[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_DST_DoubleMu3_noVtx_CaloScouting_v6[JPsiToMuMu_lep1_idx[0]] != 1) return false;
  
  // the following was used to compute SF on period A only
  //if(Muon_fired_HLT_Mu9_IP6[JPsiToMuMu_lep1_idx[0]] != 1) return false;
  //if(Muon_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2[JPsiToMuMu_lep1_idx[0]] == -1) return false;
  
  // requirement of the tag muon to fire any BParking HLT line
  //if(Muon_isTriggeringBPark[JPsiToMuMu_lep1_idx[0]] != 1) return false;
  //if(Muon_fired_HLT_Mu12_IP6[JPsiToMuMu_lep1_idx[0]] != 1) return false;
  if(Muon_fired_HLT_Mu9_IP6[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu12_IP6[JPsiToMuMu_lep1_idx[0]] != 1) return false;
  //if(Muon_fired_HLT_Mu10p5_IP3p5[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu8_IP3[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu12_IP6[JPsiToMuMu_lep1_idx[0]] != 1) return false;
  //if(Muon_fired_HLT_Mu10p5_IP3p5[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu12_IP6[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu7_IP4[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu8_IP3[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu8_IP5[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu8_IP6[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu8p5_IP3p5[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu9_IP4[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu9_IP5[JPsiToMuMu_lep1_idx[0]] != 1 && Muon_fired_HLT_Mu9_IP6[JPsiToMuMu_lep1_idx[0]] != 1) return false;
    
  // number of candidates in the event
  UInt_t nCand = *nJPsiToMuMu; 

  if(nCand > 0){

    // keep probe muon with given id
    //if(Muon_looseId[JPsiToMuMu_lep2_idx[0]] != 1) return false;

    // only keep matched candidates - although matching efficiency close to 1
    if(isMC && JPsiToMuMu_isMatched[0] != 1) return false;

    // by default, take the first candidate (possible since <permille events have more than one candidate per event)
    the_pt = JPsiToMuMu_pt[0];
    the_eta = fabs(JPsiToMuMu_eta[0]);
    the_phi = JPsiToMuMu_phi[0];
    the_mass = JPsiToMuMu_mass[0];
    the_cos2d = JPsiToMuMu_cos2D[0];
    the_deltar = JPsiToMuMu_deltaR[0];
    the_lxy = JPsiToMuMu_lxy[0];
    the_lxy_sig = JPsiToMuMu_lxy_sig[0];
    the_ismatched = fabs(JPsiToMuMu_isMatched[0]);

    the_tag_pt = JPsiToMuMu_lep1_pt[0];
    the_tag_eta = fabs(JPsiToMuMu_lep1_eta[0]);
    the_tag_phi = JPsiToMuMu_lep1_phi[0];
    the_tag_dxy = fabs(Muon_dxy[JPsiToMuMu_lep1_idx[0]]);
    the_tag_dxy_bs = fabs(Muon_dxy_BS[JPsiToMuMu_lep1_idx[0]]);
    the_tag_dz = fabs(Muon_dz[JPsiToMuMu_lep1_idx[0]]);
    the_tag_dxy_sig = fabs(Muon_dxyS[JPsiToMuMu_lep1_idx[0]]);
    the_tag_dxy_sig_bs = fabs(Muon_dxyS_BS[JPsiToMuMu_lep1_idx[0]]);
    the_tag_dxy_sig_bs_rdst = fabs(Muon_dxyS_BS_alaRdst[JPsiToMuMu_lep1_idx[0]]);
    the_tag_dz_sig = fabs(Muon_dzS[JPsiToMuMu_lep1_idx[0]]);
    the_tag_isloose = Muon_looseId[JPsiToMuMu_lep1_idx[0]];
    the_tag_issoft = Muon_softId[JPsiToMuMu_lep1_idx[0]];
    the_tag_ismedium = Muon_mediumId[JPsiToMuMu_lep1_idx[0]];
    the_tag_istight = Muon_tightId[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu7_IP4 = Muon_fired_HLT_Mu7_IP4[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu8_IP6 = Muon_fired_HLT_Mu8_IP6[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu8_IP5 = Muon_fired_HLT_Mu8_IP5[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu8_IP3 = Muon_fired_HLT_Mu8_IP3[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu8p5_IP3p5 = Muon_fired_HLT_Mu8p5_IP3p5[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu9_IP6 = Muon_fired_HLT_Mu9_IP6[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu9_IP5 = Muon_fired_HLT_Mu9_IP5[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu9_IP4 = Muon_fired_HLT_Mu9_IP4[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu10p5_IP3p5 = Muon_fired_HLT_Mu10p5_IP3p5[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu12_IP6 = Muon_fired_HLT_Mu12_IP6[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu8_v1 = Muon_fired_HLT_Mu8_v1[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu8_v12 = Muon_fired_HLT_Mu8_v12[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_Mu7p5_Track7_Jpsi_v11 = Muon_fired_HLT_Mu7p5_Track7_Jpsi_v11[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_L2Mu23NoVtx_2Cha_v1 = Muon_fired_HLT_L2Mu23NoVtx_2Cha_v1[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1 = Muon_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_DST_DoubleMu1_noVtx_CaloScouting_v2 = Muon_fired_DST_DoubleMu1_noVtx_CaloScouting_v2[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_DST_DoubleMu3_noVtx_CaloScouting_v6 = Muon_fired_DST_DoubleMu3_noVtx_CaloScouting_v6[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3 = Muon_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3[JPsiToMuMu_lep1_idx[0]];
    the_tag_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13 = Muon_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu7_IP4 = Muon_prescale_HLT_Mu7_IP4[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu8_IP6 = Muon_prescale_HLT_Mu8_IP6[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu8_IP5 = Muon_prescale_HLT_Mu8_IP5[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu8_IP3 = Muon_prescale_HLT_Mu8_IP3[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu8p5_IP3p5 = Muon_prescale_HLT_Mu8p5_IP3p5[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu9_IP6 = Muon_prescale_HLT_Mu9_IP6[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu9_IP5 = Muon_prescale_HLT_Mu9_IP5[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu9_IP4 = Muon_prescale_HLT_Mu9_IP4[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu10p5_IP3p5 = Muon_prescale_HLT_Mu10p5_IP3p5[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu12_IP6 = Muon_prescale_HLT_Mu12_IP6[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu8_v1 = Muon_prescale_HLT_Mu8_v1[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu8_v12 = Muon_prescale_HLT_Mu8_v12[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_Mu7p5_Track7_Jpsi_v11 = Muon_prescale_HLT_Mu7p5_Track7_Jpsi_v11[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_L2Mu23NoVtx_2Cha_v1 = Muon_prescale_HLT_L2Mu23NoVtx_2Cha_v1[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1 = Muon_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2 = Muon_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6 = Muon_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3 = Muon_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3[JPsiToMuMu_lep1_idx[0]];
    the_tag_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13 = Muon_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13[JPsiToMuMu_lep1_idx[0]];

    the_probe_pt = JPsiToMuMu_lep2_pt[0];
    the_probe_eta = fabs(JPsiToMuMu_lep2_eta[0]);
    the_probe_phi = JPsiToMuMu_lep2_phi[0];
    the_probe_dxy = fabs(Muon_dxy[JPsiToMuMu_lep2_idx[0]]);
    the_probe_dxy_bs = fabs(Muon_dxy_BS[JPsiToMuMu_lep2_idx[0]]);
    the_probe_dxy_bs_uncorrected = fabs(Muon_dxy_BS[JPsiToMuMu_lep2_idx[0]]);
    the_weight_dxy_bs = isMC ? getMCCorrection("mc_weight_probe_dxy_bs.root", the_probe_dxy_bs_uncorrected, 0.3) : 1.;

    //float smeared_corr_0p01 = gauss_function_0p01->GetRandom(); 
    //float smeared_corr_0p02 = gauss_function_0p02->GetRandom(); 
    float smeared_corr_0p03 = gauss_function_0p03->GetRandom(); 
    //float smeared_corr_0p04 = gauss_function_0p04->GetRandom(); 
    //float smeared_corr_0p05 = gauss_function_0p05->GetRandom(); 
    //float smeared_corr_0p06 = gauss_function_0p06->GetRandom(); 
    //float smeared_corr_0p07 = gauss_function_0p07->GetRandom(); 
    //float smeared_corr_0p08 = gauss_function_0p08->GetRandom(); 
    //float smeared_corr_0p09 = gauss_function_0p09->GetRandom(); 
    //float smeared_corr_0p1 = gauss_function_0p1->GetRandom(); 
    //float smeared_corr_0p15 = gauss_function_0p15->GetRandom(); 
    //float smeared_corr_0p2 = gauss_function_0p2->GetRandom(); 
    //float smeared_corr_0p25 = gauss_function_0p25->GetRandom(); 
    //float smeared_corr_0p3 = gauss_function_0p3->GetRandom(); 
    //float smeared_corr_0p35 = gauss_function_0p35->GetRandom(); 
    //float smeared_corr_0p4 = gauss_function_0p4->GetRandom(); 
    //float smeared_corr_0p45 = gauss_function_0p45->GetRandom(); 
    //float smeared_corr_0p5 = gauss_function_0p5->GetRandom(); 

    //the_smeared_corr_0p01 = smeared_corr_0p01;
    //the_smeared_corr_0p02 = smeared_corr_0p02;
    the_smeared_corr_0p03 = smeared_corr_0p03;
    //the_smeared_corr_0p04 = smeared_corr_0p04;
    //the_smeared_corr_0p05 = smeared_corr_0p05;
    //the_smeared_corr_0p06 = smeared_corr_0p06;
    //the_smeared_corr_0p07 = smeared_corr_0p07;
    //the_smeared_corr_0p08 = smeared_corr_0p08;
    //the_smeared_corr_0p09 = smeared_corr_0p09;
    //the_smeared_corr_0p1 = smeared_corr_0p1;
    //the_smeared_corr_0p15 = smeared_corr_0p15;
    //the_smeared_corr_0p2 = smeared_corr_0p2;
    //the_smeared_corr_0p25 = smeared_corr_0p25;
    //the_smeared_corr_0p3 = smeared_corr_0p3;
    //the_smeared_corr_0p35 = smeared_corr_0p35;
    //the_smeared_corr_0p4 = smeared_corr_0p4;
    //the_smeared_corr_0p45 = smeared_corr_0p45;
    //the_smeared_corr_0p5 = smeared_corr_0p5;

    //the_double_gauss_corr_1 = double_gauss_1->GetRandom();
    //the_double_gauss_corr_2 = double_gauss_2->GetRandom();

    the_probe_dz = fabs(Muon_dz[JPsiToMuMu_lep2_idx[0]]);
    the_probe_dxy_sig = fabs(Muon_dxyS[JPsiToMuMu_lep2_idx[0]]);

    the_probe_dxy_sig_bs_uncorrected = fabs(Muon_dxyS_BS[JPsiToMuMu_lep2_idx[0]]);
    //the_probe_dxy_sig_bs = fabs(Muon_dxyS_BS[JPsiToMuMu_lep2_idx[0]]);
    // linear correction applied
    //the_probe_dxy_sig_bs = isMC ? (1.178 -0.002968 * the_probe_dxy_sig_bs_uncorrected) * the_probe_dxy_sig_bs_uncorrected : the_probe_dxy_sig_bs_uncorrected;

    // scale plus smear correction applied
    the_probe_dxy_sig_bs = isMC ? the_probe_dxy_sig_bs_uncorrected * 1.12 * the_smeared_corr_0p03 : the_probe_dxy_sig_bs_uncorrected;

    the_weight_dxy_sig_bs = isMC ? getMCCorrection("mc_weight_probe_dxy_sig_bs.root", the_probe_dxy_sig_bs_uncorrected, 60) : 1.;
    the_probe_dxy_sig_bs_corrected_linearscale = isMC ? (1.178 -0.002968 * the_probe_dxy_sig_bs_uncorrected) * the_probe_dxy_sig_bs_uncorrected : the_probe_dxy_sig_bs_uncorrected;

    //std::cout << the_probe_dxy_sig_bs_uncorrected << " " << getMCCorrection("mc_weight_probe_dxy_sig_bs.root", the_probe_dxy_sig_bs_uncorrected, 60) << std::endl;

    //the_probe_dxy_sig_bs = corr_dxy_sig_bs * fabs(Muon_dxyS_BS[JPsiToMuMu_lep2_idx[0]]);

    the_probe_dxy_sig_bs_rdst = fabs(Muon_dxyS_BS_alaRdst[JPsiToMuMu_lep2_idx[0]]);
    the_probe_dz_sig = fabs(Muon_dzS[JPsiToMuMu_lep2_idx[0]]);
    the_probe_isloose = Muon_looseId[JPsiToMuMu_lep2_idx[0]];
    the_probe_issoft = Muon_softId[JPsiToMuMu_lep2_idx[0]];
    the_probe_ismedium = Muon_mediumId[JPsiToMuMu_lep2_idx[0]];
    the_probe_istight = Muon_tightId[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu7_IP4 = Muon_fired_HLT_Mu7_IP4[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu8_IP6 = Muon_fired_HLT_Mu8_IP6[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu8_IP5 = Muon_fired_HLT_Mu8_IP5[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu8_IP3 = Muon_fired_HLT_Mu8_IP3[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu8p5_IP3p5 = Muon_fired_HLT_Mu8p5_IP3p5[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu9_IP6 = Muon_fired_HLT_Mu9_IP6[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu9_IP5 = Muon_fired_HLT_Mu9_IP5[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu9_IP4 = Muon_fired_HLT_Mu9_IP4[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu10p5_IP3p5 = Muon_fired_HLT_Mu10p5_IP3p5[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu12_IP6 = Muon_fired_HLT_Mu12_IP6[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu8_v1 = Muon_fired_HLT_Mu8_v1[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu8_v12 = Muon_fired_HLT_Mu8_v12[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_Mu7p5_Track7_Jpsi_v11 = Muon_fired_HLT_Mu7p5_Track7_Jpsi_v11[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_L2Mu23NoVtx_2Cha_v1 = Muon_fired_HLT_L2Mu23NoVtx_2Cha_v1[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1 = Muon_fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_DST_DoubleMu1_noVtx_CaloScouting_v2 = Muon_fired_DST_DoubleMu1_noVtx_CaloScouting_v2[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_DST_DoubleMu3_noVtx_CaloScouting_v6 = Muon_fired_DST_DoubleMu3_noVtx_CaloScouting_v6[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3 = Muon_fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3[JPsiToMuMu_lep2_idx[0]];
    the_probe_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13 = Muon_fired_HLT_BTagMu_AK4DiJet40_Mu5_v13[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu7_IP4 = Muon_prescale_HLT_Mu7_IP4[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu8_IP6 = Muon_prescale_HLT_Mu8_IP6[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu8_IP5 = Muon_prescale_HLT_Mu8_IP5[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu8_IP3 = Muon_prescale_HLT_Mu8_IP3[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu8p5_IP3p5 = Muon_prescale_HLT_Mu8p5_IP3p5[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu9_IP6 = Muon_prescale_HLT_Mu9_IP6[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu9_IP5 = Muon_prescale_HLT_Mu9_IP5[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu9_IP4 = Muon_prescale_HLT_Mu9_IP4[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu10p5_IP3p5 = Muon_prescale_HLT_Mu10p5_IP3p5[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu12_IP6 = Muon_prescale_HLT_Mu12_IP6[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu8_v1 = Muon_prescale_HLT_Mu8_v1[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu8_v12 = Muon_prescale_HLT_Mu8_v12[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_Mu7p5_Track7_Jpsi_v11 = Muon_prescale_HLT_Mu7p5_Track7_Jpsi_v11[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_L2Mu23NoVtx_2Cha_v1 = Muon_prescale_HLT_L2Mu23NoVtx_2Cha_v1[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1 = Muon_prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2 = Muon_prescale_DST_DoubleMu1_noVtx_CaloScouting_v2[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6 = Muon_prescale_DST_DoubleMu3_noVtx_CaloScouting_v6[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3 = Muon_prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3[JPsiToMuMu_lep2_idx[0]];
    the_probe_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13 = Muon_prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13[JPsiToMuMu_lep2_idx[0]];

    if(the_probe_fired_HLT_Mu7_IP4==1 || the_probe_fired_HLT_Mu8_IP6==1 || the_probe_fired_HLT_Mu8_IP5==1 || the_probe_fired_HLT_Mu8_IP3==1 || the_probe_fired_HLT_Mu8p5_IP3p5==1 || the_probe_fired_HLT_Mu9_IP6==1 || the_probe_fired_HLT_Mu9_IP5==1 || the_probe_fired_HLT_Mu9_IP4 ==1 || the_probe_fired_HLT_Mu10p5_IP3p5==1 || the_probe_fired_HLT_Mu12_IP6==1){
      the_probe_fired_BParkingHLT = 1;
    }
    else{
      the_probe_fired_BParkingHLT = 0;
    }

    the_pv_npvs = *PV_npvs;

    // weights
    the_weight_hlt_A1 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_DST_DoubleMu1_A1_v1/scaleFactor_results_cat_pt_eta_fit.root", the_tag_pt, fabs(the_tag_eta)) : 1.;
    the_weight_hlt_A1_6 = isMC ? getTriggerScaleFactor("/t3home/anlyon/BHNL/BHNLNano/CMSSW_10_2_15/src/PhysicsTools/TagAndProbe/test/results/tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_DST_DoubleMu1_A1_6_v1/scaleFactor_results_cat_pt_eta_fit.root", the_tag_pt, fabs(the_tag_eta)) : 1.;
    the_weight_pu_A = isMC ? getPUWeight("pileup_weight_dataA_mcAutumn18.root", *Pileup_nTrueInt) : 1.;
    the_weight_pu_B = isMC ? getPUWeight("pileup_weight_dataB_mcAutumn18.root", *Pileup_nTrueInt) : 1.;
    the_weight_pu_C = isMC ? getPUWeight("pileup_weight_dataC_mcAutumn18.root", *Pileup_nTrueInt) : 1.;
    the_weight_pu_D = isMC ? getPUWeight("pileup_weight_dataD_mcAutumn18.root", *Pileup_nTrueInt) : 1.;
    the_weight_pu_tot = isMC ? getPUWeight("pileup_weight_datatot_mcAutumn18.root", *Pileup_nTrueInt) : 1.;

    tree->Fill();
  }

  return kTRUE;
}

void TagAndProbeDumper::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void TagAndProbeDumper::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  my_file->cd();
  tree->Write("", TObject::kOverwrite);
  my_file->Close();

  cout << "- End Tag and Probe Dumper -" << endl;
}
