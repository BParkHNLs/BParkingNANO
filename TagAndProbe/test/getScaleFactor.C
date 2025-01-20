// usage: root -l -b getScaleFactor()

#include<iostream>
#include<fstream>
#include "TFile.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"

// ------- Global Variables ------- //

//TString dataFileName = "results_tag_and_probe_v2_tag_fired_DST_DoubleMu1_data_A1_extraptbin.root";
//TString dataFileName = "results_tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_HLT_Mu9_IP6_A1_6.root";
//TString mcFileName = "results_tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_HLT_Mu9_IP6_mc.root";
//string dirLabel = "tag_and_probe_v2_BToJPsiKstar_V0_tag_fired_HLT_Mu9_IP6_A1_6";


// -------------------------------- //

void  write2DScaleFactor(TH2F* h, string name){
  string name_fortable = name + ".txt";
  ofstream outFile(name_fortable);

  TString name_forhist = name.c_str();
  TFile* root_file = TFile::Open(name_forhist + ".root", "RECREATE");
  TH2D* hist_scale_factor = (TH2D*)h->Clone("hist_scale_factor");
  TCanvas* c = new TCanvas("c", "c", 1200, 1000);

  int nX = h->GetNbinsX();
  int nY = h->GetNbinsY();

  for(int i=1; i<=nX; ++i) {
    Double_t pT0 = h->GetXaxis()->GetBinLowEdge(i);
    Double_t pT1 = h->GetXaxis()->GetBinLowEdge(i+1);

    for(int j=1; j<=nY; ++j) {
      Double_t x = h->GetBinContent(i,j);
      Double_t dx = h->GetBinError(i,j);
      Double_t err = (dx/x)*100;
      Double_t eta0 = h->GetYaxis()->GetBinLowEdge(j);
      Double_t eta1 = h->GetYaxis()->GetBinLowEdge(j+1);

      outFile << pT0 << " " << pT1 << " " << eta0 << " " << eta1 << " " << x << " " << dx << " " << err << "%"<< std::endl;

      hist_scale_factor->SetBinContent(i, j, x);
      hist_scale_factor->SetBinError(i, j, dx);
    }
  }

  hist_scale_factor->GetZaxis()->SetTitle("Scale Factor");
  hist_scale_factor->SetOption("colztexte");
  hist_scale_factor->Write();
  hist_scale_factor->Draw();
  gStyle->SetPaintTextFormat(".1f");
  c->SaveAs(name_forhist + ".pdf");
  root_file->Close();
}


void  write1DScaleFactor(TH1F* h, string name){
  string name_fortable = name + ".txt";
  ofstream outFile(name_fortable);

  TString name_forhist = name.c_str();
  TFile* root_file = TFile::Open(name_forhist + ".root", "RECREATE");
  TH1D* hist_scale_factor = (TH1D*)h->Clone("hist_scale_factor");
  TCanvas* c = new TCanvas("c", "c", 800, 700);

  int nX = h->GetNbinsX();

  for(int i=1; i<=nX; ++i) {
    Double_t pT0 = h->GetXaxis()->GetBinLowEdge(i);
    Double_t pT1 = h->GetXaxis()->GetBinLowEdge(i+1);

    Double_t x = h->GetBinContent(i);
    Double_t dx = h->GetBinError(i);
    Double_t err = (dx/x)*100;

    outFile << pT0 << " " << pT1 << " " << x << " " << dx << " " << err << "%"<< std::endl;

    hist_scale_factor->SetBinContent(i, x);
    hist_scale_factor->SetBinError(i, dx);
  }

  hist_scale_factor->GetYaxis()->SetTitle("Scale Factor");
  hist_scale_factor->Write();
  hist_scale_factor->Draw();
  c->SaveAs(name_forhist + ".png");
  root_file->Close();
}

void  produce2DPlot(TH2F* h, string name){
  TString name_forhist = name.c_str();
  TFile* root_file = TFile::Open(name_forhist + ".root", "RECREATE");
  TH2D* hist_scale_factor = new TH2D("hist_scale_factor", "hist_scale_factor", h->GetNbinsX(), 0, h->GetNbinsX(), h->GetNbinsY(), 0, h->GetNbinsY()-1);
  TCanvas* c = new TCanvas("c", "c", 1200, 1000);
  c->SetRightMargin(0.15);

  int nX = h->GetNbinsX();
  int nY = h->GetNbinsY();

  for(int i=1; i<=nX; ++i) {
    Double_t pT0 = h->GetXaxis()->GetBinLowEdge(i);
    Double_t pT1 = h->GetXaxis()->GetBinLowEdge(i+1);

    for(int j=1; j<=nY; ++j) {
      Double_t x = h->GetBinContent(i,j);
      Double_t dx = h->GetBinError(i,j);
      Double_t err = (dx/x)*100;
      Double_t eta0 = h->GetYaxis()->GetBinLowEdge(j);
      Double_t eta1 = h->GetYaxis()->GetBinLowEdge(j+1);

      hist_scale_factor->SetBinContent(i, j, x);
      hist_scale_factor->SetBinError(i, j, dx);
    }
  }

  hist_scale_factor->GetZaxis()->SetTitle("Scale Factor");
  hist_scale_factor->GetZaxis()->SetRangeUser(-1e-3, 1);
  hist_scale_factor->SetOption("colztexte");
  hist_scale_factor->SetTitle("");
  hist_scale_factor->Write();
  hist_scale_factor->Draw();
  gStyle->SetTitleFillColor(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".3f");
  c->SaveAs(name_forhist + ".pdf");
  c->SaveAs(name_forhist + ".png");
  root_file->Close();
}


void process(string dir, string subdir, string outdir, string method, TFile* fData, TFile* fMC, TString canvasname, TString extraname, string dim){
  TCanvas* c;

  TString results_dir = dir + "/" + subdir + "/" + method + "_eff_plots/";
  std::cout << results_dir << std::endl;
  fData->cd(results_dir);

  // data
  c = (TCanvas*) gDirectory->Get(canvasname);
  TH2F* hData2D = 0;
  TH1F* hData1D = 0;
  if(dim == "2D" || dim == "3D") hData2D = (TH2F*) c->FindObject(canvasname);
  else if(dim == "1D") hData1D = (TH1F*) c->FindObject(canvasname);
  hData2D->SetDefaultSumw2();

  string name_data = outdir + "efficiency_data_" + dir + "_" + subdir + "_" + method;
  if(dim == "3D") name_data += "_" + extraname;
  if(dim == "2D" || dim == "3D"){
    write2DScaleFactor(hData2D, name_data);
    produce2DPlot(hData2D, name_data);
  }
  else if(dim == "1D"){
    write1DScaleFactor(hData1D, name_data);
  }

  std::cout << name_data << ".txt created " << std::endl;
  std::cout << name_data << ".root created " << std::endl;
  std::cout << name_data << ".png created " << std::endl;

  // mc
  fMC->cd(results_dir);
  c = (TCanvas*) gDirectory->Get(canvasname);
  TH2F* hMC2D = 0;
  TH1F* hMC1D = 0;
  if(dim == "2D" || dim == "3D") hMC2D = (TH2F*) c->FindObject(canvasname);
  else if(dim == "1D") hMC1D = (TH1F*) c->FindObject(canvasname);
  hMC2D->SetDefaultSumw2();

  string name_mc = outdir + "efficiency_mc_" + dir + "_" + subdir + "_" + method;
  if(dim == "3D") name_mc += "_" + extraname;
  if(dim == "2D" || dim == "3D"){
    write2DScaleFactor(hMC2D, name_mc);
    produce2DPlot(hMC2D, name_mc);
  }
  else if(dim == "1D"){
    write1DScaleFactor(hMC1D, name_mc);
  }

  std::cout << name_mc << ".txt created " << std::endl;
  std::cout << name_mc << ".root created " << std::endl;
  std::cout << name_mc << ".png created " << std::endl;

  //hData2D->SetDefaultSumw2();
  //hMC2D->SetDefaultSumw2();

  // ratio
  string name = outdir + "scaleFactor_" + dir + "_" + subdir + "_" + method;
  if(dim == "3D") name += "_" + extraname;
  if(dim == "2D" || dim == "3D") hData2D->Divide(hMC2D);
  else if(dim == "1D") hData1D->Divide(hMC1D);
  if(dim == "2D" || dim == "3D"){
    write2DScaleFactor(hData2D, name);
    produce2DPlot(hData2D, name);
  }
  else if(dim == "1D"){
    write1DScaleFactor(hData1D, name);
  }

  std::cout << name << ".txt created " << std::endl;
  std::cout << name << ".root created " << std::endl;
  std::cout << name << ".png created " << std::endl;

}


void getScaleFactor(TString dataFileName, TString mcFileName, string dirLabel, string categorisation="pt_eta", string dir="results"){

  //TString dataFileName = dataFileName.c_str();
  //TString mcFileName = mcFileName.c_str();

  string outdir = "./results/" + dirLabel + "/";
  system(Form("mkdir -p %s", outdir.c_str()));

  TFile* fData = new TFile(dataFileName);
  TFile* fMC = new TFile(mcFileName);

  if(categorisation == "pt_eta"){
    process(dir, "cat_eff", outdir, "fit", fData, fMC, "probe_pt_probe_eta_PLOT", "", "2D"); 
    process(dir, "cat_eff", outdir, "fit", fData, fMC, "probe_pt_probe_eta_PLOT", "", "2D"); 
    //process(dir, "cat_pt_eta", outdir, "cnt", fData, fMC, "probe_pt_probe_eta_PLOT", "", "2D"); 
    //process(dir, "cat_pt_eta", outdir, "cnt", fData, fMC, "probe_pt_probe_eta_PLOT", "", "2D"); 
  }
  else if(categorisation == "pt_dxysig"){
    process(dir, "cat_eff", outdir, "fit", fData, fMC, "probe_pt_probe_dxy_sig_PLOT", "", "2D"); 
    process(dir, "cat_eff", outdir, "cnt", fData, fMC, "probe_pt_probe_dxy_sig_PLOT", "", "2D"); 
    //process(dir, "cat_pt_dxysig", outdir, "fit", fData, fMC, "probe_pt_probe_dxy_sig_PLOT", "", "2D"); 
    //process(dir, "cat_pt_dxysig", outdir, "cnt", fData, fMC, "probe_pt_probe_dxy_sig_PLOT", "", "2D"); 
  }
  else if(categorisation == "pt_eta_dxysig"){
    // make sure that the eta categories are correct
    process(dir, "cat_eff", outdir, "fit", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin0", "eta_0p00_0p50", "3D"); 
    process(dir, "cat_eff", outdir, "cnt", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin0", "eta_0p00_0p50", "3D"); 
    process(dir, "cat_eff", outdir, "fit", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin1", "eta_0p50_1p00", "3D"); 
    process(dir, "cat_eff", outdir, "cnt", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin1", "eta_0p50_1p00", "3D"); 
    process(dir, "cat_eff", outdir, "fit", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin2", "eta_1p00_1p50", "3D"); 
    process(dir, "cat_eff", outdir, "cnt", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin2", "eta_1p00_1p50", "3D"); 
    process(dir, "cat_eff", outdir, "fit", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin3", "eta_1p50_2p00", "3D"); 
    process(dir, "cat_eff", outdir, "cnt", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin3", "eta_1p50_2p00", "3D"); 
    //process(dir, "cat_pt_eta_dxysig", outdir, "fit", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin0", "eta_0p00_0p50", "3D"); 
    //process(dir, "cat_pt_eta_dxysig", outdir, "cnt", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin0", "eta_0p00_0p50", "3D"); 
    //process(dir, "cat_pt_eta_dxysig", outdir, "fit", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin1", "eta_0p50_1p00", "3D"); 
    //process(dir, "cat_pt_eta_dxysig", outdir, "cnt", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin1", "eta_0p50_1p00", "3D"); 
    //process(dir, "cat_pt_eta_dxysig", outdir, "fit", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin2", "eta_1p00_1p50", "3D"); 
    //process(dir, "cat_pt_eta_dxysig", outdir, "cnt", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin2", "eta_1p00_1p50", "3D"); 
    //process(dir, "cat_pt_eta_dxysig", outdir, "fit", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin2", "eta_1p50_2p00", "3D"); 
    //process(dir, "cat_pt_eta_dxysig", outdir, "cnt", fData, fMC, "probe_pt_probe_dxy_sig_PLOT_probe_eta_bin2", "eta_1p50_2p00", "3D"); 
  }

  fData->Close();
  fMC->Close();
  delete fData;
  delete fMC;

  std::cout << "Done " << std::endl;
}



