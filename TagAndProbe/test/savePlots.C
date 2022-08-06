// usage: root -l -b savePlots.C+
#include<iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TStyle.h"
#include "TH2F.h"


// ------- Global Variables ------- //

//TString inputFileName = "results_tag_and_probe_v2_tag_fired_DST_DoubleMu1_data_A1_6_forplots.root";
//TString inputFileName = "/work/anlyon/tag_and_probe/outfiles/test_mc_tag_fired_anyBParkHLT_ptetadxysig/results_test_mc_tag_fired_anyBParkHLT_ptetadxysig_incl.root";
//string dirLabel = "test_mc_tag_fired_anyBParkHLT_ptetadxysig";
//Bool_t isMC = true;


// -------------------------------- //


void savePlots(
         TString inputFileName,
         string dirLabel,
         string isMC,
			   string dir="results", 
			   bool isCutAndCount=false, 
			   bool isMCTruth = false ) {

  gStyle->SetPaintTextFormat(".1f");

  string outdir = "./results/" + dirLabel;
  if(isMC == "True"){
    outdir += "/mc/";
  }
  else{
    outdir += "/data/";
  }
  system(Form("mkdir -p %s", outdir.c_str()));

  TFile* f = new TFile(inputFileName);
  f->cd(dir.c_str());

  TKey* key;
  TCanvas* c;
  TIter next(gDirectory->GetListOfKeys()); 
  while ((key = (TKey*) next())) {
    TObject *obj = key->ReadObj();
    TString name = obj->GetName();

    if ( !(obj->IsA()->InheritsFrom( "TDirectory" )) ) continue;

    TString dirName = name + "/cnt_eff_plots/"; 
    if( !isCutAndCount ) dirName = name + "/fit_eff_plots/";
    gDirectory->cd(dirName);

      TKey *innerkey;
      TIter innernext(gDirectory->GetListOfKeys());
      while ((innerkey = (TKey*) innernext())){
        obj = innerkey->ReadObj();
        TString canvasname = obj->GetName();
        c = (TCanvas*) gDirectory->Get(canvasname);
        c->Draw();
        c->SaveAs(TString(outdir)+TString(canvasname)+TString(".png")); 
        c->SaveAs(TString(outdir)+TString(canvasname)+TString(".pdf")); 
        c->Draw();
      }

    if(!isCutAndCount) {
      //now make plot of the fit distributions
      gDirectory->cd("../");
      gDirectory->ls();

      TKey *innerkey;
      TIter innernext(gDirectory->GetListOfKeys());
      while ((innerkey = (TKey*) innernext())){
      obj = innerkey->ReadObj();
      auto innername = obj->GetName();
      if(!(obj->IsA()->InheritsFrom("TDirectory")) || 
      !(TString(innername).Contains("_bin")) ) continue;
      gDirectory->cd(innername);
      c = (TCanvas*) gDirectory->Get("fit_canvas");
      c->Draw();
      TString plotname = TString("fit")+TString(name)+TString("_")+
      TString(innername)+TString(".png");
      TString(innername)+TString(".pdf");
      std::cout << "plotname " << plotname << std::endl;
      c->SaveAs(TString(outdir) + plotname); 
      gDirectory->cd("../");
      } // end while innerkey
    } // end isCutAndCount

    //get back to the initial directory
    if( !isCutAndCount ) gDirectory->cd("../");
    else gDirectory->cd("../../");

    std::cout << "Plots saved in " << outdir << std::endl;
    std::cout << "Done" << std::endl;
  }
} 



