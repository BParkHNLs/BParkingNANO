#ifndef utils
#define utils

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

// ------------------------------------------ //
//  extra functions needed by the ntupliser
// ------------------------------------------ //


bool sortcansbydesc(const pair<int, float> &a1, const pair<int, float> &a2){
  return a1.second > a2.second;
}


bool sortcansbydesc_opp(const pair<int, float> &a1, const pair<int, float> &a2){
  return a1.second < a2.second;
}


bool sortcansbydesc_diff(const pair<int, float> &a1, const pair<int, float> &a2){
  return fabs(a1.second-6.274) < fabs(a2.second-6.274);
}


vector<pair<int,float>> createPairWithDesc(const UInt_t& nCand, const TTreeReaderArray<Float_t>& desc){
  vector<pair<int,float>> pair_candIdx_desc;
  for(unsigned int iCand(0); iCand < nCand; ++iCand){
    pair<int, float> pair_candIdx_desc_tmp;
    pair_candIdx_desc_tmp.first  = iCand;
    pair_candIdx_desc_tmp.second = desc[iCand] ;
    pair_candIdx_desc.push_back(pair_candIdx_desc_tmp);
  }
  return pair_candIdx_desc;
}


vector<pair<int,float>> updatePairWithDesc(vector<pair<int,float>> the_ini_pair, const TTreeReaderArray<Int_t>& quantity){
  vector<pair<int,float>> pair_candIdx_desc;
  for(unsigned int iCand(0); iCand < the_ini_pair.size(); ++iCand){
    pair<int, float> pair_candIdx_desc_tmp;
    pair_candIdx_desc_tmp.first = the_ini_pair[iCand].first;
    pair_candIdx_desc_tmp.second = fabs(quantity[the_ini_pair[iCand].first]); // taking the abs of the quantity
    pair_candIdx_desc.push_back(pair_candIdx_desc_tmp);
  }
  return pair_candIdx_desc;
}


vector<pair<int,float>> updatePairWithDesc(vector<pair<int,float>> the_ini_pair, const TTreeReaderArray<Float_t>& quantity){
  vector<pair<int,float>> pair_candIdx_desc;
  for(unsigned int iCand(0); iCand < the_ini_pair.size(); ++iCand){
    pair<int, float> pair_candIdx_desc_tmp;
    pair_candIdx_desc_tmp.first = the_ini_pair[iCand].first;
    pair_candIdx_desc_tmp.second = fabs(quantity[the_ini_pair[iCand].first]); // taking the abs of the quantity
    pair_candIdx_desc.push_back(pair_candIdx_desc_tmp);
  }
  return pair_candIdx_desc;
}


vector<pair<int,float>> updatePairWithDesc(vector<pair<int,float>> the_ini_pair, const TTreeReaderArray<Int_t>& mu_idx, const TTreeReaderArray<Int_t>& quantity){
  vector<pair<int,float>> pair_candIdx_desc;
  for(unsigned int iCand(0); iCand < the_ini_pair.size(); ++iCand){
    pair<int, float> pair_candIdx_desc_tmp;
    pair_candIdx_desc_tmp.first = the_ini_pair[iCand].first;
    pair_candIdx_desc_tmp.second = fabs(quantity[mu_idx[the_ini_pair[iCand].first]]); // taking the abs of the quantity
    pair_candIdx_desc.push_back(pair_candIdx_desc_tmp);
  }
  return pair_candIdx_desc;
}


Float_t get3Ddisp(const Float_t vx1, const Float_t vx2, const Float_t vy1, const Float_t vy2, const Float_t vz1, const Float_t vz2){
  return TMath::Sqrt( (vx1-vx2)*(vx1-vx2) + (vy1-vy2)*(vy1-vy2) + (vz1-vz2)*(vz1-vz2) );
}


Float_t get2Ddisp(const Float_t vx1, const Float_t vx2, const Float_t vy1, const Float_t vy2){
  return TMath::Sqrt( (vx1-vx2)*(vx1-vx2) + (vy1-vy2)*(vy1-vy2) );
}


float deltaR(float eta1, float eta2, float phi1, float phi2){
  float pi = 3.14159265359;
  float delta_eta = eta1 - eta2;
  float delta_phi = phi1 - phi2;
  if(fabs(delta_phi) > pi){
    delta_phi -= 2*pi;
  }
  return sqrt(pow(delta_eta, 2) + pow(delta_phi, 2));
}


bool checkLumi(string lumi_ranges, int lumi, int seed = 1){
  // get boundaries of lumi range
  string range_min = lumi_ranges.substr(lumi_ranges.find("[", seed)+1, int(lumi_ranges.find(",", seed)) - int(lumi_ranges.find("[", seed)+1));
  string range_max = lumi_ranges.substr(lumi_ranges.find(",", seed)+2, int(lumi_ranges.find("]", seed)) - int(lumi_ranges.find(",", seed)+2));
  //cout << range_min << " " << range_max << endl;

  if(lumi >= stoi(range_min) && lumi <= stoi(range_max)){
    //cout << "lumi accepted" << endl;
    return true;
  }
  else if(lumi_ranges.substr(int(lumi_ranges.find("]", seed))+1, 1) != "]"){ // there are further lumi blocks
    return checkLumi(lumi_ranges, lumi, int(lumi_ranges.find("]", seed)+2));
  }
  else{
    return false;
  }
}


bool lumiMask(int run, int lumi){
  //cout << "run " << run << " lumi " << lumi << endl;
  
  // put the content of the json in a tree
  boost::property_tree::ptree the_json_tree;
  boost::property_tree::read_json("golden_2018.json", the_json_tree);

  // checking that the run exists
  std::string lumi_ranges = the_json_tree.get<std::string> (std::to_string(run), "RUN NOT FOUND");

  if(lumi_ranges == "RUN NOT FOUND"){
    return false; 
  }
  else{
    // check that lumi is valid
    return checkLumi(lumi_ranges, lumi);
  }
}

float getMCCorrection(TString filename, double var, double max_val){
  // get file
  TFile* file = TFile::Open(filename);
  file->cd();

  // get histogram
  TH1D* hist = (TH1D*) file->Get("hist_ratio")->Clone("hist");

  // truncate the variable
  var = std::max(0., std::min(max_val-1e-03, double(var)));
  
  // get the bin
  int bin = hist->GetXaxis()->FindBin(var);
  
  // get weight
  Float_t mc_weight = hist->GetBinContent(bin);

  file->Close();

  return mc_weight;
}


float getPUWeight(TString filename, int var){
  // get file
  TFile* file = TFile::Open(filename);
  file->cd();

  // get histogram
  TH1D* hist = (TH1D*) file->Get("hist_weight")->Clone("hist");

  var = std::max(0, std::min(200, int(var)));

  // get weight
  Float_t pu_weight = hist->GetBinContent(var);
  
  file->Close();

  return pu_weight;
}


float getLeptonScaleFactor(TString filename, string ID, float pt, float eta, TString flag=""){
  // get file
  TFile* file = TFile::Open(filename);
  file->cd();

  // get histogram
  TH1D* hist = 0;
  if(ID == "softid"){
    hist = (TH1D*) file->Get("NUM_SoftID_DEN_genTracks_pt_abseta")->Clone("hist");
  }
  else if(ID == "looseid"){
    hist = (TH1D*) file->Get("NUM_LooseID_DEN_genTracks_pt_abseta")->Clone("hist");
  }

  // get bin
  pt = std::max(0., std::min(40., double(pt)));
  int bin_pt = hist->GetXaxis()->FindBin(pt);
  int bin_eta = hist->GetYaxis()->FindBin(eta);

  // get scale factor
  Float_t scale_factor;
  if(flag == "plus_one_sigma"){
    scale_factor = hist->GetBinContent(bin_pt, bin_eta) + hist->GetBinError(bin_pt, bin_eta);
  }
  else if(flag == "minus_one_sigma"){
    scale_factor = hist->GetBinContent(bin_pt, bin_eta) - hist->GetBinError(bin_pt, bin_eta);
  }
  else{
    scale_factor = hist->GetBinContent(bin_pt, bin_eta);
  }

  file->Close();

  return scale_factor;
}


float getTriggerScaleFactor(TString filename_sf, float pt, float dxysig){
  // get trigger scale factor file
  TFile* file_sf = TFile::Open(filename_sf);
  file_sf->cd();

  // get histogram
  TH2D* hist_sf = (TH2D*) file_sf->Get("hist_scale_factor")->Clone("hist_sf");

  pt = std::max(6. std::min(99.9, double(pt)));
  dxysig = std::max(0., std::min(499.9, double(dxysig)));

  // get bin
  int bin_pt = hist_sf->GetXaxis()->FindBin(pt);
  int bin_dxysig = hist_sf->GetYaxis()->FindBin(dxysig);

  // get scale factor
  Float_t scale_factor = hist_sf->GetBinContent(bin_pt, bin_dxysig);
  
  file_sf->Close();

  return scale_factor;
}


float getTriggerScaleFactor_v2(TString filename_sf, float pt, float dxysig){
  // get trigger scale factor file
  TFile* file_sf = TFile::Open(filename_sf);
  file_sf->cd();

  // get histogram
  TH2D* hist_sf = (TH2D*) file_sf->Get("hist_scale_factor")->Clone("hist_sf");

  //std::cout << "pt ini: " << pt << " dxysig ini: " << dxysig << std::endl;

  pt = std::max(0.1, std::min(99.9, double(pt)));
  dxysig = std::max(0., std::min(499.9, double(dxysig)));

  //std::cout << "pt fin: " << pt << " dxysig fin: " << dxysig << std::endl;

  // get bin
  int bin_pt = hist_sf->GetXaxis()->FindBin(pt);
  int bin_dxysig = hist_sf->GetYaxis()->FindBin(dxysig);

  //std::cout << bin_pt << " " << bin_dxysig << std::endl;

  // get scale factor
  Float_t scale_factor = hist_sf->GetBinContent(bin_pt, bin_dxysig);
  
  file_sf->Close();

  return scale_factor;
}


#endif
