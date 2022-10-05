// class to produce 2 pat::MuonCollections
// one matched to the Park triggers
// another fitered wrt the Park triggers


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "../interface/ETHMuon.h"

#include <TLorentzVector.h>
#include "../interface/helper.h"

using namespace std;

constexpr bool debug = false; //false;

class MuonTriggerSelector : public edm::EDProducer {

  public:

    explicit MuonTriggerSelector(const edm::ParameterSet &iConfig);

    ~MuonTriggerSelector() override {};


  private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<std::vector<reco::Track>> displacedStandaloneMuonSrc_;
    const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

    // added
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
    const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;

    // trigger muon matching
    const double max_deltaR_trigger_matching_;
    const double max_deltaPtRel_trigger_matching_;
    
    // add displaced standalone muons
    const bool add_dsa_;

    // dsa to slimmed muon matching
    const bool do_dsa_matching_;
    const double max_deltaR_dsaToSlimmed_matching_;
    const double max_deltaPtRel_dsaToSlimmed_matching_;

    // for the sel muon
    const double selmu_ptMin_;          // min pT in all muons for B candidates
    const double selmu_absEtaMax_;      // max eta ""
    std::vector<std::string> HLTPaths_;
    //std::vector<std::string> L1Seeds_;
};


MuonTriggerSelector::MuonTriggerSelector(const edm::ParameterSet &iConfig):
  muonSrc_(consumes<std::vector<pat::Muon>> (iConfig.getParameter<edm::InputTag>("muonCollection"))),
  displacedStandaloneMuonSrc_(consumes<std::vector<reco::Track>> (iConfig.getParameter<edm::InputTag>("displacedStandaloneMuonCollection"))),
  vertexSrc_( consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  // added
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  //
  max_deltaR_trigger_matching_(iConfig.getParameter<double>("max_deltaR_trigger_matching")),
  max_deltaPtRel_trigger_matching_(iConfig.getParameter<double>("max_deltaPtRel_trigger_matching")),
  add_dsa_(iConfig.getParameter<bool>("add_dsa")),
  do_dsa_matching_(iConfig.getParameter<bool>("do_dsa_matching")),
  max_deltaR_dsaToSlimmed_matching_(iConfig.getParameter<double>("max_deltaR_dsaToSlimmed_matching")),
  max_deltaPtRel_dsaToSlimmed_matching_(iConfig.getParameter<double>("max_deltaPtRel_dsaToSlimmed_matching")),
  selmu_ptMin_(iConfig.getParameter<double>("selmu_ptMin")),
  selmu_absEtaMax_(iConfig.getParameter<double>("selmu_absEtaMax")),
  HLTPaths_(iConfig.getParameter<std::vector<std::string>>("HLTPaths"))
  //L1Seeds_(iConfig.getParameter<std::vector<std::string>>("L1seeds"))
{
  // produce 2 collections: trgMuons (tags) and SelectedMuons (probes & tags if survive preselection cuts)
  // trigger muons are slimmedMuons only
  // selected muons are slimmedMuons and displacedStandaloneMuons
  // make sure that the indices between the selected muons and transient tracks are consistent
  produces<pat::MuonCollection>("trgMuons"); 
  produces<std::vector<pat::ETHMuon>>("SelectedMuons");
  produces<TransientTrackCollection>("SelectedTransientMuons");  
}


void MuonTriggerSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::ESHandle<MagneticField> bFieldHandle;
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

    std::unique_ptr<pat::MuonCollection>      trgmuons_out   ( new pat::MuonCollection );
    std::unique_ptr<std::vector<pat::ETHMuon>> ETHmuons_out  ( new std::vector<pat::ETHMuon> );
    std::unique_ptr<TransientTrackCollection> trans_muons_out( new TransientTrackCollection );
    
    edm::Handle<std::vector<pat::Muon>> slimmed_muons;
    iEvent.getByToken(muonSrc_, slimmed_muons);

    edm::Handle<std::vector<reco::Track>> displaced_standalone_muons;
    iEvent.getByToken(displacedStandaloneMuonSrc_, displaced_standalone_muons);

    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByToken(vertexSrc_, vertexHandle);
    const reco::VertexCollection& vertices = *vertexHandle;

    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(triggerPrescales_, triggerPrescales);

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    if ( ! beamSpotHandle.isValid() ) {
      edm::LogError("MuonTriggerSelector") << "No beam spot available from Event" ;
    }  
    const reco::BeamSpot& beamSpot = *beamSpotHandle;

    std::vector<int> muonIsTrigger(slimmed_muons->size(), 0);
    std::vector<int> muonIsTriggerBPark(slimmed_muons->size(), 0);
    std::vector<float> muonDR(slimmed_muons->size(),10000.);
    std::vector<float> muonDPT(slimmed_muons->size(),10000.);
    std::vector<int> loose_id(slimmed_muons->size(),0);

    std::vector<int> matched_reco_flag(slimmed_muons->size(),-1);
    std::vector<int> matched_trg_index(slimmed_muons->size(),-1);
    std::vector<float> matched_dr(slimmed_muons->size(),10000.);
    std::vector<float> matched_dpt(slimmed_muons->size(),-10000.);
    std::vector<std::vector<int>> fires;
    std::vector<std::vector<int>> prescales;
    std::vector<std::vector<float>> matcher; 
    std::vector<std::vector<float>> DR;
    std::vector<std::vector<float>> DPT;    

    // fetching the primary vertex
    const reco::Vertex & PV = vertexHandle->front();

    //std::cout << std::endl;

    //////////////////////////////
    /* method 1 */
    // proceeding to the trigger muon matching
    // for that, only use slimmedMuons
    for(const pat::Muon &muon : *slimmed_muons){
        if(debug) std::cout <<"Muon Pt="<< muon.pt() << " Eta=" << muon.eta() << " Phi=" << muon.phi() << endl;
        //std::cout <<"Muon Pt="<< muon.pt() << " Eta=" << muon.eta() << " Phi=" << muon.phi() << " nTriggerObjectMatch " << muon.triggerObjectMatches().size() << endl;

        std::vector<int> frs(HLTPaths_.size(),0); //path fires for each reco muon
        std::vector<int> prescale(HLTPaths_.size(),-1);
        //std::vector<int> sds(L1Seeds_.size(),0);// L1 Seeds for each L1 muon
        std::vector<float> temp_matched_to(HLTPaths_.size(),1000.);
        std::vector<float> temp_DR(HLTPaths_.size(),1000.);
        std::vector<float> temp_DPT(HLTPaths_.size(),1000.);
        int ipath=-1;

/*        int iseed=-1;
        for (const std::string seed: L1Seeds_){
            iseed++;
            char cstr[(seed+"*").size()+1];
            strcpy( cstr,(seed+"*").c_str());
            if(muon.triggerObjectMatches().size()!=0){
                for(size_t i=0;i<muon.triggerObjectMatches().size(); i++){
                    if(muon.triggerObjectMatch(i)!=0 && muon.triggerObjectMatch(i)->hasAlgorithmName(cstr,true)){
                        sds[iseed]=1;
                        std::cout << "L1 Seed="<< cstr <<" fired="<< sds[iseed] << endl;
                    }
                }
            }
        } */

        // for each muon, we study whether it fires a HLT line or not
        for (const std::string path: HLTPaths_){
            if(debug) std::cout << path << " " << muon.triggerObjectMatches().size() << std::endl;
            ipath++;

            // the following vectors are used in order to find the minimum DR between a reco muon and all the HLT objects that is matched to it
            // as a reco muon will be matched with only one HLT object every time, so there is a one-to-one correspondence between the two collections
            // DPt_rel is not used to create this one-to-one correspondence but only to create a few plots, debugging and be sure that everything is working fine. 
            std::vector<float> temp_dr(muon.triggerObjectMatches().size(),1000.);
            std::vector<float> temp_dpt(muon.triggerObjectMatches().size(),1000.);
            std::vector<float> temp_pt(muon.triggerObjectMatches().size(),1000.);

            // put the HLT path name into a wildcard
            char cstr[ (path+"*").size() + 1];
            strcpy( cstr, (path+"*").c_str() );       

            // Here we find all the HLT objects from each HLT path each time that are matched with the reco muon.
            if(muon.triggerObjectMatches().size()!=0){
                for(size_t i=0; i<muon.triggerObjectMatches().size();i++){
                  if(debug) std::cout << "triggerObjMatch " << i << " " << cstr << " " <<  muon.triggerObjectMatch(i)->hasPathName(cstr,true,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,false) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,true,false) << std::endl;
                    //std::cout << "trgObjMatch " << i << " " << cstr << " " <<  muon.triggerObjectMatch(i)->hasPathName(cstr,true,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,false) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,true,false) << std::endl;
                    // first bool is pathLastFilterAccepted, second is pathL3FilterAccepted
                    if(muon.triggerObjectMatch(i)!=0 && muon.triggerObjectMatch(i)->hasPathName(cstr,true,true)){
                    //if(muon.triggerObjectMatch(i)!=0 && muon.triggerObjectMatch(i)->hasPathName(cstr,false,false) && DR[iMuo][path]<max_deltaR_ && fabs(DPT[iMuo][path])<max_deltaPtRel_ && DR[iMuo][path]!=10000){
                        //if(abs(muon.triggerObjectMatch(i)->eta())>1.5) std::cout << "eta=" <<muon.triggerObjectMatch(i)->eta();
                        frs[ipath]=1;
                        float dr = reco::deltaR(muon.triggerObjectMatch(i)->p4(), muon.p4()); 
                        float dpt=(muon.triggerObjectMatch(i)->pt()-muon.pt())/muon.triggerObjectMatch(i)->pt();
                        temp_dr[i]=dr;
                        temp_dpt[i]=dpt;
                        temp_pt[i]=muon.triggerObjectMatch(i)->pt();                   
                        if(debug)std::cout <<"Path=" <<cstr << endl;
                        if(debug)std::cout <<"HLT  Pt="<<muon.triggerObjectMatch(i)->pt() <<" Eta="<<muon.triggerObjectMatch(i)->eta() <<" Phi="<<muon.triggerObjectMatch(i)->phi() << endl;
                        //std::cout <<"HLT " << cstr << " Pt="<<muon.triggerObjectMatch(i)->pt() <<" Eta="<<muon.triggerObjectMatch(i)->eta() <<" Phi="<<muon.triggerObjectMatch(i)->phi() << endl;
                        if(debug)std::cout <<"Muon Pt="<< muon.pt() << " Eta=" << muon.eta() << " Phi=" << muon.phi()  <<endl;
                        if(debug)std::cout <<"DR = " << temp_dr[i] <<endl;
                        //std::cout <<"DR = " << temp_dr[i] << " DPT = " << temp_dpt[i] <<endl;
                    }
                }
                // and now we find the real minimum between the reco muon and all its matched HLT objects. 
                temp_DR[ipath]=*min_element(temp_dr.begin(),temp_dr.end());
                int position=std::min_element(temp_dr.begin(),temp_dr.end()) - temp_dr.begin();
                temp_DPT[ipath]=temp_dpt[position];
                temp_matched_to[ipath]=temp_pt[position];
                }
            }
        //and now since we have found the minimum DR we save a few variables for plots       
        fires.push_back(frs); //This is used in order to see if a reco muon fired a Trigger (1) or not (0).
        //for(unsigned int i(0); i<frs.size(); ++i){
        //  std::cout << "fires " << frs[i] << " deltaR " << temp_DR[i] << " deltaPt " << temp_DPT[i] << std::endl;
        //}
        prescales.push_back(prescale); //This is used for the initialisation of the prescales.
        matcher.push_back(temp_matched_to); //This is used in order to see if a reco muon is matched with a HLT object. PT of the reco muon is saved in this vector. 
        DR.push_back(temp_DR);
        DPT.push_back(temp_DPT);

    }

    //now, check for different reco muons that are matched to the same HLTObject.
    for(unsigned int path=0; path<HLTPaths_.size(); path++){
        for(unsigned int iMuo=0; iMuo<slimmed_muons->size(); iMuo++){
            for(unsigned int im=(iMuo+1); im<slimmed_muons->size(); im++){
                if(matcher[iMuo][path]!=1000. && matcher[iMuo][path]==matcher[im][path]){
                //if(matcher[iMuo][path]!=1000. && matcher[iMuo][path]==matcher[im][path] && DR[iMuo][path]<max_deltaR_ && fabs(DPT[iMuo][path])<max_deltaPtRel_ && DR[im][path]<max_deltaR_ && fabs(DPT[im][path])<max_deltaPtRel_){
                    if(DR[iMuo][path]<DR[im][path]){ //Keep the one that has the minimum DR with the HLT object
                        fires[im][path]=0;
                        matcher[im][path]=1000.;
                        DR[im][path]=1000.;                       
                        DPT[im][path]=1000.;
                    }
                    else{
                        fires[iMuo][path]=0;
                        matcher[iMuo][path]=1000.;
                        DR[iMuo][path]=1000.;                       
                        DPT[iMuo][path]=1000.;
                    }
                }              
            }
            if(matcher[iMuo][path]!=1000. && DR[iMuo][path]<max_deltaR_trigger_matching_ && fabs(DPT[iMuo][path])<max_deltaPtRel_trigger_matching_ && DR[iMuo][path]!=10000){
              muonIsTrigger[iMuo]=1;
              // BParking trigger lines correspond to the 10 first element (0 to 9) of the HLTPaths vector
              if(path < 10) muonIsTriggerBPark[iMuo]=1;
              muonDR[iMuo]=DR[iMuo][path];
              muonDPT[iMuo]=DPT[iMuo][path];                
              
              for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
                //if(names.triggerName(i).find(HLTPaths_[path]) != std::string::npos || names.triggerName(i)==HLTPaths_[path]){
                if(names.triggerName(i).find(HLTPaths_[path]) != std::string::npos){
                  prescales[iMuo][path] = triggerPrescales->getPrescaleForIndex(i);
                }
              }
            }
            else{
              fires[iMuo][path]=0;
              //prescales[iMuo][path] = -1;
            }
        }
    }

    //////////////////////////////
    /* method 2 */

    /*
    // we first get the trigger objects
    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);


    // container for the trigger objects, that will be later on matched to reco muons
    std::vector<pat::TriggerObjectStandAlone> triggeringMuons;


    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackFilterLabels(iEvent, *triggerBits);
      obj.unpackPathNames(names);

      bool isTriggerMuon = false;
      for (unsigned h = 0; h < obj.filterIds().size(); ++h){
        if(obj.filterIds()[h] == 83){ 
          isTriggerMuon = true; 
          break;
        } 
      }

      if(!isTriggerMuon) continue; 
      for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
        std::string filterName = obj.filterLabels()[h];
        if(filterName.find("hltL3") != std::string::npos  && filterName.find("Park") != std::string::npos){
          isTriggerMuon = true;
          if(debug) std::cout << "\t   Filters:   " << filterName; 
          break;
        }
        else{ isTriggerMuon = false; }
      }

      if(!isTriggerMuon) continue;

      std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      // Print trigger object collection and type
      std::cout << "\t   Collection: " << obj.collection() << std::endl;
      std::cout << "\t   Type IDs:   ";

      for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
      std::cout << std::endl;

      for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
      std::cout << std::endl;

      std::vector pathNamesAll = obj.pathNames(false);
      std::vector pathNamesLast = obj.pathNames(true);

      std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";

      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
            bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
            bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
            bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
            bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
            std::cout << "   " << pathNamesAll[h];
            if (isBoth) std::cout << "(L,3)";
            if (isL3 && !isBoth) std::cout << "(*,3)";
            if (isLF && !isBoth) std::cout << "(L,*)";
            if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
      }
      std::cout << std::endl;

      triggeringMuons.push_back(obj);
      if(debug){ std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
        // Print trigger object collection and type
        std::cout << "\t   Collection: " << obj.collection() << std::endl;
      }
    }//trigger objects
    */


    if(debug)std::cout << "number of Muons=" <<slimmed_muons->size() << endl;
    // And now create a collection with trg muons from bParking line (10 first elements of HLTPaths)
    for(const pat::Muon & muon : *slimmed_muons){
        unsigned int iMuo(&muon -&(slimmed_muons->at(0)));
        if(muonIsTrigger[iMuo]==1 && (fires[iMuo][0]==1 || fires[iMuo][1]==1 || fires[iMuo][2]==1 || 
          fires[iMuo][3]==1 || fires[iMuo][4]==1 || fires[iMuo][5]==1 || fires[iMuo][6]==1 || 
          fires[iMuo][7]==1 || fires[iMuo][8]==1 || fires[iMuo][9]==1)){
            pat::Muon recoTriggerMuonCand(muon);
            trgmuons_out->emplace_back(recoTriggerMuonCand);
        }
    }

    // add the slimmed muons to the collection 
    for(const pat::Muon & slimmed_muon : *slimmed_muons){
      unsigned int iMuo(&slimmed_muon - &(slimmed_muons->at(0)));
      if(slimmed_muon.pt()<selmu_ptMin_) continue;
      if(fabs(slimmed_muon.eta())>selmu_absEtaMax_) continue;

      const reco::TransientTrack muonTT((*(slimmed_muon.bestTrack())),&(*bFieldHandle)); //sara:check,why not using inner track for muons? GM: What is this and why do we need this???
      if(!muonTT.isValid()) continue; // GM: and why do we skip this muon if muonTT is invalid? This seems to have no effect so I kept it.
      trans_muons_out->emplace_back(muonTT);

      //std::cout << "slimmed muon " << iMuo << " pt " << slimmed_muon.pt()  << " eta " << slimmed_muon.eta() << std::endl;
      
      pat::ETHMuon the_muon(slimmed_muon);
      ETHmuons_out->emplace_back(the_muon);

      for(unsigned int i=0; i<HLTPaths_.size(); i++){
        ETHmuons_out->back().addUserInt(HLTPaths_[i], fires[iMuo][i]);
        ETHmuons_out->back().addUserInt(HLTPaths_[i] + "_prescale", prescales[iMuo][i]);
      }
      ETHmuons_out->back().addUserInt("isTriggering", muonIsTrigger[iMuo]);
      ETHmuons_out->back().addUserInt("isTriggeringBPark", muonIsTriggerBPark[iMuo]);
      ETHmuons_out->back().addUserFloat("DR", muonDR[iMuo]);
      ETHmuons_out->back().addUserFloat("DPT" ,muonDPT[iMuo]);
      ETHmuons_out->back().addUserInt("isDSAMuon", 0);
      ETHmuons_out->back().addUserInt("isSlimmedMuon", 1);
      ETHmuons_out->back().addUserFloat("dz", slimmed_muon.dB(slimmed_muon.PVDZ));
      ETHmuons_out->back().addUserFloat("dzS", slimmed_muon.dB(slimmed_muon.PVDZ)/slimmed_muon.edB(slimmed_muon.PVDZ));
      ETHmuons_out->back().addUserFloat("dxy", slimmed_muon.dB(slimmed_muon.PV2D));
      ETHmuons_out->back().addUserFloat("dxyS", slimmed_muon.dB(slimmed_muon.PV2D)/slimmed_muon.edB(slimmed_muon.PV2D));

      // computing the IP with respect to the beamspot
      // first, a la R(D*)
      float dz_alaRdst = -99.;
      float dxy_BS_alaRdst = -99.;
      float dxyS_BS_alaRdst = -99.;
      if(!slimmed_muon.innerTrack().isNull()){
        auto tk = slimmed_muon.innerTrack();
        dz_alaRdst = tk->dz(PV.position());
        dxy_BS_alaRdst = tk->dxy(beamSpot);
        dxyS_BS_alaRdst = tk->dxy(beamSpot) / dxyError(*tk, beamSpot.position(), beamSpot.rotatedCovariance3D());
      }

      ETHmuons_out->back().addUserFloat("dz_alaRdst", dz_alaRdst);
      ETHmuons_out->back().addUserFloat("dxy_BS_alaRdst", dxy_BS_alaRdst);
      ETHmuons_out->back().addUserFloat("dxyS_BS_alaRdst", dxyS_BS_alaRdst);

      // then, by accounting for that the axis of the beamspot can be shifted
      // get the vertex the closest in dz to the muon
      const reco::Vertex& the_PV = PV; // initialise it to the first primary vertex
      float dist = -99;
      for(const reco::Vertex& vertex: vertices){
        // compute muon dz at this vertex
        float slimmed_muon_dz = abs(slimmed_muon.bestTrack()->dz(vertex.position()));
        if(dist == -99 || slimmed_muon_dz < dist){
          dist = slimmed_muon_dz;
          auto the_PV = vertex;
        }
      }
      auto bs_point = reco::Vertex::Point(beamSpot.x(the_PV.z()), beamSpot.y(the_PV.z()), beamSpot.z0());
      auto bs_error = beamSpot.covariance3D();
      float chi2 = 0.;
      float ndof = 0.;
      auto the_BS = reco::Vertex(bs_point, bs_error, chi2, ndof, 3);
        
      float mu_dxy_BS = -99.;
      float mu_dxyS_BS = -99.;
      if(slimmed_muon.bestTrack()){
        mu_dxy_BS = slimmed_muon.bestTrack()->dxy(the_BS.position());
        mu_dxyS_BS = slimmed_muon.bestTrack()->dxy(the_BS.position()) / dxyError(*slimmed_muon.bestTrack(), the_BS.position(), the_BS.error());
      }

      ETHmuons_out->back().addUserFloat("dxy_BS", mu_dxy_BS);
      ETHmuons_out->back().addUserFloat("dxyS_BS", mu_dxyS_BS);

      ETHmuons_out->back().addUserInt("isMatchedToSlimmedMuon", -99);
      ETHmuons_out->back().addUserInt("indexMatchedSlimmedMuon", -99);
      ETHmuons_out->back().addUserFloat("dsaToSlimmedMatching_deltaPtRel", -99.);
      ETHmuons_out->back().addUserFloat("dsaToSlimmedMatching_deltaR", -99.);
      ETHmuons_out->back().addUserFloat("dsaToSlimmedMatching_deltadxyRel", -99.);
      ETHmuons_out->back().addUserFloat("dsaToSlimmedMatching_deltadzRel", -99.);
      ETHmuons_out->back().addUserInt("passDSAMuonID", -99);
    }

    // add the displaced standalone muons to the collection
    if(add_dsa_){
      for(const reco::Track & dsa_muon : *displaced_standalone_muons){
        //unsigned int iDSAMuon(&dsa_muon - &(displaced_standalone_muons->at(0)));
        if(dsa_muon.pt()<selmu_ptMin_) continue;
        if(fabs(dsa_muon.eta())>selmu_absEtaMax_) continue;

        // add the muon to the transient track collection
        // one has to make sure that the indices in the muon and transient track collections match
        const reco::TransientTrack muonTT((*(&dsa_muon)),&(*bFieldHandle)); 
        if(!muonTT.isValid()) continue; 
        trans_muons_out->emplace_back(muonTT);

        pat::ETHMuon the_muon(dsa_muon);
        ETHmuons_out->emplace_back(the_muon);

        for(unsigned int i=0; i<HLTPaths_.size(); i++){
          ETHmuons_out->back().addUserInt(HLTPaths_[i], -1);
          ETHmuons_out->back().addUserInt(HLTPaths_[i] + "_prescale", -1);
        }
        ETHmuons_out->back().addUserInt("isTriggering", -1);
        ETHmuons_out->back().addUserInt("isTriggeringBPark", -1);
        ETHmuons_out->back().addUserFloat("DR", -1.);
        ETHmuons_out->back().addUserFloat("DPT", -1.);
        ETHmuons_out->back().addUserInt("isDSAMuon", 1);
        ETHmuons_out->back().addUserInt("isSlimmedMuon", 0);
        ETHmuons_out->back().addUserFloat("dz", dsa_muon.dz(PV.position()));
        ETHmuons_out->back().addUserFloat("dzS", dsa_muon.dz(PV.position())/dsa_muon.dzError());
        ETHmuons_out->back().addUserFloat("dxy", dsa_muon.dxy(PV.position()));
        ETHmuons_out->back().addUserFloat("dxyS", dsa_muon.dxy(PV.position())/dsa_muon.dxyError());

        //std::cout << "DSA muon " << iDSAMuon << " pt " << dsa_muon.pt()  << " eta " << dsa_muon.eta() << std::endl;

        // perform dR matching between DSA and slimmed muons
        // (without ambiguity resolving, i.e that a given slimmed muon can be matched to more than one DSA muon)
        bool isMatchedToSlimmedMuon = 0;
        int indexMatchedSlimmedMuon = -1;
        float dsaToSlimmedMatching_deltaR = -1;
        float dsaToSlimmedMatching_deltaPtRel = -1.;
        float dsaToSlimmedMatching_deltadxyRel = -1.;
        float dsaToSlimmedMatching_deltadzRel = -1.;

        if(do_dsa_matching_){
          std::vector<pair<int, std::array<float, 4>>> pairs_slimmedIdx_deltaPtRel;

          // for each DSA muon, get the slimmed muons that are within the deltaR_max cone
          for(const pat::Muon & slimmed_muon : *slimmed_muons){
            // for each slimmed  muon, apply same selection as above, to make sure that the indices are consistent
            if(slimmed_muon.pt()<selmu_ptMin_) continue;
            if(fabs(slimmed_muon.eta())>selmu_absEtaMax_) continue;
            const reco::TransientTrack muonTT((*(slimmed_muon.bestTrack())),&(*bFieldHandle));
            if(!muonTT.isValid()) continue;

            unsigned int iSlimmedMuon(&slimmed_muon - &(slimmed_muons->at(0)));

            float deltaR_dsa_slimmed = reco::deltaR(dsa_muon.eta(), dsa_muon.phi(), slimmed_muon.eta(), slimmed_muon.phi());
            if(deltaR_dsa_slimmed < max_deltaR_dsaToSlimmed_matching_){
              pair<int, std::array<float, 4>> pairs_slimmedIdx_deltaPtRel_tmp;
              pairs_slimmedIdx_deltaPtRel_tmp.first = iSlimmedMuon; //indexMatchedSlimmedMuon;
              pairs_slimmedIdx_deltaPtRel_tmp.second[0] = fabs(dsa_muon.pt() - slimmed_muon.pt()) / slimmed_muon.pt();
              pairs_slimmedIdx_deltaPtRel_tmp.second[1] = deltaR_dsa_slimmed;
              pairs_slimmedIdx_deltaPtRel_tmp.second[2] = fabs(dsa_muon.dxy(PV.position()) - slimmed_muon.dB(slimmed_muon.PV2D)) / fabs(slimmed_muon.dB(slimmed_muon.PV2D));
              pairs_slimmedIdx_deltaPtRel_tmp.second[3] = fabs(dsa_muon.dz(PV.position()) - slimmed_muon.dB(slimmed_muon.PVDZ)) / fabs(slimmed_muon.dB(slimmed_muon.PVDZ));
              pairs_slimmedIdx_deltaPtRel.push_back(pairs_slimmedIdx_deltaPtRel_tmp);
              //std::cout << "DSA muon " << iDSAMuon << " matched to slimmed muon " << pairs_slimmedIdx_deltaPtRel_tmp.first << " with deltaPtRel " << fabs(dsa_muon.pt() - slimmed_muon.pt()) / slimmed_muon.pt() << " and deltaR " << deltaR_dsa_slimmed << std::endl;
            }
          }
          // sort the pairs in deltaPtRel
          sort(pairs_slimmedIdx_deltaPtRel.begin(), pairs_slimmedIdx_deltaPtRel.end(), [](const pair<int, std::array<float, 4>> &pair_i, const pair<int, std::array<float, 4>> &pair_j){
            return pair_i.second[0] < pair_j.second[0];
          });

          // fetch the matching index
          if(pairs_slimmedIdx_deltaPtRel.size() > 0 && pairs_slimmedIdx_deltaPtRel[0].second[0] < max_deltaPtRel_dsaToSlimmed_matching_){
            isMatchedToSlimmedMuon = 1;
            indexMatchedSlimmedMuon = pairs_slimmedIdx_deltaPtRel[0].first;
            dsaToSlimmedMatching_deltaPtRel = pairs_slimmedIdx_deltaPtRel[0].second[0];
            dsaToSlimmedMatching_deltaR = pairs_slimmedIdx_deltaPtRel[0].second[1];
            dsaToSlimmedMatching_deltadxyRel = pairs_slimmedIdx_deltaPtRel[0].second[2];
            dsaToSlimmedMatching_deltadzRel = pairs_slimmedIdx_deltaPtRel[0].second[3];
            //std::cout << "DSA muon is matched " << isMatchedToSlimmedMuon << " to slimmed muon " << indexMatchedSlimmedMuon << " with deltaR " << dsaToSlimmedMatching_deltaR << " deltaPtRel " << dsaToSlimmedMatching_deltaPtRel << std::endl;
          }
        }

        ETHmuons_out->back().addUserInt("isMatchedToSlimmedMuon", isMatchedToSlimmedMuon);
        ETHmuons_out->back().addUserInt("indexMatchedSlimmedMuon", indexMatchedSlimmedMuon);
        ETHmuons_out->back().addUserFloat("dsaToSlimmedMatching_deltaPtRel", dsaToSlimmedMatching_deltaPtRel);
        ETHmuons_out->back().addUserFloat("dsaToSlimmedMatching_deltaR", dsaToSlimmedMatching_deltaR);
        ETHmuons_out->back().addUserFloat("dsaToSlimmedMatching_deltadxyRel", dsaToSlimmedMatching_deltadxyRel);
        ETHmuons_out->back().addUserFloat("dsaToSlimmedMatching_deltadzRel", dsaToSlimmedMatching_deltadzRel);

        // building DSA muon ID
        bool passDSAMuonID = false;
        if(
          (
           (dsa_muon.hitPattern().numberOfValidMuonDTHits() > 18 && dsa_muon.hitPattern().numberOfValidMuonCSCHits() == 0) ||
           (dsa_muon.hitPattern().numberOfValidMuonHits() > 12 && dsa_muon.hitPattern().numberOfValidMuonCSCHits() > 0)
          ) &&
          dsa_muon.normalizedChi2() < 2.5 &&
          dsa_muon.ptError() / dsa_muon.pt() < 1 && 
          dsa_muon.hitPattern().muonStationsWithValidHits() > 1
          ){
          passDSAMuonID = true;
        }
        ETHmuons_out->back().addUserInt("passDSAMuonID", passDSAMuonID);
      }
    }

    //std::cout << std::endl;

    iEvent.put(std::move(trgmuons_out), "trgMuons"); 
    iEvent.put(std::move(ETHmuons_out), "SelectedMuons");
    iEvent.put(std::move(trans_muons_out), "SelectedTransientMuons");
}



DEFINE_FWK_MODULE(MuonTriggerSelector);
