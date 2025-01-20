#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "../interface/helper.h"
#include <limits>
#include <algorithm>
#include "../interface/KinVtxFitter.h"
#include "../interface/ETHMuon.h"

class BToKLLBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BToKLLBuilder(const edm::ParameterSet &cfg):
    k_selection_{cfg.getParameter<std::string>("kaonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    isMC_{cfg.getParameter<bool>("isMC")},
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    kaons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("kaons") )},
    kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    genParticles_{consumes<reco::GenParticleCollection>( cfg.getParameter<edm::InputTag>("genParticles"))}, 
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>();
    }

  ~BToKLLBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> k_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const bool isMC_;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> kaons_;
  const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 
  
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BToKLLBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  evt.getByToken(dileptons_, dileptons);
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> kaons;
  evt.getByToken(kaons_, kaons);
  
  edm::Handle<TransientTrackCollection> kaons_ttracks;
  evt.getByToken(kaons_ttracks_, kaons_ttracks);  
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genParticles_, genParticles);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);

  std::vector<int> used_lep1_id, used_lep2_id, used_trk_id;

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx) {
    edm::Ptr<pat::CompositeCandidate> k_ptr(kaons, k_idx);
    if( !k_selection_(*k_ptr) ) continue;
    
    math::PtEtaPhiMLorentzVector k_p4(
      k_ptr->pt(), 
      k_ptr->eta(),
      k_ptr->phi(),
      K_MASS
      );

    for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
      edm::Ptr<pat::CompositeCandidate> ll_ptr(dileptons, ll_idx);
      edm::Ptr<reco::Candidate> l1_ptr = ll_ptr->userCand("l1");
      edm::Ptr<reco::Candidate> l2_ptr = ll_ptr->userCand("l2");
      int l1_idx = ll_ptr->userInt("l1_idx");
      int l2_idx = ll_ptr->userInt("l2_idx");
    
      math::PtEtaPhiMLorentzVector ll_p4(
        ll_ptr->pt(), 
        ll_ptr->eta(),
        ll_ptr->phi(),
        ll_ptr->mass()
        );

      pat::CompositeCandidate cand;
      cand.setP4(ll_p4 + k_p4);
      cand.setCharge(ll_ptr->charge() + k_ptr->charge());
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the lepton passing the corresponding selection
      cand.addUserCand("l1", l1_ptr);
      cand.addUserCand("l2", l2_ptr);
      cand.addUserCand("K", k_ptr);
      cand.addUserCand("dilepton", ll_ptr);

      cand.addUserInt("l1_idx", l1_idx);
      cand.addUserInt("l2_idx", l2_idx);
      cand.addUserInt("k_idx", k_idx);
    
      auto dr_info = min_max_dr({l1_ptr, l2_ptr, k_ptr});
      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);
      // TODO add meaningful variables
      
      if( !pre_vtx_selection_(cand) ) continue;
    
      KinVtxFitter fitter(
        {leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx), kaons_ttracks->at(k_idx)},
        {l1_ptr->mass(), l2_ptr->mass(), K_MASS},
        {LEP_SIGMA, LEP_SIGMA, K_SIGMA} //some small sigma for the lepton mass
        );
      if(!fitter.success()) continue;

      cand.setVertex( 
        reco::Candidate::Point( 
          fitter.fitted_vtx().x(),
          fitter.fitted_vtx().y(),
          fitter.fitted_vtx().z()
          )  
        );
      used_lep1_id.emplace_back(l1_idx);
      used_lep2_id.emplace_back(l2_idx);
      used_trk_id.emplace_back(k_idx);
      cand.addUserInt("sv_OK" , fitter.success());
      cand.addUserFloat("sv_chi2", fitter.chi2());
      cand.addUserFloat("sv_ndof", fitter.dof()); // float??
      cand.addUserFloat("sv_prob", fitter.prob());
      cand.addUserFloat("fitted_mll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass());
      auto fit_p4 = fitter.fitted_p4();
      cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
      cand.addUserFloat("fitted_eta" , fit_p4.eta());
      cand.addUserFloat("fitted_phi" , fit_p4.phi());
      cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());      
      cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
      cand.addUserFloat(
        "cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, cand.p4())
        );
      cand.addUserFloat(
        "fitted_cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, fit_p4)
        );
      auto lxy = l_xy(fitter, *beamspot);
      cand.addUserFloat("l_xy", lxy.value());
      cand.addUserFloat("l_xy_unc", lxy.error());
      cand.addUserFloat("vtx_x", cand.vx());
      cand.addUserFloat("vtx_y", cand.vy());
      cand.addUserFloat("vtx_z", cand.vz());
      cand.addUserFloat("vtx_ex", sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      cand.addUserFloat("vtx_ey", sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      cand.addUserFloat("vtx_ez", sqrt(fitter.fitted_vtx_uncertainty().czz()));

      cand.addUserFloat("fitted_l1_pt" , fitter.daughter_p4(0).pt()); 
      cand.addUserFloat("fitted_l1_eta", fitter.daughter_p4(0).eta());
      cand.addUserFloat("fitted_l1_phi", fitter.daughter_p4(0).phi());
      cand.addUserFloat("fitted_l2_pt" , fitter.daughter_p4(1).pt()); 
      cand.addUserFloat("fitted_l2_eta", fitter.daughter_p4(1).eta());
      cand.addUserFloat("fitted_l2_phi", fitter.daughter_p4(1).phi());
      cand.addUserFloat("fitted_k_pt"  , fitter.daughter_p4(2).pt()); 
      cand.addUserFloat("fitted_k_eta" , fitter.daughter_p4(2).eta());
      cand.addUserFloat("fitted_k_phi" , fitter.daughter_p4(2).phi());
   
      cand.addUserFloat("ll_sv_prob", ll_ptr->userFloat("sv_prob"));

      if( !post_vtx_selection_(cand) ) continue;        

      // add kaon information
      cand.addUserFloat("kaon_pt"                  , k_ptr->pt()                             );
      cand.addUserFloat("kaon_eta"                 , k_ptr->eta()                            );
      cand.addUserFloat("kaon_phi"                 , k_ptr->phi()                            );
      cand.addUserFloat("kaon_mass"                , k_ptr->mass()                           );
      cand.addUserInt("kaon_charge"                , k_ptr->charge()                         );
      cand.addUserInt("kaon_pdgId"                 , k_ptr->pdgId()                          );
      cand.addUserFloat("kaon_vx"                  , k_ptr->vx()                             );
      cand.addUserFloat("kaon_vy"                  , k_ptr->vy()                             );
      cand.addUserFloat("kaon_vz"                  , k_ptr->vz()                             );
      cand.addUserFloat("kaon_dz"                  , k_ptr->userFloat("dz")                  );
      cand.addUserFloat("kaon_dxy"                 , k_ptr->userFloat("dxy")                 );
      cand.addUserFloat("kaon_dzS"                 , k_ptr->userFloat("dzS")                 );
      cand.addUserFloat("kaon_dxyS"                , k_ptr->userFloat("dxyS")                );
      // additionnal quantities
      cand.addUserFloat("kaon_DCASig"              , k_ptr->userFloat("DCASig")              );
      cand.addUserFloat("kaon_DCASig_corr"         , k_ptr->userFloat("DCASig_corr")         );
      cand.addUserInt("kaon_ispacked"              , k_ptr->userInt("isPacked")              );
      cand.addUserInt("kaon_islost"                , k_ptr->userInt("isLostTrk")             );
      cand.addUserFloat("kaon_chi2"                , k_ptr->userFloat("chi2")                );
      cand.addUserFloat("kaon_normalisedChi2"      , k_ptr->userFloat("normalisedChi2")      );
      cand.addUserFloat("kaon_validFraction"       , k_ptr->userFloat("validFraction")       );
      cand.addUserInt("kaon_ndof"                  , k_ptr->userInt("ndof")                  );
      cand.addUserInt("kaon_numberOfValidHits"     , k_ptr->userInt("nValidHits")            );
      cand.addUserInt("kaon_numberOfLostHits"      , k_ptr->userInt("numberOfLostHits")      );
      cand.addUserInt("kaon_numberOfValidPixelHits", k_ptr->userInt("numberOfValidPixelHits"));
      cand.addUserInt("kaon_numberOfTrackerLayers" , k_ptr->userInt("numberOfTrackerLayers") );
      cand.addUserInt("kaon_numberOfPixelLayers"   , k_ptr->userInt("numberOfPixelLayers")   );
      cand.addUserInt("kaon_qualityIndex"          , k_ptr->userInt("qualityIndex")          );
      cand.addUserInt("kaon_highPurityFlag"        , k_ptr->userInt("highPurityFlag")        );

      //compute isolation
      float l1_iso03 = 0;
      float l1_iso04 = 0;
      float l2_iso03 = 0;
      float l2_iso04 = 0;
      float k_iso03  = 0;
      float k_iso04  = 0;
      float b_iso03  = 0;
      float b_iso04  = 0;
      // with conditions: best track + close to B z-vertex
      float l1_iso03_close = 0; 
      float l1_iso04_close = 0; 
      float l2_iso03_close = 0; 
      float l2_iso04_close = 0; 
      float k_iso03_close = 0; 
      float k_iso04_close = 0; 
      float b_iso03_close = 0; 
      float b_iso04_close = 0; 

      /*
      unsigned int nTracks     = iso_tracks->size();
      unsigned int totalTracks = nTracks + iso_lostTracks->size();

      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
      
        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if( !isotrk_selection_(trk) ) continue;
        // check if the track is the kaon
        if (k_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        // check if the track is one of the two leptons
        if (track_to_lepton_match(l1_ptr, iso_tracks.id(), iTrk) || 
            track_to_lepton_match(l2_ptr, iso_tracks.id(), iTrk) ) continue;

        // add to final particle iso if dR < cone
        float dr_to_l1 = deltaR(cand.userFloat("fitted_l1_eta"), cand.userFloat("fitted_l1_phi"), trk.eta(), trk.phi());
        float dr_to_l2 = deltaR(cand.userFloat("fitted_l2_eta"), cand.userFloat("fitted_l2_phi"), trk.eta(), trk.phi());
        float dr_to_k  = deltaR(cand.userFloat("fitted_k_eta") , cand.userFloat("fitted_k_phi") , trk.eta(), trk.phi());
        float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());

        if (dr_to_l1 < 0.4){
          l1_iso04 += trk.pt();
          if ( dr_to_l1 < 0.3) l1_iso03 += trk.pt();
        }
        if (dr_to_l2 < 0.4){
          l2_iso04 += trk.pt();
          if (dr_to_l2 < 0.3)  l2_iso03 += trk.pt();
        }
        if (dr_to_k < 0.4){
          k_iso04 += trk.pt();
          if (dr_to_k < 0.3) k_iso03 += trk.pt();
        }
        if (dr_to_b < 0.4){
          b_iso04 += trk.pt();
          if (dr_to_b < 0.3) b_iso03 += trk.pt();
        }
        
        // add requirement of the tracks to be close to the B
        if (!l1_ptr->bestTrack() || fabs(trk.dz() - l1_ptr->bestTrack()->dz()) > 0.4) continue;
        if (!l2_ptr->bestTrack() || fabs(trk.dz() - l2_ptr->bestTrack()->dz()) > 0.4) continue;
        if (fabs(trk.dz() - k_ptr->userFloat("dz")) > 0.4) continue;

        if (dr_to_l1 < 0.4){
          l1_iso04_close += trk.pt();
          if ( dr_to_l1 < 0.3) l1_iso03_close += trk.pt();
        }
        if (dr_to_l2 < 0.4){
          l2_iso04_close += trk.pt();
          if ( dr_to_l2 < 0.3) l2_iso03_close += trk.pt();
        }
        if (dr_to_k < 0.4){
          k_iso04_close += trk.pt();
          if ( dr_to_k < 0.3) k_iso03_close += trk.pt();
        }
        if (dr_to_b < 0.4){
          b_iso04_close += trk.pt();
          if (dr_to_b < 0.3) b_iso03_close += trk.pt();
        }
      }
      */

      cand.addUserFloat("l1_iso03", l1_iso03);
      cand.addUserFloat("l1_iso04", l1_iso04);
      cand.addUserFloat("l2_iso03", l2_iso03);
      cand.addUserFloat("l2_iso04", l2_iso04);
      cand.addUserFloat("k_iso03" , k_iso03 );
      cand.addUserFloat("k_iso04" , k_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );
      
      cand.addUserFloat("l1_iso03_close", l1_iso03_close);
      cand.addUserFloat("l1_iso04_close", l1_iso04_close);
      cand.addUserFloat("l2_iso03_close", l2_iso03_close);
      cand.addUserFloat("l2_iso04_close", l2_iso04_close);
      cand.addUserFloat("k_iso03_close", k_iso03_close);
      cand.addUserFloat("k_iso04_close", k_iso04_close);
      cand.addUserFloat("b_iso03_close", b_iso03_close);
      cand.addUserFloat("b_iso04_close", b_iso04_close);

      // gen-matching
      int isMatched = 0;
      int l1_genIdx(-1), l2_genIdx(-1), k_genIdx(-1);
      int genMuon1Mother_genPdgId(-1), genMuon2Mother_genPdgId(-1), genKaonMother_genPdgId(-1);
      float matched_b_pt(-99.), matched_b_eta(-99.), matched_b_phi(-99.), matched_b_mass(-99.);
      float matched_l1_pt(-99.), matched_l1_eta(-99.), matched_l1_phi(-99.), matched_l1_mass(-99.);
      float matched_l2_pt(-99.), matched_l2_eta(-99.), matched_l2_phi(-99.), matched_l2_mass(-99.);
      float matched_k_pt(-99.), matched_k_eta(-99.), matched_k_phi(-99.), matched_k_mass(-99.);

      // for MC only
      if(isMC_ == true){

        // pdgId of the gen particle to which the final-state particles are matched
        int l1_genPdgId = ll_ptr->userInt("l1_mcMatch");
        int l2_genPdgId = ll_ptr->userInt("l2_mcMatch");
        int k_genPdgId  = k_ptr->userInt("mcMatch");
        
        // index of the gen particle to which the final-state particles are matched
        l1_genIdx = ll_ptr->userInt("l1_mcMatchIndex"); 
        l2_genIdx = ll_ptr->userInt("l2_mcMatchIndex"); 
        k_genIdx  = k_ptr->userInt("mcMatchIndex"); 

        if(l1_genIdx != -1 && l2_genIdx != -1 && k_genIdx != -1){

          // getting the associated gen particles
          edm::Ptr<reco::GenParticle> genMuon1_ptr(genParticles, l1_genIdx);
          edm::Ptr<reco::GenParticle> genMuon2_ptr(genParticles, l2_genIdx);
          edm::Ptr<reco::GenParticle> genKaon_ptr(genParticles, k_genIdx);

          // index of the associated mother particle
          int genMuon1Mother_genIdx = -1;
          int genMuon2Mother_genIdx = -1;
          int genKaonMother_genIdx  = -1;
          if(genMuon1_ptr->numberOfMothers()>0) genMuon1Mother_genIdx = genMuon1_ptr->motherRef(0).key();
          if(genMuon2_ptr->numberOfMothers()>0) genMuon2Mother_genIdx = genMuon2_ptr->motherRef(0).key();
          if(genKaon_ptr->numberOfMothers()>0) genKaonMother_genIdx = genKaon_ptr->motherRef(0).key();

          // getting the mother particles
          edm::Ptr<reco::GenParticle> genMuon1Mother_ptr(genParticles, genMuon1Mother_genIdx);
          edm::Ptr<reco::GenParticle> genMuon2Mother_ptr(genParticles, genMuon2Mother_genIdx);
          edm::Ptr<reco::GenParticle> genKaonMother_ptr(genParticles, genKaonMother_genIdx);

          // pdgId of the mother particles
          genMuon1Mother_genPdgId = genMuon1Mother_ptr->pdgId();
          genMuon2Mother_genPdgId = genMuon2Mother_ptr->pdgId();
          genKaonMother_genPdgId  = genKaonMother_ptr->pdgId();

          // grand-mothers of the muons
          int genMuon1GMother_genIdx = (genMuon1Mother_ptr->motherRef(0).key()) ? genMuon1Mother_ptr->motherRef(0).key() : -1;
          int genMuon2GMother_genIdx = (genMuon2Mother_ptr->motherRef(0).key()) ? genMuon2Mother_ptr->motherRef(0).key() : -1;
          edm::Ptr<reco::GenParticle> genMuon1GMother_ptr(genParticles, genMuon1GMother_genIdx);
          edm::Ptr<reco::GenParticle> genMuon2GMother_ptr(genParticles, genMuon2GMother_genIdx);
          int genMuon1GMother_genPdgId = genMuon1GMother_ptr->pdgId();
          int genMuon2GMother_genPdgId = genMuon2GMother_ptr->pdgId();
             
          // matching: muons fully matched, kaon not necessarily
          if(
             fabs(l1_genPdgId) == 13 && fabs(genMuon1Mother_genPdgId) == 443 && 
             fabs(l2_genPdgId) == 13 && fabs(genMuon2Mother_genPdgId) == 443 && 
             (
               (
                 (  
                    (fabs(k_genPdgId) == 321 && fabs(genKaonMother_genPdgId) == 521) || 
                    (fabs(genMuon1GMother_genPdgId) == fabs(genMuon2GMother_genPdgId) && fabs(genMuon2GMother_genPdgId)==521) 
                 )
               ) || (
                 (  
                    (fabs(k_genPdgId) == 211 && fabs(genKaonMother_genPdgId) == 541) || 
                    (fabs(genMuon1GMother_genPdgId) == fabs(genMuon2GMother_genPdgId) && fabs(genMuon2GMother_genPdgId)==541) 
                 )
               )
             )
            ){
              isMatched = 1;
              matched_b_pt = genKaonMother_ptr->pt();
              matched_b_eta = genKaonMother_ptr->eta();
              matched_b_phi = genKaonMother_ptr->phi();
              matched_b_mass = genKaonMother_ptr->mass();
              matched_l1_pt = genMuon1_ptr->pt();
              matched_l1_eta = genMuon1_ptr->eta();
              matched_l1_phi = genMuon1_ptr->phi();
              matched_l1_mass = genMuon1_ptr->mass();
              matched_l2_pt = genMuon2_ptr->pt();
              matched_l2_eta = genMuon2_ptr->eta();
              matched_l2_phi = genMuon2_ptr->phi();
              matched_l2_mass = genMuon2_ptr->mass();
              matched_k_pt = genKaon_ptr->pt();
              matched_k_eta = genKaon_ptr->eta();
              matched_k_phi = genKaon_ptr->phi();
              matched_k_mass = genKaon_ptr->mass();
               
          }
        }
      }

      cand.addUserInt("isMatched", isMatched);
      cand.addUserInt("matching_l1_genIdx", l1_genIdx);
      cand.addUserInt("matching_l2_genIdx", l2_genIdx);
      cand.addUserInt("matching_k_genIdx", k_genIdx);
      cand.addUserInt("matching_l1_motherPdgId", genMuon1Mother_genPdgId);
      cand.addUserInt("matching_l2_motherPdgId", genMuon2Mother_genPdgId);
      cand.addUserInt("matching_k_motherPdgId", genKaonMother_genPdgId);
      cand.addUserFloat("matched_b_pt", matched_b_pt);
      cand.addUserFloat("matched_b_eta", matched_b_eta);
      cand.addUserFloat("matched_b_phi", matched_b_phi);
      cand.addUserFloat("matched_b_mass", matched_b_mass);
      cand.addUserFloat("matched_l1_pt", matched_l1_pt);
      cand.addUserFloat("matched_l1_eta", matched_l1_eta);
      cand.addUserFloat("matched_l1_phi", matched_l1_phi);
      cand.addUserFloat("matched_l1_mass", matched_l1_mass);
      cand.addUserFloat("matched_l2_pt", matched_l2_pt);
      cand.addUserFloat("matched_l2_eta", matched_l2_eta);
      cand.addUserFloat("matched_l2_phi", matched_l2_phi);
      cand.addUserFloat("matched_l2_mass", matched_l2_mass);
      cand.addUserFloat("matched_k_pt", matched_k_pt);
      cand.addUserFloat("matched_k_eta", matched_k_eta);
      cand.addUserFloat("matched_k_phi", matched_k_phi);
      cand.addUserFloat("matched_k_mass", matched_k_mass);


      ret_val->push_back(cand);

    } // for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
  } // for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx)

  for (auto & cand: *ret_val){
    cand.addUserInt("n_k_used", std::count(used_trk_id.begin(),used_trk_id.end(),cand.userInt("k_idx")));
    cand.addUserInt("n_l1_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l1_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l1_idx")));
    cand.addUserInt("n_l2_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l2_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l2_idx")));
  }

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToKLLBuilder);
