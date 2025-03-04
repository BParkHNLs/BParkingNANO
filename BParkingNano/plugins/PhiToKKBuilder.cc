#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "../interface/helper.h"
#include <limits>
#include <algorithm>
#include "../interface/KinVtxFitter.h"

class PhiToKKBuilder : public edm::global::EDProducer<> {

public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit PhiToKKBuilder(const edm::ParameterSet &cfg):
    // kaon tracks
    kaons_ {consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("kaons"))},
    kaons_ttracks_ {consumes<TransientTrackCollection>(cfg.getParameter<edm::InputTag>("kaonsTransientTracks"))},

    // selection phi candidate
    k1_selection_ {cfg.getParameter<std::string>("k1_selection")},
    k2_selection_ {cfg.getParameter<std::string>("k2_selection")},
    pre_vtx_selection_phi_ {cfg.getParameter<std::string>("pre_vtx_selection_phi")},
    post_vtx_selection_phi_ {cfg.getParameter<std::string>("post_vtx_selection_phi")},

    //// selection second phi candidate (phi_2)
    //k3_selection_ {cfg.getParameter<std::string>("k3_selection")},
    //k4_selection_ {cfg.getParameter<std::string>("k4_selection")},
    //pre_vtx_selection_phi2_ {cfg.getParameter<std::string>("pre_vtx_selection_phi2")},
    //post_vtx_selection_phi2_ {cfg.getParameter<std::string>("post_vtx_selection_phi2")},

    //// selection Bs candidate
    //pre_vtx_selection_Bs_ {cfg.getParameter<std::string>("pre_vtx_selection_Bs")},
    //post_vtx_selection_Bs_ {cfg.getParameter<std::string>("post_vtx_selection_Bs")},

    // gen-particles
    genParticles_ {consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("genParticles"))}, 

    isMC_ {cfg.getParameter<bool>("isMC")},

    // beamspot
    beamspot_{consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))}{ 
       produces<pat::CompositeCandidateCollection>();
    }

  ~PhiToKKBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  // kaon tracks 
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> kaons_;
  const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;

  // selection phi candidate
  const StringCutObjectSelector<pat::CompositeCandidate> k1_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> k2_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_phi_;
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_phi_; 

  //// selection second phi candidate (phi_2)
  //const StringCutObjectSelector<pat::CompositeCandidate> k3_selection_; 
  //const StringCutObjectSelector<pat::CompositeCandidate> k4_selection_; 
  //const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_phi2_;
  //const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_phi2_; 
  //  
  //// selection Bs candidate
  //const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_Bs_;
  //const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_Bs_; 

  // gen-particles
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;

  const bool isMC_;

  // beamspot
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
  
};

void PhiToKKBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //TODO study strategy for tag muon
  
  // input
  edm::Handle<pat::CompositeCandidateCollection> kaons;
  evt.getByToken(kaons_, kaons);
  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(kaons_ttracks_, ttracks);

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genParticles_, genParticles);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());

  // needed to sort in pt
  std::vector<pat::CompositeCandidate> vector_phi_candidates;
  
  // first, build phi candidate
  for(size_t k1_idx = 0; k1_idx < kaons->size(); ++k1_idx) {
    edm::Ptr<pat::CompositeCandidate> k1_ptr(kaons, k1_idx);

    // selection on first kaon track
    if(!k1_selection_(*k1_ptr)) continue; //TODO would be nice to have track ID there 

    math::PtEtaPhiMLorentzVector k1_p4(
      k1_ptr->pt(), 
      k1_ptr->eta(),
      k1_ptr->phi(),
      K_MASS // impose kaon mass
      );
    
    for(size_t k2_idx = k1_idx + 1; k2_idx < kaons->size(); ++k2_idx) {
      edm::Ptr<pat::CompositeCandidate> k2_ptr(kaons, k2_idx);

      // only consider kaon2 with charge opposite to that of kaon1 (phi is neutral)
      if(k2_ptr->charge() == k1_ptr->charge()) continue;

      // then apply selection on second kaon track
      if(!k2_selection_(*k2_ptr)) continue;

      math::PtEtaPhiMLorentzVector k2_p4(
        k2_ptr->pt(), 
        k2_ptr->eta(),
        k2_ptr->phi(),
        K_MASS // impose kaon mass
        );

      pat::CompositeCandidate phi_cand;
      phi_cand.setP4(k1_p4 + k2_p4);
      phi_cand.setCharge(k1_ptr->charge() + k2_ptr->charge());
      
      phi_cand.addUserFloat("deltaR_prefit", reco::deltaR(*k1_ptr, *k2_ptr));
      
      // save indices
      phi_cand.addUserInt("k1_idx", k1_idx);
      phi_cand.addUserInt("k2_idx", k2_idx);

      // save candidates
      phi_cand.addUserCand("k1", k1_ptr);
      phi_cand.addUserCand("k2", k2_ptr);

      if(isMC_){
        phi_cand.addUserInt("k1_mcMatch", k1_ptr->userInt("mcMatch"));
        phi_cand.addUserInt("k2_mcMatch", k2_ptr->userInt("mcMatch"));
        phi_cand.addUserInt("k1_mcMatchIndex", k1_ptr->userInt("mcMatchIndex"));
        phi_cand.addUserInt("k2_mcMatchIndex", k2_ptr->userInt("mcMatchIndex"));
      }

      // apply pre-fit selection on first kaon pair
      if( !pre_vtx_selection_phi_(phi_cand) ) continue;

      // fit the kaon tracks to a common vertex
      //TODO use mass-constrained fit?
      KinVtxFitter fitter_phi(
        {ttracks->at(k1_idx), ttracks->at(k2_idx)},
        {K_MASS, K_MASS}, // force kaon mass
        {K_SIGMA, K_SIGMA}
        );
        //1.0195
        //);

      if(!fitter_phi.success()) continue;

      phi_cand.setVertex( 
        reco::Candidate::Point( 
          fitter_phi.fitted_vtx().x(),
          fitter_phi.fitted_vtx().y(),
          fitter_phi.fitted_vtx().z()
          )  
        );

      phi_cand.addUserFloat("phi_vx", fitter_phi.fitted_vtx().x());
      phi_cand.addUserFloat("phi_vy", fitter_phi.fitted_vtx().y());
      phi_cand.addUserFloat("phi_vz", fitter_phi.fitted_vtx().z());

      const auto& covMatrix = fitter_phi.fitted_vtx_uncertainty();
      phi_cand.addUserFloat("phi_vtx_cxx", covMatrix.cxx());
      phi_cand.addUserFloat("phi_vtx_cyy", covMatrix.cyy());
      phi_cand.addUserFloat("phi_vtx_czz", covMatrix.czz());
      phi_cand.addUserFloat("phi_vtx_cyx", covMatrix.cyx());
      phi_cand.addUserFloat("phi_vtx_czx", covMatrix.czx());
      phi_cand.addUserFloat("phi_vtx_czy", covMatrix.czy());

      auto phi_fit_p4 = fitter_phi.fitted_p4();
      auto phi_lxy = l_xy(fitter_phi, *beamspot);

      phi_cand.addUserFloat("phi_sv_chi2", fitter_phi.chi2());
      phi_cand.addUserFloat("phi_sv_ndof", fitter_phi.dof()); 
      phi_cand.addUserFloat("phi_sv_prob", fitter_phi.prob());
      phi_cand.addUserFloat("phi_fitted_mass", fitter_phi.success() ? fitter_phi.fitted_candidate().mass() : -1);
      phi_cand.addUserFloat("phi_fitted_massErr", fitter_phi.success() ? sqrt(fitter_phi.fitted_candidate().kinematicParametersError().matrix()(6,6)) : -1);
      phi_cand.addUserFloat("phi_charge", phi_cand.charge());

      phi_cand.addUserFloat("deltaR_postfit", reco::deltaR(fitter_phi.daughter_p4(0), fitter_phi.daughter_p4(1)));
        
      phi_cand.addUserFloat("phi_lxy", phi_lxy.value());
      phi_cand.addUserFloat("phi_lxy_sig", phi_lxy.value() / phi_lxy.error());

      phi_cand.addUserFloat("phi_fitted_pt", phi_fit_p4.pt()); 
      phi_cand.addUserFloat("phi_fitted_eta", phi_fit_p4.eta());
      phi_cand.addUserFloat("phi_fitted_phi", phi_fit_p4.phi());
      phi_cand.addUserFloat("phi_cos_theta_2D", cos_theta_2D(fitter_phi, *beamspot, phi_fit_p4));
      
      phi_cand.addUserFloat("phi_fitted_k1_pt", fitter_phi.daughter_p4(0).pt()); 
      phi_cand.addUserFloat("phi_fitted_k2_pt", fitter_phi.daughter_p4(1).pt()); 

      // apply pots-fit selection on phi candidate
      //std::cout << "before postfit selection" << std::endl;
      if(!post_vtx_selection_phi_(phi_cand)) continue; //TODO move to after definition of user floats?
        
      // save information on kaons 
      //phi_cand.addUserFloat("phi_fitted_k1_pt", fitter_phi.daughter_p4(0).pt()); 
      phi_cand.addUserFloat("phi_fitted_k1_eta", fitter_phi.daughter_p4(0).eta());
      phi_cand.addUserFloat("phi_fitted_k1_phi", fitter_phi.daughter_p4(0).phi());
      phi_cand.addUserFloat("phi_fitted_k1_mass", fitter_phi.daughter_p4(0).mass());
      phi_cand.addUserFloat("phi_k1_charge", k1_ptr->charge());
      phi_cand.addUserFloat("phi_k1_pdgId", k1_ptr->pdgId());
      phi_cand.addUserFloat("phi_k1_vx", k1_ptr->vx());
      phi_cand.addUserFloat("phi_k1_vy", k1_ptr->vy());
      phi_cand.addUserFloat("phi_k1_vz", k1_ptr->vz());
      phi_cand.addUserFloat("phi_k1_dz", k1_ptr->userFloat("dz"));
      phi_cand.addUserFloat("phi_k1_dxy", k1_ptr->userFloat("dxy"));
      phi_cand.addUserFloat("phi_k1_dzS", k1_ptr->userFloat("dzS"));
      phi_cand.addUserFloat("phi_k1_dxyS", k1_ptr->userFloat("dxyS"));
      phi_cand.addUserFloat("phi_k1_DCASig", k1_ptr->userFloat("DCASig_corr"));
      phi_cand.addUserFloat("phi_k1_ptErr", k1_ptr->userFloat("ptErr"));
      phi_cand.addUserInt("phi_k1_ispacked", k1_ptr->userInt("isPacked"));
      phi_cand.addUserInt("phi_k1_islost", k1_ptr->userInt("isLostTrk"));
      phi_cand.addUserFloat("phi_k1_chi2", k1_ptr->userFloat("chi2"));
      phi_cand.addUserFloat("phi_k1_normalisedChi2", k1_ptr->userFloat("normalisedChi2"));
      phi_cand.addUserFloat("phi_k1_validFraction", k1_ptr->userFloat("validFraction"));
      phi_cand.addUserInt("phi_k1_ndof", k1_ptr->userInt("ndof"));
      phi_cand.addUserInt("phi_k1_numberOfValidHits", k1_ptr->userInt("nValidHits"));
      phi_cand.addUserInt("phi_k1_numberOfPixelHits", k1_ptr->userInt("numberOfPixelHits"));
      phi_cand.addUserInt("phi_k1_numberOfLostHits", k1_ptr->userInt("numberOfLostHits"));
      phi_cand.addUserInt("phi_k1_numberOfValidPixelHits", k1_ptr->userInt("numberOfValidPixelHits"));
      phi_cand.addUserInt("phi_k1_numberOfTrackerLayers", k1_ptr->userInt("numberOfTrackerLayers"));
      phi_cand.addUserInt("phi_k1_numberOfPixelLayers", k1_ptr->userInt("numberOfPixelLayers"));
      phi_cand.addUserInt("phi_k1_qualityIndex", k1_ptr->userInt("qualityIndex"));
      phi_cand.addUserInt("phi_k1_highPurityFlag", k1_ptr->userInt("highPurityFlag"));
      phi_cand.addUserFloat("phi_k1_covQopQop", k1_ptr->userFloat("covQopQop"));
      phi_cand.addUserFloat("phi_k1_covLamLam", k1_ptr->userFloat("covLamLam"));
      phi_cand.addUserFloat("phi_k1_covPhiPhi", k1_ptr->userFloat("covPhiPhi"));
      phi_cand.addUserFloat("phi_k1_covQopLam", k1_ptr->userFloat("covQopLam"));
      phi_cand.addUserFloat("phi_k1_covQopPhi", k1_ptr->userFloat("covQopPhi"));
      phi_cand.addUserFloat("phi_k1_covLamPhi", k1_ptr->userFloat("covLamPhi"));

      //phi_cand.addUserFloat("phi_fitted_k2_pt", fitter_phi.daughter_p4(1).pt()); 
      phi_cand.addUserFloat("phi_fitted_k2_eta", fitter_phi.daughter_p4(1).eta());
      phi_cand.addUserFloat("phi_fitted_k2_phi", fitter_phi.daughter_p4(1).phi());
      phi_cand.addUserFloat("phi_fitted_k2_mass", fitter_phi.daughter_p4(1).mass());
      phi_cand.addUserFloat("phi_k2_charge", k2_ptr->charge());
      phi_cand.addUserFloat("phi_k2_pdgId", k2_ptr->pdgId());
      phi_cand.addUserFloat("phi_k2_vx", k2_ptr->vx());
      phi_cand.addUserFloat("phi_k2_vy", k2_ptr->vy());
      phi_cand.addUserFloat("phi_k2_vz", k2_ptr->vz());
      phi_cand.addUserFloat("phi_k2_dz", k2_ptr->userFloat("dz"));
      phi_cand.addUserFloat("phi_k2_dxy", k2_ptr->userFloat("dxy"));
      phi_cand.addUserFloat("phi_k2_dzS", k2_ptr->userFloat("dzS"));
      phi_cand.addUserFloat("phi_k2_dxyS", k2_ptr->userFloat("dxyS"));
      phi_cand.addUserFloat("phi_k2_DCASig", k2_ptr->userFloat("DCASig_corr"));
      phi_cand.addUserFloat("phi_k2_ptErr", k2_ptr->userFloat("ptErr"));
      phi_cand.addUserInt("phi_k2_ispacked", k2_ptr->userInt("isPacked"));
      phi_cand.addUserInt("phi_k2_islost", k2_ptr->userInt("isLostTrk"));
      phi_cand.addUserFloat("phi_k2_chi2", k2_ptr->userFloat("chi2"));
      phi_cand.addUserFloat("phi_k2_normalisedChi2", k2_ptr->userFloat("normalisedChi2"));
      phi_cand.addUserFloat("phi_k2_validFraction", k2_ptr->userFloat("validFraction"));
      phi_cand.addUserInt("phi_k2_ndof", k2_ptr->userInt("ndof"));
      phi_cand.addUserInt("phi_k2_numberOfValidHits", k2_ptr->userInt("nValidHits"));
      phi_cand.addUserInt("phi_k2_numberOfPixelHits", k2_ptr->userInt("numberOfPixelHits"));
      phi_cand.addUserInt("phi_k2_numberOfLostHits", k2_ptr->userInt("numberOfLostHits"));
      phi_cand.addUserInt("phi_k2_numberOfValidPixelHits", k2_ptr->userInt("numberOfValidPixelHits"));
      phi_cand.addUserInt("phi_k2_numberOfTrackerLayers", k2_ptr->userInt("numberOfTrackerLayers"));
      phi_cand.addUserInt("phi_k2_numberOfPixelLayers", k2_ptr->userInt("numberOfPixelLayers"));
      phi_cand.addUserInt("phi_k2_qualityIndex", k2_ptr->userInt("qualityIndex"));
      phi_cand.addUserInt("phi_k2_highPurityFlag", k2_ptr->userInt("highPurityFlag"));
      phi_cand.addUserFloat("phi_k2_covQopQop", k2_ptr->userFloat("covQopQop"));
      phi_cand.addUserFloat("phi_k2_covLamLam", k2_ptr->userFloat("covLamLam"));
      phi_cand.addUserFloat("phi_k2_covPhiPhi", k2_ptr->userFloat("covPhiPhi"));
      phi_cand.addUserFloat("phi_k2_covQopLam", k2_ptr->userFloat("covQopLam"));
      phi_cand.addUserFloat("phi_k2_covQopPhi", k2_ptr->userFloat("covQopPhi"));
      phi_cand.addUserFloat("phi_k2_covLamPhi", k2_ptr->userFloat("covLamPhi"));

      // apply pots-fit selection on phi candidate
      //std::cout << "before postfit selection" << std::endl;
      //if(!post_vtx_selection_phi_(phi_cand)) continue; //TODO move to after definition of user floats?
      //std::cout << "after postfit phi selection" << std::endl;

      //TODO add? fetch PV
      //float hnl_lxyz = sqrt(pow(mu0_ptr->vx() - hnl_cand.vx(), 2) + pow(mu0_ptr->vy() - hnl_cand.vy(), 2) + pow(mu0_ptr->vz() - hnl_cand.vz(), 2));
      //phi_cand.addUserFloat("hnl_l_xyz", hnl_lxyz);
      //phi_cand.addUserFloat("hnl_ct", hnl_lxyz / (hnl_cand.p4().Beta() * hnl_cand.p4().Gamma()));
      //phi_cand.addUserFloat("dilepton_vzdiff", fabs(mu0_ptr->vz()-lep_ptr->vz()));
      //phi_cand.addUserFloat("dilepton_vxdiff", fabs(mu0_ptr->vx()-lep_ptr->vx()));
      //phi_cand.addUserFloat("dilepton_vydiff", fabs(mu0_ptr->vy()-lep_ptr->vy()));
      //phi_cand.addUserFloat("dilepton_Lxy"   , sqrt(pow(mu0_ptr->vx()-lep_ptr->vx(), 2) + pow(mu0_ptr->vy()-lep_ptr->vy(), 2)));
      //phi_cand.addUserFloat("dilepton_Lxyz"  , sqrt(pow(mu0_ptr->vx()-lep_ptr->vx(), 2) + pow(mu0_ptr->vy()-lep_ptr->vy(), 2) + pow(mu0_ptr->vz()-lep_ptr->vz(), 2)));
      //// difference between the z vertex position of the pion and primary muon
      //phi_cand.addUserFloat("pion_mu0_vzdiff", fabs(mu0_ptr->vz()-pi_ptr->vz()));
       
      // computation of cos(theta*), 
      // (angle between the phi's momentum direction in the lab frame and the daughter's momentum direction in the center of mass frame)
      float mass_phi_fitted = fitter_phi.fitted_candidate().mass();
      float mass_k1 = fitter_phi.daughter_p4(0).mass(); 
      float mass_k2 = fitter_phi.daughter_p4(1).mass(); 

      float energy_phi_fitted_lab = sqrt(pow(fitter_phi.fitted_candidate().globalMomentum().x(), 2) + pow(fitter_phi.fitted_candidate().globalMomentum().y(), 2) + pow(fitter_phi.fitted_candidate().globalMomentum().z(), 2) + pow(fitter_phi.fitted_candidate().mass(), 2));
      float energy_phi_fitted_cm = mass_phi_fitted;
      float energy_k1_fitted_lab = fitter_phi.daughter_p4(0).energy();
      float energy_k1_fitted_cm = (pow(mass_k1, 2) - pow(mass_k2, 2) + pow(mass_phi_fitted, 2)) / (2. * mass_phi_fitted); 
      float energy_k2_fitted_lab = fitter_phi.daughter_p4(1).energy();
      float energy_k2_fitted_cm = (pow(mass_k2, 2) - pow(mass_k1, 2) + pow(mass_phi_fitted, 2)) / (2. * mass_phi_fitted); 

      float momentum_k1_fitted_cm = sqrt(pow(energy_k1_fitted_cm, 2) - pow(mass_k1, 2));
      float momentum_k2_fitted_cm = sqrt(pow(energy_k2_fitted_cm, 2) - pow(mass_k2, 2));

      float beta = sqrt(pow(energy_phi_fitted_lab, 2) - pow(mass_phi_fitted, 2)) / energy_phi_fitted_lab;
      float gamma = 1. / sqrt(1 - pow(beta, 2));

      float cos_theta_star_k1 = (1. / (beta * momentum_k1_fitted_cm)) * (energy_k1_fitted_lab / gamma - energy_k1_fitted_cm);
      float cos_theta_star_k2 = (1. / (beta * momentum_k2_fitted_cm)) * (energy_k2_fitted_lab / gamma - energy_k2_fitted_cm);

      phi_cand.addUserFloat("cos_theta_star_k1", cos_theta_star_k1);
      phi_cand.addUserFloat("cos_theta_star_k2", cos_theta_star_k2);

      // energy-momentum conservation (lab) 
      //phi_cand.addUserFloat("energy_diff_phi_daughters_lab", energy_phi_fitted_lab - (energy_k1_fitted_lab + energy_k2_fitted_lab)); 
      phi_cand.addUserFloat("px_diff_phi_daughters_lab", fitter_phi.fitted_candidate().globalMomentum().x() - (fitter_phi.daughter_p4(0).px() + fitter_phi.daughter_p4(1).px())); 
      phi_cand.addUserFloat("py_diff_phi_daughters_lab", fitter_phi.fitted_candidate().globalMomentum().y() - (fitter_phi.daughter_p4(0).py() + fitter_phi.daughter_p4(1).py())); 
      phi_cand.addUserFloat("pz_diff_phi_daughters_lab", fitter_phi.fitted_candidate().globalMomentum().z() - (fitter_phi.daughter_p4(0).pz() + fitter_phi.daughter_p4(1).pz())); 

      // energy-momentum conservation (lab, prefit hnl), indirect way to assess fit quality 
      phi_cand.addUserFloat("energy_diff_prefitphi_daughters_lab", phi_cand.p4().energy() - (energy_k1_fitted_lab + energy_k2_fitted_lab)); 
      phi_cand.addUserFloat("px_diff_prefitphi_daughters_lab", phi_cand.p4().px() - (fitter_phi.daughter_p4(0).px() + fitter_phi.daughter_p4(1).px())); 
      phi_cand.addUserFloat("py_diff_prefitphi_daughters_lab", phi_cand.p4().py() - (fitter_phi.daughter_p4(0).py() + fitter_phi.daughter_p4(1).py())); 
      phi_cand.addUserFloat("pz_diff_prefitphi_daughters_lab", phi_cand.p4().pz() - (fitter_phi.daughter_p4(0).pz() + fitter_phi.daughter_p4(1).pz())); 

      // energy-momentum conservation (center of mass) 
      phi_cand.addUserFloat("energy_diff_phi_daughters_cm", energy_phi_fitted_cm - (energy_k1_fitted_cm + energy_k2_fitted_cm)); 
      phi_cand.addUserFloat("p_daughters_cm", momentum_k1_fitted_cm + momentum_k2_fitted_cm); 
      
      // gen-matching
      int isMatched = 0;
      int k1_isMatched(0), k2_isMatched(0);
      int k1_genIdx(-1), k2_genIdx(-1);
      int genKaon1Mother_genPdgId(-1), genKaon2Mother_genPdgId(-1);
      //int kaon1Mother_genIdx(-1), kaon2Mother_genIdx(-1);
      //float pilep_mass_reldiff(99.), lxy_reldiff(99.);
      //float gen_lxy(-1.);
      //// species of the B mother
      //int BMother_isBu = 0;
      //int BMother_isBd = 0;
      //int BMother_isBs = 0;

      // for MC only
      if(isMC_ == true){
        // pdgId of the gen particle to which the final-state particles are matched
        int k1_genPdgId = k1_ptr->userInt("mcMatch");
        int k2_genPdgId = k2_ptr->userInt("mcMatch");
        
        // index of the gen particle to which the final-state particles are matched
        k1_genIdx = k1_ptr->userInt("mcMatchIndex"); 
        k2_genIdx = k2_ptr->userInt("mcMatchIndex"); 

        //float pilep_mass_reco = fitter.fitted_candidate().mass(); // taking the fitted mass 
        //float pilep_mass_gen = 99.;

        //float mu0_vx_gen(99.), mu0_vy_gen(99.), lep_vx_gen(99.), lep_vy_gen(99.);

        if(k1_genIdx != -1){
          // getting the associated gen particles
          edm::Ptr<reco::GenParticle> genKaon1_ptr(genParticles, k1_genIdx);

          // index of the associated mother particle
          int genKaon1Mother_genIdx = -1;
          if(genKaon1_ptr->numberOfMothers()>0) genKaon1Mother_genIdx = genKaon1_ptr->motherRef(0).key();
          //kaon1Mother_genIdx = genKaon1Mother_genIdx;

          // getting the mother particle
          edm::Ptr<reco::GenParticle> genKaon1Mother_ptr(genParticles, genKaon1Mother_genIdx);

          // getting vertices
          //mu0_vx_gen = genPrimaryMuon_ptr->vx();
          //mu0_vy_gen = genPrimaryMuon_ptr->vy();

          // pdgId of the mother particles
          genKaon1Mother_genPdgId = genKaon1Mother_ptr->pdgId();

          // matching of the kaon
          if(fabs(k1_genPdgId) == 321 && fabs(genKaon1Mother_genPdgId) == 333){
            k1_isMatched = 1;
          }
        }

        if(k2_genIdx != -1){
          // getting the associated gen particles
          edm::Ptr<reco::GenParticle> genKaon2_ptr(genParticles, k2_genIdx);

          // index of the associated mother particle
          int genKaon2Mother_genIdx = -1;
          if(genKaon2_ptr->numberOfMothers()>0) genKaon2Mother_genIdx = genKaon2_ptr->motherRef(0).key();
          //kaon2Mother_genIdx = genKaon2Mother_genIdx;

          // getting the mother particle
          edm::Ptr<reco::GenParticle> genKaon2Mother_ptr(genParticles, genKaon2Mother_genIdx);

          // getting vertices
          //mu0_vx_gen = genPrimaryMuon_ptr->vx();
          //mu0_vy_gen = genPrimaryMuon_ptr->vy();

          // pdgId of the mother particles
          genKaon2Mother_genPdgId = genKaon2Mother_ptr->pdgId();

          // matching of the kaon
          if(fabs(k2_genPdgId) == 321 && fabs(genKaon2Mother_genPdgId) == 333){
            k2_isMatched = 1;
          }
        }

        // computing displacement at gen level
        //float gen_hnl_lxy = sqrt(pow(mu0_vx_gen - lep_vx_gen, 2) + pow(mu0_vy_gen - lep_vy_gen, 2));

        // computing relative difference between gen and reco quantities
        //TODO add back?
        //pilep_mass_reldiff = fabs(pilep_mass_reco - pilep_mass_gen) / pilep_mass_gen;
        //lxy_reldiff = fabs(lxy.value() - gen_hnl_lxy) / gen_hnl_lxy;

        // matching of the full candidate
        //if(k1_isMatched==1 && k2_isMatched==1 && pi_isMatched==1 && pilep_mass_reldiff<0.1 && primaryMuonMother_genIdx==hnlMother_genIdx && lep_ptr->charge()!=pi_ptr->charge()){
        if(k1_isMatched==1 && k2_isMatched==1 && k1_ptr->charge()!=k2_ptr->charge()){
          isMatched = 1;
        }
      }

      phi_cand.addUserInt("isMatched", isMatched);
      phi_cand.addUserInt("k1_isMatched", k1_isMatched);
      phi_cand.addUserInt("k2_isMatched", k2_isMatched);
      //phi_cand.addUserFloat("pilep_mass_reco_gen_reldiff", pilep_mass_reldiff);
      //phi_cand.addUserFloat("lxy_reco_gen_reldiff", lxy_reldiff);
      //phi_cand.addUserFloat("gen_lxy", gen_lxy);

      //ret_value->push_back(phi_cand);
      vector_phi_candidates.emplace_back(phi_cand);   
    }
  }

  // sort phi candidate collection in pt
  std::sort(vector_phi_candidates.begin(), vector_phi_candidates.end(), 
             [] (auto & phi_cand1, auto & phi_cand2) -> 
                  bool {return (phi_cand1.pt() > phi_cand2.pt());} 
           );

  for (auto & phi_cand: vector_phi_candidates){
    ret_value -> emplace_back(phi_cand);
  }
  
  evt.put(std::move(ret_value));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PhiToKKBuilder);
