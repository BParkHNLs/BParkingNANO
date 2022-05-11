#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "../interface/helper.h"
#include <limits>
#include <algorithm>
#include "../interface/KinVtxFitter.h"

#include "../interface/ETHMuon.h"

template<typename Lepton>
class BToMuLPiBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<Lepton> LeptonCollection;
  typedef std::vector<pat::ETHMuon> ETHMuonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BToMuLPiBuilder(const edm::ParameterSet &cfg):
    pi_selection_          {cfg.getParameter<std::string>("pionSelection")},
    isotrk_selection_      {cfg.getParameter<std::string>("isoTracksSelection")},
    trgmu_selection_       {cfg.getParameter<std::string>("trgMuonSelection")},
    lep_selection_         {cfg.getParameter<std::string>("leptonSelection")},
    pre_vtx_selection_     {cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_    {cfg.getParameter<std::string>("postVtxSelection")},
    extra_selection_       {cfg.getParameter<std::string>("extraSelection")},
    pi_selection_dsa_      {cfg.getParameter<std::string>("pionSelection_dsa")},
    trgmu_selection_dsa_   {cfg.getParameter<std::string>("trgMuonSelection_dsa")},
    lep_selection_dsa_     {cfg.getParameter<std::string>("leptonSelection_dsa")},
    post_vtx_selection_dsa_{cfg.getParameter<std::string>("postVtxSelection_dsa")},
    lepton_type_           {cfg.getParameter<std::string>("label")},
    isMC_                  {cfg.getParameter<bool>("isMC")},

    // these two collections are ideally created beforehand by MuonTriggerSelector.cc
    //    * the former are muons that pass the preselection defined there AND match one of the 
    //      BParking triggers
    //    * the latter are all muons that pass the preselection (regardless whether they 
    //      fired the trigger). It's a superset of the previous collection
    trg_muons_         {consumes<ETHMuonCollection>                ( cfg.getParameter<edm::InputTag>("trgMuons"               ) )},
    leptons_           {consumes<LeptonCollection>                 ( cfg.getParameter<edm::InputTag>("leptons"                ) )},
    leptons_ttracks_   {consumes<TransientTrackCollection>         ( cfg.getParameter<edm::InputTag>("leptonsTransientTracks" ) )},
    pions_             {consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pions"                  ) )},
    pions_ttracks_     {consumes<TransientTrackCollection>         ( cfg.getParameter<edm::InputTag>("pionsTransientTracks"   ) )},
    isotracksToken_    {consumes<pat::PackedCandidateCollection>   ( cfg.getParameter<edm::InputTag>("tracks"                 ) )},
    isolostTracksToken_{consumes<pat::PackedCandidateCollection>   ( cfg.getParameter<edm::InputTag>("lostTracks"             ) )},
    genParticles_      {consumes<reco::GenParticleCollection>      ( cfg.getParameter<edm::InputTag>("genParticles"           ) )}, 
    beamspot_          {consumes<reco::BeamSpot>                   ( cfg.getParameter<edm::InputTag>("beamSpot"               ) )} 
    {
      produces<pat::CompositeCandidateCollection>();
    }

    // added for fetching the PV
    //vertexSrc_         { consumes<reco::VertexCollection>          ( iConfig.getParameter<edm::InputTag>( "vertexCollection"  ) )},
    //vertexSrc_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) )

  ~BToMuLPiBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  // pre-fitter preselection 
  const StringCutObjectSelector<pat::CompositeCandidate> pi_selection_; 
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_;
  const StringCutObjectSelector<pat::ETHMuon> trgmu_selection_; 
  const StringCutObjectSelector<Lepton> lep_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; 
  // post-fitter preselection 
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; 
  // extra preselection
  const StringCutObjectSelector<pat::CompositeCandidate> extra_selection_; 
  // preselection on dsa candidates
  const StringCutObjectSelector<pat::CompositeCandidate> pi_selection_dsa_; 
  const StringCutObjectSelector<pat::ETHMuon> trgmu_selection_dsa_; 
  const StringCutObjectSelector<Lepton> lep_selection_dsa_; 
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_dsa_; 

  const std::string lepton_type_;
  const bool isMC_;

  const edm::EDGetTokenT<ETHMuonCollection> trg_muons_;
  const edm::EDGetTokenT<LeptonCollection> leptons_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pions_;
  const edm::EDGetTokenT<TransientTrackCollection> pions_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;

  const edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  

  // for PV
  //const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  
};

template<typename Lepton>
void BToMuLPiBuilder<Lepton>::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<ETHMuonCollection> trg_muons;
  evt.getByToken(trg_muons_, trg_muons);

  edm::Handle<LeptonCollection> leptons;
  evt.getByToken(leptons_, leptons);
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> pions;
  evt.getByToken(pions_, pions);
  
  edm::Handle<TransientTrackCollection> pions_ttracks;
  evt.getByToken(pions_ttracks_, pions_ttracks);  

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genParticles_, genParticles);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  


  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);

  //unsigned int nTracks     = iso_tracks->size();
  //unsigned int totalTracks = nTracks + iso_lostTracks->size();

  // PV fetched for getting the trigger muon id (caveat: B is long lived)
  //edm::Handle<reco::VertexCollection> vertexHandle;
  //evt.getByToken(vertexSrc_, vertexHandle);
  //const reco::Vertex & PV = vertexHandle->front();
  
  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  for(size_t trg_mu_idx = 0; trg_mu_idx < trg_muons->size(); ++trg_mu_idx) {
    edm::Ptr<pat::ETHMuon> trg_mu_ptr(trg_muons, trg_mu_idx);

    // only select the trigger muons (are slimmedMuons by construction)
    if(lepton_type_ == "electron" && trg_mu_ptr->userInt("isTriggeringBPark") != 1) continue;

    // selection on the trigger muon
    if( !trgmu_selection_(*trg_mu_ptr) ) continue;

    math::PtEtaPhiMLorentzVector trg_mu_p4(
      trg_mu_ptr->pt(), 
      trg_mu_ptr->eta(),
      trg_mu_ptr->phi(),
      trg_mu_ptr->mass()
      );

    for(size_t pi_idx = 0; pi_idx < pions->size(); ++pi_idx) {
      edm::Ptr<pat::CompositeCandidate> pi_ptr(pions, pi_idx);

      // selection on the pion
      if( !pi_selection_(*pi_ptr) ) continue;
      
      math::PtEtaPhiMLorentzVector pi_p4(
        pi_ptr->pt(), 
        pi_ptr->eta(),
        pi_ptr->phi(),
        PI_MASS
        );
  
      // loop on selected muons and for a mu-pi candidate 
      // as well as a B candidate, that is HNL + trg mu
      for(size_t lep_idx = 0; lep_idx < leptons->size(); ++lep_idx) {
        edm::Ptr<Lepton> lep_ptr(leptons, lep_idx);
        //if(lep_ptr->isDSAMuon()) continue;
       
        // the second muon must be _other_ than the trigger muon
        if(lep_ptr->pt()==trg_mu_ptr->pt()) { // lacking of any better idea for a comparison by pointer... 
            // save anyways the position in the collection
            // trigger muons are a subset of selected muons and selected muons are those that 
            // are saved in the tress eventually (see muonsBPark_cff.py), so
            // find the position of the trigger muon in the collection of selected muons
//             std::cout << __LINE__ << "]\t selected muon pt\t"     << sel_mu_ptr->pt()  << std::endl
//                                   << "    \t trigger  muon pt\t"  << trg_mu_ptr->pt()  << std::endl
//                                   << "    \t selected muon eta\t" << sel_mu_ptr->eta() << std::endl
//                                   << "    \t trigger  muon eta\t" << trg_mu_ptr->eta() << std::endl
//                                   << "    \t selected muon phi\t" << sel_mu_ptr->phi() << std::endl
//                                   << "    \t trigger  muon phi\t" << trg_mu_ptr->phi() << std::endl
//                                   << std::endl;
            continue;
        }

        // in the muon channel, ask either muon to trigger
        if(lepton_type_ == "muon" && trg_mu_ptr->userInt("isTriggeringBPark") != 1 && lep_ptr->userInt("isTriggeringBPark") != 1) continue;

        // selection on the lepton
        if( !lep_selection_(*lep_ptr) ) continue;

        //if( lep_ptr->userInt("isDSAMuon")!=1 && !trgmu_selection_(*trg_mu_ptr) ) continue;
        //if( lep_ptr->userInt("isDSAMuon")==1 && !trgmu_selection_dsa_(*trg_mu_ptr) ) continue;

        //if( lep_ptr->userInt("isDSAMuon")!=1 && !pi_selection_(*pi_ptr) ) continue;
        //if( lep_ptr->userInt("isDSAMuon")==1 && !pi_selection_dsa_(*pi_ptr) ) continue;

        //if( lep_ptr->userInt("isDSAMuon")!=1 && !lep_selection_(*lep_ptr) ) continue;
        //if( lep_ptr->userInt("isDSAMuon")==1 && !lep_selection_dsa_(*lep_ptr) ) continue;

        math::PtEtaPhiMLorentzVector lep_p4(
          lep_ptr->pt(), 
          lep_ptr->eta(),
          lep_ptr->phi(),
          lep_ptr->mass()
          );

        // HNL candidate
        pat::CompositeCandidate hnl_cand;
        hnl_cand.setP4(lep_p4 + pi_p4);
        hnl_cand.setCharge(lep_ptr->charge() + pi_ptr->charge());
        
        hnl_cand.addUserCand("lep", lep_ptr);
        hnl_cand.addUserCand("pi", pi_ptr);

        // check if pass pre vertex cut
        if( !pre_vtx_selection_(hnl_cand) ) continue;

        // fit the mu-pi vertex
        KinVtxFitter fitter(
          {leptons_ttracks->at(lep_idx), pions_ttracks->at(pi_idx)},
          {lep_ptr->mass()             , PI_MASS                  },
          {LEP_SIGMA                   , PI_SIGMA                 } //some small sigma for the lepton mass
        );
        if(!fitter.success()) continue; // hardcoded, but do we need otherwise?
        hnl_cand.setVertex( 
          reco::Candidate::Point( 
            fitter.fitted_vtx().x(),
            fitter.fitted_vtx().y(),
            fitter.fitted_vtx().z()
          )  
        );

        auto fit_p4 = fitter.fitted_p4();
        auto lxy    = l_xy(fitter, *beamspot);

        // B candidate
        pat::CompositeCandidate b_cand;
        b_cand.setP4(hnl_cand.p4() + trg_mu_p4);
        b_cand.setCharge(hnl_cand.charge() + trg_mu_ptr->charge());

//         b_cand.addUserCand("trg_mu", trg_mu_ptr);
        // https://cmssdt.cern.ch/lxr/source/DataFormats/Candidate/interface/Candidate.h
        
//         std::cout << __LINE__ << "]\t" << hnl_cand.pt() << std::endl;
//         std::cout << __LINE__ << "]\t" << hnl_cand.originalObjectRef().isNull() << std::endl;
//         edm::Ptr<pat::CompositeCandidate> hnl_cand_ptr = edm::refToPtr(hnl_cand.originalObjectRef());
//         b_cand.addUserCand("hnl", hnl_cand.originalObjectRef());
//         b_cand.addUserCand("hnl", hnl_cand.sourceCandidatePtr(0));

        b_cand.addDaughter(*trg_mu_ptr, "trg_mu");
        b_cand.addDaughter( hnl_cand  , "hnl"   );

        b_cand.addUserInt  ("hnl_vtx_OK"             , fitter.success()                                                        );
        b_cand.addUserFloat("hnl_vtx_chi2"           , fitter.chi2()                                                           );
        b_cand.addUserFloat("hnl_vtx_ndof"           , fitter.dof()                                                            ); // float??
        b_cand.addUserFloat("hnl_vtx_prob"           , fitter.prob()                                                           );
        b_cand.addUserFloat("hnl_fitted_pt"          , fit_p4.pt()                                                             ); 
        b_cand.addUserFloat("hnl_fitted_eta"         , fit_p4.eta()                                                            );
        b_cand.addUserFloat("hnl_fitted_phi"         , fit_p4.phi()                                                            );
        b_cand.addUserFloat("hnl_fitted_mass"        , fitter.fitted_candidate().mass()                                        );      
        b_cand.addUserFloat("hnl_fitted_massErr"     , sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
        b_cand.addUserFloat("hnl_cos_theta_2D"       , cos_theta_2D(fitter, *beamspot, hnl_cand.p4())                          );
        b_cand.addUserFloat("hnl_fitted_cos_theta_2D", cos_theta_2D(fitter, *beamspot, fit_p4)                                 );
        b_cand.addUserFloat("hnl_l_xy"               , lxy.value()                                                             );
        b_cand.addUserFloat("hnl_l_xy_unc"           , lxy.error()                                                             );
        b_cand.addUserFloat("hnl_ls_xy"              , lxy.value()/lxy.error()                                                 );
        b_cand.addUserFloat("hnl_charge"             , hnl_cand.charge()                                                       );
        b_cand.addUserFloat("hnl_vtx_x"              , hnl_cand.vx()                                                           );
        b_cand.addUserFloat("hnl_vtx_y"              , hnl_cand.vy()                                                           );
        b_cand.addUserFloat("hnl_vtx_z"              , hnl_cand.vz()                                                           );
        b_cand.addUserFloat("hnl_vtx_ex"             , sqrt(fitter.fitted_vtx_uncertainty().cxx())                             );
        b_cand.addUserFloat("hnl_vtx_ey"             , sqrt(fitter.fitted_vtx_uncertainty().cyy())                             );
        b_cand.addUserFloat("hnl_vtx_ez"             , sqrt(fitter.fitted_vtx_uncertainty().czz())                             );
        b_cand.addUserFloat("hnl_fitted_lep_pt"      , fitter.daughter_p4(0).pt()                                              ); 
        b_cand.addUserFloat("hnl_fitted_lep_eta"     , fitter.daughter_p4(0).eta()                                             );
        b_cand.addUserFloat("hnl_fitted_lep_phi"     , fitter.daughter_p4(0).phi()                                             );
        b_cand.addUserFloat("hnl_fitted_lep_mass"    , fitter.daughter_p4(0).mass()                                            );
        b_cand.addUserFloat("hnl_fitted_pi_pt"       , fitter.daughter_p4(1).pt()                                              ); 
        b_cand.addUserFloat("hnl_fitted_pi_eta"      , fitter.daughter_p4(1).eta()                                             );
        b_cand.addUserFloat("hnl_fitted_pi_phi"      , fitter.daughter_p4(1).phi()                                             );
        b_cand.addUserFloat("hnl_fitted_pi_mass"     , fitter.daughter_p4(1).mass()                                            );

        // post fit selection
        if( !post_vtx_selection_(b_cand) ) continue;        

        // computation of cos(theta*), 
        // (angle between the hnl's momentum direction in the lab frame and the daughter's momentum direction in the center of mass frame)
        float mass_hnl_fitted = fitter.fitted_candidate().mass();
        float mass_lep = fitter.daughter_p4(0).mass(); 
        float mass_pi = fitter.daughter_p4(1).mass(); 

        float energy_hnl_fitted_lab = sqrt(pow(fitter.fitted_candidate().globalMomentum().x(), 2) + pow(fitter.fitted_candidate().globalMomentum().y(), 2) + pow(fitter.fitted_candidate().globalMomentum().z(), 2) + pow(fitter.fitted_candidate().mass(), 2));
        float energy_hnl_fitted_cm = mass_hnl_fitted;
        float energy_pion_fitted_lab = fitter.daughter_p4(1).energy();
        float energy_pion_fitted_cm = (pow(mass_pi, 2) - pow(mass_lep, 2) + pow(mass_hnl_fitted, 2)) / (2. * mass_hnl_fitted); 
        float energy_lepton_fitted_lab = fitter.daughter_p4(0).energy();
        float energy_lepton_fitted_cm = (pow(mass_lep, 2) - pow(mass_pi, 2) + pow(mass_hnl_fitted, 2)) / (2. * mass_hnl_fitted); 

        float momentum_pion_fitted_cm = sqrt(pow(energy_pion_fitted_cm, 2) - pow(mass_pi, 2));
        float momentum_lepton_fitted_cm = sqrt(pow(energy_lepton_fitted_cm, 2) - pow(mass_lep, 2));

        float beta = sqrt(pow(energy_hnl_fitted_lab, 2) - pow(mass_hnl_fitted, 2)) / energy_hnl_fitted_lab;
        float gamma = 1. / sqrt(1 - pow(beta, 2));

        float cos_theta_star_pion = (1. / (beta * momentum_pion_fitted_cm)) * (energy_pion_fitted_lab / gamma - energy_pion_fitted_cm);
        float cos_theta_star_lepton = (1. / (beta * momentum_lepton_fitted_cm)) * (energy_lepton_fitted_lab / gamma - energy_lepton_fitted_cm);

        b_cand.addUserFloat("cos_theta_star_pion", cos_theta_star_pion);
        b_cand.addUserFloat("cos_theta_star_lepton", cos_theta_star_lepton);

        // energy-momentum conservation (lab) 
        b_cand.addUserFloat("energy_diff_hnl_daughters_lab", energy_hnl_fitted_lab - (energy_lepton_fitted_lab + energy_pion_fitted_lab)); 
        b_cand.addUserFloat("px_diff_hnl_daughters_lab", fitter.fitted_candidate().globalMomentum().x() - (fitter.daughter_p4(0).px() + fitter.daughter_p4(1).px())); 
        b_cand.addUserFloat("py_diff_hnl_daughters_lab", fitter.fitted_candidate().globalMomentum().y() - (fitter.daughter_p4(0).py() + fitter.daughter_p4(1).py())); 
        b_cand.addUserFloat("pz_diff_hnl_daughters_lab", fitter.fitted_candidate().globalMomentum().z() - (fitter.daughter_p4(0).pz() + fitter.daughter_p4(1).pz())); 

        // energy-momentum conservation (lab, prefit hnl), indirect way to assess fit quality 
        b_cand.addUserFloat("energy_diff_prefithnl_daughters_lab", hnl_cand.p4().energy() - (energy_lepton_fitted_lab + energy_pion_fitted_lab)); 
        b_cand.addUserFloat("px_diff_prefithnl_daughters_lab", hnl_cand.p4().px() - (fitter.daughter_p4(0).px() + fitter.daughter_p4(1).px())); 
        b_cand.addUserFloat("py_diff_prefithnl_daughters_lab", hnl_cand.p4().py() - (fitter.daughter_p4(0).py() + fitter.daughter_p4(1).py())); 
        b_cand.addUserFloat("pz_diff_prefithnl_daughters_lab", hnl_cand.p4().pz() - (fitter.daughter_p4(0).pz() + fitter.daughter_p4(1).pz())); 

        // energy-momentum conservation (center of mass) 
        b_cand.addUserFloat("energy_diff_hnl_daughters_cm", energy_hnl_fitted_cm - (energy_lepton_fitted_cm + energy_pion_fitted_cm)); 
        b_cand.addUserFloat("p_daughters_cm", momentum_lepton_fitted_cm + momentum_pion_fitted_cm); 

        // displacement
        float hnl_lxyz = sqrt(pow(trg_mu_ptr->vx() - hnl_cand.vx(), 2) + pow(trg_mu_ptr->vy() - hnl_cand.vy(), 2) + pow(trg_mu_ptr->vz() - hnl_cand.vz(), 2));
        b_cand.addUserFloat("hnl_l_xyz", hnl_lxyz);
        b_cand.addUserFloat("hnl_ct", hnl_lxyz / (hnl_cand.p4().Beta() * hnl_cand.p4().Gamma()));
      
        // adding trigger muon information to the b candidate
        b_cand.addUserFloat("trg_muon_pt" , trg_mu_ptr->pt());
        b_cand.addUserFloat("trg_muon_eta", trg_mu_ptr->eta());
        b_cand.addUserFloat("trg_muon_phi", trg_mu_ptr->phi());

        // difference between the z vertex position of the selected muon and tigger muon
        // computed at the prefit stage 
        b_cand.addUserFloat("dilepton_vzdiff", fabs(trg_mu_ptr->vz()-lep_ptr->vz()));
        b_cand.addUserFloat("dilepton_vxdiff", fabs(trg_mu_ptr->vx()-lep_ptr->vx()));
        b_cand.addUserFloat("dilepton_vydiff", fabs(trg_mu_ptr->vy()-lep_ptr->vy()));
        b_cand.addUserFloat("dilepton_Lxy"   , sqrt(pow(trg_mu_ptr->vx()-lep_ptr->vx(), 2) + pow(trg_mu_ptr->vy()-lep_ptr->vy(), 2)));
        b_cand.addUserFloat("dilepton_Lxyz"  , sqrt(pow(trg_mu_ptr->vx()-lep_ptr->vx(), 2) + pow(trg_mu_ptr->vy()-lep_ptr->vy(), 2) + pow(trg_mu_ptr->vz()-lep_ptr->vz(), 2)));
        
        // difference between the z vertex position of the pion and tigger muon
        b_cand.addUserFloat("pion_trgmuon_vzdiff", fabs(trg_mu_ptr->vz()-pi_ptr->vz()));

        // fetch the id of the sel muon at the secondary vertex (use instead info saved in the muonsBPark collection?)
        //if(lepton_type_ == "muon"){
        //  float sel_muon_isSoft   = lep_ptr->isSoftMuon  ((const reco::Vertex&) fitter) ? 1. : 0. ;
        //  float sel_muon_isTight  = lep_ptr->isTightMuon ((const reco::Vertex&) fitter) ? 1. : 0. ;
        //  float sel_muon_isMedium = lep_ptr->isMediumMuon()                             ? 1. : 0. ;
        //  float sel_muon_isLoose  = lep_ptr->isLooseMuon ()                             ? 1. : 0. ;

        //  b_cand.addUserFloat("sel_muon_isSoft"       , sel_muon_isSoft      );
        //  b_cand.addUserFloat("sel_muon_isTight"      , sel_muon_isTight     );
        //  b_cand.addUserFloat("sel_muon_isMedium"     , sel_muon_isMedium    );
        //  b_cand.addUserFloat("sel_muon_isLoose"      , sel_muon_isLoose     );
        //}

        // adding dR quantities (with fitted quantities)
        float dR_lep_pi    = reco::deltaR(fitter.daughter_p4(0), fitter.daughter_p4(1)); 
        float dR_trgmu_hnl = reco::deltaR((*trg_mu_ptr), hnl_cand                     ); 
        float dR_trgmu_lep = reco::deltaR((*trg_mu_ptr), fitter.daughter_p4(0)        ); 
        float dR_trgmu_pi  = reco::deltaR((*trg_mu_ptr), fitter.daughter_p4(1)        ); 
        b_cand.addUserFloat("dr_lep_pi"      , dR_lep_pi    );
        b_cand.addUserFloat("dr_trgmu_hnl"   , dR_trgmu_hnl );
        b_cand.addUserFloat("dr_trgmu_lep"   , dR_trgmu_lep  );
        b_cand.addUserFloat("dr_trgmu_pi"    , dR_trgmu_pi  );

        float dPhi_lep_pi    = reco::deltaPhi(fitter.daughter_p4(0).phi(), fitter.daughter_p4(1).phi()); 
        float dPhi_trgmu_hnl = reco::deltaPhi(trg_mu_ptr->phi(), fit_p4.phi()); 
        float dPhi_trgmu_lep = reco::deltaPhi(trg_mu_ptr->phi(), fitter.daughter_p4(0).phi()); 
        float dPhi_trgmu_pi  = reco::deltaPhi(trg_mu_ptr->phi(), fitter.daughter_p4(1).phi()); 
        b_cand.addUserFloat("dphi_lep_pi"      , dPhi_lep_pi    );
        b_cand.addUserFloat("dphi_trgmu_hnl"   , dPhi_trgmu_hnl );
        b_cand.addUserFloat("dphi_trgmu_lep"   , dPhi_trgmu_lep );
        b_cand.addUserFloat("dphi_trgmu_pi"    , dPhi_trgmu_pi  );

        float dEta_lep_pi    = fitter.daughter_p4(0).eta() - fitter.daughter_p4(1).eta(); 
        float dEta_trgmu_hnl = trg_mu_ptr->eta() - fit_p4.eta(); 
        float dEta_trgmu_lep = trg_mu_ptr->eta() - fitter.daughter_p4(0).eta();
        float dEta_trgmu_pi  = trg_mu_ptr->eta() - fitter.daughter_p4(1).eta();
        b_cand.addUserFloat("deta_lep_pi"      , dEta_lep_pi    );
        b_cand.addUserFloat("deta_trgmu_hnl"   , dEta_trgmu_hnl );
        b_cand.addUserFloat("deta_trgmu_lep"   , dEta_trgmu_lep );
        b_cand.addUserFloat("deta_trgmu_pi"    , dEta_trgmu_pi  );

        // difference of the kinematics of the objects and their fitted value
        float dE_pi_fit_pi   = pi_p4.energy() - fitter.daughter_p4(1).energy();
        float dPt_pi_fit_pi  = pi_ptr->pt() - fitter.daughter_p4(1).pt(); 
        float dPx_pi_fit_pi  = pi_p4.px() - fitter.daughter_p4(1).px();
        float dPy_pi_fit_pi  = pi_p4.py() - fitter.daughter_p4(1).py();
        float dPz_pi_fit_pi  = pi_p4.pz() - fitter.daughter_p4(1).pz();
        float dEta_pi_fit_pi = pi_ptr->eta() - fitter.daughter_p4(1).eta(); 
        float dPhi_pi_fit_pi = reco::deltaPhi(pi_ptr->phi(), fitter.daughter_p4(1).phi()); 
        b_cand.addUserFloat("de_pi_fit_pi"  , dE_pi_fit_pi);
        b_cand.addUserFloat("dpt_pi_fit_pi" , dPt_pi_fit_pi);
        b_cand.addUserFloat("dpx_pi_fit_pi" , dPx_pi_fit_pi);
        b_cand.addUserFloat("dpy_pi_fit_pi" , dPy_pi_fit_pi);
        b_cand.addUserFloat("dpz_pi_fit_pi" , dPz_pi_fit_pi);
        b_cand.addUserFloat("deta_pi_fit_pi", dEta_pi_fit_pi);
        b_cand.addUserFloat("dphi_pi_fit_pi", dPhi_pi_fit_pi);

        float dE_lep_fit_lep   = lep_p4.energy() - fitter.daughter_p4(0).energy();
        float dPt_lep_fit_lep  = lep_ptr->pt() - fitter.daughter_p4(0).pt(); 
        float dPx_lep_fit_lep  = lep_p4.px() - fitter.daughter_p4(0).px();
        float dPy_lep_fit_lep  = lep_p4.py() - fitter.daughter_p4(0).py();
        float dPz_lep_fit_lep  = lep_p4.pz() - fitter.daughter_p4(0).pz();
        float dEta_lep_fit_lep = lep_ptr->eta() - fitter.daughter_p4(0).eta(); 
        float dPhi_lep_fit_lep = reco::deltaPhi(lep_ptr->phi(), fitter.daughter_p4(0).phi()); 
        b_cand.addUserFloat("de_lep_fit_lep"   , dE_lep_fit_lep);
        b_cand.addUserFloat("dpt_lep_fit_lep"  , dPt_lep_fit_lep);
        b_cand.addUserFloat("dpx_lep_fit_lep"  , dPx_lep_fit_lep);
        b_cand.addUserFloat("dpy_lep_fit_lep"  , dPy_lep_fit_lep);
        b_cand.addUserFloat("dpz_lep_fit_lep"  , dPz_lep_fit_lep);
        b_cand.addUserFloat("deta_lep_fit_lep" , dEta_lep_fit_lep);
        b_cand.addUserFloat("dphi_lep_fit_lep" , dPhi_lep_fit_lep);

        float dE_hnl_fit_hnl  = hnl_cand.p4().energy() - energy_hnl_fitted_lab;
        float dPt_hnl_fit_hnl = hnl_cand.p4().pt() - fit_p4.pt();
        float dPx_hnl_fit_hnl = hnl_cand.p4().px() - fitter.fitted_candidate().globalMomentum().x();
        float dPy_hnl_fit_hnl = hnl_cand.p4().py() - fitter.fitted_candidate().globalMomentum().y();
        float dPz_hnl_fit_hnl = hnl_cand.p4().pz() - fitter.fitted_candidate().globalMomentum().z();
        float dEta_hnl_fit_hnl = hnl_cand.p4().eta() - fit_p4.eta(); 
        float dPhi_hnl_fit_hnl = reco::deltaPhi(hnl_cand.p4().phi(), fit_p4.phi()); 
        b_cand.addUserFloat("de_hnl_fit_hnl"  , dE_hnl_fit_hnl);
        b_cand.addUserFloat("dpt_hnl_fit_hnl" , dPt_hnl_fit_hnl);
        b_cand.addUserFloat("dpx_hnl_fit_hnl" , dPx_hnl_fit_hnl);
        b_cand.addUserFloat("dpy_hnl_fit_hnl" , dPy_hnl_fit_hnl);
        b_cand.addUserFloat("dpz_hnl_fit_hnl" , dPz_hnl_fit_hnl);
        b_cand.addUserFloat("dphi_hnl_fit_hnl", dPhi_hnl_fit_hnl);
        b_cand.addUserFloat("deta_hnl_fit_hnl", dEta_hnl_fit_hnl);

        // save information on pions 
        b_cand.addUserFloat("pion_pt"                  , pi_ptr->pt()                             );
        b_cand.addUserFloat("pion_eta"                 , pi_ptr->eta()                            );
        b_cand.addUserFloat("pion_phi"                 , pi_ptr->phi()                            );
        b_cand.addUserFloat("pion_mass"                , pi_ptr->mass()                           );
        b_cand.addUserInt("pion_charge"                , pi_ptr->charge()                         );
        b_cand.addUserInt("pion_pdgId"                 , pi_ptr->pdgId()                          );
        b_cand.addUserFloat("pion_vx"                  , pi_ptr->vx()                             );
        b_cand.addUserFloat("pion_vy"                  , pi_ptr->vy()                             );
        b_cand.addUserFloat("pion_vz"                  , pi_ptr->vz()                             );
        b_cand.addUserFloat("pion_dz"                  , pi_ptr->userFloat("dz")                  );
        b_cand.addUserFloat("pion_dxy"                 , pi_ptr->userFloat("dxy")                 );
        b_cand.addUserFloat("pion_dzS"                 , pi_ptr->userFloat("dzS")                 );
        b_cand.addUserFloat("pion_dxyS"                , pi_ptr->userFloat("dxyS")                );
        b_cand.addUserFloat("pion_DCASig"              , pi_ptr->userFloat("DCASig")              );
        b_cand.addUserInt("pion_ispacked"              , pi_ptr->userInt("isPacked")              );
        b_cand.addUserInt("pion_islost"                , pi_ptr->userInt("isLostTrk")             );
        b_cand.addUserFloat("pion_chi2"                , pi_ptr->userFloat("chi2")                );
        b_cand.addUserFloat("pion_normalisedChi2"      , pi_ptr->userFloat("normalisedChi2")      );
        b_cand.addUserFloat("pion_validFraction"       , pi_ptr->userFloat("validFraction")       );
        b_cand.addUserInt("pion_ndof"                  , pi_ptr->userInt("ndof")                  );
        b_cand.addUserInt("pion_numberOfValidHits"     , pi_ptr->userInt("nValidHits")            );
        b_cand.addUserInt("pion_numberOfLostHits"      , pi_ptr->userInt("numberOfLostHits")      );
        b_cand.addUserInt("pion_numberOfValidPixelHits", pi_ptr->userInt("numberOfValidPixelHits"));
        b_cand.addUserInt("pion_numberOfTrackerLayers" , pi_ptr->userInt("numberOfTrackerLayers") );
        b_cand.addUserInt("pion_numberOfPixelLayers"   , pi_ptr->userInt("numberOfPixelLayers")   );
        b_cand.addUserInt("pion_qualityIndex"          , pi_ptr->userInt("qualityIndex")          );
        b_cand.addUserInt("pion_highPurityFlag"        , pi_ptr->userInt("highPurityFlag")        );


        // post fit selection
        //if( !post_vtx_selection_(b_cand) ) continue;        

        //if( lep_ptr->userInt("isDSAMuon")!=1 && !post_vtx_selection_(b_cand) ) continue;
        //if( lep_ptr->userInt("isDSAMuon")==1 && !post_vtx_selection_dsa_(b_cand) ) continue;

        // isolation
        float trg_mu_iso03 = 0; 
        float trg_mu_iso04 = 0;
        float lep_iso03 = 0; 
        float lep_iso04 = 0;
        float pi_iso03  = 0; 
        float pi_iso04  = 0;
        float hnl_iso03 = 0;
        float hnl_iso04 = 0;

        // with conditions: best track + close to B 
        float lep_iso03_close = 0;
        float lep_iso04_close = 0;
        float trg_mu_iso03_close = 0; 
        float trg_mu_iso04_close = 0;
        float pi_iso03_close  = 0; 
        float pi_iso04_close  = 0;
        float hnl_iso03_close = 0;
        float hnl_iso04_close = 0;
        float lep_iso03_rel_close = 0;
        float lep_iso04_rel_close = 0;
        float trg_mu_iso03_rel_close = 0; 
        float trg_mu_iso04_rel_close = 0;
        float pi_iso03_rel_close  = 0; 
        float pi_iso04_rel_close  = 0;
        float hnl_iso03_rel_close = 0;
        float hnl_iso04_rel_close = 0;
        
        // nTracks     = iso_tracks->size();
        // totalTracks = nTracks + iso_lostTracks->size();

        /*
        for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
        
          const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];

          // same preselection as for pion tracks
          if( !isotrk_selection_(trk) ) continue;

          // check if the track is the pion
          if (pi_ptr->userCand("b_cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
         
          // check if the track is one of the two leptons 
          if (track_to_lepton_match(trg_mu_ptr, iso_tracks.id(), iTrk) || 
              track_to_lepton_match(lep_ptr, iso_tracks.id(), iTrk) ) continue;

          // add to final particle iso if dR < cone
          float dr_to_trgmu = deltaR(b_cand.userFloat("trg_muon_eta")     , b_cand.userFloat("trg_muon_phi")     , trk.eta(), trk.phi());
          float dr_to_lep = deltaR(b_cand.userFloat("hnl_fitted_lep_eta"), b_cand.userFloat("hnl_fitted_lep_phi"), trk.eta(), trk.phi());
          float dr_to_pi    = deltaR(b_cand.userFloat("hnl_fitted_pi_eta"), b_cand.userFloat("hnl_fitted_pi_phi"), trk.eta(), trk.phi());
          float dr_to_hnl   = deltaR(b_cand.userFloat("hnl_fitted_eta")   , b_cand.userFloat("hnl_fitted_phi")   , trk.eta(), trk.phi());

          if (dr_to_trgmu < 0.4){
            trg_mu_iso04 += trk.pt();
            if ( dr_to_trgmu < 0.3) trg_mu_iso03 += trk.pt();
          }
          if (dr_to_lep < 0.4){
            lep_iso04 += trk.pt();
            if (dr_to_lep < 0.3)  lep_iso03 += trk.pt();
          }
          if (dr_to_pi < 0.4){
            pi_iso04 += trk.pt();
            if (dr_to_pi < 0.3) pi_iso03 += trk.pt();
          }
          if (dr_to_hnl < 0.4){
            hnl_iso04 += trk.pt();
            if (dr_to_hnl < 0.3) hnl_iso03 += trk.pt();
          }

          // add requirement of the tracks to be close to the B
          if (!trg_mu_ptr->bestTrack() || fabs(trk.dz() - trg_mu_ptr->bestTrack()->dz()) > 0.4) continue;
          if (!lep_ptr->bestTrack() || fabs(trk.dz() - lep_ptr->bestTrack()->dz()) > 0.4) continue;
          //if (!pi_ptr->bestTrack() || fabs(trk.dz() - pi_ptr->bestTrack()->dz()) > 0.4) continue; // pion never passes bestTrack requirement
          if (fabs(trk.dz() - pi_ptr->userFloat("dz")) > 0.4) continue; //dropping requirement of best track

          if (dr_to_trgmu < 0.4){
            trg_mu_iso04_close += trk.pt();
            if ( dr_to_trgmu < 0.3) trg_mu_iso03_close += trk.pt();
          }
          if (dr_to_lep < 0.4){
            lep_iso04_close += trk.pt();
            if (dr_to_lep < 0.3)  lep_iso03_close += trk.pt();
          }
          if (dr_to_pi < 0.4){
            pi_iso04_close += trk.pt();
            if (dr_to_pi < 0.3) pi_iso03_close += trk.pt();
          }
          if (dr_to_hnl < 0.4){
            hnl_iso04_close += trk.pt();
            if (dr_to_hnl < 0.3) hnl_iso03_close += trk.pt();
          }
        }
        */

        trg_mu_iso03_rel_close = trg_mu_iso03_close / trg_mu_ptr->pt();
        trg_mu_iso04_rel_close = trg_mu_iso04_close / trg_mu_ptr->pt();
        lep_iso03_rel_close = lep_iso03_close / lep_ptr->pt();
        lep_iso04_rel_close = lep_iso04_close / lep_ptr->pt();
        pi_iso03_rel_close = pi_iso03_close / pi_ptr->pt();
        pi_iso04_rel_close = pi_iso04_close / pi_ptr->pt();
        hnl_iso03_rel_close = hnl_iso03_close / fit_p4.pt();
        hnl_iso04_rel_close = hnl_iso04_close / fit_p4.pt();

        b_cand.addUserFloat("trg_mu_iso03", trg_mu_iso03);
        b_cand.addUserFloat("trg_mu_iso04", trg_mu_iso04);
        b_cand.addUserFloat("lep_iso03", lep_iso03);
        b_cand.addUserFloat("lep_iso04", lep_iso04);
        b_cand.addUserFloat("pi_iso03" , pi_iso03 );
        b_cand.addUserFloat("pi_iso04" , pi_iso04 );
        b_cand.addUserFloat("hnl_iso03" , hnl_iso03 );
        b_cand.addUserFloat("hnl_iso04" , hnl_iso04 );

        // add requirement of the tracks to be close to the B
        b_cand.addUserFloat("trg_mu_iso03_close", trg_mu_iso03_close);
        b_cand.addUserFloat("trg_mu_iso04_close", trg_mu_iso04_close);
        b_cand.addUserFloat("sel_mu_iso03_close", lep_iso03_close);
        b_cand.addUserFloat("sel_mu_iso04_close", lep_iso04_close);
        b_cand.addUserFloat("pi_iso03_close", pi_iso03_close);
        b_cand.addUserFloat("pi_iso04_close", pi_iso04_close);
        b_cand.addUserFloat("hnl_iso03_close", hnl_iso03_close);
        b_cand.addUserFloat("hnl_iso04_close", hnl_iso04_close);
        b_cand.addUserFloat("trg_mu_iso03_rel_close", trg_mu_iso03_rel_close);
        b_cand.addUserFloat("trg_mu_iso04_rel_close", trg_mu_iso04_rel_close);
        b_cand.addUserFloat("sel_mu_iso03_rel_close", lep_iso03_rel_close);
        b_cand.addUserFloat("sel_mu_iso04_rel_close", lep_iso04_rel_close);
        b_cand.addUserFloat("pi_iso03_rel_close", pi_iso03_rel_close);
        b_cand.addUserFloat("pi_iso04_rel_close", pi_iso04_rel_close);
        b_cand.addUserFloat("hnl_iso03_rel_close", hnl_iso03_rel_close);
        b_cand.addUserFloat("hnl_iso04_rel_close", hnl_iso04_rel_close);

        // position of the muons / tracks in their own collections
        b_cand.addUserInt("trg_mu_idx", trg_mu_idx);
        b_cand.addUserInt("lep_idx", lep_idx);
        b_cand.addUserInt("pi_idx", pi_idx);

        // invariant masses
        float dilepton_mass = (fitter.daughter_p4(0) + trg_mu_p4).mass();
        float trgmu_pi_mass = (fitter.daughter_p4(1) + trg_mu_p4).mass();
        b_cand.addUserFloat("dilepton_mass", dilepton_mass);
        b_cand.addUserFloat("trgmu_pi_mass", trgmu_pi_mass);

        float dilepton_pt = (fitter.daughter_p4(0) + trg_mu_p4).pt();
        float trgmu_pi_pt = (fitter.daughter_p4(1) + trg_mu_p4).pt();
        b_cand.addUserFloat("dilepton_pt", dilepton_pt);
        b_cand.addUserFloat("trgmu_pi_pt", trgmu_pi_pt);

        // extra selection
        if( !extra_selection_(b_cand) ) continue;        

        // gen-matching
        int isMatched = 0;
        int trg_mu_isMatched(0), sel_mu_isMatched(0), pi_isMatched(0);
        int sel_mu_isFake(0), pi_isFake(0);
        int trg_mu_genIdx(-1), sel_mu_genIdx(-1), pi_genIdx(-1);
        int genTriggerMuonMother_genPdgId(-1), genMuonMother_genPdgId(-1), genPionMother_genPdgId(-1);
        int triggerMuonMother_genIdx(-1), hnlMother_genIdx(-1);
        float mupi_mass_reldiff(99.), lxy_reldiff(99.);
        float gen_lxy(-1.);

        // for MC only
        if(isMC_ == true){
          // pdgId of the gen particle to which the final-state particles are matched
          int trg_mu_genPdgId = trg_mu_ptr->userInt("mcMatch");
          int sel_mu_genPdgId = lep_ptr->userInt("mcMatch");
          int pi_genPdgId     = pi_ptr->userInt("mcMatch");
          
          // index of the gen particle to which the final-state particles are matched
          trg_mu_genIdx   = trg_mu_ptr->userInt("mcMatchIndex"); 
          sel_mu_genIdx   = lep_ptr->userInt("mcMatchIndex"); 
          pi_genIdx       = pi_ptr->userInt("mcMatchIndex"); 

          float mupi_mass_reco = fitter.fitted_candidate().mass(); // taking the fitted mass 
          float mupi_mass_gen = 99.;

          float trgmu_vx_gen(99.), trgmu_vy_gen(99.), mu_vx_gen(99.), mu_vy_gen(99.);

          if(trg_mu_genIdx != -1){
            // getting the associated gen particles
            edm::Ptr<reco::GenParticle> genTriggerMuon_ptr(genParticles, trg_mu_genIdx);

            // index of the associated mother particle
            int genTriggerMuonMother_genIdx = -1;
            if(genTriggerMuon_ptr->numberOfMothers()>0) genTriggerMuonMother_genIdx = genTriggerMuon_ptr->motherRef(0).key();
            triggerMuonMother_genIdx = genTriggerMuonMother_genIdx;

            // getting the mother particle
            edm::Ptr<reco::GenParticle> genTriggerMuonMother_ptr(genParticles, genTriggerMuonMother_genIdx);

            // getting vertices
            trgmu_vx_gen = genTriggerMuon_ptr->vx();
            trgmu_vy_gen = genTriggerMuon_ptr->vy();

            // pdgId of the mother particles
            genTriggerMuonMother_genPdgId = genTriggerMuonMother_ptr->pdgId();

            // matching of the trigger muon
            if(fabs(trg_mu_genPdgId) == 13 && (fabs(genTriggerMuonMother_genPdgId) == 511 || fabs(genTriggerMuonMother_genPdgId) == 521 
                  || fabs(genTriggerMuonMother_genPdgId) == 531 || fabs(genTriggerMuonMother_genPdgId) == 541)){
              trg_mu_isMatched = 1;
            }
          }

          if(sel_mu_genIdx != -1){
            // getting the associated gen particles
            edm::Ptr<reco::GenParticle> genMuon_ptr(genParticles, sel_mu_genIdx);

            // index of the associated mother particle
            int genMuonMother_genIdx = -1;
            if(genMuon_ptr->numberOfMothers()>0) genMuonMother_genIdx = genMuon_ptr->motherRef(0).key();

            // getting the mother particle
            edm::Ptr<reco::GenParticle> genMuonMother_ptr(genParticles, genMuonMother_genIdx);

            // index of the grand-mother particle
            if(genMuonMother_ptr->numberOfMothers()>0) hnlMother_genIdx = genMuonMother_ptr->motherRef(0).key();

            // fetching mass
            mupi_mass_gen = genMuonMother_ptr->mass();

            // getting vertices
            mu_vx_gen = genMuon_ptr->vx();
            mu_vy_gen = genMuon_ptr->vy();

            // pdgId of the mother particles
            genMuonMother_genPdgId = genMuonMother_ptr->pdgId();

            // matching of the displaced lepton
            if(fabs(sel_mu_genPdgId) == 13 && fabs(genMuonMother_genPdgId) == 9900015){
              sel_mu_isMatched = 1;
            }
            else if(fabs(sel_mu_genPdgId) == 211 && fabs(genMuonMother_genPdgId) == 9900015){
              sel_mu_isFake = 1;
              sel_mu_isMatched = 1;
            }
          }

          if(pi_genIdx != -1){
            // getting the associated gen particles
            edm::Ptr<reco::GenParticle> genPion_ptr(genParticles, pi_genIdx);

            // index of the associated mother particle
            int genPionMother_genIdx = -1;
            if(genPion_ptr->numberOfMothers()>0) genPionMother_genIdx = genPion_ptr->motherRef(0).key();

            // getting the mother particles
            edm::Ptr<reco::GenParticle> genPionMother_ptr(genParticles, genPionMother_genIdx);

            // pdgId of the mother particle
            genPionMother_genPdgId = genPionMother_ptr->pdgId();

            // matching of the displaced pion
            if(fabs(pi_genPdgId) == 211 && fabs(genPionMother_genPdgId) == 9900015){
              pi_isMatched = 1;
            }
            else if(fabs(pi_genPdgId) == 13 && fabs(genPionMother_genPdgId) == 9900015){
              pi_isFake = 1;
              pi_isMatched = 1;
            }
          }

          // computing displacement at gen level
          float gen_hnl_lxy = sqrt(pow(trgmu_vx_gen - mu_vx_gen, 2) + pow(trgmu_vy_gen - mu_vy_gen, 2));

          // computing relative difference between gen and reco quantities
          mupi_mass_reldiff = fabs(mupi_mass_reco - mupi_mass_gen) / mupi_mass_gen;
          lxy_reldiff = fabs(lxy.value() - gen_hnl_lxy) / gen_hnl_lxy;
        
          // matching of the full mulpi candidate
          if(fabs(pi_genPdgId)!=fabs(sel_mu_genPdgId) && trg_mu_isMatched==1 && sel_mu_isMatched==1 && pi_isMatched==1 && mupi_mass_reldiff<0.1 && lxy_reldiff<1 && triggerMuonMother_genIdx==hnlMother_genIdx && lep_ptr->charge()!=pi_ptr->charge()){
            isMatched = 1;
          }
        }

        b_cand.addUserInt("isMatched", isMatched);
        b_cand.addUserInt("trg_mu_isMatched", trg_mu_isMatched);
        b_cand.addUserInt("sel_mu_isMatched", sel_mu_isMatched);
        b_cand.addUserInt("pi_isMatched", pi_isMatched);
        b_cand.addUserInt("sel_mu_isFake", sel_mu_isFake);
        b_cand.addUserInt("pi_isFake", pi_isFake);
        b_cand.addUserInt("matching_trg_mu_genIdx", trg_mu_genIdx);
        b_cand.addUserInt("matching_sel_mu_genIdx", sel_mu_genIdx);
        b_cand.addUserInt("matching_pi_genIdx", pi_genIdx);
        b_cand.addUserInt("matching_trg_mu_motherPdgId", genTriggerMuonMother_genPdgId);
        b_cand.addUserInt("matching_sel_mu_motherPdgId", genMuonMother_genPdgId);
        b_cand.addUserInt("matching_pi_motherPdgId", genPionMother_genPdgId);
        b_cand.addUserFloat("mupi_mass_reco_gen_reldiff", mupi_mass_reldiff);
        b_cand.addUserFloat("lxy_reco_gen_reldiff", lxy_reldiff);
        b_cand.addUserFloat("gen_lxy", gen_lxy);

        ret_val->push_back(b_cand);
              
      } // for(size_t sel_mu_idx = 0; sel_mu_idx < sel_muons->size(); ++sel_mu_idx)
      
    } // for(size_t pi_idx = 0; pi_idx < kaons->size(); ++pi_idx)

  } // for(size_t trg_mu_idx = 0; trg_mu_idx < trg_muons->size(); ++trg_mu_idx)
  
  evt.put(std::move(ret_val));
}

typedef BToMuLPiBuilder<pat::ETHMuon> BToMuMuPiBuilder;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToMuMuPiBuilder);
