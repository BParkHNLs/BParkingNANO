import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

PhiToKK = cms.EDProducer(
    'PhiToKKBuilder',
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    genParticles = cms.InputTag("finalGenParticlesBPark"),
    isMC = cms.bool(False),
    #isMC = cms.bool(True), #FIXME
    beamSpot = cms.InputTag("offlineBeamSpot"),

    k1_selection = cms.string(''), #cms.string('pt > 1.5'),
    k2_selection = cms.string(''),
    #preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
    #                             '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'),
    #preVtxSelection = cms.string('userFloat("lep_deltaR") > 0.03'),

    #pre_vtx_selection_phi = cms.string('abs(mass-1.0195)<0.05 && userFloat("deltaR_prefit")<0.5'), #TODO add sum/product of pt?
    #post_vtx_selection_phi = cms.string('abs(userFloat("phi_fitted_mass")-1.0195)<0.015 && userFloat("phi_fitted_k1_pt")>0.8 && userFloat("phi_fitted_k2_pt")>0.7 && userFloat("deltaR_postfit")<0.25 && userFloat("phi_sv_prob")>0.01'),

    #post_vtx_selection_phi = cms.string('abs(userFloat("phi_fitted_mass")-1.0195)<0.015 && userFloat("phi_fitted_k1_pt")>0.8 && userFloat("phi_fitted_k2_pt")>0.7 && userFloat("deltaR_postfit")<0.25 && userFloat("phi_sv_prob")>0.01 && userFloat("phi_cos_theta_2D")>0.85'),

    # tight preselection
    #pre_vtx_selection_phi = cms.string('userFloat("deltaR_prefit") < 0.3 && abs(mass-1.0195)<0.015'), #TODO add charge=0? sum/product of pt? mass range, deltaR?
    #post_vtx_selection_phi = cms.string('userFloat("deltaR_postfit") < 0.2 && userFloat("phi_sv_prob") > 0.2 && userFloat("phi_fitted_k1_pt") > 1.5 && abs(userFloat("phi_fitted_k1_eta")) < 2 && userFloat("phi_k1_DCASig_corr") > 0.3 && userFloat("phi_fitted_k2_pt") > 1.2 && abs(userFloat("phi_fitted_k2_eta")) < 2 && userFloat("phi_k2_DCASig_corr") > 0.2 && userFloat("phi_fitted_pt") > 3.2 && abs(userFloat("phi_fitted_eta")) < 2 && userFloat("phi_cos_theta_2D")>0.99 && abs(userFloat("phi_fitted_mass")-1.0195)<0.008'),

    # no cuts
    pre_vtx_selection_phi = cms.string('abs(mass-1.0195)<0.1'),
    #pre_vtx_selection_phi = cms.string(''),
    #post_vtx_selection_phi = cms.string(''),
    post_vtx_selection_phi = cms.string('abs(userFloat("phi_fitted_mass")-1.0195)<0.1'),

    # loose preselection
    #pre_vtx_selection_phi = cms.string('abs(mass-1.0195)<0.03'), #TODO add charge=0? sum/product of pt? mass range, deltaR?
    #post_vtx_selection_phi = cms.string('userFloat("phi_sv_prob") > 0.01 && userFloat("phi_cos_theta_2D")>0.9 && abs(userFloat("phi_fitted_mass")-1.0195)<0.01'),

    # very loose preselection
    #pre_vtx_selection_phi = cms.string('abs(mass-1.0195)<0.2'), #TODO add charge=0? sum/product of pt? mass range, deltaR?
    ##post_vtx_selection_phi = cms.string('userFloat("phi_sv_prob") > 1e-5 && userFloat("phi_cos_theta_2D")>0.85 && abs(userFloat("phi_fitted_mass")-1.0195)<0.1'),
    #post_vtx_selection_phi = cms.string('userFloat("phi_sv_prob") > 1e-5 && abs(userFloat("phi_fitted_mass")-1.0195)<0.1'),
)

PhiToKKMC = PhiToKK.clone( 
    isMC = cms.bool(True),
)

BsToPhiPhiTo4K = cms.EDProducer(
    'BsToPhiPhiTo4KBuilder',
    # phi candidates
    phis = cms.InputTag('PhiToKK'),
    #phisTransientTracks = phis.kaonsTransientTracks,
    phisTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),

    #kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    #kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),

    genParticles = cms.InputTag("finalGenParticlesBPark"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    isMC = cms.bool(False),
    #isMC = cms.bool(True), #FIXME

    # selection Bs
    #pre_vtx_selection_Bs = cms.string('abs(mass-5.367)<0.4'),
    #pre_vtx_selection_Bs = cms.string('abs(mass-5.367)<0.6'),
    #post_vtx_selection_Bs = cms.string('userFloat("Bs_sv_prob") > 0.001 && userFloat("Bs_cos_theta_2D")>0.9 && abs(userFloat("Bs_fitted_mass")-5.367)<0.3'),

    #pre_vtx_selection_Bs = cms.string('abs(mass-5.367)<3'),
    #post_vtx_selection_Bs = cms.string(' && '.join([
    #    'abs(userFloat("phi1_fitted_mass") - 1.0195) < 0.012', 
    #    'abs(userFloat("phi2_fitted_mass") - 1.0195) < 0.012', 
    #    'userFloat("phi1_fitted_pt") > 2.5',
    #    'userFloat("phi2_fitted_pt") > 1.8', 
    #    '(userFloat("phi1_fitted_pt") * userFloat("phi1_fitted_pt")) > 6',
    #    'userFloat("Bs_fitted_pt") > 1.5',
    #    'abs(userFloat("Bs_fitted_eta")) < 2.5',
    #    'userFloat("deltaR_min") < 0.15',
    #    'userFloat("deltaR_max") < 2.5', 
    #    'userFloat("Bs_lxy_sig") > 1',
    #    'userFloat("Bs_sv_prob") > 0.001', 
    #    'userFloat("Bs_cos_theta_2D") > 0.9', 
    #    'abs(userFloat("Bs_fitted_mass")-5.367)<0.3',
    #    ])
    #),
        
    # no cuts
    pre_vtx_selection_Bs = cms.string('abs(mass-5.367)<3'),
    #pre_vtx_selection_Bs = cms.string(''),
    #post_vtx_selection_Bs = cms.string(''),
    post_vtx_selection_Bs = cms.string('abs(mass-5.367)<0.5'),

    # loose preselection
    #pre_vtx_selection_Bs = cms.string('abs(mass-5.367)<3'),
    #post_vtx_selection_Bs = cms.string('userFloat("Bs_sv_prob") > 1e-5 && abs(userFloat("Bs_fitted_mass")-5.367)<0.3'),
)


BsToPhiPhiTo4KMC = BsToPhiPhiTo4K.clone( 
    phis = cms.InputTag('PhiToKKMC'),
    isMC = cms.bool(True),
)


PhiToKKTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("PhiToKK"),
    cut = cms.string(""),
    name = cms.string("PhiToKK"),
    doc = cms.string("PhiToKK Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        phi_sv_chi2 = ufloat('phi_sv_chi2'),
        phi_sv_ndof = ufloat('phi_sv_ndof'),
        phi_sv_prob = ufloat('phi_sv_prob'),
        phi_vx = ufloat('phi_vx'),
        phi_vy = ufloat('phi_vy'),
        phi_vz = ufloat('phi_vz'),
        phi_cxx = ufloat('phi_vtx_cxx'),
        phi_cyy = ufloat('phi_vtx_cyy'),
        phi_czz = ufloat('phi_vtx_czz'),
        phi_cyx = ufloat('phi_vtx_cyx'),
        phi_czx = ufloat('phi_vtx_czx'),
        phi_czy = ufloat('phi_vtx_czy'),
        phi_mass = ufloat('phi_fitted_mass'),
        phi_masserr = ufloat('phi_fitted_massErr'),
        phi_pt = ufloat('phi_fitted_pt'),
        phi_eta = ufloat('phi_fitted_eta'),
        phi_phi = ufloat('phi_fitted_phi'),
        phi_charge = ufloat('phi_charge'),
        phi_cos2D = ufloat('phi_cos_theta_2D'),

        k1_idx = uint('k1_idx'),
        k2_idx = uint('k2_idx'),
        phi_k1_pt = ufloat('phi_fitted_k1_pt'),
        phi_k1_eta = ufloat('phi_fitted_k1_eta'),
        phi_k1_phi = ufloat('phi_fitted_k1_phi'),
        phi_k1_mass = ufloat('phi_fitted_k1_mass'),
        phi_k1_charge = ufloat('phi_k1_charge'),
        phi_k1_pdgid = ufloat('phi_k1_pdgId'),
        phi_k1_vx = ufloat('phi_k1_vx'),
        phi_k1_vy = ufloat('phi_k1_vy'),
        phi_k1_vz = ufloat('phi_k1_vz'),
        phi_k1_dz = ufloat('phi_k1_dz'),
        phi_k1_dxy = ufloat('phi_k1_dxy'),
        phi_k1_dzS = ufloat('phi_k1_dzS'),
        phi_k1_dxyS = ufloat('phi_k1_dxyS'),
        phi_k1_DCASig = ufloat('phi_k1_DCASig'),
        phi_k1_ptErr = ufloat('phi_k1_ptErr'),
        phi_k1_ispacked = uint('phi_k1_ispacked'),
        phi_k1_islost = uint('phi_k1_islost'),
        phi_k1_chi2 = ufloat('phi_k1_chi2'),
        phi_k1_normalisedChi2 = ufloat('phi_k1_normalisedChi2'),
        phi_k1_validFraction = ufloat('phi_k1_validFraction'),
        phi_k1_ndof = uint('phi_k1_ndof'),
        phi_k1_numberOfValidHits = uint('phi_k1_numberOfValidHits'),
        phi_k1_numberOfPixelHits = uint('phi_k1_numberOfPixelHits'),
        phi_k1_numberOfLostHits = uint('phi_k1_numberOfLostHits'),
        phi_k1_numberOfValidPixelHits = uint('phi_k1_numberOfValidPixelHits'),
        phi_k1_numberOfTrackerLayers = uint('phi_k1_numberOfTrackerLayers'),
        phi_k1_numberOfPixelLayers = uint('phi_k1_numberOfPixelLayers'),
        phi_k1_qualityIndex = uint('phi_k1_qualityIndex'),
        phi_k1_highPurityFlag = uint('phi_k1_highPurityFlag'),
        phi_k1_covQopQop = ufloat('phi_k1_covQopQop'),
        phi_k1_covLamLam = ufloat('phi_k1_covLamLam'),
        phi_k1_covPhiPhi = ufloat('phi_k1_covPhiPhi'),
        phi_k1_covQopLam = ufloat('phi_k1_covQopLam'),
        phi_k1_covQopPhi = ufloat('phi_k1_covQopPhi'),
        phi_k1_covLamPhi = ufloat('phi_k1_covLamPhi'),

        phi_k2_pt = ufloat('phi_fitted_k2_pt'),
        phi_k2_eta = ufloat('phi_fitted_k2_eta'),
        phi_k2_phi = ufloat('phi_fitted_k2_phi'),
        phi_k2_mass = ufloat('phi_fitted_k2_mass'),
        phi_k2_charge = ufloat('phi_k2_charge'),
        phi_k2_pdgid = ufloat('phi_k2_pdgId'),
        phi_k2_vx = ufloat('phi_k2_vx'),
        phi_k2_vy = ufloat('phi_k2_vy'),
        phi_k2_vz = ufloat('phi_k2_vz'),
        phi_k2_dz = ufloat('phi_k2_dz'),
        phi_k2_dxy = ufloat('phi_k2_dxy'),
        phi_k2_dzS = ufloat('phi_k2_dzS'),
        phi_k2_dxyS = ufloat('phi_k2_dxyS'),
        phi_k2_DCASig = ufloat('phi_k2_DCASig'),
        phi_k2_ptErr = ufloat('phi_k2_ptErr'),
        phi_k2_ispacked = uint('phi_k2_ispacked'),
        phi_k2_islost = uint('phi_k2_islost'),
        phi_k2_chi2 = ufloat('phi_k2_chi2'),
        phi_k2_normalisedChi2 = ufloat('phi_k2_normalisedChi2'),
        phi_k2_validFraction = ufloat('phi_k2_validFraction'),
        phi_k2_ndof = uint('phi_k2_ndof'),
        phi_k2_numberOfValidHits = uint('phi_k2_numberOfValidHits'),
        phi_k2_numberOfPixelHits = uint('phi_k2_numberOfPixelHits'),
        phi_k2_numberOfLostHits = uint('phi_k2_numberOfLostHits'),
        phi_k2_numberOfValidPixelHits = uint('phi_k2_numberOfValidPixelHits'),
        phi_k2_numberOfTrackerLayers = uint('phi_k2_numberOfTrackerLayers'),
        phi_k2_numberOfPixelLayers = uint('phi_k2_numberOfPixelLayers'),
        phi_k2_qualityIndex = uint('phi_k2_qualityIndex'),
        phi_k2_highPurityFlag = uint('phi_k2_highPurityFlag'),
        phi_k2_covQopQop = ufloat('phi_k2_covQopQop'),
        phi_k2_covLamLam = ufloat('phi_k2_covLamLam'),
        phi_k2_covPhiPhi = ufloat('phi_k2_covPhiPhi'),
        phi_k2_covQopLam = ufloat('phi_k2_covQopLam'),
        phi_k2_covQopPhi = ufloat('phi_k2_covQopPhi'),
        phi_k2_covLamPhi = ufloat('phi_k2_covLamPhi'),

        phi_lxy = ufloat('phi_lxy'),
        phi_lxy_sig = ufloat('phi_lxy_sig'),
        deltaR_prefit = ufloat('deltaR_prefit'),
        deltaR_postfit = ufloat('deltaR_postfit'),
        cos_theta_star_k1 = ufloat('cos_theta_star_k1'),
        cos_theta_star_k2 = ufloat('cos_theta_star_k2'),
        #energy_diff_phi_daughters_lab = ufloat('energy_diff_phi_daughters_lab'),
        px_diff_phi_daughters_lab = ufloat('px_diff_phi_daughters_lab'),
        py_diff_phi_daughters_lab = ufloat('py_diff_phi_daughters_lab'),
        pz_diff_phi_daughters_lab = ufloat('pz_diff_phi_daughters_lab'),
        energy_diff_prefitphi_daughters_lab = ufloat('energy_diff_prefitphi_daughters_lab'),
        px_diff_prefitphi_daughters_lab = ufloat('px_diff_prefitphi_daughters_lab'),
        py_diff_prefitphi_daughters_lab = ufloat('py_diff_prefitphi_daughters_lab'),
        pz_diff_prefitphi_daughters_lab = ufloat('pz_diff_prefitphi_daughters_lab'),
        energy_diff_phi_daughters_cm = ufloat('energy_diff_phi_daughters_cm'),
        p_daughters_cm = ufloat('p_daughters_cm'),
        #lep_vzdiff = ufloat('lep_vzdiff'),
        isMatched = uint('isMatched'),
    )
)

BsToPhiPhiTo4KTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BsToPhiPhiTo4K"),
    cut = cms.string(""),
    name = cms.string("BsToPhiPhiTo4K"),
    doc = cms.string("BsToPhiPhiTo4K Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        Bs_vx = ufloat('Bs_vx'),
        Bs_vy = ufloat('Bs_vy'),
        Bs_vz = ufloat('Bs_vz'),
        Bs_cxx = ufloat('Bs_vtx_cxx'),
        Bs_cyy = ufloat('Bs_vtx_cyy'),
        Bs_czz = ufloat('Bs_vtx_czz'),
        Bs_cyx = ufloat('Bs_vtx_cyx'),
        Bs_czx = ufloat('Bs_vtx_czx'),
        Bs_czy = ufloat('Bs_vtx_czy'),
        beamspot_x = ufloat('beamspot_x'),
        beamspot_y = ufloat('beamspot_y'),
        beamspot_z = ufloat('beamspot_z'),
        Bs_sv_chi2 = ufloat('Bs_sv_chi2'),
        Bs_sv_ndof = ufloat('Bs_sv_ndof'),
        Bs_sv_prob = ufloat('Bs_sv_prob'),
        Bs_mass = ufloat('Bs_fitted_mass'),
        Bs_massErr = ufloat('Bs_fitted_massErr'),
        Bs_mass_corr = ufloat('Bs_fitted_mass_corr'),
        Bs_pt = ufloat('Bs_fitted_pt'),
        Bs_eta = ufloat('Bs_fitted_eta'),
        Bs_phi = ufloat('Bs_fitted_phi'),
        Bs_charge = ufloat('Bs_charge'),
        Bs_cos2D = ufloat('Bs_cos_theta_2D'),
        Bs_lxy = ufloat('Bs_lxy'),
        Bs_lxy_sig = ufloat('Bs_lxy_sig'),

        k1_idx = uint('k1_idx'),
        k2_idx = uint('k2_idx'),
        k3_idx = uint('k3_idx'),
        k4_idx = uint('k4_idx'),
        phi1_idx = uint('phi1_idx'),
        phi2_idx = uint('phi2_idx'),

        deltaR_phi1phi2 = ufloat('deltaR_phi1phi2'), 
        deltaR_k1k3 = ufloat('deltaR_k1k3'),
        deltaR_k1k4 = ufloat('deltaR_k1k4'),
        deltaR_k2k3 = ufloat('deltaR_k2k3'),
        deltaR_k2k4 = ufloat('deltaR_k2k4'),
        deltaR_k1k2 = ufloat('deltaR_k1k2'),
        deltaR_k3k4 = ufloat('deltaR_k3k4'),
        deltaR_max = ufloat('deltaR_max'),
        deltaR_min = ufloat('deltaR_min'),

        k1k3_mass = ufloat('k1k3_mass'),
        k1k4_mass = ufloat('k1k4_mass'),
        k2k3_mass = ufloat('k2k3_mass'),
        k2k4_mass = ufloat('k2k4_mass'),
        k1k3_pt = ufloat('k1k3_pt'),
        k1k4_pt = ufloat('k1k4_pt'),
        k2k3_pt = ufloat('k2k3_pt'),
        k2k4_pt = ufloat('k2k4_pt'),

        cos_theta_star_phi1 = ufloat('cos_theta_star_phi1'),
        cos_theta_star_phi2 = ufloat('cos_theta_star_phi2'),
        Bs_beta = ufloat('Bs_beta'),

        phi1_mass = ufloat('phi1_fitted_mass'),
        phi1_pt = ufloat('phi1_fitted_pt'),
        phi1_eta = ufloat('phi1_fitted_eta'),
        phi1_phi = ufloat('phi1_fitted_phi'),
        phi2_mass = ufloat('phi2_fitted_mass'),
        phi2_pt = ufloat('phi2_fitted_pt'),
        phi2_eta = ufloat('phi2_fitted_eta'),
        phi2_phi = ufloat('phi2_fitted_phi'),

        k1_mass = ufloat('k1_fitted_mass'),
        k1_pt = ufloat('k1_fitted_pt'),
        k1_eta = ufloat('k1_fitted_eta'),
        k1_phi = ufloat('k1_fitted_phi'),
        k2_mass = ufloat('k2_fitted_mass'),
        k2_pt = ufloat('k2_fitted_pt'),
        k2_eta = ufloat('k2_fitted_eta'),
        k2_phi = ufloat('k2_fitted_phi'),
        k3_mass = ufloat('k3_fitted_mass'),
        k3_pt = ufloat('k3_fitted_pt'),
        k3_eta = ufloat('k3_fitted_eta'),
        k3_phi = ufloat('k3_fitted_phi'),
        k4_mass = ufloat('k4_fitted_mass'),
        k4_pt = ufloat('k4_fitted_pt'),
        k4_eta = ufloat('k4_fitted_eta'),
        k4_phi = ufloat('k4_fitted_phi'),
        #lep_vzdiff = ufloat('lep_vzdiff'),
        isMatched = uint('isMatched'),
    )
)

CountPhiToKK = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("PhiToKK")
) 

CountBsToPhiPhiTo4K = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BsToPhiPhiTo4K")
) 


PhiToKKSequence = cms.Sequence(PhiToKK)
PhiToKKSequenceMC = cms.Sequence(PhiToKKMC)
BsToPhiPhiTo4KSequence = cms.Sequence(BsToPhiPhiTo4K)
BsToPhiPhiTo4KSequenceMC = cms.Sequence(BsToPhiPhiTo4KMC)
BsToPhiPhiTo4KTables = cms.Sequence(BsToPhiPhiTo4KTable)

