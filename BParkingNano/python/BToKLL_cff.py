import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

# NOT USED
electronPairsForKee = cms.EDProducer(
    'DiElectronBuilder',
    src = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    transientTracksSrc = cms.InputTag('electronsForAnalysis', 'SelectedTransientElectrons'),
    lep1Selection = cms.string('pt > 1.3'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string(
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
        '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03 && userInt("nlowpt")<2'
        
    ),
    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'),
)

# NOT USED
BToKee = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('electronPairsForKee'),
    leptonTransientTracks = electronPairsForKee.transientTracksSrc,
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    kaonSelection = cms.string(''),
    isoTracksSelection = cms.string('pt > 0.5 && abs(eta)<2.5'),
    preVtxSelection = cms.string(
        'pt > 1.75 && userFloat("min_dr") > 0.03 '
        '&& mass < 7. && mass > 4.'
        ),
    postVtxSelection = cms.string(
         'userInt("sv_OK") == 1 && userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.'
    )
)

# Di-muons for JPSi reconstruction
muonPairsForKmumu = cms.EDProducer(
    'DiMuonBuilder',
    isMC = cms.bool(False),
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection = cms.string(' && '.join([
        'pt > 7',                       # we require l1 to be the triggering one
        'abs(eta) < 2.'
      ])
    ),
    lep2Selection = cms.string(' && '.join([
        'pt > 1.0',
        'abs(eta) < 2.'
      ])
    ),
    preVtxSelection = cms.string(' && '.join([
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1.',
        'mass() < 5',
        'mass() > 0',
        'charge() == 0',
        'userFloat("lep_deltaR") > 0.03',
      ])
    ),
    postVtxSelection = cms.string(' && '.join([
        #'userInt("sv_OK") == 1',        # hard-coded in builder
        #'userFloat("sv_chi2") < 998',   # arbitrary, was removed
        'userFloat("sv_prob") > 0.01',   # can be tightened
        'userFloat("fitted_mass") > 2.', # can be tightened
        'userFloat("fitted_mass") < 4.', # can be tightened
      ])
    ),
    label = cms.string('muon'),
)


muonPairsForKmumuMC = muonPairsForKmumu.clone(
    isMC = cms.bool(True),
)


BToKmumu = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('muonPairsForKmumu'),
    leptonTransientTracks = muonPairsForKmumu.transientTracksSrc,
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    genParticles = cms.InputTag("finalGenParticlesBPark"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    isMC = cms.bool(False),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    kaonSelection = cms.string(''),
    isoTracksSelection = BToKee.isoTracksSelection,
    # This in principle can be different between electrons and muons
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& mass < 9. && mass > 4.'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 9.'
    )
)


BToKmumuMC = BToKmumu.clone(
    dileptons = cms.InputTag('muonPairsForKmumuMC'),
    leptonTransientTracks = muonPairsForKmumuMC.transientTracksSrc,
    isMC = cms.bool(True),
)


BToKeeTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToKee"),
    cut = cms.string(""),
    name = cms.string("BToKEE"),
    doc = cms.string("BToKEE Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),
        kIdx = uint('k_idx'),
        minDR = ufloat('min_dr'),
        maxDR = ufloat('max_dr'),
        # fit and vtx info
        chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        sv_ndof = ufloat('sv_ndof'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
        vtx_x = ufloat('vtx_x'),
        vtx_y = ufloat('vtx_y'),
        vtx_z = ufloat('vtx_z'),
        vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float), # this might not work
        mllErr_llfit = Var('userCand("dilepton").userFloat("fitted_massErr")', float), # this might not work
        mll_fullfit = ufloat('fitted_mll'),
        ll_sv_prob = ufloat('ll_sv_prob'), # probability of dimuon fit
        # Cos(theta)
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
        fit_l1_pt = ufloat('fitted_l1_pt'),
        fit_l1_eta = ufloat('fitted_l1_eta'),
        fit_l1_phi = ufloat('fitted_l1_phi'),
        fit_l2_pt = ufloat('fitted_l2_pt'),
        fit_l2_eta = ufloat('fitted_l2_eta'),
        fit_l2_phi = ufloat('fitted_l2_phi'),
        fit_k_pt = ufloat('fitted_k_pt'),
        fit_k_eta = ufloat('fitted_k_eta'),
        fit_k_phi = ufloat('fitted_k_phi'),
        l1_iso03 = ufloat('l1_iso03'),
        l1_iso04 = ufloat('l1_iso04'),
        l2_iso03 = ufloat('l2_iso03'),
        l2_iso04 = ufloat('l2_iso04'),
        k_iso03  = ufloat('k_iso03'),
        k_iso04  = ufloat('k_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
        l1_iso03_close = ufloat('l1_iso03_close'),
        l1_iso04_close = ufloat('l1_iso04_close'),
        l2_iso03_close = ufloat('l2_iso03_close'),
        l2_iso04_close = ufloat('l2_iso04_close'),
        k_iso03_close = ufloat('k_iso03_close'),
        k_iso04_close = ufloat('k_iso04_close'),
        b_iso03_close = ufloat('b_iso03_close'),
        b_iso04_close = ufloat('b_iso04_close'),
        k_pt = ufloat('kaon_pt'), 
        k_eta = ufloat('kaon_eta'), 
        k_phi = ufloat('kaon_phi'), 
        k_mass = ufloat('kaon_mass'), 
        k_charge = uint('kaon_charge'), 
        k_pdgid = uint('kaon_pdgId'), 
        k_vx = ufloat('kaon_vx'), 
        k_vy = ufloat('kaon_vy'), 
        k_vz = ufloat('kaon_vz'), 
        k_dz = ufloat('kaon_dz'), 
        k_dxy = ufloat('kaon_dxy'), 
        k_dzS = ufloat('kaon_dzS'), 
        k_dxyS = ufloat('kaon_dxyS'), 
        k_DCASig = ufloat('kaon_DCASig'), 
        k_ispacked = uint('kaon_ispacked'), 
        k_islost = uint('kaon_islost'), 
        k_highPurityFlag = uint('kaon_highPurityFlag'), 
        n_k_used = uint('n_k_used'),
        n_l1_used = uint('n_l1_used'),
        n_l2_used = uint('n_l2_used'),
        isMatched = Var("userInt('isMatched')", int, mcOnly=True),
        matching_l1_genIdx = Var("userInt('matching_l1_genIdx')", int, mcOnly=True),
        matching_l2_genIdx = Var("userInt('matching_l2_genIdx')", int, mcOnly=True),
        matching_k_genIdx = Var("userInt('matching_k_genIdx')", int, mcOnly=True),
        matching_l1_motherPdgId = Var("userInt('matching_l1_motherPdgId')", int, mcOnly=True),
        matching_l2_motherPdgId = Var("userInt('matching_l2_motherPdgId')", int, mcOnly=True),
        matching_k_motherPdgId = Var("userInt('matching_k_motherPdgId')", int, mcOnly=True),
        matched_b_pt = Var("userFloat('matched_b_pt')", float, mcOnly=True),
        matched_b_eta = Var("userFloat('matched_b_eta')", float, mcOnly=True),
        matched_b_phi = Var("userFloat('matched_b_phi')", float, mcOnly=True),
        matched_b_mass = Var("userFloat('matched_b_mass')", float, mcOnly=True),
        matched_l1_pt = Var("userFloat('matched_l1_pt')", float, mcOnly=True),
        matched_l1_eta = Var("userFloat('matched_l1_eta')", float, mcOnly=True),
        matched_l1_phi = Var("userFloat('matched_l1_phi')", float, mcOnly=True),
        matched_l1_mass = Var("userFloat('matched_l1_mass')", float, mcOnly=True),
        matched_l2_pt = Var("userFloat('matched_l2_pt')", float, mcOnly=True),
        matched_l2_eta = Var("userFloat('matched_l2_eta')", float, mcOnly=True),
        matched_l2_phi = Var("userFloat('matched_l2_phi')", float, mcOnly=True),
        matched_l2_mass = Var("userFloat('matched_l2_mass')", float, mcOnly=True),
        matched_k_pt = Var("userFloat('matched_k_pt')", float, mcOnly=True),
        matched_k_eta = Var("userFloat('matched_k_eta')", float, mcOnly=True),
        matched_k_phi = Var("userFloat('matched_k_phi')", float, mcOnly=True),
        matched_k_mass = Var("userFloat('matched_k_mass')", float, mcOnly=True),

        
    )
)

BToKmumuTable = BToKeeTable.clone(
    src = cms.InputTag("BToKmumu"),
    name = cms.string("BToKMuMu"),
    doc = cms.string("BToKMuMu Variable")
)


CountBToKee = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToKee")
) 


CountBToKmumu = CountBToKee.clone(
    minNumber = cms.uint32(1),
    src = cms.InputTag("BToKmumu")
)


BToKMuMuSequence = cms.Sequence(
    (muonPairsForKmumu * BToKmumu)
)

BToKMuMuSequenceMC = cms.Sequence(
    (muonPairsForKmumuMC * BToKmumuMC)
)

BToKEESequence = cms.Sequence(
    (electronPairsForKee * BToKee)
)

BToKLLSequence = cms.Sequence(
    (electronPairsForKee * BToKee) +
    (muonPairsForKmumu * BToKmumu)
)

BToKLLTables = cms.Sequence(BToKeeTable + BToKmumuTable)

