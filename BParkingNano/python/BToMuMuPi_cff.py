import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import uint, ufloat, Var, CandVars

BToMuMuPi = cms.EDProducer(
    'BToMuMuPiBuilder',
    trgMuons                = cms.InputTag('muonTrgSelector', 'trgMuons'),
    selMuons                = cms.InputTag('muonTrgSelector', 'SelectedMuons'), 
    selMuonsTransientTracks = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'), 
    pions                   = cms.InputTag('tracksBPark', 'SelectedTracks'),
    pionsTransientTracks    = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks                  = cms.InputTag("packedPFCandidates"), 
    lostTracks              = cms.InputTag("lostTracks")        , 
    beamSpot                = cms.InputTag("offlineBeamSpot")   , 

    # pre-fitter preselection
    pionSelection           = cms.string('pt > 0.55 && abs(eta)<2'), # pion preselection to be modified here 
    isoTracksSelection      = cms.string('pt > 0.55 && abs(eta)<2'),
    trgMuonSelection        = cms.string('pt > 5 && abs(eta) < 1.7'),
    selMuonSelection        = cms.string('pt > 1.5 && abs(eta) < 2'),
    preVtxSelection = cms.string(' & '.join([
        'pt > 2',
        'mass > 0.2',        
        'mass < 7.0',        
        ])
    ), # applied on the HNL cand 

    # post-fitter preselection
    postVtxSelection = cms.string(' & '.join([
        'userInt("hnl_vtx_OK") == 1',
        'userFloat("hnl_vtx_prob") > 0.0001',
        'mass < 8',
        'pt > 11',
        'abs(eta) < 1.7',
        'userFloat("hnl_vtx_chi2") < 9',
        'userFloat("trg_muon_sip3d") > 0.8',
        ##'userFloat("sel_muon_ip3d") > 0.0015',
        ##'abs(userFloat("sel_muon_dxy")) > 0.0005',
        ##'abs(userFloat("pion_dz")) > 0.001',
        ##'abs(userFloat("pion_dxy")) > 0.0003',
        'abs(userFloat("pion_dzS")) > 0.5',
        'abs(userFloat("pion_dxyS")) > 0.2',
        'abs(userFloat("pion_DCASig")) > 0.3',
        ##'userFloat("muons_Lxyz") > 0.005',
        ##'userFloat("pion_muon_vzdiff") > 0.001',
        ##'userFloat("dr_mu_pi") < 1.8',
        'userFloat("dr_trgmu_hnl") < 0.8',
        'userFloat("hnl_fitted_mass") > 0.5',
        'userFloat("hnl_fitted_mass") < 6.5',
        'userFloat("hnl_fitted_pt") > 4',
        'abs(userFloat("hnl_fitted_eta")) < 1.8',
        'userFloat("hnl_fitted_cos_theta_2D") > 0.95',
        ])
    ), # applied on the B cand
)
    

BToMuMuPiTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToMuMuPi"),
    cut = cms.string(""),
    name = cms.string("b"),
    doc = cms.string("HNL Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        trg_mu_idx      = uint('trg_mu_idx'),
        sel_mu_idx      = uint('sel_mu_idx'),
        pi_idx          = uint('pi_idx'    ),
        # mu-pi fit and vtx info
        sv_chi2         = ufloat('hnl_vtx_chi2' ),
        sv_prob         = ufloat('hnl_vtx_prob' ),
        sv_lxy          = ufloat('hnl_l_xy'     ),
        sv_lxye         = ufloat('hnl_l_xy_unc' ),
        sv_lxy_sig      = ufloat('hnl_ls_xy'    ),
        sv_x            = ufloat('hnl_vtx_x'    ),
        sv_y            = ufloat('hnl_vtx_y'    ),
        sv_z            = ufloat('hnl_vtx_z'    ),
        sv_xe           = ufloat('hnl_vtx_ex'   ), ## only saving diagonal elements of the cov matrix
        sv_ye           = ufloat('hnl_vtx_ey'   ),
        sv_ze           = ufloat('hnl_vtx_ez'   ),
        # HNL 
        hnl_mass        = Var('daughter("hnl").mass()', float),
        hnl_pt          = Var('daughter("hnl").pt()'  , float),
        hnl_eta         = Var('daughter("hnl").eta()' , float),
        hnl_phi         = Var('daughter("hnl").phi()' , float),
        hnl_charge      = Var('daughter("hnl").charge()', int),
        # Cos(theta)
        hnl_cos2D       = ufloat('hnl_cos_theta_2D'       ),
        hnl_fit_cos2D   = ufloat('hnl_fitted_cos_theta_2D'),
        # post-fit momentum
        hnl_fit_mass    = ufloat('hnl_fitted_mass'   ),
        hnl_fit_masserr = ufloat('hnl_fitted_massErr'),
        hnl_fit_pt      = ufloat('hnl_fitted_pt'     ),
        hnl_fit_eta     = ufloat('hnl_fitted_eta'    ),
        hnl_fit_phi     = ufloat('hnl_fitted_phi'    ),
        fit_mu_pt       = ufloat('hnl_fitted_mu_pt'  ),
        fit_mu_eta      = ufloat('hnl_fitted_mu_eta' ),
        fit_mu_phi      = ufloat('hnl_fitted_mu_phi' ),
        fit_mu_mass     = ufloat('hnl_fitted_mu_mass'),
        fit_pi_pt       = ufloat('hnl_fitted_pi_pt'  ),
        fit_pi_eta      = ufloat('hnl_fitted_pi_eta' ),
        fit_pi_phi      = ufloat('hnl_fitted_pi_phi' ),
        fit_pi_mass     = ufloat('hnl_fitted_pi_mass'),
        # trigger muon
        trg_mu_pt       = ufloat('trg_muon_pt'       ), 
        trg_mu_eta      = ufloat('trg_muon_eta'      ), 
        trg_mu_phi      = ufloat('trg_muon_phi'      ), 
        # vertex difference between the two muons
        mu_vxdiff       = ufloat('muons_vxdiff'      ),
        mu_vydiff       = ufloat('muons_vydiff'      ),
        mu_vzdiff       = ufloat('muons_vzdiff'      ),
        mu_Lxy          = ufloat('muons_Lxy'         ),
        mu_Lxyz         = ufloat('muons_Lxyz'        ),
        # vertex difference between the trigger muon and pion
        pi_mu_vzdiff    = ufloat('pion_muon_vzdiff'  ),
        # Id WP of the selected muon
        sel_mu_isSoft   = ufloat('sel_muon_isSoft'   ),
        sel_mu_isTight  = ufloat('sel_muon_isTight'  ),
        sel_mu_isMedium = ufloat('sel_muon_isMedium' ),
        sel_mu_isLoose  = ufloat('sel_muon_isLoose'  ),
        # impact paramaters
        trg_mu_ip3d     = ufloat('trg_muon_ip3d'  ), 
        trg_mu_sip3d    = ufloat('trg_muon_sip3d' ), 
        trg_mu_dxy      = ufloat('trg_muon_dxy'   ), 
        trg_mu_dz       = ufloat('trg_muon_dz'    ), 
        sel_mu_ip3d     = ufloat('sel_muon_ip3d'  ), 
        sel_mu_sip3d    = ufloat('sel_muon_sip3d' ), 
        sel_mu_dxy      = ufloat('sel_muon_dxy'   ), 
        sel_mu_dz       = ufloat('sel_muon_dz'    ), 
        pi_dz           = ufloat('pion_dz'        ), 
        pi_dxy          = ufloat('pion_dxy'       ), 
        pi_dzS          = ufloat('pion_dzS'       ), 
        pi_dxyS         = ufloat('pion_dxyS'      ), 
        pi_DCASig       = ufloat('pion_DCASig'    ), 
        # dR quantities
        dr_mu_pi        = ufloat('dr_mu_pi'       ),
        dr_trgmu_hnl    = ufloat('dr_trgmu_hnl'   ),
    )
)

CountBToMuMuPi = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToMuMuPi")
)    

BToMuMuPiSequence = cms.Sequence( BToMuMuPi )
#BToMuMuPiSequence = cms.Sequence( BToMuMuPi * CountBToMuMuPi )
