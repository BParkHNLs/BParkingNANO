import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *


# make sure that the BPark HLT lines correspond to the 10 first elements
#Path=["HLT_Mu7_IP4","HLT_Mu8_IP6","HLT_Mu8_IP5","HLT_Mu8_IP3","HLT_Mu8p5_IP3p5","HLT_Mu9_IP6","HLT_Mu9_IP5","HLT_Mu9_IP4","HLT_Mu10p5_IP3p5","HLT_Mu12_IP6", "HLT_Mu8_v1", "HLT_Mu8_v12", "HLT_Mu7p5_Track7_Jpsi_v11", "HLT_L2Mu23NoVtx_2Cha_v1", "HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1", "DST_DoubleMu1_noVtx_CaloScouting_v2", "DST_DoubleMu3_noVtx_CaloScouting_v6", "DST_DoubleMu3_noVtx_Mass10_PFScouting_v3", "HLT_BTagMu_AK4DiJet40_Mu5_v13"]
Path=["HLT_Mu7_IP4","HLT_Mu8_IP6","HLT_Mu8_IP5","HLT_Mu8_IP3","HLT_Mu8p5_IP3p5","HLT_Mu9_IP6","HLT_Mu9_IP5","HLT_Mu9_IP4","HLT_Mu10p5_IP3p5","HLT_Mu12_IP6"]
#Path=["HLT_Mu12_IP6"]

muonTrgSelector = cms.EDProducer("MuonTriggerSelector",
                                 muonCollection = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
                                 displacedStandaloneMuonCollection = cms.InputTag("displacedStandAloneMuons"), #same collection as in NanoAOD                                                           
                                 vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'), 

                                 # added for restoring first trigger matching method
                                 bits = cms.InputTag("TriggerResults", "", "HLT"), # format: module, label, process
                                 triggerObjects = cms.InputTag("slimmedPatTrigger"),
                                 prescales = cms.InputTag("patTrigger"),
                                 beamSpot = cms.InputTag("offlineBeamSpot"),

                                 # trigger muon matching conditions
                                 max_deltaR_trigger_matching = cms.double(0.05),
                                 max_deltaPtRel_trigger_matching = cms.double(0.1),

                                 # add displaced standalone muons
                                 add_dsa = cms.bool(False),

                                 # DSA to slimmed muon matching conditions
                                 do_dsa_matching = cms.bool(False),
                                 max_deltaR_dsaToSlimmed_matching = cms.double(0.1),
                                 max_deltaPtRel_dsaToSlimmed_matching = cms.double(0.2),
                                 
                                 # selection for the selected muon
                                 selmu_ptMin = cms.double(1.),
                                 #selmu_ptMin = cms.double(0.3),
                                 selmu_absEtaMax = cms.double(2.3),
                                 #selmu_absEtaMax = cms.double(2.8),
                                 HLTPaths=cms.vstring(Path),
                                 #L1seeds=cms.vstring(Seed),
                             )
#cuts minimun number in B both mu and e, min number of trg, dz muon, dz and dr track, 
countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("muonTrgSelector", "trgMuons")
)


muonBParkTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("muonTrgSelector:SelectedMuons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Muon"),
    doc  = cms.string("slimmedMuons + DSAMuons for BPark after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(CandVars,
        isSlimmedMuon = Var("isSlimmedMuon()", int, doc="muon is slimmedMuon"), 
        isDSAMuon = Var("isDSAMuon()", int, doc="muon is displacedStandaloneMuon"), 
        isMatchedToSlimmedMuon = Var("userInt('isMatchedToSlimmedMuon')", int, doc="DSA muon is matched to slimmed muon"), 
        indexMatchedSlimmedMuon = Var("userInt('indexMatchedSlimmedMuon')", int, doc="index of the slimmed muon to which the dsa muon is matched"), 
        dsaToSlimmedMatching_deltaR = Var("userFloat('dsaToSlimmedMatching_deltaR')", float, doc="deltaR between the dsa and slimmed muon"), 
        dsaToSlimmedMatching_deltaPtRel = Var("userFloat('dsaToSlimmedMatching_deltaPtRel')", float, doc="relative difference in pt between the dsa and matched slimmed muons"), 
        dsaToSlimmedMatching_deltadxyRel = Var("userFloat('dsaToSlimmedMatching_deltadxyRel')", float, doc="relative difference in dxy between the dsa and matched slimmed muons"), 
        dsaToSlimmedMatching_deltadzRel = Var("userFloat('dsaToSlimmedMatching_deltadzRel')", float, doc="relative difference in dz between the dsa and matched slimmed muons"), 
        vx = Var("vx()", float, doc="x coordinate of vertex position, in cm", precision=6),
        vy = Var("vy()", float, doc="y coordinate of vertex position, in cm", precision=6),
        vz = Var("vz()", float, doc="z coordinate of vertex position, in cm", precision=6),
        ##ptErr = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
        dz = Var("userFloat('dz')", float, doc="dz (with sign) wrt first PV, in cm", precision=10),
        ##dzErr = Var("userFloat('dzErr')", float, doc="dz uncertainty, in cm", precision=6),
        dzS = Var("userFloat('dzS')", float, doc="dz (with sign) significance wrt first PV", precision=10),
        dz_alaRdst = Var("userFloat('dz_alaRdst')", float, doc="dz (innerTrack) (with sign) wrt first PV (a la R(D*)), in cm", precision=10),
        dxy = Var("userFloat('dxy')", float, doc="dxy (with sign) wrt first PV, in cm", precision=10),
        ##dxyErr = Var("userFloat('dxyErr')", float, doc="dxy uncertainty, in cm", precision=6),
        dxyS = Var("userFloat('dxyS')", float, doc="dxy (with sign) significance wrt first PV", precision=10),
        dxy_BS_alaRdst = Var("userFloat('dxy_BS_alaRdst')", float, doc="dxy (with sign) wrt first BS (a la R(D*)), in cm", precision=10),
        dxyS_BS_alaRdst = Var("userFloat('dxyS_BS_alaRdst')", float, doc="dxy (with sign) significance wrt BS (a la R(D*))", precision=10),
        dxy_BS = Var("userFloat('dxy_BS')", float, doc="dxy (with sign) wrt BS, in cm", precision=10),
        dxyS_BS = Var("userFloat('dxyS_BS')", float, doc="dxy (with sign) significance wrt BS", precision=10),
        ##ip3d = Var("ip3d()", float, doc="3D impact parameter wrt first PV, in cm", precision=10),
        ##sip3d = Var("sip3d()", float, doc="3D impact parameter significance wrt first PV", precision=10),
        ##pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component"),
        #pfiso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))", float, doc="PF isolation dR=0.3, total (deltaBeta corrections)"),
        pfiso03_all = Var("pfiso03_all()", float, doc="PF isolation dR=0.3, total (deltaBeta corrections)"),
        #pfiso03Rel_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt", float, doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
        pfiso03Rel_all = Var("pfiso03Rel_all()", float, doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
        #pfiso03_trk = Var("isolationR03().sumPt", float, doc = "raw isolationR03"),
        #pfiso03Rel_trk = Var("isolationR03().sumPt/pt", float, doc = "raw isolationR03 relative"),
        #pfiso03_ch = Var("pfIsolationR03().sumChargedHadronPt", float, doc = "raw pf iso charged hadrons"),
        #pfiso03Rel_ch = Var("pfIsolationR03().sumChargedHadronPt/pt", float, doc = "raw pf iso charged hadrons"),
        #pfiso03_n = Var("pfIsolationR03().sumNeutralHadronEt", float, doc = "raw iso neutral hadron"),
        #pfiso03Rel_n = Var("pfIsolationR03().sumNeutralHadronEt/pt", float, doc = "raw iso neutral hadron"),
        #pfiso03_pho = Var("pfIsolationR03().sumPhotonEt", float, doc = "raw pf iso photons"),
        #pfiso03Rel_pho = Var("pfIsolationR03().sumPhotonEt/pt", float, doc = "raw pf iso photons"),
        #pfiso03_pu = Var("pfIsolationR03().sumPUPt", float, doc = "raw pf iso PU"),
        #pfiso03Rel_pu = Var("pfIsolationR03().sumPUPt/pt", float, doc = "raw pf iso PU"),
        #pfiso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))", float, doc="PF isolation dR=0.4, total (deltaBeta corrections)"),
        pfiso04_all = Var("pfiso04_all()", float, doc="PF isolation dR=0.4, total (deltaBeta corrections)"),
        #pfiso04Rel_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt", float, doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
        pfiso04Rel_all = Var("pfiso04Rel_all()", float, doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
        ##pfiso05_trk = Var("isolationR05().sumPt",float, doc = "raw isolationR03"),
        ##pfiso05Rel_trk = Var("isolationR05().sumPt/pt",float, doc = "raw isolationR03 relative"),
        #pfiso05_trk = Var("isolationR05().sumPt", float, doc = "raw isolationR05"),
        #pfiso05Rel_trk = Var("isolationR05().sumPt/pt", float, doc = "raw relative isolationR05"),
        #pfiso04_ch = Var("pfIsolationR04().sumChargedHadronPt", float, doc = "raw pf iso charged hadrons"),
        #pfiso04Rel_ch = Var("pfIsolationR04().sumChargedHadronPt/pt", float, doc = "raw relative pf iso charged hadrons"),
        #pfiso04_n= Var("pfIsolationR04().sumNeutralHadronEt", float, doc = "raw pf iso neutral hadron"),
        #pfiso04Rel_n = Var("pfIsolationR04().sumNeutralHadronEt/pt", float, doc = "raw pf iso neutral hadron"),
        #pfiso04_pho = Var("pfIsolationR04().sumPhotonEt", float, doc = "raw pf iso photons"),
        #pfiso04Rel_pho = Var("pfIsolationR04().sumPhotonEt/pt", float, doc = "raw pf iso photons"),
        #pfiso04_pu = Var("pfIsolationR04().sumPUPt", float, doc = "raw pf iso PU"),
        #pfiso04Rel_pu = Var("pfIsolationR04().sumPUPt/pt", float, doc = "raw pf iso PU"),

        ##tightCharge = Var("?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0",int,doc="Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)"),
        isPF = Var("isPFSlimmedMuon()", int, doc="muon is PF candidate"),
        isGlobalMuon = Var("isGlobalSlimmedMuon()", int, doc="muon is global muon"),
        isTrackerMuon = Var("isTrackerSlimmedMuon()", int, doc="muon is tracker muon"),
        isGlobalOrTrackerMuon = Var("isGlobalOrTrackerSlimmedMuon()", int, doc="muon is global muon or tracker muon"),
        isGlobalNotTrackerMuon = Var("isGlobalNotTrackerSlimmedMuon()", int, doc="muon is global muon and not tracker muon"),
        isTrackerNotGlobalMuon = Var("isTrackerNotGlobalSlimmedMuon()", int, doc="muon is tracker muon and not global muon"),
        looseId = Var("looseId()", int, doc="reco muon is Loose"),
        mediumId = Var("mediumId()", int, doc="cut-based ID, medium WP"),
        ##mediumPromptId = Var("passed('CutBasedIdMediumPrompt')",bool,doc="cut-based ID, medium prompt WP"),
        tightId = Var("tightId()", int, doc="cut-based ID, tight WP"),
        softId = Var("softId()", int, doc="soft cut-based ID"),
        pfIsoId = Var("pfIsoId()", "uint8", doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"),
        tkIsoId = Var("tkIsoId()", "uint8", doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"),  
        triggerIdLoose = Var("triggerIdLoose()", int, doc="TriggerIdLoose ID"),
        passDSAMuonID = Var("userInt('passDSAMuonID')", int, doc="DSA muon passes customised DSA muon ID"), 
        ##highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"),
        ##softMvaId = Var("passed('SoftMvaId')",bool,doc="soft MVA ID"),
        ##mvaId = Var("passed('MvaLoose')+passed('MvaMedium')+passed('MvaTight')","uint8",doc="Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight)"),
        ##mvaId = Var("passed('MvaLoose')+passed('MvaMedium')+passed('MvaTight')","uint8",doc="Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight)"),
        ##miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"),
        ##multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"),

        inTimeMuon = Var("inTimeMuon()", int, doc="inTimeMuon ID"),
        segmentCompatibility = Var("segmentCompatibility()", float, doc = "muon segment compatibility: propagating the tracker tracks to the muon system and evaluate the number of matched segments and the closeness of the matching", precision=14), # keep higher precision since people have cuts with 3 digits on this
        caloCompatibility = Var("caloCompatibility()", float, doc = "calorimetric compatibility"),
        validHitFraction = Var("validHitFraction()", float, doc = "fraction of hits a tracker track uses (among inner tracker layers it traverses)"),
        kinkFinderChi2 = Var("kinkFinderChi2()", float, doc = "chi2 of kink-finding algorithm: how likely it is that a track is made of more than one single track"),
        globalNormalisedChi2 = Var("globalNormalisedChi2()", float, doc = "chi2/ndof of global fit"),
        localPositionChi2 = Var("localPositionChi2()", float, doc = "chi2 of the position match between the tracker muon and the standalone muon"),
        trackerHighPurityFlag = Var("trackerHighPurityFlag()", int, doc = "tracker high-purity flag"), # int? 
        numberOfValidMuonHits = Var("numberOfValidMuonHits()", int, doc = "number of hits in the muon stations"),
        numberOfValidPixelHits = Var("numberOfValidPixelHits()", int, doc = "number of pixel hits"),
        numberOfTrackerLayers = Var("numberOfTrackerLayers()", int, doc = "number of tracker layers with hits"),
        numberOfPixelLayers = Var("numberOfPixelLayers()", int, doc = "number of pixel layers with hits"),
        numberOfStations = Var("numberOfStations()", int, doc = "number of matched stations with default arbitration (segment & track)"),

        is_triggering = Var("userInt('is_triggering')", int, doc="flag the reco muon is also triggering"),
        isTriggering = Var("userInt('isTriggering')", int, doc="flag the reco muon is also triggering"),
        isTriggeringBPark = Var("userInt('isTriggeringBPark')", int, doc="flag the reco muon is also triggering (only for BPark lines)"),
        matched_dr = Var("userFloat('DR')", float, doc="dr with the matched triggering muon" ),
        matched_dpt = Var("userFloat('DPT')", float, doc="dpt/pt with the matched triggering muon" ),   

        #triggered_HLT_Mu7_IP4 = Var("userInt('triggered_HLT_Mu7_IP4')", int, doc="reco muon triggered this trigger"),
        #triggered_HLT_Mu8_IP6 = Var("userInt('triggered_HLT_Mu8_IP6')", int, doc="reco muon triggered this trigger"),
        #triggered_HLT_Mu8_IP5 = Var("userInt('triggered_HLT_Mu8_IP5')", int, doc="reco muon triggered this trigger"),
        #triggered_HLT_Mu8_IP3 = Var("userInt('triggered_HLT_Mu8_IP3')", int, doc="reco muon triggered this trigger"),
        #triggered_HLT_Mu8p5_IP3p5 = Var("userInt('triggered_HLT_Mu8p5_IP3p5')", int, doc="reco muon triggered this trigger"),
        #triggered_HLT_Mu9_IP6 = Var("userInt('triggered_HLT_Mu9_IP6')", int, doc="reco muon triggered this trigger"),
        #triggered_HLT_Mu9_IP5 = Var("userInt('triggered_HLT_Mu9_IP5')", int, doc="reco muon triggered this trigger"),
        #triggered_HLT_Mu9_IP4 = Var("userInt('triggered_HLT_Mu9_IP4')", int, doc="reco muon triggered this trigger"),
        #triggered_HLT_Mu10p5_IP3p5 = Var("userInt('triggered_HLT_Mu10p5_IP3p5')", int, doc="reco muon triggered this trigger"),
        triggered_HLT_Mu12_IP6 = Var("userInt('triggered_HLT_Mu12_IP6')", int, doc="reco muon fired this trigger"),

        #fired_HLT_Mu7_IP4 = Var("userInt('HLT_Mu7_IP4')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu8_IP6 = Var("userInt('HLT_Mu8_IP6')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu8_IP5 = Var("userInt('HLT_Mu8_IP5')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu8_IP3 = Var("userInt('HLT_Mu8_IP3')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu8p5_IP3p5 = Var("userInt('HLT_Mu8p5_IP3p5')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu9_IP6 = Var("userInt('HLT_Mu9_IP6')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu9_IP5 = Var("userInt('HLT_Mu9_IP5')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu9_IP4 = Var("userInt('HLT_Mu9_IP4')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu10p5_IP3p5 = Var("userInt('HLT_Mu10p5_IP3p5')", int, doc="reco muon fired this trigger"),
        fired_HLT_Mu12_IP6 = Var("userInt('HLT_Mu12_IP6')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu8_v1 = Var("userInt('HLT_Mu8_v1')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu8_v12 = Var("userInt('HLT_Mu8_v12')", int, doc="reco muon fired this trigger"),
        #fired_HLT_Mu7p5_Track7_Jpsi_v11 = Var("userInt('HLT_Mu7p5_Track7_Jpsi_v11')", int, doc="reco muon fired this trigger"),
        #fired_HLT_L2Mu23NoVtx_2Cha_v1 = Var("userInt('HLT_L2Mu23NoVtx_2Cha_v1')", int, doc="reco muon fired this trigger"),
        #fired_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1 = Var("userInt('HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1')", int, doc="reco muon fired this trigger"),
        #fired_DST_DoubleMu1_noVtx_CaloScouting_v2 = Var("userInt('DST_DoubleMu1_noVtx_CaloScouting_v2')", int, doc="reco muon fired this trigger"),
        #fired_DST_DoubleMu3_noVtx_CaloScouting_v6 = Var("userInt('DST_DoubleMu3_noVtx_CaloScouting_v6')", int, doc="reco muon fired this trigger"),
        #fired_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3 = Var("userInt('DST_DoubleMu3_noVtx_Mass10_PFScouting_v3')", int, doc="reco muon fired this trigger"),
        #fired_HLT_BTagMu_AK4DiJet40_Mu5_v13 = Var("userInt('HLT_BTagMu_AK4DiJet40_Mu5_v13')", int, doc="reco muon fired this trigger"),

        #prescale_HLT_Mu7_IP4 = Var("userInt('HLT_Mu7_IP4_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu8_IP6 = Var("userInt('HLT_Mu8_IP6_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu8_IP5 = Var("userInt('HLT_Mu8_IP5_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu8_IP3 = Var("userInt('HLT_Mu8_IP3_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu8p5_IP3p5 = Var("userInt('HLT_Mu8p5_IP3p5_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu9_IP6 = Var("userInt('HLT_Mu9_IP6_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu9_IP5 = Var("userInt('HLT_Mu9_IP5_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu9_IP4 = Var("userInt('HLT_Mu9_IP4_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu10p5_IP3p5 = Var("userInt('HLT_Mu10p5_IP3p5_prescale')", int, doc="reco muon prescale this trigger"),
        prescale_HLT_Mu12_IP6 = Var("userInt('HLT_Mu12_IP6_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu8_v1 = Var("userInt('HLT_Mu8_v1_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu8_v12 = Var("userInt('HLT_Mu8_v12_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_Mu7p5_Track7_Jpsi_v11 = Var("userInt('HLT_Mu7p5_Track7_Jpsi_v11_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_L2Mu23NoVtx_2Cha_v1 = Var("userInt('HLT_L2Mu23NoVtx_2Cha_v1_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1 = Var("userInt('HLT_L2Mu23NoVtx_2Cha_CosmicSeed_v1_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_DST_DoubleMu1_noVtx_CaloScouting_v2 = Var("userInt('DST_DoubleMu1_noVtx_CaloScouting_v2_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_DST_DoubleMu3_noVtx_CaloScouting_v6 = Var("userInt('DST_DoubleMu3_noVtx_CaloScouting_v6_prescale')", int, doc="reco muon prescale this trigger"),
        #prescale_DST_DoubleMu3_noVtx_Mass10_PFScouting_v3 = Var("userInt('DST_DoubleMu3_noVtx_Mass10_PFScouting_v3_prescale')", int, doc="reco muon fired this trigger"),
        #prescale_HLT_BTagMu_AK4DiJet40_Mu5_v13 = Var("userInt('HLT_BTagMu_AK4DiJet40_Mu5_v13_prescale')", int, doc="reco muon prescale this trigger"),
    ),
)

muonsBParkMCMatchForTable = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = muonBParkTable.src,                         # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),                             # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.25),                           # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.25),                           # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(False),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
    motherPdgId = cms.vint32(9900015, 443, 511, 521, 531, 541),
)

muonBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = muonBParkTable.src,
    mcMap   = cms.InputTag("muonsBParkMCMatchForTable"),
    objName = muonBParkTable.name,
    objType = muonBParkTable.name, 
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

selectedMuonsMCMatchEmbedded = cms.EDProducer(
    'MuonMatchEmbedder',
    src = cms.InputTag('muonTrgSelector:SelectedMuons'),
    matching = cms.InputTag('muonsBParkMCMatchForTable')
)


muonTriggerBParkTable = muonBParkTable.clone(
    src = cms.InputTag("muonTrgSelector:trgMuons"),
    name = cms.string("TriggerMuon"),
    doc  = cms.string("HLT Muons matched with reco muons"), #reco muon matched to triggering muon"),
    variables = cms.PSet(CandVars,
        vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6)####################,
#       trgMuonIndex = Var("userInt('trgMuonIndex')", int,doc="index in trigger muon collection")
   )
)

# not used in the end
muonsTriggerBParkMCMatchForTable = cms.EDProducer("MCMatcher",# cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = muonTriggerBParkTable.src,                  # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),                             # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.1),                            # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.25),                           # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
    motherPdgId = cms.vint32(511, 521, 531, 541),
)

muonTriggerBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = muonTriggerBParkTable.src,
    mcMap   = cms.InputTag("muonsTriggerBParkMCMatchForTable"),
    objName = muonTriggerBParkTable.name,
    objType = muonTriggerBParkTable.name, 
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

triggerMuonsMCMatchEmbedded = cms.EDProducer(
    'TriggerMuonMatchEmbedder',
    src = cms.InputTag('muonTrgSelector', 'trgMuons'),
    matching = cms.InputTag('muonsTriggerBParkMCMatchForTable')
)

#muonBParkSequence = cms.Sequence(muonTrgSelector * countTrgMuons)
muonBParkSequence = cms.Sequence(muonTrgSelector)

muonBParkMC = cms.Sequence(muonTrgSelector + muonsBParkMCMatchForTable + selectedMuonsMCMatchEmbedded + muonBParkMCTable)
muonBParkMCWithTriggerMuon = cms.Sequence(muonBParkMC + muonsTriggerBParkMCMatchForTable + triggerMuonsMCMatchEmbedded + muonTriggerBParkMCTable)
muonBParkTables = cms.Sequence(muonBParkTable)
muonTriggerMatchedTables = cms.Sequence(muonTriggerBParkTable)
