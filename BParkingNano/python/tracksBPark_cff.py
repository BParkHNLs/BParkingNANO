import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

tracksBPark = cms.EDProducer('TrackMerger',
                             beamSpot   = cms.InputTag("offlineBeamSpot"),
                             trgMuon    = cms.InputTag('muonTrgSelector', 'trgMuons'),
                             tracks     = cms.InputTag("packedPFCandidates"),
                             lostTracks = cms.InputTag("lostTracks"),
                             muons      = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
                             pfElectrons= cms.InputTag("slimmedElectrons"),
                             vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                             #lowPtElectrons=cms.InputTag("slimmedLowPtElectrons"),
                             #gsf2packed=cms.InputTag("lowPtGsfLinks:packedCandidates"),
                             #gsf2lost=cms.InputTag("lowPtGsfLinks:lostTracks"),

                             # keep the cuts as tight as possible without affecting any physics builder
                             trkPtCut = cms.double(0.6),
                             #trkPtCut = cms.double(0.3),
                             trkEtaCut = cms.double(2.1),
                             #trkEtaCut = cms.double(2.8),

                             # clean tracks wrt trigger muons (checked that not relevant for our study)
                             do_trgmu_cleaning = cms.bool(False),
                             dzTrg_cleaning = cms.double(-1), # initial value: 1.8
                             #drTrg_cleaning = cms.double(0.03), # deltaR requirement between track and trigger muon  
                             drTrg_cleaning = cms.double(-1), # keep track even in trgmu jet, deltaR requirement between track and trigger muon  
                             # clean tracks wrt muons (checked that not relevant for our study)
                             do_mu_cleaning = cms.bool(False),
                             # clean tracks wrt electrons (checked that not relevant for our study)
                             do_el_cleaning = cms.bool(False),
                             # request PackedCandidate to have high purity
                             do_trk_highpurity = cms.bool(False),

                             dcaSig = cms.double(-100000),
                             trkNormChiMin = cms.int32(-1),
                             trkNormChiMax = cms.int32(-1)
                            )


trackBParkTable = cms.EDProducer(
    "SimpleCompositeCandidateFlatTableProducer",
    src = cms.InputTag("tracksBPark:SelectedTracks"),
    cut = cms.string(""),
    name = cms.string("ProbeTracks"),
    doc  = cms.string("track collection probe side for BPark after basic selection"),
    singleton = cms.bool(False),
    extension = cms.bool(False), 
    variables = cms.PSet(
         CandVars,
        vx = Var("vx()", float, doc="x coordinate of vertex position, in cm", precision=10),
        vy = Var("vy()", float, doc="y coordinate of vertex position, in cm", precision=10),
        vz = Var("vz()", float, doc="z coordinate of vertex position, in cm", precision=10),
        isPacked = Var("userInt('isPacked')",int,doc="track from packedCandidate collection", precision=10),
        isLostTrk = Var("userInt('isLostTrk')",int,doc="track from lostTrack collection", precision=10),
        dz = Var("userFloat('dz')",float,doc="dz (with sign) wrt first PV, in cm", precision=10),
        dxy = Var("userFloat('dxy')",float,doc="dxy (with sign) wrt first PV, in cm", precision=10),
        dzS = Var("userFloat('dzS')", float, doc="dz/err (with sign) wrt first PV, in cm", precision=10),
        dxyS = Var("userFloat('dxyS')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
        DCASig=Var("userFloat('DCASig')", float,doc="significance of xy-distance of closest approach wrt beamspot", precision=10),
        DCASig_corr=Var("userFloat('DCASig_corr')", float,doc="significance of xy-distance of closest approach wrt beamspot", precision=10),
        drTrg = Var("userFloat('drTrg')", float,doc="deltaR between the track and one trigger muon", precision=10),
        dzTrg = Var("userFloat('dzTrg')", float,doc="dz from the corresponding trigger muon, in cm", precision=10),
        isMatchedToMuon = Var("userInt('isMatchedToMuon')",bool,doc="track was used to build a muon", precision=10),
        isMatchedToLooseMuon = Var("userInt('isMatchedToLooseMuon')",bool,doc="track was used to build a muon passing LooseID", precision=10),
        isMatchedToSoftMuon = Var("userInt('isMatchedToSoftMuon')",bool,doc="track was used to build a muon passing softID", precision=10),
        isMatchedToMediumMuon = Var("userInt('isMatchedToMediumMuon')",bool,doc="track was used to build a muon passing mediumID", precision=10),
        isMatchedToEle = Var("userInt('isMatchedToEle')",bool,doc="track was used to build a PF ele", precision=10),
        #isMatchedToLowPtEle = Var("userInt('isMatchedToLowPtEle')",bool,doc="track was used to build a low-pT ele", precision=10),
        nValidHits = Var("userInt('nValidHits')", int,doc="Number of valid hits on track", precision=10),
        #dEdXStrip=Var("userFloat('dEdXStrip')", float,doc="dE/dX from strips of associated isolated track"),
        #dEdXPixel=Var("userFloat('dEdXPixel')", float,doc="dE/dX from pixels of associated isolated track"),
        chi2 = Var("userFloat('chi2')", float, doc="chi2 of the track fit", precision=10),
        ndof = Var("userInt('ndof')", int, doc="ndof of the track fit", precision=10),
        normalisedChi2 = Var("userFloat('normalisedChi2')", float, doc="chi2/ndof of the track fit", precision=10),
        numberOfValidHits = Var("userInt('numberOfValidHits')", int, doc="number of valid hits of the track", precision=10),
        numberOfLostHits = Var("userInt('numberOfLostHits')", int, doc="number of lost hits of the track", precision=10),
        numberOfValidPixelHits = Var("userInt('numberOfValidPixelHits')", int, doc="number of valid pixel hits", precision=10),
        numberOfTrackerLayers = Var("userInt('numberOfTrackerLayers')", int, doc="number of tracker layers", precision=10),
        numberOfPixelLayers = Var("userInt('numberOfPixelLayers')", int, doc="number of pixel layers", precision=10),
        qualityIndex = Var("userInt('qualityIndex')", int, doc="quality index", precision=10),
        highPurityFlag = Var("userInt('highPurityFlag')", int, doc="high purity flag", precision=10),
        validFraction = Var("userFloat('validFraction')", float, doc="valid fraction", precision=10),
        ),
)


tracksBParkMCMatchForTable = cms.EDProducer("MCMatcher",   # cut on deltaR, deltaPt/Pt; pick best by deltaR
#tracksBParkMCMatchForTable = cms.EDProducer("MCMatcherByPt",   # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = trackBParkTable.src,                     # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"),  # final mc-truth particle collection
    mcPdgId     = cms.vint32(321,211),                     # one or more PDG ID (321 = charged kaon, 211 = charged pion); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.15), #0.03       # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(False),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
    motherPdgId = cms.vint32(9900015, 521, 541),
    #minPt = tracksBPark.trkPtCut,
    #maxEta = tracksBPark.trkEtaCut,
)

tracksBParkMCMatchEmbedded = cms.EDProducer(
    'CompositeCandidateMatchEmbedder',
    src = trackBParkTable.src,
    matching = cms.InputTag("tracksBParkMCMatchForTable")
)

tracksBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = tracksBParkMCMatchForTable.src,
    mcMap   = cms.InputTag("tracksBParkMCMatchForTable"),
    objName = trackBParkTable.name,
    objType = trackBParkTable.name,
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 kaons or pions"),
)


tracksBParkSequence = cms.Sequence(tracksBPark)
tracksBParkTables = cms.Sequence(trackBParkTable)
tracksBParkMC = cms.Sequence(tracksBParkSequence + tracksBParkMCMatchForTable + tracksBParkMCMatchEmbedded)
tracksBParkMCWithTable = cms.Sequence(tracksBParkMC + tracksBParkMCTable)


