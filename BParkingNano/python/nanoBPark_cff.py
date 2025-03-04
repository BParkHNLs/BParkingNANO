from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.BParkingNano.trgbits_cff import *



##for gen and trigger muon
from PhysicsTools.BParkingNano.genparticlesBPark_cff import finalGenParticlesBPark, genParticleBParkTable, genParticleBParkSequence, genParticleBParkTables
from PhysicsTools.BParkingNano.particlelevelBPark_cff import *
from PhysicsTools.BParkingNano.triggerObjectsBPark_cff import *
from PhysicsTools.BParkingNano.muonsBPark_cff import * 

## filtered input collections
from PhysicsTools.BParkingNano.electronsBPark_cff import * 
from PhysicsTools.BParkingNano.tracksBPark_cff import *

## B collections
from PhysicsTools.BParkingNano.BsToPhiPhiTo4K_cff import *
#from PhysicsTools.BParkingNano.PhiToKK_cff import *
from PhysicsTools.BParkingNano.BToMuLPi_cff import BToMuMuPi, BToMuMuPiMC, BToMuMuPiTable, BToMuMuPiSequence, BToMuMuPiSequenceMC, CountBToMuMuPi
from PhysicsTools.BParkingNano.BToKLL_cff import *
from PhysicsTools.BParkingNano.HNLToMuPi_cff import HNLToMuPi, HNLToMuPiMC, HNLToMuPiTable, HNLToMuPiSequence, HNLToMuPiSequenceMC, CountHNLToMuPi
from PhysicsTools.BParkingNano.BToKstarLL_cff import *
from PhysicsTools.BParkingNano.TagAndProbeJPsiToMuMu_cff import * 



nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)

nanoSequence = cms.Sequence(nanoMetadata + 
                            vertexSequence +           
                            globalTables + vertexTables + 
                            triggerObjectBParkTables + l1bits)

nanoSequenceMC = cms.Sequence(particleLevelBParkSequence + genParticleBParkSequence + 
                              globalTablesMC + genWeightsTable + genParticleBParkTables + lheInfoTable) 



def nanoAOD_customizeMuonTriggerBPark(process, addTriggerMuonCollection=False):
    if addTriggerMuonCollection:
      process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonTriggerMatchedTables + muonBParkTables)
    else:
      process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonBParkTables)
    return process

def nanoAOD_customizeTrackFilteredBPark(process, addProbeTracksCollection=False):
    if addProbeTracksCollection:
      process.nanoSequence = cms.Sequence( process.nanoSequence + tracksBParkSequence + tracksBParkTables)
    else:
      process.nanoSequence = cms.Sequence( process.nanoSequence + tracksBParkSequence)
    return process

def nanoAOD_customizeElectronFilteredBPark(process):
    process.nanoBKeeSequence     = cms.Sequence( electronsBParkSequence + electronBParkTables)
    process.nanoBKstarEESequence = cms.Sequence( electronsBParkSequence + electronBParkTables)
    return process

def nanoAOD_customizeTriggerBitsBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + trgTables)
    return process

def nanoAOD_customizeBsToPhiPhiTo4K(process, isMC=False):
    if isMC == False:
      #process.nanoBsToPhiPhiTo4KSequence = cms.Sequence( BsToPhiPhiTo4KSequence + BsToPhiPhiTo4KTable )
      process.nanoBsToPhiPhiTo4KSequence = cms.Sequence( PhiToKKSequence + BsToPhiPhiTo4KSequence + PhiToKKTable + BsToPhiPhiTo4KTable )
    else:
      #process.nanoBsToPhiPhiTo4KSequence = cms.Sequence( BsToPhiPhiTo4KSequenceMC + BsToPhiPhiTo4KTable )
      #process.nanoBsToPhiPhiTo4KSequence = cms.Sequence( PhiToKKSequenceMC + BsToPhiPhiTo4KSequenceMC + PhiToKKTable + BsToPhiPhiTo4KTable )
      #FIXME
      #process.nanoBsToPhiPhiTo4KSequence = cms.Sequence( PhiToKKSequence + BsToPhiPhiTo4KSequence + PhiToKKTable + BsToPhiPhiTo4KTable ) # worked so far
      process.nanoBsToPhiPhiTo4KSequence = cms.Sequence( PhiToKKSequenceMC + BsToPhiPhiTo4KSequenceMC + PhiToKKTable + BsToPhiPhiTo4KTable )
    return process

#def nanoAOD_customizePhiToKK(process, isMC=False):
#    if isMC == False:
#      process.nanoPhiToKKSequence = cms.Sequence( PhiToKKSequence + PhiToKKTable )
#    else:
#      process.nanoPhiToKKSequence = cms.Sequence( PhiToKKSequenceMC + PhiToKKTable )
#    return process

def nanoAOD_customizeBToMuMuPi(process, isMC=False):
    if isMC == False:
      process.nanoBMuMuPiSequence = cms.Sequence( BToMuMuPiSequence + BToMuMuPiTable )
    else:
      process.nanoBMuMuPiSequence = cms.Sequence( BToMuMuPiSequenceMC + BToMuMuPiTable )
    return process

def nanoAOD_customizeHNLToMuPi(process, isMC=False):
    if isMC == False:
      process.nanoHNLToMuPiSequence = cms.Sequence( HNLToMuPiSequence + HNLToMuPiTable )
    else:
      process.nanoHNLToMuPiSequence = cms.Sequence( HNLToMuPiSequenceMC + HNLToMuPiTable )
    return process

def nanoAOD_customizeBToKMuMu(process, isMC=False):
    if isMC == False:
      process.nanoBKMuMuSequence = cms.Sequence( BToKMuMuSequence + BToKmumuTable )
    else:
      process.nanoBKMuMuSequence = cms.Sequence( BToKMuMuSequenceMC + BToKmumuTable )
    return process

def nanoAOD_customizeBToKLL(process):
    process.nanoBKeeSequence   = cms.Sequence( process.nanoBKeeSequence + BToKEESequence    + BToKeeTable   )
    process.nanoBKMuMuSequence = cms.Sequence( BToKMuMuSequence + BToKmumuTable )
    return process

def nanoAOD_customizeTagAndProbeJPsiToMuMu(process, isMC=False):
    if isMC == False:
      process.nanoJPsiToMuMuSequence = cms.Sequence( JPsiToMuMuSequence + JPsiToMuMuTable )
    else:
      process.nanoJPsiToMuMuSequence = cms.Sequence( JPsiToMuMuSequenceMC + JPsiToMuMuTable )
    return process

#three possibilities for K*LL
def nanoAOD_customizeBToKstarLL(process):
    process.nanoBKstarLLSequence   = cms.Sequence( KstarToKPiSequence + BToKstarLLSequence + KstarToKPiTable + BToKstarLLTables )
    return process

def nanoAOD_customizeBToKstarEE(process):
    process.nanoBKstarEESequence   = cms.Sequence( process.nanoBKstarEESequence + BToKstarEESequence + BToKstarEETable + KstarToKPiTable )
    return process

def nanoAOD_customizeBToKstarMuMu(process):
    process.nanoBKstarMuMuSequence = cms.Sequence( BToKstarMuMuSequence + BToKstarMuMuTable + KstarToKPiTable )
    return process

from FWCore.ParameterSet.MassReplace import massSearchReplaceAnyInputTag
#def nanoAOD_customizeMC(process, ancestor_particles=[511, 521, 531, 541, 9900015], addTriggerMuonCollection=False, addProbeTracksCollection=False):  
def nanoAOD_customizeMC(process, ancestor_particles=[531, 333], addTriggerMuonCollection=False, addProbeTracksCollection=False):  
    for name, path in process.paths.iteritems():
        # replace all the non-match embedded inputs with the matched ones
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:SelectedMuons', 'selectedMuonsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:trgMuons', 'triggerMuonsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'electronsForAnalysis:SelectedElectrons', 'selectedElectronsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'tracksBPark:SelectedTracks', 'tracksBParkMCMatchEmbedded')

        # make the PhiToKKTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'PhiToKK', 'PhiToKKMC')

        # make the BsToPhiPhiTo4KTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'BsToPhiPhiTo4K', 'BsToPhiPhiTo4KMC')

        # make the BToMuMuPiTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'BToMuMuPi', 'BToMuMuPiMC')
        
        # make the BToKMuMuTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'BToKmumu', 'BToKmumuMC')

        # make the HNLToMuPiTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'HNLToMuPi', 'HNLToMuPiMC')

        # make the JPsiToMuMuTable/count talk to the correct producer
        massSearchReplaceAnyInputTag(path, 'JPsiToMuMu', 'JPsiToMuMuMC')

        # save the all descendants of ancestor_particles
        to_save = ' || '.join(['abs(pdgId) == %d'%ipdg for ipdg in ancestor_particles])
        to_save = 'keep++ (%s)'%to_save
        finalGenParticlesBPark.select.append(to_save)

        # modify the path to include mc-specific info
        path.insert(0, nanoSequenceMC)
        if addTriggerMuonCollection:
          path.replace(process.muonBParkSequence, process.muonBParkMCWithTriggerMuon)
        else:
          path.replace(process.muonBParkSequence, process.muonBParkMC)
        path.replace(process.electronsBParkSequence, process.electronBParkMC)
        if addProbeTracksCollection:
          path.replace(process.tracksBParkSequence, process.tracksBParkMCWithTable)
        else:
          path.replace(process.tracksBParkSequence, process.tracksBParkMC)
