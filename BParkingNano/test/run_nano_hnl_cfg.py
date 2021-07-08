# E.g.: cmsRun run_nano_hnl_cfg.py isMC=true
from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms
from glob import glob

options = VarParsing('python')

options.register('isMC'           ,  True           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run this on real data"                  )
options.register('skipDuplicated' ,  True           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Skip duplicated events. True by default")
options.register('globalTag'      , 'NOTSET'        , VarParsing.multiplicity.singleton, VarParsing.varType.string, "Set global tag"                         )
options.register('wantSummary'    ,  True           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run this on real data"                  )
options.register('wantFullRECO'   ,  False          , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run this on real data"                  )
options.register('reportEvery'    ,  1000           , VarParsing.multiplicity.singleton, VarParsing.varType.int   , "report every N events"                  )
options.register('skip'           ,  0              , VarParsing.multiplicity.singleton, VarParsing.varType.int   , "skip first N events"                    )
options.register('inputFile'      , None            , VarParsing.multiplicity.singleton, VarParsing.varType.string, "inputFile name"                         )
options.register('outFile'        , 'bparknano.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, "outputFile name"                        )


options.setDefault('maxEvents', -1)
options.parseArguments()


globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
if options._beenSet['globalTag']:
    globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['BParkNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['BParkFullEvt', extension[options.isMC], options.tag])+'.root')


if not options.inputFiles:
    #options.inputFiles = ['/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v4/110000/A1ADA5DA-4A57-7945-9B9D-6FAC167A1627.root'] if options.isMC else \
    options.inputFiles = ['/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/F7E7EF39-476F-1C48-95F7-74CB5C7A542C.root'] if not options.isMC else \
                         ['file:%s' %i for i in glob('/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V20_Bc/mass3.0_ctau184.0/step4_nj1.root')]
                         #['file:%s' %i for i in glob('/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V20_emu/mass3.0_ctau184.0/step4_nj15.root')]

annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BParkNANO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.BParkingNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles) if not options.inputFile else cms.untracked.vstring('file:{}'.format(options.inputFile)),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
    duplicateCheckMode = cms.untracked.string('checkEachFile' if options.skipDuplicated else 'noDuplicateCheck'),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileFEVT,
    outputCommands = (cms.untracked.vstring('keep *',
                                            'drop *_*_SelectedTransient*_*',
                     )),
    splitLevel = cms.untracked.int32(0)
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    #fileName = outputFileNANO,
    fileName = outputFileNANO if not options.outFile else cms.untracked.string('file:{}'.format(options.outFile)),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep nanoaodFlatTable_*Table_*_*',     # event data
        'keep nanoaodUniqueString_nanoMetadata_*_*',   # basic metadata
        'keep nanoaodMergeableCounterTable_*Table_*_*',  # includes gentables
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

from PhysicsTools.BParkingNano.nanoBPark_cff import *
process = nanoAOD_customizeMuonTriggerBPark  (process)
process = nanoAOD_customizeTrackFilteredBPark(process)
process = nanoAOD_customizeBToMuMuPi         (process, isMC=options.isMC)
process = nanoAOD_customizeBToKMuMu          (process, isMC=options.isMC) 
process = nanoAOD_customizeTriggerBitsBPark  (process)

# Path and EndPath definitions
#process.nanoAOD_general_step = cms.Path(process.nanoSequence)
process.nanoAOD_MuMuPi_step = cms.Path(process.nanoSequence + process.nanoBMuMuPiSequence + CountBToMuMuPi )
process.nanoAOD_KMuMu_step  = cms.Path(process.nanoSequence + process.nanoBKMuMuSequence + CountBToKmumu ) 

# customisation of the process.
if options.isMC:
    from PhysicsTools.BParkingNano.nanoBPark_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process, ancestor_particles=[511, 521, 531, 541, 9900015]) 

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
    #process.nanoAOD_general_step,
    process.nanoAOD_MuMuPi_step,
    process.nanoAOD_KMuMu_step, 
    process.endjob_step, 
    process.NANOAODoutput_step
)
    
if options.wantFullRECO:
    process.schedule = cms.Schedule(
        #process.nanoAOD_general_step,
        process.nanoAOD_MuMuPi_step,
        process.nanoAOD_KMuMu_step, 
        process.endjob_step, 
        process.FEVTDEBUGHLToutput_step, 
        process.NANOAODoutput_step
    )
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('nanoAOD_MuMuPi_step', 'nanoAOD_KMuMu_step') 
    #SelectEvents = cms.vstring('nanoAOD_MuMuPi_step', 'nanoAOD_KMuMu_step', 'nanoAOD_general_step') 
    #SelectEvents = cms.vstring('nanoAOD_MuMuPi_step', 'nanoAOD_general_step') 
)

### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
