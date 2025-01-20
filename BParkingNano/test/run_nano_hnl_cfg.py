# E.g.: cmsRun run_nano_hnl_cfg.py isMC=true

# info on GT: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable

from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms
from glob import glob

options = VarParsing('python')

options.register('isMC'                    , True            , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run this on real data"                  )
options.register('doSignal'                , True            , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run the BToMuMuPiBuilder"               )
options.register('doControl'               , False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run the BToKMuMuBuilder"                )
options.register('doHNL'                   , False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run the HNLToMuPiBuilder"               )
options.register('doTagAndProbe'           , False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run the TagAndProbeJpsiToMuMu"          )
options.register('doGeneral'               , False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run without builder"                    )
options.register('addTriggerMuonCollection', False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Add the TriggerMuon_* branches"         )
options.register('addProbeTracksCollection', False           , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Add the ProbeTracks_* branches"         )
options.register('skipDuplicated'          , True            , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Skip duplicated events. True by default")
options.register('globalTag'               ,'NOTSET'         , VarParsing.multiplicity.singleton, VarParsing.varType.string, "Set global tag"                         )
options.register('wantSummary'             , True            , VarParsing.multiplicity.singleton, VarParsing.varType.bool  , "Run this on real data"                  )
options.register('reportEvery'             , 1000            , VarParsing.multiplicity.singleton, VarParsing.varType.int   , "report every N events"                  )
options.register('skip'                    ,  0              , VarParsing.multiplicity.singleton, VarParsing.varType.int   , "skip first N events"                    )
options.register('inputFile'               , None            , VarParsing.multiplicity.singleton, VarParsing.varType.string, "inputFile name"                         )
options.register('outFile'                 , 'bparknano.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, "outputFile name"                        )

options.setDefault('maxEvents', -1)

options.parseArguments()

globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
if options._beenSet['globalTag']:
    globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['BParkNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['BParkFullEvt', extension[options.isMC], options.tag])+'.root')


if not options.inputFiles:
    options.inputFiles = ['file:/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/miniaod/data/copy_AM_part1/miniaod_D4_nj39.root']
    #options.inputFiles = ['/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v4/110000/A1ADA5DA-4A57-7945-9B9D-6FAC167A1627.root']
    #options.inputFiles = ['/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v4/110000/A1ADA5DA-4A57-7945-9B9D-6FAC167A1627.root'] if options.isMC else \
    #options.inputFiles = ['/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/F7E7EF39-476F-1C48-95F7-74CB5C7A542C.root'] if not options.isMC else \
    #options.inputFiles = ['/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/20000/72A32519-B235-654F-AB7A-3ABAFD0A0818.root'] if not options.isMC else \
    #options.inputFiles = ['/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/270000/8D4EC1FC-2E3D-934D-85FC-DC4B94925765.root'] if not options.isMC else \
    #options.inputFiles = ['/store/data/Run2018D/ParkingBPH1/MINIAOD/05May2019promptD-v1/270000/F85EA23D-7ACA-CC47-ABA5-3F0D8DFFE32E.root'] if not options.isMC else \
    #                     ['file:%s' %i for i in glob('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39/mass3.0_ctau100.0/step4_nj15.root')]

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
    #fileNames = cms.untracked.vstring('file:{}'.format(options.inputFiles[0])), #FIXME
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
process = nanoAOD_customizeMuonTriggerBPark      (process, addTriggerMuonCollection=options.addTriggerMuonCollection)
process = nanoAOD_customizeTrackFilteredBPark    (process, addProbeTracksCollection=options.addProbeTracksCollection)
process = nanoAOD_customizeBToMuMuPi             (process, isMC=options.isMC)
process = nanoAOD_customizeBToKMuMu              (process, isMC=options.isMC) 
process = nanoAOD_customizeHNLToMuPi             (process, isMC=options.isMC)
process = nanoAOD_customizeTagAndProbeJPsiToMuMu (process, isMC=options.isMC) 
process = nanoAOD_customizeTriggerBitsBPark      (process)

# Path and EndPath definitions
process.nanoAOD_general_step = cms.Path(process.nanoSequence)
process.nanoAOD_MuMuPi_step = cms.Path(process.nanoSequence + process.nanoBMuMuPiSequence + CountBToMuMuPi )
process.nanoAOD_KMuMu_step  = cms.Path(process.nanoSequence + process.nanoBKMuMuSequence + CountBToKmumu ) 
process.nanoAOD_HNLToMuPi_step = cms.Path(process.nanoSequence + process.nanoHNLToMuPiSequence + CountHNLToMuPi )
process.nanoAOD_JPsiToMuMu_step  = cms.Path(process.nanoSequence + process.nanoJPsiToMuMuSequence + CountJPsiToMuMu ) 

# customisation of the process.
if options.isMC:
    from PhysicsTools.BParkingNano.nanoBPark_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process, ancestor_particles=[511, 521, 531, 541, 9900015], addTriggerMuonCollection=options.addTriggerMuonCollection, addProbeTracksCollection=options.addProbeTracksCollection) 

if not options.isMC:
  # apply lumi mask, see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGoodLumiSectionsJSONFile#cmsRun
  import FWCore.PythonUtilities.LumiList as LumiList
  process.source.lumisToProcess = LumiList.LumiList(url='https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt').getVLuminosityBlockRange()

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
)
if options.doGeneral:
  process.schedule += cms.Schedule(
      process.nanoAOD_general_step,
  )
if options.doSignal:
  process.schedule += cms.Schedule(
      process.nanoAOD_MuMuPi_step,
  )
if options.doControl:
  process.schedule += cms.Schedule(
      process.nanoAOD_KMuMu_step, 
  )
if options.doHNL:
  process.schedule += cms.Schedule(
      process.nanoAOD_HNLToMuPi_step, 
  )
if options.doTagAndProbe:
  process.schedule += cms.Schedule(
      process.nanoAOD_JPsiToMuMu_step, 
  )
process.schedule += cms.Schedule(
    process.endjob_step, 
    process.NANOAODoutput_step
)
    
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process_string = cms.vstring()
if options.doGeneral: process_string.append('nanoAOD_general_step')
if options.doSignal: process_string.append('nanoAOD_MuMuPi_step')
if options.doControl: process_string.append('nanoAOD_KMuMu_step')
if options.doHNL: process_string.append('nanoAOD_HNLToMuPi_step')
if options.doTagAndProbe: process_string.append('nanoAOD_JPsiToMuMu_step')

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
    SelectEvents = process_string
)

### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
