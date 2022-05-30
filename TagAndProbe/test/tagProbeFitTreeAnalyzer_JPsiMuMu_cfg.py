import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os
import sys


options = VarParsing('analysis')
options.register('outLabel' , 'test'      , VarParsing.multiplicity.singleton, VarParsing.varType.string, "outLabel" )
options.register('outSuffix', 'chunk0_nj1', VarParsing.multiplicity.singleton, VarParsing.varType.string, "outSuffix")
options.parseArguments()

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )    

process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    #InputFileNames = cms.vstring("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/flat_bparknano_tag_and_probe_v1.root"),
    InputFileNames = cms.vstring(options.inputFiles),

    InputTreeName = cms.string("tree"),
    # output
    #OutputFileName = cms.string("results.root"),
    OutputFileName = cms.string("results_{}_{}.root".format(options.outLabel, options.outSuffix)),

    #number of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies whether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),

    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("mass", "2.9", "3.3", "GeV/c^{2}"),
        probe_pt = cms.vstring("probe_pt", "0", "100", "GeV/c"),
        probe_eta = cms.vstring("probe_eta", "0", "2.5", ""),
        probe_dxy_sig = cms.vstring("probe_dxy_sig", "0", "500", ""),
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        ismatched = cms.vstring("ismatched", "dummy[true=1,false=0]"),
        probe_istight = cms.vstring("probe_istight", "dummy[true=1,false=0]"),
        probe_fired_BParkingHLT = cms.vstring("probe_fired_BParkingHLT", "dummy[true=1,false=0]"),
    ),

    Cuts = cms.PSet(
        probe_pt = cms.vstring("probe_pt", "probe_pt", "10"),
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
        gaussPlusLinear = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.03,0.01,0.05])",
            "Chebychev::backgroundPass(mass, cPass[0,-1,1])",
            "Chebychev::backgroundFail(mass, cFail[0,-1,1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        gaussPlusQuadratic = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.03,0.01,0.05])",
            "Chebychev::backgroundPass(mass, {cPass1[0,-1,1], cPass2[0,-1,1]})",
            "Chebychev::backgroundFail(mass, {cFail1[0,-1,1], cFail2[0,-1,1]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        )
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
        # used for final scale factors
        ##the name of the parameter set becomes the name of the directory
        cat_pt_eta = cms.PSet(
            #specifies the efficiency of which category and state to measure 
            EfficiencyCategoryAndState = cms.vstring("probe_fired_BParkingHLT","true"),
            #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
            UnbinnedVariables = cms.vstring("mass"),
            #specifies the binning of parameters
            BinnedVariables = cms.PSet(
                probe_pt = cms.vdouble(6.0, 7.0, 8.0, 8.5, 9.0, 10.0, 10.5, 11.0, 12.0, 20.0, 100.0),
                probe_eta = cms.vdouble(0.0, 0.5, 1.0, 1.5, 2.0),
            ),
            #first string is the default followed by binRegExp - PDFname pairs
            #BinToPDFmap = cms.vstring("gaussPlusLinear", "*pt_bin0*", "gaussPlusQuadratic")
            BinToPDFmap = cms.vstring("gaussPlusLinear")
        ),
        # for better looking 2D plot
        #cat_pt_eta = cms.PSet(
        #    EfficiencyCategoryAndState = cms.vstring("probe_fired_BParkingHLT","true"),
        #    UnbinnedVariables = cms.vstring("mass"),
        #    BinnedVariables = cms.PSet(
        #        probe_pt = cms.vdouble(6.0, 7.0, 8.0, 8.5, 9.0, 10.0, 10.5, 11.0, 12.0, 20.0),
        #        probe_eta = cms.vdouble(0.0, 0.5, 1.0, 1.5, 2.0),
        #    ),
        #    BinToPDFmap = cms.vstring("gaussPlusLinear")
        #),
        #cat_pt_dxysig = cms.PSet(
        #    EfficiencyCategoryAndState = cms.vstring("probe_fired_BParkingHLT","true"),
        #    UnbinnedVariables = cms.vstring("mass"),
        #    BinnedVariables = cms.PSet(
        #        probe_pt = cms.vdouble(6.0, 7.0, 8.0, 8.5, 9.0, 10.0, 10.5, 11.0, 12.0, 20.0, 100.0),
        #        probe_dxy_sig = cms.vdouble(0.0, 4.0, 6.0, 8.0, 10.0, 20.0, 500.0),
        #    ),
        #    BinToPDFmap = cms.vstring("gaussPlusLinear")
        #),
        #cat_pt_eta_dxysig = cms.PSet(
        #    EfficiencyCategoryAndState = cms.vstring("probe_fired_BParkingHLT","true"),
        #    UnbinnedVariables = cms.vstring("mass"),
        #    BinnedVariables = cms.PSet(
        #        probe_pt = cms.vdouble(6.0, 7.0, 8.0, 8.5, 9.0, 10.0, 10.5, 11.0, 12.0, 20.0, 100.0),
        #        probe_eta = cms.vdouble(0.0, 0.5, 1.0, 1.5, 2.0),
        #        probe_dxy_sig = cms.vdouble(0.0, 4.0, 6.0, 8.0, 10.0, 20.0, 500.0),
        #    ),
        #    BinToPDFmap = cms.vstring("gaussPlusLinear")
        #),
    )
)

process.fitness = cms.Path(
    process.TagProbeFitTreeAnalyzer
)

