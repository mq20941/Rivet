import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('standard')
options.register('runOnly', '', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Run only specified analysis")
options.register('yodafile', 'test.yoda', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Name of yoda output file")
options.setDefault('maxEvents', 100)
if(hasattr(sys, "argv")):
    options.parseArguments()
print(options)

process = cms.Process("runRivetAnalysis")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("EmptySource")

process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    nEvents = cms.untracked.uint32(options.maxEvents),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh'),
    numberOfParameters = cms.uint32(1),
    args = cms.vstring('/afs/cern.ch/work/m/mseidel/generator/CMSSW_7_1_20/src/TT_hdamp_TuneT4_noweights_NNPDF30_13TeV_powheg_hvq.tgz')
    #args = cms.vstring('/afs/cern.ch/work/m/mseidel/generator/TT_hdamp_noweights_NNPDF31_13TeV_powheg_hvq.tgz')
)

from Configuration.Generator.Pythia8CommonSettings_cfi import *
#from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *
process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    PythiaParameters = cms.PSet(
    pythia8CommonSettingsBlock,
    pythia8PowhegEmissionVetoSettingsBlock,
    pythia8CP5Settings = cms.vstring(
        'Tune:pp 14',
        'Tune:ee 7',
        'MultipartonInteractions:ecmPow=0.03344',
        'PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118',
        'MultipartonInteractions:bProfile=2',
        'MultipartonInteractions:pT0Ref=1.41',
        'MultipartonInteractions:coreRadius=0.7634',
        'MultipartonInteractions:coreFraction=0.63',
        'ColourReconnection:range=5.176',
        'SigmaTotal:zeroAXB=off',
        'SpaceShower:alphaSorder=2',
        'SpaceShower:alphaSvalue=0.118',
        'SigmaProcess:alphaSvalue=0.118',
        'SigmaProcess:alphaSorder=2',
        'MultipartonInteractions:alphaSvalue=0.118',
        'MultipartonInteractions:alphaSorder=2',
        'TimeShower:alphaSorder=2',
        'TimeShower:alphaSvalue=0.118',
    ),
    processParameters = cms.vstring(
                'POWHEG:nFinal = 2',
                'TimeShower:mMaxGamma = 1.0',
                'Tune:pp 14',
                'Tune:ee 7',
                'MultipartonInteractions:ecmPow=0.25208',
                'SpaceShower:alphaSvalue=0.1108',
                'PDF:pSet=LHAPDF6:NNPDF30_lo_as_0130',
                'MultipartonInteractions:pT0Ref=2.197139e+00',
                'MultipartonInteractions:expPow=1.608572e+00',
                'ColourReconnection:range=6.593269e+00',
    ),
    parameterSets = cms.vstring(
      'pythia8CommonSettings',
      'pythia8PowhegEmissionVetoSettings',
      #'pythia8CP5Settings',
      'processParameters'
      )
    )
)

#from Configuration.Generator.Pythia8CommonSettings_cfi import *
#from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
#from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *
#process.generator = cms.EDFilter("Pythia8HadronizerFilter",
#    maxEventsToPrint = cms.untracked.int32(0),
#    pythiaPylistVerbosity = cms.untracked.int32(0),
#    filterEfficiency = cms.untracked.double(1.0),
#    pythiaHepMCVerbosity = cms.untracked.bool(False),
#    comEnergy = cms.double(13000.),
#    PythiaParameters = cms.PSet(
#        pythia8CommonSettingsBlock,
#        pythia8CUEP8M1SettingsBlock,
#        pythia8PowhegEmissionVetoSettingsBlock,
#        processParameters = cms.vstring(
#            'POWHEG:nFinal = 2',
#            'TimeShower:mMaxGamma = 1.0',
#            'Tune:pp 14',
#            'Tune:ee 7',
#            'MultipartonInteractions:ecmPow=0.25208',
#            'SpaceShower:alphaSvalue=0.1108',
#            'PDF:pSet=LHAPDF6:NNPDF30_lo_as_0130',
#            'MultipartonInteractions:pT0Ref=2.197139e+00',
#            'MultipartonInteractions:expPow=1.608572e+00',
#            'ColourReconnection:range=6.593269e+00',
#        ),
#        fsrParameters = cms.vstring(
#            'TimeShower:renormMultFac   = ' + str(options.scale**2),
#            'TimeShower:factorMultFac   = ' + str(options.scale**2),
#            'TimeShower:alphaSvalue     = ' + str(options.asfsr),
#            'TimeShower:MEcorrections   = ' + options.me,
#        ),
#        parameterSets = cms.vstring('pythia8CommonSettings',
#            #'pythia8CUEP8M1Settings',
#            'pythia8PowhegEmissionVetoSettings',
#            'processParameters',
#            'fsrParameters'
#        )
#    )
#)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

if options.runOnly:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(options.runOnly)
else:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(
        #'CMS_2016_I1434354', # diff xs lepton+jets
        #'MC_TTBAR', # MC analysis for lepton+jets
        #'MC_TOPMASS_LJETS', # MC analysis for lepton+jets top mass
        #'CMS_LesHouches2015', # MC analysis for dilepton
        #'MC_GENERIC', # MC generic analysis
        #'MC_XS', # MC xs analysis
        'CMS_2018_I1690148',  # jet substructure
    )
process.rivetAnalyzer.OutputFile      = options.yodafile
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection    = 831.76 # NNLO (arXiv:1303.6254)

process.lhe_step = cms.Path(process.externalLHEProducer)
process.generation_step = cms.Path(process.generator)
process.analysis_step = cms.Path(process.rivetAnalyzer)

process.schedule = cms.Schedule(process.lhe_step,process.generation_step,process.analysis_step)
