#setting up environment
import sys 
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('standard')
options.register('runOnly', '', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Run only specified analysis")
options.register('yodafile', 'MC.yoda', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Name of yoda output file") #adds option to name output yoda file in command line, default name MC.yoda 
#options.register('model', 'CP5', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Name of input model") #adds option to input model in command line, default model CP5
options.setDefault('maxEvents', 100)
if(hasattr(sys, "argv")): 
    options.parseArguments()
# print(options)

process = cms.Process("runRivetAnalysis")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

if options.runOnly:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(options.runOnly)
else: 
     process.rivetAnalyzer.AnalysisNames = cms.vstring(
        'MC_ColourFlow',  
    )

process.rivetAnalyzer.OutputFile      = options.yodafile #sets up output file
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generatorSmeared")
process.rivetAnalyzer.CrossSection    = 831.76 # NNLO (arXiv:1303.6254)

process.p = cms.Path(process.rivetAnalyzer)

process.source.fileNames = [
'file:/hdfs/user/ec6821/Jonah/inputFiles/Herwig7/002D85B6-F78F-1442-BB71-5AF71962C986.root',
]
