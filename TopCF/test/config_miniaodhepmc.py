import sys
import FWCore.ParameterSet.Config as cms # imports CMS-specific Python classes and functions 
import FWCore.ParameterSet.VarParsing as VarParsing # handles parsing of commandline objects for you

options = VarParsing.VarParsing ('standard')
options.register('runOnly', '', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Run only specified analysis") # single input object, string 
options.register('yodafile', 'MC.yoda', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Name of yoda output file")
options.setDefault('maxEvents', 100)
if(hasattr(sys, "argv")):
    options.parseArguments()
# print(options)

process = cms.Process("runRivetAnalysis") # process object aggregates all the configuration information for the cmsRun exectuable

# attach imported objects to the process 
process.load('Configuration.StandardSequences.Services_cff') # cff=configurationfragment

process.load("FWCore.MessageLogger.MessageLogger_cfi") # cfi=moduleinitialisaton - cfi files are used to initialize modules with all of their parameter settings 
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

# modules (configured with C++ type name and parameter set for each object)
# module = dynamically loadable processing component 
# EDProducer ctreates new data to placed in the Event, reconstructs data into HepMC format to be fed into Rivet routine
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"), #uncompressed data 
    inputPacked = cms.InputTag("packedGenParticles"), #compressed data
)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi") # takes merged particle collection and converts into HepMC format
process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("mergedGenParticles"),
    genEventInfo = cms.InputTag("generator", "", "SIM"),
    signalParticlePdgIds = cms.vint32()
)
#

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

if options.runOnly:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(options.runOnly)
else:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(
        'MC_ColourFlow',  
    )
process.rivetAnalyzer.OutputFile      = options.yodafile
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection    = 831.76 # NNLO (arXiv:1303.6254)

process.p = cms.Path(process.mergedGenParticles*process.generator*process.rivetAnalyzer) # specify modules to execute in order given

process.source.fileNames = [
'file:/hdfs/user/ec6821/Jonah/inputFiles/Herwig7/002D85B6-F78F-1442-BB71-5AF71962C986.root',
]
