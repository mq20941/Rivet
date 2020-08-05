import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('standard')
options.register('runOnly', '', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Run only specified analysis")
options.register('yodafile', 'test.yoda', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Name of yoda output file")
options.setDefault('maxEvents', 10000)
if(hasattr(sys, "argv")):
    options.parseArguments()
print(options)

process = cms.Process("runRivetAnalysis")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource")
process.source.fileNames = cms.untracked.vstring(
    '/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FEFC23FC-37C7-E811-97DA-0CC47AA53D86.root',
    '/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FEF89E38-48C7-E811-BF05-002590D9D8B8.root',
    '/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FEEAB1A9-FAC6-E811-991C-0025901AA5AE.root',
    '/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FEC9FD38-0EC7-E811-8C4B-00259029E81A.root',
    '/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FE984109-17C7-E811-8D1F-0025901AC0FA.root',
)

process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles"),
)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("mergedGenParticles"),
    genEventInfo = cms.InputTag("generator", "", "SIM"),
    signalParticlePdgIds = cms.vint32()
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

if options.runOnly:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(options.runOnly)
else:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(
        'CMS_2018_I1711625', 'CMS_2019_I1753680', 'MC_WEIGHTS', 'MC_XS',
    )
process.rivetAnalyzer.OutputFile      = options.yodafile
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection = 6077.22

process.p = cms.Path(process.mergedGenParticles*process.generator*process.rivetAnalyzer)


