import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("runRivetAnalysis")

from TopQuarkAnalysis.TopEventProducers.tqafInputFiles_cff import relValTTbar
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(relValTTbar)
)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("genParticles"),
    genEventInfo = cms.InputTag("generator", "", "SIM"),
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring(
    'ParticleLevelProducerTemplate', # test/template for ParticleLevelProducer
)
process.rivetAnalyzer.OutputFile      = "ParticleLevelProducer.yoda"
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection    = 252.89 # NNLO (arXiv:1303.6254)

process.p = cms.Path(process.generator*process.rivetAnalyzer)
