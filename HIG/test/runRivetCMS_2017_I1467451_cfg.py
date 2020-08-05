import FWCore.ParameterSet.Config as cms

process = cms.Process("runRivetAnalysis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",  
                fileNames = cms.untracked.vstring('/store/user/lviliani/ggHPowheg_JHUGen_Pythia8/GEN/ggHPowheg_JHUGen_Pythia8_101.root')
                # more files in /store/user/lviliani/ggHPowheg_JHUGen_Pythia8/GEN/
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2017_I1467451')
process.rivetAnalyzer.HepMCCollection = cms.InputTag('generator','')
process.rivetAnalyzer.OutputFile = cms.string('mcfile.yoda')
process.rivetAnalyzer.CrossSection = 439.8 

process.p = cms.Path(process.rivetAnalyzer)
