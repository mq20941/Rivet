#import os
import sys
import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")

# Getting enviroment virables from set by the batch submiter
SEED = 1234 #int(os.getenv('SEED'))
EVENTS_PER_JOB = 1000 #int(os.getenv('EVENTS_PER_JOB'))
OUTPUT_FILE = "FSQ_16_011_PYTHIA8_CUETP8M1_MB_13TeV.yoda" #os.getenv('OUTPUT')

process.RandomNumberGeneratorService.generator.initialSeed = SEED
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(EVENTS_PER_JOB))

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.source = cms.Source("EmptySource")

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.generator = cms.EDFilter("Pythia8GeneratorFilter",
        comEnergy = cms.double(13000.0),
        crossSection = cms.untracked.double(78300000000),
        filterEfficiency = cms.untracked.double(1),
        maxEventsToPrint = cms.untracked.int32(0),
        pythiaHepMCVerbosity = cms.untracked.bool(False),
        pythiaPylistVerbosity = cms.untracked.int32(1),
        PythiaParameters = cms.PSet(
            processParameters = cms.vstring(
            'Main:timesAllowErrors    = 10000',

            'ParticleDecays:limitTau0 = on',
	        'ParticleDecays:tauMax = 10',

            # Processes:
            'SoftQCD:inelastic = on',

            # Tune parameters (CUETP8M1):
            'Tune:ee = 7',
            'Tune:pp = 14',
            'MultipartonInteractions:pT0Ref = 2.4024',
            'MultipartonInteractions:ecmPow = 0.25208',
            'MultipartonInteractions:expPow = 1.61',

            # Section to play:

            ),
    parameterSets = cms.vstring('processParameters')
	)
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.rivetAnalyzer.OutputFile = cms.string(OUTPUT_FILE)
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.AnalysisNames   = cms.vstring('CMS_2018_I1680318')

process.p = cms.Path(process.generator*process.rivetAnalyzer)
