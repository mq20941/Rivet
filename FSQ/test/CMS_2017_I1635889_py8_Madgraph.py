import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")

process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
#from GeneratorInterface.RivetInterface.Pythia8CUEP8M2T4Settings_cfi import *
#from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *

'''process.source = cms.Source("LHESource",
        fileNames = cms.untracked.vstring(
               'file:lheFileInput'

	 )
)'''

process.source = cms.Source( "PoolSource",
                             fileNames = cms.untracked.vstring(
			     '/store/mc/RunIIWinter15wmLHE/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/LHE/MCRUN2_71_V1-v1/50000/009C282B-5CDC-E411-9499-0002C90EB9D8.root'
			     )
)

#from genAnalysis.analysis.CFIFILE import *
#from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *


#-------Produce charged jets-----------
from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

#----------------------------------------------

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
                         maxEventsToPrint = cms.untracked.int32(1),
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         filterEfficiency = cms.untracked.double(1.0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(13000),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
#	pythia8CUEP8M2T4SettingsBlock,
        JetMatchingParameters = cms.vstring(
            'JetMatching:setMad = off',
            'JetMatching:scheme = 1',
            'JetMatching:merge = on',
            'JetMatching:jetAlgorithm = 2',
            'JetMatching:etaJetMax = 5.',
            'JetMatching:coneRadius = 1.',
            'JetMatching:slowJetPower = 1',
            'JetMatching:qCut = 19.', #this is the actual merging scale
            'JetMatching:nQmatch = 5', #4 corresponds to 4-flavour scheme (no matching of b-quarks), 5 for 5-flavour scheme
            'JetMatching:nJetMax = 4', #number of partons in born matrix element for highest multiplicity
            'JetMatching:doShowerKt = off', #off for MLM matching, turn on for shower-kT matching
            'TimeShower:mMaxGamma = 4.0',
	    ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
#				    'pythia8CUEP8M2T4Settings',
                                    'JetMatchingParameters'
                                    )
        )
                         )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

process.rivetAnalyzer_CMSueDY = cms.EDAnalyzer('RivetAnalyzer',
  AnalysisNames = cms.vstring('CMS_2017_PAS_FSQ_16_008'),
  HepMCCollection = cms.InputTag('generatorSmeared'),
  UseExternalWeight = cms.bool(False),
  GenEventInfoCollection = cms.InputTag('generator'),
  CrossSection = cms.double(340),
  DoFinalize = cms.bool(True),
  ProduceDQMOutput = cms.bool(False),
  OutputFile = cms.string('cmsUEinDY_rootFileOutput.yoda')
)

process.load("Configuration.EventContent.EventContent_cff")
process.p = cms.Path(process.generator * process.pgen * process.rivetAnalyzer_CMSueDY)

