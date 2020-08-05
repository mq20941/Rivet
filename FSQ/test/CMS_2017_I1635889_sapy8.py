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
#from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from GeneratorInterface.RivetInterface.Pythia8CUEP8M2T4Settings_cfi import *
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *

'''process.source = cms.Source("LHESource",
        fileNames = cms.untracked.vstring(
               'file:lheFileInput'

	 )
)'''

'''process.source = cms.Source( "PoolSource",
                             fileNames = cms.untracked.vstring(
			     'rootFileInput'
			     )
)'''

process.source = cms.Source("EmptySource")
process.options = cms.untracked.PSet(

)

#from genAnalysis.analysis.CFIFILE import *
#from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *


#-------Produce charged jets-----------
from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

#----------------------------------------------

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
                         comEnergy = cms.double(13000.0),
                         crossSection = cms.untracked.double(6.44),
                         filterEfficiency = cms.untracked.double(1),
                         maxEventsToPrint = cms.untracked.int32(1),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
#        pythia8CUEP8M1SettingsBlock,
        pythia8CUEP8M2T4SettingsBlock,
        processParameters = cms.vstring(
            'WeakSingleBoson:ffbar2gmZ = on',
            '23:onMode = off',
            '23:onIfAny = 11',
            '23:onIfAny = 13',
            '23:onIfAny = 15',
            '23:mMin = 50.',
            ),
        parameterSets = cms.vstring('pythia8CommonSettings',
#                                    'pythia8CUEP8M1Settings',
				    'pythia8CUEP8M2T4Settings',
                                    'processParameters',
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
  OutputFile = cms.string('cmsUEinDY_rootFileOutput')
)

process.load("Configuration.EventContent.EventContent_cff")
process.p = cms.Path(process.generator * process.pgen * process.rivetAnalyzer_CMSueDY)

