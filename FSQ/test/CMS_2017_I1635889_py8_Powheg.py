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
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *

process.source = cms.Source("LHESource",
        fileNames = cms.untracked.vstring(
               'file:lheFileInput'

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
                         comEnergy = cms.double(energyS),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
#	pythia8CUEP8M2T4SettingsBlock,
        pythia8PowhegEmissionVetoSettingsBlock,
        processParameters = cms.vstring(
            'POWHEG:nFinal = 3',   ## Number of final state particles
                                   ## (BEFORE THE DECAYS) in the LHE
                                   ## other than emitted extra parton
 #           '25:onMode = off',      # turn OFF all H decays
 #           '25:onIfMatch = 23 23', # turn ON H->ZZ
 #           '23:onMode = off',      # turn OFF all Z decays
 #           '23:onIfAny = 13',      # turn ON Z->mumu 
#             'Tune:pp 21',           #ATLAS A14 central tune with NNPDF2.3LO
        #     'Tune:pp 17',           #"ATLAS Tune AZ", is tuned to the pT spectrum of the Z^/gamma^*0 boson in a set of rapidity bins
        #    'PartonLevel:MPI = off'     
        #     'ColourReconnection:reconnect = off'
       ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
#				    'pythia8CUEP8M2T4Settings',
                                    'pythia8PowhegEmissionVetoSettings',
                                    'processParameters'
                                    )
        )
                         )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
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

