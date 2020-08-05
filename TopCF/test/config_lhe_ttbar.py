import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys


options = VarParsing.VarParsing ('standard') 
options.register('outputFile', 'TOP-RunIIFall18wmLHEGS-00382_inLHE.root', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Name of output root file containing LHE events") 
options.setDefault('maxEvents', 1002)
if(hasattr(sys, "argv")): 
    options.parseArguments()

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('LHE',Run2_2018) # process variable refers to a Process object which aggregates all the configuration information for the cmsRun executable.

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeVEarly2018Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# number of events to process
# Is this section of code neccessary for LHE config file (as is not reading in any events)?
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents), #processing terminates after producing n events
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("EmptySource") # generating ME events from scratch 

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info 
# What does this section of code do?
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/GenProduction/python/TOP-RunIIFall18wmLHEGS-00382-fragment.py nevts:1002'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.LHEoutput = cms.OutputModule("PoolOutputModule", # OutputModule stores the data from the event
# PoolOutputModule  writes data to a standard CMS format ROOT file, which is CMSSW-readable 
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('LHE'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.LHEEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
# from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v11', '')

process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/hdfs/user/ec6821/Jonah/inputFiles/gridpack/TT_hdamp_NNPDF31_NNLO_ljets.tgz'),
    nEvents = cms.untracked.uint32(options.maxEvents),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)


# Path and EndPath definitions
process.lhe_step = cms.Path(process.externalLHEProducer)
process.LHEoutput_step = cms.EndPath(process.LHEoutput)

# Schedule definition
process.schedule = cms.Schedule(process.lhe_step,process.LHEoutput_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

# filter all path with the production filter sequence

for path in process.paths:
	if path in ['lhe_step']: continue
	getattr(process,path).insert(0, process.ProductionFilterSequence)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions

# Customisation from command line

process.RandomNumberGeneratorService.externalLHEProducer.initialSeed=int(39)
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
