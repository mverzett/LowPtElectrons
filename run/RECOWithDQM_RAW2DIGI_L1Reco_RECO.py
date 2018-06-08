# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RECOWithDQM -s RAW2DIGI,L1Reco,RECO --datatier RECO --runUnscheduled --nThreads 4 --data --era Run2_2017 --scenario pp --conditions 94X_dataRun2_PromptLike_v9 --eventcontent RECO --filein file:/afs/cern.ch/work/m/mverzett/public/RAW-RECO_ZElectron-94X_dataRun2_PromptLike_v5_RelVal_doubEG2017B-v1.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'file:/afs/cern.ch/work/m/mverzett/public/RAW-RECO_ZElectron-94X_dataRun2_PromptLike_v5_RelVal_doubEG2017B-v1.root'
      #'/store/data/Run2017F/SingleElectron/RAW/v1/000/306/459/00000/B6410DA6-C8C5-E711-947F-02163E01A45E.root'
      ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RECOWithDQM nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

## Output definition
## process.RECOoutput = cms.OutputModule("PoolOutputModule",
##     dataset = cms.untracked.PSet(
##         dataTier = cms.untracked.string('RECO'),
##         filterName = cms.untracked.string('')
##     ),
##     fileName = cms.untracked.string('RECOWithDQM_RAW2DIGI_L1Reco_RECO.root'),
##     outputCommands = process.RECOEventContent.outputCommands,
##     splitLevel = cms.untracked.int32(0)
## )

process.ntuples = cms.EDAnalyzer(
   'TrackerDrivenSeedFeatures',
   filename = cms.string('ntuples.root'),# options.outname),
   )
for name, val in process.trackerDrivenElectronSeeds.parameters_().iteritems():
   setattr(process.ntuples, name, val)
process.ntuples.MaxPt = cms.double(10.0)
process.ntuples.MinPt = cms.double(1.0)
process.ntuples.tracks = process.ntuples.TkColList[0]
process.reconstruction *= process.ntuples

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v2', '')
#'94X_dataRun2_PromptLike_v9', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.ntuples_step = cms.Path(process.ntuples_seq)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(4)
#process.options.numberOfStreams=cms.untracked.uint32(0)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options   = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(False), #True),
      #SkipEvent = cms.untracked.vstring('ProductNotFound'),
)
