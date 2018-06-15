# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RECOWithDQM -s RAW2DIGI,L1Reco,RECO --datatier RECO --runUnscheduled --nThreads 4 --data --era Run2_2017 --scenario pp --conditions 94X_dataRun2_PromptLike_v9 --eventcontent RECO --filein file:/afs/cern.ch/work/m/mverzett/public/RAW-RECO_ZElectron-94X_dataRun2_PromptLike_v5_RelVal_doubEG2017B-v1.root --no_exec
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('outname', 'track_features.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('ichunk', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    ""
)
options.register('nchunks', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    ""
)
options.register('data', 'RAWTest',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    ""
)
options.register('globalTag', '94X_dataRun2_ReReco_EOY17_v2',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    ""
)
options.register('isMC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    ""
)
options.setDefault('maxEvents', -1)
options.parseArguments()

from CMGTools.LowPtElectrons.samples import all_samples
#split into even chunks
input_files = all_samples[options.data]
n = len(input_files)/options.nchunks
chunks = [input_files[i:i + n] for i in xrange(0, len(input_files), n)]
leftovers = sum(chunks[options.nchunks:], [])
chunks = chunks[:options.nchunks]
for i in range(len(leftovers)):
   chunks[i].append(leftovers[i])


import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if options.isMC:
   process.load('SimGeneral.MixingModule.mixNoPU_cfi')
   process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
   process.load('Configuration.StandardSequences.MagneticField_cff')
   process.load('Configuration.StandardSequences.DigiDMPreMix_cff')
   process.load('SimGeneral.MixingModule.digi_MixPreMix_cfi')
   process.load('Configuration.StandardSequences.DataMixerPreMix_cff')
   process.load('Configuration.StandardSequences.SimL1EmulatorDM_cff')
   process.load('Configuration.StandardSequences.DigiToRawDM_cff')
   process.load('Configuration.StandardSequences.MagneticField_cff')
   process.load('Configuration.StandardSequences.RawToDigi_cff')
   process.load('Configuration.StandardSequences.Reconstruction_cff')
   process.load('Configuration.StandardSequences.RecoSim_cff')
   process.load('CommonTools.ParticleFlow.EITopPAG_cff')
else:
   process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
   process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
   process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      chunks[options.ichunk]
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

process.electronFeatures = cms.Sequence()
if options.isMC:
   process.load('SimTracker/TrackAssociation/trackingParticleRecoTrackAsssociation_cfi')
   #process.reconstruction *= process.simSiPixelDigis
   process.electronFeatures *= process.tpClusterProducer
   process.electronFeatures *= process.quickTrackAssociatorByHits
   process.quickTrackAssociatorByHits.useClusterTPAssociation = False
   process.electronFeatures *= process.trackingParticleRecoTrackAsssociation
   process.genElectronTracks = cms.EDProducer(
      'TracksFromGenParticles',
      tracks = cms.InputTag('generalTracks'),
      association = cms.InputTag('trackingParticleRecoTrackAsssociation')
      )
   process.electronFeatures *= process.genElectronTracks
   from SimGeneral.DataMixingModule.customiseForPremixingInput import customiseForPreMixingInput
   customiseForPreMixingInput(process)
else:
   process.tracksFromConversions = cms.EDProducer(
      'TracksFromConversions',
      src = cms.InputTag('allConversions'),
      beamspot = cms.InputTag('offlineBeamSpot'),
      tracks = cms.InputTag('generalTracks'),
      generalTracksOnly    = cms.bool(True),
      arbitratedMerged     = cms.bool(False),
      arbitratedEcalSeeded = cms.bool(False),
      ecalalgotracks       = cms.bool(False),
      highPurity           = cms.bool(True),
      minProb = cms.double(-1.),
      maxHitsBeforeVtx = cms.uint32(1),
      minLxy = cms.double(-9999.9),
      minVtxR = cms.double(1.5),
      maxVtxR = cms.double(3.5),
      minLeadPt = cms.double(0.),
      )
   process.electronFeatures *= process.tracksFromConversions
   
   process.looseTracksFromConversions = cms.EDProducer(
      'TracksFromConversions',
      src = cms.InputTag('allConversions'),
      beamspot = cms.InputTag('offlineBeamSpot'),
      tracks = cms.InputTag('generalTracks'),
      generalTracksOnly    = cms.bool(True),
      arbitratedMerged     = cms.bool(False),
      arbitratedEcalSeeded = cms.bool(False),
      ecalalgotracks       = cms.bool(False),
      highPurity           = cms.bool(True),
      minProb = cms.double(-1.),
      maxHitsBeforeVtx = cms.uint32(999),
      minLxy = cms.double(-9999.9),
      minVtxR = cms.double(0.),
      maxVtxR = cms.double(999.),
      minLeadPt = cms.double(0.),
      )
   process.electronFeatures *= process.looseTracksFromConversions


process.ntuples = cms.EDAnalyzer(
   'TrackerDrivenSeedFeatures',
   filename = cms.string(options.outname.replace('.root', '_signal.root')),
   beamspot = cms.InputTag('offlineBeamSpot'),
   prescale = cms.double(1.)
   )
for name, val in process.trackerDrivenElectronSeeds.parameters_().iteritems():
   setattr(process.ntuples, name, val)
process.ntuples.tracks = cms.InputTag('genElectronTracks') \
   if options.isMC else \
   cms.InputTag('tracksFromConversions', 'electrons')
process.ntuples.MaxPt = cms.double(10.0)
process.ntuples.MinPt = cms.double(1.0)
process.electronFeatures *= process.ntuples

if not options.isMC:
   process.ntuplesBackground = process.ntuples.clone(
      filename = cms.string(options.outname.replace('.root', '_background.root')),
      tracks = cms.InputTag('looseTracksFromConversions', 'NOTelectrons'),
      prescale = cms.double(0.05),
      )
   process.electronFeatures *= process.ntuplesBackground

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, '')
#'94X_dataRun2_PromptLike_v9', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
if options.isMC:
   #process.electronsFeatsMC = cms.Path(process.electronFeatures)
   ## for name in process.electronFeatures.moduleNames():
   ##    process.recosim.add(getattr(process, name))
   process.recosim_step = cms.Path(process.recosim)
   process.reconstruction_step *= process.electronFeatures
   process.eventinterpretaion_step = cms.Path(process.EIsequence)
   process.schedule = cms.Schedule(
      process.raw2digi_step, process.reconstruction_step, 
      process.recosim_step, process.eventinterpretaion_step
      )
else:
   process.reconstruction *= process.electronFeatures
   process.L1Reco_step = cms.Path(process.L1Reco)
   process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step)

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
      wantSummary = cms.untracked.bool(False),
      #SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

open('pydump.py','w').write(process.dumpPython())
