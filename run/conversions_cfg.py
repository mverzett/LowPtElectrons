import FWCore.ParameterSet.Config as cms
import math

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('outname', 'lowpt_ntuples.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('ichunk', 9,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    ""
)
options.register('nchunks', 10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    ""
)
options.setDefault('maxEvents', -1)
options.parseArguments()

from CMGTools.LowPtElectrons.samples.fill6371 import input_files
#split into even chunks
n = len(input_files)/options.nchunks
chunks = [input_files[i:i + n] for i in xrange(0, len(input_files), n)]
leftovers = sum(chunks[options.nchunks:], [])
chunks = chunks[:options.nchunks]
for i in range(len(leftovers)):
   chunks[i].append(leftovers[i])


process = cms.Process("SoftElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'

process.TFileService = cms.Service(
   "TFileService", 
   fileName = cms.string(options.outname)
   )

process.options   = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(True),
      SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

#put input file
process.source = cms.Source(
   "PoolSource",
   fileNames = cms.untracked.vstring(
      #chunks[options.ichunk]
      'file:/afs/cern.ch/work/m/mverzett/public/AOD_300Evts.root'
      ),
   )

## #APPLY JSON
## import PhysicsTools.PythonAnalysis.LumiList as LumiList
## myLumis = LumiList.LumiList(
##    filename = '/afs/cern.ch/work/m/mverzett/RK94New/src/CMGTools/LowPtElectrons/run/fill6371_JSON.txt'
##    )
## process.source.lumisToProcess = myLumis.getVLuminosityBlockRange()


from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(
   process,applyEnergyCorrections=False,
   applyVIDOnCorrectedEgamma=False,
   isMiniAOD=False
)

process.ntuples = cms.EDAnalyzer(
   'ConversionsNtuples',
   tracks = cms.InputTag('generalTracks'),
   electrons = cms.InputTag('gedGsfElectrons'),
   electronTracks = cms.InputTag('electronGsfTracks'),
   eid = cms.InputTag('egmGsfElectronIDs', 'mvaEleID-Fall17-noIso-V1-wp90'),
   trigger = cms.InputTag('TriggerResults', '','HLT'),
   triggerEvent = cms.InputTag('hltTriggerSummaryAOD'),
   #CUTS
   trkMinPt = cms.double(1.),
   trkMaxPt = cms.double(20.),
   trkMaxDTheta = cms.double(math.pi/6.),
   vtxMaxChi2 = cms.double(10.),
   vtxMaxMass = cms.double(10.),
)

process.seq = cms.Sequence(
   process.egammaPostRecoSeq
   * process.ntuples
)

## process.tsk = cms.Task()
## #Trick to make it work in 9_1_X
## for mod in process.producers_().itervalues():
##     process.tsk.add(mod)
## for mod in process.filters_().itervalues():
##     process.tsk.add(mod)

process.p = cms.Path(
   process.seq
   #process.tsk
)

## process.edmOut = cms.OutputModule(
##    "PoolOutputModule",
##    # use this in case of filter available
##    outputCommands = cms.untracked.vstring( 
##       'keep *',
##       ),
##    fileName = cms.untracked.string('edmTEST.root')
##    )
## 
## process.end = cms.EndPath(
##    process.edmOut
##    )

