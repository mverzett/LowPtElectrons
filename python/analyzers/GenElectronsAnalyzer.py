import ROOT
from itertools import product, combinations
import math

from PhysicsTools.Heppy.analyzers.core.Analyzer      import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle    import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar             import deltaR, deltaR2, bestMatch
from PhysicsTools.Heppy.physicsobjects.Muon          import Muon
from PhysicsTools.Heppy.physicsobjects.Electron      import Electron
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject

from pdb import set_trace

class GenElectronsAnalyzer(Analyzer):
   '''
   '''
   ##def __init__(self, *args, **kwargs):
   ##   super(GenElectronsAnalyzer, self).__init__(*args, **kwargs)

   def declareHandles(self):
      super(GenElectronsAnalyzer, self).declareHandles()
      self.handles['generalTracks'] = AutoHandle(
          ('generalTracks', ''),
	  'vector<reco::Track>'
	  )

      self.mchandles['genParticles'] = AutoHandle(
          ('genParticles', ''),
	  'vector<reco::GenParticle>'
	  )

      self.handles['gedGsfElectrons'] = AutoHandle(
          ('gedGsfElectrons', ''),
	  'vector<reco::GsfElectron>'
	  )
      
   def beginLoop(self, setup):
      super(GenElectronsAnalyzer, self).beginLoop(setup)
      self.counters.addCounter('GenElectronsAnalyzer')
      count = self.counters.counter('GenElectronsAnalyzer')
      count.register('all events')

   def process(self, event):
      self.readCollections(event.input)
      counter = self.counters.counter('GenElectronsAnalyzer')
      counter.inc('all events')
      set_trace()
      gen_ps = self.mchandles['genParticles'].product()
      electrons = [
	      i for i in gen_ps 
	      if abs(i.pdgId()) == 11
	      if abs(i.mother().pdgId()) == 521]
      for ele in electrons:
         pass
	      
    
   @staticmethod
   def match_tracks(obj, collection):
    '''Return the best match to object in matchCollection, which is the closest object in delta R'''
    candidates = [
	    (i, deltaR(obj, i)) for i in collection 
	    if abs(i.dxy(obj.vertex())) < 0.5
	    if abs(i.dz(obj.vertex())) < 0.5
	    if deltaR(obj, i) < 0.1
	    ]
    if not candidates: 
	    return None, float('+inf')

    best, dr = min(candidates, key=lambda x: x[1])
    return best, dr

    
    
    
    
    
    
