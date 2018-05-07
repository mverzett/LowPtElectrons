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

class ConversionsAnalyzer(Analyzer):
   '''
   '''
   def declareHandles(self):
      super(ConversionsAnalyzer, self).declareHandles()
      self.handles['generalTracks'] = AutoHandle(
          ('generalTracks', ''),
	  'vector<reco::Track>'
	  )

      self.handles['gedGsfElectrons'] = AutoHandle(
          ('gedGsfElectrons', ''),
	  'vector<reco::GsfElectron>'
	  )
      
   def beginLoop(self, setup):
      super(ConversionsAnalyzer, self).beginLoop(setup)
      self.counters.addCounter('ConversionsAnalyzer')
      count = self.counters.counter('ConversionsAnalyzer')
      count.register('all events')
      count.register('passing trigger')
      count.register('passing electron selection')

   def process(self, event):
      self.readCollections(event.input)
      counter = self.counters.counter('ConversionsAnalyzer')
      counter.inc('all events')
      ele_pt = None
      if event.HLT_HLT_Ele27:
         ele_pt = 27
      elif event.HLT_HLT_Ele32:
         ele_pt = 32
      else:
         return False
      counter.inc('passing trigger')
      
      electrons = self.handles['gedGsfElectrons'].product()
      seeds_ele = [i for i in electrons if i.pt() > ele_pt]
      set_trace()
      

    
    
    
    
    
    
