from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer  import * 
from CMGTools.BKstLL.analyzers.treeVariables import *
from pdb import set_trace

class BKstLLHLTInfoTreeProducer(AutoFillTreeProducer):
	def __init__(self, cfg_ana, cfg_comp, looperName):
		super(BKstLLHLTInfoTreeProducer, self).__init__(cfg_ana, cfg_comp, looperName)
		self.globalVariables += [
			NTupleVariable('l1_mu_pt' , lambda ev: ev.l1mu.pt()),
			NTupleVariable('l1_mu_eta', lambda ev: ev.l1mu.eta()),
			NTupleVariable('l1_mu_phi', lambda ev: ev.l1mu.phi()),
			NTupleVariable('mu_pt' , lambda ev: ev.tag_mu.pt()),
			NTupleVariable('mu_eta', lambda ev: ev.tag_mu.eta()),
			NTupleVariable('mu_phi', lambda ev: ev.tag_mu.phi()),
			NTupleVariable('mu_iso', lambda ev: ev.tag_mu.absIsoWithFSR()),
			#NTupleVariable('mu_sip2D', lambda ev: ev.tag_mu.absIsoWithFSR()),
			#NTupleVariable('mu_sip3D', lambda ev: ev.tag_mu.absIsoWithFSR()),			
			NTupleVariable('jet_pt' , lambda ev: ev.probe_jet.pt()  if ev.probe_jet is not None else -1),
			NTupleVariable('jet_eta', lambda ev: ev.probe_jet.eta() if ev.probe_jet is not None else -1),
			NTupleVariable('jet_phi', lambda ev: ev.probe_jet.phi() if ev.probe_jet is not None else -1),
			NTupleVariable('jet_m'  , lambda ev: ev.probe_jet.mass()if ev.probe_jet is not None else -1),
			NTupleVariable('jet_deepcsv', lambda ev: \
											 ev.probe_jet.bDiscriminator('pfDeepCSVJetTags:probb') + \
											 ev.probe_jet.bDiscriminator('pfDeepCSVJetTags:probbb') \
											 if ev.probe_jet is not None else -1
										 )


			]
		
		self.initDone = True

## class BKstLLHLTInfoTreeProducer(TreeAnalyzerNumpy):
## 
## 	def __init__(self, *args):
## 		super(BKstLLHLTInfoTreeProducer, self).__init__(*args)
## 		self.skimFunction = 'True'
## 		if hasattr(self.cfg_ana, 'skimFunction'):
## 			self.skimFunction = self.cfg_ana.skimFunction
## 
## 	def declareVariables(self,setup):
## 		#branches definitions
## 		self.tree.var('l1_mu_pt' , float)
## 		self.tree.var('l1_mu_eta', float)
## 		self.tree.var('l1_mu_phi', float)
## 		self.tree.var('mu_pt' , float)
## 		self.tree.var('mu_eta', float)
## 		self.tree.var('mu_phi', float)
## 		self.tree.var('mu_iso', float)
## 		self.tree.var('jet_pt', float)
## 		self.tree.var('jet_eta', float)
## 		self.tree.var('jet_phi', float)
## 		self.tree.var('jet_m', float)
## 		self.tree.var('jet_deepcsv', float)
## 
## 	def process(self, event):
## 		self.readCollections(event.input)
## 		self.tree.reset()
## 		
## 		if not eval(self.skimFunction):
## 			return False
## 
## 		self.tree.fill('l1_mu_pt' , event.l1mu.pt())
## 		self.tree.fill('l1_mu_eta', event.l1mu.eta())
## 		self.tree.fill('l1_mu_phi', event.l1mu.phi())
## 
## 		self.tree.fill('mu_pt' , event.tag_mu.pt())
## 		self.tree.fill('mu_eta', event.tag_mu.eta())
## 		self.tree.fill('mu_phi', event.tag_mu.phi())
## 		self.tree.fill('mu_iso', event.tag_mu.absIsoWithFSR())
## 		
## 		if event.probe_jet is not None:
## 			self.tree.fill('jet_pt' , event.probe_jet.pt())
## 			self.tree.fill('jet_eta', event.probe_jet.eta())
## 			self.tree.fill('jet_phi', event.probe_jet.phi())
## 			self.tree.fill('jet_m'  , event.probe_jet.mass())
## 			self.tree.fill(
## 				'jet_deepcsv', 
## 				event.probe_jet.bDiscriminator('pfDeepCSVJetTags:probb') +
## 				event.probe_jet.bDiscriminator('pfDeepCSVJetTags:probbb') 
## 			)
## 		else:
## 			self.tree.fill('jet_pt' , -1)
## 			self.tree.fill('jet_eta', -1)
## 			self.tree.fill('jet_phi', -1)
## 			self.tree.fill('jet_m'  , -1)
## 			self.tree.fill('jet_deepcsv', -1)
## 
## 		if eval(self.skimFunction):
## 			self.tree.tree.Fill()
			
