import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator

creator = ComponentCreator()

fill6371 = None
fill6371_test = cfg.DataComponent(
   name = 'fill6371_test',
   files = ['/eos/user/m/mverzett/low_pt_electrons/DATA_AOD.root'],
   intLumi=1,
   triggers = [],
   )
