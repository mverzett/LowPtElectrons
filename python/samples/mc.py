import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator

creator = ComponentCreator()

BtoKee  = creator.makeMyPrivateMCComponent(
   "BtoKee",
   "/BToKee_Pythia/tstreble-BToKee_Pythia_MINIAODSIM_18_03_21-393c5b9eb90ffa47a3c9a0f5562ec979/USER",
   "PRIVATE", ".*root", 'phys03', 1, useAAA=True
   )

BtoKee_ext= creator.makeMyPrivateMCComponent(
   "BtoKee",
   "/BToKee_Pythia/tstreble-BToKee_Pythia_AODSIM_18_03_19-69dd6d102bc2e11aebadb61c0acf0068/USER",
   "PRIVATE", ".*root", 'phys03', 1, useAAA=True
   )

BtoKee_test = creator.makeMCComponentFromLocal(
    "BtoKee_test",
    "XXX",
    path = "/eos/user/m/mverzett/low_pt_electrons/", 
    pattern=".*BtoKee_MC_AOD.root",
    xSec=1.0,
)
